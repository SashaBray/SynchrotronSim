# -*- coding: utf-8 -*-
"""EngineThread — тело режима real-time (TbtPlant + Scheduler) в отдельном потоке, за швом Link.

Движок крутит модель с ПОЛНОЙ скоростью и владеет доской; UI в главном потоке только читает снимок.
Это перенос тела `control/mode5_realtime.py` (пошаговый цикл, чтение ручек на лету) со сменой приёмника:
вместо записи HTML-дашборда — публикация СНИМКА подписанных адресов + меты через Link.

Базовый такт `sch.step()` = `stride` оборотов (plant.TREV·stride = dt); пучок трекается КАЖДЫЙ такт. FOFB
крутится СВОЕЙ частотой (rate_hz в планировщике, задаётся по модельному времени — развязано от публикации).
Снимок на доски/графики публикуется раз в `publish_every` тактов (≤ управляющего цикла) — меньше значение →
СВЕЖЕЕ данные чаще (между тактами орбита реально эволюционирует, не дубли). Между тактами движок осушает
очередь команд Link и применяет ручки ШТАТНЫМ путём (команда → доска → plant/Limits); ядро (диффур/карты/токи)
не трогается. Модельное время = Σdt (state.t); скорость = модельное/стенное время на шаги. Публикуется и на паузе.

Движок развязан от конструкции plant'а: принимает готовый Scheduler. `build_session` собирает реальную
связку режима 4 на TbtPlant (ленивый импорт ring_cpp), а лёгкие тесты шва подают фейковый plant.
"""
from __future__ import annotations

import threading
import time
from time import perf_counter
from typing import Dict, Optional

import numpy as np

from control import StateView

from . import addresses as A
from .link import (Link, Snapshot, RUN, STOP, RESTART, SHUTDOWN,
                   COMMAND, CONTROL, SET_SCREEN, WATCH, UNWATCH, RECORD)

META_NAMES = ("model_time", "wall_time", "speed", "cycle", "seq", "running", "turn", "alive")


class EngineThread(threading.Thread):
    """Поток-движок над готовым Scheduler. start() запускает цикл; control('run'/'stop') им управляет;
    stop() (или control('shutdown')) гасит поток. Публикует снимок ТОЛЬКО подписанных адресов + мету."""

    def __init__(self, link: Link, scheduler, steps_per_cycle: int = 9, publish_every: int = 3,
                 idle_sleep: float = 0.003, screen_period: float = 1.0, name: str = "engine",
                 run_root: Optional[str] = None):
        super().__init__(name=name, daemon=True)
        self.link = link
        self.sch = scheduler
        self.state = scheduler.state
        self.plant = scheduler.plant
        self.steps_per_cycle = int(steps_per_cycle)   # длина управляющего цикла (FOFB), тактов — СПРАВОЧНО
        self.publish_every = max(1, int(publish_every))  # публиковать снимок раз в N тактов (свежесть досок)
        self._idle = float(idle_sleep)
        self._running = False
        self._stop = threading.Event()
        self._cycle = 0
        self._step = 0                                # счётчик базовых тактов (для кадра раз в publish_every)
        self._seq = 0
        self._cav_cache = None                        # (ключ такта и ВЧ, состояние резонатора) — один расчёт на кадр
        self._wall_run = 0.0                 # накопленные стенные секунды шагов — знаменатель скорости
        self._screens: Dict[str, dict] = {}  # addr → {"area","res"} — настройки экрана (set_screen)
        self._screen_val: Dict[str, dict] = {}   # addr → последний посчитанный результат экрана
        self._screen_last: Dict[str, float] = {}  # addr → стенное время последнего пересчёта
        self._screen_key: Dict[str, tuple] = {}   # addr → (настройки окна, зазор источника) последнего пересчёта
        self._screen_job: Dict[str, object] = {}  # addr → Future расчёта в фоновом потоке (≤1 на экран)
        self._screen_pool = None                  # ThreadPoolExecutor(1) — лениво, при первом экране
        self.screen_period = float(screen_period)  # период пересчёта экрана, с (медленный такт)
        self._unknown: Optional[str] = None  # последняя нераспознанная команда (диагностика)
        self._err: Optional[str] = None      # последнее исключение тела потока
        self._design: Optional[dict] = None  # правки проекта, применённые в модели: отпечаток машины, число, ошибка
        self._kcur: Dict[int, list] = {}     # магнит → силы ядра (k1, k1v, k2, k3) после горячих правок
        # исполнитель опыта (библиотека): шаг, рестарт прохода, лимит оборотов, запись прогона (Watch → Run)
        from synchrotronsim.simulation import Runner
        self.run_root = run_root             # корень проекта/рабочей папки для workspace/runs (None — проект)
        self.runner = Runner(scheduler, run_root=run_root, kind="cockpit", source="cockpit")

    # ---- управление из главного потока ----
    def stop(self) -> None:
        self._stop.set()
        pool, self._screen_pool = self._screen_pool, None
        if pool is not None:                              # фоновый расчёт экрана: не ждать, новых не брать
            pool.shutdown(wait=False, cancel_futures=True)

    @property
    def running(self) -> bool:
        return self._running

    # ---- тело потока ----
    def run(self) -> None:
        try:
            while not self._stop.is_set():
                self._drain()
                if self._stop.is_set():
                    break
                if not self._running:
                    self._publish()                       # кадр и на паузе: UI живой
                    time.sleep(self._idle)
                    continue
                t0 = perf_counter()
                self.runner.step()                        # ОДИН базовый такт = stride оборотов
                self._wall_run += perf_counter() - t0
                self._step += 1
                if self._step % self.publish_every == 0:  # снимок раз в publish_every тактов (не каждый — дешевле)
                    self._cycle += 1
                    self._publish()
                    self.runner.record_tick(self._cycle, {"wall_time": float(self._wall_run)})   # запись — темпом снимка
                if self.runner.limit_reached():
                    self.runner.end_pass("done")          # проход дошёл до лимита
                    if self.runner.repeat:
                        self._restart()                   # новый проход с нуля (и новый прогон, если пишем)
                    else:
                        self._running = False
        except Exception as e:                            # поток не должен падать молча
            self._err = f"{type(e).__name__}: {e}"
            try:
                self._publish()
            except Exception:
                pass
        finally:
            self.runner.end_pass("stopped")

    # ---- очередь команд ----
    def _drain(self) -> None:
        for m in self.link.pull():                        # subscribe/unsubscribe уже применены в pull()
            if m.kind == CONTROL:
                self._control(m.value)
            elif m.kind == COMMAND:
                self._command(m.addr, m.value)
            elif m.kind == SET_SCREEN:
                self._screens[m.addr] = dict(m.value or {})
            elif m.kind == WATCH:
                if not self.runner.watch(m.addr):
                    self._mark_unknown(m.addr)
            elif m.kind == UNWATCH:
                self.runner.unwatch(m.addr)
            elif m.kind == RECORD:
                v = dict(m.value or {})
                self.runner.set_recording(bool(v.get("on")), v.get("max_turns"), bool(v.get("repeat")))

    def _control(self, action: str) -> None:
        if action == RUN:
            self._running = True
        elif action == STOP:
            self._running = False
        elif action == RESTART:
            self.runner.end_pass("stopped")               # прерванный рестартом проход
            self._restart()
        elif action == SHUTDOWN:
            self._stop.set()

    def _restart(self) -> None:
        """Новый проход с нуля (Runner.reset_model) + счётчики публикации; режим пуск/пауза сохраняется."""
        was = self._running
        self._running = False
        self._cycle = 0
        self._step = 0
        self._wall_run = 0.0
        self.runner.reset_model()
        self._running = was

    def _command(self, addr: str, value) -> None:
        """Применить ручку ШТАТНЫМ путём (команда → доска). Диспетчер по типу адреса. Неизвестное —
        записать в диагностику, не падать."""
        addr = A.canonical(addr)                          # прежние имена (drive.*, actual.*) — псевдонимы
        kind = A.kind_of(addr)
        try:
            if kind == "undulator":
                gid, _, attr = addr[len("undulator."):].partition(".")
                u = self.state.undulators.get(gid)
                if u is None:
                    return self._mark_unknown(addr)
                if attr in ("gap", "gap_target"):
                    u.command(float(value))               # уставка привода; реализует зазор plant.advance
                elif attr == "slew":
                    u.slew_mm_s = float(value)
                else:
                    self._mark_unknown(addr)
            elif kind == "loop":
                name, _, fld = addr[len("loop."):].rpartition(".")
                spec = self.sch.loops.get(name)
                if spec is None:
                    return self._mark_unknown(addr)
                if fld == "enabled":
                    self.sch.set_enabled(name, bool(value))
                elif fld == "gain" and "gain" in (getattr(spec.controller, "params", None) or {}):
                    spec.controller.params["gain"] = float(value)   # закон контура читает params["gain"]
                else:
                    self._mark_unknown(addr)
            elif kind == "design":                        # правка ПРОЕКТА (этап 4): горячая — ручкой ядра
                self._design_edit(value)
            elif kind == "cavity":                        # ручки резонатора: cavity.volt_mv | phase_deg | df_hz
                obj = addr.split(".", 1)[1]
                tp = self._design_plant()
                if obj in A.CAVITY_WRITABLE and hasattr(tp, "set_cavity"):
                    tp.set_cavity(**{obj: float(value)})
                else:
                    self._mark_unknown(addr)
            elif kind == "act":                           # привод актуатора: act.<aid>.<поле>
                p = A.parse(addr)
                aid, fld = p.obj, p.attr
                if fld == "slew":
                    self.state.limits.set(aid, slew=float(value))
                elif fld == "lag_tau":
                    self.state.limits.set(aid, lag_tau=float(value))
                elif fld == "res":
                    self.state.limits.set(aid, res=float(value))
                elif fld == "act":
                    self.state.limits.set(aid, act=bool(value))
                else:
                    self._mark_unknown(addr)
            else:
                self._mark_unknown(addr)
        except (KeyError, IndexError, ValueError, AttributeError, TypeError):
            self._mark_unknown(addr)

    def _mark_unknown(self, addr: str) -> None:
        self._unknown = addr

    # ---- правки проекта (synchrotronsim.edits)
    def _design_plant(self):
        """Plant с кольцом орбиты (у CompositePlant — его TbtPlant)."""
        return getattr(self.plant, "tbt", None) or self.plant

    def _design_edit(self, value) -> None:
        """Горячая правка проекта: машина модели пересчитывается (тем же разбором, что из таблиц), сила
        реализуется ручкой ядра базовым возмущением (переживает шаги и рестарт). Холодная — только пересборкой.
        Ошибка — в мете design.error, поток не падает."""
        from synchrotronsim.edits import Edit, apply_machine, hot_commands
        try:
            e = Edit.from_dict(value)
            if e.cost != "hot":
                raise ValueError(f"{e.addr}: правка применяется пересборкой модели, не на ходу")
            tp = self._design_plant()
            m_old = getattr(tp, "machine", None) or getattr(self.sch, "machine", None)
            if m_old is None:
                raise ValueError("у модели нет документа машины")
            m_new = apply_machine(m_old, e)
            ring = tp.ring
            cmds, cur = hot_commands(m_old, m_new, e.element, self._kcur.get(int(e.element)))
            for c, j, v in cmds:
                ring.add_base(c, j, v)
            self._kcur[int(e.element)] = cur               # силы ядра после правки (те же операции, что в C++)
            if hasattr(tp, "_last_kick"):
                tp._last_kick = None                       # следующий шаг — RESET + базовые: ровно как пересборка
            tp.machine = m_new
            if hasattr(ring, "_machine"):
                ring._machine = m_new
            self.sch.machine = m_new                       # паспорт прогона — машина, которая считается
            self.runner.design_edits.append(e.to_dict())
            self._design = {"fingerprint": m_new.fingerprint, "edits": len(self.runner.design_edits), "error": None}
        except Exception as ex:                            # noqa: BLE001 — в мету, не ронять поток
            self._design = dict(self._design or {}, error=f"{type(ex).__name__}: {ex}")

    # ---- ВЧ-резонатор (synchrotronsim.longitudinal): параметры, яма, пучок на резонаторе
    def _cavity_bucket(self):
        tp = self._design_plant()
        ring = getattr(tp, "ring", None)
        rf = getattr(ring, "rf", None)
        if not rf or not hasattr(tp, "U0_mev"):
            return None
        from synchrotronsim.longitudinal import bucket
        opt = getattr(self.state, "optics", None)
        rf0 = getattr(tp, "_rf0", None) or rf               # заданные V, φ (ядру уходят шагом позже)
        u0 = float(tp.u0_now()) if hasattr(tp, "u0_now") else float(tp.U0_mev)   # U₀ текущей машины (зазоры ИД…)
        b = bucket(float(rf0["volt_mv"]), u0, float(ring.design_mev), float(rf0["harmonic"]),
                   float(getattr(opt, "alphac", 0.0) or 0.0), float(rf0["phase_deg"]), float(ring.rf_frequency()),
                   float(tp.F0))
        b["manual_df_hz"] = float(getattr(tp, "rf_manual_df", 0.0))
        b["loop_df_hz"] = float(getattr(tp, "rf_df", 0.0)) - b["manual_df_hz"]
        return b

    def _cavity_value(self, addr: str):
        from synchrotronsim.longitudinal import particle_phase_deg, inside
        obj = addr.split(".", 1)[1]
        tp = self._design_plant()
        try:
            b = self._cavity_bucket()
            if b is None:
                return None
            if obj == "volt_mv":
                return b["V_mv"]
            if obj == "phase_deg":
                return b["phase_deg"]
            if obj == "df_hz":
                return b["manual_df_hz"]
            if obj == "freq_hz":
                return b["freq_hz"]
            L = getattr(tp, "long_tbt", None)
            ph = d = None
            if L is not None and len(L):                   # фаза — по ВЧ, с которым такт трекался
                p0, hh, crf = getattr(tp, "long_rf", None) or (b["phase_deg"], b["h"], b["c_rf"])
                ph = particle_phase_deg(L[:, :, 0], p0, hh, crf)
                d = L[:, :, 1]
            if obj == "bunch":
                if ph is None:
                    return None
                return {"phase_deg": np.asarray(ph, float), "delta": np.asarray(d, float),
                        "turn": int(getattr(tp, "turn", 0)), "alive": bool(getattr(tp, "alive", True))}
            ck = (int(getattr(tp, "turn", 0)), id(L), b["V_mv"], b["phase_deg"], b["freq_hz"], b["U0_mev"],
                  b["alphac"], id(getattr(tp, "_rad", None)))
            if self._cav_cache is not None and self._cav_cache[0] == ck:
                out = self._cav_cache[1]
            else:
                out = dict(b)
                if ph is not None:                         # пучок — последний оборот такта, живые частицы
                    p, dd = ph[:, -1], d[:, -1]
                    ok = np.isfinite(p) & np.isfinite(dd)
                    out.update(n_part=int(len(p)), n_lost=int(np.sum(~ok)))
                    if np.any(ok):
                        p, dd = p[ok], dd[ok]
                        out.update(phi_mean=float(np.mean(p)), phi_std=float(np.std(p)),
                                   delta_mean=float(np.mean(dd)), delta_std=float(np.std(dd)),
                                   n_inside=int(np.sum(inside(b, p, dd))))
                    else:
                        out["n_inside"] = 0
                from synchrotronsim.longitudinal import derived
                rad = tp.radiation_now() if hasattr(tp, "radiation_now") else None
                P = getattr(getattr(tp, "_und_sr", None), "P", None)     # ток — допущение undulator_sr.P
                derived(out, float(tp.ring.circumference), rad=rad,
                        current_a=float(P.I) if P is not None and getattr(P, "I", None) else None)
                self._cav_cache = (ck, out)
            if obj == "state":
                return out
            v = out.get(obj)                               # скаляр каталога → ряд во времени (график)
            return float(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else None
        except Exception:                                  # noqa: BLE001 — окно резонатора опционально
            return None

    def _design_value(self, addr: str):
        """Значение `design.<instance_id>.<столбец>` из машины модели (ячейка строки; NaN — пусто → None)."""
        from synchrotronsim import addresses as AD
        try:
            p = AD.parse(addr)
            m = getattr(self._design_plant(), "machine", None)
            I = m.instances
            r = int(np.nonzero(I.instance_id == p.obj)[0][0])
            v = getattr(I, p.attr)[r]
            v = float(v)
            return None if v != v else v
        except Exception:                                  # noqa: BLE001
            return None

    # ---- публикация снимка ----
    def _meta(self) -> dict:
        mt = float(self.state.t)
        speed = (mt / self._wall_run) if self._wall_run > 0 else 0.0
        m = {"model_time": mt, "wall_time": float(self._wall_run), "speed": float(speed),
             "cycle": int(self._cycle), "seq": int(self._seq), "running": bool(self._running),
             "turn": int(getattr(self.plant, "turn", 0)), "alive": bool(getattr(self.plant, "alive", True))}
        R = self.runner
        if self._design:
            m["design"] = dict(self._design)
        m["recording"] = R.run is not None
        m["run_id"] = R.run.run_id if R.run is not None else ""
        m["pass"] = int(R.pass_no)
        m["max_turns"] = int(R.max_turns or 0)
        m["watches"] = sorted(R.watches)
        if R.err and not self._err:
            self._err = R.err
        if self._unknown:
            m["unknown_command"] = self._unknown
        if self._err:
            m["error"] = self._err
        self._seq += 1
        return m

    def _publish(self) -> None:
        meta = self._meta()
        view = StateView(self.state)
        subs = self.link.subscriptions()
        values = {}
        for a in subs:
            if A.kind_of(a) == "screen":
                values[a] = self._screen_value(a)
            elif A.kind_of(a) == "diag":
                values[a] = self._diag_value(a)
            elif A.kind_of(a) == "design":
                values[a] = self._design_value(a)
            elif A.kind_of(a) == "cavity":
                values[a] = self._cavity_value(a)
            else:
                values[a] = A.resolve(view, a, loops=self.sch.loops, meta=meta)
        self.link.publish(Snapshot(values=values, meta=meta))

    # ---- стандартные диагностические графики (`diag.<name>`): статичны для сессии, считаются один раз ----
    def _diag_value(self, addr: str) -> Optional[dict]:
        """Данные стандартного графика (client/diagnostics.py). Пересчёт — только при смене ключа кэша
        (напр. rcond FOFB); ошибку отдаём окну словарём, поток не роняем."""
        name = addr.split(".", 1)[1] if "." in addr else ""
        cache = self.__dict__.setdefault("_diag_cache", {})            # addr → (ключ, значение)
        try:
            from . import diagnostics as D
            key = D.cache_key(name, self.sch.loops)
            hit = cache.get(addr)
            if hit is None or hit[0] != key:
                cache[addr] = (key, D.compute(name, self.plant, self.state, self.sch.loops, self.sch.resources))
            return cache[addr][1]
        except Exception as e:                      # noqa: BLE001 — диагностика опциональна
            return {"error": f"{type(e).__name__}: {e}"}

    # ---- экран: тяжёлый пересчёт картины в ФОНОВОМ потоке (ScreenRecorder) ----
    def _screen_value(self, addr: str) -> Optional[dict]:
        """Последний результат экрана; пересчёт — в ФОНОВОМ потоке, модель не ждёт. Картина зависит лишь от
        зазора ондулятора-источника и настроек окна (точка спектра, разрешение): новые настройки — пересчёт сразу;
        движение зазора — не чаще screen_period; ничего не менялось — без пересчёта. Пока считается — отдаём
        прежнюю картину с флагом computing (первый раз — только флаг). Вход расчёта — СНИМОК доски на момент
        запуска (зазоры, photon.x/y), поэтому поток модели и поток экрана не делят изменяемое состояние."""
        job = self._screen_job.get(addr)
        if job is not None and job.done():                 # фоновый расчёт закончился — забрать результат
            try:
                self._screen_val[addr] = job.result()
            except Exception as e:                         # noqa: BLE001 — экран опционален
                self._screen_val[addr] = {"error": f"{type(e).__name__}: {e}"}
            del self._screen_job[addr]; job = None
        cfg = self._screens.get(addr, {}) or {}
        ckey = (int(cfg.get("res") or 121), tuple(float(v) for v in (cfg.get("point") or (0.0, 0.0))),
                _quantity(cfg))
        try:
            from control.recorders import ScreenRecorder
            sid = addr.split(".", 1)[1] if "." in addr else addr
            gap = round(float(ScreenRecorder(params=sid)._live_gap(self.state)), 6)
        except Exception:                                  # noqa: BLE001 — экран опционален
            gap = None
        last = self._screen_key.get(addr)                  # ключ последнего ЗАПУЩЕННОГО расчёта
        now = perf_counter()
        due = (now - self._screen_last.get(addr, -1e9)) >= self.screen_period
        need = last is None or last[0] != ckey or (last[1] != gap and due)
        if need and job is None and not self._stop.is_set():
            if self._screen_pool is None:
                from concurrent.futures import ThreadPoolExecutor
                self._screen_pool = ThreadPoolExecutor(max_workers=1, thread_name_prefix="screen")
            self._screen_last[addr] = now
            self._screen_key[addr] = (ckey, gap)
            self._screen_job[addr] = job = self._screen_pool.submit(
                self._compute_screen, addr, dict(cfg), self._screen_inputs())
        val = self._screen_val.get(addr)
        if job is None:
            return val
        return dict(val, computing=True) if isinstance(val, dict) else {"computing": True}

    def _screen_inputs(self):
        """Снимок доски для расчёта экрана в фоне: только то, что читает ScreenRecorder (живые зазоры
        ондуляторов, положение фотонного пятна) — копии, без ссылок на изменяемое состояние модели."""
        from types import SimpleNamespace
        obs = {k: np.array(self.state.observables[k], float, copy=True)
               for k in ("photon.x", "photon.y") if self.state.observables.get(k) is not None}
        und = {gid: SimpleNamespace(gap_mm=float(u.gap_mm))
               for gid, u in (getattr(self.state, "undulators", None) or {}).items()}
        return SimpleNamespace(observables=obs, undulators=und)

    def _compute_screen(self, addr: str, cfg: Optional[dict] = None, st=None) -> dict:
        sid = addr.split(".", 1)[1] if "." in addr else addr
        cfg = (self._screens.get(addr, {}) or {}) if cfg is None else cfg
        st = self.state if st is None else st
        mesh = int(cfg.get("res") or 121)
        pt = cfg.get("point") or (0.0, 0.0)                 # точка на экране (мм) для спектра; (0,0)=на оси
        bend = self._bend_context(sid)
        if bend is not None:                                # экран у диполя: веера, поток на ε_c, сплошной спектр
            return self._compute_bend(sid, bend, mesh, pt, _quantity(cfg), st)
        if _quantity(cfg) == "power":
            return self._compute_power(sid, mesh, pt, st)
        try:
            from control.recorders import ScreenRecorder
            rec = ScreenRecorder(params=sid, mesh=mesh, point=pt, machine=self._machine())
            rec.sample(st, 0, None)                         # ScreenRecorder читает только доску (снимок)
            rec.finalize(st, None)
            art = rec.artifact()
            if not art:
                return {"error": "экран не дал артефакта"}
            data = art[1]
            heat = np.asarray(data.get("measured", data.get("ideal_end", [])), float)
            out = {"heatmap": heat.tolist() if heat.size else None,
                   "extent": np.asarray(data.get("extent", []), float).tolist(),
                   "e1_series": np.asarray(data.get("gap_series", []), float).tolist(),
                   "mesh": mesh}
            if "gap_end_mm" in data:
                out["gap"] = float(np.ravel(data["gap_end_mm"])[0])
            if "energy_end_keV" in data:
                out["e1"] = float(np.ravel(data["energy_end_keV"])[0])
            sp = data.get("spectrum")                        # настоящий спектр (flux vs энергия) в точке
            if isinstance(sp, dict) and sp.get("E_keV") is not None:
                out["spectrum"] = {"E_keV": np.asarray(sp["E_keV"], float).tolist(),
                                   "flux": np.asarray(sp["flux"], float).tolist(),
                                   "u_mm": float(sp["u_mm"]), "v_mm": float(sp["v_mm"]),
                                   "E1": float(sp["E1_keV"])}
            out["quantity"] = "flux"
            return out
        except Exception as e:                      # экран опционален/тяжёл — ошибку отдаём окну, не падаем
            return {"error": f"{type(e).__name__}: {e}"}

    def _sa(self):
        import os, sys
        und = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "undulators")
        if und not in sys.path:
            sys.path.insert(0, und)
        import screen_analysis as sa
        return sa

    def _lattice(self) -> int:
        return int(getattr(self.plant, "lattice", None) or 26)

    def _machine(self):
        """Документ машины plant'а (экраны и геометрия — из него, таблицы не читаются); нет — None."""
        return getattr(self.plant, "machine", None)

    def _lat(self):
        """Машина для screen_analysis: документ, а без него — номер решётки (прежний путь через таблицы)."""
        m = self._machine()
        return m if m is not None else self._lattice()

    def _bend_context(self, sid: str):
        """Контекст живого экрана, если его источник — ДИПОЛЬ (экраны у линеек SCR_DP…, замороженные mode=bend);
        иначе None (экран ондулятора — прежний путь)."""
        try:
            sa = self._sa()
            if not str(sid).startswith("SCR_DP") and not any(
                    r.get("screen_id") == sid and r.get("lattice_id") == str(self._lattice()) and r.get("mode") == "bend"
                    for r in sa._screen_rows(self._lat())):
                return None
            c = sa.live_context(self._lat(), sid)
            return c if c is not None and c["scr"].get("mode") == "bend" else None
        except Exception:                                  # noqa: BLE001 — нет таблиц: не диполь
            return None

    def _compute_bend(self, sid: str, c: dict, mesh: int, pt, quantity: str, st) -> dict:
        """Экран у диполя — по траектории электрона (точный режим, мягкий край, без стыков вееров): «мощность» —
        Вт/мм², «поток» — на ε_c центрального модуля, ф/с/мрад²/0.1%BW; спектр в точке — сплошной, с локальной ε_c
        участка дуги, смотрящего на точку. Карта статична (зазоры не влияют) — считается один раз на сетку."""
        try:
            sa = self._sa()
            v = sa.bend_view(self._lat(), sid, mesh, quantity)        # по траектории; кэш — считается один раз
            out = {"mesh": mesh, "quantity": quantity, "kind": "bend", "ec": v["ec"], "dev": v["dev"], "D": v["D"],
                   "heatmap": np.asarray(v["S"], float).tolist(), "extent": [float(e) for e in v["extent"]]}
            if quantity == "power":
                out.update(peak=v["peak"], total=v["total"], und_W=0.0, bend_W=v["bend_W"], n_bend=v["n_bend"],
                           centroid=[v["cu"], v["cv"]])
            ec0 = v["ec"]
            sp = sa.bend_spectrum(self._lat(), sid, float(pt[0]), float(pt[1]))
            out["spectrum"] = {"E_keV": np.asarray(sp["E_keV"], float).tolist(),
                               "flux": np.asarray(sp["flux"], float).tolist(),
                               "u_mm": float(pt[0]), "v_mm": float(pt[1]), "E1": ec0}
            return out
        except Exception as e:                      # noqa: BLE001 — экран опционален
            return {"error": f"{type(e).__name__}: {e}"}

    def _compute_power(self, sid: str, mesh: int, pt, st) -> dict:
        """Карта ПЛОТНОСТИ МОЩНОСТИ (Вт/мм², быстрый режим `screen_analysis.live_power_map`) на живых зазорах
        всех ондуляторов + веера соседних диполей; геометрия экрана — как у сцены. Спектр в точке — тот же,
        что в режиме потока (зазор источника)."""
        try:
            import os, sys
            und = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "undulators")
            if und not in sys.path:
                sys.path.insert(0, und)
            import screen_analysis as sa
            import undulator_sr as sr
            from control.recorders import ScreenRecorder
            gaps = {gid: float(u.gap_mm) for gid, u in (getattr(st, "undulators", None) or {}).items()}
            lattice = int(getattr(self.plant, "lattice", None) or 26)
            fp = sa.live_power_map(self._lat(), sid, gaps, mesh)
            if fp is None:
                return {"error": f"карта мощности: экран {sid} не построить (источник — не ондулятор решётки {lattice})"}
            rec = ScreenRecorder(params=sid, mesh=mesh, point=pt, machine=self._machine())
            gap = float(rec._live_gap(st))
            out = {"heatmap": np.asarray(fp["S"], float).tolist(), "extent": [float(e) for e in fp["extent"]],
                   "mesh": mesh, "quantity": "power", "gap": gap,
                   "e1": float(sr.harmonic_energy_keV(1, sr.K_of_B0(sr.B0_of_gap(gap)), 0.0)),
                   "peak": fp["peak"], "total": fp["total"], "und_W": fp["und_W"], "bend_W": fp["bend_W"],
                   "n_bend": fp["n_bend"], "centroid": [fp["cu"], fp["cv"]], "D": fp.get("D")}
            sp = rec.point_spectrum(gap)                     # плечо — как в режиме потока (dist_m экрана)
            if isinstance(sp, dict):
                out["spectrum"] = {"E_keV": np.asarray(sp["E_keV"], float).tolist(),
                                   "flux": np.asarray(sp["flux"], float).tolist(),
                                   "u_mm": float(sp["u_mm"]), "v_mm": float(sp["v_mm"]), "E1": float(sp["E1_keV"])}
            return out
        except Exception as e:                      # экран опционален — ошибку отдаём окну, не падаем
            return {"error": f"{type(e).__name__}: {e}"}


def _quantity(cfg) -> str:
    """Что рисует экран: "power" (мощность, Вт/мм²) или "flux" (поток на E₁, по умолчанию — как раньше)."""
    return "power" if str((cfg or {}).get("quantity") or "").lower() == "power" else "flux"


def build_session(lattice: int = 26, nsl: int = 8, settle_turns: int = 250, fofb: bool = True,
                  disturbance_spec=None, gap_mm: float = 5.0, steps_per_cycle: int = 9,
                  meas_noise: float = 1e-7, n_photon: int = 1, optics: bool = False,
                  optics_drift=(1e-3, 1e-2, 3e-2), optics_order: int = 4, optics_period_ticks: int = 90,
                  optics_nsl: int = 2, rf: bool = False, id_pole_nsl: int = 1, **plant_kw):
    """Собрать (plant, state, scheduler) режима 4. Контуры добавлены ВЫКЛ (оператор включает ручкой/панелью).
    `optics=False` — только орбита на TbtPlant (fofb_x/y). `optics=True` — CompositePlant (орбита + линеариз.
    оптика) и ПЯТЬ контуров: fofb_x/y (≈2 кГц) + tune/chroma/adts (медленный темп `optics_period_ticks`).
    Оптические контуры — самодостаточные плагины (control/controllers/optics_*.py), их схему вычисления
    собираем из `provides()`. Импорт ленивый, чтобы лёгкие тесты шва работали без ring_cpp.
    `rf=True` — + контур ВЧ-частоты `rf` (control/controllers/rf_frequency.py: энергия по дисперсионной орбите →
    частота ВЧ; актуатор `rf.freq` в схеме), темп — как у FOFB (раз в цикл; внутри усреднение ≈ T_s).
    `id_pole_nsl` — срезов у полюсов ондуляторов в трекинге (1: вдвое быстрее, физика та же до 1e-8 по тюну;
    None — как у остальных элементов, nsl)."""
    # с шага 5 плана библиотеки — опыт-объект + ЕДИНЫЙ сборщик (synchrotronsim.simulation.build_session);
    # параметры дословно прежние (Experiment.cockpit), опыт кладётся в планировщик для паспорта прогона
    from synchrotronsim.experiment import Experiment
    from synchrotronsim.simulation import build_session as _build
    exp = Experiment.cockpit(lattice=lattice, nsl=nsl, settle_turns=settle_turns, fofb=fofb,
                             disturbance_spec=disturbance_spec, gap_mm=gap_mm, steps_per_cycle=steps_per_cycle,
                             meas_noise=meas_noise, n_photon=n_photon, optics=optics, optics_drift=optics_drift,
                             optics_order=optics_order, optics_period_ticks=optics_period_ticks,
                             optics_nsl=optics_nsl, rf=rf, id_pole_nsl=id_pole_nsl, **plant_kw)
    S = _build(exp)
    plant, state, sch = S.plant, S.state, S.sch
    sch.experiment = exp
    return plant, state, sch


# ---------------------------------------------------------------------------------------------- описание сеанса
# Всё, что интерфейсу кокпита нужно знать о собранной модели, — словарём (числа/строки): так же для модели в этом
# процессе и для модели на вычислителе (synchrotronsim.live шлёт его по сети). Без Qt — годится серверу.

def loops_info(loops: dict) -> list:
    """Контуры САУ для панели: [{name, file, enabled, tip}]. `file` — имя файла, описывающего логику контура
    (модуль класса контроллера, напр. integral_svd.py); разные контуры одного типа различаются по name."""
    import inspect
    import os
    out = []
    for name, spec in (loops or {}).items():
        ctrl = getattr(spec, "controller", None)
        try:
            file = os.path.basename(inspect.getsourcefile(type(ctrl)) or inspect.getfile(type(ctrl)))
        except (TypeError, OSError):
            file = "?"
        cls = type(ctrl).__name__ if ctrl is not None else "?"
        rate = float(getattr(spec, "rate_hz", 0.0) or 0.0)
        tip = f"{cls}  ·  {rate:.0f} Гц" if rate else cls
        out.append({"name": name, "file": file, "enabled": bool(getattr(spec, "enabled", False)), "tip": tip})
    return out


OPTICS_COMPONENTS = {
    "tune":   [("obs.tune:0", "тюн νx"), ("obs.tune:1", "тюн νy")],
    "chroma": [("obs.chroma:0", "хроматизм ξx"), ("obs.chroma:1", "хроматизм ξy")],
    "adts":   [("obs.adts:0", "adts ∂νx/∂Jx"), ("obs.adts:1", "adts ∂νy/∂Jy"),
               ("obs.adts:2", "adts ∂νx/∂Jy")],
}


def optics_specs(state) -> list:
    """Строки оптики для панели: [(addr, label, fmt)] по КОМПОНЕНТАМ тюна/хроматизма/adts — только для тех
    наблюдаемых, что реально есть в состоянии (optics=True). Пусто при orbit-only (optics=False)."""
    obs = getattr(state, "observables", {}) or {}
    specs = []
    for o in ("tune", "chroma", "adts"):
        if o in obs:
            for addr, label in OPTICS_COMPONENTS[o]:
                specs.append((addr, label, "{:.4f}"))
    return specs


def init_undulators(state, slew, open_gid="und_01", open_mm=5.0, closed_mm=50.0):
    """Стартовая конфигурация ондуляторов: ВКЛЮЧЁН один (`open_gid`, челюсти сведены на `open_mm`), остальные
    ВЫКЛЮЧЕНЫ (челюсти разведены на `closed_mm` — поток ∝ K² падает ~1e-7). Поднимаем `gap_max` (деф. 20),
    иначе `command/advance` клипнут 50 обратно; ставим gap_mm==gap_target, чтобы `advance()` не уводил старт."""
    for gid, u in (getattr(state, "undulators", {}) or {}).items():
        try:
            u.slew_mm_s = float(slew)
            u.gap_max = max(float(getattr(u, "gap_max", 20.0)), float(closed_mm))
            g = float(open_mm) if gid == open_gid else float(closed_mm)
            u.gap_mm = u.gap_target = g
        except Exception:
            pass


def apply_gaps(state, gaps: Optional[dict]) -> None:
    """Перенос зазоров при пересборке: {gid: (зазор, уставка, скорость привода)} — как были до пересборки."""
    for gid, v in (gaps or {}).items():
        u = (getattr(state, "undulators", {}) or {}).get(gid)
        if u is None:
            continue
        g, gt, slew = (float(x) for x in v)
        u.slew_mm_s = slew
        u.gap_max = max(float(getattr(u, "gap_max", 20.0)), g, gt)
        u.gap_mm, u.gap_target = g, gt


def read_screen_ids(machine=None) -> list:
    """Список экранов: у документа машины — её экраны; без машины — таблица screens (через synchrotronsim.db)."""
    if machine is not None:
        return [s.screen_id for s in machine.screens]
    try:
        from synchrotronsim.db import Repository, PROJECT_ROOT
        return [v.strip() for v in Repository(PROJECT_ROOT).table("screens", required=False).column("screen_id")
                if v.strip()]
    except (OSError, ValueError):
        return []


def loop_targets(state, loops) -> list:
    """Строки «контур: текущее → цель» для панели: [{loop, label, cur, tgt, tgt_value, scale, units}]. Текущее —
    наблюдаемая-датчик контура (компоненты тюна/хроматизма/ADTS, СКО орбиты, энергия δ̂), цель — уставка контура
    на доске (`ref.<имя>`) или ноль, если уставки нет (проектная орбита / энергия)."""
    obs = getattr(state, "observables", {}) or {}
    refs = getattr(state, "refs", {}) or {}
    out = []
    for name, spec in (loops or {}).items():
        ctrl = getattr(spec, "controller", None)
        sensor = getattr(ctrl, "sensor", None)
        ref = getattr(ctrl, "ref", None)
        if not sensor or sensor not in obs:
            continue
        has_ref = bool(ref) and ref in refs
        if sensor in OPTICS_COMPONENTS:
            for i, (addr, label) in enumerate(OPTICS_COMPONENTS[sensor]):
                out.append(dict(loop=name, label=label, cur=addr, tgt=f"ref.{ref}:{i}" if has_ref else None,
                                tgt_value=None if has_ref else 0.0, scale=1.0, units="", fmt="{:.4f}"))
        elif sensor in ("ebpm.x", "ebpm.y"):
            pl = sensor[-1]
            out.append(dict(loop=name, label=f"орбита {pl}, СКО", cur=f"obs.{sensor}:rms",
                            tgt=f"ref.{ref}:rms" if has_ref else None, tgt_value=None if has_ref else 0.0,
                            scale=1e6, units="мкм", fmt="{:.3f}"))
        elif sensor == "rf.delta":
            out.append(dict(loop=name, label="энергия δ̂", cur="obs.rf.delta:0", tgt=f"ref.{ref}:0" if has_ref else None,
                            tgt_value=None if has_ref else 0.0, scale=1e6, units="ppm", fmt="{:.2f}"))
        else:
            out.append(dict(loop=name, label=f"{sensor}, СКО", cur=f"obs.{sensor}:rms",
                            tgt=f"ref.{ref}:rms" if has_ref else None, tgt_value=None if has_ref else 0.0,
                            scale=1.0, units="", fmt="{:.6g}"))
    return out


def describe_session(plant, state, sch, exp=None) -> dict:
    """Описание собранного сеанса для интерфейса: опыт, геометрия сцены, ондуляторы, контуры, оптика, наблюдаемые,
    экраны, машина. Только числа и строки — уходит по сети как есть."""
    from .scene3d import build_geometry
    geom = build_geometry(plant, state)
    m = getattr(plant, "machine", None)
    scr = [(s["sid"], s.get("kind", "")) for s in (geom.get("screens") or [])] or \
        [(sid, "") for sid in read_screen_ids(m)]
    return {"experiment": exp.to_dict() if exp is not None else None,
            "geometry": geom, "undulators": list(getattr(state, "undulators", {}) or {}),
            "loops": loops_info(sch.loops), "optics": optics_specs(state),
            "observables": sorted(getattr(state, "observables", {}) or {}), "screens": scr,
            "targets": loop_targets(state, sch.loops),
            "cavity": bool(hasattr(getattr(plant, "tbt", None) or plant, "set_cavity")),
            "machine": ({"lattice_id": int(m.lattice_id), "fingerprint": m.fingerprint, "source": str(m.source)}
                        if m is not None else None)}
