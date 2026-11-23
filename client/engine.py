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
                   COMMAND, CONTROL, SET_SCREEN)

META_NAMES = ("model_time", "wall_time", "speed", "cycle", "seq", "running", "turn", "alive")


class EngineThread(threading.Thread):
    """Поток-движок над готовым Scheduler. start() запускает цикл; control('run'/'stop') им управляет;
    stop() (или control('shutdown')) гасит поток. Публикует снимок ТОЛЬКО подписанных адресов + мету."""

    def __init__(self, link: Link, scheduler, steps_per_cycle: int = 9, publish_every: int = 3,
                 idle_sleep: float = 0.003, screen_period: float = 1.0, name: str = "engine"):
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
                self.sch.step()                           # ОДИН базовый такт = stride оборотов
                self._wall_run += perf_counter() - t0
                self._step += 1
                if self._step % self.publish_every == 0:  # снимок раз в publish_every тактов (не каждый — дешевле)
                    self._cycle += 1
                    self._publish()
        except Exception as e:                            # поток не должен падать молча
            self._err = f"{type(e).__name__}: {e}"
            try:
                self._publish()
            except Exception:
                pass

    # ---- очередь команд ----
    def _drain(self) -> None:
        for m in self.link.pull():                        # subscribe/unsubscribe уже применены в pull()
            if m.kind == CONTROL:
                self._control(m.value)
            elif m.kind == COMMAND:
                self._command(m.addr, m.value)
            elif m.kind == SET_SCREEN:
                self._screens[m.addr] = dict(m.value or {})

    def _control(self, action: str) -> None:
        if action == RUN:
            self._running = True
        elif action == STOP:
            self._running = False
        elif action == RESTART:
            was = self._running
            self._running = False
            self.plant.reset(self.state)                  # новый прогон с нуля (штатный reset плана)
            self.state.t = 0.0
            self._cycle = 0
            self._step = 0
            self._wall_run = 0.0
            self.sch.tick = 0
            for slot in self.state.target_slots.values():  # обнулить вклады контуров прошлого прогона
                slot.v[:] = 0.0
            self._running = was                           # рестарт сохраняет режим пуск/пауза
        elif action == SHUTDOWN:
            self._stop.set()

    def _command(self, addr: str, value) -> None:
        """Применить ручку ШТАТНЫМ путём (команда → доска). Диспетчер по типу адреса. Неизвестное —
        записать в диагностику, не падать."""
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
                elif fld == "gain" and hasattr(spec.controller, "gain"):
                    spec.controller.gain = float(value)
                else:
                    self._mark_unknown(addr)
            elif kind == "drive":
                aid, _, fld = addr[len("drive."):].rpartition(".")
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

    # ---- публикация снимка ----
    def _meta(self) -> dict:
        mt = float(self.state.t)
        speed = (mt / self._wall_run) if self._wall_run > 0 else 0.0
        m = {"model_time": mt, "wall_time": float(self._wall_run), "speed": float(speed),
             "cycle": int(self._cycle), "seq": int(self._seq), "running": bool(self._running),
             "turn": int(getattr(self.plant, "turn", 0)), "alive": bool(getattr(self.plant, "alive", True))}
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
        ckey = (int(cfg.get("res") or 121), tuple(float(v) for v in (cfg.get("point") or (0.0, 0.0))))
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
        try:
            from control.recorders import ScreenRecorder
            rec = ScreenRecorder(params=sid, mesh=mesh, point=pt)
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
            return out
        except Exception as e:                      # экран опционален/тяжёл — ошибку отдаём окну, не падаем
            return {"error": f"{type(e).__name__}: {e}"}


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
    from control import ExperimentState, Limits, Scheduler
    from control.controllers.integral_svd import IntegralSVD

    loops = []                                                 # объекты контуров ПЕРВЫМИ (нужны для сборки Observer'ов)
    if fofb:
        loops += [IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5, rcond=1e-2),
                  IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=0.5, rcond=1e-2)]
    if optics:
        from control.controllers.optics_tune import TuneLoop
        from control.controllers.optics_chroma import ChromaLoop
        from control.controllers.optics_adts import AdtsLoop
        loops += [TuneLoop(), ChromaLoop(), AdtsLoop()]
    if rf:
        from control.controllers.rf_frequency import RfFrequencyLoop
        loops += [RfFrequencyLoop()]
    tbt_serves = {"ebpm.x", "ebpm.y", "photon.x", "photon.y"} | ({"rf.delta"} if rf else set())

    if optics:
        from control.composite_plant import CompositePlant
        from control.obs_assembly import observer_for
        tbt_obs = observer_for(loops, serves=tbt_serves)
        optics_obs = observer_for(loops, serves={"tune", "chroma", "adts"})
        plant = CompositePlant(lattice=lattice, nsl=nsl, settle_turns=settle_turns, meas_noise=meas_noise,
                               n_photon=n_photon, gap_mm=gap_mm, disturbance_spec=disturbance_spec,
                               drift=optics_drift, order=optics_order, optics_nsl=optics_nsl,
                               tbt_observer=tbt_obs, optics_observer=optics_obs, rf_actuator=rf,
                               id_pole_nsl=id_pole_nsl, **plant_kw)
    else:
        from control.tbt_plant import TbtPlant
        obs = None
        if rf:
            from control.obs_assembly import observer_for
            obs = observer_for(loops, serves=tbt_serves)
        plant = TbtPlant(lattice=lattice, nsl=nsl, settle_turns=settle_turns, meas_noise=meas_noise,
                         n_photon=n_photon, gap_mm=gap_mm, disturbance_spec=disturbance_spec,
                         rf_actuator=rf, observer=obs, id_pole_nsl=id_pole_nsl, **plant_kw)

    dt = plant.stride * plant.TREV
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, dt)
    rate_f = 1.0 / (dt * steps_per_cycle)                      # FOFB: раз в управляющий цикл (≈2 кГц)
    rate_o = 1.0 / (dt * max(1, optics_period_ticks))         # оптика: медленный темп (квазистатика)
    for lp in loops:
        sch.add(lp, rate_hz=(rate_f if lp.name in ("fofb_x", "fofb_y", "rf") else rate_o), enabled=False)
    return plant, state, sch
