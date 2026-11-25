# -*- coding: utf-8 -*-
"""Стандартный цикл SynchrotronSim (docs/LIBRARY_ARCHITECTURE.md, раздел 5.2; шаг 5).

    exp = Experiment.from_fb(3, smoke=True)
    run = Simulation(exp).run()                 # → synchrotronsim.run.Run (папка прогона с паспортом)
    run.series("rec.orbit_rms.x")

`build_session(exp)` — ЕДИНЫЙ сборщик (plant, доска, планировщик с контурами) для всех входов: блокнот,
`control/mode4_run.py`, кокпит (`client.engine.build_session`). Раньше сборщиков было два, и они расходились в
мелочах; теперь различия — данные опыта (`Experiment.plant.kwargs`, темпы контуров), а не код.

Этапы цикла и точки настройки: (1) предварительные расчёты — внутри plant'а и `setup` контуров (матрицы отклика,
замкнутая орбита; кэш по отпечатку); (2) контуры — плагины `Controller` в планировщике; (3) общая обработка —
рекордеры подписки `exp.record` (Watch: адрес или рекордер) по окнам усреднения; (4) выход — `Run`.
"""
from __future__ import annotations

from typing import Callable, Optional

from .experiment import Experiment, LoopSpec
from .run import Run

TBT_SERVES = ("ebpm.x", "ebpm.y", "photon.x", "photon.y", "rf.delta")   # наблюдаемые с сырьём у TbtPlant
OPTICS_SERVES = ("tune", "chroma", "adts")


def make_loop(l: LoopSpec):
    """LoopSpec → объект контура (класс по виду). None-поля — значения по умолчанию класса."""
    from control.controllers.integral_svd import IntegralSVD
    from control.controllers.rf_frequency import RfFrequencyLoop
    from control.controllers.optics_tune import TuneLoop
    from control.controllers.optics_chroma import ChromaLoop
    from control.controllers.optics_adts import AdtsLoop
    kw = {"name": l.name}
    for k in ("gain", "rcond"):
        if getattr(l, k) is not None:
            kw[k] = float(getattr(l, k))
    if l.ref is not None:
        kw["ref"] = l.ref or None                          # "" — явно без уставки
    if l.kind in ("integral_svd", "rf_frequency"):
        if l.sensor is not None:
            kw["sensor"] = l.sensor
        if l.actuator is not None:
            kw["actuator_key"] = l.actuator
    cls = {"integral_svd": IntegralSVD, "rf_frequency": RfFrequencyLoop, "optics_tune": TuneLoop,
           "optics_chroma": ChromaLoop, "optics_adts": AdtsLoop}[l.kind]
    return cls(**kw)


class Session:
    """Собранная связка опыта: plant, доска (state), планировщик (sch), контуры, такт dt и окно усреднения."""

    def __init__(self, exp, plant, state, sch, loops, dt, window, machine=None):
        self.exp, self.plant, self.state, self.sch, self.loops, self.dt, self.window = \
            exp, plant, state, sch, loops, dt, window
        self.machine = machine                             # документ машины, по которому собрано кольцо

    def close(self):
        try:
            self.plant.close()
        except Exception:
            pass


def loop_rate_hz(l: LoopSpec, dt: float, window: int, run) -> float:
    if l.rate == "auto":
        return 1.0 / (dt * window)                         # раз в окно усреднения
    if l.rate == "cycle":
        return 1.0 / (dt * max(1, int(run.steps_per_cycle)))
    if l.rate == "optics":
        return 1.0 / (dt * max(1, int(run.optics_period_ticks)))
    return float(l.rate)


def resolve_machine(exp: Experiment, store=None, machine=None):
    """Машина опыта: заданный документ или из хранилища по ссылке опыта (решётка + ревизия; "" — последняя)."""
    if machine is not None:
        if int(machine.lattice_id) != int(exp.lattice_id):
            raise ValueError(f"машина решётки {machine.lattice_id}, а опыт «{exp.name}» — на решётке {exp.lattice_id}")
        if exp.machine_rev and not machine.fingerprint.startswith(exp.machine_rev):
            raise ValueError(f"опыт «{exp.name}» ссылается на ревизию {exp.machine_rev[:12]}, а передана машина "
                             f"{machine.fingerprint[:12]}")
        return machine
    if store is None:
        from .store import open_store
        store = open_store()
    ref = exp.machine_ref()
    return store.get("machine", ref.id, ref.rev or None)


def build_session(exp: Experiment, store=None, machine=None) -> Session:
    """Опыт → (plant, доска, планировщик с контурами). Машина — документ: `machine=` или из хранилища `store` по
    ссылке опыта (`exp.machine_ref()`); plant'ы собирают кольцо из него, таблиц не читают. Тяжёлые импорты
    ленивые (ring_cpp и т.п.)."""
    from control import ExperimentState, Scheduler, NoiseModel, DRIVE_PROFILE_DEFAULTS
    from control.obs_assembly import observer_for
    from control.controllers.rf_frequency import RfFrequencyLoop

    machine = resolve_machine(exp, store, machine)
    loops = [make_loop(l) for l in exp.loops]
    sensors = list(loops)
    if "rf.delta" in exp.plant.extra_sensors and not any(getattr(c, "sensor", "") == "rf.delta" for c in loops):
        sensors.append(RfFrequencyLoop())                  # датчик энергии δ̂ для записи без контура ВЧ
    kw = dict(exp.plant.kwargs)
    kind = exp.plant.kind
    if kind == "optics":
        from control.optics_plant import OpticsPlant
        plant = OpticsPlant(lattice=exp.lattice_id, machine=machine, **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
        for q, v in plant.nominal().items():               # уставки контуров = проектная оптика
            state.refs[f"{q}.target"] = v
    elif kind == "composite":
        from control.composite_plant import CompositePlant
        plant = CompositePlant(lattice=exp.lattice_id,
                               tbt_observer=observer_for(sensors, serves=set(TBT_SERVES) if kw.get("rf_actuator")
                                                         else set(TBT_SERVES) - {"rf.delta"}),
                               optics_observer=observer_for(sensors, serves=set(OPTICS_SERVES)), machine=machine,
                               **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
    else:
        from control.tbt_plant import TbtPlant
        obs = observer_for(sensors, TBT_SERVES) if exp.plant.observer == "loops" else None
        plant = TbtPlant(lattice=exp.lattice_id, observer=obs, machine=machine, **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
    if exp.noise is not None:
        state.noise = NoiseModel(exp.noise)                # закон шума — в доске (живой)
    if exp.drive_profile is not None:                      # реалистичные приводы ВСЕХ магнитов (slew+лаг+ОДЗ)
        state.limits.set_profile(DRIVE_PROFILE_DEFAULTS if exp.drive_profile == "default" else exp.drive_profile)
    dt = (plant.stride * plant.TREV) if hasattr(plant, "stride") else 1.0   # optics: абстрактный такт
    window = int(getattr(plant, "avg_chunks", 1) or 1)
    sch = Scheduler(plant, state, dt)
    sch.machine = machine                                  # паспорт прогона: машина с ревизией (Runner фиксирует)
    for spec, ctl in zip(exp.loops, loops):
        sch.add(ctl, rate_hz=loop_rate_hz(spec, dt, window, exp.run), enabled=bool(spec.enabled))
    return Session(exp, plant, state, sch, loops, dt, window, machine=machine)


def records_of(exp: Experiment):
    """Подписка опыта → строки fb_config.Record для control.recorders.build (рекордеры и адреса)."""
    from control.fb_config import Record
    out = []
    for w in exp.record:
        params = w.params
        if not w.recorder and w.every > 1 and "every=" not in (params or ""):
            params = f"every={w.every}" + (f";{params}" if params else "")
        out.append(Record(w.recorder or w.addr, params, w.cadence, True))
    return out


class Runner:
    """Исполнитель опыта по шагам — ОДИН для пакетного прогона (`Simulation.run`) и живого кокпита
    (`client.engine.EngineThread` — тонкий адаптер: поток + шов Link + публикация снимка).

    Шаг планировщика, рестарт прохода (штатный reset plant'а, модельное время и вклады контуров — в ноль), лимит
    оборотов прохода, запись прогона: Watch-адреса (`watch`) пишутся на каждом `record_tick`, каждый проход — свой
    `Run` (`begin_pass`/`end_pass`), в паспорте — опыт с текущей подпиской и лимитом. Ошибки записи не роняют
    исполнение — они в `err`."""

    def __init__(self, sch, experiment=None, run_root: Optional[str] = None, kind: str = "cockpit",
                 source: str = "cockpit"):
        self.sch, self.state, self.plant = sch, sch.state, sch.plant
        self.experiment = experiment if experiment is not None else getattr(sch, "experiment", None)
        self.run_root, self.kind, self.source = run_root, kind, source
        self.watches: set = set()             # канонические адреса, которые пишутся в прогон
        self.rec_on = False                   # запись включена: каждый проход — свой Run
        self.max_turns: Optional[int] = None  # лимит оборотов прохода (None — без лимита)
        self.repeat = False                   # после лимита — новый проход (рестарт), иначе пауза
        self.run: Optional[Run] = None        # текущий прогон или None
        self.pass_no = 0                      # номер прохода (растёт на каждом рестарте)
        self.err: Optional[str] = None        # последняя ошибка записи

    # ---- шаг и проход ----
    def step(self) -> None:
        self.sch.step()

    def turn(self) -> int:
        return int(getattr(self.plant, "turn", 0) or 0)

    def limit_reached(self) -> bool:
        return bool(self.max_turns) and self.turn() >= self.max_turns

    def reset_model(self) -> None:
        """Новый проход с нуля: штатный reset plant'а, модельное время и такт — в ноль, вклады контуров прошлого
        прохода — в ноль; при включённой записи — новый прогон."""
        self.plant.reset(self.state)
        self.state.t = 0.0
        self.sch.tick = 0
        for slot in self.state.target_slots.values():
            slot.v[:] = 0.0
        self.pass_no += 1
        if self.rec_on:
            self.begin_pass()

    # ---- подписка и запись ----
    def watch(self, addr: str) -> bool:
        from . import addresses as AD
        if not AD.is_address(addr):
            return False
        self.watches.add(AD.canonical(addr))
        return True

    def unwatch(self, addr: str) -> None:
        from . import addresses as AD
        self.watches.discard(AD.canonical(addr))

    def set_recording(self, on: bool, max_turns: Optional[int] = None, repeat: bool = False) -> None:
        self.max_turns = int(max_turns) if max_turns else None
        self.repeat = bool(repeat)
        if on and not self.rec_on:
            self.rec_on = True
            self.begin_pass()
        elif not on and self.rec_on:
            self.rec_on = False
            self.end_pass("stopped")

    def experiment_snapshot(self):
        """Опыт с ТЕКУЩЕЙ подпиской и лимитом — в паспорт прогона; без опыта — сводка исполнителя."""
        import dataclasses
        from .run import Watch
        exp = self.experiment
        if exp is None:
            return {"engine": self.kind, "plant": type(self.plant).__name__,
                    "watches": sorted(self.watches), "max_turns": self.max_turns, "repeat": self.repeat}
        return dataclasses.replace(exp, record=[Watch(a, cadence="per_step") for a in sorted(self.watches)],
                                   run=dataclasses.replace(exp.run, max_turns=self.max_turns, repeat=self.repeat))

    def _machine(self):
        """Машина прогона: та, по которой собрано кольцо (sch.machine / plant.machine), иначе — из хранилища."""
        m = getattr(self.sch, "machine", None) or getattr(self.plant, "machine", None)
        if m is not None:
            return m
        lat = getattr(self.experiment, "lattice_id", None) if self.experiment is not None else None
        lat = lat if lat is not None else getattr(self.plant, "lattice", None)
        if lat is None:
            return None
        try:
            from .store import open_store
            return open_store(root=self.run_root).machine(int(lat))
        except Exception:
            return None

    def _pin(self, machine):
        """Зафиксировать ревизию машины в хранилище (снимок документа) → текст ссылки или None (ошибка — в err)."""
        if machine is None:
            return None
        try:
            from .store import open_store
            return str(open_store(root=self.run_root).pin(machine))
        except Exception as e:                            # запись снимка не должна ронять исполнение
            self.err = f"снимок машины: {type(e).__name__}: {e}"
            return None

    def begin_pass(self, config=None, note: Optional[str] = None) -> Optional[Run]:
        """Начать прогон (прежний, если был, — закрыть как stopped). config/note — для пакетного прогона."""
        self.end_pass("stopped")
        try:
            m = self._machine()
            self.run = Run.create(self.kind, config=self.experiment_snapshot() if config is None else config,
                                  machine=m, machine_ref=self._pin(m), seed=getattr(self.experiment, "seed", None),
                                  source=self.source,
                                  note=note if note is not None else
                                  f"проход {self.pass_no}; старт с оборота {self.turn()}", root=self.run_root)
        except Exception as e:                            # запись не должна ронять исполнение
            self.run = None
            self.err = f"запись прогона: {type(e).__name__}: {e}"
        return self.run

    def end_pass(self, status: str) -> None:
        run, self.run = self.run, None
        if run is not None:
            try:
                run.finish(status)
            except Exception as e:
                self.err = f"запись прогона: {type(e).__name__}: {e}"

    def record_tick(self, step: int, meta: Optional[dict] = None) -> None:
        """Точка прогона по каждому Watch-адресу: шаг = номер тика (у кокпита — цикл публикации), время — модельное."""
        run = self.run
        if run is None or not self.watches:
            return
        from control import StateView
        from . import addresses as AD
        view = StateView(self.state)
        mt = float(self.state.t)
        m = {"model_time": mt, "cycle": int(step), "turn": self.turn()}
        m.update(meta or {})
        for a in sorted(self.watches):
            v = AD.resolve(view, a, loops=self.sch.loops, meta=m)
            if v is not None:
                run.record(a, step, v, t=mt)


class Simulation:
    """Исполнение опыта без UI: шаги планировщика, рекордеры по окнам, прогон Run с паспортом.

    legacy_fb — для опытов из БД (source «fb_exp_id=N») писать и прежние fb_results/fb_artifacts/статус."""

    def __init__(self, exp: Experiment, root: Optional[str] = None, legacy_fb: Optional[bool] = None,
                 progress: Optional[Callable[[str], None]] = print, store=None, machine=None):
        self.exp = exp
        self.root = root
        self.store, self.machine = store, machine         # хранилище и/или готовый документ машины
        self.legacy_fb = exp.source.startswith("fb_exp_id=") if legacy_fb is None else bool(legacy_fb)
        self.progress = progress or (lambda s: None)
        self.session: Optional[Session] = None

    def build(self) -> Session:
        if self.session is None:
            self.session = build_session(self.exp, store=self.store, machine=self.machine)
        return self.session

    def _fb_id(self):
        return int(self.exp.source.split("=", 1)[1]) if self.legacy_fb else None

    def run(self, steps: Optional[int] = None, close: bool = True) -> Run:
        """Пакетный прогон через Runner: один Run (паспорт — полный опыт), рекордеры подписки по окнам."""
        from control import recorders
        exp = self.exp
        S = self.build()
        n_steps = int(steps if steps is not None else (exp.run.steps or 0))
        if n_steps <= 0:
            raise ValueError("Simulation.run: не задано число шагов (exp.run.steps или steps=)")
        recs = recorders.build(records_of(exp), mesh=exp.run.record_mesh, loops=S.sch.loops, machine=S.machine)
        R = Runner(S.sch, experiment=exp, run_root=self.root, kind="fb_mode4" if self.legacy_fb else "simulation",
                   source=exp.source)
        R.max_turns = exp.run.max_turns
        run = R.begin_pass(config=exp, note=exp.note)
        if run is None:
            raise RuntimeError(R.err or "не удалось создать прогон")
        fb = None
        if self.legacy_fb:
            from control import fb_store
            fb = fb_store
            fb.set_status(self._fb_id(), "in_progress")
        t_win, status = {}, "failed"
        try:
            widx = 0
            window = S.window
            for step in range(n_steps):
                R.step()
                if (step + 1) % window == 0:
                    t_win[widx] = float(S.state.t)
                    for rec in recs:
                        rec.sample(S.state, widx, S.plant)
                    widx += 1
                    if widx % 4 == 0 or step + 1 == n_steps:
                        last = recs[0].rows[-1] if (recs and recs[0].rows) else None
                        msg = f"{last[0]}={last[2]:.4g}" if last else "—"
                        self.progress(f"  окно {widx:3d}/{max(1, n_steps // window)}: {msg}")
                if R.limit_reached():
                    self.progress(f"  лимит {R.max_turns} оборотов — стоп")
                    break
            for rec in recs:
                rec.finalize(S.state, S.plant)
            rows = []
            for rec in recs:
                rows += rec.result_rows()
                run.record_rows(rec.series_rows(), prefix=rec.prefix, t_of_step=t_win)
                art = rec.artifact()
                if art:
                    run.artifact(art[0], art[1], art[2])
                    if fb is not None:
                        rel = fb.write_artifact(self._fb_id(), art[0], art[1], art[2])
                        self.progress(f"  артефакт [{art[0]}] → {rel}")
            if fb is not None:
                n = fb.write_results(self._fb_id(), rows)
                fb.set_status(self._fb_id(), "done")
                self.progress(f"\n[OK] записано {n} строк результатов в fb_results.csv; статус эксперимента → done")
            status = "done"
        except KeyboardInterrupt:
            status = "stopped"
            raise
        finally:
            R.end_pass(status)
            if close:
                S.close()
                self.session = None
        return run
