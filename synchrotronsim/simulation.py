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

    def __init__(self, exp, plant, state, sch, loops, dt, window):
        self.exp, self.plant, self.state, self.sch, self.loops, self.dt, self.window = \
            exp, plant, state, sch, loops, dt, window

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


def build_session(exp: Experiment) -> Session:
    """Опыт → (plant, доска, планировщик с контурами). Тяжёлые импорты ленивые (ring_cpp и т.п.)."""
    from control import ExperimentState, Scheduler, NoiseModel, DRIVE_PROFILE_DEFAULTS
    from control.obs_assembly import observer_for
    from control.controllers.rf_frequency import RfFrequencyLoop

    loops = [make_loop(l) for l in exp.loops]
    sensors = list(loops)
    if "rf.delta" in exp.plant.extra_sensors and not any(getattr(c, "sensor", "") == "rf.delta" for c in loops):
        sensors.append(RfFrequencyLoop())                  # датчик энергии δ̂ для записи без контура ВЧ
    kw = dict(exp.plant.kwargs)
    kind = exp.plant.kind
    if kind == "optics":
        from control.optics_plant import OpticsPlant
        plant = OpticsPlant(lattice=exp.lattice_id, **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
        for q, v in plant.nominal().items():               # уставки контуров = проектная оптика
            state.refs[f"{q}.target"] = v
    elif kind == "composite":
        from control.composite_plant import CompositePlant
        plant = CompositePlant(lattice=exp.lattice_id,
                               tbt_observer=observer_for(sensors, serves=set(TBT_SERVES) if kw.get("rf_actuator")
                                                         else set(TBT_SERVES) - {"rf.delta"}),
                               optics_observer=observer_for(sensors, serves=set(OPTICS_SERVES)), **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
    else:
        from control.tbt_plant import TbtPlant
        obs = observer_for(sensors, TBT_SERVES) if exp.plant.observer == "loops" else None
        plant = TbtPlant(lattice=exp.lattice_id, observer=obs, **kw)
        state = ExperimentState(plant.schema(), history=plant.history())
    if exp.noise is not None:
        state.noise = NoiseModel(exp.noise)                # закон шума — в доске (живой)
    if exp.drive_profile is not None:                      # реалистичные приводы ВСЕХ магнитов (slew+лаг+ОДЗ)
        state.limits.set_profile(DRIVE_PROFILE_DEFAULTS if exp.drive_profile == "default" else exp.drive_profile)
    dt = (plant.stride * plant.TREV) if hasattr(plant, "stride") else 1.0   # optics: абстрактный такт
    window = int(getattr(plant, "avg_chunks", 1) or 1)
    sch = Scheduler(plant, state, dt)
    for spec, ctl in zip(exp.loops, loops):
        sch.add(ctl, rate_hz=loop_rate_hz(spec, dt, window, exp.run), enabled=bool(spec.enabled))
    return Session(exp, plant, state, sch, loops, dt, window)


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


class Simulation:
    """Исполнение опыта без UI: шаги планировщика, рекордеры по окнам, прогон Run с паспортом.

    legacy_fb — для опытов из БД (source «fb_exp_id=N») писать и прежние fb_results/fb_artifacts/статус."""

    def __init__(self, exp: Experiment, root: Optional[str] = None, legacy_fb: Optional[bool] = None,
                 progress: Optional[Callable[[str], None]] = print):
        self.exp = exp
        self.root = root
        self.legacy_fb = exp.source.startswith("fb_exp_id=") if legacy_fb is None else bool(legacy_fb)
        self.progress = progress or (lambda s: None)
        self.session: Optional[Session] = None

    def build(self) -> Session:
        if self.session is None:
            self.session = build_session(self.exp)
        return self.session

    def _fb_id(self):
        return int(self.exp.source.split("=", 1)[1]) if self.legacy_fb else None

    def run(self, steps: Optional[int] = None, close: bool = True) -> Run:
        from control import recorders
        from .db import open_db
        exp = self.exp
        S = self.build()
        n_steps = int(steps if steps is not None else (exp.run.steps or 0))
        if n_steps <= 0:
            raise ValueError("Simulation.run: не задано число шагов (exp.run.steps или steps=)")
        recs = recorders.build(records_of(exp), mesh=exp.run.record_mesh, loops=S.sch.loops)
        try:
            machine = open_db(self.root).machine(exp.lattice_id)
        except Exception:
            machine = None
        run = Run.create("fb_mode4" if self.legacy_fb else "simulation", config=exp, machine=machine, seed=exp.seed,
                         source=exp.source, note=exp.note, root=self.root)
        fb = None
        if self.legacy_fb:
            from control import fb_store
            fb = fb_store
            fb.set_status(self._fb_id(), "in_progress")
        t_win, status = {}, "failed"
        try:
            widx = 0
            window = S.window
            max_turns = exp.run.max_turns
            for step in range(n_steps):
                S.sch.step()
                if (step + 1) % window == 0:
                    t_win[widx] = float(S.state.t)
                    for rec in recs:
                        rec.sample(S.state, widx, S.plant)
                    widx += 1
                    if widx % 4 == 0 or step + 1 == n_steps:
                        last = recs[0].rows[-1] if (recs and recs[0].rows) else None
                        msg = f"{last[0]}={last[2]:.4g}" if last else "—"
                        self.progress(f"  окно {widx:3d}/{max(1, n_steps // window)}: {msg}")
                if max_turns and int(getattr(S.plant, "turn", 0) or 0) >= max_turns:
                    self.progress(f"  лимит {max_turns} оборотов — стоп")
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
            run.finish(status)
            if close:
                S.close()
                self.session = None
        return run
