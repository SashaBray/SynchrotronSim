# -*- coding: utf-8 -*-
"""Режим 4 — ТОЧКА ВХОДА эксперимента с обратной связью: конфиг из БД по ключу → моделирование →
результаты в БД и в прогон. Отчёт строится ОТДЕЛЬНО, по требованию (control/mode4_report.py).

С шага 5 плана библиотеки (docs/LIBRARY_ARCHITECTURE.md) это тонкая обёртка над стандартным циклом:
    exp = synchrotronsim.experiment.Experiment.from_fb(fb_exp_id, smoke)   # таблицы fb_* → опыт (объект)
    run = synchrotronsim.simulation.Simulation(exp).run()                  # сборка + шаги + рекордеры → Run
Тот же опыт можно запустить из блокнота теми же двумя строками.

Поток: опыт → TbtPlant (нелинейное 6D-кольцо СИЛА) с многоканальным возмущением из disturbance_spec (или
OpticsPlant для plant=optics) + Scheduler с контурами по fb_loops (класс — по датчику: rf.delta → контур ВЧ-частоты
RfFrequencyLoop; прочие — IntegralSVD) → по окнам усреднения рекордеры ПОДПИСКИ (fb_record: signal — имя
рекордера или любой адрес единого пространства, напр. obs.tune / undulator.und_01.gap / act.corr.000.hkick) →
результаты в fb_results.csv + артефакты (.npz + fb_artifacts) и в ПРОГОН workspace/runs/<id> (manifest.json с
паспортом, signals.npz, tables/signals.csv; индекс — workspace/tables/runs.csv).
Статус эксперимента (new→in_progress→done) — прерывание-безопасность: упал на середине → остался
in_progress. Диффур-ядро/поля заморожены; всё через apply_perturb.

ЗАПУСК: py -3.14 control/mode4_run.py --exp 1 [smoke]
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from control.controllers.integral_svd import IntegralSVD
from control.controllers.rf_frequency import RfFrequencyLoop
from synchrotronsim.experiment import Experiment
from synchrotronsim.simulation import Simulation, TBT_SERVES  # noqa: F401  (TBT_SERVES — совместимость)

# датчик контура → класс плагина (свой вход в схему вычисления наблюдаемой — provides()); прочие — IntegralSVD
LOOP_CLASSES = {"rf.delta": RfFrequencyLoop}


def make_loop(l):
    """fb_config.Loop → объект контура (класс по датчику). Совместимость: опыт строит контуры сам
    (synchrotronsim.simulation.make_loop по LoopSpec)."""
    cls = LOOP_CLASSES.get(l.sensor, IntegralSVD)
    return cls(name=l.name, sensor=l.sensor, actuator_key=l.actuator_key, gain=l.gain, rcond=l.rcond,
               ref=(l.ref or None))


def main():
    ap = argparse.ArgumentParser(description="режим 4 — эксперимент ОС из БД")
    ap.add_argument("mode", nargs="?", default="")
    ap.add_argument("--exp", type=int, required=True, help="ключ эксперимента fb_exp_id")
    a = ap.parse_args()
    smoke = a.mode == "smoke"

    exp = Experiment.from_fb(a.exp, smoke=smoke)              # конфиг ИЗ БД по ключу (таблицы fb_*) → опыт

    # --- АЛЬТЕРНАТИВА: опыт БЕЗ чтения таблиц fb_* -----------------------------------------------------
    # Закомментируйте строку выше и соберите опыт прямо в коде — например, орбитный FOFB на СИЛА-26:
    #
    # from synchrotronsim.experiment import LoopSpec, PlantSpec, RunLimits
    # from synchrotronsim.run import Watch
    # exp = Experiment(
    #     name="fofb_hardcoded", lattice_id=26, seed=1, source="fb_exp_id=1",   # source — куда писать fb_results
    #     plant=PlantSpec("tbt", dict(nsl=8, sr=True, stride=15, settle_turns=250, seed=1, dist_seed=1,
    #                                 n_photon=8, relative_orbit=False, center_bunch=None, rf_actuator=False)),
    #     loops=[LoopSpec("fofb_x", sensor="ebpm.x", actuator="hkick", gain=0.5, rcond=0.01, ref=""),
    #            LoopSpec("fofb_y", sensor="ebpm.y", actuator="vkick", gain=0.5, rcond=0.01, ref="")],
    #     noise=[{"name": "ground_vibration", "target": "quad", "plane": "both", "form": "tone",
    #             "correlation": "girder", "quantity": "position", "amplitude": 3e-7, "freq_hz": 12.0,
    #             "phase_deg": 0.0}],
    #     record=[Watch("orbit_rms", recorder="orbit_rms"), Watch("obs.ebpm.x:rms")],
    #     run=RunLimits(steps=180))
    # ------------------------------------------------------------------------------------------------

    pk = exp.plant.kwargs
    print(f"=== РЕЖИМ 4: эксперимент ОС {exp.source} ({exp.plant.kind}) ===")
    print(f"решётка {exp.lattice_id} · {exp.run.steps} шагов · nsl {pk.get('nsl')}"
          + (f" · СИ {'вкл' if pk.get('sr') else 'выкл'} · stride {pk.get('stride')} · осадка {pk.get('settle_turns')}"
             if exp.plant.kind == "tbt" else ""))
    if exp.plant.kind == "optics":
        print(f"возмущение: дрейф градиентов k1/k2/k3 = {(pk.get('drift_k1'), pk.get('drift_k2'), pk.get('drift_k3'))}")
    else:
        print(f"возмущения: {', '.join(c['name'] for c in (exp.noise or [])) or '—'}")
    print(f"контуры: {', '.join(l.name for l in exp.loops if l.enabled) or '—'}")
    print(f"подписка (запись): {', '.join(w.addr + ('[' + w.params + ']' if w.params else '') for w in exp.record) or '—'}\n",
          flush=True)

    sim = Simulation(exp, progress=lambda s: print(s, flush=True))
    S = sim.build()
    if exp.plant.kind == "tbt" and (pk.get("relative_orbit") or getattr(S.plant, "rf_ids", None)):
        print(f"BPM {'относительно опорной орбиты' if pk.get('relative_orbit') else 'абсолютно'}; "
              f"актуатор ВЧ-частоты: {'да (номинал ' + format(S.plant.F0 / 1e6, '.6f') + ' МГц)' if S.plant.rf_ids else 'нет'}")
    if exp.drive_profile is not None:
        print(f"реалистичные приводы (профиль {'дефолт' if exp.drive_profile == 'default' else 'БД'})")
    run = sim.run()
    print(f"прогон: {os.path.relpath(run.path)} (manifest.json, signals.npz, tables/signals.csv)")
    if exp.source.startswith("fb_exp_id="):
        print(f"отчёт: py -3.14 control/mode4_report.py --exp {exp.source.split('=', 1)[1]}")


if __name__ == "__main__":
    main()
