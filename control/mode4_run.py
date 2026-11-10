# -*- coding: utf-8 -*-
"""Режим 4 — ТОЧКА ВХОДА эксперимента с обратной связью: конфиг из БД по ключу → моделирование →
результаты в БД. Отчёт строится ОТДЕЛЬНО, по требованию (control/mode4_report.py).

Поток: fb_config.load(fb_exp_id) → TbtPlant (нелинейное 6D-кольцо СИЛА) с многоканальным возмущением из
disturbance_spec + Scheduler с контурами (IntegralSVD по fb_loops) → по окнам усреднения зовём рекордеры
ПОДПИСКИ (fb_record: за чем следить) → результаты пишем в fb_results.csv + артефакты (.npz + fb_artifacts).
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

from control import ExperimentState, Scheduler, NoiseModel
from control.tbt_plant import TbtPlant
from control.controllers.integral_svd import IntegralSVD
from control import fb_config, fb_store, recorders


def main():
    ap = argparse.ArgumentParser(description="режим 4 — эксперимент ОС из БД")
    ap.add_argument("mode", nargs="?", default="")
    ap.add_argument("--exp", type=int, required=True, help="ключ эксперимента fb_exp_id")
    a = ap.parse_args()
    smoke = a.mode == "smoke"

    cfg = fb_config.load(a.exp)                              # конфиг ИЗ БД по ключу (таблицы fb_*)

    # --- АЛЬТЕРНАТИВА: ЗАХАРДКОЖЕННЫЙ конфиг БЕЗ чтения БД --------------------------------------------
    # Закомментируйте строку `cfg = fb_config.load(a.exp)` выше и РАСКОММЕНТИРУЙТЕ блок ниже, чтобы
    # прогнать опыт на данных, зашитых прямо в код (эквивалент fb_exp_id=1) — ни одна таблица fb_*/
    # disturbance_spec.csv при этом не читается. Контуры (loops), возмущения (disturbance_spec) и подписка
    # (records) объявлены явно. Результаты всё равно пишутся в БД по cfg.fb_exp_id (set_status/write_* —
    # это запись, не чтение; при отсутствии fb_experiments.csv set_status просто ничего не делает).
    #
    # cfg = fb_config.FbConfig(
    #     fb_exp_id=1, lattice_id=26, steps=180, stride=15, settle=250, seed=1, sr=True, nsl=8,
    #     note="hardcoded (без БД)",
    #     loops=[
    #         #        name      sensor    actuator_key gain  rcond rate_hz ref  enabled
    #         fb_config.Loop("fofb_x", "ebpm.x", "hkick", 0.5, 0.01, "auto", "", True),
    #         fb_config.Loop("fofb_y", "ebpm.y", "vkick", 0.5, 0.01, "auto", "", True),
    #     ],
    #     disturbance_spec=[   # формат строк — как у disturbance.load_spec(): амплитуды в СИ (м / отн.доля)
    #         {"name": "ground_vibration", "target": "quad", "plane": "both", "form": "tone",
    #          "correlation": "girder", "quantity": "position", "amplitude": 3e-7, "freq_hz": 12.0, "phase_deg": 0.0},
    #         {"name": "fast_tone", "target": "quad", "plane": "x", "form": "tone",
    #          "correlation": "independent", "quantity": "position", "amplitude": 1e-7, "freq_hz": 40.0, "phase_deg": 0.0},
    #         {"name": "dip_field_ripple", "target": "dipole", "plane": "", "form": "tone",
    #          "correlation": "common", "quantity": "field:b1", "amplitude": 5e-5, "freq_hz": 50.0, "phase_deg": 0.0},
    #     ],
    #     records=[            # подписка: signal, params, cadence, enabled
    #         fb_config.Record("orbit_rms", "", "per_window", True),
    #         fb_config.Record("light", "", "per_window", True),
    #         fb_config.Record("phase", "", "final", True),
    #         fb_config.Record("screen", "SCR_U19", "final", True),
    #     ],
    # )
    # ------------------------------------------------------------------------------------------------

    if smoke:
        steps = min(cfg.steps, 8 if cfg.plant == "optics" else 36)   # optics: normal_form на шаг дорог → меньше шагов
    else:
        steps = cfg.steps
    settle = min(cfg.settle, 120) if smoke else cfg.settle
    mesh = 81 if smoke else 121

    print(f"=== РЕЖИМ 4: эксперимент ОС fb_exp_id={cfg.fb_exp_id} ({cfg.plant}) ===")
    print(f"решётка {cfg.lattice_id} · {steps} шагов · nsl {cfg.nsl}"
          + (f" · СИ {'вкл' if cfg.sr else 'выкл'} · stride {cfg.stride} · осадка {settle}" if cfg.plant == "tbt" else ""))
    if cfg.plant == "optics":
        print(f"возмущение: дрейф градиентов k1/k2/k3 = {cfg.drift}")
    else:
        print(f"возмущения: {', '.join(c['name'] for c in cfg.disturbance_spec) or '—'}")
    print(f"контуры: {', '.join(l.name for l in cfg.enabled_loops()) or '—'}")
    print(f"подписка (запись): {', '.join(r.signal + ('[' + r.params + ']' if r.params else '') for r in cfg.enabled_records()) or '—'}\n", flush=True)

    # --- plant по типу эксперимента (plant-агностичный режим 4) ---
    if cfg.plant == "optics":                                # квазистатика: оптические контуры на OpticsPlant
        from control.optics_plant import OpticsPlant
        plant = OpticsPlant(lattice=cfg.lattice_id, nsl=cfg.nsl, observe=("tune", "chroma", "adts"),
                            drift_k1=cfg.drift[0], drift_k2=cfg.drift[1], drift_k3=cfg.drift[2],
                            rho=1.0, seed=cfg.seed)
        state = ExperimentState(plant.schema(), history=plant.history())
        for q, v in plant.nominal().items():                 # уставки контуров = проектная оптика
            state.refs[f"{q}.target"] = v
    else:                                                    # tbt: орбита turn-by-turn на TbtPlant
        plant = TbtPlant(lattice=cfg.lattice_id, nsl=cfg.nsl, sr=cfg.sr, stride=cfg.stride,
                         settle_turns=settle, seed=cfg.seed, dist_seed=cfg.seed, n_photon=8)
        state = ExperimentState(plant.schema(), history=plant.history())
        state.noise = NoiseModel(cfg.disturbance_spec)       # закон шума — в доске (живой)

    dt = (plant.stride * plant.TREV) if hasattr(plant, "stride") else 1.0   # tbt: по оборотам; optics: абстрактный такт
    window = getattr(plant, "avg_chunks", 1)                 # окно усреднения (optics — без усреднения)
    sch = Scheduler(plant, state, dt)
    for l in cfg.enabled_loops():
        rate = 1.0 / (dt * window) if l.rate_hz == "auto" else float(l.rate_hz)
        sch.add(IntegralSVD(name=l.name, sensor=l.sensor, actuator_key=l.actuator_key,
                            gain=l.gain, rcond=l.rcond, ref=(l.ref or None)), rate_hz=rate)
    recs = recorders.build(cfg.enabled_records(), mesh=mesh)

    fb_store.set_status(cfg.fb_exp_id, "in_progress")
    try:
        widx = 0
        for step in range(steps):
            sch.step()
            if (step + 1) % window == 0:
                for rec in recs:
                    rec.sample(state, widx, plant)
                widx += 1
                if widx % 4 == 0 or step + 1 == steps:
                    last = recs[0].rows[-1] if (recs and recs[0].rows) else None
                    msg = f"{last[0]}={last[2]:.4g}" if last else "—"
                    print(f"  окно {widx:3d}/{max(1, steps // window)}: {msg}", flush=True)
        for rec in recs:
            rec.finalize(state, plant)
        rows = []
        for rec in recs:
            rows += rec.result_rows()
            art = rec.artifact()
            if art:
                rel = fb_store.write_artifact(cfg.fb_exp_id, art[0], art[1], art[2])
                print(f"  артефакт [{art[0]}] → {rel}", flush=True)
        n = fb_store.write_results(cfg.fb_exp_id, rows)
        fb_store.set_status(cfg.fb_exp_id, "done")
        print(f"\n[OK] записано {n} строк результатов в fb_results.csv; статус эксперимента → done")
        print(f"отчёт: py -3.14 control/mode4_report.py --exp {cfg.fb_exp_id}")
    finally:
        plant.close()


if __name__ == "__main__":
    main()
