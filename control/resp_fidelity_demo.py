# -*- coding: utf-8 -*-
"""Опция 7: горизонтальный пол ×3.0 из опции 3 — зазор МОДЕЛЬ↔РЕАЛЬНОСТЬ.

Гипотеза: горизонталь не добивается, потому что контур обращает 4D-модель отклика (δ=0, без СИ), а
реальный 6D-пучок горизонтально несёт дисперсию/энергию. Вертикаль дисперсии не имеет → 4D-модель там
точна → ×50. Проверка: снять модель контура на СОБСТВЕННОМ 6D-кольце plant'а (аналог ИЗМЕРЕННОЙ матрицы
LOCO) и сравнить подавление с 4D-моделью. Та же платформа, тот же плагин, меняется ТОЛЬКО источник
матрицы отклика (resp_dims=4 vs 6) — ровно ради этого отклик сделан РЕСУРСОМ.

ЗАПУСК: py -3.14 control/resp_fidelity_demo.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control import ExperimentState, Limits, Scheduler
from control.tbt_plant import TbtPlant
from control.controllers.integral_svd import IntegralSVD

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
N = 60 if SMOKE else 180
TACT_CHUNKS = 9


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def run(resp_dims, enabled):
    plant = TbtPlant(lattice=26, nsl=8, settle_turns=(150 if SMOKE else 300),
                     dist_amp=0.3e-6, meas_noise=1e-7, resp_dims=resp_dims)
    dt = plant.stride * plant.TREV
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, dt)
    sch.add(IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5, rcond=1e-2),
            rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    sch.add(IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=0.5, rcond=1e-2),
            rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    ox, oy = [], []
    for _ in range(N):
        sch.step()
        ox.append(np.asarray(state.observables["ebpm.x"]))
        oy.append(np.asarray(state.observables["ebpm.y"]))
    tail = slice(len(ox) // 2, None)
    alive = plant.alive
    plant.close()
    return rms(np.mean(ox[tail], axis=0)), rms(np.mean(oy[tail], axis=0)), alive


def main():
    print("=== ОПЦИЯ 7: горизонтальный пол — 4D-модель контура vs 6D (LOCO-подобная) ===")
    print(f"{N} чанков×15 об; искажение feed-down 0.3 мкм; контур 2 кГц; меняется ТОЛЬКО источник отклика\n")

    bx, by, _ = run(resp_dims=4, enabled=False)                # базовая орбита (контур выкл)
    print("снимаю 4D-модель и гоняю FOFB…", flush=True)
    x4, y4, a4 = run(resp_dims=4, enabled=True)
    print("снимаю 6D-модель (на своём кольце, dims=6 — первый раз дольше) и гоняю FOFB…", flush=True)
    x6, y6, a6 = run(resp_dims=6, enabled=True)

    print(f"\n{'режим':26} | {'орбита X':>9} | {'орбита Y':>9} | пучок")
    print(f"{'без контура':26} | {bx*1e6:8.3f} | {by*1e6:8.3f} | —")
    print(f"{'FOFB, 4D-модель':26} | {x4*1e6:8.3f} | {y4*1e6:8.3f} | {'жив' if a4 else 'ПОТЕРЯН'}")
    print(f"{'FOFB, 6D-модель (LOCO)':26} | {x6*1e6:8.3f} | {y6*1e6:8.3f} | {'жив' if a6 else 'ПОТЕРЯН'}")
    print(f"\nгоризонталь: без {bx*1e6:.3f} → 4D ×{bx/max(x4,1e-30):.1f} ({x4*1e6:.3f}) "
          f"→ 6D ×{bx/max(x6,1e-30):.1f} ({x6*1e6:.3f}) мкм")
    print(f"вертикаль:   без {by*1e6:.3f} → 4D ×{by/max(y4,1e-30):.1f} → 6D ×{by/max(y6,1e-30):.1f} мкм")
    if x6 < 0.8 * x4:                                          # значимое улучшение, не шум
        print(f"\n[вывод] 6D-модель УЛУЧШАЕТ горизонталь в {x4/max(x6,1e-30):.1f}× — пол опции 3 был зазором "
              f"модель↔реальность (дисперсия/энергия в H); лечится матрицей, снятой с реальной машины (LOCO)")
    else:
        print(f"\n[вывод] 6D-модель горизонталь НЕ улучшила ({x4*1e6:.3f}→{x6*1e6:.3f} мкм, в пределах шума) — "
              f"пол НЕ в отклике. Значит это степень свободы, недостижимая чистым рулением: энергия/дисперсия "
              f"(ВЧ-частота) либо подпространство, срезанное rcond. Следующий рычаг — ВЧ/энергия, не матрица")


if __name__ == "__main__":
    main()
