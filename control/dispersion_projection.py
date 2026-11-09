# -*- coding: utf-8 -*-
"""Опция 7a: дисперсионная проекция горизонтального пола ОС.

Опция 7 опровергла «зазор модель↔реальность»: 6D-матрица отклика (LOCO-подобная) дала тот же
горизонтальный остаток ×2.9, что и 4D. Значит горизонталь упирается в степень свободы, НЕдостижимую
рулением диполями. Гипотеза: это энергия/дисперсия — горизонтальное искажение сдвигает равновесную
энергию δ (через ВЧ и длину орбиты), и остаётся дисперсионная орбита η_x·δ; корректоры-диполи обнуляют
орбиту в своём подпространстве, но не могут снять дисперсионную компоненту. Вертикаль дисперсии не имеет
(η_y≈0) → чистые ×50.

ТЕСТ (однозначный): спроецировать горизонтальную орбиту на η_x в точках BPM.
  δ̂ = (η·r)/(η·η);  R² = 1 − ||r − η·δ̂||²/||r||²  — доля орбиты, объяснённая ОДНИМ скаляром энергии.
Сравниваем ДВА состояния на одной машине:
  • БАЗА (контур выкл): орбита — смесь корректируемого и дисперсионного → R² невысокий;
  • ОСТАТОК (FOFB): контур снял всё, что достаёт → если остаток коллапсирует на η_x (R²→1), пол
    ЭНЕРГЕТИЧЕСКИЙ, и рычаг — ВЧ-частота, а не добавление корректоров.
Контроль: вертикаль на η_y (должна не проецироваться — там дисперсии нет).

ЗАПУСК: py -3.14 control/dispersion_projection.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import ring_analysis as ra

from control import ExperimentState, Limits, Scheduler
from control.tbt_plant import TbtPlant
from control.ring_model import RingModel
from control.controllers.integral_svd import IntegralSVD

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
N = 60 if SMOKE else 180
TACT_CHUNKS = 9
LAT, NSL = 26, 8


def project(r, basis):
    """δ̂, R² и косинус угла: сколько вектора r объясняет один масштаб базиса (η·δ)."""
    r = np.asarray(r, float); basis = np.asarray(basis, float)
    bb = float(basis @ basis)
    if bb == 0.0 or float(r @ r) == 0.0:
        return 0.0, 0.0, 0.0
    d = float(basis @ r) / bb
    resid = r - basis * d
    r2 = 1.0 - float(resid @ resid) / float(r @ r)
    cos = float(basis @ r) / np.sqrt(bb * float(r @ r))
    return d, r2, cos


def orbit_tail(enabled):
    """Средняя орбита (x,y) в BPM на хвосте: контур выкл → база, вкл → остаток FOFB. Одна машина."""
    plant = TbtPlant(lattice=LAT, nsl=NSL, settle_turns=(150 if SMOKE else 300),
                     dist_amp=0.3e-6, meas_noise=1e-7, resp_dims=4)
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
    return np.mean(ox[tail], axis=0), np.mean(oy[tail], axis=0), alive


def main():
    print("=== ОПЦИЯ 7a: дисперсионная проекция горизонтального пола ===")
    print(f"{N} чанков×15 об; одна машина (feed-down 0.3 мкм), контур 2 кГц; проекция орбиты BPM на η\n")

    # --- дисперсия η в точках BPM (номинальная решётка, 4D off-energy орбита) ---
    rm = RingModel(lattice=LAT, nsl=NSL)
    tw = ra.twiss(rm.ring)
    ex = tw["etax"][rm.mon]; ey = tw["etay"][rm.mon]
    rm.close()
    print(f"η_x в BPM: СКО {np.std(ex)*1e3:.1f} мм, размах {ex.min()*1e3:.1f}…{ex.max()*1e3:.1f} мм;  "
          f"η_y: СКО {np.std(ey)*1e3:.2f} мм (контроль — должна быть ~0)\n")

    bx, by, _ = orbit_tail(enabled=False)
    rx, ry, a = orbit_tail(enabled=True)

    dbx, r2bx, _ = project(bx, ex)                            # база, горизонталь на η_x
    drx, r2rx, cosrx = project(rx, ex)                        # остаток FOFB на η_x
    _, r2ry, _ = project(ry, ey)                              # контроль: остаток верт. на η_y

    def rms(v): return float(np.sqrt(np.mean(np.asarray(v) ** 2)))
    print(f"{'состояние':24} | {'rms X':>8} | {'R²(η_x)':>8} | δ̂ (энергия)")
    print(f"{'база (контур выкл)':24} | {rms(bx)*1e6:7.2f} | {r2bx:8.3f} | {dbx:+.2e}")
    print(f"{'остаток FOFB':24} | {rms(rx)*1e6:7.2f} | {r2rx:8.3f} | {drx:+.2e}  (cos={cosrx:+.2f})")
    print(f"\nвертикаль-контроль: остаток rms {rms(ry)*1e6:.3f} мкм, R²(η_y) {r2ry:.3f} (пучок {'жив' if a else 'ПОТЕРЯН'})")

    # δ̂ совпадает в базе и остатке? → энергосдвиг НЕ снимается контуром (ключевой тест, не величина R²)
    d_match = dbx != 0.0 and abs(drx - dbx) / abs(dbx) < 0.25 and np.sign(drx) == np.sign(dbx)
    amp_frac = abs(cosrx)                                     # доля остатка по амплитуде вдоль η_x (√R²)
    print()
    if r2rx > 0.7 and r2rx > r2bx + 0.2:
        print(f"[вывод] ПОДТВЕРЖДЕНО: остаток FOFB коллапсирует на дисперсию (R² {r2bx:.2f}→{r2rx:.2f}). "
              f"Горизонтальный пол — ЭНЕРГЕТИЧЕСКАЯ степень свободы (η_x·δ, δ̂≈{drx:+.1e}), диполями не "
              f"снимается. Рычаг — ВЧ-частота/энергия как актуатор (опция 7b), не добавление корректоров.")
    elif r2rx > r2bx + 0.2 and d_match:
        print(f"[вывод] ПОДТВЕРЖДЕНО ЧАСТИЧНО (смешанный пол). Энергия/дисперсия — РЕАЛЬНАЯ несжимаемая часть: "
              f"δ̂ базы {dbx:+.2e} и остатка {drx:+.2e} совпадают (контур её не снимает), вдоль η_x лежит "
              f"{amp_frac*100:.0f}% амплитуды остатка (R² {r2bx:.2f}→{r2rx:.2f}). НО не весь пол: ортогональная "
              f"компонента {np.sqrt(max(0,1-r2rx))*100:.0f}% — подпространство rcond=1e-2 либо нелинейный feed-down.\n"
              f"        ДВА рычага: (7b) ВЧ-частота/энергия снимет дисперсионные ~{amp_frac*100:.0f}%; "
              f"снижение rcond/нелинейная поправка — ортогональный остаток.")
    else:
        print(f"[вывод] ОПРОВЕРГНУТО: остаток на η_x НЕ проецируется (R² {r2rx:.2f}, δ̂ база/остаток "
              f"{'совпали' if d_match else 'разошлись'}). Пол НЕ энергетический — подпространство rcond=1e-2 "
              f"либо нелинейный feed-down. ВЧ-рычаг не ключ; снизить rcond и/или нелинейная поправка модели.")


if __name__ == "__main__":
    main()
