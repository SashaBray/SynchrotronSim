# -*- coding: utf-8 -*-
"""Опция 7c: развёртка rcond — truncation подпространства vs нелинейный feed-down.

7a/7b установили: горизонтальный пол FOFB (~6.2 мкм, ×2.9) — НЕ энергия/дисперсия (реальный сдвиг δ
ничтожен) и НЕ фиделити матрицы (6D=4D). Остаётся гипотеза: геометрия корректируемого подпространства —
SVD при rcond=1e-2 выбрасывает слабые моды, в т.ч. крупные гладкие (η_x-подобные), и остаток = то, что
срезано. Прямой тест: снизить rcond (впустить слабые моды) и смотреть горизонтальный остаток.

ДИСКРИМИНАТОР:
  • остаток ПАДАЕТ к вертикальным ×50 → пол был TRUNCATION, моды реальны, модель совпадает с реальностью.
    Рычаг — rcond (и/или регуляризация), дёшево.
  • остаток НЕ падает / пучок дестабилизируется / кики раздуваются → слабые моды в РЕАЛЬНОСТИ не ведут
    себя как в линейной модели (нелинейный feed-down, зазор модель↔реальность на слабых модах). Рычаг —
    нелинейная поправка модели отклика, не rcond.

ЗАПУСК: py -3.14 control/rcond_sweep.py [smoke]
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
N = 50 if SMOKE else 150
TACT_CHUNKS = 9
LAT, NSL = 26, 8
RCONDS = [1e-2, 5e-3, 3e-3] if SMOKE else [2e-2, 1e-2, 5e-3, 3e-3]  # по переходу числа удержанных мод


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def proj_r2(r, basis):
    r = np.asarray(r, float); basis = np.asarray(basis, float)
    bb = float(basis @ basis); rr = float(r @ r)
    if bb == 0.0 or rr == 0.0:
        return 0.0
    d = float(basis @ r) / bb
    return 1.0 - float((r - basis * d) @ (r - basis * d)) / rr


def run(rcond, enabled, etax):
    plant = TbtPlant(lattice=LAT, nsl=NSL, settle_turns=(150 if SMOKE else 300),
                     dist_amp=0.3e-6, meas_noise=1e-7, resp_dims=4)
    dt = plant.stride * plant.TREV
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, dt)
    cx = IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5, rcond=rcond)
    cy = IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=0.5, rcond=rcond)
    sch.add(cx, rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    sch.add(cy, rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    ox, oy = [], []
    for _ in range(N):
        sch.step()
        ox.append(np.asarray(state.observables["ebpm.x"]))
        oy.append(np.asarray(state.observables["ebpm.y"]))
    tail = slice(len(ox) // 2, None)
    rx = np.mean(ox[tail], axis=0); ry = np.mean(oy[tail], axis=0)
    # макс. кик корректора (раздувание на слабых модах) из внутреннего состояния контуров
    kick_x = float(np.max(np.abs(cx.u))) if enabled and hasattr(cx, "u") else 0.0
    kick_y = float(np.max(np.abs(cy.u))) if enabled and hasattr(cy, "u") else 0.0
    alive = plant.alive
    plant.close()
    return dict(rx=rms(rx), ry=rms(ry), r2x=proj_r2(rx, etax), kick_x=kick_x, kick_y=kick_y, alive=alive)


def main():
    print("=== ОПЦИЯ 7c: развёртка rcond — truncation vs нелинейность ===")
    print(f"{N} чанков×15 об; искажение feed-down 0.3 мкм; контур 2 кГц; sweep rcond обе плоскости\n")

    # η_x в BPM + спектр сингулярных чисел модели отклика (сколько мод режет каждый rcond)
    rm = RingModel(lattice=LAT, nsl=NSL)
    etax = ra.twiss(rm.ring)["etax"][rm.mon]
    Hc = rm.response(rm.hk_ids, ["ebpm.x"])["ebpm.x"]
    rm.close()
    s = np.linalg.svd(Hc, compute_uv=False)
    print(f"модель отклика X: {Hc.shape[0]}×{Hc.shape[1]}, сингулярные σ_max={s[0]:.3e} σ_min={s[-1]:.3e}, "
          f"обусловленность {s[0]/s[-1]:.1e}")
    for rc in RCONDS:
        keep = int(np.sum(s >= rc * s[0]))
        print(f"  rcond={rc:.0e}: удержано {keep}/{len(s)} мод (порог σ≥{rc*s[0]:.2e})")
    print()

    base = run(RCONDS[0], enabled=False, etax=etax)
    print(f"база (контур выкл): X {base['rx']*1e6:.2f} мкм, Y {base['ry']*1e6:.2f} мкм\n")

    print(f"{'rcond':>7} | {'остаток X':>9} | {'×X':>5} | {'R²(η_x)':>7} | {'остаток Y':>9} | {'×Y':>5} | "
          f"{'макс кик X':>10} | пучок")
    rows = []
    for rc in RCONDS:
        r = run(rc, enabled=True, etax=etax)
        rows.append((rc, r))
        print(f"{rc:>7.0e} | {r['rx']*1e6:8.3f} | {base['rx']/max(r['rx'],1e-30):4.1f} | {r['r2x']:7.3f} | "
              f"{r['ry']*1e6:8.3f} | {base['ry']/max(r['ry'],1e-30):4.1f} | {r['kick_x']*1e6:9.2f}µ | "
              f"{'жив' if r['alive'] else 'ПОТЕРЯН'}")

    # --- вывод: два вопроса. (1) выигрыш от де-усечения; (2) остаётся ли X/Y-разрыв при ПОЛНОМ ранге ---
    cur = next(r for rc, r in rows if abs(rc - 1e-2) < 1e-12)    # текущая рабочая точка rcond=1e-2
    full = min(rows, key=lambda kr: kr[0])                        # наименьший rcond = полный ранг
    detrunc = cur['rx'] / max(full[1]['rx'], 1e-30)              # выигрыш от впуска срезанных мод
    gap = full[1]['rx'] / max(full[1]['ry'], 1e-30)             # остаточный разрыв X/Y при полном ранге
    blow = full[1]['kick_x'] > 10 * cur['kick_x'] or not full[1]['alive']
    print()
    print(f"(1) де-усечение rcond 1e-2→{full[0]:.0e} (все {len(s)} мод): горизонталь {cur['rx']*1e6:.2f}→"
          f"{full[1]['rx']*1e6:.2f} мкм (×{detrunc:.1f}), R²(η_x) {cur['r2x']:.2f}→{full[1]['r2x']:.2f}, "
          f"макс кик {cur['kick_x']*1e6:.1f}→{full[1]['kick_x']*1e6:.1f} мкрад, пучок "
          f"{'жив' if full[1]['alive'] else 'ПОТЕРЯН'}.")
    if detrunc > 1.3 and not blow:
        print(f"    → rcond=1e-2 ПЕРЕсекал (обусловленность всего {s[0]/s[-1]:.0f}, все моды здоровы). "
              f"Дешёвый выигрыш ×{detrunc:.1f} простой сменой rcond≈{full[0]:.0e}.")
    elif blow:
        print(f"    → низкий rcond раздувает кики/губит пучок — слабые моды в реальности не как в модели; "
              f"держать rcond≈1e-2.")
    else:
        print(f"    → де-усечение почти не помогает (×{detrunc:.1f}) — 1e-2 уже около оптимума.")
    print(f"(2) при ПОЛНОМ ранге горизонталь {full[1]['rx']*1e6:.2f} мкм против вертикали "
          f"{full[1]['ry']*1e6:.2f} мкм — остаточный разрыв ×{gap:.1f}.")
    if gap > 1.8:
        print(f"    → разрыв НЕ закрывается модами: горизонтальный остаток — структурный предел ЛИНЕЙНОЙ "
              f"дипольной коррекции нелинейной машины (feed-down сквозь секступоли в H, которого нет в модели "
              f"отклика). Рычаг на остаток — НЕЛИНЕЙНАЯ ПОПРАВКА модели отклика, не rcond.")
    else:
        print(f"    → при полном ранге горизонталь догнала вертикаль: пол был ЦЕЛИКОМ усечением мод. "
              f"Достаточно рабочего rcond≈{full[0]:.0e}.")
    # (3) плоскостная асимметрия: оптимальный rcond у X и Y разный (фиделити слабых мод разная)
    bx = min(rows, key=lambda kr: kr[1]['rx'] if kr[1]['alive'] else 1e9)
    by = min(rows, key=lambda kr: kr[1]['ry'] if kr[1]['alive'] else 1e9)
    if abs(bx[0] - by[0]) > 1e-12:
        print(f"(3) АСИММЕТРИЯ ПЛОСКОСТЕЙ: оптимум X при rcond={bx[0]:.0e} ({bx[1]['rx']*1e6:.2f} мкм), "
              f"Y при rcond={by[0]:.0e} ({by[1]['ry']*1e6:.3f} мкм). Низкий rcond помогает X, но ПОРТИТ Y "
              f"(слабые верт. моды возбуждают нелинейное движение) → рабочая точка ПОПЛОСКОСТНАЯ: X≈{bx[0]:.0e}, "
              f"Y≈{by[0]:.0e}, а не один общий rcond.")


if __name__ == "__main__":
    main()
