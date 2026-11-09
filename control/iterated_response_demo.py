# -*- coding: utf-8 -*-
"""Опция 7d: итеративная (нелинейная) модель отклика против остаточного горизонтального пола.

7c разложил горизонтальный пол: ×1.6 было лишним усечением rcond (снято), остаток ×2.6 к вертикали —
структурный предел ЛИНЕЙНОЙ дипольной коррекции. Механизм: H_nominal снята на ПРОЕКТНОЙ орбите (≈0), а
скорректированный пучок идёт со смещением сквозь секступоли (feed-down: секступоль при x работает как
квадруполь ∝S·x) → реальный отклик корректор→орбита ≠ H_nominal → FOFB по H_nominal оставляет остаток.

ТЕСТ (детерминированный, на замкнутой орбите — без бетатронного шума; feed-down чисто поперечный, dims=4):
  A (линейная модель, как живой FOFB): итерация θ −= pinv(H_nom)·орбита, H_nom снята на ЧИСТОМ кольце и
    ФИКСИРОВАНА. Сходится к тому же полу, что живой контур.
  B (Гаусс-Ньютон, нелинейная модель): та же итерация, но H ПЕРЕСНИМАЕТСЯ на текущей рабочей орбите на
    каждом шаге (feed-down секступолей входит в H естественно).
Если B сбивает горизонталь к вертикальному уровню, а A застревает — остаток ИМЕННО feed-down, и рычаг =
итеративная/нелинейная модель отклика. Поплоскостной rcond из 7c (X 3e-3, Y 1e-2).

ЗАПУСК: py -3.14 control/iterated_response_demo.py
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import ring_cpp
import facility_report as fr
import ring_analysis as ra

from control.ring_model import RingModel

LAT, NSL = 26, 8
PHASE_DEG = 165.0
DIST_AMP = 0.3e-6
DIST_SEED = 11
EPS = 1e-6
RCOND = {"x": 3e-3, "y": 1e-2}          # поплоскостной оптимум из 7c
N_ITER_A = 8                             # фикс-H итерация до сходимости
N_ITER_B = 4                             # Гаусс-Ньютон (пересъём H)


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    model = RingModel(lattice=LAT, nsl=NSL)
    mon = model.mon; corr = np.asarray(model.corr, int); nc = len(corr)
    quads = model.quads; K1L = np.asarray(model.K1L, float); nq = len(quads)
    etax = ra.twiss(model.ring)["etax"][mon]
    model.close()

    rp = np.random.default_rng(DIST_SEED)
    patx = rp.normal(0, 1, nq); patx *= 1.0 / np.std(patx)
    paty = rp.normal(0, 1, nq); paty *= 1.0 / np.std(paty)
    dkx = K1L * (patx * DIST_AMP); dky = -K1L * (paty * DIST_AMP)

    ring = fr._ring(inf, nsl=NSL)        # ТО ЖЕ кольцо, что у RingModel/живого FOFB: СИ/ВЧ выкл, полная
    #                                      нелинейная решётка (секступоли есть → feed-down в замкнутой орбите)

    def orbit_at(thx, thy, dist=True):
        """Замкнутая орбита (dims=4) при корректорах thx/thy (+ feed-down искажение, если dist)."""
        kk = [(int(corr[j]), float(thx[j])) for j in range(nc)]
        vv = [(int(corr[j]), float(thy[j])) for j in range(nc)]
        if dist:
            kk += [(int(quads[j]), float(dkx[j])) for j in range(nq)]
            vv += [(int(quads[j]), float(dky[j])) for j in range(nq)]
        ring.apply_perturb(kicks=kk, vkicks=vv)
        co = ra.closed_orbit(ring, dims=4)
        return ring.track_elements([co])[0]

    def resp_at(thx, thy, plane, dist=True):
        """Якобиан корректор→орбита в плоскости plane вокруг рабочей точки (thx,thy). feed-down входит,
        если орбита смещена (dist=True и/или корректоры ненулевые)."""
        col = 0 if plane == "x" else 2
        base = orbit_at(thx, thy, dist)[mon, col]
        H = np.zeros((len(mon), nc))
        for j in range(nc):
            tx, ty = thx.copy(), thy.copy()
            if plane == "x":
                tx[j] += EPS
            else:
                ty[j] += EPS
            H[:, j] = (orbit_at(tx, ty, dist)[mon, col] - base) / EPS
        return H

    thx0 = np.zeros(nc); thy0 = np.zeros(nc)
    base = orbit_at(thx0, thy0, dist=True)
    print(f"корректоров {nc}, мониторов {len(mon)}; база (искажение, корректоры=0): "
          f"X {rms(base[mon,0])*1e6:.2f} мкм, Y {rms(base[mon,2])*1e6:.2f} мкм\n")

    # --- H_nominal: снята на ЧИСТОМ кольце (как живой FOFB) ---
    print("снимаю H_nominal (чистое кольцо)…", flush=True)
    Hx_nom = resp_at(thx0, thy0, "x", dist=False)
    Hy_nom = resp_at(thx0, thy0, "y", dist=False)
    Px_nom = np.linalg.pinv(Hx_nom, rcond=RCOND["x"])
    Py_nom = np.linalg.pinv(Hy_nom, rcond=RCOND["y"])

    # === A: фикс-H итерация (линейная модель) ===
    print("A: итерация с ФИКСИРОВАННОЙ H_nominal…", flush=True)
    thx, thy = thx0.copy(), thy0.copy()
    for _ in range(N_ITER_A):
        o = orbit_at(thx, thy)
        thx -= Px_nom @ o[mon, 0]
        thy -= Py_nom @ o[mon, 2]
    oA = orbit_at(thx, thy)
    rAx, rAy = rms(oA[mon, 0]), rms(oA[mon, 2])

    # === B: Гаусс-Ньютон (H переснимается на рабочей орбите) ===
    print("B: Гаусс-Ньютон (пересъём H на рабочей орбите каждый шаг)…", flush=True)
    thx, thy = thx0.copy(), thy0.copy()
    for it in range(N_ITER_B):
        o = orbit_at(thx, thy)
        Hx = resp_at(thx, thy, "x", dist=True)          # feed-down на ТЕКУЩЕЙ орбите
        Hy = resp_at(thx, thy, "y", dist=True)
        thx -= np.linalg.pinv(Hx, rcond=RCOND["x"]) @ o[mon, 0]
        thy -= np.linalg.pinv(Hy, rcond=RCOND["y"]) @ o[mon, 2]
        oB = orbit_at(thx, thy)
        print(f"   итер {it+1}: X {rms(oB[mon,0])*1e6:.3f}  Y {rms(oB[mon,2])*1e6:.3f} мкм", flush=True)
    rBx, rBy = rms(oB[mon, 0]), rms(oB[mon, 2])
    ring.apply_perturb(); ring.close()

    print(f"\n{'метод':34} | {'остаток X':>9} | {'остаток Y':>9}")
    print(f"{'A: линейная модель (H_nom фикс)':34} | {rAx*1e6:8.3f} | {rAy*1e6:8.3f}")
    print(f"{'B: Гаусс-Ньютон (пересъём H)':34} | {rBx*1e6:8.3f} | {rBy*1e6:8.3f}")
    print(f"\nгоризонталь: A {rAx*1e6:.2f} → B {rBx*1e6:.2f} мкм "
          f"(×{rAx/max(rBx,1e-30):.1f}); вертикаль A {rAy*1e6:.3f} → B {rBy*1e6:.3f} мкм")

    gap_A = rAx / max(rAy, 1e-30); gap_B = rBx / max(rBy, 1e-30)
    print()
    if rBx < 0.6 * rAx:
        print(f"[вывод] ПОДТВЕРЖДЕНО: пересъём H на рабочей орбите сбивает горизонталь A {rAx*1e6:.2f}→"
              f"B {rBx*1e6:.2f} мкм (×{rAx/rBx:.1f}), разрыв X/Y {gap_A:.1f}→{gap_B:.1f}. Остаточный пол 7c "
              f"= feed-down секступолей, отсутствующий в H_nominal. Рычаг: итеративная/нелинейная модель "
              f"отклика (Гаусс-Ньютон / пересъём H на скорректированной орбите, аналог LOCO). Встроить в "
              f"Resources как переснимаемую матрицу → живой контур добьёт горизонталь.")
    elif rBx < 0.9 * rAx:
        print(f"[вывод] ЧАСТИЧНО: пересъём H помогает (A {rAx*1e6:.2f}→B {rBx*1e6:.2f} мкм, ×{rAx/rBx:.1f}), но "
              f"разрыв X/Y не закрыт ({gap_B:.1f}). Часть остатка — feed-down (берётся итерацией), часть — "
              f"высший порядок/связь, не линеаризуемый одной матрицей. Рычаг частичный.")
    else:
        print(f"[вывод] ОПРОВЕРГНУТО: пересъём H горизонталь не сбил (A {rAx*1e6:.2f}→B {rBx*1e6:.2f} мкм). "
              f"Остаток НЕ feed-down первого порядка — это высший порядок/нелинейная динамика, не берётся "
              f"якобианом корректор→орбита. Дипольными корректорами горизонталь на этом уровне предельна.")


if __name__ == "__main__":
    main()
