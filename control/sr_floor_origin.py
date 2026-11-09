# -*- coding: utf-8 -*-
"""Опция 7e/разбор: 7e ОПРОВЕРГЛА 7d — 6D-пол 3.84 мкм НЕ дисперсионная орбита (проекция на η_x лишь
1.36 мкм, 3.59 мкм ⊥η_x). Тогда ЧТО это? Между 4D-полом (1.05) и 6D-полом (3.84) два связанных отличия:
  • радиация (SR вкл/выкл),
  • продольная степень свободы 6D (dims 4→6, ВЧ-резонатор),
и ещё вопрос — матрица отклика: живой FOFB берёт H с 4D-модели (SR-выкл), а плант 6D+SR.

РАЗДЕЛЯЕМ (одно и то же искажение, 10 итераций коррекции):
  H_4D — отклик с ЧИСТОГО 4D-кольца (SR-выкл), как у живого FOFB;
  H_6D — отклик с ЧИСТОГО 6D-кольца (SR-вкл, dims=6), «истинный» плант.

  (1) SR-выкл, dims=4, H_4D         → 1.05 (эталон)
  (2) SR-выкл, dims=6 (ВЧ, без рад.), H_4D → энергетическая СС без радиации
  (3) SR-вкл,  dims=6, H_4D         → 3.84 (= живой FOFB)
  (4) SR-вкл,  dims=6, H_6D         → перфектная модель: если →1.05, пол = МОДЕЛЬНОЕ рассогласование
                                       H_4D↔6D-плант; если 3.84 — пол ФИЗИЧЕН и неустраним диполями.

ЗАПУСК: py -3.14 control/sr_floor_origin.py
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
RCOND = {"x": 3e-3, "y": 1e-2}
N_ITER = 10


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    volt, harm = inf["rf"]; energy = inf["energy_mev"]

    model = RingModel(lattice=LAT, nsl=NSL)
    mon = model.mon; corr = np.asarray(model.corr, int); nc = len(corr)
    quads = model.quads; K1L = np.asarray(model.K1L, float); nq = len(quads)

    rp = np.random.default_rng(DIST_SEED)
    patx = rp.normal(0, 1, nq); patx *= 1.0 / np.std(patx)
    paty = rp.normal(0, 1, nq); paty *= 1.0 / np.std(paty)
    dkx = K1L * (patx * DIST_AMP); dky = -K1L * (paty * DIST_AMP)
    z = np.zeros(nc)

    def kicks(ring, thx, thy, dist):
        kk = [(int(corr[j]), float(thx[j])) for j in range(nc)]
        vv = [(int(corr[j]), float(thy[j])) for j in range(nc)]
        if dist:
            kk += [(int(quads[j]), float(dkx[j])) for j in range(nq)]
            vv += [(int(quads[j]), float(dky[j])) for j in range(nq)]
        ring.apply_perturb(kicks=kk, vkicks=vv)

    def orbit(ring, thx, thy, dims, dist):
        kicks(ring, thx, thy, dist)
        return ring.track_elements([ra.closed_orbit(ring, dims=dims)])[0]

    def build_H(ring, dims):
        b = orbit(ring, z, z, dims, dist=False)
        bx, by = b[mon, 0], b[mon, 2]
        Hx = np.zeros((len(mon), nc)); Hy = np.zeros((len(mon), nc))
        for j in range(nc):
            e = z.copy(); e[j] = EPS
            Hx[:, j] = (orbit(ring, e, z, dims, dist=False)[mon, 0] - bx) / EPS
            Hy[:, j] = (orbit(ring, z, e, dims, dist=False)[mon, 2] - by) / EPS
        return np.linalg.pinv(Hx, rcond=RCOND["x"]), np.linalg.pinv(Hy, rcond=RCOND["y"])

    def correct(ring, dims, Px, Py):
        thx, thy = z.copy(), z.copy()
        for _ in range(N_ITER):
            kicks(ring, thx, thy, dist=True)
            o = ring.track_elements([ra.closed_orbit(ring, dims=dims)])[0]
            thx = thx - Px @ o[mon, 0]
            thy = thy - Py @ o[mon, 2]
        o = orbit(ring, thx, thy, dims, dist=True)
        return rms(o[mon, 0]), rms(o[mon, 2])

    # H_4D на чистом 4D-кольце (SR-выкл) — как живой FOFB
    r4clean = model.ring
    Px4, Py4 = build_H(r4clean, dims=4)
    model.close()

    # H_6D на чистом 6D-кольце (SR-вкл)
    def new6(sr):
        return ring_cpp.CppRing(lattice=LAT, sr=sr,
                                rf={"volt_mv": volt, "harmonic": harm, "phase_deg": PHASE_DEG},
                                design_mev=energy, nsl=NSL)
    r6clean = new6(True)
    Px6, Py6 = build_H(r6clean, dims=6)
    r6clean.close()

    # (1) SR-выкл dims=4, H_4D
    r = fr._ring(inf, nsl=NSL); x1, y1 = correct(r, 4, Px4, Py4); r.apply_perturb(); r.close()
    # (2) SR-выкл dims=6 (ВЧ, без радиации), H_4D
    r = new6(False); x2, y2 = correct(r, 6, Px4, Py4); r.apply_perturb(); r.close()
    # (3) SR-вкл dims=6, H_4D  (= живой FOFB)
    r = new6(True);  x3, y3 = correct(r, 6, Px4, Py4); r.apply_perturb(); r.close()
    # (4) SR-вкл dims=6, H_6D  (перфектная модель)
    r = new6(True);  x4, y4 = correct(r, 6, Px6, Py6); r.apply_perturb(); r.close()

    print(f"{'конфигурация коррекции':40} | {'rms X':>7} | {'rms Y':>7}")
    print(f"{'(1) SR-выкл, dims=4,  H_4D (эталон)':40} | {x1*1e6:6.3f} | {y1*1e6:6.3f}")
    print(f"{'(2) SR-выкл, dims=6,  H_4D (ВЧ, без рад.)':40} | {x2*1e6:6.3f} | {y2*1e6:6.3f}")
    print(f"{'(3) SR-вкл,  dims=6,  H_4D (= живой FOFB)':40} | {x3*1e6:6.3f} | {y3*1e6:6.3f}")
    print(f"{'(4) SR-вкл,  dims=6,  H_6D (перфект.мод.)':40} | {x4*1e6:6.3f} | {y4*1e6:6.3f}")

    print("\n[разбор]")
    # радиация vs энергетическая DOF
    if x2 < 0.5 * x3:
        print(f"  • РАДИАЦИЯ: без SR (dims=6, ВЧ) X={x2*1e6:.2f}≈эталон {x1*1e6:.2f} — энергетическая DOF сама "
              f"пол не поднимает; его создаёт СИНХРОТРОННОЕ ИЗЛУЧЕНИЕ (SR вкл → {x3*1e6:.2f}).")
    elif x2 > 0.7 * x3:
        print(f"  • ЭНЕРГО-DOF: даже без радиации 6D-коррекция X={x2*1e6:.2f}≈{x3*1e6:.2f} — пол от продольной "
              f"степени свободы 6D, не от самой радиации.")
    else:
        print(f"  • СМЕШАННО: без SR X={x2*1e6:.2f} (между {x1*1e6:.2f} и {x3*1e6:.2f}) — вклад и DOF, и радиации.")
    # модель vs физика
    if x4 < 0.5 * x3:
        print(f"  • МОДЕЛЬ: перфектная H_6D давит X {x3*1e6:.2f}→{x4*1e6:.2f} — пол был РАССОГЛАСОВАНИЕМ H_4D↔6D-"
              f"плант. Живой FOFB встаёт на 3.84 из-за H с 4D-модели; с 6D-откликом горизонталь корректируема. "
              f"Лечится ВХОДОМ: строить матрицу отклика на 6D+SR кольце.")
    elif x4 > 0.7 * x3:
        print(f"  • ФИЗИКА: даже перфектная H_6D оставляет X={x4*1e6:.2f}≈{x3*1e6:.2f} — пол ФИЗИЧЕН, диполь-"
              f"корректорами неустраним (не вопрос модели). Это радиационно-равновесная орбита вне бетатронного "
              f"подпространства; снимается лишь энергетическим/ВЧ-частотным каналом (в ядре нет, 7b-A).")
    else:
        print(f"  • ЧАСТИЧНО модель: H_6D помогает ({x3*1e6:.2f}→{x4*1e6:.2f}), но не до эталона {x1*1e6:.2f} — "
              f"есть и рассогласование, и физический остаток.")


if __name__ == "__main__":
    main()
