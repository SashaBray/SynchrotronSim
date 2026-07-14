# -*- coding: utf-8 -*-
"""Подбор тока I1 каждого DQ на СЕКТОРНОЙ модели под ПРОЕКТНЫЙ ИЗГИБ (только ∫b0, I2=I3=0).
Гипотеза: гор.увод орбиты в диффуре = DQ недоворачивает (осевая coupling дала токи ~5А → изгиб
×4 меньше проектного). Здесь I1 solo ставит проектный изгиб при МИНИМАЛЬНОМ секступоле (I2=0),
жертвуя градиентом (focusing уйдёт — потом добираем QP). Цель: ∫b0·ds = −mat_B·L_eff.
Запуск: python dq_bend_solve.py [lattice=101]"""
import os
import sys
import csv
import numpy as np
import dq_sector_profile as dsp

TB = os.path.join("workspace", "tables")


def sector_b(devs, I1, hw):
    prof = dsp.build(list(devs), [I1, 0.0, 0.0], hw)
    return float(np.trapezoid(prof["b0"], prof["Xg"]))


def solve_i1(devs, target, hw, I0, iters=12):
    I = float(I0); r = None
    for _ in range(iters):
        f = sector_b(devs, I, hw); r = f - target
        if abs(r) < 1e-9:
            break
        dI = max(0.2, abs(I) * 0.05)
        J = (sector_b(devs, I + dI, hw) - f) / dI
        if abs(J) < 1e-30:
            break
        I -= np.clip(r / J, -30, 30)
    return I, (r if r is not None else 0.0)


def main(lattice=101):
    hw = {r["device_type_id"]: r for r in csv.DictReader(
        open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";")}
    eff = {d: float(r["eff_length"]) for d, r in hw.items()}
    path = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";")); head = rows[0]
    li, di = head.index("lattice_id"), head.index("device_type_id")
    ai, mb = head.index("arg_val"), head.index("mat_B")
    xi, zi = head.index("x"), head.index("z")
    ring = [r for r in rows[1:] if r and r[li] == str(lattice)]
    i = 0; groups = []
    while i < len(ring):
        r = ring[i]
        if r[di].startswith("DQ_I"):
            g = [r]
            while i + 1 < len(ring) and ring[i + 1][xi:zi + 1] == r[xi:zi + 1]:
                g.append(ring[i + 1]); i += 1
            if len(g) == 3:
                groups.append(g)
        i += 1
    print(f"Решётка {lattice}: DQ-групп {len(groups)}. Подбор I1 под ПРОЕКТНЫЙ ИЗГИБ (I2=I3=0)…")
    cache = {}; n = 0; maxres = 0.0; maxI = 0.0
    for g in groups:
        devs = tuple(x[di] for x in g); Leff = eff[devs[0]]
        matB = float(g[0][mb] or 0); target = -matB * Leff
        key = (devs, round(matB, 10))
        if key in cache:
            I1, res = cache[key]
        else:
            I1, res = solve_i1(devs, target, hw, float(g[0][ai] or 0)); cache[key] = (I1, res)
        g[0][ai] = f"{I1:.9g}"; g[1][ai] = "0"; g[2][ai] = "0"
        n += 1; maxres = max(maxres, abs(res)); maxI = max(maxI, abs(I1))
        if n <= 3:
            b0 = sector_b(devs, I1, hw)
            print(f"  {devs[0]}: цель ∫b0={target:+.6f} → I1={I1:+.3f}А (∫b0={b0:+.6f}, изгиб={b0* -84.925*1e3:+.2f}мрад... )")
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    print(f"[OK] {n} DQ (уник {len(cache)}). Макс невязка {maxres:.2e}, макс |I1| {maxI:.1f}А.")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
