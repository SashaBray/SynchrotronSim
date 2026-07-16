# -*- coding: utf-8 -*-
"""Токи обмоток DQ в кольце — из АРОЧНОГО датасета dq_arc_i1i2.csv ОБРАТНОЙ интерполяцией (НЕ через
матрицу связей). Для каждого DQ берём проектные (θ, G) из lattice_configs (mat_B, mat_G), в плоскости
(θ_дуга, G_дуга) датасета линейно интерполируем ближайшие → (i1, i2); i3=0.

Длина: датасет для базового DQ (L≈0.45); проверено θ∝Leff, G≈const ⇒ цель по базовой длине
  θ_base = θ_проект·(L_base/L_ring),  G_base = G_проект.
Пишем arg_val обмоток DQ (nominal_arg=100 ⇒ arg_val=ток[А]). Остальные элементы — как обычно
phys_to_current. Запуск ПОСЛЕ phys_to_current: py -3.14 dq_currents_from_dataset.py [lattice=101]"""
import os, csv, sys, shutil
import numpy as np
from collections import defaultdict
from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
from bip_matrix import QP
import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
KD = os.path.join("workspace", "research", "kickdoc")
EFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}
LBASE = EFF["DQ_I1_1A"]                                       # длина базового DQ датасета (0.451)


def build_interp():
    """(θ_дуга, G_дуга) → (i1, i2) линейной интерполяцией; вне охвата — ближайший."""
    D = list(csv.DictReader(open(os.path.join(KD, "dq_arc_i1i2.csv"), encoding="utf-8-sig"), delimiter=";"))
    pts = np.array([(float(r["theta_arc_mrad"]), float(r["G_arc_Tlm"])) for r in D])
    i1 = np.array([float(r["i1_A"]) for r in D]); i2 = np.array([float(r["i2_A"]) for r in D])
    L1 = LinearNDInterpolator(pts, i1); L2 = LinearNDInterpolator(pts, i2)
    N1 = NearestNDInterpolator(pts, i1); N2 = NearestNDInterpolator(pts, i2)

    def interp(th, G):
        a = L1(th, G)
        if not np.isfinite(a):
            return float(N1(th, G)), float(N2(th, G)), "nearest"
        return float(a), float(L2(th, G)), "linear"
    return interp, len(D)


def main(lattice=101):
    interp, ndata = build_interp()
    print(f"Датасет dq_arc_i1i2.csv: {ndata} точек. Базовая длина DQ = {LBASE:.3f} м. QP={QP:.3f}")
    path = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(rows[0].keys())
    lc = [r for r in rows if r["lattice_id"] == str(lattice)]

    groups = defaultdict(list)                               # со-локованные обмотки (как build_line)
    for r in lc:
        groups[(r["x"], r["y"], r["z"])].append(r)

    n = 0; nnear = 0; by = defaultdict(list)
    print(f"\n{'DQ-тип':22}{'Leff':>6}{'θпр,мрад':>9}{'Gпр':>9}{'θбаз':>8}{'i1,А':>8}{'i2,А':>8}  метод")
    for g in groups.values():
        i1r = next((x for x in g if x["device_type_id"].startswith("DQ_I1")), None)
        if not i1r:
            continue
        i2r = next((x for x in g if x["device_type_id"].startswith("DQ_I2")), None)
        i3r = next((x for x in g if x["device_type_id"].startswith("DQ_I3")), None)
        dev = i1r["device_type_id"]; Leff = EFF[dev]
        mB = float(i1r.get("mat_B") or 0); mG = float(i1r.get("mat_G") or 0)
        th_des = QP * mB * Leff * 1e3                        # проектный изгиб, мрад
        th_base = th_des * (LBASE / Leff)                   # пересчёт к базовой длине датасета
        a, b, how = interp(th_base, mG)
        i1r["arg_val"] = f"{a:.4f}"
        if i2r: i2r["arg_val"] = f"{b:.4f}"
        if i3r: i3r["arg_val"] = "0"
        n += 1; nnear += (how == "nearest"); by[dev].append((a, b))
        if len(by[dev]) == 1:                                # печать по одному примеру на тип
            print(f"{dev:22}{Leff:>6.2f}{th_des:>9.2f}{mG:>9.4f}{th_base:>8.2f}{a:>8.3f}{b:>8.3f}  {how}")

    shutil.copy(path, path + ".bak_dqds")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows)
    print(f"\n[OK] токи {n} DQ-групп из датасета (i3=0) → lattice_configs {lattice} (бэкап .bak_dqds); "
          f"вне охвата (nearest): {nnear}")
    for dev, v in by.items():
        v = np.array(v)
        print(f"    {dev:22}: {len(v)} шт, i1∈[{v[:,0].min():.2f},{v[:,0].max():.2f}] i2∈[{v[:,1].min():.2f},{v[:,1].max():.2f}] А")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
