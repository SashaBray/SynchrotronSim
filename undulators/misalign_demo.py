# -*- coding: utf-8 -*-
"""Демо ЖЁСТКОГО перекоса ондулятора: всё устройство (полюса+фокус) смещается/поворачивается как целое.

Показывает, что система хранения/импорта это поддерживает: берём группу из element_groups.csv, применяем
один 6-мерный сдвиг (dx,dy,dz,dyaw,dpitch,droll) ко ВСЕМ её инстансам вокруг pivot (errors.apply_group_rigid),
пишем возмущённую рабочую папку, трекаем. Сравниваем тюн и замкнутую орбиту с номиналом.

ЗАПУСК: py -3.14 undulators/misalign_demo.py --lattice 26 --group und_19 --dy 200e-6 --droll 1e-3
"""
import os, sys, csv, shutil, argparse, tempfile
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import ring_cpp as rc
import ring_analysis as ra
import errors as er

TB = os.path.join("workspace", "tables")


def perturbed_group_ws(lattice, group_id, dof, out_dir):
    """Рабочая папка с ЖЁСТКО перекошенной группой group_id."""
    groups = {g["group_id"]: g for g in er.load_groups(lattice)}
    if group_id not in groups:
        raise SystemExit(f"группа {group_id} не найдена (есть: {sorted(groups)[:5]}…)")
    dst = os.path.join(out_dir, "tables"); os.makedirs(dst, exist_ok=True)
    rows = list(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    hdr = rows[0]; ci = {c: i for i, c in enumerate(hdr)}
    er.apply_group_rigid(rows[1:], ci, groups[group_id], dof)
    with open(os.path.join(dst, "lattice_configs.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(rows[1:])
    for t in ("element_connectivity.csv", "particle_groups.csv"):
        shutil.copyfile(os.path.join(TB, t), os.path.join(dst, t))
    return out_dir


def state(ring):
    tw = ra.twiss(ring)
    orb = np.asarray(tw["orbit"])                      # замкнутая орбита на выходах элементов (x,px,y,py,…)
    return tw["nux"], tw["nuy"], np.sqrt(np.mean(orb[:, 0] ** 2)), np.sqrt(np.mean(orb[:, 2] ** 2))


def main():
    ap = argparse.ArgumentParser(description="демо жёсткого перекоса ондулятора")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--group", default="und_19")
    ap.add_argument("--dx", type=float, default=0.0); ap.add_argument("--dy", type=float, default=200e-6)
    ap.add_argument("--dz", type=float, default=0.0); ap.add_argument("--dyaw", type=float, default=0.0)
    ap.add_argument("--dpitch", type=float, default=0.0); ap.add_argument("--droll", type=float, default=0.0)
    a = ap.parse_args()
    dof = dict(dx=a.dx, dy=a.dy, dz=a.dz, dyaw=a.dyaw, dpitch=a.dpitch, droll=a.droll)

    r0 = rc.CppRing(lattice=a.lattice); nx0, ny0, xr0, yr0 = state(r0); r0.close()
    tmp = tempfile.mkdtemp(prefix="und_mis_")
    try:
        perturbed_group_ws(a.lattice, a.group, dof, tmp)
        r1 = rc.CppRing(lattice=a.lattice, ws=tmp); nx1, ny1, xr1, yr1 = state(r1); r1.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print(f"  решётка {a.lattice}, группа {a.group}")
    print(f"  перекос: dx={a.dx*1e6:.0f} dy={a.dy*1e6:.0f} dz={a.dz*1e6:.0f} мкм; "
          f"dyaw={a.dyaw*1e3:.2f} dpitch={a.dpitch*1e3:.2f} droll={a.droll*1e3:.2f} мрад")
    print(f"  номинал:  νx={nx0:.5f} νy={ny0:.5f}  орбита СКО x={xr0*1e6:.3f} y={yr0*1e6:.3f} мкм")
    print(f"  перекос:  νx={nx1:.5f} νy={ny1:.5f}  орбита СКО x={xr1*1e6:.3f} y={yr1*1e6:.3f} мкм")
    print(f"  Δνx={nx1-nx0:+.6f} Δνy={ny1-ny0:+.6f};  Δорбита x={(xr1-xr0)*1e6:+.3f} y={(yr1-yr0)*1e6:+.3f} мкм")


if __name__ == "__main__":
    main()
