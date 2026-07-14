# -*- coding: utf-8 -*-
"""Построить ПОЛНОСТЬЮ ЗАМКНУТОЕ кольцо = 40× стандартный супериод (каждый гнёт ровно 9°, 40×9°=360°).

Реплицирует геометрию + НАСТРОЕННЫЕ mat_* стандартного супериода жёстким трансформом T (поворот 9° +
сдвиг, извлечён из реш. 107: SP0→SP1). Пишет решётку 140 в lattice_configs + connectivity (замкнутую) +
эксперимент 41. Прогоняет пучок, проверяет ПОТЕРИ и ЗАМЫКАНИЕ опорной. Это идеализация БЕЗ junction-
секций (их коэф. не настраивались; реальное кольцо с junction — реш. 101).

Запуск:  py -3.14 build_closed_ring.py
Затем:   py -3.14 matrix_traj_export.py 41 140   (траектория для визуализатора)
         py -3.14 visualize.py 41
"""
import os, csv, math, shutil
import numpy as np
import matrix_tracker as mt
import matrix_optimizer as mo

TB = mt.TB
NSP, LAT, EXP = 40, 140, 41


def rot(th):
    c, s = math.cos(th), math.sin(th)
    return np.array([[c, -s], [s, c]])


def _transform():
    """T (R, t) из реш. 107: SP0[0]→SP1[0], поворот РОВНО 2π/40 = 9°."""
    allrows = list(csv.DictReader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    c107 = sorted([r for r in allrows if r["lattice_id"] == "107"], key=lambda r: int(r["instance_id"]))
    dyaw = 2 * math.pi / NSP
    R = rot(dyaw)
    p0 = np.array([float(c107[0]["x"]), float(c107[0]["y"])])
    p1 = np.array([float(c107[69]["x"]), float(c107[69]["y"])])
    return R, p1 - R @ p0, dyaw, c107[:69]


def build_rows():
    path = os.path.join(TB, "lattice_configs.csv")
    allrows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(allrows[0].keys())
    R, t, dyaw, base = _transform()
    allrows = [r for r in allrows if r["lattice_id"] != str(LAT)]           # убрать старую 140
    xy = [np.array([float(r["x"]), float(r["y"])]) for r in base]
    yaw = [float(r["yaw"]) for r in base]
    new = []; iid = 1
    for k in range(NSP):
        for i, r in enumerate(base):
            nr = dict(r); nr["lattice_id"] = str(LAT); nr["instance_id"] = str(iid); iid += 1
            nr["x"] = f"{xy[i][0]:.6f}"; nr["y"] = f"{xy[i][1]:.6f}"; nr["yaw"] = f"{yaw[i]:.9f}"
            new.append(nr)
        xy = [R @ p + t for p in xy]; yaw = [y + dyaw for y in yaw]
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(path, os.path.join(TB, "backups", "lattice_configs_pre_closedring.csv"))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(allrows + new)
    return len(new)


def build_connectivity():
    path = os.path.join(TB, "element_connectivity.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(rows[0].keys())
    rows = [r for r in rows if r["lattice_id"] != str(LAT)]
    cfg = sorted([r for r in csv.DictReader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"),
                                            delimiter=";") if r["lattice_id"] == str(LAT)],
                 key=lambda r: int(r["instance_id"]))
    n = len(cfg); new = []
    for i, r in enumerate(cfg):
        new.append({"lattice_id": str(LAT), "instance_id": r["instance_id"], "device_type_id": r["device_type_id"],
                    "prev_instance_id": cfg[(i - 1) % n]["instance_id"],       # замыкание в кольцо
                    "next_instance_id": cfg[(i + 1) % n]["instance_id"]})
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows + new)


def add_experiment():
    path = os.path.join(TB, "experiments.csv")
    rows = [r for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";") if r["exp_id"] != str(EXP)]
    fields = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";").fieldnames)
    tmpl = next(r for r in rows if r["lattice_id"] == "107").copy()
    tmpl["exp_id"] = str(EXP); tmpl["lattice_id"] = str(LAT)
    rows.append(tmpl); rows.sort(key=lambda r: int(r["exp_id"]))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows)


def main():
    print(f"строю замкнутое кольцо: {NSP}× стандартный супериод (9° каждый)")
    print(f"  решётка {LAT}: {build_rows()} инстансов")
    build_connectivity(); print("  connectivity записана (замкнута)")
    add_experiment(); print(f"  эксперимент {EXP} → решётка {LAT}")
    mags = mo.build_ring(LAT); total = len(mags)
    print(f"  магнитов (сгруппировано): {total} = {NSP}×{total // NSP}")
    # выживание
    P = mo.make_beam(radii=(1e-3, 2e-3))
    L, info = mo.loss_of(mags, P, 0.03, 5)
    print(f"\nТРЕК пучка {P.shape[0]} частиц (1/2мм) через ПОЛНОЕ кольцо:")
    print(f"  горизонт={info['horizon']}/{total}  потеряно={info['nlost']}/{P.shape[0]}  "
          f"макс.смещ={info['maxoff_surv']*1e3:.1f}мм")
    # замыкание опорной
    import matrix_traj_export as mte
    ref = mte.track_global(mags, (0, 0, 0, 0, 0, 0))
    gap = float(np.linalg.norm(ref[-1] - ref[0])) if len(ref) else -1
    print(f"  опорная: {len(ref)} точек, |замыкание старт↔финиш| = {gap*1e3:.1f}мм "
          f"(старт={ref[0].round(3)}, финиш={ref[-1].round(3)})")


if __name__ == "__main__":
    main()
