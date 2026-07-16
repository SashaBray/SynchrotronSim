# -*- coding: utf-8 -*-
"""Конфигурация «1 супериод + первый элемент 2-го SP» и трек ОБЕИМИ моделями для визуального сравнения.

СИЛА, решётка 101: 1 супериод = 78 элементов (шаблон standard), первый элемент 2-го SP = инстанс 79.
Создаёт решётку 111 = инстансы 1..79 (открытая линия) и гоняет:
  • КИК (матричная модель, mat_*): инжекция группы 1 → трек через 79 элементов → exp 59;
  • ДИФФУР (C++ SynchrotronTracker по картам поля): финиш-плоскость ЗА инстансом 79 → exp 58.
Обе траектории в глобальных X,Y,Z (results/trajectories_exp{58,59}.bin) для визуализатора.

Запуск: py -3.14 track_sp1.py
Просмотр: py -3.14 viz_sp1.py   (или viz_html.py 58 / 59)"""
import os, csv, struct, shutil, subprocess
import numpy as np
import matrix_tracker as mt
import matrix_traj_export as mte
import beam_injection as bi

TB = os.path.join("workspace", "tables")
RESULTS = os.path.join("workspace", "results")
LAT_FULL, LAT_SUB, GROUP, NEL = 101, 111, 1, 79
DIFF_EXP, KICK_EXP = 58, 59


def make_sublattice():
    """Решётка LAT_SUB = инстансы 1..NEL решётки LAT_FULL (открытая линия, последний next='')."""
    # lattice_configs
    p = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";")); fld = list(rows[0].keys())
    rows = [r for r in rows if r["lattice_id"] != str(LAT_SUB)]                 # снести старую копию
    sub = [dict(r, lattice_id=str(LAT_SUB)) for r in rows
           if r["lattice_id"] == str(LAT_FULL) and 1 <= int(r["instance_id"]) <= NEL]
    if not os.path.exists(p + ".bak_sp1"):
        shutil.copy(p, p + ".bak_sp1")
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fld, delimiter=";"); w.writeheader(); w.writerows(rows + sub)
    # element_connectivity
    cp = os.path.join(TB, "element_connectivity.csv")
    cr = list(csv.DictReader(open(cp, encoding="utf-8"), delimiter=";")); cfld = list(cr[0].keys())
    cr = [r for r in cr if r["lattice_id"] != str(LAT_SUB)]
    csub = []
    for r in cr:
        if r["lattice_id"] == str(LAT_FULL) and 1 <= int(r["instance_id"]) <= NEL:
            nr = dict(r, lattice_id=str(LAT_SUB))
            if int(r["instance_id"]) == NEL:
                nr["next_instance_id"] = ""                                     # обрыв линии
            csub.append(nr)
    with open(cp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=cfld, delimiter=";"); w.writeheader(); w.writerows(cr + csub)
    print(f"[OK] решётка {LAT_SUB} = инстансы 1..{NEL} (SP1 78 эл + первый элемент SP2): {len(sub)} магнитов")


def colocated_devs(lattice, lo, hi):
    """Device первых обмоток co-located групп инстансов [lo..hi] — как их сливает build_line."""
    conf = {(r["lattice_id"], r["instance_id"]): r for r in mt._read("lattice_configs.csv")}
    conn = sorted([r for r in mt._read("element_connectivity.csv") if r["lattice_id"] == str(lattice)],
                  key=lambda r: int(r["instance_id"]))
    conn = [r for r in conn if lo <= int(r["instance_id"]) <= hi]
    devs = []; i = 0
    while i < len(conn):
        r = conf[(str(lattice), conn[i]["instance_id"])]; devs.append(r["device_type_id"])
        while i + 1 < len(conn):
            rn = conf[(str(lattice), conn[i + 1]["instance_id"])]
            if (rn["x"], rn["y"], rn["z"]) == (r["x"], r["y"], r["z"]):
                i += 1
            else:
                break
        i += 1
    return devs


def kick_track():
    """КИК: инжекция группы GROUP (source=physical, mat_*), трек через SP1 + первый элемент SP2 → bin."""
    plane, ordered, states, info = bi.inject(LAT_FULL, GROUP, source="physical")
    exp_dev = colocated_devs(LAT_FULL, 1, NEL)                    # слитые магниты инстансов 1..79
    ksub = len(exp_dev)
    seg = ordered[:ksub]
    got_dev = [m.dev for m in seg]
    assert got_dev == exp_dev, f"порядок не совпал!\nждём {exp_dev[:6]}…\nполучили {got_dev[:6]}…"
    print(f"  слитых магнитов в инстансах 1..{NEL}: {ksub} (SP1 + первый элемент SP2)")
    trajs = []
    for st in states:
        s = tuple(st); pts = [mte._g(plane, 0.0, st[0], st[2])]
        try:
            s, pp = mte.drift_global(s, plane, seg[0]); pts += pp
            s, pp = mte.body_global(seg[0], s); pts += pp
            for j in range(1, len(seg)):
                s, pp = mte.drift_global(s, seg[j - 1], seg[j]); pts += pp
                s, pp = mte.body_global(seg[j], s); pts += pp
        except mt.ParticleLost:
            pass
        trajs.append(np.asarray(pts, float))
    os.makedirs(RESULTS, exist_ok=True)
    mte.write_bin(trajs, os.path.join(RESULTS, f"trajectories_exp{KICK_EXP}.bin"))
    mte.ensure_experiment(KICK_EXP, LAT_SUB)
    lens = [len(t) for t in trajs]
    print(f"КИК (матрица, mat_*): {len(trajs)} частиц, точек {min(lens)}..{max(lens)} "
          f"(центр.частица прошла все {len(seg)} магнитов: {'да' if lens[0] > len(seg) else 'ПОТЕРЯ'})")
    return plane, seg, info, trajs


def setup_finish_plane(last, skip=0, stride=10, dt=1e-11):
    """Строка эксперимента DIFF_EXP: финиш-плоскость ЗА выходом last (инстанс 79)."""
    c = last.exit + 0.05 * last.axis                              # 5 см за выходом инстанса 79
    u, v, ax = last.uhat, last.vhat, last.axis
    p1 = c + u + v; p2 = c - u + v; p3 = c - v
    row = [str(DIFF_EXP), str(LAT_SUB), str(GROUP), f"{dt:.2e}".replace("e", "E"), str(skip), str(stride)]
    for pt in (p1, p2, p3):
        row += [f"{pt[0]:.6f}", f"{pt[1]:.6f}", f"{pt[2]:.6f}"]
    row += [f"{ax[0]:.6f}", f"{ax[1]:.6f}", f"{ax[2]:.6f}"]
    path = os.path.join(TB, "experiments.csv")
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    rows = [r for r in rows if not (r and r[0] == str(DIFF_EXP))]
    rows.append(row)
    shutil.copy(path, path + ".bak_sp1")
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    print(f"  финиш-плоскость: центр {c.round(3)}  нормаль {ax.round(3)}  (skip={skip}, stride={stride})")


def diff_track(last):
    """ДИФФУР: C++ SynchrotronTracker по решётке LAT_SUB до финиш-плоскости за инстансом 79."""
    setup_finish_plane(last)
    exe = os.path.join("build", "SynchrotronTracker.exe")
    if not os.path.exists(exe):
        exe = os.path.join("build", "Release", "SynchrotronTracker.exe")
    print(f"  запуск {exe} {DIFF_EXP} …")
    r = subprocess.run([exe, str(DIFF_EXP)], capture_output=True, text=True, encoding="utf-8", errors="replace")
    print("\n".join((r.stdout or "").splitlines()[-8:]))


if __name__ == "__main__":
    make_sublattice()
    plane, seg, info, ktr = kick_track()
    print(f"\nинжекция: '{info['after_elem']}' → '{info['before_elem']}', центр {info['center'].round(3)}")
    print("ДИФФУР:")
    diff_track(seg[-1])
    print(f"\n[ГОТОВО] кик → exp {KICK_EXP}, диффур → exp {DIFF_EXP}. Визуализация: py -3.14 viz_sp1.py")
