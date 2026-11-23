# -*- coding: utf-8 -*-
"""Миграция БД: ЯВНЫЕ старт/финиш-плоскости опыта, пучок в ЛОКАЛЬНЫХ координатах старт-плоскости, таблица
резонаторов установки.

БЫЛО. Старт-плоскость выводилась из пучка (центроид + средняя скорость → дрейф, где опорная орбита её
пересекает), частицы — в глобальных XYZ/скоростях (particle_groups.csv), финиш матричного ядра = старт (кольцо)
или выход последнего элемента (отрезок), а p1..p3/dir в experiments.csv читал только Boris (по номерам столбцов).

СТАЛО.
  experiments.csv   start_x..start_nz (точка = начало локальных координат, нормаль = направление пролёта) и
                    finish_x..finish_nz (пусто ⇒ финиш = старт, кольцо) ВМЕСТО p1x..p3z;dirx..dirz.
  particle_groups   group_id;part_id;m;q;x;xp;y;yp;ct;delta — координаты в рамке старт-плоскости (СИ; xp=dx/ds).
  rf_cavities.csv   cavity_id;lattice_id;x;y;z;nx;ny;nz;volt_share;note — резонаторы УСТАНОВКИ (пусто ⇒ ВЧ как было).

КАК. Для каждого опыта — ДОСЛОВНО прежний алгоритм инжекции (вшит ниже): старт = та же рамка (центр и ось до
нормировки, repr ⇒ ядро получает те же числа — Python бит-в-бит). Финиш отрезка: прежний финиш матричного ядра
(выход последнего элемента), если он лежит на p-плоскости Boris; иначе p-плоскость (демо: x=3,5). Частицы
группы → локальные (x=rel·u, y=rel·v, xp=(d·u)/(d·a), yp=(d·v)/(d·a), δ как раньше, ct=0); проверка: новая
каноника px=(1+δ)·xp/√(1+xp²+yp²) совпадает с прежней (1+δ)·(d·u). Общая группа, дающая разные локальные
координаты у разных опытов, клонируется; группа без опытов — относительно своей плоскости (запись в отчёте).
Резервная копия — workspace/tables/backups/explicit_planes_<время>/. Повторный запуск — отказ (уже мигрировано).

ЗАПУСК
  py -3.14 migrate_explicit_planes.py --dry-run     # отчёт без записи
  py -3.14 migrate_explicit_planes.py
"""
import os, sys, csv, math, time, shutil, argparse
sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT); os.chdir(ROOT)
import contextlib, io
import numpy as np
import matrix_tracker as mt
import bip_matrix as bm

TB = os.path.join("workspace", "tables")
EXP = os.path.join(TB, "experiments.csv")
PG = os.path.join(TB, "particle_groups.csv")
CAV = os.path.join(TB, "rf_cavities.csv")
OLD_PLANE = ["p1x", "p1y", "p1z", "p2x", "p2y", "p2z", "p3x", "p3y", "p3z", "dirx", "diry", "dirz"]
NEW_PLANE = ["start_x", "start_y", "start_z", "start_nx", "start_ny", "start_nz",
             "finish_x", "finish_y", "finish_z", "finish_nx", "finish_ny", "finish_nz"]
PG_COLS = ["group_id", "part_id", "m", "q", "x", "xp", "y", "yp", "ct", "delta"]
CAV_COLS = ["cavity_id", "lattice_id", "x", "y", "z", "nx", "ny", "nz", "volt_share", "note"]
ME, C = bm.ME, bm.C_LIGHT
P0 = bm.P_MOM


# ---------------- ПРЕЖНИЙ алгоритм инжекции (дословно beam_injection до миграции) ----------------
def _gamma(v):
    return 1.0 / math.sqrt(1.0 - (v / C) ** 2)


def _frame(axis):
    a = np.array(axis, float); a = a / np.linalg.norm(a)
    u = np.array([-a[1], a[0], 0.0])
    nu = np.linalg.norm(u); u = u / nu if nu > 1e-9 else np.array([0.0, 1.0, 0.0])
    return a, u, np.cross(a, u)


def _fit_plane(parts, tol=1e-4):
    pos = np.array([p["pos"] for p in parts]); vel = np.array([p["vel"] for p in parts])
    vmean = vel.mean(axis=0); n = vmean / np.linalg.norm(vmean)
    centroid = pos.mean(axis=0)
    if np.abs((pos - centroid) @ n).max() > tol:
        raise ValueError("частицы НЕ компланарны")
    return centroid, n


def _find_site(mags, plane_pt, normal, vel_dir):
    N = len(mags); best = None
    for i in range(N):
        A = mags[i].exit; B = mags[(i + 1) % N].entry
        AB = B - A; L = np.linalg.norm(AB)
        if L < 1e-9:
            continue
        axis = AB / L
        denom = float(normal @ AB)
        if abs(denom) < 1e-12:
            continue
        t = float(normal @ (plane_pt - A)) / denom
        if -0.02 <= t <= 1.02:
            center = A + t * AB
            if float(vel_dir @ axis) <= 0:
                continue
            dist = np.linalg.norm(center - plane_pt)
            if best is None or dist < best[-1]:
                best = (i, center, axis, dist)
    if best is None:
        raise ValueError("нет дрейфа")
    return best[0], best[1], best[2]


def old_inject(mags, parts):
    """→ (start_center, start_axis, closed, finish_center, finish_axis, canon_states) как прежний bi.inject."""
    centroid, normal = _fit_plane(parts)
    vel_dir = np.mean([p["vel"] / np.linalg.norm(p["vel"]) for p in parts], axis=0)
    vel_dir /= np.linalg.norm(vel_dir)
    try:
        i, center, axis = _find_site(mags, centroid, normal, vel_dir)
        closed, fc, fa = True, None, None
    except ValueError:
        center, axis = centroid, vel_dir
        last = mags[-1]
        closed, fc, fa = False, np.asarray(last.exit, float), last.axis
    a, u, v = _frame(axis)
    states = []
    for p in parts:
        rel = p["pos"] - center
        x = float(rel @ u); y = float(rel @ v)
        vv = np.linalg.norm(p["vel"]); d = p["vel"] / vv
        delta = (_gamma(vv) * p["m"] * vv - P0) / P0
        states.append((x, (1.0 + delta) * float(d @ u), y, (1.0 + delta) * float(d @ v), 0.0, delta))
    return center, axis, closed, fc, fa, states


# ---------------- новая запись частицы (локально) ----------------
def to_local(p, center, axis):
    """Глоб. частица → (x, xp, y, yp, ct, delta, |продольный отступ|) в рамке плоскости (center, axis)."""
    a, u, v = _frame(axis)
    rel = p["pos"] - center
    vv = np.linalg.norm(p["vel"]); d = p["vel"] / vv
    da = float(d @ a)
    delta = (_gamma(vv) * p["m"] * vv - P0) / P0
    return (float(rel @ u), float(d @ u) / da, float(rel @ v), float(d @ v) / da, 0.0, delta,
            abs(float(rel @ a)))


def local_to_canon(loc):
    """ТА ЖЕ формула, что новый beam_injection.local_to_canonical."""
    x, xp, y, yp, ct, d = loc[:6]
    k = (1.0 + d) / math.sqrt(1.0 + xp * xp + yp * yp)
    return (x, k * xp, y, k * yp, ct, d)


def _f(v):
    return repr(float(v) + 0.0)                  # −0.0 → 0.0


def _plane_of_points(r):
    """p1..p3+dir → (центр, нормаль) Boris-плоскости или None. Центр = ((p1+p2)/2 + p3)/2 (авторинг c±u+v, c−v)."""
    try:
        p = [np.array([float(r[f"p{k}{c}"]) for c in "xyz"]) for k in (1, 2, 3)]
        dvec = np.array([float(r[f"dir{c}"]) for c in "xyz"])
    except (ValueError, KeyError):
        return None
    n = np.cross(p[1] - p[0], p[2] - p[0]); n = n / np.linalg.norm(n)
    if float(n @ dvec) < 0:
        n = -n
    return 0.5 * (0.5 * (p[0] + p[1]) + p[2]), n


def migrate(dry=False):
    rows = list(csv.DictReader(open(EXP, encoding="utf-8"), delimiter=";"))
    hdr = list(csv.reader(open(EXP, encoding="utf-8"), delimiter=";"))[0]
    if "start_x" in hdr:
        raise SystemExit("[ОТКАЗ] experiments.csv уже содержит start_x — миграция была выполнена")
    groups = {}
    for r in csv.DictReader(open(PG, encoding="utf-8"), delimiter=";"):
        groups.setdefault(r["group_id"], []).append(r)
    parts_of = {g: [dict(pos=np.array([float(r["x"]), float(r["y"]), float(r["z"])]),
                         vel=np.array([float(r["vx"]), float(r["vy"]), float(r["vz"])]),
                         m=float(r["m"]), q=float(r["q"]), part_id=r["part_id"], mr=r["m"], qr=r["q"])
                    for r in rs] for g, rs in groups.items()}
    eff, phys, coup, trans = mt.load_db()
    report = []; local_of = {}; new_rows = []; worst = 0.0; worst_long = 0.0
    next_gid = max(int(g) for g in groups) + 1
    for r in rows:
        exp, lat, gid = r["exp_id"], int(r["lattice_id"]), r["group_id"]
        with contextlib.redirect_stdout(io.StringIO()):
            mags = mt.build_line(lat, eff, phys, coup, trans)
        sc, sa, closed, fc, fa, canon = old_inject(mags, parts_of[gid])
        if closed:                                                     # кольцо: старт = плоскость САМОГО пучка (центроид,
            cen, _n = _fit_plane(parts_of[gid])                        # направление) — пересечение с опорной сидит в ~пм от
            vd = np.mean([p["vel"] / np.linalg.norm(p["vel"]) for p in parts_of[gid]], axis=0)   # него (невязка
            vd /= np.linalg.norm(vd)                                   # геометрии решётки) ⇒ чистые 0;0;0 / 1;0;0
            a0, _, _ = _frame(sa)
            dev = max(float(np.linalg.norm(cen - sc)), float(np.linalg.norm(np.cross(vd, a0))))
            globals()["_SHIFT"] = max(globals().get("_SHIFT", 0.0), dev)
            sc, sa = cen, vd
            a, u, v = _frame(sa)                                       # прежняя формула каноники — в новой рамке
            canon = []
            for p in parts_of[gid]:
                rel = p["pos"] - sc; vv = np.linalg.norm(p["vel"]); d = p["vel"] / vv
                de = (_gamma(vv) * p["m"] * vv - P0) / P0
                canon.append((float(rel @ u), (1 + de) * float(d @ u), float(rel @ v), (1 + de) * float(d @ v), 0.0, de))
        loc = [to_local(p, sc, sa) for p in parts_of[gid]]
        for L, s0 in zip(loc, canon):                                  # новая каноника = прежней
            s1 = local_to_canon(L)
            worst = max(worst, max(abs(a - b) for a, b in zip(s0, s1)))
            worst_long = max(worst_long, L[6])
        out = dict(r)
        for k in OLD_PLANE:
            out.pop(k, None)
        out.update({k: "" for k in NEW_PLANE})
        out.update(start_x=_f(sc[0]), start_y=_f(sc[1]), start_z=_f(sc[2]),
                   start_nx=_f(sa[0]), start_ny=_f(sa[1]), start_nz=_f(sa[2]))
        how = "кольцо (финиш = старт)"
        if not closed:
            pb = _plane_of_points(r)
            use_old = pb is None
            if pb is not None:
                cb, nb = pb
                a, _, _ = _frame(fa)
                on = abs(float(nb @ (fc - cb))) < 1e-5 and float(np.linalg.norm(np.cross(nb, a))) < 1e-6
                use_old = on
            if use_old:
                c, n = fc, fa; how = "отрезок, финиш = выход последнего элемента (как матричное ядро)"
            else:
                c, n = cb, nb; how = f"отрезок, финиш = p-плоскость Boris {np.round(cb, 6).tolist()}"
            out.update(finish_x=_f(c[0]), finish_y=_f(c[1]), finish_z=_f(c[2]),
                       finish_nx=_f(n[0]), finish_ny=_f(n[1]), finish_nz=_f(n[2]))
        key = tuple(tuple(L[:6]) for L in loc)
        if gid in local_of and local_of[gid][0] != key:                # общая группа, но другие координаты
            ng = str(next_gid); next_gid += 1
            local_of[ng] = (key, gid); out["group_id"] = ng
            how += f"; группа {gid} → клон {ng} (другая старт-плоскость)"
        elif gid not in local_of:
            local_of[gid] = (key, gid)
        new_rows.append(out)
        report.append(f"опыт {exp:>4} (реш. {lat}, гр. {out['group_id']}): старт "
                      f"{np.round(sc, 6).tolist()} n={np.round(np.asarray(sa) / np.linalg.norm(sa), 6).tolist()}; {how}")
    for g, parts in parts_of.items():                                  # группы без опытов — своя плоскость
        if g in local_of:
            continue
        c, n = _fit_plane(parts)
        local_of[g] = (tuple(tuple(to_local(p, c, n)[:6]) for p in parts), g)
        report.append(f"группа {g} без опытов: координаты относительно своей плоскости {np.round(c, 6).tolist()} "
                      f"n={np.round(n, 6).tolist()}")
    report.append(f"новая каноника vs прежняя: max|Δ| = {worst:.3e}; продольный отступ частиц от старта ≤ "
                  f"{worst_long * 1e3:.3e} мм (отбрасывался и раньше); старт колец = плоскость пучка, отличие от "
                  f"прежней рамки (пересечение с опорной) ≤ {globals().get('_SHIFT', 0.0):.3e} (м, рад)")
    if dry:
        return report, None
    bdir = os.path.join(TB, "backups", "explicit_planes_" + time.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(bdir, exist_ok=True)
    for p in (EXP, PG):
        shutil.copy2(p, bdir)
    i = hdr.index("p1x")
    new_hdr = hdr[:i] + NEW_PLANE + [h for h in hdr[i:] if h not in OLD_PLANE]
    with open(EXP, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=new_hdr, delimiter=";", extrasaction="ignore")
        w.writeheader(); w.writerows(new_rows)
    with open(PG, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(PG_COLS)
        for g in sorted(local_of, key=int):
            key, src = local_of[g]
            for p, L in zip(parts_of[src], key):
                w.writerow([g, p["part_id"], p["mr"], p["qr"]] + [_f(v) for v in L])
    if not os.path.exists(CAV):
        with open(CAV, "w", newline="", encoding="utf-8") as f:
            csv.writer(f, delimiter=";").writerow(CAV_COLS)
    return report, bdir


def main():
    ap = argparse.ArgumentParser(description="явные плоскости опыта + локальный пучок + rf_cavities")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rep, bdir = migrate(a.dry_run)
    print("\n".join(rep))
    if bdir:
        print(f"резервная копия: {bdir}")


if __name__ == "__main__":
    main()
