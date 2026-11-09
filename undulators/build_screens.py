# -*- coding: utf-8 -*-
"""Размещение ФОТОННЫХ ЭКРАНОВ: манифест (относительно источника) → замороженный ГЛОБАЛЬНЫЙ абсолют в БД.

ЗАЧЕМ ИМЕННО ТАК. Экран/поглотитель — это железо front-end, прикрученное к полу: при перекосе магнита
оно НЕ должно ехать. Поэтому экран задаётся удобно (относительно источника), но РАЗРЕШАЕТСЯ один раз по
НОМИНАЛЬНОЙ геометрии источника и ЗАМОРАЖИВАЕТСЯ в абсолют (центр + нормаль + оси плоскости) в сайдкар
workspace/tables/screens.csv. Потом перекос источника двигает ПЯТНО по неподвижному экрану, а не сам экран
(анализ: screen_analysis/screen_spectrum для ондуляторов, bend_screen для поворотного СИ — читают ФАКТ. геометрию).

Экран — НЕ трекинговый элемент (его нет на орбите) ⇒ отдельная таблица, НЕ в lattice_configs/connectivity.

МАНИФЕСТ (undulators/screens.csv), колонки ';':
  screen_id ; source ; mode ; dist_m ; los_h_mrad ; los_v_mrad ; tilt_h_mrad ; tilt_v_mrad ; roll_mrad ;
  width_mm ; height_mm ; gx ; gy ; gz ; extract_mrad
  • mode=axis   (ОНДУЛЯТОР) — source=group_id (und_01..); центр = pivot + dist вдоль оси ±(los_h,los_v).
  • mode=global (ОНДУЛЯТОР) — центр явно (gx,gy,gz, геодезия); нормаль смотрит на source.
  • mode=bend   (ДИПОЛЬ/волноотвод) — source=instance_id диполя; центр = точка дуги под УГЛОМ ОТБОРА
                extract_mrad (от входа, по ходу; пусто=центр диполя) + dist вдоль КАСАТЕЛЬНОЙ ±(los_h,los_v).
  • los_*  — малый увод линии визирования; tilt_* — наклон САМОГО экрана; roll — крен плоскости. Всё 0 по умолч.

ЗАПУСК
  py -3.14 undulators/build_screens.py --lattice 27                       # собрать по манифесту
  py -3.14 undulators/build_screens.py --lattice 27 --list-bends          # меню диполей (instance_id для source)
"""
import os, sys, csv, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt
import bend_sr as bsr

TB = os.path.join("workspace", "tables")
EG = os.path.join(TB, "element_groups.csv")
LC = os.path.join(TB, "lattice_configs.csv")
SCR = os.path.join(TB, "screens.csv")
HERE = os.path.dirname(__file__)
WORLD_UP = np.array([0.0, 0.0, 1.0])

SCR_HDR = ["lattice_id", "screen_id", "source_ref", "source_kind",
           "cx", "cy", "cz", "nx", "ny", "nz", "ux", "uy", "uz", "vx", "vy", "vz",
           "width_mm", "height_mm", "mode", "dist_m", "extract_mrad",
           "los_h_mrad", "los_v_mrad", "tilt_h_mrad", "tilt_v_mrad", "roll_mrad"]


def _rot(vec, axis, ang):
    v = np.asarray(vec, float); axis = np.asarray(axis, float)
    c, s = math.cos(ang), math.sin(ang)
    return v * c + np.cross(axis, v) * s + axis * float(np.dot(axis, v)) * (1.0 - c)


def _unit(v):
    v = np.asarray(v, float); n = np.linalg.norm(v)
    return v / n if n else v


def screen_frame(n, roll=0.0, world_up=WORLD_UP):
    n = _unit(n); up = np.asarray(world_up, float)
    if abs(float(np.dot(n, up))) > 0.999:
        up = np.array([1.0, 0.0, 0.0])
    u = _unit(np.cross(up, n)); v = _unit(np.cross(n, u))
    if roll:
        u, v = _rot(u, n, roll), _rot(v, n, roll)
    return u, v


def load_groups(lattice):
    out = {}
    if not os.path.exists(EG):
        raise SystemExit(f"нет {EG} — сначала вставьте ондуляторы (insert_undulators.py)")
    for r in csv.DictReader(open(EG, encoding="utf-8-sig"), delimiter=";"):
        if r.get("lattice_id") != str(lattice):
            continue
        out[r["group_id"]] = dict(
            pivot=np.array([float(r["pivot_x"]), float(r["pivot_y"]), float(r["pivot_z"])]),
            yaw=float(r["yaw"]), pitch=float(r["pitch"]), roll=float(r["roll"]))
    return out


def load_lattice(lattice):
    rows = list(csv.reader(open(LC, encoding="utf-8-sig"), delimiter=";"))
    ci = {c: i for i, c in enumerate(rows[0])}
    byid = {r[ci["instance_id"]]: r for r in rows[1:] if r and r[ci["lattice_id"]] == str(lattice)}
    return ci, byid


def bend_source(iid, extract_mrad, ci, byid, E=6.0):
    """Геометрия отбора поворотного СИ: точка на дуге под углом отбора + касательная там (НОМ. или ФАКТ.)."""
    r = byid.get(str(iid))
    if r is None:
        raise SystemExit(f"диполь instance_id={iid} не найден (см. --list-bends)")
    kind = r[ci["kind"]].strip()
    C = np.array([float(r[ci[a]] or 0) for a in ("x", "y", "z")])
    axis, uhat, vhat = (np.asarray(v, float) for v in
                        mt.frame(*(float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))))
    L = float(r[ci["eff_length"]] or 0); matB = float(r[ci["mat_B"]] or 0)
    theta0 = mt.QP * matB * L                                  # угол изгиба (знак по изгибу)
    if abs(theta0) < 1e-12:
        raise SystemExit(f"диполь {iid} ({r[ci['device_type_id']]}) не изгибает (θ≈0) — не источник СИ")
    rho = L / theta0; O = C + rho * uhat                        # центр кривизны (как в matrix_tracker)
    phi = (abs(theta0) / 2.0) if extract_mrad is None else (extract_mrad * 1e-3)
    if phi < -1e-9 or phi > abs(theta0) + 1e-9:
        raise SystemExit(f"диполь {iid}: угол отбора {phi*1e3:.3f} мрад вне дуги [0,{abs(theta0)*1e3:.3f}] мрад")
    alpha = -theta0 / 2.0 + math.copysign(phi, theta0)         # от входа по ходу
    P = O + _rot(C - O, vhat, alpha)
    t = _unit(_rot(axis, vhat, alpha))
    rho_, B, ec = bsr.bend_geom(theta0, L, E)
    return dict(P=P, t=t, theta0=theta0, rho=rho, L=L, B=B, ec=ec, phi=phi,
                uhat=uhat, vhat=vhat, dev=r[ci["device_type_id"]], kind=kind)


def list_bends(lattice, E=6.0):
    ci, byid = load_lattice(lattice)
    by_dev = {}
    for iid, r in byid.items():
        if r[ci["kind"]].strip() not in ("dipole", "dq"):
            continue
        matB = float(r[ci["mat_B"]] or 0); L = float(r[ci["eff_length"]] or 0)
        th = mt.QP * matB * L
        if abs(th) < 1e-12:
            continue
        by_dev.setdefault(r[ci["device_type_id"]], []).append((int(iid), th, L))
    print(f"ДИПОЛИ-ИСТОЧНИКИ решётки {lattice} (по типам; любой instance_id типа годится как source):")
    print(f"  {'device_type':16s} {'шт':>4} {'θ,мрад':>9} {'ρ,м':>8} {'B,Т':>7} {'ε_c,кэВ':>8}  пример instance_id")
    for dev, lst in sorted(by_dev.items()):
        lst.sort(); iid0, th, L = lst[0]
        rho, B, ec = bsr.bend_geom(th, L, E)
        print(f"  {dev:16s} {len(lst):>4} {th*1e3:>+9.3f} {rho:>8.2f} {B:>7.4f} {ec:>8.3f}  {iid0}")


def build(lattice, manifest, out=SCR, E=6.0):
    groups = load_groups(lattice)
    ci, byid = load_lattice(lattice)
    rows = list(csv.DictReader(open(manifest, encoding="utf-8-sig"), delimiter=";"))
    recs = []
    for row in rows:
        sid = (row.get("screen_id") or "").strip()
        if not sid or sid.startswith("#"):
            continue
        src = row["source"].strip(); mode = (row.get("mode") or "axis").strip().lower()
        f = lambda k, d=0.0: float(row[k]) if row.get(k) not in (None, "") else d
        los_h, los_v = f("los_h_mrad") * 1e-3, f("los_v_mrad") * 1e-3
        tilt_h, tilt_v = f("tilt_h_mrad") * 1e-3, f("tilt_v_mrad") * 1e-3
        roll = f("roll_mrad") * 1e-3
        w, h = f("width_mm", 50.0), f("height_mm", 50.0)
        extract = row.get("extract_mrad")
        extract = float(extract) if extract not in (None, "") else None

        if mode == "bend":                                     # ── ПОВОРОТНОЕ СИ (волноотвод) ──
            bs = bend_source(src, extract, ci, byid, E)
            P, t = bs["P"], bs["t"]
            ut = _unit(np.cross(WORLD_UP, t))                  # горизонт-попереч. к касательной
            d = _unit(_rot(_rot(t, WORLD_UP, los_h), ut, los_v))
            dist = f("dist_m")
            if dist <= 0:
                raise SystemExit(f"[{sid}] mode=bend требует dist_m>0")
            C = P + dist * d; skind = "bend"
            extract_used = bs["phi"] * 1e3
            note = (f"диполь {bs['dev']} θ={bs['theta0']*1e3:+.2f} отбор {extract_used:.2f} мрад "
                    f"ε_c={bs['ec']:.2f} кэВ")
        else:                                                  # ── ОНДУЛЯТОР (axis/global) ──
            grp = groups.get(src)
            if grp is None:
                raise SystemExit(f"[{sid}] источник-ондулятор {src} не найден в element_groups реш.{lattice}")
            axis, uhat, vhat = (np.asarray(v, float) for v in mt.frame(grp["yaw"], grp["pitch"], grp["roll"]))
            P = grp["pivot"]; skind = "undulator"; extract_used = 0.0
            if mode == "global":
                C = np.array([f("gx"), f("gy"), f("gz")]); d = _unit(C - P); dist = float(np.linalg.norm(C - P))
            else:
                d = _unit(_rot(_rot(axis, vhat, los_h), uhat, los_v)); dist = f("dist_m")
                if dist <= 0:
                    raise SystemExit(f"[{sid}] mode=axis требует dist_m>0")
                C = P + dist * d
            note = f"ондулятор {src}"

        n0 = -d
        us, vs = screen_frame(n0)
        n = _unit(_rot(_rot(n0, vs, tilt_h), us, tilt_v))
        us, vs = screen_frame(n, roll)
        recs.append([str(lattice), sid, src, skind,
                     *(f"{x:.9f}" for x in C), *(f"{x:.9f}" for x in n),
                     *(f"{x:.9f}" for x in us), *(f"{x:.9f}" for x in vs),
                     f"{w:.4f}", f"{h:.4f}", mode, f"{dist:.6f}", f"{extract_used:.5f}",
                     f"{f('los_h_mrad'):.4f}", f"{f('los_v_mrad'):.4f}",
                     f"{f('tilt_h_mrad'):.4f}", f"{f('tilt_v_mrad'):.4f}", f"{f('roll_mrad'):.4f}"])
        print(f"  {sid:12s} ← {src:10s} {mode:6s} D={dist:7.3f} м  центр=({C[0]:.3f},{C[1]:.3f},{C[2]:.3f})  "
              f"{w:.0f}×{h:.0f} мм  [{note}]")

    keep = []
    if os.path.exists(out):
        old = list(csv.reader(open(out, encoding="utf-8"), delimiter=";"))
        if old:
            li = old[0].index("lattice_id")
            keep = [r for r in old[1:] if r and r[li] != str(lattice)]
        import shutil; shutil.copy(out, out + ".bak")
    with open(out, "w", newline="", encoding="utf-8") as fo:
        wr = csv.writer(fo, delimiter=";"); wr.writerow(SCR_HDR); wr.writerows(keep + recs)
    print(f"[OK] {out} ← решётка {lattice}: {len(recs)} экранов (чужих строк сохранено {len(keep)})")


def main():
    ap = argparse.ArgumentParser(description="разместить фотонные экраны по манифесту (заморозить в БД)")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--manifest", default=os.path.join(HERE, "screens.csv"))
    ap.add_argument("--list-bends", action="store_true", help="показать диполи-источники (для source в mode=bend)")
    a = ap.parse_args()
    if a.list_bends:
        list_bends(a.lattice)
    else:
        build(a.lattice, a.manifest)


if __name__ == "__main__":
    main()
