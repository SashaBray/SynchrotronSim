# -*- coding: utf-8 -*-
"""Off-screen СКРИНШОТЫ сцены визуализатора SynchrotronSim с камерой, наведённой на заданный
элемент/точку. Переиспользует геометрию и траектории из visualize.py. Для обследования
траектории у конкретных магнитов (например, у 4-го диполь-квадруполя).

Запуск:
  py -3.14 viz_shot.py --exp 54 --mag 29 --view top  --span 8 --out dq4_top.png
  py -3.14 viz_shot.py --exp 54 --mag 29 --view persp --span 6 --out dq4_persp.png
  py -3.14 viz_shot.py --exp 54 --focus 15.25 0.48 0 --view side --span 10 --out dq4_side.png
"""
import os, sys, argparse
import numpy as np
import pandas as pd
from collections import defaultdict
import pyvista as pv
pv.OFF_SCREEN = True
import visualize as V

WS = os.path.join(V.SCRIPT_DIR, "workspace") + os.sep
OUTDIR = os.path.join("workspace", "results", "shots")
os.makedirs(OUTDIR, exist_ok=True)


def build_scene(plotter, exp_id, hi_mag=None):
    """Собрать сцену (магниты + опорный путь + траектории) в plotter. Возврат: dict меток."""
    exp_df = pd.read_csv(WS + "tables/experiments.csv", sep=";", skipinitialspace=True)
    lattice_df = pd.read_csv(WS + "tables/lattice_configs.csv", sep=";", skipinitialspace=True)
    hw_df = pd.read_csv(WS + "tables/hardware_library.csv", sep=";", skipinitialspace=True)
    lat_id = int(exp_df[exp_df["exp_id"] == exp_id]["lattice_id"].iloc[0])
    cur = lattice_df[lattice_df["lattice_id"] == lat_id]
    hw = {}
    for row in hw_df.itertuples():
        d = str(getattr(row, "device_type_id", "")).strip()
        if d and d != "nan":
            hw[d] = (str(getattr(row, "stl_path", "")).strip(),
                     str(getattr(row, "Active", "")).strip().lower() in ("true", "1", "yes"))
    stl_cache = {}
    for stl_name, _ in hw.values():
        if stl_name in ("-", "", "nan") or stl_name in stl_cache:
            continue
        p = os.path.join(WS, "geometry", stl_name)
        if os.path.exists(p):
            m = pv.read(p); m.scale(0.001, inplace=True); stl_cache[stl_name] = m
        else:
            stl_cache[stl_name] = None
    cols = ["device_type_id", "yaw", "pitch", "roll", "x", "y", "z"]
    dev = cur["device_type_id"].astype(str).str.strip().values
    arr = {c: cur[c].astype(float).values for c in cols[1:]}
    grouped = defaultdict(list)
    for i in range(len(dev)):
        if dev[i] in hw and stl_cache.get(hw[dev[i]][0]) is not None:
            grouped[dev[i]].append(tuple(arr[c][i] for c in cols[1:]))
    op = V.VIZ_CFG.get("lattice_opacity", {})
    for d, tl in grouped.items():
        base = stl_cache[hw[d][0]]
        big = V.build_combined_mesh(base, np.asarray(tl, float))
        active = hw[d][1]
        plotter.add_mesh(big, color="silver",
                         opacity=float(op.get("EM", 0.25) if active else op.get("PM", 0.4)),
                         show_edges=True)
    # опорный путь
    rpts = V.load_reference_path_points(WS, exp_id)
    if rpts is not None and len(rpts) >= 2:
        poly = pv.PolyData(np.asarray(rpts, float)); n = len(rpts)
        c = np.full((1, n + 1), n, dtype=np.int_); c[0, 1:] = np.arange(n); poly.lines = c
        plotter.add_mesh(poly, color="dimgray", line_width=3, opacity=0.8)
    # траектории
    parts = V.load_trajectories(os.path.join(WS, "results", f"trajectories_exp{exp_id}.bin"))
    colors = V.VIZ_CFG.get("particle_colors", ["magenta", "cyan", "lime", "yellow", "red"])
    for i, pts in enumerate(parts):
        if len(pts) < 2:
            continue
        poly = pv.PolyData(pts); n = len(pts)
        c = np.full((1, n + 1), n, dtype=np.int_); c[0, 1:] = np.arange(n); poly.lines = c
        plotter.add_mesh(poly, color=colors[i % len(colors)], line_width=5)
    return {"n_parts": len(parts)}


def mag_entry(mag_idx):
    import beam_injection as bi
    plane, ordered, states, info = bi.inject(101, 1, source="current")
    m = ordered[mag_idx]
    return np.array(m.entry, float), getattr(m, "dev", "?"), getattr(m, "kind", "?")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", type=int, default=54)
    ap.add_argument("--mag", type=int, default=None, help="индекс магнита (фокус на его вход)")
    ap.add_argument("--focus", type=float, nargs=3, default=None, help="точка фокуса X Y Z")
    ap.add_argument("--view", choices=["top", "side", "persp"], default="top")
    ap.add_argument("--span", type=float, default=8.0, help="ширина кадра, м")
    ap.add_argument("--out", default="shot.png")
    a = ap.parse_args()

    if a.mag is not None:
        f, dev, kind = mag_entry(a.mag); tag = f"m{a.mag} {dev} ({kind})"
    elif a.focus:
        f = np.array(a.focus, float); tag = "точка"
    else:
        f = np.array([15.25, 0.48, 0.0]); tag = "DQ#4"

    pl = pv.Plotter(off_screen=True, window_size=[1600, 1000])
    pl.set_background("white")
    build_scene(pl, a.exp, hi_mag=a.mag)
    # ЗЕЛЁНАЯ проектная орбита (входы магнитов) + маркеры DQ для сравнения
    import beam_injection as bi
    plane, ordered, states, info = bi.inject(101, 1, source="current")
    ents = np.array([m.entry for m in ordered], float)
    poly = pv.PolyData(ents[:120]); n = 120
    cc = np.full((1, n + 1), n, dtype=np.int_); cc[0, 1:] = np.arange(n); poly.lines = cc
    pl.add_mesh(poly, color="green", line_width=3, opacity=0.7)
    for j, m in enumerate(ordered[:120]):
        if getattr(m, "kind", "") == "dq":
            pl.add_mesh(pv.Sphere(radius=0.05, center=m.entry), color="blue")
    D = a.span
    if a.view == "top":                                    # вид сверху (плоскость кольца X-Y)
        pos = (f[0], f[1], f[2] + D); up = (0, 1, 0)
    elif a.view == "side":                                 # вид сбоку (X-Z, видно вертикаль)
        pos = (f[0], f[1] + D, f[2]); up = (0, 0, 1)
    else:                                                  # перспектива сверху-сбоку
        pos = (f[0] - D * 0.4, f[1] + D * 0.5, f[2] + D * 0.8); up = (0, 1, 0)
    pl.camera_position = [pos, tuple(f), up]
    pl.enable_parallel_projection()
    pl.camera.parallel_scale = D * 0.5
    pl.add_text(f"exp{a.exp}  фокус: {tag}  вид: {a.view}", font_size=11, color="black")
    out = os.path.join(OUTDIR, a.out)
    pl.screenshot(out)
    print(f"[OK] {out}  (фокус {f.round(2)}, {tag})")


if __name__ == "__main__":
    main()
