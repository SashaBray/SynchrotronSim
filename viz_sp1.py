# -*- coding: utf-8 -*-
"""Оверлей-визуализация конфигурации SP1 (решётка 111): КИК (матрица, mat_*) vs ДИФФУР (C++ Boris по
картам поля) над магнитами 1-го супериода + первого элемента 2-го SP. Две траектории в ОДНОЙ сцене:
  • sp1_kick_vs_diff.html — самодостаточный интерактивный 3D (крутить/зумить в браузере);
  • sp1_kick_vs_diff.png  — 2D-проекции ВИД СВЕРХУ (X-Y) и СБОКУ (X-Z), где расхождение сразу видно.
КИК = зелёный (проектная орбита по mat_*), ДИФФУР = красный (реальная симуляция по полю).
Запуск: py -3.14 viz_sp1.py"""
import os
import numpy as np
import pyvista as pv
import plotly.graph_objects as go
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from collections import defaultdict
import visualize as vz
import viz_html as vh

WS = "workspace"
HTML_DIR = os.path.join(WS, "results", "html")
LAT, KICK_EXP, DIFF_EXP = 111, 59, 58
C_KICK, C_DIFF = "#2f9e44", "#e15759"


def _load(exp):
    return vz.load_trajectories(os.path.join(WS, "results", f"trajectories_exp{exp}.bin"))


def _magnets(fig, lc, hw, kind):
    """Магниты решётки как лёгкие блоки (glyph) по габаритам STL, полупрозрачные."""
    cache, grouped = {}, defaultdict(list)
    for r in lc:
        dev = r["device_type_id"]; stl = hw.get(dev, {}).get("stl_path", "-")
        if stl in ("-", "", "nan"):
            continue
        if stl not in cache:
            p = os.path.join(WS, "geometry", stl)
            m = pv.read(p).scale(0.001, inplace=False) if os.path.exists(p) else None
            cache[stl] = pv.Box(bounds=m.bounds).triangulate() if m is not None else None
        if cache[stl] is None:
            continue
        grouped[dev].append((float(r["yaw"]), float(r["pitch"]), float(r["roll"]),
                             float(r["x"]), float(r["y"]), float(r["z"])))
    for dev, tl in grouped.items():
        big = vz.build_combined_mesh(cache[hw[dev]["stl_path"]], np.asarray(tl, float))
        V = np.asarray(big.points); F = np.asarray(big.faces).reshape(-1, 4)[:, 1:]
        fig.add_trace(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                                color=vh.KINDCOL.get(kind.get(dev, ""), "silver"), opacity=0.28,
                                name=dev, hoverinfo="skip", showlegend=False))


def _add_traj(fig, exp, color, label):
    for i, p in enumerate(_load(exp)):
        if len(p) < 2:
            continue
        central = (i == 0)
        fig.add_trace(go.Scatter3d(x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                                   line=dict(width=7 if central else 2, color=color),
                                   opacity=1.0 if central else 0.30,
                                   name=(label if central else None), showlegend=central,
                                   legendgroup=label, hovertext=label, hoverinfo="text"))


def html_scene(lc, hw, kind):
    fig = go.Figure()
    _magnets(fig, lc, hw, kind)
    ref = np.array([(float(r["x"]), float(r["y"]), float(r["z"]))
                    for r in sorted(lc, key=lambda r: int(r["instance_id"]))])
    fig.add_trace(go.Scatter3d(x=ref[:, 0], y=ref[:, 1], z=ref[:, 2], mode="lines",
                               line=dict(width=2, color="dimgray", dash="dot"), name="опорный путь", opacity=0.5))
    _add_traj(fig, KICK_EXP, C_KICK, "КИК (матрица, mat_*)")
    _add_traj(fig, DIFF_EXP, C_DIFF, "ДИФФУР (Boris по полю)")
    fig.update_layout(title="SynchrotronSim — SP1 + первый элемент SP2: КИК vs ДИФФУР (решётка 111)",
                      scene=dict(aspectmode="data", xaxis_title="X, м", yaxis_title="Y, м",
                                 zaxis_title="Z, м", bgcolor="white"),
                      autosize=True, margin=dict(l=0, r=0, t=40, b=0),
                      legend=dict(x=0.99, xanchor="right", y=0.99, yanchor="top",
                                  bgcolor="rgba(255,255,255,0.75)", bordercolor="gray", borderwidth=1))
    os.makedirs(HTML_DIR, exist_ok=True)
    out = os.path.join(HTML_DIR, "sp1_kick_vs_diff.html")
    fig.write_html(out, include_plotlyjs="inline", default_width="100%", default_height="100vh",
                   config={"responsive": True})
    print(f"[OK] 3D: {out} ({os.path.getsize(out)/1e6:.1f} МБ) — крутить/зумить в браузере")


def png_projections(lc):
    K = _load(KICK_EXP); D = _load(DIFF_EXP)
    fig, ax = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    # магниты (центры) серыми точками по X
    mx = [float(r["x"]) for r in lc]; my = [float(r["y"]) for r in lc]
    for a, comp, ylab, title in [(ax[0], 1, "Y, м (гориз.)", "ВИД СВЕРХУ (X-Y): горизонтальный увод"),
                                  (ax[1], 2, "Z, м (верт.)", "ВИД СБОКУ (X-Z): вертикальный увод")]:
        for i, p in enumerate(K):
            a.plot(p[:, 0], p[:, comp], color=C_KICK, lw=2.2 if i == 0 else 0.7,
                   alpha=1 if i == 0 else 0.3, label="КИК (матрица)" if i == 0 else None, zorder=3)
        for i, p in enumerate(D):
            a.plot(p[:, 0], p[:, comp], color=C_DIFF, lw=2.2 if i == 0 else 0.7,
                   alpha=1 if i == 0 else 0.3, label="ДИФФУР (Boris)" if i == 0 else None, zorder=3)
        if comp == 1:
            a.scatter(mx, my, s=8, c="gray", alpha=0.5, marker="s", label="магниты", zorder=1)
        a.set_ylabel(ylab); a.set_title(title); a.grid(alpha=0.3); a.legend(fontsize=9, loc="upper left")
    ax[1].set_xlabel("X, м (вдоль пучка)")
    fig.suptitle("SP1 + первый элемент SP2: расхождение кик ↔ диффур (СИЛА, решётка 111)", fontsize=13)
    fig.tight_layout()
    out = os.path.join(HTML_DIR, "sp1_kick_vs_diff.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    print(f"[OK] 2D: {out} — проекции сверху/сбоку")


def main():
    hw = {r["device_type_id"]: r for r in vh._read("hardware_library.csv")}
    kind = {r["device_type_id"]: r["kind"] for r in vh._read("element_transforms.csv")}
    lc = [r for r in vh._read("lattice_configs.csv") if r["lattice_id"] == str(LAT)]
    html_scene(lc, hw, kind)
    png_projections(lc)
    print("\nКИК (зелёный) = проектная орбита по mat_*; ДИФФУР (красный) = реальная симуляция по полю.")


if __name__ == "__main__":
    main()
