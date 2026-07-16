# -*- coding: utf-8 -*-
"""Поверхности физкоэффициентов DQ (θ→B, G, S) над плоскостью токов из арочного датасета.
Каждый коэф. — отдельно: статичный PNG + интерактивный HTML (plotly, вращается). Две поверхности:
ДУГА (симметр. проезд, физичная калибровка, цветная) и ОСЬ (прямой приосевой пробник, серая) — видно
расхождение. Точки раскрашены по невязке матчинга движка. Плоскость — (I1, I2+I3): при i3=0 это (i1, i2).
Запуск: py -3.14 dq_surfaces.py
        py -3.14 dq_surfaces.py --csv dq_arc_i1i2.csv --prefix dq_i1i2_surface --ylabel "i2, А (i3=0)" """
import os, csv, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bip_matrix import QP
import matrix_tracker as mt

OUT = os.path.join("workspace", "research", "kickdoc")
_ap = argparse.ArgumentParser()
_ap.add_argument("--csv", default="dq_arc_vs_axis.csv", help="входной датасет (в workspace/research/kickdoc)")
_ap.add_argument("--prefix", default="dq_surface", help="префикс выходных файлов")
_ap.add_argument("--ylabel", default="I2+I3, А", help="подпись оси Y (плоскость тока)")
_A = _ap.parse_args()
R = list(csv.DictReader(open(os.path.join(OUT, _A.csv), encoding="utf-8-sig"), delimiter=";"))
LEFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}["DQ_I1_1A"]
plt.rcParams.update({"font.size": 10, "font.family": "DejaVu Sans"})
RELI = 1.0                                                    # поверхность — по точкам с невязкой < 1 мм


def I(r, k):
    return round(float(r[k]), 4)                              # ключ тока (float — держит нецелый шаг)


def ykey(r):
    return round(I(r, "i2_A") + I(r, "i3_A"), 4)              # ось Y = i2+i3 (при i3=0 = i2)


def fnum(x):
    try:
        v = float(x); return v if np.isfinite(v) else None
    except (ValueError, TypeError):
        return None


def col(r, k, tr=lambda v: v):
    v = fnum(r.get(k))
    return None if v is None else tr(v)


i1v = sorted(set(I(r, "i1_A") for r in R))
sums = sorted(set(ykey(r) for r in R))
tb = lambda v: v / 1e3 / (QP * LEFF)                          # θ[мрад] → B[Тл]

PAN = [
    ("B", "B — изгибное поле", "B, Тл",
     lambda r: col(r, "theta_arc_mrad", tb), lambda r: col(r, "theta_axis_mrad", tb)),
    ("G", "G — градиент (квадруполь)", "G, Тл/м",
     lambda r: col(r, "G_arc_Tlm"), lambda r: col(r, "G_axis_Tlm")),
    ("S", "S — секступоль", "S, Тл/м²",
     lambda r: col(r, "S_arc_Tlm2"), lambda r: col(r, "S_axis_Tlm2")),
    # октуполь исключён: 2-мм пучок его не определяет (∝r³) — фитим только θ,G,S
]


def grid(fn, use_reli):
    """Поверхность коэф. над (I1, I2+I3); use_reli — только надёжные (невязка<RELI)."""
    Z = np.full((len(i1v), len(sums)), np.nan)
    for a, i1 in enumerate(i1v):
        for b, s in enumerate(sums):
            vals = []
            for r in R:
                if I(r, "i1_A") != i1 or ykey(r) != s:
                    continue
                if use_reli and float(r["fit_rms_mm"]) >= RELI:
                    continue
                v = fn(r)
                if v is not None:
                    vals.append(v)
            if vals:
                Z[a, b] = np.mean(vals)
    return Z


for tag, title, zlab, arc_fn, axis_fn in PAN:
    Zarc = grid(arc_fn, True)
    Zaxis = grid(axis_fn, False) if axis_fn else None
    pts = [(I(r, "i1_A"), ykey(r), arc_fn(r), float(r["fit_rms_mm"]), r) for r in R]
    pts = [p for p in pts if p[2] is not None]
    xs = np.array([p[0] for p in pts]); ys = np.array([p[1] for p in pts])
    zs = np.array([p[2] for p in pts]); rr = np.array([p[3] for p in pts])

    # --- статичный PNG ---
    I1g, Sg = np.meshgrid(i1v, sums, indexing="ij")
    fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111, projection="3d")
    if Zaxis is not None:
        ax.plot_surface(I1g, Sg, np.ma.masked_invalid(Zaxis), color="0.6", alpha=0.28,
                        linewidth=0, rstride=1, cstride=1)
    ax.plot_surface(I1g, Sg, np.ma.masked_invalid(Zarc), cmap="viridis", alpha=0.85,
                    linewidth=0.3, edgecolor="0.3", rstride=1, cstride=1)
    sc = ax.scatter(xs, ys, zs, c=rr, cmap="Reds", vmin=0, vmax=2, s=22,
                    edgecolor="k", linewidth=0.3, depthshade=False)
    ax.set_xlabel("I1, А"); ax.set_ylabel(_A.ylabel); ax.set_zlabel(zlab)
    ax.set_title(f"{title}  (дуга — цвет, ось — серая)", fontsize=12, weight="bold")
    ax.view_init(elev=24, azim=-58)
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1); cb.set_label("невязка матчинга, мм")
    fig.savefig(os.path.join(OUT, f"{_A.prefix}_{tag}.png"), dpi=130, bbox_inches="tight")
    plt.close(fig)

    # --- интерактивный HTML (вращается) ---
    figh = go.Figure()
    figh.add_trace(go.Surface(x=i1v, y=sums, z=Zarc.T, colorscale="Viridis", opacity=0.9,
                              showscale=False, name="ДУГА (симметр.)", showlegend=True))
    if Zaxis is not None:
        figh.add_trace(go.Surface(x=i1v, y=sums, z=Zaxis.T, colorscale=[[0, "gray"], [1, "gray"]],
                                  opacity=0.35, showscale=False, name="ОСЬ (прямой пробник)", showlegend=True))
    figh.add_trace(go.Scatter3d(
        x=xs, y=ys, z=zs, mode="markers",
        marker=dict(size=3, color=rr, colorscale="Reds", cmin=0, cmax=2,
                    colorbar=dict(title="невязка,мм", x=0.97, len=0.7), line=dict(width=0.4, color="black")),
        text=[f"I=({I(p[4],'i1_A'):g},{I(p[4],'i2_A'):g},{I(p[4],'i3_A'):g})А  невязка={p[3]:.3f}мм" for p in pts],
        hovertemplate="%{text}<br>" + zlab + "=%{z:.4g}<extra></extra>", name="дуга: точки"))
    figh.update_layout(title=f"{title} — DQ (дуга vs ось), {_A.csv}",
                       scene=dict(xaxis_title="I1, А", yaxis_title=_A.ylabel, zaxis_title=zlab,
                                  domain=dict(x=[0.0, 0.92], y=[0.0, 1.0])),   # место под шкалу
                       autosize=True, legend=dict(x=0.01, y=0.99), margin=dict(l=0, r=0, t=40, b=0))
    figh.write_html(os.path.join(OUT, f"{_A.prefix}_{tag}.html"), include_plotlyjs="inline",
                    default_width="100%", default_height="100vh", config={"responsive": True})
    print(f"  {tag}: {_A.prefix}_{tag}.png + {_A.prefix}_{tag}.html  "
          f"(поверхность по {int(np.isfinite(Zarc).sum())} ячейкам)")

print(f"[OK] 3 коэф. (B,G,S) × (PNG + вращаемый HTML) из {_A.csv} → {OUT}")
