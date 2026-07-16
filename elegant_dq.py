# -*- coding: utf-8 -*-
"""Датасет DQ через КИКОВУЮ МОДЕЛЬ ELEGANT (CSBEND): для каждой комбинации токов обмоток (i1,i2,i3)
подбираем параметры elegant-поворота ANGLE, K1, K2, K3 так, чтобы траектория киковой модели elegant
повторила траекторию диффурной модели (Boris по 3D-картам), включая 5-метровый хвост за элементом.

Модель elegant CSBEND = combined-function секторный поворот: постоянные мультиполи вдоль изогнутой
опорной дуги, симплектический интегратор Форест-Рута 4-го порядка + СТАНДАРТНОЕ мультипольное ядро
Re/Im[(x+iy)ⁿ] (elegant_model.py). Кик на срезе (по отклонению du=u−y_ref от орбиты):
  Δqx = −[K0 + K1·du + ½K2(du²−v²) + (K3/6)(du³−3du·v²)]·ds
  Δqy = +[K1·v + K2·du·v + (K3/6)(3du²·v−v³)]·ds,   K0=ANGLE/L_eff.
elegant-силы: ANGLE=θ0=(q/p)B·L, K1=(q/p)G, K2=(q/p)S, K3=(q/p)O.

Всё в отдельной папке workspace/elegant_dq/: dataset CSV + наложения траекторий + поверхности
(PNG + вращаемый HTML). Запуск: py -3.14 elegant_dq.py"""
import os, sys, csv
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
sys.path.insert(0, "C:/Sasha/Work/SILA/SynchrotronSim_archive/legacy/scripts")
import elegant_model as elm
import dq_dataset as dqd                       # переиспуем карты/орбиту/diff_probe/arc_coeffs
import kickdoc_compare as kd
import matrix_tracker as mt
from bip_matrix import QP

OUT = os.path.join("workspace", "research", "elegant_dq"); os.makedirs(OUT, exist_ok=True)
HALF, LEFF, N, GRID, PROBES, DRIFT, DS = dqd.HALF, dqd.LEFF, dqd.N, dqd.GRID, dqd.PROBES, dqd.DRIFT, dqd.DS
DS_G, DS_F = dqd.DS_G, dqd.DS_F
plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans"})


def csbend(theta0, K1, K2, K3, u0, Xo, Uo, v0=0.0, want_traj=False):
    """Кик elegant CSBEND (Форест-Рут) над зоной 2·half; мультиполи на du=u−y_ref(X). + 5-м хвост.
    Возврат: (u_вых, u'_вых, U после хвоста) или полная траектория (X,u,v)."""
    K0 = theta0 / LEFF
    x = -HALF; u, v = u0, v0
    qx, qy = elm._slopes_to_momenta(0.0, 0.0)
    tr = [(x, u, v)]
    for _ in range(N):
        for st in range(4):
            up, vp = elm._momenta_to_slopes(qx, qy)
            dg = DS_G * elm.DRIFT4[st]
            x += dg; u += up * dg; v += vp * dg
            if elm.KICK4[st] == 0.0:
                break
            du = u - float(np.interp(x, Xo, Uo))
            if want_traj and (abs(du) > 0.1 or abs(v) > 0.1):
                return np.array(tr)                           # частица улетела — обрезаем рисунок
            duc = 0.03 if du > 0.03 else (-0.03 if du < -0.03 else du)   # клип кика вне физ.зоны
            vc = 0.03 if v > 0.03 else (-0.03 if v < -0.03 else v)       # (±4мм точны; стабильность фита)
            kf = DS_F * elm.KICK4[st]
            qx += -(K0 + K1 * duc + 0.5 * K2 * (duc * duc - vc * vc) + (K3 / 6.0) * (duc**3 - 3 * duc * vc * vc)) * kf
            qy += (K1 * vc + K2 * duc * vc + (K3 / 6.0) * (3 * duc * duc * vc - vc**3)) * kf
        if want_traj:
            tr.append((x, u, v))
    up, vp = elm._momenta_to_slopes(qx, qy)
    if want_traj:
        tr.append((x + DRIFT, u + up * DRIFT, v + vp * DRIFT))
        return np.array(tr)
    return u, up, u + up * DRIFT


def fit_elegant(fld, Xo, Uo):
    """Подгонка (ANGLE,K1,K2,K3) elegant под хвосты диффура. Старт — из разложения поля вдоль дуги."""
    diffs = [dqd.diff_probe(fld, u0) for u0 in PROBES]
    B, G, S, O = dqd.arc_coeffs(fld, Xo, Uo)
    x0 = [QP * B * LEFF, QP * G, QP * S, QP * O]           # θ0,K1,K2,K3 старт

    def resid(p):
        th, K1, K2, K3 = p; r = []
        for u0, (ue_d, _, U5_d) in zip(PROBES, diffs):
            try:
                ue_k, _, U5_k = csbend(th, K1, K2, K3, u0, Xo, Uo)
                r += [(ue_k - ue_d) * 1e3, (U5_k - U5_d) * 1e3]
            except Exception:
                r += [1e6, 1e6]
        return r

    sol = least_squares(resid, x0, x_scale=[1e-1, 1e1, 1e2, 1e4], method="trf")
    return sol.x, float(np.sqrt(np.mean(np.array(sol.fun) ** 2)))


def build_dataset():
    rows = []
    print(f"Elegant CSBEND датасет: {GRID} А ({len(GRID)**3} комб.), 5-м хвост, L_eff={LEFF:.4f} м")
    for i1 in GRID:
        for i2 in GRID:
            for i3 in GRID:
                fld = dqd.comb_field(i1, i2, i3)
                Xo, Uo = dqd.orbit(fld)
                (th, K1, K2, K3), rms = fit_elegant(fld, Xo, Uo)
                B, Gc, Sc, Oc = th / (QP * LEFF), K1 / QP, K2 / QP, K3 / QP
                rows.append([i1, i2, i3, f"{th*1e3:.5f}", f"{K1:.5f}", f"{K2:.4f}", f"{K3:.3f}",
                             f"{B:.6e}", f"{Gc:.6e}", f"{Sc:.6e}", f"{Oc:.6e}", f"{rms:.4f}"])
    with open(os.path.join(OUT, "elegant_dq_dataset.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["i1_A", "i2_A", "i3_A", "ANGLE_mrad", "K1_m2", "K2_m3", "K3_m4",
                    "B_Tl", "G_Tlm", "S_Tlm2", "O_Tlm3", "fit_rms_mm"])
        w.writerows(rows)
    rmsv = np.array([float(r[-1]) for r in rows])
    print(f"[OK] {len(rows)} строк → {OUT}/elegant_dq_dataset.csv")
    print(f"невязка elegant CSBEND: медиана {np.median(rmsv):.3f} мм, 90%={np.percentile(rmsv,90):.3f}, "
          f"макс {rmsv.max():.3f}; <1мм {100*np.mean(rmsv<1):.0f}%")
    return rows


def overlays():
    """Наложение траекторий: диффур (Boris, сплошн.) vs elegant CSBEND (пунктир) для 6 настроек."""
    setts = [(50, 0, 0), (-100, 0, 0), (100, 100, 100), (0, 100, 100), (50, 50, 50), (100, -100, -100)]
    pr = [0.0, 2e-3, -2e-3, 4e-3, -4e-3]; cols = ["#000", "#c0392b", "#2980b9", "#27ae60", "#e67e22"]
    fig, axs = plt.subplots(2, 3, figsize=(15, 8)); axs = axs.ravel()
    for k, (i1, i2, i3) in enumerate(setts):
        ax = axs[k]; fld = dqd.comb_field(i1, i2, i3); Xo, Uo = dqd.orbit(fld)
        (th, K1, K2, K3), rms = fit_elegant(fld, Xo, Uo)
        for c, u0 in enumerate(pr):
            Sb, Ub, Vb, upb, vpb = kd.boris(fld, HALF, u0, 0.0, DS)
            ax.plot((Sb - HALF) * 1e3, Ub * 1e3, "-", color=cols[c], lw=2.2, alpha=0.5)
            tr = csbend(th, K1, K2, K3, u0, Xo, Uo, want_traj=True)
            m = tr[:, 0] <= HALF + 1e-9                      # только зона элемента (хвост опускаем на графике)
            ax.plot(tr[m, 0] * 1e3, tr[m, 1] * 1e3, "--", color=cols[c], lw=1.0)
        ax.axvspan(-HALF * 1e3, HALF * 1e3, color="0.93", zorder=0)
        ax.set_title(f"I=({i1},{i2},{i3})  θ0={th*1e3:+.1f}мрад K1={K1:+.2f}\nневязка хвоста {rms:.2f}мм", fontsize=9)
        ax.set_xlabel("X (вдоль), мм"); ax.set_ylabel("u (гор.), мм"); ax.grid(alpha=0.3)
    fig.suptitle("DQ: диффур (сплошн.) vs киковая модель elegant CSBEND (пунктир), пучок ±4мм", fontsize=12, weight="bold")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "elegant_overlays.png"), dpi=130); plt.close(fig)
    print(f"[OK] {OUT}/elegant_overlays.png")


def surfaces(rows):
    """Поверхности elegant-сил (ANGLE,K1,K2,K3) над (I1, I2+I3): PNG + вращаемый HTML."""
    hdr = ["i1_A", "i2_A", "i3_A", "ANGLE_mrad", "K1_m2", "K2_m3", "K3_m4", "B_Tl", "G_Tlm", "S_Tlm2", "O_Tlm3", "fit_rms_mm"]
    R = [dict(zip(hdr, r)) for r in rows]
    i1v = sorted(set(int(r["i1_A"]) for r in R)); sums = sorted(set(int(r["i2_A"]) + int(r["i3_A"]) for r in R))
    I1g, Sg = np.meshgrid(i1v, sums, indexing="ij")
    xs = np.array([int(r["i1_A"]) for r in R]); ys = np.array([int(r["i2_A"]) + int(r["i3_A"]) for r in R])
    rr = np.array([float(r["fit_rms_mm"]) for r in R])

    def grid(c):
        Z = np.full((len(i1v), len(sums)), np.nan)
        for a, i1 in enumerate(i1v):
            for b, s in enumerate(sums):
                vals = [float(r[c]) for r in R if int(r["i1_A"]) == i1 and int(r["i2_A"]) + int(r["i3_A"]) == s and float(r["fit_rms_mm"]) < 5]
                if vals:
                    Z[a, b] = np.mean(vals)
        return Z

    for c, title, zl in [("ANGLE_mrad", "ANGLE (θ0)", "мрад"), ("K1_m2", "K1", "м⁻²"),
                         ("K2_m3", "K2", "м⁻³"), ("K3_m4", "K3", "м⁻⁴")]:
        Z = grid(c); zs = np.array([float(r[c]) for r in R]); tag = c.split("_")[0]
        fig = plt.figure(figsize=(9, 7)); ax = fig.add_subplot(111, projection="3d")
        ax.plot_surface(I1g, Sg, np.ma.masked_invalid(Z), cmap="viridis", alpha=0.7, edgecolor="0.3", lw=0.3)
        sc = ax.scatter(xs, ys, zs, c=rr, cmap="Reds", vmin=0, vmax=10, s=26, edgecolor="k", lw=0.3, depthshade=False)
        ax.set_xlabel("I1, А"); ax.set_ylabel("I2+I3, А"); ax.set_zlabel(f"{title}, {zl}")
        ax.set_title(f"elegant {title}", weight="bold"); ax.view_init(elev=24, azim=-58)
        cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.1); cb.set_label("невязка, мм")
        fig.savefig(os.path.join(OUT, f"elegant_surface_{tag}.png"), dpi=130, bbox_inches="tight"); plt.close(fig)
        figh = go.Figure()
        figh.add_trace(go.Surface(x=i1v, y=sums, z=Z.T, colorscale="Viridis", opacity=0.8, showscale=False))
        figh.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode="markers",
                       marker=dict(size=4, color=rr, colorscale="Reds", cmin=0, cmax=10, colorbar=dict(title="невязка,мм")),
                       text=[f"I=({r['i1_A']},{r['i2_A']},{r['i3_A']}) нев={float(r['fit_rms_mm']):.2f}" for r in R],
                       hovertemplate="%{text}<br>" + title + "=%{z:.4g}<extra></extra>"))
        figh.update_layout(title=f"elegant {title} над (I1, I2+I3)",
                           scene=dict(xaxis_title="I1, А", yaxis_title="I2+I3, А", zaxis_title=f"{title}, {zl}"), width=950, height=780)
        figh.write_html(os.path.join(OUT, f"elegant_surface_{tag}.html"), include_plotlyjs="inline")
    print(f"[OK] поверхности ANGLE,K1,K2,K3 (PNG+HTML) → {OUT}/")


if __name__ == "__main__":
    rows = build_dataset()
    overlays()
    surfaces(rows)
