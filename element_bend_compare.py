# -*- coding: utf-8 -*-
"""Диполь-квадруполи (DQ, обмотки I1/I2/I3) и диполи (DP_L, DP_sh): сравнение траекторий
Boris(карта) vs идеальная комбинированная модель (диполь b0 + градиент G + секступоль).
Поле диполь+квад линейно/константно → грубая 5-мм карта интерполируется точно.
Для DQ: коэффициент перевода обмотка(ток)→(b0, G) — это строки матрицы связи.
Изгиб при 100А сильный (низкая энергия), поэтому для трекинга поле масштабируем до ~0.1 рад."""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
from field_compare_viz import fit_coeffs
import field_line_viz as flv
from _ct_asym_track import make_field, boris_track, ME, QE, V, GAMMA0

KD = r"C:/Sasha/Work/SILA/ansis approximations/KD"
FIGS = os.path.join("workspace", "bip_compare", "figs")
P = GAMMA0 * ME * V; BRHO = P / abs(QE); NOM_A = 100.0

# (тег, карта(v2 полная 5мм), тип, назв)
ELEMS = [
    ("DQ_I1", "DQ/DQ_I1_100A_v2.fld", "dq", "Диполь-квадруполь, обмотка I1"),
    ("DQ_I2", "DQ/DQ_I2_100A_v2.fld", "dq", "Диполь-квадруполь, обмотка I2"),
    ("DQ_I3", "DQ/DQ_I3_100A_v2.fld", "dq", "Диполь-квадруполь, обмотка I3"),
    ("DP_L",  "DP_L/dp_L_PM_1A.fld",  "dp", "Диполь длинный (ПМ)"),
    ("DP_sh", "DP_sh/dp_sh_PM_1mm.fld", "dp", "Диполь короткий (ПМ)"),
]


def prof(h, B, orders, r1, r2):
    """Профили c_n(x) для набора порядков (кольцо r1..r2), свёрнуто к x≥0."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    R = max(2, int(math.ceil(r2 / h["sy"])) + 1)
    js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
    yg = (np.arange(h["ny"])[js] - j0) * h["sy"]; zg = (np.arange(h["nz"])[ks] - k0) * h["sz"]
    Z, Y = np.meshgrid(zg, yg); rr = np.sqrt(Y ** 2 + Z ** 2); mask = (rr >= r1) & (rr <= r2)
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    C = {n: [] for n in orders}
    for i in range(h["nx"]):
        d = fit_coeffs(Y, Z, B[i, js, ks, 1], B[i, js, ks, 2], orders, mask)
        for n in orders: C[n].append(d[n])
    pos = xg >= -1e-9; xa = xg[pos]; o = np.argsort(xa)
    xa = xa[o]
    return xa, {n: np.array(C[n])[pos][o] for n in orders}


def ideal_fn(xa, Cn):
    def f(pos):
        x, y, z = pos; ax = abs(x)
        if ax > xa[-1]: return np.zeros(3)
        w = y + 1j * z; fc = 0j
        for n, ca in Cn.items():
            k = np.interp(ax, xa, ca.real) + 1j * np.interp(ax, xa, ca.imag)
            fc += k * w ** (n - 1)
        return np.array([0.0, fc.imag, fc.real])
    return f


def integ(xa, ca):
    return 2.0 * np.trapezoid(ca.real, xa)         # ∫ по полной длине (свёрнуто ×2)


def analyse(tag, fld, kind, title, ax):
    h, B = load(os.path.join(KD, fld))
    orders = (1, 2, 3)
    xa, Cn = prof(h, B, orders, 3e-3, 9e-3)
    b0 = Cn[1].real; G = Cn[2].real; S = 2 * Cn[3].real
    ib = int(np.argmax(np.abs(b0))) if kind == "dp" else int(np.argmax(np.abs(G)))
    b0_body = b0[ib]; G_body = G[ib]
    Bdl = integ(xa, Cn[1])                          # ∫Bz dx  [Тл·м] диполь
    Gdl = integ(xa, Cn[2])                          # ∫G dx   [Тл]   квадруполь
    Sdl = integ(xa, Cn[3]) * 2                       # ∫ ∂²Bz dx
    Leff_d = Bdl / (b0_body + 1e-30) if abs(b0_body) > 1e-9 else float("nan")
    Leff_q = Gdl / (G_body + 1e-30) if abs(G_body) > 1e-9 else float("nan")
    angle = Bdl / BRHO                               # рад при нашем p (диагностика)

    # трекинг при слабом изгибе: угол ~0.02 рад, чтобы пучок оставался у оси (в зоне фита ≤9мм)
    scale = 0.02 * BRHO / (abs(Bdl) + 1e-30)
    scale = min(scale, 1.0)
    fr = make_field(h, B, "full")
    def fr_s(pos): return fr(pos) * scale
    Cn_s = {n: ca * scale for n, ca in Cn.items()}
    fi = ideal_fn(xa, Cn_s)
    x0, x1 = h["mnx"] + 1e-6, h["mxx"] - 1e-6
    devs = []
    for a, col in [(2e-3, "tab:red"), (-2e-3, "tab:blue")]:
        tr, _ = boris_track(a, fr_s, x0, x1); ti, _ = boris_track(a, fi, x0, x1)
        ypr = (tr[-1, 1] - tr[-2, 1]) / (tr[-1, 0] - tr[-2, 0])
        ypi = (ti[-1, 1] - ti[-2, 1]) / (ti[-1, 0] - ti[-2, 0])
        devs.append(abs(ypr - ypi) / (abs(ypi) + 1e-30) * 100)
        ax.plot(tr[:, 0] * 1e3, tr[:, 1] * 1e3, "-", color=col, lw=1.7, label=f"реал y0={a*1e3:+.0f}")
        ax.plot(ti[:, 0] * 1e3, ti[:, 1] * 1e3, "--", color="k", lw=1.0)
    dev = float(np.mean(devs))
    ax.axvline(0, color="0.6", ls=":", lw=0.8); ax.grid(True, alpha=0.3)
    ax.set_xlabel("X, мм"); ax.set_ylabel("Y, мм"); ax.legend(fontsize=6.5)
    ax.set_title(f"{tag}: {title}\nрасх. Boris↔идеал = {dev:.2f}% (изгиб ×{scale:.3f})", fontsize=8.5)

    print(f"\n=== {tag}  ({title}) ===")
    print(f"  диполь b0(тело)   = {b0_body*1e3:+.4g} мТл    ∫Bz dx = {Bdl*1e3:+.4g} мТл·м   L_eff(дип)={Leff_d*1e3:.1f} мм")
    print(f"  градиент G(тело)  = {G_body:+.4g} Тл/м    ∫G dx  = {Gdl:+.4g} Тл    L_eff(квад)={Leff_q*1e3:.1f} мм")
    print(f"  секступоль ∂²Bz(тело)={S[ib]:+.3g} Тл/м²   ∫∂²Bz dx={Sdl:+.3g} Тл/м")
    if kind == "dq":
        print(f"  на Ампер:  b0={b0_body/NOM_A*1e6:+.4g} мкТл/А   G={G_body/NOM_A*1e3:+.4g} мТл/м/А")
    print(f"  угол изгиба ∫Bz/Bρ = {angle*1e3:+.3g} мрад @наш p (100А)   расх траекторий={dev:.2f}%")
    return dict(tag=tag, kind=kind, b0=b0_body, G=G_body, Bdl=Bdl, Gdl=Gdl,
                Leff_d=Leff_d, Leff_q=Leff_q, dev=dev, angle=angle)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9)); axes = axes.ravel()
    res = [analyse(*e, axes[i]) for i, e in enumerate(ELEMS)]
    axes[-1].axis("off")
    t = "СВОДКА  DQ / диполи (100А / ПМ)\n\n"
    t += f"{'элем':6}{'b0тело/мТл':>12}{'Gтело Тл/м':>12}{'∫Bz мТл·м':>11}{'расх%':>7}\n"
    for r in res:
        t += f"{r['tag']:6}{r['b0']*1e3:+12.4g}{r['G']:+12.4g}{r['Bdl']*1e3:+11.4g}{r['dev']:6.2f}\n"
    t += f"\nBρ={BRHO*1e3:.3f} мТл·м @наш p"
    axes[-1].text(0.02, 0.98, t, va="top", ha="left", family="monospace", fontsize=8.5)
    fig.suptitle("DQ (диполь-квадруполь, обмотки I1/I2/I3) и диполи: Boris(карта) vs идеал (диполь+градиент)\n"
                 "новые карты KD. Изгиб масштабирован для мягкой траектории (расхождение линейно по току)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(FIGS, "fig_element_bend_compare.png"); fig.savefig(p, dpi=115); plt.close(fig)
    print("\nфигура:", os.path.basename(p))


if __name__ == "__main__":
    main()
