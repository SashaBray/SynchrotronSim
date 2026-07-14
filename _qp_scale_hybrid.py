# -*- coding: utf-8 -*-
"""Гибридный ресайз карты квадруполя: РАЗНЫЕ коэффициенты сжатия тела (c_body) и фринджа (c_fringe).
Кусочно-линейное отображение новой X в старую: в теле old=|X|/c_body (тело z-независимо ⇒ ТОЧНО),
во фриндже old=Xb+(|X|−Xb·c_body)/c_fringe (фриндж сжат в c_fringe). Полукарта→полный квадруполь.
Идея: c_body — свободный ТОЧНЫЙ рычаг силы ∫G (тело масштабируется без нарушения Максвелла);
c_fringe — рычаг края (ширина фринджа, краевой фокус) ценой малого нарушения Максвелла.
Проверяем: L_eff = c_body·L_body + c_fringe·L_fringe (линейно), ∇·B(c_body,c_fringe) — плоско по
c_body, растёт от c_fringe=1. Выход: workspace/qp_scale/qp_hybrid.png."""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import _qp_scale as qs

OUT = qs.OUT


def body_halfwidth(h, B, frac=0.98):
    """Полуширина ТЕЛА Xb: где |G| впервые < frac·Gmax при сканировании от центра (X≥0)."""
    xg, G, _, _ = qs.grad_profile(h, B); Gn = np.abs(G) / np.abs(G).max()
    pos = np.where(xg >= 0)[0]; xp = xg[pos]; gp = Gn[pos]
    return next((xp[i] for i in range(len(gp)) if gp[i] < frac), xp[-1])


def divB_fringe(h, B):
    """∇·B RMS в зоне ФРИНДЖА (0.1<|G|/Gmax<0.9) у оси — там нарушение Максвелла при c_fringe≠1."""
    xg, G, _, _ = qs.grad_profile(h, B); Gn = np.abs(G) / np.abs(G).max()
    fr = (Gn > 0.1) & (Gn < 0.9)
    yg = h["mny"] + np.arange(h["ny"]) * h["sy"]; zg = h["mnz"] + np.arange(h["nz"]) * h["sz"]
    jm = [j for j in range(h["ny"]) if abs(yg[j]) <= 6e-3]; km = [k for k in range(h["nz"]) if abs(zg[k]) <= 6e-3]
    dBx = np.gradient(B[:, :, :, 0], h["sx"], axis=0)
    dBy = np.gradient(B[:, :, :, 1], h["sy"], axis=1)
    dBz = np.gradient(B[:, :, :, 2], h["sz"], axis=2)
    d = (dBx + dBy + dBz)[fr][:, jm][:, :, km]
    return np.sqrt(np.mean(d ** 2)) if d.size else 0.0


def resize(h, B, cb, cf, Xb):
    """Гибридный ресайз полной карты: тело ×cb, фриндж ×cf. Кусочно-линейное old(new), odd Bx учтён картой."""
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]; sx = h["sx"]; Xmax = xg[-1]
    Xbp = Xb * cb                                              # новая полуширина тела
    Xmax_new = Xbp + (Xmax - Xb) * cf
    nh = int(round(Xmax_new / sx)); xn = np.arange(-nh, nh + 1) * sx
    a = np.abs(xn)
    oldabs = np.where(a <= Xbp, a / cb, Xb + (a - Xbp) / cf)   # тело: a/cb; фриндж: Xb+(a−Xbp)/cf
    oldX = np.sign(xn) * oldabs
    fx = (oldX - xg[0]) / sx; i0 = np.clip(fx.astype(int), 0, h["nx"] - 2); fr = (fx - i0)[:, None, None, None]
    Bn = (1 - fr) * B[i0] + fr * B[i0 + 1]                     # линейная интерполяция вдоль X
    hn = dict(h); hn["nx"] = len(xn); hn["mnx"] = xn[0]
    return hn, Bn


def main():
    h0, B0 = qs.load(qs.MAP); h, B = qs.mirror_full(h0, B0)
    Xb = body_halfwidth(h, B)
    Gm = np.abs(qs.grad_profile(h, B)[1]).max()
    iG0 = qs.intG(h, B); Leff0 = iG0 / (-Gm)                   # знак: G<0
    Lbody = 2 * Xb; Lfr = Leff0 - Lbody
    print(f"Тело: полуширина Xb={Xb*1e3:.0f}мм (полное тело {Lbody*1e3:.0f}мм); "
          f"L_eff={Leff0*1e3:.1f}мм ⇒ вклад фринджа L_fr={Lfr*1e3:.1f}мм")
    print("Проверка тождества L_eff(cb,cf) = cb·L_body + cf·L_fr:")
    print(f"{'cb':>5}{'cf':>6}{'∫G,мТл·м':>11}{'L_eff,мм':>10}{'cb·Lb+cf·Lf':>13}{'фриндж,мм':>11}{'∇·B-фрдж':>11}")
    for cb, cf in [(1.0, 1.0), (0.5, 1.0), (2.0, 1.0), (1.0, 0.6), (1.0, 1.4), (1.5, 0.7), (0.7, 1.5)]:
        hh, BB = resize(h, B, cb, cf, Xb)
        iG = qs.intG(hh, BB); Leff = iG / (-Gm); pred = cb * Lbody + cf * Lfr
        fw = qs.fringe_width(hh, BB); dv = divB_fringe(hh, BB)
        print(f"{cb:>5.1f}{cf:>6.1f}{iG*1e3:>11.2f}{Leff*1e3:>10.1f}{pred*1e3:>13.1f}{fw*1e3:>11.1f}{dv*1e3:>11.3f}")

    # 2D-скан ∇·B и L_eff по (cb, cf)
    CB = np.linspace(0.5, 2.0, 7); CF = np.linspace(0.6, 1.6, 7)
    DV = np.zeros((len(CF), len(CB))); LE = np.zeros_like(DV)
    for a, cf in enumerate(CF):
        for b, cb in enumerate(CB):
            hh, BB = resize(h, B, cb, cf, Xb)
            DV[a, b] = divB_fringe(hh, BB) * 1e3; LE[a, b] = qs.intG(hh, BB) / (-Gm) * 1e3
    base = divB_fringe(h, B) * 1e3
    print(f"\n2D: ∇·B минимально по всей оси c_body (тело точно), растёт только с |c_fringe−1|.")
    print(f"    базлайн ∇·B={base:.2f} мТл/м; при c_fringe=1 (любой c_body) ∇·B≈{DV[np.argmin(np.abs(CF-1))].mean():.2f}; "
          f"при c_fringe=0.6 ∇·B≈{DV[0].mean():.2f} мТл/м")
    visualize(CB, CF, DV, LE, base)
    print("Фигура:", os.path.join(OUT, "qp_hybrid.png"))


def visualize(CB, CF, DV, LE, base):
    fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
    im0 = ax[0].imshow(DV, origin="lower", aspect="auto", cmap="RdYlGn_r",
                       extent=[CB[0], CB[-1], CF[0], CF[-1]], vmin=base * 0.9, vmax=base * 1.6)
    ax[0].axhline(1.0, color="k", ls="--", lw=1.2, label="c_fringe=1 (физичный фриндж)")
    ax[0].set_xlabel("c_body (масштаб тела)"); ax[0].set_ylabel("c_fringe (масштаб фринджа)")
    ax[0].set_title("Нарушение Максвелла ∇·B RMS, мТл/м\n(плоско по c_body — тело точно; растёт от c_fringe=1)")
    ax[0].legend(fontsize=8); fig.colorbar(im0, ax=ax[0])
    im1 = ax[1].imshow(LE, origin="lower", aspect="auto", cmap="viridis",
                       extent=[CB[0], CB[-1], CF[0], CF[-1]])
    ax[1].set_xlabel("c_body"); ax[1].set_ylabel("c_fringe")
    ax[1].set_title("Эффективная длина L_eff, мм\n= c_body·L_body + c_fringe·L_fringe (два рычага)")
    fig.colorbar(im1, ax=ax[1])
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "qp_hybrid.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
