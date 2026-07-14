# -*- coding: utf-8 -*-
"""Сравнение траекторий: интегральный трекинг по картам полей (Boris) vs матричная/идеальная
модель (аналитический мультиполь с фит-коэффициентом). Для КАЖДОГО элемента:
  * физкоэффициент из фита поля  ∂ⁿ⁻¹Bz/∂yⁿ⁻¹ = (n−1)!·Re(Kₙ),
  * эффективная длина L_eff = ∫c(x)dx / c_тело,
  * коэффициент перевода ток→коэф (на Ампер, поле ∝ ток),
  * расхождение траекторий Boris(карта) vs идеал, на пучке ±2мм.
Линейные поля (квадруполь) — по полной 5-мм карте (трилин. точна); нелинейные (SP ∝r², OP ∝r³)
— по мелкой v3 (2мм, зеркало X<0), фит по кольцу r≥3мм (у оси трилин. на грубой сетке врёт)."""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
from field_compare_viz import fit_coeffs
from _ct_asym_track import make_field, boris_track, C, ME, QE, V, GAMMA0

KD = r"C:/Sasha/Work/SILA/ansis approximations/KD"
FIGS = os.path.join("workspace", "bip_compare", "figs")
P = GAMMA0 * ME * V            # импульс
QP = QE / P                    # q/p  (= −1/Bρ по модулю)
BRHO = P / abs(QE)
NOM_A = 100.0

# (тег, карта_фит, карта_трек, N, halfmap_track, (r1,r2) кольца фита, L_lte, назв)
ELEMS = [
    ("QP_L",  "QP_L/QP_L_100A_full_v2.fld",  "QP_L/QP_L_100A_full_v2.fld",  2, False, (3e-3, 9e-3), 0.295, "Квадруполь длинный"),
    ("QP_sh", "QP_sh/QP_sh_100A_full_v3.fld", "QP_sh/QP_sh_100A_full_v2.fld", 2, False, (3e-3, 9e-3), 0.162, "Квадруполь короткий"),
    ("SP_I1", "SP/SP_I1_100A_v3.fld",         "SP/SP_I1_100A_v3.fld",         3, True,  (3e-3, 9e-3), 0.166, "Секступоль (обмотка I1)"),
    ("SP_I2", "SP/SP_I2_100A_v3.fld",         "SP/SP_I2_100A_v3.fld",         3, True,  (3e-3, 9e-3), 0.166, "Секступоль (обмотка I2)"),
    ("OP",    "OP/OP_100A_v3.fld",            "OP/OP_100A_v3.fld",            4, True,  (5e-3, 15e-3), 0.09, "Октуполь"),
]
FACT = {2: 1, 3: 2, 4: 6}       # (n−1)!
UNIT = {2: "Тл/м", 3: "Тл/м²", 4: "Тл/м³"}


def coeff_profile(h, B, N, r1, r2):
    """Комплексный c(x) порядка N по срезам (кольцо r1..r2), свёрнутый к x≥0."""
    j0, k0 = h["ny"] // 2, h["nz"] // 2
    R = max(2, int(math.ceil(r2 / h["sy"])) + 1)
    js = slice(j0 - R, j0 + R + 1); ks = slice(k0 - R, k0 + R + 1)
    yg = (np.arange(h["ny"])[js] - j0) * h["sy"]; zg = (np.arange(h["nz"])[ks] - k0) * h["sz"]
    Z, Y = np.meshgrid(zg, yg); rr = np.sqrt(Y ** 2 + Z ** 2); mask = (rr >= r1) & (rr <= r2)
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    c = np.array([fit_coeffs(Y, Z, B[i, js, ks, 1], B[i, js, ks, 2], (N,), mask)[N] for i in range(h["nx"])])
    pos = xg >= -1e-9
    xa = xg[pos]; ca = c[pos]
    o = np.argsort(xa); return xa[o], ca[o]


def ideal_fn(xa, ca, N):
    def f(pos):
        x, y, z = pos
        ax = abs(x)
        if ax > xa[-1]: return np.zeros(3)
        k = np.interp(ax, xa, ca.real) + 1j * np.interp(ax, xa, ca.imag)
        fc = k * (y + 1j * z) ** (N - 1)
        return np.array([0.0, fc.imag, fc.real])
    return f


def analyse(tag, fld_fit, fld_trk, N, halftrk, ring, L_lte, title, ax):
    hf, Bf = load(os.path.join(KD, fld_fit))
    ht, Bt = load(os.path.join(KD, fld_trk))
    xa, ca = coeff_profile(hf, Bf, N, *ring)
    reC = ca.real
    ib = int(np.argmax(np.abs(reC))); c_body = reC[ib]                 # тело = макс |Re c|
    deriv_body = FACT[N] * c_body                                      # ∂ⁿ⁻¹Bz/∂yⁿ⁻¹ на теле
    Sint = 2.0 * np.trapezoid(FACT[N] * reC, xa)                        # ∫ deriv dx (интегр. сила, знаковый)
    Leff = Sint / (deriv_body + 1e-30)                                 # L_eff = ∫G dx / G_тело (шум сокращается)
    Knorm = QP * deriv_body                                            # нормир. K на теле (1/м^N) при нашем p

    f_real = make_field(ht, Bt, "mirror" if halftrk else "full")
    f_ideal = ideal_fn(xa, ca, N)
    x0, x1 = ht["mnx"] + 1e-6, ht["mxx"] - 1e-6
    if halftrk: x0 = -ht["mxx"] + 1e-6                                 # полная длина через зеркало
    amps = [2e-3, 4e-3, 6e-3]; devs = []
    cols = plt.cm.viridis(np.linspace(0, 0.8, len(amps)))
    rows = []
    for a, col in zip(amps, cols):
        tr, _ = boris_track(a, f_real, x0, x1); ti, _ = boris_track(a, f_ideal, x0, x1)
        ypr = (tr[-1, 1] - tr[-2, 1]) / (tr[-1, 0] - tr[-2, 0])
        ypi = (ti[-1, 1] - ti[-2, 1]) / (ti[-1, 0] - ti[-2, 0])
        dev = abs(ypr - ypi) / (abs(ypi) + 1e-30) * 100; devs.append(dev)
        rows.append((a, ypr * 1e3, ypi * 1e3, dev))
        ax.plot(tr[:, 0] * 1e3, tr[:, 1] * 1e3, "-", color=col, lw=1.7, label=f"реал y0={a*1e3:.0f}")
        ax.plot(ti[:, 0] * 1e3, ti[:, 1] * 1e3, "--", color=col, lw=1.1)
    dev2 = rows[0][3]
    ax.axvline(0, color="0.6", ls=":", lw=0.8); ax.grid(True, alpha=0.3)
    ax.set_xlabel("X, мм"); ax.set_ylabel("Y, мм"); ax.legend(fontsize=6.5, ncol=1)
    ax.set_title(f"{tag}: {title}\nрасх. Boris↔идеал @2мм = {dev2:.1f}%  (сплошн.=карта, штрих=идеал)", fontsize=9)

    print(f"\n=== {tag}  ({title}) ===  N={N}")
    print(f"  физкоэф на теле ∂^{N-1}Bz/∂y^{N-1} = {deriv_body:+.5g} {UNIT[N]}   (при 100А)")
    print(f"  на Ампер: {deriv_body/NOM_A:+.5g} {UNIT[N]}/А     нормир. K{N-1}={Knorm:+.4g} 1/м^{N} @наш p")
    print(f"  L_eff = {Leff*1e3:.1f} мм   (L из .lte = {L_lte*1e3:.0f} мм)   ∫сила = {Sint:+.4g} {UNIT[N]}·м")
    print(f"  расхождение траекторий Boris(карта) vs идеал:")
    for a, r, i, d in rows:
        print(f"     y0={a*1e3:.0f}мм:  реал Y'={r:+8.4f}  идеал Y'={i:+8.4f} мрад   расх {d:4.1f}%")
    return dict(tag=tag, N=N, deriv=deriv_body, perA=deriv_body / NOM_A, Knorm=Knorm,
                Leff=Leff, Sint=Sint, dev2=dev2, L_lte=L_lte)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(16, 9))
    axes = axes.ravel()
    res = []
    for i, e in enumerate(ELEMS):
        res.append(analyse(*e[:8], axes[i]))
    axes[-1].axis("off")
    txt = "СВОДКА  (карта 100А)\n\n" + f"{'элемент':7} {'коэф/тело':>12} {'L_eff/мм':>9} {'расх@2мм':>9}\n"
    for r in res:
        txt += f"{r['tag']:7} {r['deriv']:+12.4g} {r['Leff']*1e3:9.1f} {r['dev2']:8.1f}%\n"
    txt += f"\nBρ={BRHO*1e3:.3f} мТл·м,  p={P:.4e} кг·м/с\nq/p={QP:.3f} (1/Тл·м)"
    axes[-1].text(0.02, 0.98, txt, va="top", ha="left", family="monospace", fontsize=9)
    fig.suptitle("Сравнение: интегральный трекинг по картам (Boris) vs идеальный мультиполь (матричная модель)\n"
                 "новые карты KD, 100А. Секступоль/октуполь — по мелкой карте v3 (2мм)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    p = os.path.join(FIGS, "fig_element_matrix_compare.png"); fig.savefig(p, dpi=115); plt.close(fig)
    print("\nфигура:", os.path.basename(p))


if __name__ == "__main__":
    main()
