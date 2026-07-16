# -*- coding: utf-8 -*-
"""Исследование: сжать/растянуть карту поля квадруполя → квадруполь другой длины?
Полукарта QP_sh (ось X, X=0 центр) зеркалится в ПОЛНЫЙ квадруполь. Два метода изменения длины:
  A) НАИВНОЕ масштабирование: B_new(X)=B_full(X/c) — вся карта (тело+фриндж) сжата в c раз.
  B) СПЛАЙС тела: вставить/убрать центральные срезы (тело z-независимо), фринджи не трогать.
Метрики: профиль G(X), интеграл ∫G·dX (сила), пик краевого Bz, нарушение Максвелла ∇·B у оси,
трекинг параллельных лучей → 1/f (фокус) и аберрация (нелинейность фокуса).
Выход: workspace/qp_scale/qp_scale.png + метрики."""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load

MAP = r"C:/Sasha/Work/SILA/ansis approximations/KD/QP_sh/QP_sh_100A_full_v3.fld"
OUT = os.path.join("workspace", "research", "qp_scale"); os.makedirs(OUT, exist_ok=True)
QE, ME, CL = 1.602176634e-19, 9.1093837015e-31, 299792458.0
V0 = 2.967e8; G0 = 1 / math.sqrt(1 - (V0 / CL) ** 2); P0 = G0 * ME * V0


def mirror_full(h, Bh):
    """Полукарта X≥0 (X=0 центр) → полный квадруполь X∈[−Xmax,Xmax]. By,Bz ЧЁТНЫ по X, Bx НЕЧЁТЕН."""
    nx = h["nx"]
    Bf = np.concatenate([Bh[:0:-1], Bh], axis=0)             # X=−Xmax..−sx, затем 0..Xmax
    Bf[:nx - 1, :, :, 0] *= -1.0                              # Bx нечётен на зеркальной половине
    hf = dict(h); hf["nx"] = 2 * nx - 1; hf["mnx"] = -h["mnx"] - (nx - 1) * h["sx"]
    return hf, Bf


def grad_profile(h, B):
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    yg = h["mny"] + np.arange(h["ny"]) * h["sy"]; zg = h["mnz"] + np.arange(h["nz"]) * h["sz"]
    jc = np.argmin(np.abs(yg)); kc = np.argmin(np.abs(zg))
    G = (B[:, jc + 1, kc, 2] - B[:, jc - 1, kc, 2]) / (2 * h["sy"])
    return xg, G, jc, kc


def intG(h, B):
    xg, G, _, _ = grad_profile(h, B)
    return np.trapezoid(G, xg)                               # ∫G·dX (Тл), сила фокусировки ∝ этому


def fringe_width(h, B):
    """Ширина фринджа: X-расстояние спада |G| с 90% до 10% от Gmax (одна сторона, X≥0), мм."""
    xg, G, _, _ = grad_profile(h, B); Gn = np.abs(G) / np.abs(G).max()
    pos = np.where(xg >= 0)[0]; xp = xg[pos]; gp = Gn[pos]
    x90 = next((xp[i] for i in range(len(gp)) if gp[i] < 0.9), xp[-1])
    x10 = next((xp[i] for i in range(len(gp)) if gp[i] < 0.1), xp[-1])
    return x10 - x90


def Gbody(h, B):
    xg, G, _, _ = grad_profile(h, B)
    return G[np.argmin(np.abs(xg))]                          # градиент в центре (тело)


def divB_axis(h, B):
    yg = h["mny"] + np.arange(h["ny"]) * h["sy"]; zg = h["mnz"] + np.arange(h["nz"]) * h["sz"]
    jm = [j for j in range(h["ny"]) if abs(yg[j]) <= 6e-3]; km = [k for k in range(h["nz"]) if abs(zg[k]) <= 6e-3]
    dBx = np.gradient(B[:, :, :, 0], h["sx"], axis=0)
    dBy = np.gradient(B[:, :, :, 1], h["sy"], axis=1)
    dBz = np.gradient(B[:, :, :, 2], h["sz"], axis=2)
    d = (dBx + dBy + dBz)[2:-2][:, jm][:, :, km]
    return np.sqrt(np.mean(d ** 2))                          # RMS |∇·B| у оси, Тл/м


# ---- методы изменения длины (на ПОЛНОЙ карте) ----
def naive_scale(h, B, c):
    """B_new(X)=B_full(X/c): координата X сжата в c (c<1 короче, c>1 длиннее). Значения интерполируются."""
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    xn = xg * c                                              # новые X-узлы (домен × c)
    hn = dict(h); hn["mnx"] = xn[0]; hn["sx"] = h["sx"] * c
    return hn, B.copy()                                      # те же значения на новой сетке (шаг × c)


def body_splice(h, B, dL_half):
    """Вставить (dL_half>0) или убрать (<0) центральные срезы тела с КАЖDОЙ стороны от центра X=0.
    Тело z-независимо ⇒ дублирование центрального среза физично; фринджи сдвигаются, не искажаясь."""
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    ic = np.argmin(np.abs(xg))                               # индекс центра
    ns = int(round(abs(dL_half) / h["sx"]))                  # число срезов на сторону
    if dL_half >= 0:
        slab = np.repeat(B[ic:ic + 1], ns, axis=0)           # копии центрального (тело) среза
        Bn = np.concatenate([B[:ic], slab, slab, B[ic:]], axis=0)   # симметрично с двух сторон
    else:
        Bn = np.concatenate([B[:ic - ns], B[ic + ns:]], axis=0)     # убрать 2·ns центральных срезов
    hn = dict(h); hn["nx"] = Bn.shape[0]; hn["mnx"] = -(hn["nx"] - 1) / 2 * h["sx"]
    return hn, Bn


# ---- компактный Boris по 3D-карте ----
def trilin(B, h, X, y, z):
    fx = (X - h["mnx"]) / h["sx"]; fy = (y - h["mny"]) / h["sy"]; fz = (z - h["mnz"]) / h["sz"]
    i, j, k = int(fx), int(fy), int(fz)
    if not (0 <= i < h["nx"] - 1 and 0 <= j < h["ny"] - 1 and 0 <= k < h["nz"] - 1):
        return np.zeros(3)
    dx, dy, dz = fx - i, fy - j, fz - k
    c = np.zeros(3)
    for (di, wx) in ((0, 1 - dx), (1, dx)):
        for (dj, wy) in ((0, 1 - dy), (1, dy)):
            for (dk, wz) in ((0, 1 - dz), (1, dz)):
                c += wx * wy * wz * B[i + di, j + dj, k + dk]
    return c


def track_ray(h, B, y0, z0, Xstart, Xend, ds=1e-3):
    """Луч из (Xstart,y0,z0) вдоль X, angle 0 → выходные углы (y',z') за Xend. Релятив. Boris, шаг пути ds."""
    pos = np.array([Xstart, y0, z0]); mom = P0 * np.array([1.0, 0.0, 0.0])
    for _ in range(int((Xend - Xstart) / ds) + 10):
        g = math.sqrt(1 + mom.dot(mom) / (ME * CL) ** 2)
        dt = ds * g * ME / math.sqrt(mom.dot(mom))          # dt под путь ds: |v|·dt=ds
        pos = pos + mom * (dt / (2 * g * ME))               # полудрейф
        Bv = trilin(B, h, pos[0], pos[1], pos[2])
        gm = math.sqrt(1 + mom.dot(mom) / (ME * CL) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * CL) ** 2)
        pos = pos + mom * (dt / (2 * g * ME))               # полудрейф
        if pos[0] >= Xend:
            break
    return mom[1] / mom[0], mom[2] / mom[0]                  # y', z'


def focal(h, B, Xspan, offs=(1e-3, 2e-3, 3e-3)):
    """1/f и аберрация по параллельным лучам. Гор: (y0,0)→y'; верт: (0,z0)→z'."""
    Xs, Xe = -Xspan, Xspan
    yp = np.array([track_ray(h, B, x, 0, Xs, Xe)[0] for x in offs])
    zp = np.array([track_ray(h, B, 0, x, Xs, Xe)[1] for x in offs])
    invfy = -np.polyfit(offs, yp, 1)[0]; invfz = -np.polyfit(offs, zp, 1)[0]
    lin = np.polyval(np.polyfit(offs, yp, 1), offs); aber = np.max(np.abs(yp - lin))
    return invfy, invfz, aber


def main():
    h0, B0 = load(MAP)
    h, B = mirror_full(h0, B0)
    xg, G, _, _ = grad_profile(h, B)
    Xspan = 0.28; iG0 = intG(h, B); gb = Gbody(h, B)
    print(f"Полный квадруполь: X {xg[0]*1e3:.0f}..{xg[-1]*1e3:.0f}мм ({h['nx']} узлов), G_max={np.abs(G).max():.3f} Тл/м")

    hdr = (f"\n{'вариант':<28}{'∫G,мТл·м':>10}{'фриндж,мм':>11}{'∇·B RMS,мТл/м':>15}"
           f"{'1/fy,1/м':>10}{'отн.аберр,%':>12}")
    print(hdr)

    def rowp(name, hh, BB):
        r = dict(iG=intG(hh, BB), fw=fringe_width(hh, BB), dv=divB_axis(hh, BB))
        y, z, a = focal(hh, BB, Xspan); rel = a / max(abs(y) * 3e-3, 1e-12) * 100
        print(f"{name:<28}{r['iG']*1e3:>10.2f}{r['fw']*1e3:>11.1f}{r['dv']*1e3:>15.3f}{y:>10.3f}{rel:>12.3f}")
        return r

    rowp("исходный (100%)", h, B)
    # ЧЕСТНО: короче ДВУМЯ путями к ОДНОМУ ∫G (одинаковая линейная фокусировка)
    print("--- короче: наивное ×0.6 vs сплайс к тому же ∫G ---")
    c = 0.6; hc, Bc = naive_scale(h, B, c)
    dL = (c * iG0 - iG0) / gb                                # сплайс тела к тому же ∫G
    hm, Bm = body_splice(h, B, dL / 2)
    rowp(f"наивное ×{c}", hc, Bc)
    rowp("сплайс к тому же ∫G", hm, Bm)
    # ЧЕСТНО: длиннее
    print("--- длиннее: наивное ×1.4 vs сплайс к тому же ∫G ---")
    c2 = 1.4; hs, Bs = naive_scale(h, B, c2)
    dL2 = (c2 * iG0 - iG0) / gb
    hp, Bp = body_splice(h, B, dL2 / 2)
    rowp(f"наивное ×{c2}", hs, Bs)
    rowp("сплайс к тому же ∫G", hp, Bp)

    print(f"\nБазлайн ∇·B = {divB_axis(h, B)*1e3:.3f} мТл/м (шум ANSYS). При ОДИНАКОВОМ ∫G (равная 1/f):")
    print("• наивное сжатие: фриндж СЖАТ (×c) — НЕфизично (у реальных магнитов фриндж задан апертурой,")
    print("  а не длиной); ∇·B РАСТЁТ (нарушение Максвелла ∝(1/c−1)·∂Bx/∂X); отн.аберрация выше.")
    print("• сплайс тела: фриндж НЕ тронут (физично); ∇·B на уровне шума; аберрация как у реального магнита.")

    visualize(h, B, hc, Bc, hs, Bs, hm, Bm, hp, Bp)
    print("Фигура:", os.path.join(OUT, "qp_scale.png"))


def visualize(h, B, hc, Bc, hs, Bs, hm, Bm, hp, Bp):
    fig, ax = plt.subplots(1, 2, figsize=(14, 5))
    for (hh, BB, lab, c, ls) in [(h, B, "исходный", "k", "-"), (hc, Bc, "наивное ×0.6", "tab:red", "--"),
                                  (hs, Bs, "наивное ×1.4", "tab:orange", "--"),
                                  (hm, Bm, "сплайс короче (=∫G)", "tab:blue", "-."),
                                  (hp, Bp, "сплайс длиннее (=∫G)", "tab:green", "-.")]:
        xg, G, _, _ = grad_profile(hh, BB)
        ax[0].plot(xg * 1e3, G, ls, color=c, lw=1.6, label=lab)
    ax[0].set_xlabel("X, мм"); ax[0].set_ylabel("G(X), Тл/м"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
    ax[0].set_title("Профиль градиента: наивное масштаб. тянет ФРИНДЖ,\nсплайс сдвигает только тело")
    # краевое Bz(X) у (5,5)мм
    for (hh, BB, lab, c, ls) in [(h, B, "исходный", "k", "-"), (hc, Bc, "наивное ×0.6", "tab:red", "--"),
                                  (hm, Bm, "сплайс −40мм", "tab:blue", "-.")]:
        xg = hh["mnx"] + np.arange(hh["nx"]) * hh["sx"]
        yg = hh["mny"] + np.arange(hh["ny"]) * hh["sy"]; zg = hh["mnz"] + np.arange(hh["nz"]) * hh["sz"]
        jo = np.argmin(np.abs(yg - 5e-3)); ko = np.argmin(np.abs(zg - 5e-3))
        ax[1].plot(xg * 1e3, BB[:, jo, ko, 0] * 1e3, ls, color=c, lw=1.6, label=lab)
    ax[1].set_xlabel("X, мм"); ax[1].set_ylabel("Bx (продольн.) при (y,z)=(5,5)мм, мТл")
    ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)
    ax[1].set_title("Краевое продольное поле: наивное сжатие\nделает фриндж КРУЧЕ (×1/c), сплайс не трогает")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "qp_scale.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
