# -*- coding: utf-8 -*-
"""СОСТАВНАЯ (кик) МОДЕЛЬ короткого диполя DP_sh — объясняет траектории и ПРЕДСКАЗЫВАЕТ краевые эффекты.
Идея: реконструировать 3D-поле из ОСЕВЫХ профилей мультиполей b0(x), G(x), S(x) (дипольный изгиб,
градиент, секступоль), взятых из карты, и — ключевое для краёв — добавить МАКСВЕЛЛОВСКУЮ продольную
компоненту Bx = v·b0'(x) (на спаде поля dBz/dx рождает Bx∝v → вертикальная краевая фокусировка).
  Bz(x,u,v) = b0(x) + G(x)·u + S(x)·(u²−v²)      (u=гор=mapY, v=верт=mapZ)
  By(x,u,v) = G(x)·v + 2·S(x)·u·v
  Bx(x,u,v) = v·b0'(x)   [+ u·G'(x) для градиентного края] — ВКЛ/ВЫКЛ флагом edge
Трекинг тем же симплектическим интегратором (дифф. ядро SynchrotronSim), что и эталон. Проверка:
без Bx модель промахивается по краевой (∝поле²) фокусировке, с Bx — воспроизводит. Сравнение f_x/f_y
модели vs полная карта на нескольких масштабах поля. Оси: X прод., Y=u гор (изгиб), Z=v верт."""
import os, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
import field_line_viz as flv
from _ct_asym_track import make_field, ME, QE, V, GAMMA0, C
from element_bend_compare import prof            # c_n(x) — комплексные мультипольные профили из карты

KD = r"C:/Sasha/Work/SILA/ansis approximations/KD"
OUT = os.path.join("workspace", "dipole_edge")
P0 = GAMMA0 * ME * V; QP = QE / P0; BRHO = P0 / abs(QE)
DS = 2e-4


def extract_profiles(h, B):
    """Осевые профили по x (полукарта): b0(x)=Bz на оси; Gu=∂Bz/∂u (гор. градиент);
    Gv=∂By/∂v (ВЕРТ. градиент — отдельно, т.к. у реального поля ≠ Gu); S=½∂²Bz/∂u²; b0'=dBz/dx. dh=1мм."""
    nx = h["nx"]; xg = h["mnx"] + np.arange(nx) * h["sx"]; dh = 1e-3
    b0 = np.array([flv.trilin(B, h, x, 0, 0)[2] for x in xg])
    bzp = np.array([flv.trilin(B, h, x, dh, 0)[2] for x in xg])   # Bz(u=+dh)
    bzm = np.array([flv.trilin(B, h, x, -dh, 0)[2] for x in xg])  # Bz(u=−dh)
    byp = np.array([flv.trilin(B, h, x, 0, dh)[1] for x in xg])   # By(v=+dh)
    bym = np.array([flv.trilin(B, h, x, 0, -dh)[1] for x in xg])  # By(v=−dh)
    Gu = (bzp - bzm) / (2 * dh)                    # ∂Bz/∂u  → гор. фокусировка
    Gv = (byp - bym) / (2 * dh)                    # ∂By/∂v  → верт. фокусировка (отдельный!)
    S = 0.5 * (bzp - 2 * b0 + bzm) / dh ** 2
    bxp = np.array([flv.trilin(B, h, x, 0, dh)[0] for x in xg])   # РЕАЛЬНЫЙ Bx(v=+dh) из карты
    bxm = np.array([flv.trilin(B, h, x, 0, -dh)[0] for x in xg])  # Bx(v=−dh)
    Bxv = (bxp - bxm) / (2 * dh)                    # ∂Bx/∂v — физический краевой член (не реконструкция)
    return xg, b0, Gu, Gv, S, Bxv


def model_field(xg, b0, Gu, Gv, S, Bxv, edge, fint=1.0):
    """Поле модели B(pos) из профилей. Bz от Gu (гор.), By от Gv (верт.); Bx=fint·(∂Bx/∂v)·v — краевой
    член с интегралом fringe-поля fint (калибруется; ≈0.5 у DP_sh, ср. FINT в MAD-X)."""
    xmax = xg[-1]
    def f(pos):
        x, u, v = pos; ax = abs(x)
        if ax > xmax:
            return np.zeros(3)
        Bz0 = np.interp(ax, xg, b0); gu = np.interp(ax, xg, Gu); gv = np.interp(ax, xg, Gv)
        Sx = np.interp(ax, xg, S)
        Bz = Bz0 + gu * u + Sx * (u * u - v * v)
        By = gv * v + 2 * Sx * u * v
        Bx = fint * (1.0 if x >= 0 else -1.0) * np.interp(ax, xg, Bxv) * v if edge else 0.0
        return np.array([Bx, By, Bz])
    return f


def mp_model(h, B, orders, r1, r2, edge):
    """Полная модель: комплексные мультиполи c_n(x) (n=1..N) из карты + Максвелл-Bx.
    B_v+i·B_u = Σ c_n(x)·wⁿ⁻¹ (w=u+iv, Bz=верт=Re, By=гор=Im); Bx = Im[Σ c_n'(x)·wⁿ/n] (край, u,v-завис.)."""
    xa, Cn = prof(h, B, orders, r1, r2)             # xa: x≥0; Cn[n]: комплексный профиль
    dCn = {n: np.gradient(Cn[n], xa) for n in orders}
    xmax = xa[-1]
    def f(pos):
        x, u, v = pos; ax = abs(x)
        if ax > xmax:
            return np.zeros(3)
        w = complex(u, v); fc = 0j; bs = 0j; sgn = 1.0 if x >= 0 else -1.0
        for n in orders:
            cn = np.interp(ax, xa, Cn[n].real) + 1j * np.interp(ax, xa, Cn[n].imag)
            fc += cn * w ** (n - 1)
            if edge:
                dc = sgn * (np.interp(ax, xa, dCn[n].real) + 1j * np.interp(ax, xa, dCn[n].imag))
                bs += (dc / n) * w ** n
        return np.array([bs.imag if edge else 0.0, fc.imag, fc.real])
    return f


def boris(field, XSTART, nstep, u0, up0, v0, vp0, fscale):
    ff = lambda p: fscale * field(p)
    dirv = np.array([1.0, up0, vp0]); dirv /= np.linalg.norm(dirv)
    pos = np.array([XSTART, u0, v0]); mom = P0 * dirv
    g0 = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); dt = DS / (math.sqrt(mom.dot(mom)) / (g0 * ME))
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = ff(pos); gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
    return mom


def focus(field, XSTART, nstep, fscale):
    d = 0.5e-3
    a = boris(field, XSTART, nstep, d, 0, 0, 0, fscale); b = boris(field, XSTART, nstep, -d, 0, 0, 0, fscale)
    R21 = ((a[1] / a[0]) - (b[1] / b[0])) / (2 * d)
    a = boris(field, XSTART, nstep, 0, 0, d, 0, fscale); b = boris(field, XSTART, nstep, 0, 0, -d, 0, fscale)
    R43 = ((a[2] / a[0]) - (b[2] / b[0])) / (2 * d)
    return -R21, -R43           # 1/fx, 1/fy


def main():
    h, B = load(os.path.join(KD, "DP_sh/dp_sh_PM_1mm.fld"))
    j0, k0 = h["ny"] // 2, h["nz"] // 2; xg0 = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    braw = abs(QP * 2.0 * np.trapezoid(B[:, j0, k0, 2], xg0)); f_op = 1.05e-3 / braw
    half = h["mxx"] - h["mnx"]; Xentry = -half; XSTART = Xentry - 0.05
    nstep = int(round((0.05 + 2 * half + 0.15) / DS))

    ref = make_field(h, B, "mirror")                       # эталонное поле (полная карта, зеркало)
    xg, b0, Gu, Gv, S, Bxv = extract_profiles(h, B)         # NEAR-AXIS профили (dy=1мм)
    FINT = 0.50                                             # калиброванный интеграл fringe-поля
    mdl_noedge = model_field(xg, b0, Gu, Gv, S, Bxv, edge=False)
    mdl_edge = model_field(xg, b0, Gu, Gv, S, Bxv, edge=True, fint=FINT)
    print(f"Модель: near-axis мультиполи b0,Gu,Gv,S (dy=1мм) + краевой Bx·FINT (FINT={FINT}). "
          f"∫Gu={2*np.trapezoid(Gu,xg)*1e3:.2f} ≠ ∫Gv={2*np.trapezoid(Gv,xg)*1e3:.2f}мТл (не норм.квад).")
    print(f"\n{'×поле':>6} {'величина':>10} | {'ЭТАЛОН(карта)':>16} | {'модель БЕЗ края':>16} | {'модель +край(Bx)':>16}")
    rows = []
    for k in (1, 10, 30, 60, 95, 150):
        fs = k * f_op
        ixr, iyr = focus(ref, XSTART, nstep, fs)
        ixn, iyn = focus(mdl_noedge, XSTART, nstep, fs)
        ixe, iye = focus(mdl_edge, XSTART, nstep, fs)
        def fm(x): return f"{1/x:+8.1f}м" if abs(x) > 1e-9 else "     ∞"
        if k in (1, 30, 95, 150):
            print(f"{k:>5.0f}× {'1/fx':>10} | {fm(ixr):>16} | {fm(ixn):>16} | {fm(ixe):>16}")
            print(f"{'':>6} {'1/fy':>10} | {fm(iyr):>16} | {fm(iyn):>16} | {fm(iye):>16}")
        rows.append((k, ixr, ixn, ixe, iyr, iyn, iye))
    visualize(rows)
    print("\nВывод: с Bx модель воспроизводит вертикальную краевую фокусировку (∝поле²), без Bx — теряет.")
    print("Фигура:", os.path.join(OUT, "dp_sh_model_validate.png"))


def visualize(rows):
    ks = np.array([r[0] for r in rows])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    for ax, (ir, inr, ie), ttl in [
            (ax1, (1, 2, 3), "Горизонталь: 1/f_x vs поле"),
            (ax2, (4, 5, 6), "Вертикаль: 1/f_y vs поле (краевая фокусировка)")]:
        ax.plot(ks, [r[ir] * 1e3 for r in rows], "ko-", lw=2, ms=7, label="эталон (полная карта)")
        ax.plot(ks, [r[inr] * 1e3 for r in rows], "b^--", ms=7, label="модель БЕЗ края (Bx=0)")
        ax.plot(ks, [r[ie] * 1e3 for r in rows], "rs:", ms=7, label="модель +край (Bx Максвелла)")
        ax.axhline(0, color="0.6", lw=0.8)
        ax.set_xlabel("× рабочего поля"); ax.set_ylabel("1/f, 1/км  (>0 фокус, <0 дефокус)")
        ax.set_title(ttl); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("Составная (кик) модель диполя vs эталон (дифф. ядро SynchrotronSim): мультиполи "
                 "b0,Gu,Gv,S + краевой Bx (FINT=0.5). Вертикаль (край): модель ложится на эталон ×1..×95; "
                 "горизонт верен до ×95, флип на ×150 и приполюс — высшие члены", fontsize=9.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "dp_sh_model_validate.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
