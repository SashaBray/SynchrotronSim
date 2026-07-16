# -*- coding: utf-8 -*-
"""Датасет DQ по СИММЕТРИЧНО-АРОЧНОМУ алгоритму (7 шагов):
  1) при токах (i1,i2,i3) пускаем пучок сквозь элемент (Boris по суперпозиции полей обмоток);
  2) смотрим, на сколько повернул (изгиб θ);
  3) ПОВОРАЧИВАЕМ ЭЛЕМЕНТ на φ, чтобы угол влёта = угол вылета (симметрия дуги);
  4-6) итерируем φ до максимальной симметрии (φ→θ/2, но θ зависит от пути ⇒ несколько итераций);
  7) ВОСПРОИЗВОДИМ опыт на матричном движке (секторный кик) и ПОДБИРАЕМ физкоэф. (θ,G,S,O)
     траекторным матчингом под диффур — 5 пробных частиц (центр + ±2/±4 мм) + 5-м хвост.
Коэф. из ТРАЕКТОРИИ на проектной дуге, а не из поля вокруг прямой оси (как recalibrate) — это и есть
рабочая точка изгибающего элемента. Для контраста считаем и ОСЕВОЙ отбор (recalibrate: прямые пробники).

«Поворот элемента на φ» ≡ влёт пучка под −φ в системе элемента (поле сэмплируется вдоль наклонённого
пути). Симметрия: наклон на выходе = −наклон на входе.
Выход: workspace/research/kickdoc/dq_arc_vs_axis.csv + figs/dq_arc_vs_axis.png
Запуск: py -3.14 dq_arc_dataset.py"""
import os, csv, math, sys
import numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matrix_tracker as mt
import dq_sector_profile as dsp
from bip_matrix import QP

HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
EFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}
# выбор длины DQ: --dq L1000 / L400 / (пусто=база). Peek ДО argparse — карты грузятся на уровне модуля.
_SUF = next((("_" + sys.argv[i + 1].lstrip("_")) for i, a in enumerate(sys.argv) if a == "--dq" and i + 1 < len(sys.argv)), "")
DEV = tuple(f"DQ_I{k}_1A{_SUF}" for k in (1, 2, 3))
DESIGN = dict(theta=0.0247680586728496, K1=-1.960219550894432) if "1000" in _SUF \
    else dict(theta=0.00523006, K1=-1.558437304101693)     # .lte: длинный DQ1B (L=1) или короткий DQ2C (L≈0.4/база)
LEFF = EFF[DEV[0]]
DS = 1e-3                                          # шаг Boris
A = 2e-3                                           # ±плечо ОСЕВОГО прямого пробника (recalibrate, гориз.)
N = 40                                            # срезов секторного кика
RBEAM = 1e-3                                        # радиус пучка (диаметр 2 мм)
PROBES = [(0.0, 0.0)] + [(RBEAM * math.cos(k * math.pi / 6.0), RBEAM * math.sin(k * math.pi / 6.0))
                         for k in range(12)]         # 13 частиц: центр + 12 по кругу через 30° в плоскости (u,v)
DRIFT = 5.0                                        # 5-м хвост за элементом
OUT = os.path.join("workspace", "research", "kickdoc", "dq_arc_vs_axis.csv")
FIG = os.path.join("workspace", "research", "kickdoc", "figs", "dq_arc_vs_axis.png")

_M = {}
for d in DEV:
    h, B, nom = dsp._load_full(d, HW[d])
    _M[d] = (dsp._field_fn(h, B), nom, float(h["mxx"]))
HALF = _M[DEV[0]][2]
_h1 = dsp._load_full(DEV[0], HW[DEV[0]])[0]
YMAX = min(abs(_h1["mny"]), abs(_h1["mny"] + (_h1["ny"] - 1) * _h1["sy"]))   # гор. предел карты (поле=0 вне)
ZMAX = min(abs(_h1["mnz"]), abs(_h1["mnz"] + (_h1["nz"] - 1) * _h1["sz"]))   # верт. предел карты
RMAX = min(YMAX, ZMAX)                                                        # предел по |попереч| (пучок 2D)
(f1, n1, _), (f2, n2, _), (f3, n3, _) = _M[DEV[0]], _M[DEV[1]], _M[DEV[2]]
GDES = DESIGN["K1"] / QP
LPHYS = {r["device_type_id"]: float(r["phys_length"]) for r in mt._read("hardware_library.csv")}[DEV[0]]
LEAD = max(0.0, (2.0 * HALF - LPHYS) / 2.0)       # поле-свободный отступ: физ.элемент L_phys по центру зоны 2·half
DS_G = LPHYS / N; DS_F = LEFF / N                  # ДРЕЙФ по L_phys (как ring-ядро element_bench), КИК по Leff — НЕ 2·half!


class OutOfMap(Exception):
    """Пучок вышел за поперечную границу карты (поле там 0 — данные недостоверны)."""


def comb(i1, i2, i3):
    w1, w2, w3 = i1 / n1, i2 / n2, i3 / n3
    return lambda p: w1 * f1(p) + w2 * f2(p) + w3 * f3(p)


def boris(fld, u0=0.0, v0=0.0, up0=0.0, vp0=0.0, ds=DS):
    """Boris по полю, влёт (u0,v0,up0,vp0). → X[], u[], v[], u'вых, v'вых, макс|попереч| (для guard)."""
    x = -HALF; u, v, up, vp = u0, v0, up0, vp0
    n = int(round(2 * HALF / ds)); Xs = [x]; Us = [u]; Vs = [v]; rmax = max(abs(u0), abs(v0))
    for _ in range(n):
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        B = fld(np.array([x, u, v])); up += -QP * B[2] * ds; vp += QP * B[1] * ds
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        Xs.append(x); Us.append(u); Vs.append(v)
        r = abs(u) if abs(u) > abs(v) else abs(v)
        if r > rmax:
            rmax = r
    return np.array(Xs), np.array(Us), np.array(Vs), up, vp, rmax


def bend_axis_integral(fld):
    """θ по интегралу поля вдоль ПРЯМОЙ оси u=0 (всегда в карте) — робастный seed для симметрии при больших токах."""
    Xg = np.linspace(-HALF + 1e-4, HALF - 1e-4, int(2 * HALF / DS))
    bz = np.array([fld(np.array([X, 0.0, 0.0]))[2] for X in Xg])
    return -QP * np.trapezoid(bz, Xg)


def boris_probe(fld, u0, v0, up0, vp0):
    """Пробная частица: влёт (u0,v0,up0,vp0). → uвых,vвых, u',v', хвосты u(+5м),v(+5м), макс|попереч|."""
    _, Us, Vs, up, vp, rmax = boris(fld, u0, v0, up0, vp0)
    uend, vend = Us[-1], Vs[-1]
    return uend, vend, up, vp, uend + up * DRIFT, vend + vp * DRIFT, rmax


def axis_coeffs(fld):
    """ОСЬ (recalibrate): прямые ГОРИЗ. пробники наклон=0 на u0=0,±A. → θ, G, S (NaN, если вышел за карту).
    По построению гориз.-пробниковый (симметрия u'(±a)) — как recalibrate_coupling."""
    _, _, _, up0, _, m0 = boris(fld, 0.0, 0.0, 0.0, 0.0)
    _, _, _, upp, _, mp = boris(fld, +A, 0.0, 0.0, 0.0); _, _, _, upm, _, mm = boris(fld, -A, 0.0, 0.0, 0.0)
    if max(m0, mp, mm) > RMAX:
        return float("nan"), float("nan"), float("nan")
    odd = (upp - upm) / 2.0; even = (upp + upm) / 2.0
    return up0, -odd / (QP * A * LEFF), -2.0 * (even - up0) / (QP * A * A * LEFF)


def symmetrize(fld, verbose=False):
    """Шаги 3-6: итеративно крутим элемент на φ, пока наклон вых = −наклон вх. → φ, θ_дуги, история.
    Старт: 0 если прямой пучок в карте (наглядно), иначе seed по интегралу поля (без блужданий за картой)."""
    _, _, _, _, _, m0 = boris(fld, 0.0, 0.0, 0.0, 0.0)
    phi = 0.0 if m0 < RMAX else bend_axis_integral(fld) / 2.0
    hist = []; theta = 0.0
    for k in range(8):
        _, _, _, up, _, _ = boris(fld, 0.0, 0.0, -phi, 0.0)   # центр.частица, влёт −φ
        theta = up + phi                                # полный изгиб = вых − вх = up − (−φ)
        asym = up - phi                                 # симметрия: вых(up) = −вх(−φ) = +φ  ⇒ asym→0
        hist.append((phi, theta, asym))
        if abs(asym) < 1e-7:
            break
        phi = theta / 2.0                               # новый угол поворота
    return phi, theta, hist


def kick_sym(theta, G, S, O, u0, v0, xref, phi):
    """МАТРИЧНЫЙ движок (секторный кик) на симметр. дуге: влёт (u0,v0) под −φ; мультиполи на du=x−xref и y=v.
    → uвых,vвых, u',v', хвосты u,v."""
    k1 = QP * G; k2 = QP * S; k3 = QP * O
    st = mt._drift((u0, -phi, v0, 0.0, 0.0, 0.0), LEAD)     # поле-свободный отступ (влёт −φ, верт. v0)
    for i in range(N):
        st = mt._drift(st, DS_G / 2.0)
        x, px, y, py, t, d = st
        du = x - xref[i]
        px += theta / N - k1 * du * DS_F - 0.5 * k2 * (du * du - y * y) * DS_F - (k3 / 6.0) * (du**3 - 3 * du * y * y) * DS_F
        py += k1 * y * DS_F + k2 * du * y * DS_F + (k3 / 6.0) * (3 * du * du * y - y**3) * DS_F
        st = mt._drift((x, px, y, py, t, d), DS_G / 2.0)
    st = mt._drift(st, LEAD)                                 # поле-свободный отступ на выходе
    uend, vend, up, vp = st[0], st[2], st[1], st[3]
    return uend, vend, up, vp, uend + up * DRIFT, vend + vp * DRIFT


def fit_matrix(fld, phi, x0):
    """Шаг 7: подгоняем (θ,G,S,O) секторного кика под траектории диффура 13 частиц (центр+круг Ø2мм) на дуге."""
    Xs, Us, Vs, _, _, orb_m = boris(fld, 0.0, 0.0, -phi, 0.0)   # симметр. опорная орбита (центр.частица)
    if orb_m > RMAX:
        raise OutOfMap                                          # дуга вышла за карту — данные недостоверны
    xk = np.array([-LPHYS / 2.0 + (i + 0.5) * DS_G for i in range(N)])   # точки кика — внутри центрального L_phys
    xref = np.interp(xk, Xs, Us)
    refs = [boris_probe(fld, u0, v0, -phi, 0.0) for (u0, v0) in PROBES]  # 13 частиц влёт под −φ
    if max(r[6] for r in refs) > RMAX:                          # какой-то пробник вышел за карту
        raise OutOfMap

    def resid(p):
        theta, G, S = p; r = []                             # октуполь 2-мм пучком не определяется (∝r³) — не фитим
        for (u0, v0), (ue, ve, upe, vpe, ut, vt, _) in zip(PROBES, refs):
            uk, vk, upk, vpk, utk, vtk = kick_sym(theta, G, S, 0.0, u0, v0, xref, phi)
            r += [(uk - ue) * 1e3, (vk - ve) * 1e3, (utk - ut) * 1e3, (vtk - vt) * 1e3]   # гор.+верт. выход и хвост
        return r

    sol = least_squares(resid, x0[:3], x_scale=[1e-3, 1e-2, 1e0], method="trf")
    rms = float(np.sqrt(np.mean(np.array(sol.fun) ** 2)))
    return np.append(sol.x, 0.0), rms, sagitta(Xs, Us)      # (θ,G,S,O=0) — O не определяется 2-мм пучком


def sagitta(Xo, Uo):
    chord = Uo[0] + (Uo[-1] - Uo[0]) * (Xo - Xo[0]) / (Xo[-1] - Xo[0])
    return float(np.max(np.abs(Uo - chord)))


def analyse(i1, i2, i3, verbose=False):
    fld = comb(i1, i2, i3)
    th_ax, G_ax, S_ax = axis_coeffs(fld)                # осевой отбор (NaN, если прямой вне карты)
    phi, th_sym, hist = symmetrize(fld, verbose)        # шаги 1-6
    g0 = G_ax if np.isfinite(G_ax) else 0.0             # init фита не должен быть NaN
    s0 = S_ax if np.isfinite(S_ax) else 0.0
    (th_m, G_m, S_m, O_m), rms, sag = fit_matrix(fld, phi, [th_sym, g0, s0, 0.0])  # шаг 7 (raise OutOfMap)
    return dict(axis=(th_ax, G_ax, S_ax), mat=(th_m, G_m, S_m, O_m),
                phi=phi, rms=rms, sag=sag, hist=hist)


def main(imax=4.0, npts=5, fix_i3=None, out_base=None):
    GRID = [round(float(x), 4) for x in np.linspace(-imax, imax, npts)]
    twod = fix_i3 is not None
    i3_grid = [round(float(fix_i3), 4)] if twod else GRID
    i3ref = round(float(fix_i3), 4) if twod else 1.0
    base = out_base or "dq_arc_vs_axis"
    out_csv = os.path.join("workspace", "research", "kickdoc", base + ".csv")
    out_fig = os.path.join("workspace", "research", "kickdoc", "figs", base + ".png")
    print(f"DQ симметрично-арочный датасет. L_eff={LEFF:.4f} м, полузона={HALF:.3f} м. "
          f"Проектно ({DEV[0]}): θ={DESIGN['theta']*1e3:.2f} мрад, K1={DESIGN['K1']:.3f} ⇒ G={GDES:+.4f} Тл/м\n")

    print(f"=== ТОЧКА (1,1,{i3ref:g}) А — итерации симметрии (шаги 1-6) ===")
    print(f"{'итер':>5}{'φ, мрад':>10}{'θ дуги, мрад':>14}{'асимметрия, мкрад':>18}")
    r = analyse(1.0, 1.0, i3ref, verbose=True)
    for k, (phi, th, asym) in enumerate(r["hist"]):
        print(f"{k:>5}{phi*1e3:>10.4f}{th*1e3:>14.4f}{asym*1e6:>18.2f}")
    print(f"→ поворот элемента φ={r['phi']*1e3:.4f} мрад")

    print(f"\n=== КОЭФФИЦИЕНТЫ в точке (1,1,{i3ref:g}) А ===")
    print(f"{'метод':30}{'θ, мрад':>10}{'G, Тл/м':>10}{'S, Тл/м²':>11}")
    th, G, S = r["axis"]
    print(f"{'ОСЬ (прямой пробник)':30}{th*1e3:>10.3f}{G:>10.4f}{S:>11.3f}")
    th, G, S, O = r["mat"]
    print(f"{'ДУГА-СИММ (матч. движка)':30}{th*1e3:>10.3f}{G:>10.4f}{S:>11.3f}")
    print(f"невязка матчинга {r['rms']:.4f} мм\n")

    mode = f"2D-СВИП i1,i2 (i3={fix_i3:g} А фикс.)" if twod else "3D-СВИП по обмоткам"
    print(f"=== {mode} ∈ {GRID} А ({len(GRID) * len(GRID) * len(i3_grid)} комб.) ===")
    rows = []; lost = 0
    for i1 in GRID:
        for i2 in GRID:
            for i3 in i3_grid:
                try:
                    a = analyse(i1, i2, i3)                             # OutOfMap если дуга за картой
                    if not np.isfinite(a["rms"]):
                        raise ValueError
                except Exception:
                    lost += 1; continue
                th_ax, G_ax, S_ax = a["axis"]; th_m, G_m, S_m, O_m = a["mat"]
                rows.append([f"{i1:.3f}", f"{i2:.3f}", f"{i3:.3f}",
                             f"{th_ax*1e3:.5f}", f"{th_m*1e3:.5f}", f"{G_ax:.6f}", f"{G_m:.6f}",
                             f"{S_ax:.5f}", f"{S_m:.5f}", f"{O_m:.4f}",
                             f"{a['phi']*1e3:.5f}", f"{a['sag']*1e3:.4f}", f"{a['rms']:.5f}"])

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    with open(out_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["i1_A", "i2_A", "i3_A", "theta_axis_mrad", "theta_arc_mrad",
                    "G_axis_Tlm", "G_arc_Tlm", "S_axis_Tlm2", "S_arc_Tlm2", "O_arc_Tlm3",
                    "phi_mrad", "sagitta_mm", "fit_rms_mm"])
        w.writerows(rows)

    rmsv = np.array([float(r[-1]) for r in rows])
    print(f"\nневязка матчинга движок↔диффур: медиана {np.median(rmsv)*1e3:.1f} мкм, "
          f"90%={np.percentile(rmsv,90)*1e3:.1f}, макс {rmsv.max()*1e3:.1f} мкм")
    print(f"надёжность (θ,G,S воспроизводят диффур): <0.1мм — {100*np.mean(rmsv<0.1):.0f}% строк, "
          f"<1мм — {100*np.mean(rmsv<1.0):.0f}%")

    os.makedirs(os.path.dirname(out_fig), exist_ok=True)
    if twod:                                                    # 2D тепловые карты θ_дуга, G_дуга над (i1, i2)
        gi = sorted(set(float(x[0]) for x in rows)); gj = sorted(set(float(x[1]) for x in rows))
        ii = {v: k for k, v in enumerate(gi)}; jj = {v: k for k, v in enumerate(gj)}
        Zth = np.full((len(gj), len(gi)), np.nan); ZG = np.full((len(gj), len(gi)), np.nan)
        for x in rows:
            Zth[jj[float(x[1])], ii[float(x[0])]] = float(x[4]); ZG[jj[float(x[1])], ii[float(x[0])]] = float(x[6])
        fig, ax = plt.subplots(1, 2, figsize=(10.5, 4.3))
        c0 = ax[0].pcolormesh(gi, gj, Zth, cmap="RdBu_r", shading="nearest"); fig.colorbar(c0, ax=ax[0]).set_label("θ_дуга, мрад")
        ax[0].set_title(f"Изгиб θ (дуга), i3={fix_i3:g} А"); ax[0].set_xlabel("i1, А"); ax[0].set_ylabel("i2, А")
        c1 = ax[1].pcolormesh(gi, gj, ZG, cmap="RdBu_r", shading="nearest"); fig.colorbar(c1, ax=ax[1]).set_label("G_дуга, Тл/м")
        ax[1].set_title("Градиент G (дуга)"); ax[1].set_xlabel("i1, А"); ax[1].set_ylabel("i2, А")
        fig.tight_layout(); fig.savefig(out_fig, dpi=130); plt.close(fig)
        print(f"\n[OK] {len(rows)} строк ({lost} потеряно) → {out_csv}")
        print(f"     2D-карты θ,G(i1,i2) → {out_fig}")
    else:                                                       # 3D: диагональ i1=i2=i3
        diag = sorted([[float(x[0]), float(x[3]), float(x[4]), float(x[5]), float(x[6])]
                       for x in rows if x[0] == x[1] == x[2]], key=lambda z: z[0])
        print(f"\nдиагональ i1=i2=i3:  {'ток':>6}{'θ_ось':>9}{'θ_дуга':>9}{'Δθ%':>7}{'G_ось':>9}{'G_дуга':>9}")
        for s, tax, tarc, gax, garc in diag:
            ok = math.isfinite(tax) and tax != 0
            dpc = (abs(tarc) - abs(tax)) / abs(tax) * 100 if ok else float("nan")
            sax = f"{tax:>9.3f}" if math.isfinite(tax) else f"{'—(вне)':>9}"
            sgax = f"{gax:>9.4f}" if math.isfinite(gax) else f"{'—':>9}"
            sdpc = f"{dpc:>7.1f}" if math.isfinite(dpc) else f"{'—':>7}"
            print(f"{'':18}{s:>6.2f}{sax}{tarc:>9.3f}{sdpc}{sgax}{garc:>9.4f}")
        if diag:
            d = np.array(diag)
            fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
            ax[0].plot(d[:, 0], d[:, 1], "o-", color="#e15759", label="ОСЬ"); ax[0].plot(d[:, 0], d[:, 2], "s-", color="#2f9e44", label="ДУГА")
            ax[0].set_xlabel("i1=i2=i3, А"); ax[0].set_ylabel("θ, мрад"); ax[0].set_title("Изгиб"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
            ax[1].plot(d[:, 0], d[:, 3], "o-", color="#e15759", label="ОСЬ"); ax[1].plot(d[:, 0], d[:, 4], "s-", color="#2f9e44", label="ДУГА")
            ax[1].axhline(GDES, ls="--", color="#555", lw=1, label="проект")
            ax[1].set_xlabel("i1=i2=i3, А"); ax[1].set_ylabel("G, Тл/м"); ax[1].set_title("Квадруполь"); ax[1].grid(alpha=.3); ax[1].legend(fontsize=8)
            fig.tight_layout(); fig.savefig(out_fig, dpi=130); plt.close(fig)
        print(f"\n[OK] {len(rows)} строк ({lost} потеряно) → {out_csv}")
        print(f"     фигура (диагональ) → {out_fig}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="DQ симметрично-арочный датасет: ось vs дуга, 3D или 2D свип обмоток")
    ap.add_argument("--imax", type=float, default=4.0, help="макс. |ток| обмотки, А (сетка −imax..+imax)")
    ap.add_argument("--n", type=int, default=5, help="точек на обмотку (нечётное — включает 0)")
    ap.add_argument("--i3", type=float, default=None, help="зафиксировать i3 → 2D-свип только i1,i2 (без флага — 3D)")
    ap.add_argument("--out", default=None, help="имя выходного файла без расширения (по умолч. dq_arc_vs_axis)")
    ap.add_argument("--dq", default="", help="суффикс длины DQ: L1000 / L400 (пусто — базовый DQ_I1_1A)")
    a = ap.parse_args()
    main(a.imax, a.n, a.i3, a.out)
