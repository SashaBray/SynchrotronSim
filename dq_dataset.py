# -*- coding: utf-8 -*-
"""Датасет combined-function DQ: перебор токов трёх обмоток (i1,i2,i3) по ОДЗ с шагом 50 А и
физкоэффициенты (B, G, S, O), дающие ТУ ЖЕ траекторию, что реальное поле при этих токах.

РАЗЛОЖЕНИЕ ВДОЛЬ ДУГИ (не осевое): физкоэффициенты извлекаются из ПОЛЯ вдоль РЕАЛЬНОЙ изогнутой
орбиты, а не подгонкой (подгонка октуполя вырождена — переоценивает O). Рецепт:
  1) суперпозиция полей обмоток B(r)=Σ_w (I_w/nom_w)·map_w(r);
  2) опорная дуга y_ref(X) — трек осевой частицы этим полем (Boris);
  3) в каждой точке дуги раскладываем Bz по отклонению du поперёк орбиты:
       Bz = b0 + G·du + ½S·du² + (1/6)O·du³  → b0,G=∂Bz,S=∂²Bz,O=∂³Bz (полином. фит по ±4мм);
  4) эффективные коэф. = интеграл профиля вдоль дуги / L_eff: B=∫b0/L, G=∫G/L, S=∫S/L, O=∫O/L.
Это физичные коэф. вдоль орбиты (октуполь не раздувается). Октуполь O добавлен.

Контроль: за элементом 5-м пустой промежуток; секторным киком (мультиполи на du от РЕАЛЬНОЙ орбиты)
считаем невязку хвостов vs диффур — насколько 4 коэф. воспроизводят траекторию.
Выход: workspace/kickdoc/dq_dataset.csv. Запуск: py -3.14 dq_dataset.py"""
import os, csv
import numpy as np
from scipy.optimize import least_squares
import matrix_tracker as mt
import dq_sector_profile as dsp
import kickdoc_compare as kd
from bip_matrix import QP

HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
EFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}
LEFF = EFF["DQ_I1_1A"]
DS = 1e-3                                     # шаг Boris
DRIFT = 5.0                                   # 5-м пустой промежуток за элементом
STEP = 50
GRID = list(range(-100, 101, STEP))          # -100,-50,0,50,100
N = 40                                        # киков
NXS = 200                                     # X-срезов вдоль дуги для профиля
UFIT = np.linspace(-4e-3, 4e-3, 9)           # ±4мм, полином. фит поля поперёк орбиты
PROBES = [0.0, +2e-3, -2e-3, +4e-3, -4e-3]   # probe для невязки хвостов
OUT = os.path.join("workspace", "research", "kickdoc", "dq_dataset.csv")

_M = {}
for d in ("DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A"):
    h, B, nom = dsp._load_full(d, HW[d])
    _M[d] = (dsp._field_fn(h, B), nom, float(h["mxx"]))
HALF = _M["DQ_I1_1A"][2]
(f1, n1, _), (f2, n2, _), (f3, n3, _) = _M["DQ_I1_1A"], _M["DQ_I2_1A"], _M["DQ_I3_1A"]
DS_G = 2.0 * HALF / N; DS_F = LEFF / N


def comb_field(i1, i2, i3):
    w1, w2, w3 = i1 / n1, i2 / n2, i3 / n3
    return lambda p: w1 * f1(p) + w2 * f2(p) + w3 * f3(p)


def orbit(fld):
    """Опорная дуга: на оси Boris → (X, u_orbit) продольно."""
    S, U, V, up, vp = kd.boris(fld, HALF, 0.0, 0.0, DS)
    return S - HALF, U                                       # X (продольно), u_ref(X)


def arc_coeffs(fld, Xo, Uo):
    """Разложение поля ВДОЛЬ дуги: b0,G,S,O в каждой точке (полином. фит по du), интеграл/L_eff."""
    Xg = np.linspace(-HALF + 1e-4, HALF - 1e-4, NXS)
    Uc = np.interp(Xg, Xo, Uo)
    b0 = np.empty(NXS); G = np.empty(NXS); Sf = np.empty(NXS); O = np.empty(NXS)
    for i in range(NXS):
        X, yc = Xg[i], Uc[i]
        bz = np.array([fld(np.array([X, yc + u, 0.0]))[2] for u in UFIT])
        c = np.polyfit(UFIT, bz, 4)[::-1]                    # c0..c4 (по возр. степени)
        b0[i] = c[0]; G[i] = c[1]; Sf[i] = 2 * c[2]; O[i] = 6 * c[3]
    # B со ЗНАКОМ «−»: кик даёт +θ0=QP·B·L, а изгиб поля up0=−QP·∫b0 ⇒ B=−∫b0/L (у G,S,O знака нет,
    # т.к. их кики уже с минусом −k·x). См. recalibrate_coupling: per_B=+u'0/(QP·I·L).
    return (-np.trapezoid(b0, Xg) / LEFF, np.trapezoid(G, Xg) / LEFF,
            np.trapezoid(Sf, Xg) / LEFF, np.trapezoid(O, Xg) / LEFF)


def diff_probe(fld, u0):
    _, U, V, up, vp = kd.boris(fld, HALF, u0, 0.0, DS)
    st5 = mt._drift((U[-1], up, V[-1], vp, 0.0, 0.0), DRIFT)
    return U[-1], up, st5[0]


def kick_sector(B, G, S, O, u0, xref):
    """Секторный кик: мультиполи на du=x−xref(s) от РЕАЛЬНОЙ орбиты + 5-м дрейф."""
    th = QP * B * LEFF; k1 = QP * G; k2 = QP * S; k3 = QP * O
    st = (u0, 0.0, 0.0, 0.0, 0.0, 0.0)
    for i in range(N):
        st = mt._drift(st, DS_G / 2.0)
        x, px, y, py, t, d = st
        du = x - xref[i]
        px += th / N - k1 * du * DS_F - 0.5 * k2 * (du * du - y * y) * DS_F - (k3 / 6.0) * (du**3 - 3 * du * y * y) * DS_F
        py += k1 * y * DS_F + k2 * du * y * DS_F + (k3 / 6.0) * (3 * du * du * y - y**3) * DS_F
        st = mt._drift((x, px, y, py, t, d), DS_G / 2.0)
    st5 = mt._drift(st, DRIFT)
    return st[0], st5[0]


def fit_arc(fld, Xo, Uo, x0):
    """Подгонка (B,G,S,O) в СЕКТОРНОМ кике (мультиполи на du от РЕАЛЬНОЙ орбиты) под хвосты диффура.
    Старт x0 — из полевого разложения вдоль дуги (физичная инициализация, октуполь не раздувается)."""
    xk = np.array([-HALF + (i + 0.5) * DS_G for i in range(N)])
    xref = np.interp(xk, Xo, Uo)                             # реальная орбита в точках кика
    diffs = [diff_probe(fld, u0) for u0 in PROBES]

    def resid(p):
        B, G, S, O = p; r = []
        for u0, (ue_d, _, U5_d) in zip(PROBES, diffs):
            try:
                ue_k, U5_k = kick_sector(B, G, S, O, u0, xref)
                r += [(ue_k - ue_d) * 1e3, (U5_k - U5_d) * 1e3]
            except mt.ParticleLost:
                r += [1e6, 1e6]
        return r

    sol = least_squares(resid, x0, x_scale=[1e-3, 1e-1, 1e0, 1e2], method="trf")
    rms = float(np.sqrt(np.mean(np.array(sol.fun) ** 2)))
    return sol.x, rms


def main():
    rows = []
    print(f"DQ-датасет (разложение ВДОЛЬ ДУГИ + октуполь): {GRID} А ({len(GRID)**3} комб.), 5-м хвост, "
          f"L_eff={LEFF:.4f} м")
    print(f"{'i1':>5}{'i2':>5}{'i3':>5} | {'B[Тл]':>11}{'G[Тл/м]':>11}{'S[Тл/м²]':>11}{'O[Тл/м³]':>11} | "
          f"{'θ0мрад':>8}{'k1':>7}{'k2':>8}{'k3':>8} {'нев[мм]':>8}")
    for i1 in GRID:
        for i2 in GRID:
            for i3 in GRID:
                fld = comb_field(i1, i2, i3)
                Xo, Uo = orbit(fld)
                x0 = arc_coeffs(fld, Xo, Uo)                 # полевое разложение вдоль дуги — старт
                (B, G, S, O), rms = fit_arc(fld, Xo, Uo, list(x0))
                th0 = QP * B * LEFF * 1e3; k1 = QP * G; k2 = QP * S; k3 = QP * O
                rows.append([i1, i2, i3, f"{B:.6e}", f"{G:.6e}", f"{S:.6e}", f"{O:.6e}",
                             f"{th0:.5f}", f"{k1:.5f}", f"{k2:.4f}", f"{k3:.3f}", f"{rms:.4f}"])
                if (i2 == 0 and i3 == 0) or (i1 in (-100, 100) and i2 == i3) or (i1 == 100 and i2 == -100 and i3 == -100):
                    print(f"{i1:5d}{i2:5d}{i3:5d} | {B:+.4e}{G:+.4e}{S:+.4e}{O:+.4e} | "
                          f"{th0:+8.2f}{k1:+7.2f}{k2:+8.2f}{k3:+8.1f} {rms:8.3f}")
    with open(OUT, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["i1_A", "i2_A", "i3_A", "B_Tl", "G_Tlm", "S_Tlm2", "O_Tlm3",
                    "theta0_mrad", "k1_m2", "k2_m3", "k3_m4", "fit_rms_mm"])
        w.writerows(rows)
    rmsv = np.array([float(r[-1]) for r in rows])
    print(f"\n[OK] {len(rows)} строк → {OUT}")
    print(f"невязка хвостов: медиана {np.median(rmsv):.3f} мм, 90%={np.percentile(rmsv,90):.3f}, "
          f"макс {rmsv.max():.3f}; доля <1мм {100*np.mean(rmsv<1):.0f}%")


if __name__ == "__main__":
    main()
