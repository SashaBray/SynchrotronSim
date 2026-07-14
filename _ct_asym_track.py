# -*- coding: utf-8 -*-
"""Сравнение траектории частицы сквозь АСИММЕТРИЧНОЕ полное поле QP_L (100А) с матричной
моделью (идеальный симметричный квадруполь). Дифф. ядро SynchrotronSim (интегратор Бориса, порт C++), 3D трилин. поле."""
import os, math, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
import field_line_viz as flv

C = 299792458.0
ME = 9.1093837e-31
QE = -1.602176634e-19
V = 2.967e8
GAMMA0 = 1.0 / math.sqrt(1 - (V / C) ** 2)
FIGS = os.path.join("workspace", "bip_compare", "figs")
FLD = r"C:/Sasha/Work/SILA/ansis approximations/KD/QP_L/QP_L_100A_full_v2.fld"


def make_field(h, B, kind="full"):
    """Возвращает функцию B(pos). kind: 'full' — как есть; 'sym' — симметризованное (среднее
    половин); 'mirror' — зеркало половины x≥0."""
    ymax = h["mny"] + (h["ny"] - 1) * h["sy"]; zmax = h["mnz"] + (h["nz"] - 1) * h["sz"]
    def inside(x, y, z):
        return h["mnx"] < x < h["mxx"] and h["mny"] < y < ymax and h["mnz"] < z < zmax
    def f(pos):
        x, y, z = pos
        if not inside(x, y, z):
            if kind in ("full", "mirror", "sym") and inside(-x, y, z):
                pass
            else:
                return np.zeros(3)
        if kind == "full":
            return flv.trilin(B, h, x, y, z) if inside(x, y, z) else np.zeros(3)
        # для sym/mirror используем зеркало
        bp = flv.trilin(B, h, x, y, z) if inside(x, y, z) else None
        bm = flv.trilin(B, h, -x, y, z) if inside(-x, y, z) else None
        bm = np.array([-bm[0], bm[1], bm[2]]) if bm is not None else None  # зеркало: Bx→−Bx
        if kind == "mirror":
            return (bp if x >= 0 and bp is not None else (bm if bm is not None else np.zeros(3)))
        # sym: среднее доступных половин
        vals = [v for v in (bp, bm) if v is not None]
        return np.mean(vals, axis=0) if vals else np.zeros(3)
    return f


def boris_track(y0, field, x0=-0.40, x1=0.40, dx=5e-4):
    """Прогон одной частицы через field(pos). Старт x0, конец x1, шаг по длине dx."""
    pos = np.array([x0, y0, 0.0]); mom = GAMMA0 * ME * np.array([V, 0.0, 0.0])
    dt = dx / V
    traj = [pos.copy()]
    nstep = int((x1 - x0) / dx)
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        pos = pos + mom * (dt / (2 * g * ME))
        Bv = field(pos)
        gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv
        s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t)
        mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        pos = pos + mom * (dt / (2 * g * ME))
        traj.append(pos.copy())
    return np.array(traj), mom


def matrix_quad(y0, h, B, x0=-0.40, x1=0.40, dx=5e-4):
    """Идеальный симметричный квадруполь: интегр. градиент из поля, аналит. фокусировка."""
    # градиент на оси по x: G(x)=∂Bz/∂y, интеграл по телу
    dy = 4e-3
    xs = np.arange(h["mnx"], h["mxx"], h["sx"])
    Gx = np.array([(flv.trilin(B, h, x, dy, 0)[2] - flv.trilin(B, h, x, -dy, 0)[2]) / (2 * dy) for x in xs])
    p = GAMMA0 * ME * V; qp = QE / p
    kappa = qp * Gx                                  # κ(x) = (q/p)·G(x), 1/м²  (симметризуем)
    xg = np.arange(x0, x1, dx)
    ks = np.interp(xg, xs, kappa)
    ksym = 0.5 * (ks + ks[::-1])                     # СИММЕТРИЧНЫЙ профиль (среднее ±x)
    y, yp = y0, 0.0
    traj = [(x0, y)]
    for i in range(len(xg) - 1):
        y += yp * dx / 2                             # leapfrog (drift-kick-drift), симплектик
        yp += -ksym[i] * y * dx                      # y'' = −κ·y  (κ=(q/p)G; фокусировка при κ>0)
        y += yp * dx / 2
        traj.append((xg[i + 1], y))
    return np.array(traj)


def main():
    os.makedirs(FIGS, exist_ok=True)
    fld = sys.argv[1] if len(sys.argv) > 1 else FLD
    name = os.path.splitext(os.path.basename(fld))[0]
    h, B = load(fld)
    x0, x1 = h["mnx"] + 1e-6, h["mxx"] - 1e-6            # диапазон из заголовка карты
    f_full = make_field(h, B, "full"); f_sym = make_field(h, B, "sym")
    fig, ax = plt.subplots(figsize=(10, 5.6))
    print(f"[{name}]  Финальное смещение Y и угол Y' (мм / мрад):")
    devs = []
    for y0, col in [(2e-3, "tab:red"), (-2e-3, "tab:blue")]:
        tf, _ = boris_track(y0, f_full, x0, x1); ts, _ = boris_track(y0, f_sym, x0, x1)
        tm = matrix_quad(y0, h, B, x0, x1)
        def fin(tr):
            yf = tr[-1, 1]; yp = (tr[-1, 1] - tr[-2, 1]) / (tr[-1, 0] - tr[-2, 0]); return yf, yp
        for lab, tr in [("Boris полное поле", tf), ("Boris симметризов.", ts)]:
            yf, yp = fin(tr); print(f"  y0={y0*1e3:+.0f}мм {lab:22}: Y={yf*1e3:+.4f}  Y'={yp*1e3:+.4f}")
        print(f"  y0={y0*1e3:+.0f}мм {'матрица (идеал симм.)':22}: Y={tm[-1,1]*1e3:+.4f}")
        devs.append(abs(tf[-1, 1] - tm[-1, 1]) / (abs(tm[-1, 1]) + 1e-30))
        ax.plot(tf[:, 0] * 1e3, tf[:, 1] * 1e3, "-", color=col, lw=1.8, label=f"Boris полное поле (y0={y0*1e3:+.0f}мм)")
        ax.plot(ts[:, 0] * 1e3, ts[:, 1] * 1e3, "--", color=col, lw=1.2, label=f"Boris симметризов. (y0={y0*1e3:+.0f}мм)")
        ax.plot(tm[:, 0] * 1e3, tm[:, 1] * 1e3, ":", color="k", lw=1.3, label=f"матрица (идеал, y0={y0*1e3:+.0f}мм)")
    bL = h["mxx"] * 1e3 * 0.5
    ax.axvspan(-bL, bL, color="0.9"); ax.axvline(0, color="0.6", ls=":", lw=1)
    ax.set_xlabel("X, мм"); ax.set_ylabel("Y, мм"); ax.grid(True, alpha=0.3); ax.legend(fontsize=7.5)
    ax.set_title(f"{name}: траектория Boris по реальному полю (сплошная) vs матричная модель "
                 f"(идеальный симм. квадруполь — точки)\nрасхождение Boris↔матрица ≈ {np.mean(devs)*100:.2f}%")
    fig.tight_layout(); p = os.path.join(FIGS, f"fig_asym_track_{name}.png"); fig.savefig(p, dpi=120); plt.close(fig)
    print(f"\nрасхождение Boris↔матрица ≈ {np.mean(devs)*100:.2f}%   фигура: {os.path.basename(p)}")


if __name__ == "__main__":
    main()
