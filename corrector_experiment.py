# -*- coding: utf-8 -*-
"""ОПЫТ по КОРРЕКТОРУ (5 обмоток): эталонные траектории дифференциального ядра SynchrotronSim по
реальному 3D-полю для широкого скана комбинаций пяти токов. Обмотки: main (гор. диполь-корректор,
∫Bz≠0) + correct_1..4 (квадрупольные, градиенты ±3.3мТл, слабые верт. компоненты). Поле —
СУПЕРПОЗИЦИЯ: B(I) = Σ (Iⱼ/100)·Bⱼ(100A). Карты v3 (2мм, полукарта X≥0, зеркалим по X). Компоновка:
0.5м дрейф → корректор (~0.4м с зеркалом) → 1.0м дрейф. Токи ±100А (поле корректора мало — пучок в
апертуре). Сканы токов: опорная (0), одиночные обмотки (±50/±100), main+correct пары, случайные
5-мерные комбинации (seeded). Каждой комбинации — пучок разного размера, смещений, углов.
Оси: X прод., Y=u гор, Z=v верт. Выход: workspace/corrector_edge/corr_traj.pkl (+ figs)."""
import os, math, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
from _ct_asym_track import make_field, ME, QE, V, GAMMA0, C

KD = r"C:/Sasha/Work/SILA/ansis approximations/KD/corrector"
WINDINGS = ["main", "correct_1", "correct_2", "correct_3", "correct_4"]
MAPS = [f"corrector_{w}_100A_v3.fld" for w in WINDINGS]
OUT = os.path.join("workspace", "corrector_edge")
os.makedirs(OUT, exist_ok=True)

P0 = GAMMA0 * ME * V; BRHO = P0 / abs(QE); QP = QE / P0
DRIFT1, DRIFT2 = 0.5, 1.0; DS = 1e-3; SUB = 5


def load_cached(name):
    npy = os.path.join(OUT, name.replace(".fld", ".npy")); hpk = npy.replace(".npy", "_hdr.pkl")
    if os.path.exists(npy) and os.path.exists(hpk):
        return pickle.load(open(hpk, "rb")), np.load(npy)
    h, B = load(os.path.join(KD, name)); np.save(npy, B); pickle.dump(h, open(hpk, "wb"))
    return h, B


def boris(pos0, mom0, field, nstep, ds):
    pos = pos0.astype(float).copy(); mom = mom0.astype(float).copy()
    g0 = math.sqrt(1.0 + mom.dot(mom) / (ME * C) ** 2)
    dt = ds / (math.sqrt(mom.dot(mom)) / (g0 * ME)); traj = [(0.0, pos[0], pos[1], pos[2])]; s = 0.0
    for i in range(nstep):
        g = math.sqrt(1.0 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = field(pos); gm = math.sqrt(1.0 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; sfac = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, sfac)
        g = math.sqrt(1.0 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME)); s += ds
        if (i + 1) % SUB == 0:
            traj.append((s, pos[0], pos[1], pos[2]))
    return np.array(traj), pos.copy(), mom.copy()


def beam_groups():
    g = {}
    size = []
    for r in (1e-3, 3e-3):
        for j in range(8):
            a = 2 * math.pi * j / 8
            size.append((r * math.cos(a), 0.0, r * math.sin(a), 0.0))
    g["size"] = size
    g["u_offset"] = [(uo, 0, 0, 0) for uo in (-5e-3, -2e-3, 2e-3, 5e-3)]
    g["v_offset"] = [(0, 0, vo, 0) for vo in (-5e-3, -2e-3, 2e-3, 5e-3)]
    g["u_angle"] = [(0, ap, 0, 0) for ap in (-5e-3, 5e-3)]
    g["v_angle"] = [(0, 0, 0, ap) for ap in (-5e-3, 5e-3)]
    g["ref"] = [(0, 0, 0, 0)]
    return g


def gen_currents():
    """Широкий скан 5-мерного пространства токов (I_main, I1, I2, I3, I4)."""
    combos = [("ref", [0, 0, 0, 0, 0])]
    for i, w in enumerate(WINDINGS):                       # одиночные обмотки
        for I in (50, 100, -50, -100):
            c = [0] * 5; c[i] = I; combos.append((f"{w}_{I:+d}", c))
    for i in range(1, 5):                                  # main + каждый correct
        c = [70, 0, 0, 0, 0]; c[i] = 70; combos.append((f"main+c{i}", c))
    combos.append(("H_kick", [80, 0, 0, 0, 0]))            # структурные
    combos.append(("quad_norm", [0, 60, -60, 60, -60]))    # попытка нормального квада
    combos.append(("quad_skew", [0, 60, 60, -60, -60]))    # попытка скью
    rng = np.random.default_rng(42)                        # случайные (seeded, воспроизводимо)
    for k in range(28):
        combos.append((f"rand{k:02d}", [int(x) for x in rng.integers(-100, 101, 5)]))
    return combos


def main():
    print("Загрузка 5 обмоток корректора (v3, кэш)…")
    hs, Bs = zip(*[load_cached(m) for m in MAPS]); h = hs[0]
    # проверка X-симметрии (валидность зеркала)
    j0, k0 = h["ny"] // 2, h["nz"] // 2; nx = h["nx"]
    sym = np.sqrt(np.mean((Bs[0][:, j0 + 2, k0, 2] - Bs[0][::-1, j0 + 2, k0, 2]) ** 2)) if h["mnx"] < -1e-9 else 0.0
    half = h["mxx"] - h["mnx"]; Xentry = -half; XSTART = Xentry - DRIFT1
    Lfull = 2.0 * half; total = DRIFT1 + Lfull + DRIFT2; nstep = int(round(total / DS))
    print(f"Магнит X∈[{Xentry*1e3:.0f},{-Xentry*1e3:.0f}]мм (зеркало), путь {total:.1f}м, {nstep} шагов")

    def field_for(cur):
        Bsup = sum((cur[j] / 100.0) * Bs[j] for j in range(5))
        return make_field(h, Bsup, "mirror")

    def init(u, up, v, vp):
        dirv = np.array([1.0, up, vp]); dirv /= np.linalg.norm(dirv)
        return np.array([XSTART, u, v]), P0 * dirv

    groups = beam_groups(); currents = gen_currents()
    xg = h["mnx"] + np.arange(nx) * h["sx"]
    data = {"meta": {"windings": WINDINGS, "maps": MAPS, "drift1": DRIFT1, "drift2": DRIFT2,
                     "Lfull": Lfull, "Xentry": Xentry, "XSTART": XSTART, "ds": DS, "P0": P0,
                     "BRHO": BRHO, "QP": QP, "aperture_mm": 50}, "settings": {}}
    print(f"Комбинаций токов: {len(currents)}, частиц/комб: {sum(len(v) for v in groups.values())}")
    for lab, cur in currents:
        field = field_for(cur)
        # справочный кик по полю на оси (гор Bz, верт By)
        Bz = sum((cur[j] / 100.0) * Bs[j][:, j0, k0, 2] for j in range(5))
        By = sum((cur[j] / 100.0) * Bs[j][:, j0, k0, 1] for j in range(5))
        hk = QP * 2.0 * np.trapezoid(Bz, xg); vk = QP * 2.0 * np.trapezoid(By, xg)
        sett = {"currents": cur, "hkick_mrad": hk * 1e3, "vkick_mrad": vk * 1e3, "groups": {}, "traj": {}}
        gid = 0; exc = 0.0
        for gname, gl in groups.items():
            inits, fpos, fdir = [], [], []
            for (u, up, v, vp) in gl:
                pos0, mom0 = init(u, up, v, vp)
                tr, pf, mf = boris(pos0, mom0, field, nstep, DS)
                inits.append([u, up, v, vp]); fpos.append(pf); fdir.append(mf / np.linalg.norm(mf))
                sett["traj"][gid] = tr; gid += 1
                ins = np.abs(tr[:, 1]) < half
                if ins.any(): exc = max(exc, np.abs(tr[ins][:, 2:4]).max())
            sett["groups"][gname] = {"init": np.array(inits), "final_pos": np.array(fpos),
                                     "final_dir": np.array(fdir), "gid0": gid - len(gl)}
        rd = sett["groups"]["ref"]["final_dir"][0]
        sett["ref_hkick_mrad"] = math.atan2(rd[1], rd[0]) * 1e3
        sett["ref_vkick_mrad"] = math.atan2(rd[2], rd[0]) * 1e3
        sett["max_exc_mm"] = exc * 1e3
        data["settings"][lab] = sett
    print(f"X-симметрия (валидность зеркала): RMS {sym*1e3:.3f}мТл (≈0 ⇒ ок)")
    print(f"{'настройка':>10} {'токи (main,c1..c4)':>22} {'H-кик':>8} {'V-кик':>8} {'экскур':>7}")
    for lab, sett in list(data["settings"].items())[:14]:
        print(f"{lab:>10} {str(sett['currents']):>22} {sett['ref_hkick_mrad']:+7.1f}м {sett['ref_vkick_mrad']:+7.1f}м {sett['max_exc_mm']:6.0f}мм")
    print(f"  … всего {len(data['settings'])} настроек")
    with open(os.path.join(OUT, "corr_traj.pkl"), "wb") as f:
        pickle.dump(data, f)
    visualize(data, half)
    print("Сохранено:", os.path.join(OUT, "corr_traj.pkl"))


def visualize(data, half):
    show = ["main_+100", "correct_1_+100", "quad_norm", "rand00", "rand05", "main+c1"]
    show = [s for s in show if s in data["settings"]][:6]
    fig, axs = plt.subplots(len(show), 2, figsize=(13, 2.6 * len(show))); axs = np.atleast_2d(axs)
    for row, lab in enumerate(show):
        s = data["settings"][lab]
        for col, comp in enumerate((2, 3)):
            ax = axs[row, col]; ax.axvspan(-half * 1e3, half * 1e3, color="0.9", zorder=0); ax.grid(alpha=0.3)
            for gn, cmap in [("size", plt.cm.Blues), ("u_offset" if col == 0 else "v_offset", plt.cm.Reds)]:
                g = s["groups"][gn]; n = len(g["init"])
                for i in range(n):
                    tr = s["traj"][g["gid0"] + i]
                    ax.plot(tr[:, 1] * 1e3, tr[:, comp] * 1e3, lw=0.7, color=cmap(0.4 + 0.5 * i / n))
            tr0 = s["traj"][s["groups"]["ref"]["gid0"]]
            ax.plot(tr0[:, 1] * 1e3, tr0[:, comp] * 1e3, "k-", lw=1.2)
            ax.set_ylabel(("Y гориз" if col == 0 else "Z верт") + ", мм")
            if row == 0: ax.set_title("плоскость X–Y (гор.)" if col == 0 else "плоскость X–Z (верт.)")
            if row == len(show) - 1: ax.set_xlabel("X, мм")
            if col == 0:
                ax.text(0.02, 0.95, f"{lab}\nI={s['currents']}\nH{s['ref_hkick_mrad']:+.0f} V{s['ref_vkick_mrad']:+.0f}мрад",
                        transform=ax.transAxes, va="top", fontsize=6.5, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    fig.suptitle("Корректор (5 обмоток): эталонные траектории дифф. ядра SynchrotronSim. Размер — синие, "
                 "смещения — красные, опорная — чёрная", fontsize=10.5)
    fig.tight_layout(); p = os.path.join(OUT, "corr_traj.png"); fig.savefig(p, dpi=130); plt.close(fig)
    print("Фигура:", p)


if __name__ == "__main__":
    main()
