# -*- coding: utf-8 -*-
"""ШАГ 3: разложение поля DQ вокруг ИЗОГНУТОЙ опорной орбиты (секторный подход).
Проблема шагов 1–2: осевые мультиполи неверны на большой сагитте (пучок уходит на ~40мм от оси).
Решение: (1) найти опорную дугу y_ref(X) (трек осевой частицы полным полем); (2) снять мультиполи
ВДОЛЬ ОРБИТЫ — b0(X),Gu(X),Gv(X),S(X) в точке (X, y_ref(X), 0); (3) реконструировать поле как
разложение по ОТКЛОНЕНИЮ от орбиты du=y−y_ref: Bz=b0+Gu·du+S(du²−z²), By=Gv·z+2S·du·z. Пучок
идёт вдоль орбиты ⇒ разложение точно там, где он есть. Трек тем же интегратором, сравнение с эталоном
(dq_traj.pkl) и с шагами 1–2. Оси: X прод., Y=u гор, Z=v верт.
Выход: workspace/dq_edge/dq_sector.png + RMS."""
import os, math, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dq_experiment as dqe
from _ct_asym_track import make_field

OUT = dqe.OUT


def orbit_multipoles(field_full, Xo, Yo, half, dy=1e-3, nX=500):
    """Мультиполи вдоль орбиты y_ref(X): b0,Gu,Gv,S на сетке X."""
    Xg = np.linspace(-half + 1e-4, half - 1e-4, nX)
    b0 = np.zeros(nX); Gu = np.zeros(nX); Gv = np.zeros(nX); S = np.zeros(nX)
    for i, X in enumerate(Xg):
        yc = float(np.interp(X, Xo, Yo))
        bz0 = field_full(np.array([X, yc, 0.0]))[2]
        bzp = field_full(np.array([X, yc + dy, 0.0]))[2]; bzm = field_full(np.array([X, yc - dy, 0.0]))[2]
        byp = field_full(np.array([X, yc, dy]))[1]; bym = field_full(np.array([X, yc, -dy]))[1]
        b0[i] = bz0; Gu[i] = (bzp - bzm) / (2 * dy); S[i] = 0.5 * (bzp - 2 * bz0 + bzm) / dy ** 2
        Gv[i] = (byp - bym) / (2 * dy)
    return Xg, b0, Gu, Gv, S


def sector_field(Xg, b0, Gu, Gv, S, Xo, Yo, half):
    """Поле-разложение вокруг орбиты: du=y−y_ref(X)."""
    def f(pos):
        X, y, z = pos
        if abs(X) > half:
            return np.zeros(3)
        yc = float(np.interp(X, Xo, Yo)); du = y - yc
        b = np.interp(X, Xg, b0); gu = np.interp(X, Xg, Gu); gv = np.interp(X, Xg, Gv); s = np.interp(X, Xg, S)
        Bz = b + gu * du + s * (du * du - z * z); By = gv * z + 2 * s * du * z
        return np.array([0.0, By, Bz])
    return f


def main():
    ref = pickle.load(open(os.path.join(OUT, "dq_traj.pkl"), "rb"))
    meta = ref["meta"]; XSTART = meta["XSTART"]; P0 = meta["P0"]; ds = meta["ds"]
    hs, Bs = zip(*[dqe.load_cached(m) for m in dqe.MAPS]); h = hs[0]
    half = h["mxx"] - h["mnx"]; nstep = int(round((meta["drift1"] + 2 * half + meta["drift2"]) / ds))
    print(f"Секторное разложение вокруг орбиты (диполь+квад+секступоль вдоль дуги), {nstep} шагов")
    print(f"{'настройка':>10} {'токи':>14} {'сагитта':>8} {'RMS,мм':>8} {'макс,мм':>8}")
    rows = []; viz = {}
    for lab, sett in ref["settings"].items():
        Bsup = sum((c / 100.0) * B for c, B in zip(sett["currents"], Bs))
        field_full = make_field(h, Bsup, "mirror")
        rt = sett["traj"][sett["groups"]["ref"]["gid0"]]        # опорная орбита [s,X,Y,Z]
        o = np.argsort(rt[:, 1]); Xo = rt[o, 1]; Yo = rt[o, 2]
        Xg, b0, Gu, Gv, S = orbit_multipoles(field_full, Xo, Yo, half)
        field = sector_field(Xg, b0, Gu, Gv, S, Xo, Yo, half)
        sagitta = np.abs(Yo[np.abs(Xo) < half]).max() * 1e3
        errs = []; trs_v = []
        for gname, g in sett["groups"].items():
            for i in range(len(g["init"])):
                u, up, v, vp = g["init"][i]
                dirv = np.array([1.0, up, vp]); dirv /= np.linalg.norm(dirv)
                tr, pf, mf = dqe.boris(np.array([XSTART, u, v]), P0 * dirv, field, nstep, ds)
                rp = g["final_pos"][i]; errs.append(math.hypot(pf[1] - rp[1], pf[2] - rp[2]))
                if lab in ("op", "s_I1") and gname == "size":
                    trs_v.append(tr)
        rms = np.sqrt(np.mean(np.square(errs))) * 1e3; mx = np.max(errs) * 1e3
        rows.append((lab, rms, mx))
        print(f"{lab:>10} {str(sett['currents']):>14} {sagitta:7.0f}м {rms:8.2f} {mx:8.2f}")
        if trs_v:
            viz[lab] = (sett, trs_v)
    print(f"\nСредний RMS (шаг 3, секторное разложение): {np.mean([r[1] for r in rows]):.2f}мм")
    print("(шаг 1 осевые мультиполи: 4.3мм; шаг 2 near-axis матрица связи: 6.5мм)")
    visualize(viz, half)
    print("Фигура:", os.path.join(OUT, "dq_sector.png"))


def visualize(viz, half):
    labs = list(viz.keys())
    fig, axs = plt.subplots(len(labs), 2, figsize=(13, 4.2 * len(labs))); axs = np.atleast_2d(axs)
    for row, lab in enumerate(labs):
        sett, trs = viz[lab]; g = sett["groups"]["size"]
        for col, comp in enumerate((2, 3)):
            ax = axs[row, col]; ax.axvspan(-half * 1e3, half * 1e3, color="0.9", zorder=0); ax.grid(alpha=0.3)
            for i in range(len(g["init"])):
                rt = sett["traj"][g["gid0"] + i]
                ax.plot(rt[:, 1] * 1e3, rt[:, comp] * 1e3, "-", lw=1.3, color="k", alpha=0.5)
                ax.plot(trs[i][:, 1] * 1e3, trs[i][:, comp] * 1e3, "--", lw=0.9, color="tab:purple")
            ax.set_ylabel(("Y гориз" if col == 0 else "Z верт") + ", мм")
            if row == 0: ax.set_title("плоскость изгиба X–Y" if col == 0 else "вертикаль X–Z")
            if row == len(labs) - 1: ax.set_xlabel("X, мм")
            if col == 0:
                ax.text(0.02, 0.95, f"{lab}\nI={sett['currents']}", transform=ax.transAxes, va="top",
                        fontsize=8, bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    fig.suptitle("DQ секторное разложение вокруг орбиты (фиолет. штрих) vs ЭТАЛОН (чёрный) — пучок 5мм",
                 fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "dq_sector.png"), dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
