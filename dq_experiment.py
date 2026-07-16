# -*- coding: utf-8 -*-
"""ОПЫТ по диполь-квадруполю DQ (combined-function): эталонные траектории дифференциального ядра
SynchrotronSim по реальному 3D-полю для набора пучков и НАБОРА ТОКОВ трёх обмоток. У DQ три обмотки
(DQ_I1/I2/I3), поле = СУПЕРПОЗИЦИЯ: B(I1,I2,I3) = Σ (Iᵢ/100)·Bᵢ(100A) (постулат суперпозиции
подтверждён). Обмотки: I1 — диполь+квад (отриц. градиент), I2/I3 — диполь+квад (полож., I2≈I3),
плюс секступоль ⇒ три тока задают (изгиб B, градиент G, секступоль S).
Карты v3 (мелкая 2мм, полукарта X≥0, зеркалим по X — симметрия 0.1%). Компоновка:
0.5 м дрейф → элемент (~0.8 м с фринджем) → 1.5 м дрейф. Токи РЕАЛЬНЫЕ (при 100А изгиб ~205 мрад,
масштаб не нужен). Сканы: размер пучка, величина тока, КОМБИНАЦИИ (I1,I2,I3), смещения, углы влёта.
Оси: X — продольная, Y=u — горизонт (плоскость изгиба), Z=v — вертикаль.
Выход: workspace/dq_edge/dq_traj.pkl (+ figs)."""
import os, math, pickle
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from field_quality import load
from _ct_asym_track import make_field, ME, QE, V, GAMMA0, C

KD = r"C:/Sasha/Work/SILA/ansis approximations/KD/DQ"
MAPS = ["DQ_I1_100A_v3.fld", "DQ_I2_100A_v3.fld", "DQ_I3_100A_v3.fld"]
OUT = os.path.join("workspace", "research", "dq_edge")
os.makedirs(OUT, exist_ok=True)

P0 = GAMMA0 * ME * V; BRHO = P0 / abs(QE); QP = QE / P0
DRIFT1, DRIFT2 = 0.5, 0.5                            # м, пустые промежутки до/после
DS = 1e-3; SUB = 5                                  # шаг интегрирования, прореживание


def load_cached(name):
    """Загрузка карты с .npy-кэшем (парсинг 74МБ .fld — один раз)."""
    npy = os.path.join(OUT, name.replace(".fld", ".npy")); hpk = npy.replace(".npy", "_hdr.pkl")
    if os.path.exists(npy) and os.path.exists(hpk):
        return pickle.load(open(hpk, "rb")), np.load(npy)
    h, B = load(os.path.join(KD, name)); np.save(npy, B); pickle.dump(h, open(hpk, "wb"))
    return h, B


def boris(pos0, mom0, field, nstep, ds):
    """6D интегратор Бориса (движок дифф. ядра). Возврат (траектория [s,X,Y,Z], финал pos, mom)."""
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
    """Группы пучков (u,u',v,v') в СИ. δ=0 (монохром)."""
    g = {}
    size = []
    for r in (1e-3, 2e-3, 5e-3):
        for j in range(8):
            a = 2 * math.pi * j / 8
            size.append((r * math.cos(a), 0.0, r * math.sin(a), 0.0))
    g["size"] = size
    g["u_offset"] = [(uo, 0, 0, 0) for uo in (-5e-3, -3e-3, -1e-3, 1e-3, 3e-3, 5e-3)]
    g["v_offset"] = [(0, 0, vo, 0) for vo in (-5e-3, -3e-3, -1e-3, 1e-3, 3e-3, 5e-3)]
    g["u_angle"] = [(0, ap, 0, 0) for ap in (-10e-3, -5e-3, 5e-3, 10e-3)]
    g["v_angle"] = [(0, 0, 0, ap) for ap in (-10e-3, -5e-3, 5e-3, 10e-3)]
    g["ref"] = [(0, 0, 0, 0)]
    return g


def main():
    print("Загрузка трёх обмоток DQ (v3, кэш)…")
    hs, Bs = zip(*[load_cached(m) for m in MAPS]); h = hs[0]
    half = h["mxx"] - h["mnx"]; Xentry = -half; XSTART = Xentry - DRIFT1
    Lfull = 2.0 * half; total = DRIFT1 + Lfull + DRIFT2; nstep = int(round(total / DS))
    print(f"Магнит X∈[{Xentry*1e3:.0f},{-Xentry*1e3:.0f}]мм (зеркало), путь {total:.1f}м, {nstep} шагов ×{DS*1e3:.0f}мм")

    def field_for(I1, I2, I3):
        Bsup = (I1 / 100.0) * Bs[0] + (I2 / 100.0) * Bs[1] + (I3 / 100.0) * Bs[2]
        return make_field(h, Bsup, "mirror")

    def init(u, up, v, vp):
        dirv = np.array([1.0, up, vp]); dirv /= np.linalg.norm(dirv)
        return np.array([XSTART, u, v]), P0 * dirv

    # НАБОР ТОКОВ (I1,I2,I3): рабочее СООТНОШЕНИЕ I1:I2≈27:44 масштабировано ВНИЗ до апертуры
    # (при полном рабочем ±(27,44)А квад DQ за кроссовером, пучок вылетает из ±50мм). Одиночные
    # обмотки — малые (поодиночке квад сильнее: рабочая комбинация I1/I2 частично гасит друг друга).
    CURRENTS = [
        ("op", -14, -22, 0), ("op_half", -7, -11, 0),                 # рабочее соотношение (в апертуре)
        ("I1_a", -18, -22, 0), ("I1_b", -8, -22, 0),                  # варьируем I1 (I2 фикс.)
        ("I2_a", -14, -14, 0), ("I2_b", -14, -28, 0),                 # варьируем I2 (I1 фикс.)
        ("I3_p", -14, -22, 10), ("I3_m", -14, -22, -10),              # добавляем I3
        ("s_I1", -12, 0, 0), ("s_I1p", 12, 0, 0),                     # изолируем обмотку (малый ток)
        ("s_I2", 0, -15, 0), ("s_I3", 0, 0, -15),
    ]
    groups = beam_groups()
    data = {"meta": {"maps": MAPS, "drift1": DRIFT1, "drift2": DRIFT2, "Lfull": Lfull,
                     "Xentry": Xentry, "XSTART": XSTART, "ds": DS, "P0": P0, "BRHO": BRHO,
                     "QP": QP, "GAMMA0": GAMMA0, "aperture_mm": 50}, "settings": {}}
    for (lab, I1, I2, I3) in CURRENTS:
        field = field_for(I1, I2, I3)
        # эталонный изгиб/градиент этой комбинации (для справки)
        xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]; j0, k0 = h["ny"] // 2, h["nz"] // 2
        Bsup = (I1 / 100.0) * Bs[0][:, j0, k0, 2] + (I2 / 100.0) * Bs[1][:, j0, k0, 2] + (I3 / 100.0) * Bs[2][:, j0, k0, 2]
        bend = QP * 2.0 * np.trapezoid(Bsup, xg)
        sett = {"currents": (I1, I2, I3), "bend_mrad": bend * 1e3, "groups": {}, "traj": {}}
        gid = 0
        for gname, gl in groups.items():
            inits, fpos, fdir, fmom = [], [], [], []
            for (u, up, v, vp) in gl:
                pos0, mom0 = init(u, up, v, vp)
                tr, pf, mf = boris(pos0, mom0, field, nstep, DS)
                inits.append([u, up, v, vp]); fpos.append(pf)
                fdir.append(mf / np.linalg.norm(mf)); fmom.append(np.linalg.norm(mf))
                sett["traj"][gid] = tr; gid += 1
            sett["groups"][gname] = {"init": np.array(inits), "final_pos": np.array(fpos),
                                     "final_dir": np.array(fdir), "gid0": gid - len(gl)}
        rd = sett["groups"]["ref"]["final_dir"][0]
        sett["ref_bend_mrad"] = math.atan2(rd[1], rd[0]) * 1e3
        # макс поперечная экскурсия ВНУТРИ магнита (|X|<half) — контроль апертуры ±50мм
        exc = 0.0
        for g in sett["traj"].values():
            ins = np.abs(g[:, 1]) < half
            if ins.any():
                exc = max(exc, np.abs(g[ins][:, 2:4]).max())
        sett["max_exc_mm"] = exc * 1e3
        data["settings"][lab] = sett
        flag = "ВЫЛЕТ!" if exc * 1e3 > 50 else "ок"
        print(f"  [{lab:12}] I=({I1:+4d},{I2:+4d},{I3:+4d})А  изгиб опорной {sett['ref_bend_mrad']:+7.1f}мрад  "
              f"экскурсия в магните {exc*1e3:4.0f}мм [{flag}]  ({gid} частиц)")

    with open(os.path.join(OUT, "dq_traj.pkl"), "wb") as f:
        pickle.dump(data, f)
    visualize(data)
    print("Сохранено:", os.path.join(OUT, "dq_traj.pkl"))


def visualize(data):
    m = data["meta"]; Xe = m["Xentry"]
    show = ["op", "s_I1", "s_I2", "I3_p"]
    fig, axs = plt.subplots(len(show), 2, figsize=(13, 3.0 * len(show)))
    for row, lab in enumerate(show):
        s = data["settings"][lab]
        for col, (comp, ylab, ttl) in enumerate([(2, "Y гориз (изгиб), мм", "плоскость изгиба X–Y"),
                                                  (3, "Z верт, мм", "вертикаль X–Z")]):
            ax = axs[row, col]; ax.axvspan(Xe * 1e3, -Xe * 1e3, color="0.9", zorder=0); ax.grid(alpha=0.3)
            for gn, cmap in [("size", plt.cm.Blues), ("u_offset" if col == 0 else "v_offset", plt.cm.Reds)]:
                g = s["groups"][gn]; n = len(g["init"])
                for i in range(n):
                    tr = s["traj"][g["gid0"] + i]
                    ax.plot(tr[:, 1] * 1e3, tr[:, comp] * 1e3, lw=0.7, color=cmap(0.4 + 0.5 * i / n))
            tr0 = s["traj"][s["groups"]["ref"]["gid0"]]
            ax.plot(tr0[:, 1] * 1e3, tr0[:, comp] * 1e3, "k-", lw=1.2)
            ax.set_ylabel(ylab)
            if row == 0: ax.set_title(ttl)
            if row == len(show) - 1: ax.set_xlabel("X, мм")
            if col == 0:
                ax.text(0.02, 0.95, f"{lab}\nI={s['currents']}А\nизгиб {s['ref_bend_mrad']:+.0f}мрад",
                        transform=ax.transAxes, va="top", fontsize=7.5,
                        bbox=dict(boxstyle="round", fc="white", alpha=0.7))
    fig.suptitle("DQ (диполь-квадруполь): эталонные траектории дифф. ядра SynchrotronSim по 3D-полю "
                 "(суперпозиция трёх обмоток). Размер пучка — синие, смещения — красные, опорная — чёрная",
                 fontsize=10.5)
    fig.tight_layout(); p = os.path.join(OUT, "dq_traj.png"); fig.savefig(p, dpi=130); plt.close(fig)
    print("Фигура:", p)


if __name__ == "__main__":
    main()
