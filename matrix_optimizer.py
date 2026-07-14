# -*- coding: utf-8 -*-
"""Оптимизатор матричной модели SynchrotronSim — аналог диффурного (Boris) оптимизатора, но
настраивает ФИЗКОЭФФИЦИЕНТЫ матричных элементов методом последовательных приближений.

Цель (как в Boris): ВЫЖИВАНИЕ ПУЧКА на полном кольце. Трекаем пучок через 10 супериодов (реш. 107,
520 магнитов), минимизируем максимальную бетатронную амплитуду (откл. от опорной орбиты) + штраф за
потерянные частицы. Ручки — ОБЩИЕ на супериод физкоэф. (G квадруполей, S секступолей, O октуполей,
G/S у DQ); повторяются на все 10 SP (принцип общих токов). Изгиб B диполей/DQ ФИКСИРОВАН геометрией.

Старт — проектные K из .lte (lte_warmstart): физически устойчивая точка. Метод — координатный спуск
с line-search по каждой ручке (последовательное приближение). Трекинг векторизован по частицам (numpy).

Зависимости: matrix_tracker (геометрия кольца, рамки, модифиц. дрейф), lte_warmstart (проектные K).
"""
import os, csv, sys, time, shutil
import numpy as np
import matrix_tracker as mt
import lte_warmstart as ws

QP = mt.QP
PERIOD = 52          # сгруппированных магнитов в стандартном супериоде
COEFF_IDX = {"G": 0, "S": 1, "O": 2}


# ==================== кольцо и ручки ====================
def build_ring(lattice=107):
    eff, phys, coup, trans = mt.load_db()
    return mt.build_line(lattice, eff, phys, coup, trans, source="physical")


def build_knobspec():
    """Warm-start карта (52 группы) → список ручек [(group_idx, coeff, init_value)] и стартовый вектор.
    Квадруполь→G, секступоль→S, октуполь→O, DQ→(G,S). Диполи/корректоры без ручек."""
    mapping = ws.build_mapping(105, verbose=False)
    knobs = []
    for r in mapping:
        if r["kind"] == "quad":
            knobs.append((r["idx"], "G", r["G"]))
        elif r["kind"] == "sext":
            knobs.append((r["idx"], "S", r["S"]))
        elif r["kind"] == "oct":
            knobs.append((r["idx"], "O", r["O"]))
        elif r["kind"] == "dq":
            knobs.append((r["idx"], "G", r["G"]))
            knobs.append((r["idx"], "S", r["S"]))
    x0 = np.array([k[2] for k in knobs], float)
    return mapping, knobs, x0


def apply_knobs(mags, mapping, knobs, xv):
    """Проставить k1/k2/k3 всем 520 магнитам из вектора ручек xv (общий на супериод). Изгиб theta0
    диполей И DQ — ПРОЕКТНЫЙ (ANGLE из .lte, геометрия), а не из устаревших токов БД (у DQ mat_B в БД
    завышен впятеро от старого тюна)."""
    gso = {r["idx"]: [r["G"], r["S"], r["O"]] for r in mapping}   # база = проект (для диполей нули)
    ang = {r["idx"]: r["angle"] for r in mapping}
    kind = {r["idx"]: r["kind"] for r in mapping}
    for (gi, c, _), v in zip(knobs, xv):
        gso[gi][COEFF_IDX[c]] = v
    for k, m in enumerate(mags):
        gi = k % PERIOD
        G, S, O = gso[gi]
        m.k1 = QP * G; m.k2 = QP * S; m.k3 = QP * O
        if kind[gi] in ("dipole", "dq") and ang[gi]:                 # изгиб — проектный (θ = ANGLE)
            m.theta0 = ang[gi]
        m.nonzero = any(abs(v) > 0 for v in (m.theta0, m.k1, m.k2, m.k3))


# ==================== векторный трекинг пучка ====================
def _drift_vec(S, alive, ell):
    d, px, py = S["d"], S["px"], S["py"]
    arg = (1.0 + d) ** 2 - px * px - py * py
    newly = alive & (arg <= 1e-12)
    alive &= ~newly
    ps = np.sqrt(np.where(arg > 1e-12, arg, 1.0))
    S["x"] = np.where(alive, S["x"] + ell * px / ps, S["x"])
    S["y"] = np.where(alive, S["y"] + ell * py / ps, S["y"])
    S["t"] = np.where(alive, S["t"] + ell * (1.0 + d) / ps, S["t"])


def _body_vec(m, S, alive, nsl):
    if not m.nonzero:
        _drift_vec(S, alive, m.L); return
    N = nsl; dsg = m.L / N; dsf = m.Leff / N
    th, k1, k2, k3 = m.theta0, m.k1, m.k2, m.k3
    for _ in range(N):
        _drift_vec(S, alive, dsg / 2.0)
        x, y = S["x"], S["y"]
        S["px"] = np.where(alive, S["px"] + th / N - k1 * x * dsf
                           - 0.5 * k2 * (x * x - y * y) * dsf
                           - (k3 / 6.0) * (x ** 3 - 3 * x * y * y) * dsf, S["px"])
        S["py"] = np.where(alive, S["py"] + k1 * y * dsf + k2 * x * y * dsf
                           + (k3 / 6.0) * (3 * x * x * y - y ** 3) * dsf, S["py"])
        _drift_vec(S, alive, dsg / 2.0)


def _mdrift_vec(S, alive, A, B):
    """Геометрический пролёт плоскость выхода A → плоскость входа B, векторно по частицам."""
    x, y, px, py, d = S["x"], S["y"], S["px"], S["py"], S["d"]
    P = A.exit[None, :] + x[:, None] * A.uhat[None, :] + y[:, None] * A.vhat[None, :]   # (Np,3)
    arg = (1.0 + d) ** 2 - px * px - py * py
    newly = alive & (arg <= 1e-12); alive &= ~newly
    ps = np.sqrt(np.where(arg > 1e-12, arg, 1.0))
    dir_ = (px[:, None] * A.uhat + py[:, None] * A.vhat + ps[:, None] * A.axis) / (1.0 + d)[:, None]
    denom = dir_ @ B.axis
    newly = alive & (denom <= 1e-9); alive &= ~newly
    denom_s = np.where(np.abs(denom) > 1e-12, denom, 1.0)
    tt = ((B.entry[None, :] - P) @ B.axis) / denom_s
    newly = alive & (tt < 0); alive &= ~newly
    Q = P + tt[:, None] * dir_
    rel = Q - B.entry[None, :]
    S["x"] = np.where(alive, rel @ B.uhat, S["x"])
    S["y"] = np.where(alive, rel @ B.vhat, S["y"])
    S["px"] = np.where(alive, (1.0 + d) * (dir_ @ B.uhat), S["px"])
    S["py"] = np.where(alive, (1.0 + d) * (dir_ @ B.vhat), S["py"])
    S["t"] = np.where(alive, S["t"] + tt, S["t"])


def track_beam(mags, P0, aper=0.03, nsl=5):
    """Пучок P0 (Np,6) через кольцо. Потеря = |смещение от оси бокса|>aper или p⊥>p. Смещение —
    АБСОЛЮТНОЕ в локальной рамке (не относит. опорной): устойчиво к нестабильной опорной орбите.
    Возврат: (maxoff_per_particle, reached_fraction, lost_at_magnet, alive_mask)."""
    Np = P0.shape[0]
    S = {k: P0[:, i].astype(float).copy() for i, k in enumerate(("x", "px", "y", "py", "t", "d"))}
    alive = np.ones(Np, bool)
    maxoff = np.zeros(Np)
    lost_at = np.full(Np, -1)
    total = len(mags)

    def stamp(j):
        nl = (~alive) & (lost_at < 0); lost_at[nl] = j

    _body_vec(mags[0], S, alive, nsl); stamp(0)
    for j in range(1, total):
        _mdrift_vec(S, alive, mags[j - 1], mags[j]); stamp(j)
        off = np.sqrt(S["x"] * S["x"] + S["y"] * S["y"])
        maxoff = np.maximum(maxoff, np.where(alive, off, maxoff))
        hit = alive & (off > aper)
        lost_at[hit & (lost_at < 0)] = j
        alive &= ~hit
        _body_vec(mags[j], S, alive, nsl); stamp(j)
    reached = np.where(lost_at < 0, 1.0, lost_at / total)
    maxoff = np.where(alive, maxoff, aper)                      # потерянные — на апертуре
    return maxoff, reached, lost_at, alive


def loss_of(mags, P0, aper=0.03, nsl=5):
    """Скаляр loss (м) = среднее по частицам: выжившие → макс. смещение; потерянные → aper·(2−доля
    пройденного кольца) (раньше потеряна ⇒ больше). Двигает к ВЫЖИВАНИЮ (горизонт) и малому пучку."""
    maxoff, reached, lost_at, alive = track_beam(mags, P0, aper, nsl)
    cost = np.where(alive, maxoff, aper * (2.0 - reached))
    loss = float(np.mean(cost))
    horizon = int(np.min(np.where(lost_at < 0, len(mags), lost_at)))
    info = dict(nlost=int(np.sum(~alive)), horizon=horizon, lost_at=lost_at, alive=alive,
                maxoff_surv=float(np.max(maxoff[alive])) if alive.any() else 0.0)
    return loss, info


# ==================== пучок ====================
def make_beam(radii=(1e-3,), deltas=(0.0,), with_dispersion=True):
    """Рабочий пучок ±1мм: опорная (на оси) + полное кольцо в (x,y) радиуса 1мм + дисперсионные (δ)
    частицы. Частицы с амплитудой >2мм считаем неустойчивыми — оптимизируем только рабочую область."""
    P = [(0, 0, 0, 0, 0, 0)]
    for dlt in deltas:
        for r in radii:
            P += [(r, 0, 0, 0, 0, dlt), (-r, 0, 0, 0, 0, dlt), (0, 0, r, 0, 0, dlt), (0, 0, -r, 0, 0, dlt),
                  (r, 0, r, 0, 0, dlt), (-r, 0, r, 0, 0, dlt), (r, 0, -r, 0, 0, dlt), (-r, 0, -r, 0, 0, dlt)]
    if with_dispersion:
        P += [(0, 0, 0, 0, 0, +0.005), (0, 0, 0, 0, 0, -0.005)]
    return np.array(P, float)


# ==================== оптимизатор (координатный спуск / line-search) ====================
def optimize(P0, aper=0.02, nsl=5, rounds=12, lattice=107, verbose=True):
    mapping, knobs, x0 = build_knobspec()
    mags = build_ring(lattice)

    def f(xv):
        apply_knobs(mags, mapping, knobs, xv)
        return loss_of(mags, P0, aper, nsl)

    x = x0.copy()
    best, info = f(x)
    step = np.abs(x0) * 0.25 + 3e-4                      # относит.+абс. шаг на ручку
    hist = [best]
    total = len(mags)
    if verbose:
        print(f"Кольцо {lattice}: {total} магнитов, {len(knobs)} ручек, пучок {P0.shape[0]} частиц, "
              f"апертура {aper*1e3:.0f}мм, срезов {nsl}.")
        print(f"Старт (проект .lte): loss={best*1e3:.2f}мм  горизонт={info['horizon']}/{total}  "
              f"потеряно={info['nlost']}  макс.смещ={info['maxoff_surv']*1e3:.1f}мм")
    t0 = time.time()
    for rnd in range(rounds):
        improved = 0
        for i in range(len(x)):
            for s in (step[i], -step[i]):
                xt = x.copy(); xt[i] += s
                lt, it = f(xt)
                if lt < best - 1e-9:
                    best, info, x = lt, it, xt; improved += 1
                    break
        if verbose:
            print(f"  раунд {rnd + 1:2}: loss={best*1e3:7.2f}мм  горизонт={info['horizon']:3}/{total}  "
                  f"потеряно={info['nlost']:2}  макс.смещ={info['maxoff_surv']*1e3:5.1f}мм  "
                  f"улучшено ручек={improved:2}  [{time.time()-t0:.0f}с]")
        hist.append(best)
        if improved == 0:
            step *= 0.5
            if np.all(step < 1e-6):
                break
    if verbose:
        print(f"Финал: loss={best*1e3:.2f}мм  горизонт={info['horizon']}/{total}  "
              f"потеряно={info['nlost']}  за {time.time()-t0:.0f}с")
    return x, best, info, mapping, knobs, hist


# ==================== сохранение настроенных коэффициентов ====================
def final_coeffs(x, mapping, knobs):
    """Итоговые физкоэф. на 52 группы супериода (проект + перекрытия ручек)."""
    gso = {r["idx"]: [r["G"], r["S"], r["O"]] for r in mapping}
    for (gi, c, _), v in zip(knobs, x):
        gso[gi][COEFF_IDX[c]] = v
    return [dict(idx=r["idx"], devs=r["devs"], lte=r["lte"], kind=r["kind"], angle=r["angle"],
                 G=gso[r["idx"]][0], S=gso[r["idx"]][1], O=gso[r["idx"]][2]) for r in mapping]


def save_tune_csv(rows, path=None):
    """Настроенный супериод → dedicated таблица (авторитетный результат тюнинга)."""
    path = path or os.path.join(mt.TB, "matrix_tune.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["group_idx", "devices", "lte_name", "kind", "G_Tlm", "S_Tlm2", "O_Tlm3", "bend_angle_rad"])
        for r in rows:
            w.writerow([r["idx"], r["devs"], r["lte"] or "", r["kind"],
                        f"{r['G']:.9g}", f"{r['S']:.9g}", f"{r['O']:.9g}", f"{r['angle']:.9g}"])
    print(f"  сохранено: {os.path.relpath(path)} ({len(rows)} групп)")


def write_lattice_configs(rows, eff, lattices=(105, 107), period=PERIOD, backup=True):
    """Записать mat_B/G/S/O в lattice_configs для period-52 решёток (реплика супериода на все SP).
    Многообмоточные: коэф. в ПЕРВУЮ обмотку, остальные 0 (сумма трекера = коэф). mat_B = θ/(QP·Leff)."""
    coeff = {r["idx"]: r for r in rows}
    path = os.path.join(mt.TB, "lattice_configs.csv")
    allrows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(allrows[0].keys())
    if backup:
        bdir = os.path.join(mt.TB, "backups"); os.makedirs(bdir, exist_ok=True)
        bpath = os.path.join(bdir, "lattice_configs_pre_matrixtune.csv")
        shutil.copy(path, bpath); print(f"  бэкап: {os.path.relpath(bpath)}")
    for lat in (str(l) for l in lattices):
        rws = sorted([r for r in allrows if r["lattice_id"] == lat], key=lambda r: int(r["instance_id"]))
        i = 0; gpos = 0
        while i < len(rws):
            r0 = rws[i]; grp = [r0]
            while i + 1 < len(rws) and (rws[i + 1]["x"], rws[i + 1]["y"], rws[i + 1]["z"]) == (r0["x"], r0["y"], r0["z"]):
                grp.append(rws[i + 1]); i += 1
            c = coeff[gpos % period]; Leff = eff.get(grp[0]["device_type_id"], 1.0)
            matB = c["angle"] / (QP * Leff) if (c["kind"] in ("dipole", "dq") and c["angle"]) else 0.0
            for k, gr in enumerate(grp):
                gr["mat_B"], gr["mat_G"], gr["mat_S"], gr["mat_O"] = (
                    (f"{matB:.9g}", f"{c['G']:.9g}", f"{c['S']:.9g}", f"{c['O']:.9g}") if k == 0
                    else ("0", "0", "0", "0"))
            gpos += 1; i += 1
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(allrows)
    print(f"  обновлены mat_* в lattice_configs: решётки {list(lattices)}")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "warmstart"
    aper = float(sys.argv[2]) if len(sys.argv) > 2 else 0.02
    P0 = make_beam()
    if mode == "warmstart":
        mapping, knobs, x0 = build_knobspec()
        mags = build_ring(107)
        apply_knobs(mags, mapping, knobs, x0)
        t0 = time.time()
        L, info = loss_of(mags, P0, aper)
        total = len(mags)
        print(f"WARM-START (проектные K из .lte), кольцо 107, {total} магнитов, {P0.shape[0]} частиц:")
        print(f"  loss={L*1e3:.2f}мм  горизонт={info['horizon']}/{total}  потеряно={info['nlost']}/"
              f"{P0.shape[0]}  макс.смещ={info['maxoff_surv']*1e3:.1f}мм  [{time.time()-t0:.1f}с/трек]")
        if info["nlost"]:
            la = info["lost_at"]
            print(f"  потеряны на магнитах: {sorted(int(v) for v in la if v>=0)[:20]}")
    elif mode == "optimize":
        x, best, info, mapping, knobs, hist = optimize(P0, aper=aper)
        print("Сохранение настроенных коэффициентов:")
        rows = final_coeffs(x, mapping, knobs)
        save_tune_csv(rows)
        write_lattice_configs(rows, mt.load_db()[0], lattices=(105, 107))
    elif mode == "verify":
        # round-trip: трекер читает НАСТРОЕННЫЕ mat_* из БД (без warm-start/apply_knobs)
        mags = build_ring(107); total = len(mags)
        Pbig = make_beam(radii=(1e-3, 3e-3, 5e-3), with_dispersion=True)
        print(f"ВЕРИФИКАЦИЯ настроенного кольца из БД (source=physical), {total} магнитов, {Pbig.shape[0]} частиц "
              f"(кольца 1/3/5мм + δ=±0.5%):")
        for nsl in (5, 40):
            L, info = loss_of(mags, Pbig, aper, nsl)
            print(f"  срезов={nsl:2}: loss={L*1e3:6.2f}мм  горизонт={info['horizon']}/{total}  "
                  f"потеряно={info['nlost']}/{Pbig.shape[0]}  макс.смещ={info['maxoff_surv']*1e3:.1f}мм")
