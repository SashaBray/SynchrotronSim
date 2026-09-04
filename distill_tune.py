# -*- coding: utf-8 -*-
"""Двухэтапная настройка неявной сети: плотный учитель — только для СТАРТА, доводка на РЕАЛЬНЫХ картах.

Зачем два этапа. У неявной сети внутри неподвижная точка, ландшафт задачи трудный, и на редких
реальных картах оптимизатор срывается (DQ при двух скрытых узлах: 9.0% вместо 2%). Плотные точки от
уже настроенного учителя ландшафт сглаживают и выводят в нужную область параметров.

Но обучаться НА учителе нельзя: между реальными точками он рисует кривую своей формы, и неявная
модель, повторяя её, съезжает с физически верных параметров. На QP_sh это стоило переносимости
смены материала (4.3% → 6.3%). Поэтому синтетика работает как стартовая точка и на этом её роль
кончается — доводка идёт по настоящим измерениям.

Учитель выбирается по элементу: у DQ лучше персептрон, у корректора — сетевая модель (МЭС).

    py -3.14 distill_tune.py corrector --teacher mec --hidden 2
    py -3.14 distill_tune.py DQ --teacher perceptron --hidden 2
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import implicit_magnet_net as im
import magnet_perceptron as mp

ORDER_BY_KEY = {"B": 1, "G": 2, "S": 3, "O": 4}


def rays(single, winds, n_per):
    """Сетка токов вдоль лучей одиночных обмоток — там, где учителя проверяли, то есть ИНТЕРПОЛЯЦИЯ.

    Комбинации не трогаем: они остаются независимой проверкой, а учитель на них и не обучался."""
    out = []
    for j, w in enumerate(winds):
        vals = [c.get(w, 0.0) for c, _, _, _ in single if w in c]
        if not vals:
            continue
        lo, hi = min(vals + [0.0]), max(vals + [0.0])
        for a in np.linspace(lo, hi, n_per):
            v = np.zeros(len(winds)); v[j] = a
            out.append(v)
    return out


def teacher_perceptron(el, winds, keys, K):
    p, w_, act, pts, dims = mp.fit(el, K, seeds=16, verbose=False)
    order = [w_.index(w) for w in winds]

    def f(I):
        return mp.predict(p, *dims, act, np.asarray(I)[order])
    return f, "персептрон"


def teacher_mec(el, winds, keys):
    """Сетевая модель: своя подгонка на каждую гармонику, знак берётся проекцией на опорную."""
    import branch_network as bn
    models = []
    for k in keys:
        cir, p, w = bn.fit(el, ORDER_BY_KEY[k], verbose=False)
        ref = bn.predict(cir, {winds[0]: 5000.0}, p, w)
        models.append((cir, p, w, ref))

    def f(I):
        cur = {w: float(v) for w, v in zip(winds, I) if v}
        out = []
        for cir, p, w, ref in models:
            if not cur:
                out.append(0.0); continue
            v = bn.predict(cir, cur, p, w)
            out.append(float(np.real(v * np.conj(ref))) / abs(ref))
        return np.array(out)
    return f, "МЭС (сетевая модель)"


def calib(teach, single, X, keys):
    """Учитель может жить в своих единицах (гармоника вместо физкоэффициента) — один масштаб на ключ."""
    k = []
    for i, key in enumerate(keys):
        num = den = 0.0
        for t in single:
            pr = teach(X[id(t)])[i]
            num += pr * t[1][key]; den += pr * pr
        k.append(num / den if den else 1.0)
    return np.array(k)


def make_resid(rows, Y, n_in, K, n_out, Hf, floor):
    def resid(p):
        out = []
        for I, y in zip(rows, Y):
            pr = im.predict(p, n_in, K, n_out, Hf, I)
            out.extend((pr - y) / (np.abs(y) + floor))
        return np.array(out)
    return resid


def run(el, teacher="perceptron", K=2, n_per=25, seeds=20, t_hidden=2):
    keys = mp.KEYS_BY_EL.get(el, ("B", "G"))
    pts = mp.load_points(el)
    winds = sorted({w for c, _, _, _ in pts for w in c})
    single = [t for t in pts if not t[2]]
    combo = [t for t in pts if t[2]]
    n_in, n_out = len(winds), len(keys)
    Hf = im.h_curve(im.sm.MATERIAL)
    X = {id(t): np.array([t[0].get(w, 0.0) for w in winds]) for t in pts}
    scale = {k: max(abs(m[k]) for _, m, _, _ in pts) for k in keys}
    floor = np.array([0.02 * scale[k] for k in keys])
    imax = max(np.max(np.abs(v)) for v in X.values()) or 1.0
    ymax = max(scale.values())

    Xr = [X[id(t)] for t in single]
    Yr = [np.array([t[1][k] for k in keys]) for t in single]
    res_real = make_resid(Xr, Yr, n_in, K, n_out, Hf, floor)

    teach, tname = (teacher_perceptron(el, winds, keys, t_hidden) if teacher == "perceptron"
                    else teacher_mec(el, winds, keys))
    kcal = calib(teach, single, X, keys)
    Xs = rays(single, winds, n_per)
    Ys = [kcal * teach(I) for I in Xs]
    res_syn = make_resid(Xs, Ys, n_in, K, n_out, Hf, floor)

    def starts():
        for seed in range(seeds):
            rng = np.random.default_rng(seed)
            yield np.concatenate([
                rng.lognormal(0, 0.5, K * n_in),
                np.log(rng.lognormal(np.log(imax / ymax), 1.0, K)),
                np.log(rng.lognormal(np.log(imax / ymax * 1e-3), 1.0, K)),
                np.log(rng.lognormal(np.log(20.0), 1.0, K)),
                rng.normal(0, 1.0, n_out * K),
                rng.normal(0, 1e-3, n_out * n_in) / imax])

    def best_of(res, p0s):
        b = None
        for p0 in p0s:
            try:
                s = least_squares(res, p0, xtol=1e-14, ftol=1e-14, max_nfev=20000)
            except Exception:
                continue
            if b is None or s.cost < b.cost:
                b = s
        return b

    def err(p, rows):
        e = []
        for t in rows:
            pr = im.predict(p, n_in, K, n_out, Hf, X[id(t)])
            for i, k in enumerate(keys):
                e.append(100 * abs(pr[i] - t[1][k]) / abs(t[1][k]))
        return (float(np.mean(e)), float(np.max(e))) if e else (float("nan"),) * 2

    print(f"=== {el}: неявная сеть, скрытых узлов {K}; учитель — {tname} ===")
    print(f"  обмоток {n_in} {winds}, карт {len(pts)} "
          f"(одиночных {len(single)}, комбинаций {len(combo)}), синтетики {len(Xs)}")
    rows = []
    a = best_of(res_real, list(starts()))
    rows.append(("только реальные карты", a))
    b = best_of(res_syn, list(starts()))
    rows.append(("обучение НА синтетике учителя", b))
    c = best_of(res_real, [b.x])
    rows.append(("синтетика ДЛЯ СТАРТА + доводка на реальных", c))
    print(f"  {'постановка':44} {'одиночные':>18} {'КОМБИНАЦИИ':>18}")
    for lbl, s in rows:
        es, xs = err(s.x, single); ec, xc = err(s.x, combo)
        print(f"  {lbl:44} {es:7.1f}% (макс{xs:6.1f}%) {ec:7.1f}% (макс{xc:6.1f}%)")
    return rows[-1][1].x


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="corrector")
    ap.add_argument("--teacher", choices=("perceptron", "mec"), default="perceptron")
    ap.add_argument("--hidden", type=int, default=2, help="скрытых узлов у неявной сети")
    ap.add_argument("--teacher-hidden", type=int, default=2)
    ap.add_argument("--points", type=int, default=25, help="синтетических точек на обмотку")
    a = ap.parse_args(argv)
    run(a.element, a.teacher, a.hidden, a.points, t_hidden=a.teacher_hidden)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
