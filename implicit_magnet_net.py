# -*- coding: utf-8 -*-
"""Неявная сеть магнита: персептрон с ОБРАТНОЙ СВЯЗЬЮ, кривая намагничивания внутри уравнения.

Отличие от `magnet_perceptron.py` — одно, но решающее. Там кривая применяется к линейной функции
ТОКА за один проход вперёд:

    h = σ(W₁·T·I)                       рабочая точка кривой фиксируется при обучении

здесь скрытый слой ищется как НЕПОДВИЖНАЯ ТОЧКА, а кривая применяется к ИНДУКЦИИ в железе:

    W₁⊙h + W₂⊙H(V⊙h) = T·I              баланс ампер-витков против падения на железе
    k = G·h + M₃·I                      кики: свёртка откликов плюс линейный обход (воздух)

Это ровно структура уравнения магнитной цепи (узловой метод), записанная как обучаемая сеть:
матрица инцидентности и размещение обмоток заменены обучаемыми матрицами, а условие баланса в узлах
превращается в поиск неподвижной точки. Если развернуть итерации — получится глубокая сеть с общими
весами и обходными связями.

ЗАЧЕМ ЭТО НУЖНО: только такая форма переносит СМЕНУ МАТЕРИАЛА. Проверено на QP_sh (две серии карт,
сталь 1008 против кастомной кривой с насыщением 1.198 Тл вместо 2.085):

    персептрон     K=1: подгонка 3.0%, перенос 15.8%   K=2: подгонка 0.2%, перенос 25.6%
    неявная сеть   K=1: подгонка 0.6%, перенос  4.3%   K=2: подгонка 0.0%, перенос  4.3%

У персептрона улучшение подгонки ПОРТИТ перенос; у неявной сети ширина слоя безопасна.

    py -3.14 implicit_magnet_net.py DQ 2
    py -3.14 implicit_magnet_net.py QP_sh 1 --swap workspace/materials/lowsat_custom_bh.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm
from magnet_perceptron import KEYS_BY_EL, load_points

MAXIT, TOL = 200, 1e-14


def h_curve(path):
    """Обратная кривая H(B): аргумент — ИНДУКЦИЯ, как и должно быть в уравнении цепи."""
    bh = np.loadtxt(path, delimiter=";", skiprows=1)
    B, H = bh[:, 1], bh[:, 0]
    mu0 = 4e-7 * np.pi

    def f(b):
        a = np.abs(b)
        return np.sign(b) * np.where(a <= B[-1], np.interp(a, B, H), H[-1] + (a - B[-1]) / mu0)
    return f


def unpack(p, n_in, K, n_out):
    i = 0
    T = p[i:i + K * n_in].reshape(K, n_in); i += K * n_in
    W1 = np.exp(p[i:i + K]); i += K
    W2 = np.exp(p[i:i + K]); i += K
    V = np.exp(p[i:i + K]); i += K
    G = p[i:i + n_out * K].reshape(n_out, K); i += n_out * K
    M3 = p[i:i + n_out * n_in].reshape(n_out, n_in)
    return T, W1, W2, V, G, M3


def solve_h(T, W1, W2, V, Hf, I):
    """Неподвижная точка: W1⊙h + W2⊙H(V⊙h) = T·I. Ньютон с численной производной."""
    F = T @ I
    h = np.zeros_like(F)
    for _ in range(MAXIT):
        r = W1 * h + W2 * Hf(V * h) - F
        eps = np.maximum(np.abs(h), 1e-9) * 1e-6
        d = W1 + W2 * (Hf(V * (h + eps)) - Hf(V * h)) / eps
        hn = h - r / np.maximum(d, 1e-12)
        if np.max(np.abs(hn - h)) < TOL * max(1.0, np.max(np.abs(hn))):
            return hn
        h = hn
    return h


def predict(p, n_in, K, n_out, Hf, I):
    T, W1, W2, V, G, M3 = unpack(p, n_in, K, n_out)
    return G @ solve_h(T, W1, W2, V, Hf, I) + M3 @ I


def fit(el, K=1, material=None, seeds=20, verbose=True):
    keys = KEYS_BY_EL.get(el, ("B", "G"))
    pts = load_points(el)
    winds = sorted({w for c, _, _, _ in pts for w in c})
    n_in, n_out = len(winds), len(keys)
    Hf = h_curve(material or sm.MATERIAL)
    single = [t for t in pts if not t[2]]
    combo = [t for t in pts if t[2]]
    X = {id(t): np.array([t[0].get(w, 0.0) for w in winds]) for t in pts}
    scale = {k: max(abs(m[k]) for _, m, _, _ in pts) for k in keys}
    imax = max(np.max(np.abs(v)) for v in X.values()) or 1.0
    ymax = max(scale.values())

    def resid(p):
        out = []
        for t in single:
            pr = predict(p, n_in, K, n_out, Hf, X[id(t)])
            for i, k in enumerate(keys):
                out.append((pr[i] - t[1][k]) / (abs(t[1][k]) + 0.02 * scale[k]))
        return np.array(out)

    best = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        p0 = np.concatenate([
            rng.lognormal(0, 0.5, K * n_in),
            np.log(rng.lognormal(np.log(imax / ymax), 1.0, K)),
            np.log(rng.lognormal(np.log(imax / ymax * 1e-3), 1.0, K)),
            np.log(rng.lognormal(np.log(20.0), 1.0, K)),
            rng.normal(0, 1.0, n_out * K),
            rng.normal(0, 1e-3, n_out * n_in) / imax])
        try:
            s = least_squares(resid, p0, xtol=1e-14, ftol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    p = best.x

    def err(rows, hf=None):
        e = []
        for t in rows:
            pr = predict(p, n_in, K, n_out, hf or Hf, X[id(t)])
            for i, k in enumerate(keys):
                e.append(100 * abs(pr[i] - t[1][k]) / abs(t[1][k]))
        return float(np.mean(e)) if e else float("nan")

    if verbose:
        n_par = K * n_in + 3 * K + n_out * K + n_out * n_in
        print(f"=== {el}: неявная сеть, скрытых узлов {K}, свободных чисел {n_par} ===")
        print(f"  подгонка на одиночных {err(single):.1f}%   "
              f"ПРОВЕРКА на комбинациях {err(combo):.1f}%")
    return p, (n_in, K, n_out), Hf, pts, X, keys, err


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="DQ")
    ap.add_argument("hidden", nargs="?", type=int, default=1)
    ap.add_argument("--swap", help="кривая другого материала: предсказать её БЕЗ переподгонки")
    a = ap.parse_args(argv)
    p, dims, Hf, pts, X, keys, err = fit(a.element, a.hidden)
    if a.swap:
        single = [t for t in pts if not t[2]]
        print(f"  ПЕРЕНОС на {os.path.basename(a.swap)} (ноль новых параметров): "
              f"{err(single, h_curve(a.swap)):.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
