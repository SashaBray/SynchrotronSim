# -*- coding: utf-8 -*-
"""Модель магнита как ПЕРСЕПТРОН с кривой намагничивания в роли функции активации.

    F = M1 · I                   — токи обмоток переводятся матрицей в аргументы кривой
    y = B(F)                     — сама кривая намагничивания, применяется поэлементно
    k = M2 · y + M3 · I + b      — кики: свёртка откликов ПЛЮС линейный обход и смещение

Вся нелинейность магнита сидит в одной функции — кривой стали, остальное линейно. Линейный обход
M3·I физически необходим: часть потока идёт по ВОЗДУХУ и не насыщается никогда, поэтому у отклика
всегда есть составляющая, пропорциональная току. Идея в том,
что поверхность момента над плоскостью токов выглядит как одна и та же кривая насыщения, у которой
вход и выход линейно растягиваются вдоль осей.

Наблюдаемая — коэффициенты, подогнанные по ТРАЕКТОРИЯМ (element_compare), то есть у изгибающих
элементов на симметричной дуге, а не гармоника у прямой оси: это разные величины.

    py -3.14 magnet_perceptron.py DQ 2          # элемент и число скрытых узлов
    py -3.14 magnet_perceptron.py corrector 3
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

# Какие коэффициенты описываем. У DQ кик двухкомпонентный (изгиб и градиент), у корректора
# осмыслен только сам кик B: G у него побочный и на малых значениях даёт бессмысленные проценты.
KEYS_BY_EL = {"DQ": ("B", "G"), "corrector": ("B",), "QP_sh": ("G",), "QP_L": ("G",),
              "SP": ("S",), "OP": ("O",)}
KEYS = ("B", "G")


def curve():
    """Кривая намагничивания как функция активации: B(H), нечётная, за таблицей — по воздуху."""
    bh = np.loadtxt(sm.MATERIAL, delimiter=";", skiprows=1)
    H, B = bh[:, 0], bh[:, 1]
    mu0 = 4e-7 * np.pi

    def f(x):
        a = np.abs(x)
        return np.sign(x) * np.where(a <= H[-1], np.interp(a, H, B),
                                     B[-1] + (a - H[-1]) * mu0)
    return f


def load_points(el):
    global KEYS
    KEYS = KEYS_BY_EL.get(el, KEYS)
    """Точки карт: (вектор токов по обмоткам в ампер-витках, {коэффициенты}, комбинация?, имя)."""
    d = sm.load_dataset(el)
    out = []
    for r in d["records"]:
        if not r["coeffs"] or any(k not in r["coeffs"] for k in KEYS):
            continue
        cur = {}
        for part in r["excitation"].replace("+", ";").split(";"):
            if "=" in part:
                k, v = part.split("=", 1)
                cur[k.strip()] = float(v)
        if cur:
            out.append((cur, {k: r["coeffs"][k] for k in KEYS}, len(cur) > 1, r["map_id"]))
    return out


def unpack(p, n_in, n_hid, n_out):
    i = 0
    M1 = p[i:i + n_hid * n_in].reshape(n_hid, n_in); i += n_hid * n_in
    M2 = p[i:i + n_out * n_hid].reshape(n_out, n_hid); i += n_out * n_hid
    M3 = p[i:i + n_out * n_in].reshape(n_out, n_in); i += n_out * n_in
    b = p[i:i + n_out]
    return M1, M2, M3, b


def predict(p, n_in, n_hid, n_out, act, I):
    M1, M2, M3, b = unpack(p, n_in, n_hid, n_out)
    return M2 @ act(M1 @ I) + M3 @ I + b


def fit(el, n_hid=2, seeds=12, verbose=True):
    pts = load_points(el)
    if not pts:
        raise SystemExit(f"нет данных по {el} — сначала saturation_model.py --element {el}")
    winds = sorted({w for c, _, _, _ in pts for w in c})
    act = curve()
    single = [t for t in pts if not t[2]]
    combo = [t for t in pts if t[2]]
    n_in, n_out = len(winds), len(KEYS)
    X = {id(t): np.array([t[0].get(w, 0.0) for w in winds]) for t in pts}
    scale = {k: max(abs(m[k]) for _, m, _, _ in pts) for k in KEYS}

    # Нормировка ОТНОСИТЕЛЬНАЯ, по каждой карте: при общем масштабе слабые карты в подгонке не
    # участвуют, и смещение на выходе становится для них больше самого сигнала (ошибки в сотни
    # процентов при формально малой невязке). Мягкий пол — чтобы не делить на почти ноль.
    def resid(p):
        out = []
        for t in single:
            pr = predict(p, n_in, n_hid, n_out, act, X[id(t)])
            for i, k in enumerate(KEYS):
                out.append((pr[i] - t[1][k]) / (abs(t[1][k]) + 0.02 * scale[k]))
        return np.array(out)

    # Масштаб первой матрицы: аргумент кривой должен попадать в её рабочую область, а токи
    # заданы в ампер-витках (тысячи) — без нормировки старт уходит в насыщенный хвост.
    imax = max(np.max(np.abs(v)) for v in X.values()) or 1.0
    best = None
    for seed in range(seeds):
        rng = np.random.default_rng(seed)
        p0 = np.concatenate([rng.normal(0, 1.0, n_hid * n_in) * (300.0 / imax),
                             rng.normal(0, 1e-3, n_out * n_hid),
                             rng.normal(0, 1e-3, n_out * n_in) / imax,
                             np.zeros(n_out)])
        try:
            s = least_squares(resid, p0, xtol=1e-14, ftol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if best is None or s.cost < best.cost:
            best = s
    p = best.x

    def err(rows):
        e = []
        for t in rows:
            pr = predict(p, n_in, n_hid, n_out, act, X[id(t)])
            for i, k in enumerate(KEYS):
                e.append(100 * abs(pr[i] - t[1][k]) / abs(t[1][k]))
        return (float(np.mean(e)), float(np.max(e))) if e else (float("nan"),) * 2

    if verbose:
        es, xs = err(single); ec, xc = err(combo)
        n_par = n_hid * n_in + n_out * n_hid + n_out * n_in + n_out
        print(f"=== {el}: персептрон, скрытых узлов {n_hid}, свободных чисел {n_par} ===")
        print(f"  обмоток {n_in} {winds}, карт {len(pts)} (одиночных {len(single)}, "
              f"комбинаций {len(combo)})")
        print(f"  подгонка на одиночных {es:.1f}% (макс {xs:.1f}%)   "
              f"ПРОВЕРКА на комбинациях {ec:.1f}% (макс {xc:.1f}%)")
        M1, M2, M3, b = unpack(p, n_in, n_hid, n_out)
        print("  M1 — токи в аргументы кривой:")
        for k in range(n_hid):
            print(f"     узел {k+1}: " + ", ".join(f"{w} {M1[k][j]:+.5g}"
                                                   for j, w in enumerate(winds)))
        print("  M2 — отклики кривой в кики, M3 — линейный обход, b — смещение:")
        for i, key in enumerate(KEYS):
            print(f"     {key}: " + ", ".join(f"узел{k+1} {M2[i][k]:+.5g}" for k in range(n_hid))
                  + " | обход " + ", ".join(f"{w} {M3[i][j]:+.4g}" for j, w in enumerate(winds))
                  + f" | смещение {b[i]:+.5g}")
    return p, winds, act, pts, (n_in, n_hid, n_out)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="DQ")
    ap.add_argument("hidden", nargs="?", type=int, default=2)
    a = ap.parse_args(argv)
    fit(a.element, a.hidden)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
