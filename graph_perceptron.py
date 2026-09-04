# -*- coding: utf-8 -*-
"""Граф-персептрон: матрица графа берётся из чертежей, свободна только свёртка потоков в кик.

Идея (предложена пользователем): вместо свободных матриц персептрона взять матрицу, которую задаёт
САМ ГРАФ магнита, — токи в потоки, потоки через кривую намагничивания, потом в кики.

При буквальной реализации схема вырождается: линейное решение графа даёт потоки, ЛИНЕЙНЫЕ по току, и
кривая, применённая к потокам зазоров, никакого насыщения не создаст. Нелинейность появляется лишь
когда пересчитанная кривой проводимость железа возвращается ОБРАТНО в решение графа. Поэтому здесь
схема развёрнута в N шагов Пикара:

    g⁰ = μ_r·μ₀·S/ℓ                              линейная проводимость железа (одно число μ_r)
    повторить N раз:  Φ = решение графа(g)       линейная задача, матрица из геометрии
                      H = напор / ℓ              напряжённость в железном ребре
                      g = S·B(H) / (ℓ·H)         секущая проводимость по кривой материала
    кик = Σ w_класс · (Σ Φ_ребра · направление)  направления — из координат наконечников

N = 1 — это в точности персептрон (один проход вперёд), N → ∞ — сошедшаяся магнитная цепь (МЭС).
Ручка N меряет, насколько для данного магнита важно ПЕРЕРАСПРЕДЕЛЕНИЕ потока при насыщении.

Свободных чисел мало: μ_r плюс по одному весу на класс воздушных ветвей (корректор 1+6, DQ 1+2), и
веса находятся точным МНК, а не поиском.

    py -3.14 graph_perceptron.py corrector
    py -3.14 graph_perceptron.py DQ --order 2
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import branch_network as bn

STEPS = (1, 2, 3, 5, 10, 40)


def linear_solve(cir, F, g):
    """Один линейный проход по графу: A·Φ = 0 при Φ = g·(напор). Матрица — из геометрии."""
    ua = cir.A[:, :]
    G = ua * g                                        # узел × ребро, взвешенное проводимостью
    M = G @ cir.A.T
    rhs = G @ F                                       # A·diag(g)·(F − Aᵀu) = 0
    try:
        u = np.linalg.solve(M + 1e-18 * np.eye(M.shape[0]), rhs)
    except np.linalg.LinAlgError:
        u = np.linalg.lstsq(M, rhs, rcond=None)[0]
    d = np.where(cir.ia >= 0, u[np.maximum(cir.ia, 0)], 0.0) \
        - np.where(cir.ic >= 0, u[np.maximum(cir.ic, 0)], 0.0) + F
    return g * d, d


def fluxes(cir, currents, mu_r, steps):
    """Потоки ветвей после N шагов Пикара. steps=1 — линейный проход (персептрон)."""
    cir._prep()
    lens, perms, secs = cir.geom
    Fd = cir.mmf(currents)
    F = np.array([Fd[b] for b in cir.branches])
    L = np.array([lens[cir.sym[b]] if cir.kind[b] != "gap" else 1.0 for b in cir.branches])
    S = np.array([secs.get(cir.sym[b], 1.0) if cir.kind[b] != "gap" else 1.0
                  for b in cir.branches])
    P = np.array([perms[cir.sym[b]] if cir.kind[b] == "gap" else 0.0 for b in cir.branches])
    g = np.where(cir.is_gap, P, mu_r * cir.mu0 * S / L)
    phi = None
    for _ in range(steps):
        phi, d = linear_solve(cir, F, g)
        H = np.where(cir.is_gap, 0.0, d / L)
        # секущая проводимость: g = S·B(H)/(ℓ·H); у нуля берём начальную (иначе деление на ноль)
        with np.errstate(divide="ignore", invalid="ignore"):
            gs = S * cir.B_of_H(H) / (L * H)
        g = np.where(cir.is_gap, P, np.where(np.abs(H) > 1e-12, gs, mu_r * cir.mu0 * S / L))
    return {b: float(v) for b, v in zip(cir.branches, phi)}


def basis(cir, currents, mu_r, steps):
    """Базис «поток → кик»: на класс воздушных ветвей одно число, направление из геометрии."""
    phi = fluxes(cir, currents, mu_r, steps)
    return np.array([sum(phi[b] * cir.direction(b) for b in cir.gaps if cir.sym[b] == c)
                     for c in cir.kick_classes], dtype=complex)


def weights(cir, rows, mu_r, steps, wts):
    A = np.array([basis(cir, cur, mu_r, steps) for cur, _, _, _ in rows])
    y = np.array([m for _, m, _, _ in rows], dtype=complex)
    den = (np.abs(y) + 0.05 * float(np.max(np.abs(y)))) \
        / np.array([wts[tuple(sorted(c))] for c, _, _, _ in rows])
    M = np.vstack([np.real(A) / den[:, None], np.imag(A) / den[:, None]])
    b = np.concatenate([np.real(y) / den, np.imag(y) / den])
    lam = 1e-3 * float(np.sqrt((M ** 2).sum() / M.shape[0]))
    M = np.vstack([M, lam * np.eye(M.shape[1])])
    b = np.concatenate([b, np.zeros(M.shape[1])])
    return np.linalg.lstsq(M, b, rcond=None)[0]


def run(element, order=1, steps=STEPS, verbose=True):
    from collections import Counter
    cir = bn.Circuit(element)
    if cir.geom is None:
        raise SystemExit(f"{element}: геометрия задана не у всех ветвей — матрицу графа не построить")
    maps = bn.read_maps(element, order)
    single = [t for t in maps if len(t[0]) == 1]
    combo = [t for t in maps if len(t[0]) > 1]
    cnt = Counter(tuple(sorted(t[0])) for t in single)
    wts = {k: 1.0 / np.sqrt(v) for k, v in cnt.items()}
    scale = max(abs(t[1]) for t in single)
    bn.calibrate_senses(cir, single, np.zeros(1))       # знаки намотки — из данных, как в МЭС

    if verbose:
        print(f"=== {element}, гармоника {order}: граф-персептрон "
              f"(матрица графа из чертежей, свободных чисел {1 + len(cir.kick_classes)}) ===")
        print(f"  карт {len(maps)} (одиночных {len(single)}, комбинаций {len(combo)}); "
              f"классов воздуха для кика {len(cir.kick_classes)}")
        print(f"  {'шагов Пикара':>14} {'одиночные':>20} {'КОМБИНАЦИИ':>20}")
    out = {}
    for n in steps:
        def resid(p):
            mu = float(np.exp(p[0]))
            w = weights(cir, single, mu, n, wts)
            r = []
            for cur, meas, _, _ in single:
                k = wts[tuple(sorted(cur))]
                d = k * (complex(basis(cir, cur, mu, n) @ w) - meas) / (abs(meas) + 0.05 * scale)
                r += [float(np.real(d)), float(np.imag(d))]
            return np.array(r)

        best = None
        for mu0 in (1e2, 1e3, 5e3):
            try:
                s = least_squares(resid, [np.log(mu0)], xtol=1e-12, ftol=1e-12, max_nfev=200)
            except Exception:
                continue
            if best is None or s.cost < best.cost:
                best = s
        mu = float(np.exp(best.x[0]))
        w = weights(cir, single, mu, n, wts)

        def err(rows):
            e = [100 * abs(complex(basis(cir, c, mu, n) @ w) - m) / abs(m) for c, m, _, _ in rows]
            return (float(np.mean(e)), float(np.max(e))) if e else (float("nan"),) * 2
        es, xs = err(single); ec, xc = err(combo)
        out[n] = (mu, es, ec)
        if verbose:
            tag = "  (= персептрон)" if n == 1 else ("  (= МЭС)" if n == max(steps) else "")
            print(f"  {n:14d} {es:9.1f}% (макс{xs:6.1f}%) {ec:9.1f}% (макс{xc:6.1f}%)"
                  f"   μr={mu:.0f}{tag}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="corrector")
    ap.add_argument("--order", type=int, default=1, help="порядок гармоники: 1 диполь, 2 квадруполь")
    ap.add_argument("--steps", type=int, nargs="*", default=list(STEPS))
    a = ap.parse_args(argv)
    run(a.element, a.order, tuple(a.steps))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
