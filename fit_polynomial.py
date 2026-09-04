# -*- coding: utf-8 -*-
"""Полиномиальная калибровочная кривая против закона насыщения — на одних данных.

Стандартная практика магнитных измерений: интегральную силу приближают НЕЧЁТНЫМ полиномом тока

    ∫G·dl(I) = a₁·I + a₃·I³ + a₅·I⁵ + …

Нечётность здесь не условность: поле магнита строго нечётно по току (проверено на паре карт
корректора при ±1808 А·в, |A+B| = 1e−12 Тл), поэтому чётные степени запрещены физически.

Наш закон насыщения того же назначения, но другой формы:

    NI = A·y + L·H(C·y)                A — воздух, L·H(C·y) — падение на железе

Сравнение ведётся ПО ОШИБКЕ ТОКА и перекрёстной проверкой с исключением по одной карте: подгонка
на всех точках ничего не говорит о том, что модель предскажет в промежутке.

ВАЖНОЕ РАЗЛИЧИЕ ФОРМ. Полином не умеет выходить на полку: при a₃ < 0 он рано или поздно
поворачивает ВНИЗ, и за пределами подгонки ответ теряет смысл. Закон насыщения асимптотику
воспроизводит по построению, потому что она заложена в кривую материала.

    py -3.14 fit_polynomial.py
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
from scipy.optimize import brentq, least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

KEY = {"QP_sh": "G", "QP_L": "G", "OP": "O", "SP": "S", "DQ": "G"}


def points(element, key, winding):
    d = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))
    out = []
    for r in d["records"]:
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not r["coeffs"] or key not in r["coeffs"]:
            continue
        if r.get("winding") != winding or "+" in (r.get("winding") or ""):
            continue
        out.append((ni, abs(r["coeffs"][key])))
    return sorted(out)


def poly_fit(ni, y, nterms):
    """Нечётный полином y = Σ a_k · I^(2k+1). Линеен по коэффициентам — точный МНК.

    Подгонка по ОТНОСИТЕЛЬНОЙ невязке: у слабых карт коэффициент на два порядка меньше, и при
    абсолютной метрике они в подгонке не участвовали бы вовсе."""
    M = np.vstack([ni ** (2 * k + 1) / y for k in range(nterms)]).T
    a, *_ = np.linalg.lstsq(M, np.ones_like(y), rcond=None)
    return a


def poly_eval(a, ni):
    return sum(a[k] * ni ** (2 * k + 1) for k in range(len(a)))


def poly_invert(a, target, hi):
    """Ток по заданной силе. Полином немонотонен при a₃<0 — берём ПЕРВЫЙ корень."""
    f = lambda x: poly_eval(a, x) - target
    lo = 1e-9
    if f(lo) * f(hi) > 0:
        # за максимумом полинома корня нет — ищем на восходящей ветви
        xs = np.linspace(lo, hi, 2000)
        v = poly_eval(a, xs) - target
        s = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
        if not len(s):
            return float("nan")
        lo, hi = xs[s[0]], xs[s[0] + 1]
    return float(brentq(f, lo, hi, xtol=1e-10, rtol=1e-13))


def law_fit(ni, y, Hof, seeds=30):
    def res(q):
        A, C, L = q[0], np.exp(q[1]), np.exp(q[2])
        return np.array([(A * v + L * float(np.atleast_1d(Hof(C * v))[0]) - n) / n
                         for n, v in zip(ni, y)])
    best = None
    sc = y.mean()
    for s in range(seeds):
        g = np.random.default_rng(s)
        q0 = [np.mean(ni / y) * g.lognormal(0, .3), np.log(1 / sc * g.lognormal(0, 1.2)),
              np.log(0.05 * g.lognormal(0, 1.2))]
        try:
            t = least_squares(res, q0, xtol=1e-14, ftol=1e-14, max_nfev=8000)
        except Exception:
            continue
        if best is None or t.cost < best.cost:
            best = t
    return best.x[0], np.exp(best.x[1]), np.exp(best.x[2])


def loo(pts, kind, Hof=None, nterms=2):
    """Ошибка ТОКА с исключением по одной карте."""
    err = []
    for i in range(len(pts)):
        tr = [q for j, q in enumerate(pts) if j != i]
        need = nterms if kind == "poly" else 3
        if len(tr) < need + 1:
            continue
        ni = np.array([q[0] for q in tr], float)
        y = np.array([q[1] for q in tr], float)
        n0, y0 = pts[i]
        if kind == "poly":
            a = poly_fit(ni, y, nterms)
            pred = poly_invert(a, y0, ni.max() * 3)
        else:
            A, C, L = law_fit(ni, y, Hof)
            pred = A * y0 + L * float(np.atleast_1d(Hof(C * y0))[0])
        if np.isfinite(pred):
            err.append(100 * abs(pred - n0) / n0)
    return float(np.mean(err)) if err else float("nan"), len(err)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--terms", type=int, nargs="*", default=[1, 2, 3, 4])
    ap.add_argument("--window", type=float, default=0.0,
                    help="сузить до рабочей области: во сколько раз ниже максимального тока "
                         "брать нижнюю границу (напр. 4 = от I_max/4 до I_max). "
                         "На практике калибровочную кривую строят именно так, а не на всём "
                         "диапазоне в 190 раз по току")
    a = ap.parse_args(argv)
    Hof = sm._material()
    print(f"  {'элемент':10} {'обм':4} {'карт':>5} "
          + "  ".join(f"полином {2 * t - 1}".rjust(12) for t in a.terms)
          + f" {'закон насыщ.':>14}")
    for el in ("QP_sh", "QP_L", "OP", "SP", "DQ"):
        fp = os.path.join(sm.OUTDIR, f"{el}.json")
        if not os.path.exists(fp):
            continue
        d = json.load(open(fp, encoding="utf-8"))
        key = KEY[el]
        for wd in sorted({r["winding"] for r in d["records"] if "+" not in (r["winding"] or "")}):
            pts = points(el, key, wd)
            if a.window > 1:
                top = max(q[0] for q in pts)
                pts = [q for q in pts if q[0] >= top / a.window]
            if len(pts) < 6:
                continue
            cells = [loo(pts, "poly", nterms=t)[0] for t in a.terms]
            law, _ = loo(pts, "law", Hof)
            print(f"  {el:10} {wd:4} {len(pts):5d} "
                  + "  ".join(f"{c:11.1f}%" for c in cells) + f" {law:13.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
