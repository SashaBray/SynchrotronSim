# -*- coding: utf-8 -*-
"""Точность модели В ТОКЕ: на сколько процентов ошибётся уставка источника.

Проценты по коэффициенту кика обманчивы: в насыщении кривая пологая, и промах в кике на 10% может
означать промах в токе втрое. Практический вопрос обратный — какой ток выставлять, чтобы получить
нужный кик. Поэтому метрика здесь такая: для каждой карты ищется множитель λ, при котором модель на
токе λ·I даёт ИЗМЕРЕННЫЙ кик; |λ − 1| и есть относительная ошибка тока.

Проверяются две формы связи «потоки → кик»:
    база          кик = Σ w_класс · Φ_класс                     (нынешняя, веса постоянны)
    +насыщение    кик = Σ (w_класс + σ·v_класс) · Φ_класс       (веса зависят от состояния железа)

σ — потоко-взвешенная доля насыщения по железным ветвям, от 0 (линейный режим) до 1 (все ветви за
коленом). Смысл поправки: с ростом возбуждения меняется не только величина, но и ФОРМА поля —
у OP это намерено до 17%, у DQ даёт перекос диполь/квадруполь на 12%. Веса остаются линейными,
поэтому оба варианта по-прежнему решаются точным МНК, а не поиском.

    py -3.14 current_accuracy.py corrector
    py -3.14 current_accuracy.py DQ --order 2
"""
from __future__ import annotations

import argparse
import os
import sys
from collections import Counter

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import branch_network as bn
import mec_sigma as ms
import saturation_model as sm


# Базис и веса берутся ИЗ mec_sigma — единственного источника модели. Своя копия здесь уже
# один раз разошлась с продакшном: разделение кик-классов по столбцу полюса попало в mec_sigma,
# а тут осталась старая версия, и сводка молча показывала прежние цифры.
def state(cir, currents, p):
    return ms.state(cir, currents, p)


def basis(cir, currents, p, augment):
    return ms.basis(cir, currents, augment, p)


def fit_weights(cir, rows, p, augment, wts):
    A = np.array([basis(cir, c, p, augment) for c, _, _, _ in rows])
    y = np.array([m for _, m, _, _ in rows], dtype=complex)
    den = (np.abs(y) + 0.05 * float(np.max(np.abs(y)))) \
        / np.array([wts[tuple(sorted(c))] for c, _, _, _ in rows])
    M = np.vstack([np.real(A) / den[:, None], np.imag(A) / den[:, None]])
    v = np.concatenate([np.real(y) / den, np.imag(y) / den])
    lam = 1e-3 * float(np.sqrt((M ** 2).sum() / M.shape[0]))
    M = np.vstack([M, lam * np.eye(M.shape[1])])
    v = np.concatenate([v, np.zeros(M.shape[1])])
    return np.linalg.lstsq(M, v, rcond=None)[0]


def current_error(cir, cur, meas, p, augment, w):
    """|λ − 1| в процентах: во сколько раз надо изменить ТОК, чтобы модель дала измеренный кик."""
    u = meas / abs(meas)

    def f(lam):
        pr = complex(basis(cir, {k: v * lam for k, v in cur.items()}, p, augment) @ w)
        return float(np.real(pr * np.conj(u))) - abs(meas)

    lo, hi = 1e-3, 50.0
    try:
        if f(lo) * f(hi) > 0:
            return float("nan")
        return 100.0 * abs(brentq(f, lo, hi, xtol=1e-10, rtol=1e-12) - 1.0)
    except (ValueError, RuntimeError):
        return float("nan")


def run(element, order=1):
    cir = bn.Circuit(element)
    if cir.geom is None:
        raise SystemExit(f"{element}: геометрия задана не полностью")
    maps = bn.read_maps(element, order)
    single = [t for t in maps if len(t[0]) == 1]
    combo = [t for t in maps if len(t[0]) > 1]
    cnt = Counter(tuple(sorted(t[0])) for t in single)
    wts = {k: 1.0 / np.sqrt(v) for k, v in cnt.items()}
    p = np.zeros(1)
    bn.calibrate_senses(cir, single, p)

    print(f"=== {element}, гармоника {order}: ошибка В ТОКЕ (|λ−1|, %) ===")
    print(f"  карт {len(maps)} (одиночных {len(single)}, комбинаций {len(combo)})")
    res = {}
    for augment in (False, True):
        w = fit_weights(cir, single, p, augment, wts)
        res[augment] = w
        rows = []
        for lbl, grp in (("одиночные", single), ("КОМБИНАЦИИ", combo)):
            e = [current_error(cir, c, m, p, augment, w) for c, m, _, _ in grp]
            e = [x for x in e if np.isfinite(x)]
            rows.append((lbl, np.mean(e), np.max(e), len(e), len(grp)))
        name = "+насыщение" if augment else "база      "
        print(f"  {name}  свободных чисел {len(w):2d}: " + "   ".join(
            f"{lbl} {m:5.1f}% (макс {x:5.1f}%, n={n}/{t})" for lbl, m, x, n, t in rows))

    w = res[True]
    # Разбивка ПО ОБМОТКАМ: ошибка в токе резко различается у главной и у боковых, и общее среднее
    # это скрывает. Комбинации сведены в отдельную строку.
    groups = {}
    for cur, meas, _, mid in maps:
        key = "КОМБИНАЦИИ" if len(cur) > 1 else next(iter(cur))
        groups.setdefault(key, []).append((cur, meas))
    print()
    print(f"  {'обмотка':12} {'карт':>5} {'база: сред':>11} {'макс':>7} {'+насыщ: сред':>13} {'макс':>7} {'усиление':>9}")
    for key in sorted(groups, key=lambda k: (k == "КОМБИНАЦИИ", k)):
        e0 = [current_error(cir, c, m, p, False, res[False]) for c, m in groups[key]]
        e1 = [current_error(cir, c, m, p, True, res[True]) for c, m in groups[key]]
        e0 = [x for x in e0 if np.isfinite(x)]; e1 = [x for x in e1 if np.isfinite(x)]
        # Обусловленность на ВЕРХНЕМ конце тока: d ln|ток| / d ln|кик| по самим картам. Если
        # отклик пологий (усиление 5x), то ошибка в кике на 1% даёт 5% ошибки тока у ЛЮБОЙ модели —
        # это свойство задачи, а не промах подгонки.
        amp = float("nan")
        if len(groups[key]) > 1:
            xs = np.log([max(abs(v) for v in c.values()) for c, _ in groups[key]])
            ys = np.log([abs(m) for _, m in groups[key]])
            o = np.argsort(xs)
            xs, ys = xs[o], ys[o]
            # наклон между ДВУМЯ ВЕРХНИМИ картами: глобальная подгонка усредняет линейный участок
            # с насыщенным и занижает усиление там, где оно и важно
            dx, dy = xs[-1] - xs[-2], ys[-1] - ys[-2]
            amp = abs(dx / dy) if abs(dy) > 1e-9 else float("inf")
        print(f"  {key:12} {len(groups[key]):5d} {np.mean(e0):10.1f}% {np.max(e0):6.1f}% "
              f"{np.mean(e1):12.1f}% {np.max(e1):6.1f}% {amp:9.1f}x")


    print(f"\n  {'карта':32} {'ток, А·в':>10} {'база':>8} {'+насыщение':>11}")
    for cur, meas, _, mid in sorted(maps, key=lambda t: t[3]):
        amp = max(abs(v) for v in cur.values())
        e0 = current_error(cir, cur, meas, p, False, res[False])
        e1 = current_error(cir, cur, meas, p, True, w)
        mark = "" if not np.isfinite(e0) or not np.isfinite(e1) else (" <" if e1 < e0 - 0.5 else "")
        print(f"  {mid:32} {amp:10.0f} {e0:7.1f}% {e1:10.1f}%{mark}")
    return res


def single_branch(element):
    """Элементы БЕЗ графа (одно вырожденное ребро): закон NI = A·y + L·H(C·y).

    Обратная задача здесь ЯВНАЯ — по измеренному коэффициенту ампер-витки считаются напрямую,
    без поиска корня. Поэтому ошибка тока получается сразу и не зависит от обусловленности поиска."""
    import json
    d = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))
    Hof = sm._material()
    out = {}
    for wd, f in d["fits"].items():
        key, A, C, L = f["coeff"], f["A"], f["C"], f["L"]
        e = []
        for r in d["records"]:
            if r["winding"] != wd or not r["coeffs"] or key not in r["coeffs"]:
                continue
            ni = abs(float(r["amp_turns"]))
            if ni <= 0:
                continue
            y = abs(r["coeffs"][key])
            pred = A * y + L * float(np.atleast_1d(Hof(C * y))[0])
            e.append(100.0 * abs(pred - ni) / ni)
        if e:
            out[wd] = (float(np.mean(e)), float(np.max(e)), len(e), key)
    return out


def library():
    """Точность В ТОКЕ по всей библиотеке: МЭС-σ там, где есть граф, однобранчевый закон иначе."""
    import mec_sigma as ms
    print("=== точность В ТОКЕ по библиотеке (|λ−1|, %) ===")
    print(f"  {'элемент':10} {'модель':12} {'обмотка':10} {'карт':>5} {'среднее':>9} {'макс':>8}")
    M = ms.load()
    for (el, order), (cir, w) in sorted(M.items()):
        maps = bn.read_maps(el, order)
        groups = {}
        for cur, meas, _, mid in maps:
            groups.setdefault("КОМБИНАЦИИ" if len(cur) > 1 else next(iter(cur)), []).append((cur, meas))
        for k in sorted(groups, key=lambda k: (k == "КОМБИНАЦИИ", k)):
            e = [100 * abs(ms.current_scale(el, order, c, m, M) - 1.0) for c, m in groups[k]]
            e = [x for x in e if np.isfinite(x)]
            print(f"  {el:10} {'МЭС-σ n=' + str(order):12} {k:10} {len(e):5d} "
                  f"{np.mean(e):8.1f}% {np.max(e):7.1f}%")
    for el in ("QP_sh", "QP_L", "SP", "OP"):
        try:
            res = single_branch(el)
        except (OSError, KeyError, ValueError) as exc:
            print(f"  {el:10} — нет данных ({exc})")
            continue
        for wd, (m, x, n, key) in sorted(res.items()):
            print(f"  {el:10} {'одна ветвь':12} {wd + ' (' + key + ')':10} {n:5d} "
                  f"{m:8.1f}% {x:7.1f}%")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="corrector")
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--library", action="store_true", help="сводка по всей библиотеке")
    a = ap.parse_args(argv)
    if a.library:
        library()
        return 0
    run(a.element, a.order)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
