# -*- coding: utf-8 -*-
"""МЭС-σ — рабочая модель связи «токи обмоток ↔ кик» с насыщением.

Магнитная эквивалентная схема, целиком заданная чертежами (граф ветвей, сечения и длины из
`saturation_branches.csv`), плюс единственное отличие от обычной МЭС — базис кика зависит от
состояния железа:

    Φ = решение цепи(I)                    нелинейное, кривая намагничивания в железных ветвях
    σ = Σ|Φ|·min(|Φ|/S/B_max, 1) / Σ|Φ|    потоко-взвешенная доля насыщения, 0 … 1
    кик = Σ_класс (w + σ·v) · Φ_класс      направление трубки — из координат наконечников

Смысл поправки: с ростом возбуждения меняется не только величина поля, но и его ФОРМА, поэтому
свёртка потоков в кик не может иметь постоянные веса. Веса при этом остаются ЛИНЕЙНЫМИ, так что
обе половины (w и v) находятся точным МНК, а не поиском: свободных чисел 12 у корректора, 4 у DQ,
и ни одного нелинейного — вся геометрия взята из чертежей.

Точность в ТОКЕ (|λ−1|, см. `current_accuracy.py`), база → МЭС-σ:
    корректор main         13.1% → 5.2%      комбинации 5.1% → 3.7%
    DQ квадруполь           4.6% → 4.3%      DQ диполь 36.4% → 18.3%

Поправка работает там, где железо насыщается, и не трогает то, где насыщения нет.

    py -3.14 mec_sigma.py --fit              # подогнать и записать веса
    py -3.14 mec_sigma.py --show             # показать записанные веса
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import Counter

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import branch_network as bn

WEIGHTS = os.path.join("workspace", "tables", "mec_sigma_weights.csv")
MODELS = (("corrector", 1), ("DQ", 1), ("DQ", 2))
_CACHE = {}


def state(cir, currents, p=None):
    """Базис классов зазоров и скалярная доля насыщения σ — за одно решение цепи."""
    p = np.zeros(1) if p is None else p
    lens, perms, secs = bn.unpack(p, cir)
    phi, _ = cir.solve(currents, lens, perms, secs=secs)
    b = np.array([sum(phi[x] * cir.direction(x) for x in cir.gaps if cir.sym[x] == c)
                  for c in cir.kick_classes], dtype=complex)
    num = den = 0.0
    for x in cir.iron:
        S = secs.get(cir.sym[x], 1.0)
        f = abs(phi[x])
        num += f * min(f / S / cir.Barr[-1], 1.0)
        den += f
    return b, (num / den if den else 0.0)


def basis(cir, currents, augment=True, p=None):
    b, s = state(cir, currents, p)
    return np.concatenate([b, s * b]) if augment else b


def fit_weights(cir, rows, augment=True, wts=None):
    """Веса точным МНК: они линейны, искать их поиском незачем."""
    if wts is None:
        cnt = Counter(tuple(sorted(t[0])) for t in rows)
        wts = {k: 1.0 / np.sqrt(v) for k, v in cnt.items()}
    A = np.array([basis(cir, c, augment) for c, _, _, _ in rows])
    y = np.array([m for _, m, _, _ in rows], dtype=complex)
    den = (np.abs(y) + 0.05 * float(np.max(np.abs(y)))) \
        / np.array([wts[tuple(sorted(c))] for c, _, _, _ in rows])
    M = np.vstack([np.real(A) / den[:, None], np.imag(A) / den[:, None]])
    v = np.concatenate([np.real(y) / den, np.imag(y) / den])
    # Регуляризация: базисные векторы почти коллинеарны, без неё МНК находит гигантские
    # взаимно сокращающиеся веса, при которых главная обмотка предсказывается нулём.
    lam = 1e-3 * float(np.sqrt((M ** 2).sum() / M.shape[0]))
    M = np.vstack([M, lam * np.eye(M.shape[1])])
    v = np.concatenate([v, np.zeros(M.shape[1])])
    return np.linalg.lstsq(M, v, rcond=None)[0]


def train(element, order):
    """Подгонка ТОЛЬКО на одиночных обмотках: комбинации остаются независимой проверкой."""
    cir = bn.Circuit(element)
    if cir.geom is None:
        raise SystemExit(f"{element}: геометрия задана не полностью — МЭС-σ требует чертежей")
    maps = bn.read_maps(element, order)
    single = [t for t in maps if len(t[0]) == 1]
    bn.calibrate_senses(cir, single, np.zeros(1))
    return cir, fit_weights(cir, single), maps


def save(rows, path=WEIGHTS):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["element", "order", "kick_class", "w", "w_sigma", "sense"])
        for r in rows:
            w.writerow(r)


def load(path=WEIGHTS):
    """Готовая модель: (цепь, веса) по элементу и гармонике, из записанной таблицы."""
    if not os.path.exists(path):
        raise SystemExit(f"нет {path} — сначала `py -3.14 mec_sigma.py --fit`")
    out = {}
    for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
        out.setdefault((r["element"], int(r["order"])), []).append(r)
    models = {}
    for key, rows in out.items():
        cir = bn.Circuit(key[0])
        by = {r["kick_class"]: r for r in rows}
        w = np.array([float(by[c]["w"]) for c in cir.kick_classes]
                     + [float(by[c]["w_sigma"]) for c in cir.kick_classes])
        for r in rows:
            if r["sense"]:
                for kv in r["sense"].split("|"):
                    if kv:
                        wd, s = kv.split("=")
                        cir.sense[wd] = float(s)
        models[key] = (cir, w)
    return models


def kick(element, order, currents, models=None):
    """Кик как комплексное число (вертикальная и горизонтальная составляющие сразу)."""
    if models is None:
        if not _CACHE:
            _CACHE.update(load())
        models = _CACHE
    cir, w = models[(element, order)]
    return complex(basis(cir, currents) @ w)


def current_scale(element, order, currents, target, models=None):
    """Во сколько раз изменить ТОК, чтобы получить заданный кик. Основа обратной задачи."""
    u = target / abs(target)

    def f(lam):
        pr = kick(element, order, {k: v * lam for k, v in currents.items()}, models)
        return float(np.real(pr * np.conj(u))) - abs(target)

    lo, hi = 1e-3, 50.0
    if f(lo) * f(hi) > 0:
        return float("nan")
    return brentq(f, lo, hi, xtol=1e-10, rtol=1e-12)


def cmd_fit(path=WEIGHTS):
    rows = []
    for el, order in MODELS:
        cir, w, maps = train(el, order)
        n = len(cir.kick_classes)
        sense = "|".join(f"{k}={v:g}" for k, v in sorted(cir.sense.items()))
        for i, c in enumerate(cir.kick_classes):
            rows.append([el, order, c, f"{w[i]:.10g}", f"{w[n + i]:.10g}",
                         sense if i == 0 else ""])
        single = [t for t in maps if len(t[0]) == 1]
        combo = [t for t in maps if len(t[0]) > 1]

        def err(g):
            if not g:
                return float("nan")
            return float(np.mean([100 * abs(complex(basis(cir, c) @ w) - m) / abs(m)
                                  for c, m, _, _ in g]))
        print(f"  {el} гармоника {order}: свободных чисел {len(w):2d}, "
              f"по кику одиночные {err(single):.1f}%  комбинации {err(combo):.1f}%")
    save(rows, path)
    print(f"веса -> {path}")


def cmd_show(path=WEIGHTS):
    for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
        print(f"  {r['element']:10} n={r['order']} {r['kick_class']:16} "
              f"w={float(r['w']):+12.5g}  w·σ={float(r['w_sigma']):+12.5g}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fit", action="store_true", help="подогнать и записать веса")
    ap.add_argument("--show", action="store_true", help="показать записанные веса")
    a = ap.parse_args(argv)
    if a.fit:
        cmd_fit()
    elif a.show:
        cmd_show()
    else:
        ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
