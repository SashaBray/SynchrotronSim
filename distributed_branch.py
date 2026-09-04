# -*- coding: utf-8 -*-
"""Распределённое железное ребро: сходимость модели при измельчении ноги полюса.

Сейчас железное ребро — СОСРЕДОТОЧЕННЫЙ элемент: одно H на всю длину, падение ℓ·H(Φ/S). Ребро либо
на колене кривой, либо нет. В настоящем магните насыщение начинается в самом узком месте и
расползается вдоль ноги, а поток по пути утекает вбок. Предельная форма — магнитная линия передачи:

    dΦ/ds = −p(s)·U(s)                утечка на единицу длины
    dU/ds = −H(Φ(s)/S(s)) + n(s)·I    МДС витков, размазанная по высоте катушки

ГОТЧА: дробить ребро с ПОСТОЯННЫМ сечением и без утечки бессмысленно — последовательные одинаковые
нелинейные элементы складываются в один тождественно (m·(ℓ/m)·H(Φ/S) = ℓ·H(Φ/S)). Разница возникает
только там, где вдоль ребра меняется сечение либо распределены МДС и утечка. У корректора верно и
то и другое: катушка сидит на ноге, а междуполюсная утечка привешена дискретными рёбрами.

Здесь нога дробится на m участков, МДС катушки делится поровну, а утечка привешивается на каждом из
m уровней с сечением, поделённым так, чтобы СУММАРНАЯ проводимость утечки не менялась. При m = 2
получается ровно нынешний граф, поэтому ряд по m — честная проверка сходимости к континууму.

    py -3.14 distributed_branch.py --steps 2 3 4 6 8 12
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import branch_network as bn
import saturation_model as sm

LEGS = ("T0", "T1", "T2", "B0", "B1", "B2")
ROOT_TOTAL = 0.170          # полная длина ноги (два участка по 85 мм в исходном графе)
BASE_LEVELS = 2             # столько уровней утечки в исходном графе


def refine(m, alpha=0.0, taper=False):
    """Граф с ногой из m участков и утечкой на m уровнях. m = BASE_LEVELS, alpha = 0 даёт исходный граф.

    alpha задаёт ПРОФИЛЬ утечки вдоль ноги: вес уровня ∝ (i/m)^alpha при сохранении суммарной
    проводимости. alpha = 0 — равномерно, alpha > 0 — прижато к наконечнику, где ноги ближе всего
    и разность потенциалов наибольшая. Каждому уровню нужен свой класс симметрии: проводимость
    хранится по классу, и при общем классе разные сечения затёрли бы друг друга."""
    g = sm.BRANCH_GRAPH["corrector"]
    new = {k: (dict(v) if isinstance(v, dict) else list(v)) for k, v in g.items()}
    new["coils"] = {w: list(bs) for w, bs in g["coils"].items()}
    new["share"] = dict(g["share"])

    root = [b for b in g["branches"] if b[-2] == "r"]                  # T1r1, T1r2, ...
    leak = [b for b in g["branches"] if b.startswith("lk_")]
    drop = set(root) | set(leak)
    new["branches"] = [b for b in g["branches"] if b not in drop]
    for d in ("edge", "kind", "sym", "geom", "pos", "mirror", "material"):
        new[d] = {k: v for k, v in g[d].items() if k not in drop}
    for w in new["coils"]:
        new["coils"][w] = [b for b in new["coils"][w] if b not in drop]
    new["share"] = {k: v for k, v in new["share"].items() if k[1] not in drop}
    new["unwound"] = [b for b in g["unwound"] if b not in drop]

    coil_of = {b: w for w, bs in g["coils"].items() for b in bs}       # ребро → обмотка
    def node(leg, i):                                                  # i = 0 низ … m верх
        return f"n{leg}" if i == 0 else (f"m{leg}" if i == m else f"q{i}{leg}")

    # Сужение ноги: сечение линейно спадает от плиты к наконечнику. Насыщение тогда начинается в
    # самом узком месте и расползается вниз, чего сосредоточенное ребро с одним S дать не может.
    S_tip = g["geom"][f"{LEGS[0]}t"][0]
    for leg in LEGS:
        proto = f"{leg}r1"
        S, _ = g["geom"][proto]
        w = coil_of.get(proto)
        for i in range(1, m + 1):
            b = f"{leg}r{i}"
            new["branches"].append(b)
            new["edge"][b] = (node(leg, i - 1), node(leg, i))
            new["kind"][b] = "iron"
            new["sym"][b] = f"{g['sym'][proto]}_S{i}" if taper else g["sym"][proto]
            new["material"][b] = g["material"][proto]
            Si = S + (S_tip - S) * (i - 0.5) / m if taper else S
            new["geom"][b] = (Si, ROOT_TOTAL / m)
            new["pos"][b] = g["pos"][proto]
            new["mirror"][b] = g["mirror"][proto]
            if w:
                new["coils"][w].append(b)
                new["share"][(w, b)] = 1.0 / m                          # катушка размазана по ноге
            else:
                new["unwound"].append(b)

    # Утечка: те же пары ног, что и раньше, но на каждом из m уровней. Сечение делится так, чтобы
    # суммарная проводимость междуполюсного рассеяния осталась прежней.
    pairs = []
    for b in leak:
        if not b.startswith("lk_m"):
            continue
        a_, c_ = g["edge"][b]
        pairs.append((a_[1:], c_[1:], g["sym"][b], g["geom"][b], g["pos"][b], g["mirror"][b]))
    wl = np.array([(i / m) ** alpha for i in range(1, m + 1)])
    wl = wl / wl.sum() * BASE_LEVELS
    for i in range(1, m + 1):
        for a_, c_, symc, (S, L), pos, mir in pairs:
            b = f"lk{i}_{a_}_{c_}"
            new["branches"].append(b)
            new["edge"][b] = (node(a_, i), node(c_, i))
            new["kind"][b] = "gap"
            new["sym"][b] = f"{symc}_L{i}"
            new["material"][b] = g["material"].get(f"lk_m{a_}_{c_}", "")
            new["geom"][b] = (S * float(wl[i - 1]), L)
            new["pos"][b] = pos
            new["mirror"][b] = mir
            new["unwound"].append(b)
    return new


def evaluate(m, order=1, alpha=0.0, taper=False, verbose=True):
    name = f"corrector_m{m}_a{alpha:.3f}_t{int(taper)}"
    sm.BRANCH_GRAPH[name] = refine(m, alpha, taper)
    cir = bn.Circuit(name)
    if cir.geom is None:
        raise SystemExit("геометрия неполна после измельчения")
    maps = bn.read_maps("corrector", order)
    single = [t for t in maps if len(t[0]) == 1]
    combo = [t for t in maps if len(t[0]) > 1]
    from collections import Counter
    cnt = Counter(tuple(sorted(t[0])) for t in single)
    wts = {k: 1.0 / np.sqrt(v) for k, v in cnt.items()}
    p = np.zeros(1)                                    # геометрия задана — свободных нелинейных нет
    bn.calibrate_senses(cir, single, p)
    w = bn.solve_weights(cir, single, p, wts)

    def err(rows):
        e = [100 * abs(bn.predict(cir, c, p, w) - y) / abs(y) for c, y, _, _ in rows]
        return (float(np.mean(e)), float(np.max(e))) if e else (float("nan"),) * 2
    es, xs = err(single); ec, xc = err(combo)
    if verbose:
        tag = "  (= нынешний граф)" if m == BASE_LEVELS and alpha == 0 else ""
        print(f"  {m:6d} {alpha:7.2f} {'да' if taper else 'нет':>7} {es:9.1f}% (макс{xs:6.1f}%) "
              f"{ec:9.1f}% (макс{xc:6.1f}%){tag}")
    return es, ec


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--steps", type=int, nargs="*", default=[2, 3, 4, 6, 8, 12])
    ap.add_argument("--order", type=int, default=1)
    ap.add_argument("--taper", action="store_true", help="сужать ногу к наконечнику")
    ap.add_argument("--alpha", type=float, nargs="*", default=None,
                    help="профили утечки; по умолчанию перебор")
    a = ap.parse_args(argv)
    print(f"=== корректор, гармоника {a.order}: измельчение ноги полюса ===")
    print("  МДС катушки и междуполюсная утечка размазываются по m участкам;")
    print("  свободных чисел всегда 6 (веса «поток → кик»), геометрия из чертежей")
    print(f"  {'m':>6} {'alpha':>7} {'сужение':>7} {'одиночные':>20} {'КОМБИНАЦИИ':>20}")
    alphas = a.alpha if a.alpha is not None else [0.0, 1.0, 2.0, 4.0, 8.0]
    for m in a.steps:
        for al in alphas:
            for tp in ((False, True) if a.taper else (False,)):
                evaluate(m, a.order, al, tp)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
