# -*- coding: utf-8 -*-
"""Обратная задача: ТОКИ ОБМОТОК по заданным кикам. Это то, ради чего строилась вся модель.

Проектные токи никто не задаёт — известны проектные КИКИ и геометрия, а токи надо посчитать.
Прямая задача (`mec_sigma.kick`) даёт кик по токам; здесь решается обратная.

ГЛАВНОЕ — ПОСТАНОВКА, а не вычисление. Число обмоток почти нигде не равно числу требуемых
величин, и задача вырождена в обе стороны:

    квадруполи, OP   1 обмотка → 1 величина      определена однозначно
    DQ               3 обмотки → 2 момента (B,G) НЕДООПРЕДЕЛЕНА: континуум решений
    корректор        5 обмоток → 1 кик           НЕДООПРЕДЕЛЕНА сильно
    SP               2 обмотки → 1 момент        НЕДООПРЕДЕЛЕНА

Без явно заявленной цели решатель молча вернёт произвольную точку из континуума — формально
верную, физически бессмысленную (например, гигантские взаимно гасящиеся токи). Поэтому цель
задаётся явно, вторым слагаемым:

    минимизировать  Σ_n |кик_n(I) − цель_n|² / масштаб_n²  +  λ · Σ_w (I_w / ном_w)²

Первое слагаемое требует попасть в кики, второе выбирает СРЕДИ попадающих решение с наименьшим
током. λ мал (по умолчанию 1e-3), поэтому на точность попадания он практически не влияет, но
вырождение снимает. `--fix` замораживает обмотку на заданном токе (например, трим в нуле).

    py -3.14 inverse_currents.py DQ --target B=2.2e-4 --target G=0.0207
    py -3.14 inverse_currents.py corrector --target B=1.2e-4 --fix c2=0 --fix c3=0 --fix c4=0
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mec_sigma as ms

# Порядок гармоники по имени физкоэффициента: кик диполя это n=1, градиент n=2 и так далее.
ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}


def windings(element, models=None):
    cir = (models or ms.load())[(element, 1)][0] if (element, 1) in (models or ms.load()) \
        else None
    if cir is None:
        import branch_network as bn
        cir = bn.Circuit(element)
    return sorted(cir.coil_branches)


def nominal(element, models=None):
    """Масштаб тока обмотки — по самой сильной карте, на которой модель училась.

    Нужен только чтобы штраф за ток был безразмерным и одинаково весил разные обмотки: у
    корректора main и боковые различаются на порядок, и без нормировки штраф давил бы одну."""
    import branch_network as bn
    out = {}
    for order in (1, 2):
        try:
            maps = bn.read_maps(element, order)
        except SystemExit:
            continue
        for cur, _, _, _ in maps:
            for w, v in cur.items():
                out[w] = max(out.get(w, 0.0), abs(v))
        if out:
            break
    return {w: (v or 1.0) for w, v in out.items()}


_CAL = {}


def calibration(element, key):
    """Множитель «гармоника → физкоэффициент»: mat_* = k · |C_n|.

    ГОТЧА, стоившая порядка величины: `mec_sigma.kick` возвращает ГАРМОНИКУ карты, а проектные
    величины заданы физкоэффициентами (mat_B в Тл, mat_G в Тл/м). Единицы разные, и без перевода
    обратная задача молча выдаёт токи в разы меньше нужных, попадая при этом в «цель» идеально.

    Множитель берётся по картам: у каждой есть и гармоника (branch_network), и физкоэффициент,
    подогнанный по траекториям (saturation_model). Разброс отношения печатается — если он велик,
    значит наблюдаемые не пропорциональны и переводу верить нельзя."""
    if (element, key) in _CAL:
        return _CAL[(element, key)]
    import branch_network as bn
    import saturation_model as sm
    harm = {t[3]: t[1] for t in bn.read_maps(element, ORDER[key])}
    d = sm.load_dataset(element)
    ks = []
    for r in (d or {}).get("records", []):
        c = harm.get(r["map_id"])
        if c is None or not r["coeffs"] or key not in r["coeffs"] or not abs(c):
            continue
        ks.append(abs(r["coeffs"][key]) / abs(c))
    if not ks:
        raise SystemExit(f"{element}/{key}: не с чем калибровать перевод гармоники в {key}")
    ks = np.array(ks)
    _CAL[(element, key)] = (float(ks.mean()), float(ks.std() / ks.mean() * 100), len(ks))
    return _CAL[(element, key)]


def calibration_by_winding(element, key):
    """То же отношение, но РАЗДЕЛЬНО по обмоткам — проверка, законен ли единый множитель.

    У DQ множитель один для всех обмоток (разброс ~2%), у корректора он различается В 2.8 РАЗА
    между main (0.00453) и боковыми (0.0016–0.0018). Значит у корректора «кик» — не одна величина:
    гармоника на радиусе и коэффициент по траекториям по-разному взвешивают обмотки, потому что у
    них разные эффективные длины. Единый перевод там применять НЕЛЬЗЯ."""
    import collections

    import branch_network as bn
    import saturation_model as sm
    harm = {t[3]: t[1] for t in bn.read_maps(element, ORDER[key])}
    d = sm.load_dataset(element)
    g = collections.defaultdict(list)
    for r in (d or {}).get("records", []):
        c = harm.get(r["map_id"])
        if c is None or not r["coeffs"] or key not in r["coeffs"] or not abs(c):
            continue
        g[r.get("winding") or "?"].append(abs(r["coeffs"][key]) / abs(c))
    return {w: (float(np.mean(v)), len(v)) for w, v in g.items()}


def solve(element, targets, fix=None, lam=1e-3, models=None, seeds=8):
    """Токи обмоток под заданные кики. targets: {'B': значение, 'G': значение, …}.

    Возвращает (токи в ампер-витках, достигнутые кики, невязка в процентах по каждой цели)."""
    models = models or ms.load()
    fix = dict(fix or {})
    ws = windings(element, models)
    free = [w for w in ws if w not in fix]
    if not free:
        raise SystemExit("все обмотки зафиксированы — решать нечего")
    nom = nominal(element, models)
    keys = list(targets)
    for k in keys:
        if (element, ORDER[k]) not in models:
            raise SystemExit(f"{element}: нет модели для гармоники {k} (n={ORDER[k]}) — "
                             f"сначала mec_sigma.py --fit")
    scale = {k: abs(targets[k]) or 1.0 for k in keys}
    cal = {k: calibration(element, k) for k in keys}

    def currents_of(x):
        cur = dict(fix)
        cur.update({w: x[i] * nom.get(w, 1.0) for i, w in enumerate(free)})
        return cur

    def resid(x):
        cur = currents_of(x)
        out = []
        for k in keys:
            pr = cal[k][0] * abs(ms.kick(element, ORDER[k], cur, models))
            # Сравниваем МОДУЛИ: знак кика задаётся намоткой и в проектной величине не участвует.
            out.append((pr - abs(targets[k])) / scale[k])
        # Штраф за ток: выбирает среди попадающих решений наименьшее по току.
        out += [np.sqrt(lam) * v for v in x]
        return np.array(out)

    best = None
    for s in range(seeds):
        rng = np.random.default_rng(s)
        x0 = rng.lognormal(0, 0.7, len(free)) * (0.3 if s else 1.0)
        try:
            r = least_squares(resid, x0, xtol=1e-14, ftol=1e-14, max_nfev=4000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    cur = currents_of(best.x)
    got = {k: cal[k][0] * abs(ms.kick(element, ORDER[k], cur, models)) for k in keys}
    err = {k: 100.0 * (got[k] - abs(targets[k])) / abs(targets[k]) for k in keys}
    return cur, got, err, cal


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element")
    ap.add_argument("--target", action="append", default=[], metavar="ИМЯ=ЗНАЧЕНИЕ",
                    help="требуемый физкоэффициент, напр. B=2.2e-4 или G=0.0207")
    ap.add_argument("--fix", action="append", default=[], metavar="ОБМОТКА=ТОК",
                    help="заморозить обмотку на этом токе (ампер-витки)")
    ap.add_argument("--lam", type=float, default=1e-3, help="вес штрафа за ток")
    a = ap.parse_args(argv)
    if not a.target:
        raise SystemExit("нужен хотя бы один --target")
    targets = {}
    for t in a.target:
        k, v = t.split("=", 1)
        targets[k.strip()] = float(v)
    fix = {}
    for t in a.fix:
        k, v = t.split("=", 1)
        fix[k.strip()] = float(v)

    models = ms.load()
    ws = windings(a.element, models)
    print(f"=== {a.element}: обмоток {len(ws)} {ws}, требуемых величин {len(targets)}")
    n_free = len(ws) - len(fix)
    if n_free > len(targets):
        print(f"  задача НЕДООПРЕДЕЛЕНА ({n_free} свободных обмоток на {len(targets)} величин): "
              f"из континуума решений берётся наименьшее по току")
    elif n_free < len(targets):
        print(f"  задача ПЕРЕОПРЕДЕЛЕНА ({n_free} обмоток на {len(targets)} величин): "
              f"точное попадание невозможно, ищется наилучшее приближение")
    cur, got, err, cal = solve(a.element, targets, fix, a.lam, models)
    for k, (v, sp, n) in cal.items():
        print(f"  перевод гармоники в {k}: множитель {v:.5g} (разброс по {n} картам {sp:.1f}%)")
        by = calibration_by_winding(a.element, k)
        singles = {w: t for w, t in by.items() if "+" not in w}
        if len(singles) > 1:
            lo = min(singles.values())[0]; hi = max(singles.values())[0]
            if hi > 1.15 * lo:
                print(f"    ВНИМАНИЕ: множитель РАЗНЫЙ у обмоток — в {hi / lo:.1f} раза "
                      f"(от {lo:.4g} до {hi:.4g}).")
                print("    Значит гармоника и траекторный коэффициент по-разному взвешивают")
                print("    обмотки (разная эффективная длина), и единый перевод НЕЗАКОНЕН.")
                print("    Ответ ниже ориентировочный; модель надо подгонять сразу на")
                print("    траекторных коэффициентах, а не переводить из гармоник.")
                for w, (t, n_) in sorted(singles.items()):
                    print(f"      {w:8} {t:.5g}  (карт {n_})")
    turns = {}
    try:
        import ansys_maps as am
        turns = am.DESIGN.get(a.element, {}).get("turns", {})
    except Exception:
        pass
    print(f"\n  {'обмотка':10} {'ампер-витков':>14} {'ампер':>10}")
    for w in ws:
        n = turns.get(w)
        amps = f"{cur.get(w, 0.0) / n:10.2f}" if n else " " * 10
        mark = "  (заморожена)" if w in fix else ""
        print(f"  {w:10} {cur.get(w, 0.0):14.1f} {amps}{mark}")
    print(f"\n  {'величина':10} {'требуется':>14} {'получено':>14} {'невязка':>10}")
    for k in targets:
        print(f"  {k:10} {abs(targets[k]):14.6g} {got[k]:14.6g} {err[k]:9.2f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
