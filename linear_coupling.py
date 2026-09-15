# -*- coding: utf-8 -*-
"""Линейная матрица «ток → кик» многообмоточного магнита в рабочей области — модель корректора.

Корректор — тонкий прибор: штатно он работает при малых токах (номинал 2 А) и насыщения не
достигает. Поэтому связь «токи обмоток → физкоэффициенты» описывается у него линейной матрицей,
подогнанной на картах рабочей области — с токами до 2.5 А на катушку. Модель насыщения для него не
нужна: за коленом (5–10 А) корректор не работает.

Матрица хранится как ИНТЕГРАЛ момента на ампер-виток: это ровно то, что даёт карта, и от принятой
эффективной длины оно не зависит. Физкоэффициент получается делением на l_eff исполнения. Знаковое
соглашение и множитель (n−1)! — те же, что у решателя `saturation_solve`; радиус пробы каждого
порядка — `branch_network.PROBE_R`.

Подгонка ведётся ТОЛЬКО по одиночным картам: комбинации в неё не входят и проверяют суперпозицию.
Компоненты, запрещённые симметрией обмотки (скошенные у main), закреплены нулём.

    py -3.14 linear_coupling.py corrector          # подогнать, записать в БД, показать точность
"""
from __future__ import annotations

import argparse
import csv
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

TABLE = os.path.join("workspace", "tables", "linear_coupling.csv")
RANGE_A = {"corrector": 2.5}          # верх рабочей области, А на катушку
SYM = {(1, 0): "B", (1, 1): "Bv", (2, 0): "G", (2, 1): "Gv", (3, 0): "S", (3, 1): "Sv"}
ORDER_OF = {s: nc for nc, s in SYM.items()}
UNITS = {1: "Тл·м/(А·в)", 2: "Тл/(А·в)", 3: "Тл/м/(А·в)"}
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0}


def turns(element):
    """Витки на катушку по обмоткам — из hardware_library, как у решателя."""
    import saturation_solve as ss
    t, _ = ss._hw(element)
    return {w: v for w, v in t.items() if w}


def _integral(n, part):
    """Часть комплексной гармоники карты → интеграл физкоэффициента (соглашение решателя)."""
    import branch_network as bn
    return MAT_SIGN[n] * math.factorial(n - 1) / bn.probe_r(n) ** (n - 1) * part


def _allowed(element, winding, comp):
    import coupling_model as cm
    return cm.allowed(element, winding, comp)


def maps_in_range(element, n):
    """Карты с токами всех обмоток не выше рабочей границы: [(ампер-витки, гармоника, split, id)]."""
    import branch_network as bn
    tp = turns(element)
    lim = RANGE_A[element]
    return [t for t in bn.read_maps(element, n)
            if all(abs(v) / tp.get(w, 1.0) <= lim + 1e-9 for w, v in t[0].items())]


def fit(element):
    """Матрица {символ: {обмотка: интеграл на ампер-виток}} по одиночным картам рабочей области."""
    ws = sorted(turns(element))
    K, used = {}, 0
    for n in (1, 2, 3):
        rows = [t for t in maps_in_range(element, n) if len(t[0]) == 1]
        used = max(used, len(rows))
        for comp in (0, 1):
            cols = [w for w in ws if _allowed(element, w, comp)]
            A = np.array([[t[0].get(w, 0.0) for w in cols] for t in rows], float)
            y = np.array([_integral(n, complex(t[1]).real if comp else complex(t[1]).imag)
                          for t in rows], float)
            sol = np.linalg.lstsq(A, y, rcond=None)[0]
            K[SYM[(n, comp)]] = {w: (float(sol[cols.index(w)]) if w in cols else 0.0) for w in ws}
    return K, used


def write(element, K, n_maps, path=TABLE):
    """Строки таблицы: элемент; выход; обмотка; значение; единицы; граница, А; радиус, мм; карт."""
    import branch_network as bn
    keep = []
    if os.path.exists(path):
        keep = [r for r in csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";")
                if r["element"] != element]
    head = ["element", "output", "winding", "per_amp_turn", "units", "range_a", "probe_r_mm",
            "n_single_maps"]
    rows = list(keep)
    for sym in SYM.values():
        n = ORDER_OF[sym][0]
        for w, v in K[sym].items():
            rows.append(dict(element=element, output=sym, winding=w, per_amp_turn=f"{v:.10g}",
                             units=UNITS[n], range_a=f"{RANGE_A[element]:g}",
                             probe_r_mm=f"{bn.probe_r(n) * 1e3:g}", n_single_maps=n_maps))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=head, delimiter=";")
        wr.writeheader()
        wr.writerows(rows)


def load(element, path=TABLE):
    """{символ: {обмотка: интеграл на ампер-виток}} из БД."""
    if not os.path.exists(path):
        raise SystemExit(f"нет таблицы {path} — запустите linear_coupling.py {element}")
    K = {}
    for r in csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"):
        if r["element"] == element:
            K.setdefault(r["output"], {})[r["winding"]] = float(r["per_amp_turn"])
    if not K:
        raise SystemExit(f"в {path} нет матрицы для {element}")
    return K


def integrals(element, currents_ni, K=None):
    """Ампер-витки обмоток → интегралы физкоэффициентов {символ: значение}."""
    K = K or load(element)
    return {sym: sum(row.get(w, 0.0) * ni for w, ni in currents_ni.items())
            for sym, row in K.items()}


def accuracy(element, K=None):
    """Расхождение по кику Σ|Δ|/Σ|изм|, % — по порядку (комплексно: нормальная и скошенная вместе).

    Возврат {(порядок, «одиночные»|«комбинации»): (%, число карт)}."""
    K = K or load(element)
    out = {}
    for n in (1, 2, 3):
        rows = maps_in_range(element, n)
        for group, sel in (("одиночные", [t for t in rows if len(t[0]) == 1]),
                           ("комбинации", [t for t in rows if len(t[0]) > 1])):
            if not sel:
                continue
            num = den = 0.0
            for cur, h, _, _ in sel:
                I = integrals(element, cur, K)
                pred = complex(I[SYM[(n, 1)]], I[SYM[(n, 0)]])
                h = complex(h)
                meas = complex(_integral(n, h.real), _integral(n, h.imag))
                num += abs(pred - meas)
                den += abs(meas)
            out[(n, group)] = (100 * num / den if den else float("nan"), len(sel))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="corrector")
    a = ap.parse_args(argv)
    K, n_maps = fit(a.element)
    write(a.element, K, n_maps)
    print(f"{a.element}: матрица по {n_maps} одиночным картам (токи ≤ {RANGE_A[a.element]:g} А) "
          f"-> {TABLE}")
    for (n, g), (e, k) in sorted(accuracy(a.element).items()):
        print(f"   {'BGS'[n - 1]}  {g:11s} {e:6.2f}%   ({k} карт)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
