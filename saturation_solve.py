# -*- coding: utf-8 -*-
"""Единый расчёт «ток ↔ физкоэффициент» по моделям насыщения из БД, прямая и обратная задача.

Покрывает шесть железных электромагнитов и все их исполнения (длины). Ядро модели даёт физический
коэффициент g (Тл/м, Тл/м², …), который от длины не зависит; эффективная длина и число витков
берутся у исполнения из `hardware_library`. Постоянные магниты и внешние s2_* не поддержаны —
у них нет токовой модели.

    py -3.14 saturation_solve.py --forward QP_sh --currents 60
    py -3.14 saturation_solve.py --inverse DQ --targets B=0.262,G=-31.2 --fix I3=0
    py -3.14 saturation_solve.py --forward corrector --device corrector_c1_100A_L100 --currents c1=7

Как функции: `forward(element, currents_A)` и `inverse(element, targets, fixed)`.
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

MODELS = os.path.join("workspace", "tables", "saturation_models.csv")
HWLIB = os.path.join("workspace", "tables", "hardware_library.csv")
R_PROBE = 8e-3
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0, 4: 1.0}
ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}
_H = None


def material_H(b):
    global _H
    if _H is None:
        import saturation_model as sm
        Hof = sm._material()
        _H = lambda x: float(np.atleast_1d(Hof(x))[0])
    return _H(b)


def _models():
    return {r["element"]: r for r in csv.DictReader(io.open(MODELS, encoding="utf-8-sig"),
                                                    delimiter=";")}


def _hw(element):
    """Число витков обмоток и базовая l_eff элемента из hardware_library.

    N постоянно по длинам (свойство катушки), поэтому берётся у любого исполнения, где задано."""
    turns, leff = {}, {}
    import re
    for r in csv.DictReader(io.open(HWLIB, encoding="utf-8-sig"), delimiter=";"):
        d = r["device_type_id"]
        if not d.startswith(element):
            continue
        m = re.search(r"_(I\d+|c\d+|main)", d)
        w = m.group(1) if m else "I1"
        if (r.get("turns_per_coil") or "").strip():
            turns.setdefault(w, float(r["turns_per_coil"]))
        if (r.get("eff_length") or "").strip():
            leff.setdefault(d, float(r["eff_length"]))
    return turns, leff


def device_leff(element, device):
    _, leff = _hw(element)
    if device and device in leff:
        return leff[device]
    return float(_models()[element]["leff_ref"])


# ───────────────────────── однообмоточные ─────────────────────────

def _single_ni_of_g(m, g):
    """Ампер-витки для заданного физкоэффициента: NI = a_g·g + l·H(c_g·g). Прямая формула."""
    a, c, l = float(m["a_g"]), float(m["c_g"]), float(m["l"])
    return a * g + l * material_H(c * g)


def _single_g_of_ni(m, ni):
    """Физкоэффициент по ампер-виткам: обращение монотонного уравнения делением отрезка."""
    if ni == 0:
        return 0.0
    s = 1.0 if ni > 0 else -1.0
    ni = abs(ni)
    lo, hi = 0.0, 1.0
    while _single_ni_of_g(m, hi) < ni and hi < 1e12:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _single_ni_of_g(m, mid) < ni:
            lo = mid
        else:
            hi = mid
    return s * 0.5 * (lo + hi)


# ───────────────────────── многообмоточные ─────────────────────────

def _graph_coef(element, currents_ni, order, leff_ref):
    import coupling_model as cm
    k = cm.predict(element, order, currents_ni)
    return MAT_SIGN[order] * math.factorial(order - 1) / R_PROBE ** (order - 1) * k.imag / leff_ref


# ───────────────────────── публичный интерфейс ─────────────────────────

def forward(element, currents_A, device=None):
    """Токи обмоток (А) → физкоэффициенты. Возврат: {символ: (g, интеграл g·l_eff)}."""
    M = _models()
    m = M[element]
    turns, _ = _hw(element)
    leff = device_leff(element, device)
    if m["model"] == "single":
        sym = m["coefficient"]
        w = next(iter(currents_A))
        ni = currents_A[w] * turns.get(w, 1.0)
        g = _single_g_of_ni(m, ni)
        return {sym: (g, g * leff)}
    ni = {w: currents_A[w] * turns.get(w, 1.0) for w in currents_A}
    out = {}
    for sym, n in (("B", 1), ("G", 2), ("S", 3)):
        g = _graph_coef(element, ni, n, float(m["leff_ref"]))
        out[sym] = (g, g * leff)
    return out


def inverse(element, targets, fixed=None, device=None):
    """Физкоэффициенты → токи обмоток (А). `targets` и `fixed` — по физкоэффициентам/токам."""
    M = _models()
    m = M[element]
    turns, _ = _hw(element)
    if m["model"] == "single":
        sym = m["coefficient"]
        g = targets[sym]
        w = next(w for w in turns) if turns else "I1"
        return {w: _single_ni_of_g(m, g) / turns.get(w, 1.0)}

    from scipy.optimize import fsolve
    import coupling_model as cm
    cm.train(element)
    fixed = fixed or {}
    free = [w for w in sorted(turns) if w not in fixed]
    tgt = [(ORDER[s], v) for s, v in targets.items()]
    if len(tgt) != len(free):
        raise SystemExit(f"{element}: целей {len(tgt)}, свободных обмоток {len(free)} — "
                         "должно совпадать (закрепите лишние через fix)")
    lr = float(m["leff_ref"])

    def resid(x):
        cur = dict(fixed)
        cur.update({w: x[i] for i, w in enumerate(free)})
        ni = {w: cur[w] * turns[w] for w in cur}
        return [_graph_coef(element, ni, n, lr) - v for n, v in tgt]

    x = fsolve(resid, [50.0] * len(free))
    cur = dict(fixed)
    cur.update({w: float(x[i]) for i, w in enumerate(free)})
    return cur


def _fmt(d):
    return ", ".join(f"{k}={v:g}" for k, v in d.items())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--forward")
    ap.add_argument("--inverse")
    ap.add_argument("--currents", help="ток(и): «60» или «c1=7,c2=3»")
    ap.add_argument("--targets", help="цели: «G=-31.2» или «B=0.26,G=-31»")
    ap.add_argument("--fix", help="закреплённые токи: «I3=0»")
    ap.add_argument("--device", help="исполнение (для l_eff), напр. QP_sh_1A_L162")
    a = ap.parse_args(argv)

    def parse(s, one_key="I1"):
        if s is None:
            return {}
        if "=" not in s:
            return {one_key: float(s)}
        return {k: float(v) for k, v in (p.split("=") for p in s.split(","))}

    if a.forward:
        cur = parse(a.currents)
        res = forward(a.forward, cur, device=a.device)
        print(f"{a.forward}  токи {_fmt(cur)} А" + (f"  ({a.device})" if a.device else ""))
        for sym, (g, integ) in res.items():
            print(f"   {sym} = {g:+.4g}  (интеграл {integ:+.4g})")
    elif a.inverse:
        tgt = {k: float(v) for k, v in (p.split("=") for p in a.targets.split(","))}
        fixed = parse(a.fix) if a.fix else {}
        res = inverse(a.inverse, tgt, fixed=fixed, device=a.device)
        print(f"{a.inverse}  цели {_fmt(tgt)}  фикс {_fmt(fixed)}")
        print("   токи, А: " + _fmt({k: round(v, 3) for k, v in res.items()}))
    else:
        ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
