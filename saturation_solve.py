# -*- coding: utf-8 -*-
"""Единый расчёт «ток ↔ физкоэффициент» по моделям насыщения из БД, прямая и обратная задача.

Покрывает шесть железных электромагнитов и все их исполнения (длины). Однообмоточные считаются
графом магнитной цепи с утечкой (`single_winding_model`), многообмоточные — графом ветвей плюс
матрицей связи «потоки → кик» (`coupling_model`). Постоянные магниты и внешние s2_* не поддержаны —
у них нет токовой модели.

Параметры графа подобраны отдельно для каждого исполнения: растяжение магнита меняет соотношение
железа и утечки, и на полке одна подгонка на все длины не переносится. Поэтому `device` — не
косметика: без него берётся базовая строка, снятая на картах прототипа.

    py -3.14 saturation_solve.py --forward QP_sh --currents 60 --device QP_sh_1A_L212
    py -3.14 saturation_solve.py --inverse DQ --targets B=0.262,G=-31.2 --fix I3=0
    py -3.14 saturation_solve.py --forward SP --device SP_I2_1A_L200 --currents I2=7

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
BH = os.path.join("workspace", "materials", "steel_1008_bh.csv")
R_PROBE = 8e-3
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0, 4: 1.0}
ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}
_H = None


def _hof():
    """Обратная кривая намагничивания H(B) стали — та же, на которой подбирались параметры."""
    global _H
    if _H is None:
        import single_winding_model as sw
        _H = sw.material(BH)
    return _H


def _rows():
    return [r for r in csv.DictReader(io.open(MODELS, encoding="utf-8-sig"), delimiter=";")]


def _models():
    """Базовая строка каждого элемента: {элемент: строка}. Исполнение выбирается через `_row`."""
    out = {}
    for r in _rows():
        if not r["device"]:
            out.setdefault(r["element"], r)
    return out


def _row(element, winding=None, device=None):
    """Строка модели под исполнение и обмотку: точное совпадение, иначе базовая строка обмотки."""
    rows = [r for r in _rows() if r["element"] == element]
    if not rows:
        raise SystemExit(f"нет модели насыщения для {element} — запустите build_saturation_models.py")
    for r in rows:
        if device and r["device"] == device and (winding is None or r["winding"] == winding):
            return r
    for r in rows:
        if not r["device"] and (winding is None or r["winding"] == winding):
            return r
    return rows[0]


def _par(row):
    par = np.array([float(v) for v in row["par"].split()], float)
    return par, int(row["nseries"])


def _hw(element):
    """Число витков обмоток и l_eff исполнений из hardware_library.

    Ключ обмотки берётся из идентификатора исполнения; у однообмоточных его там нет, и обмотка
    зовётся так же, как в реестре карт (её имя приходит из строки модели)."""
    import re
    turns, leff = {}, {}
    for r in csv.DictReader(io.open(HWLIB, encoding="utf-8-sig"), delimiter=";"):
        d = r["device_type_id"]
        if not d.startswith(element + "_"):
            continue
        m = re.search(r"_(I\d+|c\d+|main)_", d)
        w = m.group(1) if m else None
        if (r.get("turns_per_coil") or "").strip():
            turns.setdefault(w, float(r["turns_per_coil"]))
        if (r.get("eff_length") or "").strip():
            leff.setdefault(d, float(r["eff_length"]))
    return turns, leff


def _turns(element, winding):
    """Витки на катушку обмотки: у однообмоточных имя обмотки в идентификатор не входит."""
    turns, _ = _hw(element)
    if winding in turns:
        return turns[winding]
    return turns.get(None, 1.0)


def device_leff(element, device):
    _, leff = _hw(element)
    if device and device in leff:
        return leff[device]
    return float(_row(element)["leff"])


# ───────────────────────── однообмоточные ─────────────────────────

def _single_ni_of_y(row, y):
    """Ампер-витки для интеграла момента y — замкнутой формулой графа с утечкой."""
    import single_winding_model as sw
    par, ns = _par(row)
    return float(sw.graph_mmf(par, np.atleast_1d(abs(y)), _hof(), ns)[0])


def _single_y_of_ni(row, ni):
    """Интеграл момента по ампер-виткам: обращение монотонной зависимости делением отрезка."""
    if ni == 0:
        return 0.0
    s = 1.0 if ni > 0 else -1.0
    ni = abs(ni)
    lo, hi = 0.0, 1.0
    while _single_ni_of_y(row, hi) < ni and hi < 1e9:
        hi *= 2.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if _single_ni_of_y(row, mid) < ni:
            lo = mid
        else:
            hi = mid
    return s * 0.5 * (lo + hi)


# ───────────────────────── многообмоточные ─────────────────────────

def _graph_coef(element, currents_ni, order, leff_ref, comp=0, mat_sign=True):
    """Физкоэффициент по токам в ампер-витках: comp=0 — нормальная гармоника, 1 — скошенная.

    `mat_sign=False` отдаёт величину в сырой конвенции карт, без знакового соглашения столбцов
    mat_* — в ней сравнивают модель с гармониками, снятыми прямо с карт."""
    import coupling_model as cm
    k = cm.predict(element, order, currents_ni)
    part = k.real if comp else k.imag
    sgn = MAT_SIGN[order] if mat_sign else 1.0
    return sgn * math.factorial(order - 1) / R_PROBE ** (order - 1) * part / leff_ref


def graph_moment(element, currents_A, order, comp=0, device=None, mat_sign=True):
    """Один момент многообмоточного магнита: (физкоэффициент, интеграл по длине исполнения).

    Отдельно от `forward`, потому что тому нужны только нормальные B, G, S, а корректор описан и
    скошенными компонентами — их и берут кривые отчёта."""
    m = _row(element, device=device)
    leff = device_leff(element, device)
    ni = {w: currents_A[w] * _turns(element, w) for w in currents_A}
    g = _graph_coef(element, ni, order, float(m["leff"]), comp, mat_sign)
    return g, g * leff


# ───────────────────────── публичный интерфейс ─────────────────────────

def forward(element, currents_A, device=None):
    """Токи обмоток (А) → физкоэффициенты. Возврат: {символ: (g, интеграл g·l_eff)}.

    У однообмоточных модель даёт сам интеграл момента, а физкоэффициент получается делением на
    l_eff той подгонки, из которой взяты параметры, — так обе величины остаются согласованными."""
    m = _row(element, device=device)
    if m["model"].startswith("single"):
        w = next(iter(currents_A))
        r = _row(element, winding=w, device=device)
        ni = currents_A[w] * _turns(element, w)
        y = _single_y_of_ni(r, ni)
        leff = float(r["leff"])
        return {r["coefficient"]: (y / leff, y)}
    leff = device_leff(element, device)
    ni = {w: currents_A[w] * _turns(element, w) for w in currents_A}
    out = {}
    for sym, n in (("B", 1), ("G", 2), ("S", 3)):
        g = _graph_coef(element, ni, n, float(m["leff"]))
        out[sym] = (g, g * leff)
    return out


def inverse(element, targets, fixed=None, device=None, winding=None):
    """Физкоэффициенты → токи обмоток (А). `targets` и `fixed` — по физкоэффициентам/токам."""
    m = _row(element, device=device)
    if m["model"].startswith("single"):
        r = _row(element, winding=winding, device=device)
        sym = r["coefficient"]
        w = r["winding"]
        ni = _single_ni_of_y(r, targets[sym] * float(r["leff"]))
        return {w: math.copysign(ni, targets[sym]) / _turns(element, w)}

    from scipy.optimize import fsolve
    import coupling_model as cm
    cm.train(element)
    turns, _ = _hw(element)
    fixed = fixed or {}
    free = [w for w in sorted(k for k in turns if k) if w not in fixed]
    tgt = [(ORDER[s], v) for s, v in targets.items()]
    if len(tgt) != len(free):
        raise SystemExit(f"{element}: целей {len(tgt)}, свободных обмоток {len(free)} — "
                         "должно совпадать (закрепите лишние через fix)")
    lr = float(m["leff"])

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
    ap.add_argument("--device", help="исполнение (для l_eff и параметров), напр. QP_sh_1A_L162")
    a = ap.parse_args(argv)

    def parse(s, element=None):
        if s is None:
            return {}
        if "=" not in s:
            return {_row(element, device=a.device)["winding"] or "I1": float(s)}
        return {k: float(v) for k, v in (p.split("=") for p in s.split(","))}

    if a.forward:
        cur = parse(a.currents, a.forward)
        res = forward(a.forward, cur, device=a.device)
        print(f"{a.forward}  токи {_fmt(cur)} А" + (f"  ({a.device})" if a.device else ""))
        for sym, (g, integ) in res.items():
            print(f"   {sym} = {g:+.4g}  (интеграл {integ:+.4g})")
    elif a.inverse:
        tgt = {k: float(v) for k, v in (p.split("=") for p in a.targets.split(","))}
        fixed = parse(a.fix, a.inverse) if a.fix else {}
        res = inverse(a.inverse, tgt, fixed=fixed, device=a.device)
        print(f"{a.inverse}  цели {_fmt(tgt)}  фикс {_fmt(fixed)}"
              + (f"  ({a.device})" if a.device else ""))
        print("   токи, А: " + _fmt({k: round(v, 3) for k, v in res.items()}))
    else:
        ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
