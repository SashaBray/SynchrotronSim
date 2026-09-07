# -*- coding: utf-8 -*-
"""Сборка параметров моделей насыщения ВСЕХ железных электромагнитов в одну таблицу БД.

Модель насыщения даёт связь «ток обмотки ↔ физкоэффициент» с учётом насыщения железа. Две формы
ядра:

  однообмоточные (QP_sh, QP_L, SP, OP) — граф магнитной цепи с утечкой (`single_winding_model`):
      Q(y) = Rз·y + Σₖ mₖ·ℓₖ·H(y/Sₖ),   p(y) = y + (P/2)·Q(y),   NI = ℓп·H(p/Sп) + Q(y)/2,
      где y = ∫коэф·dl — интеграл ведущего момента, H(·) — обратная кривая намагничивания.
      Параметры подбираются НА КАРТАХ КАЖДОГО ИСПОЛНЕНИЯ (карты снятой длины пересчитываются на
      длину исполнения синтезом): растяжение магнита меняет соотношение железа и утечки, поэтому
      на полке одна и та же модель на разные длины не переносится.

  многообмоточные (DQ, corrector) — граф магнитной цепи (`saturation_branches`) плюс матрица
      связи «потоки → кик» (`mec_sigma_weights`); физкоэф считается моделью `coupling_model`.

Таблица `saturation_models.csv` хранит по строке на исполнение и обмотку: вектор параметров графа
(в логарифмах, как их держит `single_winding_model`), число последовательных классов ветвей и
эффективную длину, при которой подгонка сделана. Строка с пустым `device` — базовая, на снятых
картах: её берут, когда исполнение не указано. Постоянные магниты (DP) и внешние s2_* сюда не
входят: у первых нет тока управления, у вторых нет свипа по току.

    py -3.14 build_saturation_models.py
    py -3.14 build_saturation_models.py --only QP_sh
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join("workspace", "tables", "saturation_models.csv")
HWLIB = os.path.join("workspace", "tables", "hardware_library.csv")
SINGLE = {"QP_sh": "G", "QP_L": "G", "SP": "S", "OP": "O"}
MULTI = ("DQ", "corrector")
ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}


def _devices(element, windings):
    """Исполнения элемента: [(device_type_id, обмотка, l_eff)] — только те, где задано число витков.

    Число витков есть у исполнений, стоящих в решётке; строки-прототипы (без витков) описывают
    снятые карты и в таблицу исполнений не идут — для них пишется базовая строка. Имя обмотки
    берётся из идентификатора исполнения, а у однообмоточных — единственное имя из реестра карт
    (в идентификаторе его нет)."""
    import re
    out = []
    for r in csv.DictReader(io.open(HWLIB, encoding="utf-8-sig"), delimiter=";"):
        d = r["device_type_id"]
        if not d.startswith(element + "_"):
            continue
        if not (r.get("turns_per_coil") or "").strip() or not (r.get("eff_length") or "").strip():
            continue
        m = re.search(r"_(I\d+|c\d+|main)_", d)
        w = m.group(1) if m and m.group(1) in windings else (
            windings[0] if len(windings) == 1 else None)
        if w is None:
            print(f"  {d}: обмотка не опознана — пропуск")
            continue
        out.append((d, w, float(r["eff_length"])))
    return out


def _base(element):
    """Снятые карты элемента: ({обмотка: (ампер-витки, интеграл момента)}, длина отнесения)."""
    import single_winding_model as sw
    ser, _ = sw.series(element)
    return ser, sw.coef_scale(element)[1]


def _variant_series(element, winding, L, order):
    """Карты исполнения длины L: (ампер-витки, интеграл момента) синтезом из снятых карт."""
    import library_report as lr
    pts = lr.variant_kicks(element, L, order).get(winding)
    if not pts:
        return None
    ni = np.array([p[0] for p in pts], float)
    y = np.array([abs(p[1][order][0]) for p in pts], float)
    o = np.argsort(ni)
    return ni[o], y[o]


def _fit(ni, y):
    """Подгонка графа с утечкой: (вектор параметров, СКО по току %)."""
    import single_winding_model as sw
    Hof = sw.material(os.path.join("workspace", "materials", "steel_1008_bh.csv"))
    good = ni > 0
    return sw.best_fit(ni[good], y[good], Hof, nseries=2, tries=3, seeds=40)


def single_rows(element, sym):
    """Строки таблицы для однообмоточного элемента: базовая плюс по строке на исполнение."""
    order = ORDER[sym]
    rows, seen = [], set()
    base, Lb = _base(element)
    for winding in sorted(base):
        par, err = _fit(*base[winding])
        rows.append([element, "", winding, "single_graph", sym,
                     " ".join(f"{v:.10g}" for v in par), (len(par) - 4) // 2,
                     f"{Lb:.6g}", f"{err:.2f}", "граф с утечкой на снятых картах"])
    for dev, winding, L in _devices(element, sorted(base)):
        if (winding, round(L, 6)) in seen:
            continue
        seen.add((winding, round(L, 6)))
        ser = _variant_series(element, winding, L, order)
        if ser is None:
            print(f"  {dev}: нет карт обмотки {winding} — пропуск")
            continue
        par, err = _fit(*ser)
        rows.append([element, dev, winding, "single_graph", sym,
                     " ".join(f"{v:.10g}" for v in par), (len(par) - 4) // 2,
                     f"{L:.6g}", f"{err:.2f}", "граф с утечкой на картах исполнения"])
    return rows


def _leff_ref(element):
    import json
    import saturation_model as sm
    recs = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))["records"]
    xs = [float(r["Leff"]) for r in recs if r.get("Leff")]
    return sum(xs) / len(xs)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", help="пересобрать только этот элемент, остальные строки сохранить")
    a = ap.parse_args(argv)

    keep = []
    if a.only and os.path.exists(OUT):
        keep = [r for r in csv.reader(io.open(OUT, encoding="utf-8-sig"), delimiter=";")][1:]
        keep = [r for r in keep if r and r[0] != a.only]

    head = ["element", "device", "winding", "model", "coefficient", "par", "nseries", "leff",
            "fit_err_pct", "note"]
    out = list(keep)
    for el, sym in SINGLE.items():
        if a.only and el != a.only:
            continue
        print(el)
        out += single_rows(el, sym)
    for el in MULTI:
        if a.only and el != a.only:
            continue
        out.append([el, "", "", "graph", "B,G,S", "", "", f"{_leff_ref(el):.6g}", "",
                    "граф saturation_branches + веса mec_sigma_weights; расчёт coupling_model"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with io.open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(head)
        w.writerows(out)
    print("параметры моделей ->", OUT)
    for r in out:
        print("  ", r[0], r[1] or "(база)", r[2], r[3], "l_eff=" + str(r[7]),
              ("СКО " + r[8] + "%") if r[8] else "")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
