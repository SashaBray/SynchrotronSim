# -*- coding: utf-8 -*-
"""Сборка параметров моделей насыщения ВСЕХ железных электромагнитов в одну таблицу БД.

Модель насыщения даёт связь «ток обмотки ↔ физкоэффициент» с учётом насыщения железа. Физкоэф
(G, S, O …) зависит от сечения и ампер-витков, но НЕ от длины магнита, поэтому ядро модели одно на
элемент, а эффективная длина и число витков берутся у конкретного исполнения из `hardware_library`.
Так одна модель покрывает элемент и все его исполнения (длины).

Две формы ядра:

  однообмоточные (QP_sh, QP_L, SP, OP) — скалярное уравнение в шкале ФИЗКОЭФФИЦИЕНТА
      NI = a·g + l·H(c·g),      g — физкоэффициент, H(·) — обратная кривая намагничивания;
      параметры a, c, l от длины не зависят (сняты на базовом исполнении).

  многообмоточные (DQ, corrector) — граф магнитной цепи (`saturation_branches`) плюс матрица
      связи «потоки → кик» (`mec_sigma_weights`); физкоэф считается моделью `coupling_model`.

Таблица `saturation_models.csv` хранит всё, что нужно для прямой и обратной задачи, ссылаясь на
`hardware_library` за l_eff и N исполнения. Постоянные магниты (DP) и внешние s2_* сюда не входят:
у первых нет тока управления, у вторых нет свипа по току.

    py -3.14 build_saturation_models.py
"""
from __future__ import annotations

import csv
import io
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join("workspace", "tables", "saturation_models.csv")
FITS = os.path.join("workspace", "materials", "saturation_fits.csv")
SINGLE = {"QP_sh": ("G", "int_G"), "QP_L": ("G", "int_G"),
          "SP": ("S", "int_S"), "OP": ("O", "int_O")}
MULTI = ("DQ", "corrector")
R_PROBE = 8e-3


def _material_H():
    import saturation_model as sm
    Hof = sm._material()
    return lambda b: float(np.atleast_1d(Hof(b))[0])


def _leff(element):
    import json
    import saturation_model as sm
    recs = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))["records"]
    xs = [float(r["Leff"]) for r in recs if r.get("Leff")]
    return sum(xs) / len(xs)


def single_params(element, coef_col):
    """Параметры однообмоточного ядра в шкале ФИЗКОЭФФИЦИЕНТА: a_g, c_g, l, плюс l_eff базы.

    Из `saturation_fits` берём фит на интеграле y=∫coef·dl (c, l, A_over_naive), множитель a
    восстанавливаем по контрольной точке при малом токе. Переводим в шкалу g = y/l_eff:
        a_g = a·l_eff,  c_g = c·l_eff,  l = l (длина железа в уравнении неизменна)."""
    fit = {r["section"]: r for r in csv.DictReader(
        (l for l in io.open(FITS, encoding="utf-8-sig") if not l.startswith("#")), delimiter=";")}
    q = fit[element]
    C, L, aon = float(q["C"]), float(q["L_m"]), float(q["A_over_naive"])
    rows = sorted(csv.DictReader(io.open(os.path.join("workspace", "datasets", "currents_to_kicks",
                  f"{element}.csv"), encoding="utf-8-sig"), delimiter=";"),
                  key=lambda r: float(r[[c for c in r if c.startswith("AT_")][0]]))
    at = [c for c in rows[0] if c.startswith("AT_")][0]
    ni0, y0 = float(rows[0][at]), abs(float(rows[0][coef_col]))
    a = aon * ni0 / y0
    leff = _leff(element)
    return dict(a_g=a * leff, c_g=C * leff, l=L, leff_ref=leff, err=float(q["max_err_pct"]))


def main():
    head = ["element", "model", "coefficient", "a_g", "c_g", "l", "leff_ref", "fit_err_pct",
            "note"]
    out = []
    for el, (sym, col) in SINGLE.items():
        p = single_params(el, col)
        out.append([el, "single", sym, f'{p["a_g"]:.6g}', f'{p["c_g"]:.6g}', f'{p["l"]:.6g}',
                    f'{p["leff_ref"]:.6g}', f'{p["err"]:.2f}',
                    "NI = a_g·g + l·H(c_g·g); g — физкоэф, длинонезависимо"])
    for el in MULTI:
        out.append([el, "graph", "B,G,S", "", "", "", f'{_leff(el):.6g}', "",
                    "граф saturation_branches + веса mec_sigma_weights; расчёт coupling_model"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with io.open(OUT, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(head)
        w.writerows(out)
    print("параметры моделей ->", OUT)
    for r in out:
        print("  ", r[0], r[1], r[2], "a_g=" + r[3], "c_g=" + r[4], "l=" + r[5], "l_eff=" + r[6])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
