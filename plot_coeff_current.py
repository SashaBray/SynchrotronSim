# -*- coding: utf-8 -*-
"""График «физкоэффициент от тока» по снятым картам, с разметкой проектных уровней.

Эффективная длина элемента — ПРОЕКТНАЯ величина и не меняется, поэтому коэффициент определён
однозначно: это интеграл поля вдоль пучка, делённый на L_eff. Именно так он и считается в
датасетах `saturation_model` (подгонка по траекториям), и именно эта величина стоит в решётке
как `mat_G`.

На график наносятся:
  • точки — снятые карты (то, что посчитано, без сглаживания и интерполяции);
  • горизонтальные уровни — проектные градиенты по ТИПАМ УСТРОЙСТВ решётки, с их длинами и
    количеством в кольце;
  • потолок — максимум достижимого коэффициента по снятым картам.

    py -3.14 plot_coeff_current.py                    # квадруполи
    py -3.14 plot_coeff_current.py --element SP --key S
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

TURNS = {"QP_sh": 68, "QP_L": 85, "OP": 34, "SP": 41, "DQ": 53}
KEY = {"QP_sh": "G", "QP_L": "G", "OP": "O", "SP": "S", "DQ": "G"}
UNITS = {"B": "Тл", "G": "Тл/м", "S": "Тл/м²", "O": "Тл/м³"}


def curve(element, key, winding=None):
    """Точки (ампер-витки, коэффициент) по одиночным картам обмотки."""
    d = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))
    out = []
    for r in d["records"]:
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not r["coeffs"] or key not in r["coeffs"]:
            continue
        if winding and r.get("winding") != winding:
            continue
        if "+" in (r.get("winding") or ""):
            continue                                   # комбинации на кривую обмотки не ложатся
        out.append((ni, abs(r["coeffs"][key])))
    return sorted(out), d["records"][0].get("Leff")


def first_crossing(x, y, level):
    """Ток, при котором коэффициент ВПЕРВЫЕ достигает заданного уровня.

    Именно первое пересечение и есть рабочая точка. Простая интерполяция здесь непригодна: у QP_sh
    кривая НЕМОНОТОННА — выше 235 А градиент идёт вниз, и `np.interp` по неупорядоченному массиву
    выдавала точку на нисходящей ветви (279 А вместо ~190 А). Ток за максимумом рабочим быть не
    может: тот же градиент там достигается меньшим током."""
    for i in range(1, len(y)):
        if (y[i - 1] - level) * (y[i] - level) <= 0 and y[i] != y[i - 1]:
            t = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(x[i - 1] + t * (x[i] - x[i - 1]))
    return None


def design_levels(element, key):
    """Проектные уровни по типам устройств решётки 1: {(тип, L_eff): (значение, количество)}."""
    col = {"G": "mat_G", "B": "mat_B", "S": "mat_S", "O": "mat_O"}[key]
    rows = list(csv.DictReader(open(os.path.join("workspace", "tables", "lattice_configs.csv"),
                                    encoding="utf-8-sig"), delimiter=";"))
    c = collections.Counter()
    for r in rows:
        if str(r.get("lattice_id")) != "1" or not r.get("device_type_id", "").startswith(element):
            continue
        try:
            v = abs(float(r.get(col) or 0))
            L = float(r.get("eff_length") or 0)
        except ValueError:
            continue
        if v > 1e-9:
            c[(r["device_type_id"], round(L, 4), round(v, 5))] += 1
    # оставляем только уровни, встречающиеся заметное число раз: единичные это отстройки оптики
    tot = sum(c.values())
    return {k: n for k, n in c.items() if n >= max(3, 0.02 * tot)}


def main(argv=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", action="append", default=None)
    ap.add_argument("--key", default=None, help="физкоэффициент: G, B, S, O")
    ap.add_argument("--out", default=os.path.join("workspace", "results", "gradient_vs_current.png"))
    a = ap.parse_args(argv)
    els = a.element or ["QP_sh", "QP_L"]

    fig, axes = plt.subplots(1, len(els), figsize=(7.2 * len(els), 5.4))
    axes = np.atleast_1d(axes)
    for ax, el in zip(axes, els):
        key = a.key or KEY[el]
        pts, leff = curve(el, key)
        if not pts:
            ax.set_title(f"{el}: нет данных"); continue
        n = TURNS.get(el, 1)
        x = np.array([p[0] for p in pts]) / n
        y = np.array([p[1] for p in pts])
        ax.plot(x, y, "-o", color="#1f4e79", ms=5, lw=1.6, zorder=3,
                label=f"снятые карты ({len(pts)} шт)")
        ceil = y.max()
        ax.axhline(ceil, color="#c0392b", ls="--", lw=1.2, zorder=2)
        falls = y.argmax() < len(y) - 1
        ceil_txt = (f"потолок {ceil:.4f} при {x[y.argmax()]:.0f} А"
                    + (" — дальше кривая идёт ВНИЗ" if falls else ""))
        colors = ["#e8a33d", "#4a9d5f", "#8e5ea8", "#1a7f8e", "#c0392b"]
        levels = sorted(design_levels(el, key).items(), key=lambda t: t[0][2])
        rows_txt = []
        for i, ((dev, L, v), cnt) in enumerate(levels):
            c = colors[i % len(colors)]
            ax.axhline(v, color=c, lw=1.1, alpha=0.8, zorder=1)
            lab = f"{dev.replace(el + '_1A_', ''):5} L={L:.3f} м ×{cnt:3d}  G={v:.4f}"
            xi = first_crossing(x, y, v)
            if xi is None:
                lab += "   ВЫШЕ ПОТОЛКА"
            else:
                ax.plot([xi], [v], "D", color=c, ms=7, zorder=4, mec="white", mew=0.8)
                near = 100.0 * (1.0 - v / ceil)
                lab += f"   →  {xi:.0f} А" + (f"   (до потолка {near:.1f}%)" if near < 8 else "")
            rows_txt.append((lab, c))
        # Подписи уровней — отдельным блоком, а не у самих линий: у QP_sh три уровня L212 идут
        # вплотную и надписи налезали друг на друга.
        rows_txt.append((ceil_txt, "#c0392b"))
        for j, (lab, c) in enumerate(rows_txt):
            ax.annotate(lab, xy=(0.02, 0.97 - 0.045 * j), xycoords="axes fraction",
                        fontsize=8.5, color=c, family="monospace", va="top")
        ax.set_xlabel(f"ток обмотки, А  ({n} витков на катушку)")
        ax.set_ylabel(f"{key}, {UNITS.get(key, '')}")
        ax.set_title(f"{el}: {key} от тока" + (f"   (L_eff карты {leff:.3f} м)" if leff else ""))
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=9)
        ax.set_ylim(0, max(ceil, max(v for _, _, v in design_levels(el, key))) * 1.12)
    fig.suptitle("Градиент квадруполей в зависимости от тока: снятые карты и проектные уровни "
                 "решётки", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=140)
    print("график ->", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
