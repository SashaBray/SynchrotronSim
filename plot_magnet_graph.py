# -*- coding: utf-8 -*-
"""Рисунок магнитного графа элемента по таблице ветвей БД — чтобы схему можно было СВЕРИТЬ глазами.

Граф магнита живёт в `saturation_branches.csv` числами: вершины, рёбра, сечения, длины, классы
симметрии и координаты концов. Числа сверять глазами тяжело, а ошибка в топологии портит всё
дальнейшее молча, поэтому здесь он выводится картинкой в РЕАЛЬНОЙ геометрии (координаты из той же
таблицы), с катушками на рёбрах и сечениями в подписях.

    py -3.14 plot_magnet_graph.py DQ
    py -3.14 plot_magnet_graph.py corrector --out ...
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import sys

BRANCHES = os.path.join("workspace", "tables", "saturation_branches.csv")
MEMBERS = os.path.join("workspace", "tables", "saturation_branch_members.csv")
KIND_STYLE = {"iron": dict(color="#1f4e79", lw=2.6, ls="-"),
              "gap": dict(color="#c0392b", lw=1.3, ls="--")}


def load(element, branches=BRANCHES, members=MEMBERS):
    """Ветви элемента и катушки на них: (список ветвей, {ребро: [(обмотка, витки)]})."""
    rows = [r for r in csv.DictReader(open(branches, encoding="utf-8-sig"), delimiter=";")
            if r["element"] == element]
    coils = collections.defaultdict(list)
    for r in csv.DictReader(open(members, encoding="utf-8-sig"), delimiter=";"):
        if r["element"] == element and r.get("branch"):
            coils[r["branch"]].append((r["winding"], r.get("turns_per_coil") or "?"))
    return rows, coils


def node_positions(rows):
    """Координаты вершин из концов рёбер. У ветвей без координат вершины расставляются по кругу —
    иначе граф вырождается в точку и рисунок ничего не показывает."""
    import math
    pos = {}
    for r in rows:
        for n, x, z in ((r["node_from"], r.get("pos_from"), r.get("z_from")),
                        (r["node_to"], r.get("pos_to"), r.get("z_to"))):
            if n and n not in pos and (x or "").strip() and (z or "").strip():
                try:
                    pos[n] = (float(x), float(z))
                except ValueError:
                    pass
    missing = [n for r in rows for n in (r["node_from"], r["node_to"]) if n and n not in pos]
    for i, n in enumerate(dict.fromkeys(missing)):
        a = 2 * math.pi * i / max(len(set(missing)), 1)
        pos[n] = (60 * math.cos(a), 60 * math.sin(a))
    return pos


def draw(element, out=None, annotate=True):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows, coils = load(element)
    if not rows:
        raise SystemExit(f"нет ветвей элемента {element} в {BRANCHES}")
    pos = node_positions(rows)

    fig, ax = plt.subplots(figsize=(9.5, 8.4))
    for r in rows:
        a, b = pos[r["node_from"]], pos[r["node_to"]]
        st = KIND_STYLE.get(r["kind"], dict(color="#888", lw=1.0, ls=":"))
        ax.plot([a[0], b[0]], [a[1], b[1]], zorder=2, **st)
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        if annotate:
            lab = r["branch"]
            if r["branch"] in coils:
                lab += "\n" + ", ".join(f"{w}({t})" for w, t in coils[r["branch"]])
            if (r.get("S_m2") or "").strip():
                lab += f"\nS={float(r['S_m2']) * 1e4:.0f} см²"
            ax.annotate(lab, xy=(mx, my), fontsize=7.5, ha="center", va="center", zorder=4,
                        color=st["color"],
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.2))
    for n, (x, z) in pos.items():
        ax.plot([x], [z], "o", ms=11, color="white", mec="#111", mew=1.6, zorder=5)
        ax.annotate(n, xy=(x, z), fontsize=8, ha="center", va="center", zorder=6)

    nir = sum(1 for r in rows if r["kind"] == "iron")
    ngap = sum(1 for r in rows if r["kind"] == "gap")
    ax.set_title(f"Магнитный граф {element}: вершин {len(pos)}, ветвей {len(rows)} "
                 f"(железо {nir}, зазоры {ngap})", fontsize=12)
    ax.set_xlabel("поперёк, мм")
    ax.set_ylabel("вертикаль, мм")
    ax.set_aspect("equal")
    ax.grid(alpha=0.2)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], label="железо", **KIND_STYLE["iron"]),
                       Line2D([], [], label="зазор", **KIND_STYLE["gap"])],
              fontsize=9, loc="upper right")
    fig.tight_layout()
    out = out or os.path.join("workspace", "results", f"graph_{element}.png")
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element")
    ap.add_argument("--out")
    ap.add_argument("--no-labels", action="store_true")
    a = ap.parse_args(argv)
    print("граф ->", draw(a.element, a.out, annotate=not a.no_labels))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
