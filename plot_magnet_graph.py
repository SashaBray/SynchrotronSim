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
    """Ветви элемента и катушки на них: (список ветвей, {ребро: [(обмотка, витки, знак)]}).

    Знак МДС катушки (`mmf_share`) выводится на рисунок: у соседних полюсов квадруполя полярность
    противоположна, значит и намотка встречная, и перепутать её ничего не стоит."""
    rows = [r for r in csv.DictReader(open(branches, encoding="utf-8-sig"), delimiter=";")
            if r["element"] == element]
    coils = collections.defaultdict(list)
    for r in csv.DictReader(open(members, encoding="utf-8-sig"), delimiter=";"):
        if r["element"] == element and r.get("branch"):
            try:
                sign = float(r.get("mmf_share") or 1)
            except ValueError:
                sign = 1.0
            coils[r["branch"]].append((r["winding"], r.get("turns_per_coil") or "?", sign))
    return rows, coils


def rule_positions(rows):
    """СХЕМАТИЧЕСКАЯ укладка по именам узлов — для элементов с регулярной структурой (корректор).

    В таблице корректора столбцы `pos_*` заняты не координатами, а номером столбца полюса (так их
    читает `branch_network`), а `z_*` пустые — координат для рисунка там нет, и вершины ложились по
    кругу, превращая схему в паутину. В реальном масштабе она не лучше: шаг полюсов 30 мм против
    413 мм ярма, наконечники слипаются в точку. Поэтому здесь укладка схематическая, как на
    эскизе от руки: ярмо слева отдельным столбцом, полюса тремя столбцами по порядку вдоль
    пластины (ближний к ярму — первый), уровни вдоль полюса равномерные.

    Возврат: {узел: (x, y)} в условных единицах. Пусто, если имена узлов правилу не отвечают."""
    colx = {"1": 1.0, "0": 2.0, "2": 3.0}          # порядок вдоль пластины, а не координата
    lvl = {"n": 3.4, "s": 2.5, "m": 1.65, "p": 0.62}
    pos = {}
    for r in rows:
        for n in (r["node_from"], r["node_to"]):
            if not n or n in pos:
                continue
            row = 1.0 if "T" in n[:3] else -1.0
            if n[0] in "pmn" and len(n) >= 3 and n[2] in colx:
                pos[n] = (colx[n[2]], row * lvl[n[0]])
            elif n.startswith("s1") and len(n) >= 4 and n[3] in colx:
                pos[n] = (colx[n[3]], row * lvl["s"])
            elif n in ("cT", "cB"):
                pos[n] = (0.0, row * 3.4)
            elif n in ("yM1", "yM2"):
                pos[n] = (0.0, 1.3 if n == "yM1" else -1.3)
    return pos


def node_positions(rows):
    """Координаты вершин и признак «схема»: из БД, где заданы; иначе по правилу имён; остаток — по кругу."""
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
    real = bool(pos)
    for n, xy in rule_positions(rows).items():
        pos.setdefault(n, xy)
    missing = [n for r in rows for n in (r["node_from"], r["node_to"]) if n and n not in pos]
    for i, n in enumerate(dict.fromkeys(missing)):
        a = 2 * math.pi * i / max(len(set(missing)), 1)
        pos[n] = (60 * math.cos(a), 60 * math.sin(a))
    return pos, real


def draw(element, out=None, annotate=True, branches=None, members=None):
    """Рисунок графа. Пути к таблицам можно передать явно — например, чтобы посмотреть вариант
    топологии, не трогая базу. Подменять модульные константы для этого НЕЛЬЗЯ: они попадают в
    значения по умолчанию при определении `load`, и подмена молча не срабатывает."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rows, coils = load(element, branches or BRANCHES, members or MEMBERS)
    if not rows:
        raise SystemExit(f"нет ветвей элемента {element} в {BRANCHES}")
    pos, real = node_positions(rows)

    fig, ax = plt.subplots(figsize=(13.5, 10.5))

    def curvature(a, b, kind):
        """Дуга вместо прямой там, где прямая пройдёт ПО ЧУЖОМУ узлу.

        Ветвь через столбец (1↔3) идёт ровно сквозь средний полюс, а диагонали через апертуру
        пересекаются все в одной точке — прямыми их подписи садятся друг на друга и на соседей."""
        (xa, ya), (xb, yb) = a, b
        dx, dy = abs(xb - xa), abs(yb - ya)
        if dx > 1.5:                                   # через столбец
            return 0.5 if ya > 0 else -0.5
        if dx < 0.05 and dy > 0.3 and kind == "gap":   # зазор вдоль железа — иначе лежит на нём
            return 0.34 if dy > 3.0 else (0.5 if dy > 1.0 else 1.1)
        if dx > 0.05 and ya * yb < 0:                  # диагональ через апертуру
            return 0.22 * (1.0 if xb > xa else -1.0) * (1.0 + 0.35 * dx)
        return 0.0

    # Параллельные рёбра (одни и те же концы) — например, воздушный шунт полюса поверх железного
    # полюса: без дуги они рисуются друг на друге, и шунт не виден.
    import collections
    pairn = collections.Counter(frozenset((r["node_from"], r["node_to"])) for r in rows)

    from matplotlib.patches import FancyArrowPatch
    for r in rows:
        a, b = pos[r["node_from"]], pos[r["node_to"]]
        st = KIND_STYLE.get(r["kind"], dict(color="#888", lw=1.0, ls=":"))
        rad = curvature(a, b, r["kind"])
        if pairn[frozenset((r["node_from"], r["node_to"]))] > 1 and r["kind"] == "gap":
            # шунт идёт теми же концами, что и железный полюс; выгибаем его в ДРУГУЮ сторону
            # (и сильнее), иначе он ложится на полюс и не виден
            rad = -1.8 * (rad if rad else 0.4)
        if rad:
            ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", shrinkA=0, shrinkB=0, zorder=2,
                                         connectionstyle=f"arc3,rad={rad}",
                                         color=st["color"], lw=st["lw"], ls=st["ls"]))
        else:
            ax.plot([a[0], b[0]], [a[1], b[1]], zorder=2, **st)
        # Середина квадратичной дуги: M + 0.5·rad·перпендикуляр(конец−начало). Знак
        # перпендикуляра именно такой, как у matplotlib в `arc3` — с обратным подписи уезжают
        # на противоположную от дуги сторону.
        mx, my = (a[0] + b[0]) / 2, (a[1] + b[1]) / 2
        if rad:
            mx += 0.5 * rad * (b[1] - a[1])
            my += 0.5 * rad * -(b[0] - a[0])
        if annotate:
            # У зазоров подпись короткая: их вдвое больше железных, и полные подписи с сечениями
            # перекрывают друг друга так, что схему не сверить — а сверяют по ней топологию.
            lab = r["branch"]
            if r["branch"] in coils:
                # Доля МДС обязательна в подписи: обмотка, разложенная по двум участкам полюса,
                # выводится дважды, и без доли это читается как ДВЕ катушки вместо одной.
                lab += "\n" + ", ".join(
                    f"{'−' if sg < 0 else '+'}{w} {t}в" + (f"×{abs(sg):g}" if abs(sg) != 1 else "")
                    for w, t, sg in coils[r["branch"]])
            if r["kind"] == "iron" and (r.get("S_m2") or "").strip():
                lab += f"\nS={float(r['S_m2']) * 1e4:.0f} см²"
            ax.annotate(lab, xy=(mx, my), ha="center", va="center", zorder=4,
                        color=st["color"], alpha=1.0 if r["kind"] == "iron" else 0.8,
                        fontsize=8.0 if r["kind"] == "iron" else 6.6,
                        bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.0))
    for n, (x, z) in pos.items():
        ax.plot([x], [z], "o", ms=13, color="white", mec="#111", mew=1.5, zorder=5)
        ax.annotate(n, xy=(x, z), fontsize=7.0, ha="center", va="center", zorder=6)

    nir = sum(1 for r in rows if r["kind"] == "iron")
    ngap = sum(1 for r in rows if r["kind"] == "gap")
    ax.set_title(f"Магнитный граф {element}: вершин {len(pos)}, ветвей {len(rows)} "
                 f"(железо {nir}, зазоры {ngap})", fontsize=12)
    if real:
        ax.set_xlabel("поперёк, мм")
        ax.set_ylabel("вертикаль, мм")
        ax.set_aspect("equal")
    else:
        ax.set_xlabel("схема (не масштаб): ярмо — левый столбец, полюса — по порядку вдоль пластины")
        ax.set_ylabel("вверх по полюсу от апертуры")
        ax.set_xticks([])
        ax.set_yticks([])
    ax.grid(alpha=0.2)
    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], label="железо", **KIND_STYLE["iron"]),
                       Line2D([], [], label="зазор", **KIND_STYLE["gap"])],
              fontsize=9, loc="upper left", bbox_to_anchor=(0.01, 0.99))
    ax.margins(x=0.09, y=0.05)
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
