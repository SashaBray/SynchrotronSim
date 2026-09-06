# -*- coding: utf-8 -*-
"""Рёбра утечки в графе DQ: шунт полюса и обход ярма между катушками.

Две разные утечки, обе воздушные и обе дают поле у пучка:

  `lk_i`  — шунт полюса: воздух от угла ярма n_i прямо к наконечнику p_i, в обход железа полюса.
            Часть потока катушки не идёт по железу, поэтому полюс нагружается сверх полезного
            потока и, насытившись, забирает добавку себе. На однообмоточных магнитах именно это
            давало плато.

  `by_ij` — обход ярма: воздух между углами n_i и n_j, минуя железо ярма. Когда ярмо насыщается,
            поток идёт мимо него — в том числе к соседним магнитам.

Обе получают класс, начинающийся с `gap_`, поэтому их поток входит в базис кика наравне с
зазорами между наконечниками. Это существенно: до правки координат в `branch_network` вклад таких
ветвей был тождественно нулевым, и гипотезу нельзя было проверить.

Сечения из чертежей неизвестны, поэтому подбираются по данным; длины берутся геометрические —
у шунта длина полюса, у обхода длина участка ярма. Скрипт ИДЕМПОТЕНТЕН: прежние рёбра утечки
удаляются перед добавлением, так что его можно звать многократно с разными сечениями.

    py -3.14 dq_leakage_edges.py --shunt 0.03 --bypass 0.01
    py -3.14 dq_leakage_edges.py --shunt 0 --bypass 0        # убрать утечки
"""
from __future__ import annotations

import argparse
import csv
import io
import itertools
import os

BR = os.path.join("workspace", "tables", "saturation_branches.csv")
ELEMENT = "DQ"


def rewrite(shunt, bypass, path=BR, element=ELEMENT):
    """Переписать рёбра утечки элемента. Возврат: число ветвей элемента после правки."""
    rows = list(csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"))
    head = list(rows[0])
    rows = [r for r in rows
            if not (r["element"] == element and r["branch"].startswith(("lk_", "by_")))]

    poles = {r["branch"]: r for r in rows
             if r["element"] == element and r["branch"].startswith("pole")}
    if not poles:
        raise SystemExit(f"{element}: ветвей полюсов не найдено — нечего шунтировать")
    corner = {p["node_from"]: (p["pos_from"], p["z_from"]) for p in poles.values()}
    yoke_len = next((r["len_m"] for r in rows
                     if r["element"] == element and r["sym_class"] == "yoke"), "0.14")

    new = []
    for name, p in sorted(poles.items()):
        if shunt <= 0:
            break
        new.append({c: "" for c in head} | {
            "element": element, "branch": f"lk_{p['node_from']}", "material": "air",
            "S_m2": f"{shunt:.6g}", "len_m": p["len_m"], "origin": "fit",
            "comment": "шунт полюса: воздух от ярма к наконечнику мимо железа",
            "node_from": p["node_from"], "node_to": p["node_to"], "kind": "gap",
            "sym_class": "gap_poleleak", "mirror": "0",
            "pos_from": p["pos_from"], "pos_to": p["pos_to"],
            "z_from": p["z_from"], "z_to": p["z_to"]})
    for a, b in itertools.combinations(sorted(corner), 2):
        if bypass <= 0:
            break
        new.append({c: "" for c in head} | {
            "element": element, "branch": f"by_{a}{b}", "material": "air",
            "S_m2": f"{bypass:.6g}", "len_m": yoke_len, "origin": "fit",
            "comment": "обход ярма: воздух между катушками мимо железа",
            "node_from": a, "node_to": b, "kind": "gap", "sym_class": "gap_bypass",
            "mirror": "0", "pos_from": corner[a][0], "pos_to": corner[b][0],
            "z_from": corner[a][1], "z_to": corner[b][1]})

    i0 = next(k for k, r in enumerate(rows) if r["element"] == element)
    rows = rows[:i0] + new + rows[i0:]
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, head, delimiter=";")
        w.writeheader()
        w.writerows(rows)
    return sum(1 for r in rows if r["element"] == element)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shunt", type=float, default=0.0, help="сечение шунта полюса, кв.м")
    ap.add_argument("--bypass", type=float, default=0.0, help="сечение обхода ярма, кв.м")
    ap.add_argument("--element", default=ELEMENT)
    a = ap.parse_args(argv)
    n = rewrite(a.shunt, a.bypass, element=a.element)
    print(f"{a.element}: шунт {a.shunt:g} кв.м, обход {a.bypass:g} кв.м; ветвей стало {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
