# -*- coding: utf-8 -*-
"""Правка графа корректора: одна катушка на обмотку и только СОБСТВЕННЫЕ утечки катушек.

Две правки, обе по указанию с эскиза:

  1. Обмотки `c1…c4` стояли в таблице ДВАЖДЫ — по половине МДС на каждом из двух участков корня
     полюса (`T1r1` и `T1r2` и т.д.). Суммарная МДС при этом верна (548 витков), но читается это
     как две катушки, а физически катушка одна. Здесь она сводится на один участок целиком.

     Признак «одна катушка, размазанная по участкам» — доли МДС, дающие в сумме единицу. Если у
     записей доля 1.0 каждая, это РАЗНЫЕ катушки одной обмотки, и сводить их нельзя: у `main`
     таких две, по одной на каждую половину ярма, и это не ошибка.

  2. Корень каждого полюса стоял ДВУМЯ одинаковыми участками (`T1r1` + `T1r2` и так далее) с
     промежуточным узлом `s1*`. Делить его было нужно, только пока по участкам была размазана
     катушка; теперь это лишние рёбра, нагружающие цепь без физического смысла, и они сводятся
     в одну ветвь на полюс. Длина складывается, сечение у участков и так одно.

  3. После этого лишними становятся и узлы степени 2 между двумя железными ветвями: в них ничего
     не ветвится, значит и узла там нет. У боковых полюсов такого не выходит — в узел `m*`
     приходит утечка катушки, — а у среднего полюса наконечник (`T0t`, `B0t`) сводится с корнем.

  4. Из утечек остаются только собственные утечки катушек: воздушный путь в обход ТОЙ ЖЕ ветви, на
     которой сидит катушка, — ни больше, ни меньше. Прежние двенадцать МЕЖПОЛЮСНЫХ рёбер (`lk_*`
     между столбцами) удаляются, и вместе с ними окно `w_cT_cB`: оно шло через всё ярмо, от угла до
     угла, а значит захватывало и средний участок `yokeM` МЕЖДУ катушками `main`, который к их
     собственной утечке отношения не имеет. Вместо него у `main` теперь две утечки, по одной на
     катушку.

Класс новых рёбер — `leak_coil`, и он намеренно НЕ начинается с `gap_`: в базис кика входят только
классы `gap_*` и `tip_row*`, а собственная утечка идёт вдоль ноги полюса, далеко от пучка. Если бы
класс начинался с `gap_`, ребро попало бы в базис — и молча дало бы там тождественный ноль, потому
что у узлов `n*`/`m*` координат нет и `direction()` вернул бы 0j.

Скрипт ИДЕМПОТЕНТЕН: прежние рёбра утечки удаляются перед добавлением.

    py -3.14 corrector_graph_edit.py --leak 0.005
    py -3.14 corrector_graph_edit.py --leak 0          # убрать утечки совсем
"""
from __future__ import annotations

import argparse
import csv
import io
import os

BR = os.path.join("workspace", "tables", "saturation_branches.csv")
MB = os.path.join("workspace", "tables", "saturation_branch_members.csv")
ELEMENT = "corrector"


def merge_coils(path=MB, element=ELEMENT):
    """Свести РАЗМАЗАННУЮ катушку на один участок; обмотки из нескольких катушек оставить как есть.

    Размазанная катушка узнаётся по долям МДС, дающим в сумме единицу (0.5 + 0.5). Доля 1.0 у
    каждой записи означает отдельные катушки одной обмотки — так намотана `main`, по катушке на
    каждую половину ярма."""
    rows = list(csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"))
    head = list(rows[0])
    total = {}
    for r in rows:
        if r["element"] == element:
            total[r["winding"]] = total.get(r["winding"], 0.0) + abs(float(r["mmf_share"] or 1))
    out, seen = [], {}
    for r in rows:
        if r["element"] != element or total[r["winding"]] > 1.0001:
            out.append(r)                   # чужой элемент или обмотка из нескольких катушек
            continue
        key = r["winding"]
        if key in seen:
            continue                        # вторая половина той же катушки — она больше не нужна
        seen[key] = r["branch"]
        out.append(r | {"mmf_share": "1", "comment": "катушка на полюсе, одна на обмотку"})
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, head, delimiter=";")
        w.writeheader()
        w.writerows(out)
    return seen


def merge_pole_roots(path=BR, element=ELEMENT):
    """Свести участки корня каждого полюса в ОДНУ ветвь; промежуточные узлы `s1*` исчезают."""
    rows = list(csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"))
    head = list(rows[0])
    groups = {}
    for r in rows:
        if r["element"] == element and r["sym_class"].startswith("pole_root"):
            groups.setdefault(r["branch"][:2], []).append(r)

    out, done, dropped = [], set(), []
    for r in rows:
        pole = r["branch"][:2]
        if r["element"] != element or not r["sym_class"].startswith("pole_root"):
            out.append(r)
            continue
        if pole in done:
            dropped.append(r["branch"])
            continue
        done.add(pole)
        seg = sorted(groups[pole], key=lambda x: x["branch"])
        out.append(seg[0] | {
            "len_m": f"{sum(float(x['len_m']) for x in seg):.6g}",
            "node_to": seg[-1]["node_to"], "pos_to": seg[-1]["pos_to"],
            "z_to": seg[-1]["z_to"],
            "comment": "корень полюса одной ветвью"})
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, head, delimiter=";")
        w.writeheader()
        w.writerows(out)
    return dropped


def merge_degree2(path=BR, element=ELEMENT):
    """Слить пары железных ветвей, сходящихся в узле степени 2: такой узел ничего не разветвляет.

    Сечения у сливаемых участков разные, поэтому у общей ветви берётся ЭКВИВАЛЕНТНОЕ сечение
    L/Σ(lᵢ/Sᵢ) — то, что даёт ту же магнитную проводимость. В насыщении это приближение: узкий
    наконечник насыщается раньше корня, и одной парой (S, L) эта разница не передаётся."""
    rows = list(csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"))
    head = list(rows[0])
    mine = [r for r in rows if r["element"] == element]
    deg = {}
    for r in mine:
        for n in (r["node_from"], r["node_to"]):
            deg[n] = deg.get(n, 0) + 1

    merged, dropped = {}, []
    for n, d in sorted(deg.items()):
        if d != 2:
            continue
        inc = [r for r in mine if n in (r["node_from"], r["node_to"])
               and r["branch"] not in dropped]
        if len(inc) != 2 or any(r["kind"] != "iron" for r in inc):
            continue
        a_, b_ = sorted(inc, key=lambda r: r["branch"])
        ends = [x for r in (a_, b_) for x in (r["node_from"], r["node_to"]) if x != n]
        L = sum(float(r["len_m"]) for r in (a_, b_))
        S = L / sum(float(r["len_m"]) / float(r["S_m2"]) for r in (a_, b_))
        merged[a_["branch"]] = {"len_m": f"{L:.6g}", "S_m2": f"{S:.6g}",
                                "node_from": ends[0], "node_to": ends[1],
                                "comment": f"{a_['branch']}+{b_['branch']} одной ветвью"}
        dropped.append(b_["branch"])

    out = []
    for r in rows:
        if r["element"] == element and r["branch"] in dropped:
            continue
        out.append(r | merged.get(r["branch"], {}) if r["element"] == element else r)
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, head, delimiter=";")
        w.writeheader()
        w.writerows(out)
    return dropped


def rewrite_leaks(leak, path=BR, members=MB, element=ELEMENT):
    """Убрать прежние утечки и поставить по одной на КАЖДУЮ катушку, в обход её собственной ветви.

    Правило одно на все обмотки: катушка сидит на ветви — её утечка идёт мимо этой ветви, теми же
    концами. Поэтому у `main` их две (по катушке на половину ярма), а средний участок `yokeM` между
    ними ни в одну не входит."""
    rows = list(csv.DictReader(io.open(path, encoding="utf-8-sig"), delimiter=";"))
    head = list(rows[0])
    rows = [r for r in rows
            if not (r["element"] == element and r["branch"].startswith(("lk_", "w_")))]
    by_branch = {r["branch"]: r for r in rows if r["element"] == element}

    coiled = []
    for r in csv.DictReader(io.open(members, encoding="utf-8-sig"), delimiter=";"):
        if r["element"] == element and r.get("branch") in by_branch:
            coiled.append((r["branch"], r["winding"]))

    new = []
    for br, wd in coiled if leak > 0 else []:
        src = by_branch[br]
        new.append({c: "" for c in head} | {
            "element": element, "branch": f"lk_{br}", "material": "air",
            "S_m2": f"{leak:.6g}", "len_m": src["len_m"], "origin": "fit",
            "comment": f"собственная утечка катушки {wd}: воздух в обход ветви {br}",
            "node_from": src["node_from"], "node_to": src["node_to"], "kind": "gap",
            "sym_class": "leak_coil", "mirror": "0",
            "pos_from": src["pos_from"], "pos_to": src["pos_to"],
            "z_from": src["z_from"], "z_to": src["z_to"]})

    i0 = next(k for k, r in enumerate(rows) if r["element"] == element)
    rows = rows[:i0] + new + rows[i0:]
    with io.open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, head, delimiter=";")
        w.writeheader()
        w.writerows(rows)
    return sum(1 for r in rows if r["element"] == element), [r["branch"] for r in new]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--leak", type=float, default=0.005,
                    help="сечение собственной утечки катушки, кв.м")
    a = ap.parse_args(argv)
    seen = merge_coils()
    gone = merge_pole_roots()
    _, made = rewrite_leaks(a.leak)
    plain = merge_degree2()
    n = sum(1 for r in csv.DictReader(io.open(BR, encoding="utf-8-sig"), delimiter=";")
            if r["element"] == ELEMENT)
    print("обмотки сведены на один участок:", ", ".join(f"{w}->{b}" for w, b in sorted(seen.items())))
    print("корни полюсов сведены, убраны ветви:", ", ".join(gone))
    print("узлы степени 2 сведены, убраны ветви:", ", ".join(plain) or "нет")
    print(f"{ELEMENT}: собственные утечки {', '.join(made)}; всего ветвей {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
