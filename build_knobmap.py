# -*- coding: utf-8 -*-
"""Сгенерировать knob-map для C++-оптимизатора: какие сгруппированные магниты кольца контролирует
каждая ручка настройки. ПРАВИЛО: все стандартные супериоды одинаковы (общие ручки по позиции в SP);
junction/инжектор — уникальные ручки на каждый магнит. Корректоры/диполи заморожены (без ручек).

Индексы групп совпадают с порядком buildLine в C++ (группировка со-локованных по instance_id).
Формат matrix_knobs.csv:  knob_id;coeff;init;n_groups;group_indices(comma)
coeff ∈ G|S|O ; init — проектное значение из warm-start.

Запуск: py -3.14 build_knobmap.py [lattice=101]
"""
import os, sys, csv
import matrix_tracker as mt
import lte_warmstart as ws

TB = mt.TB
COEFF_OF_KIND = {"quad": ["G"], "sext": ["S"], "oct": ["O"], "dq": ["G", "S"]}


def _load_tune():
    """Настроенный СТАНДАРТНЫЙ супериод (matrix_tune.csv) по lte_name → {lte: {G,S,O}}. Для старта
    стандартных ручек с УЖЕ настроенного 10-SP тюна (проектный K вертикально нестабилен на 40 SP)."""
    p = os.path.join(TB, "matrix_tune.csv")
    if not os.path.exists(p):
        return {}
    out = {}
    for r in csv.DictReader(open(p, encoding="utf-8"), delimiter=";"):
        if r["lte_name"]:
            out[r["lte_name"]] = {"G": float(r["G_Tlm"]), "S": float(r["S_Tlm2"]), "O": float(r["O_Tlm3"])}
    return out


def build_knobmap(lattice=101):
    grp = ws.grouped_magnets(lattice)
    mapping = ws.build_mapping_full(lattice, verbose=False)
    tune = _load_tune()
    # инстанс → (sp_type, sp_element_id) из lattice_layout (global_id == instance_id решётки 101)
    layout = {r["global_id"]: (r["sp_type"], r["sp_element_id"])
              for r in csv.DictReader(open(os.path.join(TB, "lattice_layout.csv"), encoding="utf-8"), delimiter=";")}

    # knob key → dict(coeff, init, indices)
    knobs = {}
    order = []
    for i, g in enumerate(grp):
        m = mapping[i]; kind = m["kind"]
        if kind not in COEFF_OF_KIND:
            continue
        iid = g[0]["instance_id"]
        sp_type, sp_elem = layout.get(iid, ("standard", iid))
        base = ("std", sp_elem) if sp_type == "standard" else ("junc", iid)   # общий/уникальный
        for c in COEFF_OF_KIND[kind]:
            key = (base, c)
            if key not in knobs:
                # стандартные — с настроенного тюна (по lte_name), junction — проектный K
                init = tune.get(m["lte"], {}).get(c, m[c]) if sp_type == "standard" else m[c]
                knobs[key] = dict(coeff=c, init=init, indices=[], sp_type=sp_type)
                order.append(key)
            knobs[key]["indices"].append(i)

    # запись
    path = os.path.join(TB, "matrix_knobs.csv")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["knob_id", "coeff", "init", "sp_type", "n_groups", "group_indices"])
        for kid, key in enumerate(order):
            k = knobs[key]
            # разделитель индексов — '|', НЕ ',' (C++ readCSV конвертирует ',' → '.' для EU-десятичных)
            w.writerow([kid, k["coeff"], f"{k['init']:.9g}", k["sp_type"],
                        len(k["indices"]), "|".join(str(x) for x in k["indices"])])

    nstd = sum(1 for key in order if key[0][0] == "std")
    njunc = len(order) - nstd
    ng_std = sum(len(knobs[key]["indices"]) for key in order if key[0][0] == "std")
    ng_junc = sum(len(knobs[key]["indices"]) for key in order if key[0][0] == "junc")
    print(f"knob-map → {os.path.relpath(path)}: {len(order)} ручек "
          f"({nstd} общих standard управляют {ng_std} магнитами, {njunc} junction — {ng_junc} магнитами)")
    return order, knobs


if __name__ == "__main__":
    build_knobmap(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
