# -*- coding: utf-8 -*-
"""Записать warm-start проектными K из .lte в mat_* ВСЕГО кольца (реш. 101: 38 стандартных супериодов +
junction/инжектор). Многообмоточные — коэф. в первую обмотку, остальные 0. mat_B = θ/(QP·Leff).
Бэкап сохраняется. Запуск: py -3.14 warmstart_ring.py [lattice=1]"""
import os, sys, csv, shutil
import matrix_tracker as mt
import lte_warmstart as ws

TB = mt.TB
QP = mt.QP


def write_warmstart(lattice=1, lte_dir=None):
    if lte_dir:
        ws.set_lte_dir(lte_dir)                              # читать .lte из папки конфигурации (без копий в корень)
    mapping = ws.build_mapping_full(lattice, verbose=False)
    grp = ws.grouped_magnets(lattice)
    eff = mt.load_db()[0]
    path = os.path.join(TB, "lattice_configs.csv")
    allrows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(allrows[0].keys())
    byid = {r["instance_id"]: r for r in allrows if r["lattice_id"] == str(lattice)}

    nbend = nfoc = 0
    for i, g in enumerate(grp):
        c = mapping[i]
        dev = g[0]["device_type_id"]; Leff = eff.get(dev, 1.0)
        # ЗНАК: положительный ANGLE из .lte = поворот ВПРАВО, т.е. θ0 = −ANGLE (см. lte_to_lattice.step).
        # Геометрия уложена так же ⇒ кик и ящик гнут в одну сторону. Без минуса диполь гнул бы ПРОТИВ
        # своей же геометрии: до разворота укладки было theta+angle и θ0=+ANGLE — согласовано, теперь нет.
        matB = -c["angle"] / (QP * Leff) if (c["kind"] in ("dipole", "dq") and c["angle"]) else 0.0
        if matB: nbend += 1
        if c["G"] or c["S"] or c["O"]: nfoc += 1
        for k, inst in enumerate(g):
            row = byid[inst["instance_id"]]
            if k == 0:
                row["mat_B"] = f"{matB:.9g}"; row["mat_G"] = f"{c['G']:.9g}"
                row["mat_S"] = f"{c['S']:.9g}"; row["mat_O"] = f"{c['O']:.9g}"
            else:
                row["mat_B"] = row["mat_G"] = row["mat_S"] = row["mat_O"] = "0"

    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(path, os.path.join(TB, "backups", "lattice_configs_pre_ringwarmstart.csv"))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(allrows)
    print(f"warm-start записан в решётку {lattice}: {len(grp)} групп ({nbend} изгибающих, {nfoc} фокусирующих)")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="warm-start проектными K из .lte в mat_* кольца")
    ap.add_argument("lattice", nargs="?", type=int, default=1)
    ap.add_argument("--lte-dir", default=None, help="папка с *_desc.lte/*_line.lte (по умолч. workspace/lte)")
    a = ap.parse_args()
    write_warmstart(a.lattice, a.lte_dir)
