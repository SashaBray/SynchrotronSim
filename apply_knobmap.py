# -*- coding: utf-8 -*-
"""Применить настроенные C++-оптимизатором ручки (matrix_tune_result.csv) к mat_* решётки в
lattice_configs. Читает knob-map (matrix_knobs.csv: coeff + group_indices) и значения; для каждой
ручки пишет mat_G/S/O в ПЕРВУЮ обмотку каждой управляемой группы. Изгибы (mat_B) не трогаются
(заданы warm-start геометрией). Бэкап сохраняется. Запуск: py -3.14 apply_knobmap.py [lattice=101]"""
import os, sys, csv, shutil
import matrix_tracker as mt
import lte_warmstart as ws

TB = mt.TB
MATCOL = {"G": "mat_G", "S": "mat_S", "O": "mat_O"}


def apply(lattice=101):
    knobrows = list(csv.DictReader(open(os.path.join(TB, "matrix_knobs.csv"), encoding="utf-8"), delimiter=";"))
    tuned = {r["knob_id"]: float(r["value"])
             for r in csv.DictReader(open(os.path.join(TB, "matrix_tune_result.csv"), encoding="utf-8"), delimiter=";")}
    grp = ws.grouped_magnets(lattice)

    path = os.path.join(TB, "lattice_configs.csv")
    allrows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(allrows[0].keys())
    byid = {r["instance_id"]: r for r in allrows if r["lattice_id"] == str(lattice)}

    n = 0
    for kr in knobrows:
        col = MATCOL[kr["coeff"]]; v = tuned[kr["knob_id"]]
        for gi in kr["group_indices"].split("|"):
            first = grp[int(gi)][0]["instance_id"]
            byid[first][col] = f"{v:.9g}"; n += 1

    shutil.copy(path, os.path.join(TB, "backups", "lattice_configs_pre_ringtune.csv"))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(allrows)
    print(f"применено {len(tuned)} ручек к {n} группам решётки {lattice} (mat_G/S/O обновлены)")


if __name__ == "__main__":
    apply(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
