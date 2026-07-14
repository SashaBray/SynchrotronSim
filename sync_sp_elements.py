# -*- coding: utf-8 -*-
"""Синхронизация sp_elements.csv из настроенного lattice_configs.csv (ОБРАЩЕНИЕ assemble).

assemble_lattice_configs.py разворачивает sp_elements (шаблон супериода) в полное кольцо по
lattice_layout: ключ (sp_type, sp_element_id) → (device_type_id, arg_val, correct_coeff). Здесь —
обратно: берём настроенные arg_val/correct_coeff из lattice_configs и переносим в sp_elements, чтобы
live_visualizer стартовал с актуальных (близких к проектным) токов, а не заглушек 1.0.

Токи ОБЩИЕ на все SP ⇒ у всех инстансов одного (sp_type, sp_element_id) значения должны совпадать
(проверяем; при расхождении берём медиану по arg_val и предупреждаем). Запуск: python sync_sp_elements.py [lattice=101]"""
import os
import sys
import csv
import shutil
from collections import defaultdict

TB = os.path.join("workspace", "tables")


def main(lattice=101):
    layout = os.path.join(TB, "lattice_layout.csv")
    lc = os.path.join(TB, "lattice_configs.csv")
    sp = os.path.join(TB, "sp_elements.csv")

    # global_id -> (sp_type, sp_element_id)
    gid2key = {}
    for r in csv.DictReader(open(layout, encoding="utf-8"), delimiter=";"):
        try:
            gid2key[int(r["global_id"])] = (r["sp_type"].strip(), int(r["sp_element_id"]))
        except (KeyError, ValueError):
            continue

    # instance_id -> (arg_val, correct_coeff) из настроенной решётки
    key2vals = defaultdict(list)
    for r in csv.DictReader(open(lc, encoding="utf-8"), delimiter=";"):
        if r["lattice_id"] != str(lattice):
            continue
        try:
            gid = int(r["instance_id"])
        except ValueError:
            continue
        key = gid2key.get(gid)
        if key is None:
            continue
        key2vals[key].append((r.get("arg_val", "").strip(), r.get("correct_coeff", "").strip()))

    # представитель на ключ + проверка согласованности (общие токи)
    key2rep = {}
    inconsistent = 0
    for key, vals in key2vals.items():
        av = [v[0] for v in vals]; cc = [v[1] for v in vals]
        if len(set(av)) > 1 or len(set(cc)) > 1:
            inconsistent += 1
            # берём самый частый (моду)
            av_rep = max(set(av), key=av.count); cc_rep = max(set(cc), key=cc.count)
        else:
            av_rep, cc_rep = av[0], cc[0]
        key2rep[key] = (av_rep, cc_rep)

    # обновить sp_elements
    rows = list(csv.DictReader(open(sp, encoding="utf-8"), delimiter=";"))
    fields = list(rows[0].keys())
    updated = 0; missing = 0
    for r in rows:
        try:
            key = (r["sp_type"].strip(), int(r["sp_element_id"]))
        except (KeyError, ValueError):
            continue
        if key in key2rep:
            av, cc = key2rep[key]
            r["arg_val"] = av if av != "" else "-"
            if "correct_coeff" in r:
                r["correct_coeff"] = cc if cc != "" else "1.0"
            updated += 1
        else:
            missing += 1

    shutil.copy(sp, sp + ".bak")
    with open(sp, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows)

    print(f"sp_elements синхронизирован из lattice {lattice}: обновлено {updated} строк, "
          f"без пары в решётке {missing}. Несогласованных ключей (разные токи у инстансов): {inconsistent}.")
    print(f"Бэкап: {os.path.basename(sp)}.bak")
    # показать DQ
    print("токи DQ в sp_elements теперь:")
    for r in rows:
        if "DQ_I" in r.get("device_type_id", ""):
            print(f"  {r['device_type_id']:16} arg_val={r['arg_val']:>12} correct_coeff={r.get('correct_coeff','')}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
