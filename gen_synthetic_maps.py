# -*- coding: utf-8 -*-
"""Генерация синтетических КАРТ ПОЛЯ (и STL) для уже созданных синтетических элементов библиотеки.

Синтетические элементы (synthetic=True) добавляются импортёром LTE под нужные эффективные длины, но
по умолчанию идут БЕЗ карты (map_bin_path='-') — матричная модель работает по L_eff. Для диффурного
ядра (Boris по картам) карта нужна: этот драйвер строит её ГИБРИДНЫМ ресайзом карты эталона
(scaled_from) под целевую L_eff (тело ×c_body, фриндж сохранён — из исследования масштабирования карт).

Идемпотентно: пропускает элементы, у которых карта уже есть. Флаг --device X — только один элемент
(для проверки). Обновляет hardware_library.csv (map_bin_path и, если построилась, stl_path).
Синтетические .bin помечены путём synthetic/<id>.bin (карта), STL — synthetic/<id>.STL."""
import os
import sys
import csv
import synthetic_assets as sa

TB = os.path.join("workspace", "tables")
HW = os.path.join(TB, "hardware_library.csv")


def main(only=None, do_stl=True):
    rows = list(csv.reader(open(HW, encoding="utf-8"), delimiter=";"))
    head, data = rows[0], rows[1:]
    col = {c: k for k, c in enumerate(head)}
    hw = {r[0]: r for r in data}
    for req in ("synthetic", "scaled_from", "map_bin_path", "stl_path", "eff_length"):
        if req not in col:
            raise RuntimeError(f"hardware_library.csv нет столбца {req}")

    todo = [r for r in data if r[col["synthetic"]] == "True"
            and r[col["map_bin_path"]] in ("-", "")
            and (only is None or r[0] == only)]
    if not todo:
        print("Нечего генерировать (все синтетики уже с картами или фильтр пуст).")
        return

    print(f"К генерации: {len(todo)} синтетич. элементов.")
    made = 0
    for r in todo:
        nid = r[0]; ref = r[col["scaled_from"]]; new_eff = float(r[col["eff_length"]])
        if ref not in hw:
            print(f"  [skip] {nid}: эталон {ref} не в библиотеке"); continue
        ref_row = hw[ref]
        binp = sa.gen_field_map(ref, nid, new_eff, ref_row, col)
        if binp:
            r[col["map_bin_path"]] = binp
            sz = os.path.getsize(os.path.join("workspace", "fields_bin", binp))
            print(f"  + {nid}: {binp}  ({sz/1e6:.1f} МБ)  <- {ref}, L_eff={new_eff*1e3:.1f}мм")
            made += 1
        else:
            print(f"  [FAIL] {nid}: карта не построилась (эталон {ref})")
        if do_stl and r[col["stl_path"]] in ("-", ""):
            stlp = sa.gen_stl(ref, nid, new_eff, ref_row, col)
            if stlp:
                r[col["stl_path"]] = stlp

    with open(HW, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(head); w.writerows(data)
    print(f"\n[OK] построено {made}/{len(todo)} карт; hardware_library.csv обновлён.")


if __name__ == "__main__":
    only = None; do_stl = True
    for a in sys.argv[1:]:
        if a == "--no-stl":
            do_stl = False
        elif a.startswith("--device="):
            only = a.split("=", 1)[1]
        elif a == "--device" and False:
            pass
    if "--device" in sys.argv:
        i = sys.argv.index("--device"); only = sys.argv[i + 1]
    main(only, do_stl)
