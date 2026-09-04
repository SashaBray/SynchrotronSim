# -*- coding: utf-8 -*-
"""Оставить в БД опыт с траекториями по одной карте поля — чтобы посмотреть его визуализатором.

`saturation_model` гоняет диффур на ВРЕМЕННОМ опыте 900 и стирает его за собой вместе с файлом
`trajectories_exp900.bin`, поэтому посмотреть те траектории постфактум нельзя — их уже нет. Этот
скрипт делает тот же прогон, но НЕ убирает за собой: остаются строки опыта в БД, временное
устройство в `hardware_library` и файл траекторий.

    py -3.14 make_viz_exp.py QP_sh 19000        # верхний ток, где градиент убывает
    py -3.14 live_visualizer.py 900             # смотреть

ПОСЛЕ ПРОСМОТРА уберите за собой — иначе временные строки останутся в БД навсегда:

    py -3.14 make_viz_exp.py --clean

Тот же результат даёт любой следующий прогон `saturation_model.py --element ...`: он вычищает
хвосты временных устройств перед стартом.
"""
from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import element_compare as ec
import saturation_model as sm

# HWCSV в saturation_model — ЛОКАЛЬНАЯ переменная measure_element, не модульная константа.
HWCSV = os.path.join("workspace", "tables", "hardware_library.csv")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", help="элемент, напр. QP_sh")
    ap.add_argument("amp_turns", nargs="?", type=int, help="ампер-витки карты")
    ap.add_argument("--clean", action="store_true", help="убрать временные строки и выйти")
    a = ap.parse_args(argv)

    if a.clean:
        ec.cleanup()
        sm._purge_temp_devices(HWCSV)
        print("временные строки убраны")
        return 0
    if not a.element or not a.amp_turns:
        ap.error("нужны элемент и ампер-витки, либо --clean")

    rows = [r for r in sm.read_registry() if r["element"] == a.element]
    rec = next((r for r in rows
                if r["amp_turns"] and int(float(r["amp_turns"])) == a.amp_turns), None)
    if rec is None:
        have = sorted({int(float(r["amp_turns"])) for r in rows if r["amp_turns"]})
        raise SystemExit(f"нет карты {a.element} при {a.amp_turns} А·в; есть: {have}")

    kind = sm.KIND[a.element]
    radius = sm.PROBE_R.get(kind, 1e-3)
    base = next(d for d in ec.HW if d.startswith(a.element))

    sm._purge_temp_devices(HWCSV)          # хвосты прошлых прогонов
    name, row = sm._temp_device(ec, base, rec)
    hw = ec.mt._read(os.path.basename(HWCSV))
    header = list(hw[0].keys())
    ec.HW[name] = row
    ec._write(os.path.basename(HWCSV), header,
              [[r.get(c, "") for c in header] for r in hw + [row]])
    sm.convert_map(rec)

    # Ток пересчитывается на СУРРОГАТНУЮ энергию: диффур летит на V_BEAM, а поле карты делится
    # на ENERGY_FACTOR, чтобы траектория совпала с 6-ГэВ-ной. См. element_compare.
    arg = float(row["nominal_arg"]) / ec.ENERGY_FACTOR
    trajs, geom = sm.run_diffur_safe(ec, name, arg, radius=radius)
    ok, r_in, rev = sm.screen_trajectories(ec, trajs)

    print(f"  карта {rec['map_id']}  ({a.amp_turns} А·в)")
    print(f"  устройство {name}, опыт {ec.TMP_EXP}, частиц {len(trajs)}")
    print(f"  радиус пробы {radius * 1e3:.1f} мм, |r| внутри магнита {r_in:.2f} мм"
          + ("" if ok else "  ← ОТБРАКОВКА: " + ("разворот" if rev else "вылет за хорошее поле")))
    p = os.path.join("workspace", "results", f"trajectories_exp{ec.TMP_EXP}.bin")
    print(f"  траектории: {p} ({os.path.getsize(p) / 1e6:.1f} МБ)"
          if os.path.exists(p) else "  ФАЙЛ ТРАЕКТОРИЙ НЕ ПОЯВИЛСЯ")
    print(f"\n  смотреть:   py -3.14 live_visualizer.py {ec.TMP_EXP}")
    print(f"  убрать:     py -3.14 make_viz_exp.py --clean")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
