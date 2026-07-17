"""
Генератор БД сегментов опорной (design) орбиты из .lte.

Что делает:
    1. Парсит workspace/lte/SYLA_desc.lte и SYLA_line.lte (через функции из
       lte_to_lattice.py).
    2. Применяет «замыкание кольца» — масштабирование ANGLE чистых диполей,
       чтобы суммарный bend стал ровно 2π. ТА ЖЕ логика, что в lte_to_lattice.
    3. Делит line.lte на супериоды (по sp_start_pattern из lte_mapping.json) и
       присваивает каждому sp_type (по первому магниту, через
       sp_type_by_first_element).
    4. Для ПЕРВОГО встреченного экземпляра каждого sp_type выписывает упорядоченный
       список геометрических сегментов:
         - kind=line: EDRIFT / KQUAD / KSEXT / KOCT / BPM (точка с L=0) / CSBEND с ANGLE=0
         - kind=arc:  CSBEND с ANGLE≠0 (DL/JL/SBM/DQ)
    5. Пишет workspace/tables/design_orbit_segments.csv:
         sp_type;seg_idx;kind;length;angle

Это БД-источник; для конкретного lattice_id орбита собирается из неё
скриптом design_orbit.py (модуль build_orbit_for_lattice).

Запуск:
    python lte_to_design_orbit.py

Перезапускайте, если поменялся .lte или lte_mapping.json. Файл-результат
коммитится в репозиторий и используется visualize.py / optimizer.cpp.
"""
from __future__ import annotations
import csv
import os
import math
import sys
from pathlib import Path

from _common import project_root, workspace_dir
import lte_to_lattice as lte  # переиспользуем parse_desc, parse_line, is_sp_start

EPS_ANGLE = 1e-12   # ниже этого CSBEND считается прямым


def main() -> None:
    root = project_root()
    ws   = workspace_dir(root)
    desc_path = ws / "lte" / "SYLA_desc.lte"
    line_path = ws / "lte" / "SYLA_line.lte"
    out_path  = ws / "tables" / "design_orbit_segments.csv"
    lattice_id = int(getattr(lte, "_CFG", {}).get("lattice_id", 1))   # опорная орбита ПРИВЯЗАНА к установке

    print(f"Reading {desc_path}")
    devices = lte.parse_desc(str(desc_path))
    print(f"  parsed {len(devices)} device definitions")

    print(f"Reading {line_path}")
    line = lte.parse_line(str(line_path))
    print(f"  parsed line with {len(line)} element references")

    # --- Замыкание кольца УДАЛЕНО (был скрытый closure-хак) ---------------
    # Раньше здесь ANGLE чистых диполей масштабировался до Σ=2π. Теперь замыкание чинится
    # В ИСТОЧНИКЕ (.lte) через lte_ring.py; здесь — только предупреждение, чтобы
    # незамкнутость была ВИДНА, а не молча размазана в опорную орбиту.
    raw_pure_bend = raw_dq_bend = 0.0
    for name in line:
        d = devices.get(name)
        if not d:
            continue
        angle = d.get("ANGLE", 0.0)
        if not isinstance(angle, (int, float)) or angle == 0.0:
            continue
        k1 = d.get("K1", 0.0)
        if isinstance(k1, (int, float)) and k1 != 0.0:
            raw_dq_bend += angle
        else:
            raw_pure_bend += angle
    raw_total = raw_pure_bend + raw_dq_bend
    deficit = 2.0 * math.pi - raw_total
    if abs(deficit) > 1e-6:
        print(f"\n  [ВНИМАНИЕ] Кольцо НЕ замкнуто: Σ ANGLE = {math.degrees(raw_total):.4f}° "
              f"(недобор {deficit * 1e3:+.2f} мрад). Опорная орбита будет ОТКРЫТА. "
              f"Замкнуть: lte_ring.py <файл> --out <файл>")

    # --- Проход по line.lte, эмит сегментов по sp_type --------------------
    # Для КАЖДОГО уникального sp_type сохраняем сегменты ПЕРВОГО встреченного
    # экземпляра. Все экземпляры одного sp_type считаем геометрически идентичными.
    #
    # Тонкость: SP начинается с дрифтов (DR_01_2P651, BPM, DR_02_001), а
    # sp_type определяется по первому МАГНИТУ. Поэтому буферизуем сегменты в
    # `pending` до момента определения sp_type, затем либо сливаем буфер в БД
    # (если sp_type новый), либо отбрасываем (если SP такого типа уже сохранён).
    # Сегмент: (kind, length, angle, lte_name) — lte_name нужен для связки
    # с lattice_layout через map_to_hw (некоторые LTE-имена эмитят несколько
    # строк в lattice, например KSEXT → SP_I1 + SP_I2).
    sp_segments: dict[str, list[tuple[str, float, float, str]]] = {}
    sp_seen_count: dict[str, int] = {}
    current_sp_type: str | None = None
    collecting_active: bool = False  # пишем напрямую в sp_segments[current_sp_type]
    collecting_pending: bool = False  # пишем в буфер pending до определения sp_type
    pending: list[tuple[str, float, float, str]] = []
    sp_number = 0

    def seg_for(name: str, dev) -> tuple[str, float, float, str]:
        L     = float(dev.get("L", 0.0))
        ANGLE = float(dev.get("ANGLE", 0.0))
        dtype = dev.get("type", "")
        if dtype == "CSBEND" and abs(ANGLE) > EPS_ANGLE:
            return ("arc", L, ANGLE, name)
        return ("line", L, 0.0, name)

    for i, name in enumerate(line):
        if name not in devices:
            continue
        dev = devices[name]

        # Граница нового SP -> начинаем буферизацию
        if lte.is_sp_start(line, i):
            sp_number += 1
            current_sp_type = None
            collecting_active = False
            collecting_pending = True
            pending = []

        # Определение sp_type по первому магниту в этом SP
        if current_sp_type is None and lte.map_to_hw(dev):
            current_sp_type = lte.FIRST_ELEMENT_TO_SP_TYPE.get(name, name)
            if current_sp_type not in sp_segments:
                # Новый sp_type: сливаем буфер и продолжаем писать напрямую
                sp_segments[current_sp_type] = list(pending)
                sp_seen_count[current_sp_type] = 1
                collecting_active = True
            else:
                # Этот sp_type уже сохранён — пропускаем остаток текущего SP
                sp_seen_count[current_sp_type] += 1
                collecting_active = False
            collecting_pending = False
            pending = []

        # Эмитим сегмент в соответствующий буфер/таблицу
        seg = seg_for(name, dev)
        if collecting_active:
            sp_segments[current_sp_type].append(seg)
        elif collecting_pending:
            pending.append(seg)
        # иначе элемент находится в SP, чей sp_type уже сохранён — игнор

    # --- Сводка ---------------------------------------------------------
    print(f"\nSuperperiods seen: {sp_number}, "
          f"distinct sp_type(s): {len(sp_segments)}")
    for t, segs in sp_segments.items():
        total_L     = sum(s[1] for s in segs)
        total_angle = sum(s[2] for s in segs)
        n_arc       = sum(1 for s in segs if s[0] == "arc")
        n_line      = len(segs) - n_arc
        print(f"  sp_type='{t}':  {len(segs)} segments "
              f"({n_line} line + {n_arc} arc),  "
              f"sum_L = {total_L:.4f} m,  sum_ANGLE = {math.degrees(total_angle):.4f} deg,  "
              f"instances in ring = {sp_seen_count[t]}")

    # --- Запись CSV ------------------------------------------------------
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator="\n")
        w.writerow(["lattice_id", "sp_type", "seg_idx", "kind", "length", "angle", "lte_name"])
        for sp_type, segs in sp_segments.items():
            for idx, (kind, L, A, lte_name) in enumerate(segs):
                w.writerow([lattice_id, sp_type, idx, kind, f"{L:.10f}", f"{A:.10f}", lte_name])
    print(f"\nWrote {out_path}")


if __name__ == "__main__":
    main()
