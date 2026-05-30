"""
Сборка lattice_configs.csv из двух исходных таблиц:

  workspace/tables/sp_elements.csv    -- шаблоны супериодов (sp_type, sp_element_id,
                                          device_type_id, arg_val). РЕДАКТИРУЕТСЯ РУКАМИ.
  workspace/tables/lattice_layout.csv -- полная развёртка кольца с координатами
                                          (global_id, sp_number, sp_type, sp_element_id,
                                          x, y, z, yaw, pitch, roll). Генерируется
                                          скриптом lte_to_lattice.py.

Выход:
  workspace/tables/lattice_configs.csv -- ПЕРЕЗАПИСЫВАЕТСЯ! Итоговая таблица для
                                          C++ симулятора, с подставленными токами
                                          из sp_elements.

Логика: для каждой строки lattice_layout ищется пара (sp_type, sp_element_id) в
sp_elements, копируется device_type_id и arg_val. Если пара не найдена -- строка
пропускается с предупреждением.
"""
import csv
import json
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WS = os.path.join(SCRIPT_DIR, "workspace")
SP_ELEMENTS    = os.path.join(WS, "tables", "sp_elements.csv")
LATTICE_LAYOUT = os.path.join(WS, "tables", "lattice_layout.csv")
LATTICE_CONFIGS = os.path.join(WS, "tables", "lattice_configs.csv")
CONFIG_PATH    = os.path.join(WS, "configs", "lte_mapping.json")

def _load_lattice_id():
    if not os.path.exists(CONFIG_PATH):
        return 101  # fallback для совместимости
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return int(json.load(f).get("lattice_id", 101))

LATTICE_ID = _load_lattice_id()

# Sub-lattices загружаются из workspace/configs/sub_lattices.json.
# Это укороченные копии lattice 101 для экспериментов на части кольца.
SUB_LATTICES_CFG = os.path.join(WS, "configs", "sub_lattices.json")

def _load_sub_lattices():
    if not os.path.exists(SUB_LATTICES_CFG):
        print(f"  (no {SUB_LATTICES_CFG}, skipping sub-lattices)")
        return []
    try:
        with open(SUB_LATTICES_CFG, encoding="utf-8") as f:
            cfg = json.load(f)
        out = []
        for entry in cfg.get("sub_lattices", []):
            out.append((int(entry["lattice_id"]), int(entry["max_gid"])))
        return out
    except Exception as e:
        print(f"  WARN: cannot read {SUB_LATTICES_CFG}: {e}")
        return []

SUB_LATTICES = _load_sub_lattices()


def load_sp_elements(path):
    """sp_elements.csv -> {(sp_type, sp_element_id): (device_type_id, arg_val, correct_coeff)}.

    Колонка correct_coeff -- per-instance множитель поля (опциональная для
    обратной совместимости; default 1.0). Эффективный коэффициент в симуляторе
    = hardware_library.correct_coeff * sp_elements.correct_coeff.
    """
    if not os.path.exists(path):
        print(f"Error: {path} not found. Run lte_to_lattice.py first.")
        sys.exit(1)
    sp_dict = {}
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            try:
                sp_type  = row["sp_type"].strip()
                sp_eid   = int(row["sp_element_id"])
                dev_type = row["device_type_id"].strip()
                arg_val  = row["arg_val"].strip()
                coeff    = row.get("correct_coeff", "1.0")
                coeff    = coeff.strip() if coeff else "1.0"
                if not coeff:
                    coeff = "1.0"
            except (KeyError, ValueError):
                continue
            sp_dict[(sp_type, sp_eid)] = (dev_type, arg_val, coeff)
    return sp_dict


def assemble():
    sp_dict = load_sp_elements(SP_ELEMENTS)
    print(f"Loaded sp_elements: {len(sp_dict)} entries, "
          f"{len(set(k[0] for k in sp_dict))} SP type(s)")

    if not os.path.exists(LATTICE_LAYOUT):
        print(f"Error: {LATTICE_LAYOUT} not found. Run lte_to_lattice.py first.")
        sys.exit(1)

    output_rows = []
    missing = set()
    with open(LATTICE_LAYOUT, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            try:
                sp_type = row["sp_type"].strip()
                sp_eid = int(row["sp_element_id"])
                global_id = int(row["global_id"])
            except (KeyError, ValueError):
                continue
            key = (sp_type, sp_eid)
            if key not in sp_dict:
                missing.add(key)
                continue
            dev_type, arg_val, coeff = sp_dict[key]
            output_rows.append({
                "lattice_id":     LATTICE_ID,
                "instance_id":    global_id,
                "device_type_id": dev_type,
                "x":              row["x"],
                "y":              row["y"],
                "z":              row["z"],
                "yaw":            row["yaw"],
                "pitch":          row["pitch"],
                "roll":           row["roll"],
                "arg_val":        arg_val,
                "correct_coeff":  coeff,
            })

    if missing:
        sample = ", ".join(f"({t}, {i})" for t, i in sorted(missing)[:5])
        print(f"WARN: {len(missing)} (sp_type, sp_element_id) keys present in layout "
              f"but missing from sp_elements: {sample}"
              + (", ..." if len(missing) > 5 else ""))

    # Авто-генерация суб-lattice'ов: копии первых N строк lattice 101 с
    # изменённым lattice_id. Нужно, чтобы эксперименты на укороченных lattice
    # (102/103/104) могли работать с актуальными значениями токов.
    sub_rows = []
    for sub_id, max_gid in SUB_LATTICES:
        for r in output_rows:
            if int(r["instance_id"]) <= max_gid:
                sub_rows.append({**r, "lattice_id": sub_id})
    output_rows.extend(sub_rows)

    # Запись с retry: Windows может временно заблокировать файл (Excel/антивирус)
    import time as _time
    last_err = None
    for k in range(12):
        try:
            with open(LATTICE_CONFIGS, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f, delimiter=";", lineterminator="\n")
                w.writerow(["lattice_id", "instance_id", "device_type_id",
                            "x", "y", "z", "yaw", "pitch", "roll",
                            "arg_val", "correct_coeff"])
                for r in output_rows:
                    w.writerow([r["lattice_id"], r["instance_id"], r["device_type_id"],
                                r["x"], r["y"], r["z"],
                                r["yaw"], r["pitch"], r["roll"],
                                r["arg_val"], r["correct_coeff"]])
            break
        except PermissionError as e:
            last_err = e
            if k == 0:
                print(f"  [retry] {os.path.basename(LATTICE_CONFIGS)} locked, waiting...")
            _time.sleep(0.5)
    else:
        raise last_err
    main_count = sum(1 for r in output_rows if int(r["lattice_id"]) == LATTICE_ID)
    print(f"Wrote {LATTICE_CONFIGS}: {len(output_rows)} rows "
          f"(lattice {LATTICE_ID}: {main_count}, sub-lattices: {len(output_rows) - main_count})")


if __name__ == "__main__":
    assemble()
