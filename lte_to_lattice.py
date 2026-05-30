"""
Конвертация Elegant LTE -> две таблицы для удобного ручного редактирования.

Вход:
  workspace/lte/SYLA_desc.lte    -- описание устройств LTE
  workspace/lte/SYLA_line.lte    -- последовательность элементов (плоский список через запятую)

Выход:
  workspace/lte/devices.csv             -- таблица всех LTE-устройств с параметрами
  workspace/lte/groups.csv              -- карта LTE-имя -> hardware_library типы
  workspace/tables/sp_elements.csv      -- ШАБЛОНЫ супериодов (sp_type, sp_element_id,
                                           device_type_id, arg_val) -- РЕДАКТИРУЕТСЯ РУКАМИ
  workspace/tables/lattice_layout.csv   -- полная развёртка кольца с координатами,
                                           без токов (global_id, sp_number, sp_type,
                                           sp_element_id, x, y, z, yaw, pitch, roll)

Чтобы получить итоговый lattice_configs.csv -- запустите assemble_lattice_configs.py.

При повторной генерации sp_elements.csv существующие значения arg_val сохраняются
(идентификация по (sp_type, sp_element_id, device_type_id)), чтобы ручные правки
токов не терялись.

Правила соответствия (LTE -> hardware_library.csv):
  EDRIFT, MONI                -- пропуск (промежуток / диагностика)
  KQUAD, L >= 0.295           -> QP_L_1A
  KQUAD, L <  0.295           -> QP_sh_1A
  KSEXT, K2 == 0              -> corrector_I1_1A + corrector_I2_1A (стек)
  KSEXT, K2 != 0              -> SP_I1_1A + SP_I2_1A (стек)
  KOCT                        -> OP_1A
  CSBEND, K1 != 0             -> DQ_I1_1A + DQ_I2_1A + DQ_I3_1A (стек)
  CSBEND, K1 == 0, L <  0.1   -> DP_sh_PM
  CSBEND, K1 == 0, L >= 0.1   -> DP_L_PM (по слайсам семейств DL/JL)
"""
import os
import re
import csv
import json
import math
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
WS = os.path.join(SCRIPT_DIR, "workspace")

LTE_DESC = os.path.join(WS, "lte", "SYLA_desc.lte")
LTE_LINE = os.path.join(WS, "lte", "SYLA_line.lte")

OUT_DEVICES        = os.path.join(WS, "lte",    "devices.csv")
OUT_GROUPS         = os.path.join(WS, "lte",    "groups.csv")
OUT_SP_ELEMENTS    = os.path.join(WS, "tables", "sp_elements.csv")
OUT_LATTICE_LAYOUT = os.path.join(WS, "tables", "lattice_layout.csv")

# ----------------------------------------------------------------------------
# Конфиг -- все маппинги, пороги и паттерны парсинга вынесены в workspace/configs
# ----------------------------------------------------------------------------
CONFIG_PATH = os.path.join(WS, "configs", "lte_mapping.json")

def _load_config():
    if not os.path.exists(CONFIG_PATH):
        raise FileNotFoundError(
            f"LTE mapping config not found at {CONFIG_PATH}. "
            "Этот файл описывает, как LTE-элементы отображаются на типы hardware_library "
            "(см. README.md). Создайте его перед запуском.")
    with open(CONFIG_PATH, encoding="utf-8") as f:
        return json.load(f)

_CFG               = _load_config()
LATTICE_ID         = _CFG["lattice_id"]
FIRST_ELEMENT_TO_SP_TYPE = _CFG["sp_type_by_first_element"]
SP_START_PATTERN   = _CFG["sp_start_pattern"]["names"]
DLJL_FAMILY_RE     = re.compile(_CFG["dljl_family_regex"]["pattern"])
EPS                = float(_CFG["thresholds"]["epsilon"])
_KQUAD_LONG_MIN_L  = float(_CFG["thresholds"]["kquad_long_min_L"])
_CSBEND_SHORT_MAX_L = float(_CFG["thresholds"]["csbend_short_max_L"])
_HW_MAP            = _CFG["lte_to_hw_map"]

def _hw_entries(key):
    """Возвращает [(device_type_id, default_arg_val), ...] из конфига."""
    return [(h["id"], h["default"]) for h in _HW_MAP[key]["hw"]]


# ---------------------------------------------------------------------------
# Этап 1. Парсинг описания элементов
# ---------------------------------------------------------------------------

# Регэксп для пар key=value: ключ слово, значение -- знаковое число с дробью/экспонентой.
_KV = re.compile(r"([A-Za-z_]\w*)\s*=\s*(-?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)")

def parse_desc(path):
    """desc.lte -> {name: {'name','type', 'L', 'ANGLE', 'K1', 'K2', ...}}."""
    devices = {}
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.split("!", 1)[0].strip()
            if not line:
                continue
            m = re.match(r"([A-Za-z_][\w\-]*)\s*:\s*(\w+)\s*,?\s*(.*)", line)
            if not m:
                continue
            name, dtype, params = m.group(1), m.group(2), m.group(3)
            dev = {"name": name, "type": dtype}
            for key, val in _KV.findall(params):
                dev[key.upper()] = float(val)
            devices[name] = dev
    return devices


# ---------------------------------------------------------------------------
# Этап 2. Парсинг последовательности элементов
# ---------------------------------------------------------------------------

def parse_line(path):
    """line.lte -> [name1, name2, ...] (плоский список, разрывы строк и пробелы игнорируются)."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    text = re.sub(r"!.*", "", text)          # убрать комментарии
    names = [n.strip() for n in re.split(r"[,\s]+", text) if n.strip()]
    return names


# ---------------------------------------------------------------------------
# Этап 3. Сопоставление одного LTE-элемента с hardware_library
# ---------------------------------------------------------------------------

def map_to_hw(dev):
    """Возвращает список (hw_device_type_id, arg_val) для одного LTE-устройства.
    Пустой список -- элемент пропускается (промежуток / BPM / неизвестный).

    Правила берутся из workspace/configs/lte_mapping.json (раздел lte_to_hw_map).
    Условия (пороги по L, знак K1/K2) -- из раздела thresholds того же файла.
    """
    t = dev["type"]
    L = dev.get("L", 0.0)

    if t in _HW_MAP["skip_types"]:
        return []

    if t == "KQUAD":
        return _hw_entries("KQUAD_long" if L >= _KQUAD_LONG_MIN_L else "KQUAD_short")

    if t == "KSEXT":
        K2 = dev.get("K2", 0.0)
        return _hw_entries("KSEXT_zero_K2" if abs(K2) < EPS else "KSEXT_nonzero")

    if t == "KOCT":
        return _hw_entries("KOCT")

    if t == "CSBEND":
        K1 = dev.get("K1", 0.0)
        if abs(K1) > EPS:
            return _hw_entries("CSBEND_combined")
        if L < _CSBEND_SHORT_MAX_L:
            return _hw_entries("CSBEND_short")
        # Длинные постоянные диполи (DL, JL семейства). Все слайсы используют
        # одну field-map (DP_L_PM.bin) и один device_type. Per-instance
        # подстройка поля -- через колонку correct_coeff в sp_elements.csv.
        return _hw_entries("CSBEND_long")

    return []  # неизвестный тип


# ---------------------------------------------------------------------------
# Этап 4. Геометрия: один шаг вдоль траектории
# ---------------------------------------------------------------------------

def step(x, y, theta, L, angle):
    """Один проход через элемент.
    Возвращает (центр_x, центр_y, yaw_центра, новые x, y, theta после элемента).

    Точная геометрия дуги: радиус R = L/angle, дуга в локальной системе
    (вход в начале, касательная вдоль локального x).
    """
    if abs(angle) > 1e-12 and L > 0:
        R = L / angle
        s_half = R * math.sin(angle / 2.0)
        t_half = R * (1.0 - math.cos(angle / 2.0))
        s_full = R * math.sin(angle)
        t_full = R * (1.0 - math.cos(angle))
        ct, st = math.cos(theta), math.sin(theta)
        xc     = x + s_half * ct - t_half * st
        yc     = y + s_half * st + t_half * ct
        yaw_c  = theta + angle / 2.0
        x_new  = x + s_full * ct - t_full * st
        y_new  = y + s_full * st + t_full * ct
        th_new = theta + angle
    else:
        xc     = x + (L / 2.0) * math.cos(theta)
        yc     = y + (L / 2.0) * math.sin(theta)
        yaw_c  = theta
        x_new  = x + L * math.cos(theta)
        y_new  = y + L * math.sin(theta)
        th_new = theta
    return xc, yc, yaw_c, x_new, y_new, th_new


# ---------------------------------------------------------------------------
# Выходы
# ---------------------------------------------------------------------------

def write_devices_csv(devices, path):
    cols = ["name", "type", "L", "ANGLE", "K1", "K2", "K3", "E1", "E2", "N_SLICES"]
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cols)
        for name in sorted(devices.keys()):
            d = devices[name]
            row = [d.get(c, "") if c not in ("name", "type") else d[c.lower()]
                   if c.lower() in d else d.get(c, "") for c in cols]
            # Аккуратно: для имён/типа берём из ключей name/type
            row[0] = d["name"]
            row[1] = d["type"]
            w.writerow(row)


def write_groups_csv(devices, path):
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["lte_name", "lte_type", "L", "K_or_angle",
                    "hw_types", "n_entries", "comment"])
        for name in sorted(devices.keys()):
            d = devices[name]
            hw = map_to_hw(d)
            # Краткая колонка K_or_angle для удобства проверки
            key = ""
            if d["type"] == "KQUAD":
                key = f"K1={d.get('K1', 0)}"
            elif d["type"] == "KSEXT":
                key = f"K2={d.get('K2', 0)}"
            elif d["type"] == "KOCT":
                key = f"K3={d.get('K3', 0)}"
            elif d["type"] == "CSBEND":
                key = f"angle={d.get('ANGLE', 0)}, K1={d.get('K1', 0)}"

            if not hw:
                w.writerow([name, d["type"], d.get("L", ""), key, "(skip)", 0, ""])
            else:
                hw_str = "|".join(h[0] for h in hw)
                w.writerow([name, d["type"], d.get("L", ""), key, hw_str, len(hw), ""])


def load_existing_sp_elements(path):
    """Читает существующий sp_elements.csv (если есть) и возвращает словарь
    {(sp_type, sp_element_id, device_type_id): (arg_val, correct_coeff)}.
    Используется, чтобы при перегенерации сохранить уже отредактированные руками
    значения. Если в файле нет колонки correct_coeff -- ставим "1.0"."""
    existing = {}
    if not os.path.exists(path):
        return existing
    try:
        with open(path, encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=";")
            for row in reader:
                try:
                    key = (row["sp_type"].strip(),
                           int(row["sp_element_id"]),
                           row["device_type_id"].strip())
                    arg_val = row["arg_val"].strip()
                    coeff   = row.get("correct_coeff", "1.0").strip() or "1.0"
                    existing[key] = (arg_val, coeff)
                except (KeyError, ValueError):
                    continue
    except Exception as e:
        print(f"  (warning: couldn't load existing {path}: {e})")
    return existing


def write_sp_elements_csv(sp_templates, path):
    """sp_templates: dict {sp_type: [list of {sp_element_id, device_type_id, arg_val, correct_coeff}]}.
    Колонка correct_coeff -- per-instance множитель поля, который домножается
    с device-level correct_coeff из hardware_library.csv при загрузке симулятором."""
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator="\n")
        w.writerow(["sp_type", "sp_element_id", "device_type_id", "arg_val", "correct_coeff"])
        for sp_type in sp_templates:
            for entry in sp_templates[sp_type]:
                w.writerow([
                    sp_type,
                    entry["sp_element_id"],
                    entry["device_type_id"],
                    entry["arg_val"],
                    entry.get("correct_coeff", "1.0"),
                ])


def write_lattice_layout_csv(rows, path):
    """rows: list of {global_id, sp_number, sp_type, sp_element_id, x, y, z, yaw, pitch, roll}."""
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator="\n")
        w.writerow(["global_id", "sp_number", "sp_type", "sp_element_id",
                    "x", "y", "z", "yaw", "pitch", "roll"])
        for r in rows:
            w.writerow([
                r["global_id"],
                r["sp_number"],
                r["sp_type"],
                r["sp_element_id"],
                f"{r['x']:.6f}",
                f"{r['y']:.6f}",
                f"{r['z']:.6f}",
                f"{r['yaw']:.8f}",
                f"{r['pitch']:.6f}",
                f"{r['roll']:.6f}",
            ])


def straighten_dljl_families(layout_rows, devices):
    """Пост-обработка: для каждой группы 5 последовательных слайсов одного DL/JL
    семейства (_1.._5) перепозиционирует элементы так, чтобы они лежали НА ОДНОЙ
    ПРЯМОЙ (касательная в положении 3-го слайса), вместо арочного расположения.

    Центральный (3-й) слайс остаётся на дуге (на опорной траектории).
    Слайсы 1, 2, 4, 5 размещаются по касательной к ней, с расстоянием по длине
    самих магнитов: L_i/2 + L_(i+1)/2 между соседями.

    yaw всех 5 слайсов становится yaw'ом 3-го (центрального).
    Возвращает количество модифицированных групп.
    """
    groups_modified = 0
    i = 0
    while i + 4 < len(layout_rows):
        nm0 = layout_rows[i].get("_lte_name", "")
        m0 = DLJL_FAMILY_RE.match(nm0)
        if not m0:
            i += 1
            continue
        base = m0.group(1)
        # Собираем slice-индексы 5 последовательных строк.
        # Они могут идти в прямом (_1.._5) или обратном (_5.._1) порядке -- LTE
        # допускает оба.
        slice_idx = [int(m0.group(2))]
        ok = True
        for k in range(1, 5):
            nm_k = layout_rows[i + k].get("_lte_name", "")
            m_k = DLJL_FAMILY_RE.match(nm_k)
            if not m_k or m_k.group(1) != base:
                ok = False
                break
            slice_idx.append(int(m_k.group(2)))
        if not ok or sorted(slice_idx) != [1, 2, 3, 4, 5]:
            i += 1
            continue

        # Центральный = слайс _3 (физический центр), независимо от порядка в layout.
        central_pos = slice_idx.index(3)        # 0..4 -- позиция _3 в layout
        central = layout_rows[i + central_pos]
        cx, cy, cyaw = central["x"], central["y"], central["yaw"]
        tx, ty = math.cos(cyaw), math.sin(cyaw)

        # Длины слайсов В ПОРЯДКЕ их следования в layout.
        L = [float(devices[layout_rows[i + k]["_lte_name"]].get("L", 0.0))
             for k in range(5)]

        # Смещения центров от центра central_pos (signed, по касательной).
        d = [0.0] * 5
        # Вперёд от центра (k > central_pos в layout)
        for k in range(central_pos + 1, 5):
            d[k] = d[k - 1] + L[k - 1] / 2.0 + L[k] / 2.0
        # Назад от центра (k < central_pos)
        for k in range(central_pos - 1, -1, -1):
            d[k] = d[k + 1] - (L[k] / 2.0 + L[k + 1] / 2.0)

        # Применяем: все слайсы вдоль касательной центрального, с его yaw.
        for k in range(5):
            r = layout_rows[i + k]
            r["x"] = cx + d[k] * tx
            r["y"] = cy + d[k] * ty
            r["z"] = 0.0
            r["yaw"] = cyaw
            r["pitch"] = 0.0
            r["roll"] = 0.0

        groups_modified += 1
        i += 5
    return groups_modified


def is_sp_start(line_names, i):
    """SP начинается с последовательности из SP_START_PATTERN (из конфига)."""
    pat = SP_START_PATTERN
    if i + len(pat) > len(line_names):
        return False
    return all(line_names[i + k] == pat[k] for k in range(len(pat)))


# ---------------------------------------------------------------------------
# Главная программа
# ---------------------------------------------------------------------------

def main():
    print(f"Reading {LTE_DESC}")
    devices = parse_desc(LTE_DESC)
    print(f"  parsed {len(devices)} device definitions")

    print(f"Reading {LTE_LINE}")
    line = parse_line(LTE_LINE)
    print(f"  parsed line with {len(line)} element references")

    write_devices_csv(devices, OUT_DEVICES)
    print(f"Wrote {OUT_DEVICES}")

    write_groups_csv(devices, OUT_GROUPS)
    print(f"Wrote {OUT_GROUPS}")

    # Сводка по группам
    counts = {}
    for d in devices.values():
        hw = map_to_hw(d)
        if not hw:
            counts.setdefault("(skipped)", 0)
            counts["(skipped)"] += 1
            continue
        for h in hw:
            counts[h[0]] = counts.get(h[0], 0) + 0  # счёт по уникальным LTE-устройствам не информативен
    # лучше посчитаем сколько LTE-устройств каждого hw-целевого типа
    hw_totals = {}
    for d in devices.values():
        for h in map_to_hw(d):
            hw_totals[h[0]] = hw_totals.get(h[0], 0) + 1
    print("\nLTE definitions -> hardware library:")
    for k in sorted(hw_totals):
        print(f"  {k:18s}  {hw_totals[k]:4d}  (LTE definitions mapping to this)")

    # Загружаем существующие arg_val (чтобы сохранить ручные правки)
    existing_args = load_existing_sp_elements(OUT_SP_ELEMENTS)
    if existing_args:
        print(f"  Loaded {len(existing_args)} existing arg_val entries from {OUT_SP_ELEMENTS}")

    # =============== ЗАМЫКАНИЕ КОЛЬЦА ===============
    # LTE-файл даёт суммарный bend ~357.59° (deficit ~2.41°). Чтобы геометрия
    # замкнулась в 360°, недостающие градусы распределяем ТОЛЬКО по ЧИСТЫМ
    # ДИПОЛЯМ (DL/JL/SBM семейства, K1=0). Combined dipole-quadrupole'ы (DQ_*,
    # K1!=0) НЕ трогаем, т.к. их bend связан с фокусировкой через единое поле,
    # и подкручивать только bend, не меняя K1 -- физически невозможно.
    raw_pure_bend = 0.0   # Σ ANGLE по чистым диполям (K1=0)
    raw_dq_bend   = 0.0   # Σ ANGLE по combined dipole-quadrupole (K1!=0)
    for name in line:
        if name not in devices:
            continue
        d = devices[name]
        angle = d.get("ANGLE", 0.0)
        if not isinstance(angle, (int, float)) or angle == 0.0:
            continue
        k1 = d.get("K1", 0.0)
        if isinstance(k1, (int, float)) and k1 != 0.0:
            raw_dq_bend += angle
        else:
            raw_pure_bend += angle
    raw_total_bend = raw_pure_bend + raw_dq_bend
    if abs(raw_pure_bend) > 1e-9:
        target_pure = (2.0 * math.pi) - raw_dq_bend
        close_scale = target_pure / raw_pure_bend
        deficit_deg = 360.0 - math.degrees(raw_total_bend)
        print(f"\n  Closing the ring (pure dipoles only, DQ unchanged):")
        print(f"  Raw sum(ANGLE): {math.degrees(raw_total_bend):.4f} deg "
              f"({math.degrees(raw_pure_bend):.4f} pure + {math.degrees(raw_dq_bend):.4f} DQ)")
        print(f"  Deficit = {deficit_deg:+.4f} deg")
        print(f"  Scaling PURE dipole ANGLE x {close_scale:.8f} (DQ K1!=0 unchanged)")
        # Скалируем angle только у чистых диполей (K1=0).
        n_scaled = 0; n_skipped = 0
        for d in devices.values():
            angle = d.get("ANGLE", 0.0)
            if not isinstance(angle, (int, float)) or angle == 0.0:
                continue
            k1 = d.get("K1", 0.0)
            if isinstance(k1, (int, float)) and k1 != 0.0:
                n_skipped += 1
            else:
                d["ANGLE"] = angle * close_scale
                n_scaled += 1
        print(f"  Scaled {n_scaled} pure-dipole device definitions, "
              f"skipped {n_skipped} DQ-type definitions")

    # Обход линии: одновременно собираем шаблоны SP и развёртку
    layout_rows = []
    sp_templates = {}        # sp_type -> [{sp_element_id, device_type_id, arg_val}, ...]
    x, y, theta = 0.0, 0.0, 0.0
    global_id = 0
    sp_number = 0            # 1..40 (наращивается при is_sp_start)
    sp_element_id = 0        # сбрасывается на каждой границе SP
    current_sp_type = None   # определяется по первому магниту в SP
    unknown = {}
    total_angle = 0.0
    total_length = 0.0
    hw_uses = {}

    for i, name in enumerate(line):
        if name not in devices:
            unknown[name] = unknown.get(name, 0) + 1
            continue
        dev = devices[name]

        # Граница нового супериода?
        if is_sp_start(line, i):
            sp_number += 1
            sp_element_id = 0
            current_sp_type = None

        # Обновляем геометрию для любого элемента (включая дрифты и BPM)
        L = dev.get("L", 0.0)
        angle = dev.get("ANGLE", 0.0)
        xc, yc, yaw_c, x, y, theta = step(x, y, theta, L, angle)
        total_length += L
        total_angle += angle

        hw_list = map_to_hw(dev)
        if not hw_list:
            continue  # дрифт / BPM / skip

        # Определяем тип SP по первому магниту в нём
        if current_sp_type is None:
            current_sp_type = FIRST_ELEMENT_TO_SP_TYPE.get(name, name)
            if current_sp_type not in sp_templates:
                sp_templates[current_sp_type] = []

        # Эмитим строки -- по одной на каждую hardware-запись
        for hw_type, default_arg in hw_list:
            sp_element_id += 1
            global_id += 1
            hw_uses[hw_type] = hw_uses.get(hw_type, 0) + 1

            # Дополнить шаблон SP, если это новая (sp_type, sp_element_id) пара
            template = sp_templates[current_sp_type]
            if sp_element_id > len(template):
                # Сохраняем существующий arg_val + correct_coeff, если были
                old = existing_args.get(
                    (current_sp_type, sp_element_id, hw_type),
                    (default_arg, "1.0")
                )
                arg_val, correct_coeff = old
                template.append({
                    "sp_element_id":  sp_element_id,
                    "device_type_id": hw_type,
                    "arg_val":        arg_val,
                    "correct_coeff":  correct_coeff,
                })

            # Развёртка: только геометрия, без токов.
            # Сохраняем lte_name для пост-обработки DL/JL групп.
            layout_rows.append({
                "global_id": global_id,
                "sp_number": sp_number,
                "sp_type": current_sp_type,
                "sp_element_id": sp_element_id,
                "x": xc, "y": yc, "z": 0.0,
                "yaw": yaw_c, "pitch": 0.0, "roll": 0.0,
                "_lte_name": name,
            })

    if unknown:
        print(f"\nWARN: unknown elements in line: {dict(sorted(unknown.items()))}")

    print(f"\nLine summary:")
    print(f"  total path length:   {total_length:.4f} m")
    print(f"  total bend angle:    {total_angle:.6f} rad  "
          f"({math.degrees(total_angle):+.3f} deg, "
          f"{total_angle / (2*math.pi):+.4f} full turns)")
    print(f"  closing position:    x = {x:.4f}, y = {y:.4f}")
    print(f"  closing direction:   {math.degrees(theta):+.3f} deg")
    print(f"  superperiods:        {sp_number} ({', '.join(f'{t}:{len(v)}' for t, v in sp_templates.items())})")
    print(f"  layout rows:         {len(layout_rows)}")

    print("\nUsage of hardware types in lattice:")
    for k in sorted(hw_uses):
        print(f"  {k:18s}  {hw_uses[k]:5d}")

    # Пост-обработка: выпрямление DL/JL 5-слайсовых семейств.
    n_straightened = straighten_dljl_families(layout_rows, devices)
    print(f"\nStraightened {n_straightened} DL/JL families "
          f"(5 slices each, colinear at central slice position)")

    # Удаляем технический ключ _lte_name перед записью.
    for r in layout_rows:
        r.pop("_lte_name", None)

    write_sp_elements_csv(sp_templates, OUT_SP_ELEMENTS)
    print(f"\nWrote {OUT_SP_ELEMENTS}  ({sum(len(v) for v in sp_templates.values())} rows, "
          f"{len(sp_templates)} SP type(s))")
    write_lattice_layout_csv(layout_rows, OUT_LATTICE_LAYOUT)
    print(f"Wrote {OUT_LATTICE_LAYOUT}  ({len(layout_rows)} rows)")
    print(f"\nNext: run assemble_lattice_configs.py to produce lattice_configs.csv")


if __name__ == "__main__":
    main()
