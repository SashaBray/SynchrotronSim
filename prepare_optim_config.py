"""
Конвертер sp1_endpoint.json -> flat-текстовый конфиг для SynchrotronOptimizer.exe.

Запуск:
    python prepare_optim_config.py workspace/optim_inputs/sp1_endpoint.json
        -> создаёт workspace/optim_inputs/sp1_endpoint.optimconfig

Формат вывода (line-based, пробелы как разделители, ';' как разделитель полей в строке GROUP):

  # Комментарий
  SCALAR <key> <value>
  GROUP;<name>;<step_abs>;<step_percent>;<target>;<file>;<value_column>;<n_matchers>;<col=val>;<col=val>;...

target: CURRENT (модифицирует MagneticElement::currentScale = arg_val / nominal_arg)
        COEFF   (модифицирует MagneticElement::correctCoeff = correct_coeff)

Матчеры: каждый "<col>=<val>" -- ограничение на (sp_type, sp_element_id, device_type_id).
Можно комбинировать sp_type, sp_element_id, device_type_id для per-position групп.
Для DP_correct_coeff: матчер по device_type_id (одно или несколько значений; в этом
случае запись повторяется -- одна GROUP-строка с несколькими альтернативами для одного матчера).

Для поддержки key_values как массива (несколько device_type_id) используем
запись "col=val1|val2|val3" (через '|' -- OR в значении).
"""
import csv
import json
import os
import sys
from pathlib import Path


def expand_per_position_groups(cfg_pp, sp_elements_path):
    """Авто-генерирует группы per-position из sp_elements.csv.
    cfg_pp может быть либо одним dict, либо списком dict-ов.
    Возвращает список dict аналогично parameter_groups."""
    if cfg_pp is None:
        return []
    # Допускаем список блоков (например один для arg_val, один для correct_coeff)
    if isinstance(cfg_pp, list):
        all_groups = []
        for block in cfg_pp:
            all_groups.extend(expand_per_position_groups(block, sp_elements_path))
        return all_groups

    if not cfg_pp.get("enabled", True):
        return []
    sp_types_filter = set(cfg_pp.get("sp_types") or [])
    exclude         = set(cfg_pp.get("exclude_device_types") or [])
    include         = set(cfg_pp.get("include_device_types") or [])  # whitelist (если задан)
    value_col       = cfg_pp.get("value_column", "arg_val")
    step_abs        = float(cfg_pp.get("step_abs", 0)) if "step_abs" in cfg_pp else 0.0
    step_percent    = float(cfg_pp.get("step_percent", 0)) if "step_percent" in cfg_pp else 0.0
    name_suffix     = cfg_pp.get("name_suffix", "")  # для различения групп с разными value_column
    groups = []
    with open(sp_elements_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            sp_type = row["sp_type"].strip()
            sp_eid  = row["sp_element_id"].strip()
            dev     = row["device_type_id"].strip()
            if sp_types_filter and sp_type not in sp_types_filter:
                continue
            if include and dev not in include:
                continue
            if dev in exclude:
                continue
            name = f"{sp_type[:3]}_pos{int(sp_eid):03d}_{dev}{name_suffix}"
            groups.append({
                "name": name,
                "enabled": True,
                "file": "sp_elements.csv",
                "value_column": value_col,
                "step_abs": step_abs if step_abs > 0 else None,
                "step_percent": step_percent if step_percent > 0 else None,
                "_compound_selectors": [
                    {"column": "sp_type",         "values": [sp_type]},
                    {"column": "sp_element_id",   "values": [sp_eid]},
                    {"column": "device_type_id",  "values": [dev]},
                ],
            })
    return groups


def normalize_group(g):
    """Унифицирует структуру: возвращает (name, file, value_col, step_abs, step_percent, selectors).
    selectors -- список dict {column, values}."""
    name      = g["name"]
    file_     = g["file"]
    val_col   = g["value_column"]
    step_abs  = g.get("step_abs")
    step_pct  = g.get("step_percent")
    enabled   = g.get("enabled", True)

    # Compound селектор (для per-position)
    if "_compound_selectors" in g:
        sels = g["_compound_selectors"]
    else:
        # классический одиночный селектор: key_column + key_values
        sels = [{"column": g["key_column"], "values": list(g["key_values"])}]

    return {
        "name": name,
        "file": file_,
        "value_col": val_col,
        "step_abs": float(step_abs) if step_abs is not None else 0.0,
        "step_percent": float(step_pct) if step_pct is not None else 0.0,
        "enabled": bool(enabled),
        "selectors": sels,
    }


def resolve_extra_targets(extra_target_magnets, layout_path):
    """Преобразует extra_target_magnets -> [(x, y, z, weight), ...]

    Поддерживаются два формата (вес опционален, по умолчанию 1.0):
        {"sp_number": N, "sp_element_id": E, "weight": W}  -- lookup в layout
        {"x": X, "y": Y, "z": Z, "weight": W}              -- raw координаты
    Можно смешивать в одном списке.
    """
    if not extra_target_magnets:
        return []
    targets = []
    rows_cache = None  # lazy load lattice_layout
    for et in extra_target_magnets:
        weight = float(et.get("weight", 1.0))
        # Raw coords -- shortcut
        if "x" in et and "y" in et and "z" in et:
            targets.append((float(et["x"]), float(et["y"]), float(et["z"]), weight))
            continue
        # Если нет ни x/y/z, ни sp_number -- это комментарий (_doc), пропускаем.
        if "sp_number" not in et or "sp_element_id" not in et:
            continue
        # Lookup by (sp_number, sp_element_id)
        if rows_cache is None:
            if not layout_path.exists():
                print(f"WARN: lattice_layout.csv not found at {layout_path}; skipping name-based targets")
                continue
            with open(layout_path, encoding="utf-8") as f:
                rows_cache = list(csv.DictReader(f, delimiter=";"))
        spn = int(et["sp_number"])
        seid = int(et["sp_element_id"])
        match = None
        for r in rows_cache:
            try:
                if int(r["sp_number"]) == spn and int(r["sp_element_id"]) == seid:
                    match = r
                    break
            except (KeyError, ValueError):
                continue
        if match is None:
            print(f"WARN: extra target (sp_number={spn}, sp_element_id={seid}) not found")
            continue
        targets.append((float(match["x"]), float(match["y"]), float(match["z"]), weight))
    return targets


def write_flat_config(cfg, out_path, sp_elements_path, history_output):
    explicit = list(cfg.get("parameter_groups", []))
    auto = expand_per_position_groups(cfg.get("parameter_groups_per_position"), sp_elements_path)
    all_groups = [normalize_group(g) for g in (explicit + auto)]
    # Сохраняем только enabled
    groups = [g for g in all_groups if g["enabled"]]

    loss = cfg.get("loss", {})
    sp_nums = loss.get("magnet_sp_numbers", [1])
    sp_nums_str = " ".join(str(int(s)) for s in sp_nums)

    # Список device_type_id, исключаемых из набора магнитов для loss
    # (например, ["DP_L_PM", "DP_sh_PM"] -- постоянные магниты не подстраиваются).
    exclude_devs = loss.get("exclude_device_types") or []
    exclude_devs_str = " ".join(str(d) for d in exclude_devs)

    # Резолвинг "extra_target_magnets" -> координаты из lattice_layout.csv
    layout_path = Path(sp_elements_path).parent / "lattice_layout.csv"
    extra_targets = resolve_extra_targets(loss.get("extra_target_magnets"), layout_path)

    lines = []
    lines.append("# Flat config for SynchrotronOptimizer.exe (auto-generated)")
    lines.append("# Do NOT edit manually -- regenerate via prepare_optim_config.py")
    lines.append("")
    lines.append(f"SCALAR EXP_ID {int(cfg['experiment_id'])}")
    lines.append(f"SCALAR N_ITER {int(cfg['n_iterations'])}")
    lines.append(f"SCALAR SEED {int(cfg.get('random_seed', 0))}")
    lines.append(f"SCALAR KICK_EVERY {int(cfg.get('kick_every', 0))}")
    lines.append(f"SCALAR KICK_SCALE {float(cfg.get('kick_scale', 0.0))}")
    lines.append(f"SCALAR REVERT_AFTER_KICKS {int(cfg.get('revert_after_kicks', 3))}")
    lines.append(f"SCALAR LOSS_TYPE {loss.get('type', 'orbit_avg')}")
    lines.append(f"SCALAR LOSS_SQUARED {1 if loss.get('squared', False) else 0}")
    lines.append(f"SCALAR LOSS_SP {sp_nums_str}")
    if exclude_devs_str:
        lines.append(f"SCALAR LOSS_EXCLUDE_DEVICES {exclude_devs_str}")
    # Опциональная "обрезанная" lattice: оптимизатор использует только
    # элементы lattice_layout с global_id <= этого значения.
    max_gid = cfg.get("lattice_max_gid")
    if max_gid is not None and int(max_gid) > 0:
        lines.append(f"SCALAR LATTICE_MAX_GID {int(max_gid)}")
    lines.append(f"SCALAR HISTORY_OUT {history_output}")
    if extra_targets:
        lines.append("")
        lines.append("# Extra explicit targets. Format: TARGET <x> <y> <z> [weight].")
        lines.append("# Weight 1.0 = same influence as one lattice magnet in loss.")
        for x, y, z, w in extra_targets:
            if abs(w - 1.0) < 1e-12:
                lines.append(f"TARGET {x:.10g} {y:.10g} {z:.10g}")
            else:
                lines.append(f"TARGET {x:.10g} {y:.10g} {z:.10g} {w:.10g}")

    # Угловая компонента loss: финальный импульс должен быть параллелен
    # эталонному направлению (= касательной опорной орбиты в TARGET'е).
    dir_loss = loss.get("direction_loss") or {}
    ref = dir_loss.get("reference")
    dir_w = float(dir_loss.get("weight", 0.0))
    if ref and isinstance(ref, list) and len(ref) == 3 and dir_w > 0:
        lines.append("")
        lines.append("# Direction loss term. Penalizes angle between final momentum")
        lines.append("# and reference orbit direction. Loss_dir = weight × Σ sin²(θ_i).")
        lines.append(f"SCALAR DIRECTION_REF {float(ref[0]):.10g} {float(ref[1]):.10g} {float(ref[2]):.10g}")
        lines.append(f"SCALAR DIRECTION_WEIGHT {dir_w:.10g}")

    # Envelope-loss: штраф за рост амплитуды бетатронных колебаний.
    env_loss = loss.get("envelope_loss") or {}
    env_w = float(env_loss.get("weight", 0.0))
    if env_w > 0:
        lines.append("")
        lines.append("# Envelope-growth loss term. Penalizes growth of betatron")
        lines.append("# oscillation amplitude (only growth, not damping).")
        lines.append("# Loss_env = weight × Σ_satellites max(0, A_final - A_init)².")
        lines.append(f"SCALAR ENVELOPE_WEIGHT {env_w:.10g}")
    lines.append("")
    lines.append("# GROUP;<name>;<step_abs>;<step_pct>;<target>;<file>;<value_col>;<n_matchers>;<col=val>;...")
    # target:
    #   CURRENT     -- arg_val в sp_elements (ток электромагнитов; el.currentScale)
    #   COEFF_HW    -- correct_coeff в hardware_library (device-level калибровка)
    #   COEFF_INST  -- correct_coeff в sp_elements (per-instance калибровка;
    #                   умножается с COEFF_HW при загрузке)
    for g in groups:
        val_col = g["value_col"]
        if val_col == "correct_coeff":
            target = "COEFF_HW" if g["file"] == "hardware_library.csv" else "COEFF_INST"
        else:
            target = "CURRENT"
        # Каждый selector кодируем как "col=val1|val2|..."
        match_strs = []
        for s in g["selectors"]:
            vals = "|".join(s["values"])
            match_strs.append(f"{s['column']}={vals}")
        n_matchers = len(match_strs)
        parts = ["GROUP",
                 g["name"],
                 f"{g['step_abs']:.10g}",
                 f"{g['step_percent']:.10g}",
                 target,
                 g["file"],
                 g["value_col"],
                 str(n_matchers)] + match_strs
        lines.append(";".join(parts))

    Path(out_path).write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {out_path}")
    print(f"  {len(groups)} groups (out of {len(all_groups)} total, {len(all_groups) - len(groups)} disabled)")
    print(f"  experiment_id={cfg['experiment_id']}  n_iter={cfg['n_iterations']}")


def main():
    if len(sys.argv) < 2:
        cfg_path = "workspace/optim_inputs/sp1_endpoint.json"
    else:
        cfg_path = sys.argv[1]
    cfg_path = Path(cfg_path).resolve()
    out_path = cfg_path.with_suffix(".optimconfig")
    history_output = str(cfg_path.with_suffix(".history.json").name)
    with open(cfg_path, encoding="utf-8") as f:
        cfg = json.load(f)
    sp_elements_path = Path(cfg_path).parent.parent.parent / "workspace" / "tables" / "sp_elements.csv"
    # ^ предполагает workspace/optim_inputs/*.json -> workspace/tables/sp_elements.csv
    # На случай нестандартного расположения позволим override через окружение
    sp_env = os.environ.get("SP_ELEMENTS_PATH")
    if sp_env:
        sp_elements_path = Path(sp_env)
    if not sp_elements_path.exists():
        # ещё одна попытка: относительно скрипта
        guess = Path(__file__).resolve().parent / "workspace" / "tables" / "sp_elements.csv"
        if guess.exists():
            sp_elements_path = guess
    write_flat_config(cfg, out_path, sp_elements_path, history_output)


if __name__ == "__main__":
    main()
