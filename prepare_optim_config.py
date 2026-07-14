"""
Конвертер sp1_endpoint.json -> flat-текстовый конфиг для SynchrotronOptimizer.exe.

Запуск:
    python prepare_optim_config.py workspace/optim_inputs/sp1_endpoint.json
        -> создаёт workspace/optim_inputs/sp1_endpoint.optimconfig
        -> создаёт workspace/results/design_orbit_lat{lattice_id}.csv
           (опорная орбита для envelope_loss; собирается из БД сегментов
           design_orbit_segments.csv, см. design_orbit.py)

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


def resolve_lattice_id(cfg, ws):
    """Возвращает lattice_id для эксперимента.

    Источники по приоритету:
      1. cfg['lattice_id'] (явный override в JSON)
      2. experiments.csv: lookup по experiment_id
    Если оба отсутствуют -- возвращает None (envelope-loss отключится).
    """
    if "lattice_id" in cfg:
        return int(cfg["lattice_id"])
    exp_id = int(cfg["experiment_id"])
    exp_csv = ws / "tables" / "experiments.csv"
    if not exp_csv.exists():
        return None
    try:
        with open(exp_csv, encoding="utf-8") as f:
            for row in csv.DictReader(f, delimiter=";", skipinitialspace=True):
                try:
                    if int(row["exp_id"]) == exp_id:
                        return int(row["lattice_id"])
                except (KeyError, ValueError):
                    continue
    except Exception as e:
        print(f"  WARN: cannot read {exp_csv}: {e}")
    return None


def generate_magnet_endpoints_csv(lattice_id, ws, max_gid=None, min_gid=None):
    """Собирает вход/выход всех магнитов lattice на эталонной орбите и пишет в
    workspace/results/magnet_endpoints_lat{id}.csv. Используется per_element-
    режимом C++ для постройки start_pose и finish_plane на каждом окне.

    min_gid: нижняя граница global_id (для участков, начинающихся НЕ с gid 1 --
    напр. junction в середине кольца: пучок стартует со входа супериода, и опорный
    путь / окна должны охватывать только этот участок, а не префикс SP1..предыдущие).

    Возвращает (path, n_rows) либо (None, 0) при неудаче.
    """
    try:
        from design_orbit import compute_magnet_endpoints
    except Exception as e:
        print(f"  WARN: модуль design_orbit недоступен: {e}")
        return None, 0
    out_dir = ws / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"magnet_endpoints_lat{lattice_id}.csv"
    try:
        eps, _final_yaw = compute_magnet_endpoints(lattice_id, workspace=ws)
    except Exception as e:
        print(f"  WARN: не удалось собрать magnet endpoints для lattice "
              f"{lattice_id}: {e}")
        return None, 0
    # Фильтрация по max_gid (sub_lattice) и min_gid (участок не с начала кольца)
    if max_gid is not None and max_gid > 0:
        eps = [e for e in eps if e.global_id <= max_gid]
    if min_gid is not None and min_gid > 0:
        eps = [e for e in eps if e.global_id >= min_gid]
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator="\n")
        w.writerow(["global_id", "device_type_id", "lte_name",
                    "x_in", "y_in", "yaw_in",
                    "x_out", "y_out", "yaw_out"])
        for e in eps:
            w.writerow([
                e.global_id, e.device_type_id, e.lte_name,
                f"{e.entry[0]:.10f}", f"{e.entry[1]:.10f}", f"{e.entry_yaw:.10f}",
                f"{e.exit [0]:.10f}", f"{e.exit [1]:.10f}", f"{e.exit_yaw :.10f}",
            ])
    print(f"  Wrote {out_csv.name}  ({len(eps)} magnets)")
    return out_csv, len(eps)


def generate_design_orbit_csv(lattice_id, ws, points_per_segment=10):
    """Собирает design orbit для lattice_id и пишет в
    workspace/results/design_orbit_lat{id}.csv.

    Возвращает (path, n_points) либо (None, 0) при неудаче.
    """
    try:
        from design_orbit import build_orbit_for_lattice
    except Exception as e:
        print(f"  WARN: модуль design_orbit недоступен: {e}")
        return None, 0
    out_dir = ws / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_csv = out_dir / f"design_orbit_lat{lattice_id}.csv"
    try:
        pts, final_yaw = build_orbit_for_lattice(
            lattice_id, workspace=ws, points_per_segment=points_per_segment)
    except Exception as e:
        print(f"  WARN: не удалось собрать design orbit для lattice {lattice_id}: {e}")
        return None, 0
    with open(out_csv, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator="\n")
        w.writerow(["x", "y", "z"])
        for p in pts:
            w.writerow([f"{p[0]:.8f}", f"{p[1]:.8f}", f"{p[2]:.8f}"])
    print(f"  Wrote {out_csv.name}  ({len(pts)} pts, "
          f"final_yaw={final_yaw:.6f} rad)")
    return out_csv, len(pts)


def write_flat_config(cfg, out_path, sp_elements_path, history_output, ws=None):
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
    # max_weight: штраф за МАКСИМАЛЬНОЕ отклонение (loss += w*max(d^2)). 0 = выкл.
    lines.append(f"SCALAR LOSS_MAX_WEIGHT {float(loss.get('max_weight', 0.0)):.10g}")
    # interparticle_weight: штраф за расхождение пучка (Σ по парам частиц |p_i-p_j|²
    # в контрольных точках). Масштаб ~как pos. 0 = выкл.
    lines.append(f"SCALAR INTERPARTICLE_WEIGHT {float(loss.get('interparticle_weight', 0.0)):.10g}")
    # envelope_decay_weight: штраф за РОСТ амплитуды взаимного разлёта частиц
    # S_m = Σ_{i<j}|p_i-p_j|² выше стартовой S_init -> огибающая должна угасать.
    # decay = Σ_m max(0, S_m - S_init). Активен в ОКОННОМ режиме. 0 = выкл.
    lines.append(f"SCALAR ENVELOPE_DECAY_WEIGHT {float(loss.get('envelope_decay_weight', 0.0)):.10g}")
    # parallelism_weight: штраф за УГЛОВУЮ расходимость пучка -- Σ_{i<j}|p̂_i-p̂_j|²
    # (направления импульса частиц должны быть параллельны). Активен в оконном режиме.
    lines.append(f"SCALAR PARALLELISM_WEIGHT {float(loss.get('parallelism_weight', 0.0)):.10g}")
    # svd_angle_weight: вес ПЕРИОДИЧЕСКОГО члена в невязке SVD (выход супериода ∥ орбите).
    lines.append(f"SCALAR SVD_ANGLE_WEIGHT {float(loss.get('svd_angle_weight', 0.0)):.10g}")
    # global_search_mode: МЕТОД поиска в ГЛОБАЛЬНОМ режиме.
    # "monte_carlo" (алиас "random") -- случайные greedy-пертурбации (Монте-Карло);
    # "line_search" (default) -- последовательное приближение: покоординатный спуск
    # с детерминированным 1D-поиском по всей конфигурации.
    lines.append(f"SCALAR GLOBAL_SEARCH_MODE {cfg.get('global_search_mode', 'line_search')}")
    # position_reference: "iron" (default, центры железа из lattice_layout.csv)
    # или "orbit" (середина дуги магнита из magnet_endpoints CSV). Для не-bending
    # магнитов (квадруполи, секступоли, октуполи) разницы нет; для диполей -- до 17 мм.
    pos_ref = str(loss.get("position_reference", "iron")).lower()
    lines.append(f"SCALAR LOSS_POSITION_REFERENCE {pos_ref}")
    lines.append(f"SCALAR LOSS_SP {sp_nums_str}")
    if exclude_devs_str:
        lines.append(f"SCALAR LOSS_EXCLUDE_DEVICES {exclude_devs_str}")
    # Опциональная "обрезанная" lattice: оптимизатор использует только
    # элементы lattice_layout с global_id <= этого значения.
    max_gid = cfg.get("lattice_max_gid")
    if max_gid is not None and int(max_gid) > 0:
        lines.append(f"SCALAR LATTICE_MAX_GID {int(max_gid)}")
    # max_particles: моделировать только первые N частиц (фаза поворотных = 1,
    # центральная). 0/отсутствует = все частицы группы.
    max_parts = cfg.get("max_particles")
    if max_parts is not None and int(max_parts) > 0:
        lines.append(f"SCALAR MAX_PARTICLES {int(max_parts)}")
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

    # ---------------- design orbit + magnet endpoints ----------------
    # Эти два CSV нужны нескольким режимам: envelope_loss (design_orbit),
    # per_element (magnet_endpoints). Поэтому генерим всегда, если lattice_id
    # резолвится. Чёткий разрыв с прежней логикой: даже без envelope_loss
    # C++ получит путь к опорной орбите (это не повредит — она используется
    # только если ENVELOPE_WEIGHT > 0 или включён per_element).
    if ws is not None:
        env_loss = loss.get("envelope_loss") or {}
        lat_id = resolve_lattice_id(cfg, ws)
        if lat_id is None:
            print("  WARN: lattice_id не резолвится; design orbit / endpoints "
                  "не сгенерированы.")
        else:
            pts_per_seg = int(env_loss.get("points_per_segment", 10))
            max_gid = cfg.get("lattice_max_gid")
            max_gid = int(max_gid) if (max_gid is not None and int(max_gid) > 0) else None
            # lattice_min_gid: нижняя граница (участок не с начала кольца, напр. junction).
            # Опорный путь и окна охватывают [min_gid, max_gid]; пучок стартует со
            # входа этого участка (particle_groups).
            min_gid = cfg.get("lattice_min_gid")
            min_gid = int(min_gid) if (min_gid is not None and int(min_gid) > 0) else None
            orbit_csv, _n_pts = generate_design_orbit_csv(
                lat_id, ws, points_per_segment=pts_per_seg)
            eps_csv,  _n_eps = generate_magnet_endpoints_csv(
                lat_id, ws, max_gid=max_gid, min_gid=min_gid)
            lines.append("")
            lines.append("# Design orbit + magnet endpoints (paths to CSVs):")
            if orbit_csv is not None:
                lines.append(f"STRING DESIGN_ORBIT_CSV {orbit_csv.as_posix()}")
            if eps_csv is not None:
                lines.append(f"STRING MAGNET_ENDPOINTS_CSV {eps_csv.as_posix()}")

            # Reference path (PATHPT) для composite_next: ПЛОТНЫЙ опорный путь --
            # центры элементов + уплотнённые прямые участки (включая 50м-дрейф).
            # Loss в composite_next считает отклонение от ЭТИХ точек (а не от
            # центров магнитов): частые точки на прямых задают НАПРАВЛЕНИЕ, поэтому
            # "выстрелить" одним магнитом нельзя -- это сразу уводит с линии.
            pe_cfg = cfg.get("per_element") or {}
            if pe_cfg.get("loss_mode") == "composite_next" and eps_csv is not None:
                try:
                    from build_reference_path import load_endpoints, build_path
                    p_els = load_endpoints(eps_csv)
                    p_spacing = float(pe_cfg.get("path_spacing", 0.25))
                    p_dend = (extra_targets[-1][0], extra_targets[-1][1]) if extra_targets else None
                    p_dgid = p_els[-1]["gid"] if p_els else None
                    # Прогрессивные gid'ы дрейфа -- только при включённой фазе
                    # расширения (drift_extension_m>0). Иначе дрейф не уплотняем.
                    p_chunk = float(pe_cfg.get("drift_chunk_m", 2.0)) \
                        if float(pe_cfg.get("drift_extension_m", 0.0)) > 0 else 0.0
                    p_pts = build_path(p_els, spacing=p_spacing,
                                       drift_end=p_dend, drift_gid=p_dgid,
                                       drift_chunk_m=p_chunk,
                                       drift_point_spacing=float(pe_cfg.get("drift_point_spacing", 1.0)))
                    lines.append("")
                    lines.append(f"# Reference path (composite_next): {len(p_pts)} pts. "
                                 "Format: PATHPT <gid> <x> <y> <z>")
                    for pg, px, py, pz in p_pts:
                        lines.append(f"PATHPT {pg} {px:.10g} {py:.10g} {pz:.10g}")
                    print(f"  Reference path: {len(p_pts)} PATHPT-точек (composite_next)")
                except Exception as e:
                    print(f"  WARN: reference path не сгенерирован: {e}")

    # Envelope-loss: штраф за рост амплитуды бетатронных колебаний.
    env_loss = loss.get("envelope_loss") or {}
    env_w = float(env_loss.get("weight", 0.0))
    if env_w > 0:
        lines.append("")
        lines.append("# Envelope-growth loss term. Penalizes growth of betatron")
        lines.append("# oscillation amplitude (only growth, not damping).")
        lines.append("# Loss_env = weight × Σ_satellites max(0, A_final - A_init)².")
        lines.append(f"SCALAR ENVELOPE_WEIGHT {env_w:.10g}")

    # РЕЖИМ оптимизации: "global" (default) | "windowed" | "gradient".
    #   global   -- настраивается вся конфигурация сразу;
    #   windowed -- скользящее окно по магнитам (оконный режим);
    #   gradient -- FD + Adam.
    # Поле JSON: optimization_mode (предпочтительно) или optimization_method (legacy).
    # CLI-флаг --method=/--mode= имеет приоритет над JSON. Историч. алиасы:
    # montecarlo/mc -> global, per_element/pe -> windowed.
    # ВНИМАНИЕ: Монте-Карло -- это МЕТОД (случайные кики), а не режим; задаётся
    # через global_search_mode / PER_ELEM_SEARCH_MODE.
    _mode_raw = str(cfg.get("optimization_mode",
                            cfg.get("optimization_method", "global"))).lower()
    _MODE_ALIASES = {
        "global": "global", "montecarlo": "global", "mc": "global",
        "windowed": "windowed", "window": "windowed",
        "per_element": "windowed", "perel": "windowed", "pe": "windowed",
        "lens": "lens", "lenses": "lens",
        "gradient": "gradient", "grad": "gradient",
    }
    mode = _MODE_ALIASES.get(_mode_raw, _mode_raw)
    lines.append("")
    lines.append("# Optimization MODE: 'global' | 'windowed' | 'gradient'.")
    lines.append(f"SCALAR METHOD {mode}")

    # Режим loss: "composite" (default; position + direction + envelope) или
    # "orbit" (только position; для грубой начальной настройки).
    # JSON поле loss.mode или alias loss_mode. CLI-флаг --loss=composite|orbit
    # имеет приоритет.
    loss_mode = str(loss.get("mode") or cfg.get("loss_mode") or "composite").lower()
    lines.append(f"# Loss mode: 'composite' (position+direction+envelope) "
                 f"или 'orbit' (только position от опорной орбиты)")
    lines.append(f"SCALAR LOSS_MODE {loss_mode}")

    # Параметры dynamic loss reweighting (только MC-режим).
    # Веса компонентов composite (position, direction, envelope) вращаются по
    # синусоиде с равномерно распределёнными фазами, что помогает MC выбираться
    # из локальных минимумов "статически-сбалансированной" loss-функции.
    dyn = cfg.get("dynamic_loss") or {}
    if dyn:
        lines.append("")
        lines.append("# === DYNAMIC LOSS REWEIGHTING (MC only) ===")
        lines.append("# enabled:   1 = включить, 0 = классический static composite (как раньше).")
        lines.append("# amplitude: A в формуле w_i = base × [1 + A·cos(2πt/T - φ_i)]. 0..1.")
        lines.append("# period:    T в итерациях (0 = auto 5×N_params).")
        lines.append(f"SCALAR DYN_LOSS_ENABLED    {1 if dyn.get('enabled', False) else 0}")
        lines.append(f"SCALAR DYN_LOSS_AMPLITUDE  {float(dyn.get('amplitude', 0.7)):.10g}")
        lines.append(f"SCALAR DYN_LOSS_PERIOD     {int(dyn.get('period_iterations', 0))}")

    # Параметры ОКОННОГО (windowed) режима (скользящее окно настраиваемых магнитов).
    # JSON-блок исторически называется "per_element"; параметры line_search_* /
    # диапазоны переиспользуются и ГЛОБАЛЬНЫМ режимом. Активны при method=='windowed'.
    # CLI-флаг --method=windowed имеет приоритет над JSON.
    pe = cfg.get("per_element") or {}
    if pe:
        lines.append("")
        lines.append("# === PER-ELEMENT MODE (sliding window) ===")
        lines.append("# window_size:      сколько магнитов в одном окне (W).")
        lines.append("# window_stride:    на сколько сдвигать окно за шаг.")
        lines.append("# iters_per_window: m -- число MC-итераций на каждом окне.")
        lines.append("# passes:           сколько раз пройти всю решётку.")
        lines.append(f"SCALAR PER_ELEM_WINDOW_SIZE       {int(pe.get('window_size', 5))}")
        lines.append(f"SCALAR PER_ELEM_WINDOW_STRIDE     {int(pe.get('window_stride', 1))}")
        lines.append(f"SCALAR PER_ELEM_ITERS_PER_WINDOW  {int(pe.get('iters_per_window', 200))}")
        lines.append(f"SCALAR PER_ELEM_PASSES            {int(pe.get('passes', 1))}")
        # start_from_sp_start: false (default) -> идеальный старт с entry окна.
        # true -> частицы стартуют с начала супериода, к которому принадлежит
        # окно (накопленный пучок). Симуляция дольше, но физически реалистичнее.
        lines.append(f"SCALAR PER_ELEM_START_FROM_SP_START  "
                     f"{1 if pe.get('start_from_sp_start', False) else 0}")
        # path_only=true -> per-window loss = ТОЛЬКО Σ d²(P1, design_polyline).
        # Чистая сумма квадратов отклонений центра пучка от опорной орбиты.
        # path_weight игнорируется в этом режиме (берётся = 1).
        # При path_only=false -- комбинированный loss (см. C++ runWindowSim).
        lines.append(f"SCALAR PER_ELEM_PATH_ONLY            "
                     f"{1 if pe.get('path_only', False) else 0}")
        lines.append(f"SCALAR PER_ELEM_PATH_WEIGHT          "
                     f"{float(pe.get('path_weight', 0.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_PATH_SAMPLE_STRIDE   "
                     f"{int(pe.get('path_sample_stride', 10))}")
        # loss_mode: "endpoint" (default), "path_only" или "composite"
        # composite = position + dirWeight·dir + envWeight·env -- но для магнитов
        # фильтруется по gid окна (extra TARGET'ы всегда учитываются).
        lines.append(f"SCALAR PER_ELEM_LOSS_MODE            "
                     f"{str(pe.get('loss_mode', 'endpoint'))}")
        # accumulate_groups: если true, набор активных ручек РАСТЁТ внутри
        # pass'а -- следующее окно видит свои + все предыдущие. Между pass'ами
        # сбрасывается. iters_per_param: если >0 -- m_iter = iters_per_param ×
        # N_active_groups (динамически растёт с накоплением).
        lines.append(f"SCALAR PER_ELEM_ACCUMULATE_GROUPS    "
                     f"{1 if pe.get('accumulate_groups', False) else 0}")
        lines.append(f"SCALAR PER_ELEM_ITERS_PER_PARAM      "
                     f"{int(pe.get('iters_per_param', 0))}")
        # total_iters: общий бюджет симуляций на проход, РАВНОМЕРНО по окнам
        # (m_iter = total_iters / num_windows). Приоритет над iters_per_param. 0=выкл.
        lines.append(f"SCALAR PER_ELEM_TOTAL_ITERS          "
                     f"{int(pe.get('total_iters', 0))}")
        # iters_per_new: новый магнит окна получает фикс. число итераций, старые --
        # по остатку (m = N_new*iters_per_new + остаток/окно). Нужен total_iters>0.
        lines.append(f"SCALAR PER_ELEM_ITERS_PER_NEW        "
                     f"{int(pe.get('iters_per_new', 0))}")
        # group_dipole_arrays: окно сразу захватывает всю линейку диполей целиком
        # (слайсы одной системы тюнятся совместно), не по одному.
        lines.append(f"SCALAR PER_ELEM_GROUP_DIPOLE_ARRAYS  "
                     f"{1 if pe.get('group_dipole_arrays', False) else 0}")
        # advance_mean_threshold_mm: окно не сдвигается, пока среднее отклонение
        # траекторий (мм) не упадёт <= порога. 0 = гейт выкл. advance_max_iter_mult:
        # защитный потолок итераций = mult × базовый бюджет окна.
        lines.append(f"SCALAR PER_ELEM_ADVANCE_MEAN_MM      "
                     f"{float(pe.get('advance_mean_threshold_mm', 0.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_ADVANCE_MAX_ITER_MULT "
                     f"{int(pe.get('advance_max_iter_mult', 10))}")
        # lagged-lens: линзы (lens_lagged_device_types) включаются в тюнинг на
        # lens_lag_windows окон позже поворотных. 0 = выкл (всё сразу).
        lines.append(f"SCALAR PER_ELEM_LENS_LAG_WINDOWS     "
                     f"{int(pe.get('lens_lag_windows', 0))}")
        lens_dts = pe.get('lens_lagged_device_types') or cfg.get('lens_lagged_device_types') or []
        if lens_dts:
            lines.append(f"STRING LENS_LAGGED_DEVICE_TYPES {','.join(lens_dts)}")
        # LENS-режим: window_unit='lens' -> окно скользит по ЛИНЗАМ (lens_device_types,
        # дефолт QP,DQ), lens_window_size линз в окне, в начало окна инжектируется
        # свежий прямой пучок ±inject_spread_m. Активируется и при method='lens'.
        lines.append(f"SCALAR PER_ELEM_WINDOW_UNIT          "
                     f"{pe.get('window_unit', 'magnet')}")
        lines.append(f"SCALAR PER_ELEM_LENS_WINDOW_SIZE     "
                     f"{int(pe.get('lens_window_size', 4))}")
        lines.append(f"SCALAR PER_ELEM_INJECT_SPREAD_M      "
                     f"{float(pe.get('inject_spread_m', 0.001)):.10g}")
        lens_device_types = pe.get('lens_device_types') or cfg.get('lens_device_types') or []
        if lens_device_types:
            lines.append(f"STRING LENS_DEVICE_TYPES {','.join(lens_device_types)}")
        # search_mode: МЕТОД поиска в окне. ДЕФОЛТ 'line_search' (последовательное
        # приближение -- основной метод). 'monte_carlo' (алиас 'random'; случайные
        # кики/свипы) -- опциональная добавка (задать явно).
        lines.append(f"SCALAR PER_ELEM_SEARCH_MODE          "
                     f"{pe.get('search_mode', 'line_search')}")
        # line_search_points: K точек на раунд (>=3). line_search_max_rounds: потолок
        # раундов зума на элемент за свип. current_range/coeff_range: |ток|<=range А,
        # |coeff|<=range -- весь допустимый промежуток, по которому идёт поиск.
        lines.append(f"SCALAR PER_ELEM_LS_POINTS            "
                     f"{int(pe.get('line_search_points', 3))}")
        lines.append(f"SCALAR PER_ELEM_LS_MAX_ROUNDS        "
                     f"{int(pe.get('line_search_max_rounds', 12))}")
        # line_search_min_rounds: ГАРАНТИРОВАННЫЙ нижний порог раундов зума на элемент
        # за визит (даже при дефиците бюджета). По умолчанию 1 (старое поведение).
        lines.append(f"SCALAR PER_ELEM_LS_MIN_ROUNDS        "
                     f"{int(pe.get('line_search_min_rounds', 1))}")
        lines.append(f"SCALAR PER_ELEM_CURRENT_RANGE        "
                     f"{float(pe.get('current_range', 100.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_COEFF_RANGE          "
                     f"{float(pe.get('coeff_range', 3.0)):.10g}")
        # line_search_local_current/coeff: если >0, line-search ищет в ОКРЕСТНОСТИ
        # вокруг текущего значения ([cur±radius], зажатой в диапазон), а не во всём
        # диапазоне. Лекарство от railing'а: токи двигаются инкрементально. 0=весь.
        lines.append(f"SCALAR PER_ELEM_LS_LOCAL_CURRENT     "
                     f"{float(pe.get('line_search_local_current', 0.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_LS_LOCAL_COEFF       "
                     f"{float(pe.get('line_search_local_coeff', 0.0)):.10g}")
        # line_search_kick_every_sweeps/amplitude: кик в глоб. line-search — каждые N
        # свипов встряхнуть все группы на ±amplitude×|значение| (0.01=1%), потом доводка.
        # 0 свипов = выкл (как раньше).
        lines.append(f"SCALAR PER_ELEM_LS_KICK_EVERY_SWEEPS "
                     f"{int(pe.get('line_search_kick_every_sweeps', 0))}")
        lines.append(f"SCALAR PER_ELEM_LS_KICK_AMPLITUDE    "
                     f"{float(pe.get('line_search_kick_amplitude', 0.0)):.10g}")
        # line_search_block_size: блочный спуск — n групп тюнятся вместе (градиент+LS
        # вдоль -grad). 1 = по одной группе (как раньше).
        lines.append(f"SCALAR PER_ELEM_LS_BLOCK_SIZE        "
                     f"{int(pe.get('line_search_block_size', 1))}")
        # line_search_accept_final: 1 -> сохранить ФИНАЛЬНОЕ состояние (откат разрешён,
        # шаг назад ради двух вперёд), а не best. 0 = best (как раньше).
        lines.append(f"SCALAR PER_ELEM_LS_ACCEPT_FINAL      "
                     f"{1 if pe.get('line_search_accept_final', False) else 0}")
        # line_search_first_improve: 1 -> принять ПЕРВОЕ улучшающее значение группы (без дозума).
        lines.append(f"SCALAR PER_ELEM_LS_FIRST_IMPROVE     "
                     f"{1 if pe.get('line_search_first_improve', False) else 0}")
        # line_search_spline: 1 -> сплайн-шаг (K проб + минимум кубич. сплайна) вместо зума.
        lines.append(f"SCALAR PER_ELEM_LS_SPLINE            "
                     f"{1 if pe.get('line_search_spline', False) else 0}")
        # line_search_parabola: 1 -> параболический шаг (Брент: вершина параболы по тройке).
        lines.append(f"SCALAR PER_ELEM_LS_PARABOLA          "
                     f"{1 if pe.get('line_search_parabola', False) else 0}")
        # drift_extension_m: >0 -> после оконной укладки SP1 фаза "вползания" в дрейф
        # по drift_chunk_m метров (до drift_extension_m). drift_iters_per_chunk -- бюджет
        # оценок line-search на рубеж.
        lines.append(f"SCALAR PER_ELEM_DRIFT_EXTENSION_M    "
                     f"{float(pe.get('drift_extension_m', 0.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_DRIFT_CHUNK_M        "
                     f"{float(pe.get('drift_chunk_m', 2.0)):.10g}")
        lines.append(f"SCALAR PER_ELEM_DRIFT_ITERS_PER_CHUNK "
                     f"{int(pe.get('drift_iters_per_chunk', 480))}")

    # Параметры градиентного режима (валидны только при method == 'gradient').
    g = cfg.get("gradient") or {}
    if g:
        lines.append("")
        lines.append("# === GRADIENT MODE ===")
        lines.append("# FD: 'central' (точнее, 2M вызовов) или 'forward' (M+1, шумнее).")
        lines.append("# update_rule: 'normalized_gd' (Δθ = -lr·g/max|g|) или 'adam'.")
        lines.append(f"SCALAR GRAD_FD_METHOD       {g.get('fd_method', 'central')}")
        lines.append(f"SCALAR GRAD_UPDATE_RULE     {g.get('update_rule', 'normalized_gd')}")
        lines.append(f"SCALAR GRAD_N_ITER          {int(g.get('n_iterations', 500))}")
        lines.append(f"SCALAR GRAD_LEARNING_RATE   {float(g.get('learning_rate', 0.05)):.10g}")
        lines.append(f"SCALAR GRAD_BETA1           {float(g.get('adam_beta1', 0.9)):.10g}")
        lines.append(f"SCALAR GRAD_BETA2           {float(g.get('adam_beta2', 0.999)):.10g}")
        lines.append(f"SCALAR GRAD_EPS             {float(g.get('adam_eps', 1e-8)):.10g}")
        lines.append(f"SCALAR GRAD_FD_STEP_ABS     {float(g.get('fd_step_abs', 0.01)):.10g}")
        lines.append(f"SCALAR GRAD_FD_STEP_PERCENT {float(g.get('fd_step_percent', 0.001)):.10g}")
        lines.append(f"SCALAR GRAD_MAX_HALVINGS    {int(g.get('max_halvings', 5))}")

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
    from _common import workspace_dir
    ws = workspace_dir()
    if len(sys.argv) < 2:
        cfg_path = ws / "optim_inputs" / "sp1_endpoint.json"
    else:
        cfg_path = Path(sys.argv[1])
    cfg_path = cfg_path.resolve()
    out_path = cfg_path.with_suffix(".optimconfig")
    history_output = str(cfg_path.with_suffix(".history.json").name)
    with open(cfg_path, encoding="utf-8") as f:
        cfg = json.load(f)
    # sp_elements.csv ищем в workspace проекта (через _common). Может быть
    # переопределён переменной окружения SP_ELEMENTS_PATH для нестандартных
    # раскладок.
    sp_elements_path = ws / "tables" / "sp_elements.csv"
    sp_env = os.environ.get("SP_ELEMENTS_PATH")
    if sp_env:
        sp_elements_path = Path(sp_env)
    write_flat_config(cfg, out_path, sp_elements_path, history_output, ws=ws)


if __name__ == "__main__":
    main()
