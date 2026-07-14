#!/usr/bin/env python3
"""continuation_tuner.py — автоматический режим-продолжатель (continuation tuner).

Автоматизирует ручной цикл настройки кольца методом продолжения:
  1. Прогон пучка по ВСЕМУ кольцу (10 SP, lattice 107) текущим тюном.
  2. Поиск рубежа: ПОСЛЕДНИЙ магнит ДО первого превышения порога (12мм) по той же
     по-магнитной метрике ближайшего подхода, что минимизирует оптимизатор => окно
     [старт..рубеж] держит max ≤ порога ПО ПОСТРОЕНИЮ.
  3. Финиш-плоскость ДИНАМИЧЕСКАЯ: finish = round(fraction × рубеж), зажат [1, конец кольца].
     fraction=1 — ровно у последнего хорошего магнита (max окна ≤ порога); <1 — отступ
     назад; >1 — намеренно ЗА точку потери (max>порога, толчок онсета).
  4. Глобальная полировка магнитов до finish_gid (lattice_max_gid), magnet_sp_numbers=[1..N].
  5. Повтор N циклов; рубеж динамический, откат рубежа допускается.
В конце — снова прогон всего кольца, отчёт, насколько далеко дошёл пучок.

Запуск:
  python continuation_tuner.py --cycles N --budget B --fraction F
                               [--threshold-mm 12] [--search line_search|monte_carlo]
"""
from __future__ import annotations
import argparse
import json
import math
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

try:   # UTF-8 вывод (Cyrillic/символы в Windows-консоли иначе падают)
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

import analyze_trajectory_stability as ats          # noqa: E402
from find_orbit_deviation_peaks import find_nearest_element  # noqa: E402
from prepare_optim_config import generate_magnet_endpoints_csv  # noqa: E402
from _common import tracker_exe, optimizer_exe, workspace_dir   # noqa: E402

WS            = workspace_dir()
FULL_RING_EXP = 10      # lattice 107, финиш на КОНЦЕ кольца (полный прогон)
POLISH_EXP    = 15      # выделенный exp для полировки (финиш на рубеже finish_gid)
RING_LAT      = 107
RING_END_GID  = 690     # последний gid 10-SP кольца
SP_PERIOD     = 69      # gid на один супериод (SP1=1..69, SP2=70..138, ...)
EXP_DT        = "6,74E-13"   # шаг интегрирования (EU-запятая, verbatim)

# Парсинг ИТОГ-блока оптимизатора (формат optimizer.cpp ~4209): "mean:  X mm  ->  Y mm"
_RE_MEAN = re.compile(r"mean:\s+([\d.]+)\s*mm\s*->\s*([\d.]+)\s*mm")
_RE_MAX  = re.compile(r"max:\s+([\d.]+)\s*mm\s*->\s*([\d.]+)\s*mm")


def gid_to_sp(gid):
    """gid -> метка 'SPs.ee' (69 gid/супериод). None -> 'всё кольцо'; 0 -> 'потерян на старте'."""
    if gid is None:
        return "всё кольцо"
    if gid == 0:
        return "потерян на старте"
    sp = (gid - 1) // SP_PERIOD + 1
    el = (gid - 1) % SP_PERIOD + 1
    return f"SP{sp}.{el:02d}"


def _gid_lbl(gid):
    """'gid N (SPs.ee, X.XX СП)' — и элементная метка, и ДРОБНЫЕ супериоды. Нужно, чтобы метку
    'SP2.26' (супериод 2, элемент 26) не путали с '2.26 СП' (на деле gid 95 = 1.38 СП)."""
    if gid is None:
        return "ВСЁ КОЛЬЦО"
    return f"gid {gid} ({gid_to_sp(gid)}, {gid / SP_PERIOD:.2f} СП)"


def _sp_elements_path() -> Path:
    return WS / "tables" / "sp_elements.csv"


def backup_sp_elements(tag) -> Path:
    """Сохранить текущий sp_elements (правило: ВСЕГДА бэкап перед оптимизатором)."""
    dst = WS / "tables" / "backups" / f"sp_elements_pre_{tag}.csv"
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(_sp_elements_path(), dst)
    return dst


def restore_sp_elements(path):
    shutil.copy2(path, _sp_elements_path())


# ----------------------------- superperiod boundaries ----------------------
def load_sp_starts():
    """gid начал супериодов (разрешённые границы финиш-плоскости) из конфига."""
    cfg = json.loads((WS / "configs" / "superperiod_starts.json").read_text(encoding="utf-8"))
    return sorted(int(g) for g in cfg["superperiod_start_gids"])


def snap_finish_to_sp(loss_gid, sp_starts):
    """Снап финиша на ГРАНИЦУ супериода: финиш = (начало SP, содержащего потерю) − 1 =
    последний gid последнего ПОЛНОГО супериода до потери. Так настраивается целое число
    периодов (loss отражает периодичность). Возвращает finish_gid (>=0)."""
    below = [s for s in sp_starts if s <= loss_gid]
    containing_start = max(below) if below else sp_starts[0]
    return containing_start - 1


# ----------------------------- experiments.csv -----------------------------
def _exp_path() -> Path:
    return WS / "tables" / "experiments.csv"


def set_finish_plane(exp_id, p1, p2, p3, direction, *, max_turns=0, stride=100):
    """Создать/обновить строку experiments.csv: финиш-плоскость (3 точки) +
    вектор направления пересечения. lattice 107, group 1, dt verbatim."""
    lines = [ln.rstrip("\n") for ln in _exp_path().read_text(encoding="utf-8").splitlines()
             if ln.strip()]
    hdr, body = lines[0], lines[1:]
    cells = [str(exp_id), str(RING_LAT), "1", EXP_DT, str(max_turns), str(stride),
             f"{p1[0]:.6f}", f"{p1[1]:.6f}", f"{p1[2]:.1f}",
             f"{p2[0]:.6f}", f"{p2[1]:.6f}", f"{p2[2]:.1f}",
             f"{p3[0]:.6f}", f"{p3[1]:.6f}", f"{p3[2]:.1f}",
             f"{direction[0]:.6f}", f"{direction[1]:.6f}", f"{direction[2]:.6f}"]
    row = ";".join(cells)
    out, found = [hdr], False
    for r in body:
        if r.split(";")[0] == str(exp_id):
            out.append(row); found = True
        else:
            out.append(r)
    if not found:
        out.append(row)
    _exp_path().write_text("\n".join(out) + "\n", encoding="utf-8", newline="\n")


# ----------------------------- geometry ------------------------------------
def load_endpoints() -> dict:
    """Полная геометрия магнитов кольца (gid -> exit/entry/yaw). Геометрия фиксирована
    (бенды — PM), поэтому считаем один раз. Перегенерируем ПОЛНЫЙ файл (без max_gid),
    т.к. полировка перезаписывает его обрезанной версией."""
    generate_magnet_endpoints_csv(RING_LAT, WS)   # full, no max_gid
    csv = WS / "results" / f"magnet_endpoints_lat{RING_LAT}.csv"
    df = pd.read_csv(csv, sep=";", skipinitialspace=True)
    d = {}
    for _, r in df.iterrows():
        d[int(r["global_id"])] = {
            "x_in": float(r["x_in"]),   "y_in": float(r["y_in"]),
            "x_out": float(r["x_out"]), "y_out": float(r["y_out"]),
            "yaw_out": float(r["yaw_out"]),
        }
    return d


def finish_plane_from_gid(endpoints, gid):
    """3 точки финиш-плоскости + вектор направления у выхода магнита gid.
    Конвенция как в experiments.csv: плоскость ⊥ касательной, нормаль = касательная."""
    e = endpoints[gid]
    x, y, yaw = e["x_out"], e["y_out"], e["yaw_out"]
    s, c = math.sin(yaw), math.cos(yaw)
    p1 = (x + s, y - c, -1.0)
    p2 = (x - s, y + c, -1.0)
    p3 = (x, y, 1.0)
    direction = (math.cos(yaw), math.sin(yaw), 0.0)
    return p1, p2, p3, direction


def snap_gid(endpoints, target):
    """Сдвинуть finish_gid к КОНЦУ co-located группы (последняя обмотка магнита,
    чтобы не обрезать DQ/корректор по середине)."""
    g = min(int(target), RING_END_GID)
    while g + 1 in endpoints \
            and abs(endpoints[g + 1]["x_in"] - endpoints[g]["x_in"]) < 1e-6 \
            and abs(endpoints[g + 1]["y_in"] - endpoints[g]["y_in"]) < 1e-6:
        g += 1
    return g


# ----------------------------- track + analyse -----------------------------
def track_full_ring(endpoints):
    """Финиш exp 10 = конец кольца + прямое направление; прогон трекером."""
    p1, p2, p3, dr = finish_plane_from_gid(endpoints, RING_END_GID)
    set_finish_plane(FULL_RING_EXP, p1, p2, p3, dr, max_turns=0)
    subprocess.run([str(tracker_exe()), str(FULL_RING_EXP)],
                   cwd=str(SCRIPT_DIR), check=True)


def _monotonic_dev(part, starts, dirs, seg_len2, fwd_window=300):
    """Отклонение от орбиты с МОНОТОННЫМ прогрессом: ближайший сегмент орбиты ищется в
    окне [last, last+fwd_window] вперёд (не прыгает назад на другую часть кривой).
    Проверено: даёт правильный рубеж на эталонах (4-SP -> 302). Малое окно НЕ годится —
    лагает и ложно завышает отклонение (W=12 давало 27 вместо 302). Метрика мягковата к
    сильно разлетевшимся состояниям (пучок может «прятаться» у соседнего сегмента), но
    строгие варианты (длина пути / малое окно) ломают хорошие состояния. Опт. метрика
    (по-магнитно) надёжнее, но не извлекается из сохранённой траектории просто.
    Возвращает (dev[], closest_pts[] на орбите)."""
    nseg = len(starts)
    dev = np.empty(len(part))
    cp = np.empty((len(part), 3))
    last = 0
    for t in range(len(part)):
        lo, hi = last, min(nseg, last + fwd_window)
        rel = part[t] - starts[lo:hi]
        tt = np.clip((rel * dirs[lo:hi]).sum(1) / seg_len2[lo:hi], 0.0, 1.0)
        c = starts[lo:hi] + tt[:, None] * dirs[lo:hi]
        d = np.linalg.norm(part[t] - c, axis=1)
        j = int(np.argmin(d))
        dev[t] = d[j]; cp[t] = c[j]; last = lo + j
    return dev, cp


def find_loss_gid(threshold_m):
    """Рубеж = ПОСЛЕДНИЙ магнит ДО первого превышения порога, по метрике ближайшего подхода
    к центру железа магнита — ТОЙ ЖЕ, что МИНИМИЗИРУЕТ оптимизатор (исключая DP_L_PM, как его
    exclude_device_types). Это даёт СТРОГУЮ гарантию: окно [старт..рубеж] держит max ≤ порога
    ПО ПОСТРОЕНИЮ (все магниты окна в пороге), значит Baseline max оптимизатора ≤ порога без
    всякого fraction<1. Финиш = round(fraction × рубеж): fraction=1 — ровно до последнего
    хорошего магнита, fraction<1 — отступ назад, fraction>1 — намеренно ЗА точку потери
    (тогда max>порога — это сознательный «толчок онсета»).

    Раньше метрика была «расстояние до ЛИНИИ орбиты»: она иногда ОТСТАЁТ от по-магнитной
    (пучок скользит вдоль линии не на своём s-месте) и завышает рубеж — из-за чего finish
    залезал в уже-потерянную зону, а окно показывало max>порога.

    Возвращает (loss_gid, max_line_dev_mm):
      loss_gid = None  -> пучок держит ВСЁ кольцо (ни один магнит не за порогом);
      loss_gid = 0     -> потерян на первом же магните (нет ни одного магнита в пороге);
      loss_gid > 0     -> последний магнит, до которого max ≤ порога.
    max_line_dev_mm — масштаб разлёта (макс. отклонение от ЛИНИИ орбиты, только для отчёта)."""
    positions, yaws, first_tan, last_tan, lat_id = ats.load_design_orbit(FULL_RING_EXP)
    poly_ext = ats.extend_polyline(ats.refine_polyline_with_arcs(positions, yaws, 10),
                                   first_tan, last_tan, 10.0, 20.0)
    starts = poly_ext[:-1]
    dirs = poly_ext[1:] - poly_ext[:-1]
    seg_len2 = np.maximum((dirs * dirs).sum(1), 1e-30)
    parts = ats.load_trajectories(WS / "results" / f"trajectories_exp{FULL_RING_EXP}.bin")
    lat_df = pd.read_csv(WS / "tables" / "lattice_configs.csv", sep=";", skipinitialspace=True)
    lat = lat_df[lat_df["lattice_id"] == lat_id].sort_values("instance_id").reset_index(drop=True)
    keep = lat["device_type_id"].astype(str).str.strip() != "DP_L_PM"   # как exclude_device_types
    cen = lat.loc[keep, ["x", "y", "z"]].to_numpy(float)
    gids = lat.loc[keep, "instance_id"].to_numpy(int)
    pts = np.vstack(parts)                                  # все точки всех частиц
    # ближайший подход к центру каждого магнита В ПОРЯДКЕ gid. Рубеж = ПОСЛЕДНИЙ магнит ДО
    # первого превышения => окно [старт..рубеж] держит max ≤ порога ПО ПОСТРОЕНИЮ.
    loss_gid = None                                        # None => держит всё кольцо
    thr2 = threshold_m * threshold_m
    for k in range(len(cen)):
        if ((pts - cen[k]) ** 2).sum(1).min() > thr2:
            loss_gid = int(gids[k - 1]) if k > 0 else 0    # 0 => потерян на первом магните
            break
    # масштаб разлёта — макс. отклонение от ЛИНИИ орбиты (только для отчёта). _monotonic_dev
    # это питон-цикл по точкам; на разлетевшихся пучках их сотни тысяч (тут — главный расход
    # времени). Для МАСШТАБА прорежаем до ~4000 точек/частицу: max меняется на доли %, а время
    # падает в ~stride раз. Рубеж считается отдельно (по-магнитно, выше) — там прореживания нет.
    max_dev = 0.0
    for part in parts:
        stride = max(1, len(part) // 4000)
        dev, _ = _monotonic_dev(part[::stride], starts, dirs, seg_len2)
        max_dev = max(max_dev, float(dev.max()))
    return loss_gid, max_dev * 1000.0


# ----------------------------- polish --------------------------------------
def polish_config(finish_gid, budget, search="line_search",
                  kick_every=0, kick_scale=0.0, n_sp=10,
                  ls_kick_sweeps=0, ls_kick_amp=0.0, ls_max_rounds=12,
                  env_w=1.0, par_w=100.0, block_size=1, accept_final=False,
                  svd_angle_w=0.0, first_improve=False, spline=False, parabola=False,
                  full_range=False) -> dict:
    return {
        "experiment_id": POLISH_EXP,
        "optimization_mode": "global",
        "lattice_id": RING_LAT,
        "lattice_max_gid": finish_gid,
        "n_iterations": budget,
        "random_seed": 9001,
        "kick_every": kick_every, "kick_scale": kick_scale, "revert_after_kicks": 0,
        "global_search_mode": search,
        "loss": {
            "type": "orbit_avg", "mode": "composite",
            "magnet_sp_numbers": list(range(1, n_sp + 1)),
            "position_reference": "orbit",
            "exclude_device_types": ["DP_L_PM"],
            "extra_target_magnets": [], "squared": True,
            "comment": f"continuation tuner: position + parallelism({par_w}) + envelope_decay({env_w})",
            "direction_loss": {"reference": [1.0, 0.0, 0.0], "weight": 0.0},
            "envelope_loss": {"weight": 0.0, "polyline_mode": "arcs", "arc_subdivide": 10},
            "max_weight": 0.0, "interparticle_weight": 0.0,
            "envelope_decay_weight": env_w, "parallelism_weight": par_w,
            "svd_angle_weight": svd_angle_w,
        },
        "parameter_groups": [],
        "parameter_groups_per_position": [
            {"enabled": True, "file": "sp_elements.csv", "sp_types": ["standard"],
             "exclude_device_types": ["DP_L_PM", "DP_sh_PM", "SP_I2_1A", "DQ_I3_1A",
                                      "corrector_I1_1A", "corrector_I2_1A"],
             "value_column": "arg_val", "step_abs": 0.05},
            {"enabled": True, "file": "sp_elements.csv", "sp_types": ["standard"],
             "include_device_types": ["DP_L_PM", "DP_sh_PM"],
             "value_column": "correct_coeff", "step_percent": 0.025, "name_suffix": "_coeff"},
        ],
        "plot_every": 1,
        "per_element": {
            "loss_mode": "composite", "line_search_points": 3, "line_search_max_rounds": ls_max_rounds,
            "current_range": 100.0, "coeff_range": 3.0,
            "line_search_local_current": 0.0 if full_range else 10.0,
            "line_search_local_coeff": 0.0 if full_range else 0.3,
            "line_search_kick_every_sweeps": ls_kick_sweeps,
            "line_search_kick_amplitude": ls_kick_amp,
            "line_search_block_size": block_size,
            "line_search_accept_final": accept_final,
            "line_search_first_improve": first_improve,
            "line_search_spline": spline,
            "line_search_parabola": parabola,
            "total_iters": budget,
        },
        "dynamic_loss": {"enabled": False, "amplitude": 0.7, "period_iterations": 0},
        "max_particles": 0,
    }


def run_polish(endpoints, finish_gid, n_sp, budget, search="line_search",
               kick_every=0, kick_scale=0.0, ls_kick_sweeps=0, ls_kick_amp=0.0, ls_max_rounds=12,
               env_w=1.0, par_w=100.0, block_size=1, accept_final=False, svd_angle_w=0.0,
               first_improve=False, spline=False, parabola=False, full_range=False):
    """Финиш exp 15 = finish_gid (n_sp супериодов); конфиг lattice_max_gid=finish_gid;
    запуск оптимизатора.

    Стримит ВЕСЬ вывод оптимизатора (по-элементные строки) и парсит mean/max до→после.
    search: line_search | monte_carlo (с киками kick_every/kick_scale).
    ls_kick_sweeps/amp: кик в line-search каждые N свипов на ±amp×|value|."""
    p1, p2, p3, dr = finish_plane_from_gid(endpoints, finish_gid)   # finish_gid уже на границе SP
    set_finish_plane(POLISH_EXP, p1, p2, p3, dr, max_turns=0)
    cfg_path = WS / "optim_inputs" / "continuation_polish.json"
    cfg_path.write_text(json.dumps(polish_config(finish_gid, budget, search, kick_every, kick_scale,
                                                  n_sp, ls_kick_sweeps, ls_kick_amp, ls_max_rounds,
                                                  env_w, par_w, block_size, accept_final, svd_angle_w,
                                                  first_improve, spline, parabola, full_range),
                                   indent=2, ensure_ascii=False), encoding="utf-8")
    stats = {"finish_gid": finish_gid, "mean_before": None, "mean_after": None,
             "max_before": None, "max_after": None}
    proc = subprocess.Popen(
        [sys.executable, "optimize.py", str(cfg_path), "--method=global", f"--budget={budget}"],
        cwd=str(SCRIPT_DIR), stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        encoding="utf-8", errors="replace", bufsize=1)
    for line in proc.stdout:
        line = line.rstrip("\n")
        print(line, flush=True)                    # пасс-throughput: весь вывод оптимизатора,
        m = _RE_MEAN.search(line)                  # включая по-элементные [LS ...] строки
        if m:
            stats["mean_before"], stats["mean_after"] = float(m.group(1)), float(m.group(2))
        m = _RE_MAX.search(line)
        if m:
            stats["max_before"], stats["max_after"] = float(m.group(1)), float(m.group(2))
    rc = proc.wait()
    if rc != 0:
        print(f"  ! оптимизатор завершился с кодом {rc}", flush=True)
    return stats


# ----------------------------- main loop -----------------------------------
def _assemble_track_reach(endpoints, thr):
    """Пересобрать lattice_configs.csv, прогнать всё кольцо, вернуть (loss_gid|None, max_dev_mm)."""
    subprocess.run([sys.executable, "assemble_lattice_configs.py"],
                   cwd=str(SCRIPT_DIR), check=True)
    track_full_ring(endpoints)
    return find_loss_gid(thr)


def _cap_reach(reach, cap_gid):
    """Применить ПОТОЛОК рубежа (cap_gid в gid, None=без потолка). reach=None (держит всё
    кольцо) трактуется как конец кольца. Возвращает эффективный рубеж в gid для расчёта финиша."""
    eff = RING_END_GID if reach is None else reach
    return min(eff, cap_gid) if cap_gid is not None else eff


def _fmt_adv(adv):
    """Строка продвижения рубежа: '+12' / '-3' / '0' / '—'."""
    if adv is None:
        return "—"
    return f"+{adv}" if adv >= 0 else str(adv)


def print_cycle_summary(cyc, reach_before, dev_before, reach_post, dev_post, stats, dt):
    """Сводка по одному циклу: что было, как отполировали, что стало (reach_post —
    динамический замер ПОСЛЕ полировки, принимается как есть, в т.ч. откат)."""
    adv = (reach_post - reach_before) if (reach_post is not None and reach_before is not None) else None
    reach_str = _gid_lbl(reach_post)
    print(f"  ----- СВОДКА ЦИКЛА {cyc} -----", flush=True)
    if stats["mean_after"] is not None and stats["max_after"] is not None:
        print(f"    полировка участка: mean {stats['mean_before']:.3f} -> {stats['mean_after']:.3f} мм, "
              f"max {stats['max_before']:.3f} -> {stats['max_after']:.3f} мм", flush=True)
    else:
        print("    полировка участка: метрики ИТОГ не распознаны (см. лог выше)", flush=True)
    print(f"    рубеж потери:   {_gid_lbl(reach_before)} -> {reach_str}"
          f"  ({_fmt_adv(adv)} gid)", flush=True)
    print(f"    max откл. по кольцу: {dev_before:.1f} -> {dev_post:.1f} мм", flush=True)
    print(f"    финиш-плоскость: {_gid_lbl(stats['finish_gid'])}", flush=True)
    print(f"    время цикла: {dt:.0f} с ({dt/60:.1f} мин)", flush=True)
    if adv is not None and adv < 0:
        print("    (рубеж откатился — принято, динамический режим)", flush=True)


def print_final_summary(rows, start_reach, start_dev, reach, max_dev, total_sec):
    """Большая итоговая сводка: таблица по циклам + общее продвижение."""
    print("\n==================== ИТОГОВАЯ СВОДКА ====================", flush=True)
    print(f"  циклов выполнено: {len(rows)}    общее время: {total_sec/60:.1f} мин", flush=True)
    print(f"  старт: рубеж gid {start_reach} ({gid_to_sp(start_reach)}), max откл. {start_dev:.1f} мм",
          flush=True)
    if rows:
        print(f"  {'цикл':>4} | {'рубеж до':>11} | {'финиш':>10} | {'mean,мм':>8} | "
              f"{'max,мм':>7} | {'рубеж после':>13} | {'+gid':>6}", flush=True)
        print("  " + "-" * 80, flush=True)
        for r in rows:
            after_s = gid_to_sp(r["after"]) if r["after"] is not None else "всё кольцо"
            ma = f"{r['mean_after']:.2f}" if r["mean_after"] is not None else "—"
            mx = f"{r['max_after']:.2f}" if r["max_after"] is not None else "—"
            print(f"  {r['cyc']:>4} | {gid_to_sp(r['before']):>11} | {gid_to_sp(r['finish']):>10} | "
                  f"{ma:>8} | {mx:>7} | {after_s:>13} | {_fmt_adv(r['adv']):>6}", flush=True)
    final_str = f"gid {reach} ({gid_to_sp(reach)})" if reach is not None else "ВСЁ КОЛЬЦО (10 SP)"
    total_adv = (reach - start_reach) if (reach is not None and start_reach is not None) else None
    print("  " + "-" * 80, flush=True)
    tail = ""
    if total_adv is not None:
        tail = (f"  (старт gid {start_reach}, продвижение {_fmt_adv(total_adv)} gid "
                f"~ {total_adv/SP_PERIOD:.1f} SP)")
    print(f"  ИТОГ: пучок доходит до {final_str}, max откл. {max_dev:.1f} мм{tail}", flush=True)
    if reach is None:
        print("  [OK] Кольцо настроено целиком!", flush=True)
    print("  =======================================================", flush=True)


def main():
    ap = argparse.ArgumentParser(description="Continuation tuner: ползущий рубеж по кольцу.")
    ap.add_argument("--cycles", type=int, required=True, help="число циклов (поиск+полировка)")
    ap.add_argument("--budget", type=int, required=True, help="бюджет симуляций на полировку")
    ap.add_argument("--fraction", type=float, required=True,
                    help="finish = round(доля × рубеж); МОЖЕТ быть >1 (финиш ЗА точкой потери), "
                         "зажат концом кольца (gid 690)")
    ap.add_argument("--max-sp", type=float, default=None,
                    help="ПОТОЛОК рубежа в супериодах (можно дробный, напр. 2.1): программа не "
                         "сдвигает финиш дальше этого, даже если пучок держит (max<порог) за ним. "
                         "По умолчанию без потолка (рубеж до конца кольца).")
    ap.add_argument("--threshold-mm", type=float, default=12.0, help="порог потери (мм)")
    ap.add_argument("--search", choices=["line_search", "monte_carlo", "svd"], default="line_search",
                    help="метод глобальной полировки: line_search (покоординатно) | "
                         "monte_carlo (кики+hill-climber) | svd (матрица отклика + Гаусс-Ньютон, "
                         "коррекция орбиты — все ручки сразу)")
    ap.add_argument("--kick-every", type=int, default=0, help="MC: кик каждые N итераций (0=выкл, по умолчанию)")
    ap.add_argument("--kick-scale", type=float, default=25.0,
                    help="MC: масштаб кика (×шаг группы; 25 ~ ±1.25A на ток)")
    ap.add_argument("--ls-kick-every-sweeps", type=int, default=0,
                    help="line-search: кик каждые N свипов (0=выкл; 2 = каждый 2-й свип)")
    ap.add_argument("--ls-kick-amplitude", type=float, default=0.01,
                    help="line-search: амплитуда кика, доля от |значения| (0.01 = 1%%)")
    ap.add_argument("--ls-max-rounds", type=int, default=12,
                    help="раундов зума на группу за свип (меньше = дешевле свип = больше свипов/киков)")
    ap.add_argument("--envelope-weight", type=float, default=1.0,
                    help="вес НАКЛОНА огибающей (новый член ~метры наклона; масштаб ~×100-1000 "
                         "больше старого, потому дефолт 1, не 1000)")
    ap.add_argument("--parallelism-weight", type=float, default=100.0,
                    help="вес параллельности parallelism в loss")
    ap.add_argument("--svd-angle-weight", type=float, default=0.0,
                    help="вес ПЕРИОДИЧЕСКОГО члена: выход ЦЕНТРОИДА супериода ∥ орбите "
                         "(замкнутая орбита). Работает в svd (в невязке) И в line_search "
                         "(в loss — устойчив к джиттеру). 0=выкл. ~1-5 — рычаг против барьера "
                         "периодичности. С line_search обычно эффективнее: коарс-пробинг видит "
                         "малый угловой градиент выше дискретизац. шума")
    ap.add_argument("--block-size", type=int, default=1,
                    help="блочный спуск: n групп тюнятся ВМЕСТЕ (градиент+LS вдоль -grad); "
                         "1 = по одной (как раньше), 4-5 — против координатных застреваний")
    ap.add_argument("--ls-first-improvement", action="store_true",
                    help="line_search: принимать ПЕРВОЕ значение группы, улучшившее тюн, без "
                         "дозума к минимуму (экономит бюджет; грубее => реже застревает в "
                         "остром экстремуме)")
    ap.add_argument("--ls-spline", action="store_true",
                    help="line_search: СПЛАЙН-шаг — K(>=5) проб в окрестности, кубический сплайн "
                         "loss(значение) и его минимум на отрезке как кандидат (сглаживает "
                         "дискретизац. шум, находит минимум МЕЖДУ пробами). Кандидат проверяется "
                         "симуляцией. Несовместим с --ls-first-improvement (сплайн имеет приоритет)")
    ap.add_argument("--ls-parabola", action="store_true",
                    help="line_search: ПАРАБОЛА (МНК по всем K>=5 пробам) — вершина least-squares "
                         "параболы как кандидат (с защитой: минимум и в диапазоне, иначе откат на "
                         "пробу). Усредняет шум, не ныряет в минус как сплайн. Приоритет над "
                         "--ls-spline. Точки печатаются в лог [PARAB ...] для анализа")
    ap.add_argument("--ls-full-range", action="store_true",
                    help="line_search: искать/аппроксимировать loss во ВСЁМ ОДЗ (|ток|<=100A, "
                         "|coeff|<=3), а не в окрестности ±10A. С параболой/сплайном число проб "
                         "автоматически растёт (11). Глобальный покоординатный поиск — против "
                         "локальных экстремумов. Узкую яму грубая сетка может промахнуть (зум надёжнее)")
    ap.add_argument("--ramp-stages", type=int, default=0,
                    help="ПРОГРЕССИВНЫЙ рубеж: прежде чем тюнить на полный finish, прогнать N "
                         "стадий с финишем k/N×finish (k=1..N) — окно ползёт от начала к рубежу, "
                         "ранние участки кондиционируются до расширения. Каждая стадия = полная "
                         "полировка (бюджет --budget). 0=выкл (сразу полный finish). Стоимость ~N×")
    ap.add_argument("--allow-rollback", action="store_true",
                    help="сохранять ФИНАЛЬНОЕ состояние (даже хуже best) — откат разрешён "
                         "(шаг назад ради двух вперёд); кап безопасности от разлёта")
    args = ap.parse_args()
    thr = args.threshold_mm / 1000.0
    ke = args.kick_every if args.search == "monte_carlo" else 0
    ks = args.kick_scale if args.search == "monte_carlo" else 0.0
    cap_gid = round(args.max_sp * SP_PERIOD) if args.max_sp else None   # потолок рубежа (gid)

    t_start = time.time()
    endpoints = load_endpoints()

    # Стартовый рубеж — где пучок теряется ДО настройки.
    reach, max_dev = _assemble_track_reach(endpoints, thr)
    start_reach, start_dev = reach, max_dev
    print(f"\n>>> СТАРТ: пучок доходит до {_gid_lbl(reach)}, "
          f"max откл. {max_dev:.1f} мм", flush=True)
    if cap_gid is not None:
        print(f">>> ПОТОЛОК рубежа: {args.max_sp:.2f} СП = gid {cap_gid} "
              f"— финиш не пойдёт дальше этого gid", flush=True)

    rows = []

    for cyc in range(1, args.cycles + 1):
        if reach is None and cap_gid is None:
            print("Пучок держит ВСЁ кольцо — настраивать больше нечего.", flush=True)
            break
        t_cyc = time.time()
        reach_before, dev_before = reach, max_dev
        # Эффективный рубеж = реальный, зажатый ПОТОЛКОМ (--max-sp). reach=None (держит всё
        # кольцо) -> конец кольца, затем потолок. Финиш = доля × эфф.рубеж. fraction>1 ставит
        # финиш ЗА точкой потери (зажат концом кольца). Откат рубежа допускается каждый раунд.
        reach_eff = _cap_reach(reach_before, cap_gid)
        finish_gid = round(args.fraction * reach_eff)
        if cap_gid is not None:
            finish_gid = min(finish_gid, cap_gid)   # ПОТОЛОК жёстко ограничивает финиш (в т.ч. fraction>1)
        finish_gid = max(1, min(finish_gid, RING_END_GID))
        capped = cap_gid is not None and finish_gid >= cap_gid    # финиш уперся в потолок
        n_sp = min(10, (finish_gid - 1) // SP_PERIOD + 1)        # для magnet_sp_numbers
        past = "  ЗА точкой потери!" if (reach_before is not None and finish_gid > reach_before) else ""
        cap_note = f"  [потолок {args.max_sp} SP]" if capped else ""

        print(f"\n========== ЦИКЛ {cyc}/{args.cycles} ==========", flush=True)
        print(f"  рубеж до:  {_gid_lbl(reach_before)}, "
              f"max откл. {dev_before:.1f} мм", flush=True)
        print(f"  финиш:     {_gid_lbl(finish_gid)} "
              f"= round({args.fraction}*{reach_eff}){past}{cap_note}, ~{n_sp} SP, "
              f"бюджет {args.budget}, метод {args.search}"
              + (f" (кики {ke}×{ks})" if (args.search == "monte_carlo" and ke > 0) else "")
              + (f" +LS-кик/{args.ls_kick_every_sweeps}св×{args.ls_kick_amplitude*100:.0f}%"
                 if args.ls_kick_every_sweeps > 0 else "")
              + (f" +блок×{args.block_size}" if args.block_size > 1 else ""), flush=True)

        backup_sp_elements(f"cyc{cyc}")                   # бэкап ПЕРЕД полировкой (всегда)
        # ПРОГРЕССИВНЫЙ рубеж (--ramp-stages N): финиш ползёт k/N×finish, k=1..N. Ранние
        # (узкие) окна кондиционируются ДО расширения — гомотопия по окну. N<=1 => как раньше.
        stages = args.ramp_stages if args.ramp_stages > 0 else 1
        for st in range(1, stages + 1):
            stage_finish = finish_gid if st == stages else max(1, round(st * finish_gid / stages))
            stage_nsp = min(10, (stage_finish - 1) // SP_PERIOD + 1)
            if stages > 1:
                print(f"  [рамп {st}/{stages}] финиш {_gid_lbl(stage_finish)}", flush=True)
            stats = run_polish(endpoints, stage_finish, stage_nsp, args.budget, args.search, ke, ks,
                               args.ls_kick_every_sweeps, args.ls_kick_amplitude, args.ls_max_rounds,
                               args.envelope_weight, args.parallelism_weight, args.block_size,
                               args.allow_rollback, args.svd_angle_weight, args.ls_first_improvement,
                               args.ls_spline, args.ls_parabola, args.ls_full_range)

        # Пере-замер. Рубеж принимается КАК ЕСТЬ, в т.ч. откат (динамический режим).
        reach, max_dev = _assemble_track_reach(endpoints, thr)
        dt = time.time() - t_cyc

        print_cycle_summary(cyc, reach_before, dev_before, reach, max_dev, stats, dt)
        rows.append({"cyc": cyc, "before": reach_before, "finish": stats["finish_gid"],
                     "after": reach,
                     "adv": (reach - reach_before)
                            if (reach is not None and reach_before is not None) else None,
                     "mean_after": stats["mean_after"], "max_after": stats["max_after"],
                     "dev_after": max_dev, "dt": dt})

    print_final_summary(rows, start_reach, start_dev, reach, max_dev, time.time() - t_start)


if __name__ == "__main__":
    main()
