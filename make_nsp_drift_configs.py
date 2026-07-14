#!/usr/bin/env python3
"""make_nsp_drift_configs.py — генерирует «N супериодов + 50м drift» конфиги по образцу
config 8 (3 SP + drift). Для каждого N:
  - суб-lattice (max_gid = N*69) в sub_lattices.json (N=10 переиспользует lattice 107);
  - строка в experiments.csv с финиш-плоскостью 50м ЗА выходом N-го SP (перпенд. касательной);
  - optim JSON (global, pos-loss, 5 целей вдоль дрейфа через 10м, weight 30).
Метод проверен: на config 8 даёт ровно его финиш и цели. Идемпотентно (существующие exp/lattice
пропускаются). Запуск: python make_nsp_drift_configs.py
"""
import json, math, copy
from pathlib import Path
import design_orbit as do
from _common import workspace_dir

WS = workspace_dir()
SP, DRIFT_M, N_TARGETS = 69, 50.0, 5
# N -> (exp_id, lattice_id). N=10 переиспользует lattice 107 (max_gid 690 уже есть).
PLAN = {5: (16, 112), 6: (17, 113), 7: (18, 114), 8: (19, 115), 9: (20, 116), 10: (21, 107)}

eps, _ = do.compute_magnet_endpoints(107)
by = {e.global_id: e for e in eps}
tmpl = json.loads((WS / "optim_inputs" / "sp_3sp_drift.json").read_text(encoding="utf-8"))

exp_rows, sub_entries = [], []
for N, (exp_id, lat_id) in sorted(PLAN.items()):
    g = N * SP
    e = by[g]
    x, y, yaw = e.exit[0], e.exit[1], e.exit_yaw
    cx, cy = math.cos(yaw), math.sin(yaw)            # касательная орбиты (направление дрейфа)
    sx, sy = math.sin(yaw), -math.cos(yaw)           # поперечная (размах финиш-плоскости)
    targets = [(x + k * 10 * cx, y + k * 10 * cy) for k in range(1, N_TARGETS + 1)]
    fx, fy = targets[-1]                              # финиш = точка 50м
    # финиш-плоскость (15 колонок, как config 8): p1=f+perp, p2=f-perp (z=-1), p3=f (z=+1)
    coords = (fx + sx, fy + sy, -1.0, fx - sx, fy - sy, -1.0, fx, fy, 1.0)
    row = f"{exp_id};{lat_id};1;6,74E-13;1;100;" + ";".join(f"{v:.6f}" for v in coords)
    exp_rows.append((exp_id, row))
    if lat_id != 107:
        sub_entries.append({"lattice_id": lat_id, "max_gid": g,
                            "description": f"{N} стандартных супериодов (gid 1..{g}) + 50м drift. "
                                           f"Финиш-плоскость 50м за выходом SP{N}. exp {exp_id}."})
    # optim JSON
    cfg = copy.deepcopy(tmpl)
    cfg["_comment"] = [f"{N} СУПЕРИОДОВ (lattice {lat_id}, gid 1..{g}) + 50м drift.",
                       "Глобальная оптимизация (line_search по умолчанию), pos-loss (dir=0, env=0).",
                       f"loss по магнитам {N} SP + {N_TARGETS} целей вдоль drift (шаг 10м, weight 30).",
                       "Гипотеза: длинный дрейф = усилитель непараллельности/непериодичности на выходе N SP."]
    cfg["experiment_id"] = exp_id
    cfg["lattice_id"] = lat_id
    cfg["lattice_max_gid"] = g
    cfg["loss"]["magnet_sp_numbers"] = list(range(1, N + 1))
    cfg["loss"]["extra_target_magnets"] = (
        [{"_doc": f"{N_TARGETS} целей вдоль 50м drift после {N} SP, шаг 10м, weight 30"}]
        + [{"x": round(tx, 6), "y": round(ty, 6), "z": 0.0, "weight": 30.0} for tx, ty in targets])
    (WS / "optim_inputs" / f"sp_{N}sp_drift.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"N={N:2d}: exp {exp_id}, lattice {lat_id} (gid 1..{g}, yaw={math.degrees(yaw):.0f}°), "
          f"finish=({fx:.3f},{fy:.3f}) -> sp_{N}sp_drift.json")

# --- append experiments.csv (idempotent) ---
exp_path = WS / "tables" / "experiments.csv"
lines = exp_path.read_text(encoding="utf-8").splitlines()
existing = {ln.split(";", 1)[0] for ln in lines[1:] if ln.strip()}
added = 0
for exp_id, row in exp_rows:
    if str(exp_id) not in existing:
        lines.append(row); added += 1
exp_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"experiments.csv: +{added} строк")

# --- add sub_lattices.json (idempotent) ---
sl_path = WS / "configs" / "sub_lattices.json"
sl = json.loads(sl_path.read_text(encoding="utf-8"))
have = {s["lattice_id"] for s in sl["sub_lattices"]}
addl = 0
for ent in sub_entries:
    if ent["lattice_id"] not in have:
        sl["sub_lattices"].append(ent); addl += 1
sl_path.write_text(json.dumps(sl, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"sub_lattices.json: +{addl} суб-lattice")
