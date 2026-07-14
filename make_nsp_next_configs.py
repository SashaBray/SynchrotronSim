#!/usr/bin/env python3
"""make_nsp_next_configs.py — «N супериодов + начало SP(N+1)» конфиги (АЛЬТЕРНАТИВА 50м дрейфу).
Финиш-плоскость ставится ЗА первым магнитом следующего супериода (gid N*69+1, QP_L) — это
напрямую форсирует корректный переход SP->SP (периодичность), без пустого дрейфа и без
дрейф-целей. Для каждого N=1..10: суб-lattice (max_gid=N*69+1; для N=2,3,4 переиспользуются
существующие 108/109/110), строка experiments.csv, optim JSON (global, pos+par+env, веса 1/1/10,
magnet_sp_numbers=[1..N+1]). Идемпотентно. Запуск: python make_nsp_next_configs.py
"""
import json, math, copy
import design_orbit as do
from _common import workspace_dir

WS = workspace_dir()
SP = 69
# N -> (exp_id, lattice_id, is_new_lattice). N=2,3,4 переиспользуют 108/109/110.
PLAN = {1: (23, 118, True),  2: (24, 108, False), 3: (25, 109, False), 4: (26, 110, False),
        5: (27, 119, True),  6: (28, 120, True),  7: (29, 121, True),
        8: (30, 122, True),  9: (31, 123, True),  10: (32, 124, True)}

eps, _ = do.compute_magnet_endpoints(101)
by = {e.global_id: e for e in eps}
tmpl = json.loads((WS / "optim_inputs" / "sp_3sp_drift.json").read_text(encoding="utf-8"))

exp_rows, sub_entries = [], []
for N, (exp_id, lat_id, is_new) in sorted(PLAN.items()):
    g = N * SP + 1                                  # первый магнит SP(N+1)
    e = by[g]
    fx, fy, yaw = e.exit[0], e.exit[1], e.exit_yaw
    sx, sy = math.sin(yaw), -math.cos(yaw)          # поперечная финиш-плоскости
    coords = (fx + sx, fy + sy, -1.0, fx - sx, fy - sy, -1.0, fx, fy, 1.0)
    exp_rows.append((exp_id, f"{exp_id};{lat_id};1;6,74E-13;1;100;" + ";".join(f"{v:.6f}" for v in coords)))
    if is_new:
        sub_entries.append({"lattice_id": lat_id, "max_gid": g,
            "description": f"{N} стандартных SP + первый магнит SP{N+1} (gid 1..{g}, QP_L). "
                           f"Финиш ЗА входом SP{N+1} — тест перехода SP->SP (альтернатива дрейфу). exp {exp_id}."})
    cfg = copy.deepcopy(tmpl)
    cfg["_comment"] = [f"{N} СУПЕРИОДОВ + начало SP{N+1} (lattice {lat_id}, gid 1..{g}).",
                       "Альтернатива 50м дрейфу: финиш ЗА первым магнитом следующего SP форсирует",
                       "корректный переход SP->SP (периодичность). Global, pos+parallelism+envelope_decay,",
                       "веса 1/1/10. Дрейф-цели НЕ нужны — констрейнт даёт сама финиш-плоскость на SP(N+1)."]
    cfg["experiment_id"] = exp_id
    cfg["lattice_id"] = lat_id
    cfg["lattice_max_gid"] = g
    cfg["loss"]["magnet_sp_numbers"] = list(range(1, N + 2))   # SP1..SP(N+1)
    cfg["loss"]["extra_target_magnets"] = []                   # без дрейф-целей
    (WS / "optim_inputs" / f"sp_{N}sp_next.json").write_text(
        json.dumps(cfg, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"N={N:2d}: exp {exp_id}, lattice {lat_id} (gid 1..{g}{' [reuse]' if not is_new else '       '}), "
          f"finish=({fx:.3f},{fy:.3f}) -> sp_{N}sp_next.json")

# --- experiments.csv (idempotent) ---
ep = WS / "tables" / "experiments.csv"
lines = ep.read_text(encoding="utf-8").splitlines()
have = {ln.split(";", 1)[0] for ln in lines[1:] if ln.strip()}
n = 0
for exp_id, row in exp_rows:
    if str(exp_id) not in have:
        lines.append(row); n += 1
ep.write_text("\n".join(lines) + "\n", encoding="utf-8")
print(f"experiments.csv: +{n} строк")

# --- sub_lattices.json (idempotent) ---
sl = WS / "configs" / "sub_lattices.json"
d = json.loads(sl.read_text(encoding="utf-8"))
havel = {s["lattice_id"] for s in d["sub_lattices"]}
m = 0
for ent in sub_entries:
    if ent["lattice_id"] not in havel:
        d["sub_lattices"].append(ent); m += 1
sl.write_text(json.dumps(d, indent=2, ensure_ascii=False), encoding="utf-8")
print(f"sub_lattices.json: +{m} суб-lattice")
