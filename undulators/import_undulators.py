# -*- coding: utf-8 -*-
"""Импорт составной решётки с ондуляторами в БД + фокус по инстансу + реестр жёстких групп.

ШАГИ
  1. s2lattice_to_db.run — плоский .lte → БД как решётка LAT (устройства `und_*`, геометрия проходом
     кольца). База переимпортируется вместе с ондуляторами ⇒ решётка самодостаточна, выверенную реш.1
     не трогает (сравнение «вкл/выкл» делается ВНУТРИ этой решётки).
  2. ПАТЧ по ИНСТАНСУ (носитель — строка lattice_configs): фокус-элемент UFC (kind=quad → body_kick)
     получает mat_Gv (вертикальный k1v — естеств. фокус k_y) и mat_G (горизонт. k1 — гасит побочную
     1/ρ² секторных полюсов) ⇒ нетто: горизонталь дрейф, вертикаль k_y. Полюса остаются dipole (изгиб-
     орбита + увод при перекосе), их mat_B×=k по зазору. k1/k1v в body_kick ядра НЕЗАВИСИМЫ (matrix_core.h).
     Полюса одного устройства — непрерывный прогон инстансов типов und_UPP/UPM/UE0/UE1; прогоны ↔ слоты
     по порядку s. Параметры (k, k1v) — из спутника undulator_knobs.csv.
  3. element_groups.csv — реестр ЖЁСТКИХ групп: group_id → инстансы + опорная точка (центр устройства) и
     рамка. На него опирается errors.py (жёсткий сдвиг/поворот всего ондулятора как целого).

ЗАПУСК
  py -3.14 undulators/import_undulators.py --config configuration_v6_sila2026_und --lattice 26
"""
import os, sys, csv, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import s2lattice_to_db as s2
import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
LC = os.path.join(TB, "lattice_configs.csv")
ET = os.path.join(TB, "element_transforms.csv")
POLE_DEVS = None                                      # заполняется префиксом: {und_UPP, und_UPM, und_UE0, und_UE1}


def _read(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    return rows[0], rows[1:]


def _write(path, hdr, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(rows)


def patch_focus(lat, prefix, knobs, manifest):
    """Патч: фокус-элемент UFC (kind=quad) ← mat_G/mat_Gv устройства; полюса (dipole) ← mat_B×k.
    Прогон = непрерывные инстансы полюсов+UFC; прогоны ↔ слоты по порядку s. Группа = весь прогон."""
    pole_devs = {prefix + n for n in ("UPP", "UPM", "UE0", "UE1")}
    fc_dev = prefix + "UFC"
    member_devs = pole_devs | {fc_dev}
    hdr, rows = _read(LC)
    ci = {c: i for i, c in enumerate(hdr)}
    li, di, ii = ci["lattice_id"], ci["device_type_id"], ci["instance_id"]

    s_of = {r["slot_id"]: 0.5 * (float(r["s_start"]) + float(r["s_end"])) for r in manifest}
    kn = sorted(knobs, key=lambda r: s_of.get(r["slot_id"], 0.0))   # слоты по порядку s

    runs, cur = [], []                                 # непрерывные прогоны (полюса+UFC) → устройства
    for r in rows:
        if r[li] != str(lat):
            continue
        if r[di] in member_devs:
            cur.append(r)
        elif cur:
            runs.append(cur); cur = []
    if cur:
        runs.append(cur)
    if len(runs) != len(kn):
        raise SystemExit(f"[ОТКАЗ] прогонов ондуляторов {len(runs)} ≠ слотов в спутнике {len(kn)}")

    groups = []
    for run, knb in zip(runs, kn):
        k = float(knb["k"]); fc_G = float(knb["fc_mat_G"]); fc_Gv = float(knb["fc_mat_Gv"])
        inst = []; xs = ys = zs = 0.0
        for r in run:
            if r[di] == fc_dev:                        # фокус-элемент: вертикаль k1v + гашение горизонт. k1
                r[ci["mat_G"]] = f"{fc_G:.15g}"; r[ci["mat_Gv"]] = f"{fc_Gv:.15g}"
            elif k != 1.0 and r[ci["mat_B"]] not in ("", "0"):   # полюса: сила изгиба по зазору
                r[ci["mat_B"]] = f"{float(r[ci['mat_B']]) * k:.15g}"
            inst.append(r[ii]); xs += float(r[ci["x"]]); ys += float(r[ci["y"]]); zs += float(r[ci["z"]])
        n = len(run); mid = run[n // 2]
        groups.append(dict(slot=knb["slot_id"], lat=str(lat), insts=inst,
                           px=xs / n, py=ys / n, pz=zs / n,
                           yaw=mid[ci["yaw"]], pitch=mid[ci["pitch"]], roll=mid[ci["roll"]]))
    _write(LC, hdr, rows)
    return groups


def write_groups(groups):
    """element_groups.csv — реестр жёстких групп (один ондулятор = одно тело)."""
    path = os.path.join(TB, "element_groups.csv")
    hdr = ["group_id", "lattice_id", "kind", "pivot_x", "pivot_y", "pivot_z", "yaw", "pitch", "roll", "instance_ids"]
    # сохранить чужие группы (другие решётки), переписать свои
    keep = []
    if os.path.exists(path):
        old_hdr, old = _read(path)
        lats = {g["lat"] for g in groups}
        keep = [r for r in old if r and r[old_hdr.index("lattice_id")] not in lats]
    rows = []
    for i, g in enumerate(sorted(groups, key=lambda g: g["slot"]), 1):
        rows.append([f"und_{i:02d}", g["lat"], "undulator", f"{g['px']:.12f}", f"{g['py']:.12f}",
                     f"{g['pz']:.12f}", g["yaw"], g["pitch"], g["roll"], "|".join(g["insts"])])
    _write(path, hdr, keep + rows)
    return path, len(rows)


def main():
    ap = argparse.ArgumentParser(description="импорт решётки с ондуляторами + фокус + группы")
    ap.add_argument("--config", default="configuration_v6_sila2026_und", help="папка конфига в workspace/lte/")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--prefix", default="und_", help="префикс устройств (уникальный, без коллизий)")
    ap.add_argument("--n-slices", type=int, default=5)
    a = ap.parse_args()

    cfg_dir = os.path.join("workspace", "lte", a.config)
    lte = os.path.join(cfg_dir, "sila2026_und.lte")
    knobs = list(csv.DictReader(open(os.path.join(cfg_dir, "undulator_knobs.csv"), encoding="utf-8"), delimiter=";"))
    manifest = list(csv.DictReader(open(os.path.join(os.path.dirname(__file__), "placements.csv"), encoding="utf-8"), delimiter=";"))

    print(f"[1/3] импорт {lte} → решётка {a.lattice} (prefix {a.prefix})")
    s2.run(lte_path=lte, lat_id=a.lattice, exp_id=a.lattice, prefix=a.prefix, n_slices=a.n_slices,
           monitors=True)

    print(f"[2/3] патч фокуса/силы по инстансам ({len(knobs)} устройств)")
    groups = patch_focus(a.lattice, a.prefix, knobs, manifest)

    print(f"[3/3] реестр жёстких групп")
    path, n = write_groups(groups)
    print(f"[OK] фокус проставлен {len(groups)} устройствам; {path}: {n} групп")
    print(f"     проверка: py -3.14 ring_analysis.py twiss  (на решётке {a.lattice} — через ura_twiss)")


if __name__ == "__main__":
    main()
