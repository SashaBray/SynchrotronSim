# -*- coding: utf-8 -*-
"""Импорт стороннего синхротрона s2lattice (configuration_v3_syberia/s2lattice.lte, PyAT/elegant)
в БД проекта как решётку 2 для МАТРИЧНОЙ модели (kick-phys). Другая машина: 124 м, 12 SP,
KQUAD/KSEXT/CSBEND. Физкоэф. из .lte: mat_G=K1/QP, mat_S=K2/QP, mat_O=K3/QP, mat_B=ANGLE/(QP·Leff)
(⇒ трекер-k1=K1, θ=ANGLE — оптика верна независимо от энергии пучка). Геометрия — проход кольца
(дрейф вдоль курса, CSBEND поворачивает курс на ANGLE). STL — плейсхолдеры SILA. Устройства s2_*.
Запуск: py -3.14 s2lattice_to_db.py"""
import re, math, csv, os
from bip_matrix import QP
import synthetic_assets as sa                                  # gen_stl: масштаб STL под длину

TB = os.path.join("workspace", "tables")
LTE = os.path.join("workspace", "lte", "configuration_v3_syberia", "s2lattice.lte")
LAT_ID, EXP_ID, GROUP = "2", "2", "1"
PFX = "s2_"
MAGKIND = {"KQUAD": "quad", "KSEXT": "sext", "KOCT": "oct", "CSBEND": "dipole"}
ACTIVE = {"quad": "k1", "sext": "k2", "oct": "k3", "dipole": "theta_y"}
STL = {"quad": "QPL_asm.stl", "sext": "SP_asm.stl", "oct": "OP_asm.stl", "dipole": "DPL_asm.stl"}


def parse():
    txt = open(LTE, encoding="utf-8").read()
    elems = {}
    for m in re.finditer(r"^([A-Za-z0-9_]+)\s*:\s*([A-Za-z]+)\s*(.*)$", txt, re.M):
        name, typ, rest = m.group(1), m.group(2).upper(), m.group(3)
        d = {"type": typ}
        for kv in re.finditer(r"([A-Za-z0-9]+)\s*=\s*([-\d.eE]+)", rest):
            d[kv.group(1).upper()] = float(kv.group(2))
        if typ in ("DRIF", "KQUAD", "KSEXT", "KOCT", "CSBEND", "MARK", "MONI"):
            elems[name] = d
    mline = re.search(r"ring\s*:\s*LINE\s*=\s*\((.*?)\)", txt, re.S | re.I)
    seq = [s.strip() for s in re.sub(r"[&\s]", " ", mline.group(1)).replace(",", " ").split() if s.strip()]
    return elems, seq


def mat_of(kind, e, L):
    m = dict(mat_B=0.0, mat_G=0.0, mat_S=0.0, mat_O=0.0)
    if kind == "quad":   m["mat_G"] = e.get("K1", 0.0) / QP
    elif kind == "sext": m["mat_S"] = e.get("K2", 0.0) / QP
    elif kind == "oct":  m["mat_O"] = e.get("K3", 0.0) / QP
    elif kind == "dipole":
        m["mat_B"] = e.get("ANGLE", 0.0) / (QP * L) if L else 0.0
        if e.get("K1"): m["mat_G"] = e["K1"] / QP           # комбинир.-функц. диполь
    return m


def main():
    elems, seq = parse()
    # уникальные магниты-типы
    devs = {nm: e for nm, e in elems.items() if e["type"] in MAGKIND}
    # 1) hardware_library + element_transforms
    hw_all = list(csv.reader(open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";"))
    hw_hdr = hw_all[0]; hw_col = {c: i for i, c in enumerate(hw_hdr)}
    hw_by_dev = {r[hw_col["device_type_id"]]: r for r in hw_all[1:] if r}
    et_hdr = next(csv.reader(open(os.path.join(TB, "element_transforms.csv"), encoding="utf-8"), delimiter=";"))
    REF = {"quad": "QP_L_1A", "sext": "SP_I1_1A", "oct": "OP_1A", "dipole": "DP_L_PM"}   # SILA-эталон STL по типу
    hw_rows, et_rows = [], []
    for nm, e in devs.items():
        kind = MAGKIND[e["type"]]; L = e.get("L", 0.0); dev = PFX + nm
        refrow = hw_by_dev.get(REF[kind])                      # синтетич. STL под длину L из SILA-эталона
        stl = (sa.gen_stl(REF[kind], dev, L, refrow, hw_col) if refrow else None) or STL[kind]
        hw_rows.append({"device_type_id": dev, "Active": "True", "field_type": "B",
                        "map_raw_path": "-", "map_bin_path": "-", "stl_path": stl,
                        "use_symmetry": "1", "nominal_arg": "-", "field_direction": "1",
                        "correct_coeff": "1", "eff_length": f"{L:.5f}",
                        "edge_ff": "0", "synthetic": "False", "scaled_from": ""})
        et_rows.append({"device_type_id": dev, "kind": kind, "active_coeffs": ACTIVE[kind],
                        "n_slices": "40", "notes": f"s2lattice {e['type']}"})
    upsert_by(os.path.join(TB, "hardware_library.csv"), "device_type_id", hw_hdr, hw_rows)
    upsert_by(os.path.join(TB, "element_transforms.csv"), "device_type_id", et_hdr, et_rows)

    # 2) lattice_configs — проход кольца
    lc_hdr = next(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    x = y = th = 0.0; iid = 0; lc_rows = []
    for nm in seq:
        e = elems.get(nm)
        if not e: continue
        L = e.get("L", 0.0); ang = e.get("ANGLE", 0.0); mid = th + ang / 2
        xc, yc = x + L / 2 * math.cos(mid), y + L / 2 * math.sin(mid)
        if nm in devs:
            iid += 1; kind = MAGKIND[e["type"]]
            r = {c: "" for c in lc_hdr}
            r.update(lattice_id=LAT_ID, instance_id=str(iid), device_type_id=PFX + nm,
                     x=f"{xc:.6f}", y=f"{yc:.6f}", z="0", yaw=f"{mid:.6f}", pitch="0", roll="0",
                     arg_val="0", correct_coeff="1")
            for k, v in mat_of(kind, e, L).items():
                r[k] = f"{v:.9g}"
            lc_rows.append(r)
        x += L * math.cos(mid); y += L * math.sin(mid); th += ang
    upsert_lattice(os.path.join(TB, "lattice_configs.csv"), lc_hdr, LAT_ID, lc_rows)

    # element_connectivity — ПОРЯДОК инстансов кольца (build_line читает ЕГО, не lattice_configs напрямую)
    ec_hdr = ["lattice_id", "instance_id", "device_type_id", "prev_instance_id", "next_instance_id"]
    ec_rows = []
    for k, r in enumerate(lc_rows):
        iid = int(r["instance_id"])
        ec_rows.append({"lattice_id": LAT_ID, "instance_id": r["instance_id"], "device_type_id": r["device_type_id"],
                        "prev_instance_id": str(iid - 1) if k > 0 else "",
                        "next_instance_id": str(iid + 1) if k < len(lc_rows) - 1 else ""})
    upsert_lattice(os.path.join(TB, "element_connectivity.csv"), ec_hdr, LAT_ID, ec_rows)

    # 3) эксперимент 2 (1 оборот) + пучок group 1
    ex_hdr = next(csv.reader(open(os.path.join(TB, "experiments.csv"), encoding="utf-8"), delimiter=";"))
    ex = {c: "" for c in ex_hdr}
    ex.update(exp_id=EXP_ID, lattice_id=LAT_ID, group_id=GROUP, dt="6,74E-13", max_turns="1", traj_stride="1")
    upsert_by(os.path.join(TB, "experiments.csv"), "exp_id", ex_hdr, [ex])
    print(f"[OK] s2lattice → БД: решётка {LAT_ID} ({iid} магнитов), эксперимент {EXP_ID}, группа {GROUP}")
    print(f"     устройств s2_*: {len(devs)}; периметр {x:.1f}м close ({x:.3f},{y:.3f}); изгиб {math.degrees(th):.1f}°")


def upsert_by(path, key, hdr, newrows):
    import shutil
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    head = rows[0]; ci = head.index(key); keys = {r[key] for r in newrows}
    keep = [r for r in rows[1:] if r and r[ci] not in keys]
    shutil.copy(path, path + ".bak_s2")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(head)
        w.writerows(keep); w.writerows([[r.get(c, "") for c in head] for r in newrows])


def upsert_lattice(path, hdr, lat_id, newrows):
    import shutil
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";")); head = rows[0]
    li = head.index("lattice_id")
    keep = [r for r in rows[1:] if r and r[li] != lat_id]
    shutil.copy(path, path + ".bak_s2")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(head)
        w.writerows(keep); w.writerows([[r.get(c, "") for c in head] for r in newrows])


def run(lte_path=None, lat_id=None, exp_id=None, group=None):
    """Параметризованный запуск (для load_config.py): переопределяет источник .lte и id решётки/эксперимента."""
    global LTE, LAT_ID, EXP_ID, GROUP
    if lte_path:        LTE = lte_path
    if lat_id is not None:  LAT_ID = str(lat_id)
    if exp_id is not None:  EXP_ID = str(exp_id)
    if group is not None:   GROUP = str(group)
    main()


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="импорт стороннего .lte (PyAT/elegant) → решётка БД")
    ap.add_argument("--lte", default=None, help="путь к .lte (по умолч. configuration_v3_syberia/s2lattice.lte)")
    ap.add_argument("--lattice", default=None, help="id решётки (по умолч. 2)")
    ap.add_argument("--exp", default=None, help="id эксперимента (по умолч. = id решётки)")
    a = ap.parse_args()
    run(lte_path=a.lte, lat_id=a.lattice, exp_id=a.exp)
