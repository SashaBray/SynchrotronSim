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
import lte_to_lattice as _geom                                  # ЕДИНСТВЕННЫЙ источник конвенции знака (step)

TB = os.path.join("workspace", "tables")
LTE = os.path.join("workspace", "lte", "configuration_v3_syberia", "s2lattice.lte")
LAT_ID, EXP_ID, GROUP = "2", "2", "1"
PFX = "s2_"                                                    # префикс имён устройств (см. run())
NSL = "40"                                                     # срезов на элемент в строке решётки
MONITORS = False                                               # заводить ли MONI как датчики положения
MONI_DEV = "BPM"                                               # готовое устройство kind=monitor из БД
MAGKIND = {"KQUAD": "quad", "KSEXT": "sext", "KOCT": "oct", "CSBEND": "dipole"}
ACTIVE = {"quad": "k1", "sext": "k2", "oct": "k3", "dipole": "theta_y", "dq": "theta_y|k1|k2",
          "kicker": "kick_y|kick_z"}
STL = {"quad": "QPL_asm.stl", "sext": "SP_asm.stl", "oct": "OP_asm.stl", "dipole": "DPL_asm.stl",
       "dq": "DPL_asm.stl", "kicker": "SP_asm.stl"}


def kind_of(e):
    """Вид элемента для ядра. Изгиб с градиентом — «dq», а НЕ «dipole»: у «dipole» ядро считает
    mat_G относительным градиентом (умножает на mat_B), у «dq» — абсолютным, как его и пишет mat_of.
    Секступоль с НУЛЕВОЙ K2 — это слот КОРРЕКТОРА (SH у СИЛА, SX/SZ у Сибири): kind=kicker (кик орбиты),
    как в решётке 1. Мат-коэф. нулевые, кик задаётся коррекцией/опытами через KICK/VKICK."""
    k = MAGKIND[e["type"]]
    if k == "dipole" and e.get("K1"):
        return "dq"
    if k == "sext" and not e.get("K2"):
        return "kicker"
    return k


def parse():
    txt = open(LTE, encoding="utf-8").read()
    elems = {}
    # ВАЖНО: после типа стоит [ 	]*, а не \s* — \s захватывает перевод строки, и определение БЕЗ
    # параметров («BPM: MONI») съедало следующую строку: BPM получал чужую длину, а тот элемент терялся.
    for m in re.finditer(r"^([A-Za-z0-9_]+)\s*:\s*([A-Za-z]+)[ 	]*(.*)$", txt, re.M):
        name, typ, rest = m.group(1), m.group(2).upper(), m.group(3)
        d = {"type": typ}
        for kv in re.finditer(r"([A-Za-z0-9_]+)\s*=\s*([-+\d.eE]+)", rest):
            d[kv.group(1).upper()] = float(kv.group(2))
        if typ in ("DRIF", "KQUAD", "KSEXT", "KOCT", "CSBEND", "MARK", "MONI"):
            elems[name] = d
    mline = re.search(r"ring\s*:\s*LINE\s*=\s*\((.*?)\)", txt, re.S | re.I)
    seq = [s.strip() for s in re.sub(r"[&\s]", " ", mline.group(1)).replace(",", " ").split() if s.strip()]
    return elems, seq


def mat_of(kind, e, L):
    """K из .lte → физкоэф. БД. ЗНАК ИЗГИБА: положительный ANGLE = поворот ВПРАВО (общая конвенция
    проекта, см. lte_to_lattice.step) ⇒ θ0 = −ANGLE, значит mat_B = −ANGLE/(QP·L). Мультиполи знака
    не меняют: кик ∝xⁿ переживает отражение при нечётном n (квадруполь x¹, октуполь x³), а чётные
    (диполь x⁰, секступоль x²) — нет; секступолю переворот не нужен именно потому, что геометрия
    уложена НЕ зеркально. Раньше здесь было +ANGLE и th+=ang — решётка выходила зеркальной, и
    секступоли получали неверный знак (та же ошибка, что чинилась у СИЛЫ, коммит c97742a)."""
    m = dict(mat_B=0.0, mat_G=0.0, mat_S=0.0, mat_O=0.0)
    if kind == "quad":   m["mat_G"] = e.get("K1", 0.0) / QP
    elif kind == "sext": m["mat_S"] = e.get("K2", 0.0) / QP
    elif kind == "oct":  m["mat_O"] = e.get("K3", 0.0) / QP
    elif kind in ("dipole", "dq"):
        m["mat_B"] = -e.get("ANGLE", 0.0) / (QP * L) if L else 0.0
        # у «dq» ядро берёт mat_G/mat_S абсолютными — так их и пишем (у «dipole» их не бывает)
        if e.get("K1"): m["mat_G"] = e["K1"] / QP
        if e.get("K2"): m["mat_S"] = e["K2"] / QP
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
    REF = {"quad": "QP_L_1A", "sext": "SP_I1_1A", "oct": "OP_1A", "dipole": "DP_L_PM",
           "dq": "DQ_I1_1A", "kicker": "SP_I1_1A"}                  # SILA-эталон STL по виду
    hw_rows, et_rows = [], []
    for nm, e in devs.items():
        kind = kind_of(e); L = e.get("L", 0.0); dev = PFX + nm
        refrow = hw_by_dev.get(REF[kind])                      # синтетич. STL под длину L из SILA-эталона
        stl = (sa.gen_stl(REF[kind], dev, L, refrow, hw_col) if refrow else None) or STL[kind]
        hw_rows.append({"device_type_id": dev, "Active": "True", "field_type": "B",
                        "map_raw_path": "-", "map_bin_path": "-", "stl_path": stl,
                        "use_symmetry": "1", "nominal_arg": "-", "field_direction": "1",
                        "correct_coeff": "1", "eff_length": f"{L:.12f}",
                        "edge_ff": "0", "synthetic": "False", "scaled_from": ""})
        et_rows.append({"device_type_id": dev, "kind": kind, "active_coeffs": ACTIVE[kind],
                        "n_slices": NSL, "notes": f"импорт .lte, {e['type']}"})
    upsert_by(os.path.join(TB, "hardware_library.csv"), "device_type_id", hw_hdr, hw_rows)
    upsert_by(os.path.join(TB, "element_transforms.csv"), "device_type_id", et_hdr, et_rows)

    # 2) lattice_configs — проход кольца
    lc_hdr = next(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    x = y = th = 0.0; iid = 0; lc_rows = []
    for nm in seq:
        e = elems.get(nm)
        if not e: continue
        L = e.get("L", 0.0); ang = e.get("ANGLE", 0.0)
        # ГЕОМЕТРИЯ — через lte_to_lattice.step: ЕДИНСТВЕННЫЙ источник конвенции знака (ANGLE>0 =
        # поворот вправо) и точной дуги. Своя формула здесь клала хорду по СРЕДНЕМУ направлению
        # (приближение) и крутила th += ang (зеркально) — решётка выходила несогласованной с СИЛОЙ.
        xc, yc, yaw_c, x, y, th = _geom.step(x, y, th, L, ang)
        if MONITORS and e["type"] in ("MONI", "MARK") and nm not in devs:
            # датчик положения: поля нет (дрейф). Длина из .lte (обычно 0 → плоскость); при L>0 — тонкий
            # ящик, чтение по центральной плоскости. Съём центроида, см. [[bpm-monitor-element]].
            iid += 1
            r = {c: "" for c in lc_hdr}
            r.update(lattice_id=LAT_ID, instance_id=str(iid), device_type_id=MONI_DEV,
                     x=f"{xc:.12f}", y=f"{yc:.12f}", z="0", yaw=f"{yaw_c:.12f}", pitch="0", roll="0",
                     arg_val="0", correct_coeff="1", eff_length=f"{L:.12f}", kind="monitor", n_slices="1",
                     mat_B="0", mat_G="0", mat_S="0", mat_O="0")
            lc_rows.append(r)
            continue
        if nm in devs:
            iid += 1; kind = kind_of(e)
            r = {c: "" for c in lc_hdr}
            r.update(lattice_id=LAT_ID, instance_id=str(iid), device_type_id=PFX + nm,
                     # геометрия пишется с полной точностью double: округление до микрона давало
                     # остаточную замкнутую орбиту ~2 мкм и мешало сверке аналитики с elegant
                     x=f"{xc:.12f}", y=f"{yc:.12f}", z="0", yaw=f"{yaw_c:.12f}", pitch="0", roll="0",
                     arg_val="0", correct_coeff="1",
                     # C++-ядро берёт длину, вид и число срезов ИЗ СТРОКИ решётки (Python умеет падать
                     # на hardware_library/element_transforms, C++ — нет): заполняем явно
                     eff_length=f"{L:.12f}", kind=kind, n_slices=NSL,
                     # края: знак тот же, что у угла (θ0 = −ANGLE) ⇒ edge = −E, как в решётке 1
                     edge1=f"{-e.get('E1', 0.0):.15g}", edge2=f"{-e.get('E2', 0.0):.15g}")
            for k, v in mat_of(kind, e, L).items():
                r[k] = f"{v:.15g}"
            lc_rows.append(r)
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
    n_mon = sum(1 for r in lc_rows if r["kind"] == "monitor")
    print(f"[OK] .lte → БД: решётка {LAT_ID} ({len(lc_rows) - n_mon} магнитов, {n_mon} датчиков), "
          f"эксперимент {EXP_ID}, группа {GROUP}")
    print(f"     устройств {PFX}*: {len(devs)}; периметр {sum(elems[n].get('L', 0) for n in seq if n in elems):.3f} м; "
          f"зазор кольца ({x:.6f},{y:.6f}) м; изгиб {math.degrees(th):.1f}°")


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


def run(lte_path=None, lat_id=None, exp_id=None, group=None, prefix=None, n_slices=None,
        monitors=None):
    """Параметризованный запуск (для load_config.py): источник .lte, id решётки/эксперимента,
    префикс имён устройств, число срезов и надо ли заводить датчики положения."""
    global LTE, LAT_ID, EXP_ID, GROUP, PFX, NSL, MONITORS
    if lte_path:        LTE = lte_path
    if lat_id is not None:  LAT_ID = str(lat_id)
    if exp_id is not None:  EXP_ID = str(exp_id)
    if group is not None:   GROUP = str(group)
    if prefix is not None:  PFX = prefix
    if n_slices is not None: NSL = str(n_slices)
    if monitors is not None: MONITORS = bool(monitors)
    main()


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="импорт стороннего .lte (PyAT/elegant) → решётка БД")
    ap.add_argument("--lte", default=None, help="путь к .lte (по умолч. configuration_v3_syberia/s2lattice.lte)")
    ap.add_argument("--lattice", default=None, help="id решётки (по умолч. 2)")
    ap.add_argument("--exp", default=None, help="id эксперимента (по умолч. = id решётки)")
    ap.add_argument("--prefix", default=None, help="префикс имён устройств (по умолч. s2_)")
    ap.add_argument("--n-slices", default=None, help="срезов на элемент в строке решётки (по умолч. 40)")
    ap.add_argument("--monitors", action="store_true", help="завести MONI как датчики положения (BPM)")
    a = ap.parse_args()
    run(lte_path=a.lte, lat_id=a.lattice, exp_id=a.exp, prefix=a.prefix, n_slices=a.n_slices,
        monitors=a.monitors or None)
