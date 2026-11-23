# -*- coding: utf-8 -*-
"""Импорт стороннего синхротрона s2lattice (configuration_v3_syberia/s2lattice.lte, PyAT/elegant)
в БД проекта как решётку 2 для МАТРИЧНОЙ модели (kick-phys). Другая машина: 124 м, 12 SP,
KQUAD/KQUSE/KSEXT/KOCT/CSBEND. Физкоэф. из .lte у ЛЮБОГО типа: mat_G=K1/QP, mat_S=K2/QP, mat_O=K3/QP
(абсолютные), у изгибов ещё mat_B=−ANGLE/(QP·Leff) (⇒ трекер-k1=K1, θ=ANGLE — оптика верна независимо
от энергии пучка). Неучтённые параметры и типы — предупреждение, не молча. Геометрия — проход кольца
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
MAGKIND = {"KQUAD": "quad", "KQUSE": "quad", "KSEXT": "sext", "KOCT": "oct", "CSBEND": "dipole"}
ACTIVE = {"quad": "k1", "sext": "k2", "oct": "k3", "dipole": "theta_y", "dq": "theta_y|k1|k2",
          "kicker": "kick_y|kick_z"}
STL = {"quad": "QPL_asm.stl", "sext": "SP_asm.stl", "oct": "OP_asm.stl", "dipole": "DPL_asm.stl",
       "dq": "DPL_asm.stl", "kicker": "SP_asm.stl"}
# Параметры .lte, которые импорт ПЕРЕНОСИТ в БД; и те, что управляют интегратором/излучением elegant
# (SynchrotronSim задаёт их сам — пропуск безвреден). Всё прочее НЕНУЛЕВОЕ у магнита — предупреждение.
USED = {"L", "ANGLE", "K1", "K2", "K3", "F1", "F2", "F3", "XREFERENCE", "E1", "E2", "HGAP", "FINT",
        "FINT1", "FINT2", "EDGE1_EFFECTS", "EDGE2_EFFECTS", "EDGE_ORDER"}
HARMLESS = {"N_SLICES", "N_KICKS", "INTEGRATION_ORDER", "SYNCH_RAD", "ISR", "ISR1PART", "USE_RAD_DIST",
            "ADD_OPENING_ANGLE", "SQRT_ORDER", "NONLINEAR", "EXPANSION_ORDER", "GROUP"}


def lte_multipoles(e):
    """(K1, K2, K3) нормальных мультиполей элемента .lte — у ЛЮБОГО типа (KQUAD/KQUSE/KSEXT/KOCT/CSBEND).
    У CSBEND к ним добавляется дробная форма F1…F3 при XREFERENCE > 0 (мультиполи в долях дипольного поля на
    радиусе XREF — формат «гармоники карты на опорном радиусе»). Как elegant csbend.c: b_n += F_n·n!/XREFⁿ,
    b_n = K_n·ρ, ρ = L/ANGLE ⇒ K_n += F_n·n!/XREFⁿ·ANGLE/L."""
    k = [e.get("K1", 0.0), e.get("K2", 0.0), e.get("K3", 0.0)]
    xr = e.get("XREFERENCE", 0.0)
    if e["type"] == "CSBEND" and xr > 0 and e.get("ANGLE") and e.get("L"):
        inv_rho, term = e["ANGLE"] / e["L"], 1.0 / xr
        for i in range(3):
            k[i] += e.get(f"F{i + 1}", 0.0) * term * inv_rho
            term *= (i + 2) / xr
    return tuple(k)


def kind_of(e):
    """Вид элемента для ядра = ГЕОМЕТРИЯ тела, а не состав поля: мультиполи (K1/K2/K3) допустимы у любого
    вида и пишутся абсолютными (mat_of). Изгиб с градиентом помечается «dq» (метка; физика у dipole и dq одна —
    секторное тело). KSEXT БЕЗ поля (все K = 0) — слот КОРРЕКТОРА (SH у СИЛА, SX/SZ у Сибири): kind=kicker,
    кик задаётся коррекцией/опытами через KICK/VKICK. KSEXT с одним K1 остаётся секступолем (прямое тело с
    градиентом)."""
    k = MAGKIND[e["type"]]
    if k == "dipole" and lte_multipoles(e)[0]:
        return "dq"
    if k == "sext" and not any(lte_multipoles(e)):
        return "kicker"
    return k


def active_of(kind, e):
    """active_coeffs устройства: собственные моменты вида + присутствующие в .lte мультиполи других порядков."""
    acts = ACTIVE[kind].split("|")
    acts += [f"k{n}" for n, v in zip((1, 2, 3), lte_multipoles(e)) if v and f"k{n}" not in acts]
    return "|".join(acts)


def unused_params(e):
    """Ненулевые параметры магнита .lte, которых SynchrotronSim не учитывает (K4+, B1…B8, F4+, G1…G8, J1,
    TILT, DX/DY, FSE, HKICK/VKICK…) + строковые (SYSTEMATIC_MULTIPOLES=… и т.п.) — для предупреждения.
    F1…F3 без XREFERENCE или без изгиба тоже не учитываются."""
    out = [k for k, v in e.items() if k not in USED | HARMLESS | {"type", "_str"} and v]
    if not (e["type"] == "CSBEND" and e.get("XREFERENCE", 0.0) > 0 and e.get("ANGLE")):
        out += [k for k in ("F1", "F2", "F3") if e.get(k)]
    return sorted(out) + sorted(e.get("_str", ()))


def parse():
    """(elems, seq, skipped): определения поддержанных типов, порядок кольца и {имя: тип} определений
    НЕПОДДЕРЖАННЫХ типов (MULT, QUAD, SBEN… — их длина в геометрию не попадёт; см. _warn_skipped)."""
    txt = open(LTE, encoding="utf-8").read()
    elems, skipped = {}, {}
    # ВАЖНО: после типа стоит [ 	]*, а не \s* — \s захватывает перевод строки, и определение БЕЗ
    # параметров («BPM: MONI») съедало следующую строку: BPM получал чужую длину, а тот элемент терялся.
    for m in re.finditer(r"^([A-Za-z0-9_]+)\s*:\s*([A-Za-z]+)[ 	]*(.*)$", txt, re.M):
        name, typ, rest = m.group(1), m.group(2).upper(), m.group(3).split("!", 1)[0]   # «!» — комментарий
        d = {"type": typ}
        for kv in re.finditer(r"([A-Za-z0-9_]+)\s*=\s*([-+\d.eE]+)", rest):
            d[kv.group(1).upper()] = float(kv.group(2))
        strs = [kv.group(1).upper() for kv in re.finditer(r'([A-Za-z0-9_]+)\s*=\s*"', rest)]
        if strs:
            d["_str"] = strs                                   # строковые параметры (файлы мультиполей и т.п.)
        if typ in ("DRIF", "KQUAD", "KQUSE", "KSEXT", "KOCT", "CSBEND", "MARK", "MONI"):
            elems[name] = d
        elif typ != "LINE":
            skipped[name] = (typ, d.get("L", 0.0))
    mline = re.search(r"ring\s*:\s*LINE\s*=\s*\((.*?)\)", txt, re.S | re.I)
    seq = [s.strip() for s in re.sub(r"[&\s]", " ", mline.group(1)).replace(",", " ").split() if s.strip()]
    return elems, seq, skipped


def mat_of(kind, e, L):
    """K из .lte → физкоэф. БД. ЗНАК ИЗГИБА: положительный ANGLE = поворот ВПРАВО (общая конвенция
    проекта, см. lte_to_lattice.step) ⇒ θ0 = −ANGLE, значит mat_B = −ANGLE/(QP·L). Мультиполи знака
    не меняют: кик ∝xⁿ переживает отражение при нечётном n (квадруполь x¹, октуполь x³), а чётные
    (диполь x⁰, секступоль x²) — нет; секступолю переворот не нужен именно потому, что геометрия
    уложена НЕ зеркально. Раньше здесь было +ANGLE и th+=ang — решётка выходила зеркальной, и
    секступоли получали неверный знак (та же ошибка, что чинилась у СИЛЫ, коммит c97742a)."""
    # Мультиполи — у ЛЮБОГО вида и всегда АБСОЛЮТНЫЕ (единая конвенция ядра; у «dipole» раньше ядро умножало
    # их на mat_B — CSBEND с K2 без K1 из-за этого получал неверный секступоль). Поле от формы полюсов
    # (секступоль квадруполя, градиент секступоля, октуполь изгиба…) доходит до ядра без потерь.
    K1, K2, K3 = lte_multipoles(e)
    m = dict(mat_B=0.0, mat_G=K1 / QP, mat_S=K2 / QP, mat_O=K3 / QP)
    if kind in ("dipole", "dq"):
        m["mat_B"] = -e.get("ANGLE", 0.0) / (QP * L) if L else 0.0
    return m


def edges_of(e):
    """Края CSBEND по правилам elegant: (E1, E2, L_f1, L_f2). EDGE1_EFFECTS/EDGE2_EFFECTS = 0 выключают
    грань целиком (угол и фриндж); L_f = 2·HGAP·FINTi (HGAP по умолч. 0, FINT 0.5, FINT1/2 < 0 ⇒ FINT) —
    см. matrix_tracker.lte_fringe_lengths, единственный источник этой конвенции."""
    import matrix_tracker as mt
    lf1, lf2 = mt.lte_fringe_lengths(e)
    E1 = e.get("E1", 0.0) if e.get("EDGE1_EFFECTS", 1.0) != 0.0 else 0.0
    E2 = e.get("E2", 0.0) if e.get("EDGE2_EFFECTS", 1.0) != 0.0 else 0.0
    return E1, E2, lf1, lf2


def _warn_edge_methods(devs):
    """SynchrotronSim реализует край elegant по умолчанию (EDGE_EFFECTS=1, EDGE_ORDER=1: линейная линза,
    tg E/ρ и tg(E−ψ)/ρ). Другие методы elegant (2–5, EDGE_ORDER≥2) дают ещё нелинейные/нулевого порядка
    члены — предупредить, что в ядре их нет."""
    odd = sorted(nm for nm, e in devs.items() if e["type"] == "CSBEND" and (
        e.get("EDGE1_EFFECTS", 1.0) not in (0.0, 1.0) or e.get("EDGE2_EFFECTS", 1.0) not in (0.0, 1.0)
        or e.get("EDGE_ORDER", 1.0) >= 2.0))
    if odd:
        print(f"[!] у {len(odd)} изгибов метод края elegant ≠ 1 или EDGE_ORDER≥2 ({', '.join(odd[:5])}"
              f"{'…' if len(odd) > 5 else ''}): SynchrotronSim считает край методом 1 (линейная линза "
              f"tg E/ρ, по вертикали tg(E−ψ)/ρ) — нелинейных членов этих методов в ядре нет")


def _warn_unused(devs):
    """Сводка параметров .lte, которые импорт НЕ переносит (раньше выбрасывались молча)."""
    by_par = {}
    for nm, e in devs.items():
        for p in unused_params(e):
            by_par.setdefault(p, []).append(nm)
    for p, names in sorted(by_par.items()):
        print(f"[!] параметр {p} не учитывается SynchrotronSim ({len(names)} эл.: {', '.join(sorted(names)[:5])}"
              f"{'…' if len(names) > 5 else ''}) — ядро знает нормальные K1…K3 (и F1…F3 при XREFERENCE у CSBEND)")


def _warn_skipped(seq, skipped):
    """Элементы кольца неподдержанных типов: их длина НЕ попадёт в геометрию (предупредить, а не молча)."""
    cnt = {}
    for nm in seq:
        if nm in skipped:
            t, L = skipped[nm]
            c = cnt.setdefault(t, [0, 0.0]); c[0] += 1; c[1] += L
    for t, (n, L) in sorted(cnt.items()):
        print(f"[!] тип {t} не импортируется ({n} в кольце, суммарная длина {L:.6g} м выпадает из геометрии) — "
              f"сначала lte_flatten.py (QUAD/SEXT/SBEN… → K-типы); MULT (тонкий мультиполь) ядро не поддерживает")


def main():
    elems, seq, skipped = parse()
    # уникальные магниты-типы
    devs = {nm: e for nm, e in elems.items() if e["type"] in MAGKIND}
    _warn_unused(devs)
    _warn_skipped(seq, skipped)
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
        et_rows.append({"device_type_id": dev, "kind": kind, "active_coeffs": active_of(kind, e),
                        "n_slices": NSL, "notes": f"импорт .lte, {e['type']}"})
    upsert_by(os.path.join(TB, "hardware_library.csv"), "device_type_id", hw_hdr, hw_rows)
    upsert_by(os.path.join(TB, "element_transforms.csv"), "device_type_id", et_hdr, et_rows)

    # 2) lattice_configs — проход кольца
    lc_hdr = next(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    fringe = {nm: edges_of(e) for nm, e in devs.items() if e["type"] == "CSBEND"}
    if any(f[2] or f[3] for f in fringe.values()):
        lc_hdr = lc_hdr + [c for c in ("lf1", "lf2") if c not in lc_hdr]   # столбцы L_f — в конец схемы
    _warn_edge_methods(devs)
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
            E1, E2, lf1, lf2 = fringe.get(nm, (e.get("E1", 0.0), e.get("E2", 0.0), 0.0, 0.0))
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
                     edge1=f"{-E1:.15g}", edge2=f"{-E2:.15g}")
            if lf1 or lf2:                                     # длина зоны спада торцов L_f = 2·HGAP·FINT, м
                r.update(lf1=f"{lf1:.15g}" if lf1 else "", lf2=f"{lf2:.15g}" if lf2 else "")
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
    import exp_planes as ep                                    # явные плоскости: кольцо — старт s = 0 (шов), финиш пуст
    st, fin = ep.default_planes(int(LAT_ID))
    ep.write_experiment(EXP_ID, start=st, finish=fin)
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
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    head = rows[0] + [c for c in hdr if c not in rows[0]]      # новые столбцы (напр. lf1/lf2) — в конец схемы
    li = head.index("lattice_id")
    keep = [r + [""] * (len(head) - len(r)) for r in rows[1:] if r and r[li] != lat_id]
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
