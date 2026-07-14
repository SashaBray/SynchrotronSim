# -*- coding: utf-8 -*-
"""Расширение БД SynchrotronSim под матричную модель. Добавляет/создаёт:
  1) hardware_library.csv  — колонка eff_length (эфф. длина, м) каждому элементу;
  2) element_transforms.csv — какие матричные преобразования делает элемент библиотеки
     (kind + активные физкоэффициенты + модель canonical|composite + число срезов);
  3) coupling_matrix.csv   — матрица связи «входной ток/коэф. → ФИЗИЧЕСКИЕ коэффициенты»
     (на единицу входа; суперпозиция: у со-локованных обмоток вклады складываются);
  4) element_connectivity.csv — для каждого экземпляра решётки: предшественник и преемник.

Источники: composite_fullrange.json (L_eff, первичные коэф.), bip_matrix (матрица связи DQ, q/p, Bρ).
Физика: k1=(q/p)G, k2=(q/p)S, k3=(q/p)O, θy=(q/p)·B·L_eff ⇒ обратно G=k1/(q/p) и т.д."""
import os, csv, json, shutil
import bip_matrix as bm

TB = os.path.join("workspace", "tables")
QP, BRHO = bm.QP, bm.BRHO
J = json.load(open("workspace/calibration/composite_fullrange.json", encoding="utf-8"))


def prim(dev):
    import numpy as np
    d = J["single"][dev]; wv = d["working_vmax"]
    kpu = [x["k_per_unit"] for x in d["linearity"] if abs(x["v"]) <= wv + 1e-9]
    return float(np.mean(kpu)), d["L_eff"], d["is_coeff"]


# --- эффективные длины: АВТОРИТЕТНЫЙ источник — matrix_eff_length/eff_length.csv (L_eff_fit),
#     подогнанные по совпадению траекторий. НЕ размер карты 2×maxX (composite_fullrange для диполей/
#     октуполя ошибочно брал размер карты). Корректоры — номинал (ждут характеризации). ---
def _load_eff():
    tag2dev = {"QP_L": "QP_L_1A", "QP_sh": "QP_sh_1A", "SP_I1": "SP_I1_1A", "SP_I2": "SP_I2_1A",
               "OP": "OP_1A", "DQ_I1": "DQ_I1_1A", "DQ_I2": "DQ_I2_1A", "DQ_I3": "DQ_I3_1A",
               "DP_L": "DP_L_PM", "DP_sh": "DP_sh_PM"}
    rows = csv.DictReader(open("workspace/calibration/eff_length.csv", encoding="utf-8-sig"), delimiter=";")
    e = {tag2dev[r["tag"]]: float(r["L_eff_fit_mm"]) / 1000.0 for r in rows if r["tag"] in tag2dev}
    e["corrector_I1_1A"] = 0.15; e["corrector_I2_1A"] = 0.15
    return e


EFF = _load_eff()

# --- физическая длина ПО МЕТАЛЛУ (габарит, для геометрии) — от пользователя, на магнит ---
PHYS = {
    "QP_L_1A": 0.3745, "QP_sh_1A": 0.2120,
    "SP_I1_1A": 0.1660, "SP_I2_1A": 0.1660,
    "OP_1A": 0.0900,
    "DQ_I1_1A": 0.4000, "DQ_I2_1A": 0.4000, "DQ_I3_1A": 0.4000,
    "DP_L_PM": 0.4100, "DP_sh_PM": 0.0400,
    "corrector_I1_1A": 0.1080, "corrector_I2_1A": 0.1080,
}

# --- какие преобразования делает каждый элемент библиотеки ---
# (kind, активные силы матрицы, модель, срезов, примечание)
TRANSFORMS = {
    "QP_L_1A":  ("quad", "k1", "canonical", 1, "квадруполь длинный"),
    "QP_sh_1A": ("quad", "k1", "canonical", 1, "квадруполь короткий"),
    "SP_I1_1A": ("sext", "k2", "canonical", 1, "секступоль, обмотка I1"),
    "SP_I2_1A": ("sext", "k2", "canonical", 1, "секступоль, обмотка I2"),
    "OP_1A":    ("oct",  "k3", "canonical", 1, "октуполь"),
    "DQ_I1_1A": ("dq", "theta_y|k1|k2", "composite", 40, "диполь-квадруполь combined, обмотка I1"),
    "DQ_I2_1A": ("dq", "theta_y|k1|k2", "composite", 40, "диполь-квадруполь combined, обмотка I2"),
    "DQ_I3_1A": ("dq", "theta_y|k1|k2", "composite", 40, "диполь-квадруполь combined, обмотка I3"),
    "DP_L_PM":  ("dipole", "theta_y", "canonical", 1, "поворотный ПМ длинный"),
    "DP_sh_PM": ("dipole", "theta_y", "canonical", 1, "поворотный ПМ короткий"),
    "corrector_I1_1A": ("kicker", "kick_y", "canonical", 1, "корректор гориз. (ЖДЁТ характеризации)"),
    "corrector_I2_1A": ("kicker", "kick_z", "canonical", 1, "корректор верт. (ЖДЁТ характеризации)"),
}


def coupling_rows():
    """device → физкоэф. на единицу входа: per_B[Тл], per_G[Тл/м], per_S[Тл/м²], per_O[Тл/м³]."""
    R = {}
    for dev in ("QP_L_1A", "QP_sh_1A"):
        k1, L, isc = prim(dev); R[dev] = dict(B=0, G=k1 / QP, S=0, O=0, itype="current")
    for dev in ("SP_I1_1A", "SP_I2_1A"):
        k2, L, isc = prim(dev); R[dev] = dict(B=0, G=0, S=k2 / QP, O=0, itype="current")
    k3, L, isc = prim("OP_1A"); R["OP_1A"] = dict(B=0, G=0, S=0, O=k3 / QP, itype="current")
    # DQ — из 3-входовой матрицы связи (вклад каждой обмотки в B,G,S; суперпозиция)
    for j, dev in enumerate(("DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A")):
        R[dev] = dict(B=bm.DQ_COUPLING_B[j], G=bm.DQ_COUPLING_G[j], S=bm.DQ_COUPLING_S[j], O=0, itype="current")
    for dev in ("DP_L_PM", "DP_sh_PM"):
        ty, _, isc = prim(dev); R[dev] = dict(B=ty / (QP * EFF[dev]), G=0, S=0, O=0, itype="coeff")  # B=θy/((q/p)L)
    for dev in ("corrector_I1_1A", "corrector_I2_1A"):
        R[dev] = dict(B=0, G=0, S=0, O=0, itype="current")                    # ЖДЁТ характеризации
    return R


def add_lengths():
    """Добавляет в hardware_library две длины: eff_length (полевая, для силы) и phys_length
    (по металлу, для геометрии). Идемпотентно."""
    p = os.path.join(TB, "hardware_library.csv")
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter=";"))
    hdr = rows[0]
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    if "eff_length" not in hdr:
        shutil.copy(p, os.path.join(TB, "backups", "hardware_library_pre_lengths.csv"))
    for col, data in (("eff_length", EFF), ("phys_length", PHYS)):
        if col not in hdr:
            hdr.append(col)
            for r in rows[1:]:
                r.append("")
        i = hdr.index(col)
        for r in rows[1:]:
            if r[0] in data:                       # обновляем только БАЗОВЫЕ элементы;
                r[i] = f"{data[r[0]]:.4f}"          # синтетику/корректор/BPM НЕ трогаем (сохраняем их eff/phys)
    with open(p, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    return p


def write_transforms():
    """СОХРАНЯЮЩЕ: обновляет/добавляет базовые (TRANSFORMS), но НЕ трогает прочие строки
    (корректор, BPM, синтетика, добавленные вне шаблона)."""
    p = os.path.join(TB, "element_transforms.csv")
    hdr = ["device_type_id", "kind", "active_coeffs", "model", "n_slices", "notes"]
    existing = {}
    if os.path.exists(p):
        rd = list(csv.reader(open(p, encoding="utf-8"), delimiter=";"))
        for r in rd[1:]:
            if r:
                existing[r[0]] = r
    for dev, (kind, coeffs, model, ns, note) in TRANSFORMS.items():
        existing[dev] = [dev, kind, coeffs, model, str(ns), note]
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(existing.values())
    return p


def write_coupling():
    """СОХРАНЯЮЩЕ: обновляет/добавляет базовые (TRANSFORMS), сохраняет прочие (корректор/синтетика)."""
    p = os.path.join(TB, "coupling_matrix.csv"); R = coupling_rows()
    hdr = ["device_type_id", "input_type", "per_B_Tl", "per_G_Tlm", "per_S_Tlm2", "per_O_Tlm3",
           "per_Bs_Tl", "per_Gs_Tlm", "per_Ss_Tlm2", "per_Os_Tlm3", "notes"]
    existing = {}
    if os.path.exists(p):
        rd = list(csv.reader(open(p, encoding="utf-8"), delimiter=";"))
        for r in rd[1:]:
            if r:
                existing[r[0]] = r
    for dev in TRANSFORMS:
        c = R.get(dev, dict(B=0, G=0, S=0, O=0, itype="current"))
        Bs, Gs, Ss, Os = (c.get(k, 0.0) for k in ("Bs", "Gs", "Ss", "Os"))
        note = "ждёт характеризации" if dev.startswith("corrector") else ""
        existing[dev] = [dev, c["itype"], f"{c['B']:.6e}", f"{c['G']:.6e}", f"{c['S']:.6e}", f"{c['O']:.6e}",
                         f"{Bs:.6e}", f"{Gs:.6e}", f"{Ss:.6e}", f"{Os:.6e}", note]
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(existing.values())
    return p


def write_connectivity():
    """Для каждого экземпляра решётки — предшественник/преемник по порядку instance_id внутри lattice."""
    p = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    by_lat = {}
    for r in rows:
        by_lat.setdefault(r["lattice_id"], []).append(r)
    out = os.path.join(TB, "element_connectivity.csv")
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["lattice_id", "instance_id", "device_type_id", "prev_instance_id", "next_instance_id"])
        for lat, lst in by_lat.items():
            lst.sort(key=lambda r: int(r["instance_id"]))
            for i, r in enumerate(lst):
                prev = lst[i - 1]["instance_id"] if i > 0 else ""          # торцы пусты (не кольцо; при нужде замкнуть)
                nxt = lst[i + 1]["instance_id"] if i < len(lst) - 1 else ""
                w.writerow([lat, r["instance_id"], r["device_type_id"], prev, nxt])
    return out, len(by_lat), len(rows)


def add_matrix_columns():
    """Добавляет в lattice_configs физ-параметры для матриц (mat_B/G/S/O) — ВТОРОЙ вход рядом с
    токами/коэффициентами. По умолчанию = coupling×вход; диполям mat_B = из ПРОЕКТНОГО изгиба
    (θ из compute_magnet_endpoints), т.к. инжекционный B ≠ coeff×per_B. Трекер берёт либо эти
    физкоэф. напрямую, либо считает из токов."""
    import design_orbit as do
    from collections import OrderedDict
    p = os.path.join(TB, "lattice_configs.csv")
    data = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    base = [c for c in data[0].keys() if not c.startswith("mat_")]
    R = coupling_rows()
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(p, os.path.join(TB, "backups", "lattice_configs_pre_matcols.csv"))
    by_lat = OrderedDict()
    for r in data:
        by_lat.setdefault(r["lattice_id"], []).append(r)
    ndip = 0
    for lat, insts in by_lat.items():
        insts.sort(key=lambda r: int(r["instance_id"]))
        try:
            eps, _ = do.compute_magnet_endpoints(int(lat))
        except Exception:
            eps = None
        aligned = eps if (eps and len(eps) == len(insts)) else None
        def _f(x):                                          # безопасный парс ('-'/пусто → 0)
            try:
                return float(x)
            except (TypeError, ValueError):
                return 0.0
        for j, r in enumerate(insts):
            dev = r["device_type_id"]; c = R.get(dev, dict(B=0, G=0, S=0, O=0, itype="current"))
            inp = _f(r.get("arg_val")) if c["itype"] == "current" else _f(r.get("correct_coeff"))
            B, G, S, O = c["B"] * inp, c["G"] * inp, c["S"] * inp, c["O"] * inp
            if dev.endswith("_PM") and aligned is not None:          # диполь: B из проектного изгиба
                bend = aligned[j].exit_yaw - aligned[j].entry_yaw
                B = bend / (QP * EFF[dev]); ndip += 1
            r["mat_B"] = f"{B:.6e}"; r["mat_G"] = f"{G:.6e}"; r["mat_S"] = f"{S:.6e}"; r["mat_O"] = f"{O:.6e}"
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=base + ["mat_B", "mat_G", "mat_S", "mat_O",
                           "mat_Gv", "mat_Bs", "mat_Gs", "mat_Ss", "mat_Os"], delimiter=";")
        w.writeheader(); w.writerows(data)
    return p, ndip


def main():
    print("hardware_library + eff_length/phys_length:", add_lengths())
    print("element_transforms:", write_transforms())
    print("coupling_matrix:", write_coupling())
    conn, nlat, nrows = write_connectivity()
    print(f"element_connectivity: {conn}  ({nlat} решёток, {nrows} экземпляров)")
    mp, ndip = add_matrix_columns()
    print(f"lattice_configs + mat_B/G/S/O: {mp}  (диполям B из проектного изгиба: {ndip} шт.)")


if __name__ == "__main__":
    main()
