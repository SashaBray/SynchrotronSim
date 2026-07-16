# -*- coding: utf-8 -*-
"""Расширение БД SynchrotronSim под матричную модель. Добавляет/создаёт:
  1) hardware_library.csv  — колонка eff_length (эфф. длина, м) каждому элементу;
  2) element_transforms.csv — какие матричные преобразования делает элемент библиотеки
     (kind + активные физкоэффициенты + число срезов тела);
  3) element_connectivity.csv — для каждого экземпляра решётки: предшественник и преемник.

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

# --- какие преобразования делает каждый элемент библиотеки ---
# (kind, активные силы матрицы, срезов тела, примечание)  [модель canonical/composite убрана]
TRANSFORMS = {
    "QP_L_1A":  ("quad", "k1", 40, "квадруполь длинный"),
    "QP_sh_1A": ("quad", "k1", 40, "квадруполь короткий"),
    "SP_I1_1A": ("sext", "k2", 40, "секступоль, обмотка I1"),
    "SP_I2_1A": ("sext", "k2", 40, "секступоль, обмотка I2"),
    "OP_1A":    ("oct",  "k3", 40, "октуполь"),
    "DQ_I1_1A": ("dq", "theta_y|k1|k2", 40, "диполь-квадруполь combined, обмотка I1"),
    "DQ_I2_1A": ("dq", "theta_y|k1|k2", 40, "диполь-квадруполь combined, обмотка I2"),
    "DQ_I3_1A": ("dq", "theta_y|k1|k2", 40, "диполь-квадруполь combined, обмотка I3"),
    "DP_L_PM":  ("dipole", "theta_y", 40, "поворотный ПМ длинный"),
    "DP_sh_PM": ("dipole", "theta_y", 40, "поворотный ПМ короткий"),
    # 5-обмоточный корректор (main + 4 корректирующие): каждая обмотка — со-локованный кикер,
    # H+V-эффект берётся из матрицы связи per_B/per_Bs (фит карт v3); суперпозиция как у DQ.
    "corrector_main_100A": ("kicker", "kick_y|kick_z", 40, "корректор 5-обм.: главная"),
    "corrector_c1_100A": ("kicker", "kick_y|kick_z", 40, "корректор 5-обм.: корр.1"),
    "corrector_c2_100A": ("kicker", "kick_y|kick_z", 40, "корректор 5-обм.: корр.2"),
    "corrector_c3_100A": ("kicker", "kick_y|kick_z", 40, "корректор 5-обм.: корр.3"),
    "corrector_c4_100A": ("kicker", "kick_y|kick_z", 40, "корректор 5-обм.: корр.4"),
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
    # 5-обмоточный корректор (corrector_main/c1..c4_100A) — per_* уже в coupling_matrix (фит карт v3),
    # write_coupling их СОХРАНЯЕТ (не трогаем здесь). Старые corrector_I1/I2 удалены из библиотеки.
    return R


def add_lengths():
    """Добавляет в hardware_library эффективную длину eff_length (полевая, для силы И геометрии).
    phys_length (по металлу) УБРАН — кик работает только по L_eff. Идемпотентно."""
    p = os.path.join(TB, "hardware_library.csv")
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter=";"))
    hdr = rows[0]
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    if "eff_length" not in hdr:
        shutil.copy(p, os.path.join(TB, "backups", "hardware_library_pre_lengths.csv"))
    for col, data in (("eff_length", EFF),):
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
    (корректор, BPM, синтетика). Столбец model (canonical/composite) УБРАН; n_slices=40 у всех."""
    p = os.path.join(TB, "element_transforms.csv")
    hdr = ["device_type_id", "kind", "active_coeffs", "n_slices", "notes"]
    existing = {}
    if os.path.exists(p):
        for r in csv.DictReader(open(p, encoding="utf-8"), delimiter=";"):   # по имени — устойчиво к старой схеме с model
            existing[r["device_type_id"]] = [r["device_type_id"], r.get("kind") or "", r.get("active_coeffs") or "",
                                             "40", r.get("notes") or ""]     # n_slices=40 всем
    for dev, (kind, coeffs, ns, note) in TRANSFORMS.items():
        existing[dev] = [dev, kind, coeffs, str(ns), note]
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
    # СТРУКТУРНЫЕ столбцы (eff_length/kind/n_slices/edge1/edge2) держим В КОНЦЕ схемы, ПОСЛЕ mat_*
    # — их заполняет add_structural_columns(). Исключаем из base, чтобы не задвоить и не переставить.
    STRUCT = ["eff_length", "kind", "n_slices", "edge1", "edge2"]
    base = [c for c in data[0].keys() if not c.startswith("mat_") and c not in STRUCT]
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
        w = csv.DictWriter(f, fieldnames=base + ["mat_B", "mat_G", "mat_S", "mat_O", "mat_Gv"],
                           delimiter=";", extrasaction="ignore")   # дропает скью mat_*s и структурные (их добавит add_structural_columns)
        w.writeheader(); w.writerows(data)
    return p, ndip


def add_structural_columns():
    """ЕДИНЫЙ ИСТОЧНИК: дописывает в КОНЕЦ схемы lattice_configs 6 СТРУКТУРНЫХ столбцов
    (eff_length, kind, n_slices, edge1, edge2) — из них матричное ядро берёт ВСЮ per-элементную
    геометрию/тип. Значения по device_type_id: eff_length из hardware_library, kind/n_slices из
    element_transforms; edge1/edge2 — краевые углы E1/E2 из .lte (Фаза 3, пока пусто/сохраняются).
    Идемпотентна: mat_* и прочие столбцы НЕ трогает, только эти 6 (фикс-порядок после mat_*).
    C++ ядро читает их по индексам [20..25]."""
    p = os.path.join(TB, "lattice_configs.csv")
    data = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    if not data:
        return p, 0
    def _tbl(name):
        return {r["device_type_id"]: r for r in
                csv.DictReader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";")}
    hw, tr = _tbl("hardware_library.csv"), _tbl("element_transforms.csv")
    STRUCT = ["eff_length", "kind", "n_slices", "edge1", "edge2"]
    DROP = ["mat_Bs", "mat_Gs", "mat_Ss", "mat_Os", "model"]   # СКЬЮ + model (canonical/composite) удалены — дропаем из схемы
    for r in data:
        dev = r["device_type_id"]; h = hw.get(dev, {}); t = tr.get(dev, {})
        r["eff_length"] = h.get("eff_length", "")
        r["kind"] = t.get("kind", "")
        r["n_slices"] = t.get("n_slices") or "40"              # число срезов тела (единое 40; модель убрана)
        r["edge1"] = r.get("edge1", "")            # E1 из .lte — Фаза 3 (сохранить, если уже задан)
        r["edge2"] = r.get("edge2", "")            # E2 из .lte
    base = [c for c in data[0].keys() if c not in STRUCT and c not in DROP]
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(p, os.path.join(TB, "backups", "lattice_configs_pre_structcols.csv"))
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=base + STRUCT, delimiter=";", extrasaction="ignore")   # дропает скью mat_*s
        w.writeheader(); w.writerows(data)
    return p, len(data)


def add_edge_angles():
    """Заполняет edge1/edge2 в lattice_configs краевыми углами E1/E2 из .lte (per-instance, радианы).
    Мост instance→.lte: design_orbit.compute_magnet_endpoints(lat) даёт список в порядке global_id
    (== instance_id) с полем lte_name, выровненный 1:1 с insts (как в add_matrix_columns); E1/E2 берём
    из lte_to_lattice.parse_desc(SYLA_desc.lte)[lte_name] (там уже float). У НЕизгибов (квад/секст/окт/
    корректор/BPM) E1/E2 нет → пусто. Стойка DQ: все обмотки — один магнит, один lte_name, один E1/E2.
    Ядро эти столбцы ПОКА игнорирует (Фаза 3a). close_scale их не масштабирует — значения дословно.
    Идемпотентна, трогает только edge1/edge2 (после add_structural_columns, где они уже заведены)."""
    import design_orbit as do
    import lte_to_lattice as lte
    from collections import OrderedDict
    lte_path = os.path.join("workspace", "lte", "configuration_v1", "SYLA_desc.lte")
    if not os.path.exists(lte_path):
        lte_path = os.path.join("workspace", "lte", "SYLA_desc.lte")
    devices = lte.parse_desc(lte_path) if os.path.exists(lte_path) else {}
    p = os.path.join(TB, "lattice_configs.csv")
    data = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    if not data or "edge1" not in data[0]:
        return p, 0                                          # столбцов нет — сначала add_structural_columns()
    by_lat = OrderedDict()
    for r in data:
        by_lat.setdefault(r["lattice_id"], []).append(r)
    nfill = 0
    for lat, insts in by_lat.items():
        insts.sort(key=lambda r: int(r["instance_id"]))
        try:
            eps, _ = do.compute_magnet_endpoints(int(lat))
        except Exception:
            eps = None
        aligned = eps if (eps and len(eps) == len(insts)) else None
        for j, r in enumerate(insts):
            e1 = e2 = ""
            if aligned is not None:
                d = devices.get(getattr(aligned[j], "lte_name", "") or "", {})
                if d.get("E1") is not None:
                    e1 = f"{float(d['E1']):.12g}"
                if d.get("E2") is not None:
                    e2 = f"{float(d['E2']):.12g}"
            r["edge1"], r["edge2"] = e1, e2
            if e1 or e2:
                nfill += 1
    fields = list(data[0].keys())
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(p, os.path.join(TB, "backups", "lattice_configs_pre_edges.csv"))
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";")
        w.writeheader(); w.writerows(data)
    return p, nfill


def main():
    print("hardware_library + eff_length:", add_lengths())
    print("element_transforms:", write_transforms())
    conn, nlat, nrows = write_connectivity()
    print(f"element_connectivity: {conn}  ({nlat} решёток, {nrows} экземпляров)")
    mp, ndip = add_matrix_columns()
    print(f"lattice_configs + mat_B/G/S/O: {mp}  (диполям B из проектного изгиба: {ndip} шт.)")
    sp, nr = add_structural_columns()
    print(f"lattice_configs + структурные (eff_length/kind/n_slices/edge1/edge2): {sp}  ({nr} строк)")
    ep, ne = add_edge_angles()
    print(f"lattice_configs + краевые углы E1/E2 из .lte: {ep}  (заполнено {ne} изгибающих экземпляров)")


if __name__ == "__main__":
    main()
