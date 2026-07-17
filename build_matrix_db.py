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

# Папка конфигурации .lte, из которой берутся краевые углы E1/E2 и δ (add_edge_angles / add_edge_delta).
# ПАРАМЕТР, а не константа: иначе при загрузке другой конфигурации (напр. замкнутой configuration_v1_closed)
# края молча считались бы по .lte ПРЕЖНЕЙ конфигурации. Переопределяется --lte-dir (load_config передаёт).
LTE_DIR = os.path.join("workspace", "lte", "configuration_v1")


def _desc_lte():
    """Путь к *_desc.lte текущей конфигурации (fallback — старое место в корне workspace/lte)."""
    p = os.path.join(LTE_DIR, "SYLA_desc.lte")
    if os.path.exists(p):
        return p
    import glob as _glob
    got = _glob.glob(os.path.join(LTE_DIR, "*_desc.lte"))
    return got[0] if got else os.path.join("workspace", "lte", "SYLA_desc.lte")


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
    """Для каждого экземпляра решётки — предшественник/преемник по порядку instance_id внутри lattice.

    СТОРОННИЕ решётки СОХРАНЯЮТСЯ. Импорт чужого синхротрона (s2lattice_to_db) пишет связность НАПРЯМУЮ,
    минуя lattice_configs, поэтому слепая перезапись «всё из lattice_configs» их стирала — так пропадала
    Сибирь (реш. 2) при каждой перегенерации СИЛЫ. Решётки, которых нет в lattice_configs, переносим как есть."""
    p = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    by_lat = {}
    for r in rows:
        by_lat.setdefault(r["lattice_id"], []).append(r)
    out = os.path.join(TB, "element_connectivity.csv")
    keep = []                                                              # строки решёток не из lattice_configs
    if os.path.exists(out):
        for r in csv.reader(open(out, encoding="utf-8"), delimiter=";"):
            if r and r[0] != "lattice_id" and r[0] not in by_lat:
                keep.append(r)
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["lattice_id", "instance_id", "device_type_id", "prev_instance_id", "next_instance_id"])
        for lat, lst in by_lat.items():
            lst.sort(key=lambda r: int(r["instance_id"]))
            for i, r in enumerate(lst):
                prev = lst[i - 1]["instance_id"] if i > 0 else ""          # торцы пусты (не кольцо; при нужде замкнуть)
                nxt = lst[i + 1]["instance_id"] if i < len(lst) - 1 else ""
                w.writerow([lat, r["instance_id"], r["device_type_id"], prev, nxt])
        w.writerows(keep)                                                  # сторонние решётки — без изменений
    nkeep = len({r[0] for r in keep})
    return out, len(by_lat) + nkeep, len(rows) + len(keep)


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
    """Заполняет edge1/edge2 в lattice_configs КАК ЯЩИК: краевой угол = (θ0_полн + δ_полн)/2 на ВНЕШНИХ
    гранях сгруппированного диполя (внутренние грани срезов = 0, они гасятся), знак — по .lte-краю.
    θ0 = ANGLE из .lte; δ из hardware_library.edge_delta (по device_type). Группа = подряд идущие изгибы
    одного device (срезы длинного диполя); обмотки DQ (одна позиция, разные device, один lte_name) —
    каждая своя группа, но θ0/δ/край одинаковы → Magnet читает первую. Асимметричные (DQ2C: расчёт ≠
    .lte-край) → сохраняем per-slice .lte. Мост instance→.lte: compute_magnet_endpoints(lat).lte_name.
    Ядро ПОКА игнорирует (Фаза 3a→3b). Идемпотентна. Возврат (path, n_заполнено)."""
    import re
    import design_orbit as do
    from collections import OrderedDict
    lte_path = _desc_lte()
    p = os.path.join(TB, "lattice_configs.csv")
    data = list(csv.DictReader(open(p, encoding="utf-8"), delimiter=";"))
    if not data or "edge1" not in data[0]:
        return p, 0                                          # столбцов нет — сначала add_structural_columns()
    lte = {}                                                 # lte_name → (ANGLE, E1, E2), рад
    if os.path.exists(lte_path):
        desc = open(lte_path, encoding="utf-8").read()
        for m in re.finditer(r"([A-Za-z0-9_]+)\s*:\s*CSBEND\s*,([^\n]+)", desc):
            b = m.group(2)
            def _g(k, b=b):
                mm = re.search(rf"\b{k}=(-?[0-9.eE+-]+)", b); return float(mm.group(1)) if mm else 0.0
            lte[m.group(1)] = (_g("ANGLE"), _g("E1"), _g("E2"))
    dlt = {}                                                 # device_type → δ
    for r in csv.DictReader(open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";"):
        v = (r.get("edge_delta") or "").strip()
        if v:
            try: dlt[r["device_type_id"]] = float(v)
            except ValueError: pass
    by_lat = OrderedDict()
    for r in data:
        by_lat.setdefault(r["lattice_id"], []).append(r)
    nfill = 0
    for lat, insts in by_lat.items():
        insts.sort(key=lambda r: int(r["instance_id"]))
        for r in insts:
            r["edge1"] = r["edge2"] = ""                     # сброс (идемпотентность)
        try:
            eps, _ = do.compute_magnet_endpoints(int(lat))
        except Exception:
            eps = None
        if not (eps and len(eps) == len(insts)):
            continue
        names = [getattr(eps[j], "lte_name", "") or "" for j in range(len(insts))]
        j = 0
        while j < len(insts):
            if names[j] not in lte:
                j += 1; continue
            dev = insts[j]["device_type_id"]; k = j            # группа подряд идущих изгибов одного device
            while k + 1 < len(insts) and names[k + 1] in lte and insts[k + 1]["device_type_id"] == dev:
                k += 1
            grp = list(range(j, k + 1))
            th = sum(lte[names[i]][0] for i in grp)            # суммарный изгиб группы
            dl = sum(dlt.get(insts[i]["device_type_id"], 0.0) for i in grp)
            # ЗНАК: положительный ANGLE .lte = поворот ВПРАВО (θ0<0, см. lte_to_lattice.step), а E меняет
            # знак вместе с ANGLE ⇒ край отрицателен. Тогда kEdge=tan(E)·θ0/L чётен по изгибу: верт. фокус
            # края не зависит от того, в какую сторону гнёт диполь (иначе зеркальный диполь ДЕфокусировал бы).
            box = -(th + dl) / 2.0                             # край ящика на внешней грани
            e1o = lte[names[grp[0]]][1]; e2o = lte[names[grp[-1]]][2]  # .lte внешние (первый.E1, последний.E2)
            if abs(abs(box) - abs(e1o)) < 2e-5 and abs(abs(box) - abs(e2o)) < 2e-5:
                insts[grp[0]]["edge1"] = f"{(box if e1o >= 0 else -box):.12g}"   # ЯЩИК: только внешние грани
                insts[grp[-1]]["edge2"] = f"{(box if e2o >= 0 else -box):.12g}"
                nfill += 2
            else:                                              # асимметрия (DQ2C) — per-slice .lte (со знаком)
                for i in grp:
                    _, e1, e2 = lte[names[i]]
                    insts[i]["edge1"] = f"{-e1:.12g}"; insts[i]["edge2"] = f"{-e2:.12g}"; nfill += 2
            j = k + 1
    fields = list(data[0].keys())
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    shutil.copy(p, os.path.join(TB, "backups", "lattice_configs_pre_edges.csv"))
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";")
        w.writeheader(); w.writerows(data)
    return p, nfill


def add_edge_delta():
    """Добавляет в hardware_library колонку edge_delta: δ = E1+E2−ANGLE из .lte (радианы), на
    ДЛИНО-специфичный device_type. δ — фиксированный доворот полюсных граней сверх геометрического
    клина, КОНСТАНТА типа магнита (не зависит от ρ/изгиба). Используется при расчёте E1/E2 = геом.+δ/2.
    Мост lte_name→device: те же compute_magnet_endpoints + lattice_configs, что в add_edge_angles.
    Идемпотентна, трогает только edge_delta. Возврат (path, n_заполнено, {device: δ})."""
    import re
    import design_orbit as do
    lte_path = _desc_lte()
    p = os.path.join(TB, "hardware_library.csv")
    if not os.path.exists(lte_path):
        return p, 0, {}
    desc = open(lte_path, encoding="utf-8").read()
    dl = {}                                                    # lte_name → δ (рад)
    for m in re.finditer(r"([A-Za-z0-9_]+)\s*:\s*CSBEND\s*,([^\n]+)", desc):
        b = m.group(2)
        def g(k, b=b):
            mm = re.search(rf"\b{k}=(-?[0-9.eE+-]+)", b)
            return float(mm.group(1)) if mm else 0.0
        dl[m.group(1)] = g("E1") + g("E2") - g("ANGLE")
    lc = {r["instance_id"]: r["device_type_id"]
          for r in csv.DictReader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";")
          if r["lattice_id"] == "1"}
    try:
        eps, _ = do.compute_magnet_endpoints(1)
    except Exception:
        eps = []
    dev_delta = {}                                            # длино-специфичный device → δ (проверяем единственность)
    for e in eps:
        if e.lte_name in dl:
            dev = lc.get(str(e.global_id), e.device_type_id)
            d = round(dl[e.lte_name], 12)
            if dev in dev_delta and abs(dev_delta[dev] - d) > 1e-9:
                raise ValueError(f"δ неоднозначна для {dev}: {dev_delta[dev]} vs {d} — проверь суффикс длины")
            dev_delta[dev] = d
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter=";"))
    hdr = rows[0]
    os.makedirs(os.path.join(TB, "backups"), exist_ok=True)
    if "edge_delta" not in hdr:
        shutil.copy(p, os.path.join(TB, "backups", "hardware_library_pre_delta.csv"))
        hdr.append("edge_delta")
        for r in rows[1:]:
            r.append("")
    i = hdr.index("edge_delta")
    n = 0
    for r in rows[1:]:
        if r[0] in dev_delta:
            r[i] = f"{dev_delta[r[0]]:.12g}"; n += 1
    with open(p, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    return p, n, dev_delta


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
    dp, nd, dd = add_edge_delta()
    print(f"hardware_library + edge_delta (δ из .lte): {dp}  ({nd} типов: {dict((k, round(v*1e3,4)) for k,v in dd.items())})")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Расширение БД SynchrotronSim под матричную модель.")
    ap.add_argument("--lte-dir", default=None,
                    help="папка конфигурации .lte для краевых углов E1/E2 и δ "
                         "(по умолчанию workspace/lte/configuration_v1)")
    a = ap.parse_args()
    if a.lte_dir:
        LTE_DIR = a.lte_dir
    main()
