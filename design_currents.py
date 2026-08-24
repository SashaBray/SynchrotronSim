# -*- coding: utf-8 -*-
"""Проектные токи обмоток по МАТРИЦАМ СВЯЗИ из БД.

Читает постоянную таблицу `workspace/tables/coupling_matrices.csv` (её пишет `element_compare.py library`)
и проектные мультипольные коэффициенты `mat_*` решётки (`lattice_configs.csv`), решает обратную задачу
«проектные моменты → токи обмоток» и выдаёт таблицу проектных токов:

    сырой ток (arg карты, вход модели) → физический ток питания (А при 6 ГэВ) → % номинала.

Пересчёта траекторий НЕТ — только линейная алгебра над готовыми матрицами, поэтому работает мгновенно
(в отличие от `element_compare.py matrices`, который матрицы заново измеряет диффуром).

    py -3.14 design_currents.py                 # таблица на экран + workspace/tables/design_currents.csv
    py -3.14 design_currents.py --lattice 1     # решётка, откуда берутся проектные mat_* (по умолч. 1)
    py -3.14 design_currents.py --no-save       # только показать
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
TB = os.path.join(ROOT, "workspace", "tables")
COUPLING_TABLE = "coupling_matrices.csv"
STAMP_FILE = "coupling_matrices.inputs.json"     # отпечатки входных данных на момент измерения матриц
OUT_TABLE = "design_currents.csv"
MEASURE_LATTICE = 1                              # решётка, на геометрии которой меряются матрицы
# Поля, реально влияющие на измеренную матрицу: карта, её масштаб/ориентация и геометрия элемента
HW_FIELDS = ("map_bin_path", "nominal_arg", "use_symmetry", "field_direction", "correct_coeff",
             "rot_x_deg", "rot_y_deg", "rot_z_deg")
LC_FIELDS = ("eff_length", "kind")

_BASE_RE = re.compile(r"_I\d+|_c\d+|_main")            # суффикс обмотки (как в element_compare.element_windings)
_VARIANT_RE = re.compile(r"_sext$")                    # вариант МОДЕЛИ того же магнита (DQ с секступольным киком)


def _base_name(device: str) -> str:
    """Имя элемента без суффикса обмотки: SP_I1_1A_L166 → SP_1A_L166."""
    return _BASE_RE.sub("", device)


def _strip_variant(element: str) -> str:
    """Убрать пометку варианта модели: `DQ_I1_1A_L400_sext` → `DQ_I1_1A_L400`.

    Один и тот же магнит может стоять в таблице несколькими моделями (у DQ — с секступольным киком и
    без). Обмотки у них ОДНИ И ТЕ ЖЕ, поэтому устройства ищем по имени без пометки.
    """
    return _VARIANT_RE.sub("", element)


def _winding_label(device: str) -> str:
    """Метка обмотки устройства: SP_I1_1A_L166 → I1; одиночная обмотка → I."""
    m = re.search(r"_(I\d+|c\d+|main)", device)
    return m.group(1) if m else "I"


def _norm_label(label: str) -> str:
    """Метка тока из таблицы связи → метка обмотки устройства: `I_c1` → `c1`, `I1`/`I` — как есть.
    (в таблице токи названы по `element_compare._curname`: c1/main получают префикс `I_`)."""
    return label[2:] if (label.startswith("I_") and len(label) > 2) else label


def load_coupling(tables_dir: str = TB) -> list[dict]:
    """Матрицы связи из БД → список dict в формате element_compare (name/kind/labels/outputs/M/devs...).

    Прямые (`forward`) строки собираются в матрицу M [моменты × обмотки]; `fixed` (PM — поле задано
    конструкцией, обмоток нет) отдаётся отдельной веткой с pm=True.
    """
    path = os.path.join(tables_dir, COUPLING_TABLE)
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"нет таблицы матриц связи: {path}\n"
            f"построить: py -3.14 element_compare.py library")
    with open(path, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.DictReader(f, delimiter=";"))

    by_el: dict[str, list[dict]] = {}
    for r in rows:
        by_el.setdefault(r["element"], []).append(r)

    mats = []
    for el, rs in by_el.items():
        kind = rs[0]["kind"]
        leff = rs[0].get("eff_length") or ""
        fixed = [r for r in rs if r["direction"] == "fixed"]
        if fixed:                                            # PM: обмоток нет, поле фиксировано
            mats.append(dict(name=el, kind=kind, pm=True, L=leff,
                             outputs=[r["out_label"] for r in fixed],
                             fixed={r["out_label"]: float(r["value"]) for r in fixed}))
            continue
        fwd = [r for r in rs if r["direction"] == "forward"]
        if not fwd:
            continue
        outputs, labels = [], []                             # порядок как в таблице (он же порядок матрицы)
        for r in fwd:
            if r["out_label"] not in outputs: outputs.append(r["out_label"])
            if r["in_label"] not in labels: labels.append(r["in_label"])
        M = np.zeros((len(outputs), len(labels)))
        for r in fwd:
            M[outputs.index(r["out_label"]), labels.index(r["in_label"])] = float(r["value"])
        mats.append(dict(name=el, kind=kind, pm=False, L=leff, outputs=outputs, labels=labels, M=M,
                         i_op=rs[0].get("i_op") or "", imax=rs[0].get("imax") or "", r2=rs[0].get("r2") or ""))
    return sorted(mats, key=lambda m: m["name"])


def attach_devices(mats: list[dict], lcrow: dict) -> None:
    """Восстановить устройства обмоток (`devs`) сопоставлением с решёткой по базовому имени и метке.

    Обратный разбор имени ненадёжен (у DQ имя элемента — это имя ГЛАВНОЙ обмотки `DQ_I1_1A_L400`, у
    секступоля — база `SP_1A_L166`), поэтому и имя элемента, и device_type_id решётки приводятся к
    ОДНОЙ базе (без суффикса обмотки), а обмотка выбирается по метке.
    """
    by_base: dict[str, dict[str, str]] = {}
    for dev in lcrow:
        by_base.setdefault(_base_name(dev), {})[_winding_label(dev)] = dev
    for m in mats:
        cand = by_base.get(_base_name(_strip_variant(m["name"])), {})
        if m.get("pm"):
            m["devs"] = [cand.get("I") or next(iter(cand.values()), m["name"])]
            continue
        m["devs"] = [cand.get(_norm_label(lb), m["name"]) for lb in m["labels"]]


REBUILD_HINT = "пересобрать: py -3.14 element_compare.py library --no-report --no-demo --no-html"


def _read_table(name: str, tables_dir: str = TB) -> list[dict]:
    path = os.path.join(tables_dir, name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=";"))


def inputs_hash(devs, hw: dict, lcrow: dict) -> str:
    """Отпечаток ВХОДНЫХ ДАННЫХ, из которых мерилась матрица элемента.

    Берём только то, что действительно меняет измеренную матрицу: карту поля и её масштаб/ориентацию
    (`hardware_library`) плюс геометрию элемента (`lattice_configs`). Витки и номинал сюда НЕ входят —
    они пересчитывают arg в амперы уже после и читаются всегда свежими.
    """
    parts: list[str] = []
    for d in sorted(devs):
        h, l = hw.get(d, {}), lcrow.get(d, {})
        parts.append(d)
        parts += [f"{k}={(h.get(k) or '').strip()}" for k in HW_FIELDS]
        parts += [f"{k}={(l.get(k) or '').strip()}" for k in LC_FIELDS]
    return hashlib.sha1("|".join(parts).encode("utf-8")).hexdigest()[:16]


def current_hashes(tables_dir: str = TB) -> dict[str, str]:
    """Отпечатки входных данных СЕЙЧАС, по элементам таблицы матриц (база измерения — решётка 1)."""
    hw = {r["device_type_id"]: r for r in _read_table("hardware_library.csv", tables_dir)}
    lcrow = {r["device_type_id"]: r for r in _read_table("lattice_configs.csv", tables_dir)
             if str(r.get("lattice_id")) == str(MEASURE_LATTICE)}
    if not hw or not lcrow:
        return {}
    mats = load_coupling(tables_dir)
    attach_devices(mats, lcrow)
    return {m["name"]: inputs_hash(m["devs"], hw, lcrow) for m in mats}


def read_stamps(tables_dir: str = TB) -> dict[str, str]:
    """Отпечатки, записанные в момент ИЗМЕРЕНИЯ матриц (пишет element_compare при обновлении таблицы)."""
    path = os.path.join(tables_dir, STAMP_FILE)
    if not os.path.exists(path):
        return {}
    try:
        with open(path, encoding="utf-8") as f:
            return dict(json.load(f).get("elements", {}))
    except (OSError, json.JSONDecodeError, AttributeError):
        return {}


def write_stamps(stamps: dict[str, str], tables_dir: str = TB) -> str:
    """Дописать отпечатки обновлённых элементов (остальные не трогаем — таблица обновляется поэлементно)."""
    path = os.path.join(tables_dir, STAMP_FILE)
    data = read_stamps(tables_dir)
    data.update(stamps)
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"_doc": "Отпечатки входных данных (карта поля, ориентация, геометрия) на момент "
                           "измерения матриц связи. Расхождение с текущими ⇒ матрица элемента устарела.",
                   "elements": data}, f, ensure_ascii=False, indent=1, sort_keys=True)
    return path


def staleness(tables_dir: str = TB) -> str:
    """Предупреждение, если матрицы в базе измерены по ДРУГИМ входным данным, чем лежат сейчас.

    Сравниваются не даты файлов, а сами данные: дата ненадёжна — `element_compare` переписывает
    `lattice_configs.csv` по ходу работы (временные одно-элементные конфигурации, демо-конфигурации),
    и таблица оказывалась «старее» своих же входов сразу после честной пересборки.
    Пустая строка ⇒ матрицы соответствуют текущей машине.
    """
    if not os.path.exists(os.path.join(tables_dir, COUPLING_TABLE)):
        return ""
    stamps, now = read_stamps(tables_dir), current_hashes(tables_dir)
    if not now:
        return ""
    if not stamps:                                            # таблица снята до появления отпечатков
        return (f"ВНИМАНИЕ: у матриц связи в базе нет отпечатка входных данных (таблица снята прежней "
                f"версией) — соответствие текущим картам поля и геометрии не проверить.\n           {REBUILD_HINT}")
    stale = sorted(el for el, h in now.items() if stamps.get(el) != h)
    if not stale:
        return ""
    head = ", ".join(stale[:4]) + (f" и ещё {len(stale) - 4}" if len(stale) > 4 else "")
    return (f"ВНИМАНИЕ: карты поля или геометрия менялись после измерения матриц — устарели "
            f"{len(stale)} из {len(now)} элементов ({head}); их проектные токи неверны.\n"
            f"           {REBUILD_HINT}")


def compute(lattice: int = 1, tables_dir: str = TB):
    """Проектные токи всех элементов: (матрицы связи из БД) + (проектные mat_* решётки) → токи обмоток.

    Возврат (rows, mats): rows — строки итоговой таблицы (по обмоткам), mats — матрицы с добавленными
    полями `Ivec` (сырой arg по обмоткам), `ctrl` (управляемые моменты), `design`/`achieved` (проверка).
    """
    import element_compare as ec                          # хелперы перевода arg↔А и решатель проектных токов

    lcrow = {r["device_type_id"]: r for r in ec.mt._read("lattice_configs.csv")
             if str(r.get("lattice_id")) == str(lattice)}
    if not lcrow:
        raise SystemExit(f"в lattice_configs.csv нет строк решётки {lattice}")
    ec.LCROW = lcrow                                       # решать против ВЫБРАННОЙ решётки

    mats = load_coupling(tables_dir)
    attach_devices(mats, lcrow)

    rows = []
    for m in mats:
        if m.get("pm"):
            rows.append(dict(element=m["name"], kind=m["kind"], winding="—", device=m["devs"][0],
                             arg="", current_a="", turns="", nominal_a="", percent="", status="PM (фикс. поле)"))
            continue
        if any(d not in lcrow for d in m["devs"]):         # матрица есть, а элемента в этой решётке нет
            rows.append(dict(element=m["name"], kind=m["kind"], winding="—", device=",".join(m["devs"]),
                             arg="", current_a="", turns="", nominal_a="", percent="",
                             status=f"нет в решётке {lattice}"))
            continue
        Ivec, ctrl, design = ec._design_ivec(m)
        m["Ivec"], m["ctrl"], m["design"] = Ivec, ctrl, design
        m["achieved"] = np.asarray(m["M"], float) @ Ivec
        for j, lb in enumerate(m["labels"]):
            dev = m["devs"][j]
            amps = ec.arg_to_amps(dev, float(Ivec[j])) if ctrl else None
            nom = ec.dev_nominal(dev)
            tpc = ec.dev_tpc(dev)
            rows.append(dict(
                element=m["name"], kind=m["kind"], winding=lb, device=dev,
                arg=f"{Ivec[j]:.6f}" if ctrl else "",
                current_a=f"{amps:.3f}" if amps is not None else "",
                turns=str(tpc) if tpc else "",
                nominal_a=f"{nom:.1f}" if nom else "",
                percent=f"{100.0 * abs(amps) / nom:.1f}" if (amps is not None and nom) else "",
                status=("" if ctrl else "нет проектных mat_* (не управляется)")
                       + ("" if (amps is not None or not ctrl) else "нет витков в БД — только сырой arg")))
    return rows, mats


FIELDS = ["element", "kind", "winding", "device", "arg", "current_a", "turns", "nominal_a", "percent", "status"]
HEADERS = {"element": "элемент", "kind": "вид", "winding": "обмотка", "device": "устройство",
           "arg": "сырой ток (arg)", "current_a": "проектный ток, А", "turns": "витк./катушку",
           "nominal_a": "номинал, А", "percent": "% номинала", "status": "примечание"}


def save(rows, tables_dir: str = TB, name: str = OUT_TABLE) -> str:
    """Записать таблицу проектных токов (`;`, точка-десятичная — как остальные машинные таблицы БД)."""
    path = os.path.join(tables_dir, name)
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(FIELDS)
        for r in rows:
            w.writerow([r[k] for k in FIELDS])
    return path


def _print_table(rows) -> None:
    cols = [c for c in FIELDS]
    width = {c: max(len(HEADERS[c]), *(len(str(r[c])) for r in rows)) for c in cols}
    print("  ".join(HEADERS[c].ljust(width[c]) for c in cols))
    print("  ".join("-" * width[c] for c in cols))
    for r in rows:
        print("  ".join(str(r[c]).ljust(width[c]) for c in cols))


def _print_moments(mats) -> None:
    """Проверка: даёт ли решение проектные моменты (расхождение — признак вырожденной матрицы)."""
    bad = []
    for m in mats:
        if m.get("pm") or "achieved" not in m or not m.get("ctrl"):
            continue
        for i, mp in enumerate(m["outputs"]):
            d = m["design"].get(mp, 0.0)
            if abs(d) < 1e-30:
                continue
            rel = abs(m["achieved"][i] - d) / abs(d) * 100.0
            if rel > 1.0:
                bad.append(f"  {m['name']}: {mp} проект {d:+.4e} ≠ достигнуто {m['achieved'][i]:+.4e} ({rel:.1f}%)")
    if bad:
        print("\nмоменты, не воспроизведённые решением (проверьте матрицу связи):")
        print("\n".join(bad))
    else:
        print("\nвсе проектные моменты воспроизведены решением (расхождение < 1%).")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Проектные токи обмоток по матрицам связи из БД")
    ap.add_argument("--lattice", type=int, default=1, help="решётка, откуда берутся проектные mat_* (по умолч. 1)")
    ap.add_argument("--no-save", action="store_true", help="не записывать workspace/tables/design_currents.csv")
    ap.add_argument("--quiet", action="store_true", help="без таблицы на экран (только путь к файлу)")
    a = ap.parse_args(argv)

    warn = staleness()
    rows, mats = compute(a.lattice)
    if not a.quiet:
        _print_table(rows)
        _print_moments(mats)
    if warn:
        print("\n" + warn)
    if not a.no_save:
        print(f"\nтаблица проектных токов ({len(rows)} строк) -> {save(rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
