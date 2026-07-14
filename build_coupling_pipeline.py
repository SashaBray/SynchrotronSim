# -*- coding: utf-8 -*-
"""АВТО-ПАЙПЛАЙН матриц связи «ток → физкоэффициенты» для ВСЕХ элементов библиотеки — единая команда.

ПОЧЕМУ НЕ ОСЕВОЙ ФИТ. Коэф. мультиполей per_G/per_S/per_O НЕЛЬЗЯ надёжно снять осевым мультипольным
разложением карты: у оси секступоль (∝r²) и октуполь (∝r³) тонут в шуме карты и артефакте грубой
сетки (проверено этим же пайплайном: осевой фит расходится с эталоном ×2.7…×19). Достоверный способ —
ТРЕКИНГ-МАТЧИНГ: дифференциальное ядро (Boris) гоняет пробный пучок через карту при свипе тока, а
коэффициент матричной модели подгоняется так, чтобы воспроизвести траекторию. Это делает
match_coefficients.py → composite_fullrange.json (наша работа BIP). Диполи — исключение: крупное
низкопорядковое поле, осевой интеграл per_B=−(∫Bz)/L_eff надёжен (совпал с эталоном 0.0–0.2%),
поэтому берётся прямо И служит независимой кросс-проверкой пайплайна.

ИСТОЧНИКИ per-амперных коэф. по типам устройств:
  QP / SP / OP / DQ (база)  — трекинг-матчинг (composite_fullrange.json → build_matrix_db.coupling_rows)
  DP (диполи)               — осевой интеграл карты (валидирован), input_type=coeff
  корректор 5-обмоточный    — мультипольный фит карт v3 (build_coupling_table.corrector_matrix)
  СИНТЕТИЧЕСКИЕ элементы     — копия эталона scaled_from (per-амперные коэф. длино-независимы: тело
                              поля идентично, синтетик отличается только длиной тела/L_eff)
  legacy-корректоры I1/I2   — 0 (не входят в кольцо; ждут характеризации)

ВЫХОД: coupling_matrix.csv (плоский вид для трекера) + element_coupling.csv (структурный: обмотка=строка).
ОТЧЁТ: провенанс каждого значения + осевой фит как кросс-проверка (Δ% к записанному).

ЗАПУСК:
  python build_coupling_pipeline.py             # собрать таблицы из достоверных источников + отчёт
  python build_coupling_pipeline.py --dry-run    # только отчёт, без записи
  python build_coupling_pipeline.py --recompute   # СПЕРВА пересчитать composite трекинг-матчингом (долго)
"""
import os
import sys
import csv
import subprocess
import numpy as np
from field_quality import load
from dipole_composite_model import prof                       # комплексный мультипольный фит (валидация)
from build_matrix_db import coupling_rows                     # база: трекинг-матчинг (composite_fullrange)
from build_coupling_table import corrector_matrix, CORR_DEV   # корректор 5-обм.: фит карт v3

TB = os.path.join("workspace", "tables")
FB = os.path.join("workspace", "fields_bin")
FR = os.path.join("workspace", "fields_raw")
R1, R2 = 2e-3, 6e-3
FACT = {1: 1.0, 2: 1.0, 3: 2.0, 4: 6.0}
CM_HEAD = ["device_type_id", "input_type", "per_B_Tl", "per_G_Tlm", "per_S_Tlm2", "per_O_Tlm3",
           "per_Bs_Tl", "per_Gs_Tlm", "per_Ss_Tlm2", "per_Os_Tlm3", "notes"]
MULT8 = ["per_B_Tl", "per_Bs_Tl", "per_G_Tlm", "per_Gs_Tlm",
         "per_S_Tlm2", "per_Ss_Tlm2", "per_O_Tlm3", "per_Os_Tlm3"]
EC_HEAD = ["element_id", "winding", "device_type_id", "input_type"] + MULT8 + ["notes"]
NCOL = {1: ("per_B_Tl", "per_Bs_Tl"), 2: ("per_G_Tlm", "per_Gs_Tlm"),
        3: ("per_S_Tlm2", "per_Ss_Tlm2"), 4: ("per_O_Tlm3", "per_Os_Tlm3")}


def read_rows(name):
    return list(csv.DictReader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";"))


def blank_cm():
    return {c: "0.000000e+00" for c in CM_HEAD[2:-1]}


def locate_map(dev, hw_row):
    mb = (hw_row.get("map_bin_path") or "-").strip()
    if mb not in ("-", "") and os.path.exists(os.path.join(FB, mb)):
        return os.path.join(FB, mb)
    for cand in (os.path.join(FR, dev + ".fld"), os.path.join(FB, dev + ".bin")):
        if os.path.exists(cand):
            return cand
    return None


def axis_fit(path, L_eff, cur):
    """Осевой мультипольный фит карты → 8 per-амперных коэф. (надёжен для диполей; для мультиполей — кросс-проверка)."""
    h, B = load(path)
    B = np.nan_to_num(np.asarray(B, dtype=float))
    mir = 2.0 if h["mnx"] >= -1e-9 else 1.0                    # полукарта X≥0 → ×2
    xa, Cn = prof(h, B, [1, 2, 3, 4], R1, R2)
    per = {c: 0.0 for c in MULT8}
    for n in (1, 2, 3, 4):
        k = FACT[n] * mir / L_eff / cur
        re = float(np.trapezoid(np.nan_to_num(Cn[n].real), xa)) * k
        im = float(np.trapezoid(np.nan_to_num(Cn[n].imag), xa)) * k
        if n == 1:
            re = -re
        nc, sc = NCOL[n]; per[nc] = re; per[sc] = im
    return per


def assemble():
    """Собрать per-амперные коэф. всех устройств из достоверных источников. Возврат dict dev→(cm_row, provenance)."""
    hw = {r["device_type_id"]: r for r in read_rows("hardware_library.csv")}
    eff = {d: float(r["eff_length"]) for d, r in hw.items()}
    trans = {r["device_type_id"]: (r.get("kind") or "") for r in read_rows("element_transforms.csv")}
    out = {}                                                   # dev → {"row": {cm cols}, "src": str}

    # 1) База (QP/SP/OP/DQ/DP/legacy-корр.) — трекинг-матчинг + осевой интеграл диполей (coupling_rows)
    base = coupling_rows()
    for dev, d in base.items():
        row = blank_cm()
        row["per_B_Tl"] = f"{d['B']:.6e}"; row["per_G_Tlm"] = f"{d['G']:.6e}"
        row["per_S_Tlm2"] = f"{d['S']:.6e}"; row["per_O_Tlm3"] = f"{d['O']:.6e}"
        k = trans.get(dev, "")
        src = ("осевой интеграл карты" if k == "dipole" else
               "0 (ждёт характ.)" if dev.startswith("corrector_I") else
               "трекинг-матчинг (composite)")
        out[dev] = {"row": row, "itype": d["itype"], "src": src, "kind": k}

    # 2) Корректор 5-обмоточный — мультипольный фит карт v3
    corr = corrector_matrix()                                  # winding → [per_B,per_Bs,per_G,per_Gs,per_S,per_Ss,per_O,per_Os]
    for wind, dev in CORR_DEV.items():
        v = corr[wind]; row = blank_cm()
        for i, col in enumerate(MULT8):
            row[col] = f"{v[i]:.6e}"
        out[dev] = {"row": row, "itype": "current", "src": "фит карт v3 (5-обм.)",
                    "kind": trans.get(dev, "corrector")}

    # 3) Синтетические элементы — копия эталона scaled_from (per-амперные коэф. длино-независимы)
    for dev, r in hw.items():
        if r.get("synthetic") != "True":
            continue
        ref = r.get("scaled_from")
        if ref in out:
            row = dict(out[ref]["row"])
            out[dev] = {"row": row, "itype": out[ref]["itype"],
                        "src": f"копия эталона {ref}", "kind": out[ref]["kind"]}
    return out, hw, eff, trans


def report(out, hw, eff):
    """Отчёт: провенанс + осевой фит как кросс-проверка (для базовых устройств с картой)."""
    prim = {"quad": "per_G_Tlm", "sext": "per_S_Tlm2", "oct": "per_O_Tlm3",
            "dipole": "per_B_Tl", "dq": "per_G_Tlm", "corrector": "per_G_Tlm", "kicker": "per_G_Tlm"}
    print(f"{'устройство':22} {'kind':8} {'источник':26} {'гл.мультип':11} "
          f"{'значение':>13} {'осев.фит':>13} {'Δ%':>7}")
    for dev in sorted(out):
        if hw.get(dev, {}).get("synthetic") == "True":
            continue                                           # синтетики — в сводке ниже
        c = out[dev]; k = c["kind"]; pc = prim.get(k, "per_G_Tlm")
        val = float(c["row"][pc])
        path = locate_map(dev, hw[dev]); af = None
        if path:
            cur = 1.0 if c["itype"] == "coeff" else (100.0 if "_100A" in dev else 1.0)
            try:
                af = axis_fit(path, eff[dev], cur)[pc]
            except Exception:
                af = None
        if af is None:
            afs, dps = "—", ""
        else:
            dp = abs(val - af) / abs(val) * 100 if abs(val) > 1e-12 else (0 if abs(af) < 1e-12 else 999)
            afs = f"{af:+.5e}"; dps = f"{dp:7.1f}" + ("  ✓" if dp < 5 else "  (треков.)" if dp > 20 else "")
        print(f"{dev:22} {k:8} {c['src']:26} {pc:11} {val:+13.5e} {afs:>13} {dps}")
    nsyn = sum(1 for d in out if hw.get(d, {}).get("synthetic") == "True")
    print(f"\nСинтетических собрано (копия эталона): {nsyn}. Диполи с ✓ — осевой фит подтверждает "
          f"трекинг-матчинг; «(треков.)» — мультиполь, осевой фит недостоверен (ожидаемо).")


def write_tables(out):
    # coupling_matrix.csv (плоский)
    with open(os.path.join(TB, "coupling_matrix.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(CM_HEAD)
        for dev in _order(out):
            c = out[dev]; r = c["row"]
            w.writerow([dev, c["itype"]] + [r[col] for col in CM_HEAD[2:-1]] + [c["src"]])
    # element_coupling.csv (структурный)
    ec = []
    for dev in _order(out):
        c = out[dev]; elem, wind = _elem_winding(dev)
        ec.append([elem, wind, dev, c["itype"]] + [c["row"][col] for col in MULT8] + [c["src"]])
    ec.sort(key=lambda r: (r[0], r[1]))
    with open(os.path.join(TB, "element_coupling.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(EC_HEAD); w.writerows(ec)
    print(f"\n[OK] coupling_matrix.csv: {len(out)} устройств; "
          f"element_coupling.csv: {len(ec)} обмоток, {len(set(r[0] for r in ec))} элементов.")


def _order(out):
    base = ["QP_L_1A", "QP_sh_1A", "SP_I1_1A", "SP_I2_1A", "OP_1A", "DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A",
            "DP_L_PM", "DP_sh_PM", "corrector_I1_1A", "corrector_I2_1A"]
    corr = list(CORR_DEV.values())
    head = [d for d in base + corr if d in out]
    syn = [d for d in out if d not in head]
    return head + syn


def _elem_winding(dev):
    for w in ("_I1", "_I2", "_I3"):
        if w in dev:
            return dev.replace(w, ""), w[1:]
    if dev.startswith("corrector_") and "_100A" in dev:
        return "corrector", dev.split("_")[1]
    return dev, "main"


def main(dry=False, recompute=False):
    if recompute:
        print("[recompute] пересчёт composite_fullrange.json трекинг-матчингом (match_coefficients.py)…")
        subprocess.run([sys.executable, "match_coefficients.py"], check=False)
    out, hw, eff, trans = assemble()
    report(out, hw, eff)
    if not dry:
        write_tables(out)
    else:
        print("\n[dry-run] таблицы не записаны.")


if __name__ == "__main__":
    main(dry="--dry-run" in sys.argv, recompute="--recompute" in sys.argv)
