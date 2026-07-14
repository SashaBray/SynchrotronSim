# -*- coding: utf-8 -*-
"""Генератор ОТДЕЛЬНОЙ таблицы БД матриц связи элементов библиотеки: workspace/tables/element_coupling.csv.
Строка = ОБМОТКА элемента (столбец матрицы связи), колонки = 8 мультипольных выходов на 1А:
  per_B(θx-диполь) per_Bs(θy-верт.диполь) per_G(норм.квад) per_Gs(скью-квад)
  per_S(норм.секступоль) per_Ss(скью-секступоль) per_O(норм.октуполь) per_Os(скью-октуполь).
Полная матрица связи элемента = набор его строк (обмоток) × 8 мультиполей.

Многообмоточные (КОРРЕКТОР 5 обмоток) считаются из карт: комплексные мультиполи c_n(x)=b_n+i·a_n
(prof), интеграл по элементу (зеркало ×2) → эффективный коэф. поля / L_eff / 100А. Множитель (n−1)!:
mat_S=∂²Bz/∂u²=2·b3, mat_O=∂³Bz/∂u³=6·b4. Одиночные (QP/SP/OP/DQ/DP) — из coupling_matrix.csv (скью=0).
Также ДОПИСЫВАЕТ новые обмотки корректора в coupling_matrix.csv (плоский вид для трекера).
Выход: element_coupling.csv (+ обновл. coupling_matrix.csv)."""
import os, csv
import numpy as np
import corrector_experiment as ce
from dipole_composite_model import prof

TB = os.path.join("workspace", "tables")
QP = ce.QP
R1, R2 = 2e-3, 6e-3                                        # радиус подгонки мультиполей (зона пучка ±5мм)
LEFF_CORR = 0.1615                                        # эфф. длина корректора (hardware_library)
FACT = {1: 1.0, 2: 1.0, 3: 2.0, 4: 6.0}                  # (n−1)!: ∂ⁿ⁻¹Bz/∂uⁿ⁻¹ = (n−1)!·bₙ
MULT = ["per_B_Tl", "per_Bs_Tl", "per_G_Tlm", "per_Gs_Tlm",
        "per_S_Tlm2", "per_Ss_Tlm2", "per_O_Tlm3", "per_Os_Tlm3"]
HEAD = ["element_id", "winding", "device_type_id", "input_type"] + MULT + ["notes"]

# device_type → (element_id, winding). Многообмоточные группируются под одним element_id.
EW = {
    "QP_L_1A": ("QP_L", "main"), "QP_sh_1A": ("QP_sh", "main"),
    "SP_I1_1A": ("SP", "I1"), "SP_I2_1A": ("SP", "I2"), "OP_1A": ("OP", "main"),
    "DQ_I1_1A": ("DQ", "I1"), "DQ_I2_1A": ("DQ", "I2"), "DQ_I3_1A": ("DQ", "I3"),
    "DP_L_PM": ("DP_L", "main"), "DP_sh_PM": ("DP_sh", "main"),
    "corrector_I1_1A": ("corrector_legacy", "I1"), "corrector_I2_1A": ("corrector_legacy", "I2"),
    "corrector_main_100A": ("corrector", "main"), "corrector_c1_100A": ("corrector", "c1"),
    "corrector_c2_100A": ("corrector", "c2"), "corrector_c3_100A": ("corrector", "c3"),
    "corrector_c4_100A": ("corrector", "c4"),
}
CORR_DEV = {"main": "corrector_main_100A", "c1": "corrector_c1_100A", "c2": "corrector_c2_100A",
            "c3": "corrector_c3_100A", "c4": "corrector_c4_100A"}
SHORT = {"main": "main", "correct_1": "c1", "correct_2": "c2", "correct_3": "c3", "correct_4": "c4"}


def corrector_matrix():
    """5×8 матрица связи корректора из карт обмоток: per-амперные мультипольные коэф. поля."""
    hs, Bs = zip(*[ce.load_cached(m) for m in ce.MAPS]); h = hs[0]
    rows = {}
    for j, w in enumerate(ce.WINDINGS):
        xa, Cn = prof(h, Bs[j], [1, 2, 3, 4], R1, R2)      # c_n(x)=b_n+i·a_n, x≥0
        vals = []
        for n in (1, 2, 3, 4):
            k = FACT[n] * 2.0 / LEFF_CORR / 100.0          # (n−1)!·(зеркало ×2)/L_eff/100А
            re = np.trapezoid(Cn[n].real, xa) * k          # нормальный b_n
            im = np.trapezoid(Cn[n].imag, xa) * k          # скью a_n
            if n == 1:
                re = -re                                   # mat_B=−(∫Bz)/L_eff: матр.конвенция px+=+QP·B, но физ. du'=−QP·Bz
            vals.append(re); vals.append(im)
        rows[SHORT[w]] = vals                              # [per_B,per_Bs,per_G,per_Gs,per_S,per_Ss,per_O,per_Os]
    return rows


def load_existing():
    """Одиночные элементы из coupling_matrix.csv → per_* (скью читаются, если есть, иначе 0)."""
    p = os.path.join(TB, "coupling_matrix.csv")
    out = {}
    for r in csv.DictReader(open(p, encoding="utf-8"), delimiter=";"):
        out[r["device_type_id"]] = r
    return out


def main():
    corr = corrector_matrix()
    existing = load_existing()

    # sanity: θx/θy/k1 корректора должны совпасть с corrector_model (88.3мрад main и т.д.)
    print("Матрица связи корректора (проверка, на 100А):")
    print(f"{'обмотка':>8} {'θx,мрад':>9} {'θy,мрад':>9} {'k1·L':>8} {'k2·L':>9} {'k3·L':>9}")
    for w in ce.WINDINGS:
        v = corr[SHORT[w]]
        thx = QP * v[0] * 100 * LEFF_CORR * 1e3; thy = QP * v[1] * 100 * LEFF_CORR * 1e3
        k1L = QP * v[2] * 100 * LEFF_CORR; k2L = QP * v[4] * 100 * LEFF_CORR; k3L = QP * v[6] * 100 * LEFF_CORR
        print(f"{w:>8} {thx:>9.2f} {thy:>9.2f} {k1L:>8.3f} {k2L:>9.4f} {k3L:>9.4f}")

    # element_coupling.csv: все элементы, обмотка = строка
    rows = []
    order = ["QP_L", "QP_sh", "SP", "OP", "DQ", "DP_L", "DP_sh", "corrector", "corrector_legacy"]
    seen = set()
    for dev, (elem, wind) in EW.items():
        if elem == "corrector":                            # НОВЫЙ корректор — из карт
            v = corr[wind]; itype = "current"
            note = "5-обмоточный, характеризован из карт v3"
            per = {MULT[i]: f"{v[i]:.6e}" for i in range(8)}
        else:                                              # из coupling_matrix.csv (скью 0)
            c = existing.get(dev, {})
            itype = c.get("input_type", "current")
            note = c.get("notes", "") or ""
            per = {m: (c.get(m) or "0.000000e+00") for m in MULT}
        row = {"element_id": elem, "winding": wind, "device_type_id": dev, "input_type": itype, "notes": note}
        row.update(per); rows.append(row); seen.add(dev)
    rows.sort(key=lambda r: (order.index(r["element_id"]) if r["element_id"] in order else 99, r["winding"]))
    p = os.path.join(TB, "element_coupling.csv")
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=HEAD, delimiter=";"); w.writeheader(); w.writerows(rows)
    print(f"\n[OK] {p}: {len(rows)} строк (обмоток), {len(set(r['element_id'] for r in rows))} элементов")

    # coupling_matrix.csv: дописать 5 обмоток корректора (плоский вид для трекера)
    cp = os.path.join(TB, "coupling_matrix.csv")
    cm = list(csv.reader(open(cp, encoding="utf-8"), delimiter=";")); head = cm[0]
    have = {r[0] for r in cm[1:]}
    idx = {name: head.index(name) for name in MULT}
    for wind, dev in CORR_DEV.items():
        if dev in have:
            continue
        v = corr[wind]; r = ["0.000000e+00"] * len(head)
        r[0] = dev; r[1] = "current"; r[head.index("notes")] = "корректор 5-обм."
        for i, m in enumerate(MULT):
            r[idx[m]] = f"{v[i]:.6e}"
        cm.append(r)
    with open(cp, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(cm)
    print(f"[OK] {cp}: +{len(CORR_DEV)} обмоток корректора (плоский вид для трекера)")


if __name__ == "__main__":
    main()
