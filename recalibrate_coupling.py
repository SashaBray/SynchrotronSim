# -*- coding: utf-8 -*-
"""Пересчёт матрицы связи coupling_matrix.csv по ТЕКУЩИМ картам fields_bin методом
ТРАЕКТОРНОГО матчинга (как исходная калибровка, но карты KD → fields_bin).

Для каждого элемента гоняем Boris-probe через поле и по СИММЕТРИИ отклонения выделяем
ТОЛЬКО активный мультиполь (element_transforms.active_coeffs), без паразитных кросс-членов:
  quad  k1:  per_G = -(u'(+a)-u'(-a)) / (2·QP·I·a·Leff)                 (нечёт)
  sext  k2:  per_S = -2·((u'(+a)+u'(-a))/2 - u'(0)) / (QP·I·a²·Leff)     (чёт − изгиб)
  oct   k3:  per_O = -6·((u'(+a)-u'(-a))/2) / (QP·I·a³·Leff)             (нечёт, без квада)
  dq  theta_y|k1|k2: per_B=-u'(0)/(QP·I·Leff), per_G и per_S как выше    (изгиб+квад+секступ)
  kicker kick_y: per_B=-u'(0)/(QP·I·Leff);  kick_z: per_Bs=+v'(0)/(QP·I·Leff)
Диполи DP (input_type=coeff) НЕ трогаем — они уже совпадают (осевой интеграл карты, ~10 мкрад).
Синтетические копии («копия эталона X») получают активные столбцы эталона X.
Запуск: py -3.14 recalibrate_coupling.py"""
import os, csv, shutil
import matrix_tracker as mt
import kickdoc_compare as kd
from bip_matrix import QP

TB = os.path.join("workspace", "tables")
HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
EFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}
CM = {r["device_type_id"]: r for r in mt._read("coupling_matrix.csv")}
TR = {r["device_type_id"]: r for r in mt._read("element_transforms.csv")}
# ток probe, А. КВАДРУПОЛИ пробуем при РАБОЧЕМ ~1А (новые карты 100А сильные: при 10А линза толстая,
# пучок проходит фокус → тонколинзовая формула per_G смещена вниз ~30%). Секступоль/октуполь не фокусируют
# линейно (нет k1) → толстолинзового смещения нет, оставляем 10А (нужен сигнал ∝a²,a³).
IPROBE = {"quad": 1.0, "sext": 10.0, "oct": 10.0, "dq": 10.0, "kicker": 20.0}
# полуразмах probe для выделения мультиполя: у октуполя сигнал ∝a³ тонет у оси → нужен a≥5мм
AMPL = {"quad": 2e-3, "sext": 3e-3, "oct": 6e-3, "dq": 2e-3, "kicker": 2e-3}


def defl(dev, val, u0, v0=0.0):
    field, half = kd.scaled_field(dev, val)
    *_ , up, vp = kd.boris(field, half, u0, v0)
    return up, vp


def calibrate(dev):
    """→ dict {столбец: значение} только для активных мультиполей элемента dev."""
    kind = TR[dev]["kind"]; act = TR[dev]["active_coeffs"]; Leff = EFF[dev]
    I = IPROBE.get(kind, 10.0); A = AMPL.get(kind, 2e-3); out = {}
    up0, vp0 = defl(dev, I, 0.0)
    upp, _ = defl(dev, I, +A); upm, _ = defl(dev, I, -A)
    odd = (upp - upm) / 2.0; even = (upp + upm) / 2.0
    if "k1" in act:
        out["per_G_Tlm"] = -odd / (QP * I * A * Leff)
    if "k2" in act:
        out["per_S_Tlm2"] = -2.0 * (even - up0) / (QP * I * A * A * Leff)
    if "k3" in act:
        out["per_O_Tlm3"] = -6.0 * odd / (QP * I * A * A * A * Leff)
    # изгиб: кик даёт +θ0 (без минуса, в отличие от −k1·x), поэтому per_B=+u'0/(QP·I·Leff)
    if "theta_y" in act and kind == "dq":
        out["per_B_Tl"] = up0 / (QP * I * Leff)
    if "kick_y" in act:
        out["per_B_Tl"] = up0 / (QP * I * Leff)
    if "kick_z" in act:
        # карта корр.I2 несёт то же Bz-поле (изгиб горизонтальный, u'0), в кольце развёрнут roll'ом
        # в вертикаль → сила вертикального кика = |ток→изгиб| из карты
        out["per_Bs_Tl"] = up0 / (QP * I * Leff)
    return out, dict(kind=kind, I=I, up0=up0, vp0=vp0)


def main():
    # корректоры (5-обмоточные) НЕ рекалибруем траекторным матчингом — их per_* из фита карт v3
    # (многообмоточные, отдельная калибровка). Старые corrector_I1/I2 удалены из библиотеки.
    base = ["QP_L_1A", "QP_sh_1A", "SP_I1_1A", "SP_I2_1A", "OP_1A",
            "DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A"]
    p = os.path.join(TB, "coupling_matrix.csv")
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter=";")); hdr = rows[0]
    idx = {c: hdr.index(c) for c in hdr}
    shutil.copy(p, p + ".bak")
    vals = {}
    print(f"{'элемент':16} {'старое':>13} {'новое':>13}  диагностика")
    for dev in base:
        newcols, diag = calibrate(dev)
        vals[dev] = newcols
        col = next(iter(newcols))                                   # главный активный столбец
        old = float(CM[dev].get(col) or 0)
        print(f"  {dev:14} {old:+.5e} {newcols[col]:+.5e}  {col} kind={diag['kind']} "
              f"u'0={diag['up0']*1e6:+.1f} v'0={diag['vp0']*1e6:+.1f} мкрад" +
              (f"  [+ещё {list(newcols)[1:]}]" if len(newcols) > 1 else ""))
    # применить к эталонам + синтетическим копиям
    note_i = idx["notes"]
    for r in rows[1:]:
        d = r[idx["device_type_id"]]; src = None
        if d in vals:
            src = d
        else:
            note = r[note_i]
            if "копия эталона" in note:
                cand = note.replace("копия эталона", "").strip()
                if cand in vals:
                    src = cand
        if src:
            for col, v in vals[src].items():
                r[idx[col]] = f"{v:.6e}"
    with open(p, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    print(f"[OK] coupling_matrix.csv пересчитан траекторным матчингом (бэкап .bak). DP не тронуты.")


if __name__ == "__main__":
    main()
