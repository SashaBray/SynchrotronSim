# -*- coding: utf-8 -*-
"""Драйвер: сглаживание + симметризация базовых карт поля и раскладка по папкам в fields_bin/<тип>/
+ обновление ссылок в БД (hardware_library.map_bin_path).

СИММЕТРИЗАЦИЯ (проекция на симметрию мультиполя, symmetrize_map) — ТОЛЬКО квадруполи/секступоли/
октуполи, с учётом ЧИСЛА ПОЛЮСОВ и порядка симметрии: quad n=2 (4 полюса), sext n=3 (6), oct n=4 (8);
quad/oct разделяются грид-точным поворотом 90°, sext — отражениями. Комбинир.-функц. DQ мультипольно
НЕ симметризуется (только сглаживание). Постоянные диполи DP — только перенос (уже симметризованы).
Оригиналы карт остаются на месте (не трогаем); пишем ОЧИЩЕННЫЕ копии в подпапки.
Запуск: py -3.14 symmetrize_maps.py"""
import os, csv, shutil
import numpy as np
import symmetrize_map as sm
from field_quality import load

FB = os.path.join("workspace", "fields_bin")
TB = os.path.join("workspace", "tables")
SIGMA = 1.0
# (device, kind|None, подпапка, обрабатывать): kind=None → только сглаживание; обрабатывать=False → только перенос
PLAN = [
    ("QP_L_1A", "quad", "quadrupole", True), ("QP_sh_1A", "quad", "quadrupole", True),
    ("SP_I1_1A", "sext", "sextupole", True), ("SP_I2_1A", "sext", "sextupole", True),
    ("OP_1A", "oct", "octupole", True),
    ("DQ_I1_1A", None, "dipole_quadrupole", True), ("DQ_I2_1A", None, "dipole_quadrupole", True),
    ("DQ_I3_1A", None, "dipole_quadrupole", True),
    ("DP_L_PM", None, "dipole", False), ("DP_sh_PM", None, "dipole", False),
]


def midplane_by(h, B):
    """Паразитный медианный By на оси (должен быть ~0 у чистого мультиполя) — маркер артефакта."""
    jy = h["ny"] // 2; jz = h["nz"] // 2
    return float(np.mean(np.abs(B[:, jy, jz, 1]))) * 1e3   # мТл


def main():
    # ИСТОЧНИК — новый сырой .fld из fields_raw (map_raw_path), а не старый .bin; nominalArg из nominal_arg.
    hw = {r["device_type_id"]: r for r in csv.DictReader(open(os.path.join(TB, "hardware_library.csv"),
                                                              encoding="utf-8"), delimiter=";")}
    FR = os.path.join("workspace", "fields_raw")
    newpaths = {}
    print(f"{'элемент':12} {'тип':6} {'сетка':>12} {'By_ось до':>10} {'после':>9} мТл  → папка")
    for dev, kind, sub, proc in PLAN:
        raw = hw[dev]["map_raw_path"]; na = hw[dev]["nominal_arg"]
        na = float(na) if na not in ("", "-") else 1.0
        h, B = load(os.path.join(FR, raw))                 # ← новый .fld
        by0 = midplane_by(h, B)
        if not proc:                                       # только перенос (DP)
            SS = B
        elif kind is None:                                 # только сглаживание (DQ, комбинир.-функц.)
            SS = sm.smooth(B, SIGMA)
        else:                                              # сглаживание + симметризация (quad/sext/oct)
            S = sm.symmetrize(h, B.copy(), kind)
            SS = sm.symmetrize(h, sm.smooth(S, SIGMA), kind)
        by1 = midplane_by(h, SS)
        outdir = os.path.join(FB, sub); os.makedirs(outdir, exist_ok=True)
        outrel = f"{sub}/{os.path.splitext(os.path.basename(raw))[0]}.bin"   # имя .bin = имя исходного .fld
        sm.write_bin(h, SS, os.path.join(FB, outrel), nominal_arg=na)   # ← nominalArg=100
        newpaths[dev] = outrel
        grid = f"{h['nx']}x{h['ny']}x{h['nz']}"
        print(f"  {dev:12} {str(kind):6} {grid:>12} {by0:10.4f} {by1:9.4f}   {outrel}")
    # обновить hardware_library.map_bin_path
    p = os.path.join(TB, "hardware_library.csv")
    rows = list(csv.reader(open(p, encoding="utf-8"), delimiter=";")); hdr = rows[0]
    mi = hdr.index("map_bin_path"); di = hdr.index("device_type_id")
    shutil.copy(p, p + ".bak")
    for r in rows[1:]:
        if r[di] in newpaths:
            r[mi] = newpaths[r[di]]
    with open(p, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    print(f"\n[OK] {len(newpaths)} карт обработано → подпапки fields_bin/; map_bin_path обновлён (бэкап .bak)")


if __name__ == "__main__":
    main()
