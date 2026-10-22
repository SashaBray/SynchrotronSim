# -*- coding: utf-8 -*-
"""Правка геометрии диполей решётки 26 под РЕАЛЬНОЕ ЖЕЛЕЗО (3D-модель), идемпотентно.

Зачем. Линейки длинных диполей (DL/JL) были заведены как 5 слайсов по 417.83 мм ВПЛОТНУЮ
(zero-gap «сплошное железо»). Физически каждый модуль короче — 410 мм, между ними зазор.
Короткий диполь (SBM) был двумя слайсами по 24.43 мм, а это один магнит 59 мм.

Что делает (только решётка 26, по месту, с бэкапами):
  • ДЛИННЫЙ модуль (kind=dipole, L≈0.41783): eff_length → 0.410, mat_B ×= old_L/0.410.
    Центр (x,y,z,yaw) НЕ трогается → зазор 417.83−410 = 7.83 мм открывается сам.
  • КОРОТКИЙ (sila26_SBM_1 + следующий sila26_SBM_2): сливаются в ОДИН 0.059-магнит,
    центр = середина пары, mat_B = суммарный изгиб / (QP·0.059), края = внешние (E1 первого,
    E2 второго). Строка SBM_2 удаляется, element_connectivity сшивается локально.
  • hardware_library.csv синхронизируется (eff_length длинных → 0.410, sila26_SBM_1 → 0.059),
    чтобы Python-тракт не расходился с C++ (C++ берёт длину из строки решётки).

Поле считается в БД как mat_B = −ANGLE/(QP·L), поэтому масштаб mat_B ∝ 1/L сохраняет
∫B·dl = угол изгиба → частица выходит на опорную орбиту. Скрипт печатает инвариант
Σ|mat_B·L| (∝ полному изгибу кольца) до/после: обязан совпасть.

Запуск: py -3.14 apply_dipole_iron_geometry.py   (идемпотентно: повтор — no-op)
"""
import csv, os, shutil, sys

sys.stdout.reconfigure(encoding="utf-8")
TB = os.path.join("workspace", "tables")
LC = os.path.join(TB, "lattice_configs.csv")
EC = os.path.join(TB, "element_connectivity.csv")
HW = os.path.join(TB, "hardware_library.csv")
LAT = "26"
LONG_IRON = 0.410           # м, железо длинного модуля DP_L
SHORT_IRON = 0.059          # м, железо короткого диполя DP_sh (один магнит вместо двух слайсов)
LONG_MIN = 0.30             # дискриминатор: dipole с L>0.30 — длинный модуль DL/JL
SHORT_MAX = 0.10            # dipole с L<0.10 — короткий (SBM)


def _read(path):
    rows = list(csv.reader(open(path, encoding="utf-8-sig"), delimiter=";"))
    return rows[0], rows[1:]


def _write(path, head, data):
    shutil.copy(path, path + ".bak_irongeom")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(head); w.writerows(data)


def _sum_bend(data, ci):
    """Σ|mat_B·L| по изгибающим (dipole/dq) решётки 26 — инвариант, ∝ полному углу."""
    s = 0.0
    for r in data:
        if r and r[ci["lattice_id"]] == LAT and r[ci["kind"]].strip() in ("dipole", "dq"):
            try:
                s += abs(float(r[ci["mat_B"]] or 0) * float(r[ci["eff_length"]] or 0))
            except ValueError:
                pass
    return s


def main():
    lh, ld = _read(LC); ci = {c: i for i, c in enumerate(lh)}
    before = _sum_bend(ld, ci)

    # идемпотентность: если длинных слайсов 0.41783 уже нет и SBM_2 уже нет — выходим
    has_long_old = any(r and r[ci["lattice_id"]] == LAT and r[ci["kind"]].strip() == "dipole"
                       and LONG_MIN < float(r[ci["eff_length"]] or 0) and abs(float(r[ci["eff_length"]]) - LONG_IRON) > 1e-4
                       for r in ld)
    has_sbm2 = any(r and r[ci["lattice_id"]] == LAT and r[ci["device_type_id"]].endswith("SBM_2") for r in ld)
    if not has_long_old and not has_sbm2:
        print("[no-op] правка уже применена (длинные 0.410, SBM слит)."); return

    # индекс по instance_id для сшивки пар SBM
    by_iid = {r[ci["instance_id"]]: r for r in ld if r and r[ci["lattice_id"]] == LAT}

    n_long = 0; n_merge = 0; del_iids = set()
    long_devs = set(); short_devs = set()
    out = []
    for r in ld:
        if not r or r[ci["lattice_id"]] != LAT:
            out.append(r); continue
        kind = r[ci["kind"]].strip()
        dev = r[ci["device_type_id"]]
        L = float(r[ci["eff_length"]] or 0)

        # ---- длинные модули DL/JL ----
        if kind == "dipole" and L > LONG_MIN:
            if abs(L - LONG_IRON) > 1e-6:
                scale = L / LONG_IRON
                r[ci["mat_B"]] = f"{float(r[ci['mat_B']] or 0) * scale:.15g}"
                r[ci["eff_length"]] = f"{LONG_IRON:.12f}"
                n_long += 1
            long_devs.add(dev); out.append(r); continue

        # ---- короткий: SBM_1 -> сливаем с SBM_2 ----
        if kind == "dipole" and L < SHORT_MAX and dev.endswith("SBM_1"):
            # партнёр SBM_2 = следующий инстанс
            iid2 = str(int(r[ci["instance_id"]]) + 1)
            r2 = by_iid.get(iid2)
            if r2 is None or not r2[ci["device_type_id"]].endswith("SBM_2"):
                print(f"  [warn] SBM_1 iid={r[ci['instance_id']]}: партнёр SBM_2 не найден — пропуск")
                out.append(r); continue
            b1 = float(r[ci["mat_B"]] or 0) * float(r[ci["eff_length"]] or 0)
            b2 = float(r2[ci["mat_B"]] or 0) * float(r2[ci["eff_length"]] or 0)
            matB_new = (b1 + b2) / SHORT_IRON                     # суммарный изгиб / новая длина
            # центр = середина пары (слайсы равной длины ⇒ среднее центров)
            for c in ("x", "y", "yaw"):
                r[ci[c]] = f"{0.5 * (float(r[ci[c]] or 0) + float(r2[ci[c]] or 0)):.12f}"
            r[ci["eff_length"]] = f"{SHORT_IRON:.12f}"
            r[ci["mat_B"]] = f"{matB_new:.15g}"
            r[ci["edge2"]] = r2[ci["edge2"]]                      # внешний край = E2 второго слайса
            del_iids.add(iid2)                                   # SBM_2 удаляем
            short_devs.add(dev); n_merge += 1
            out.append(r); continue

        # SBM_2 — удаляем (слит в SBM_1)
        if kind == "dipole" and L < SHORT_MAX and dev.endswith("SBM_2"):
            continue

        out.append(r)

    # удалить строки SBM_2 (по instance_id)
    out = [r for r in out if not (r and r[ci["lattice_id"]] == LAT and r[ci["instance_id"]] in del_iids)]
    after = _sum_bend(out, ci)
    _write(LC, lh, out)

    # ---- element_connectivity: сшить вокруг удалённых SBM_2 ----
    eh, ed = _read(EC); eci = {c: i for i, c in enumerate(eh)}
    ec_by_iid = {r[eci["instance_id"]]: r for r in ed if r and r[eci["lattice_id"]] == LAT}
    ec_out = []
    for r in ed:
        if r and r[eci["lattice_id"]] == LAT and r[eci["instance_id"]] in del_iids:
            continue                                             # строку SBM_2 убираем
        ec_out.append(r)
    for iid2 in del_iids:
        r2 = ec_by_iid.get(iid2)
        if not r2:
            continue
        prev, nxt = r2[eci["prev_instance_id"]], r2[eci["next_instance_id"]]  # prev=SBM_1, nxt=следующий
        if prev in ec_by_iid: ec_by_iid[prev][eci["next_instance_id"]] = nxt  # SBM_1.next -> следующий
        if nxt in ec_by_iid:  ec_by_iid[nxt][eci["prev_instance_id"]] = prev  # следующий.prev -> SBM_1
    _write(EC, eh, ec_out)

    # ---- hardware_library: синхронизировать eff_length ----
    hh, hd = _read(HW); hci = {c: i for i, c in enumerate(hh)}
    n_hw = 0
    for r in hd:
        if not r:
            continue
        dev = r[hci["device_type_id"]]
        if dev in long_devs and abs(float(r[hci["eff_length"]] or 0) - LONG_IRON) > 1e-6:
            r[hci["eff_length"]] = f"{LONG_IRON:.12f}"; n_hw += 1
        elif dev in short_devs:                                  # sila26_SBM_1 -> 0.059
            r[hci["eff_length"]] = f"{SHORT_IRON:.12f}"; n_hw += 1
    _write(HW, hh, hd)

    print(f"[OK] решётка 26:")
    print(f"  длинных модулей укорочено до {LONG_IRON*1e3:.0f} мм: {n_long}  (зазор {(0.41783-LONG_IRON)*1e3:.2f} мм)")
    print(f"  коротких диполей слито 2→1 ({SHORT_IRON*1e3:.0f} мм): {n_merge}  (удалено строк SBM_2: {len(del_iids)})")
    print(f"  hardware_library строк обновлено: {n_hw}")
    print(f"  ИНВАРИАНТ Σ|mat_B·L| (∝ полный изгиб): до={before:.9f}  после={after:.9f}  "
          f"Δ={abs(after-before):.2e}  {'OK' if abs(after-before) < 1e-9 else '!!! РАСХОЖДЕНИЕ'}")
    print(f"  бэкапы: *.bak_irongeom")


if __name__ == "__main__":
    main()
