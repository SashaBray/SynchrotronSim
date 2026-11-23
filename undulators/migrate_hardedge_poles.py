# -*- coding: utf-8 -*-
"""Миграция БД: полюса ондуляторов → эквивалент «жёсткий край», отдельная линза UFC убрана.

ЗАЧЕМ. Раньше полюс = ступенька во весь полупериод (λu/2, поле (2/π)·B0): угол верный, а ∫B² — 81% синусоиды,
поэтому естественную вертикальную фокусировку ½(B0/Bρ)² нёс отдельный элемент UFC (kind=quad) на выходе устройства
(он же гасил горизонтальную слабую фокусировку секторных полюсов). Теперь фокусировка — свойство самих полюсов:

  длина полюса  L' = (∫B)²/∫B² = (8/π²)·λu/2   (полупериод 9 мм → 7,30 мм, промежутки 1,70 мм),
  поле          B' = (π/4)·B0                    (угол θ за полупериод — прежний: mat_B ×= L/L'),
  грани         E1 = E2 = θ/2                    (прямоугольный полюс: горизонталь — дрейф, вертикаль — ∫B²/(Bρ)²).

Совпадают сразу ∫B (орбита, K, замыкание I1=I2=0 — центры полюсов не двигаются), ∫B² (вертикальная фокусировка
½(B0/Bρ)² и мощность СИ) и нулевая горизонтальная фокусировка. Знак граней: в ядре край — py −= tan(E)·(θ0/L)·y,
вертикальный ФОКУС при tan(E)·θ0 > 0 ⇒ E того же знака, что θ0 = QP·mat_B·L.

ЧТО ДЕЛАЕТ (решётки 26 и 27 — те, где стоят ондуляторы): строки полюсов und_UPP/UPM/UE0/UE1 — eff_length, mat_B,
edge1/edge2; строки und_UFC удаляются, инстансы перенумеровываются 1..M (build_line упорядочивает по instance_id),
element_connectivity перестраивается цепочкой, element_groups — без UFC, опорная точка = среднее центров полюсов
(центр устройства; раньше UFC смещал её ~4,6 мм к выходу). Перед записью — резервная копия трёх таблиц в
workspace/tables/backups/hardedge_<время>/. Повторный запуск ничего не меняет (полюса уже короткие, UFC нет).

Точечная миграция, а не переимпорт: полный переимпорт решётки 26 из .lte затёр бы правки БД, сделанные после него.
Генераторы (build_undulator_lte/assemble_undulators/import_undulators/insert_undulators) дают ту же модель; одно
отличие — при новой сборке устройство центрируется в промежутке без бывших 2 см UFC, т.е. на 1 см дальше по s,
чем в мигрированной БД (там центры полюсов сохранены как были). На оптику и орбиту это не влияет.

ЗАПУСК
  py -3.14 undulators/migrate_hardedge_poles.py              # решётки 26 и 27
  py -3.14 undulators/migrate_hardedge_poles.py --dry-run    # только отчёт
"""
import os, sys, csv, math, time, shutil, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
LC = os.path.join(TB, "lattice_configs.csv")
CN = os.path.join(TB, "element_connectivity.csv")
EG = os.path.join(TB, "element_groups.csv")
POLE_DEVS = {"und_UPP", "und_UPM", "und_UE0", "und_UE1"}
UFC_DEV = "und_UFC"
FRAC = 8.0 / math.pi ** 2                       # L'/(λu/2) эквивалента «жёсткий край»


def _read(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    return rows[0], rows[1:]


def _write(path, hdr, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(rows)


def migrate(lats=(26, 27), dry=False):
    hdr, rows = _read(LC)
    ci = {c: i for i, c in enumerate(hdr)}
    li, ii, di = ci["lattice_id"], ci["instance_id"], ci["device_type_id"]
    report, idmap = {}, {}
    new_lat_rows = {}
    for lat in map(str, lats):
        mine = sorted((r for r in rows if r[li] == lat), key=lambda r: int(r[ii]))
        poles = [r for r in mine if r[di] in POLE_DEVS]
        n_ufc = sum(r[di] == UFC_DEV for r in mine)
        if not poles:
            report[lat] = "ондуляторов нет"; continue
        L0s = {round(float(r[ci["eff_length"]]), 9) for r in poles}
        done = all(abs(L - FRAC * 0.009) < 1e-9 for L in L0s)
        if done and not n_ufc:
            report[lat] = "уже мигрирована"; continue
        if done or len(L0s) != 1 or abs(next(iter(L0s)) - 0.009) > 1e-9:
            raise SystemExit(f"[ОТКАЗ] решётка {lat}: длины полюсов {sorted(L0s)} — ожидались 0,009 м (λu/2)")
        out, m = [], {}
        for r in mine:
            if r[di] == UFC_DEV:
                continue
            r = list(r)
            if r[di] in POLE_DEVS:
                L0 = float(r[ci["eff_length"]]); mB = float(r[ci["mat_B"]] or 0.0)
                th = mt.QP * mB * L0                               # угол полюса (сохраняется)
                Lp = FRAC * L0
                r[ci["eff_length"]] = f"{Lp:.12f}"
                r[ci["mat_B"]] = f"{th / (mt.QP * Lp):.15g}"
                r[ci["edge1"]] = r[ci["edge2"]] = f"{0.5 * th:.15g}"
            m[r[ii]] = str(len(out) + 1)
            r[ii] = m[r[ii]]
            out.append(r)
        idmap[lat] = m; new_lat_rows[lat] = out
        report[lat] = f"полюсов {len(poles)} (L {L0s.pop()*1e3:.3f} → {FRAC*9:.4f} мм, грани θ/2), UFC удалено {n_ufc}, инстансов {len(mine)} → {len(out)}"
    if dry or not new_lat_rows:
        return report, None
    # --- резервная копия
    bdir = os.path.join(TB, "backups", "hardedge_" + time.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(bdir, exist_ok=True)
    for p in (LC, CN, EG):
        shutil.copy2(p, bdir)
    # --- lattice_configs: строки мигрируемых решёток на прежнее место (по первой строке решётки)
    res, placed = [], set()
    for r in rows:
        lat = r[li]
        if lat in new_lat_rows:
            if lat not in placed:
                res.extend(new_lat_rows[lat]); placed.add(lat)
            continue
        res.append(r)
    _write(LC, hdr, res)
    # --- connectivity: цепочка по новым номерам
    chdr, crows = _read(CN)
    cci = {c: i for i, c in enumerate(chdr)}
    cres, placed = [], set()
    for r in crows:
        lat = r[cci["lattice_id"]]
        if lat in new_lat_rows:
            if lat not in placed:
                lr = new_lat_rows[lat]; n = len(lr)
                for k, x in enumerate(lr, 1):
                    row = [""] * len(chdr)
                    row[cci["lattice_id"]] = lat; row[cci["instance_id"]] = str(k)
                    row[cci["device_type_id"]] = x[di]
                    row[cci["prev_instance_id"]] = str(k - 1) if k > 1 else ""
                    row[cci["next_instance_id"]] = str(k + 1) if k < n else ""
                    cres.append(row)
                placed.add(lat)
            continue
        cres.append(r)
    _write(CN, chdr, cres)
    # --- element_groups: без UFC, новые номера, опорная точка = среднее центров полюсов
    ghdr, grows = _read(EG)
    gci = {c: i for i, c in enumerate(ghdr)}
    for r in grows:
        lat = r[gci["lattice_id"]]
        if lat not in idmap:
            continue
        old = r[gci["instance_ids"]].split("|")
        new = [idmap[lat][o] for o in old if o in idmap[lat]]
        byid = {x[ii]: x for x in new_lat_rows[lat]}
        mem = [byid[n] for n in new if byid[n][di] in POLE_DEVS]
        r[gci["instance_ids"]] = "|".join(new)
        for ax, col in (("x", "pivot_x"), ("y", "pivot_y"), ("z", "pivot_z")):
            r[gci[col]] = f"{sum(float(x[ci[ax]]) for x in mem) / len(mem):.12f}"
        mid = mem[len(mem) // 2]
        r[gci["yaw"]], r[gci["pitch"]], r[gci["roll"]] = mid[ci["yaw"]], mid[ci["pitch"]], mid[ci["roll"]]
    _write(EG, ghdr, grows)
    return report, bdir


def main():
    ap = argparse.ArgumentParser(description="полюса ИД → эквивалент «жёсткий край», UFC убран")
    ap.add_argument("--lattice", type=int, nargs="*", default=[26, 27])
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    rep, bdir = migrate(a.lattice, a.dry_run)
    for lat, s in rep.items():
        print(f"решётка {lat}: {s}")
    if bdir:
        print(f"резервная копия: {bdir}")


if __name__ == "__main__":
    main()
