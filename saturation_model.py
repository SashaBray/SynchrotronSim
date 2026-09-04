# -*- coding: utf-8 -*-
"""Графовая модель насыщения по элементам библиотеки: подгонка по траекториям + Word-отчёт.

Построена на `element_compare.py`: оттуда берутся диффурный прогон, киковая модель и подгонка
физкоэффициентов СРАВНЕНИЕМ ТРАЕКТОРИЙ (`fit_element` минимизирует расхождение траекторий диффура и
кика). Отличие от `element_compare`: там строится ЛИНЕЙНАЯ матрица связи «ток → физкоэф», здесь —
нелинейная модель насыщения на графе ветвей.

Что делает:

1. Читает реестр карт `workspace/tables/field_maps.csv` (колонка `split`: `train` — по три карты на
   обмотку для подгонки, `test` — всё остальное, в подгонке НЕ участвуют).
2. На каждую карту заводит ВРЕМЕННОЕ устройство в `hardware_library` (карта берётся как есть, без
   масштабирования по току — в этом и смысл: каждая карта решена ANSYS при своём возбуждении),
   конвертирует `.fld` → `.bin` и гоняет диффур.
3. Подгоняет физкоэффициенты по траекториям — получает точки «ампер-витки → коэффициент».
4. По этим точкам подгоняет модель насыщения: на ветвь два числа (сечение и длина железа),
   `NI = A·y + L·H(C·y)`; для многообмоточных — граф ветвей.
5. Проверяет на комбинационных картах и собирает Word-отчёт.

    py -3.14 saturation_model.py --list                 # что есть в реестре, чего не хватает
    py -3.14 saturation_model.py --element QP_sh        # прогнать один элемент
    py -3.14 saturation_model.py --all                  # прогнать всё (долго)
    py -3.14 saturation_model.py --report               # собрать отчёт по готовым датасетам

Датасеты кешируются в `workspace/results/saturation_model/`, поэтому `--report` не пересчитывает
физику. Таблицы БД правятся только на время прогона и восстанавливаются (как в `element_compare`).
"""
from __future__ import annotations

import argparse
import csv
import json
import os

REGISTRY = os.path.join("workspace", "tables", "field_maps.csv")
OUTDIR = os.path.join("workspace", "results", "saturation_model")
MATERIAL = os.path.join("workspace", "materials", "steel_1008_bh.csv")

# сколько карт на обмотку нужно для подгонки: две точки кривой железа + линейный наклон
MAPS_PER_WINDING = 3

# вид магнита для element_compare (field_type в БД — направление поля, не вид)
KIND = {"QP_sh": "quad", "QP_L": "quad", "SP": "sext", "OP": "oct",
        "DQ": "dq", "corrector": "kicker"}

# Радиус пучка/проб по типу элемента. Секступоль даёт отклонение ∝r², октуполь ∝r³ — на круге 1 мм
# (годном для квадруполя) они тонут в шуме карты, и подгонка S/O скачет. Радиус берётся с запасом
# внутри области хорошего поля (GOOD_FIELD_MM).
PROBE_R = {"quad": 1e-3, "dq": 1e-3, "kicker": 1e-3, "sext": 5e-3, "oct": 8e-3}

BRANCHES = os.path.join("workspace", "tables", "saturation_branches.csv")
MEMBERS = os.path.join("workspace", "tables", "saturation_branch_members.csv")


def read_graph(branches=BRANCHES, members=MEMBERS):
    """Граф ветвей из БД: ветви элемента, размещение обмоток, материал и геометрия железа.

    Возврат {элемент: {branches, unwound, coils, turns, material, geom, edge}}. `geom` — заданные по
    чертежу сечение и длина железа ветви; пусто ⇒ подбираются. `edge[ветвь] = (узел, узел)` — топология
    магнитной цепи. `coils[обмотка]` — СПИСОК ветвей: одна обмотка может быть парой катушек на разных
    полюсах, соединённых последовательно (один ток, два ребра) — так устроен DQ."""
    def _rows(p):
        if not os.path.exists(p):
            raise SystemExit(f"нет таблицы {p} — граф ветвей задаётся в БД")
        with open(p, encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter=";"))
    g = {}
    for r in _rows(branches):
        d = g.setdefault(r["element"], dict(branches=[], unwound=[], coils={}, turns={},
                                            material={}, geom={}, edge={}, kind={}, sym={},
                                            mirror={}, pos={}, share={}, xy={}))
        d["branches"].append(r["branch"])
        d["material"][r["branch"]] = r["material"] or "steel_1008"
        d["edge"][r["branch"]] = ((r.get("node_from") or "").strip(), (r.get("node_to") or "").strip())
        d["kind"][r["branch"]] = (r.get("kind") or "iron").strip() or "iron"
        d["sym"][r["branch"]] = (r.get("sym_class") or r["branch"]).strip() or r["branch"]
        d["mirror"][r["branch"]] = float((r.get("mirror") or "0").strip() or 0)
        d["pos"][r["branch"]] = (float((r.get("pos_from") or "0").strip() or 0),
                                 float((r.get("pos_to") or "0").strip() or 0))
        # Координаты концов ребра в поперечной плоскости (мм от оси пучка), если заданы в БД:
        # pos_* — горизонталь, z_* — вертикаль. Нужны для направлений трубок потока и мультиполей.
        za, zb = (r.get("z_from") or "").strip(), (r.get("z_to") or "").strip()
        if za or zb:
            d["xy"][r["branch"]] = ((d["pos"][r["branch"]][0], float(za or 0)),
                                    (d["pos"][r["branch"]][1], float(zb or 0)))
        S, L = (r.get("S_m2") or "").strip(), (r.get("len_m") or "").strip()
        if S and L:
            d["geom"][r["branch"]] = (float(S.replace(",", ".")), float(L.replace(",", ".")))
    for r in _rows(members):
        d = g.get(r["element"])
        if d is None:
            raise SystemExit(f"обмотка {r['winding']}: элемента {r['element']} нет в {branches}")
        if r["branch"] not in d["branches"]:
            raise SystemExit(f"{r['element']}: ветви {r['branch']} нет в {branches}")
        d["coils"].setdefault(r["winding"], []).append(r["branch"])
        # Доля МДС обмотки на этом ребре. У DQ и main это РАЗНЫЕ катушки — по 1.0 на каждой (МДС
        # складывается). У обмотки, размазанной по участкам одного корня, доля 1/N: катушка одна.
        sh = (r.get("mmf_share") or "").strip()
        d["share"][(r["winding"], r["branch"])] = float(sh.replace(",", ".")) if sh else 1.0
        t = (r.get("turns_per_coil") or "").strip()
        if t:
            d["turns"][r["winding"]] = int(float(t))
    for el, d in g.items():                                # ветви без обмоток — производная величина
        wound = {b for bs in d["coils"].values() for b in bs}
        d["unwound"] = [b for b in d["branches"] if b not in wound]
    return g


BRANCH_GRAPH = read_graph() if os.path.exists(BRANCHES) else {}


# ─────────────────────────────── реестр карт ────────────────────────────────

def read_registry(path=REGISTRY):
    """Реестр карт модели: список словарей."""
    if not os.path.exists(path):
        raise SystemExit(f"нет реестра карт {path} — сначала разложите карты по папкам")
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def inventory(rows):
    """Свод по элементам: сколько карт на обмотку и чего не хватает."""
    inv = {}
    for r in rows:
        if r["split"] != "train":
            continue
        inv.setdefault(r["element"], {}).setdefault(r["winding"], []).append(int(r["amp_turns"]))
    for el in inv:
        for w in inv[el]:
            inv[el][w] = sorted(set(inv[el][w]))
    return inv


def missing(inv):
    """Чего не хватает до трёх карт на обмотку."""
    out = []
    for el in sorted(inv):
        for w in sorted(inv[el]):
            lv = inv[el][w]
            if len(lv) < MAPS_PER_WINDING:
                out.append((el, w, MAPS_PER_WINDING - len(lv), lv))
    return out


def cmd_list(args):
    rows = read_registry()
    inv = inventory(rows)
    tests = {}
    for r in rows:
        if r["split"] == "test":
            tests.setdefault(r["element"], set()).add(r["excitation"])
    print(f"{'элемент':11} {'обмотка':8} {'карт':>5}  уровни, ампер-витки")
    for el in sorted(inv):
        for w in sorted(inv[el]):
            lv = inv[el][w]
            flag = "" if len(lv) >= MAPS_PER_WINDING else f"   ← нужно ещё {MAPS_PER_WINDING - len(lv)}"
            print(f"{el:11} {w:8} {len(lv):5}  {lv}{flag}")
        g = BRANCH_GRAPH.get(el)
        extra = f"  ветвей {len(g['branches'])}, без обмоток {len(g['unwound'])}" if g else "  ГРАФ НЕ ЗАДАН"
        print(f"{'':11} {'проверочных комбинаций':30} {len(tests.get(el, ())):3}{extra}\n")
    bad = []                                        # сверка графа из БД с реестром карт
    for el, winds in inv.items():
        g = BRANCH_GRAPH.get(el)
        if not g:
            bad.append(f"{el}: нет в {BRANCHES}")
            continue
        for w in winds:
            if w not in g["coils"]:
                bad.append(f"{el}/{w}: карты есть, а в {MEMBERS} обмотки нет")
        for w in g["coils"]:
            if w not in winds:
                bad.append(f"{el}/{w}: в графе есть, а карт нет")
    if bad:
        print("РАСХОЖДЕНИЯ ГРАФА И РЕЕСТРА:")
        for b in bad:
            print("   " + b)
        print()
    miss = missing(inv)
    if miss:
        print("НЕ ХВАТАЕТ КАРТ:")
        for el, w, k, lv in miss:
            print(f"   {el} / {w}: есть {lv}, снять ещё {k}")
    else:
        print("карт хватает на все обмотки")
    return 0


# ─────────────────────── физика: диффур на каждой карте ──────────────────────

def _temp_device(ec, base_dev, rec):
    """Временное устройство под одну карту: та же геометрия, своя карта, nominal_arg = ампер-витки.

    Карта берётся как есть — масштабирования по току нет, поэтому насыщение в ней уже учтено.
    Пути в `hardware_library` ОТНОСИТЕЛЬНЫЕ: ядро ищет их от `workspace/fields_raw/` и
    `workspace/fields_bin/` (см. src/main.cpp)."""
    name = "_satmodel_" + rec["map_id"].replace("+", "p").replace(".", "_")
    row = dict(ec.HW[base_dev])
    row["device_type_id"] = name
    row["map_raw_path"] = rec["raw_path"]
    row["map_bin_path"] = rec["bin_path"]
    row["nominal_arg"] = str(abs(int(rec["amp_turns"]))) if rec["amp_turns"] else "1"
    row["synthetic"] = "False"
    return name, row


def convert_map(rec):
    """`.fld` → `.bin` ТОЛЬКО для этой карты (штатный `--convert` пересобирает всю библиотеку)."""
    import numpy as np
    import pandas as pd
    from synthetic_assets import _write_bin
    dst = os.path.join("workspace", "fields_bin", rec["bin_path"])
    if os.path.exists(dst):
        return dst
    src = os.path.join("workspace", "fields_raw", rec["raw_path"])
    d = pd.read_csv(src, sep=r"\s+", skiprows=2, header=None,
                    names=list("xyzabc"), engine="c").to_numpy(float)
    xs, ys, zs = (np.unique(d[:, i]) for i in range(3))
    B = d[np.lexsort((d[:, 2], d[:, 1], d[:, 0])), 3:].reshape(len(xs), len(ys), len(zs), 3)
    h = dict(mnx=xs[0], mny=ys[0], mnz=zs[0], mxx=xs[-1],
             sx=xs[1] - xs[0], sy=ys[1] - ys[0], sz=zs[1] - zs[0],
             nx=len(xs), ny=len(ys), nz=len(zs))
    _write_bin(dst, h, B, nominal_arg=float(row_nominal(rec)))
    return dst


def row_nominal(rec):
    """Ток снятия карты в тех же единицах, в каких потом задаётся arg_val."""
    return abs(int(rec["amp_turns"])) if rec["amp_turns"] else 1


DIFFUR_TIMEOUT = 180        # с. Жёсткий предел на один прогон ядра: зависшее ядро съедает всю память
GOOD_FIELD_MM = 20.0        # мм. Радиус области хорошего поля: дальше разложение поля не описывает
MAGNET_HALF_M = 0.15        # м. Полудлина области магнита — где и меряется отклонение


def screen_trajectories(ec, trajs):
    """Годится ли карта для подгонки ПО ТРАЕКТОРИЯМ.

    Киковая модель считает элемент тонким и не знает, что частица покинула апертуру. Если при
    возбуждении карты пучок разлетается или разворачивается, две модели несравнимы в принципе —
    оптимизатору нечего сводить. Возврат (годится, макс |r| внутри магнита в мм, есть ли разворот)."""
    import numpy as np
    r_in, reversed_ = 0.0, False
    for t in trajs:
        m = np.abs(t[:, 0] - ec.CENTER) < MAGNET_HALF_M
        if m.any():
            r_in = max(r_in, 1000.0 * float(np.hypot(t[m, 1], t[m, 2]).max()))
        if np.any(np.diff(t[:, 0]) < 0):
            reversed_ = True
    return (not reversed_ and r_in <= GOOD_FIELD_MM), r_in, reversed_
MAP_TIMEOUT = 600           # с. Предел на обработку ОДНОЙ карты в дочернем процессе


def scale_scan_range(ec, kind, amp_turns):
    """Расширить диапазон грубого скана под возбуждение карты. Возврат прежнего значения.

    `element_compare` рассчитан на карты при 100 ампер-витках (`SCAN_RANGE` для квадруполя ±3).
    Карты серии сняты до 15000 ампер-витков, где коэффициент в сотню раз больше: без расширения
    скан упирается в границу, Нелдер-Мид не сходится и впустую перебирает 12000 вычислений."""
    key, span = ec.SCAN_RANGE[kind]
    factor = max(1.0, abs(float(amp_turns)) / 100.0)
    ec.SCAN_RANGE[kind] = (key, span * factor)
    return key, span


def probe_ring(radius, n=12):
    """Пробы кика = стартовые точки пучка диффура: центр + n по кругу. Радиусы обязаны совпадать,
    иначе траектории сравниваются в разных точках поля."""
    import math
    return [(0.0, 0.0)] + [(radius * math.cos(2 * math.pi * k / n),
                            radius * math.sin(2 * math.pi * k / n)) for k in range(n)]


def run_diffur_safe(ec, dev, cur, yaw=0.0, timeout=DIFFUR_TIMEOUT, radius=1e-3):
    """`element_compare.run_diffur` с ЖЁСТКИМ таймаутом: при зависании ядро убивается.

    Штатный `run_diffur` ждёт ядро бесконечно. Один зависший процесс занимал 5.7 ГБ и вешал машину,
    поэтому здесь — своя обёртка с `timeout`, которая гарантированно завершает потомка."""
    import subprocess
    from _common import no_window
    geom = ec.build_config(dev, cur, yaw, radius=radius)
    try:
        r = subprocess.run([ec.EXE, str(ec.TMP_EXP)], capture_output=True, text=True,
                           cwd=ec.ROOT, timeout=timeout, **no_window())
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"ядро не уложилось в {timeout} с — прогон прерван, процесс убит")
    if r.returncode != 0:
        raise RuntimeError(f"ядро упало: {r.stderr[-400:]}")
    return ec.read_bin(ec.TMP_EXP), geom


def _purge_temp_devices(path):
    """Убрать временные устройства `_satmodel_*`, если прошлый прогон был убит и не прибрал за собой.

    Бэкап снимается ПОСЛЕ этой чистки: иначе хвосты попадают в бэкап и восстановление их же и вернёт."""
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    keep = [ln for ln in lines if "_satmodel_" not in ln]
    if len(keep) != len(lines):
        print(f"  убрано хвостов прошлого прогона: {len(lines) - len(keep)}")
        with open(path, "w", encoding="utf-8", newline="") as f:
            f.writelines(keep)


def measure_element(element, rows, base_dev, kind, probes=None, limit=None):
    """Прогнать диффур на всех картах элемента и подогнать физкоэф ПО ТРАЕКТОРИЯМ.

    Возврат: список записей {winding, amp_turns, role, coeffs, Leff, rms_mm}."""
    import shutil
    import element_compare as ec
    radius = PROBE_R.get("dq" if kind == "kicker" else kind, 1e-3)
    if probes is None:
        probes = probe_ring(radius)      # пробы кика ровно на стартовых точках пучка диффура
    HWCSV = os.path.join("workspace", "tables", "hardware_library.csv")
    BACKUP = HWCSV + ".satmodel_bak"
    recs = []
    mine = [x for x in rows if x["element"] == element]
    if limit:
        mine = mine[:limit]
    _purge_temp_devices(HWCSV)     # хвосты убитого прогона: бэкап поверх них сделал бы грязь вечной
    shutil.copy2(HWCSV, BACKUP)                           # СВОЙ бэкап: guard_tables эту таблицу не хранит
    ec.guard_tables()
    try:
        hw_rows = ec.mt._read(os.path.basename(HWCSV))
        header = list(hw_rows[0].keys())
        added = []
        for r in mine:
            name, row = _temp_device(ec, base_dev, r)
            added.append(row)
            ec.HW[name] = row
        # ВАЖНО: _write пишет обычным csv.writer — строки должны быть СПИСКАМИ значений, не словарями
        ec._write(os.path.basename(HWCSV), header,
                  [[r.get(c, "") for c in header] for r in hw_rows + added])
        for i, r in enumerate(mine, 1):                   # .bin ДО трекинга: ядро не должно искать сырьё
            print(f"  [{i}/{len(mine)}] конвертация {r['map_id']}", flush=True)
            convert_map(r)
        # ── проход 1: диффур на каждой карте. Подгонка НЕ здесь: L_eff — геометрическая константа
        #    элемента, одна на все карты, поэтому коэффициенты фитятся общим проходом ниже.
        bending = kind in ("dipole", "dq", "kicker")
        keff = "dq" if kind == "kicker" else kind        # корректор кикается сектором, как dq
        runs = []
        for r, row in zip(mine, added):
            dev = row["device_type_id"]
            try:
                arg = float(row["nominal_arg"]) / ec.ENERGY_FACTOR
                trajs, geom = run_diffur_safe(ec, dev, arg, radius=radius)
                ok, r_in, rev = screen_trajectories(ec, trajs)
                if not ok:                                # карта непригодна для подгонки по траекториям
                    why = "разворот частиц" if rev else f"вылет за {GOOD_FIELD_MM:.0f} мм"
                    recs.append(dict(map_id=r["map_id"], winding=r["winding"],
                                     amp_turns=r["amp_turns"], role="rejected",
                                     excitation=r["excitation"], coeffs={}, Leff=0.0, rms_mm=0.0,
                                     r_in_mm=r_in, reason=why))
                    print(f"  {r['map_id']:24} ОТБРАКОВКА: {why}, |r| внутри {r_in:.1f} мм")
                    continue
                runs.append(dict(cur=arg, diff_su=ec.diff_local(trajs, geom["xstart"], bending),
                                 geom=geom, theta0=ec._full_exit_angle(trajs[0]) if bending else 0.0,
                                 rec=r, r_in=r_in))
            except Exception as e:
                print(f"  {r['map_id']:24} {r['split']:5} ПРОПУСК: {type(e).__name__}: {e}")
            finally:
                ec.cleanup()
        # ── проход 2: ОДИН общий L_eff на элемент (проектная величина) + физкоэф на каждую карту
        if runs:
            print(f"  подбор общего L_eff по {len(runs)} картам...", flush=True)
            Lshared = ec.fit_shared(base_dev, keff, list(probes), runs)
            print(f"  L_eff(общий) = {Lshared:.4f} м")
            for q in runs:
                r = q["rec"]
                recs.append(dict(map_id=r["map_id"], winding=r["winding"], amp_turns=r["amp_turns"],
                                 role=r["split"], excitation=r["excitation"], coeffs=q["coeffs"],
                                 Leff=q["Leff"], rms_mm=q["rms_mm"], r_in_mm=q["r_in"]))
                print(f"  {r['map_id']:24} {r['split']:5} " +
                      " ".join(f"{k}={v:+.4e}" for k, v in q["coeffs"].items()) +
                      f"  RMS={q['rms_mm']:.3f}мм")
    finally:
        shutil.copy2(BACKUP, HWCSV)                       # восстановление в ЛЮБОМ случае
        os.remove(BACKUP)
        ec.release_tables()
    return recs


# ─────────────────────── подгонка модели насыщения ──────────────────────────

def _material():
    import numpy as np
    import pandas as pd
    bh = pd.read_csv(MATERIAL, sep=";")
    H, B = bh.iloc[:, 0].to_numpy(float), bh.iloc[:, 1].to_numpy(float)
    mu0 = 4e-7 * np.pi

    def Hof(b):
        a = np.abs(np.atleast_1d(b))
        return np.sign(np.atleast_1d(b)) * np.where(a <= B[-1], np.interp(a, B, H),
                                                    H[-1] + (a - B[-1]) / mu0)
    return Hof


def fit_saturation(points, Hof):
    """Однованный случай: подобрать (A, C, L) в NI = A·y + L·H(C·y) по точкам (NI, y).

    Возврат (A, C, L, средняя ошибка %, максимальная %)."""
    import numpy as np
    from scipy.optimize import least_squares
    # Модель НЕЧЁТНА по знаку возбуждения (сменить знак тока = сменить знак кика), поэтому подгоняем
    # по МОДУЛЯМ. Иначе карта с отрицательным током (у корректора такая есть) даёт A0 < 0, а он
    # логарифмируется — least_squares падает на «Initial guess is outside of provided bounds».
    NI = np.abs(np.array([p[0] for p in points], float))
    y = np.abs(np.array([p[1] for p in points], float))
    order = np.argsort(NI)
    NI, y = NI[order], y[order]
    A0 = NI[0] / y[0]

    def res(p):
        A, C, L = np.exp(p)
        return (A * y + L * Hof(C * y) - NI) / NI
    best = None
    for c0 in (0.003, 0.03, 0.3, 3.0, 30.0):
        f = least_squares(res, np.log([A0 * 0.97, c0 / max(y), 0.05]), max_nfev=20000)
        if best is None or (f.fun ** 2).sum() < (best.fun ** 2).sum():
            best = f
    A, C, L = (float(v) for v in np.exp(best.x))
    e = 100.0 * np.abs(best.fun)
    return A, C, L, float(e.mean()), float(e.max())


def check_holdout(fit, Hof, recs, key, wind):
    """Проверка НА КОНТРОЛЬНЫХ картах (role=test), НЕ участвовавших в подгонке: по измеренному
    физкоэф y предсказать NI = A·y + L·H(C·y) и сравнить с фактическими ампер-витками карты.

    Возврат (средняя ошибка %, максимальная %, точек) или (None, None, 0), если контроля нет."""
    import numpy as np
    A, C, L = fit["A"], fit["C"], fit["L"]
    pts = [(float(r["amp_turns"]), r["coeffs"][key]) for r in recs
           if r["role"] == "test" and r["winding"] == wind and key in r["coeffs"]]
    if not pts:
        return None, None, 0
    # модель нечётна по знаку: сравниваем МОДУЛИ возбуждения, иначе карта с отрицательным током
    # (у корректора есть −1808 А·в) даёт отрицательный знаменатель и отрицательную «погрешность»
    NI = np.abs(np.array([p[0] for p in pts], float)); y = np.abs(np.array([p[1] for p in pts], float))
    NI_pred = A * y + L * Hof(C * y)
    e = 100.0 * np.abs(NI_pred - NI) / NI
    return float(e.mean()), float(e.max()), len(pts)


def save_dataset(element, recs, fits):
    os.makedirs(OUTDIR, exist_ok=True)
    path = os.path.join(OUTDIR, f"{element}.json")
    json.dump(dict(element=element, records=recs, fits=fits), open(path, "w", encoding="utf-8"),
              indent=1, ensure_ascii=False)
    return path


def load_dataset(element):
    path = os.path.join(OUTDIR, f"{element}.json")
    return json.load(open(path, encoding="utf-8")) if os.path.exists(path) else None


# ─────────────────────────────── команды ────────────────────────────────

def cmd_element(args):
    rows = read_registry()
    el = args.element
    if el not in {r["element"] for r in rows}:
        raise SystemExit(f"элемента {el} нет в реестре; см. --list")
    import element_compare as ec
    base = args.device or next((d for d in ec.HW if d.startswith(el)), None)
    if not base:
        raise SystemExit(f"не найдено базовое устройство для {el} — укажите --device")
    kind = args.kind or KIND.get(el) or "quad"
    print(f"элемент {el}: базовое устройство {base}, тип {kind}")
    recs = measure_element(el, rows, base, kind, limit=args.limit)
    if not recs:
        raise SystemExit("ни одной карты не удалось обработать")
    Hof = _material()
    good = [r for r in recs if r["coeffs"]]
    bad = [r for r in recs if r["role"] == "rejected"]
    if bad:
        print(f"  отбраковано карт: {len(bad)} из {len(recs)} "
              f"(поле слишком сильное — киковая модель неприменима)")
    if not good:
        raise SystemExit("все карты элемента отбракованы — подгонять нечего")
    fits = fit_windings(good, Hof, el)
    print("сохранено:", save_dataset(el, recs, fits))
    return 0


# Какой физкоэффициент несёт СМЫСЛ элемента. Выбирать «самый большой по модулю» нельзя: величины
# разной размерности (B в Тл, G в Тл/м) несопоставимы, и у корректора побеждал побочный G вместо
# кика B — контроль по полю 32.5% вместо 0.5%.
KEY_BY_KIND = {"quad": "G", "sext": "S", "oct": "O", "dq": "G", "kicker": "B", "dipole": "B"}


def target_coeff(el, good):
    """Коэффициент, по которому подгоняется насыщение: по НАЗНАЧЕНИЮ элемента, а не по величине."""
    want = KEY_BY_KIND.get(KIND.get(el, ""), "")
    if want and all(want in r["coeffs"] for r in good):
        return want
    return sorted(good[0]["coeffs"], key=lambda k: -abs(good[0]["coeffs"][k]))[0]


def fit_windings(good, Hof, el=""):
    """Подгонка насыщения по обмоткам: точки train → (A, C, L), контроль на картах test.

    ВАЖНО про «невязку подгонки»: на ТРЁХ точках она тождественно нулевая (три параметра проходят
    через три точки) и ничего не доказывает. На ЧЕТЫРЁХ и более она становится содержательной —
    большая невязка означает, что сама форма модели не описывает этот магнит."""
    fits = {}
    key = target_coeff(el, good)
    for wind in sorted({r["winding"] for r in good if r["role"] == "train"}):
        pts = [(float(r["amp_turns"]), r["coeffs"][key]) for r in good
               if r["role"] == "train" and r["winding"] == wind]
        if len(pts) < 3:
            print(f"  обмотка {wind}: годных точек {len(pts)} из трёх — подгонка невозможна, "
                  f"нужна карта на меньшем возбуждении")
            continue
        A, C, L, em, ex = fit_saturation(pts, Hof)
        fits[wind] = dict(coeff=key, A=A, C=C, L=L, err_mean=em, err_max=ex, n=len(pts))
        tm, tx, tn = check_holdout(fits[wind], Hof, good, key, wind)
        note = "" if len(pts) > 3 else "  (на 3 точках невязка всегда 0)"
        if tn:
            fits[wind].update(test_err_mean=tm, test_err_max=tx, test_n=tn)
            print(f"  обмотка {wind}: A={A:.5g} C={C:.5g} L={1000*L:.1f}мм  "
                  f"подгонка {em:.2f}% / {ex:.2f}% (n={len(pts)}){note}  "
                  f"КОНТРОЛЬ {tm:.2f}% / {tx:.2f}% (n={tn})")
        else:
            print(f"  обмотка {wind}: A={A:.5g} C={C:.5g} L={1000*L:.1f}мм  "
                  f"подгонка {em:.2f}% / {ex:.2f}% (n={len(pts)}){note}  КОНТРОЛЬ: карт нет")
    return fits


def cmd_refit(args):
    """Пересчитать подгонку по КЭШИРОВАННЫМ датасетам, взяв роли train/test из реестра заново.

    Физику (диффур + подгонка траекторий) не трогаем: коэффициент карты от разбиения не зависит,
    поэтому смена train/test — это только повторное решение (A, C, L). Секунды вместо получаса."""
    rows = read_registry()
    roles = {r["map_id"]: r["split"] for r in rows}
    Hof = _material()
    for el in sorted({r["element"] for r in rows}):
        d = load_dataset(el)
        if not d:
            continue
        print(f"\n=== {el} ===")
        for rec in d["records"]:
            if rec["role"] != "rejected" and rec["map_id"] in roles:
                rec["role"] = roles[rec["map_id"]]
        good = [r for r in d["records"] if r["coeffs"]]
        if not good:
            print("  нет годных карт"); continue
        print("сохранено:", save_dataset(el, d["records"], fit_windings(good, Hof, el)))
    return 0


def cmd_all(args):
    rows = read_registry()
    for el in sorted({r["element"] for r in rows}):
        print(f"\n=== {el} ===")
        args.element = el
        try:
            cmd_element(args)
        except SystemExit as e:
            print(f"  пропуск: {e}")
    return 0


def _graph_layout(nodes, edges):
    """Координаты узлов: узлы-развилки (степень ≥ 2) в ряд по горизонтали, свободные концы (степень 1,
    полюсные наконечники у зазора) веером от своей развилки наружу."""
    deg = {n: 0 for n in nodes}
    for a, b in edges.values():
        deg[a] += 1; deg[b] += 1
    junc = [n for n in nodes if deg[n] > 1] or list(nodes)[:1]
    leaf = [n for n in nodes if n not in junc]
    pos = {}
    for i, n in enumerate(junc):
        pos[n] = (2.2 * i, 0.0)
    if len(junc) == 1 and len(leaf) == 1:                 # простейший случай: одно ребро
        pos[leaf[0]] = (2.2, 0.0)
        return pos
    own = {j: [] for j in junc}                            # каждый свободный конец — к своей развилке
    for b, (a, c) in edges.items():
        for x, y in ((a, c), (c, a)):
            if x in own and y in leaf:
                own[x].append(y)
    for k, (j, kids) in enumerate(own.items()):
        outward = -1.0 if k == 0 else 1.0                  # первая развилка ветвится влево, прочие вправо
        jx, jy = pos[j]
        for i, n in enumerate(kids):
            off = i - (len(kids) - 1) / 2.0
            pos[n] = (jx + outward * 1.9, jy - 1.25 * off)   # первая обмотка списка — сверху
    return pos


def _branch_figure(el, g, path):
    """Граф магнитной цепи: ВЕРШИНЫ — узлы (развилки железа и полюсные наконечники), РЁБРА — ветви
    магнитопровода, катушки — штрихами прямо на своём ребре.

    Топология берётся из БД (`node_from`/`node_to` в `saturation_branches.csv`). Свободные концы — это
    наконечники полюсов у зазора; замыкание потока через зазор на рисунке не показано."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib.patches import Circle
    branches = g.get("branches", [])
    edges = {b: e for b, e in (g.get("edge") or {}).items() if all(e)}
    if not branches or len(edges) != len(branches):
        return None
    coils_by_branch = {}
    for w, bs in sorted(g.get("coils", {}).items()):
        for b in bs:
            coils_by_branch.setdefault(b, []).append(w)
    nodes = []
    for a, c in edges.values():
        for n in (a, c):
            if n not in nodes:
                nodes.append(n)
    pos = _graph_layout(nodes, edges)
    xs = [p[0] for p in pos.values()]; ys = [p[1] for p in pos.values()]
    fig, ax = plt.subplots(figsize=(max(4.5, (max(xs) - min(xs)) * 1.7 + 2.4),
                                    max(2.6, (max(ys) - min(ys)) * 1.15 + 1.8)))
    deg = {n: 0 for n in nodes}
    for a, c in edges.values():
        deg[a] += 1; deg[c] += 1
    for b in branches:
        (na, nc) = edges[b]
        (x0, y0), (x1, y1) = pos[na], pos[nc]
        wound = coils_by_branch.get(b)
        is_gap = (g.get("kind") or {}).get(b) == "gap"
        ax.plot([x0, x1], [y0, y1], lw=1.2 if is_gap else 1.8,
                color="#b03030" if is_gap else "#33475b", ls=(0, (5, 4)) if is_gap else "-",
                zorder=1, solid_capstyle="round")
        mx, my = (x0 + x1) / 2.0, (y0 + y1) / 2.0
        dx, dy = x1 - x0, y1 - y0
        ln = max(1e-9, (dx * dx + dy * dy) ** 0.5)
        ux, uy = dx / ln, dy / ln                          # вдоль ребра
        px, py = -uy, ux                                   # поперёк ребра
        if wound:                                          # катушка — три штриха поперёк СВОЕГО ребра
            hx, hy = x0 + dx * 0.36, y0 + dy * 0.36
            for k in (-1, 0, 1):
                cx, cy = hx + ux * 0.15 * k, hy + uy * 0.15 * k
                ax.plot([cx - px * 0.19, cx + px * 0.19], [cy - py * 0.19, cy + py * 0.19],
                        lw=2.4, color="#8a5a12", zorder=3, solid_capstyle="round")
        # одна подпись на ребро, прямо НА нём (иначе у веера подпись уезжает к соседнему ребру)
        lx, ly = x0 + dx * 0.70, y0 + dy * 0.70
        ax.text(lx, ly, f"{b}\n[{'+'.join(wound)}]" if wound else b, ha="center", va="center",
                fontsize=8.5, linespacing=1.25, zorder=4,
                color="#5a3d0c" if wound else "#8894a2",
                bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="none", alpha=0.95))
    for n, (x, y) in pos.items():
        junction = deg[n] > 1
        ax.add_patch(Circle((x, y), 0.10 if junction else 0.075,
                            fc="#33475b" if junction else "white",
                            ec="#33475b", lw=1.4, zorder=5))
    n_wound = len(coils_by_branch)
    ax.set_title(f"{el}: граф магнитной цепи — {len(branches)} ветв"
                 f"{'ь' if len(branches) == 1 else ('и' if len(branches) < 5 else 'ей')} (рёбра), "
                 f"из них с обмотками {n_wound};\nполые кружки — наконечники полюсов у зазора",
                 fontsize=9.5)
    ax.set_xlim(min(xs) - 1.0, max(xs) + 1.0); ax.set_ylim(min(ys) - 0.9, max(ys) + 0.9)
    ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout()
    fig.savefig(path, dpi=160); plt.close(fig)
    return path


def build_report(data, rows, path):
    """Word-отчёт по готовым датасетам: на элемент — граф, формула, матрицы, проверка."""
    from docx import Document
    import latex2omml as lo
    import saturation_report as sr                       # общие помощники оформления
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    sr._tcap.n = sr._fcap.n = 0
    doc = Document()
    doc.add_heading("Модель насыщения по элементам библиотеки СИЛА", level=0)
    doc.add_paragraph(
        "Для каждого элемента: граф ветвей, расчётная формула, подобранные числа и проверка на "
        "картах, которые в подгонке не участвовали. Физкоэффициенты получены сравнением траекторий "
        "диффурной и киковой моделей; насыщение подогнано по сериям карт при разных возбуждениях.")
    tests = {}
    for r in rows:
        if r["split"] == "test":
            tests.setdefault(r["element"], []).append(r["excitation"])
    for el in sorted(data):
        d = data[el]
        g = BRANCH_GRAPH.get(el, {})
        doc.add_heading(el, level=1)
        sr._table(doc, ["Свойство", "Значение"], [
            ("Ветвей магнитопровода", str(len(g.get("branches", [])))),
            ("Из них без обмоток", str(len(g.get("unwound", [])))),
            ("Обмоток", str(len(g.get("coils", {})))),
            ("Карт для подгонки", str(sum(1 for r in d["records"] if r["role"] == "train"))),
            ("Проверочных комбинаций", str(len(tests.get(el, ())))),
        ], caption=f"{el} — состав данных")
        fig_path = _branch_figure(el, g, os.path.join(OUTDIR, f"graph_{el}.png"))
        if fig_path:
            sr._pic(doc, fig_path, cm=15)
            sr._fcap(doc, f"{el} — граф магнитной цепи: рёбра между узлами = ветви магнитопровода "
                          f"(тонкие серые — ветви без обмоток)")
        if g.get("coils"):
            sr._table(doc, ["Обмотка", "Ребро (ветвь) графа", "витков/катушку"],
                      [[w, " + ".join(g["coils"][w]), str(g["turns"].get(w, "—"))] for w in sorted(g["coils"])],
                      caption=f"{el} — какая обмотка на каком ребре графа (матрица T)")
        doc.add_heading("Расчётная формула", level=2)
        sr._formula(doc, r"\Phi = Q\left( T\,\mathbf{I} - \Lambda\, H(\Sigma^{-1}\Phi) \right),"
                         r" \qquad \mathbf{y} = G\,\Phi", lo)
        sr._symtab(doc, [
            ["Φ", "поток ветви", "вектор", "Вб"],
            ["Q", "проводимости ветвей", "матрица", "Вб/А"],
            ["T", "размещение обмоток по рёбрам", "матрица", "—"],
            ["I", "токи обмоток", "вектор", "А"],
            ["Λ", "длина пути по железу", "диаг. матрица", "м"],
            ["Σ", "сечение сердечника", "диаг. матрица", "м²"],
            ["H(·)", "обратная кривая намагничивания", "функция", "А/м"],
            ["y", "физкоэффициенты магнита", "вектор", "Тл/м^(n−1)"],
            ["G", "поток → физкоэффициент", "матрица", "Тл/(м^(n−1)·Вб)"],
        ], caption="обозначения")
        if len(g.get("branches", [])) == 1:
            doc.add_paragraph("Одна ветвь — система сворачивается в одно уравнение:")
            sr._formula(doc, r"N I = A\,y + L\,H\!\left( C\,y \right)", lo)
            sr._symtab(doc, [
                ["NI", "ампер-витки обмотки", "скаляр", "А"],
                ["y", "физкоэффициент (G, S или O)", "скаляр", "Тл/м^(n−1)"],
                ["A", "линейный коэффициент цепи", "скаляр", "А/[y]"],
                ["L", "длина пути по железу", "скаляр", "м"],
                ["C", "y → индукция сердечника", "скаляр", "Тл/[y]"],
                ["H(·)", "обратная кривая намагничивания", "функция", "А/м"],
            ], caption="обозначения")
        doc.add_heading("Подобранные числа", level=2)
        if d["fits"]:
            def _ctrl(f):
                return (f"{f['test_err_mean']:.1f}% / {f['test_err_max']:.1f}% (n={f['test_n']})"
                        if f.get("test_n") else "карт нет")
            sr._table(doc, ["Обмотка", "Коэфф.", "A", "C", "L, мм", "точек", "невязка подгонки",
                            "КОНТРОЛЬ на неучтённых картах"],
                      [[w, f["coeff"], f"{f['A']:.5g}", f"{f['C']:.5g}", f"{1000*f['L']:.1f}",
                        str(f["n"]), f"{f['err_mean']:.2f}% / {f['err_max']:.2f}%", _ctrl(f)]
                       for w, f in sorted(d["fits"].items())],
                      caption=f"{el} — параметры насыщения по обмоткам. Невязка подгонки на ТРЁХ "
                              f"точках тождественно нулевая (три параметра проходят через три точки) "
                              f"и ничего не доказывает; на четырёх и более она содержательна. "
                              f"Контроль — карты, не участвовавшие в подгонке")
        else:
            doc.add_paragraph("Параметры не подобраны: не хватает карт (нужно три на обмотку).")
        doc.add_heading("Измеренные точки", level=2)
        keys = sorted({k for r in d["records"] for k in r["coeffs"]})
        sr._table(doc, ["Возбуждение", "Роль"] + keys + ["L_eff, м", "RMS, мм"],
                  [[r["excitation"], r["role"]] + [f"{r['coeffs'].get(k, 0):+.4e}" for k in keys]
                   + [f"{r['Leff']:.4f}", f"{r['rms_mm']:.3f}"] for r in d["records"]],
                  caption=f"{el} — физкоэффициенты, подогнанные по траекториям", pt=8)
    doc.add_heading("Кривая намагничивания", level=1)
    F = sr._figures()
    sr._pic(doc, F["bh"], cm=16)
    sr._fcap(doc, "Кривая намагничивания материала сердечника")
    import pandas as pd
    bh = pd.read_csv(MATERIAL, sep=";")
    sr._table(doc, ["H, А/м", "B, Тл"], [[f"{r[0]:g}", f"{r[1]:g}"] for r in bh.to_numpy()],
              caption="Точки таблицы материала")
    sr._body_format(doc)
    doc.save(path)
    return path


def cmd_report(args):
    rows = read_registry()
    have = {el: d for el in sorted({r["element"] for r in rows}) if (d := load_dataset(el))}
    if not have:
        raise SystemExit(f"нет готовых датасетов в {OUTDIR} — сначала --element или --all")
    print("отчёт ->", build_report(have, rows, args.out))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description="Графовая модель насыщения: подгонка и Word-отчёт")
    ap.add_argument("--list", action="store_true", help="что есть в реестре карт и чего не хватает")
    ap.add_argument("--element", help="прогнать один элемент")
    ap.add_argument("--all", action="store_true", help="прогнать все элементы реестра")
    ap.add_argument("--report", action="store_true", help="собрать Word-отчёт по готовым датасетам")
    ap.add_argument("--refit", action="store_true",
                    help="пересчитать (A,C,L) по кэшу с новыми ролями train/test из реестра (без физики)")
    ap.add_argument("--limit", type=int, help="обработать только первые N карт элемента")
    ap.add_argument("--device", help="базовое устройство (если имя не выводится из элемента)")
    ap.add_argument("--kind", help="тип элемента (quad/sext/oct/dq/kicker)")
    ap.add_argument("--out", default=os.path.join(OUTDIR, "saturation_model_report.docx"))
    a = ap.parse_args(argv)
    if a.list:
        return cmd_list(a)
    if a.element:
        return cmd_element(a)
    if a.all:
        return cmd_all(a)
    if a.refit:
        return cmd_refit(a)
    if a.report:
        return cmd_report(a)
    ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
