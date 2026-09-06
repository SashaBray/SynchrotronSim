# -*- coding: utf-8 -*-
"""Датасеты «токи обмоток → киковые моменты» по всем элементам библиотеки — под разработку теории.

Собирает в отдельную папку по одному CSV на элемент: строка = карта поля, столбцы = токи обмоток и
киковые моменты, извлечённые из этой карты. Комбинационные карты (несколько обмоток под током
одновременно) входят наравне с одиночными — на них теория и проверяется.

Моменты СО ЗНАКОМ и в двух видах: интеграл гармоники вдоль элемента (физические единицы при 6 ГэВ,
длина не участвует) и физкоэффициент в шкале `mat_*` решётки. Рядом — коэффициенты из подгонки по
траекториям (`saturation_model`), где они есть: это независимый способ, и расхождение с
разложением карт видно прямо в таблице.

    py -3.14 build_kick_datasets.py            # собрать всё
    py -3.14 build_kick_datasets.py --element DQ
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

OUTDIR = os.path.join("workspace", "datasets", "currents_to_kicks")
CACHE = os.path.join(OUTDIR, "_harmonics_cache.json")

# Порядок мультиполя → имя нормального и скью момента (обозначения `element_compare.MULTIPOLES`).
ORDER_NAMES = {1: ("B", "Bv"), 2: ("G", "Gv"), 3: ("S", "Ss"), 4: ("O", "Os")}
ORDERS = (1, 2, 3, 4)
# Радиус пробы для каждого порядка. Поле мультиполя ∝ r^(n−1), а шум карты аддитивен, поэтому
# высокий момент у самой оси лежит под шумовым полом: октуполь на 1 мм мерился бы в шуме.
PROBE_R = {1: 1e-3, 2: 5e-3, 3: 8e-3, 4: 14e-3}
# Радиусы для проверки разрешимости — свои у каждого порядка, те, где момент вообще измерим.
RES_R = {1: (1e-3, 3e-3, 5e-3, 8e-3), 2: (1e-3, 3e-3, 5e-3, 8e-3),
         3: (5e-3, 8e-3, 11e-3, 14e-3), 4: (8e-3, 11e-3, 14e-3)}
RESOLVED = 25.0        # % разброса по радиусу, выше которого момент считается шумом
# Ведущий момент элемента — тот, ради которого магнит и сделан. По нему определяется длина.
# Выбирать его сравнением величин НЕЛЬЗЯ: у моментов разные единицы (октуполь в Тл/м² численно
# больше квадруполя в Тл на четыре порядка), и «самым сильным» всегда оказывался бы октуполь.
KIND_ORDER = {"quad": 2, "sext": 3, "oct": 4, "dq": 2, "kicker": 1, "dipole": 1}
# Знак при переводе гармоники в шкалу `mat_*`: у ДИПОЛЯ он противоположен, у остальных порядков
# совпадает. Это не догадка: сверено с подгонкой по траекториям на диполе-квадруполе и корректоре —
# после снятия длинового множителя отношение mat/fit равно −1.0000 у B и +0.997 у G. У дипольного
# кика знак задаётся зарядом частицы и стороной изгиба, у производных поля — нет.
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0, 4: 1.0}


def load_full(path):
    """Карта из .bin, приведённая к ПОЛНОЙ: библиотека хранит половинки (X≥0)."""
    from field_quality import load
    import _qp_scale as qs
    h, B = load(path)
    if h["mnx"] >= -1e-9:
        h, B = qs.mirror_full(h, B)
    return h, B


def harmonics(path, radii, rot_z=0.0):
    """{радиус: {порядок: (Re, Im)}} — КОМПЛЕКСНЫЕ интегралы ∫b_n·dl, со знаком.

    Гармоника считается как (By + i·Bz) на окружности радиуса r, разложенная в ряд Фурье и
    проинтегрированная вдоль элемента. Мнимая часть — нормальный мультиполь, вещественная —
    скошенный (повёрнутый на пол-периода). Модуль брать нельзя: он их смешивает.

    Поворот установленного магнита (`rot_z_deg` в hardware_library) применяется АНАЛИТИЧЕСКИ:
    поворот на 180° вокруг поперечной оси переводит C_n → (−1)ⁿ·conj(C_n). Проверено численно
    против поворота самой сетки. У всей библиотеки, кроме диполя-квадруполя, поворота нет."""
    from field_quality import _circle
    h, B = load_full(path)
    xs = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    out = {}
    for r in radii:
        acc = {n: [] for n in ORDERS}
        for ix in range(h["nx"]):
            by, bz = _circle(h, B, ix, r, nth=64)
            sp = np.fft.fft(by + 1j * bz) / 64
            for n in ORDERS:
                acc[n].append(sp[0] if n == 1 else sp[-(n - 1)])
        out[r] = {}
        for n in ORDERS:
            v = np.array(acc[n])
            c = complex(np.trapezoid(v.real, xs), np.trapezoid(v.imag, xs)) / r ** (n - 1)
            if abs(rot_z - 180.0) < 1e-9:
                c = ((-1) ** n) * c.conjugate()
            out[r][n] = (c.real, c.imag)
    return out


def cached_harmonics(rows, rot_by_el):
    """Гармоники всех карт с кэшем на диске: пересчёт занимает минуты, а нужен он редко."""
    cache = {}
    if os.path.exists(CACHE):
        try:
            cache = json.load(open(CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    dirty = False
    radii = sorted({r for rs in RES_R.values() for r in rs} | set(PROBE_R.values()))
    for i, r in enumerate(rows, 1):
        p = os.path.join("workspace", "fields_bin", r["bin_path"])
        if not os.path.exists(p):
            continue
        have = cache.setdefault(r["map_id"], {})
        need = [rr for rr in radii if f"{rr:.6f}" not in have]
        if need:
            print(f"  [{i}/{len(rows)}] {r['map_id']}", flush=True)
            s = harmonics(p, need, rot_by_el.get(r["element"], 0.0))
            have.update({f"{rr:.6f}": {str(n): list(v) for n, v in d.items()}
                         for rr, d in s.items()})
            dirty = True
    if dirty:
        os.makedirs(os.path.dirname(CACHE) or ".", exist_ok=True)
        json.dump(cache, open(CACHE, "w", encoding="utf-8"))
    return {k: {float(rr): {int(n): tuple(v) for n, v in d.items()} for rr, d in v0.items()}
            for k, v0 in cache.items()}


def parse_excitation(text):
    """`I1=1891;I2=-800` → {'I1': 1891.0, 'I2': -800.0} (ампер-витки, со знаком)."""
    out = {}
    for part in (text or "").split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            try:
                out[k.strip()] = float(v)
            except ValueError:
                pass
    return out


def resolution(series, order):
    """Разброс момента по радиусу пробы, % от медианы — по верхней половине возбуждений.

    Настоящий момент, поделённый на r^(n−1), от радиуса не зависит; шум зависит."""
    top = sorted(series, key=lambda t: t[0])[len(series) // 2:]
    sc = []
    for _, hh in top:
        v = np.array([abs(hh[r][order][1]) for r in RES_R[order] if r in hh])
        m = np.median(v)
        if m > 0:
            sc.append(100 * (v.max() - v.min()) / m)
    return float(np.mean(sc)) if sc else float("nan")


def winding_turns():
    """{элемент: {обмотка: витков на катушку}} из таблицы ветвей насыщения.

    Не из `hardware_library`: там витки записаны у ИСПОЛНЕНИЙ (`SP_I1_1A_L166`), а карты привязаны
    к базовым устройствам (`SP_I1_1A`), у которых поля витков пустые."""
    out = collections.defaultdict(dict)
    with open(os.path.join("workspace", "tables", "saturation_branch_members.csv"),
              encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter=";"):
            if r.get("turns_per_coil"):
                out[r["element"]][r["winding"]] = int(r["turns_per_coil"])
    return out


def nominal_current(ec, element, winding):
    """Номинальный ток обмотки — берётся у исполнения элемента, где он заполнен."""
    for d, hw in ec.HW.items():
        if not d.startswith(element) or not (hw.get("nominal_current_a") or "").strip():
            continue
        tok = d[len(element) + 1:].split("_")[0]
        if tok == winding or tok.endswith("A"):
            try:
                return float(hw["nominal_current_a"].split()[0])
            except (ValueError, IndexError):
                pass
    return None


def fitted_coeffs(element):
    """{map_id: {момент: коэффициент}} из подгонки по траекториям, если датасет есть."""
    p = os.path.join(sm.OUTDIR, f"{element}.json")
    if not os.path.exists(p):
        return {}, None
    d = json.load(open(p, encoding="utf-8"))
    out = {r["map_id"]: (r.get("coeffs") or {}) for r in d["records"]}
    leff = next((r.get("Leff") for r in d["records"] if r.get("Leff")), None)
    return out, leff


def native_length(element, rows):
    """Эффективная длина ИЗМЕРЕННОГО магнита по форме поля: ∫профиля / максимум. Одна на элемент.

    Нужна, чтобы перевести интегралы в шкалу `mat_*`. Паспортные L_eff для этого не годятся: они
    получены подгонкой по траекториям, а подгонка определяет только ПРОИЗВЕДЕНИЕ «коэффициент ×
    длина» и не определяет, как разделить его на сомножители, — оттого у обмоток одного и того же
    магнита паспортные длины разные (у диполя-квадруполя 0.4509, 0.4549 и 0.4608 м). Длина, взятая
    из формы поля, одна на элемент и совпадает с проектной: у длинного квадруполя 0.374 против
    проектных 0.3745, у диполя-квадруполя 0.403 против 0.400.

    Берётся по карте среднего возбуждения: у слабых всё тонет в шуме, у самых сильных профиль
    заострён насыщением."""
    from field_quality import _circle
    singles = [r for r in rows if ";" not in r["excitation"] and r["amp_turns"]]
    if not singles:
        return None, 2
    lead = KIND_ORDER.get(sm.KIND.get(element, ""), 2)
    m = sorted(singles, key=lambda x: abs(float(x["amp_turns"])))[len(singles) // 2]
    h, B = load_full(os.path.join("workspace", "fields_bin", m["bin_path"]))
    xs = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    v = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, PROBE_R[lead], nth=64)
        sp = np.fft.fft(by + 1j * bz) / 64
        v.append(abs((sp[0] if lead == 1 else sp[-(lead - 1)]).imag))
    v = np.array(v)
    return float(np.trapezoid(v, xs) / v.max()), lead


def build_element(element, rows, harm, ec, outdir):
    """CSV датасета элемента + JSON описания. Возврат: строка сводки."""
    wnames = []
    for r in rows:
        for w in parse_excitation(r["excitation"]):
            if w not in wnames:
                wnames.append(w)
    wnames.sort()

    # Устройство обмотки берётся из реестра карт (столбец device_type_id) — это надёжная связь
    # карта→устройство; разбирать имя обратно нельзя, у секступоля и DQ схемы имён разные.
    devs, turns, nominal, Ldev = {}, {}, {}, {}
    wt = winding_turns().get(element, {})
    for w in wnames:
        dev = next((r["device_type_id"] for r in rows
                    if list(parse_excitation(r["excitation"])) == [w]), None)
        devs[w] = dev
        turns[w] = wt.get(w) or (ec.dev_tpc(dev) if dev else None) or 1
        nominal[w] = nominal_current(ec, element, w)
        Le = (ec.HW.get(dev, {}) or {}).get("eff_length") if dev else None
        Ldev[w] = float(Le) if Le else None

    Lm, lead = native_length(element, rows)
    Lconv = Lm                               # длина для перевода интеграла в шкалу mat_*
    ef = ec.ENERGY_FACTOR

    fit, fit_leff = fitted_coeffs(element)

    head = (["map_id", "split", "excitation"]
            + [f"AT_{w}" for w in wnames] + [f"I_{w}" for w in wnames]
            + [f"int_{nm}" for n in ORDERS for nm in ORDER_NAMES[n]]
            + [f"mat_{nm}" for n in ORDERS for nm in ORDER_NAMES[n]]
            + [f"fit_{nm}" for n in ORDERS for nm in ORDER_NAMES[n] if nm in ("B", "G", "S", "O")])
    out = []
    for r in sorted(rows, key=lambda x: x["map_id"]):
        if r["map_id"] not in harm:
            continue
        exc = parse_excitation(r["excitation"])
        hh = harm[r["map_id"]]
        row = [r["map_id"], r["split"], r["excitation"]]
        row += [f"{exc.get(w, 0.0):.6g}" for w in wnames]
        row += [f"{exc.get(w, 0.0) / turns[w]:.6g}" for w in wnames]
        ints, mats = [], []
        for n in ORDERS:
            re_, im_ = hh[PROBE_R[n]][n]
            k = MAT_SIGN[n] * math.factorial(n - 1) / (Lconv * ef) if Lconv else float("nan")
            for val in (im_, re_):                       # нормальный, затем скошенный
                ints.append(f"{val:.6e}")
                mats.append(f"{val * k:.6e}")
        row += ints + mats
        # Подгонка приведена к ТОЙ ЖЕ длине, что и mat_*: сама она делит интеграл на свою L_eff,
        # и без приведения столбцы отличались бы ровно на отношение длин, а не по существу.
        fc = fit.get(r["map_id"], {})
        kL = (fit_leff / Lconv) if (fit_leff and Lconv) else 1.0
        row += [(f"{fc[nm] * kL:.6e}" if nm in fc else "") for n in ORDERS
                for nm in ORDER_NAMES[n] if nm in ("B", "G", "S", "O")]
        out.append(row)

    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, f"{element}.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(head)
        w.writerows(out)

    res = {}
    for w in wnames:
        ser = [(abs(parse_excitation(r["excitation"])[w]), harm[r["map_id"]])
               for r in rows
               if list(parse_excitation(r["excitation"])) == [w] and r["map_id"] in harm]
        if ser:
            res[w] = {ORDER_NAMES[n][0]: round(resolution(ser, n), 1) for n in ORDERS}

    meta = dict(
        element=element, kind=sm.KIND.get(element, ""),
        maps=len(out),
        singles=sum(1 for r in out if ";" not in r[2]),
        combos=sum(1 for r in out if ";" in r[2]),
        windings=wnames, devices=devs, turns_per_coil=turns, nominal_current_a=nominal,
        L_eff_from_field_m=Lm, L_used_for_mat_m=Lconv, lead_order=lead,
        L_eff_device_m=Ldev, L_eff_trajectory_fit_m=fit_leff,
        energy_factor=ef, probe_radius_m={str(n): PROBE_R[n] for n in ORDERS},
        rot_z_deg=float((ec.HW.get(devs[wnames[0]], {}) or {}).get("rot_z_deg") or 0.0),
        resolution_pct=res, resolved_threshold_pct=RESOLVED,
        columns=dict(
            AT="возбуждение обмотки, ампер-витки (со знаком)",
            I="ток питания обмотки, А = AT / витков на катушку",
            int="интеграл гармоники вдоль элемента, физические единицы при 6 ГэВ, со знаком; "
                "нормальная компонента (B, G, S, O) и скошенная (Bv, Gv, Ss, Os)",
            mat="то же в шкале mat_* решётки: int × знак × (n−1)! / (L_used × energy_factor); "
                "знак у диполя противоположный, у остальных порядков совпадает",
            fit="коэффициент из подгонки по траекториям (saturation_model), приведённый к той же "
                "длине L_used — сравним с mat_* напрямую"),
        mat_sign={str(n): MAT_SIGN[n] for n in ORDERS})
    with open(os.path.join(outdir, f"{element}.meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)
    return meta


README = """# Датасеты «токи обмоток → киковые моменты»

Собрано скриптом `build_kick_datasets.py` из реестра карт поля (`workspace/tables/field_maps.csv`).
На каждый элемент — `<элемент>.csv` (строка = карта поля) и `<элемент>.meta.json` (паспорт: обмотки,
витки, номиналы, длины, радиусы пробы, разрешимость моментов). Сводка по всем — `index.csv`.

## Столбцы

* `map_id`, `split` (train/test), `excitation` — исходная карта.
* `AT_<обмотка>` — возбуждение, ампер-витки со знаком. Ноль означает обесточенную обмотку.
* `I_<обмотка>` — ток питания, А: ампер-витки, делённые на число витков в катушке.
* `int_<момент>` — интеграл гармоники вдоль элемента, физические единицы при 6 ГэВ, СО ЗНАКОМ.
  Длина в него не входит: это то, что чувствует пучок.
* `mat_<момент>` — то же в шкале `mat_*` решётки: `int × (n−1)! / (L_used × energy_factor)`.
* `fit_<момент>` — коэффициент из подгонки по траекториям (`saturation_model`), где он есть.
  Независимый способ: на сильных картах сходится с разложением в пределах процента.

## Моменты

`B`/`G`/`S`/`O` — нормальные диполь, квадруполь, секступоль, октуполь; `Bv`/`Gv`/`Ss`/`Os` — их
скошенные (повёрнутые на пол-периода) пары. Гармоника считается как `(By + i·Bz)` на окружности
радиуса пробы: мнимая часть даёт нормальный момент, вещественная — скошенный.

## На что смотреть при разработке теории

* Комбинационные карты (`;` в `excitation`) — главный материал: одиночные карты линейную модель
  не опровергают, нелинейность видна только на комбинациях.
* Радиус пробы у каждого порядка свой: поле мультиполя растёт как `r^(n−1)`, а шум карты от радиуса
  не зависит, поэтому высокий момент у самой оси лежит под шумовым полом. Разброс момента по
  радиусу приведён в `meta.json` — момент с разбросом выше 25% в картах не разрешается, и его
  величину по ним назвать нельзя.
* Поворот установленного магнита учтён: у диполя-квадруполя `rot_z_deg = 180`, и он применён
  аналитически (`C_n → (−1)ⁿ·conj(C_n)`). Без него знак квадруполя был бы противоположен проекту.
"""


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", help="только один элемент")
    ap.add_argument("--out", default=OUTDIR)
    a = ap.parse_args(argv)

    import element_compare as ec
    rows = [r for r in sm.read_registry()
            if os.path.exists(os.path.join("workspace", "fields_bin", r["bin_path"]))]
    els = sorted({r["element"] for r in rows})
    if a.element:
        els = [a.element]
        rows = [r for r in rows if r["element"] == a.element]

    rot = {}
    for el in els:
        dev = next((d for d in ec.HW if d.startswith(el)), None)
        rot[el] = float((ec.HW.get(dev, {}) or {}).get("rot_z_deg") or 0.0)
    bad = {e: v for e, v in rot.items() if v not in (0.0, 180.0)}
    if bad:
        raise SystemExit(f"поворот, не сводимый к аналитике: {bad}")

    print(f"карт: {len(rows)}, элементов: {len(els)}")
    harm = cached_harmonics(rows, rot)

    os.makedirs(a.out, exist_ok=True)
    metas = []
    for el in els:
        meta = build_element(el, [r for r in rows if r["element"] == el], harm, ec, a.out)
        metas.append(meta)
        print(f"  {el:10} карт {meta['maps']:3d} (одиночных {meta['singles']:3d}, "
              f"комбинаций {meta['combos']:3d})  обмоток {len(meta['windings'])}  "
              f"L={meta['L_used_for_mat_m']:.4f} м  rot_z={meta['rot_z_deg']:.0f}")

    with open(os.path.join(a.out, "index.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["element", "kind", "maps", "singles", "combos", "windings",
                    "L_used_for_mat_m", "rot_z_deg"])
        for m in metas:
            w.writerow([m["element"], m["kind"], m["maps"], m["singles"], m["combos"],
                        ",".join(m["windings"]), f"{m['L_used_for_mat_m']:.6g}",
                        f"{m['rot_z_deg']:.0f}"])
    with open(os.path.join(a.out, "README.md"), "w", encoding="utf-8") as f:
        f.write(README)
    print(f"\nдатасеты -> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
