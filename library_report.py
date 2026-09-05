# -*- coding: utf-8 -*-
"""Word-отчёт по библиотеке элементов: паспорт, кривые «кик от тока», проектные рабочие точки.

На каждый элемент — таблица параметров из БД (витки, число катушек, номинальный ток, эффективная
длина) и график зависимости физкоэффициента от тока обмотки, построенный по СНЯТЫМ КАРТАМ. На
графике отмечены проектные уровни решётки и ток, при котором каждый из них достигается.

Коэффициенты берутся из датасетов `saturation_model`, то есть получены сопоставлением ТРАЕКТОРИЙ
(диффурное ядро по карте поля), а не гармоническим разложением.

У корректора отдельный раздел: он многообмоточный и его обмотки создают не одну величину, а смесь
гармоник, поэтому на каждую обмотку строится свой график со всеми гармониками сразу.

    py -3.14 library_report.py
"""
from __future__ import annotations

import argparse
import collections
import csv
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

OUT = os.path.join("workspace", "results", "library_report.docx")
PNG = os.path.join("workspace", "results", "library_report_figs")

ELEMENTS = ["QP_sh", "QP_L", "SP", "OP", "DQ"]
KEY = {"QP_sh": "G", "QP_L": "G", "SP": "S", "OP": "O", "DQ": "G", "corrector": "B"}
UNITS = {"B": "Тл", "G": "Тл/м", "S": "Тл/м²", "O": "Тл/м³"}
COLUMN = {"B": "mat_B", "G": "mat_G", "S": "mat_S", "O": "mat_O"}
NAMES = {"QP_sh": "Квадруполь короткий", "QP_L": "Квадруполь длинный",
         "SP": "Секступоль", "OP": "Октуполь", "DQ": "Диполь-квадруполь",
         "corrector": "Корректор"}
ORDER_NAME = {1: "диполь B", 2: "квадруполь G", 3: "секступоль S", 4: "октуполь O"}
KICK_KEY = {1: "B", 2: "G", 3: "S", 4: "O"}
# Единицы ИНТЕГРАЛА коэффициента по длине — это и есть проектный кик: ∫G·dl, ∫S·dl, ∫O·dl.
INT_UNITS = {"B": "Тл·м", "G": "Тл", "S": "Тл/м", "O": "Тл/м²"}


# ───────────────────────────── данные ─────────────────────────────

def hardware():
    import matrix_tracker as mt
    return mt._read("hardware_library.csv")


def devices(element):
    """Семейства устройств элемента с паспортом: витки, катушки, номинальный ток, длина."""
    out = []
    for r in hardware():
        d = r["device_type_id"]
        if not d.startswith(element) or not (r.get("nominal_current_a") or "").strip():
            continue
        out.append(dict(dev=d, turns=r.get("turns_per_coil", ""), coils=r.get("num_coils", ""),
                        total=r.get("total_turns", ""), inom=r.get("nominal_current_a", ""),
                        leff=r.get("eff_length", "")))
    return sorted(out, key=lambda x: x["dev"])


def nominal_by_winding(element):
    """{обмотка: номинальный ток, А} — из имён устройств вида `SP_I1_1A_L166`.

    У однообмоточных элементов в имени обмотки нет, и номинал возвращается под ключом None."""
    out = {}
    for d in devices(element):
        tok = d["dev"][len(element) + 1:].split("_")[0]
        try:
            out[tok if not tok.endswith("A") else None] = float(d["inom"].split()[0])
        except (ValueError, IndexError):
            pass
    return out


def design_levels(element, key):
    """Проектные значения из решётки 1: {(устройство, L_eff, значение): количество}."""
    col = COLUMN[key]
    rows = list(csv.DictReader(open(os.path.join("workspace", "tables", "lattice_configs.csv"),
                                    encoding="utf-8-sig"), delimiter=";"))
    c = collections.Counter()
    for r in rows:
        if str(r.get("lattice_id")) != "1" or not r.get("device_type_id", "").startswith(element):
            continue
        try:
            v, L = abs(float(r.get(col) or 0)), float(r.get("eff_length") or 0)
        except ValueError:
            continue
        if v > 1e-9:
            c[(r["device_type_id"], round(L, 4), round(v, 5))] += 1
    tot = sum(c.values()) or 1
    return {k: n for k, n in c.items() if n >= max(3, 0.02 * tot)}


def curve(element, key, winding=None):
    """(ампер-витки, коэффициент) по одиночным картам; плюс общий L_eff подгонки."""
    p = os.path.join(sm.OUTDIR, f"{element}.json")
    if not os.path.exists(p):
        return [], None
    d = json.load(open(p, encoding="utf-8"))
    out, leff = [], None
    for r in d["records"]:
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not r["coeffs"] or key not in r["coeffs"]:
            continue
        w = r.get("winding") or ""
        if "+" in w or (winding and w != winding):
            continue
        out.append((ni, abs(r["coeffs"][key])))
        leff = r.get("Leff") or leff
    return sorted(out), leff


def first_crossing(x, y, level):
    """Ток, при котором коэффициент ВПЕРВЫЕ достигает уровня. Кривая может быть немонотонна на
    плато, поэтому простая интерполяция непригодна — берём первое пересечение."""
    for i in range(1, len(y)):
        if (y[i - 1] - level) * (y[i] - level) <= 0 and y[i] != y[i - 1]:
            t = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(x[i - 1] + t * (x[i] - x[i - 1]))
    return None


PROBE_SCAN = (1e-3, 3e-3, 5e-3, 8e-3, 14e-3)
SPEC_CACHE = os.path.join("workspace", "results", "library_report_spectra.json")


def spectrum(path, orders, radii, nth=64):
    """Гармоники на НЕСКОЛЬКИХ радиусах пробы за одну загрузку карты: {радиус: {порядок: |кик|}}.

    Несколько радиусов нужны не для полноты. Поле мультиполя ∝ r^(n−1), а шум карты аддитивен,
    поэтому у слабой высокой гармоники сигнал у оси лежит ПОД шумовым полом. Настоящую гармонику
    от шума отличает то, что после деления на r^(n−1) она НЕ ЗАВИСИТ от радиуса пробы."""
    from field_quality import load, _circle
    h, B = load(path)
    xs = np.linspace(h["mnx"], h["mxx"], h["nx"])
    out = {}
    for r in radii:
        acc = {n: [] for n in orders}
        for ix in range(h["nx"]):
            by, bz = _circle(h, B, ix, r, nth=nth)
            sp = np.fft.fft(by + 1j * bz) / nth
            for n in orders:
                acc[n].append(sp[0] if n == 1 else sp[-(n - 1)])
        out[r] = {}
        for n in orders:
            v = np.array(acc[n])
            c = complex(2 * np.trapezoid(v.real, xs), 2 * np.trapezoid(v.imag, xs))
            out[r][n] = abs(c) / r ** (n - 1)
    return out


def corrector_spectra(orders=(1, 2, 3, 4), radii=PROBE_SCAN):
    """{обмотка: [(ампер-витки, {радиус: {порядок: значение}})]} по одиночным картам, с кэшем."""
    cache = {}
    if os.path.exists(SPEC_CACHE):
        try:
            cache = json.load(open(SPEC_CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    dirty = False
    data = collections.defaultdict(list)
    for r in sm.read_registry():
        if r["element"] != "corrector" or ";" in r["excitation"] or "+" in r["excitation"]:
            continue
        p = os.path.join("workspace", "fields_bin", r["bin_path"])
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not os.path.exists(p):
            continue
        key = r["map_id"]
        if key not in cache:
            s = spectrum(p, orders, radii)
            cache[key] = {f"{rr:.6f}": {str(n): v for n, v in d.items()} for rr, d in s.items()}
            dirty = True
        s = {float(rr): {int(n): v for n, v in d.items()} for rr, d in cache[key].items()}
        data[r["winding"]].append((ni, s))
    if dirty:
        os.makedirs(os.path.dirname(SPEC_CACHE) or ".", exist_ok=True)
        json.dump(cache, open(SPEC_CACHE, "w", encoding="utf-8"))
    return {w: sorted(v, key=lambda t: t[0]) for w, v in data.items()}


def resolution(pts, order, radii=PROBE_SCAN):
    """Разброс гармоники по радиусу пробы, % от медианы — по верхней половине токов.

    Настоящая гармоника даёт единицы процентов: поделённая на r^(n−1), она от радиуса не зависит.
    Шум даёт десятки и сотни — у него такой зависимости нет."""
    top = pts[len(pts) // 2:]
    sc = []
    for ni, s in top:
        v = np.array([s[r][order] for r in radii if r in s])
        m = np.median(v)
        if m > 0:
            sc.append(100 * (v.max() - v.min()) / m)
    return float(np.mean(sc)) if sc else float("nan")


# ───────────────────────────── графики ─────────────────────────────

def winding_turns():
    """{элемент: {обмотка: витков на катушку}} — из таблицы ветвей насыщения."""
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join("workspace", "tables",
                                              "saturation_branch_members.csv"),
                                 encoding="utf-8-sig"), delimiter=";"):
        if r.get("turns_per_coil"):
            out[r["element"]][r["winding"]] = int(r["turns_per_coil"])
    return out


def plot_element(element, key, outdir):
    """Кривая «коэффициент от тока» с проектными уровнями решётки.

    У многообмоточных элементов (SP, DQ) на каждую обмотку своя кривая: витков на катушку у них
    разное (у DQ 53, 13 и 90), и общая ось тока в амперах для всех сразу бессмысленна — карты
    разных обмоток слились бы в одну ломаную с зубцами."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    turns = winding_turns().get(element, {})
    wds = sorted(turns) or [None]
    series = []
    leff = None
    for w in wds:
        pts, L = curve(element, key, w)
        leff = L or leff
        if pts:
            series.append((w, np.array([p[0] for p in pts]) / (turns.get(w) or 1),
                           np.array([p[1] for p in pts])))
    if not series:
        return None, None

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    wcol = ["#1f4e79", "#c0392b", "#4a9d5f", "#8e5ea8"]
    ymax = max(y.max() for _, _, y in series)
    for i, (w, x, y) in enumerate(series):
        ax.plot(x, y, "-o", color=wcol[i % len(wcol)], ms=5, lw=1.6, zorder=3,
                label=(f"обмотка {w} ({turns.get(w)} витков)" if w and len(series) > 1
                       else f"снятые карты ({len(x)} шт)"))

    colors = ["#e8a33d", "#4a9d5f", "#8e5ea8", "#1a7f8e", "#c0392b"]
    rows = []
    lv = sorted(design_levels(element, key).items(), key=lambda t: t[0][2])
    for i, ((dev, L, v), cnt) in enumerate(lv):
        c = colors[i % len(colors)]
        ax.axhline(v, color=c, lw=1.1, alpha=0.85, zorder=1)
        cross = []
        for w, x, y in series:
            xi = first_crossing(x, y, v)
            if xi is not None:
                ax.plot([xi], [v], "D", color=c, ms=8, zorder=4, mec="white", mew=0.9)
                ax.plot([xi, xi], [0, v], color=c, lw=1.0, ls="--", alpha=0.7, zorder=1)
            cross.append((w, xi))
        rows.append((dev, L, v, cnt, cross))
        txt = "  ".join(f"{w}:{xi:.1f}А" if xi is not None else f"{w}:—"
                        for w, xi in cross) if len(series) > 1 else (
            f"{cross[0][1]:.1f} А" if cross[0][1] is not None else "выше достижимого")
        ax.annotate(f"{dev.replace(element + '_', '')}  {v:.4g}  →  {txt}",
                    xy=(0.015, 0.97 - 0.05 * i), xycoords="axes fraction",
                    fontsize=8.5, color=c, family="monospace", va="top",
                    bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.2))
        ymax = max(ymax, v)

    ax.set_xlabel("ток обмотки, А")
    ax.set_ylabel(f"{key}, {UNITS[key]}")
    ax.set_title(f"{NAMES[element]} ({element}): {key} от тока обмотки")
    ax.grid(alpha=0.25)
    ax.legend(loc="lower right", fontsize=9)
    ax.set_ylim(0, ymax * 1.15)
    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{element}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png, dict(rows=rows, leff=leff, npts=sum(len(x) for _, x, _ in series),
                     nseries=len(series))


CORR_TURNS = {"main": 622, "c1": 548, "c2": 548, "c3": 548, "c4": 548}
CORR_RADIUS = {1: 1e-3, 2: 5e-3, 3: 8e-3, 4: 14e-3}
RESOLVED = 25.0        # % разброса по радиусу пробы, выше которого гармоника считается шумом


def plot_corrector(outdir):
    """По графику на обмотку: все РАЗРЕШИМЫЕ гармоники сразу, каждая в своих единицах.

    Что гармоника разрешима, решает не наше желание её увидеть, а `resolution`: поделённая на
    r^(n−1), настоящая гармоника от радиуса пробы не зависит. У корректора этот отбор проходят
    диполь у всех обмоток и квадруполь у боковых; у `main` квадруполь и все высшие гармоники
    в картах не разрешаются, и рисовать их значило бы рисовать шум."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    orders = (1, 2, 3, 4)
    spec = corrector_spectra(orders)
    cols = {1: "#1f4e79", 2: "#c0392b", 3: "#4a9d5f", 4: "#e8a33d"}
    out = []
    for w in ("main", "c1", "c2", "c3", "c4"):
        pts = spec.get(w)
        if not pts:
            continue
        res = {n: resolution(pts, n) for n in orders}
        show = [n for n in orders if res[n] <= RESOLVED]
        x = np.array([p[0] for p in pts]) / CORR_TURNS[w]
        fig, ax = plt.subplots(figsize=(8.6, 5.0))
        axes = [ax] + [ax.twinx() for _ in show[1:]]
        for k, (n, a) in enumerate(zip(show, axes)):
            y = np.array([p[1][CORR_RADIUS[n]][n] for p in pts])
            a.plot(x, y, "-o", color=cols[n], ms=4.5, lw=1.5)
            a.set_ylabel(f"{ORDER_NAME[n]}, {UNITS[KICK_KEY[n]]}"
                         f"   (проба {CORR_RADIUS[n] * 1e3:.0f} мм)", color=cols[n], fontsize=9)
            a.tick_params(axis="y", labelcolor=cols[n], labelsize=8)
            if k > 1:
                a.spines["right"].set_position(("outward", 52 * (k - 1)))
            a.set_ylim(0, max(y.max() * 1.12, 1e-12))
        ax.set_xlabel(f"ток обмотки, А   ({CORR_TURNS[w]} витков на катушку)")
        ax.set_title(f"Корректор, обмотка {w}: "
                     + ", ".join(ORDER_NAME[n] for n in show) + " от тока")
        ax.grid(alpha=0.25)
        fig.tight_layout()
        os.makedirs(outdir, exist_ok=True)
        png = os.path.join(outdir, f"corrector_{w}.png")
        fig.savefig(png, dpi=150)
        plt.close(fig)
        out.append(dict(w=w, png=png, npts=len(pts), show=show, res=res))
    return out


# ───────────────────────────── документ ─────────────────────────────

def build(out=OUT, figs=PNG):
    from docx import Document
    from docx.shared import Cm, Pt

    doc = Document()
    doc.add_heading("Библиотека магнитных элементов СИЛА: связь токов и киков", 0)
    p = doc.add_paragraph()
    p.add_run("Что здесь. ").bold = True
    p.add_run(
        "На каждый элемент библиотеки приведён паспорт из базы (число витков, число катушек, "
        "номинальный ток, эффективная длина) и кривая зависимости физического коэффициента от "
        "тока обмотки, построенная по снятым картам поля. Горизонтальными линиями отмечены "
        "проектные значения из решётки, ромбом и вертикальной чертой — ток, при котором каждое "
        "из них достигается.")
    p = doc.add_paragraph()
    p.add_run("Как получены коэффициенты. ").bold = True
    p.add_run(
        "Через каждую карту поля прогоняется пучок дифференциальным ядром, и по траекториям "
        "подбирается физический коэффициент при общей для элемента эффективной длине. Это не "
        "гармоническое разложение поля, а именно то, что чувствует пучок.")

    doc.add_heading("Сводка по библиотеке", level=1)
    reg = collections.Counter(r["element"] for r in sm.read_registry())
    t = doc.add_table(rows=1, cols=4); t.style = "Light Grid Accent 1"
    for i, h in enumerate(("Элемент", "Обозначение", "Карт поля", "Обмоток")):
        t.rows[0].cells[i].text = h
    wcount = {"QP_sh": 1, "QP_L": 1, "SP": 2, "OP": 1, "DQ": 3, "corrector": 5}
    for el in ELEMENTS + ["corrector"]:
        c = t.add_row().cells
        c[0].text, c[1].text = NAMES[el], el
        c[2].text, c[3].text = str(reg.get(el, 0)), str(wcount[el])

    for el in ELEMENTS:
        key = KEY[el]
        doc.add_page_break()
        doc.add_heading(f"{NAMES[el]} ({el})", level=1)

        dv = devices(el)
        if dv:
            t = doc.add_table(rows=1, cols=6); t.style = "Light Grid Accent 1"
            for i, h in enumerate(("Исполнение", "Витков на катушку", "Катушек",
                                   "Всего витков", "Номинальный ток, А", "L_eff, м")):
                t.rows[0].cells[i].text = h
            for d in dv:
                c = t.add_row().cells
                for i, v in enumerate((d["dev"], d["turns"], d["coils"], d["total"],
                                       d["inom"], d["leff"])):
                    c[i].text = str(v)

        png, info = plot_element(el, key, figs)
        if png:
            doc.add_picture(png, width=Cm(16.5))
            if info["leff"]:
                q = doc.add_paragraph()
                q.add_run("Эффективная длина подгонки: ").bold = True
                q.add_run(f"{info['leff']:.4f} м (одна на все {info['npts']} карт элемента). ")
                fam = {round(float(d["leff"]), 4) for d in dv if d["leff"]}
                if fam and all(abs(info["leff"] - f) > 0.005 for f in fam):
                    q.add_run(
                        f"Внимание: она не совпадает с проектными длинами исполнений "
                        f"({', '.join(f'{f:.3f}' for f in sorted(fam))} м). Коэффициент равен "
                        f"интегралу поля, делённому на длину, поэтому при разных длинах "
                        f"величины прямо не сопоставимы.").italic = True
            if info["rows"]:
                nom_w = nominal_by_winding(el)
                t = doc.add_table(rows=1, cols=6); t.style = "Light Grid Accent 1"
                for i, h in enumerate(("Исполнение", "L_eff, м", f"Проектный {key}, {UNITS[key]}",
                                       f"Проектный кик ∫{key}·dl, {INT_UNITS[key]}",
                                       "Штук в кольце", "Ток обмотки, А")):
                    t.rows[0].cells[i].text = h
                for dev, L, v, cnt, cross in info["rows"]:
                    c = t.add_row().cells
                    c[0].text = dev.replace(el + "_", "")
                    c[1].text = f"{L:.4f}"
                    c[2].text = f"{v:.5g}"
                    c[3].text = f"{v * L:.5g}"
                    c[4].text = str(cnt)
                    parts = []
                    for w, xi in cross:
                        pre = f"{w}: " if w and info["nseries"] > 1 else ""
                        if xi is None:
                            parts.append(pre + "не достигается")
                            continue
                        nom = nom_w.get(w) or nom_w.get(None)
                        parts.append(pre + f"{xi:.1f}"
                                     + (f"  ({100 * xi / nom:.0f}% ном.)" if nom else ""))
                    c[5].text = ";  ".join(parts)
                if info["nseries"] > 1:
                    doc.add_paragraph(
                        "Ток указан для случая, когда проектное значение обеспечивает ОДНА "
                        "обмотка; в кольце они работают вместе. Поэтому у подстроечных обмоток "
                        "(SP I2, DQ I3) процент от номинала выходит много больше ста: в одиночку "
                        "нести всю силу они и не рассчитаны.").italic = True

    doc.add_page_break()
    doc.add_heading(f"{NAMES['corrector']} (corrector)", level=1)
    p = doc.add_paragraph()
    p.add_run(
        "Корректор пятиобмоточный, и его обмотки создают не одну величину, а смесь гармоник: "
        "главная обмотка — почти чистый диполь, а каждая боковая даёт и диполь, и квадруполь. "
        "Поэтому на каждую обмотку приведён отдельный график со всеми гармониками сразу — у "
        "каждой своя шкала, так как единицы разные.")
    p = doc.add_paragraph()
    p = doc.add_paragraph()
    p.add_run("Проектного кика у корректора в решётке нет: ").bold = True
    p.add_run(
        "все 120 корректоров решётки 1 стоят с нулевыми коэффициентами — они не оттюнены, и "
        "рабочая точка им ещё не назначена. Поэтому здесь нет ни проектного уровня, ни тока под "
        "него; приведена сама зависимость кика от тока, по которой ток и будет определяться, "
        "когда рабочая точка появится.")
    p = doc.add_paragraph()
    p.add_run("Какие гармоники показаны. ").bold = True
    p.add_run(
        "Поле мультиполя растёт как r в степени n−1, а шум карты от радиуса не зависит, поэтому "
        "слабая высокая гармоника у оси попросту тонет. Проверка простая: настоящая гармоника, "
        "поделённая на r в степени n−1, от радиуса пробы НЕ зависит. Ниже приведён разброс "
        "каждой гармоники по пробам от 1 до 14 мм; на графики вынесены те, у кого он не "
        "превышает "
        f"{RESOLVED:.0f}%. Остальные в имеющихся картах не разрешаются, и рисовать их значило бы "
        "рисовать шум.")
    dv = devices("corrector")
    if dv:
        t = doc.add_table(rows=1, cols=6); t.style = "Light Grid Accent 1"
        for i, h in enumerate(("Обмотка", "Витков на катушку", "Катушек", "Всего витков",
                               "Номинальный ток, А", "L_eff, м")):
            t.rows[0].cells[i].text = h
        for d in dv:
            c = t.add_row().cells
            w = d["dev"].split("_")[1]
            for i, v in enumerate((w, d["turns"], d["coils"], d["total"], d["inom"], d["leff"])):
                c[i].text = str(v)
    cr = plot_corrector(figs)
    t = doc.add_table(rows=1, cols=6); t.style = "Light Grid Accent 1"
    for i, h in enumerate(("Обмотка", "Карт") + tuple(ORDER_NAME[n] for n in (1, 2, 3, 4))):
        t.rows[0].cells[i].text = h
    for r in cr:
        c = t.add_row().cells
        c[0].text, c[1].text = r["w"], str(r["npts"])
        for i, n in enumerate((1, 2, 3, 4)):
            c[2 + i].text = (f"{r['res'][n]:.0f}%"
                             + (" — на графике" if n in r["show"] else " — не разрешается"))
    doc.add_paragraph("Разброс гармоники по радиусу пробы, % от медианы (по верхней половине "
                      "токов). Единицы: B — Тл, G — Тл/м, S — Тл/м², O — Тл/м³.").italic = True
    for r in cr:
        doc.add_heading(f"Обмотка {r['w']}", level=2)
        doc.add_picture(r["png"], width=Cm(16.5))
        doc.add_paragraph(
            f"Карт: {r['npts']}. Показаны: "
            + ", ".join(ORDER_NAME[n] for n in r["show"])
            + ("." if len(r["show"]) == 4 else
               ". Не разрешаются: " + ", ".join(ORDER_NAME[n] for n in (1, 2, 3, 4)
                                                if n not in r["show"]) + "."))

    doc.add_page_break()
    doc.add_heading("Оговорки", level=1)
    for s in (
        "Эффективная длина у большинства элементов ПОДБИРАЕТСЯ по траекториям, одна на элемент, и "
        "не обязана совпадать с проектной длиной исполнения. Коэффициент равен интегралу поля, "
        "делённому на эту длину, поэтому при расхождении длин сравнение с проектом требует "
        "приведения.",
        "Ток на графиках — ток ОДНОЙ катушки в амперах: ампер-витки карты делятся на число витков "
        "катушки. Возбуждение карты задаётся в ампер-витках.",
        "Кривые построены по снятым картам без сглаживания и интерполяции; точки — то, что "
        "посчитано.",
        "У короткого квадруполя выше ~13000 ампер-витков кривая выходит на плато: ток растёт "
        "втрое, а градиент меняется в пределах полупроцента. Ток по заданному градиенту в этой "
        "области определяется плохо.",
    ):
        doc.add_paragraph(s, style="List Bullet")

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    doc.save(out)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    print("отчёт ->", build(a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
