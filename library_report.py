# -*- coding: utf-8 -*-
"""Word-отчёт по библиотеке элементов: паспорт, кики каждой обмотки, проектные рабочие точки.

Устройство отчёта: на КАЖДУЮ ОБМОТКУ свой график, и на нём ВСЕ кики, которые эта обмотка создаёт —
диполь, квадруполь, секступоль, октуполь. Обмотка редко создаёт одну чистую гармонику: боковая
обмотка корректора даёт и диполь, и квадруполь, а диполь-квадруполь на то и combined-function.

По вертикали отложен КИК — интеграл гармоники вдоль элемента, а не коэффициент. Пучок чувствует
именно интеграл; коэффициент получается делением его на эффективную длину, и это соглашение, а не
измерение (см. `map_length`). Единицы физические, при проектной энергии 6 ГэВ.

Кики берутся гармоническим разложением карт. Сверено с независимым способом — подгонкой по
траекториям через диффурное ядро (`saturation_model`): совпадение 0.1–0.8%, см. `crosscheck`.

    py -3.14 library_report.py
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

OUT = os.path.join("workspace", "results", "library_report.docx")
PNG = os.path.join("workspace", "results", "library_report_figs")

ELEMENTS = ["QP_sh", "QP_L", "SP", "OP", "DQ", "corrector"]
KEY = {"QP_sh": "G", "QP_L": "G", "SP": "S", "OP": "O", "DQ": "G", "corrector": "B"}
COLUMN = {"B": "mat_B", "G": "mat_G", "S": "mat_S", "O": "mat_O"}
NAMES = {"QP_sh": "Квадруполь короткий", "QP_L": "Квадруполь длинный",
         "SP": "Секступоль", "OP": "Октуполь", "DQ": "Диполь-квадруполь",
         "corrector": "Корректор"}
ORDERS = (1, 2, 3, 4)
ORDER_NAME = {1: "диполь", 2: "квадруполь", 3: "секступоль", 4: "октуполь"}
KICK_KEY = {1: "B", 2: "G", 3: "S", 4: "O"}
ORDER_OF = {v: k for k, v in KICK_KEY.items()}
# Единицы ИНТЕГРАЛА гармоники по длине — это и есть кик. By ∝ b_n·r^(n−1), отсюда степени метра.
INT_UNITS = {"B": "Тл·м", "G": "Тл", "S": "Тл/м", "O": "Тл/м²"}
ORDER_COLOR = {1: "#1f4e79", 2: "#c0392b", 3: "#4a9d5f", 4: "#e8a33d"}

PROBE_SCAN = (1e-3, 3e-3, 5e-3, 8e-3, 11e-3, 14e-3)  # карты все ±50 мм, 14 мм внутри хор. поля
PLOT_RADIUS = {1: 1e-3, 2: 5e-3, 3: 8e-3, 4: 14e-3}
# Радиусы, НА КОТОРЫХ проверяется устойчивость данного порядка. Проверять октуполь на 1 мм
# бессмысленно: там он заведомо под шумовым полом, и настоящая гармоника отбраковалась бы вместе
# с шумом — у секступоля SP, его ГЛАВНОЙ гармоники, так и выходило 56%.
RES_RADII = {1: (1e-3, 3e-3, 5e-3, 8e-3), 2: (1e-3, 3e-3, 5e-3, 8e-3),
             3: (5e-3, 8e-3, 11e-3, 14e-3), 4: (8e-3, 11e-3, 14e-3)}
RESOLVED = 25.0        # % разброса по радиусу пробы, выше которого гармоника считается шумом
SKEW_SHOW = 5.0        # % от нормальной, начиная с которых скошенную компоненту стоит рисовать

SPEC_CACHE = os.path.join("workspace", "results", "library_report_spectra_ri.json")
MAPLEN_CACHE = os.path.join("workspace", "results", "library_report_maplen.json")


# ───────────────────────────── паспорт из БД ─────────────────────────────

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


def winding_turns():
    """{элемент: {обмотка: витков на катушку}} — из таблицы ветвей насыщения."""
    out = collections.defaultdict(dict)
    for r in csv.DictReader(open(os.path.join("workspace", "tables",
                                              "saturation_branch_members.csv"),
                                 encoding="utf-8-sig"), delimiter=";"):
        if r.get("turns_per_coil"):
            out[r["element"]][r["winding"]] = int(r["turns_per_coil"])
    return out


def nominal_by_winding(element):
    """{обмотка: номинальный ток, А} — из имён устройств вида `SP_I1_1A_L166`.

    У однообмоточных элементов обмотки в имени нет, и номинал ложится под ключ None."""
    out = {}
    for d in devices(element):
        tok = d["dev"][len(element) + 1:].split("_")[0]
        try:
            out[tok if not tok.endswith("A") else None] = float(d["inom"].split()[0])
        except (ValueError, IndexError):
            pass
    return out


def energy_factor():
    """Множитель «шкала трекера → физическая»: трекинг идёт на суррогатной энергии γ≈7, а карты
    и паспорта — на проектных 6 ГэВ."""
    import element_compare as ec
    return ec.ENERGY_FACTOR


def design_levels(element, key):
    """Проектные кики решётки 1 в ФИЗИЧЕСКИХ единицах: {(устройство, L, ∫кик·dl): количество}.

    В `lattice_configs` лежит коэффициент в шкале трекера. Физический кик получается умножением на
    длину, на энергетический множитель и на (n−1)! — таково соглашение столбцов mat_*: у них
    By ∝ mat_S·x²/2!, а у гармонического разложения By ∝ b₃·x²."""
    col, n = COLUMN[key], ORDER_OF[key]
    scale = energy_factor() / math.factorial(n - 1)
    c = collections.Counter()
    for r in csv.DictReader(open(os.path.join("workspace", "tables", "lattice_configs.csv"),
                                 encoding="utf-8-sig"), delimiter=";"):
        if str(r.get("lattice_id")) != "1" or not r.get("device_type_id", "").startswith(element):
            continue
        try:
            v, L = abs(float(r.get(col) or 0)), float(r.get("eff_length") or 0)
        except ValueError:
            continue
        if v > 1e-9 and L > 0:
            c[(r["device_type_id"], round(L, 4), v * L * scale)] += 1
    tot = sum(c.values()) or 1
    return {k: m for k, m in c.items() if m >= max(3, 0.02 * tot)}


def first_crossing(x, y, level):
    """Ток, при котором кик ВПЕРВЫЕ достигает уровня. В насыщении кривая может быть немонотонна,
    поэтому обычная интерполяция непригодна — берём первое пересечение."""
    for i in range(1, len(y)):
        if (y[i - 1] - level) * (y[i] - level) <= 0 and y[i] != y[i - 1]:
            t = (level - y[i - 1]) / (y[i] - y[i - 1])
            return float(x[i - 1] + t * (x[i] - x[i - 1]))
    return None


# ───────────────────────────── кики из карт ─────────────────────────────

def spectrum(path, orders, radii, nth=64):
    """Кики на НЕСКОЛЬКИХ радиусах пробы за одну загрузку карты: {радиус: {порядок: (норм, скью)}}.

    Гармоника комплексна: (By + i·Bz) даёт нормальный мультиполь и повёрнутый на пол-периода
    (скошенный). Брать МОДУЛЬ нельзя — он их смешивает: у боковых обмоток корректора скошенная
    часть сравнима с нормальной, и модуль завышал дипольный кик на 59%. У чистых элементов
    скошенная часть — доли процента, и разложение это подтверждает.

    Множитель 2 при интеграле — потому что карты половинные: поле максимально при x=0 (середина
    магнита) и спадает к дальнему концу, вторую половину даёт симметрия.

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
            out[r][n] = (abs(c.imag) / r ** (n - 1), abs(c.real) / r ** (n - 1))
    return out


def element_spectra(element, orders=ORDERS, radii=PROBE_SCAN):
    """{обмотка: [(ампер-витки, {радиус: {порядок: кик}})]} по ОДИНОЧНЫМ картам, с кэшем на диске.

    Карты с несколькими одновременно запитанными обмотками пропускаются: по ним нельзя сказать,
    какая обмотка какой кик создала, а вопрос отчёта именно такой."""
    cache = {}
    if os.path.exists(SPEC_CACHE):
        try:
            cache = json.load(open(SPEC_CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    dirty = False
    data = collections.defaultdict(list)
    for r in sm.read_registry():
        if r["element"] != element or ";" in r["excitation"] or "+" in r["excitation"]:
            continue
        p = os.path.join("workspace", "fields_bin", r["bin_path"])
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not os.path.exists(p):
            continue
        key = r["map_id"]
        have = cache.setdefault(key, {})
        # Досчитываем ТОЛЬКО недостающие радиусы: иначе добавление радиуса в скан либо молча
        # оставляло бы старый набор, либо заставляло пересчитывать всю библиотеку заново.
        need = [rr for rr in radii if f"{rr:.6f}" not in have]
        if need:
            s = spectrum(p, orders, need)
            have.update({f"{rr:.6f}": {str(n): list(v) for n, v in d.items()}
                         for rr, d in s.items()})
            dirty = True
        s = {float(rr): {int(n): tuple(v) for n, v in d.items()} for rr, d in have.items()}
        data[r["winding"]].append((ni, s))
    if dirty:
        os.makedirs(os.path.dirname(SPEC_CACHE) or ".", exist_ok=True)
        json.dump(cache, open(SPEC_CACHE, "w", encoding="utf-8"))
    return {w: sorted(v, key=lambda t: t[0]) for w, v in data.items()}


def resolution(pts, order, radii=None):
    """Разброс кика по радиусу пробы, % от медианы — по верхней половине токов.

    Настоящая гармоника даёт единицы процентов: поделённая на r^(n−1), она от радиуса не зависит.
    Шум даёт десятки и сотни — у него такой зависимости нет. Радиусы берутся из RES_RADII: для
    каждого порядка свои, те, на которых гармоника этого порядка вообще может быть измерена."""
    radii = radii or RES_RADII[order]
    top = pts[len(pts) // 2:]
    sc = []
    for ni, s in top:
        v = np.array([s[r][order][0] for r in radii if r in s])
        m = np.median(v)
        if m > 0:
            sc.append(100 * (v.max() - v.min()) / m)
    return float(np.mean(sc)) if sc else float("nan")


def map_length(element, key):
    """Жёстко-краевая эффективная длина ИЗМЕРЕННОГО магнита: 2·∫c·dx / max c (карта половинная).

    Это НЕ подгоночная L_eff из датасетов. Подгонка определяет только ПРОИЗВЕДЕНИЕ «коэффициент ×
    длина» (оно сходится с интегралом карты до 0.8%), а как разделить его на сомножители — не
    определяет вовсе: при смене сетки подгоночная длина QP_sh уехала на 13%, у OP она ниже
    настоящей на 39%. Длина же, взятая из ФОРМЫ поля, попадает в проектную: QP_L 0.374 против
    0.3745, DQ 0.367 против 0.400, SP 0.164 против 0.166."""
    cache = {}
    if os.path.exists(MAPLEN_CACHE):
        try:
            cache = json.load(open(MAPLEN_CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    ck = f"{element}|{key}"
    if ck in cache:
        return cache[ck]
    from field_quality import load, _circle
    n = ORDER_OF[key]
    r = PLOT_RADIUS[n]
    reg = [x for x in sm.read_registry()
           if x["element"] == element and ";" not in x["excitation"] and "+" not in x["excitation"]]
    if not reg:
        return None
    # Средний ток: у самых слабых карт всё тонет в шуме, у самых сильных профиль уже искажён
    # насыщением, и жёстко-краевая длина заметно расползается.
    m = sorted(reg, key=lambda x: abs(float(x["amp_turns"] or 0)))[len(reg) // 2]
    h, B = load(os.path.join("workspace", "fields_bin", m["bin_path"]))
    xs = np.linspace(h["mnx"], h["mxx"], h["nx"])
    v = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, r)
        sp = np.fft.fft(by + 1j * bz) / 64
        v.append(abs((sp[0] if n == 1 else sp[-(n - 1)]).imag))
    v = np.array(v)
    L = float(2 * np.trapezoid(v, xs) / v.max())
    cache[ck] = L
    os.makedirs(os.path.dirname(MAPLEN_CACHE) or ".", exist_ok=True)
    json.dump(cache, open(MAPLEN_CACHE, "w", encoding="utf-8"))
    return L


def crosscheck(element, key):
    """Сверка кика из гармоник с киком из подгонки траекторий на самой сильной карте элемента.

    Два независимых способа: гармоническое разложение поля и прогон пучка диффурным ядром через ту
    же карту. Совпадают — значит ни в разложении, ни в подгонке нет систематики."""
    p = os.path.join(sm.OUTDIR, f"{element}.json")
    if not os.path.exists(p):
        return None
    d = json.load(open(p, encoding="utf-8"))
    rec = [x for x in d["records"]
           if x["coeffs"] and key in x["coeffs"] and "+" not in (x["winding"] or "")]
    if not rec:
        return None
    x = sorted(rec, key=lambda z: abs(float(z["amp_turns"])))[-1]
    reg = {r["map_id"]: r for r in sm.read_registry() if r["element"] == element}
    m = reg.get(x["map_id"])
    if m is None:
        return None
    n = ORDER_OF[key]
    r = PLOT_RADIUS[n]
    harm = spectrum(os.path.join("workspace", "fields_bin", m["bin_path"]), (n,), (r,))[r][n][0]
    traj = abs(x["coeffs"][key]) * x["Leff"] * energy_factor() / math.factorial(n - 1)
    return dict(ni=abs(float(x["amp_turns"])), key=key, harm=harm, traj=traj,
                diff=100 * (traj / harm - 1))


# ───────────────────────────── графики ─────────────────────────────

def plot_winding(element, winding, pts, turns, outdir):
    """Один график = одна обмотка. На нём ВСЕ кики, которые она создаёт, каждый в своей шкале.

    Что кик реален, решает не наше желание его увидеть, а `resolution`. Нарисовать шум под именем
    октуполя хуже, чем не рисовать ничего: по такой кривой будут считать ток."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    res = {n: resolution(pts, n) for n in ORDERS}
    show = [n for n in ORDERS if res[n] <= RESOLVED]
    if not show:
        return None
    x = np.array([p[0] for p in pts]) / (turns or 1)
    Lmap = map_length(element, KEY[element])

    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    axes = {n: (ax if k == 0 else ax.twinx()) for k, n in enumerate(show)}
    rows, notes = [], []
    skew = {}
    for k, n in enumerate(show):
        a, col = axes[n], ORDER_COLOR[n]
        y = np.array([p[1][PLOT_RADIUS[n]][n][0] for p in pts])
        ys = np.array([p[1][PLOT_RADIUS[n]][n][1] for p in pts])
        a.plot(x, y, "-o", color=col, ms=4.5, lw=1.6, zorder=3)
        # Скошенная компонента — это ТОТ ЖЕ мультиполь, повёрнутый на пол-периода: он связывает
        # плоскости, а не фокусирует в них. У чистых элементов её доли процента, у боковых обмоток
        # корректора она сравнима с нормальной, и умолчать о ней значило бы соврать про кик.
        skew[n] = 100 * float(np.mean(ys[len(ys) // 2:] / np.maximum(y[len(y) // 2:], 1e-30)))
        if skew[n] >= SKEW_SHOW:
            a.plot(x, ys, "--s", color=col, ms=3.5, lw=1.3, alpha=0.75, zorder=2)
        a.set_ylabel(f"∫{KICK_KEY[n]}·dl, {INT_UNITS[KICK_KEY[n]]}"
                     f"   (проба {PLOT_RADIUS[n] * 1e3:.0f} мм)", color=col, fontsize=9)
        a.tick_params(axis="y", labelcolor=col, labelsize=8)
        if k > 1:
            a.spines["right"].set_position(("outward", 58 * (k - 1)))
        ymax = y.max()

        for (dev, L, v), cnt in sorted(design_levels(element, KICK_KEY[n]).items(),
                                       key=lambda t: t[0][2]):
            same = Lmap and abs(L - Lmap) / Lmap <= 0.10
            a.axhline(v, color=col, lw=1.0, alpha=0.75, ls="-" if same else ":", zorder=1)
            xi = first_crossing(x, y, v)
            if xi is not None:
                a.plot([xi], [v], "D", color=col, ms=7.5, zorder=4, mec="white", mew=0.9)
                a.plot([xi, xi], [0, v], color=col, lw=0.9, ls="--", alpha=0.6, zorder=1)
            rows.append(dict(order=n, dev=dev, L=L, level=v, count=cnt, cur=xi, same=same))
            notes.append((col, f"{ORDER_NAME[n]}  {dev.replace(element + '_', '')}  {v:.4g}  →  "
                               + (f"{xi:.1f} А" if xi is not None else "не достигается")
                               + ("" if same else "  (магнит другой длины)")))
            # Уровень чужой длины НЕ растягивает шкалу: он недостижим по построению, и из-за него
            # настоящая кривая сплющивалась бы в нижнюю треть графика.
            if same:
                ymax = max(ymax, v)
        a.set_ylim(0, max(ymax, ys.max() if skew[n] >= SKEW_SHOW else 0) * 1.15)

    for i, (col, txt) in enumerate(notes):
        ax.annotate(txt, xy=(0.015, 0.97 - 0.052 * i), xycoords="axes fraction", fontsize=8.5,
                    color=col, family="monospace", va="top",
                    bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.2))

    ax.set_xlabel(f"ток обмотки, А   ({turns} витков на катушку)")
    sub = [f"длина магнита в картах {Lmap:.3f} м"] if Lmap else []
    if any(v >= SKEW_SHOW for v in skew.values()):
        sub.append("пунктир — скошенная компонента")
    ax.set_title(f"{NAMES[element]} ({element}), обмотка {winding}: "
                 + ", ".join(ORDER_NAME[n] for n in show)
                 + ("\n" + ";  ".join(sub) if sub else ""), fontsize=11)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{element}_{winding}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return dict(w=winding, png=png, npts=len(pts), show=show, res=res, rows=rows, Lmap=Lmap,
                skew=skew)


def plot_element(element, outdir):
    """Все обмотки элемента: по графику на каждую."""
    spec = element_spectra(element)
    turns = winding_turns().get(element, {})
    out = []
    for w in sorted(spec, key=lambda z: (z != "main", z)):
        r = plot_winding(element, w, spec[w], turns.get(w), outdir)
        if r:
            out.append(r)
    return out


# ───────────────────────────── документ ─────────────────────────────

def _table(doc, headers):
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Light Grid Accent 1"
    for i, h in enumerate(headers):
        t.rows[0].cells[i].text = h
    return t


def build(out=OUT, figs=PNG):
    from docx import Document
    from docx.shared import Cm

    doc = Document()
    doc.add_heading("Библиотека магнитных элементов СИЛА: связь токов и киков", 0)

    p = doc.add_paragraph()
    p.add_run("Что здесь. ").bold = True
    p.add_run(
        "На каждую ОБМОТКУ каждого элемента приведён свой график, и на нём все кики, которые эта "
        "обмотка создаёт: диполь, квадруполь, секступоль, октуполь. Обмотка редко создаёт одну "
        "чистую гармонику — боковая обмотка корректора даёт и диполь, и квадруполь, а "
        "диполь-квадруполь на то и combined-function. Горизонтальными линиями отмечены проектные "
        "кики из решётки, ромбом и вертикальной чертой — ток, при котором каждый достигается.")

    p = doc.add_paragraph()
    p.add_run("Почему кик, а не коэффициент. ").bold = True
    p.add_run(
        "Пучок чувствует ИНТЕГРАЛ поля вдоль элемента — им и определяется угол поворота. "
        "Коэффициент получается делением интеграла на эффективную длину, и это соглашение, а не "
        "измерение. Подгонка траекторий определяет произведение «коэффициент × длина» с точностью "
        "меньше процента, а как разделить его на сомножители — не определяет: при смене сетки "
        "подгоночная длина короткого квадруполя сместилась на 13% при неизменном произведении. "
        "Поэтому по вертикали отложен интеграл — он измеряется, и длина в нём не участвует.")

    p = doc.add_paragraph()
    p.add_run("Единицы. ").bold = True
    p.add_run(
        f"Физические, при проектной энергии 6 ГэВ. Кик порядка n — это ∫b_n·dl, где поле "
        f"раскладывается как By ∝ b_n·r^(n−1): диполь в Тл·м, квадруполь в Тл, секступоль в Тл/м, "
        f"октуполь в Тл/м². Коэффициенты решётки хранятся в шкале трекера (суррогатная энергия "
        f"γ≈7) и здесь переведены в физические умножением на {energy_factor():.1f}. Ток по "
        f"горизонтали — ток ОДНОЙ катушки: ампер-витки карты, делённые на число витков катушки.")

    doc.add_heading("Сводка по библиотеке", level=1)
    reg = collections.Counter(r["element"] for r in sm.read_registry())
    wt = winding_turns()
    t = _table(doc, ("Элемент", "Обозначение", "Карт поля", "Обмоток", "Обмотки"))
    for el in ELEMENTS:
        c = t.add_row().cells
        c[0].text, c[1].text = NAMES[el], el
        c[2].text, c[3].text = str(reg.get(el, 0)), str(len(wt.get(el, {})))
        c[4].text = ", ".join(sorted(wt.get(el, {})))

    doc.add_heading("Сверка двух независимых способов", level=1)
    doc.add_paragraph(
        "Кики на графиках получены гармоническим разложением карт. Тот же кик можно получить "
        "иначе — прогнать через карту пучок диффурным ядром и подобрать силу по траекториям. Ниже "
        "сравнение обоих способов на самой сильной карте каждого элемента. Совпадение означает, "
        "что систематики нет ни в разложении, ни в подгонке. Сравнивается НОРМАЛЬНАЯ компонента: "
        "гармоника комплексна, и её модуль смешал бы нормальную часть со скошенной — у боковых "
        "обмоток корректора это давало ложные 59% расхождения.")
    t = _table(doc, ("Элемент", "Карта, А·в", "Кик", "Из гармоник", "Из траекторий", "Расхождение"))
    for el in ELEMENTS:
        cc = crosscheck(el, KEY[el])
        if not cc:
            continue
        c = t.add_row().cells
        for i, v in enumerate((NAMES[el], f"{cc['ni']:.0f}",
                               f"∫{cc['key']}·dl, {INT_UNITS[cc['key']]}",
                               f"{cc['harm']:.6g}", f"{cc['traj']:.6g}", f"{cc['diff']:+.2f}%")):
            c[i].text = v

    for el in ELEMENTS:
        doc.add_page_break()
        doc.add_heading(f"{NAMES[el]} ({el})", level=1)

        dv = devices(el)
        if dv:
            t = _table(doc, ("Исполнение", "Витков на катушку", "Катушек", "Всего витков",
                             "Номинальный ток, А", "Проектная L_eff, м"))
            for d in dv:
                c = t.add_row().cells
                for i, v in enumerate((d["dev"], d["turns"], d["coils"], d["total"],
                                       d["inom"], d["leff"])):
                    c[i].text = str(v)

        res = plot_element(el, figs)
        if not res:
            continue

        Lmap = res[0]["Lmap"]
        if Lmap:
            q = doc.add_paragraph()
            q.add_run("Длина магнита в картах: ").bold = True
            q.add_run(f"{Lmap:.3f} м (по форме поля: удвоенный интеграл продольного профиля, "
                      f"делённый на его максимум). ")
            fam = sorted({round(float(d["leff"]), 4) for d in dv if d["leff"]})
            near = [f for f in fam if abs(f - Lmap) / Lmap <= 0.10]
            if near:
                q.add_run("Совпадает с проектной длиной исполнения "
                          + ", ".join(f"{f:.3f}" for f in near) + " м.")
            elif fam:
                q.add_run("Проектные длины исполнений: " + ", ".join(f"{f:.3f}" for f in fam)
                          + " м — ни одна не совпала, рабочие точки ниже условны.").italic = True
            if len(fam) > 1:
                doc.add_paragraph(
                    "В решётке этот магнит стоит в нескольких длинах, а карты сняты только с "
                    "одной. Уровни исполнений другой длины показаны пунктиром и помечены в "
                    "таблице: их кик наш магнит не обязан достигать — он просто другой длины.",
                    style="List Bullet")

        t = _table(doc, ("Обмотка", "Карт") + tuple(ORDER_NAME[n] for n in ORDERS))
        for r in res:
            c = t.add_row().cells
            c[0].text, c[1].text = r["w"], str(r["npts"])
            for i, n in enumerate(ORDERS):
                c[2 + i].text = (f"{r['res'][n]:.0f}%"
                                 + (" — на графике" if n in r["show"] else " — не разрешается"))
        doc.add_paragraph(
            "Разброс кика по радиусу пробы, % от медианы (по верхней половине токов). Настоящая "
            "гармоника, поделённая на r^(n−1), от радиуса не зависит; шум зависит. На графики "
            f"вынесены кики с разбросом не более {RESOLVED:.0f}%.", style="List Bullet")

        nom_w = nominal_by_winding(el)
        for r in res:
            doc.add_heading(f"Обмотка {r['w']}", level=2)
            doc.add_picture(r["png"], width=Cm(16.5))
            doc.add_paragraph(
                f"Карт: {r['npts']}. Создаёт: " + ", ".join(ORDER_NAME[n] for n in r["show"])
                + ("." if len(r["show"]) == len(ORDERS) else
                   ". Не разрешаются в картах: "
                   + ", ".join(ORDER_NAME[n] for n in ORDERS if n not in r["show"]) + "."))
            sk = [n for n in r["show"] if r["skew"][n] >= SKEW_SHOW]
            if sk:
                doc.add_paragraph(
                    "Скошенная компонента (тот же мультиполь, повёрнутый на пол-периода — он "
                    "связывает плоскости, а не фокусирует в них): "
                    + ", ".join(
                        f"{ORDER_NAME[n]} — "
                        + (f"в {r['skew'][n] / 100:.1f} раза БОЛЬШЕ нормальной"
                           if r["skew"][n] > 100 else f"{r['skew'][n]:.0f}% от нормальной")
                        for n in sk)
                    + ". На графике пунктиром.")
            if not r["rows"]:
                doc.add_paragraph("Проектного кика для этой обмотки в решётке нет.")
                continue
            t = _table(doc, ("Кик", "Исполнение", "Проектная L_eff, м", "Проектный кик",
                             "Штук в кольце", "Карты этой длины", "Ток обмотки, А"))
            nom = nom_w.get(r["w"]) or nom_w.get(None)
            for q in r["rows"]:
                c = t.add_row().cells
                cur = ("не достигается" if q["cur"] is None else
                       f"{q['cur']:.1f}" + (f"  ({100 * q['cur'] / nom:.0f}% ном.)" if nom else ""))
                for i, v in enumerate((
                        f"∫{KICK_KEY[q['order']]}·dl, {INT_UNITS[KICK_KEY[q['order']]]}",
                        q["dev"].replace(el + "_", ""), f"{q['L']:.4f}", f"{q['level']:.5g}",
                        str(q["count"]), "да" if q["same"] else "нет", cur)):
                    c[i].text = v

    doc.add_page_break()
    doc.add_heading("Оговорки", level=1)
    for s in (
        "Длина магнита взята не из подгонки, а из ФОРМЫ поля: удвоенный интеграл продольного "
        "профиля, делённый на его максимум (карты половинные — поле максимально в середине "
        "магнита при x=0, вторую половину даёт симметрия). Так измеренная длина сходится с "
        "проектной: длинный квадруполь 0.374 против проектных 0.3745, секступоль 0.164 против "
        "0.166. Подгоночная эффективная длина из датасетов не использована: она определена плохо, "
        "у октуполя оказалась на 39% ниже настоящей.",
        "Показаны только кики, разрешимые в имеющихся картах. Поле мультиполя растёт как r^(n−1), "
        "а шум карты от радиуса не зависит, поэтому слабая высокая гармоника у оси тонет. "
        "Неразрешимый кик не значит нулевой — значит, что по этим картам его величину назвать "
        "нельзя. Разброс по радиусу приведён таблицей, вывод проверяем.",
        "Карты с несколькими одновременно запитанными обмотками в графики не входят: по ним "
        "нельзя сказать, какая обмотка какой кик создала. У корректора таких карт больше половины "
        "набора, и они используются отдельно — для проверки модели на комбинациях.",
        "Ток указан для случая, когда проектный кик обеспечивает ОДНА обмотка; в кольце они "
        "работают вместе. Поэтому у подстроечных обмоток (SP I2, DQ I3) процент от номинала "
        "выходит много больше ста: в одиночку нести всю силу они не рассчитаны.",
        "Кривые построены по снятым картам без сглаживания и интерполяции; точки — то, что "
        "посчитано.",
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
