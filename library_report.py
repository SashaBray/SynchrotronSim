# -*- coding: utf-8 -*-
"""Достоверная копия `element_compare_report`, дополненная двумя вещами и только ими.

Исходный отчёт НЕ пересобирается: он открывается как есть, и в конец каждой главы элемента
добавляются кривые намагниченности обмоток (киковые моменты в зависимости от тока питания, по
трёхмерным картам поля) и оценка проектных токов С УЧЁТОМ НАСЫЩЕНИЯ. Всё остальное — вводный
раздел, матрицы связи, таблицы свипа, траектории, проектные токи, сводная таблица диапазонов —
остаётся нетронутым.

    py -3.14 library_report.py
"""
from __future__ import annotations

import argparse
import collections
import json
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

SRC = os.path.join("workspace", "results", "element_compare", "element_compare_report.docx")
OUT = os.path.join("workspace", "results", "library_report.docx")
PNG = os.path.join("workspace", "results", "library_report_figs")

# Киковый момент → (порядок мультиполя, компонента): 0 — нормальная, 1 — скью. Обозначения те же,
# что в `element_compare.MULTIPOLES`.
MOMENT = {"B": (1, 0), "Bv": (1, 1), "G": (2, 0), "Gv": (2, 1), "S": (3, 0), "O": (4, 0)}
# Единицы ИНТЕГРАЛА момента по длине: By ∝ b_n·r^(n−1), отсюда степени метра.
MOMENT_UNITS = {"B": "Тл·м", "Bv": "Тл·м", "G": "Тл", "Gv": "Тл", "S": "Тл/м", "O": "Тл/м²"}
MOMENT_COLOR = {"B": "#1f4e79", "Bv": "#5b9bd5", "G": "#c0392b", "Gv": "#e07b6a",
                "S": "#4a9d5f", "O": "#e8a33d"}
ELEMENTS = ("QP_sh", "QP_L", "SP", "OP", "DQ", "corrector")   # элементы, у которых есть карты полей
# Моменты, рассматриваемые у элементов данного вида. У диполя-квадруполя это B, G и S: секступоль
# там несущий (через feed-down он и объясняет связь G с током), а октуполь в модель не входит.
KIND_MOMENTS = {"dq": ("B", "G", "S")}

PROBE_SCAN = (1e-3, 3e-3, 5e-3, 8e-3, 11e-3, 14e-3)   # карты все ±50 мм, 14 мм внутри хорошего поля
PLOT_RADIUS = {1: 1e-3, 2: 5e-3, 3: 8e-3, 4: 14e-3}
# Радиусы, НА КОТОРЫХ проверяется устойчивость данного порядка. Проверять октуполь на 1 мм
# бессмысленно: там он заведомо под шумовым полом, и настоящий момент отбраковался бы вместе с
# шумом — у секступоля SP, его ГЛАВНОГО момента, так и выходило 56%.
RES_RADII = {1: (1e-3, 3e-3, 5e-3, 8e-3), 2: (1e-3, 3e-3, 5e-3, 8e-3),
             3: (5e-3, 8e-3, 11e-3, 14e-3), 4: (8e-3, 11e-3, 14e-3)}
RESOLVED = 25.0        # % разброса по радиусу пробы, выше которого момент считается шумом

SPEC_CACHE = os.path.join("workspace", "results", "library_report_spectra_ri.json")
VAR_CACHE = os.path.join("workspace", "results", "library_report_variants.json")
MAPLEN_CACHE = os.path.join("workspace", "results", "library_report_maplen.json")


# ───────────────────────────── моменты из карт ─────────────────────────────

def _load_full(path):
    """Карта из .bin, приведённая к ПОЛНОЙ: библиотека хранит половинки (X≥0)."""
    from field_quality import load
    import _qp_scale as qs
    h, B = load(path)
    if h["mnx"] >= -1e-9:
        h, B = qs.mirror_full(h, B)
    return h, B


def _profile(h, B, order, radius, nth=64):
    """Продольный профиль момента: (x, |b_n(x)|), нормальная компонента."""
    from field_quality import _circle
    xs = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    v = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, radius, nth=nth)
        sp = np.fft.fft(by + 1j * bz) / nth
        v.append((sp[0] if order == 1 else sp[-(order - 1)]).imag)
    return xs, np.abs(np.array(v)) / radius ** (order - 1)


def _kicks(h, B, orders, nth=64):
    """{порядок: (нормальный, скью)} на ПОЛНОЙ карте, каждый порядок на своём радиусе пробы.

    Момент комплексен: (By + i·Bz) даёт нормальный мультиполь и повёрнутый на пол-периода (скью).
    Брать МОДУЛЬ нельзя — он их смешивает: у боковых обмоток корректора скью-часть больше
    нормальной, и модуль завышал дипольный момент на 59%."""
    from field_quality import _circle
    xs = h["mnx"] + np.arange(h["nx"]) * h["sx"]
    out = {}
    for r in sorted({PLOT_RADIUS[n] for n in orders}):
        acc = {n: [] for n in orders if PLOT_RADIUS[n] == r}
        if not acc:
            continue
        for ix in range(h["nx"]):
            by, bz = _circle(h, B, ix, r, nth=nth)
            sp = np.fft.fft(by + 1j * bz) / nth
            for n in acc:
                acc[n].append(sp[0] if n == 1 else sp[-(n - 1)])
        for n, a in acc.items():
            v = np.array(a)
            c = complex(np.trapezoid(v.real, xs), np.trapezoid(v.imag, xs))
            out[n] = (abs(c.imag) / r ** (n - 1), abs(c.real) / r ** (n - 1))
    return out


def spectrum(path, orders, radii, nth=64):
    """Моменты на НЕСКОЛЬКИХ радиусах пробы за одну загрузку карты: {радиус: {порядок: (норм, скью)}}.

    Несколько радиусов нужны не для полноты. Поле мультиполя ∝ r^(n−1), а шум карты аддитивен,
    поэтому у слабого высокого момента сигнал у оси лежит ПОД шумовым полом. Настоящий момент от
    шума отличает то, что после деления на r^(n−1) он НЕ ЗАВИСИТ от радиуса пробы."""
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


def element_spectra(element, orders=(1, 2, 3, 4), radii=PROBE_SCAN):
    """{обмотка: [(ампер-витки, {радиус: {порядок: (норм, скью)}})]} по ОДИНОЧНЫМ картам, с кэшем.

    Карты с несколькими одновременно запитанными обмотками пропускаются: по ним нельзя сказать,
    какая обмотка какой момент создала."""
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
        have = cache.setdefault(r["map_id"], {})
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
    """Разброс момента по радиусу пробы, % от медианы — по верхней половине токов.

    Настоящий момент даёт единицы процентов: поделённый на r^(n−1), он от радиуса не зависит. Шум
    даёт десятки и сотни. Радиусы берутся из RES_RADII: для каждого порядка свои, те, на которых
    момент этого порядка вообще может быть измерен."""
    radii = radii or RES_RADII[order]
    top = pts[len(pts) // 2:]
    sc = []
    for ni, s in top:
        v = np.array([s[r][order][0] for r in radii if r in s])
        m = np.median(v)
        if m > 0:
            sc.append(100 * (v.max() - v.min()) / m)
    return float(np.mean(sc)) if sc else float("nan")


def map_length(element, order):
    """Жёстко-краевая эффективная длина ИЗМЕРЕННОГО магнита: ∫c·dx / max c по полной карте.

    Запасной вариант на случай, когда у элемента нет L_eff в базе. Это НЕ подгоночная L_eff из
    датасетов: подгонка определяет произведение «коэффициент × длина», но не его разложение на
    сомножители."""
    cache = {}
    if os.path.exists(MAPLEN_CACHE):
        try:
            cache = json.load(open(MAPLEN_CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    ck = f"{element}|{order}"
    if ck in cache:
        return cache[ck]
    reg = [x for x in sm.read_registry()
           if x["element"] == element and ";" not in x["excitation"] and "+" not in x["excitation"]]
    if not reg:
        return None
    # Средний ток: у самых слабых карт всё тонет в шуме, у самых сильных профиль искажён насыщением.
    m = sorted(reg, key=lambda x: abs(float(x["amp_turns"] or 0)))[len(reg) // 2]
    h, B = _load_full(os.path.join("workspace", "fields_bin", m["bin_path"]))
    xs, p = _profile(h, B, order, PLOT_RADIUS[order])
    L = float(np.trapezoid(p, xs) / p.max())
    cache[ck] = L
    os.makedirs(os.path.dirname(MAPLEN_CACHE) or ".", exist_ok=True)
    json.dump(cache, open(MAPLEN_CACHE, "w", encoding="utf-8"))
    return L


def stretch(h, B, order, radius, L_target):
    """Карта магнита ДРУГОЙ ДЛИНЫ: тело растягивается, краевое поле сохраняется.

    Штатный гибридный ресайз проекта (`_qp_scale_hybrid.resize`, тот же, что у
    `gen_synthetic_maps`): удлинить магнит — значит добавить прямого участка, а не растянуть его
    целиком, потому что краевое поле определяется поперечником и от длины не зависит.

    Тело выделяется по профилю ВЕДУЩЕГО момента, а не по |B| в стороне от оси, как в
    продакшн-эвристике: у наших карт та срывалась в запасной вариант «тело = вся длина». Разбиение
    считается ДЛЯ КАЖДОЙ КАРТЫ отдельно: в насыщении профиль заостряется."""
    import _qp_scale_hybrid as hy
    xs, p = _profile(h, B, order, radius)
    mx = p.max()
    if mx <= 0:
        return None, None, None
    Xb = float(xs[p >= 0.9 * mx].max())
    L_body, L0 = 2 * Xb, float(np.trapezoid(p, xs) / mx)
    if L_body <= 1e-6:
        return None, None, None
    cb = min(max((L_target - (L0 - L_body)) / L_body, 0.05), 20.0)
    h2, B2 = hy.resize(h, B, cb, 1.0, Xb)
    return h2, B2, cb


def variant_kicks(element, L_target, order_lead, orders=(1, 2, 3, 4)):
    """{обмотка: [(ампер-витки, {порядок: (норм, скью)}, c_body)]} для магнита длины L_target.

    Карты сняты с ОДНОЙ длины семейства, а в библиотеке магнит стоит в нескольких исполнениях.
    Кривые для остальных длин получаются синтезом (`stretch`), иначе кривая намагниченности
    относилась бы к магниту другой длины."""
    cache = {}
    if os.path.exists(VAR_CACHE):
        try:
            cache = json.load(open(VAR_CACHE, encoding="utf-8"))
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
        ck = f"{r['map_id']}|{L_target:.4f}"
        if ck not in cache:
            h, B = _load_full(p)
            h2, B2, cb = stretch(h, B, order_lead, PLOT_RADIUS[order_lead], L_target)
            if h2 is None:
                continue
            k = _kicks(h2, B2, orders)
            cache[ck] = {"cb": cb, "k": {str(n): list(v) for n, v in k.items()}}
            dirty = True
        e = cache[ck]
        data[r["winding"]].append((ni, {int(n): tuple(v) for n, v in e["k"].items()}, e["cb"]))
    if dirty:
        os.makedirs(os.path.dirname(VAR_CACHE) or ".", exist_ok=True)
        json.dump(cache, open(VAR_CACHE, "w", encoding="utf-8"))
    return {w: sorted(v, key=lambda t: t[0]) for w, v in data.items()}


# ───────────────────── привязка к главам исходного отчёта ─────────────────────

def element_of(name):
    """Имя элемента отчёта (`QP_L_1A_L295`) → элемент библиотеки карт (`QP_L`)."""
    cand = [e for e in ELEMENTS if name.startswith(e)]
    return max(cand, key=len) if cand else None


def lead_order(m):
    """Порядок ВЕДУЩЕГО момента элемента — по нему выделяется тело магнита при синтезе длины."""
    return max(MOMENT[mp][0] for mp in m["outputs"]) if m["outputs"] else 2


def chapter_anchor(doc, name):
    """_p следующего заголовка после главы элемента — точка вставки в КОНЕЦ главы.

    Имя главы сравнивается вместе с открывающей скобкой: `DQ_I1_1A_L400` иначе совпало бы и с
    главой `DQ_I1_1A_L400_sext`."""
    ps = doc.paragraphs
    head = name + "  ("
    idx = next((i for i, p in enumerate(ps)
                if p.style.name == "Heading 1" and p.text.strip().startswith(head)), None)
    if idx is None:
        return None
    return next((ps[j]._p for j in range(idx + 1, len(ps))
                 if ps[j].style.name.startswith("Heading")), None)


def saturation_currents(m, ec, pts_by_w, L):
    """Проектные токи С УЧЁТОМ НАСЫЩЕНИЯ: (λ, {обмотка: ток, А}, управляемые моменты).

    Матрица связи линейна — она снята свипом в рабочей области и продолжена прямой. Кривые
    намагниченности содержат насыщение железа, поэтому линейный вектор токов на проектные моменты
    вообще говоря не выводит. Ищется ОДИН множитель λ ко всему вектору: направление задаёт матрица
    (какая обмотка какую долю несёт), а величину — измеренные кривые. λ = 1 значит, что проектная
    точка лежит в линейной области.

    Моменты обмоток складываются: постулат суперпозиции проверен на секступоле (16 комбинаций
    ANSYS против суммы, 1.5%). Знак вклада берётся у линейного предсказания матрицы — сами кривые
    приведены по модулю."""
    from scipy.optimize import brentq
    Ivec, ctrl, design = ec._design_ivec(m)
    if not ctrl:
        return None, {}, ctrl
    tpc = [ec.dev_tpc(d) or 1 for d in m["devs"]]
    Ilin = np.array([Ivec[j] * ec.ENERGY_FACTOR / tpc[j] for j in range(len(tpc))])
    M = np.asarray(m["M"], float)
    curves = {}
    for j, w in enumerate(m["labels"]):
        p = pts_by_w.get(w) or pts_by_w.get("I1" if w == "I" else "I")
        if p:
            x = np.array([q[0] for q in p]) / tpc[j]
            curves[j] = (x, {n: np.array([q[1][n][0] for q in p]) for n in (1, 2, 3, 4)},
                         {n: np.array([q[1][n][1] for q in p]) for n in (1, 2, 3, 4)})
    if not curves:
        return None, {}, ctrl

    def target(mp):                                   # проектный момент в физических единицах
        return abs(design[mp]) * L * ec.ENERGY_FACTOR / math.factorial(MOMENT[mp][0] - 1)

    def resid(lam):
        tot = 0.0
        for mp in ctrl:
            n, comp = MOMENT[mp]
            k = m["outputs"].index(mp)
            s = 0.0
            for j, (x, nrm, skw) in curves.items():
                y = (nrm if comp == 0 else skw)[n]
                sgn = np.sign(M[k, j] * Ivec[j]) or 1.0
                s += sgn * np.interp(abs(lam * Ilin[j]), x, y, right=y[-1])
            tot += (abs(s) - target(mp)) / target(mp)
        return tot

    try:
        lam = float(brentq(resid, 0.05, 20.0, xtol=1e-6))
    except (ValueError, RuntimeError):
        return None, {}, ctrl
    return lam, {w: lam * Ilin[j] for j, w in enumerate(m["labels"])}, ctrl


# ───────────────────────── кривые намагниченности ─────────────────────────

def to_coef(value, moment, L, ef):
    """Интеграл момента по длине (физические единицы) → физкоэффициент в шкале отчёта (mat_*).

    Обратное преобразование к тому, которым проектный mat_* переводится в измеримую величину:
    делим на длину, на энергетический множитель и умножаем на (n−1)! — таково соглашение столбцов
    mat_*. После него ось графика и столбец таблицы — одно и то же число."""
    return value * math.factorial(MOMENT[moment][0] - 1) / (L * ef)


def plot_winding(m, label, pts, turns, L, moments, design, cur, ef, outdir):
    """Кривая намагниченности обмотки: физкоэффициенты в зависимости от тока питания.

    Проектная точка обозначена линиями: горизонталь — проектный физкоэффициент из .lte, вертикаль —
    проектный ток обмотки с учётом насыщения, ромб — их пересечение."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    x = np.array([p[0] for p in pts]) / (turns or 1)
    cb = float(np.mean([p[2] for p in pts]))
    fig, ax = plt.subplots(figsize=(9.0, 5.2))
    axes = {mp: (ax if k == 0 else ax.twinx()) for k, mp in enumerate(moments)}
    for k, mp in enumerate(moments):
        n, comp = MOMENT[mp]
        a, col = axes[mp], MOMENT_COLOR[mp]
        y = np.array([to_coef(p[1][n][comp], mp, L, ef) for p in pts])
        a.plot(x, y, "-o", color=col, ms=4.5, lw=1.6, zorder=3)
        a.set_ylabel(mp, color=col, fontsize=10)
        a.tick_params(axis="y", labelcolor=col, labelsize=8)
        if k > 1:
            a.spines["right"].set_position(("outward", 46 * (k - 1)))
        ymax = y.max()
        d = abs(design.get(mp, 0.0))
        if d > 1e-30:
            a.axhline(d, color=col, lw=1.2, ls="-", alpha=0.85, zorder=1)
            a.annotate(f"{mp} = {d:.4e}", xy=(0.015, 0.965 - 0.055 * k), xycoords="axes fraction",
                       fontsize=8.5, color=col, family="monospace", va="top",
                       bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
            if cur is not None and abs(cur) > 1e-9:
                a.plot([abs(cur)], [d], "D", color=col, ms=8, mec="white", mew=1.0, zorder=5)
            ymax = max(ymax, d)
        a.set_ylim(0, ymax * 1.15)
    if cur is not None and abs(cur) > 1e-9:
        ax.axvline(abs(cur), color="#444", lw=1.3, ls="--", zorder=2)
        ax.annotate(f"{abs(cur):.1f} А", xy=(abs(cur), 0.02), xycoords=("data", "axes fraction"),
                    fontsize=9, color="#222", ha="left", va="bottom",
                    bbox=dict(fc="white", ec="none", alpha=0.85, pad=1.5))
    ax.set_xlabel(f"ток обмотки, А   ({turns} витк./катушку)")
    ax.set_title(f"{m['name']} — обмотка {label}: " + ", ".join(moments)
                 + f";  L_eff = {L:.3f} м"
                 + ("" if abs(cb - 1) < 0.03 else f" (тело карты ×{cb:.2f})"), fontsize=10.5)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{m['name']}_{label}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# ───────────────────────── дополнение отчёта ─────────────────────────

def build(out=OUT, src=SRC, figs=PNG):
    """Копия исходного отчёта, дополненная кривыми намагниченности и оценкой проектных токов с
    учётом насыщения. Документ не пересобирается: открывается как есть, и в конец каждой главы
    элемента вставляются только эти два дополнения."""
    from docx import Document
    from docx.shared import Inches
    import design_currents as dc
    import element_compare as ec

    if not os.path.exists(src):
        raise SystemExit(f"нет исходного отчёта: {src}\n"
                         f"построить: py -3.14 element_compare.py library")
    doc = Document(src)
    mats = dc.load_coupling()
    dc.attach_devices(mats, ec.LCROW)
    for m in mats:
        if m.get("labels"):
            m["labels"] = [dc._norm_label(l) for l in m["labels"]]

    added = 0
    for m in mats:
        el = None if m.get("pm") else element_of(m["name"])
        if not el:
            continue
        anchor = chapter_anchor(doc, m["name"])
        if anchor is None:
            print(f"  главы {m['name']} в отчёте нет — пропуск")
            continue
        lead = lead_order(m)
        L = float(m["L"]) if m["L"] else map_length(el, lead)
        pts_by_w = variant_kicks(el, L, lead)
        native = element_spectra(el)
        allow = KIND_MOMENTS.get(m["kind"], tuple(MOMENT))
        lam, cur, ctrl = saturation_currents(m, ec, pts_by_w, L)
        design = ec._design_ivec(m)[2]

        new = []
        p = doc.add_paragraph(f"Кривые намагниченности {m['name']}: киковые моменты обмоток в "
                              f"зависимости от тока питания, по трёхмерным картам поля.")
        new.append(p._p)
        for j, lb in enumerate(m["labels"]):
            pts = pts_by_w.get(lb) or pts_by_w.get("I1" if lb == "I" else "I")
            nat = native.get(lb) or native.get("I1" if lb == "I" else "I")
            if not pts or not nat:
                continue
            shown = [mp for mp in MOMENT
                     if mp in allow and resolution(nat, MOMENT[mp][0]) <= RESOLVED
                     and (MOMENT[mp][1] == 0 or mp in m["outputs"])]
            if not shown:
                continue
            png = plot_winding(m, lb, pts, ec.dev_tpc(m["devs"][j]), L, shown,
                               design, cur.get(lb), ec.ENERGY_FACTOR, figs)
            pic = doc.add_paragraph()
            pic.add_run().add_picture(png, width=Inches(6.3))
            new.append(pic._p)

        cap = doc.add_paragraph()
        cap.add_run(f"{m['name']} — проектные токи с учётом насыщения").bold = True
        new.append(cap._p)
        t = doc.add_table(rows=1, cols=5)
        t.style = "Light Grid Accent 1"
        for k, h in enumerate(("элемент", "обмотка", "проектные физкоэффициенты (.lte)",
                               "проектный ток, А", "% номинала")):
            t.rows[0].cells[k].text = h
        coefs = ";  ".join(f"{mp} = {design[mp]:+.4e}" for mp in ctrl) if ctrl else "—"
        for j, lb in enumerate(m["labels"]):
            c = t.add_row().cells
            c[0].text, c[1].text, c[2].text = m["name"], lb, coefs
            I = cur.get(lb)
            if I is None:
                c[3].text = c[4].text = "—"
                continue
            nom = ec.dev_nominal(m["devs"][j])
            c[3].text = f"{I:+.1f}"
            c[4].text = f"{100 * abs(I) / nom:.0f}%" if nom else "—"
        new.append(t._tbl)
        new.append(doc.add_paragraph()._p)

        for node in new:                                   # перенос созданных блоков в конец главы
            anchor.addprevious(node)
        added += 1

    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    try:
        doc.save(out)
    except PermissionError:
        out = os.path.splitext(out)[0] + "_new.docx"
        doc.save(out)
        print(f"  (отчёт заблокирован Word — записан {os.path.basename(out)})")
    print(f"дополнено глав: {added}")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--src", default=SRC, help="исходный element_compare_report.docx")
    a = ap.parse_args(argv)
    print("отчёт ->", build(a.out, a.src))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
