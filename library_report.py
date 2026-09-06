# -*- coding: utf-8 -*-
"""Word-отчёт: киковые моменты B/G/S/O от тока обмоток по элементам библиотеки.

Устройство отчёта повторяет `element_compare_report`: вводный методический раздел, однородный раздел
на КАЖДЫЙ элемент библиотеки (то есть на исполнение — `QP_L_1A_L295`, `QP_L_1A_L374`, … отдельными
разделами), расшифровка моментов, проектные токи обмоток и сверка проекта с фактом. Терминология и
таблицы берутся штатными функциями `element_compare`, поэтому совпадают дословно.

Отличие — содержание разделов: вместо таблиц свипа и траекторий приведены ГРАФИКИ зависимости
киковых моментов от тока обмотки, снятые по картам поля, и на них отмечен проектный ток обмотки.

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

OUT = os.path.join("workspace", "results", "library_report.docx")
PNG = os.path.join("workspace", "results", "library_report_figs")

# Киковый момент → (порядок мультиполя, компонента): 0 — нормальная, 1 — скью. Обозначения и порядок
# те же, что в `element_compare.MULTIPOLES`.
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
    """Жёстко-краевая эффективная длина ИЗМЕРЕННОГО магнита: 2·∫c·dx / max c (карта половинная).

    Это НЕ подгоночная L_eff из датасетов: подгонка определяет только ПРОИЗВЕДЕНИЕ «коэффициент ×
    длина», а как разделить его на сомножители — не определяет (при смене сетки подгоночная длина
    QP_sh уехала на 13% при неизменном произведении). Длина, взятая из ФОРМЫ поля, попадает в
    проектную: QP_L 0.374 против 0.3745, DQ 0.403 против 0.400."""
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
    целиком, потому что краевое поле определяется поперечником и от длины не зависит. Отсюда
    тождество L_eff = c_body·L_тела + L_краевого.

    Тело выделяется по профилю ВЕДУЩЕГО момента (его интеграл и есть кик), а не по |B| в стороне от
    оси, как в продакшн-эвристике: у наших карт та срывалась в запасной вариант «тело = вся длина».
    Разбиение считается ДЛЯ КАЖДОЙ КАРТЫ отдельно: в насыщении профиль заостряется."""
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

    Карты сняты с ОДНОЙ длины из семейства, а в библиотеке магнит стоит в нескольких исполнениях.
    Кривые для остальных длин получаются синтезом (`stretch`), иначе проектный момент выглядел бы
    недостижимым просто потому, что сравнивался с более коротким магнитом."""
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


# ───────────────────────── привязка к элементам отчёта ─────────────────────────

def element_of(name):
    """Имя элемента отчёта (`QP_L_1A_L295`) → элемент библиотеки карт (`QP_L`)."""
    cand = [e for e in ELEMENTS if name.startswith(e)]
    return max(cand, key=len) if cand else None


def lead_order(m):
    """Порядок ВЕДУЩЕГО момента элемента — по нему выделяется тело магнита при синтезе длины."""
    return max(MOMENT[mp][0] for mp in m["outputs"]) if m["outputs"] else 2


def design_currents(m, ec):
    """{обмотка: (проектный ток, А; % номинала)} — те же числа, что в таблице «проектные токи
    обмоток». Считаются обращением матрицы связи под проектные mat_* решётки (`_design_ivec`),
    поэтому у многообмоточных элементов это КОМБИНАЦИЯ, а не ток одной обмотки."""
    Ivec, ctrl, design = ec._design_ivec(m)
    out = {}
    if not ctrl:
        return out, design
    for j, lb in enumerate(m["labels"]):
        dev = m["devs"][j]
        tpc = ec.dev_tpc(dev)
        if not tpc:
            continue
        nom = ec.dev_nominal(dev)
        cur = Ivec[j] * ec.ENERGY_FACTOR / tpc
        out[lb] = (cur, (100.0 * abs(cur) / nom) if nom else None)
    return out, design


def curve_current(m, el, L, lead, design, pts_by_w, ec):
    """{момент: (ток по матрице, ток по кривой карт)} — для ОДНООБМОТОЧНЫХ элементов.

    Матрица связи линейна: она снята свипом в рабочей области и продолжена прямой. Кривая, снятая
    по картам, содержит насыщение железа. Там, где проектная точка лежит на изгибе кривой, эти два
    тока расходятся, и расхождение — прямая мера того, насколько линейная матрица недооценивает
    нужный ток. У многообмоточных элементов вопрос не ставится: проектную точку они берут
    комбинацией, и по кривой одной обмотки её не прочитать."""
    if len(m["labels"]) != 1:
        return {}
    w = m["labels"][0]
    pts = pts_by_w.get(w) or pts_by_w.get("I1" if w == "I" else "I")
    tpc = ec.dev_tpc(m["devs"][0])
    if not pts or not tpc:
        return {}
    x = np.array([q[0] for q in pts]) / tpc
    out = {}
    for mp in m["outputs"]:
        d = abs(design.get(mp, 0.0))
        if d < 1e-30:
            continue
        n, comp = MOMENT[mp]
        y = np.array([q[1][n][comp] for q in pts])
        tgt = d * L * ec.ENERGY_FACTOR / math.factorial(n - 1)
        xi = None
        for i in range(1, len(y)):
            if (y[i - 1] - tgt) * (y[i] - tgt) <= 0 and y[i] != y[i - 1]:
                xi = float(x[i - 1] + (tgt - y[i - 1]) / (y[i] - y[i - 1]) * (x[i] - x[i - 1]))
                break
        out[mp] = xi
    return out


# ───────────────────────────── графики ─────────────────────────────

def plot_winding(m, label, pts, turns, L, moments, design, outdir):
    """График одной обмотки: киковые моменты от тока, вертикаль — проектный ток обмотки.

    Показываются только моменты, РАЗРЕШИМЫЕ в снятых картах (`resolution`): нарисовать шум под
    именем октуполя хуже, чем не рисовать ничего, — по такой кривой будут считать ток."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    x = np.array([p[0] for p in pts]) / (turns or 1)
    cb = float(np.mean([p[2] for p in pts]))
    fig, ax = plt.subplots(figsize=(9.0, 5.4))
    axes = {mp: (ax if k == 0 else ax.twinx()) for k, mp in enumerate(moments)}
    for k, mp in enumerate(moments):
        n, comp = MOMENT[mp]
        a, col = axes[mp], MOMENT_COLOR[mp]
        y = np.array([p[1][n][comp] for p in pts])
        a.plot(x, y, "-o", color=col, ms=4.5, lw=1.6, zorder=3)
        a.set_ylabel(f"{mp}, {MOMENT_UNITS[mp]}   (проба {PLOT_RADIUS[n] * 1e3:.0f} мм)",
                     color=col, fontsize=9)
        a.tick_params(axis="y", labelcolor=col, labelsize=8)
        if k > 1:
            a.spines["right"].set_position(("outward", 58 * (k - 1)))
        a.set_ylim(0, y.max() * 1.15)

    cur = design.get(label)
    if cur and abs(cur[0]) > 1e-9:
        ax.axvline(abs(cur[0]), color="#444", lw=1.4, ls="--", zorder=2)
        ax.annotate(f"проектный ток {abs(cur[0]):.1f} А"
                    + (f"  ({cur[1]:.0f}% номинала)" if cur[1] is not None else ""),
                    xy=(0.015, 0.965), xycoords="axes fraction", fontsize=9, va="top",
                    color="#222", bbox=dict(fc="white", ec="#bbb", alpha=0.9, pad=2.5))
    sub = [f"L_eff = {L:.3f} м" + ("" if abs(cb - 1) < 0.03 else
                                   f" (карта растянута ×{cb:.2f} по телу)")]
    ax.set_xlabel(f"ток обмотки, А   ({turns} витк./катушку)")
    ax.set_title(f"{m['name']} — обмотка {label}: моменты " + ", ".join(moments) + " от тока\n"
                 + ";  ".join(sub), fontsize=11)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{m['name']}_{label}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


# ───────────────────────────── документ ─────────────────────────────

def _intro(doc):
    """Вводный методический раздел — в терминах и стиле `element_compare_report`."""
    doc.add_heading("Вводная информация по расчетам", level=1)
    doc.add_paragraph(
        "Отчёт показывает, как киковые моменты элементов библиотеки зависят от токов их обмоток. "
        "Моменты получены гармоническим разложением трёхмерных карт магнитного поля, снятых для "
        "каждой обмотки: поле раскладывается на окружности вокруг оси пучка, а гармоника порядка n "
        "интегрируется вдоль элемента. Отложен именно ИНТЕГРАЛ момента по длине — его и чувствует "
        "пучок; сам коэффициент получается делением интеграла на эффективную длину, и это "
        "соглашение, а не измерение.")
    doc.add_paragraph(
        "Разложение сверено с независимым способом — подбором моментов по траекториям, тем самым, "
        "которым получены матрицы связи: дифференциальное ядро SynchrotronSim интегрирует силу "
        "Лоренца по карте, а моменты кика подбираются так, чтобы киковая траектория её "
        "воспроизвела. На самых сильных картах шести элементов оба способа сходятся в пределах "
        "0.72%.")
    doc.add_paragraph(
        "Разделы отчёта соответствуют элементам библиотеки, то есть исполнениям: длинный квадруполь "
        "стоит в решётке в трёх длинах, и каждой отведён свой раздел. Карты полей сняты с одной "
        "длины семейства, поэтому кривые остальных исполнений получены штатным синтезом карт "
        "(`gen_synthetic_maps`): тело магнита растягивается или сжимается, а краевое поле "
        "сохраняется, потому что определяется поперечником, а не длиной. Коэффициент растяжения "
        "тела указан в заголовке графика; синтез проверен на попадание в целевую L_eff — "
        "расхождение 0.2–1.2%.")
    doc.add_paragraph(
        "Токи и коэффициенты. Все токи в отчёте — физические токи питания обмоток (амперы при "
        "проектной энергии 6 ГэВ). «Проектные токи обмоток» (таблица в конце каждого раздела) — это "
        "токи питания, при которых элемент даёт ПРОЕКТНЫЕ мультипольные коэффициенты, заданные в "
        "файле решётки .lte; они получены обращением матрицы связи под проектные коэффициенты. "
        "Рядом приведены «ампер-витки на катушку» — возбуждение магнита, равное току обмотки, "
        "умноженному на число витков в катушке. Именно эта величина задаётся при расчёте поля и "
        "определяет состояние железа, поэтому она приведена явно: ток обмотки сам по себе о "
        "насыщении ничего не говорит. У многообмоточных элементов проектный ток — это КОМБИНАЦИЯ: "
        "элемент выходит на проектные моменты всеми обмотками сразу, а не каждой по отдельности.")
    doc.add_paragraph(
        "На графиках показаны только моменты, РАЗРЕШИМЫЕ в снятых картах. Поле мультиполя растёт "
        "как r^(n−1), а шум карты от радиуса не зависит, поэтому слабый высокий момент у оси тонет. "
        "Настоящий момент, поделённый на r^(n−1), от радиуса пробы не зависит; по этому признаку "
        f"моменты и отбираются (допуск {RESOLVED:.0f}% разброса). Неразрешимый момент не значит "
        "нулевой — значит, что по этим картам его величину назвать нельзя.")


def build(out=OUT, figs=PNG):
    from docx import Document
    from docx.shared import Inches
    import design_currents as dc
    import element_compare as ec

    mats = dc.load_coupling()
    dc.attach_devices(mats, ec.LCROW)
    for m in mats:                       # `I_c1` → `c1`: в таблицах исходного отчёта обмотки без префикса
        if m.get("labels"):
            m["labels"] = [dc._norm_label(l) for l in m["labels"]]
    ec._tcap.n = 0                                   # сквозная нумерация таблиц — с начала

    doc = Document()
    doc.add_heading("SynchrotronSim — киковые моменты от тока обмоток: расчёт по элементам", level=0)
    _intro(doc)

    stale = dc.staleness()
    if stale:
        doc.add_paragraph("Внимание: " + stale)

    for m in mats:
        doc.add_heading(f"{m['name']}  ({m['kind']}, "
                        + ("PM, фикс. поле" if m.get("pm") else f"обмоток: {len(m['labels'])}")
                        + ")", level=1)
        if m["kind"] == "dq":
            doc.add_paragraph("Ядро кика: диполь + квадруполь + секступоль (моменты "
                              + ", ".join(m["outputs"]) + ")." if "S" in m["outputs"] else
                              "Ядро кика: диполь + квадруполь (моменты "
                              + ", ".join(m["outputs"]) + ").")
        ec._moment_legend(doc, m["outputs"])

        if m.get("pm"):
            ec._tcap(doc, m["name"])
            t = doc.add_table(rows=1, cols=len(m["outputs"])); t.style = "Light Grid Accent 1"
            for j, mp in enumerate(m["outputs"]):
                t.rows[0].cells[j].text = f"{mp}, {MOMENT_UNITS[mp]}"
            c = t.add_row().cells
            for j, mp in enumerate(m["outputs"]):
                c[j].text = f"{m['fixed'][mp]:+.4e}"
            doc.add_paragraph()
            doc.add_paragraph("Постоянный магнит — обмоток и токов нет, карт поля по обмоткам не "
                              "снимают; состав поля приведён в таблице выше.")
            continue

        el = element_of(m["name"])
        cur_by_w, design = design_currents(m, ec)
        by_curve = {}
        if el:
            lead = lead_order(m)
            L = float(m["L"]) if m["L"] else map_length(el, lead)
            pts_by_w = variant_kicks(el, L, lead)
            by_curve = curve_current(m, el, L, lead, design, pts_by_w, ec)
            native = element_spectra(el)
            for j, lb in enumerate(m["labels"]):
                w = lb
                pts = pts_by_w.get(w) or pts_by_w.get("I1" if w == "I" else "I")
                nat = native.get(w) or native.get("I1" if w == "I" else "I")
                if not pts or not nat:
                    continue
                allow = KIND_MOMENTS.get(m["kind"], tuple(MOMENT))
                shown = [mp for mp in MOMENT
                         if mp in allow and resolution(nat, MOMENT[mp][0]) <= RESOLVED
                         and (MOMENT[mp][1] == 0 or mp in m["outputs"])]
                if not shown:
                    continue
                png = plot_winding(m, lb, pts, ec.dev_tpc(m["devs"][j]), L, shown,
                                   cur_by_w, figs)
                doc.add_picture(png, width=Inches(6.3))
            Lm = map_length(el, lead)
            if Lm and abs(Lm - L) / L > 0.10:
                doc.add_paragraph(
                    f"Замечание по длине: L_eff элемента в базе {L:.3f} м, а продольный профиль "
                    f"снятой карты даёт {Lm:.3f} м (интеграл профиля, делённый на максимум). Кривые "
                    f"выше приведены к базовой L_eff синтезом, но само расхождение стоит проверить.")
            doc.add_paragraph()

        ec._working_point_table(doc, m)

        if cur_by_w:
            over = [f"{lb} — {abs(v[0]):.1f} А ({v[1]:.0f}% номинала)"
                    for lb, v in cur_by_w.items() if v[1] and v[1] > 100]
            p = doc.add_paragraph()
            p.add_run("Проектные токи. ").bold = True
            p.add_run("Ток — свойство элемента, а не отдельной обмотки: приведённые значения "
                      "получены обращением матрицы связи под проектные коэффициенты решётки, то "
                      "есть это набор токов, который магнит держит одновременно. "
                      + (("За номинал выходят: " + "; ".join(over) + ".") if over else
                         "Все обмотки укладываются в номинальный ток."))
            if by_curve:
                lb = m["labels"][0]
                Im = abs(cur_by_w.get(lb, (0.0, None))[0])
                parts = []
                for mp, xi in by_curve.items():
                    if xi is None:
                        parts.append(f"{mp} — проектное значение кривой не достигается")
                    elif Im > 0:
                        parts.append(f"{mp} — {xi:.1f} А ({100 * (xi / Im - 1):+.0f}% к матричному)")
                q = doc.add_paragraph()
                q.add_run("Матрица и кривая. ").bold = True
                q.add_run(
                    "Матрица связи линейна: она снята свипом в рабочей области и продолжена "
                    "прямой. Кривая, снятая по картам, содержит насыщение железа, поэтому ток, при "
                    "котором элемент реально выходит на проектный момент, может отличаться: "
                    + "; ".join(parts) + ". Расхождение — прямая мера того, насколько проектная "
                    "точка ушла на изгиб кривой намагничивания.")
        else:
            p = doc.add_paragraph()
            p.add_run("Проектные токи. ").bold = True
            p.add_run("Проектных коэффициентов у этого элемента в решётке нет — он не оттюнен, и "
                      "рабочая точка ему ещё не назначена. Графики выше показывают саму "
                      "зависимость моментов от тока, по которой ток и будет определяться, когда "
                      "рабочая точка появится.")

    ec._apply_body_formatting(doc)
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    try:
        doc.save(out)
    except PermissionError:
        out = os.path.splitext(out)[0] + "_new.docx"
        doc.save(out)
        print(f"  (отчёт заблокирован Word — записан {os.path.basename(out)})")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    print("отчёт ->", build(a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
