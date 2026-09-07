# -*- coding: utf-8 -*-
"""Достоверная копия `element_compare_report`, дополненная двумя вещами и только ими.

Исходный отчёт НЕ пересобирается: он открывается как есть, и в конец каждой главы элемента
добавляются кривые намагниченности обмоток (киковые моменты в зависимости от тока питания, по
трёхмерным картам поля) и таблица проектных токов С УЧЁТОМ НАСЫЩЕНИЯ, а в самый конец — раздел с
кривой намагничивания материала сердечника. Всё остальное — вводный
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

SMODELS = os.path.join("workspace", "tables", "saturation_models.csv")
LEAD_MOMENT = {"quad": "G", "sext": "S", "oct": "O"}   # ведущий момент однообмоточного элемента
SINGLE_KINDS = ("quad", "sext", "oct")
SPEC_CACHE = os.path.join("workspace", "results", "library_report_spectra_ri.json")
VAR_CACHE = os.path.join("workspace", "results", "library_report_variants.json")
MAPLEN_CACHE = os.path.join("workspace", "results", "library_report_maplen.json")
MATERIAL = os.path.join("workspace", "materials", "steel_1008_bh.csv")


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

def _smodels():
    """Параметры ядра моделей насыщения из БД: {элемент: строка saturation_models.csv}."""
    import csv
    if not os.path.exists(SMODELS):
        return {}
    return {r["element"]: r for r in csv.DictReader(open(SMODELS, encoding="utf-8-sig"),
                                                    delimiter=";")}


def _g_design(design, mp, ef):
    """Проектный физкоэффициент (Тл/м …) из проектного mat_* отчёта: g = mat·EF/(n−1)!."""
    n = MOMENT[mp][0]
    return design[mp] * ef / math.factorial(n - 1)


def sat_currents_solve(el, kind, m, design, ec, pts_by_w=None, L=None):
    """Проектные токи С УЧЁТОМ НАСЫЩЕНИЯ. Возврат {обмотка: ток, А}.

    Однообмоточные считаются граф-моделью с утечкой, подогнанной на картах ЭТОГО исполнения
    (`single_current`) — той же, что рисует линию на графике. Многообмоточные — графом магнитной
    цепи через `saturation_solve`. Токи обмоток, которые модель не описывает (тримы), берутся
    линейными: на малом токе насыщения нет, и разницы с моделью там тоже нет."""
    import saturation_solve as ss
    ef = ec.ENERGY_FACTOR
    tpc = [ec.dev_tpc(d) or 1 for d in m["devs"]]
    lin = {w: float(np.asarray(ec._design_ivec(m)[0])[j]) * ef / tpc[j]
           for j, w in enumerate(m["labels"])}
    # Знак цели берётся из модели на ЛИНЕЙНЫХ токах: проектный mat_* хранится по модулю в своей
    # конвенции, а физический знак (например, дефокусирующий квадруполь G < 0) даёт сама модель.
    def signed(el_, cur_lin, mps):
        fwd = ss.forward(el_, cur_lin)
        out = {}
        for mp in mps:
            s = 1.0 if fwd.get(mp, (0.0,))[0] >= 0 else -1.0
            out[mp] = math.copysign(abs(_g_design(design, mp, ef)), s)
        return out

    try:
        if kind == "dq":
            mps = [mp for mp in ("B", "G") if mp in design and abs(design[mp]) > 1e-30]
            tgt = signed("DQ", {w: lin[w] for w in lin}, mps)
            cur = ss.inverse("DQ", tgt, fixed={"I3": 0.0})
            return {w: cur.get(w, lin.get(w, 0.0)) for w in m["labels"]}
        lead = LEAD_MOMENT[kind]
        if lead not in design or abs(design[lead]) < 1e-30:
            return lin
        w0 = m["labels"][0]
        n = MOMENT[lead][0]
        pts = (pts_by_w or {}).get(w0) or (pts_by_w or {}).get("I1") or (pts_by_w or {}).get("I")
        if not pts or not L:
            return lin
        y_design = _g_design(design, lead, ef) * L          # проектный интеграл момента
        out = dict(lin)
        out[w0] = single_current(el, L, pts, n, y_design,
                                 ec.dev_tpc(m["devs"][0]) or 1,
                                 lin[w0] if lin[w0] else 1.0, w0)
        return out
    except Exception:
        return lin


def formula_paragraph(doc, el, kind, par, L, lo):
    """Модель насыщения элемента с подобранными числами — граф магнитной цепи с утечкой.

    Та же модель, что рисует линию на графике и даёт проектный ток: последовательные ветви (ярмо и
    участки полюс–ярмо) плюс утечка, шунтирующая катушку. Параметры подобраны на картах поля этого
    исполнения."""
    import single_winding_model as sw
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    lead = LEAD_MOMENT[kind]
    v = np.exp(np.asarray(par, float))
    lp, Sp, Rg, P = v[0], v[1], v[2], v[3]
    terms = []
    for k in range(_ns(par)):
        mult = sw.SERIES_MULT[k] if k < len(sw.SERIES_MULT) else 1.0
        pre = "" if mult == 1.0 else "%g " % mult
        terms.append(r"%s%.4g H\left(\frac{y}{%.4g}\right)" % (pre, v[4 + 2 * k], v[5 + 2 * k]))

    doc.add_paragraph(
        f"Модель насыщения: ампер-витки обмотки связаны с интегралом момента y = ∫{lead}·dl "
        f"уравнениями магнитной цепи (ярмо, участки полюс–ярмо, зазор и утечка, шунтирующая "
        f"катушку); подобранные параметры подставлены.")
    for tex in (r"Q(y) = %.4g y + %s" % (Rg, " + ".join(terms)),
                r"p(y) = y + %.4g Q(y)" % (0.5 * P),
                r"NI = %.4g H\left(\frac{p(y)}{%.4g}\right) + \frac{Q(y)}{2}" % (lp, Sp)):
        q = doc.add_paragraph()
        q.alignment = WD_ALIGN_PARAGRAPH.CENTER
        lo.add_formula(q, tex)
    doc.add_paragraph(
        f"здесь y = ∫{lead}·dl — интеграл момента вдоль элемента (физкоэффициент "
        f"{lead} = y/lэф, lэф = {L:.4g} м), Q(y) — сумма магнитных напряжений "
        f"на последовательных ветвях и зазоре, p(y) — поток полюса (полезный плюс утечка), H(·) — "
        f"обратная кривая намагничивания стали, NI = N·I — ампер-витки. Проектный ток получается "
        f"обращением этой цепочки под проектный {lead} и делением на число витков N.")


def to_coef(value, moment, L, ef):
    """Интеграл момента по длине (физические единицы) → физкоэффициент в шкале отчёта (mat_*).

    Обратное преобразование к тому, которым проектный mat_* переводится в измеримую величину:
    делим на длину, на энергетический множитель и умножаем на (n−1)! — таково соглашение столбцов
    mat_*. После него ось графика и столбец таблицы — одно и то же число."""
    return value * math.factorial(MOMENT[moment][0] - 1) / (L * ef)


_HOF = None
_GRAPH = {}


def _hof():
    import single_winding_model as sw
    global _HOF
    if _HOF is None:
        _HOF = sw.material(os.path.join("workspace", "materials", "steel_1008_bh.csv"))
    return _HOF


def _ns(par):
    """Сколько последовательных классов ветвей в подобранном графе: par = [ℓп, Sп, Rз, P, (ℓ, S)…].

    Подгонка вложенная по числу участков и может остановиться на меньшем, если оно не хуже, —
    поэтому число берётся из самого вектора, а не задаётся снаружи."""
    return (len(par) - 4) // 2


def single_graph_fit(el, L, pts, order, tag=""):
    """Граф-модель с утечкой (single_winding_model) на картах ЭТОГО исполнения: (par, y_max).

    Именно эта модель строилась для однообмоточных магнитов: сердечник + катушка + утечка + участки
    полюс–ярмо. Подгоняется на синтетических картах конкретной длины, поэтому масштаб по длине не
    экстраполируется, а берётся из самих карт исполнения."""
    import single_winding_model as sw
    key = (el, round(float(L), 4), order, tag)    # обмотка в ключе: у SP их две, цепи разные
    if key in _GRAPH:
        return _GRAPH[key]
    ni = np.array([p[0] for p in pts], float)
    y = np.array([abs(p[1][order][0]) for p in pts], float)
    o = np.argsort(ni)
    ni, y = ni[o], y[o]
    good = ni > 0
    par, _ = sw.best_fit(ni[good], y[good], _hof(), nseries=2, tries=3, seeds=40)
    _GRAPH[key] = (par, float(y.max()))
    return _GRAPH[key]


def single_current(el, L, pts, order, y_design, tpc, sign, tag=""):
    """Проектный ток однообмоточного магнита: обращение граф-модели NI(y) при проектном кике y."""
    import single_winding_model as sw
    par, _ = single_graph_fit(el, L, pts, order, tag)
    ni = float(sw.graph_mmf(par, np.array([abs(y_design)]), _hof(), _ns(par))[0])
    return math.copysign(ni / (tpc or 1), sign)


def single_model_curve(el, L, pts, order, tpc, ef, tag=""):
    """Кривая граф-модели для показа: (ток_А, значение_в_шкале_графика) свипом по кику y."""
    import single_winding_model as sw
    par, ymax = single_graph_fit(el, L, pts, order, tag)
    yg = np.linspace(0.0, ymax * 1.05, 400)
    ni = sw.graph_mmf(par, yg, _hof(), _ns(par))
    x = ni / (tpc or 1)
    yv = yg * math.factorial(order - 1) / (L * ef)          # то же, что to_coef
    # На полке ток растёт быстрее любой полезной шкалы: за пределами снятых карт линия только
    # растянула бы ось, поэтому обрезаем её по измеренному диапазону тока.
    xmax = max(p[0] for p in pts) / (tpc or 1) * 1.05
    k = x <= xmax
    return x[k], yv[k]


def moment_possible(el, kind, winding, comp):
    """Разрешена ли компонента симметрией обмотки (у многообмоточных — по маске модели связи).

    У обмотки `main` корректора скошенных компонент быть не может: в картах на их месте остаётся
    численный шум порядка 1e-7, который на собственной оси выглядит как сигнал."""
    if comp == 0 or kind not in ("dq", "kicker"):
        return True
    try:
        import coupling_model as cm
        return cm.allowed(el, winding, comp)
    except Exception:
        return True


def model_curves(el, kind, label, moments, pts, turns, ef, L):
    """Кривые МОДЕЛИ насыщения для показываемых моментов: {момент: (ток_А, значение_в_шкале_графика)}.

    Модель считает физкоэффициент g (Тл/м …); шкала графика — та же, что у измерений (`to_coef`),
    в ней значение равно g·(n−1)!/EF (длина сокращается). Обмотку однообмоточного элемента модель
    знает как I1, у многообмоточных метка совпадает."""
    import saturation_solve as ss
    wkey = "I1" if kind in SINGLE_KINDS else label
    have = {"quad": ("G",), "sext": ("S",), "oct": ("O",),
            "dq": ("B", "G", "S"), "kicker": ("B", "Bv", "G", "Gv", "S")}.get(kind, ())
    xs = np.array([p[0] for p in pts]) / (turns or 1)
    grid = np.linspace(0.0, float(xs.max()) * 1.02, 80)
    out = {}
    for mp in moments:
        if mp not in have:
            continue
        n, comp = MOMENT[mp]
        if kind in SINGLE_KINDS:
            # Однообмоточные: граф-модель с утечкой, подогнанная на картах ЭТОГО исполнения.
            try:
                out[mp] = single_model_curve(el, L, pts, n, turns, ef, label)
            except Exception:
                pass
            continue
        try:
            # Модель даёт ИНТЕГРАЛ момента (её физкоэф отнесён к своей опорной длине). Шкала
            # графика — та же, что у измерений: интеграл, делённый на длину ЭТОГО исполнения.
            out[mp] = (grid, np.array(
                [abs(ss.graph_moment(el, {wkey: float(I)}, n, comp)[1]) / (L * ef)
                 for I in grid]))
        except Exception:
            continue
    return out


def plot_winding(m, label, pts, turns, L, moments, design, cur, ef, outdir, model=None):
    """Кривая намагниченности обмотки: физкоэффициенты в зависимости от тока питания.

    Линия — модель насыщения (`model`: {момент: (ток_А, значение)}), точки — измерения по картам.
    Проектная точка обозначена линиями: горизонталь — проектный физкоэффициент из .lte, вертикаль —
    проектный ток обмотки с учётом насыщения, ромб — их пересечение."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    model = model or {}
    x = np.array([p[0] for p in pts]) / (turns or 1)
    cb = float(np.mean([p[2] for p in pts]))
    # Каждый момент живёт на своей оси; при пяти-шести моментах оси съедают поле, поэтому
    # холст расширяется под их число, а заголовок ставится на ФИГУРУ, иначе его срезает.
    fig, ax = plt.subplots(figsize=(9.0 + 0.85 * max(0, len(moments) - 3), 5.2))
    axes = {mp: (ax if k == 0 else ax.twinx()) for k, mp in enumerate(moments)}
    for k, mp in enumerate(moments):
        n, comp = MOMENT[mp]
        a, col = axes[mp], MOMENT_COLOR[mp]
        y = np.array([to_coef(p[1][n][comp], mp, L, ef) for p in pts])
        a.plot(x, y, "o", color=col, ms=5, mec="white", mew=0.7, zorder=3)   # измерения — точки
        if mp in model:                                                       # модель — линия
            mx, my = model[mp]
            a.plot(mx, my, "-", color=col, lw=1.8, zorder=2)
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
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.suptitle(f"{m['name']} — обмотка {label}: " + ", ".join(moments)
                 + f";  L_eff = {L:.3f} м"
                 + ("" if abs(cb - 1) < 0.03 else f" (тело карты ×{cb:.2f})"), fontsize=10.5)
    fig.subplots_adjust(top=0.90)
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, f"{m['name']}_{label}.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png


def plot_material(path, outdir):
    """Кривая намагничивания материала B(H): линейный и логарифмический масштаб по H.

    Логарифм нужен, потому что таблица покрывает четыре порядка по напряжённости: в линейном
    масштабе весь начальный участок, где сталь и работает, сжимается в отрезок у нуля."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import pandas as pd
    d = pd.read_csv(path, sep=";")
    H, B = d.iloc[:, 0].to_numpy(float), d.iloc[:, 1].to_numpy(float)
    fig, ax = plt.subplots(1, 2, figsize=(9.6, 4.2))
    for a, log in zip(ax, (False, True)):
        a.plot(H, B, "-o", color="#1f4e79", ms=4, lw=1.6)
        if log:
            a.set_xscale("log")
            a.set_xlim(max(H[H > 0].min() * 0.6, 1e-1), H.max() * 1.4)
        else:
            a.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
        a.set_xlabel("H, А/м")
        a.set_ylabel("B, Тл")
        a.grid(alpha=0.25, which="both")
        a.set_title("логарифмический масштаб по H" if log else "линейный масштаб",
                    fontsize=10)
    fig.suptitle(f"Кривая намагничивания материала: {os.path.basename(path)}", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    os.makedirs(outdir, exist_ok=True)
    png = os.path.join(outdir, "material_bh.png")
    fig.savefig(png, dpi=150)
    plt.close(fig)
    return png, H, B


def material_section(doc, figs, path=MATERIAL):
    """Заключительный раздел: кривая намагничивания материала и таблица её точек."""
    from docx.shared import Inches
    if not os.path.exists(path):
        return False
    doc.add_heading("Кривая намагничивания материала", level=1)
    png, H, B = plot_material(path, figs)
    doc.add_picture(png, width=Inches(6.3))
    p = doc.add_paragraph()
    p.add_run(f"Материал сердечника: {os.path.basename(path)}; точек в таблице: {len(H)}; "
              f"насыщение при B = {B.max():.4f} Тл.")
    t = doc.add_table(rows=1, cols=3)
    t.style = "Light Grid Accent 1"
    for k, h in enumerate(("№", "H, А/м", "B, Тл")):
        t.rows[0].cells[k].text = h
    for i, (h_, b_) in enumerate(zip(H, B), 1):
        c = t.add_row().cells
        c[0].text, c[1].text, c[2].text = str(i), f"{h_:.1f}", f"{b_:.4f}"
    doc.add_paragraph()
    return True


# ───────────────────────── дополнение отчёта ─────────────────────────

def mark_linear_tables(doc):
    """Пометить существующие таблицы проектных токов как ЛИНЕЙНОЕ ПРИБЛИЖЕНИЕ.

    Эти таблицы в исходном отчёте получены обращением линейной матрицы связи; ниже в каждой главе
    добавляется таблица того же формата, но с учётом насыщения. Помечаем подпись, чтобы читатель
    не путал два расчёта."""
    n = 0
    for p in doc.paragraphs:
        t = p.text.strip()
        if "проектные токи обмоток" in t and "линейн" not in t.lower():
            for r in p.runs:
                if "проектные токи обмоток" in r.text:
                    r.text = r.text.rstrip() + " (линейное приближение)"
                    n += 1
                    break
    return n


def _cell_text(tc, text):
    """Заменить текст ячейки, сохранив её оформление (ячейка приходит копией соседней)."""
    from docx.table import _Cell
    c = _Cell(tc, None)
    for extra in c.paragraphs[1:]:
        extra._p.getparent().remove(extra._p)
    par = c.paragraphs[0]
    for extra in par.runs[1:]:
        extra._r.getparent().remove(extra._r)
    if par.runs:
        par.runs[0].text = text
    else:
        par.add_run(text)


def insert_column(table, after, header, values, key_col=0):
    """Вставить столбец сразу за столбцом `after`: заголовок и значения по имени в `key_col`.

    Ячейки создаются копией соседних, поэтому рамки, заливка и шрифт совпадают со столбцом, рядом
    с которым столбец встаёт; python-docx умеет добавлять столбец только в конец таблицы."""
    import copy
    from docx.oxml.ns import qn
    grid = table._tbl.find(qn("w:tblGrid"))
    cols = grid.findall(qn("w:gridCol"))
    grid.insert(list(grid).index(cols[after]) + 1, copy.deepcopy(cols[after]))
    for i, tr in enumerate(table._tbl.findall(qn("w:tr"))):
        tcs = tr.findall(qn("w:tc"))
        new = copy.deepcopy(tcs[after])
        tcs[after].addnext(new)
        if i == 0:
            _cell_text(new, header)
        else:
            name = _Cell_text(tcs[key_col])
            _cell_text(new, values.get(name, "—"))


def _Cell_text(tc):
    from docx.table import _Cell
    return _Cell(tc, None).text.strip()


def summary_table(doc):
    """Сводная таблица отчёта: та, где в шапке есть столбец проектных токов."""
    for t in doc.tables:
        head = [_Cell_text(c._tc) for c in t.rows[0].cells]
        if any(h.startswith("проектные токи") for h in head):
            return t, [h.startswith("проектные токи") for h in head].index(True)
    return None, None


def page_break_chapters(doc):
    """Каждая глава — с новой страницы (кроме первой: перед ней и так начало документа)."""
    n = 0
    for i, par in enumerate(p for p in doc.paragraphs if p.style.name == "Heading 1"):
        if i:
            par.paragraph_format.page_break_before = True
            n += 1
    return n


def build(out=OUT, src=SRC, figs=PNG):
    """Достоверная копия element_compare_report, дополненная по элементам: кривые намагниченности,
    формула модели (у квадруполей/секступолей/октуполей) и таблица проектных токов С УЧЁТОМ
    НАСЫЩЕНИЯ (той же формы, что линейная). Существующие линейные таблицы помечаются."""
    from docx import Document
    from docx.shared import Inches
    import design_currents as dc
    import element_compare as ec
    import latex2omml as lo

    if not os.path.exists(src):
        raise SystemExit(f"нет исходного отчёта: {src}\n"
                         f"построить: py -3.14 element_compare.py library")
    doc = Document(src)
    mark_linear_tables(doc)
    sat_by_name = {}
    sm = _smodels()
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
        kind = m["kind"]
        lead = lead_order(m)
        L = float(m["L"]) if m["L"] else map_length(el, lead)
        pts_by_w = variant_kicks(el, L, lead)
        native = element_spectra(el)
        allow = KIND_MOMENTS.get(kind, tuple(MOMENT))
        design = ec._design_ivec(m)[2]
        cur = sat_currents_solve(el, kind, m, design, ec, pts_by_w, L)
        if kind != "kicker":
            sat_by_name[m["name"]] = ", ".join(f"{lb}={cur.get(lb, 0.0):+.1f}"
                                               for lb in m["labels"])

        new = []

        def add(node):
            new.append(node)

        # Кривые физкоэффициентов от тока (у всех элементов с картами; corrector — только они).
        add(doc.add_paragraph(f"Кривые физкоэффициентов {m['name']} от тока обмоток "
                              f"(по трёхмерным картам поля; линии — модель, точки — измерения)")._p)
        for j, lb in enumerate(m["labels"]):
            pts = pts_by_w.get(lb) or pts_by_w.get("I1" if lb == "I" else "I")
            nat = native.get(lb) or native.get("I1" if lb == "I" else "I")
            if not pts or not nat:
                continue
            shown = [mp for mp in MOMENT
                     if mp in allow and resolution(nat, MOMENT[mp][0]) <= RESOLVED
                     and (MOMENT[mp][1] == 0 or mp in m["outputs"])
                     and moment_possible(el, kind, lb, MOMENT[mp][1])]
            if not shown:
                continue
            tpc = ec.dev_tpc(m["devs"][j])
            model = model_curves(el, kind, lb, shown, pts, tpc, ec.ENERGY_FACTOR, L)
            png = plot_winding(m, lb, pts, tpc, L, shown, design, cur.get(lb),
                               ec.ENERGY_FACTOR, figs, model=model)
            pic = doc.add_paragraph()
            pic.add_run().add_picture(png, width=Inches(6.3))
            add(pic._p)

        # Формула модели — только у однообмоточных (квадруполь, секступоль, октуполь).
        par_lead = None
        if kind in SINGLE_KINDS:
            w0 = m["labels"][0]
            pts0 = pts_by_w.get(w0) or pts_by_w.get("I1") or pts_by_w.get("I")
            if pts0:
                try:
                    par_lead = single_graph_fit(el, L, pts0, MOMENT[LEAD_MOMENT[kind]][0], w0)[0]
                except Exception:
                    par_lead = None
        if par_lead is not None:
            marker = doc.add_paragraph()._p            # сдвиг: docx кладёт абзацы в конец, переносим
            formula_paragraph(doc, el, kind, par_lead, L, lo)
            body = doc.paragraphs
            i0 = body.index(next(p for p in body if p._p is marker))
            for p in body[i0 + 1:]:
                add(p._p)

        # Таблица проектных токов с учётом насыщения — у всех, КРОМЕ корректора (там только кривые).
        if kind != "kicker":
            ctrl = [mp for mp in MOMENT if mp in design and abs(design[mp]) > 1e-30
                    and (kind != "dq" or mp in ("B", "G", "S"))]
            if kind == "dq":
                add(doc.add_paragraph(
                    "Диполь-квадруполь несимметричен: направление влёта пучка существенно для "
                    "токов. Магнит установлен развёрнутым на 180° вокруг вертикали (rot_z = 180°) — "
                    "это переворачивает знак квадрупольного градиента при сохранении дипольного "
                    "поля, благодаря чему проектные B и G достигаются токами, близкими к "
                    "номинальным. Токи ниже приведены для этой установленной ориентации.")._p)
            cap = doc.add_paragraph()
            cap.add_run(f"{m['name']} — проектные токи с учётом насыщения").bold = True
            add(cap._p)
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
                nom = ec.dev_nominal(m["devs"][j])
                c[3].text = f"{I:+.1f}" if I is not None else "—"
                c[4].text = f"{100 * abs(I) / nom:.0f}%" if (I is not None and nom) else "—"
            add(t._tbl)
            add(doc.add_paragraph()._p)

        for node in new:                                   # перенос созданных блоков в конец главы
            anchor.addprevious(node)
        added += 1

    if material_section(doc, figs):                 # заключительный раздел — в САМЫЙ конец
        added += 1
    tbl, col = summary_table(doc)
    if tbl is not None:                     # сводка: старый столбец пометить, новый вставить рядом
        _cell_text(tbl.rows[0].cells[col]._tc, "проектные токи, А (линейное приближение)")
        insert_column(tbl, col, "проектные токи, А (с насыщением)", sat_by_name)
    print("глав с новой страницы:", page_break_chapters(doc))
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
