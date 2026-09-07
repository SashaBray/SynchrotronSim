# -*- coding: utf-8 -*-
"""Word-методичка: модель «токи обмоток → кики» как расчётная функция, пригодная для счёта руками.

Плотное изложение всей математики цепочки от вектора токов к вектору киков: матрица намотки,
уравнения магнитной цепи с явным шагом Ньютона и проводимостями ветвей, свёртка потоков киком через
матрицу связи, перевод в единицы решётки. Обозначения собраны в один справочный список. Завершает
методичку РАЗОБРАННЫЙ ЧИСЛОВОЙ ПРИМЕР для однообмоточного элемента — его можно повторить на бумаге.

Все числа примера и таблиц считаются живой моделью при сборке.

    py -3.14 model_report.py
    py -3.14 model_report.py --out путь/отчёт.docx
"""
from __future__ import annotations

import argparse
import csv
import io
import os
import sys
import warnings

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join("workspace", "results", "currents_to_kicks_model.docx")
FITS = os.path.join("workspace", "materials", "saturation_fits.csv")
ELEMENTS = ("QP_sh", "QP_L", "SP", "OP", "DQ", "corrector")

IO = {
    "QP_sh": ("I1", "G"), "QP_L": ("I1", "G"), "OP": ("I1", "O"),
    "SP": ("I1, I2", "S"), "DQ": ("I1, I2, I3", "B, G, S"),
    "corrector": ("I1…I5", "B, G, S, Bs, Gs, Ss"),
}
# столбец интеграла гармоники в датасете, по которому подгоняется однообмоточное уравнение
SINGLE_COEF = {"QP_sh": "int_G", "QP_L": "int_G", "SP": "int_S", "OP": "int_O"}


# ─────────────────────────────── оформление ────────────────────────────────

def _formula(doc, latex, lo):
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    lo.add_formula(p, latex)


def _where(doc, lo, triples):
    """Пояснения переменных под формулой: «Где:», затем на каждую строку символ НАБРАН ФОРМУЛОЙ
    (нативный OMML, как в самой формуле) и пояснение мелким шрифтом."""
    from docx.shared import Cm, Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    h = doc.add_paragraph()
    h.paragraph_format.space_after = Pt(1)
    h.add_run("Где:").font.size = Pt(9)
    for sym, kind, txt in triples:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(1.0)
        p.paragraph_format.first_line_indent = Cm(-1.0)
        p.paragraph_format.space_after = Pt(1)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.add_run("— ").font.size = Pt(9)
        lo.add_formula(p, sym)
        p.add_run(f" — {txt}; {kind}.").font.size = Pt(9)
    doc.add_paragraph()          # отступ строкой после пояснений


def _table(doc, header, rows, caption=None, pt=9):
    from docx.shared import Pt
    if caption:
        _table.n = getattr(_table, "n", 0) + 1
        c = doc.add_paragraph()
        r = c.add_run(f"Таблица {_table.n}. {caption}")
        r.bold = True
        r.font.size = Pt(9)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Light Grid Accent 1"
    for i, h in enumerate(header):
        t.rows[0].cells[i].text = h
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = str(v)
    for row in t.rows:
        for cell in row.cells:
            for par in cell.paragraphs:
                for run in par.runs:
                    run.font.size = Pt(pt)


def _pic(doc, path, cm=15.0, caption=None):
    from docx.shared import Cm, Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    if not os.path.exists(path):
        return
    doc.add_picture(path, width=Cm(cm))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER
    if caption:
        _pic.n = getattr(_pic, "n", 0) + 1
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        r = p.add_run(f"Рисунок {_pic.n}. {caption}")
        r.bold = True
        r.font.size = Pt(9)


def _body_format(doc):
    from docx.shared import Cm, Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    for p in doc.paragraphs:
        p.paragraph_format.space_after = Pt(4)
        if not p.text.strip() or p.style.name.startswith("Heading") or p.style.name == "Title":
            continue
        if p.alignment == WD_ALIGN_PARAGRAPH.CENTER:
            continue
        if p.text.startswith("— ") or p.text.strip() == "Где:":   # пояснения — свой отступ
            continue
        p.paragraph_format.first_line_indent = Cm(1.0)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY


# ─────────────────────────────── живые числа ────────────────────────────────

COEF = {"QP_sh": ("G", "Тл/м", "int_G"), "QP_L": ("G", "Тл/м", "int_G"),
        "SP": ("S", "Тл/м²", "int_S"), "OP": ("O", "Тл/м³", "int_O")}
NTURNS = {"QP_sh": "68", "QP_L": "85", "SP": "41+390", "OP": "34"}


def _leff(element):
    """Средняя эффективная длина элемента по записям датасета (проектная, фиксированная)."""
    import json
    import saturation_model as sm
    recs = json.load(open(os.path.join(sm.OUTDIR, f"{element}.json"), encoding="utf-8"))["records"]
    xs = [float(r["Leff"]) for r in recs if r.get("Leff")]
    return sum(xs) / len(xs)


def facts():
    import numpy as np
    import saturation_model as sm

    F = {"fits": {}, "single": [], "example": None}
    for r in csv.DictReader((l for l in io.open(FITS, encoding="utf-8-sig")
                             if not l.startswith("#")), delimiter=";"):
        F["fits"].setdefault(r["section"], []).append(r)

    fit = {r["section"]: r for rs in F["fits"].values() for r in rs}
    for el in ("QP_sh", "QP_L", "SP", "OP"):
        sym, unit, _ = COEF[el]
        r = fit[el]
        F["single"].append(dict(el=el, sym=sym, unit=unit, N=NTURNS[el], leff=_leff(el),
                                c=float(r["C"]), l=float(r["L_m"]),
                                err=float(r["max_err_pct"]), n=r["n_points"]))

    # Разобранный пример: QP_sh. Цепочка I → NI=N·I → χ=∫G·dl → G=χ/l_eff.
    Hof = sm._material()
    H = lambda b: float(np.atleast_1d(Hof(b))[0])
    q = fit["QP_sh"]
    C, L, aon = float(q["C"]), float(q["L_m"]), float(q["A_over_naive"])
    N, leff = 68.0, _leff("QP_sh")
    rows = sorted(csv.DictReader(io.open(os.path.join("workspace", "datasets",
                  "currents_to_kicks", "QP_sh.csv"), encoding="utf-8-sig"), delimiter=";"),
                  key=lambda r: float(r["AT_I"]))
    ni0, chi0 = float(rows[0]["AT_I"]), abs(float(rows[0]["int_G"]))
    A = aon * ni0 / chi0
    picks, seen = [], set()
    for r in rows:
        ni, chi = float(r["AT_I"]), abs(float(r["int_G"]))
        band = round(ni / 1600)
        if ni > 6000 or band in seen:
            continue
        seen.add(band)
        nim = A * chi + L * H(C * chi)
        picks.append((ni / N, ni, chi, chi / leff, nim, 100 * (nim - ni) / ni))
    F["example"] = dict(A=A, C=C, L=L, N=N, leff=leff, ni0=ni0, chi0=chi0, rows=picks[:4])

    F["dq"] = dq_example()
    return F


DQ_N = {"I1": 53, "I2": 13, "I3": 90}
DQ_LEFF = 0.4081
DQ_R = 8e-3
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0}


def _dq_coef(kappa, n):
    """Физкоэффициент DQ гармоники n из кика (мнимая часть — нормальная составляющая)."""
    import math
    return MAT_SIGN[n] * math.factorial(n - 1) / DQ_R ** (n - 1) * kappa.imag / DQ_LEFF


def dq_example(currents_A=(("I1", 150.0), ("I2", 100.0), ("I3", 40.0))):
    """Пример расчёта DQ: заданные токи → σ → B, G, S. Числа считаются моделью."""
    import numpy as np
    import branch_network as bn
    import coupling_model as cm
    A = dict(currents_A)
    NI = {k: A[k] * DQ_N[k] for k in A}
    cir, _, _ = cm.train("DQ")
    lens, perms, secs = bn.unpack(np.zeros(1), cir)
    phi, _ = cir.solve(NI, lens, perms, secs=secs)
    num = den = 0.0
    for x in cir.iron:
        S = secs.get(cir.sym[x], 1.0)
        a = abs(phi[x])
        num += a * min(a / S / cir.Barr[-1], 1.0)
        den += a
    sigma = num / den if den else 0.0
    coef = {n: _dq_coef(cm.predict("DQ", n, NI), n) for n in (1, 2, 3)}
    return dict(A=A, N=DQ_N, NI=NI, sigma=sigma, leff=DQ_LEFF,
                B=coef[1], G=coef[2], S=coef[3])


# ─────────────────────────────── разделы ────────────────────────────────

def ch_task(doc, F, lo):
    doc.add_heading("1. Задача и обозначения", level=1)
    doc.add_paragraph(
        "Модель: на входе вектор возбуждений обмоток (в ампер-витках, число витков уже учтено), "
        "на выходе вектор кик-коэффициентов магнита.")
    _formula(doc, r"\vec{k} = \mathcal{K}(\vec{I})", lo)
    _where(doc, lo, [
        (r"\vec{I}", "вектор, ампер-витки", "возбуждения обмоток (витки уже учтены)"),
        (r"\vec{k}", "вектор", "кик-коэффициенты магнита: нормальные и скошенные, по гармоникам"),
        (r"\mathcal{K}", "вектор-функция", "модель — композиция четырёх звеньев (раздел 2)"),
    ])
    doc.add_paragraph("Размерности входа и выхода заданы назначением элемента:")
    _table(doc, ["Элемент", "Вход", "Выход"], [[el, IO[el][0], IO[el][1]] for el in ELEMENTS],
           caption="Вход и выход элементов библиотеки")


def ch_math(doc, F, lo):
    doc.add_heading("2. Уравнения расчёта", level=1)
    doc.add_paragraph("Функция есть композиция четырёх звеньев; нелинейно только второе.")

    p = doc.add_paragraph()
    p.add_run("Намотка. ").bold = True
    p.add_run("Ток обмотки создаёт МДС на охватываемых ею рёбрах, линейно:")
    _formula(doc, r"\vec{F} = M\,\vec{I}", lo)
    _formula(doc, r"M_{bw} = s_w\,\alpha_{wb}", lo)
    _where(doc, lo, [
        (r"\vec{F}", "вектор, ампер-витки", "МДС на каждом ребре"),
        (r"M", "матрица", "матрица намотки «рёбра × обмотки», постоянная"),
        (r"s_w", "скаляр, ±1", "знак намотки обмотки"),
        (r"\alpha_{wb}", "скаляр", "доля МДС обмотки w на ребре b (0, если не охватывает)"),
    ])

    p = doc.add_paragraph()
    p.add_run("Магнитная цепь. ").bold = True
    p.add_run("Напор на ребре, поток ребра и закон сохранения потока в каждой вершине дают "
              "нелинейную систему относительно потенциалов вершин:")
    _formula(doc, r"\vec{d} = \vec{F} - A^{\mathsf{T}}\vec{u}", lo)
    _formula(doc, r"\varphi_b = \begin{cases} p_b\,d_b, & \text{воздух},\\ "
                  r"s_b\,B(d_b/\ell_b), & \text{железо},\end{cases}", lo)
    _formula(doc, r"A\,\vec{\varphi}(\vec{u}) = 0", lo)
    _where(doc, lo, [
        (r"\vec{d}", "вектор, ампер-витки", "магнитный напор на каждом ребре"),
        (r"A", "матрица, ±1 и 0", "инцидентная матрица графа «вершины × рёбра»"),
        (r"\vec{u}", "вектор, ампер-витки", "потенциалы вершин (одна опорная, с нулём)"),
        (r"\vec{\varphi}", "вектор, вебер", "потоки рёбер — искомое (компонента φ_b)"),
        (r"p_b", "скаляр, генри", "проводимость воздушного ребра, p_b = μ₀·s_b/ℓ_b"),
        (r"s_b,\ \ell_b", "скаляры, кв. м и м", "сечение и длина железного ребра"),
        (r"B(\cdot)", "функция, тесла", "кривая намагничивания материала"),
    ])
    doc.add_paragraph(
        "Система решается методом Ньютона по потенциалам вершин: на каждой итерации по невязке "
        "A·φ(u) и якобиану из проводимостей ветвей потенциалы поправляются")
    _formula(doc, r"\vec{u} \leftarrow \vec{u} - J^{-1}\,A\,\vec{\varphi}(\vec{u})", lo)
    _formula(doc, r"J = -A\,\operatorname{diag}(\vec{g})\,A^{\mathsf{T}}", lo)
    _formula(doc, r"g_b = \begin{cases} p_b, & \text{воздух},\\ "
                  r"(s_b/\ell_b)\,B'(d_b/\ell_b), & \text{железо}.\end{cases}", lo)
    _where(doc, lo, [
        (r"J", "матрица", "якобиан системы по потенциалам"),
        (r"\vec{g}", "вектор", "проводимости ветвей g_b (компоненты)"),
        (r"g_b", "скаляр", "проводимость ребра ∂φ_b/∂d_b: у воздуха постоянна, у железа зависит "
         "от наклона кривой B′"),
    ])
    p = doc.add_paragraph()
    p.add_run("Признаки и насыщение. ").bold = True
    p.add_run("Поле у пучка создают зазоры между наконечниками и потоки катушек (в графе предусматривается "
              "воздушные зазоры между катушками). "
              "Признак — поток узла (его компонента в чистом потоке через A); форма поля меняется "
              "с насыщением, мера которого — потоко-взвешенная доля:")
    _formula(doc, r"q_t = (A\,\vec{\varphi})_t", lo)
    _formula(doc, r"\sigma = \frac{\sum |\varphi_b|\,\min(|\varphi_b|/s_b b_{\max},\,1)}"
                  r"{\sum |\varphi_b|}", lo)
    _where(doc, lo, [
        (r"\vec{q}", "вектор, вебер", "потоки наконечников — признаки кика"),
        (r"q_t", "скаляр, вебер", "поток наконечника t (компонента q̄)"),
        (r"\sigma", "скаляр, 0…1", "доля насыщения железа"),
        (r"b_{\max}", "скаляр, тесла", "предел кривой намагничивания"),
    ])

    p = doc.add_paragraph()
    p.add_run("Кик и единицы. ").bold = True
    p.add_run("Кик — свёртка признаков матрицей связи с поправкой на насыщение; мнимая часть даёт "
              "нормальную составляющую, вещественная — скошенную. Итог приводится к "
              "физкоэффициентам при фиксированной эффективной длине:")
    _formula(doc, r"\vec{k} = W^{\mathsf{T}}\begin{bmatrix}\vec{q}\\ "
                  r"\sigma\,\vec{q}\end{bmatrix}", lo)
    _formula(doc, r"c_n = \varsigma_n\,\frac{(n-1)!}{r^{\,n-1}\,l_{\text{эфф}}}\;\kappa_n", lo)
    _where(doc, lo, [
        (r"W", "матрица", "матрица связи «признаки × компоненты кика», подбираемая"),
        (r"\vec{w}_n,\ \vec{v}_n", "векторы", "строки W: веса при ненасыщенном железе и "
         "поправка, набирающая силу с σ"),
        (r"\kappa_n", "скаляр, комплексный", "кик гармоники n (компонента k̄)"),
        (r"c_n", "скаляр, комплексный", "физкоэффициент гармоники n (G, S, O — нормальный и "
         "скошенный)"),
        (r"\varsigma_n", "скаляр, ±1", "знаковая конвенция гармоники"),
        (r"r,\ l_{\text{эфф}}", "скаляры, м", "радиус пробы и эффективная длина элемента"),
    ])
    doc.add_paragraph(
        "Матрицы M и A и величины s_b, ℓ_b взяты из чертежей. Подбираются по картам поля только "
        "матрица связи W (линейна, находится точным МНК) и сечения воздушных рёбер утечки; знаки "
        "намотки выбираются из конечного набора. Один прогон функции — это умножение на M, одно "
        "решение цепи Ньютоном, умножение на A, подсчёт σ и две свёртки; итерации есть только в "
        "решении цепи.")


def ch_single(doc, F, lo):
    doc.add_heading("3. Однообмоточный элемент", level=1)
    doc.add_paragraph(
        "У квадруполей, секступоля и октуполя обмотка одна, поток идёт единственным путём, граф "
        "вырождается в одно ребро. Расчёт от тока к физкоэффициенту идёт в три шага. Сначала ток "
        "переводится в возбуждение:")
    _formula(doc, r"NI = N\,I", lo)
    _where(doc, lo, [
        (r"I", "скаляр, ампер", "ток обмотки"),
        (r"N", "скаляр", "число витков катушки"),
        (r"NI", "скаляр, ампер-витки", "возбуждение"),
    ])
    doc.add_paragraph(
        "Возбуждение связано с интегральной силой гармоники χ (интеграл коэффициента вдоль "
        "пучка) скалярным уравнением магнитной цепи:")
    _formula(doc, r"NI = a\,\chi + l\,H(c\,\chi)", lo)
    _where(doc, lo, [
        (r"\chi", "скаляр", "интегральная сила гармоники, χ = ∫G·dl (или ∫S·dl, ∫O·dl)"),
        (r"a", "скаляр", "линейный член: воздух, не насыщается"),
        (r"l", "скаляр, м", "длина пути по железу"),
        (r"c", "скаляр", "пересчёт χ в индукцию железа"),
        (r"H(\cdot)", "функция, А/м", "обратная кривая намагничивания"),
    ])
    _formula(doc, r"G = \chi / l_{\text{эфф}}", lo)
    _where(doc, lo, [
        (r"G", "скаляр, Тл/м", "физкоэффициент (для секступоля S, Тл/м²; октуполя O, Тл/м³)"),
        (r"l_{\text{эфф}}", "скаляр, м", "эффективная длина элемента (проектная, фиксированная)"),
    ])
    doc.add_paragraph("Постоянные, число витков и эффективная длина по элементам:")
    rows = [[s["el"], f'{s["sym"]}, {s["unit"]}', s["N"], f'{s["leff"]:.4f}',
             f'{s["c"]:.4g}', f'{s["l"]:.4g}', f'{s["err"]:.1f}', s["n"]]
            for s in F["single"]]
    _table(doc, ["Элемент", "Коэффициент", "N", "l_эфф, м", "c", "l, м", "макс. ошибка, %",
                 "Точек"], rows, caption="Постоянные однообмоточного уравнения", pt=8)

    doc.add_heading("Порядок расчёта", level=2)
    doc.add_paragraph("Дано ток I, требуется физкоэффициент G. Порядок такой:")
    for t in ("возбуждение NI = N·I;",
              "взять a, c, l элемента из таблицы 2 (a — через контрольную точку, см. пример);",
              "решить уравнение цепи относительно χ: задаться пробным χ, вычислить "
              "a·χ + l·H(c·χ), взяв H по таблице материала при индукции c·χ, сравнить с NI и "
              "поправить χ (уравнение монотонно, годится деление отрезка пополам; 2–3 приближения);",
              "физкоэффициент G = χ / l_эфф."):
        doc.add_paragraph(t, style="List Number")

    ex = F["example"]
    doc.add_heading("Числовой пример (QP_sh)", level=2)
    doc.add_paragraph(
        f"Дано: N = {ex['N']:.0f} витков, l_эфф = {ex['leff']:.4f} м, постоянные c = {ex['C']:.4g}, "
        f"l = {ex['L']:.4g} м. Множитель a по контрольной точке NI₀ = {ex['ni0']:.0f} А·в, "
        f"χ₀ = {ex['chi0']:.3g}: a = {ex['ni0']:.0f}/{ex['chi0']:.3g} ≈ {ex['A']:.0f}. Для нескольких "
        f"токов проходим цепочку I → NI → χ → G и сверяем расчётное возбуждение с картой:")
    rows = [[f"{i:.1f}", f"{ni:.0f}", f"{chi:.3f}", f"{g:.2f}", f"{nim:.0f}", f"{err:+.1f}"]
            for i, ni, chi, g, nim, err in ex["rows"]]
    _table(doc, ["I, А", "NI = N·I", "χ = ∫G·dl", "G = χ/l_эфф, Тл/м", "NI расч.", "ошибка, %"],
           rows, caption="Проверка цепочки на рабочих картах QP_sh", pt=8)



def _dq_figures(figdir):
    """Готовит две иллюстрации для главы 4: граф DQ и графики моментов B, G, S от токов."""
    import json
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import branch_network as bn
    import coupling_model as cm
    import plot_magnet_graph as pmg

    os.makedirs(figdir, exist_ok=True)
    graph = os.path.join(figdir, "dq_graph.png")
    pmg.draw("DQ", out=graph)

    cm.train("DQ")
    turns = json.load(open(os.path.join("workspace", "datasets", "currents_to_kicks",
                                         "DQ.meta.json"), encoding="utf-8"))["turns_per_coil"]
    ORD = ((1, "B", "Тл"), (2, "G", "Тл/м"), (3, "S", "Тл/м²"))
    maps = {n: bn.read_maps("DQ", n) for n, _, _ in ORD}
    wnames = ["I1", "I2", "I3"]
    fig, ax = plt.subplots(3, 3, figsize=(11, 8), squeeze=False)
    for i, (n, nm, u) in enumerate(ORD):
        for k, w in enumerate(wnames):
            a = ax[i][k]
            pts = sorted((c[w], v) for c, v, *_ in maps[n] if list(c) == [w])
            if pts:
                at = np.array([p[0] for p in pts])
                meas = np.array([_dq_coef(complex(p[1]), n) for p in pts])
                grid = np.linspace(0, at.max() * 1.02, 80)
                mod = np.array([_dq_coef(cm.predict("DQ", n, {w: g}), n) for g in grid])
                a.plot(grid / turns[w], mod, "-", color="#1f4e79", lw=1.7)
                a.plot(at / turns[w], meas, "o", color="#1f4e79", ms=5, mec="white", mew=.7)
            a.grid(alpha=.25); a.tick_params(labelsize=8); a.axhline(0, color="#999", lw=.6)
            if i == 0:
                a.set_title(f"обмотка {w}", fontsize=10)
            if i == 2:
                a.set_xlabel("ток, А", fontsize=9)
            if k == 0:
                a.set_ylabel(f"{nm}, {u}", fontsize=10)
    from matplotlib.lines import Line2D
    fig.legend(handles=[Line2D([], [], color="#1f4e79", lw=1.7, label="модель"),
                        Line2D([], [], color="#1f4e79", marker="o", ls="", label="карты")],
               loc="lower center", ncol=2, fontsize=10, frameon=False)
    fig.suptitle("Диполь-квадруполь: B, G, S от тока каждой обмотки", fontsize=12)
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    moments = os.path.join(figdir, "dq_moments.png")
    fig.savefig(moments, dpi=140)
    plt.close(fig)
    return graph, moments


def ch_dq_example(doc, F, lo):
    dq = F["dq"]
    doc.add_heading("4. Пример расчёта: диполь-квадруполь", level=1)
    doc.add_paragraph(
        "Диполь-квадруполь несёт три обмотки и даёт три коэффициента — B, G, S. Здесь цепочка "
        "раздела 2 пройдена с числами.")
    doc.add_paragraph(
        f"Дано: токи обмоток I1 = {dq['A']['I1']:.0f} А, I2 = {dq['A']['I2']:.0f} А, "
        f"I3 = {dq['A']['I3']:.0f} А. Числа витков N1 = {dq['N']['I1']}, N2 = {dq['N']['I2']}, "
        f"N3 = {dq['N']['I3']}; эффективная длина l_эфф = {dq['leff']:.4f} м.")

    doc.add_paragraph("Шаг 1. Возбуждение ветвей — токи через матрицу намотки:")
    _formula(doc, r"\vec{F} = M\,\vec{I}", lo)
    doc.add_paragraph(
        f"ампер-витки обмоток: NI1 = {dq['NI']['I1']:.0f}, NI2 = {dq['NI']['I2']:.0f}, "
        f"NI3 = {dq['NI']['I3']:.0f}.")

    doc.add_paragraph("Шаг 2. Решение цепи методом Ньютона даёт потоки всех ветвей:")
    _formula(doc, r"A\,\vec{\varphi}(\vec{u}) = 0", lo)
    doc.add_paragraph(
        f"по потокам считается доля насыщения σ = {dq['sigma']:.2f} (железо на этом режиме "
        "насыщено примерно наполовину).")

    doc.add_paragraph("Шаг 3. Кик — свёртка потоков наконечников матрицей связи, затем перевод в "
                      "коэффициенты:")
    _formula(doc, r"\vec{k} = W^{\mathsf{T}}\begin{bmatrix}\vec{q}\\ \sigma\,\vec{q}\end{bmatrix}", lo)
    _formula(doc, r"c_n = \varsigma_n\,\frac{(n-1)!}{r^{\,n-1}\,l_{\text{эфф}}}\;\kappa_n", lo)
    doc.add_paragraph("Результат:")
    _table(doc, ["Коэффициент", "Значение", "Единицы"],
           [["B (диполь)", f"{dq['B']:+.3g}", "Тл"],
            ["G (квадруполь)", f"{dq['G']:+.4g}", "Тл/м"],
            ["S (секступоль)", f"{dq['S']:+.4g}", "Тл/м²"]],
           caption="Физкоэффициенты диполя-квадруполя при заданных токах")

    graph, moments = F["dq_figs"]
    _pic(doc, graph, 14.0, caption="Магнитный граф диполя-квадруполя: четыре полюса, ярмо, "
         "шунты полюсов и зазоры между наконечниками")
    _pic(doc, moments, 16.5, caption="Зависимости коэффициентов B, G, S от тока каждой обмотки: "
         "линии — модель, точки — карты поля")


def build(out=OUT):
    import docx
    import latex2omml as lo
    from docx.shared import Pt

    F = facts()
    F["dq_figs"] = _dq_figures(os.path.join(os.path.dirname(out) or ".", "model_report_figs"))
    doc = docx.Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(11)

    doc.add_heading("Модель «токи обмоток → кики»: методика расчёта", level=0)
    p = doc.add_paragraph()
    p.add_run("Уравнения всех звеньев функции ток→кик и разобранные числовые примеры для расчёта "
              "вручную. Связь потоков с киком — матрицей связи, без координат полюсов.").italic = True

    for ch in (ch_task, ch_math, ch_single, ch_dq_example):
        ch(doc, F, lo)

    _body_format(doc)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
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
