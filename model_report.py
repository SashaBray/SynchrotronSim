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

def facts():
    import numpy as np
    import saturation_model as sm

    F = {"fits": {}, "example": None}
    for r in csv.DictReader((l for l in io.open(FITS, encoding="utf-8-sig")
                             if not l.startswith("#")), delimiter=";"):
        F["fits"].setdefault(r["section"], []).append(r)

    # Разобранный пример: QP_sh, уравнение NI = A·y + L·H(C·y) на рабочих картах.
    Hof = sm._material()
    H = lambda b: float(np.atleast_1d(Hof(b))[0])
    fit = {r["section"]: r for rs in F["fits"].values() for r in rs}
    q = fit["QP_sh"]
    C, L, aon = float(q["C"]), float(q["L_m"]), float(q["A_over_naive"])
    rows = sorted(csv.DictReader(io.open(os.path.join("workspace", "datasets",
                  "currents_to_kicks", "QP_sh.csv"), encoding="utf-8-sig"), delimiter=";"),
                  key=lambda r: float(r["AT_I"]))
    ni0, y0 = float(rows[0]["AT_I"]), abs(float(rows[0]["int_G"]))
    A = aon * ni0 / y0
    picks, seen = [], set()
    for r in rows:
        ni, y = float(r["AT_I"]), abs(float(r["int_G"]))
        band = round(ni / 1600)
        if ni > 6000 or band in seen:
            continue
        seen.add(band)
        air, iron = A * y, L * H(C * y)
        picks.append((ni, y, C * y, H(C * y), air, iron, air + iron,
                      100 * (air + iron - ni) / ni))
    F["example"] = dict(A=A, C=C, L=L, ni0=ni0, y0=y0, rows=picks[:4])
    return F


# ─────────────────────────────── разделы ────────────────────────────────

def ch_task(doc, F, lo):
    doc.add_heading("1. Задача и обозначения", level=1)
    doc.add_paragraph(
        "Модель — расчётная функция одного аргумента: на входе вектор возбуждений обмоток "
        "(в ампер-витках, число витков уже учтено), на выходе вектор кик-коэффициентов магнита.")
    _formula(doc, r"\mathbf{k} = \mathcal{K}(\mathbf{I})", lo)
    _where(doc, lo, [
        (r"\mathbf{I}", "вектор, ампер-витки", "возбуждения обмоток (витки уже учтены)"),
        (r"\mathbf{k}", "вектор", "кик-коэффициенты магнита: нормальные и скошенные, по гармоникам"),
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
    _formula(doc, r"\mathbf{F} = M\mathbf{I}", lo)
    _formula(doc, r"M_{bw} = s_w\,\alpha_{wb}", lo)
    _where(doc, lo, [
        (r"\mathbf{F}", "вектор, ампер-витки", "МДС на каждом ребре"),
        (r"M", "матрица", "матрица намотки «рёбра × обмотки», постоянная"),
        (r"s_w", "±1", "знак намотки обмотки"),
        (r"\alpha_{wb}", "скаляр", "доля МДС обмотки w на ребре b (0, если не охватывает)"),
    ])

    p = doc.add_paragraph()
    p.add_run("Магнитная цепь. ").bold = True
    p.add_run("Напор на ребре, поток ребра и закон сохранения потока в каждой вершине дают "
              "нелинейную систему относительно потенциалов вершин:")
    _formula(doc, r"\mathbf{d} = \mathbf{F} - A^{\mathsf{T}}\mathbf{u}", lo)
    _formula(doc, r"\Phi_b = \begin{cases} P_b\,d_b, & \text{воздух},\\[2pt] "
                  r"S_b\,B(d_b/\ell_b), & \text{железо},\end{cases}", lo)
    _formula(doc, r"A\boldsymbol{\Phi}(\mathbf{u}) = \mathbf{0}", lo)
    _where(doc, lo, [
        (r"\mathbf{d}", "вектор, ампер-витки", "магнитный напор на каждом ребре"),
        (r"A", "матрица, ±1 и 0", "инцидентная матрица графа «вершины × рёбра»"),
        (r"\mathbf{u}", "вектор, ампер-витки", "потенциалы вершин (одна опорная, с нулём)"),
        (r"\boldsymbol{\Phi}", "вектор, вебер", "потоки рёбер — искомое"),
        (r"P_b", "скаляр, генри", "проводимость воздушного ребра, P_b = μ₀·S_b/ℓ_b"),
        (r"S_b,\ \ell_b", "скаляры, кв. м и м", "сечение и длина железного ребра"),
        (r"B(\cdot)", "функция, тесла", "кривая намагничивания материала"),
    ])
    doc.add_paragraph(
        "Система решается методом Ньютона по потенциалам вершин: на каждой итерации по невязке "
        "A·Φ(u) и якобиану из проводимостей ветвей потенциалы поправляются")
    _formula(doc, r"\mathbf{u} \leftarrow \mathbf{u} - J^{-1}\,A\boldsymbol{\Phi}(\mathbf{u})", lo)
    _formula(doc, r"J = -A\,\operatorname{diag}(\mathbf{g})\,A^{\mathsf{T}}", lo)
    _formula(doc, r"g_b = \begin{cases} P_b, & \text{воздух},\\[2pt] "
                  r"(S_b/\ell_b)\,B'(d_b/\ell_b), & \text{железо}.\end{cases}", lo)
    _where(doc, lo, [
        (r"J", "матрица", "якобиан системы по потенциалам"),
        (r"g_b", "скаляр", "проводимость ребра ∂Φ_b/∂d_b: у воздуха постоянна, у железа зависит "
         "от наклона кривой B′"),
    ])
    doc.add_paragraph(
        "Итерации сходятся за единицы шагов; шаг при необходимости дробится пополам, пока невязка "
        "убывает. Решение не зависит от гармоники — потоки одни, гармоники различаются лишь "
        "свёрткой ниже.")

    p = doc.add_paragraph()
    p.add_run("Признаки и насыщение. ").bold = True
    p.add_run("Поле у пучка создают зазоры между наконечниками и, где есть, шунты полюсов. "
              "Признак — поток узла (его компонента в чистом потоке через A); форма поля меняется "
              "с насыщением, мера которого — потоко-взвешенная доля:")
    _formula(doc, r"q_t = (A\boldsymbol{\Phi})_t", lo)
    _formula(doc, r"\sigma = \frac{\sum |\Phi_b|\,\min(|\Phi_b|/S_bB_{\max},\,1)}{\sum |\Phi_b|}", lo)
    _where(doc, lo, [
        (r"\mathbf{q}", "вектор, вебер", "потоки наконечников — признаки кика"),
        (r"\sigma", "скаляр, 0…1", "доля насыщения железа"),
        (r"B_{\max}", "скаляр, тесла", "предел кривой намагничивания"),
    ])

    p = doc.add_paragraph()
    p.add_run("Кик и единицы. ").bold = True
    p.add_run("Кик — свёртка признаков матрицей связи с поправкой на насыщение; мнимая часть даёт "
              "нормальную составляющую, вещественная — скошенную. Итог приводится к "
              "физкоэффициентам решётки:")
    _formula(doc, r"\mathbf{k} = W^{\mathsf{T}}\begin{bmatrix}\mathbf{q}\\ \sigma\mathbf{q}"
                  r"\end{bmatrix}", lo)
    _formula(doc, r"c_n = \varsigma_n\,\frac{(n-1)!}{r^{\,n-1}\,L_{\text{эфф}}\,E}\;\kappa_n", lo)
    _where(doc, lo, [
        (r"W", "матрица", "матрица связи «признаки × компоненты кика», подбираемая"),
        (r"\mathbf{w}_n,\ \mathbf{v}_n", "векторы", "строки W: веса при ненасыщенном железе и "
         "поправка, набирающая силу с σ"),
        (r"c_n", "комплексное", "физкоэффициент гармоники n в конвенции решётки"),
        (r"\varsigma_n", "±1", "знаковая конвенция гармоники"),
        (r"r,\ L_{\text{эфф}}", "скаляры, м", "радиус пробы и эффективная длина элемента"),
        (r"E", "скаляр", "энергетический множитель"),
    ])
    doc.add_paragraph(
        "Матрицы M и A и величины S_b, ℓ_b взяты из чертежей. Подбираются по картам поля только "
        "матрица связи W (линейна, находится точным МНК) и сечения воздушных рёбер утечки; знаки "
        "намотки выбираются из конечного набора. Один прогон функции — это умножение на M, одно "
        "решение цепи Ньютоном, умножение на A, подсчёт σ и две свёртки; итерации есть только в "
        "решении цепи.")


def ch_single(doc, F, lo):
    doc.add_heading("3. Однообмоточный элемент: расчёт на бумаге", level=1)
    doc.add_paragraph(
        "У квадруполей, секступоля и октуполя обмотка одна, поток идёт единственным путём, граф "
        "вырождается в одно ребро, и цепь сводится к скалярному уравнению, связывающему "
        "возбуждение NI с физкоэффициентом y:")
    _formula(doc, r"NI = A\,y + L\,H(C\,y)", lo)
    _where(doc, lo, [
        (r"NI", "скаляр, ампер-витки", "возбуждение обмотки"),
        (r"y", "скаляр", "физкоэффициент элемента — искомое"),
        (r"A", "скаляр", "линейный член: воздух, не насыщается"),
        (r"L", "скаляр, м", "длина пути по железу"),
        (r"C", "скаляр", "пересчёт коэффициента в индукцию железа"),
        (r"H(\cdot)", "функция, А/м", "обратная кривая намагничивания"),
    ])
    doc.add_paragraph(
        "Первое слагаемое — воздух (линеен всегда), второе — железо (почти ничего до насыщения и "
        "резкий рост после). Свободных чисел два, C и L; "
        "множитель A берётся из контрольной точки при малом токе, где железо ещё не мешает "
        "(NI ≈ A·y, откуда A = NI/y). Обе подбираемые величины геометрические — смена материала "
        "означает лишь замену таблицы H(·) без переподгонки.")
    doc.add_paragraph("Подобранные по картам постоянные:")
    rows = []
    for el in ("QP_sh", "QP_L", "SP", "OP"):
        for r in F["fits"].get(el, []):
            rows.append([el, r["coefficient"], f'{float(r["C"]):.4g}', f'{float(r["L_m"]):.4g}',
                         f'{float(r["max_err_pct"]):.1f}', r["n_points"]])
    _table(doc, ["Элемент", "Коэффициент", "C", "L, м", "макс. ошибка, %", "Точек"], rows,
           caption="Постоянные однообмоточного уравнения")

    doc.add_heading("Порядок расчёта", level=2)
    doc.add_paragraph(
        "Дано возбуждение NI, требуется коэффициент y (прямая задача — наоборот, тоже одним "
        "уравнением). Порядок такой:")
    for t in ("взять A, C, L элемента из таблицы 3 (A — через контрольную точку, см. пример);",
              "задаться пробным y и вычислить правую часть A·y + L·H(C·y), взяв H по таблице "
              "материала при индукции C·y;",
              "сравнить с NI и поправить y (уравнение монотонно по y, годится деление отрезка "
              "пополам); двух-трёх приближений достаточно."):
        doc.add_paragraph(t, style="List Number")

    ex = F["example"]
    doc.add_heading("Числовой пример (QP_sh)", level=2)
    doc.add_paragraph(
        f"Постоянные: C = {ex['C']:.4g}, L = {ex['L']:.4g} м. Множитель A по контрольной точке "
        f"NI₀ = {ex['ni0']:.0f} А·в, y₀ = {ex['y0']:.3g} Тл: A = {ex['ni0']:.0f}/{ex['y0']:.3g} "
        f"≈ {ex['A']:.0f}. Подставляя измеренный коэффициент y = int_G в уравнение, получаем "
        f"возбуждение и сверяем с картой:")
    rows = [[f"{ni:.0f}", f"{y:.3f}", f"{cy:.3f}", f"{h:.0f}", f"{air:.0f}", f"{iron:.0f}",
             f"{tot:.0f}", f"{err:+.1f}"]
            for ni, y, cy, h, air, iron, tot, err in ex["rows"]]
    _table(doc, ["NI изм., А·в", "y = int_G, Тл", "C·y, Тл", "H(C·y), А/м",
                 "A·y (возд.)", "L·H (жел.)", "NI расч.", "ошибка, %"], rows,
           caption="Проверка уравнения на рабочих картах QP_sh", pt=8)
    doc.add_paragraph(
        "Видно, как с ростом возбуждения вклад железа (L·H) растёт быстрее воздушного и на верхних "
        "токах становится заметным — это и есть насыщение. Совпадение с картами в пределах "
        "нескольких процентов подтверждает, что расчёт воспроизводим вручную.")


def build(out=OUT):
    import docx
    import latex2omml as lo
    from docx.shared import Pt

    F = facts()
    doc = docx.Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(11)

    doc.add_heading("Модель «токи обмоток → кики»: методика расчёта", level=0)
    p = doc.add_paragraph()
    p.add_run("Уравнения всех звеньев функции ток→кик и разобранный числовой пример для расчёта "
              "вручную. Связь потоков с киком — матрицей связи, без координат полюсов.").italic = True

    for ch in (ch_task, ch_math, ch_single):
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
