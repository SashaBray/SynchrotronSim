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
    return F


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
    doc.add_heading("3. Однообмоточный элемент: расчёт на бумаге", level=1)
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
        "Возбуждение связано со СТРОГО интегральной силой гармоники χ (интеграл коэффициента вдоль "
        "пучка) скалярным уравнением магнитной цепи:")
    _formula(doc, r"NI = a\,\chi + l\,H(c\,\chi)", lo)
    _where(doc, lo, [
        (r"\chi", "скаляр", "интегральная сила гармоники, χ = ∫G·dl (или ∫S·dl, ∫O·dl)"),
        (r"a", "скаляр", "линейный член: воздух, не насыщается"),
        (r"l", "скаляр, м", "длина пути по железу"),
        (r"c", "скаляр", "пересчёт χ в индукцию железа"),
        (r"H(\cdot)", "функция, А/м", "обратная кривая намагничивания"),
    ])
    doc.add_paragraph(
        "Первое слагаемое — воздух (линеен всегда), второе — железо (почти ничего до насыщения и "
        "резкий рост после). Наконец, при ФИКСИРОВАННОЙ эффективной длине физкоэффициент — это "
        "интеграл, делённый на неё:")
    _formula(doc, r"G = \chi / l_{\text{эфф}}", lo)
    _where(doc, lo, [
        (r"G", "скаляр, Тл/м", "физкоэффициент (для секступоля S, Тл/м²; октуполя O, Тл/м³)"),
        (r"l_{\text{эфф}}", "скаляр, м", "эффективная длина элемента (проектная, фиксированная)"),
    ])
    doc.add_paragraph(
        "Свободных в уравнении цепи два числа, c и l; множитель a берётся из контрольной точки при "
        "малом токе, где железо ещё не мешает (NI ≈ a·χ, откуда a = NI/χ). Величины c и l "
        "геометрические — смена материала означает лишь замену таблицы H(·) без перенастройки. "
        "Модель везде останавливается на физкоэффициенте (Тл/м, Тл/м², Тл/м³); перевод в "
        "нормированную силу решётки — деление на жёсткость пучка — тут не рассматривается.")
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
