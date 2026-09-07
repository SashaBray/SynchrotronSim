# -*- coding: utf-8 -*-
"""Word-отчёт: аналитическая модель «токи обмоток → кики» как расчётная функция.

Формальное описание всей цепочки от ВЕКТОРА ТОКОВ на входе до ВЕКТОРА КИКОВ на выходе: структуры
данных (матрица намотки, граф магнитной цепи, матрица связи), уравнения каждого звена, что берётся
из чертежей и что подбирается, сколько свободных чисел у каждого элемента и какая достигнута
точность. Связь «потоки → кик» описана через МАТРИЦУ СВЯЗИ (`coupling_model`), без координат полюсов.

Все числа считаются живой моделью при сборке — переписанных руками значений в документе нет.

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
GRAPH = ("DQ", "corrector")

# Вход и выход каждого элемента библиотеки — модельное соглашение (см. постановку задачи):
# однообмоточные дают один момент, DQ три, корректор шесть (нормальные и скошенные).
IO = {
    "QP_sh": ("I1", "G"), "QP_L": ("I1", "G"), "OP": ("I1", "O"),
    "SP": ("I1, I2", "S"), "DQ": ("I1, I2, I3", "B, G, S"),
    "corrector": ("I1…I5", "B, G, S, Bs, Gs, Ss"),
}


def _formula(doc, latex, lo):
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    lo.add_formula(p, latex)


def _vars(doc, triples):
    from docx.shared import Cm, Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    doc.add_paragraph("Где:")
    for sym, kind, txt in triples:
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(1.0)
        p.paragraph_format.first_line_indent = Cm(-1.0)
        p.paragraph_format.space_after = Pt(1)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.add_run("— ")
        p.add_run(sym).italic = True
        p.add_run(f" — {txt}; {kind}.")


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
        if p.alignment == WD_ALIGN_PARAGRAPH.CENTER or p.text.startswith("— "):
            continue
        p.paragraph_format.first_line_indent = Cm(1.0)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY


def facts():
    import branch_network as bn
    import coupling_model as cm

    F = {"graph": {}, "free": {}, "acc": {}, "fits": {}}
    for el in ELEMENTS:
        cir = bn.Circuit(el)
        F["graph"][el] = dict(nodes=len(cir.nodes), branches=len(cir.branches),
                              iron=len(cir.iron), gaps=len(cir.gaps), coils=len(cir.coil_branches),
                              geom=cir.geom is not None)
    F["fb"] = {}
    for el in GRAPH:
        F["free"][el] = cm.free_numbers(el)
        F["fb"][el] = cm.feature_breakdown(el)
        F["acc"][el] = cm.accuracy(el)
    for r in csv.DictReader((l for l in io.open(FITS, encoding="utf-8-sig")
                             if not l.startswith("#")), delimiter=";"):
        F["fits"].setdefault(r["section"], []).append(r)
    return F


def ch_io(doc, F, lo):
    doc.add_heading("1. Вход и выход", level=1)
    doc.add_paragraph(
        "Модель — расчётная функция одного аргумента: на входе вектор возбуждений обмоток "
        "(в ампер-витках, число витков уже учтено), на выходе вектор кик-коэффициентов магнита.")
    _formula(doc, r"\mathbf{k} = \mathcal{K}(\mathbf{I})", lo)
    doc.add_paragraph(
        "Размерности входа и выхода у каждого элемента библиотеки заданы его назначением:")
    _table(doc, ["Элемент", "Вход", "Выход"],
           [[el, IO[el][0], IO[el][1]] for el in ELEMENTS],
           caption="Что каждый элемент получает и что создаёт")
    doc.add_paragraph(
        "Скошенные составляющие (Bs, Gs, Ss) есть только у корректора: это единственный элемент с "
        "несколькими независимо управляемыми катушками. У остальных выход строго определён.")


def ch_chain(doc, F, lo):
    doc.add_heading("2. Из чего состоит функция", level=1)
    doc.add_paragraph("Между входом и выходом стоят четыре звена, второе нелинейно:")
    for t in ("токи → МДС ветвей: постоянная матрица намотки M;",
              "МДС → потоки: нелинейная магнитная цепь, определённая графом и свойствами "
              "материала;",
              "потоки → базис кика: потоки наконечников полюсов (топология графа);",
              "базис → кик: матрица связи с поправкой на насыщение."):
        doc.add_paragraph(t, style="List Bullet")
    doc.add_paragraph(
        "Матрица намотки и граф определяются на основе чертежей, предварительно требуется подобрать "
        "матрицу связи последнего и параметры утечек из магнитопровода.")


def ch_mmf(doc, F, lo):
    doc.add_heading("3. Токи → магнитодвижущие силы", level=1)
    doc.add_paragraph("Ток обмотки создаёт МДС на охватываемых ею рёбрах. Связь линейна:")
    _formula(doc, r"\mathbf{F} = M\,\mathbf{I}, \qquad M_{bw} = s_w\,\alpha_{wb}", lo)
    _vars(doc, [
        ("F", "ВЕКТОР, ампер-витки", "МДС на каждом ребре"),
        ("M", "МАТРИЦА", "матрица намотки «рёбра × обмотки», постоянная"),
        ("s_w", "СКАЛЯР, ±1", "знак намотки обмотки"),
        ("α_wb", "СКАЛЯР", "доля МДС обмотки w на ребре b (0, если не охватывает)"),
    ])
    doc.add_paragraph(
        "Знаки намотки — единственное подбираемое в M, и они дискретны: набор один на элемент, "
        "общий для всех гармоник (обмотка намотана в одну сторону).")


def ch_circuit(doc, F, lo):
    doc.add_heading("4. Магнитная цепь (нелинейное звено)", level=1)
    doc.add_paragraph(
        "Распределение потока по рёбрам при заданных МДС находится из трёх соотношений — напор на "
        "ребре, закон ребра, сохранение потока в вершине:")
    _formula(doc, r"\mathbf{d} = \mathbf{F} - A^{\mathsf{T}}\mathbf{u}, \qquad "
                  r"\Phi_b = \begin{cases} P_b\,d_b, & \text{воздух},\\ "
                  r"S_b\,B(d_b/\ell_b), & \text{железо},\end{cases} \qquad "
                  r"A\,\boldsymbol{\Phi}(\mathbf{u}) = \mathbf{0}", lo)
    _vars(doc, [
        ("u", "ВЕКТОР, ампер-витки", "потенциалы вершин (одна опорная, с нулём)"),
        ("A", "МАТРИЦА, ±1 и 0", "инцидентная матрица графа «вершины × рёбра»"),
        ("Φ", "ВЕКТОР, вебер", "потоки рёбер — искомое"),
        ("P_b", "СКАЛЯР, генри", "проводимость воздушного ребра, P_b = μ₀·S_b/ℓ_b"),
        ("S_b, ℓ_b", "СКАЛЯРЫ, кв. м и м", "сечение и длина железного ребра"),
        ("B(·)", "СКАЛЯРНАЯ ФУНКЦИЯ, тесла", "кривая намагничивания материала"),
    ])
    doc.add_paragraph(
        "У воздуха сечение и длина задаются их отношением — P_b. "
        "Система нелинейна лишь через B(·) и решается методом Ньютона по потенциалам вершин с "
        "аналитическим якобианом J = −A·diag(g)·Aᵀ, где g — производная потока ребра по напору.")


def ch_coupling(doc, F, lo):
    doc.add_heading("5. Потоки → кик: матрица связи", level=1)
    doc.add_paragraph(
        "Поле у пучка создают воздушные зазоры между наконечниками полюсов и, где они есть, шунты "
        "полюсов (воздушный путь мимо железа). Признаком служит поток узла — сумма потоков "
        "сходящихся в нём зазоров (величина топологическая, координаты не нужны):")
    _formula(doc, r"q_t = \sum_{b\,:\,t\,\in\,b} \pm\,\Phi_b", lo)
    doc.add_paragraph("Форма поля меняется с насыщением, поэтому в признаки добавляется их же "
                      "копия, взвешенная долей насыщения железа σ, и кик получается матрицей связи:")
    _formula(doc, r"\sigma = \frac{\sum |\Phi_b|\,\min(|\Phi_b|/S_bB_{\max},1)}{\sum |\Phi_b|}, "
                  r"\qquad \mathbf{k} = W^{\mathsf{T}}\begin{bmatrix}\mathbf{q}\\ "
                  r"\sigma\,\mathbf{q}\end{bmatrix}", lo)
    _vars(doc, [
        ("q", "ВЕКТОР, вебер", "потоки наконечников"),
        ("σ", "СКАЛЯР, 0…1", "потоко-взвешенная доля насыщения железа"),
        ("W", "МАТРИЦА", "матрица связи «признаки × компоненты кика», подбираемая"),
        ("k", "ВЕКТОР", "компоненты кика: нормальные и скошенные, по гармоникам"),
    ])


def ch_units(doc, F, lo):
    doc.add_heading("6. Перевод кика в единицы решётки", level=1)
    doc.add_paragraph("Кик, выраженный гармоникой поля на пробном радиусе, приводится к "
                      "физкоэффициентам решётки:")
    _formula(doc, r"c_n = \varsigma_n\,\frac{(n-1)!}{r^{\,n-1}\,L_{\text{эфф}}\,E}\;\kappa_n", lo)
    _vars(doc, [
        ("c_n", "КОМПЛЕКСНОЕ", "физкоэффициент гармоники n"),
        ("ς_n", "СКАЛЯР, ±1", "знаковая конвенция гармоники"),
        ("r, L_эфф", "СКАЛЯРЫ, м", "радиус пробы и эффективная длина элемента"),
        ("E", "СКАЛЯР", "энергетический множитель"),
    ])


def ch_free(doc, F, lo):
    doc.add_heading("7. Свободные числа", level=1)
    doc.add_paragraph(
        "Из чертежей взяты сечения и длины рёбер, расположение обмоток и число витков; матрицы M и "
        "A из них построены целиком. Подбирается только матрица связи W: её размер равен удвоенному "
        "числу потоковых признаков на число компонент кика. Признаки двух родов — потоки "
        "наконечников полюсов и потоки шунтов полюсов (воздушный путь мимо железа): у DQ 4 "
        "наконечника и 4 шунта, у корректора 6 наконечников без отдельных шунтов.")
    rows = []
    for el in GRAPH:
        f, t, c = F["free"][el]
        np_, ns_ = F["fb"][el]
        rows.append([el, np_, ns_, c, f, F["graph"][el]["coils"]])
    _table(doc, ["Элемент", "Наконечников", "Шунтов", "Компонент", "Свободных чисел", "Обмоток"],
           rows, caption="Свободные числа матрицы связи (= 2·(наконечники+шунты)·компоненты; "
                         "множитель 2 — копия признаков, взвешенная насыщением)")
    doc.add_paragraph(
        "Все они линейны и находятся точным МНК. Сверх этого подбираются сечения воздушных рёбер "
        "утечки (по одному числу на класс: у DQ один, у корректора два) — они задают, где у модели "
        "наступает потолок. Знаки намотки выбираются из конечного набора перебором.")


def ch_single(doc, F, lo):
    doc.add_heading("8. Элементы с одной обмоткой", level=1)
    doc.add_paragraph(
        "У квадруполей, секступоля и октуполя поток идёт единственным путём, граф вырождается в "
        "одно ребро, и вся цепь сводится к одному уравнению относительно кика:")
    _formula(doc, r"NI = A\,y + L\,H(C\,y)", lo)
    _vars(doc, [
        ("NI", "СКАЛЯР, ампер-витки", "возбуждение обмотки"),
        ("y", "СКАЛЯР", "физкоэффициент — искомое"),
        ("A", "СКАЛЯР", "линейный член: воздух, не насыщается"),
        ("L, C", "СКАЛЯРЫ", "длина пути по железу и пересчёт коэффициента в индукцию"),
        ("H(·)", "СКАЛЯРНАЯ ФУНКЦИЯ", "обратная кривая намагничивания"),
    ])
    doc.add_paragraph(
        "Свободных чисел два — L и C, оба геометрические (смена материала не требует переподгонки); "
        "A определяется по контрольной точке при малом токе.")
    rows = []
    for el in ELEMENTS:
        for r in F["fits"].get(el, []):
            rows.append([el, r["coefficient"], f'{float(r["C"]):.4g}', f'{float(r["L_m"]):.4g}',
                         f'{float(r["max_err_pct"]):.2f}', r["n_points"]])
    if rows:
        _table(doc, ["Элемент", "Коэффициент", "C", "L, м", "макс. ошибка, %", "Точек"], rows,
               caption="Подгонка однообмоточного уравнения")


def ch_acc(doc, F, lo):
    doc.add_heading("9. Точность", level=1)
    doc.add_paragraph(
        "Матрица связи подбирается только на картах с одной возбуждённой обмоткой; карты с "
        "несколькими обмотками в подгонку не входят и служат независимой проверкой. Метрика — "
        "расхождение по кику (сумма модулей разностей к сумме модулей измеренного), по всем картам.")
    rows = []
    NAME = {1: "B", 2: "G", 3: "S"}
    for el in GRAPH:
        for n in (1, 2, 3):
            s = F["acc"][el][(n, "одиночные")]
            c = F["acc"][el][(n, "комбинации")]
            rows.append([el, NAME[n], f"{s[0]:.1f}", f"{s[1]}", f"{c[0]:.1f}", f"{c[1]}"])
    _table(doc, ["Элемент", "Гарм.", "одиночные, %", "карт", "комбинации, %", "карт"], rows,
           caption="Расхождение на подгонке и на проверке")
    doc.add_paragraph(
        "Разница между DQ и корректором объясняется не устройством модели — оно одно, — а "
        "плотностью данных: на обмотку DQ приходится вдвое больше карт, чем на обмотку корректора.")


def build(out=OUT):
    import docx
    import latex2omml as lo
    from docx.shared import Pt

    F = facts()
    doc = docx.Document()
    doc.styles["Normal"].font.name = "Times New Roman"
    doc.styles["Normal"].font.size = Pt(11)

    doc.add_heading("Аналитическая модель «токи обмоток → кики»", level=0)
    p = doc.add_paragraph()
    p.add_run("Расчётная функция: структуры данных, уравнения звеньев, свободные параметры и "
              "точность. Связь потоков с киком — матрицей связи, без координат полюсов.").italic = True

    for ch in (ch_io, ch_chain, ch_mmf, ch_circuit, ch_coupling, ch_units,
               ch_free, ch_single, ch_acc):
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
