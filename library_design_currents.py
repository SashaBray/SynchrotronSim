# -*- coding: utf-8 -*-
"""Проектные токи по КАЖДОМУ элементу библиотеки СИЛА — через МОДЕЛЬ НАСЫЩЕНИЯ (не линейную матрицу связи).

Отдельный файл (не путать с `saturation_model.py`, который проверяет ТОЧНОСТЬ подогнанной модели по
типам элементов, и с `design_currents.py` — тем старым расчётом по ЛИНЕЙНОЙ матрице связи, который
насыщение не учитывает и который просили игнорировать).

Для каждого реального устройства библиотеки (конкретная длина, конкретное место в решётке —
`lattice_configs.csv`, решётка 1): расчётная формула → числа модели именно для этого устройства →
проектный ток из проектного кика (`mat_G`/`mat_S`/`mat_O`/`mat_B`) → обратная задача для проверки.
Корректор — только формула и граф (проектные токи нулевые, см. `saturation_branches.csv`).

Насыщение подогнано ПО ТИПАМ (`saturation_model.py --all`), не по каждой длине — сечение и путь
железа (A, C, L) не зависят от осевой длины магнита (это свойство поперечного сечения), а физкоэф
(G/S/O) — интенсивная величина (на единицу длины), поэтому один и тот же фит обслуживает все длины
одного типа; для каждой длины меняется только целевое значение проектного кика.

DQ (2 обмотки, I1+I2, участвуют в проекте; I3=0) — только ОДНА подогнана насыщением по G (I2 сама по
себе слишком слабая, чтобы качество подгонки было надёжным — обмотки заведомо на разных ветвях, поле
складывается). Решение: взять СООТНОШЕНИЕ токов I1:I2 из старой линейной модели (оно уже отражает,
как обмотки распределяют работу между собой) и подобрать ОБЩИЙ множитель так, чтобы СУММА кика двух
насыщающихся обмоток дала проектный G. Изгиб (B) отдельно моделью насыщения не откалиброван (в
подгонке B берётся из угла КОНКРЕТНОЙ карты, не из независимого измерения) — показан достигнутым по
той же пропорции, БЕЗ отдельной коррекции, это явно оговорено в отчёте.

    py -3.14 library_design_currents.py               # собрать workspace/results/library_design_currents.docx
"""
from __future__ import annotations

import argparse
import json
import os

import numpy as np
from scipy.optimize import brentq

import element_compare as ec
import saturation_model as sm
import design_currents as dc

OUT = os.path.join("workspace", "results", "library_design_currents.docx")
ASSET_DIR = os.path.join("workspace", "results", "element_compare")   # старый отчёт: 3D-сцены + траектории
LATTICE = 1
TARGET_COL = {"QP_sh": "mat_G", "QP_L": "mat_G", "SP": "mat_S", "OP": "mat_O", "DQ": "mat_G"}
TARGET_KEY = {"QP_sh": "G", "QP_L": "G", "SP": "S", "OP": "O", "DQ": "G"}
SOLE_WINDING = {"QP_sh": "I", "QP_L": "I", "SP": "I1", "OP": "I1"}   # у SP реально работает только I1


def _asset(prefix, name, el):
    """Картинка старого отчёта `element_compare` (3D-сцена магнита `shot_`, траектории `mtraj_`/`mxy_`)
    для этого устройства, если она снята. Имена файлов там — БЕЗ суффикса обмотки (кроме DQ, где
    обмотку в имени элемента исторически оставляют — см. `design_currents._base_name`)."""
    base = name if el == "DQ" else dc._base_name(name)
    p = os.path.join(ASSET_DIR, f"{prefix}_{base}.png")
    return p if os.path.exists(p) else None


def _lcrows(lattice=LATTICE):
    return {r["device_type_id"]: r for r in ec.mt._read("lattice_configs.csv")
            if str(r.get("lattice_id")) == str(lattice)}


def _variants(el, lcrow):
    """Устройства решётки этого типа: {имя_варианта → строка lattice_configs}, отсортировано по длине."""
    pref = el + "_"
    out = {n: r for n, r in lcrow.items() if n.startswith(pref) or n == el}
    return dict(sorted(out.items()))


def _fval(row, col):
    s = (row.get(col) or "").strip()
    try:
        v = float(s)
    except ValueError:
        return 0.0
    return 0.0 if abs(v) < 1e-30 else v


def load_fits():
    """Подогнанные (A,C,L) по типам — из датасетов saturation_model.py (--all уже прогнан)."""
    out = {}
    for el in TARGET_COL:
        d = sm.load_dataset(el)
        if not d:
            raise SystemExit(f"нет датасета {el} — сначала py -3.14 saturation_model.py --element {el}")
        out[el] = d
    return out


def ni_of_y(fit, Hof, y):
    """Прямая формула: NI = A·y + L·H(C·y). Проектный ток — ровно эта формула при y=проектный кик."""
    return fit["A"] * y + fit["L"] * Hof(fit["C"] * y)[0]


def y_of_ni(fit, Hof, ni, y_scale):
    """Обратная задача: по NI найти y (числа нет в явном виде — корень NI(y)=ni методом Брента).

    y_scale — масштаб для скобки поиска корня (обычно |y_design|·3, с запасом)."""
    if ni == 0:
        return 0.0
    f = lambda y: ni_of_y(fit, Hof, y) - ni
    lo, hi = sorted((0.0, np.sign(ni) * max(3.0 * abs(y_scale), 1e-9)))
    flo, fhi = f(lo), f(hi)
    tries = 0
    while flo * fhi > 0 and tries < 40:                  # корень не в скобке — раздвигаем
        hi *= 1.7; lo = min(lo, -0.1 * abs(hi)); flo, fhi = f(lo), f(hi); tries += 1
    return brentq(f, lo, hi, xtol=1e-12, rtol=1e-12)


def single_winding_case(el, dev_name, row, fits, Hof):
    """QP_sh/QP_L/OP/SP: одна работающая обмотка, прямая формула + обратная проверка."""
    key = TARGET_KEY[el]
    y_design = _fval(row, TARGET_COL[el])
    if y_design == 0.0:
        return None
    wind = SOLE_WINDING[el]
    fit = fits[el]["fits"].get(wind)
    if not fit or fit["coeff"] != key:
        return None
    turns = sm.BRANCH_GRAPH[el]["turns"].get(wind)
    NI = ni_of_y(fit, Hof, y_design)
    current = NI / turns if turns else None
    y_check = y_of_ni(fit, Hof, NI, y_design)
    return dict(device=dev_name, winding=wind, key=key, y_design=y_design, NI=NI, turns=turns,
               current_a=current, y_check=y_check,
               rel_err=100.0 * (y_check - y_design) / y_design if y_design else 0.0)


def dq_case(dev_name, row, fits, Hof, lin_rows):
    """DQ: I1+I2 по проекту (I3=0). Пропорция I1:I2 — из линейной модели; общий множитель x подбирается
    так, чтобы G1(x·NI1_лин) + G2(x·NI2_лин) = проектный G (поля разных ветвей складываются)."""
    G_design = _fval(row, "mat_G"); B_design = _fval(row, "mat_B")
    if G_design == 0.0:
        return None
    lin = {r["winding"]: r for r in lin_rows if r["element"] == dev_name}
    NI1_lin = float(lin["I1"]["arg"]) * ec.ENERGY_FACTOR if lin.get("I1", {}).get("arg") else 0.0
    NI2_lin = float(lin["I2"]["arg"]) * ec.ENERGY_FACTOR if lin.get("I2", {}).get("arg") else 0.0
    f1, f2 = fits["DQ"]["fits"]["I1"], fits["DQ"]["fits"]["I2"]

    def G_total(x):
        y1 = y_of_ni(f1, Hof, x * NI1_lin, G_design) if NI1_lin else 0.0
        y2 = y_of_ni(f2, Hof, x * NI2_lin, G_design) if NI2_lin else 0.0
        return y1 + y2

    lo, hi = 0.2, 5.0
    glo, ghi = G_total(lo) - G_design, G_total(hi) - G_design
    tries = 0
    while glo * ghi > 0 and tries < 30:
        hi *= 1.5; ghi = G_total(hi) - G_design; tries += 1
    x = brentq(lambda x: G_total(x) - G_design, lo, hi, xtol=1e-10, rtol=1e-10)
    NI1, NI2 = x * NI1_lin, x * NI2_lin
    t1, t2 = sm.BRANCH_GRAPH["DQ"]["turns"]["I1"], sm.BRANCH_GRAPH["DQ"]["turns"]["I2"]
    y1 = y_of_ni(f1, Hof, NI1, G_design) if NI1 else 0.0
    y2 = y_of_ni(f2, Hof, NI2, G_design) if NI2 else 0.0
    return dict(device=dev_name, x=x, NI1=NI1, NI2=NI2, cur1=NI1 / t1, cur2=NI2 / t2,
               y1=y1, y2=y2, G_check=y1 + y2, G_design=G_design, B_design=B_design,
               B_note="B отдельно насыщением не откалиброван — оценка ×масштаб I1, см. текст")


def build(fits, Hof, lo, path=OUT):
    from docx import Document
    import latex2omml as lex
    import saturation_report as sr
    doc = Document()
    sr._tcap.n = sr._fcap.n = 0
    doc.add_heading("Проектные токи элементов библиотеки СИЛА (по модели насыщения)", level=0)
    doc.add_paragraph(
        "Для каждого устройства решётки (своя длина, своё место): расчётная формула модели насыщения, "
        "числа именно для этого устройства, проектный ток из проектного кика (mat_G/mat_S/mat_O — "
        "решётка 1, закрытое кольцо) и обратная задача для проверки. Корректор — без проектных токов "
        "(они нулевые по проекту), только формула и граф ветвей.")
    lcrow = _lcrows()
    dcurrent_rows, _ = dc.compute(LATTICE)

    for el in ("QP_sh", "QP_L", "SP", "OP", "DQ", "corrector"):
        g = sm.BRANCH_GRAPH.get(el, {})
        doc.add_heading(el, level=1)
        fig = sm._branch_figure(el, g, os.path.join(sm.OUTDIR, f"graph_lib_{el}.png"))
        if fig:
            sr._pic(doc, fig, cm=15)
            sr._fcap(doc, f"{el} — граф магнитной цепи: рёбра между узлами = ветви магнитопровода "
                          f"(тонкие серые — ветви без обмоток)")
        if g.get("coils"):
            sr._table(doc, ["Обмотка", "Ребро (ветвь) графа", "витков/катушку"],
                      [[w, " + ".join(g["coils"][w]), str(g["turns"].get(w, "—"))] for w in sorted(g["coils"])],
                      caption=f"{el} — какая обмотка на каком ребре графа (матрица T)")
        if el == "corrector":
            shot = _asset("shot", "corrector_100A_L100", el)
            if shot:
                sr._pic(doc, shot, cm=14)
                sr._fcap(doc, f"{el} — трёхмерная сцена магнита")
        doc.add_heading("Расчётная формула", level=2)
        if el == "DQ":
            sr._formula(doc, r"N I = A\,y + L\,H(C\,y),\quad y_{\Sigma} = y_1(NI_1) + y_2(NI_2)", lex)
            sr._symtab(doc, [
                ["y_Σ", "суммарный физкоэффициент элемента", "скаляр", "Тл/м"],
                ["y₁, y₂", "вклад обмоток I1, I2 (разные ветви — вклады складываются)", "скаляры", "Тл/м"],
                ["NI₁, NI₂", "ампер-витки обмоток I1, I2", "скаляры", "А"],
                ["A", "линейный коэффициент цепи (свой у каждой обмотки)", "скаляр", "А/[y]"],
                ["L", "длина пути по железу", "скаляр", "м"],
                ["C", "y → индукция сердечника", "скаляр", "Тл/[y]"],
                ["H(·)", "обратная кривая намагничивания", "функция", "А/м"],
            ], caption="обозначения")
            doc.add_paragraph("Ветви I1 и I2 разные — их вклад в общий градиент складывается. Проектный "
                               "ток решается через пропорцию токов из линейной модели связи и общий "
                               "множитель на насыщение (подробности — во врезке ниже).")
        elif el == "corrector":
            sr._formula(doc, r"\Phi = Q\left(T\,\mathbf{I} - \Lambda\,H(\Sigma^{-1}\Phi)\right)", lex)
            sr._symtab(doc, [
                ["Φ", "поток ветви", "вектор", "Вб"],
                ["Q", "проводимости ветвей", "матрица", "Вб/А"],
                ["T", "размещение обмоток по рёбрам", "матрица", "—"],
                ["I", "токи обмоток", "вектор", "А"],
                ["Λ", "длина пути по железу", "диаг. матрица", "м"],
                ["Σ", "сечение сердечника", "диаг. матрица", "м²"],
                ["H(·)", "обратная кривая намагничивания", "функция", "А/м"],
            ], caption="обозначения")
            doc.add_paragraph("Проектные токи корректора по проекту нулевые (корректор — элемент "
                               "коррекции орбиты, не элемент структуры кольца) — расчёт токов не "
                               "выполняется, приведена только модель.")
        else:
            sr._formula(doc, r"N I = A\,y + L\,H(C\,y)", lex)
            sr._symtab(doc, [
                ["NI", "ампер-витки обмотки", "скаляр", "А"],
                ["y", "физкоэффициент (G, S или O)", "скаляр", "Тл/м^(n−1)"],
                ["A", "линейный коэффициент цепи", "скаляр", "А/[y]"],
                ["L", "длина пути по железу", "скаляр", "м"],
                ["C", "y → индукция сердечника", "скаляр", "Тл/[y]"],
                ["H(·)", "обратная кривая намагничивания", "функция", "А/м"],
            ], caption="обозначения")
        if el in fits and fits[el]["fits"]:
            wtab = [[w, f["coeff"], f"{f['A']:.5g}", f"{f['C']:.5g}", f"{1000*f['L']:.1f}"]
                    for w, f in sorted(fits[el]["fits"].items())]
            sr._table(doc, ["Обмотка", "Коэфф.", "A", "C", "L, мм"], wtab,
                      caption=f"{el} — числа модели (по типу, общие для всех длин)")
        if el == "corrector":
            continue
        doc.add_heading("Проектные токи устройств", level=2)
        variants = _variants(el, lcrow)
        any_row = False
        if el == "DQ":
            for name, row in variants.items():
                if not name.startswith("DQ_I1"):          # проектный кик хранится на строке I1
                    continue
                c = dq_case(name, row, fits, Hof, dcurrent_rows)
                if not c:
                    continue
                any_row = True
                doc.add_heading(name, level=3)
                shot = _asset("shot", name, el)
                if shot:
                    sr._pic(doc, shot, cm=13)
                    sr._fcap(doc, f"{name} — трёхмерная сцена магнита")
                sr._table(doc, ["проект. G", "проект. B", "x", "NI₁", "NI₂",
                                "ток I1, А", "ток I2, А", "проверка G(I1)+G(I2)"],
                          [[f"{c['G_design']:+.4e}", f"{c['B_design']:+.4e}", f"{c['x']:.4f}",
                            f"{c['NI1']:.4g}", f"{c['NI2']:.4g}", f"{c['cur1']:.3f}", f"{c['cur2']:.3f}",
                            f"{c['G_check']:+.4e}"]],
                          caption=f"{name} — проектные токи (модель насыщения; B — по проекту, "
                                  f"насыщением отдельно не корректировался)")
                mtraj, mxy = _asset("mtraj", name, el), _asset("mxy", name, el)
                if mtraj:
                    sr._pic(doc, mtraj, cm=15)
                    sr._fcap(doc, f"{name} — форма траектории (иллюстрация из характеризации элемента; "
                                  f"снята на своём наборе токов, не на проектном)")
                if mxy:
                    sr._pic(doc, mxy, cm=13)
                    sr._fcap(doc, f"{name} — сечение пучка в плоскости x-y до/после элемента")
        else:
            for name, row in variants.items():
                c = single_winding_case(el, name, row, fits, Hof)
                if not c:
                    continue
                any_row = True
                doc.add_heading(name, level=3)
                shot = _asset("shot", name, el)
                if shot:
                    sr._pic(doc, shot, cm=13)
                    sr._fcap(doc, f"{name} — трёхмерная сцена магнита")
                sr._table(doc, ["обмотка", f"проект. {c['key']}", "NI, А", "витк./кат.",
                                "проектный ток, А", f"обратно {c['key']}", "расхожд."],
                          [[c["winding"], f"{c['y_design']:+.5e}", f"{c['NI']:.5g}",
                            str(c["turns"] or "—"),
                            f"{c['current_a']:.4f}" if c["current_a"] is not None else "—",
                            f"{c['y_check']:+.5e}", f"{c['rel_err']:+.4f}%"]],
                          caption=f"{name} — проектный ток (прямая формула → обратная проверка)")
                mtraj, mxy = _asset("mtraj", name, el), _asset("mxy", name, el)
                if mtraj:
                    sr._pic(doc, mtraj, cm=15)
                    sr._fcap(doc, f"{name} — форма траектории (иллюстрация из характеризации элемента; "
                                  f"снята на своём наборе токов, не на проектном)")
                if mxy:
                    sr._pic(doc, mxy, cm=13)
                    sr._fcap(doc, f"{name} — сечение пучка в плоскости x-y до/после элемента")
        if not any_row:
            doc.add_paragraph("В решётке 1 нет устройств этого типа с ненулевым проектным киком.")

    doc.add_heading("Кривая намагничивания", level=1)
    F = sr._figures()
    sr._pic(doc, F["bh"], cm=16)
    sr._fcap(doc, "Кривая намагничивания материала сердечника")
    import pandas as pd
    bh = pd.read_csv(sm.MATERIAL, sep=";")
    sr._table(doc, ["H, А/м", "B, Тл"], [[f"{r[0]:g}", f"{r[1]:g}"] for r in bh.to_numpy()],
              caption="Точки таблицы материала")
    sr._body_format(doc)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    doc.save(path)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    fits = load_fits()
    Hof = sm._material()
    lo = None
    print("отчёт ->", build(fits, Hof, lo, a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
