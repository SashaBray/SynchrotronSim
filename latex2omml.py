# -*- coding: utf-8 -*-
"""LaTeX → OMML — НАТИВНЫЕ формулы Word (не картинки).

Пайплайн: `latex2mathml` (LaTeX → MathML presentation) → Office `MML2OMML.XSL` (MathML → OMML) через lxml.
Вставка в python-docx: `add_formula(paragraph, latex)` дописывает `<m:oMath>` в конец абзаца — Word
отрисует настоящую редактируемую формулу.

Сам документ собирается чистым Python (`python-docx`), поэтому переносим. Единственная привязка к
Windows — XSL из Microsoft Office: файл её частью и остаётся, скопировать нельзя. Без него формулы
вставляются ТЕКСТОМ (курсивом), а отчёт собирается целиком: потерять оформление формул лучше, чем
потерять весь документ. Путь к XSL переопределяется env `MML2OMML_XSL` — если файл скопирован с
Windows-машины, нативные формулы работают и на Linux.
"""
import os
import sys
from lxml import etree
import latex2mathml.converter

# XSL Office (MathML→OMML). Путь можно переопределить env MML2OMML_XSL.
_XSL = os.environ.get("MML2OMML_XSL") or r"C:\Program Files\Microsoft Office\root\Office16\MML2OMML.XSL"
_transform = None
_warned = False


def available() -> bool:
    """Доступны ли НАТИВНЫЕ формулы (есть ли XSL Office)."""
    return os.path.exists(_XSL)


def _warn_once() -> None:
    global _warned
    if not _warned:
        _warned = True
        print(f"[формулы] {_XSL} не найден — формулы вставляются текстом. "
              f"Нативные формулы Word: задайте env MML2OMML_XSL на копию MML2OMML.XSL.",
              file=sys.stderr)


def _get_transform():
    global _transform
    if _transform is None:
        if not os.path.exists(_XSL):
            raise FileNotFoundError(f"MML2OMML.XSL не найден: {_XSL} (задай env MML2OMML_XSL)")
        _transform = etree.XSLT(etree.parse(_XSL))
    return _transform


M_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/math}"


def _fill_empty_slots(root):
    """Заткнуть ПУСТЫЕ слоты OMML — Word рисует на их месте пунктирный квадратик.

    Так получается у операторов с пределами (∑, ∏, ∫, ∮): XSL кладёт знак с пределами в <m:nary>, а
    подынтегральное выражение оставляет соседними ранами, и слот операнда <m:e/> остаётся пустым.
    Правильнее писать операнд в фигурных скобках (sum_i{f_i}) — тогда слот заполняется сам; это
    защита на случай, если так не написали."""
    for e in root.iter(M_NS + "e"):
        if len(e) == 0 and not (e.text or "").strip():
            r = etree.SubElement(e, M_NS + "r")
            etree.SubElement(r, M_NS + "t").text = "​"        # нулевой ширины: слот не пуст, но не виден
    return root


def latex_to_omml(latex):
    """LaTeX-строка → lxml-элемент OMML (`<m:oMath>`) для вставки в docx."""
    mathml = latex2mathml.converter.convert(latex)
    omml_root = _get_transform()(etree.fromstring(mathml)).getroot()
    # XSL отдаёт <m:oMathPara>…<m:oMath>…; для инлайна берём внутренний <m:oMath>, иначе корень
    tag = M_NS + "oMath"
    return _fill_empty_slots(omml_root.find(tag) if omml_root.tag != tag else omml_root)


def set_size(omml, size_pt):
    """Задать кегль ранам формулы (расшифровки под формулой идут мельче основного текста)."""
    W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    half = str(int(round(size_pt * 2)))
    for r in omml.iter(M_NS + "r"):
        pr = r.find(W + "rPr")
        if pr is None:
            pr = etree.SubElement(r, W + "rPr")
            r.insert(0, pr)
        for tag in ("sz", "szCs"):
            el = pr.find(W + tag)
            if el is None:
                el = etree.SubElement(pr, W + tag)
            el.set(W + "val", half)
    return omml


def add_formula(paragraph, latex, size_pt=None):
    """Вставить формулу LaTeX в конец абзаца python-docx. Возврат paragraph.

    Нет XSL Office (обычная ситуация вне Windows) ⇒ формула идёт курсивным текстом, сборка отчёта
    не прерывается."""
    try:
        node = latex_to_omml(latex)
    except (FileNotFoundError, OSError, etree.XMLSyntaxError, etree.XSLTError):
        _warn_once()
        r = paragraph.add_run(latex); r.italic = True
        if size_pt:
            from docx.shared import Pt
            r.font.size = Pt(size_pt)
        return paragraph
    if node is not None:
        if size_pt:
            set_size(node, size_pt)
        paragraph._p.append(node)
    return paragraph
