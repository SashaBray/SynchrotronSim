# -*- coding: utf-8 -*-
"""LaTeX → OMML — НАТИВНЫЕ формулы Word (не картинки).

Пайплайн: `latex2mathml` (LaTeX → MathML presentation) → Office `MML2OMML.XSL` (MathML → OMML) через lxml.
Вставка в python-docx: `add_formula(paragraph, latex)` дописывает `<m:oMath>` в конец абзаца — Word
отрисует настоящую редактируемую формулу.
"""
import os
from lxml import etree
import latex2mathml.converter

# XSL Office (MathML→OMML). Путь можно переопределить env MML2OMML_XSL.
_XSL = os.environ.get("MML2OMML_XSL") or r"C:\Program Files\Microsoft Office\root\Office16\MML2OMML.XSL"
_transform = None


def _get_transform():
    global _transform
    if _transform is None:
        if not os.path.exists(_XSL):
            raise FileNotFoundError(f"MML2OMML.XSL не найден: {_XSL} (задай env MML2OMML_XSL)")
        _transform = etree.XSLT(etree.parse(_XSL))
    return _transform


def latex_to_omml(latex):
    """LaTeX-строка → lxml-элемент OMML (`<m:oMath>`) для вставки в docx."""
    mathml = latex2mathml.converter.convert(latex)
    omml_root = _get_transform()(etree.fromstring(mathml)).getroot()
    # XSL отдаёт <m:oMathPara>…<m:oMath>…; для инлайна берём внутренний <m:oMath>, иначе корень
    tag = "{http://schemas.openxmlformats.org/officeDocument/2006/math}oMath"
    return omml_root.find(tag) if omml_root.tag != tag else omml_root


def add_formula(paragraph, latex):
    """Вставить формулу LaTeX как нативную формулу Word в конец абзаца python-docx. Возврат paragraph."""
    node = latex_to_omml(latex)
    if node is not None:
        paragraph._p.append(node)
    return paragraph
