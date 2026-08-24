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


def latex_to_omml(latex):
    """LaTeX-строка → lxml-элемент OMML (`<m:oMath>`) для вставки в docx."""
    mathml = latex2mathml.converter.convert(latex)
    omml_root = _get_transform()(etree.fromstring(mathml)).getroot()
    # XSL отдаёт <m:oMathPara>…<m:oMath>…; для инлайна берём внутренний <m:oMath>, иначе корень
    tag = "{http://schemas.openxmlformats.org/officeDocument/2006/math}oMath"
    return omml_root.find(tag) if omml_root.tag != tag else omml_root


def add_formula(paragraph, latex):
    """Вставить формулу LaTeX в конец абзаца python-docx. Возврат paragraph.

    Нет XSL Office (обычная ситуация вне Windows) ⇒ формула идёт курсивным текстом, сборка отчёта
    не прерывается."""
    try:
        node = latex_to_omml(latex)
    except (FileNotFoundError, OSError, etree.XMLSyntaxError, etree.XSLTError):
        _warn_once()
        paragraph.add_run(latex).italic = True
        return paragraph
    if node is not None:
        paragraph._p.append(node)
    return paragraph
