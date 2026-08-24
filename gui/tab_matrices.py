# -*- coding: utf-8 -*-
"""Вкладка «Матрицы связи»: измерение матриц ток→физкоэф, Word-отчёт и проектные токи обмоток.

Два разных по цене действия разведены явно:
  • измерение матриц (`element_compare.py matrices` / `library`) гоняет диффур по картам поля —
    это десятки минут;
  • расчёт проектных токов (`design_currents.py`) — линейная алгебра над уже измеренными матрицами
    из базы, мгновенно.
Поэтому у второй кнопки показывается возраст таблицы матриц: считать токи по устаревшим матрицам
можно, но результат будет описывать прежнюю машину.
"""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QGroupBox, QHBoxLayout, QLabel, QPushButton,
                               QVBoxLayout)

from . import paths
from .csvtable import CsvEditor
from .runner import py_job
from .tabbase import TabBase
from .widgets import SectionLabel, open_in_os

REPORT = paths.RESULTS / "element_compare" / "element_compare_report.docx"


class MatricesTab(TabBase):
    def __init__(self, runner, log, parent=None):
        super().__init__(runner, log, parent)

        # --- измерение матриц ------------------------------------------------
        self.btn_matrices = QPushButton("Матрицы + отчёт (matrices)")
        self.btn_matrices.setToolTip("element_compare.py matrices — измерить матрицы связи всех элементов\n"
                                     "и собрать Word-отчёт. Таблицу базы и HTML не трогает.")
        self.btn_library = QPushButton("Полный конвейер (library)")
        self.btn_library.setToolTip("element_compare.py library — матрицы + запись в таблицу базы\n"
                                    "coupling_matrices.csv + отчёт + демо-конфигурации + HTML.")
        self.btn_one = QPushButton("Один элемент")
        self.btn_one.setToolTip("element_compare.py library <элемент> — пересчитать один магнит\n"
                                "и обновить только его запись в таблице базы.")
        self.btn_compare = QPushButton("Сравнение моделей (all)")
        self.btn_compare.setToolTip("element_compare.py all — отчёт сопоставления киковой и диффурной\n"
                                    "моделей по элементам (свип токов, невязки траекторий).")
        self.elem = QComboBox()
        self.elem.setMinimumWidth(220)

        self.ck_report = QCheckBox("отчёт")
        self.ck_demo = QCheckBox("демо-конфигурации")
        self.ck_html = QCheckBox("HTML-отчёты элементов")
        for c in (self.ck_report, self.ck_demo):
            c.setChecked(True)
        self.ck_html.setChecked(False)
        self.ck_html.setToolTip("HTML-отчёты элементов заметно дольше остального конвейера")

        self.btn_matrices.clicked.connect(self._run_matrices)
        self.btn_library.clicked.connect(self._run_library)
        self.btn_one.clicked.connect(self._run_one)
        self.btn_compare.clicked.connect(self._run_all)

        box1 = QGroupBox("Измерение матриц связи (диффур по картам поля — долго)")
        v1 = QVBoxLayout(box1)
        r1 = QHBoxLayout()
        r1.addWidget(self.btn_matrices)
        r1.addWidget(self.btn_compare)
        r1.addStretch(1)
        r2 = QHBoxLayout()
        r2.addWidget(self.btn_library)
        for c in (self.ck_report, self.ck_demo, self.ck_html):
            r2.addWidget(c)
        r2.addStretch(1)
        r3 = QHBoxLayout()
        r3.addWidget(QLabel("Элемент:"))
        r3.addWidget(self.elem)
        r3.addWidget(self.btn_one)
        r3.addStretch(1)
        for r in (r1, r2, r3):
            v1.addLayout(r)

        # --- проектные токи ---------------------------------------------------
        self.btn_currents = QPushButton("Рассчитать проектные токи по матрицам")
        self.btn_currents.setToolTip("design_currents.py — обратная задача: проектные mat_* из .lte →\n"
                                     "токи обмоток (сырой arg и физические амперы, % номинала).")
        self.btn_open_report = QPushButton("Открыть отчёт (.docx)")
        self.btn_open_report.clicked.connect(lambda: self._open(REPORT))
        # .docx собирается чистым Python и открывается и в Word, и в LibreOffice; конвертация — для
        # другого формата (правка в LibreOffice, рассылка PDF) и делается самой LibreOffice.
        self.conv_fmt = QComboBox()
        for f in ("odt", "pdf"):
            self.conv_fmt.addItem(f"в {f.upper()}", f)
        self.btn_convert = QPushButton("Конвертировать")
        self.btn_convert.setToolTip("office_convert.py — конвертация отчёта силами LibreOffice "
                                    "(soffice --headless). Нужна установленная LibreOffice; путь к ней "
                                    "можно задать переменной SS_SOFFICE.")
        self.btn_convert.clicked.connect(self._convert_report)
        self.btn_currents.clicked.connect(self._run_currents)
        self.stale = QLabel("")
        self.stale.setWordWrap(True)

        box2 = QGroupBox("Проектные токи обмоток (по уже измеренным матрицам — мгновенно)")
        v2 = QVBoxLayout(box2)
        rr = QHBoxLayout()
        rr.addWidget(self.btn_currents)
        rr.addWidget(self.btn_open_report)
        rr.addWidget(self.conv_fmt)
        rr.addWidget(self.btn_convert)
        rr.addStretch(1)
        v2.addLayout(rr)
        v2.addWidget(self.stale)

        self.table = CsvEditor(read_only=True)
        self.table.message.connect(lambda t, e: self.say(t, e))
        v2.addWidget(self.table, 1)

        self.guard(self.btn_matrices, self.btn_library, self.btn_one, self.btn_compare,
                   self.btn_currents, self.btn_convert)

        lay = QVBoxLayout(self)
        lay.addWidget(SectionLabel("Матрицы связи ток → физические коэффициенты"))
        lay.addWidget(box1)
        lay.addWidget(box2, 1)

    # ---- обновление ------------------------------------------------------
    def refresh(self) -> None:
        self._fill_elements()
        self._update_stale()
        p = paths.table("design_currents.csv")
        if p.exists():
            self.table.load(p)

    def _fill_elements(self) -> None:
        """Список элементов — из таблицы матриц (что уже измерено), иначе из решётки проекта."""
        cur = self.elem.currentText()
        names: list[str] = []
        tab = paths.table("coupling_matrices.csv")
        if tab.exists():
            import csv
            with open(tab, encoding="utf-8-sig", newline="") as f:
                for r in csv.DictReader(f, delimiter=";"):
                    if r.get("element") and r["element"] not in names:
                        names.append(r["element"])
        self.elem.clear()
        self.elem.addItems(names)
        if cur in names:
            self.elem.setCurrentText(cur)

    def _update_stale(self) -> None:
        try:
            import sys
            sys.path.insert(0, str(paths.ROOT))
            import design_currents as dc
            msg = dc.staleness()
        except Exception as e:                       # инструмент недоступен — просто молчим
            self.stale.setText("")
            self.say(f"проверка свежести матриц недоступна: {type(e).__name__}: {e}", True)
            return
        if msg:
            self.stale.setText("⚠ " + msg.splitlines()[0].replace("ВНИМАНИЕ: ", ""))
            self.stale.setStyleSheet("color:#b06000;")
        else:
            self.stale.setText("матрицы в базе свежее исходных данных — токи считаются по актуальным матрицам")
            self.stale.setStyleSheet("color:#666;")

    def _open(self, p) -> None:
        err = open_in_os(p)
        if err:
            self.say(err, True)

    # ---- команды ---------------------------------------------------------
    def _run_matrices(self) -> None:
        self.runner.run(py_job("матрицы связи + Word-отчёт", "element_compare.py", "matrices",
                               on_done=lambda c: self.refresh()))

    def _run_library(self) -> None:
        flags = []
        flags.append("--report" if self.ck_report.isChecked() else "--no-report")
        flags.append("--demo" if self.ck_demo.isChecked() else "--no-demo")
        flags.append("--html" if self.ck_html.isChecked() else "--no-html")
        self.runner.run(py_job("полный конвейер матриц (library)", "element_compare.py", "library", *flags,
                               on_done=lambda c: self.refresh()))

    def _run_one(self) -> None:
        name = self.elem.currentText().strip()
        if not name:
            self.say("не выбран элемент", True)
            return
        self.runner.run(py_job(f"матрица связи: {name}", "element_compare.py", "library", name,
                               on_done=lambda c: self.refresh()))

    def _run_all(self) -> None:
        self.runner.run(py_job("сравнение моделей по элементам", "element_compare.py", "all",
                               "quad,sext,oct,dipole,dq", on_done=lambda c: self.refresh()))

    def _convert_report(self) -> None:
        """Конвертировать готовый отчёт в выбранный формат (делает LibreOffice)."""
        if not REPORT.exists():
            self.say(f"нет отчёта {REPORT} — сначала соберите его", True)
            return
        fmt = self.conv_fmt.currentData()
        self.runner.run(py_job(f"конвертация отчёта в {fmt.upper()}", "office_convert.py",
                               str(REPORT), "--to", fmt))

    def _run_currents(self) -> None:
        lat = self.project.lattice if self.project else 1
        self.runner.run(py_job("проектные токи по матрицам связи", "design_currents.py", "--lattice", lat,
                               on_done=lambda c: self.refresh()))
