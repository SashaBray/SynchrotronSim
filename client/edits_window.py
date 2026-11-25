# -*- coding: utf-8 -*-
"""Окно «Правки проекта» кокпита (docs/DATA_PROTOCOL.md, этап 4): журнал правок параметров проекта
(synchrotronsim.edits.Journal) — что изменено, было → стало, применено ли в модели; отмена и повтор;
пересборка модели с правками, которые применяются только так; сохранение в хранилище (новая ревизия машины
через рабочую копию таблиц); сброс всех правок. Само окно ничего не меняет — шлёт сигналы кокпиту."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (QHBoxLayout, QHeaderView, QInputDialog, QLabel, QPushButton, QTableWidget,
                               QTableWidgetItem, QVBoxLayout, QWidget, QAbstractItemView)


class EditsWindow(QWidget):
    undo = Signal()
    redo = Signal()
    rebuild = Signal()
    save = Signal(str)          # заметка к ревизии
    discard = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Правки проекта")
        self.head = QLabel("журнал правок пуст")
        self.head.setWordWrap(True)
        self.head.setTextFormat(Qt.RichText)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["элемент", "поле", "было", "стало", "в модели"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.NoSelection)
        self.b_undo = QPushButton("Отменить"); self.b_redo = QPushButton("Повторить")
        self.b_rebuild = QPushButton("Пересобрать модель"); self.b_save = QPushButton("Сохранить ревизией…")
        self.b_discard = QPushButton("Сбросить все")
        self.b_rebuild.setToolTip("применить правки, которые не применяются на ходу (положение, длина, поле диполя…)")
        self.b_save.setToolTip("наложить правки на рабочую копию таблиц и зафиксировать новую ревизию машины")
        self.b_undo.clicked.connect(self.undo.emit); self.b_redo.clicked.connect(self.redo.emit)
        self.b_rebuild.clicked.connect(self.rebuild.emit); self.b_discard.clicked.connect(self.discard.emit)
        self.b_save.clicked.connect(self._ask_save)
        row = QHBoxLayout()
        for b in (self.b_undo, self.b_redo, self.b_rebuild, self.b_save, self.b_discard):
            row.addWidget(b)
        row.addStretch(1)
        lay = QVBoxLayout(self)
        lay.addWidget(self.head)
        lay.addWidget(self.table, 1)
        lay.addLayout(row)
        self.refresh(None)

    def _ask_save(self):
        note, ok = QInputDialog.getText(self, "Сохранить правки", "Заметка к ревизии:")
        if ok:
            self.save.emit(note.strip())

    def refresh(self, journal, reason: str = "", engine_fp: str = "") -> None:
        """Показать журнал (None — правка недоступна: reason)."""
        if journal is None:
            self.table.setRowCount(0)
            self.head.setText(f"<i>правка проекта недоступна</i>{': ' + reason if reason else ''}")
            for b in (self.b_undo, self.b_redo, self.b_rebuild, self.b_save, self.b_discard):
                b.setEnabled(False)
            return
        es = journal.entries
        self.table.setRowCount(len(es))
        from synchrotronsim.edits import TITLES, UNITS
        for n, e in enumerate(es):
            ed = e.edit
            u = UNITS.get(ed.column, "")
            state = "да" if e.applied else "ждёт пересборки"
            vals = [ed.name or ed.keys[0], TITLES.get(ed.column, ed.column) + (f", {u}" if u else ""),
                    ed.old[0] or "(пусто)", ed.new[0] or "(пусто)",
                    state + (" · на ходу" if ed.cost == "hot" else " · пересборкой")]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(str(v))
                if c == 4 and not e.applied:
                    it.setForeground(Qt.darkRed)
                self.table.setItem(n, c, it)
        txt = (f"машина <b>{journal.machine_id}</b> · база {journal.base_rev[:12]} · правок <b>{len(es)}</b>"
               + (" · <span style='color:#a00'>модель отстаёт — нужна пересборка</span>" if journal.rebuild_needed else ""))
        if engine_fp and not journal.rebuild_needed:
            ok = engine_fp == journal.machine().fingerprint
            txt += " · модель = журнал" if ok else " · <span style='color:#a00'>модель и журнал разошлись</span>"
        self.head.setText(txt)
        self.b_undo.setEnabled(bool(es)); self.b_redo.setEnabled(bool(journal.redo_stack))
        self.b_rebuild.setEnabled(journal.rebuild_needed); self.b_save.setEnabled(bool(es))
        self.b_discard.setEnabled(bool(es))
