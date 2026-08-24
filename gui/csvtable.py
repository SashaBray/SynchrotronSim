# -*- coding: utf-8 -*-
"""Редактор CSV-таблиц базы данных проекта.

Таблицы БД — текстовые CSV с разделителем `;` и ТОЧКОЙ-десятичной (машинный формат, его читают оба
ядра и C++ по индексам столбцов). Редактор сохраняет файл ровно в том же виде: тот же разделитель,
та же кодировка (BOM сохраняется), порядок и состав столбцов не меняются — иначе ядра перестанут
находить свои колонки. Перед записью делается резервная копия в `workspace/tables/backups/`.
"""
from __future__ import annotations

import csv
import datetime as dt
import shutil
from pathlib import Path

from PySide6.QtCore import QAbstractTableModel, QModelIndex, QSortFilterProxyModel, Qt, Signal
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (QAbstractItemView, QHBoxLayout, QHeaderView, QLabel, QLineEdit,
                               QMessageBox, QPushButton, QTableView, QVBoxLayout, QWidget)

from . import paths
from .widgets import MONO


def read_csv(path: Path) -> tuple[list[str], list[list[str]], str]:
    """Прочитать `;`-CSV → (заголовок, строки, кодировка-для-записи с учётом BOM)."""
    raw = path.read_bytes()
    enc = "utf-8-sig" if raw.startswith(b"\xef\xbb\xbf") else "utf-8"
    text = raw.decode(enc, errors="replace")
    rows = list(csv.reader(text.splitlines(), delimiter=";"))
    if not rows:
        return [], [], enc
    return rows[0], [r for r in rows[1:] if r], enc


def write_csv(path: Path, header: list[str], rows: list[list[str]], enc: str = "utf-8") -> Path:
    """Записать таблицу, предварительно сняв резервную копию. Возврат — путь копии.

    Переводы строк сохраняются такими, какими были в файле: в базе соседствуют таблицы с CRLF и с LF,
    и смена концов строк давала бы «изменённый» файл с тем же содержимым (лишний шум в истории)."""
    backups = paths.TABLES / "backups"
    backups.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = backups / f"{path.stem}.{stamp}.csv"
    nl = "\r\n"
    if path.exists():
        shutil.copy2(path, bak)
        head = path.read_bytes()[:65536]
        nl = "\r\n" if b"\r\n" in head else ("\n" if b"\n" in head else "\r\n")
    with open(path, "w", encoding=enc, newline="") as f:
        w = csv.writer(f, delimiter=";", lineterminator=nl)
        w.writerow(header)
        w.writerows(rows)
    return bak


class CsvTableModel(QAbstractTableModel):
    """Модель `;`-CSV: значения хранятся строками (никаких приведений типов — формат чисел не портится)."""

    dirty_changed = Signal(bool)

    def __init__(self, parent=None):
        super().__init__(parent)
        self.path: Path | None = None
        self.header: list[str] = []
        self.rows: list[list[str]] = []
        self.enc = "utf-8"
        self._dirty = False

    # ---- загрузка/сохранение --------------------------------------------
    def load(self, path: Path) -> None:
        self.beginResetModel()
        self.path = path
        self.header, self.rows, self.enc = read_csv(path)
        self.endResetModel()
        self._set_dirty(False)

    def save(self) -> Path:
        if self.path is None:
            raise RuntimeError("таблица не загружена")
        bak = write_csv(self.path, self.header, self.rows, self.enc)
        self._set_dirty(False)
        return bak

    @property
    def dirty(self) -> bool:
        return self._dirty

    def _set_dirty(self, v: bool) -> None:
        if v != self._dirty:
            self._dirty = v
            self.dirty_changed.emit(v)

    # ---- Qt --------------------------------------------------------------
    def rowCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self.rows)

    def columnCount(self, parent=QModelIndex()) -> int:
        return 0 if parent.isValid() else len(self.header)

    def data(self, index, role=Qt.DisplayRole):
        if not index.isValid() or role not in (Qt.DisplayRole, Qt.EditRole):
            return None
        row = self.rows[index.row()]
        return row[index.column()] if index.column() < len(row) else ""

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return None
        return self.header[section] if orientation == Qt.Horizontal else section + 1

    def flags(self, index):
        return Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsEditable

    def setData(self, index, value, role=Qt.EditRole) -> bool:
        if role != Qt.EditRole or not index.isValid():
            return False
        row = self.rows[index.row()]
        while len(row) < len(self.header):                 # выровнять «рваную» строку под схему
            row.append("")
        if row[index.column()] == str(value):
            return False
        row[index.column()] = str(value)
        self.dataChanged.emit(index, index, [Qt.DisplayRole, Qt.EditRole])
        self._set_dirty(True)
        return True

    # ---- правка строк ----------------------------------------------------
    def add_row(self, after: int | None = None) -> int:
        """Новая строка — копия выделенной (удобно: id меняешь, остальное готово) либо пустая."""
        src = self.rows[after] if (after is not None and 0 <= after < len(self.rows)) else None
        new = list(src) if src else [""] * len(self.header)
        at = (after + 1) if after is not None else len(self.rows)
        self.beginInsertRows(QModelIndex(), at, at)
        self.rows.insert(at, new)
        self.endInsertRows()
        self._set_dirty(True)
        return at

    def remove_rows(self, indexes: list[int]) -> None:
        for r in sorted(set(indexes), reverse=True):
            if 0 <= r < len(self.rows):
                self.beginRemoveRows(QModelIndex(), r, r)
                self.rows.pop(r)
                self.endRemoveRows()
        self._set_dirty(True)


class NumericSortProxy(QSortFilterProxyModel):
    """Сортировка по клику на заголовок: числа — как числа, остальное — как текст.

    Значения хранятся строками, а обычное строковое сравнение ставит «16» перед «2» — для таблиц,
    где столбцы это номера опытов, токи и координаты, это бесполезно.
    """

    @staticmethod
    def _num(v: str):
        try:
            return float(str(v).replace(",", "."))            # запятая-десятичная тоже читается
        except (TypeError, ValueError):
            return None

    def lessThan(self, left, right) -> bool:
        a, b = self.sourceModel().data(left), self.sourceModel().data(right)
        na, nb = self._num(a), self._num(b)
        if na is not None and nb is not None:
            return na < nb
        if (na is None) != (nb is None):                      # числа выше пустых/текстовых ячеек
            return na is not None
        return str(a) < str(b)


class CsvEditor(QWidget):
    """Просмотр и правка одной таблицы БД: фильтр, добавление/удаление строк, сохранение с бэкапом."""

    saved = Signal(str)         # путь сохранённой таблицы
    message = Signal(str, bool)  # (текст, ошибка) — в общий журнал

    def __init__(self, parent=None, read_only: bool = False):
        super().__init__(parent)
        self.read_only = read_only
        self.model = CsvTableModel(self)
        self.proxy = NumericSortProxy(self)
        self.proxy.setSourceModel(self.model)
        self.proxy.setFilterKeyColumn(-1)                 # фильтр по всем столбцам
        self.proxy.setFilterCaseSensitivity(Qt.CaseInsensitive)

        self.view = QTableView()
        self.view.setModel(self.proxy)
        self.view.setSortingEnabled(True)
        self.view.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)   # по умолчанию — порядок файла
        self.view.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.view.setEditTriggers(QAbstractItemView.DoubleClicked | QAbstractItemView.EditKeyPressed |
                                  QAbstractItemView.AnyKeyPressed)
        self.view.setFont(QFont(MONO, 9))
        self.view.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.view.horizontalHeader().setStretchLastSection(True)
        self.view.verticalHeader().setDefaultSectionSize(20)

        self.filter = QLineEdit()
        self.filter.setPlaceholderText("фильтр по любому столбцу…")
        self.filter.textChanged.connect(self.proxy.setFilterFixedString)

        self.btn_add = QPushButton("Добавить строку")
        self.btn_del = QPushButton("Удалить строки")
        self.btn_save = QPushButton("Сохранить")
        self.btn_reload = QPushButton("Перечитать")
        self.btn_add.clicked.connect(self._add)
        self.btn_del.clicked.connect(self._del)
        self.btn_save.clicked.connect(self.save)
        self.btn_reload.clicked.connect(self.reload)
        self.btn_save.setEnabled(False)
        self.model.dirty_changed.connect(self.btn_save.setEnabled)

        self.info = QLabel("таблица не выбрана")
        self.info.setStyleSheet("color: #666;")

        if read_only:                                      # производная таблица — смотреть, но не править
            self.view.setEditTriggers(QAbstractItemView.NoEditTriggers)
            for b in (self.btn_add, self.btn_del, self.btn_save):
                b.hide()

        top = QHBoxLayout()
        top.addWidget(self.filter, 1)
        top.addWidget(self.btn_add)
        top.addWidget(self.btn_del)
        top.addWidget(self.btn_reload)
        top.addWidget(self.btn_save)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addLayout(top)
        lay.addWidget(self.view, 1)
        lay.addWidget(self.info)

    # ---- API -------------------------------------------------------------
    def load(self, path: Path) -> None:
        if self.model.dirty and not self._confirm_discard():
            return
        if not path.exists():
            self.message.emit(f"нет файла: {path}", True)
            self.info.setText(f"нет файла: {path.name}")
            return
        self.model.load(path)
        self.proxy.sort(-1)                                   # новая таблица — в порядке файла, без сортировки
        self.view.horizontalHeader().setSortIndicator(-1, Qt.AscendingOrder)
        self.view.resizeColumnsToContents()
        self.info.setText(f"{path.name}: строк {len(self.model.rows)}, столбцов {len(self.model.header)}  —  {path}")

    def reload(self) -> None:
        if self.model.path is not None:
            if self.model.dirty and not self._confirm_discard():
                return
            self.model.load(self.model.path)
            self.info.setText(f"перечитано: {self.model.path.name}")

    def save(self) -> None:
        if self.model.path is None:
            return
        try:
            bak = self.model.save()
        except Exception as e:
            QMessageBox.critical(self, "Ошибка записи", f"{type(e).__name__}: {e}")
            self.message.emit(f"не сохранено: {type(e).__name__}: {e}", True)
            return
        self.message.emit(f"сохранено: {self.model.path.name} (копия прежней → {bak.name})", False)
        self.saved.emit(str(self.model.path))

    # ---- внутреннее ------------------------------------------------------
    def _current_source_rows(self) -> list[int]:
        return sorted({self.proxy.mapToSource(i).row() for i in self.view.selectionModel().selectedRows()})

    def _add(self) -> None:
        if self.model.path is None:
            return
        sel = self._current_source_rows()
        at = self.model.add_row(sel[-1] if sel else None)
        self.message.emit(f"добавлена строка {at + 1}" + (" (копия выделенной)" if sel else ""), False)

    def _del(self) -> None:
        rows = self._current_source_rows()
        if not rows:
            return
        if QMessageBox.question(self, "Удалить строки", f"Удалить выделенные строки ({len(rows)})?") \
                != QMessageBox.Yes:
            return
        self.model.remove_rows(rows)

    def _confirm_discard(self) -> bool:
        return QMessageBox.question(self, "Несохранённые правки",
                                    "В таблице есть несохранённые правки. Отбросить их?") == QMessageBox.Yes
