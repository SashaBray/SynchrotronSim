# -*- coding: utf-8 -*-
"""Мелкие общие виджеты оболочки: панель логов, заголовки секций, выбор пути."""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QFont, QTextCharFormat, QTextCursor
from PySide6.QtWidgets import (QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
                               QPushButton, QVBoxLayout, QWidget)

MONO = "Consolas" if sys.platform == "win32" else "monospace"


def open_in_os(path: str | Path) -> str:
    """Открыть файл/папку системным приложением. Возврат — сообщение об ошибке ('' если удалось)."""
    p = Path(path)
    if not p.exists():
        return f"не найдено: {p}"
    try:
        if sys.platform == "win32":
            os.startfile(str(p))                       # noqa: S606 — штатный способ Windows
        elif sys.platform == "darwin":
            subprocess.Popen(["open", str(p)])
        else:
            subprocess.Popen(["xdg-open", str(p)])
    except Exception as e:
        return f"{type(e).__name__}: {e}"
    return ""


def plural(n: int, one: str, few: str, many: str) -> str:
    """Русская форма слова при числе: 1 оборот, 2 оборота, 5 оборотов, 11 оборотов."""
    n = abs(int(n))
    if n % 10 == 1 and n % 100 != 11:
        return one
    if 2 <= n % 10 <= 4 and not (12 <= n % 100 <= 14):
        return few
    return many


class SectionLabel(QLabel):
    """Подзаголовок группы элементов управления."""

    def __init__(self, text: str, parent=None):
        super().__init__(text, parent)
        f = self.font()
        f.setBold(True)
        self.setFont(f)
        self.setContentsMargins(0, 8, 0, 2)


def hline() -> QFrame:
    f = QFrame()
    f.setFrameShape(QFrame.HLine)
    f.setFrameShadow(QFrame.Sunken)
    return f


class PathPicker(QWidget):
    """Поле пути + кнопка «Обзор» (файл или папка)."""

    def __init__(self, placeholder: str = "", directory: bool = False, filt: str = "", parent=None):
        super().__init__(parent)
        self._dir, self._filt = directory, filt
        self.edit = QLineEdit()
        self.edit.setPlaceholderText(placeholder)
        btn = QPushButton("Обзор…")
        btn.clicked.connect(self._browse)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(self.edit, 1)
        lay.addWidget(btn)

    def _browse(self) -> None:
        start = self.edit.text().strip() or str(Path.cwd())
        if self._dir:
            p = QFileDialog.getExistingDirectory(self, "Выберите папку", start)
        else:
            p, _ = QFileDialog.getOpenFileName(self, "Выберите файл", start, self._filt or "Все файлы (*)")
        if p:
            self.edit.setText(os.path.normpath(p))

    def text(self) -> str:
        return self.edit.text().strip()

    def setText(self, t: str) -> None:
        self.edit.setText(t)


class LogPanel(QWidget):
    """Общее окно логов: поток вывода дочерних процессов, ошибки — красным, заголовки задач — жирным."""

    MAX_LINES = 20000

    def __init__(self, parent=None):
        super().__init__(parent)
        self.view = QPlainTextEdit()
        self.view.setReadOnly(True)
        self.view.setMaximumBlockCount(self.MAX_LINES)
        self.view.setFont(QFont(MONO, 9))
        self.view.setLineWrapMode(QPlainTextEdit.NoWrap)

        self.btn_clear = QPushButton("Очистить")
        self.btn_save = QPushButton("Сохранить…")
        self.btn_stop = QPushButton("Остановить")     # общая для всей очереди: долгие прогоны идут с любой вкладки
        self.btn_stop.setEnabled(False)
        self.btn_stop.setToolTip("Прервать текущую задачу.\nТаблицы базы, которые она временно правит, "
                                 "вернутся из копии при следующем запуске инструмента.")
        self.btn_clear.clicked.connect(self.view.clear)
        self.btn_save.clicked.connect(self._save)
        self.status = QLabel("готов")
        self.status.setStyleSheet("color: #666;")

        bar = QHBoxLayout()
        bar.addWidget(QLabel("Журнал"))
        bar.addWidget(self.status, 1)
        bar.addWidget(self.btn_stop)
        bar.addWidget(self.btn_clear)
        bar.addWidget(self.btn_save)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 2, 6, 6)
        lay.addLayout(bar)
        lay.addWidget(self.view, 1)

    # ---- вывод -----------------------------------------------------------
    def append(self, text: str, err: bool = False) -> None:
        self._write(text, QColor("#c0392b") if err else None)

    def append_title(self, text: str) -> None:
        self._write(f"\n=== {text} ===", QColor("#1a5fb4"), bold=True)

    def _write(self, text: str, color: QColor | None, bold: bool = False) -> None:
        cur = self.view.textCursor()
        cur.movePosition(QTextCursor.End)
        fmt = QTextCharFormat()
        if color is not None:
            fmt.setForeground(color)
        if bold:
            fmt.setFontWeight(QFont.Bold)
        cur.insertText(text + "\n", fmt)
        self.view.setTextCursor(cur)
        self.view.ensureCursorVisible()

    def set_status(self, text: str, busy: bool = False) -> None:
        self.status.setText(("⏳ " if busy else "") + text)
        self.status.setStyleSheet("color: #b06000;" if busy else "color: #666;")

    def _save(self) -> None:
        p, _ = QFileDialog.getSaveFileName(self, "Сохранить журнал", "synchrotronsim_log.txt",
                                           "Текст (*.txt);;Все файлы (*)")
        if p:
            Path(p).write_text(self.view.toPlainText(), encoding="utf-8")
            self.append(f"[журнал сохранён: {p}]")
