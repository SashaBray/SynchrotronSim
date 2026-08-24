# -*- coding: utf-8 -*-
"""Общий предок вкладок: доступ к очереди команд, журналу и текущему проекту."""
from __future__ import annotations

from PySide6.QtWidgets import QWidget

from .projects import Project
from .runner import Runner
from .widgets import LogPanel


class TabBase(QWidget):
    """Вкладка знает очередь (`runner`), журнал (`log`) и текущий проект."""

    def __init__(self, runner: Runner, log: LogPanel, parent=None):
        super().__init__(parent)
        self.runner = runner
        self.log = log
        self.project: Project | None = None
        self._buttons: list[QWidget] = []

    # ---- переопределяемое -------------------------------------------------
    def set_project(self, project: Project | None) -> None:
        """Проект сменился — обновить содержимое вкладки."""
        self.project = project
        self.refresh()

    def refresh(self) -> None:
        """Перечитать состояние из базы (после смены проекта или завершения команды)."""

    # ---- занятость --------------------------------------------------------
    def guard(self, *widgets: QWidget) -> None:
        """Пометить кнопки, которые нужно гасить на время работы очереди."""
        self._buttons.extend(widgets)

    def set_busy(self, busy: bool) -> None:
        for w in self._buttons:
            w.setEnabled(not busy)

    # ---- помощь -----------------------------------------------------------
    def say(self, text: str, err: bool = False) -> None:
        self.log.append(text, err)
