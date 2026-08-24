# -*- coding: utf-8 -*-
"""Главное окно оболочки и выбор проекта."""
from __future__ import annotations

import sys

from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor
from PySide6.QtWidgets import (QApplication, QDialog, QDialogButtonBox, QHBoxLayout, QLabel,
                               QListWidget, QListWidgetItem, QMainWindow, QMessageBox, QPushButton,
                               QSplitter, QTabWidget, QVBoxLayout, QWidget)

from . import paths, projects
from .projects import Project, Registry
from .runner import Runner
from .tab_matrices import MatricesTab
from .tab_project import ProjectTab
from .tab_simulation import SimulationTab
from .tab_viz import VizTab
from .widgets import LogPanel, open_in_os
from .wizard import NewProjectWizard


class ProjectDialog(QDialog):
    """Стартовое окно: открыть проект из реестра или создать новый."""

    def __init__(self, registry: Registry, parent=None):
        super().__init__(parent)
        self.setWindowTitle("SynchrotronSim — выбор проекта")
        self.resize(680, 420)
        self.registry = registry
        self.chosen: Project | None = None
        self.create_new = False
        self._found: list[Project] = []            # решётки из базы, ещё не заведённые проектами

        self.list = QListWidget()
        self.list.itemDoubleClicked.connect(lambda _i: self._open())
        self.list.currentItemChanged.connect(self._show)
        self.details = QLabel("выберите проект слева")
        self.details.setWordWrap(True)
        self.details.setAlignment(Qt.AlignTop)
        self.details.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.details.setStyleSheet("color:#333;")

        self.btn_new = QPushButton("Создать новый проект…")
        self.btn_del = QPushButton("Убрать из списка")
        self.btn_new.clicked.connect(self._new)
        self.btn_del.clicked.connect(self._remove)

        buttons = QDialogButtonBox(QDialogButtonBox.Open | QDialogButtonBox.Close)
        buttons.button(QDialogButtonBox.Open).setText("Открыть проект")
        buttons.button(QDialogButtonBox.Close).setText("Выход")
        buttons.accepted.connect(self._open)
        buttons.rejected.connect(self.reject)
        self.buttons = buttons

        cols = QHBoxLayout()
        cols.addWidget(self.list, 1)
        cols.addWidget(self.details, 1)
        side = QHBoxLayout()
        side.addWidget(self.btn_new)
        side.addWidget(self.btn_del)
        side.addStretch(1)

        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Проект связывает файл решётки .lte с решёткой и опытом в базе данных."))
        lay.addLayout(cols, 1)
        lay.addLayout(side)
        lay.addWidget(buttons)
        self._fill()

    def _fill(self) -> None:
        """Список: сначала заведённые проекты, затем решётки, уже лежащие в базе."""
        self.list.clear()
        self._found = projects.discover_projects(self.registry.items)
        for p in self.registry.items:
            it = QListWidgetItem(p.label())
            it.setData(Qt.UserRole, ("reg", p.name))
            self.list.addItem(it)
        if self._found:
            head = QListWidgetItem("— найдено в базе —")
            head.setFlags(Qt.NoItemFlags)
            head.setForeground(QColor("#888"))
            self.list.addItem(head)
            for i, p in enumerate(self._found):
                it = QListWidgetItem(p.label())
                it.setData(Qt.UserRole, ("found", i))
                it.setForeground(QColor("#1a5fb4"))
                self.list.addItem(it)
        first = next((i for i in range(self.list.count())
                      if self.list.item(i).data(Qt.UserRole) is not None), None)
        if first is not None:
            self.list.setCurrentRow(first)
        else:
            self.details.setText("Проектов пока нет, и база пуста.\n\nНажмите «Создать новый проект…» — "
                                 "мастер спросит файл решётки, проверит замкнутость кольца, поможет "
                                 "задать пучок и загрузит решётку в базу.")
        self.buttons.button(QDialogButtonBox.Open).setEnabled(first is not None)

    def _current(self) -> Project | None:
        it = self.list.currentItem()
        d = it.data(Qt.UserRole) if it else None
        if not d:
            return None
        kind, key = d
        return self.registry.get(key) if kind == "reg" else self._found[key]

    def _is_found(self) -> bool:
        it = self.list.currentItem()
        d = it.data(Qt.UserRole) if it else None
        return bool(d) and d[0] == "found"

    def _show(self) -> None:
        p = self._current()
        if p is None:
            return
        extra = ("\n\nЭта решётка уже есть в базе, но ещё не заведена проектом. "
                 "Откройте её — она попадёт в список проектов.") if self._is_found() else ""
        self.details.setText(p.describe() + extra)
        self.btn_del.setEnabled(not self._is_found())

    def _open(self) -> None:
        p = self._current()
        if p is None:
            return
        if self._is_found():                       # найденную в базе решётку заводим проектом при открытии
            self.registry.add(p)
        self.chosen = p
        self.accept()

    def _new(self) -> None:
        self.create_new = True
        self.accept()

    def _remove(self) -> None:
        p = self._current()
        if p is None:
            return
        if QMessageBox.question(self, "Убрать проект",
                                f"Убрать «{p.name}» из списка?\nФайлы решётки и база не изменятся.") \
                == QMessageBox.Yes:
            self.registry.remove(p.name)
            self._fill()


class MainWindow(QMainWindow):
    """Вкладки конвейера + общий журнал внизу."""

    def __init__(self, registry: Registry, project: Project | None):
        super().__init__()
        self.registry = registry
        self.project = project
        self.resize(1280, 860)

        self.log = LogPanel()
        self.runner = Runner(self)
        self.log.btn_stop.clicked.connect(self._stop)
        self.runner.line.connect(self.log.append)
        self.runner.job_started.connect(self._job_started)
        self.runner.job_finished.connect(self._job_finished)
        self.runner.busy_changed.connect(self._busy)

        self.tabs = QTabWidget()
        self.tab_project = ProjectTab(self.runner, self.log)
        self.tab_matrices = MatricesTab(self.runner, self.log)
        self.tab_sim = SimulationTab(self.runner, self.log)
        self.tab_viz = VizTab(self.runner, self.log)
        self.tabs.addTab(self.tab_project, "Проект и база")
        self.tabs.addTab(self.tab_matrices, "Матрицы связи")
        self.tabs.addTab(self.tab_sim, "Симуляция")
        self.tabs.addTab(self.tab_viz, "Визуализация")
        self.tabs.currentChanged.connect(lambda _i: self._current_tab().refresh())

        split = QSplitter(Qt.Vertical)
        split.addWidget(self.tabs)
        split.addWidget(self.log)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        central = QWidget()
        v = QVBoxLayout(central)
        v.setContentsMargins(6, 6, 6, 0)
        v.addWidget(split)
        self.setCentralWidget(central)

        self._menu()
        self.set_project(project)
        self.log.append(f"рабочая папка: {paths.ROOT}")
        self.log.append(f"интерпретатор: {paths.python_exe()}")

    # ---- меню -------------------------------------------------------------
    def _menu(self) -> None:
        m = self.menuBar().addMenu("Проект")
        a_switch = QAction("Сменить проект…", self)
        a_switch.triggered.connect(self._switch)
        a_new = QAction("Новый проект…", self)
        a_new.triggered.connect(self._new_project)
        m.addAction(a_switch)
        m.addAction(a_new)
        m.addSeparator()
        for title, path in (("Папка результатов", paths.RESULTS), ("Папка таблиц", paths.TABLES),
                            ("Папка HTML-отчётов", paths.HTML_DIR)):
            act = QAction(title, self)
            act.triggered.connect(lambda _c=False, p=path: open_in_os(p))
            m.addAction(act)
        m.addSeparator()
        a_quit = QAction("Выход", self)
        a_quit.triggered.connect(self.close)
        m.addAction(a_quit)

        h = self.menuBar().addMenu("Справка")
        a_manual = QAction("Мануал (docs/MANUAL.md)", self)
        a_manual.triggered.connect(lambda: open_in_os(paths.ROOT / "docs" / "MANUAL.md"))
        a_about = QAction("О программе", self)
        a_about.triggered.connect(self._about)
        h.addAction(a_manual)
        h.addAction(a_about)

    def _about(self) -> None:
        QMessageBox.information(
            self, "О программе",
            "Оболочка SynchrotronSim.\n\n"
            "Все расчёты выполняют штатные инструменты проекта — оболочка их запускает и показывает "
            "вывод. Поэтому результат из окна и из командной строки совпадает.\n\n"
            f"Корень проекта: {paths.ROOT}")

    # ---- проект -----------------------------------------------------------
    def _current_tab(self):
        return self.tabs.currentWidget()

    def set_project(self, project: Project | None) -> None:
        self.project = project
        for t in (self.tab_project, self.tab_matrices, self.tab_sim, self.tab_viz):
            t.set_project(project)
        name = project.name if project else "без проекта"
        self.setWindowTitle(f"SynchrotronSim — {name}")
        if project:
            self.log.append_title(f"проект: {project.name} (решётка {project.lattice}, опыт {project.exp})")

    def _switch(self) -> None:
        if self._busy_guard():
            return
        dlg = ProjectDialog(self.registry, self)
        if dlg.exec() != QDialog.Accepted:
            return
        if dlg.create_new:
            self._new_project()
        elif dlg.chosen is not None:
            self.set_project(dlg.chosen)

    def _new_project(self) -> None:
        if self._busy_guard():
            return
        wiz = NewProjectWizard(self)
        if wiz.exec() != QDialog.Accepted:
            return
        proj = wiz.result_project()
        msg = wiz.write_beam(proj.group)
        if msg:
            self.log.append(msg)
        self.registry.add(proj)
        self.set_project(proj)
        self.runner.run(*wiz.load_jobs(proj))

    # ---- состояние --------------------------------------------------------
    def _job_started(self, title: str) -> None:
        self.log.append_title(title)
        self.log.set_status(title, busy=True)

    def _job_finished(self, title: str, code: int) -> None:
        if code == 0:
            self.statusBar().showMessage(f"{title}: готово", 8000)
        else:
            self.statusBar().showMessage(f"{title}: ошибка (код {code})", 12000)
        for t in (self.tab_project, self.tab_matrices, self.tab_sim, self.tab_viz):
            t.refresh()

    def _stop(self) -> None:
        """Прервать очередь. Инструменты, которые временно правят таблицы базы, при убийстве процесса
        не успевают их вернуть — восстановление произойдёт при следующем их запуске, о чём и сообщаем."""
        self.runner.stop()
        self.log.append("прерванный инструмент вернёт временно изменённые таблицы базы "
                        "при следующем запуске (проверка метки прогона)")

    def _busy(self, busy: bool) -> None:
        for t in (self.tab_project, self.tab_matrices, self.tab_sim, self.tab_viz):
            t.set_busy(busy)
        self.log.btn_stop.setEnabled(busy)
        if not busy:
            self.log.set_status("готов")

    def _busy_guard(self) -> bool:
        if self.runner.busy:
            QMessageBox.information(self, "Идёт расчёт",
                                    "Дождитесь окончания текущей задачи или остановите её "
                                    "на вкладке «Симуляция».")
            return True
        return False

    def closeEvent(self, event) -> None:
        if self.runner.busy:
            if QMessageBox.question(self, "Выход", "Идёт расчёт. Прервать его и выйти?") != QMessageBox.Yes:
                event.ignore()
                return
            self.runner.stop()
        event.accept()


def run() -> int:
    """Точка входа: выбор проекта, затем главное окно."""
    app = QApplication(sys.argv)
    app.setApplicationName("SynchrotronSim")
    registry = Registry()

    dlg = ProjectDialog(registry)
    if dlg.exec() != QDialog.Accepted:
        return 0

    win = MainWindow(registry, None if dlg.create_new else dlg.chosen)
    win.show()
    if dlg.create_new:
        win._new_project()
        if win.project is None:                       # мастер отменили — предложить выбор снова
            win._switch()
            if win.project is None:
                return 0
    return app.exec()
