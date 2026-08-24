# -*- coding: utf-8 -*-
"""Вкладка «Проект»: сводка по текущей решётке и правка таблиц базы данных.

Правка таблиц — это правка входных данных обоих ядер, поэтому рядом стоят кнопки пересборки: после
изменения ШАБЛОНА супериода (`sp_elements.csv`) нужно пересобрать `lattice_configs.csv`, иначе ядра
продолжат читать прежние значения.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (QComboBox, QGroupBox, QHBoxLayout, QLabel, QPushButton, QSplitter,
                               QVBoxLayout)
from PySide6.QtCore import Qt

from . import paths, projects
from .csvtable import CsvEditor
from .runner import py_job
from .tabbase import TabBase
from .widgets import SectionLabel, open_in_os

# Таблицы базы, доступные к правке: (файл, пояснение)
TABLES = [
    ("experiments.csv", "опыты: решётка, пучок, шаг, финиш-плоскость, параметры ядра"),
    ("particle_groups.csv", "пучки: масса, заряд, старт и скорость каждой частицы"),
    ("radiation_params.csv", "излучение по опытам: проектная энергия, ISR, ВЧ"),
    ("bpm_config.csv", "датчики положения пучка (BPM): несовершенства"),
    ("sp_elements.csv", "ШАБЛОН супериода: токи и калибровки (после правки — пересборка)"),
    ("lattice_configs.csv", "большая таблица: экземпляры магнитов, геометрия и коэффициенты (её читают оба ядра)"),
    ("hardware_library.csv", "библиотека типов: карты полей, витки, номиналы, повороты"),
    ("design_currents.csv", "проектные токи обмоток (результат расчёта по матрицам связи)"),
    ("coupling_matrices.csv", "матрицы связи ток→физкоэф (результат element_compare library)"),
]


class ProjectTab(TabBase):
    def __init__(self, runner, log, parent=None):
        super().__init__(runner, log, parent)

        # --- сводка ---------------------------------------------------------
        self.summary = QLabel("проект не выбран")
        self.summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.summary.setWordWrap(True)
        self.db_info = QLabel("")
        self.db_info.setStyleSheet("color:#666;")
        self.db_info.setWordWrap(True)

        self.btn_assemble = QPushButton("Пересобрать решётку")
        self.btn_assemble.setToolTip("assemble_lattice_configs.py — шаблон + геометрия → lattice_configs.csv.\n"
                                     "Нужен после правки sp_elements.csv.")
        self.btn_warm = QPushButton("Проектные mat_* из .lte")
        self.btn_warm.setToolTip("warmstart_ring.py — записать в решётку проектные коэффициенты K из .lte")
        self.btn_snap = QPushButton("Снимок базы")
        self.btn_snap.setToolTip("db_snapshot.py save — зафиксировать состояние (таблицы не хранятся в git)")
        self.btn_check = QPushButton("Проверить снимок")
        self.btn_reload_cfg = QPushButton("Перезагрузить .lte в базу")
        self.btn_reload_cfg.setToolTip("load_config.py --force — полная перезагрузка решётки из .lte "
                                       "(mat_* сбрасываются к проектным)")
        self.btn_assemble.clicked.connect(self._assemble)
        self.btn_warm.clicked.connect(self._warmstart)
        self.btn_snap.clicked.connect(self._snapshot)
        self.btn_check.clicked.connect(self._check_snapshot)
        self.btn_reload_cfg.clicked.connect(self._reload_config)
        self.guard(self.btn_assemble, self.btn_warm, self.btn_snap, self.btn_check, self.btn_reload_cfg)

        actions = QHBoxLayout()
        for b in (self.btn_assemble, self.btn_warm, self.btn_reload_cfg, self.btn_snap, self.btn_check):
            actions.addWidget(b)
        actions.addStretch(1)

        box = QGroupBox("Текущий проект")
        v = QVBoxLayout(box)
        v.addWidget(self.summary)
        v.addWidget(self.db_info)
        v.addLayout(actions)

        # --- таблицы --------------------------------------------------------
        self.picker = QComboBox()
        for name, hint in TABLES:
            self.picker.addItem(f"{name} — {hint}", name)
        self.picker.currentIndexChanged.connect(self._load_table)
        self.btn_open_dir = QPushButton("Открыть папку таблиц")
        self.btn_open_dir.clicked.connect(lambda: self._open(paths.TABLES))

        self.editor = CsvEditor()
        self.editor.message.connect(lambda t, e: self.say(t, e))
        self.editor.saved.connect(self._after_save)

        head = QHBoxLayout()
        head.addWidget(QLabel("Таблица:"))
        head.addWidget(self.picker, 1)
        head.addWidget(self.btn_open_dir)

        tbox = QGroupBox("Таблицы базы данных")
        tv = QVBoxLayout(tbox)
        tv.addLayout(head)
        tv.addWidget(self.editor, 1)

        split = QSplitter(Qt.Vertical)
        top = QVBoxLayout()
        top.addWidget(box)
        from PySide6.QtWidgets import QWidget
        w = QWidget()
        w.setLayout(top)
        split.addWidget(w)
        split.addWidget(tbox)
        split.setStretchFactor(1, 1)

        lay = QVBoxLayout(self)
        lay.addWidget(SectionLabel("Проект и база данных"))
        lay.addWidget(split, 1)

        self._load_table()

    # ---- обновление ------------------------------------------------------
    def refresh(self) -> None:
        if self.project is None:
            self.summary.setText("проект не выбран")
            self.db_info.setText("")
            return
        self.summary.setText("<pre style='margin:0'>" + self.project.describe() + "</pre>")
        in_db = projects.lattices_in_db()
        n = in_db.get(self.project.lattice, 0)
        exps = [e["exp_id"] for e in projects.experiments()]
        self.db_info.setText(
            f"В базе: решётка {self.project.lattice} — {n} элементов; "
            f"опыты: {', '.join(exps) if exps else '—'}; пучков: {len(projects.particle_groups())}."
            + ("" if n else "  ⚠ решётка пуста — загрузите .lte в базу."))

    def _open(self, p) -> None:
        err = open_in_os(p)
        if err:
            self.say(err, True)

    def _load_table(self) -> None:
        name = self.picker.currentData()
        if name:
            self.editor.load(paths.table(name))

    def _after_save(self, path: str) -> None:
        if Path(path).name == "sp_elements.csv":
            self.say("шаблон супериода изменён — нажмите «Пересобрать решётку», иначе ядра прочтут прежние значения")

    # ---- команды ---------------------------------------------------------
    def _assemble(self) -> None:
        self.runner.run(py_job("пересборка lattice_configs", "assemble_lattice_configs.py"),
                        py_job("синхронизация суб-решёток", "assemble_lattice_configs.py", "--sync-subs"))

    def _warmstart(self) -> None:
        if self.project is None:
            return
        self.runner.run(py_job(f"проектные mat_* решётки {self.project.lattice}",
                               "warmstart_ring.py", self.project.lattice))

    def _snapshot(self) -> None:
        note = self.project.name if self.project else "снимок из оболочки"
        self.runner.run(py_job("снимок базы", "db_snapshot.py", "save", "--note", note))

    def _check_snapshot(self) -> None:
        self.runner.run(py_job("сверка со снимком", "db_snapshot.py", "check"))

    def _reload_config(self) -> None:
        p = self.project
        if p is None:
            return
        args = [p.path]
        if p.line:
            args += ["--line", p.line]
        args += ["--lattice", p.lattice, "--exp", p.exp, "--force"]
        self.runner.run(py_job(f"перезагрузка решётки {p.lattice} из .lte", "load_config.py", *args,
                               on_done=lambda code: self.refresh()))
