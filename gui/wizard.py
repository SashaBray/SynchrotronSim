# -*- coding: utf-8 -*-
"""Мастер нового проекта: файл решётки → проверка кольца → пучок → загрузка в базу.

Мастер только СОБИРАЕТ решение (какие файлы, какой пучок, какие идентификаторы) и отдаёт его главному
окну списком команд — они выполняются общей очередью, а их вывод идёт в общий журнал. Исключение —
проверка кольца (`lte_ring.py`): её результат нужен прямо в мастере, чтобы решить, замыкать кольцо или
нет, поэтому она запускается тут же и ждёт завершения (разбор `.lte`, единицы секунд).
"""
from __future__ import annotations

import csv
import math
import subprocess
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QGuiApplication
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QMessageBox, QPlainTextEdit,
                               QPushButton, QRadioButton, QSpinBox, QVBoxLayout, QWizard, QWizardPage)

from . import paths, projects
from .runner import Job, py_job
from .widgets import MONO, PathPicker

ELECTRON_M = "9.1093837e-31"
ELECTRON_Q = "-1.602176634e-19"
V_SURROGATE = 2.967e8            # рабочая скорость пучка проекта (суррогатная энергия γ≈7)


# --------------------------------------------------------------------------- 1. решётка
class LatticePage(QWizardPage):
    """Имя проекта и файл(ы) решётки; формат определяется автоматически (load_config)."""

    def __init__(self):
        super().__init__()
        self.setTitle("Решётка")
        self.setSubTitle("Укажите файл или папку с описанием решётки (.lte). Формат определится сам.")
        self.cfg: dict | None = None

        self.name = QLineEdit()
        self.name.setPlaceholderText("например: СИЛА замкнутая")
        self.path = PathPicker("путь к .lte или к папке конфигурации", filt="Решётка (*.lte);;Все файлы (*)")
        self.line = PathPicker("файл последовательности *_line.lte (если он отдельный)",
                               filt="Последовательность (*.lte);;Все файлы (*)")
        self.btn_dir = QPushButton("Выбрать папку…")
        self.btn_dir.clicked.connect(self._pick_dir)

        self.known = QComboBox()
        self.known.addItem("— выбрать из workspace/lte/ —", None)
        for c in projects.discover_configs():
            if c["format"] != "none":
                self.known.addItem(f"{c['name']}  ({projects.lc.FMT_RU.get(c['format'], c['format'])})", c["dir"])
        self.known.currentIndexChanged.connect(self._pick_known)

        self.status = QLabel("формат: —")
        self.status.setWordWrap(True)
        self.path.edit.textChanged.connect(self._detect)
        self.line.edit.textChanged.connect(self._detect)

        form = QFormLayout(self)
        form.addRow("Название проекта:", self.name)
        form.addRow("Готовая конфигурация:", self.known)
        form.addRow("Путь к решётке:", self.path)
        row = QHBoxLayout()
        row.addWidget(self.btn_dir)
        row.addStretch(1)
        form.addRow("", row)
        form.addRow("Последовательность:", self.line)
        form.addRow("", self.status)

        self.registerField("proj_name*", self.name)
        self.registerField("lte_path*", self.path.edit)
        self.registerField("lte_line", self.line.edit)

    def _pick_dir(self) -> None:
        from PySide6.QtWidgets import QFileDialog
        d = QFileDialog.getExistingDirectory(self, "Папка конфигурации", str(paths.LTE_DIR))
        if d:
            self.path.setText(d)

    def _pick_known(self) -> None:
        d = self.known.currentData()
        if d:
            self.path.setText(d)
            self.line.setText("")
            if not self.name.text().strip():
                self.name.setText(Path(d).name)

    def _detect(self) -> None:
        p = self.path.text()
        self.cfg = None
        if not p:
            self.status.setText("формат: —")
        else:
            try:
                self.cfg = projects.config_from_path(p, self.line.text())
            except ValueError as e:
                self.status.setText(f"<span style='color:#c0392b'>{e}</span>")
            else:
                fmt = projects.lc.FMT_RU.get(self.cfg["format"], self.cfg["format"])
                files = self.cfg["descs"] + self.cfg["lines"] if self.cfg["format"] == "sila" else self.cfg["flat"]
                self.status.setText(f"формат: <b>{fmt}</b><br>" +
                                    "<br>".join(Path(f).name for f in files[:4]))
        self.completeChanged.emit()

    def isComplete(self) -> bool:
        return bool(self.name.text().strip()) and self.cfg is not None and self.cfg["format"] != "none"


# --------------------------------------------------------------------------- 2. проверка кольца
class RingPage(QWizardPage):
    """Отчёт `lte_ring.py`: замкнуто ли кольцо по углу и по положению; при желании — замкнуть копию."""

    def __init__(self):
        super().__init__()
        self.setTitle("Проверка кольца")
        self.setSubTitle("Замкнутость по углу (Σ ANGLE = 2π) и по положению. Проверка не обязательна.")
        self.closed_path = ""

        self.out = QPlainTextEdit()
        self.out.setReadOnly(True)
        self.out.setFont(QFont(MONO, 9))
        self.btn = QPushButton("Проверить кольцо")
        self.btn.clicked.connect(self._check)
        self.btn_close = QPushButton("Замкнуть → сохранить копию")
        self.btn_close.clicked.connect(self._close_ring)
        self.btn_close.setEnabled(False)
        self.note = QLabel("Замыкание масштабирует ANGLE и края E1/E2 чистых диполей; "
                           "combined-function элементы не трогаются.")
        self.note.setWordWrap(True)
        self.note.setStyleSheet("color:#666;")

        bar = QHBoxLayout()
        bar.addWidget(self.btn)
        bar.addWidget(self.btn_close)
        bar.addStretch(1)
        lay = QVBoxLayout(self)
        lay.addLayout(bar)
        lay.addWidget(self.out, 1)
        lay.addWidget(self.note)

    # --- вспомогательное
    def _lte_args(self) -> list[str] | None:
        page: LatticePage = self.wizard().page(0)
        cfg = page.cfg
        if not cfg:
            return None
        if cfg["format"] == "sila":
            if not (cfg["descs"] and cfg["lines"]):
                return None
            return [cfg["descs"][0], "--line", cfg["lines"][0]]
        return [cfg["flat"][0]] if cfg["flat"] else None

    def _run(self, extra: list[str]) -> tuple[int, str]:
        args = self._lte_args()
        if args is None:
            return 1, "не определены файлы решётки"
        cmd = [paths.python_exe(), str(paths.ROOT / "lte_ring.py")] + args + extra
        QGuiApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            r = subprocess.run(cmd, cwd=str(paths.ROOT), capture_output=True, text=True,
                               encoding="utf-8", errors="replace", env=paths.child_env(), timeout=600)
            return r.returncode, (r.stdout or "") + (("\n" + r.stderr) if r.stderr else "")
        except subprocess.TimeoutExpired:
            return 1, "проверка не уложилась в 10 минут"
        except OSError as e:
            return 1, f"{type(e).__name__}: {e}"
        finally:
            QGuiApplication.restoreOverrideCursor()

    def _check(self) -> None:
        code, text = self._run([])
        self.out.setPlainText(text.strip() or "(пусто)")
        self.btn_close.setEnabled(code == 0)

    def _close_ring(self) -> None:
        page: LatticePage = self.wizard().page(0)
        cfg = page.cfg
        src = Path((cfg["descs"] or cfg["flat"])[0])
        out = src.with_name(src.stem + "_closed.lte")
        code, text = self._run(["--out", str(out)])
        self.out.setPlainText(text.strip())
        if code == 0 and out.exists():
            self.closed_path = str(out)
            self.out.appendPlainText(f"\nзамкнутая копия: {out}\nона и будет загружена в базу.")
        else:
            QMessageBox.warning(self, "Замыкание", "Не удалось записать замкнутую копию — см. вывод.")


# --------------------------------------------------------------------------- 3. пучок
class BeamPage(QWizardPage):
    """Выбор готового пучка из базы или создание нового (центр + спутники)."""

    def __init__(self):
        super().__init__()
        self.setTitle("Пучок")
        self.setSubTitle("Начальные условия частиц: взять готовую группу из базы или задать новую.")

        groups = projects.particle_groups()
        self.rb_exist = QRadioButton("Использовать готовый пучок")
        self.rb_new = QRadioButton("Создать новый пучок")
        self.rb_exist.setChecked(bool(groups))
        self.rb_new.setChecked(not groups)

        self.combo = QComboBox()
        for gid, parts in sorted(groups.items()):
            p0 = parts[0]
            self.combo.addItem(f"группа {gid}: {len(parts)} частиц, старт "
                               f"({p0['x']}, {p0['y']}, {p0['z']})", gid)
        self.combo.setEnabled(bool(groups))

        self.n_part = QComboBox()
        self.n_part.addItem("5 — центр + 4 спутника (±r по y и z)", 5)
        self.n_part.addItem("13 — центр + 12 по окружности", 13)
        self.radius = QDoubleSpinBox()
        self.radius.setRange(0.0, 100.0)
        self.radius.setDecimals(3)
        self.radius.setValue(0.25)
        self.radius.setSuffix(" мм")
        self.speed = QDoubleSpinBox()
        self.speed.setRange(1e3, 2.998e8)
        self.speed.setDecimals(0)
        self.speed.setValue(V_SURROGATE)
        self.speed.setSuffix(" м/с")
        self.x0 = QDoubleSpinBox(); self.x0.setRange(-1e4, 1e4); self.x0.setDecimals(6); self.x0.setSuffix(" м")
        self.y0 = QDoubleSpinBox(); self.y0.setRange(-1e4, 1e4); self.y0.setDecimals(6); self.y0.setSuffix(" м")
        self.z0 = QDoubleSpinBox(); self.z0.setRange(-1e4, 1e4); self.z0.setDecimals(6); self.z0.setSuffix(" м")
        self.mass = QLineEdit(ELECTRON_M)
        self.charge = QLineEdit(ELECTRON_Q)

        box_new = QGroupBox("Новый пучок")
        f = QFormLayout(box_new)
        f.addRow("Состав:", self.n_part)
        f.addRow("Радиус пучка:", self.radius)
        f.addRow("Скорость вдоль X:", self.speed)
        f.addRow("Старт x:", self.x0)
        f.addRow("Старт y:", self.y0)
        f.addRow("Старт z:", self.z0)
        f.addRow("Масса частицы, кг:", self.mass)
        f.addRow("Заряд частицы, Кл:", self.charge)
        self.box_new = box_new

        hint = QLabel("Скорость по умолчанию — рабочая для проекта (суррогатная энергия γ≈7): "
                      "трекинг ведётся на ней, а проектные 6 ГэВ используются для излучения.")
        hint.setWordWrap(True)
        hint.setStyleSheet("color:#666;")

        lay = QVBoxLayout(self)
        lay.addWidget(self.rb_exist)
        lay.addWidget(self.combo)
        lay.addWidget(self.rb_new)
        lay.addWidget(box_new, 1)
        lay.addWidget(hint)

        self.rb_new.toggled.connect(self._toggle)
        self._toggle()

    def _toggle(self) -> None:
        self.box_new.setEnabled(self.rb_new.isChecked())
        self.combo.setEnabled(self.rb_exist.isChecked())

    # --- результат
    def group_id(self) -> int:
        if self.rb_exist.isChecked() and self.combo.currentData() is not None:
            return int(self.combo.currentData())
        return projects.next_free(projects.particle_groups().keys(), 1)

    def build_rows(self, gid: int) -> list[list[str]]:
        """Строки particle_groups.csv для нового пучка (пусто, если выбран готовый)."""
        if not self.rb_new.isChecked():
            return []
        r = self.radius.value() * 1e-3
        x0, y0, z0 = self.x0.value(), self.y0.value(), self.z0.value()
        v = self.speed.value()
        m, q = self.mass.text().strip(), self.charge.text().strip()
        n = int(self.n_part.currentData())
        if n == 5:
            offs = [(0.0, 0.0), (-r, 0.0), (r, 0.0), (0.0, r), (0.0, -r)]
        else:
            offs = [(0.0, 0.0)] + [(r * math.cos(2 * math.pi * k / 12), r * math.sin(2 * math.pi * k / 12))
                                   for k in range(12)]
        return [[str(gid), str(i + 1), m, q, repr(x0), repr(y0 + dy), repr(z0 + dz), repr(v), "0.0", "0.0"]
                for i, (dy, dz) in enumerate(offs)]


# --------------------------------------------------------------------------- 4. загрузка в базу
class LoadPage(QWizardPage):
    """Идентификаторы решётки/опыта и флаги загрузки."""

    def __init__(self):
        super().__init__()
        self.setTitle("Загрузка в базу")
        self.setSubTitle("Решётка получит номер в общей базе; конфигурации сосуществуют по этому номеру.")

        in_db = projects.lattices_in_db()
        self.lattice = QSpinBox()
        self.lattice.setRange(1, 999)
        self.lattice.setValue(projects.next_free(in_db.keys(), 1))
        self.exp = QSpinBox()
        self.exp.setRange(1, 999)
        self.exp.setValue(projects.next_free([e["exp_id"] for e in projects.experiments()], 1))
        self.force = QCheckBox("Перезаписать решётку, если номер занят (--force, с резервными копиями)")
        self.assets = QCheckBox("Сгенерировать синтетические карты и STL (--gen-assets)")
        self.note = QLineEdit()
        self.note.setPlaceholderText("примечание к проекту (необязательно)")

        occupied = ", ".join(f"{k} ({v})" for k, v in sorted(in_db.items())) or "—"
        info = QLabel(f"Заняты решётки: {occupied}.\nПолная загрузка сбрасывает mat_* к проектным K из .lte.")
        info.setWordWrap(True)
        info.setStyleSheet("color:#666;")

        form = QFormLayout(self)
        form.addRow("Номер решётки:", self.lattice)
        form.addRow("Номер опыта:", self.exp)
        form.addRow("", self.force)
        form.addRow("", self.assets)
        form.addRow("Примечание:", self.note)
        form.addRow("", info)

        self.lattice.valueChanged.connect(self._sync_force)
        self._sync_force()

    def _sync_force(self) -> None:
        busy = self.lattice.value() in projects.lattices_in_db()
        self.force.setChecked(busy)


class NewProjectWizard(QWizard):
    """Собирает проект и список команд для его создания."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Новый проект")
        self.setWizardStyle(QWizard.ModernStyle)
        self.setOption(QWizard.NoBackButtonOnStartPage, True)
        self.resize(760, 560)
        self.p_lat = LatticePage()
        self.p_ring = RingPage()
        self.p_beam = BeamPage()
        self.p_load = LoadPage()
        for p in (self.p_lat, self.p_ring, self.p_beam, self.p_load):
            self.addPage(p)
        self.setButtonText(QWizard.FinishButton, "Создать проект")
        self.setButtonText(QWizard.NextButton, "Далее")
        self.setButtonText(QWizard.BackButton, "Назад")
        self.setButtonText(QWizard.CancelButton, "Отмена")

    # --- результат мастера ---------------------------------------------------
    def result_project(self) -> projects.Project:
        cfg = self.p_lat.cfg or {}
        path = self.p_ring.closed_path or self.p_lat.path.text()
        line = "" if self.p_ring.closed_path else self.p_lat.line.text()
        if not line and cfg.get("format") == "sila" and cfg.get("lines") and Path(path).is_file():
            line = cfg["lines"][0]
        return projects.Project(
            name=self.p_lat.name.text().strip(),
            path=path,
            line=line,
            fmt=cfg.get("format", "sila"),
            lattice=self.p_load.lattice.value(),
            exp=self.p_load.exp.value(),
            group=self.p_beam.group_id(),
            note=self.p_load.note.text().strip(),
        )

    def write_beam(self, gid: int) -> str:
        """Дописать новый пучок в particle_groups.csv. Возврат — сообщение для журнала."""
        rows = self.p_beam.build_rows(gid)
        if not rows:
            return ""
        path = paths.table("particle_groups.csv")
        with open(path, "a", encoding="utf-8", newline="") as f:
            csv.writer(f, delimiter=";").writerows(rows)
        return f"в particle_groups.csv добавлен пучок {gid}: {len(rows)} частиц"

    def load_jobs(self, proj: projects.Project) -> list[Job]:
        """Команды загрузки решётки в базу."""
        args = [proj.path]
        if proj.line:
            args += ["--line", proj.line]
        args += ["--lattice", proj.lattice, "--exp", proj.exp]
        if self.p_load.force.isChecked():
            args.append("--force")
        if self.p_load.assets.isChecked():
            args.append("--gen-assets")
        return [py_job(f"загрузка решётки {proj.lattice} в базу", "load_config.py", *args)]
