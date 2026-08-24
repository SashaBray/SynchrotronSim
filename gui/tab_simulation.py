# -*- coding: utf-8 -*-
"""Вкладка «Симуляция»: запуск трекинга разными ядрами, с излучением и без.

Три ядра проекта:
  • **C++ диффур** (`SynchrotronTracker`) — интегрирование силы Лоренца по картам поля, эталон.
    Излучение не моделирует — переключатели гасятся.
  • **C++ кик** (`MatrixTracker`) — автономное матричное ядро: само читает базу, трекает кольцо и
    пишет `.npz`. Излучение берёт ПРАВИЛОМ из `radiation_params.csv` по номеру опыта, флагов у него
    нет, поэтому переключатель излучения здесь правит именно эту строку базы (с резервной копией).
  • **Python кик** — эталонная реализация того же кикового ядра; принимает флаги излучения напрямую.
"""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QButtonGroup, QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout,
                               QGroupBox, QHBoxLayout, QLabel, QPushButton, QRadioButton, QSpinBox,
                               QVBoxLayout)

from . import paths, projects
from .csvtable import read_csv, write_csv
from .runner import exe_job, py_job
from .tabbase import TabBase
from .widgets import SectionLabel, plural

RAD_TABLE = "radiation_params.csv"
REC_GLOBAL, REC_LOCAL, REC_LONG = "global", "local", "long"   # что записывает трек


class SimulationTab(TabBase):
    def __init__(self, runner, log, parent=None):
        super().__init__(runner, log, parent)

        # --- опыт -------------------------------------------------------------
        self.exp = QComboBox()
        self.exp.setMinimumWidth(320)
        self.exp.currentIndexChanged.connect(self._exp_changed)
        self.btn_reload = QPushButton("Обновить список")
        self.btn_reload.clicked.connect(self.refresh)

        row_exp = QHBoxLayout()
        row_exp.addWidget(QLabel("Опыт:"))
        row_exp.addWidget(self.exp, 1)
        row_exp.addWidget(self.btn_reload)

        # Установка и пучок — из строки опыта, но их можно выбрать самому: конфигурация опыта
        # (включая финиш-плоскость) будет переписана под выбранную пару перед запуском.
        self.lattice = QComboBox()
        self.group = QComboBox()
        self.lattice.currentIndexChanged.connect(self._config_changed)
        self.group.currentIndexChanged.connect(self._config_changed)
        self.cfg_note = QLabel("")
        self.cfg_note.setWordWrap(True)
        self.cfg_note.setStyleSheet("color:#666;")

        row_cfg = QHBoxLayout()
        row_cfg.addWidget(QLabel("Установка:"))
        row_cfg.addWidget(self.lattice, 1)
        row_cfg.addWidget(QLabel("Пучок:"))
        row_cfg.addWidget(self.group, 1)

        # --- ядро -------------------------------------------------------------
        self.rb_diff = QRadioButton("C++ диффурное ядро (Boris по картам поля, эталон)")
        self.rb_cpp = QRadioButton("C++ киковое ядро (MatrixTracker, автономное)")
        self.rb_py = QRadioButton("Python киковое ядро (эталон кика, track3)")
        self.rb_diff.setChecked(True)
        self.engines = QButtonGroup(self)
        for i, rb in enumerate((self.rb_diff, self.rb_cpp, self.rb_py)):
            self.engines.addButton(rb, i)
            rb.toggled.connect(self._engine_changed)
        self.exe_info = QLabel("")
        self.exe_info.setStyleSheet("color:#666;")

        box_eng = QGroupBox("Ядро")
        v = QVBoxLayout(box_eng)
        for rb in (self.rb_diff, self.rb_cpp, self.rb_py):
            v.addWidget(rb)
        v.addWidget(self.exe_info)

        # --- параметры трека ---------------------------------------------------
        self.turns = QSpinBox()
        self.turns.setRange(1, 10_000_000)
        self.turns.setValue(1)
        self.record = QComboBox()
        self.record.addItem("глобальные траектории — для 3D-визуализации", REC_GLOBAL)
        self.record.addItem("локальные координаты — спектр, тюны, x(S)", REC_LOCAL)
        self.record.addItem("длинный турн-by-turn — затухание, огибающая", REC_LONG)
        self.record.setToolTip("Что записывает трек:\n"
                               "• глобальные — trajectories_exp<N>.bin (его читает visualize.py и HTML-сцена)\n"
                               "• локальные — trajectories_local_exp<N>.npz на границах элементов (traj_analyzer)\n"
                               "• длинный — состояние раз за оборот у плоскости; только тут работает ВЧ")
        self.record.currentIndexChanged.connect(lambda _i: self._engine_changed())
        self.x0 = QDoubleSpinBox(); self.x0.setRange(-1000, 1000); self.x0.setDecimals(3); self.x0.setSuffix(" мм")
        self.y0 = QDoubleSpinBox(); self.y0.setRange(-1000, 1000); self.y0.setDecimals(3); self.y0.setSuffix(" мм")
        self.single = QCheckBox("одна частица со смещением (вместо пучка группы)")
        self.single.toggled.connect(lambda on: (self.x0.setEnabled(on), self.y0.setEnabled(on)))
        self.x0.setEnabled(False)
        self.y0.setEnabled(False)

        box_par = QGroupBox("Трек")
        f = QFormLayout(box_par)
        f.addRow("Оборотов:", self.turns)
        f.addRow("Запись:", self.record)
        f.addRow("", self.single)
        f.addRow("Смещение x:", self.x0)
        f.addRow("Смещение y:", self.y0)

        # --- излучение ---------------------------------------------------------
        self.rad_on = QCheckBox("синхротронное излучение (потери энергии в изгибах)")
        self.rad_isr = QCheckBox("квантовое возбуждение ISR (равновесный эмиттанс)")
        self.rad_rf = QCheckBox("идеальный ВЧ-резонатор (компенсация потерь ⇒ затухание)")
        self.rad_mev = QDoubleSpinBox()
        self.rad_mev.setRange(1.0, 1e6)
        self.rad_mev.setDecimals(0)
        self.rad_mev.setValue(6000)
        self.rad_mev.setSuffix(" МэВ")
        self.rad_note = QLabel("")
        self.rad_note.setWordWrap(True)
        self.rad_note.setStyleSheet("color:#666;")
        self.rad_on.toggled.connect(self._rad_changed)
        self.rad_isr.setEnabled(False)
        self.rad_isr.setToolTip("требует включённого излучения")

        box_rad = QGroupBox("Излучение")
        fr = QFormLayout(box_rad)
        fr.addRow("", self.rad_on)
        fr.addRow("", self.rad_isr)
        fr.addRow("", self.rad_rf)
        fr.addRow("Проектная энергия:", self.rad_mev)
        fr.addRow("", self.rad_note)

        # --- запуск -------------------------------------------------------------
        self.btn_run = QPushButton("Запустить трекинг")
        self.btn_stop = QPushButton("Остановить")
        self.btn_stop.setEnabled(False)
        self.btn_run.clicked.connect(self._run)
        self.btn_stop.clicked.connect(lambda: self.window().log.btn_stop.click())   # общий обработчик остановки
        self.guard(self.btn_run)

        row_run = QHBoxLayout()
        row_run.addWidget(self.btn_run)
        row_run.addWidget(self.btn_stop)
        row_run.addStretch(1)

        cols = QHBoxLayout()
        cols.addWidget(box_par, 1)
        cols.addWidget(box_rad, 1)

        lay = QVBoxLayout(self)
        lay.addWidget(SectionLabel("Запуск симуляции"))
        lay.addLayout(row_exp)
        lay.addLayout(row_cfg)
        lay.addWidget(self.cfg_note)
        lay.addWidget(box_eng)
        lay.addLayout(cols)
        lay.addLayout(row_run)
        lay.addStretch(1)

        self._engine_changed()

    # ---- состояние -------------------------------------------------------
    def set_busy(self, busy: bool) -> None:
        super().set_busy(busy)
        self.btn_stop.setEnabled(busy)

    def refresh(self) -> None:
        cur = self.exp.currentData()
        self.exp.clear()
        for e in projects.experiments():
            self.exp.addItem(f"опыт {e['exp_id']} — решётка {e.get('lattice_id', '?')}, "
                             f"пучок {e.get('group_id', '?')}", int(e["exp_id"]))
        self.lattice.blockSignals(True)
        self.lattice.clear()
        for lid, n in sorted(projects.lattices_in_db().items()):
            self.lattice.addItem(f"решётка {lid} — {n} "
                                 + plural(n, "элемент", "элемента", "элементов")
                                 + (" (демо-конфигурация)" if n == 1 else ""), lid)
        self.lattice.blockSignals(False)
        self.group.blockSignals(True)
        self.group.clear()
        for gid, parts in sorted(projects.particle_groups().items()):
            p0 = parts[0]
            self.group.addItem(f"пучок {gid} — {len(parts)} частиц, старт "
                               f"({p0.get('x', '?')}, {p0.get('y', '?')}, {p0.get('z', '?')})", gid)
        self.group.blockSignals(False)
        if self.project is not None:
            idx = self.exp.findData(self.project.exp)
            if idx >= 0:
                self.exp.setCurrentIndex(idx)
        elif cur is not None:
            idx = self.exp.findData(cur)
            if idx >= 0:
                self.exp.setCurrentIndex(idx)
        self._exp_changed()
        self._check_exes()

    def _check_exes(self) -> None:
        miss = [s for s in ("SynchrotronTracker", "MatrixTracker") if not paths.exe_path(s).exists()]
        self.exe_info.setText("⚠ не собраны бинари: " + ", ".join(miss) +
                              "  (cmake --build build --config Release)" if miss
                              else "бинари на месте: SynchrotronTracker, MatrixTracker")
        self.exe_info.setStyleSheet("color:#b06000;" if miss else "color:#666;")

    def _current_exp(self) -> int | None:
        d = self.exp.currentData()
        return int(d) if d is not None else None

    # ---- излучение --------------------------------------------------------
    def _rad_row(self, exp: int) -> dict | None:
        p = paths.table(RAD_TABLE)
        if not p.exists():
            return None
        header, rows, _ = read_csv(p)
        for r in rows:
            d = dict(zip(header, r))
            if (d.get("exp_id") or "").strip() == str(exp):
                return d
        return None

    def _stored_config(self, exp: int) -> tuple[int | None, int | None]:
        """(установка, пучок) из строки опыта в experiments.csv."""
        return projects.experiment_config(exp)

    def _config_changed(self) -> None:
        """Подсказка: выбранная пара расходится со строкой опыта ⇒ она будет переписана."""
        exp = self._current_exp()
        if exp is None:
            return
        lat, grp = self._stored_config(exp)
        want_l, want_g = self.lattice.currentData(), self.group.currentData()
        if want_l is None or want_g is None:
            return
        if (want_l, want_g) == (lat, grp):
            self.cfg_note.setText("установка и пучок — как в строке опыта")
            self.cfg_note.setStyleSheet("color:#666;")
        else:
            self.cfg_note.setText(
                f"опыт {exp} будет ПЕРЕПИСАН под установку {want_l} и пучок {want_g} "
                f"(было: {lat}/{grp}); финиш-плоскость пересчитается по плоскости инжекции этого пучка. "
                f"Чтобы не трогать существующий опыт, выберите свободный номер в experiments.csv.")
            self.cfg_note.setStyleSheet("color:#b06000;")

    def _exp_changed(self) -> None:
        exp = self._current_exp()
        if exp is None:
            return
        lat, grp = self._stored_config(exp)               # установка/пучок — из строки выбранного опыта
        for combo, val in ((self.lattice, lat), (self.group, grp)):
            combo.blockSignals(True)
            i = combo.findData(val)
            if i >= 0:
                combo.setCurrentIndex(i)
            combo.blockSignals(False)
        self._config_changed()
        row = self._rad_row(exp)
        if row:
            self.rad_on.setChecked(True)
            self.rad_isr.setChecked((row.get("isr") or "0").strip() == "1")
            self.rad_rf.setChecked((row.get("rf") or "none").strip() == "ideal")
            try:
                self.rad_mev.setValue(float(row.get("design_energy_mev") or 6000))
            except ValueError:
                pass
        else:
            self.rad_on.setChecked(False)
        self._rad_changed()

    def _rad_changed(self) -> None:
        on = self.rad_on.isChecked()
        self.rad_isr.setEnabled(on and not self.rb_diff.isChecked())
        if not on:
            self.rad_isr.setChecked(False)
        self._engine_changed(update_rad_only=True)

    def _engine_changed(self, *_a, update_rad_only: bool = False) -> None:
        diff = self.rb_diff.isChecked()
        cpp = self.rb_cpp.isChecked()

        # доступность режимов записи: диффур пишет только глобальные, C++ кик — только свои два
        if diff:
            allowed = {REC_GLOBAL}
        elif cpp:
            allowed = {REC_LOCAL, REC_LONG}                          # C++ ядро пишет .npz, не глобальный .bin
        else:
            allowed = {REC_GLOBAL, REC_LOCAL, REC_LONG}
        model = self.record.model()
        for i in range(self.record.count()):
            model.item(i).setEnabled(self.record.itemData(i) in allowed)
        if self.record.currentData() not in allowed:                 # текущий стал недоступен → первый по списку
            first = next(i for i in range(self.record.count()) if self.record.itemData(i) in allowed)
            self.record.setCurrentIndex(first)

        long_mode = self.record.currentData() == REC_LONG
        for w in (self.rad_on, self.rad_mev):
            w.setEnabled(not diff)
        self.rad_rf.setEnabled(not diff and long_mode and self.rad_on.isChecked())
        if not self.rad_rf.isEnabled():
            self.rad_rf.setChecked(False)
        self.rad_isr.setEnabled(not diff and self.rad_on.isChecked())
        self.turns.setEnabled(not diff)
        self.record.setEnabled(not diff)
        self.single.setEnabled(not diff)
        if diff:
            self.rad_note.setText("Диффурное ядро излучение не моделирует — оно интегрирует движение "
                                  "по картам поля за один проход и пишет trajectories_exp<N>.bin.")
        elif cpp:
            self.rad_note.setText(f"C++ киковое ядро берёт излучение ПРАВИЛОМ из {RAD_TABLE} по номеру опыта. "
                                  f"Переключатели применятся к строке этого опыта (прежняя таблица уйдёт в копию).")
        else:
            self.rad_note.setText("Python киковое ядро принимает флаги излучения напрямую "
                                  "(--sr / --isr / --rf ideal / --design-mev); таблица базы не меняется.")

    def _apply_rad_row(self, exp: int) -> None:
        """Привести строку `radiation_params.csv` опыта в соответствие переключателям (для C++ ядра)."""
        p = paths.table(RAD_TABLE)
        if not p.exists():
            self.say(f"нет таблицы {RAD_TABLE} — излучение для C++ ядра задать нечем", True)
            return
        header, rows, enc = read_csv(p)
        idx = {h: i for i, h in enumerate(header)}
        keep = [r for r in rows if (r[idx["exp_id"]] if idx["exp_id"] < len(r) else "") != str(exp)]
        if self.rad_on.isChecked():
            new = [""] * len(header)
            new[idx["exp_id"]] = str(exp)
            new[idx["design_energy_mev"]] = str(int(self.rad_mev.value()))
            new[idx["isr"]] = "1" if self.rad_isr.isChecked() else "0"
            new[idx["rf"]] = "ideal" if self.rad_rf.isChecked() else "none"
            keep.append(new)
            what = (f"энергия {int(self.rad_mev.value())} МэВ, ISR "
                    f"{'вкл' if self.rad_isr.isChecked() else 'выкл'}, "
                    f"ВЧ {'ideal' if self.rad_rf.isChecked() else 'none'}")
        else:
            what = "строка удалена (излучение не учитывается)"
        keep.sort(key=lambda r: int(r[idx["exp_id"]]) if (r[idx["exp_id"]] or "").isdigit() else 0)
        bak = write_csv(p, header, keep, enc)
        self.say(f"{RAD_TABLE}, опыт {exp}: {what}  (копия прежней → {bak.name})")

    # ---- запуск -----------------------------------------------------------
    def _author_job(self, exp: int):
        """Если выбранная пара «установка+пучок» расходится со строкой опыта — переписать её ДО запуска.

        Конфигурация опыта (включая финиш-плоскость) — единый источник для всех трёх ядер; автономный
        C++ MatrixTracker читает только её, флагов установки/пучка у него нет. Строку собирает сам
        track3 (`--author-only`): финиш-плоскость считается по плоскости инжекции выбранного пучка.
        """
        want_l, want_g = self.lattice.currentData(), self.group.currentData()
        if want_l is None or want_g is None or (want_l, want_g) == self._stored_config(exp):
            return None
        return py_job(f"конфигурация опыта {exp}: установка {want_l}, пучок {want_g}",
                      "track3.py", "--author-only", "--exp", exp, "--lattice", want_l, "--group", want_g)

    def _run(self) -> None:
        exp = self._current_exp()
        if exp is None:
            self.say("не выбран опыт", True)
            return
        pre = self._author_job(exp)                        # привести строку опыта к выбранной конфигурации
        turns = self.turns.value()
        x0 = self.x0.value() if self.single.isChecked() else 0.0
        y0 = self.y0.value() if self.single.isChecked() else 0.0

        if self.rb_diff.isChecked():
            self.runner.run(*[j for j in (pre, py_job(f"диффурное ядро, опыт {exp}",
                                                      "track3.py", "--mode", "diff", "--exp", exp)) if j])
            return

        rec = self.record.currentData()
        if self.rb_cpp.isChecked():
            if not paths.exe_path("MatrixTracker").exists():
                self.say("нет build/MatrixTracker — соберите: cmake --build build --target MatrixTracker", True)
                return
            self._apply_rad_row(exp)                       # излучение C++ ядра задаётся только таблицей базы
            mode = "--full" if rec == REC_LOCAL else "--track"
            args = [mode, exp, turns] + ([x0, y0] if self.single.isChecked() else [])
            self.runner.run(*[j for j in (pre, exe_job(f"C++ киковое ядро, опыт {exp} ({mode})",
                                                      "MatrixTracker", *args)) if j])
            return

        args = ["--mode", "kick-phys", "--exp", exp, "--turns", turns]
        if rec == REC_LOCAL:
            args.append("--local")
        elif rec == REC_LONG:
            args.append("--long")
        if self.single.isChecked():
            args += ["--x0-mm", x0, "--y0-mm", y0]
        if self.rad_on.isChecked():
            args.append("--sr")
            if self.rad_isr.isChecked():
                args.append("--isr")
            args += ["--design-mev", int(self.rad_mev.value())]
        if self.rad_rf.isChecked():
            args += ["--rf", "ideal"]
        self.runner.run(*[j for j in (pre, py_job(f"Python киковое ядро, опыт {exp}",
                                                 "track3.py", *args)) if j])
