# -*- coding: utf-8 -*-
"""Окно «Опыт» кокпита: опыт (synchrotronsim.experiment.Experiment), который сейчас исполняет модель, —
просмотр и правка параметров точности/модели, включение контуров на старте, открыть/сохранить опыт (JSON).
«Применить и перестроить» возвращает НОВЫЙ опыт — кокпит пересобирает по нему модель (окна, подписки и запись
прогона сохраняются).

Опыт кокпита (source = cockpit) пересобирается через `Experiment.cockpit(**cockpit_params())` — поэтому правятся
и вид модели (оптические контуры, контур ВЧ); у опыта из файла/БД правятся аргументы plant'а и флаги контуров.
"""
from __future__ import annotations

import dataclasses
from typing import Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox, QFileDialog, QFormLayout,
                               QGroupBox, QLabel, QMessageBox, QPushButton, QSpinBox, QVBoxLayout, QHBoxLayout)

SLICES = (("keep", "как в опыте"), ("tight", "tight — набор по умолчанию (workspace/configs)"),
          ("loose", "loose — быстрее (nsl_convergence)"), ("none", "нет — только «срезов», полюса ИД 1"))


def _spin(lo, hi, val, step=1, suffix=""):
    s = QSpinBox()
    s.setRange(lo, hi)
    s.setSingleStep(step)
    s.setValue(int(val))
    if suffix:
        s.setSuffix(suffix)
    return s


def plant_default(key: str, fallback):
    """Значение аргумента TbtPlant по умолчанию (то, с чем plant собирается, если в опыте аргумента нет)."""
    try:
        import inspect
        from control.tbt_plant import TbtPlant
        v = inspect.signature(TbtPlant.__init__).parameters[key].default
        return fallback if v is inspect.Parameter.empty or v is None else v
    except Exception:
        return fallback


def slices_for(choice: str, lattice: int, current, same_lattice: bool = True):
    """Набор срезов по типам по выбору в окне: keep — как в опыте (при смене решётки — конфиг новой решётки:
    набор прежней к ней не относится); tight — конфиг кокпита; loose — из скана; none — без набора.
    Возврат (nsl_dev | None, пояснение)."""
    if choice == "keep" and not same_lattice:
        choice = "tight"
    if choice == "keep":
        return current, "как в опыте"
    if choice == "none":
        return None, "без набора"
    from synchrotronsim.convergence import cockpit_slices, load_recommendation
    if choice == "tight":
        return cockpit_slices(lattice)
    return load_recommendation(lattice, "loose"), "loose"


class ExperimentDialog(QDialog):
    """Опыт текущей модели. result_experiment() — опыт для пересборки (после «Применить и перестроить»)."""

    def __init__(self, experiment, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Опыт")
        self.base = experiment
        self._result = None
        e = experiment
        kw = e.plant.kwargs
        self.cockpit = (e.source == "cockpit")

        head = QLabel(f"<b>{e.name}</b> · модель <b>{e.plant.kind}</b> · решётка {e.lattice_id}"
                      + (f" · источник {e.source}" if e.source else "")
                      + (f"<br><i>{e.note}</i>" if e.note else ""))
        head.setTextFormat(Qt.RichText)
        head.setWordWrap(True)

        form = QFormLayout()
        self.sp_lat = _spin(1, 9999, e.lattice_id)
        self.sp_nsl = _spin(1, 200, kw.get("nsl", 8), suffix=" срезов")
        self.cb_slices = QComboBox()
        for key, title in SLICES:
            self.cb_slices.addItem(title, key)
        cur = kw.get("nsl_dev")
        self.lbl_slices = QLabel(f"сейчас: {len(cur)} типов" if cur else "сейчас: нет")
        # показанные значения: из опыта, иначе — умолчание plant'а (у опыта из файла/БД пишутся только правки)
        self._shown = {k: int(kw[k]) if k in kw else int(plant_default(k, d))
                       for k, d in (("nsl", 8), ("n_part", 12), ("stride", 15), ("settle_turns", 300))}
        self.sp_nsl.setValue(self._shown["nsl"])
        self.sp_part = _spin(1, 100000, self._shown["n_part"], suffix=" частиц")
        self.sp_stride = _spin(1, 100000, self._shown["stride"], suffix=" об./такт")
        self.sp_settle = _spin(0, 10_000_000, self._shown["settle_turns"], step=50, suffix=" об.")
        form.addRow("решётка", self.sp_lat)
        form.addRow("срезов (сектор/кик-тело)", self.sp_nsl)
        form.addRow("срезы по типам", self.cb_slices)
        form.addRow("", self.lbl_slices)
        form.addRow("частиц в пучке", self.sp_part)
        form.addRow("оборотов на такт", self.sp_stride)
        form.addRow("осадка при рестарте", self.sp_settle)
        self.chk_optics = QCheckBox("оптические контуры (тюн, хроматизм, ADTS) — модель composite")
        self.chk_rf = QCheckBox("контур частоты ВЧ")
        self.chk_optics.setChecked(e.plant.kind == "composite")
        self.chk_rf.setChecked(any(l.name == "rf" for l in e.loops))
        for c in (self.chk_optics, self.chk_rf):
            c.setEnabled(self.cockpit)
            if not self.cockpit:
                c.setToolTip("вид модели меняется только у опыта кокпита")
        form.addRow(self.chk_optics)
        form.addRow(self.chk_rf)
        box = QGroupBox("модель и точность")
        box.setLayout(form)

        self.loop_chk = {}
        lv = QVBoxLayout()
        for l in e.loops:
            c = QCheckBox(f"{l.name}  ({l.kind})")
            c.setChecked(bool(l.enabled))
            self.loop_chk[l.name] = c
            lv.addWidget(c)
        if not e.loops:
            lv.addWidget(QLabel("контуров нет"))
        box_loops = QGroupBox("контуры: включены на старте")
        box_loops.setLayout(lv)

        btns = QDialogButtonBox()
        self.btn_apply = btns.addButton("Применить и перестроить", QDialogButtonBox.AcceptRole)
        btns.addButton("Отмена", QDialogButtonBox.RejectRole)
        self.btn_apply.clicked.connect(self._apply)
        btns.rejected.connect(self.reject)
        io = QHBoxLayout()
        b_open = QPushButton("Открыть опыт…")
        b_save = QPushButton("Сохранить опыт…")
        b_open.clicked.connect(self._open)
        b_save.clicked.connect(self._save)
        io.addWidget(b_open); io.addWidget(b_save); io.addStretch(1)

        lay = QVBoxLayout(self)
        lay.addWidget(head)
        lay.addWidget(box)
        lay.addWidget(box_loops)
        lay.addLayout(io)
        lay.addWidget(btns)

    # ---- опыт из формы ----
    def experiment_from_form(self):
        """Новый опыт по значениям окна (опыт кокпита — через Experiment.cockpit, иначе — правка аргументов)."""
        from synchrotronsim.experiment import Experiment
        e = self.base
        lat = int(self.sp_lat.value())
        kw0 = e.plant.kwargs
        nsl_dev, _ = slices_for(self.cb_slices.currentData(), lat, kw0.get("nsl_dev"), lat == e.lattice_id)
        vals = {"nsl": self.sp_nsl, "n_part": self.sp_part, "stride": self.sp_stride, "settle_turns": self.sp_settle}
        # в опыт — то, что в нём уже было, и то, что правили (умолчания plant'а не вписываются)
        edits = {k: int(w.value()) for k, w in vals.items() if k in kw0 or int(w.value()) != self._shown[k]}
        if self.cockpit:
            p = e.cockpit_params()
            p.update(lattice=lat, optics=self.chk_optics.isChecked(), rf=self.chk_rf.isChecked(), **edits)
            p.pop("nsl_dev", None)
            if nsl_dev:
                p["nsl_dev"] = dict(nsl_dev)
            new = Experiment.cockpit(**p)
            new = dataclasses.replace(new, name=e.name, record=list(e.record), run=dataclasses.replace(e.run),
                                      seed=e.seed, note=e.note)
        else:
            kw = dict(kw0)
            kw.update(edits)
            kw.pop("nsl_dev", None)
            if nsl_dev:
                kw["nsl_dev"] = dict(nsl_dev)
            new = dataclasses.replace(e, lattice_id=lat, plant=dataclasses.replace(e.plant, kwargs=kw),
                                      loops=[dataclasses.replace(l) for l in e.loops])
        for l in new.loops:                                    # флаги «включён на старте» из окна
            if l.name in self.loop_chk:
                l.enabled = self.loop_chk[l.name].isChecked()
        return new

    def _apply(self):
        try:
            self._result = self.experiment_from_form()
        except Exception as ex:                                # напр. нет скана для loose
            QMessageBox.warning(self, "Опыт", f"{type(ex).__name__}: {ex}")
            return
        self.accept()

    def _open(self):
        path, _ = QFileDialog.getOpenFileName(self, "Открыть опыт", "", "Опыт (*.json)")
        if not path:
            return
        from synchrotronsim.experiment import Experiment
        try:
            self._result = Experiment.from_json(path)
        except Exception as ex:
            QMessageBox.warning(self, "Опыт", f"не удалось прочитать {path}:\n{type(ex).__name__}: {ex}")
            return
        self.accept()

    def _save(self):
        path, _ = QFileDialog.getSaveFileName(self, "Сохранить опыт", f"{self.base.name}.json", "Опыт (*.json)")
        if not path:
            return
        try:
            self.experiment_from_form().to_json(path)
        except Exception as ex:
            QMessageBox.warning(self, "Опыт", f"{type(ex).__name__}: {ex}")

    def result_experiment(self) -> Optional[object]:
        return self._result
