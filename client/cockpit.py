# -*- coding: utf-8 -*-
"""Кокпит бесконечного режима (режим 4): QMainWindow с 3D-сценой в ЦЕНТРЕ и ПЛАВАЮЩИМИ окнами-панелями.

Почему плавающие окна, а не MDI/доки: 3D-сцена на VTK — НАТИВНОЕ окно, оно перебивает z-порядок дочерних
Qt-виджетов (всплывает и прячет панели), а доки по краям вытесняют сцену. Поэтому сцена = центральный
виджет (всегда во всю область, доки её не теснят), а панель управления, дерево, инфо/графики/экраны —
ОТДЕЛЬНЫЕ верхнеуровневые окна (`QWidget` с `Qt.Window`, родитель — кокпит): их можно двигать куда угодно,
у каждого нативный крестик, и как самостоятельные ОС-окна они НЕ прячутся за сценой. UI в главном потоке
только ЧИТАЕТ снимок через Link (~30 к/с) и шлёт ручки; движок крутит модель в своём потоке; тяжёлый
TbtPlant строится в ФОНОВОМ потоке.

`format_meta` — чистая (без Qt), тестируется отдельно.
"""
from __future__ import annotations

import threading
from collections import Counter
from typing import Dict, List, Optional

from PySide6.QtCore import Qt, QTimer, QObject, Signal
from PySide6.QtGui import QCursor, QKeySequence, QShortcut
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QGridLayout, QGroupBox, QHBoxLayout, QSpinBox,
                               QInputDialog, QLabel, QMainWindow, QMessageBox, QPushButton, QToolBar,
                               QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget)

from .link import InProcLink, RUN, STOP, RESTART
from . import addresses as A
from . import layout as layoutmod

ENGINE_META = ("meta.model_time", "meta.wall_time", "meta.speed", "meta.cycle",
               "meta.running", "meta.turn", "meta.alive")
# размеры окон-подписчиков при открытии
_WIN_SIZES = {"info": (340, 330), "graph": (480, 320), "point": (400, 380), "screen": (560, 460),
              "orbit_modes": (820, 880)}


class _ToolWindow(QWidget):
    """Плавающее окно-инструмент (панель/инфо/график) — отдельное ВЕРХНЕУРОВНЕВОЕ окно над кокпитом:
    двигается куда угодно, с нативным крестиком, не прячется за сценой. closeEvent синхронно зовёт
    on_close (снять подписки/учёт сразу, не дожидаясь отложенного destroyed)."""

    def __init__(self, inner: QWidget, title: str, on_close=None, parent=None, on_top: bool = False):
        # Обычное окно (min/max/close, ЕДИНЫЙ стиль со всеми окнами); on_top добавляет «поверх всех»
        # (панель управления не уходит на задний план) — кнопка «свернуть» при этом сохраняется.
        flags = (Qt.Window | Qt.WindowStaysOnTopHint) if on_top else Qt.Window
        super().__init__(parent, flags)
        self.setWindowTitle(title)
        self.setAttribute(Qt.WA_DeleteOnClose, True)
        self._on_close = on_close
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(inner)

    def closeEvent(self, event):
        if self._on_close is not None:
            try:
                self._on_close()
            except Exception:
                pass
            self._on_close = None
        super().closeEvent(event)


# ------------------------------------------------------------------ чистое форматирование панели
def _fmt_time(s: float) -> str:
    s = float(s or 0.0)
    if s < 1e-3:
        return f"{s * 1e6:.1f} мкс"
    if s < 1.0:
        return f"{s * 1e3:.2f} мс"
    return f"{s:.3f} с"


def format_meta(meta: dict) -> dict:
    """Мета движка → строки для панели (чисто, без Qt — тестируется напрямую)."""
    m = meta or {}
    speed = float(m.get("speed", 0.0) or 0.0)
    return {
        "running": "БЕЖИТ" if m.get("running") else "ПАУЗА",
        "model_time": _fmt_time(m.get("model_time", 0.0)),
        "wall_time": _fmt_time(m.get("wall_time", 0.0)),
        "speed": (f"{speed:.3e} ×реал" if 0 < speed < 1e-2 else f"{speed:.3f} ×реал"),
        "cycle": str(int(m.get("cycle", 0))),
        "turn": str(int(m.get("turn", 0))),
        "alive": "жив" if m.get("alive", True) else "ПОТЕРЯН",
    }


# ------------------------------------------------------------------ панель управления
class ControlPanel(QWidget):
    """Пуск/стоп/рестарт + модельное/стенное время, скорость, цикл; ручка зазора ондулятора."""

    def __init__(self, link, parent=None):
        super().__init__(parent)
        self.link = link
        self._und: List[str] = []

        bar = QHBoxLayout()
        self.btn_run = QPushButton("▶ Пуск")
        self.btn_stop = QPushButton("⏸ Стоп")
        self.btn_restart = QPushButton("⟲ Рестарт")
        self.btn_run.clicked.connect(lambda: self.link.control(RUN))
        self.btn_stop.clicked.connect(lambda: self.link.control(STOP))
        self.btn_restart.clicked.connect(lambda: self.link.control(RESTART))
        for b in (self.btn_run, self.btn_stop, self.btn_restart):
            bar.addWidget(b)
        bar.addStretch(1)

        grid = QGridLayout()
        self.lbl = {}
        rows = [("running", "состояние"), ("model_time", "модельное время"), ("speed", "скорость"),
                ("wall_time", "стенное время"), ("cycle", "цикл"), ("turn", "оборот"), ("alive", "пучок")]
        for r, (key, title) in enumerate(rows):
            grid.addWidget(QLabel(title), r, 0)
            v = QLabel("—")
            v.setStyleSheet("font-weight:600;color:#185fa5;")
            v.setTextInteractionFlags(Qt.TextSelectableByMouse)
            self.lbl[key] = v
            grid.addWidget(v, r, 1)
        box_meta = QGroupBox("движок")
        box_meta.setLayout(grid)

        gap = QHBoxLayout()
        self.cmb = QComboBox()
        self.spin = QDoubleSpinBox()
        self.spin.setRange(1.0, 50.0)
        self.spin.setSingleStep(0.5)
        self.spin.setValue(5.0)
        self.spin.setSuffix(" мм")
        self.btn_gap = QPushButton("Задать зазор")
        self.btn_gap.clicked.connect(self._send_gap)
        gap.addWidget(QLabel("ондулятор"))
        gap.addWidget(self.cmb, 1)
        gap.addWidget(self.spin)
        gap.addWidget(self.btn_gap)
        box_gap = QGroupBox("ручка зазора")
        box_gap.setLayout(gap)

        # контуры САУ: мастер-переключатель + по контуру (имя файла логики), вкл/выкл
        self._loop_chk: dict = {}                      # имя контура → QCheckBox
        self.chk_all = QCheckBox("Все контуры")
        fa = self.chk_all.font(); fa.setBold(True); self.chk_all.setFont(fa)
        self.chk_all.setTristate(True)                 # показывает частичное состояние (часть вкл)
        self.chk_all.clicked.connect(self._toggle_all)
        self._loops_inner = QVBoxLayout()
        loops_lay = QVBoxLayout()
        loops_lay.addWidget(self.chk_all)
        loops_lay.addLayout(self._loops_inner)
        self.box_loops = QGroupBox("контуры САУ")
        self.box_loops.setLayout(loops_lay)
        self.box_loops.setVisible(False)               # показать, когда контуры заданы (set_loops)

        # оптика кольца: ЖИВЫЕ тюн/хроматизм/adts — клик по значению берёт его на слежение (график)
        self._opt_links: dict = {}                     # addr → ValueLink
        self._opt_grid = QGridLayout()
        self.box_optics = QGroupBox("оптика кольца")
        self.box_optics.setLayout(self._opt_grid)
        self.box_optics.setVisible(False)              # показать, когда оптика доступна (set_optics)

        # запись прогона (synchrotronsim.run.Run): каждый проход — от рестарта до лимита — своя папка
        # workspace/runs/<id> (паспорт + ряды); ЧТО писать — кнопка «● запись» в окнах графиков
        self.chk_rec = QCheckBox("записывать прогон")
        self.spin_turns = QSpinBox()
        self.spin_turns.setRange(0, 2_000_000_000)
        self.spin_turns.setSingleStep(1000)
        self.spin_turns.setSpecialValueText("без лимита")
        self.spin_turns.setSuffix(" об.")
        self.chk_repeat = QCheckBox("по лимиту — новый проход (иначе пауза)")
        self.lbl_run = QLabel("не пишется")
        self.lbl_run.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.lbl_run.setWordWrap(True)
        self.chk_rec.clicked.connect(self._send_record)
        self.chk_repeat.clicked.connect(self._send_record)
        self.spin_turns.editingFinished.connect(self._send_record)
        rec = QGridLayout()
        rec.addWidget(self.chk_rec, 0, 0, 1, 2)
        rec.addWidget(QLabel("лимит прохода"), 1, 0)
        rec.addWidget(self.spin_turns, 1, 1)
        rec.addWidget(self.chk_repeat, 2, 0, 1, 2)
        rec.addWidget(self.lbl_run, 3, 0, 1, 2)
        box_rec = QGroupBox("запись прогона")
        box_rec.setLayout(rec)
        box_rec.setToolTip("что писать — кнопка «● запись» в окнах графиков; папки — workspace/runs, индекс — "
                           "workspace/tables/runs.csv")

        lay = QVBoxLayout(self)
        lay.addLayout(bar)
        lay.addWidget(box_meta)
        lay.addWidget(box_rec)
        lay.addWidget(box_gap)
        lay.addWidget(self.box_optics)
        lay.addWidget(self.box_loops)
        lay.addStretch(1)

    def set_undulators(self, gids: List[str]):
        self._und = list(gids)
        self.cmb.clear()
        self.cmb.addItems(self._und)

    def _send_gap(self):
        gid = self.cmb.currentText()
        if gid:
            self.link.command(f"undulator.{gid}.gap", float(self.spin.value()))

    def update_meta(self, meta: dict):
        for key, text in format_meta(meta).items():
            if key in self.lbl:
                self.lbl[key].setText(text)
        n = len(meta.get("watches") or [])
        rid = meta.get("run_id") or ""
        if rid:
            self.lbl_run.setText(f"пишется: {rid} · величин {n} · проход {meta.get('pass', 0)}")
        else:
            self.lbl_run.setText(f"не пишется · отмечено величин {n}" if n else
                                 "не пишется · отметьте величины кнопкой «● запись» в окнах")

    def _send_record(self, *_):
        self.link.record(self.chk_rec.isChecked(), int(self.spin_turns.value()) or None, self.chk_repeat.isChecked())

    def set_record(self, on: bool, max_turns: Optional[int], repeat: bool):
        """Показать настройки записи (лимит/повтор из опыта, перенос при пересборке) — без отправки в движок."""
        for w in (self.chk_rec, self.spin_turns, self.chk_repeat):
            w.blockSignals(True)
        self.chk_rec.setChecked(bool(on))
        self.spin_turns.setValue(int(max_turns or 0))
        self.chk_repeat.setChecked(bool(repeat))
        for w in (self.chk_rec, self.spin_turns, self.chk_repeat):
            w.blockSignals(False)

    # ---- контуры САУ: список (имя файла логики) + переключатели вкл/выкл ----
    def set_loops(self, loops: List[dict]):
        """Заполнить панель контуров. `loops` — [{name, file, enabled, tip?}] из Scheduler.loops."""
        for c in self._loop_chk.values():              # снять прежние (при пересборке модели)
            self._loops_inner.removeWidget(c); c.setParent(None)
        self._loop_chk.clear()
        for info in loops:
            name = str(info["name"])
            chk = QCheckBox(f"{info.get('file', '?')}   ·   {name}")   # первично — имя файла логики контура
            chk.setChecked(bool(info.get("enabled", False)))
            chk.setToolTip(str(info.get("tip", name)))
            chk.clicked.connect(lambda checked, n=name: self._toggle_loop(n, checked))
            self._loop_chk[name] = chk
            self._loops_inner.addWidget(chk)
        self.box_loops.setVisible(bool(loops))
        self._refresh_all_state()

    def _toggle_loop(self, name: str, checked: bool):
        self.link.command(f"loop.{name}.enabled", bool(checked))
        self._refresh_all_state()

    def _toggle_all(self, *_):
        on = [c.isChecked() for c in self._loop_chk.values()]
        want = not (bool(on) and all(on))              # не все вкл → включить все; все вкл → выключить все
        for name, c in self._loop_chk.items():
            c.blockSignals(True); c.setChecked(want); c.blockSignals(False)
            self.link.command(f"loop.{name}.enabled", want)
        self._refresh_all_state()

    def _refresh_all_state(self):
        on = [c.isChecked() for c in self._loop_chk.values()]
        self.chk_all.blockSignals(True)
        self.chk_all.setCheckState(Qt.Checked if on and all(on)
                                   else Qt.Unchecked if not any(on) else Qt.PartiallyChecked)
        self.chk_all.blockSignals(False)

    def update_loops(self, values: dict):
        """Синхронизировать переключатели с живым состоянием контуров из снимка (без обратной отправки команд)."""
        if not self._loop_chk:
            return
        changed = False
        for name, c in self._loop_chk.items():
            v = values.get(f"loop.{name}.enabled")
            if v is None:
                continue
            if c.isChecked() != bool(v):
                c.blockSignals(True); c.setChecked(bool(v)); c.blockSignals(False)
                changed = True
        if changed:
            self._refresh_all_state()

    # ---- оптика кольца: тюн/хроматизм/adts на мониторинг ----
    def set_optics(self, specs: List[tuple], emit):
        """Заполнить бокс оптики. specs — [(addr, label, fmt)]; emit(kind, addr) — открыть график. Пусто → скрыть."""
        from .windows import ValueLink
        for i in reversed(range(self._opt_grid.count())):      # очистить прежние виджеты (пересборка модели)
            w = self._opt_grid.itemAt(i).widget()
            if w is not None:
                w.setParent(None)
        self._opt_links = {}
        for r, (addr, label, fmt) in enumerate(specs):
            self._opt_grid.addWidget(QLabel(label), r, 0)
            vl = ValueLink(addr, fmt=fmt, emit=emit)
            self._opt_links[addr] = vl
            self._opt_grid.addWidget(vl, r, 1)
        self.box_optics.setVisible(bool(specs))

    def update_optics(self, values: dict):
        """Обновить живые значения оптики из снимка (адрес есть в снимке только если подписан)."""
        for addr, vl in self._opt_links.items():
            if addr in values:
                vl.set_value(values[addr])


# ------------------------------------------------------------------ кокпит
class Cockpit(QMainWindow):
    """Кокпит: 3D-сцена в ЦЕНТРЕ + панель управления, дерево и окна-подписчики — ПЛАВАЮЩИМИ окнами поверх.
    Окна двигаются/закрываются, не прячутся за сценой. Подписки с ref-count: адрес жив, пока нужен окну."""

    def __init__(self, link, project: str = "sila26", parent=None, poll_ms: int = 33):
        super().__init__(parent)
        self.setWindowTitle("SynchrotronSim — режим 4 (живой кокпит)")
        self.resize(1360, 880)
        self.link = link
        self.project = project
        self.scene = None
        self._engine = None
        self._plant = None
        self._state = None
        self._builder = None
        self._last_seq = -1
        self._held_snap = None                 # снимок, отложенный для окон на время осмотра сцены
        self._subs: Counter = Counter()        # ref-count подписок на Link
        self._windows: List = []               # открытые SubscriberWindow
        self._ready = False
        self._panels_shown = False
        self.experiment = None                 # опыт (synchrotronsim.experiment.Experiment), который исполняет модель
        self._session_subs: List[str] = []     # подписки панели/сцены текущей сборки (снимаются при пересборке)
        self._carry: Optional[dict] = None     # что перенести в новую сборку (запись, отметки, зазоры, пуск)
        self._opts = dict(gap_scale=1.0, gap_slew_mm_s=500.0, publish_every=3, autorun=False, restore_last=True,
                          note="")

        # центр — 3D-сцена (всегда во всю область; доки её не теснят). До сборки — заглушка.
        self._placeholder = QLabel("строится модель кольца…")
        self._placeholder.setAlignment(Qt.AlignCenter)
        self._placeholder.setStyleSheet("color:#777;font-size:15px;")
        self.setCentralWidget(self._placeholder)

        # панель управления и дерево — ПЛАВАЮЩИЕ окна (показываем при первом showEvent)
        self.control = ControlPanel(self.link)
        self.win_ctl = _ToolWindow(self.control, "Управление", parent=self, on_top=True)  # всегда поверх сцены
        self.tree = QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.itemDoubleClicked.connect(self._tree_open)
        self.win_tree = _ToolWindow(self.tree, "Дерево проекта", parent=self)

        self._make_toolbar()
        self.statusBar().showMessage("инициализация…")
        # навигация камеры по элементам: Ctrl+←/→ — соседний, Shift (или Ctrl+Shift)+←/→ — соседний того же типа;
        # → — по движению пучка, ← — навстречу. На уровне приложения: работает, даже если фокус в плавающем окне
        # (поля ввода свои Ctrl/Shift+стрелки перехватывают сами — правка текста не ломается).
        self._nav_keys = []
        for seq, step, same in (("Ctrl+Right", 1, False), ("Ctrl+Left", -1, False),
                                ("Shift+Right", 1, True), ("Shift+Left", -1, True),
                                ("Ctrl+Shift+Right", 1, True), ("Ctrl+Shift+Left", -1, True)):
            sc = QShortcut(QKeySequence(seq), self)
            sc.setContext(Qt.ApplicationShortcut)
            sc.activated.connect(lambda st=step, sm=same: self._nav(st, sm))
            self._nav_keys.append(sc)

        self.timer = QTimer(self)
        self.timer.timeout.connect(self._tick)
        self.timer.start(int(poll_ms))

    def _make_toolbar(self):
        tb = QToolBar("окна")
        tb.setObjectName("toolbar_windows")
        self.addToolBar(tb)
        tb.addAction("Панель управления", lambda: self._toggle_win(self.win_ctl))
        tb.addAction("Дерево", lambda: self._toggle_win(self.win_tree))
        tb.addAction("Вид сверху", self._reset_view)
        tb.addSeparator()
        self.act_experiment = tb.addAction("Опыт…", self._experiment_dialog)   # точность/модель/контуры + пересборка
        self.act_experiment.setToolTip("опыт, который исполняет модель: правка точности и контуров с пересборкой, "
                                       "открыть/сохранить опыт (JSON)")
        tb.addSeparator()
        tb.addAction("Закрыть окна", self._close_subscriber_windows)
        tb.addSeparator()
        tb.addAction("Сохранить раскладку…", self._save_layout_dialog)
        tb.addAction("Загрузить раскладку…", self._load_layout_dialog)

    # ---- показ/размещение плавающих окон ----
    def showEvent(self, event):
        super().showEvent(event)
        if not self._panels_shown:                       # панели показываем один раз, у краёв поверх сцены
            self._panels_shown = True
            self.win_ctl.resize(330, 370)
            self._place_over(self.win_ctl, 12, 44)
            self.win_ctl.show()
            self.win_tree.resize(300, 340)
            self._place_over(self.win_tree, 12, 430)
            self.win_tree.show()

    def _place_over(self, w, dx, dy):
        g = self.frameGeometry()
        w.move(g.x() + int(dx), g.y() + int(dy))

    def _place_near_cursor(self, w):
        """Поставить плавающее окно у курсора (место клика), не вылезая за экран."""
        p = QCursor.pos()
        x, y = p.x() + 14, p.y() + 14
        scr = self.screen().availableGeometry() if self.screen() is not None else None
        if scr is not None:
            x = max(scr.x(), min(x, scr.right() - w.width()))
            y = max(scr.y(), min(y, scr.bottom() - w.height()))
        w.move(x, y)

    def _toggle_win(self, w):
        if w.isVisible():
            w.hide()
        else:
            w.show()
            w.raise_()

    def set_status(self, text: str):
        self.statusBar().showMessage(text)

    def bind_scene(self, scene):
        self.scene = scene
        self.setCentralWidget(scene)                      # сцена — ЦЕНТР: всегда во всю область
        self._placeholder = None
        try:
            scene.fit_view()
        except Exception:
            pass
        # окна открываем через QueuedConnection — ОТЛОЖЕННО, вне колбэка мыши VTK (иначе интерактор «залипал»)
        try:
            scene.picked.connect(self._open_info, Qt.QueuedConnection)            # ЛКМ по магниту → инфо
        except Exception:
            pass
        try:
            scene.hovered.connect(self._on_hover)         # наведение → имя магнита в статус-баре
        except Exception:
            pass
        try:
            scene.focused.connect(self.set_status)        # перелёт камеры (Ctrl/Shift+←/→) → элемент в статусе
        except Exception:
            pass
        try:
            scene.picked_screen.connect(self._open_screen_info, Qt.QueuedConnection)       # ЛКМ по экрану
        except Exception:
            pass
        try:
            scene.picked_undulator.connect(self._open_undulator_info, Qt.QueuedConnection)  # ЛКМ по ондулятору
        except Exception:
            pass
        try:
            scene.picked_sensor.connect(self._open_sensor_info, Qt.QueuedConnection)         # ЛКМ по сенсору
        except Exception:
            pass

    def _nav(self, step: int, same_type: bool):
        if self.scene is not None:
            try:
                self.scene.nav_step(step, same_type)
            except Exception as e:                       # noqa: BLE001 — навигация опциональна
                self.set_status(f"навигация: {type(e).__name__}: {e}")

    def _on_hover(self, name: str):
        if name:
            self.set_status(f"наведён: {name}")

    def _open_screen_info(self, screen_id: str):
        from .windows import ScreenInfoWindow
        self.set_status(f"экран: {screen_id}")
        self.open_window(ScreenInfoWindow(screen_id))

    def _open_undulator_info(self, gid: str):
        from .windows import UndulatorInfoWindow
        self.set_status(f"ондулятор: {gid}")
        self.open_window(UndulatorInfoWindow(gid))

    def _open_sensor_info(self, idx: int):
        from .windows import SensorInfoWindow
        s = None
        try:
            for sen in self.scene.model.sensors:
                if int(sen.get("idx", -1)) == int(idx):
                    s = sen
                    break
        except Exception:
            s = None
        if s is None:
            return
        pos = [s["xy"][0], s["xy"][1], float(s.get("z", 0.0))]
        self.set_status(f"сенсор: {s.get('name', '')}")
        self.open_window(SensorInfoWindow(int(idx), name=s.get("name", ""), pos=pos, sid=s.get("sid", ""),
                                          meas_noise=s.get("meas_noise")))

    def _reset_view(self):
        if self.scene is not None:
            self.scene.fit_view()

    # ---- менеджер подписок (ref-count) ----
    def _subscribe(self, addr: str):
        if self._subs[addr] == 0:
            self.link.subscribe(addr)
        self._subs[addr] += 1

    def _unsubscribe(self, addr: str):
        if self._subs[addr] <= 0:
            return
        self._subs[addr] -= 1
        if self._subs[addr] == 0:
            self.link.unsubscribe(addr)
            del self._subs[addr]

    # ---- окна-подписчики (плавающие) ----
    def open_window(self, win, geom=None, restore: bool = False):
        """Открыть окно-подписчик плавающим окном у курсора, подписать адреса, связать сигналы. Дубль того
        же окна не плодим — поднимаем уже открытое. restore=True — восстановление раскладки (geom задан)."""
        from .windows import ScreenWindow, UndulatorInfoWindow
        d0 = win.to_dict()
        if geom is None and not restore and self._already_open(d0.get("kind"), d0.get("params")):
            return None
        for a in win.addresses():
            self._subscribe(a)
        try:
            win.request_open.connect(self._open_from_request)
            win.request_watch.connect(self._watch_from_window)
        except Exception:
            pass
        if isinstance(win, (ScreenWindow, UndulatorInfoWindow)):
            win.bind_link(self.link)
        tw = _ToolWindow(win, win.title(), on_close=lambda w=win: self._on_window_closed(w), parent=self)
        win._toolwin = tw
        self._windows.append(win)
        if geom:
            tw.setGeometry(*geom)
        else:
            tw.resize(*_WIN_SIZES.get(getattr(win, "KIND", ""), (420, 320)))
            self._place_near_cursor(tw)
        tw.show()
        tw.raise_()
        return tw

    def _already_open(self, kind, params) -> bool:
        for w in list(self._windows):
            d = w.to_dict()
            if d.get("kind") == kind and d.get("params") == params:
                tw = getattr(w, "_toolwin", None)
                if tw is not None:
                    tw.show()
                    tw.raise_()
                return True
        return False

    def _on_window_closed(self, win):
        if win in self._windows:
            self._windows.remove(win)
            for a in win.addresses():
                self._unsubscribe(a)

    def _close_subscriber_windows(self):
        for win in list(self._windows):
            tw = getattr(win, "_toolwin", None)
            if tw is not None:
                tw.close()

    def _watch_from_window(self, addr: str, on: bool):
        """Кнопка «● запись» окна → писать/не писать адрес в прогон (Watch движка)."""
        (self.link.watch if on else self.link.unwatch)(addr)

    def _open_from_request(self, kind: str, addr: str):
        if kind == "graph":
            self._open_graph(addr)
        elif kind == "screen":
            self._open_screen(addr)                       # кнопка «тепловая карта» → окно-экран
        elif kind == "point":
            ax, _, ay = addr.partition("|")               # кнопка «график x×y» сенсора → окно-точка
            self._open_point(ax, ay)
        elif kind == "diag":
            self._open_diag(addr)

    def _open_diag(self, name: str):
        """Стандартный график по имени (`diag.<name>` или `<name>`) → его окно из DIAG_WINDOWS."""
        from .windows import DIAG_WINDOWS
        cls = DIAG_WINDOWS.get(name.split(".", 1)[1] if name.startswith("diag.") else name)
        if cls is not None:
            self.open_window(cls())

    def _open_graph(self, addr: str):
        from .windows import GraphWindow
        self.open_window(GraphWindow(addr))

    def _open_info(self, addr: str, name: str = ""):
        from .windows import InfoWindow
        element = name or (addr[len("element."):] if addr.startswith("element.") else addr)
        kind, aids = "", None                              # вид и актуаторы объекта — из геометрии сцены
        try:
            for n in (self.scene.model.named or []):
                if n.get("name") == element:
                    kind, aids = n.get("kind", ""), n.get("aids")
                    break
        except Exception:
            pass
        self.set_status(f"элемент: {element}")
        self.open_window(InfoWindow(element, kind=kind, aids=aids))

    def _open_point(self, x_addr: str, y_addr: str):
        from .windows import PointWindow
        self.open_window(PointWindow(x_addr, y_addr))

    def _open_screen(self, screen_id: str):
        from .windows import ScreenWindow
        self.open_window(ScreenWindow(screen_id))

    # ---- дерево проекта ----
    def build_tree(self, state, loops, screens: List[str]):
        self.tree.clear()

        def node(parent, text, payload=None):
            it = QTreeWidgetItem([text])
            if payload is not None:
                it.setData(0, Qt.UserRole, payload)
            (parent.addChild if isinstance(parent, QTreeWidgetItem) else self.tree.addTopLevelItem)(it)
            return it

        n_obs = node(None, "наблюдаемые")
        for oid in sorted(state.observables):
            node(n_obs, oid, ("graph", f"obs.{oid}:rms"))
        n_loop = node(None, "контуры")
        for name, spec in (loops or {}).items():
            node(n_loop, f"{name} [{'вкл' if spec.enabled else 'выкл'}]", ("graph", f"loop.{name}.enabled"))
        n_und = node(None, "ондуляторы")
        for gid in sorted(state.undulators):
            node(n_und, gid, ("graph", f"undulator.{gid}.gap"))
        n_scr = node(None, "экраны")
        sub = {}
        for it in screens:                                 # sid или (sid, вид): ондуляторы / линейки диполей
            sid, kind = (it, "") if isinstance(it, str) else (it[0], it[1])
            parent = n_scr
            if kind:
                title = {"undulator": "у ондуляторов", "bend": "у линеек диполей (модуль №3)"}.get(kind, kind)
                if kind not in sub:
                    sub[kind] = node(n_scr, title)
                parent = sub[kind]
            node(parent, sid, ("screen", sid))
        n_sens = node(None, "сенсоры")
        if "ebpm.x" in state.observables:
            node(n_sens, "e-BPM 0 (x×y)", ("point", "obs.ebpm.x:0", "obs.ebpm.y:0"))
        n_eng = node(None, "движок")
        for a, t in (("meta.speed", "скорость"), ("meta.model_time", "модельное время"),
                     ("meta.cycle", "цикл")):
            node(n_eng, t, ("graph", a))
        n_std = node(None, "стандартные графики")
        if "fofb_x" in (loops or {}):                      # η-узор и моды FOFB — нужна матрица отклика FOFB X
            node(n_std, A.DIAG_TITLES["orbit_modes"], ("diag", "orbit_modes"))
        self.tree.expandToDepth(0)

    def _tree_open(self, item, _col):
        payload = item.data(0, Qt.UserRole)
        if not payload:
            return
        kind = payload[0]
        if kind == "graph":
            self._open_graph(payload[1])
        elif kind == "info":
            self._open_info(f"element.{payload[1]}", payload[1])
        elif kind == "screen":
            self._open_screen(payload[1])
        elif kind == "point":
            self._open_point(payload[1], payload[2])
        elif kind == "diag":
            self._open_diag(payload[1])

    # ---- такт опроса ----
    def _scene_busy(self) -> bool:
        try:
            return bool(self.scene is not None and self.scene.is_interacting())
        except Exception:
            return False

    def _feed_windows(self, snap):
        for win in list(self._windows):
            try:
                win.on_snapshot(snap)
            except Exception:
                pass

    def _tick(self):
        if self.scene is not None:
            try:
                self.scene.animate_jaws()         # плавная анимация челюстей ондулятора — каждый кадр
            except Exception:
                pass
        busy = self._scene_busy()
        if not busy and self._held_snap is not None:   # осмотр закончился — окна догоняют последним снимком
            held, self._held_snap = self._held_snap, None
            self._feed_windows(held)
        snap = self.link.poll_snapshot()
        if snap is None:
            return
        seq = int(snap.meta.get("seq", -1))
        if seq == self._last_seq:
            return
        self._last_seq = seq
        self.control.update_meta(snap.meta)
        self.control.update_loops(snap.values)        # синхрон переключателей контуров с движком
        self.control.update_optics(snap.values)       # живые тюн/хроматизм/adts на панели
        if self.scene is not None:
            self.scene.update_snapshot(snap)
        if busy:                                      # крутят камеру: перерисовка графиков (десятки мс в этом же
            self._held_snap = snap                    # потоке) дала бы рывки — откладываем до конца осмотра
            return
        self._feed_windows(snap)

    # ---- раскладка (этап 6): окна + панели по глоб. геометрии + камера ----
    @staticmethod
    def _wgeom(w) -> list:
        return [w.x(), w.y(), w.width(), w.height()]

    def collect_layout(self) -> dict:
        wins = []
        for win in self._windows:
            tw = getattr(win, "_toolwin", None)
            d = dict(win.to_dict())
            d["geom"] = (self._wgeom(tw) if tw is not None else None)
            wins.append(d)
        cam = None
        if self.scene is not None:
            try:
                cam = [list(map(float, row)) for row in self.scene.camera_position]
            except Exception:
                cam = None
        panels = {"control": {"geom": self._wgeom(self.win_ctl), "visible": self.win_ctl.isVisible()},
                  "tree": {"geom": self._wgeom(self.win_tree), "visible": self.win_tree.isVisible()}}
        return layoutmod.make_layout(wins, camera=cam, docks=panels)

    def restore_layout(self, data: dict, with_camera: bool = True):
        if not data:
            return
        from .windows import make_window
        for win in list(self._windows):       # закрыть прежние окна перед восстановлением
            tw = getattr(win, "_toolwin", None)
            if tw is not None:
                tw.close()
        for d in data.get("windows", []):
            win = make_window(d)
            if win is not None:
                self.open_window(win, geom=d.get("geom"), restore=True)
        panels = data.get("docks") or {}
        for key, w in (("control", self.win_ctl), ("tree", self.win_tree)):
            info = panels.get(key)
            if isinstance(info, dict):
                g = info.get("geom")
                if g:
                    w.setGeometry(*g)
                w.setVisible(bool(info.get("visible", True)))
        cam = data.get("camera")
        if with_camera and cam and self.scene is not None:
            try:
                self.scene.camera_position = cam
                self.scene.render()
            except Exception:
                pass

    def _save_layout_dialog(self):
        name, ok = QInputDialog.getText(self, "Сохранить раскладку", "Имя раскладки:", text="моя")
        if ok and name:
            p = layoutmod.save_layout(self.project, name, self.collect_layout())
            self.set_status(f"раскладка сохранена: {p}")

    def _load_layout_dialog(self):
        names = layoutmod.list_layouts(self.project)
        if not names:
            QMessageBox.information(self, "Раскладки", "Сохранённых раскладок нет.")
            return
        name, ok = QInputDialog.getItem(self, "Загрузить раскладку", "Раскладка:", names, 0, False)
        if ok and name:
            data = layoutmod.load_layout(self.project, name)
            self.restore_layout(data)
            self.set_status(f"раскладка загружена: {name}")

    # ---- опыт: сборка модели по опыту, пересборка (окна, подписки и запись сохраняются) ----
    def build(self, experiment):
        """Собрать модель по опыту в ФОНЕ (первый запуск или пересборка). Готовая сборка — `_attach`."""
        from synchrotronsim.experiment import Experiment
        exp = load_experiment(experiment) if not isinstance(experiment, Experiment) else experiment
        if exp.plant.kind not in COCKPIT_PLANTS:
            raise ValueError(f"опыт «{exp.name}»: модель {exp.plant.kind} без кольца — кокпит показывает только "
                             f"{'/'.join(COCKPIT_PLANTS)}; такой опыт исполняется Simulation(exp).run()")
        note = self._opts.get("note") or ""
        self.set_status(f"строится модель кольца: опыт «{exp.name}», решётка {exp.lattice_id} — это ~15–20 с…"
                        + (f" · {note}" if note else ""))
        if self.scene is None and self._placeholder is not None:
            self._placeholder.setText(f"строится модель кольца…\nопыт «{exp.name}», решётка {exp.lattice_id}")
        b = _Builder(exp)
        self._builder = b
        b.done.connect(self._attach)
        b.start()
        return b

    def rebuild(self, experiment):
        """Пересобрать модель по новому опыту: движок останавливается (проход записи закрывается), plant
        закрывается, новая сборка — в фоне; окна-подписчики, отметки «● запись», настройки записи, зазоры
        ондуляторов и режим пуск/пауза переносятся."""
        if self._builder is not None:
            self.set_status("модель ещё строится — дождитесь готовности")
            return None
        from synchrotronsim.experiment import Experiment
        exp = load_experiment(experiment) if not isinstance(experiment, Experiment) else experiment
        if exp.plant.kind not in COCKPIT_PLANTS:                 # проверить ДО остановки текущей модели
            QMessageBox.warning(self, "Опыт", f"модель {exp.plant.kind} без кольца — кокпит её не показывает")
            return None
        self._carry = self._detach()
        self._opts["note"] = ""
        return self.build(exp)

    def _detach(self) -> dict:
        """Остановить движок и закрыть plant текущей сборки; вернуть то, что переносится в следующую."""
        eng, self._engine = self._engine, None
        carry: dict = dict(self._carry or {})
        if eng is not None:
            R = eng.runner
            carry.update(watches=set(R.watches), rec_on=bool(R.rec_on), max_turns=R.max_turns, repeat=bool(R.repeat),
                         running=bool(eng.running), screens=dict(eng._screens))
            eng.stop()
            eng.join(timeout=15)
        if self._state is not None:
            carry["gaps"] = {gid: (float(u.gap_mm), float(u.gap_target), float(u.slew_mm_s))
                             for gid, u in (getattr(self._state, "undulators", {}) or {}).items()}
        if self._plant is not None:
            try:
                self._plant.close()
            except Exception:
                pass
        self._plant = self._state = None
        for a in self._session_subs:
            self._unsubscribe(a)
        self._session_subs = []
        self._ready = False
        return carry

    def _sub_session(self, addr: str):
        self._subscribe(addr)
        self._session_subs.append(addr)

    def _attach(self, result):
        """Готовая сборка (plant, state, sch) → сцена, панель, дерево, подписки, движок; перенос из прежней."""
        self._builder = None
        if isinstance(result, Exception):
            self.set_status(f"ошибка сборки модели: {type(result).__name__}: {result}")
            if self._placeholder is not None:
                self._placeholder.setText(f"не удалось построить модель:\n{result}")
            return
        plant, state, sch = result
        exp = getattr(sch, "experiment", None)
        first = self.scene is None
        carry, self._carry = (self._carry or {}), None
        self.experiment = exp
        self._plant, self._state = plant, state
        o = self._opts
        init_undulators(state, o["gap_slew_mm_s"])      # один ондулятор сведён (5 мм), остальные разведены (50 мм)
        for gid, (g, gt, slew) in (carry.get("gaps") or {}).items():   # пересборка: зазоры как были
            u = state.undulators.get(gid)
            if u is not None:
                u.slew_mm_s = slew
                u.gap_max = max(float(getattr(u, "gap_max", 20.0)), g, gt)
                u.gap_mm, u.gap_target = g, gt
        from .scene3d import make_ring_scene, build_geometry
        scene = make_ring_scene(gap_scale=o["gap_scale"]) if first else self.scene
        scene.set_geometry(build_geometry(plant, state))
        if first:
            self.bind_scene(scene)

        gids = list(state.undulators)
        self.control.set_undulators(gids)
        self.control.set_loops(_loops_info(sch.loops))   # панель контуров САУ (имя файла логики + вкл/выкл)
        opt_specs = _optics_specs(state)                 # оптика кольца (тюн/хроматизм/adts) на панель + мониторинг
        self.control.set_optics(opt_specs, self._open_from_request)
        for addr, _, _ in opt_specs:
            self._sub_session(addr)
        for a in ENGINE_META:
            self._sub_session(a)
        for gid in gids:
            self._sub_session(f"undulator.{gid}.gap")
            self._sub_session(f"undulator.{gid}.gap_target")
        for name in sch.loops:                           # живое состояние контуров → переключатели панели
            self._sub_session(f"loop.{name}.enabled")

        scr_list = [(s["sid"], s.get("kind", "")) for s in (scene.model.screens or [])] or _read_screen_ids()
        self.build_tree(state, sch.loops, scr_list)

        from .engine import EngineThread
        spc = int(getattr(getattr(exp, "run", None), "steps_per_cycle", 9) or 9)
        eng = EngineThread(self.link, sch, steps_per_cycle=spc, publish_every=o["publish_every"])
        R = eng.runner
        for w in (getattr(exp, "record", None) or []):   # подписка опыта (адреса) → отмеченные величины
            if not w.recorder:
                R.watch(w.addr)
        for a in carry.get("watches") or ():
            R.watch(a)
        eng._screens.update(carry.get("screens") or {})
        if carry:
            rec_on, max_turns, repeat = carry.get("rec_on", False), carry.get("max_turns"), carry.get("repeat", False)
        else:
            run = getattr(exp, "run", None)
            rec_on, max_turns, repeat = False, getattr(run, "max_turns", None), bool(getattr(run, "repeat", False))
        R.set_recording(rec_on, max_turns, repeat)       # до старта потока: запись нового опыта — новый прогон
        self.control.set_record(rec_on, max_turns, repeat)
        eng._running = bool(carry.get("running", False))
        self._engine = eng
        eng.start()
        self._ready = True
        if first:
            if o["restore_last"]:
                data = layoutmod.load_last(self.project)
                if data:
                    self.restore_layout(data, with_camera=False)   # не тащить устаревшую камеру
            try:
                self.scene.fit_view()                              # всегда стартуем с чистого вида сверху
            except Exception:
                pass
        name = f"опыт «{exp.name}»" if exp is not None else "модель"
        if first and o["autorun"]:
            self.link.control(RUN)
            self.set_status(f"{name} готов — идёт прогон")
        elif eng._running:
            self.set_status(f"{name} пересобран — прогон продолжается с нуля")
        else:
            self.set_status(f"{name} {'готов' if first else 'пересобран'} — нажмите «Пуск»")

    def current_experiment(self):
        """Опыт с ЖИВЫМИ флагами контуров, отмеченными величинами и лимитом записи — то, что сохранится в JSON
        и с чего начнётся пересборка."""
        import dataclasses
        from synchrotronsim.run import Watch
        exp = self.experiment
        eng = self._engine
        if exp is None or eng is None:
            return exp
        live = eng.sch.loops
        loops = [dataclasses.replace(l, enabled=bool(live[l.name].enabled) if l.name in live else l.enabled)
                 for l in exp.loops]
        R = eng.runner
        record = [w for w in exp.record if w.recorder] + [Watch(a, cadence="per_step") for a in sorted(R.watches)]
        run = dataclasses.replace(exp.run, max_turns=R.max_turns, repeat=bool(R.repeat))
        return dataclasses.replace(exp, loops=loops, record=record, run=run)

    def _experiment_dialog(self):
        if self.experiment is None:
            self.set_status("модель ещё строится — опыт появится после сборки")
            return
        from .experiment_dialog import ExperimentDialog
        dlg = ExperimentDialog(self.current_experiment(), parent=self)
        if dlg.exec() and dlg.result_experiment() is not None:
            self.rebuild(dlg.result_experiment())

    def closeEvent(self, event):
        try:
            if self._ready:
                layoutmod.save_last(self.project, self.collect_layout())   # авто-`last` при выходе
        except Exception:
            pass
        try:
            if self._engine is not None:
                self._engine.stop()
                self._engine.join(timeout=5)
            if self._plant is not None:
                self._plant.close()
        except Exception:
            pass
        event.accept()


class _Builder(QObject):
    """Строит сессию опыта (plant + доска + планировщик — `synchrotronsim.simulation.build_session`) в фоновом
    потоке; `done` эмитится в поток-владелец. Опыт кладётся в планировщик (`sch.experiment`) — паспорт прогона."""
    done = Signal(object)   # (plant, state, sch) | Exception

    def __init__(self, experiment):
        super().__init__()
        self.experiment = experiment

    def start(self):
        threading.Thread(target=self._run, name="session-build", daemon=True).start()

    def _run(self):
        try:
            from synchrotronsim.simulation import build_session
            S = build_session(self.experiment)
            S.sch.experiment = self.experiment
            self.done.emit((S.plant, S.state, S.sch))
        except Exception as e:                      # noqa: BLE001 — ошибку доставим в UI
            self.done.emit(e)


COCKPIT_PLANTS = ("tbt", "composite")            # модели с кольцом (сцена, BPM, пучок); optics — без кольца
START_ON = ("fofb_x", "fofb_y", "tune", "chroma")  # контуры, включённые при старте кокпита (adts ВЫКЛ: октуполи ранг-1)


def load_experiment(spec):
    """Опыт для кокпита: Experiment | путь к JSON (`Experiment.to_json`) | «fb:N» — опыт ОС из таблиц fb_*."""
    from synchrotronsim.experiment import Experiment
    if isinstance(spec, Experiment):
        return spec
    s = str(spec).strip()
    if s.lower().startswith("fb:"):
        return Experiment.from_fb(int(s[3:]))
    return Experiment.from_json(s)


def cockpit_experiment(lattice: int = 26, settle_turns: int = 250, steps_per_cycle: int = 9, nsl: int = 8,
                       n_part: int = 12, stride: int = 15, fofb_on: bool = True, optics: bool = True,
                       nsl_dev="auto"):
    """Опыт кокпита по параметрам запуска → (Experiment, пояснение к срезам). Те же параметры, что прежде
    передавались сборщику: BPM — отклонение от опорной орбиты, контур ВЧ (стартует ВЫКЛ), оптика по флагу;
    `fofb_on` — контуры START_ON включены на старте. nsl_dev="auto" — набор workspace/configs/slices_<lat>.json."""
    from synchrotronsim.experiment import Experiment
    note = ""
    if isinstance(nsl_dev, str) and nsl_dev == "auto":
        from synchrotronsim.convergence import cockpit_slices
        nsl_dev, note = cockpit_slices(lattice)
    exp = Experiment.cockpit(lattice=lattice, settle_turns=settle_turns, steps_per_cycle=steps_per_cycle,
                             nsl=nsl, n_part=n_part, stride=stride,       # nsl — главный рычаг модельной скорости
                             **({"nsl_dev": dict(nsl_dev)} if nsl_dev else {}),   # срезы по типам (--slices)
                             relative_orbit=True,                         # BPM — отклонение от опорной орбиты
                             optics=optics,                               # +оптические контуры (tune/chroma/adts)
                             rf=True)                                     # +контур ВЧ-частоты (энергия), ВЫКЛ
    if fofb_on:
        for l in exp.loops:
            if l.name in START_ON:
                l.enabled = True
    return exp, note


def _loops_info(loops: dict) -> List[dict]:
    """Контуры САУ для панели: [{name, file, enabled, tip}]. `file` — имя файла, описывающего логику контура
    (модуль класса контроллера, напр. integral_svd.py); разные контуры одного типа различаются по name."""
    import inspect
    import os
    out = []
    for name, spec in (loops or {}).items():
        ctrl = getattr(spec, "controller", None)
        try:
            file = os.path.basename(inspect.getsourcefile(type(ctrl)) or inspect.getfile(type(ctrl)))
        except (TypeError, OSError):
            file = "?"
        cls = type(ctrl).__name__ if ctrl is not None else "?"
        rate = float(getattr(spec, "rate_hz", 0.0) or 0.0)
        tip = f"{cls}  ·  {rate:.0f} Гц" if rate else cls
        out.append({"name": name, "file": file, "enabled": bool(getattr(spec, "enabled", False)), "tip": tip})
    return out


_OPTICS_COMPONENTS = {
    "tune":   [("obs.tune:0", "тюн νx"), ("obs.tune:1", "тюн νy")],
    "chroma": [("obs.chroma:0", "хроматизм ξx"), ("obs.chroma:1", "хроматизм ξy")],
    "adts":   [("obs.adts:0", "adts ∂νx/∂Jx"), ("obs.adts:1", "adts ∂νy/∂Jy"),
               ("obs.adts:2", "adts ∂νx/∂Jy")],
}


def _optics_specs(state) -> List[tuple]:
    """Строки оптики для панели: [(addr, label, fmt)] по КОМПОНЕНТАМ тюна/хроматизма/adts — только для тех
    наблюдаемых, что реально есть в состоянии (optics=True). Пусто при orbit-only (optics=False)."""
    obs = getattr(state, "observables", {}) or {}
    specs = []
    for o in ("tune", "chroma", "adts"):
        if o in obs:
            for addr, label in _OPTICS_COMPONENTS[o]:
                specs.append((addr, label, "{:.4f}"))
    return specs


def init_undulators(state, slew, open_gid="und_01", open_mm=5.0, closed_mm=50.0):
    """Стартовая конфигурация ондуляторов: ВКЛЮЧЁН один (`open_gid`, челюсти сведены на `open_mm`), остальные
    ВЫКЛЮЧЕНЫ (челюсти разведены на `closed_mm` — поток ∝ K² падает ~1e-7). Поднимаем `gap_max` (деф. 20),
    иначе `command/advance` клипнут 50 обратно; ставим gap_mm==gap_target, чтобы `advance()` не уводил старт."""
    for gid, u in (getattr(state, "undulators", {}) or {}).items():
        try:
            u.slew_mm_s = float(slew)
            u.gap_max = max(float(getattr(u, "gap_max", 20.0)), float(closed_mm))
            g = float(open_mm) if gid == open_gid else float(closed_mm)
            u.gap_mm = u.gap_target = g
        except Exception:
            pass


def _read_screen_ids() -> List[str]:
    """Список экранов из таблицы screens (колонка screen_id; через synchrotronsim.db); пусто, если таблицы нет."""
    try:
        from synchrotronsim.db import Repository
        return [v.strip() for v in Repository(layoutmod._ROOT).table("screens", required=False).column("screen_id")
                if v.strip()]
    except (OSError, ValueError):
        return []


def launch_cockpit(lattice: int = 26, exp=None, parent=None, gap_scale: float = 1.0,
                   steps_per_cycle: int = 9, publish_every: int = 3, settle_turns: int = 250,
                   nsl: int = 8, n_part: int = 12, stride: int = 15, fofb_on: bool = True,
                   optics: bool = True, autorun: bool = False, gap_slew_mm_s: float = 500.0,
                   project: Optional[str] = None, restore_last: bool = True,
                   nsl_dev="auto", experiment=None) -> Cockpit:
    """Открыть кокпит: создать Link и окно сразу, модель кольца собрать в фоне ПО ОПЫТУ, затем привязать сцену,
    запустить EngineThread, подписать панель/сцену, построить дерево и (опц.) восстановить раскладку `last`.
    Возвращает окно (держите ссылку, иначе соберётся GC).

    `experiment` — опыт (Experiment | путь к JSON | «fb:N»); без него опыт собирается из параметров запуска
    (`cockpit_experiment`: lattice, nsl, срезы nsl_dev, n_part, stride, settle_turns, optics, fofb_on). Опыт
    правится и пересобирается окном «Опыт…». `exp` — прежний хук оболочки (номер опыта проекта), не используется.

    `gap_slew_mm_s` — скорость привода зазора (мм/с МОДЕЛЬНОГО времени). Модельное время сжато ~1e-4× реала
    (6D-кольцо), поэтому физический привод (0.5 мм/с) за демо-время незаметен, а слишком быстрый пучок видит
    неадиабатично. Дефолт 500 — КОМПРОМИСС: ~540 оборотов/мм (в 10× адиабатичнее прежних 5000), видно за ~13 с
    стенных. --gap-slew: меньше → ближе к реальному адиабатич. приводу (можно 0.5); больше → быстрее/нагляднее.
    Сцена отражает РЕАЛИЗОВАННЫЙ зазор; FOFB парирует возмущение орбиты от смены зазора."""
    note = ""
    if experiment is not None:
        e = load_experiment(experiment)
    else:
        e, note = cockpit_experiment(lattice=lattice, settle_turns=settle_turns, steps_per_cycle=steps_per_cycle,
                                     nsl=nsl, n_part=n_part, stride=stride, fofb_on=fofb_on, optics=optics,
                                     nsl_dev=nsl_dev)
    project = project or f"sila{e.lattice_id}"
    link = InProcLink()
    win = Cockpit(link, project=project, parent=parent)
    win._opts.update(gap_scale=gap_scale, gap_slew_mm_s=gap_slew_mm_s, publish_every=publish_every,
                     autorun=autorun, restore_last=restore_last, note=note)
    win.show()
    try:
        win.build(e)
    except Exception as ex:                          # опыт без кольца и т.п. — в окно, не падать
        win.set_status(f"опыт не запущен: {type(ex).__name__}: {ex}")
        if win._placeholder is not None:
            win._placeholder.setText(f"опыт не запущен:\n{ex}")
    return win
