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
_WIN_SIZES = {"info": (340, 330), "graph": (480, 320), "point": (400, 380), "screen": (560, 460), "cavity": (1400, 900),
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

    def set_targets(self, rows: List[dict], emit):
        """Контуры: ТЕКУЩЕЕ → ЦЕЛЬ (client.engine.loop_targets): по строке на величину, которую держит контур
        (компоненты тюна/хроматизма/ADTS, СКО орбиты, энергия). Текущее и цель кликаются → график; «⇄» —
        один график: текущее и цель пунктиром. Пусто → скрыть."""
        from .windows import ValueLink, _LINK_CSS
        for i in reversed(range(self._opt_grid.count())):      # очистить прежние виджеты (пересборка модели)
            w = self._opt_grid.itemAt(i).widget()
            if w is not None:
                w.setParent(None)
        self._opt_links = {}
        self.box_optics.setTitle("контуры: текущее → цель")
        for c, t in enumerate(("", "текущее", "цель", "")):
            h = QLabel(f"<span style='color:#888'>{t}</span>")
            self._opt_grid.addWidget(h, 0, c)
        for r, row in enumerate(rows, start=1):
            un = row.get("units", "")
            sc = float(row.get("scale", 1.0))
            fmt = row.get("fmt", "{:.6g}")
            lab = QLabel(f"{row['label']} <span style='color:#888'>· {row['loop']}</span>")
            lab.setTextFormat(Qt.RichText)
            self._opt_grid.addWidget(lab, r, 0)
            cur = ValueLink(row["cur"], fmt=fmt, units=un, scale=sc, emit=emit)
            self._opt_links[row["cur"]] = cur
            self._opt_grid.addWidget(cur, r, 1)
            if row.get("tgt"):
                tg = ValueLink(row["tgt"], fmt=fmt, units=un, scale=sc, emit=emit)
                self._opt_links[row["tgt"]] = tg
                self._opt_grid.addWidget(tg, r, 2)
                ref = row["tgt"]
            else:
                v = float(row.get("tgt_value") or 0.0)
                tl = QLabel(fmt.format(v * sc) + (f" {un}" if un else "") + " <span style='color:#888'>(проект)</span>")
                tl.setTextFormat(Qt.RichText)
                tl.setToolTip("у контура нет уставки на доске — цель: проектная величина")
                self._opt_grid.addWidget(tl, r, 2)
                ref = repr(v)
            both = QLabel(f'<a href="graph_ref::{row["cur"]}|{ref}|{sc!r}|{un}" style="{_LINK_CSS}">⇄</a>')
            both.setTextFormat(Qt.RichText)
            both.setTextInteractionFlags(Qt.TextBrowserInteraction)
            both.setToolTip("график: текущее и цель")
            both.linkActivated.connect(lambda href: emit(*href.split("::", 1)))
            self._opt_grid.addWidget(both, r, 3)
        self.box_optics.setVisible(bool(rows))

    def update_optics(self, values: dict):
        """Обновить живые значения оптики из снимка (адрес есть в снимке только если подписан)."""
        for addr, vl in self._opt_links.items():
            if addr in values:
                vl.set_value(values[addr])


# ------------------------------------------------------------------ кокпит
class Cockpit(QMainWindow):
    """Кокпит: 3D-сцена в ЦЕНТРЕ + панель управления, дерево и окна-подписчики — ПЛАВАЮЩИМИ окнами поверх.
    Окна двигаются/закрываются, не прячутся за сценой. Подписки с ref-count: адрес жив, пока нужен окну."""

    def __init__(self, link, project: str = "sila26", parent=None, poll_ms: int = 33, site=None):
        super().__init__(parent)
        self.site = site                       # площадка (synchrotronsim.Site): модель на её вычислителе (NetLink)
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
        self._machine_doc = None               # документ машины из пакета опыта (держится для пересборок той же решётки)
        self._desc = None                      # описание сеанса (client.engine.describe_session / с вычислителя)
        self.journal = None                    # журнал правок проекта (synchrotronsim.edits.Journal) — лениво
        self._journal_reason = ""              # почему правка недоступна
        self._edits_win = None                 # окно «Правки проекта» (_ToolWindow)
        self._engine_design_fp = ""            # отпечаток машины модели после горячих правок (мета design)
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
        self.act_cavity = tb.addAction("Резонатор ВЧ", self._open_cavity)  # окно резонатора: яма, фазовый портрет
        self.act_edits = tb.addAction("Правки…", self._edits_dialog)       # журнал правок проекта (этап 4)
        self.act_edits.setToolTip("правки параметров проекта: список, отмена/повтор, пересборка, сохранение ревизией")
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
        if hasattr(win, "request_edit"):                       # окно элемента: правка проекта
            win.request_edit.connect(self.design_edit)
            if getattr(win, "index", -1) >= 0:
                self._fill_design(win)
        if isinstance(win, (ScreenWindow, UndulatorInfoWindow)) or hasattr(win, "bind_link"):
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
        if kind == "graph" and "|" in addr:               # «адрес|множитель|единицы» — таблица окна резонатора
            from .windows import GraphWindow
            a, sc, un = (addr.split("|") + ["1.0", ""])[:3]
            self.open_window(GraphWindow(a, scale=float(sc or 1.0), units=un or None))
        elif kind == "graph":
            self._open_graph(addr)
        elif kind == "screen":
            self._open_screen(addr)                       # кнопка «тепловая карта» → окно-экран
        elif kind == "point":
            ax, _, ay = addr.partition("|")               # кнопка «график x×y» сенсора → окно-точка
            self._open_point(ax, ay)
        elif kind == "diag":
            self._open_diag(addr)
        elif kind == "graph_ref":                         # «⇄» панели: текущее и цель на одном графике
            cur, ref, sc, un = (addr.split("|") + ["", "1.0", ""])[:4]
            self._open_graph_ref(cur, ref, float(sc or 1.0), un)
        elif kind == "cavity":
            self._open_cavity()

    def _open_graph_ref(self, cur: str, ref: str, scale: float = 1.0, units: str = ""):
        from .windows import GraphWindow
        try:
            val = float(ref)
            self.open_window(GraphWindow(cur, ref_value=val, scale=scale, units=units or None))
        except ValueError:
            self.open_window(GraphWindow(cur, ref_addr=ref, scale=scale, units=units or None))

    def _open_cavity(self):
        from .windows import CavityWindow
        if (self._desc or {}).get("cavity") is False:
            self.set_status("у модели нет ВЧ-резонатора")
            return
        self.open_window(CavityWindow())

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
        kind, aids, idx = "", None, -1                     # вид, актуаторы, индекс магнита — из геометрии сцены
        try:
            for n in (self.scene.model.named or []):
                if n.get("name") == element:
                    kind, aids, idx = n.get("kind", ""), n.get("aids"), int(n.get("s_idx", -1))
                    break
        except Exception:
            pass
        self.set_status(f"элемент: {element}")
        self.open_window(InfoWindow(element, kind=kind, aids=aids, index=idx))

    def _open_point(self, x_addr: str, y_addr: str):
        from .windows import PointWindow
        self.open_window(PointWindow(x_addr, y_addr))

    def _open_screen(self, screen_id: str):
        from .windows import ScreenWindow
        self.open_window(ScreenWindow(screen_id))

    # ---- дерево проекта ----
    def build_tree(self, state, loops, screens: List[str]):
        """Дерево по объектам модели (состояние, Scheduler.loops) — обёртка над build_tree_desc."""
        return self.build_tree_desc(sorted(getattr(state, "observables", {}) or {}),
                                    sorted(getattr(state, "undulators", {}) or {}),
                                    [{"name": n, "enabled": bool(getattr(sp, "enabled", False))}
                                     for n, sp in (loops or {}).items()], screens)

    def build_tree_desc(self, observables, undulators, loops, screens):
        """Дерево проекта по описанию сеанса (одинаково для модели здесь и на вычислителе): наблюдаемые, контуры
        [{name, enabled}], ондуляторы, экраны (sid | (sid, вид)), сенсоры, движок, стандартные графики."""
        self.tree.clear()
        loop_names = [l["name"] for l in (loops or [])]

        def node(parent, text, payload=None):
            it = QTreeWidgetItem([text])
            if payload is not None:
                it.setData(0, Qt.UserRole, payload)
            (parent.addChild if isinstance(parent, QTreeWidgetItem) else self.tree.addTopLevelItem)(it)
            return it

        n_obs = node(None, "наблюдаемые")
        for oid in sorted(observables):
            node(n_obs, oid, ("graph", f"obs.{oid}:rms"))
        n_loop = node(None, "контуры")
        for l in (loops or []):
            node(n_loop, f"{l['name']} [{'вкл' if l.get('enabled') else 'выкл'}]", ("graph", f"loop.{l['name']}.enabled"))
        n_und = node(None, "ондуляторы")
        for gid in sorted(undulators):
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
        if "ebpm.x" in observables:
            node(n_sens, "e-BPM 0 (x×y)", ("point", "obs.ebpm.x:0", "obs.ebpm.y:0"))
        n_rf = node(None, "ВЧ-резонатор")
        node(n_rf, "резонатор: параметры, яма, фазовый портрет", ("cavity",))
        n_eng = node(None, "движок")
        for a, t in (("meta.speed", "скорость"), ("meta.model_time", "модельное время"),
                     ("meta.cycle", "цикл")):
            node(n_eng, t, ("graph", a))
        n_std = node(None, "стандартные графики")
        if "fofb_x" in loop_names:                         # η-узор и моды FOFB — нужна матрица отклика FOFB X
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
        elif kind == "cavity":
            self._open_cavity()

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
        err = getattr(self.link, "error", None) if self.remote else None
        if err and err != getattr(self, "_shown_link_err", None):   # обрыв связи с вычислителем — в статус
            self._shown_link_err = err
            self.set_status(f"{err} — пересоберите модель («Опыт…» → «Применить и перестроить»)")
        snap = self.link.poll_snapshot()
        if snap is None:
            return
        seq = int(snap.meta.get("seq", -1))
        if seq == self._last_seq:
            return
        self._last_seq = seq
        self.control.update_meta(snap.meta)
        dfp = (snap.meta.get("design") or {}).get("fingerprint") or ""
        if dfp != self._engine_design_fp:                        # модель применила правку — сверить с журналом
            self._engine_design_fp = dfp
            derr = (snap.meta.get("design") or {}).get("error")
            if derr:
                self.set_status(f"модель не приняла правку: {derr}")
            self._refresh_edits()
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
    @property
    def remote(self) -> bool:
        """Модель на вычислителе площадки (шов Link — сеть), а не в этом процессе."""
        from .netlink import NetLink
        return isinstance(self.link, NetLink)

    def build(self, experiment, machine=None):
        """Собрать модель по опыту в ФОНЕ (первый запуск или пересборка). Готовая сборка — `_attach`.
        machine — документ машины (пакет опыта); None — документ, с которым уже работали (та же решётка), иначе
        из хранилища по ссылке опыта (последняя ревизия: правки таблиц подхватываются пересборкой).
        Удалённый режим (площадка с вычислителем): модель строится на вычислителе, сюда приходит её описание."""
        exp, m = load_experiment_doc(experiment)
        machine = machine or m
        if machine is not None:                                  # документ машины задан явно (пакет опыта)
            self._machine_doc = machine
        elif self.journal is not None and int(self.journal.machine_id) == int(exp.lattice_id) and \
                (self.journal.entries or self.journal.rebuild_needed):
            machine = self.journal.machine()                     # правки проекта — в пересборку
        elif self._machine_doc is not None and int(self._machine_doc.lattice_id) == int(exp.lattice_id):
            machine = self._machine_doc
        else:
            self._machine_doc = None
        if exp.plant.kind not in COCKPIT_PLANTS:
            raise ValueError(f"опыт «{exp.name}»: модель {exp.plant.kind} без кольца — кокпит показывает только "
                             f"{'/'.join(COCKPIT_PLANTS)}; такой опыт исполняется Simulation(exp).run()")
        note = self._opts.get("note") or ""
        where = f" на вычислителе {self.link.url}" if self.remote else ""
        self.set_status(f"строится модель кольца{where}: опыт «{exp.name}», решётка {exp.lattice_id} — это ~15–20 с…"
                        + (f" · {note}" if note else ""))
        if self.scene is None and self._placeholder is not None:
            self._placeholder.setText(f"строится модель кольца{where}…\nопыт «{exp.name}», решётка {exp.lattice_id}")
        if self.remote:
            o = self._opts
            b = _RemoteBuilder(self.link, self.site, exp, machine,
                               {"gap_slew_mm_s": o["gap_slew_mm_s"], "publish_every": o["publish_every"],
                                "gaps": (self._carry or {}).get("gaps")})
        else:
            b = _Builder(exp, machine)
        self._builder = b
        b.done.connect(self._attach)
        b.start()
        return b

    def rebuild(self, experiment, machine=None):
        """Пересобрать модель по новому опыту: движок останавливается (проход записи закрывается), plant
        закрывается, новая сборка — в фоне; окна-подписчики, отметки «● запись», настройки записи, зазоры
        ондуляторов и режим пуск/пауза переносятся. machine — документ машины (пакет опыта)."""
        if self._builder is not None:
            self.set_status("модель ещё строится — дождитесь готовности")
            return None
        exp, m = load_experiment_doc(experiment)
        machine = machine or m
        if exp.plant.kind not in COCKPIT_PLANTS:                 # проверить ДО остановки текущей модели
            QMessageBox.warning(self, "Опыт", f"модель {exp.plant.kind} без кольца — кокпит её не показывает")
            return None
        self._carry = self._detach()
        self._opts["note"] = ""
        return self.build(exp, machine)

    def _detach(self) -> dict:
        """Остановить движок и закрыть plant текущей сборки; вернуть то, что переносится в следующую.
        Удалённый режим: сеанс на вычислителе закрывается; подписки, отметки, экраны, запись и пуск NetLink
        повторит сам, зазоры — из последнего снимка."""
        eng, self._engine = self._engine, None
        carry: dict = dict(self._carry or {})
        if self.remote:
            snap = self.link.poll_snapshot()
            if snap is not None:
                v = snap.values
                slew = float(self._opts.get("gap_slew_mm_s", 500.0))
                gaps = {}
                for gid in (self._desc or {}).get("undulators", []):
                    g, gt = v.get(f"undulator.{gid}.gap"), v.get(f"undulator.{gid}.gap_target")
                    if g is not None and gt is not None:
                        gaps[gid] = (float(g), float(gt), slew)
                if gaps:
                    carry["gaps"] = gaps
                self.link.remember_running(bool(snap.meta.get("running")))
            self.link.close()
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
        """Готовая сборка → сцена, панель, дерево, подписки, движок; перенос из прежней. result — (plant, state,
        sch) модели в этом процессе или ("remote", опыт, описание сеанса) модели на вычислителе."""
        self._builder = None
        if isinstance(result, Exception):
            self.set_status(f"ошибка сборки модели: {type(result).__name__}: {result}")
            if self._placeholder is not None:
                self._placeholder.setText(f"не удалось построить модель:\n{result}")
            return
        first = self.scene is None
        carry, self._carry = (self._carry or {}), None
        o = self._opts
        if isinstance(result, tuple) and len(result) == 3 and result[0] == "remote":
            _, exp, D = result
            plant = state = sch = None
        else:
            plant, state, sch = result
            exp = getattr(sch, "experiment", None)
            init_undulators(state, o["gap_slew_mm_s"])  # один ондулятор сведён (5 мм), остальные разведены (50 мм)
            apply_gaps(state, carry.get("gaps"))         # пересборка: зазоры как были
            D = describe_session(plant, state, sch, exp)
        self.experiment = exp
        self._desc = D
        self._plant, self._state = plant, state
        self._engine_design_fp = ""
        if self.journal is not None:                             # пересборка по машине журнала — правки в модели
            fp = (D.get("machine") or {}).get("fingerprint")
            if fp and self.journal.machine().fingerprint == fp:
                self.journal.mark_rebuilt()
            else:
                n = len(self.journal.entries)
                self.journal = None
                if n:
                    self.set_status(f"машина модели сменилась — правки журнала ({n}) к ней не относятся и сброшены")
        self._refresh_edits()
        from .scene3d import make_ring_scene
        scene = make_ring_scene(gap_scale=o["gap_scale"]) if first else self.scene
        scene.set_geometry(D["geometry"])
        if first:
            self.bind_scene(scene)

        gids = list(D["undulators"])
        self.control.set_undulators(gids)
        self.control.set_loops(D["loops"])               # панель контуров САУ (имя файла логики + вкл/выкл)
        targets = D.get("targets")                       # контуры: текущее → цель (панель + мониторинг)
        if targets is not None:
            self.control.set_targets(targets, self._open_from_request)
            for row in targets:
                self._sub_session(row["cur"])
                if row.get("tgt"):
                    self._sub_session(row["tgt"])
        else:                                            # описание без целей (старый вычислитель) — как раньше
            opt_specs = [tuple(x) for x in D["optics"]]
            self.control.set_optics(opt_specs, self._open_from_request)
            for addr, _, _ in opt_specs:
                self._sub_session(addr)
        for a in ENGINE_META:
            self._sub_session(a)
        for gid in gids:
            self._sub_session(f"undulator.{gid}.gap")
            self._sub_session(f"undulator.{gid}.gap_target")
        for l in D["loops"]:                             # живое состояние контуров → переключатели панели
            self._sub_session(f"loop.{l['name']}.enabled")
        self.build_tree_desc(D["observables"], gids, D["loops"], [tuple(x) if isinstance(x, list) else x
                                                                   for x in D["screens"]])

        if self.remote:                                  # движок — на вычислителе; запись/пуск NetLink повторил
            rec = self.link._record or {}
            run = getattr(exp, "run", None)
            self.control.set_record(bool(rec.get("on")), rec.get("max_turns", getattr(run, "max_turns", None)),
                                    bool(rec.get("repeat", getattr(run, "repeat", False))))
            running = bool(self.link._running)
            self.link.resume()                           # пуск — как был до пересборки, теперь, когда всё готово
        else:
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
            running = eng._running
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
        where = f" (вычислитель {self.link.url})" if self.remote else ""
        if first and o["autorun"]:
            self.link.control(RUN)
            self.set_status(f"{name} готов{where} — идёт прогон")
        elif running:
            self.set_status(f"{name} пересобран{where} — прогон продолжается с нуля")
        else:
            self.set_status(f"{name} {'готов' if first else 'пересобран'}{where} — нажмите «Пуск»")

    # ---- правки проекта (synchrotronsim.edits): журнал у клиента, горячие — в модель, холодные — пересборкой
    def _edit_store(self):
        if self.site is not None and getattr(self.site, "store_url", ""):
            return self.site.store
        from synchrotronsim.store import open_store
        return open_store()

    def _ensure_journal(self):
        """Журнал правок для машины текущей модели: срез-база — рабочая копия таблиц или ревизия хранилища, чья
        машина совпала с машиной модели по отпечатку. None — правка недоступна (причина — _journal_reason)."""
        if self.journal is not None:
            return self.journal
        D = self._desc or {}
        mi = D.get("machine")
        if not mi:
            self._journal_reason = "модель ещё не собрана"
            return None
        fp, lat = mi["fingerprint"], int(mi["lattice_id"])
        from synchrotronsim.edits import Journal
        from synchrotronsim import tableset as T
        cands = []
        try:
            st = self._edit_store()
        except Exception as e:                                   # noqa: BLE001
            st = None
            self._journal_reason = f"хранилище недоступно: {e}"
        try:
            if st is not None and getattr(st, "working", None) and hasattr(st, "working_slices"):
                lib, sls = st.working_slices([lat])
                cands.append((sls[str(lat)], lib, ""))
            else:
                from synchrotronsim.db import PROJECT_ROOT
                lib = T.library_slice(PROJECT_ROOT)
                sl = T.machine_slice(PROJECT_ROOT, lat)
                sl["library"] = T.content_rev(lib)
                cands.append((sl, lib, ""))
        except Exception:                                        # noqa: BLE001
            pass
        try:
            if st is not None and hasattr(st, "machine_slice_at"):
                row = st.revision("machine", str(lat), fp)
                sl = st.machine_slice_at(str(lat), row["rev"])
                cands.append((sl, st.library_content(sl["library"]), row["rev"]))
        except Exception:                                        # noqa: BLE001
            pass
        for sl, lib, rev in cands:
            J = Journal(sl, lib, str(lat), base_rev=rev)
            if J.machine().fingerprint == fp:
                self.journal = J
                self._journal_reason = ""
                return J
        self._journal_reason = ("машины модели нет ни в рабочей копии таблиц, ни в хранилище (модель из пакета "
                                "опыта или таблицы изменились после сборки) — правка проекта недоступна")
        return None

    def _element_name(self, j: int) -> str:
        for n in ((self._desc or {}).get("geometry") or {}).get("named", []):
            if int(n.get("s_idx", -1)) == int(j):
                return n.get("name", "")
        return ""

    def design_edit(self, j: int, column: str, value: str):
        """Правка параметра проекта магнита j: в журнал; горячая — сразу в модель, холодная — ждёт пересборки."""
        J = self._ensure_journal()
        if J is None:
            self.set_status(self._journal_reason)
            return None
        try:
            e = J.make(int(j), column, value, name=self._element_name(j))
            J.push(e)
        except Exception as ex:                                  # noqa: BLE001 — правка не принята
            self.set_status(f"правка не принята: {ex}")
            return None
        if e.cost == "hot":
            self.link.command(e.addr, e.to_dict())
            self.set_status(f"правка в модели: {e.describe()}")
        else:
            self.set_status(f"правка в журнале: {e.describe()} — применится пересборкой («Правки…» → «Пересобрать»)")
        self._refresh_edits()
        return e

    def _undo_edit(self):
        J = self.journal
        r = J.undo() if J is not None else None
        if r is None:
            return
        inv, send = r
        if send:
            self.link.command(inv.addr, inv.to_dict())
        self.set_status(f"отменено: {inv.inverse().describe()}" + (" — модель отстаёт, нужна пересборка"
                                                                    if J.rebuild_needed else ""))
        self._refresh_edits()

    def _redo_edit(self):
        J = self.journal
        e = J.redo() if J is not None else None
        if e is None:
            return
        if e.cost == "hot":
            self.link.command(e.addr, e.to_dict())
        self.set_status(f"повторено: {e.describe()}")
        self._refresh_edits()

    def _rebuild_edits(self):
        if self.journal is None or self.experiment is None:
            return
        self.rebuild(self.current_experiment(), self.journal.machine())

    def _save_edits(self, note: str = ""):
        J = self.journal
        if J is None or not J.entries:
            return None
        try:
            refs = J.save(self._edit_store(), note)
        except Exception as ex:                                  # noqa: BLE001
            self.set_status(f"правки не сохранены: {ex}")
            QMessageBox.warning(self, "Правки проекта", f"правки не сохранены:\n{ex}")
            return None
        self.set_status("правки сохранены: " + (", ".join(r.short for r in refs) or "ревизия не изменилась")
                        + (" · модель отстаёт — нужна пересборка" if J.rebuild_needed else ""))
        self._refresh_edits()
        return refs

    def _discard_edits(self):
        J = self.journal
        if J is None:
            return
        n = 0
        while J.entries:
            inv, send = J.undo()
            if send:
                self.link.command(inv.addr, inv.to_dict())
            n += 1
        J.redo_stack.clear()
        self.set_status(f"правки сброшены: {n}" + (" — модель отстаёт, нужна пересборка" if J.rebuild_needed else ""))
        self._refresh_edits()

    def _fill_design(self, win):
        J = self._ensure_journal()
        if J is None or getattr(win, "index", -1) < 0:
            return
        try:
            from synchrotronsim.edits import element_rows
            keys = element_rows(J.machine(), int(win.index))[:1]
            vals = {c: J.cell(keys, c) for c in getattr(win, "_design", {})}
            win.set_design_values(vals)
        except Exception:                                        # noqa: BLE001
            pass

    def _refresh_edits(self):
        if self._edits_win is not None:
            try:
                self._edits_win.inner.refresh(self.journal, self._journal_reason, self._engine_design_fp)
            except RuntimeError:                                 # окно закрыли
                self._edits_win = None
        for w in list(self._windows):
            if hasattr(w, "set_design_values"):
                self._fill_design(w)

    def _edits_dialog(self):
        from .edits_window import EditsWindow
        if self._edits_win is not None:
            try:
                self._edits_win.show(); self._edits_win.raise_()
                self._ensure_journal(); self._refresh_edits()
                return
            except RuntimeError:
                self._edits_win = None
        ew = EditsWindow()
        ew.undo.connect(self._undo_edit); ew.redo.connect(self._redo_edit)
        ew.rebuild.connect(self._rebuild_edits); ew.save.connect(self._save_edits)
        ew.discard.connect(self._discard_edits)
        tw = _ToolWindow(ew, "Правки проекта", on_close=lambda: setattr(self, "_edits_win", None), parent=self)
        tw.inner = ew
        tw.resize(720, 320)
        self._place_near_cursor(tw)
        tw.show()
        self._edits_win = tw
        self._ensure_journal()
        self._refresh_edits()

    def current_experiment(self):
        """Опыт с ЖИВЫМИ флагами контуров, отмеченными величинами и лимитом записи — то, что сохранится в JSON
        и с чего начнётся пересборка."""
        import dataclasses
        from synchrotronsim.run import Watch
        exp = self.experiment
        eng = self._engine
        if exp is not None and self.remote:
            snap = self.link.poll_snapshot()
            vals = snap.values if snap is not None else {}
            loops = [dataclasses.replace(l, enabled=bool(vals.get(f"loop.{l.name}.enabled", l.enabled)))
                     for l in exp.loops]
            rec = self.link._record or {}
            record = [w for w in exp.record if w.recorder] + [Watch(a, cadence="per_step")
                                                              for a in sorted(self.link._watches)]
            run = dataclasses.replace(exp.run, max_turns=rec.get("max_turns", exp.run.max_turns),
                                      repeat=bool(rec.get("repeat", exp.run.repeat)))
            return dataclasses.replace(exp, loops=loops, record=record, run=run)
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
            self.rebuild(dlg.result_experiment(), dlg.result_machine())

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
            if self.remote:
                self.link.close()                        # сеанс на вычислителе закрывается
        except Exception:
            pass
        event.accept()


class _Builder(QObject):
    """Строит сессию опыта (plant + доска + планировщик — `synchrotronsim.simulation.build_session`) в фоновом
    потоке; `done` эмитится в поток-владелец. Опыт кладётся в планировщик (`sch.experiment`) — паспорт прогона."""
    done = Signal(object)   # (plant, state, sch) | Exception

    def __init__(self, experiment, machine=None):
        super().__init__()
        self.experiment = experiment
        self.machine = machine                  # документ машины (пакет опыта); None — из хранилища по ссылке опыта

    def start(self):
        threading.Thread(target=self._run, name="session-build", daemon=True).start()

    def _run(self):
        try:
            from synchrotronsim.simulation import build_session
            S = build_session(self.experiment, machine=self.machine)
            S.sch.experiment = self.experiment
            self.done.emit((S.plant, S.state, S.sch))
        except Exception as e:                      # noqa: BLE001 — ошибку доставим в UI
            self.done.emit(e)


class _RemoteBuilder(QObject):
    """Живой сеанс на вычислителе площадки: опыт + ссылка на ревизию машины → (модель строится там) → описание
    сеанса. Машина: у площадки есть сервер хранилища — фиксируется там (черновик, если рабочая копия не
    зафиксирована), уходит ссылка; нет — документ машины целиком."""
    done = Signal(object)   # ("remote", опыт, описание) | Exception

    def __init__(self, link, site, experiment, machine, hello):
        super().__init__()
        self.link, self.site, self.experiment, self.machine, self.hello = link, site, experiment, machine, dict(hello)

    def start(self):
        threading.Thread(target=self._run, name="session-build-remote", daemon=True).start()

    def _run(self):
        try:
            exp = self.experiment
            h = dict(self.hello, experiment=exp.to_dict())
            if self.site is not None and getattr(self.site, "store_url", ""):
                st = self.site.store
                m = self.machine if self.machine is not None else st.get("machine", str(exp.lattice_id),
                                                                         exp.machine_rev or None)
                h["machine_ref"] = str(st.pin(m))
            else:
                from synchrotronsim.simulation import resolve_machine
                from synchrotronsim.documents import machine_to_bytes
                m = self.machine if self.machine is not None else resolve_machine(exp)
                h["machine_doc"] = machine_to_bytes(m)
            desc = self.link.open(h, replay_run=False)          # пуск — когда интерфейс готов (_attach)
            self.done.emit(("remote", exp, desc))
        except Exception as e:                      # noqa: BLE001 — ошибку доставим в UI
            self.done.emit(e)


COCKPIT_PLANTS = ("tbt", "composite")            # модели с кольцом (сцена, BPM, пучок); optics — без кольца
START_ON = ("fofb_x", "fofb_y", "tune", "chroma")  # контуры, включённые при старте кокпита (adts ВЫКЛ: октуполи ранг-1)


def load_experiment_doc(spec):
    """Опыт и (если есть) документ машины: Experiment | путь к JSON | пакет .zip (`Experiment.save_bundle` — опыт
    + документ машины, модель считает без БД) | «fb:N» — опыт ОС из таблиц fb_*. → (Experiment, Machine | None)."""
    from synchrotronsim.experiment import Experiment
    if isinstance(spec, Experiment):
        return spec, None
    s = str(spec).strip()
    if s.lower().startswith("fb:"):
        return Experiment.from_fb(int(s[3:])), None
    if s.lower().endswith(".zip"):
        return Experiment.open_bundle(s)
    return Experiment.from_json(s), None


def load_experiment(spec):
    """Опыт для кокпита (см. load_experiment_doc; документ машины пакета отбрасывается)."""
    return load_experiment_doc(spec)[0]


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


# описание сеанса (без Qt — общее с сервером вычислителя): client.engine
from .engine import (loops_info as _loops_info, OPTICS_COMPONENTS as _OPTICS_COMPONENTS,  # noqa: E402
                     optics_specs as _optics_specs, init_undulators, apply_gaps, read_screen_ids as _read_screen_ids,
                     describe_session)


def launch_cockpit(lattice: int = 26, exp=None, parent=None, gap_scale: float = 1.0,
                   steps_per_cycle: int = 9, publish_every: int = 3, settle_turns: int = 250,
                   nsl: int = 8, n_part: int = 12, stride: int = 15, fofb_on: bool = True,
                   optics: bool = True, autorun: bool = False, gap_slew_mm_s: float = 500.0,
                   project: Optional[str] = None, restore_last: bool = True,
                   nsl_dev="auto", experiment=None, site=None) -> Cockpit:
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
    note, m = "", None
    if experiment is not None:
        e, m = load_experiment_doc(experiment)       # пакет .zip — опыт + документ машины
    else:
        e, note = cockpit_experiment(lattice=lattice, settle_turns=settle_turns, steps_per_cycle=steps_per_cycle,
                                     nsl=nsl, n_part=n_part, stride=stride, fofb_on=fofb_on, optics=optics,
                                     nsl_dev=nsl_dev)
    project = project or f"sila{e.lattice_id}"
    if isinstance(site, str):
        from synchrotronsim import Site
        site = Site.load(site)
    if site is not None and getattr(site, "compute_url", ""):  # модель на вычислителе площадки
        from .netlink import NetLink
        link = NetLink(site.compute_url, site.token())
    else:
        link = InProcLink()
    win = Cockpit(link, project=project, parent=parent, site=site)
    win._opts.update(gap_scale=gap_scale, gap_slew_mm_s=gap_slew_mm_s, publish_every=publish_every,
                     autorun=autorun, restore_last=restore_last, note=note)
    win.show()
    try:
        win.build(e, m)
    except Exception as ex:                          # опыт без кольца и т.п. — в окно, не падать
        win.set_status(f"опыт не запущен: {type(ex).__name__}: {ex}")
        if win._placeholder is not None:
            win._placeholder.setText(f"опыт не запущен:\n{ex}")
    return win
