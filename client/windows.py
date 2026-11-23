# -*- coding: utf-8 -*-
"""Окна-подписчики кокпита — один паттерн «виджет на адрес(а)».

База `SubscriberWindow`: объявляет, какие адреса ей нужны (`addresses()`), принимает свежий снимок
(`on_snapshot`), сериализуется для раскладки (`to_dict`/`from_dict`). Кокпит подписывает адреса окна на
Link, на каждом кадре раздаёт снимок всем окнам и переоткрывает окна из сохранённой раскладки. Окна НЕ
держат доску — только адреса и значения из снимка (через шов Link ходят лишь числа).

Типы: `GraphWindow` (скаляр во времени, pyqtgraph), `PointWindow` (точка (x,y) со следом), `InfoWindow`
(атрибуты элемента из снимка; клик по атрибуту просит кокпит открыть график), `ScreenWindow` (heatmap +
спектр; тяжёлый пересчёт делает ДВИЖОК на медленном такте и шлёт как значение адреса `screen.<id>`).
Стандартные графики (`DIAG_WINDOWS`, адрес `diag.<name>`): `OrbitModesWindow` — η-узор, бетатронная волна
и моды FOFB с порогом rcond.
"""
from __future__ import annotations

from collections import deque
from typing import Dict, List, Optional

import numpy as np
import pyqtgraph as pg
from PySide6.QtCore import Qt, Signal, QRectF
from PySide6.QtWidgets import (QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout, QGridLayout, QGroupBox,
                               QHBoxLayout, QLabel, QPushButton, QSpinBox, QVBoxLayout, QWidget)

from . import addresses as A

pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "#222")
pg.setConfigOptions(antialias=True)                     # сглаживание кривых — чище линии/точки

_PALETTE = ["#1a5fb4", "#e15759", "#59a14f", "#b07aa1", "#ff9d00", "#2b9099"]
_LINK_CSS = "color:#185fa5;font-weight:600;text-decoration:none;"


def _scalar(v):
    """Снимковое значение → скаляр (float) или None. Вектор сворачиваем в rms (устойчиво)."""
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    a = np.asarray(v, float).ravel()
    return float(np.sqrt(np.mean(a ** 2))) if a.size else None


def _colormap(name="inferno"):
    """Цветовая карта как в HTML-отчётах (плотли "Inferno"). Fallback — встроенная CET-L9 (≈inferno)."""
    for src in ("matplotlib", None):
        try:
            return pg.colormap.get(name, source=src) if src else pg.colormap.get("CET-L9")
        except Exception:
            continue
    return None


def _export_png(plotwidget, suggested):
    """Сохранить график/карту в PNG через диалог. Возврат — путь или ''."""
    path, _ = QFileDialog.getSaveFileName(plotwidget, "Сохранить график", suggested, "PNG (*.png)")
    if path:
        try:
            from pyqtgraph.exporters import ImageExporter
            item = plotwidget.getPlotItem() if hasattr(plotwidget, "getPlotItem") else plotwidget.scene()
            ImageExporter(item).export(path)
        except Exception:
            return ""
    return path or ""


def _save_button(plotwidget, name):
    b = QPushButton("Сохранить…")
    b.setMaximumWidth(110)
    b.clicked.connect(lambda _c=False: _export_png(plotwidget, name))
    return b


class ValueLink(QLabel):
    """Значение как ГИПЕРССЫЛКУ: клик открывает график этого адреса. Разделитель href — '::', чтобы
    point-адрес 'ax|ay' не ломался. Используется вместо кнопок «график» (значения сами кликабельны)."""

    def __init__(self, addr, *, kind="graph", units="", fmt="{:.6g}", scale=1.0, emit=None):
        super().__init__()
        self.addr, self.kind, self.units = addr, kind, units
        self._fmt, self._scale, self._emit = fmt, float(scale), emit
        self.setTextFormat(Qt.RichText)
        self.setTextInteractionFlags(Qt.TextBrowserInteraction)     # кликаются ссылки
        self.setToolTip(f"{addr} — клик: график")
        self.linkActivated.connect(self._on_link)
        self.set_value(None)

    def set_value(self, v):
        if v is None:
            txt = "—"
        else:
            txt = self._fmt.format(v * self._scale) + (f" {self.units}" if self.units else "")
        self.setText(f'<a href="{self.kind}::{self.addr}" style="{_LINK_CSS}">{txt}</a>')

    def _on_link(self, href):
        kind, _, addr = href.partition("::")
        if self._emit is not None:
            self._emit(kind, addr)


class SubscriberWindow(QWidget):
    """База окна-подписчика. Наследник задаёт KIND, addresses(), on_snapshot(), params()/from_params()."""
    KIND = "base"
    request_open = Signal(str, str)        # (kind, addr) — попросить кокпит открыть новое окно

    def __init__(self, title: str = ""):
        super().__init__()
        self._title = title or self.KIND
        self.setWindowTitle(self._title)

    def title(self) -> str:
        return self._title

    def addresses(self) -> List[str]:
        return []

    def on_snapshot(self, snap) -> None:
        pass

    def params(self) -> dict:
        """Параметры конструкции (без геометрии окна — её хранит раскладка). Переопределяется."""
        return {}

    def to_dict(self) -> dict:
        return {"kind": self.KIND, "title": self._title, "params": self.params()}

    @classmethod
    def from_params(cls, params: dict) -> "SubscriberWindow":
        return cls(**params)


# ------------------------------------------------------------------ график (скаляр во времени)
class GraphWindow(SubscriberWindow):
    KIND = "graph"

    def __init__(self, addr: str, history: int = 600, color: str = _PALETTE[0], title: str = ""):
        m = A.meta(addr)
        super().__init__(title or f"{m.name}")
        self.addr = addr
        self.history = int(history)
        self.color = color
        self._t = deque(maxlen=self.history)
        self._y = deque(maxlen=self.history)

        self.plot = pg.PlotWidget()
        self.plot.setTitle(m.name)
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setLabel("left", m.name, units=m.units or None)
        self.plot.setLabel("bottom", "модельное время", units="с")
        self.plot.setDownsampling(mode="peak")          # лёгкий рендер длинной истории
        self.plot.setClipToView(True)
        self.curve = self.plot.plot(pen=pg.mkPen(self.color, width=2))

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("история"))
        self.spin = QSpinBox()
        self.spin.setRange(20, 20000)
        self.spin.setValue(self.history)
        self.spin.valueChanged.connect(self._set_history)
        ctl.addWidget(self.spin)
        ctl.addStretch(1)
        ctl.addWidget(_save_button(self.plot, f"{self.addr}.png"))

        lay = QVBoxLayout(self)
        lay.addWidget(self.plot, 1)
        lay.addLayout(ctl)

    def _set_history(self, n):
        self.history = int(n)
        self._t = deque(self._t, maxlen=self.history)
        self._y = deque(self._y, maxlen=self.history)

    def addresses(self) -> List[str]:
        return [self.addr]

    def on_snapshot(self, snap) -> None:
        y = _scalar(snap.values.get(self.addr))
        if y is None:
            return
        self._t.append(float(snap.meta.get("model_time", len(self._t))))
        self._y.append(y)
        self.curve.setData(list(self._t), list(self._y))

    def params(self) -> dict:
        return {"addr": self.addr, "history": self.history, "color": self.color}


# ------------------------------------------------------------------ точка (x,y) со следом
class PointWindow(SubscriberWindow):
    KIND = "point"

    def __init__(self, x_addr: str, y_addr: str, trail: int = 200, title: str = ""):
        super().__init__(title or f"точка {x_addr} × {y_addr}")
        self.x_addr = x_addr
        self.y_addr = y_addr
        self.trail = int(trail)
        self._x = deque(maxlen=self.trail)
        self._y = deque(maxlen=self.trail)

        self.plot = pg.PlotWidget()
        self.plot.setTitle(f"{A.meta(x_addr).name} × {A.meta(y_addr).name}")
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setLabel("left", A.meta(y_addr).name, units=A.meta(y_addr).units or None)
        self.plot.setLabel("bottom", A.meta(x_addr).name, units=A.meta(x_addr).units or None)
        self.plot.setAspectLocked(True)                 # орбита x/y изотропна — равный масштаб осей
        self.trace = self.plot.plot(pen=pg.mkPen("#bbb", width=1))
        self.dot = pg.ScatterPlotItem(size=11, pen=pg.mkPen("#b01c22", width=1), brush=pg.mkBrush("#e15759"))
        self.plot.addItem(self.dot)

        ctl = QHBoxLayout()
        ctl.addStretch(1)
        ctl.addWidget(_save_button(self.plot, f"{x_addr}_{y_addr}.png".replace(":", "")))
        lay = QVBoxLayout(self)
        lay.addWidget(self.plot, 1)
        lay.addLayout(ctl)

    def addresses(self) -> List[str]:
        return [self.x_addr, self.y_addr]

    def on_snapshot(self, snap) -> None:
        x = _scalar(snap.values.get(self.x_addr))
        y = _scalar(snap.values.get(self.y_addr))
        if x is None or y is None:
            return
        self._x.append(x)
        self._y.append(y)
        self.trace.setData(list(self._x), list(self._y))
        self.dot.setData([x], [y])

    def params(self) -> dict:
        return {"x_addr": self.x_addr, "y_addr": self.y_addr, "trail": self.trail}


# ------------------------------------------------------------------ инфо-панель элемента (объекто-специфичная)
_KIND_TITLE = {"corrector": "корректор", "kicker": "кикер", "quad": "квадруполь", "sext": "секступоль",
               "oct": "октуполь", "dipole": "диполь", "dq": "диполь-квадруполь", "monitor": "BPM"}
# СОБСТВЕННЫЙ параметр по роли магнита: (суффикс element-адреса, подпись, единицы)
_KIND_FIELDS = {"quad": [("k1", "градиент k1", "1/м²")], "sext": [("k2", "секступоль k2", "1/м³")],
                "oct": [("k3", "октуполь k3", "1/м⁴")], "dipole": [("theta0", "угол θ₀", "рад")],
                "dq": [("theta0", "угол θ₀", "рад"), ("k1", "градиент k1", "1/м²")]}
_OPTICS_FIELDS = [("beta_x", "β_x", "м"), ("beta_y", "β_y", "м"), ("eta_x", "η_x", "м")]
# живая Δk, которую РЕАЛИЗУЕТ оптический контур на семействе магнита (из state.actual[<famid>]) — двигается
# при вкл/выкл ОС. famid приходит в aids из геометрии (scene3d). Номинал k1/k2/k3 — проектный, статичный.
_FAM_DK = {"quad": ("Δk1 семейства (контур)", "1/м²"), "sext": ("Δk2 семейства (контур)", "1/м³"),
           "oct": ("Δk3 семейства (контур)", "1/м⁴")}
# углы ориентации прибора (рад в ядре → показываем в мрад): yaw=азимут, pitch/roll=наклоны (0 у идеала)
_ORIENT_FIELDS = [("yaw", "рысканье (yaw)", "мрад", 1e3), ("pitch", "тангаж (pitch)", "мрад", 1e3),
                  ("roll", "крен (roll)", "мрад", 1e3)]
_DEVICE_KINDS = set(_KIND_FIELDS) | {"corrector", "kicker"}   # приборы с телом → есть ориентация


def info_fields(kind, element, aids):
    """[(addr, label, units[, scale])] — ТОЛЬКО параметры, которыми объект реально обладает (по роли/виду)."""
    rows = [(f"element.{element}.x", "положение X", "м"), (f"element.{element}.y", "положение Y", "м")]
    if kind in ("corrector", "kicker") and aids:            # корректор: ЖИВЫЕ H/V дипольные кики
        rows += [(f"actual.{aids[0]}", "кик H", "рад"), (f"actual.{aids[1]}", "кик V", "рад")]
    else:
        for suf, label, un in _KIND_FIELDS.get(kind, []):   # проектные (статичные) силы магнита
            rows.append((f"element.{element}.{suf}", label, un))
        if kind in _FAM_DK and aids:                        # + ЖИВАЯ Δk семейства, реализованная контуром
            lbl, un = _FAM_DK[kind]
            rows.append((f"actual.{aids[0]}", lbl, un))
    if kind in _DEVICE_KINDS:                               # углы ориентации прибора (yaw/pitch/roll)
        rows += [(f"element.{element}.{suf}", label, un, sc) for suf, label, un, sc in _ORIENT_FIELDS]
    if kind in _KIND_FIELDS:                                 # оптика в точке — отдельной группой (для магнитов)
        rows += [(f"optics.{suf}@{element}", label, un) for suf, label, un in _OPTICS_FIELDS]
    return rows


class InfoWindow(SubscriberWindow):
    """Справка по элементу, СПЕЦИФИЧНАЯ ВИДУ: корректор → кик H/V, квадруполь → k1, секступоль → k2 и т.д.
    Значения — кликабельные гиперссылки (клик → график этого адреса)."""
    KIND = "info"

    def __init__(self, element: str, kind: str = "", aids=None, title: str = ""):
        super().__init__(title or f"{_KIND_TITLE.get(kind, 'элемент')} {element}")
        self.element = element
        self.kind = kind
        self.aids = list(aids or [])
        self._rows = {}            # addr → ValueLink

        grid = QGridLayout()
        head = f"<b>{element}</b>"
        if kind:
            head += f"  <span style='color:#888'>· {_KIND_TITLE.get(kind, kind)}</span>"
        grid.addWidget(QLabel(head), 0, 0, 1, 2)
        r = 1
        for row in info_fields(kind, element, self.aids):
            addr, label, units = row[0], row[1], row[2]
            scale = float(row[3]) if len(row) > 3 else 1.0       # опц. масштаб показа (напр. рад→мрад)
            grid.addWidget(QLabel(label), r, 0)
            vl = ValueLink(addr, units=units, scale=scale, emit=self.request_open.emit)
            self._rows[addr] = vl
            grid.addWidget(vl, r, 1)
            r += 1
        lay = QVBoxLayout(self)
        lay.addLayout(grid)
        lay.addStretch(1)

    def addresses(self) -> List[str]:
        return list(self._rows)

    def on_snapshot(self, snap) -> None:
        for addr, vl in self._rows.items():
            vl.set_value(_scalar(snap.values.get(addr)))

    def params(self) -> dict:
        return {"element": self.element, "kind": self.kind, "aids": self.aids}


# ------------------------------------------------------------------ экран (цветная карта + спектр в точке)
class ScreenWindow(SubscriberWindow):
    """Картина на экране за ондулятором. Тяжёлый расчёт делает ДВИЖОК на медленном такте (ScreenRecorder)
    и шлёт как значение `screen.<id>` = {heatmap, extent(мм), spectrum, gap, e1, error?}. Окно рисует:
    ЦВЕТНУЮ карту (inferno, как HTML-отчёты) с мм-осями и цветовой шкалой, и СПЕКТР излучения (поток vs
    энергия) в ВЫБРАННОЙ ТОЧКЕ — клик по карте выбирает точку (движок пересчитает спектр), команда
    set_screen(addr, area, res, point). Разрешение (mesh) — ручка скорость↔точность."""
    KIND = "screen"

    def __init__(self, screen_id: str, mesh: int = 121, point=None, title: str = ""):
        super().__init__(title or f"экран {screen_id}")
        self.screen_id = screen_id
        self.mesh = int(mesh)
        self._point = tuple(float(p) for p in (point or (0.0, 0.0)))   # точка на экране (мм); (0,0)=на оси
        self._link = None

        # --- тепловая карта: PlotWidget + ImageItem(inferno) + мм-оси + colorbar ---
        self.heat = pg.PlotWidget()
        self.heat.setAspectLocked(True)
        self.heat.setLabel("bottom", "u (горизонт.)", units="мм")
        self.heat.setLabel("left", "v (верт.)", units="мм")
        self.item = pg.ImageItem()
        self._cmap = _colormap("inferno")
        self.cbar = None
        if self._cmap is not None:
            self.item.setColorMap(self._cmap)
        self.heat.addItem(self.item)
        self.marker = pg.ScatterPlotItem(size=22, symbol="+", pen=pg.mkPen("#17c3d6", width=2),
                                         brush=pg.mkBrush("#17c3d6"))      # перекрестье точки спектра (не круг)
        self.marker.setData([self._point[0]], [self._point[1]])
        self.heat.addItem(self.marker)
        if self._cmap is not None:
            self.cbar = pg.ColorBarItem(colorMap=self._cmap, label="поток, ф/с/мрад²/0.1%BW")
            self.cbar.setImageItem(self.item, insert_in=self.heat.getPlotItem())
        self.heat.getPlotItem().scene().sigMouseClicked.connect(self._on_heat_click)

        # --- спектр излучения в точке: поток vs энергия (log y) ---
        self.spec = pg.PlotWidget()
        self.spec.showGrid(x=True, y=True, alpha=0.25)
        self.spec.setLogMode(x=False, y=True)
        self.spec.setLabel("left", "поток", units="ф/с/мрад²/0.1%BW")
        self.spec.setLabel("bottom", "энергия фотона", units="кэВ")
        self.specline = self.spec.plot(pen=pg.mkPen("#e15759", width=2))
        self.info = QLabel("ожидание расчёта экрана…  (клик по карте — спектр в точке)")

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("разрешение (mesh)"))
        self.cmb = QComboBox()
        for v in (61, 121, 201, 301):
            self.cmb.addItem(str(v), v)
        self.cmb.setCurrentText(str(self.mesh))
        self.cmb.currentIndexChanged.connect(self._set_mesh)
        ctl.addWidget(self.cmb)
        ctl.addStretch(1)
        ctl.addWidget(_save_button(self.heat, f"{screen_id}_heat.png"))
        ctl.addWidget(_save_button(self.spec, f"{screen_id}_spec.png"))

        lay = QVBoxLayout(self)
        lay.addWidget(self.info)
        lay.addWidget(self.heat, 3)
        lay.addWidget(self.spec, 2)
        lay.addLayout(ctl)

    def bind_link(self, link):
        self._link = link
        self._push_screen_cfg()

    def _set_mesh(self):
        self.mesh = int(self.cmb.currentData())
        self._push_screen_cfg()

    def _on_heat_click(self, ev):
        vb = self.heat.getPlotItem().vb
        if not vb.sceneBoundingRect().contains(ev.scenePos()):     # клик вне поля карты — игнор
            return
        p = vb.mapSceneToView(ev.scenePos())
        self._point = (float(p.x()), float(p.y()))                 # данные = мм на экране
        self.marker.setData([p.x()], [p.y()])
        self._push_screen_cfg()                                    # спектр пересчитается на медленном такте

    def _push_screen_cfg(self):
        if self._link is not None:
            self._link.set_screen(f"screen.{self.screen_id}", area=None, res=self.mesh, point=self._point)

    def addresses(self) -> List[str]:
        return [f"screen.{self.screen_id}"]

    def on_snapshot(self, snap) -> None:
        d = snap.values.get(f"screen.{self.screen_id}")
        if not isinstance(d, dict):
            return
        if d.get("error"):
            self.info.setText(f"экран: {d['error']}")
            return
        hm = d.get("heatmap")
        if hm is not None:
            arr = np.asarray(hm, float)
            if arr.ndim == 2 and arr.size:
                self.item.setImage(arr.T, autoLevels=False)        # .T: pyqtgraph ждёт (x,y)
                lo, hi = float(np.nanmin(arr)), float(np.nanmax(arr))
                if hi > lo:
                    self.item.setLevels((lo, hi))
                    if self.cbar is not None:
                        try:
                            self.cbar.setLevels((lo, hi))
                        except Exception:
                            pass
                ext = d.get("extent")
                if ext and len(ext) == 4:                          # [x0,x1,y0,y1] в мм — раньше ИГНОРИРОВАЛОСЬ
                    x0, x1, y0, y1 = (float(e) for e in ext)
                    self.item.setRect(QRectF(x0, y0, x1 - x0, y1 - y0))
        sp = d.get("spectrum")
        if isinstance(sp, dict) and sp.get("E_keV"):               # НАСТОЯЩИЙ спектр (flux vs энергия) в точке
            E = np.asarray(sp["E_keV"], float)
            F = np.clip(np.asarray(sp["flux"], float), 1e-30, None)   # log-ось: без нулей/отрицательных
            self.specline.setData(E, F)
        gap = d.get("gap")
        e1 = d.get("e1")
        bits = []
        if gap is not None:
            bits.append(f"зазор {float(gap):.2f} мм")
        if e1 is not None:
            bits.append(f"E₁ {float(e1):.2f} кэВ")
        bits.append(f"точка ({self._point[0]:.2f}, {self._point[1]:.2f}) мм")
        bits.append(f"mesh {self.mesh}")
        if d.get("computing"):                                     # движок пересчитывает картину в фоне
            bits.append("пересчёт…")
        self.info.setText(" · ".join(bits))

    def params(self) -> dict:
        return {"screen_id": self.screen_id, "mesh": self.mesh, "point": list(self._point)}


# ------------------------------------------------------------------ справка по экрану (→ тепловая карта)
class ScreenInfoWindow(SubscriberWindow):
    """Справочное окно фотонного экрана: id/источник + кнопка открыть тепловую карту света (ScreenWindow)."""
    KIND = "screen_info"

    def __init__(self, screen_id: str, title: str = ""):
        super().__init__(title or f"экран {screen_id}")
        self.screen_id = screen_id
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f"<b>Экран {screen_id}</b>"))
        lay.addWidget(QLabel("Фотонный экран за ондулятором — проекция светового потока."))
        btn = QPushButton("Открыть тепловую карту света")
        btn.clicked.connect(lambda: self.request_open.emit("screen", self.screen_id))
        lay.addWidget(btn)
        lay.addStretch(1)

    def params(self) -> dict:
        return {"screen_id": self.screen_id}


# ------------------------------------------------------------------ окно ондулятора (зазор)
class UndulatorInfoWindow(SubscriberWindow):
    """Окно ондулятора: живой зазор (подписка undulator.<gid>.gap) + ручка «задать зазор» (команда движку)."""
    KIND = "undulator"

    def __init__(self, gid: str, title: str = ""):
        super().__init__(title or f"ондулятор {gid}")
        self.gid = gid
        self._addr = f"undulator.{gid}.gap"
        self._link = None
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f"<b>Ондулятор {gid}</b>"))
        row = QHBoxLayout()
        row.addWidget(QLabel("зазор (живой)"))
        self.val = QLabel("—")
        self.val.setStyleSheet("font-weight:600;color:#185fa5;")
        row.addWidget(self.val)
        row.addStretch(1)
        lay.addLayout(row)
        ctl = QHBoxLayout()
        self.spin = QDoubleSpinBox()
        self.spin.setRange(1.0, 50.0)
        self.spin.setSingleStep(0.5)
        self.spin.setValue(5.0)
        self.spin.setSuffix(" мм")
        btn = QPushButton("Задать зазор")
        btn.clicked.connect(self._send)
        ctl.addWidget(self.spin)
        ctl.addWidget(btn)
        ctl.addStretch(1)
        lay.addLayout(ctl)
        lay.addStretch(1)

    def bind_link(self, link):
        self._link = link

    def _send(self):
        if self._link is not None:
            self._link.command(self._addr, float(self.spin.value()))

    def addresses(self) -> List[str]:
        return [self._addr]

    def on_snapshot(self, snap) -> None:
        v = snap.values.get(self._addr)
        if v is not None:
            self.val.setText(f"{float(v):.2f} мм")

    def params(self) -> dict:
        return {"gid": self.gid}


# ------------------------------------------------------------------ окно сенсора (BPM)
class SensorInfoWindow(SubscriberWindow):
    """Окно сенсора (BPM): параметры из структуры (имя, sensor_id, позиция, углы ориентации) + ЖИВОЕ
    показание x/y (клик → график), портрет x×y, И СВОЙСТВА ШУМА датчика (закон + СКО σ). Любой параметр
    берётся на слежение кликом как обычный график; график x×y и шум — дополнения для этого класса приборов."""
    KIND = "sensor"

    def __init__(self, idx, name: str = "", pos=None, sid: str = "", meas_noise=None, title: str = ""):
        super().__init__(title or f"BPM {name or idx}")
        self.idx = int(idx)
        self.sname = name
        self.pos = [float(p) for p in (pos or [0.0, 0.0, 0.0])]
        self.sid = sid
        self.meas_noise = None if meas_noise is None else float(meas_noise)
        self._ax = f"obs.ebpm.x:{self.idx}"
        self._ay = f"obs.ebpm.y:{self.idx}"
        self._rows = {}            # addr → ValueLink (живые поля: показания + углы ориентации)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f"<b>BPM {name}</b>"))
        grid = QGridLayout()
        r = 0
        for k, v in (("sensor_id", sid), ("индекс", str(self.idx)), ("X, м", f"{self.pos[0]:.3f}"),
                     ("Y, м", f"{self.pos[1]:.3f}"), ("Z, м", f"{self.pos[2]:.3f}")):
            grid.addWidget(QLabel(k), r, 0)
            grid.addWidget(QLabel(str(v)), r, 1)
            r += 1
        # живые показания — КЛИКАБЕЛЬНЫЕ (клик → график); для BPM график x×y — лишь ОДНО из дополнений
        for addr, label in ((self._ax, "показание x"), (self._ay, "показание y")):
            grid.addWidget(QLabel(label), r, 0)
            vl = ValueLink(addr, units="мкм", scale=1e6, fmt="{:.2f}", emit=self.request_open.emit)
            self._rows[addr] = vl
            grid.addWidget(vl, r, 1)
            r += 1
        self.vx, self.vy = self._rows[self._ax], self._rows[self._ay]   # именованные ссылки (совместимость)
        grid.addWidget(QLabel("портрет"), r, 0)
        xy = ValueLink(f"{self._ax}|{self._ay}", kind="point", emit=self.request_open.emit)
        xy.setText(f'<a href="point::{self._ax}|{self._ay}" style="{_LINK_CSS}">x × y →</a>')
        grid.addWidget(xy, r, 1)
        r += 1
        if name:                                                 # углы ориентации прибора (yaw/pitch/roll)
            for suf, label in (("yaw", "рысканье (yaw)"), ("pitch", "тангаж (pitch)"), ("roll", "крен (roll)")):
                a = f"element.{name}.{suf}"
                grid.addWidget(QLabel(label), r, 0)
                vl = ValueLink(a, units="мрад", scale=1e3, fmt="{:.3f}", emit=self.request_open.emit)
                self._rows[a] = vl
                grid.addWidget(vl, r, 1)
                r += 1
        lay.addLayout(grid)
        # --- шум датчика: закон + параметр (СКО σ). Шум АДДИТИВНЫЙ ГАУССОВ БЕЛЫЙ, σ=meas_noise (м→мкм) ---
        noise = QGridLayout()
        noise.addWidget(QLabel("закон шума"), 0, 0)
        noise.addWidget(QLabel("гауссов белый (аддитивный, независимый по тактам)"), 0, 1)
        noise.addWidget(QLabel("σ шума"), 1, 0)
        sig = "—" if self.meas_noise is None else f"{self.meas_noise * 1e6:.3f} мкм"
        noise.addWidget(QLabel(sig), 1, 1)
        box = QGroupBox("шум измерения")
        box.setLayout(noise)
        lay.addWidget(box)
        lay.addStretch(1)

    def addresses(self) -> List[str]:
        return list(self._rows)

    def on_snapshot(self, snap) -> None:
        for addr, vl in self._rows.items():
            vl.set_value(_scalar(snap.values.get(addr)))

    def params(self) -> dict:
        return {"idx": self.idx, "name": self.sname, "pos": self.pos, "sid": self.sid,
                "meas_noise": self.meas_noise}


# ------------------------------------------------------------------ стандартный график: η-узор и моды FOFB
_C_ETA, _C_BET, _C_INK2, _C_ZERO = "#2a78d6", "#eb6834", "#52514e", "#9a9994"


class OrbitModesWindow(SubscriberWindow):
    """Стандартный график «η-узор и моды FOFB» (адрес `diag.orbit_modes`, данные считает движок один раз):
    (1) η-узор — орбита в BPM при сдвиге энергии δ (x = η·δ; ручка δ, ppm); (2) бетатронная волна — орбита от
    кика 1 мкрад выбранного корректора (столбец матрицы отклика FOFB X; ручка — номер корректора);
    (3) моды SVD матрицы отклика σ/σ_max (лог) с порогом rcond FOFB: сколько мод он берёт и какую долю
    η-узора может нарисовать корректорами (эту часть он «отнимает» у контура ВЧ-частоты)."""
    KIND = "orbit_modes"
    ADDR = "diag.orbit_modes"

    def __init__(self, delta_ppm: float = 100.0, corr: int = 60, title: str = ""):
        super().__init__(title or A.DIAG_TITLES.get("orbit_modes", "η-узор и моды FOFB"))
        self.delta_ppm = float(delta_ppm)
        self.corr = int(corr)
        self._d = None
        self._sig = None

        def mkplot(title, left, units, bottom):
            p = pg.PlotWidget()
            p.setTitle(title)
            p.showGrid(x=False, y=True, alpha=0.25)
            p.setLabel("left", left, units=units)
            p.setLabel("bottom", bottom)
            return p

        self.p_eta = mkplot("η-узор: орбита при сдвиге энергии пучка (x = η·δ)", "x в BPM", "мкм",
                            "номер BPM по кольцу")
        self.p_eta.addItem(pg.InfiniteLine(pos=0, angle=0, pen=pg.mkPen(_C_ZERO, width=1)))
        self.c_eta = self.p_eta.plot(pen=pg.mkPen(_C_ETA, width=1.5))
        self.p_bet = mkplot("Бетатронная волна: орбита от кика 1 мкрад одного корректора (энергия та же)",
                            "x в BPM", "мкм", "номер BPM по кольцу")
        self.p_bet.addItem(pg.InfiniteLine(pos=0, angle=0, pen=pg.mkPen(_C_ZERO, width=1)))
        self.c_bet = self.p_bet.plot(pen=pg.mkPen(_C_BET, width=1.5))
        self.p_sv = mkplot("Моды корректоров FOFB X (SVD матрицы отклика): σ/σ_max — орбита на единицу кика",
                           "σ / σ_max", None, "номер моды (от сильной к слабой)")
        self.p_sv.setLogMode(x=False, y=True)
        self.c_sv = self.p_sv.plot(pen=pg.mkPen(_C_ETA, width=1.5), symbol="o", symbolSize=5,
                                   symbolBrush=_C_ETA, symbolPen=None)
        self.thr = pg.InfiniteLine(angle=0, pen=pg.mkPen(_C_INK2, width=1, style=Qt.DashLine))
        self.p_sv.addItem(self.thr)

        self.info = QLabel("ожидание данных от движка…")
        self.info.setWordWrap(True)

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("δ"))
        self.sp_delta = QDoubleSpinBox()
        self.sp_delta.setRange(-5000.0, 5000.0)
        self.sp_delta.setDecimals(0)
        self.sp_delta.setSingleStep(10.0)
        self.sp_delta.setSuffix(" ppm")
        self.sp_delta.setValue(self.delta_ppm)
        self.sp_delta.valueChanged.connect(self._set_delta)
        ctl.addWidget(self.sp_delta)
        ctl.addSpacing(16)
        ctl.addWidget(QLabel("корректор"))
        self.sp_corr = QSpinBox()
        self.sp_corr.setRange(0, 9999)
        self.sp_corr.setValue(self.corr)
        self.sp_corr.valueChanged.connect(self._set_corr)
        ctl.addWidget(self.sp_corr)
        self.lb_corr = QLabel("")
        self.lb_corr.setStyleSheet(f"color:{_C_INK2};")
        ctl.addWidget(self.lb_corr)
        ctl.addStretch(1)
        ctl.addWidget(_save_button(self.p_eta, "eta_pattern.png"))
        ctl.addWidget(_save_button(self.p_bet, "betatron_wave.png"))
        ctl.addWidget(_save_button(self.p_sv, "fofb_modes.png"))

        lay = QVBoxLayout(self)
        lay.addWidget(self.p_eta, 1)
        lay.addWidget(self.p_bet, 1)
        lay.addWidget(self.info)
        lay.addWidget(self.p_sv, 1)
        lay.addLayout(ctl)

    # ---- ручки ----
    def _set_delta(self, v):
        self.delta_ppm = float(v)
        self._draw_eta()

    def _set_corr(self, v):
        self.corr = int(v)
        self._draw_bet()

    # ---- отрисовка из кэша данных ----
    def _draw_eta(self):
        d = self._d
        if not d or not d.get("eta_bpm"):
            return
        eta = np.asarray(d["eta_bpm"], float)
        self.c_eta.setData(np.arange(eta.size), eta * self.delta_ppm * 1e-6 * 1e6)     # м·δ → мкм

    def _draw_bet(self):
        d = self._d
        cols = (d or {}).get("cols") or []
        if not cols:
            return
        j = max(0, min(self.corr, len(cols) - 1))
        col = np.asarray(cols[j], float)                     # м/рад; ×1 мкрад → м; ×1e6 → мкм ⇒ численно = R
        self.c_bet.setData(np.arange(col.size), col)
        ids = d.get("corr_ids") or []
        self.lb_corr.setText(f"{ids[j]}  (из {len(cols)})" if j < len(ids) else "")

    def _draw_modes(self):
        d = self._d
        sv = np.asarray(d.get("sv") or [], float)
        if not sv.size:
            return
        self.c_sv.setData(np.arange(1, sv.size + 1), sv)
        rc = float(d.get("rcond", 0.0) or 0.0)
        if rc > 0:
            self.thr.setPos(np.log10(rc))                    # лог-ось: позиция в log10
        self.info.setText(
            f"Порог FOFB X rcond = {rc:.0e}: берёт <b>{d.get('n_kept', 0)}</b> из {d.get('n_modes', sv.size)} мод. "
            f"η-узор в его подпространстве: <b>{100 * float(d.get('eta_reach', 0)):.0f}%</b> нормы "
            f"(без усечения — {100 * float(d.get('eta_reach_full', 0)):.0f}%). Эту часть энергетической орбиты "
            f"FOFB может «рисовать» киками вместо контура ВЧ-частоты.")

    def addresses(self) -> List[str]:
        return [self.ADDR]

    def on_snapshot(self, snap) -> None:
        d = snap.values.get(self.ADDR)
        if not isinstance(d, dict):
            return
        if d.get("error"):
            self.info.setText(f"диагностика: {d['error']}")
            return
        sig = (id(d), d.get("rcond"), d.get("n_kept"), len(d.get("eta_bpm") or []))
        if sig == self._sig:                                   # данные статичны — перерисовка только по смене
            return
        self._sig = sig
        self._d = d
        self.sp_corr.setMaximum(max(0, len(d.get("cols") or []) - 1))
        self._draw_eta()
        self._draw_bet()
        self._draw_modes()

    def params(self) -> dict:
        return {"delta_ppm": self.delta_ppm, "corr": self.corr}


# стандартные графики (дерево «стандартные графики» → окно по имени diag.<name>)
DIAG_WINDOWS = {"orbit_modes": OrbitModesWindow}

WINDOW_REGISTRY = {c.KIND: c for c in (GraphWindow, PointWindow, InfoWindow, ScreenWindow,
                                       ScreenInfoWindow, UndulatorInfoWindow, SensorInfoWindow,
                                       OrbitModesWindow)}


def make_window(d: dict) -> Optional[SubscriberWindow]:
    """Пересоздать окно из сериализованного описания {kind, params} (для восстановления раскладки)."""
    cls = WINDOW_REGISTRY.get(d.get("kind"))
    if cls is None:
        return None
    win = cls.from_params(dict(d.get("params", {})))
    return win
