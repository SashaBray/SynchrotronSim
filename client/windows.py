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
from PySide6.QtGui import QPainterPath
from PySide6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox, QDoubleSpinBox, QFileDialog, QFormLayout,
                               QGridLayout, QGroupBox, QHBoxLayout, QHeaderView, QLabel, QLineEdit, QPushButton,
                               QSpinBox, QSplitter, QTableWidget, QTableWidgetItem, QTreeWidget, QTreeWidgetItem,
                               QVBoxLayout, QWidget)

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
    """Сохранить график/карту через диалог: PNG (растр) или SVG (вектор — кривые и подписи остаются линиями и
    текстом; тепловая карта внутри SVG — картинкой). Формат — по выбранному фильтру/расширению. Возврат — путь или ''."""
    path, flt = QFileDialog.getSaveFileName(plotwidget, "Сохранить график", suggested,
                                            "PNG (*.png);;SVG — вектор (*.svg)")
    if path:
        if "svg" in (flt or "").lower() and not path.lower().endswith(".svg"):
            path = path.rsplit(".", 1)[0] + ".svg"
        try:
            from pyqtgraph.exporters import ImageExporter, SVGExporter
            item = plotwidget.getPlotItem() if hasattr(plotwidget, "getPlotItem") else plotwidget.scene()
            (SVGExporter if path.lower().endswith(".svg") else ImageExporter)(item).export(path)
        except Exception:
            return ""
    return path or ""


def _thin_cross() -> QPainterPath:
    """Символ «тонкое перекрестье»: две линии без заливки (штатный «+» pyqtgraph — толстый залитый крест)."""
    path = QPainterPath()
    path.moveTo(-0.5, 0.0); path.lineTo(0.5, 0.0)
    path.moveTo(0.0, -0.5); path.lineTo(0.0, 0.5)
    return path


_SI_POW = {"y": -24, "z": -21, "a": -18, "f": -15, "p": -12, "n": -9, "µ": -6, "u": -6, "m": -3,
           "k": 3, "M": 6, "G": 9, "T": 12, "P": 15, "E": 18, "Z": 21, "Y": 24}
_RU_PREFIX = {-12: "п", -9: "н", -6: "мк", -3: "м", 3: "к", 6: "М", 9: "Г", 12: "Т"}
_BASE_UNITS = {"м", "с", "рад", "Гц", "Вт", "А", "В", "Т", "эВ", "Дж", "Па", "К", "Гс"}   # можно приставку
_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


class SciAxis(pg.AxisItem):
    """Ось без латинской SI-приставки перед русскими единицами (pyqtgraph писал «m1/м²», «Еф/с/…», «µм» — читается
    как другие единицы). Простая единица (м, с, рад…) — русская приставка: мс, мкм, мрад; составная или уже с
    приставкой — множитель: «×10⁻³ 1/м²». Лог-ось — без множителя (деления подписаны абсолютными степенями)."""

    def updateAutoSIPrefix(self):
        super().updateAutoSIPrefix()
        if self.logMode:
            self.autoSIPrefixScale, self.labelUnitPrefix = 1.0, ""
        else:
            p = _SI_POW.get(self.labelUnitPrefix)
            if p:
                if self.labelUnits in _BASE_UNITS and p in _RU_PREFIX:
                    self.labelUnitPrefix = _RU_PREFIX[p]
                else:
                    self.labelUnitPrefix = "×10" + str(p).translate(_SUP) + (" " if self.labelUnits else "")
        self._updateLabel()

    def labelString(self) -> str:
        s = super().labelString()
        if self.labelUnits == "" and self.autoSIPrefix and self.autoSIPrefixScale not in (0.0, 1.0):
            p = int(round(-np.log10(self.autoSIPrefixScale)))       # безразмерная: «(×10⁻⁶)» вместо «(x1e-06)»
            s = s.replace(f"(x{1.0 / self.autoSIPrefixScale:g})", "(×10" + str(p).translate(_SUP) + ")")
        return s


def _plot() -> "pg.PlotWidget":
    """PlotWidget клиента с осями SciAxis (подписи единиц без латинских приставок)."""
    return pg.PlotWidget(axisItems={"left": SciAxis("left"), "bottom": SciAxis("bottom")})


def _place_tip(tip, vb, x, y) -> None:
    """Подсказка у точки (x, y) в координатах вида, развёрнутая внутрь графика (не обрезается краем)."""
    (x0, x1), (y0, y1) = vb.viewRange()
    tip.setAnchor((1 if x > 0.5 * (x0 + x1) else 0, 0 if y > 0.5 * (y0 + y1) else 1))
    tip.setPos(x, y)


def _fmt_num(v) -> str:
    v = float(v)
    return f"{v:.4g}" if 1e-3 <= abs(v) < 1e5 or v == 0 else f"{v:.3e}"


class CurveHover:
    """Значение под курсором для кривых pyqtgraph: ближайшая точка данных (mode="x" — по оси x, для рядов
    во времени и спектров; "xy" — по расстоянию на экране, для траекторий) подсвечивается кружком, рядом —
    подсказка с подписью и значениями. Лог-оси учитываются (данные берутся исходные, не логарифмы).
    curves — функция → [(PlotDataItem, подпись, fmt_x, fmt_y)], чтобы набор кривых мог меняться."""

    def __init__(self, plot, curves, mode="x"):
        self.pi = plot.getPlotItem(); self.vb = self.pi.vb
        self.curves = curves if callable(curves) else (lambda: curves)
        self.mode = mode
        self.text = pg.TextItem(anchor=(0, 1), color="#111", fill=pg.mkBrush(255, 255, 255, 225),
                                border=pg.mkPen("#888"))
        self.dot = pg.PlotDataItem([], [], pen=None, symbol="o", symbolSize=9, symbolBrush=None,
                                   symbolPen=pg.mkPen("#111", width=1.5))
        for it, z in ((self.dot, 900), (self.text, 1000)):
            it.setZValue(z)
            self.pi.addItem(it, ignoreBounds=True)
        self.text.hide()
        plot.scene().sigMouseMoved.connect(self._move)

    def _log(self):
        try:
            return bool(self.pi.ctrl.logXCheck.isChecked()), bool(self.pi.ctrl.logYCheck.isChecked())
        except Exception:
            return False, False

    def hide(self):
        self.text.hide(); self.dot.setData([], [])

    def _move(self, pos):
        if not self.vb.sceneBoundingRect().contains(pos):
            self.hide()
            return
        p = self.vb.mapSceneToView(pos)
        pw, ph = self.vb.viewPixelSize()
        lx, ly = self._log()
        best = None
        for curve, label, fx, fy in self.curves():
            x, y = curve.getOriginalDataset()
            if x is None or y is None or len(x) == 0:
                continue
            x = np.asarray(x, float); y = np.asarray(y, float)
            with np.errstate(divide="ignore", invalid="ignore"):
                xd = np.log10(x) if lx else x
                yd = np.log10(np.clip(y, 1e-300, None)) if ly else y
            ok = np.isfinite(xd) & np.isfinite(yd)
            if not ok.any():
                continue
            dx = (xd - p.x()) / (pw or 1.0); dy = (yd - p.y()) / (ph or 1.0)
            if self.mode == "x":
                i = int(np.nanargmin(np.where(ok, np.abs(dx), np.inf)))
                d = abs(dx[i]) * 1e3 + abs(dy[i])                     # сначала по x, потом ближайшая по y кривая
            else:
                d2 = np.where(ok, dx * dx + dy * dy, np.inf); i = int(np.argmin(d2)); d = float(d2[i])
            if best is None or d < best[0]:
                best = (d, x[i], y[i], xd[i], yd[i], label, fx, fy)
        if best is None:
            self.hide()
            return
        _, x, y, xd, yd, label, fx, fy = best
        self.dot.setData([x], [y])
        head = f"<b>{label}</b><br>" if label else ""
        self.text.setHtml(f'<div style="font-size:9pt">{head}{fx(x)}<br>{fy(y)}</div>')
        _place_tip(self.text, self.vb, xd, yd)
        self.text.show()


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


def _watch_button(win, addrs):
    """Кнопка «● запись»: включённая — адреса окна пишутся в прогон (Watch) при включённой записи на панели."""
    b = QPushButton("● запись")
    b.setCheckable(True)
    b.setMaximumWidth(110)
    b.setToolTip("писать эту величину в прогон (workspace/runs) — запись включается на панели управления")
    b.setStyleSheet("QPushButton:checked{color:#b01c22;font-weight:600;}")
    b.toggled.connect(lambda on: [win.request_watch.emit(a, bool(on)) for a in addrs])
    return b


class SubscriberWindow(QWidget):
    """База окна-подписчика. Наследник задаёт KIND, addresses(), on_snapshot(), params()/from_params()."""
    KIND = "base"
    request_open = Signal(str, str)        # (kind, addr) — попросить кокпит открыть новое окно
    request_watch = Signal(str, bool)      # (addr, писать?) — попросить кокпит писать адрес в прогон

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

    def __init__(self, addr: str, history: int = 600, color: str = _PALETTE[0], title: str = "",
                 ref_addr: Optional[str] = None, ref_value: Optional[float] = None, scale: float = 1.0,
                 units: Optional[str] = None):
        m = A.meta(addr)
        super().__init__(title or (f"{m.name} и цель" if (ref_addr or ref_value is not None) else f"{m.name}"))
        self.addr = addr
        self.ref_addr = ref_addr                 # цель контура (адрес уставки) — пунктиром на том же графике
        self.ref_value = None if ref_value is None else float(ref_value)   # или постоянная цель (0 — проектная)
        self.scale = float(scale)                # масштаб показа (м → мкм и т.п.)
        self.units = units
        self.history = int(history)
        self.color = color
        self._t = deque(maxlen=self.history)
        self._y = deque(maxlen=self.history)
        self._r = deque(maxlen=self.history)

        self.plot = _plot()
        self.plot.setTitle(m.name)
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        un = self.units if self.units is not None else (m.units if self.scale == 1.0 else "")
        self.plot.setLabel("left", m.name, units=un or None)
        self.plot.setLabel("bottom", "модельное время", units="с")
        self.plot.setDownsampling(mode="peak")          # лёгкий рендер длинной истории
        self.plot.setClipToView(True)
        self.curve = self.plot.plot(pen=pg.mkPen(self.color, width=2), antialias=False)   # без сглаживания: ×2 дешевле
        curves = [(self.curve, m.name, lambda t: f"t = {_fmt_num(t)} с", lambda v, u=un: f"{_fmt_num(v)} {u or ''}")]
        self.ref_curve = None
        if self.ref_addr or self.ref_value is not None:
            self.ref_curve = self.plot.plot(pen=pg.mkPen(_PALETTE[1], width=2, style=Qt.DashLine), antialias=False)
            curves.append((self.ref_curve, "цель", lambda t: f"t = {_fmt_num(t)} с",
                           lambda v, u=un: f"цель {_fmt_num(v)} {u or ''}"))
            leg = self.plot.addLegend(offset=(10, 10))
            leg.addItem(self.curve, "текущее")
            leg.addItem(self.ref_curve, "цель")
        self._hover = CurveHover(self.plot, curves, mode="x")

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("история"))
        self.spin = QSpinBox()
        self.spin.setRange(20, 20000)
        self.spin.setValue(self.history)
        self.spin.valueChanged.connect(self._set_history)
        ctl.addWidget(self.spin)
        ctl.addStretch(1)
        self.btn_watch = _watch_button(self, [self.addr])
        ctl.addWidget(self.btn_watch)
        ctl.addWidget(_save_button(self.plot, f"{self.addr}.png"))

        lay = QVBoxLayout(self)
        lay.addWidget(self.plot, 1)
        lay.addLayout(ctl)

    def _set_history(self, n):
        self.history = int(n)
        self._t = deque(self._t, maxlen=self.history)
        self._y = deque(self._y, maxlen=self.history)

    def addresses(self) -> List[str]:
        return [self.addr] + ([self.ref_addr] if self.ref_addr else [])

    def on_snapshot(self, snap) -> None:
        y = _scalar(snap.values.get(self.addr))
        if y is None:
            return
        self._t.append(float(snap.meta.get("model_time", len(self._t))))
        self._y.append(y * self.scale)
        self.curve.setData(list(self._t), list(self._y))
        if self.ref_curve is not None:
            r = _scalar(snap.values.get(self.ref_addr)) if self.ref_addr else self.ref_value
            self._r.append(float("nan") if r is None else r * self.scale)
            self.ref_curve.setData(list(self._t), list(self._r), connect="finite")

    def params(self) -> dict:
        d = {"addr": self.addr, "history": self.history, "color": self.color}
        if self.ref_addr or self.ref_value is not None:
            d.update(ref_addr=self.ref_addr, ref_value=self.ref_value)
        if self.scale != 1.0 or self.units is not None:
            d.update(scale=self.scale, units=self.units)
        return d


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

        self.plot = _plot()
        self.plot.setTitle(f"{A.meta(x_addr).name} × {A.meta(y_addr).name}")
        self.plot.showGrid(x=True, y=True, alpha=0.25)
        self.plot.setLabel("left", A.meta(y_addr).name, units=A.meta(y_addr).units or None)
        self.plot.setLabel("bottom", A.meta(x_addr).name, units=A.meta(x_addr).units or None)
        self.plot.setAspectLocked(True)                 # орбита x/y изотропна — равный масштаб осей
        self.trace = self.plot.plot(pen=pg.mkPen("#bbb", width=1), antialias=False)
        self.dot = pg.ScatterPlotItem(size=11, pen=pg.mkPen("#b01c22", width=1), brush=pg.mkBrush("#e15759"))
        self.plot.addItem(self.dot)
        mx, my = A.meta(x_addr), A.meta(y_addr)
        self._hover = CurveHover(self.plot, [(self.trace, "", lambda v: f"{mx.name} = {_fmt_num(v)} {mx.units or ''}",
                                              lambda v: f"{my.name} = {_fmt_num(v)} {my.units or ''}")], mode="xy")

        ctl = QHBoxLayout()
        ctl.addStretch(1)
        self.btn_watch = _watch_button(self, [x_addr, y_addr])
        ctl.addWidget(self.btn_watch)
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
        rows += [(f"act.{aids[0]}", "кик H", "рад"), (f"act.{aids[1]}", "кик V", "рад")]
    else:
        for suf, label, un in _KIND_FIELDS.get(kind, []):   # проектные (статичные) силы магнита
            rows.append((f"element.{element}.{suf}", label, un))
        if kind in _FAM_DK and aids:                        # + ЖИВАЯ Δk семейства, реализованная контуром
            lbl, un = _FAM_DK[kind]
            rows.append((f"act.{aids[0]}", lbl, un))
    if kind in _DEVICE_KINDS:                               # углы ориентации прибора (yaw/pitch/roll)
        rows += [(f"element.{element}.{suf}", label, un, sc) for suf, label, un, sc in _ORIENT_FIELDS]
    if kind in _KIND_FIELDS:                                 # оптика в точке — отдельной группой (для магнитов)
        rows += [(f"optics.{suf}@{element}", label, un) for suf, label, un in _OPTICS_FIELDS]
    return rows


class InfoWindow(SubscriberWindow):
    """Справка по элементу, СПЕЦИФИЧНАЯ ВИДУ: корректор → кик H/V, квадруполь → k1, секступоль → k2 и т.д.
    Значения — кликабельные гиперссылки (клик → график этого адреса)."""
    KIND = "info"
    request_edit = Signal(int, str, str)   # (индекс магнита, столбец, значение) — правка ПРОЕКТА (этап 4)

    def __init__(self, element: str, kind: str = "", aids=None, title: str = "", index: int = -1):
        super().__init__(title or f"{_KIND_TITLE.get(kind, 'элемент')} {element}")
        self.element = element
        self.kind = kind
        self.aids = list(aids or [])
        self.index = int(index)    # индекс магнита (порядок пучка); ≥ 0 — можно править параметры проекта
        self._rows = {}            # addr → ValueLink
        self._design = {}          # столбец → (QLabel текущее, QLineEdit новое)

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
        if self.index >= 0:
            lay.addWidget(self._design_box())
        lay.addStretch(1)

    def _design_box(self):
        """Параметры ПРОЕКТА элемента (ячейки таблицы решётки): текущее проектное значение и правка. Сила
        (G, S, O) — модель применяет сразу; остальное — пересборкой («Правки…»)."""
        from synchrotronsim.edits import FIELDS_BY_KIND, GEOMETRY, TITLES, UNITS, EDITABLE
        box = QGroupBox("проект (правка)")
        g = QGridLayout(box)
        cols = list(dict.fromkeys(list(FIELDS_BY_KIND.get(self.kind, ())) + list(GEOMETRY)))
        for n, c in enumerate(cols):
            u = UNITS.get(c, "")
            hot = EDITABLE[c][1] == "hot"
            lab = QLabel(f"{TITLES.get(c, c)}{', ' + u if u else ''}")
            lab.setToolTip(f"столбец {c}; {'модель применит сразу' if hot else 'применится пересборкой модели'}")
            cur = QLabel("—")
            cur.setTextInteractionFlags(Qt.TextSelectableByMouse)
            ed = QLineEdit()
            ed.setPlaceholderText("на ходу" if hot else "пересборкой")
            ed.returnPressed.connect(lambda col=c: self._apply_design(col))
            btn = QPushButton("→")
            btn.setFixedWidth(28)
            btn.setToolTip("применить правку")
            btn.clicked.connect(lambda _=False, col=c: self._apply_design(col))
            g.addWidget(lab, n, 0); g.addWidget(cur, n, 1); g.addWidget(ed, n, 2); g.addWidget(btn, n, 3)
            self._design[c] = (cur, ed)
        return box

    def _apply_design(self, col: str):
        cur, ed = self._design[col]
        t = ed.text().strip().replace(",", ".")
        if t:
            self.request_edit.emit(self.index, col, t)
            ed.clear()

    def set_design_values(self, values: dict) -> None:
        """Проектные значения (текст ячеек журнала правок) — показать в блоке «проект»."""
        for c, (cur, _) in self._design.items():
            v = values.get(c)
            cur.setText("—" if v is None else (v if v != "" else "(пусто)"))

    def addresses(self) -> List[str]:
        return list(self._rows)

    def on_snapshot(self, snap) -> None:
        for addr, vl in self._rows.items():
            vl.set_value(_scalar(snap.values.get(addr)))

    def params(self) -> dict:
        return {"element": self.element, "kind": self.kind, "aids": self.aids, "index": self.index}


# ------------------------------------------------------------------ экран (цветная карта + спектр в точке)
SCREEN_QUANTITIES = (("flux", "поток (E₁ / ε_c)", "поток, ф/с/мрад²/0.1%BW"),
                     ("power", "мощность, Вт/мм²", "плотность мощности, Вт/мм²"))


class ScreenWindow(SubscriberWindow):
    """Картина на экране за ондулятором. Тяжёлый расчёт делает ДВИЖОК на медленном такте и шлёт как значение
    `screen.<id>` = {heatmap, extent(мм), spectrum, gap, e1, quantity, error?}. Окно рисует ЦВЕТНУЮ карту
    (inferno, как HTML-отчёты) с мм-осями и цветовой шкалой и СПЕКТР излучения (поток vs энергия) в ВЫБРАННОЙ
    ТОЧКЕ — клик по карте выбирает точку (движок пересчитает спектр). Что на карте — переключатель «величина»:
      поток на E₁ — интерференционная картина на энергии первой гармоники в узком окне вокруг оси (центральный
                    конус, круглое пятно ~±0,3 мм; ScreenRecorder);
      мощность    — плотность мощности СИ, Вт/мм², по Льенару на ВСЁМ экране (вытянута по горизонтали ±K/γ,
                    + веера соседних диполей; screen_analysis.live_power_map, быстрый режим, живые зазоры).
    «лог. шкала» — цвет по log₁₀ (LOG_DECADES декад от пика): видны слабые веера диполей и хвосты. Наведение
    курсора — значение в точке (карта и спектр). Метки спектра: двойной клик по спектру ставит метку (с «к пику» —
    на ближайший максимум), её вертикаль можно тащить; в таблице имя и энергия редактируются, метки удаляются;
    легенда на графике повторяет таблицу. Команда set_screen(addr, area, res, point, quantity)."""
    KIND = "screen"
    LOG_DECADES = 5

    def __init__(self, screen_id: str, mesh: int = 121, point=None, title: str = "", quantity: str = "flux",
                 log: bool = False, marks=None, snap: bool = True):
        super().__init__(title or f"экран {screen_id}")
        self.screen_id = screen_id
        self.mesh = int(mesh)
        self.quantity = "power" if quantity == "power" else "flux"
        self.log = bool(log)
        self._last = None                                            # последнее значение screen.<id>
        self._hm_ref = self._sp_ref = None                           # уже нарисованные карта/спектр (по ссылке)
        self._arr = self._ext = None                                 # карта [u, v] и её extent (мм) — для подсказки
        self._zoom_pending = self.quantity == "power"                # мощность: первый кадр — зум на пятно ИД
        self._point = tuple(float(p) for p in (point or (0.0, 0.0)))   # точка на экране (мм); (0,0)=на оси
        self._link = None
        self._marks: List[dict] = []
        self._mark_seq = 0
        self._updating = False

        # --- тепловая карта: PlotWidget + ImageItem(inferno) + мм-оси + colorbar ---
        self.heat = pg.PlotWidget()
        self.heat.setAspectLocked(True)
        self.heat.setLabel("bottom", "u (горизонт.)", units="мм")
        self.heat.setLabel("left", "v (верт.)", units="мм")
        for ax in ("bottom", "left"):                               # мм без авто-приставок («ммм» = мкм путало)
            self.heat.getAxis(ax).enableAutoSIPrefix(False)
        self.item = pg.ImageItem()
        self._cmap = _colormap("inferno")
        self.cbar = None
        if self._cmap is not None:
            self.item.setColorMap(self._cmap)
        self.heat.addItem(self.item)
        self.marker = pg.ScatterPlotItem(size=16, symbol=_thin_cross(), pen=pg.mkPen("#17c3d6", width=1.2),
                                         brush=None)               # тонкое перекрестье точки спектра
        self.marker.setData([self._point[0]], [self._point[1]])
        self.heat.addItem(self.marker)
        if self._cmap is not None:
            self.cbar = pg.ColorBarItem(colorMap=self._cmap, label=self._cbar_label())
            self.cbar.setImageItem(self.item, insert_in=self.heat.getPlotItem())
        self.heat_tip = pg.TextItem(anchor=(0, 1), color="#111", fill=pg.mkBrush(255, 255, 255, 225),
                                    border=pg.mkPen("#888"))
        self.heat_tip.setZValue(1000)
        self.heat.getPlotItem().addItem(self.heat_tip, ignoreBounds=True)
        self.heat_tip.hide()
        self.heat.getPlotItem().scene().sigMouseClicked.connect(self._on_heat_click)
        self.heat.scene().sigMouseMoved.connect(self._heat_move)

        # --- спектр излучения в точке: поток vs энергия (log y) + метки с легендой ---
        self.spec = _plot()
        self.spec.showGrid(x=True, y=True, alpha=0.25)
        self.spec.setLogMode(x=False, y=True)
        self.spec.setLabel("left", "поток", units="ф/с/мрад²/0.1%BW")
        self.spec.setLabel("bottom", "энергия фотона", units="кэВ")
        self.specline = self.spec.plot(pen=pg.mkPen("#e15759", width=2))
        self.legend = pg.LegendItem(offset=(-10, 10), labelTextSize="8pt", brush=pg.mkBrush(255, 255, 255, 200),
                                    pen=pg.mkPen("#bbb"))
        self.legend.setParentItem(self.spec.getPlotItem().vb)
        self.legend.hide()
        self._spec_hover = CurveHover(self.spec, [(self.specline, "спектр", lambda e: f"E = {e:.3f} кэВ",
                                                   lambda f: f"{_fmt_num(f)} ф/с/мрад²/0.1%BW")])
        self.spec.scene().sigMouseClicked.connect(self._on_spec_click)
        self.info = QLabel("ожидание расчёта экрана…  (клик по карте — спектр в точке; двойной клик по спектру — метка)")

        # --- таблица меток спектра ---
        self.marks_tbl = QTableWidget(0, 3)
        self.marks_tbl.setHorizontalHeaderLabels(["метка", "E, кэВ", "поток, ф/с/мрад²/0.1%BW"])
        self.marks_tbl.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.marks_tbl.verticalHeader().setVisible(False)
        self.marks_tbl.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.marks_tbl.setMaximumHeight(120)
        self.marks_tbl.itemChanged.connect(self._mark_edited)
        self.marks_tbl.hide()
        mrow = QHBoxLayout()
        mrow.addWidget(QLabel("метки спектра: двойной клик по спектру — добавить"))
        self.chk_snap = QCheckBox("к пику")
        self.chk_snap.setToolTip("новая метка встаёт на ближайший максимум спектра (±1% по энергии)")
        self.chk_snap.setChecked(bool(snap))
        mrow.addWidget(self.chk_snap)
        mrow.addStretch(1)
        self.btn_del = QPushButton("удалить метку")
        self.btn_del.clicked.connect(self._delete_selected)
        self.btn_clr = QPushButton("очистить")
        self.btn_clr.clicked.connect(self._clear_marks)
        mrow.addWidget(self.btn_del)
        mrow.addWidget(self.btn_clr)
        try:
            from PySide6.QtGui import QKeySequence, QShortcut
            QShortcut(QKeySequence.Delete, self.marks_tbl, activated=self._delete_selected)
        except Exception:
            pass

        ctl = QHBoxLayout()
        ctl.addWidget(QLabel("величина"))
        self.qcmb = QComboBox()
        for key, name, _ in SCREEN_QUANTITIES:
            self.qcmb.addItem(name, key)
        self.qcmb.setCurrentIndex(0 if self.quantity == "flux" else 1)
        self.qcmb.currentIndexChanged.connect(self._set_quantity)
        ctl.addWidget(self.qcmb)
        self.chk_log = QCheckBox("лог. шкала")
        self.chk_log.setToolTip(f"цвет по log₁₀ — {self.LOG_DECADES} декад от пика: видны слабые веера диполей")
        self.chk_log.setChecked(self.log)
        self.chk_log.toggled.connect(self._set_log)
        ctl.addWidget(self.chk_log)
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
        lay.addLayout(mrow)
        lay.addWidget(self.marks_tbl)
        lay.addLayout(ctl)

        for m in (marks or []):                                      # метки из раскладки
            try:
                self._add_mark(float(m["E"]), m.get("name"), snap=False)
            except (KeyError, TypeError, ValueError):
                continue

    # ---- ручки ----
    def bind_link(self, link):
        self._link = link
        self._push_screen_cfg()

    def _set_mesh(self):
        self.mesh = int(self.cmb.currentData())
        self._push_screen_cfg()

    def _qlabel(self) -> str:
        return next(lab for key, _, lab in SCREEN_QUANTITIES if key == self.quantity)

    def _cbar_label(self) -> str:
        return ("log₁₀ " if self.log else "") + self._qlabel()

    def _set_cbar_label(self):
        if self.cbar is not None:
            try:
                self.cbar.getAxis("left").setLabel(self._cbar_label())
            except Exception:
                pass

    def _set_quantity(self):
        self.quantity = str(self.qcmb.currentData())
        self._set_cbar_label()
        self.info.setText("пересчёт: " + self.qcmb.currentText() + "…")
        self._zoom_pending = True                                    # новый кадр — свой масштаб вида
        self._hm_ref = None
        self._push_screen_cfg()

    def _set_log(self, on):
        self.log = bool(on)
        self._set_cbar_label()
        self._draw_heat()

    def _zoom_view(self, arr, ext) -> None:
        """Вид на новый кадр. Поток — вся картина (она и так в узком окне у оси). Мощность — пятно ИД: область
        выше 10% пика ×3 (не меньше ±1 мм), иначе на экране 50 мм пятно в пару мм теряется; колесо/«View All»
        показывают весь экран с веерами диполей."""
        x0, x1, y0, y1 = ext
        if self.quantity != "power" or arr.max() <= 0:
            self.heat.setRange(xRange=(x0, x1), yRange=(y0, y1), padding=0)
            return
        nu, nv = arr.shape
        m = arr >= 0.1 * float(arr.max())
        iu = np.where(m.any(axis=1))[0]; iv = np.where(m.any(axis=0))[0]
        du = (x1 - x0) / max(nu - 1, 1); dv = (y1 - y0) / max(nv - 1, 1)
        ulo, uhi = x0 + iu[0] * du, x0 + iu[-1] * du
        vlo, vhi = y0 + iv[0] * dv, y0 + iv[-1] * dv
        hu = max(1.5 * (uhi - ulo), 1.0); hv = max(1.5 * (vhi - vlo), 1.0)
        cu, cv = 0.5 * (ulo + uhi), 0.5 * (vlo + vhi)
        self.heat.setRange(xRange=(cu - hu, cu + hu), yRange=(cv - hv, cv + hv), padding=0)

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
            self._link.set_screen(f"screen.{self.screen_id}", area=None, res=self.mesh, point=self._point,
                                  quantity=self.quantity)

    # ---- тепловая карта: отрисовка (лин./лог.) и значение под курсором ----
    def _draw_heat(self):
        arr = self._arr
        if arr is None or not arr.size:
            return
        hi = float(np.nanmax(arr))
        if self.log and hi > 0:
            lo = hi * 10.0 ** (-self.LOG_DECADES)
            img = np.log10(np.clip(arr, lo, None)); lv = (np.log10(lo), np.log10(hi))
        else:
            img = arr; lv = (float(np.nanmin(arr)), hi)
        self.item.setImage(img, autoLevels=False)                  # карты [u, v] = (x, y), как ждёт col-major pyqtgraph
        if lv[1] > lv[0]:
            self.item.setLevels(lv)
            if self.cbar is not None:
                try:
                    self.cbar.setLevels(lv)
                except Exception:
                    pass
        if self._ext is not None:                                  # пиксель i — с центром в узле сетки i
            x0, x1, y0, y1 = self._ext
            nu, nv = arr.shape
            du = (x1 - x0) / max(nu - 1, 1); dv = (y1 - y0) / max(nv - 1, 1)
            self.item.setRect(QRectF(x0 - du / 2, y0 - dv / 2, (x1 - x0) + du, (y1 - y0) + dv))

    def _heat_value(self, u, v):
        """Значение карты в узле, ближайшем к (u, v) мм; None — вне карты."""
        if self._arr is None or self._ext is None:
            return None
        x0, x1, y0, y1 = self._ext
        nu, nv = self._arr.shape
        iu = int(round((u - x0) / ((x1 - x0) / max(nu - 1, 1)))); iv = int(round((v - y0) / ((y1 - y0) / max(nv - 1, 1))))
        if 0 <= iu < nu and 0 <= iv < nv:
            return float(self._arr[iu, iv])
        return None

    def _heat_move(self, pos):
        vb = self.heat.getPlotItem().vb
        if not vb.sceneBoundingRect().contains(pos):
            self.heat_tip.hide()
            return
        p = vb.mapSceneToView(pos)
        val = self._heat_value(p.x(), p.y())
        if val is None:
            self.heat_tip.hide()
            return
        unit = "Вт/мм²" if self.quantity == "power" else "ф/с/мрад²/0.1%BW"
        self.heat_tip.setHtml(f'<div style="font-size:9pt">u = {p.x():.3f} мм, v = {p.y():.3f} мм<br>'
                              f'<b>{_fmt_num(val)}</b> {unit}</div>')
        _place_tip(self.heat_tip, vb, p.x(), p.y())
        self.heat_tip.show()

    # ---- метки спектра ----
    def _spec_data(self):
        x, y = self.specline.getOriginalDataset()
        if x is None or y is None or len(x) == 0:
            return None, None
        return np.asarray(x, float), np.asarray(y, float)

    def _snap(self, E):
        """Ближайший максимум спектра в окне ±1% по энергии (не уже 2 шагов сетки)."""
        Es, Fs = self._spec_data()
        if Es is None:
            return E
        w = max(0.01 * abs(E), 2.0 * float(np.median(np.diff(Es))) if Es.size > 1 else 0.0)
        m = np.abs(Es - E) <= w
        if not m.any():
            return E
        return float(Es[m][int(np.argmax(Fs[m]))])

    def _flux_at(self, E):
        Es, Fs = self._spec_data()
        if Es is None or not (Es[0] <= E <= Es[-1]):
            return None
        return float(np.exp(np.interp(E, Es, np.log(np.clip(Fs, 1e-300, None)))))

    def _on_spec_click(self, ev):
        if not ev.double() or ev.button() != Qt.LeftButton:
            return
        vb = self.spec.getPlotItem().vb
        if not vb.sceneBoundingRect().contains(ev.scenePos()):
            return
        self._add_mark(float(vb.mapSceneToView(ev.scenePos()).x()), snap=self.chk_snap.isChecked())
        ev.accept()

    def _add_mark(self, E, name=None, snap=True):
        if snap:
            E = self._snap(E)
        self._mark_seq += 1
        color = _PALETTE[(self._mark_seq - 1) % len(_PALETTE)]
        line = pg.InfiniteLine(pos=E, angle=90, movable=True, pen=pg.mkPen(color, width=1, style=Qt.DashLine),
                               hoverPen=pg.mkPen(color, width=2.5))
        dot = pg.PlotDataItem([], [], pen=None, symbol="o", symbolSize=8, symbolBrush=color, symbolPen=None)
        mk = dict(name=str(name or f"метка {self._mark_seq}"), E=float(E), color=color, line=line, dot=dot)
        line.setToolTip("тащите — сдвинуть метку")
        line.sigPositionChanged.connect(lambda _l, mk=mk: self._mark_moved(mk))
        self.spec.getPlotItem().addItem(line, ignoreBounds=True)
        self.spec.getPlotItem().addItem(dot, ignoreBounds=True)
        self._marks.append(mk)
        self._refresh_marks()
        return mk

    def _mark_moved(self, mk):
        if self._updating:
            return
        mk["E"] = float(mk["line"].value())
        self._refresh_marks()

    def _remove_mark(self, mk):
        for it in (mk["line"], mk["dot"]):
            self.spec.getPlotItem().removeItem(it)
        self._marks.remove(mk)

    def _delete_selected(self):
        rows = sorted({i.row() for i in self.marks_tbl.selectedIndexes()} or
                      ({self.marks_tbl.currentRow()} if self.marks_tbl.currentRow() >= 0 else set()), reverse=True)
        for r in rows:
            if 0 <= r < len(self._marks):
                self._remove_mark(self._marks[r])
        self._refresh_marks()

    def _clear_marks(self):
        for mk in list(self._marks):
            self._remove_mark(mk)
        self._refresh_marks()

    def _mark_edited(self, item):
        if self._updating or not (0 <= item.row() < len(self._marks)):
            return
        mk = self._marks[item.row()]
        if item.column() == 0:
            mk["name"] = item.text().strip() or mk["name"]
        elif item.column() == 1:
            try:
                mk["E"] = float(item.text().replace(",", ".").strip())
            except ValueError:
                pass
            self._updating = True
            mk["line"].setValue(mk["E"])
            self._updating = False
        self._refresh_marks()

    def _refresh_marks(self):
        """Точки меток на кривой, легенда на графике и таблица — по текущему спектру."""
        self._updating = True
        try:
            self.legend.clear()
            self.marks_tbl.setRowCount(len(self._marks))
            for r, mk in enumerate(self._marks):
                F = self._flux_at(mk["E"])
                mk["dot"].setData([mk["E"]], [F]) if F is not None else mk["dot"].setData([], [])
                self.legend.addItem(mk["dot"], f"{mk['name']} — {mk['E']:.3f} кэВ")
                cells = (mk["name"], f"{mk['E']:.4f}", "—" if F is None else f"{F:.3e}")
                for c, txt in enumerate(cells):
                    it = self.marks_tbl.item(r, c)
                    if it is None:
                        it = QTableWidgetItem()
                        self.marks_tbl.setItem(r, c, it)
                    it.setText(txt)
                    flags = it.flags() | Qt.ItemIsEnabled | Qt.ItemIsSelectable
                    it.setFlags(flags | Qt.ItemIsEditable if c < 2 else flags & ~Qt.ItemIsEditable)
                self.marks_tbl.item(r, 0).setForeground(pg.mkBrush(mk["color"]))
            self.legend.setVisible(bool(self._marks))
            self.marks_tbl.setVisible(bool(self._marks))
        finally:
            self._updating = False

    # ---- кадр ----
    def addresses(self) -> List[str]:
        return [f"screen.{self.screen_id}"]

    def on_snapshot(self, snap) -> None:
        d = snap.values.get(f"screen.{self.screen_id}")
        if not isinstance(d, dict):
            return
        if d.get("error"):
            self.info.setText(f"экран: {d['error']}")
            return
        if d.get("quantity", "flux") != self.quantity and not d.get("computing"):
            return                                                 # кадр другой величины (до пересчёта) — не рисуем
        self._last = d
        hm = d.get("heatmap") if d.get("quantity", "flux") == self.quantity else None
        if hm is not None and hm is not self._hm_ref:              # новая карта (тот же кадр не перерисовываем)
            arr = np.asarray(hm, float)
            if arr.ndim == 2 and arr.size:
                self._hm_ref = hm
                self._arr = arr
                ext = d.get("extent")
                self._ext = tuple(float(e) for e in ext) if ext and len(ext) == 4 else None
                self._draw_heat()
                if self._ext is not None and self._zoom_pending and not d.get("computing"):
                    self._zoom_pending = False
                    self._zoom_view(arr, self._ext)
        sp = d.get("spectrum")
        if isinstance(sp, dict) and sp.get("E_keV") and sp.get("flux") is not self._sp_ref:
            self._sp_ref = sp.get("flux")                          # НАСТОЯЩИЙ спектр (flux vs энергия) в точке
            E = np.asarray(sp["E_keV"], float)
            F = np.clip(np.asarray(sp["flux"], float), 1e-30, None)   # log-ось: без нулей/отрицательных
            self.specline.setData(E, F)
            if self._marks:
                self._refresh_marks()
        gap = d.get("gap")
        e1 = d.get("e1")
        bits = []
        if d.get("kind") == "bend":                                # экран у диполя: источник и ε_c
            bits.append(f"диполь {d.get('dev', '')}".strip())
            if d.get("ec"):
                bits.append(f"поток на ε_c {float(d['ec']):.2f} кэВ" if d.get("quantity") == "flux"
                            else f"ε_c {float(d['ec']):.2f} кэВ")
        else:
            if gap is not None:
                bits.append(f"зазор {float(gap):.2f} мм")
            if e1 is not None:
                bits.append(f"E₁ {float(e1):.2f} кэВ")
        if d.get("quantity") == "power" and self.quantity == "power" and d.get("peak") is not None:
            bits.append(f"пик {float(d['peak']):.1f} Вт/мм²")
            bits.append(f"перехвачено {float(d.get('total') or 0):.0f} Вт"
                        + (f" (ИД {float(d.get('und_W') or 0):.0f} + веера {int(d.get('n_bend') or 0)} диполей "
                           f"{float(d.get('bend_W') or 0):.0f})" if d.get("n_bend") else ""))
        bits.append(f"точка ({self._point[0]:.2f}, {self._point[1]:.2f}) мм")
        bits.append(f"mesh {self.mesh}")
        if d.get("computing"):                                     # движок пересчитывает картину в фоне
            bits.append("пересчёт…")
        self.info.setText(" · ".join(bits))

    def params(self) -> dict:
        return {"screen_id": self.screen_id, "mesh": self.mesh, "point": list(self._point), "quantity": self.quantity,
                "log": self.log, "snap": self.chk_snap.isChecked(),
                "marks": [{"name": mk["name"], "E": mk["E"]} for mk in self._marks]}


# ------------------------------------------------------------------ справка по экрану (→ тепловая карта)
class ScreenInfoWindow(SubscriberWindow):
    """Справочное окно фотонного экрана: id/источник + кнопка открыть тепловую карту света (ScreenWindow)."""
    KIND = "screen_info"

    def __init__(self, screen_id: str, title: str = ""):
        super().__init__(title or f"экран {screen_id}")
        self.screen_id = screen_id
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel(f"<b>Экран {screen_id}</b>"))
        lay.addWidget(QLabel("Фотонный экран: карта мощности/потока и спектр света источника (ондулятор или диполь)."))
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
            p = _plot()
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
        bpm = lambda v: f"BPM {int(round(v))}"
        self._hovers = [CurveHover(self.p_eta, [(self.c_eta, "η-узор", bpm, lambda v: f"x = {_fmt_num(v)} мкм")]),
                        CurveHover(self.p_bet, [(self.c_bet, "бетатронная волна", bpm,
                                                 lambda v: f"x = {_fmt_num(v)} мкм")]),
                        CurveHover(self.p_sv, [(self.c_sv, "мода SVD", lambda v: f"мода {int(round(v))}",
                                                lambda v: f"σ/σ_max = {_fmt_num(v)}")])]

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


# ------------------------------------------------------------------ ВЧ-резонатор: яма, синусоида, фазовый портрет
class CavityWindow(SubscriberWindow):
    """Резонатор ВЧ (synchrotronsim.longitudinal): НЕИДЕАЛЬНЫЙ резонатор — амплитуда V, фаза φ, частота f
    (номинал h·c/C + сдвиг контура ВЧ + ручной Δf). Ручки меняют резонатор на ходу (команды cavity.*), пучок
    отвечает синхротронными колебаниями.
      синусоида   — V·sin φ (МВ) на периоде вокруг синхронной фазы, потеря U₀ за оборот, яма [φu, φ2] — где
                    частицы устойчивы по фазе; точки — фазы, на которые частицы приходят на резонатор;
      портрет     — (φ, δ): сепаратриса ямы, границы по импульсу δs ± δ_max и по фазе, частицы последнего такта и
                    след последних тактов.
    Яма и сепаратриса — гамильтониан продольного движения (без затухания: ~4·10⁻⁴ за оборот)."""
    KIND = "cavity"

    def __init__(self, trail: int = 40, title: str = ""):
        super().__init__(title or "Резонатор ВЧ")
        self.trail = int(trail)
        self._link = None
        self._b = None
        self._design = None                   # первое состояние — проектный резонатор (кнопка «проектные»)
        self._trail = deque(maxlen=self.trail)
        self._edited = False

        # ---- ручки
        self.sp_v = QDoubleSpinBox(); self.sp_v.setRange(0.0, 200.0); self.sp_v.setDecimals(4)
        self.sp_v.setSingleStep(0.1); self.sp_v.setSuffix(" МВ")
        self.sp_phi = QDoubleSpinBox(); self.sp_phi.setRange(-360.0, 720.0); self.sp_phi.setDecimals(3)
        self.sp_phi.setSingleStep(0.5); self.sp_phi.setSuffix(" °")
        self.sp_df = QDoubleSpinBox(); self.sp_df.setRange(-1e6, 1e6); self.sp_df.setDecimals(2)
        self.sp_df.setSingleStep(10.0); self.sp_df.setSuffix(" Гц")
        self.sp_nb = QSpinBox(); self.sp_nb.setRange(1, 100000)      # верх — h (с первым состоянием)
        self.sp_nb.setToolTip("сгустков в кольце (заполнение): не больше числа сепаратрис h; от него — заряд сгустка "
                              "и средний шаг. На динамику модели не влияет: трекается один сгусток")
        for sp in (self.sp_v, self.sp_phi, self.sp_df, self.sp_nb):
            sp.valueChanged.connect(self._mark_edited)
            sp.textChanged.connect(self._mark_edited)       # набор с клавиатуры — уже правка (до Enter)
            sp.setKeyboardTracking(False)
        self.b_apply = QPushButton("Применить"); self.b_apply.clicked.connect(self._apply)
        self.b_design = QPushButton("Проектные"); self.b_design.clicked.connect(self._to_design)
        for b in (self.b_apply, self.b_design):              # до первого состояния модели в полях нули —
            b.setEnabled(False)                              # «Применить» отправило бы V = 0
        self.b_design.setToolTip("вернуть амплитуду и фазу, с которыми модель собрана, ручной Δf = 0, "
                                 "заполнение — все h мест")
        form = QFormLayout()
        form.addRow("амплитуда V", self.sp_v)
        form.addRow("фаза φ (программная)", self.sp_phi)
        form.addRow("сдвиг частоты Δf (ручной)", self.sp_df)
        form.addRow("сгустков в кольце (заполнение)", self.sp_nb)
        self.chk_zoom = QCheckBox("портрет: сгусток крупно")
        self.chk_zoom.setToolTip("фазовый портрет — вокруг синхронной точки: равновесный сгусток (эллипсы 1σ и 3σ) "
                                 "и частицы модели; снять — вся яма")
        self.chk_zoom.toggled.connect(self._zoom_toggled)
        row = QHBoxLayout(); row.addWidget(self.b_apply); row.addWidget(self.b_design); row.addStretch(1)
        row.addWidget(self.chk_zoom)
        form.addRow(row)
        box_ctl = QGroupBox("резонатор (на ходу)")
        box_ctl.setLayout(form)

        # ---- сведения
        self.info = QLabel("—")
        self.info.setTextFormat(Qt.RichText)
        self.info.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.info.setWordWrap(True)
        box_info = QGroupBox("яма и пучок")
        li = QVBoxLayout(box_info); li.addWidget(self.info)

        # ---- синусоида
        self.p_sin = _plot()
        self.p_sin.setTitle("V·sin φ: куда приходят частицы")
        self.p_sin.showGrid(x=True, y=True, alpha=0.25)
        self.p_sin.setLabel("bottom", "фаза φ", units="°"); self.p_sin.setLabel("left", "прибавка энергии", units="МэВ")
        self.c_sin = self.p_sin.plot(pen=pg.mkPen(_PALETTE[0], width=2))
        self.l_u0 = pg.InfiniteLine(angle=0, pen=pg.mkPen(_PALETTE[1], width=1.5, style=Qt.DashLine),
                                    label="U₀", labelOpts={"position": 0.05, "color": _PALETTE[1]})
        self.r_sin = pg.LinearRegionItem(movable=False, brush=pg.mkBrush(89, 161, 79, 40))
        self.l_phis = pg.InfiniteLine(angle=90, pen=pg.mkPen("#555", width=1, style=Qt.DotLine),
                                      label="φs", labelOpts={"position": 0.92, "color": "#555"})
        self.s_sin = pg.ScatterPlotItem(size=8, pen=pg.mkPen("#111", width=0.8), brush=pg.mkBrush(_PALETTE[4]))
        for it in (self.r_sin, self.l_u0, self.l_phis, self.s_sin):
            self.p_sin.addItem(it)
        self._hs = CurveHover(self.p_sin, [(self.c_sin, "V·sin φ", lambda x: f"φ = {x:.2f}°",
                                            lambda y: f"{_fmt_num(y)} МэВ")], mode="x")

        # ---- фазовый портрет
        self.p_ph = _plot()
        self.p_ph.setTitle("фазовый портрет (φ, δ) и сепаратриса")
        self.p_ph.showGrid(x=True, y=True, alpha=0.25)
        self.p_ph.setLabel("bottom", "фаза φ", units="°"); self.p_ph.setLabel("left", "δ = Δp/p, 10⁻³")
        self.c_up = self.p_ph.plot(pen=pg.mkPen(_PALETTE[1], width=2))
        self.c_lo = self.p_ph.plot(pen=pg.mkPen(_PALETTE[1], width=2))
        pen_b = pg.mkPen(_PALETTE[3], width=1.2, style=Qt.DashLine)
        self.l_dp = pg.InfiniteLine(angle=0, pen=pen_b, label="δs+δmax", labelOpts={"position": 0.08, "color": _PALETTE[3]})
        self.l_dm = pg.InfiniteLine(angle=0, pen=pen_b, label="δs−δmax", labelOpts={"position": 0.08, "color": _PALETTE[3]})
        pen_f = pg.mkPen(_PALETTE[2], width=1.2, style=Qt.DashLine)
        self.l_f1 = pg.InfiniteLine(angle=90, pen=pen_f, label="φu", labelOpts={"position": 0.95, "color": _PALETTE[2]})
        self.l_f2 = pg.InfiniteLine(angle=90, pen=pen_f, label="φ2", labelOpts={"position": 0.95, "color": _PALETTE[2]})
        self.s_tr = pg.ScatterPlotItem(size=3, pen=None, brush=pg.mkBrush(120, 120, 120, 90))
        self.s_now = pg.ScatterPlotItem(size=8, pen=pg.mkPen("#111", width=0.8), brush=pg.mkBrush(_PALETTE[0]))
        self.s_sync = pg.ScatterPlotItem(size=12, symbol="+", pen=pg.mkPen("#111", width=2))
        # равновесный сгусток при нулевом токе (σφ₀, σδ₀): эллипсы 1σ и 3σ вокруг синхронной точки
        self.c_eq1 = pg.PlotDataItem(pen=pg.mkPen(_PALETTE[4], width=1.8))
        self.c_eq3 = pg.PlotDataItem(pen=pg.mkPen(_PALETTE[4], width=1.2, style=Qt.DashLine))
        for it in (self.l_dp, self.l_dm, self.l_f1, self.l_f2, self.s_tr, self.c_eq3, self.c_eq1, self.s_now,
                   self.s_sync):
            self.p_ph.addItem(it)
        self._hp = CurveHover(self.p_ph, [(self.c_up, "сепаратриса", lambda x: f"φ = {x:.2f}°",
                                           lambda y: f"δ = {_fmt_num(y)}·10⁻³"),
                                          (self.c_lo, "сепаратриса", lambda x: f"φ = {x:.2f}°",
                                           lambda y: f"δ = {_fmt_num(y)}·10⁻³")], mode="x")

        # ---- все параметры системы «резонатор — пучок» (synchrotronsim.longitudinal.PARAMS)
        from synchrotronsim.longitudinal import PARAMS
        self.table = QTreeWidget()
        self.table.setColumnCount(3)
        self.table.setHeaderLabels(["параметр", "значение", "откуда / формула"])
        self.table.setRootIsDecorated(True); self.table.setAlternatingRowColors(True)
        self.table.setToolTip("двойной клик по параметру — график во времени")
        self._rows = {}                           # ключ → (элемент дерева, множитель, единицы)
        for group, rows in PARAMS:
            g = QTreeWidgetItem([group]); f = g.font(0); f.setBold(True); g.setFont(0, f)
            self.table.addTopLevelItem(g)
            for key, name, un, sc, _fmt, src, _x in rows:
                it = QTreeWidgetItem([name, "—", src]); it.setData(0, Qt.UserRole, key)
                it.setToolTip(2, src)
                g.addChild(it)
                self._rows[key] = (it, sc, un)
            g.setExpanded(True)
        hd = self.table.header()
        hd.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hd.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hd.setStretchLastSection(True)
        self.table.itemDoubleClicked.connect(self._graph_row)

        left = QWidget(); ll = QVBoxLayout(left); ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(box_ctl); ll.addWidget(box_info); ll.addWidget(self.table, 1)
        right = QWidget(); rl = QVBoxLayout(right); rl.setContentsMargins(0, 0, 0, 0)
        rl.addWidget(self.p_sin, 1); rl.addWidget(self.p_ph, 2)
        split = QSplitter(Qt.Horizontal); split.addWidget(left); split.addWidget(right)
        split.setStretchFactor(0, 1); split.setStretchFactor(1, 1); split.setSizes([700, 640])
        lay = QVBoxLayout(self)
        lay.addWidget(split)

    def _zoom_toggled(self, on):
        if not on:
            self.p_ph.enableAutoRange()
        elif self._b is not None:
            self._zoom_bunch(self._b)

    def _zoom_bunch(self, b):
        """Вид портрета на сгусток: 3σ равновесного (σφ₀, σδ₀) и все частицы модели последнего кадра, с полями."""
        if not b.get("ok"):
            return
        c, ds = float(b["phis_deg"]), float(b["delta_s"]) * 1e3
        sp0, sd0 = float(b.get("sphi0_deg") or 0.0), float(b.get("sdelta0") or 0.0) * 1e3
        xs = [c - 3 * sp0, c + 3 * sp0]; ys = [ds - 3 * sd0, ds + 3 * sd0]
        px, py = self.s_now.getData()
        if px is not None and len(px):
            xs += [float(np.nanmin(px)), float(np.nanmax(px))]; ys += [float(np.nanmin(py)), float(np.nanmax(py))]
        x0, x1, y0, y1 = min(xs), max(xs), min(ys), max(ys)
        if x1 > x0 and y1 > y0:
            self.p_ph.setRange(xRange=(x0, x1), yRange=(y0, y1), padding=0.15)

    def _graph_row(self, item, _col=0):
        """Двойной клик по строке таблицы → график параметра во времени (адрес cavity.<ключ>)."""
        key = item.data(0, Qt.UserRole)
        if not key or key not in self._rows:
            return
        _it, sc, un = self._rows[key]
        self.request_open.emit("graph", f"cavity.{key}|{sc!r}|{un}")

    # ---- связь с движком
    def bind_link(self, link):
        self._link = link

    def _mark_edited(self, *_):
        self._edited = True

    def _send(self, attr, value):
        if self._link is not None:
            self._link.command(f"cavity.{attr}", float(value))

    def _apply(self):
        self._send("volt_mv", self.sp_v.value())
        self._send("phase_deg", self.sp_phi.value())
        self._send("df_hz", self.sp_df.value())
        self._send("n_bunches", self.sp_nb.value())
        self._edited = False
        self._trail.clear()                      # след старого режима не мешает видеть новый

    def _to_design(self):
        if not self._design:
            return
        for sp, v in ((self.sp_v, self._design["V_mv"]), (self.sp_phi, self._design["phase_deg"]), (self.sp_df, 0.0)):
            sp.blockSignals(True); sp.setValue(float(v)); sp.blockSignals(False)
        self.sp_nb.blockSignals(True); self.sp_nb.setValue(int(self._design["h"])); self.sp_nb.blockSignals(False)
        self._apply()

    def addresses(self) -> List[str]:
        return ["cavity.state", "cavity.bunch"]

    def params(self) -> dict:
        return {"trail": self.trail}

    # ---- кадр
    @staticmethod
    def _wrap(p, c):
        return c + np.mod(np.asarray(p, float) - c + 180.0, 360.0) - 180.0

    def on_snapshot(self, snap) -> None:
        b = snap.values.get("cavity.state")
        bunch = snap.values.get("cavity.bunch")
        if not isinstance(b, dict):
            return
        if self._design is None:
            self._design = dict(b)
        self._b = b
        if not self._edited:                     # ручки показывают резонатор модели, пока их не тронули
            # (фокус — не правка: первое поле получает фокус при открытии окна и иначе осталось бы нулём)
            for sp, v in ((self.sp_v, b["V_mv"]), (self.sp_phi, b["phase_deg"]), (self.sp_df, b.get("manual_df_hz", 0.0))):
                sp.blockSignals(True); sp.setValue(float(v)); sp.blockSignals(False)
            self.sp_nb.blockSignals(True)
            self.sp_nb.setMaximum(int(b["h"]))                  # сгустков не больше мест (сепаратрис)
            self.sp_nb.setValue(int(b.get("n_bunches") or b["h"]))
            self.sp_nb.blockSignals(False)
        self.b_apply.setEnabled(True); self.b_design.setEnabled(True)
        ok = bool(b.get("ok"))
        c = float(b["phis_deg"]) if ok else float(b["phase_deg"])
        V, U0 = float(b["V_mv"]), float(b["U0_mev"])
        x = np.linspace(c - 180.0, c + 180.0, 721)
        self.c_sin.setData(x, V * np.sin(np.radians(x)))
        self.l_u0.setValue(U0)
        if ok:
            f1, f2 = sorted((float(b["phiu_deg"]), float(b["phi2_deg"])))
            self.r_sin.setRegion((f1, f2)); self.r_sin.show()
            self.l_phis.setValue(c); self.l_phis.show()
            from synchrotronsim.longitudinal import separatrix
            sx, su = separatrix(b)
            ds = float(b["delta_s"])
            self.c_up.setData(sx, su * 1e3)
            self.c_lo.setData(sx, (2 * ds - su) * 1e3)
            dm = float(b["delta_max"])
            self.l_dp.setValue((ds + dm) * 1e3); self.l_dm.setValue((ds - dm) * 1e3)
            self.l_f1.setValue(f1); self.l_f2.setValue(f2)
            self.s_sync.setData([c], [ds * 1e3])
            for it in (self.l_dp, self.l_dm, self.l_f1, self.l_f2, self.s_sync):
                it.show()
            sp0, sd0 = b.get("sphi0_deg"), b.get("sdelta0")
            if sp0 and sd0:
                t = np.linspace(0.0, 2.0 * np.pi, 121)
                self.c_eq1.setData(c + sp0 * np.cos(t), (ds + sd0 * np.sin(t)) * 1e3)
                self.c_eq3.setData(c + 3 * sp0 * np.cos(t), (ds + 3 * sd0 * np.sin(t)) * 1e3)
            else:
                self.c_eq1.setData([], []); self.c_eq3.setData([], [])
        else:
            self.r_sin.hide(); self.l_phis.hide()
            self.c_up.setData([], []); self.c_lo.setData([], [])
            self.c_eq1.setData([], []); self.c_eq3.setData([], [])
            for it in (self.l_dp, self.l_dm, self.l_f1, self.l_f2, self.s_sync):
                it.hide()
        if isinstance(bunch, dict) and bunch.get("phase_deg") is not None:
            ph = np.asarray(bunch["phase_deg"], float); dl = np.asarray(bunch["delta"], float)
            good = np.isfinite(ph) & np.isfinite(dl)
            ph, dl = np.where(good, ph, np.nan), np.where(good, dl, np.nan)
            pw = self._wrap(ph, c)
            self._trail.append((pw[:, :-1].ravel(), dl[:, :-1].ravel()))
            last_p, last_d = pw[:, -1], dl[:, -1]
            m = np.isfinite(last_p)
            self.s_sin.setData(last_p[m], V * np.sin(np.radians(last_p[m])))
            self.s_now.setData(last_p[m], last_d[m] * 1e3)
            if self._trail:
                tp = np.concatenate([t[0] for t in self._trail]); td = np.concatenate([t[1] for t in self._trail])
                k = np.isfinite(tp) & np.isfinite(td)
                self.s_tr.setData(tp[k], td[k] * 1e3)
        self.info.setText(self._info_html(b))
        if self.chk_zoom.isChecked():
            self._zoom_bunch(b)
        from synchrotronsim.longitudinal import format_param
        for key, (it, _sc, _un) in self._rows.items():
            it.setText(1, format_param(key, b))

    @staticmethod
    def _info_html(b) -> str:
        """Короткая сводка над таблицей: главное — есть ли яма, синхронная фаза, νs, пучок в яме."""
        g = lambda k, f="{:.4g}", s=1.0: (f.format(float(b[k]) * s) if b.get(k) is not None else "—")
        rows = [f"f = <b>{g('freq_hz', '{:.6f}', 1e-6)}</b> МГц · V = <b>{g('V_mv')}</b> МВ · φ<sub>прог</sub> = "
                f"{g('phase_deg', '{:.3f}')}° · U₀ = <b>{g('U0_mev', '{:.2f}', 1e3)}</b> кэВ/об"]
        if b.get("ok"):
            rows.append(f"φs = <b>{g('phis_deg', '{:.3f}')}°</b> · νs = {g('nus', '{:.4e}')} "
                        f"(период {g('Ts_turns', '{:.0f}')} об.) · аксептанс ±{g('delta_max', '{:.3f}', 100)} %")
        else:
            rows.append(f"<span style='color:#a00'><b>ямы нет:</b> {b.get('reason', '')}</span>")
        if b.get("sz0_m") is not None:
            rows.append(f"сгусток (нулевой ток): σz₀ = <b>{g('sz0_m', '{:.3f}', 1e3)}</b> мм = "
                        f"<b>{g('sz0_m', '{:.2f}', 1e12 / 299792458.0)}</b> пс · σφ₀ = {g('sphi0_deg', '{:.3f}')}° · "
                        f"σδ₀ = {g('sdelta0', '{:.4f}', 1e3)}·10⁻³")
        rows.append(f"сгустков в кольце: <b>{g('n_bunches', '{:.0f}')}</b> из {g('n_buckets', '{:.0f}')} мест "
                    f"(сепаратрис, шаг {g('lambda_rf_m', '{:.4f}')} м = {g('lambda_rf_m', '{:.4f}', 1e9 / 299792458.0)} нс)"
                    f" · средний шаг между сгустками {g('bunch_spacing_m', '{:.4f}')} м")
        if b.get("n_part"):
            lost = int(b.get("n_lost", 0) or 0)
            rows.append(f"пучок: в яме <b>{int(b.get('n_inside', 0) or 0)}</b> из {int(b['n_part'])}"
                        + (f" · <span style='color:#a00'>потеряно {lost}</span>" if lost else "")
                        + (f" · ⟨φ⟩ − φs = {g('dphi_deg', '{:+.3f}')}° · амплитуда {g('amp_deg', '{:.3f}')}°"
                           if b.get("dphi_deg") is not None else ""))
        return "<br>".join(rows)


# стандартные графики (дерево «стандартные графики» → окно по имени diag.<name>)
DIAG_WINDOWS = {"orbit_modes": OrbitModesWindow}

WINDOW_REGISTRY = {c.KIND: c for c in (GraphWindow, PointWindow, InfoWindow, ScreenWindow,
                                       ScreenInfoWindow, UndulatorInfoWindow, SensorInfoWindow,
                                       OrbitModesWindow, CavityWindow)}


_ADDR_PARAMS = ("addr", "x_addr", "y_addr")                # параметры окон, в которых лежат адреса


def make_window(d: dict) -> Optional[SubscriberWindow]:
    """Пересоздать окно из сериализованного описания {kind, params} (для восстановления раскладки). Адреса
    приводятся к канонической форме: раскладки со старыми именами (drive.*, actual.*) открываются как прежде."""
    cls = WINDOW_REGISTRY.get(d.get("kind"))
    if cls is None:
        return None
    params = dict(d.get("params", {}))
    for k in _ADDR_PARAMS:
        if isinstance(params.get(k), str):
            params[k] = A.canonical(params[k])
    win = cls.from_params(params)
    return win
