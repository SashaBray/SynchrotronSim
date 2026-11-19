# -*- coding: utf-8 -*-
"""Контракт адресов сигналов — единственная адресная схема клиента.

Адрес — короткая СТАБИЛЬНАЯ между запусками строка, однозначно указывающая на величину в `StateView`.
Стабильность критична: раскладка окон (этап 6) сериализует адреса и переподписывается на НОВЫЙ прогон
по ним — поэтому адрес строится из устойчивых идентификаторов (id элемента/актуатора, имя контура,
group_id ондулятора, имя наблюдаемой), а не из порядкового номера или объекта.

Схема (разбор по префиксу до первой точки → тип):
  obs.<id>[:<red>]          наблюдаемая сенсора (id может содержать точку: `obs.ebpm.x`). Вектор →
                            список; <red> ∈ {rms,mean,max,min,<int>} сворачивает вектор в скаляр.
  drive.<aid>.<field>       привод актуатора (DriveView над Limits). aid содержит точки
                            (`drive.corr.000.hkick.slew`); field ∈ {slew,lag_tau,res,act,odz_lo,odz_hi}.
  undulator.<gid>.<attr>    ондулятор (UndulatorConfig). attr ∈ {gap,gap_target,K,B0,orbit_kick,slew}.
  optics.<attr>[@<dev>]     линейная оптика (OpticsView). Без @ — скаляры кольца (nux,nuy,alphac);
                            с @ — attr ∈ {beta_x,beta_y,alpha_x,alpha_y,eta_x,eta_xp} в месте магнита.
  element.<dev>.<attr>      магнит (MachineView). dev содержит точки (`element.SP03.COR05.x`);
                            attr ∈ {x,y,z,theta0,k1,k2,k3}.
  loop.<name>.<field>       состояние контура (из Scheduler.loops, передаётся явно). field ∈
                            {enabled,rate_hz,period}.
  meta.<name>               мета прогона движка (model_time,wall_time,speed,cycle,seq) — передаётся явно.
  diag.<name>               стандартный диагностический график (client/diagnostics.py). Как `screen.<id>`,
                            считает ДВИЖОК (статично для сессии, кэш) и отдаёт словарь чисел; `resolve` → None.

`resolve` ВСЕГДА возвращает сериализуемое (float|bool|list[float]|None): через шов Link не ходят ни
numpy-массивы, ни живые объекты доски.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np

# --- типы виджетов (предпочтительный способ показать адрес) -------------------------------------
W_GRAPH = "graph"     # скаляр во времени
W_VECTOR = "vector"   # массив (график-профиль или свернуть редукцией)
W_POINT = "point"     # точка (x,y)
W_SCREEN = "screen"   # heatmap + спектр
W_SCALAR = "scalar"   # одно число / флаг
W_INFO = "info"       # составной (инфо-панель)

DRIVE_FIELDS = ("slew", "lag_tau", "res", "act", "odz_lo", "odz_hi")
UND_ATTRS = ("gap", "gap_target", "K", "B0", "orbit_kick", "slew")
OPTICS_AT = {"beta_x": "betax", "beta_y": "betay", "alpha_x": "alphax",
             "alpha_y": "alphay", "eta_x": "etax", "eta_xp": "etaxp"}
OPTICS_SCALARS = ("nux", "nuy", "alphac")
ELEMENT_ATTRS = ("x", "y", "z", "theta0", "k1", "k2", "k3", "yaw", "pitch", "roll")
LOOP_FIELDS = ("enabled", "rate_hz", "period")
META_NAMES = ("model_time", "wall_time", "speed", "cycle", "seq")
DIAG_TITLES = {"orbit_modes": "η-узор и моды FOFB"}   # стандартные графики (diag.<name>) — заголовки для UI

_REDUCERS = {
    "rms": lambda a: float(np.sqrt(np.mean(a ** 2))),
    "mean": lambda a: float(np.mean(a)),
    "max": lambda a: float(np.max(a)),
    "min": lambda a: float(np.min(a)),
}


@dataclass(frozen=True)
class AddrMeta:
    """Метаданные адреса для UI: человеко-читаемое имя, единицы, предпочтительный виджет, тип (kind)."""
    addr: str
    name: str
    units: str
    widget: str
    kind: str


# ------------------------------------------------------------------ сериализация
def _ser(x: Any) -> Any:
    """Привести к сериализуемому: numpy-скаляр→float, numpy-массив/кортеж→list[float], bool/None как есть."""
    if x is None or isinstance(x, (bool, str)):
        return x
    if isinstance(x, (int, float, np.floating, np.integer)):
        return float(x)
    a = np.asarray(x, float)
    return [float(v) for v in a.ravel()]


# ------------------------------------------------------------------ разбор
def kind_of(addr: str) -> str:
    """Тип адреса по префиксу до первой точки ('obs'/'drive'/'undulator'/'optics'/'element'/'actual'/'loop'/'meta')."""
    return addr.split(".", 1)[0]


# ------------------------------------------------------------------ резолв
def resolve(view, addr: str, *, loops: Optional[Dict] = None, meta: Optional[Dict] = None) -> Any:
    """Значение адреса по StateView (сериализуемое) или None, если недоступно/не найдено.
    `loops` — Scheduler.loops (имя→LoopSpec) для `loop.*`; `meta` — мета движка для `meta.*`."""
    kind = kind_of(addr)
    try:
        if kind == "obs":
            return _resolve_obs(view, addr)
        if kind == "drive":
            return _resolve_drive(view, addr)
        if kind == "undulator":
            return _resolve_undulator(view, addr)
        if kind == "optics":
            return _resolve_optics(view, addr)
        if kind == "element":
            return _resolve_element(view, addr)
        if kind == "actual":
            return _resolve_actual(view, addr)
        if kind == "loop":
            return _resolve_loop(loops, addr)
        if kind == "meta":
            return _ser(meta.get(addr.split(".", 1)[1])) if meta else None
    except (KeyError, IndexError, ValueError, AttributeError, TypeError):
        return None
    return None


def _resolve_obs(view, addr):
    body = addr[len("obs."):]
    red = None
    if ":" in body:
        body, red = body.rsplit(":", 1)
    arr = view.read(body)
    if arr is None:
        return None
    a = np.asarray(arr, float).ravel()
    if red is None:
        return _ser(a)
    if red in _REDUCERS:
        return _REDUCERS[red](a)
    return float(a[int(red)])                      # числовой суффикс — индекс элемента


def _resolve_drive(view, addr):
    aid, _, fld = addr[len("drive."):].rpartition(".")
    if fld not in DRIVE_FIELDS:
        return None
    d = view.drives.of(aid)                        # {slew,lag_tau,odz:(lo,hi),res,act}
    if fld == "odz_lo":
        return float(d["odz"][0])
    if fld == "odz_hi":
        return float(d["odz"][1])
    if fld == "act":
        return bool(d["act"])
    return float(d[fld])


def _resolve_undulator(view, addr):
    gid, _, attr = addr[len("undulator."):].partition(".")
    u = view.undulators.get(gid)
    if u is None:
        return None
    if attr == "gap":
        return float(u.gap_mm)
    if attr == "gap_target":
        return float(u.gap_target)
    if attr == "slew":
        return float(u.slew_mm_s)
    if attr == "orbit_kick":
        return float(u.orbit_kick_rad())
    if attr == "K":
        return _ser(u.K())                         # None, если модуль undulator_sr недоступен
    if attr == "B0":
        return _ser(u.B0())
    return None


def _resolve_optics(view, addr):
    opt = view.optics
    if opt is None:
        return None
    body = addr[len("optics."):]
    if "@" in body:
        attr, dev = body.split("@", 1)
        key = OPTICS_AT.get(attr)
        if key is None:
            return None
        row = opt.at_magnet(dev)
        return _ser(row.get(key))
    if body in OPTICS_SCALARS:
        return float(getattr(opt, body))
    return None


def _resolve_element(view, addr):
    mv = view.machine
    if mv is None:
        return None
    dev, _, attr = addr[len("element."):].rpartition(".")
    if attr in ("x", "y", "z"):
        return float(mv.position(dev)[{"x": 0, "y": 1, "z": 2}[attr]])
    if attr in ("theta0", "k1", "k2", "k3"):
        return float(mv.field(dev)[attr])
    if attr in ("yaw", "pitch", "roll"):               # углы ориентации прибора (рад): азимут + наклоны
        return float(mv.orientation(dev)[attr])
    return None


def _resolve_actual(view, addr):
    """Реализованное значение актуатора по id схемы: `actual.corr.NNN.hkick` → `state.actual[<id>]`.
    Это ЖИВОЙ дипольный кик корректора (уставка БП, которую ведёт FOFB). Нет такого id → KeyError → None."""
    sv = getattr(view, "actual", None)
    if sv is None:
        return None
    return float(sv[addr[len("actual."):]])          # StateVector.__getitem__ по id


def _resolve_loop(loops, addr):
    if not loops:
        return None
    name, _, fld = addr[len("loop."):].rpartition(".")
    spec = loops.get(name)
    if spec is None or fld not in LOOP_FIELDS:
        return None
    if fld == "enabled":
        return bool(spec.enabled)
    if fld == "period":
        return float(spec.period)
    return float(spec.rate_hz)


# ------------------------------------------------------------------ перечисление
def enumerate_addresses(view, *, loops: Optional[Dict] = None, meta: Optional[Dict] = None) -> List[str]:
    """Все адреса, разрешимые из текущего состояния и достойные списка (для дерева проекта/подписки).
    Поэлементные `element.*` и `optics.*@dev` НЕ перечисляются поголовно (их тысячи) — они резолвятся по
    запросу через структуру (3D/дерево); здесь — наблюдаемые, приводы, ондуляторы, скаляры оптики,
    контуры, мета."""
    out: List[str] = []
    s = getattr(view, "_s", view)                      # StateView(._s) или сырой ExperimentState
    # наблюдаемые
    for oid, val in getattr(s, "observables", {}).items():
        out.append(f"obs.{oid}")
        if val is not None and np.asarray(val).size > 1:
            out.append(f"obs.{oid}:rms")               # удобный скалярный по умолчанию
    # приводы (по актуаторам схемы)
    schema = getattr(s, "schema", None)
    if schema is not None:
        for aid in schema.ids:
            for fld in DRIVE_FIELDS:
                out.append(f"drive.{aid}.{fld}")
    # ондуляторы
    for gid in getattr(s, "undulators", {}):
        for attr in UND_ATTRS:
            out.append(f"undulator.{gid}.{attr}")
    # оптика — только скаляры кольца
    if getattr(s, "optics", None) is not None:
        out += [f"optics.{k}" for k in OPTICS_SCALARS]
    # контуры
    for name in (loops or {}):
        for fld in LOOP_FIELDS:
            out.append(f"loop.{name}.{fld}")
    # мета
    for k in (meta or {}):
        out.append(f"meta.{k}")
    return out


# ------------------------------------------------------------------ метаданные для UI
_OBS_NAMES = {
    "ebpm.x": ("орбита X (e-BPM)", "м"), "ebpm.y": ("орбита Y (e-BPM)", "м"),
    "photon.x": ("свет X (photon)", "м"), "photon.y": ("свет Y (photon)", "м"),
    "tune": ("рабочая точка", ""), "chroma": ("хроматизм", ""),
}
_DRIVE_UNITS = {"slew": "ед/с", "lag_tau": "с", "res": "ед", "act": "", "odz_lo": "ед", "odz_hi": "ед"}
_UND_UNITS = {"gap": "мм", "gap_target": "мм", "slew": "мм/с", "orbit_kick": "рад", "K": "", "B0": "Тл"}
_OPT_UNITS = {"beta_x": "м", "beta_y": "м", "alpha_x": "", "alpha_y": "",
              "eta_x": "м", "eta_xp": "", "nux": "", "nuy": "", "alphac": ""}
_ELEM_UNITS = {"x": "м", "y": "м", "z": "м", "theta0": "рад", "k1": "1/м²", "k2": "1/м³", "k3": "1/м⁴",
               "yaw": "рад", "pitch": "рад", "roll": "рад"}
_META_UNITS = {"model_time": "с", "wall_time": "с", "speed": "×реал", "cycle": "такт", "seq": ""}


def meta(addr: str) -> AddrMeta:
    """Метаданные адреса (имя/единицы/виджет/kind) — чисто из разбора строки (без состояния)."""
    kind = kind_of(addr)
    if kind == "obs":
        body = addr[len("obs."):]
        red = None
        if ":" in body:
            body, red = body.rsplit(":", 1)
        nm, un = _OBS_NAMES.get(body, (f"obs {body}", "м"))
        if red:
            return AddrMeta(addr, f"{nm} · {red}", un, W_GRAPH, kind)
        return AddrMeta(addr, nm, un, W_VECTOR, kind)
    if kind == "drive":
        aid, _, fld = addr[len("drive."):].rpartition(".")
        return AddrMeta(addr, f"привод {aid} · {fld}", _DRIVE_UNITS.get(fld, ""),
                        W_SCALAR, kind)
    if kind == "undulator":
        gid, _, attr = addr[len("undulator."):].partition(".")
        return AddrMeta(addr, f"ондулятор {gid} · {attr}", _UND_UNITS.get(attr, ""), W_GRAPH, kind)
    if kind == "optics":
        body = addr[len("optics."):]
        attr = body.split("@", 1)[0]
        return AddrMeta(addr, f"оптика {body}", _OPT_UNITS.get(attr, ""), W_GRAPH, kind)
    if kind == "element":
        dev, _, attr = addr[len("element."):].rpartition(".")
        return AddrMeta(addr, f"{dev} · {attr}", _ELEM_UNITS.get(attr, ""), W_GRAPH, kind)
    if kind == "actual":
        aid = addr[len("actual."):]
        suf = aid.rsplit(".", 1)[-1]
        un = {"hkick": "рад", "vkick": "рад", "dk": "1/м²", "dk2": "1/м³", "dk3": "1/м⁴"}.get(suf, "")
        return AddrMeta(addr, f"{aid} · реализ.", un, W_GRAPH, kind)
    if kind == "loop":
        name, _, fld = addr[len("loop."):].rpartition(".")
        return AddrMeta(addr, f"контур {name} · {fld}", "", W_SCALAR, kind)
    if kind == "meta":
        nm = addr.split(".", 1)[1] if "." in addr else addr
        return AddrMeta(addr, f"движок · {nm}", _META_UNITS.get(nm, ""), W_GRAPH, kind)
    if kind == "diag":
        nm = addr.split(".", 1)[1] if "." in addr else addr
        return AddrMeta(addr, DIAG_TITLES.get(nm, f"диагностика · {nm}"), "", W_INFO, kind)
    return AddrMeta(addr, addr, "", W_SCALAR, kind)
