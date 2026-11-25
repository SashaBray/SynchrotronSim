# -*- coding: utf-8 -*-
"""Единое адресное пространство SynchrotronSim (docs/LIBRARY_ARCHITECTURE.md, раздел 3.7; шаг 4).

Адрес — короткая СТАБИЛЬНАЯ между запусками строка, однозначно указывающая на величину: её показывает клиент,
записывает `Watch` в прогон `Run`, двигает команда. Грамматика:

    <область>.<объект>[.<атрибут>][@<место>][:<редукция>]

  obs.<id>[:<red>]          наблюдаемая (выход Observer): `obs.ebpm.x`, `obs.ebpm.x:rms`, `obs.tune:0`
  act.<aid>                 РЕАЛИЗОВАННОЕ значение актуатора (state.actual): `act.corr.000.hkick`, `act.rf.freq`
  act.<aid>.<поле привода>  привод актуатора (Limits): поле ∈ {slew, lag_tau, res, act, odz_lo, odz_hi}
  undulator.<gid>.<attr>    ондулятор: attr ∈ {gap, gap_target, K, B0, orbit_kick, slew}
  optics.<attr>[@<dev>]     линейная оптика: скаляры кольца (nux, nuy, alphac) или функции Твисса у магнита
  element.<dev>.<attr>      магнит: x, y, z, theta0, k1, k2, k3, yaw, pitch, roll
  loop.<name>.<field>       контур: enabled, rate_hz, period, gain
  meta.<name>               мета прогона: model_time, wall_time, speed, cycle, seq, turn
  screen.<id>               экран (тепловая карта + спектр; ДОРОГОЙ — считает движок по запросу)
  diag.<name>               стандартный график (ДОРОГОЙ, по запросу)
  rec.<series>              ряд, выданный рекордером (обработка данных контуров): `rec.orbit_rms.x`

<red> ∈ {rms, mean, max, min, <индекс>} сворачивает вектор в скаляр.

ПСЕВДОНИМЫ прежних имён (раскладки окон, подписки, пользовательский код продолжают работать):
  drive.<aid>.<f>  → act.<aid>.<f>      (привод — бывшая область drive)
  actual.<aid>     → act.<aid>          (реализованное значение — бывшая область actual)
  <aid> без области (corr.NNN.hkick, quad.<fam>.dk, sext.<fam>.dk2, oct.<fam>.dk3, rf.freq) → act.<aid>
`canonical(addr)` приводит любой адрес к канонической форме; `resolve`/`meta` принимают и старые имена.

`resolve` ВСЕГДА возвращает сериализуемое (float | bool | list[float] | None): через шов Link не ходят ни
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

AREAS = {
    "obs": "наблюдаемая", "act": "актуатор и его привод", "undulator": "ондулятор", "optics": "линейная оптика",
    "element": "магнит", "loop": "контур", "meta": "мета прогона", "screen": "экран", "diag": "стандартный график",
    "rec": "ряд рекордера", "design": "параметр проекта", "ref": "уставка контура", "cavity": "ВЧ-резонатор",
}
CAVITY_WRITABLE = ("volt_mv", "phase_deg", "df_hz", "n_bunches")   # ручки окна «Резонатор ВЧ» (+ заполнение)
CAVITY_UNITS = {"volt_mv": "МВ", "phase_deg": "°", "df_hz": "Гц", "freq_hz": "Гц", "state": "", "bunch": "",
                "n_bunches": ""}
DRIVE_FIELDS = ("slew", "lag_tau", "res", "act", "odz_lo", "odz_hi")
UND_ATTRS = ("gap", "gap_target", "K", "B0", "orbit_kick", "slew")
OPTICS_AT = {"beta_x": "betax", "beta_y": "betay", "alpha_x": "alphax",
             "alpha_y": "alphay", "eta_x": "etax", "eta_xp": "etaxp"}
OPTICS_SCALARS = ("nux", "nuy", "alphac")
ELEMENT_ATTRS = ("x", "y", "z", "theta0", "k1", "k2", "k3", "yaw", "pitch", "roll")
LOOP_FIELDS = ("enabled", "rate_hz", "period", "gain")
META_NAMES = ("model_time", "wall_time", "speed", "cycle", "seq", "turn")
DIAG_TITLES = {"orbit_modes": "η-узор и моды FOFB"}   # стандартные графики (diag.<name>) — заголовки для UI
ACTUATOR_TYPES = ("corr", "quad", "sext", "oct")       # первый токен id схемы актуатора (+ rf.freq)

# ручки, которые принимает движок клиента (команда → доска); прочие адреса — только чтение
WRITABLE_UNDULATOR = ("gap", "gap_target", "slew")
WRITABLE_LOOP = ("enabled", "gain")
WRITABLE_DRIVE = ("slew", "lag_tau", "res", "act")

_REDUCERS = {
    "rms": lambda a: float(np.sqrt(np.mean(a ** 2))),
    "mean": lambda a: float(np.mean(a)),
    "max": lambda a: float(np.max(a)),
    "min": lambda a: float(np.min(a)),
}


@dataclass(frozen=True)
class AddrMeta:
    """Метаданные адреса: имя для человека, единицы, предпочтительный виджет, область (kind), можно ли писать,
    форма значения (scalar | vector | map | info), кто хранит (board — доска; engine — считает движок по
    запросу; recorder — выдаёт рекордер; run — мета прогона), цена (cheap | heavy)."""
    addr: str
    name: str
    units: str
    widget: str
    kind: str
    writable: bool = False
    shape: str = "scalar"
    owner: str = "board"
    cost: str = "cheap"
    scope: str = "state"      # state — живое состояние (ток, зазор, усиление); project — параметр проекта (машина)
    edit: str = ""            # для параметров проекта: hot — модель применяет сразу; cold — пересборкой


@dataclass(frozen=True)
class Parsed:
    """Разбор канонического адреса: область, объект, атрибут/поле, место (@), редукция (:)."""
    area: str
    obj: str
    attr: str = ""
    at: str = ""
    red: str = ""


# ------------------------------------------------------------------ сериализация
def _ser(x: Any) -> Any:
    """Привести к сериализуемому: numpy-скаляр→float, numpy-массив/кортеж→list[float], bool/None как есть."""
    if x is None or isinstance(x, (bool, str)):
        return x
    if isinstance(x, (int, float, np.floating, np.integer)):
        return float(x)
    a = np.asarray(x, float)
    return [float(v) for v in a.ravel()]


# ------------------------------------------------------------------ псевдонимы и разбор
def canonical(addr: str) -> str:
    """Каноническая форма адреса (псевдонимы прежних имён → новые). Канонический адрес не меняется."""
    if not addr:
        return addr
    head, _, rest = addr.partition(".")
    if head == "drive" and rest:
        return "act." + rest
    if head == "actual" and rest:
        return "act." + rest
    if head in ACTUATOR_TYPES and rest:
        return "act." + addr
    if addr == "rf.freq" or (head == "rf" and rest.startswith("freq.") and rest[5:] in DRIVE_FIELDS):
        return "act." + addr
    return addr


def kind_of(addr: str) -> str:
    """Область адреса (после приведения псевдонимов): 'obs' / 'act' / 'undulator' / … ."""
    return canonical(addr).split(".", 1)[0]


def split_reducer(body: str):
    """'ebpm.x:rms' → ('ebpm.x', 'rms'); без редукции → (body, '')."""
    if ":" in body:
        b, red = body.rsplit(":", 1)
        return b, red
    return body, ""


def parse(addr: str) -> Parsed:
    """Разбор адреса по грамматике (псевдонимы приводятся). Неизвестная область — ValueError."""
    a = canonical(addr)
    area, _, body = a.partition(".")
    if area not in AREAS:
        raise ValueError(f"адрес {addr!r}: неизвестная область {area!r} (допустимы {sorted(AREAS)})")
    if not body:
        raise ValueError(f"адрес {addr!r}: после области нет объекта")
    body, red = split_reducer(body)
    if area == "obs":
        return Parsed(area, body, red=red)
    if area == "act":
        aid, _, fld = body.rpartition(".")
        if fld in DRIVE_FIELDS and aid:
            return Parsed(area, aid, fld, red=red)
        return Parsed(area, body, red=red)
    if area == "optics":
        attr, _, dev = body.partition("@")
        return Parsed(area, "", attr, at=dev, red=red)
    if area in ("element", "loop", "design"):
        obj, _, attr = body.rpartition(".")
        return Parsed(area, obj, attr, red=red)
    if area == "undulator":
        gid, _, attr = body.partition(".")
        return Parsed(area, gid, attr, red=red)
    return Parsed(area, body, red=red)                     # meta / screen / diag / rec


def is_address(s: str) -> bool:
    """Строка — адрес известной области (в т.ч. по псевдониму)?"""
    try:
        parse(s)
        return True
    except ValueError:
        return False


def _reduce(val, red: str):
    if not red or val is None:
        return _ser(val)
    a = np.asarray(val, float).ravel()
    if red in _REDUCERS:
        return _REDUCERS[red](a)
    return float(a[int(red)])                              # числовой суффикс — индекс элемента


# ------------------------------------------------------------------ резолв
def resolve(view, addr: str, *, loops: Optional[Dict] = None, meta: Optional[Dict] = None) -> Any:
    """Значение адреса по StateView (сериализуемое) или None, если недоступно/не найдено.
    `loops` — Scheduler.loops (имя→LoopSpec) для `loop.*`; `meta` — мета движка для `meta.*`.
    Принимает и старые имена (псевдонимы)."""
    try:
        p = parse(addr)
    except ValueError:
        return None
    try:
        if p.area == "obs":
            arr = view.read(p.obj)
            return None if arr is None else _reduce(arr, p.red)
        if p.area == "act":
            return _resolve_act(view, p)
        if p.area == "undulator":
            return _resolve_undulator(view, p)
        if p.area == "optics":
            return _resolve_optics(view, p)
        if p.area == "element":
            return _resolve_element(view, p)
        if p.area == "loop":
            return _resolve_loop(loops, p)
        if p.area == "meta":
            return _ser(meta.get(p.obj)) if meta else None
        if p.area == "ref":                                # уставка контура (state.refs): массив или скаляр
            arr = view.ref(p.obj)
            return None if arr is None else _reduce(arr, p.red)
    except (KeyError, IndexError, ValueError, AttributeError, TypeError):
        return None
    return None                                            # screen / diag / rec — не с доски


def _resolve_act(view, p: Parsed):
    if p.attr:                                             # привод актуатора (Limits)
        d = view.drives.of(p.obj)                          # {slew,lag_tau,odz:(lo,hi),res,act}
        if p.attr == "odz_lo":
            return float(d["odz"][0])
        if p.attr == "odz_hi":
            return float(d["odz"][1])
        if p.attr == "act":
            return bool(d["act"])
        return float(d[p.attr])
    sv = getattr(view, "actual", None)                     # реализованное значение (StateVector по id)
    if sv is None:
        return None
    return float(sv[p.obj])


def _resolve_undulator(view, p: Parsed):
    u = view.undulators.get(p.obj)
    if u is None:
        return None
    a = p.attr
    if a == "gap":
        return float(u.gap_mm)
    if a == "gap_target":
        return float(u.gap_target)
    if a == "slew":
        return float(u.slew_mm_s)
    if a == "orbit_kick":
        return float(u.orbit_kick_rad())
    if a == "K":
        return _ser(u.K())                                 # None, если модуль undulator_sr недоступен
    if a == "B0":
        return _ser(u.B0())
    return None


def _resolve_optics(view, p: Parsed):
    opt = view.optics
    if opt is None:
        return None
    if p.at:
        key = OPTICS_AT.get(p.attr)
        if key is None:
            return None
        return _reduce(opt.at_magnet(p.at).get(key), p.red)
    if p.attr in OPTICS_SCALARS:
        return float(getattr(opt, p.attr))
    return None


def _resolve_element(view, p: Parsed):
    mv = view.machine
    if mv is None:
        return None
    if p.attr in ("x", "y", "z"):
        return float(mv.position(p.obj)[{"x": 0, "y": 1, "z": 2}[p.attr]])
    if p.attr in ("theta0", "k1", "k2", "k3"):
        return float(mv.field(p.obj)[p.attr])
    if p.attr in ("yaw", "pitch", "roll"):                 # углы ориентации прибора (рад): азимут + наклоны
        return float(mv.orientation(p.obj)[p.attr])
    return None


def _resolve_loop(loops, p: Parsed):
    if not loops:
        return None
    spec = loops.get(p.obj)
    if spec is None or p.attr not in LOOP_FIELDS:
        return None
    if p.attr == "enabled":
        return bool(spec.enabled)
    if p.attr == "period":
        return float(spec.period)
    if p.attr == "gain":                                   # усиление контура — Controller.params["gain"]
        g = (getattr(spec.controller, "params", None) or {}).get("gain")
        return None if g is None else float(g)
    return float(spec.rate_hz)


# ------------------------------------------------------------------ перечисление
def enumerate_addresses(view, *, loops: Optional[Dict] = None, meta: Optional[Dict] = None) -> List[str]:
    """Все КАНОНИЧЕСКИЕ адреса, разрешимые из текущего состояния и достойные списка (дерево проекта, подписка).
    Поэлементные `element.*` и `optics.*@dev` НЕ перечисляются поголовно (их тысячи) — они резолвятся по
    запросу через структуру (3D/дерево); здесь — наблюдаемые, приводы, ондуляторы, скаляры оптики, контуры, мета."""
    out: List[str] = []
    s = getattr(view, "_s", view)                          # StateView(._s) или сырой ExperimentState
    for oid, val in getattr(s, "observables", {}).items():
        out.append(f"obs.{oid}")
        if val is not None and np.asarray(val).size > 1:
            out.append(f"obs.{oid}:rms")                   # удобный скалярный по умолчанию
    schema = getattr(s, "schema", None)
    if schema is not None:
        for aid in schema.ids:
            for fld in DRIVE_FIELDS:
                out.append(f"act.{aid}.{fld}")
    for gid in getattr(s, "undulators", {}):
        for attr in UND_ATTRS:
            out.append(f"undulator.{gid}.{attr}")
    if getattr(s, "optics", None) is not None:
        out += [f"optics.{k}" for k in OPTICS_SCALARS]
    for name in (loops or {}):
        for fld in LOOP_FIELDS:
            out.append(f"loop.{name}.{fld}")
    for k in (meta or {}):
        out.append(f"meta.{k}")
    return out


# ------------------------------------------------------------------ метаданные для UI и записи
_OBS_NAMES = {
    "ebpm.x": ("орбита X (e-BPM)", "м"), "ebpm.y": ("орбита Y (e-BPM)", "м"),
    "photon.x": ("свет X (photon)", "м"), "photon.y": ("свет Y (photon)", "м"),
    "tune": ("рабочая точка", ""), "chroma": ("хроматизм", ""), "adts": ("ADTS", "1/м"),
    "rf.delta": ("энергия δ̂ (оценка по орбите)", ""),
}
_DRIVE_UNITS = {"slew": "ед/с", "lag_tau": "с", "res": "ед", "act": "", "odz_lo": "ед", "odz_hi": "ед"}
_ACT_UNITS = {"hkick": "рад", "vkick": "рад", "dk": "1/м²", "dk2": "1/м³", "dk3": "1/м⁴", "freq": "Гц"}
_UND_UNITS = {"gap": "мм", "gap_target": "мм", "slew": "мм/с", "orbit_kick": "рад", "K": "", "B0": "Тл"}
_OPT_UNITS = {"beta_x": "м", "beta_y": "м", "alpha_x": "", "alpha_y": "",
              "eta_x": "м", "eta_xp": "", "nux": "", "nuy": "", "alphac": ""}
_ELEM_UNITS = {"x": "м", "y": "м", "z": "м", "theta0": "рад", "k1": "1/м²", "k2": "1/м³", "k3": "1/м⁴",
               "yaw": "рад", "pitch": "рад", "roll": "рад"}
_META_UNITS = {"model_time": "с", "wall_time": "с", "speed": "×реал", "cycle": "такт", "seq": "", "turn": "об"}
# ряды рекордеров (control/recorders.py): имя ряда → (подпись, единицы)
_REC_SERIES = {
    "orbit_rms": ("СКО орбиты", "мкм"), "orbit_dc": ("СКО статической орбиты", "мкм"),
    "light": ("СКО центроида света", "мкм"), "rf.df_hz": ("отклонение частоты ВЧ", "Гц"),
    "rf.circ_um": ("сдвиг длины кольца", "мкм"), "rf.delta_ppm": ("энергия δ̂", "ppm"),
    "rf.delta_true_ppm": ("энергия δ (истинная)", "ppm"), "tune": ("рабочая точка", ""),
    "chroma": ("хроматизм", ""), "adts": ("ADTS", "1/м"), "screen": ("экран", ""),
}


def _rec_meta(series: str):
    for k in sorted(_REC_SERIES, key=len, reverse=True):
        if series == k or series.startswith(k + "."):
            nm, un = _REC_SERIES[k]
            tail = series[len(k) + 1:] if series != k else ""
            if k == "screen":                              # screen.<величина>_<единица>
                return f"экран · {tail}", (tail.rsplit("_", 1)[-1] if "_" in tail else "")
            return (f"{nm} · {tail}" if tail else nm), un
    return f"ряд {series}", ""


def meta(addr: str) -> AddrMeta:
    """Метаданные адреса — чисто из разбора строки (без состояния). Принимает и старые имена; поле addr —
    КАНОНИЧЕСКАЯ форма."""
    a = canonical(addr)
    try:
        p = parse(a)
    except ValueError:
        return AddrMeta(a, a, "", W_SCALAR, a.split(".", 1)[0])
    k = p.area
    if k == "obs":
        nm, un = _OBS_NAMES.get(p.obj, (f"obs {p.obj}", "м"))
        if p.red:
            return AddrMeta(a, f"{nm} · {p.red}", un, W_GRAPH, k)
        return AddrMeta(a, nm, un, W_VECTOR, k, shape="vector")
    if k == "act":
        if p.attr:
            return AddrMeta(a, f"привод {p.obj} · {p.attr}", _DRIVE_UNITS.get(p.attr, ""), W_SCALAR, k,
                            writable=p.attr in WRITABLE_DRIVE)
        suf = p.obj.rsplit(".", 1)[-1]
        return AddrMeta(a, f"{p.obj} · реализ.", _ACT_UNITS.get(suf, ""), W_GRAPH, k)
    if k == "undulator":
        return AddrMeta(a, f"ондулятор {p.obj} · {p.attr}", _UND_UNITS.get(p.attr, ""), W_GRAPH, k,
                        writable=p.attr in WRITABLE_UNDULATOR)
    if k == "optics":
        body = p.attr + (f"@{p.at}" if p.at else "")
        return AddrMeta(a, f"оптика {body}", _OPT_UNITS.get(p.attr, ""), W_GRAPH, k)
    if k == "element":
        return AddrMeta(a, f"{p.obj} · {p.attr}", _ELEM_UNITS.get(p.attr, ""), W_GRAPH, k)
    if k == "loop":
        return AddrMeta(a, f"контур {p.obj} · {p.attr}", "", W_SCALAR, k, writable=p.attr in WRITABLE_LOOP)
    if k == "meta":
        return AddrMeta(a, f"движок · {p.obj}", _META_UNITS.get(p.obj, ""), W_GRAPH, k, owner="run")
    if k == "diag":
        return AddrMeta(a, DIAG_TITLES.get(p.obj, f"диагностика · {p.obj}"), "", W_INFO, k,
                        shape="info", owner="engine", cost="heavy")
    if k == "screen":
        return AddrMeta(a, f"экран {p.obj}", "", W_SCREEN, k, shape="map", owner="engine", cost="heavy")
    if k == "rec":
        nm, un = _rec_meta(p.obj)
        return AddrMeta(a, nm, un, W_GRAPH, k, owner="recorder")
    if k == "ref":
        base = p.obj.rsplit(".", 1)[0]
        nm, un = _OBS_NAMES.get(base, (base, ""))
        return AddrMeta(a, f"цель · {nm}" + (f" · {p.red}" if p.red else ""), un, W_GRAPH, k)
    if k == "cavity":
        title = {"volt_mv": "амплитуда V", "phase_deg": "фаза φ", "df_hz": "сдвиг частоты Δf (ручной)",
                 "freq_hz": "частота f", "state": "состояние резонатора", "bunch": "пучок на резонаторе"}
        from .longitudinal import PARAM_INDEX            # параметры окна «Резонатор ВЧ» (единицы — СИ, как в state)
        if p.obj not in title and p.obj in PARAM_INDEX:
            title[p.obj] = PARAM_INDEX[p.obj][1]
        return AddrMeta(a, f"резонатор · {title.get(p.obj, p.obj)}", CAVITY_UNITS.get(p.obj, ""),
                        W_INFO if p.obj in ("state", "bunch") else W_GRAPH, k,
                        writable=p.obj in CAVITY_WRITABLE, owner="engine",
                        shape="info" if p.obj in ("state", "bunch") else "scalar")
    if k == "design":
        from .edits import UNITS, TITLES, EDITABLE
        return AddrMeta(a, f"{p.obj} · {TITLES.get(p.attr, p.attr)} (проект)", UNITS.get(p.attr, ""), W_SCALAR, k,
                        writable=p.attr in EDITABLE, owner="machine", scope="project",
                        edit=EDITABLE.get(p.attr, ("", "cold"))[1])
    return AddrMeta(a, a, "", W_SCALAR, k)
