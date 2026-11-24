# -*- coding: utf-8 -*-
"""Модель данных SynchrotronSim — слой 1 «Установка» (docs/LIBRARY_ARCHITECTURE.md, раздел 3.3).

Объекты НЕИЗМЕНЯЕМЫ (frozen dataclass, массивы только для чтения): правка машины даёт новую машину с новым
отпечатком. Таблицы элементов хранятся ПО СТОЛБЦАМ (numpy), чтобы почти без преобразования уходить в C++.

Единицы в памяти — СИ: м, рад, Тл (mat_B), Тл/м (mat_G, mat_Gv), Тл/м² (mat_S), Тл/м³ (mat_O), с; доли вместо
процентов. Исключения по смыслу, с единицей в имени: энергия `*_mev` (МэВ — конвенция проекта), напряжение ВЧ
`rf_volt_mv` (МВ). Перевод из единиц таблиц (суффиксы `_mm`, `_mrad`, `_deg`, `_pct`) — только в `synchrotronsim.db`.

Два уровня описания магнитов:
  * `InstanceTable` — строки `lattice_configs` решётки в порядке обхода (как в БД: обмотка = строка);
  * `ElementTable`  — физические магниты: подряд идущие строки с одинаковыми x, y, z — обмотки ОДНОГО
    магнита, их mat_* складываются (так же, как `matrix_tracker.build_line` и C++ `buildLine`).
    Связь — `member_start`/`member_count` (элемент i = строки instances[member_start[i] : +member_count[i]]).

Отпечаток (`fingerprint`) — sha256 канонической кодировки содержимого (порядок важен: это порядок пучка).
В него входят строки решётки, использованные типы устройств, группы и устройства на кольце; производная
ElementTable в отпечаток не входит (она однозначно следует из строк).
"""
from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass, field, fields
from typing import Dict, Iterable, Optional, Tuple

import numpy as np

QP = -84.925  # q/p₀ опорной жёсткости, 1/(Тл·м) — справочно (как matrix_tracker.QP); модель хранит ПОЛЯ, не силы


# ---------------------------------------------------------------------------------------------- отпечатки
_NAN = struct.pack("<Q", 0x7FF8000000000000)


class _Hasher:
    """Каноническая кодировка для отпечатка: тип + длина + данные; NaN приводится к одному битовому образу."""

    def __init__(self, tag: str):
        self.h = hashlib.sha256()
        self.str(tag)

    def _raw(self, b: bytes):
        self.h.update(struct.pack("<Q", len(b))); self.h.update(b)

    def str(self, s: str):
        self.h.update(b"s"); self._raw(s.encode("utf-8"))

    def num(self, x: float):
        self.h.update(b"f")
        x = float(x)
        self._raw(_NAN if x != x else struct.pack("<d", x))

    def arr(self, a: np.ndarray):
        a = np.asarray(a)
        if a.dtype.kind in "fc":
            v = np.ascontiguousarray(a, "<f8").copy()
            v[np.isnan(v)] = np.nan                      # единый образ NaN (знак/полезная нагрузка не важны)
            self.h.update(b"F"); self.h.update(struct.pack("<Q", v.size)); self._raw(v.tobytes())
        elif a.dtype.kind in "iub":
            self.h.update(b"I"); self._raw(np.ascontiguousarray(a, "<i8").tobytes())
        else:                                            # строки / объекты
            self.h.update(b"S"); self.h.update(struct.pack("<Q", a.size))
            for s in a.ravel():
                self.str("" if s is None else str(s))

    def obj(self, o):
        """dataclass / кортеж / словарь / скаляр — рекурсивно, поля в порядке объявления."""
        if o is None:
            self.h.update(b"N")
        elif isinstance(o, bool):
            self.h.update(b"B1" if o else b"B0")
        elif isinstance(o, (int, np.integer)):
            self.h.update(b"i"); self._raw(struct.pack("<q", int(o)))
        elif isinstance(o, (float, np.floating)):
            self.num(o)
        elif isinstance(o, str):
            self.str(o)
        elif isinstance(o, np.ndarray):
            self.arr(o)
        elif isinstance(o, dict):
            self.h.update(b"D"); self.h.update(struct.pack("<Q", len(o)))
            for k in sorted(o):
                self.str(str(k)); self.obj(o[k])
        elif isinstance(o, (tuple, list)):
            self.h.update(b"T"); self.h.update(struct.pack("<Q", len(o)))
            for v in o:
                self.obj(v)
        elif hasattr(o, "__dataclass_fields__"):
            self.str(type(o).__name__)
            for f in fields(o):
                if f.metadata.get("fingerprint", True):
                    self.str(f.name); self.obj(getattr(o, f.name))
        else:
            raise TypeError(f"отпечаток: неподдержанный тип {type(o).__name__}")

    def hexdigest(self) -> str:
        return self.h.hexdigest()


def fingerprint_of(tag: str, *parts) -> str:
    h = _Hasher(tag)
    for p in parts:
        h.obj(p)
    return h.hexdigest()


def _ro(a: np.ndarray) -> np.ndarray:
    a = np.asarray(a)
    a.flags.writeable = False
    return a


class _Columns:
    """Основа столбцовой таблицы: поля — numpy-массивы одной длины, только для чтения."""
    _COLUMNS: Tuple[str, ...] = ()

    def __post_init__(self):
        n = None
        for c in self._COLUMNS:
            a = _ro(getattr(self, c))
            object.__setattr__(self, c, a)
            if n is None:
                n = len(a)
            elif len(a) != n:
                raise ValueError(f"{type(self).__name__}: столбец {c} длины {len(a)}, ожидалось {n}")

    def __len__(self):
        return len(getattr(self, self._COLUMNS[0])) if self._COLUMNS else 0

    def column_names(self):
        return self._COLUMNS

    def as_dict(self) -> Dict[str, np.ndarray]:
        return {c: getattr(self, c) for c in self._COLUMNS}


# ---------------------------------------------------------------------------------------------- решётка
@dataclass(frozen=True, eq=False)
class InstanceTable(_Columns):
    """Строки `lattice_configs` одной решётки в порядке обхода (сортировка `instance_id` как целого).

    Пустая ячейка числового столбца — NaN (в т.ч. `mat_Gv`: пусто ⇒ берётся `mat_G`; это правило применяет
    сборка ElementTable). `n_slices` пусто ⇒ 0."""
    instance_id: np.ndarray        # str
    device_type: np.ndarray        # str
    kind: np.ndarray               # str (как в таблице; пусто — тип решает element_transforms)
    x: np.ndarray                  # м
    y: np.ndarray
    z: np.ndarray
    yaw: np.ndarray                # рад
    pitch: np.ndarray
    roll: np.ndarray
    arg_val: np.ndarray            # ток/аргумент карты (для диффур-ядра)
    correct_coeff: np.ndarray
    mat_B: np.ndarray              # Тл
    mat_G: np.ndarray              # Тл/м
    mat_Gv: np.ndarray             # Тл/м (NaN ⇒ = mat_G)
    mat_S: np.ndarray              # Тл/м²
    mat_O: np.ndarray              # Тл/м³
    eff_length: np.ndarray         # м (NaN ⇒ из hardware_library)
    n_slices: np.ndarray           # int (0 ⇒ из element_transforms / настроек опыта)
    edge1: np.ndarray              # рад (как в БД: E хранится как −E .lte)
    edge2: np.ndarray
    lf1: np.ndarray                # м, длина спада поля торца L_f = 2·HGAP·FINT (0 — жёсткий край)
    lf2: np.ndarray

    _COLUMNS = ("instance_id", "device_type", "kind", "x", "y", "z", "yaw", "pitch", "roll", "arg_val",
                "correct_coeff", "mat_B", "mat_G", "mat_Gv", "mat_S", "mat_O", "eff_length", "n_slices",
                "edge1", "edge2", "lf1", "lf2")


@dataclass(frozen=True, eq=False)
class ElementTable(_Columns):
    """Физические магниты в порядке пучка (обмотки одного магнита сложены). Поля — суммы mat_* обмоток,
    геометрия и структурные метаданные — по первой обмотке, с фолбэком на тип устройства (как build_line)."""
    instance_id: np.ndarray        # str — id первой обмотки (у датчика — его sensor_id)
    device_type: np.ndarray        # str — тип первой обмотки
    kind: np.ndarray               # str — dipole/dq/quad/sext/oct/kicker/monitor/…
    x: np.ndarray                  # м, центр «ящика» (у изгиба — середина дуги)
    y: np.ndarray
    z: np.ndarray
    yaw: np.ndarray                # рад, Tait–Bryan (каждый поворот вокруг уже повёрнутой оси)
    pitch: np.ndarray
    roll: np.ndarray
    eff_length: np.ndarray         # м
    n_slices: np.ndarray           # int, срезы box-тела; 0 ⇒ настройка опыта (Settings.n_slices)
    mat_B: np.ndarray              # Тл (сумма обмоток)
    mat_G: np.ndarray              # Тл/м
    mat_Gv: np.ndarray             # Тл/м
    mat_S: np.ndarray              # Тл/м²
    mat_O: np.ndarray              # Тл/м³
    edge1: np.ndarray              # рад
    edge2: np.ndarray
    lf1: np.ndarray                # м
    lf2: np.ndarray
    member_start: np.ndarray       # int — первая строка InstanceTable
    member_count: np.ndarray       # int — число обмоток (строк)

    _COLUMNS = ("instance_id", "device_type", "kind", "x", "y", "z", "yaw", "pitch", "roll", "eff_length",
                "n_slices", "mat_B", "mat_G", "mat_Gv", "mat_S", "mat_O", "edge1", "edge2", "lf1", "lf2",
                "member_start", "member_count")

    @property
    def centers(self) -> np.ndarray:
        return np.stack([self.x, self.y, self.z], axis=1)

    def members(self, i: int) -> slice:
        s = int(self.member_start[i])
        return slice(s, s + int(self.member_count[i]))

    def index_of(self, instance_id: str) -> int:
        hit = np.nonzero(self.instance_id == str(instance_id))[0]
        if not len(hit):
            raise KeyError(f"нет элемента с instance_id={instance_id}")
        return int(hit[0])

    def where(self, kind: str) -> np.ndarray:
        """Индексы элементов данного вида."""
        return np.nonzero(self.kind == kind)[0]


# ---------------------------------------------------------------------------------------------- устройства
@dataclass(frozen=True)
class FieldMapRef:
    """Строка реестра карт поля `field_maps` (несколько карт на устройство)."""
    map_id: str
    element: str
    winding: str
    amp_turns: float               # А·витки (NaN — не задано)
    current_a: str                 # как в таблице (может быть составным: «I1=219.25»)
    excitation: str
    split: str                     # train / test
    grid: str
    material_curve: str
    raw_path: str
    bin_path: str
    source_fld: str
    note: str


@dataclass(frozen=True)
class DeviceType:
    """Тип устройства: строка `hardware_library` + `element_transforms` (+ карты `field_maps`)."""
    device_type_id: str
    active: bool
    field_type: str
    kind: str                      # из element_transforms (пусто — не задан)
    active_coeffs: str
    n_slices: int                  # из element_transforms (0 — не задан)
    eff_length: float              # м
    nominal_arg: float
    field_direction: float
    correct_coeff: float
    edge_ff: float                 # форм-фактор краевого фринджа (0 — выключен)
    edge_delta: float
    use_symmetry: bool
    synthetic: bool
    scaled_from: str
    map_raw_path: str
    map_bin_path: str
    stl_path: str
    turns_per_coil: float
    num_coils: float
    total_turns: float
    nominal_current: float         # А
    rot_x: float                   # рад (в таблице — градусы, rot_*_deg)
    rot_y: float
    rot_z: float
    notes: str = ""
    maps: Tuple[FieldMapRef, ...] = ()


# ---------------------------------------------------------------------------------------------- устройства на кольце
@dataclass(frozen=True)
class ElementGroup:
    """Жёсткая группа элементов (`element_groups`): ондулятор und_NN и т.п."""
    group_id: str
    kind: str
    pivot: Tuple[float, float, float]      # м
    yaw: float                              # рад
    pitch: float
    roll: float
    instance_ids: Tuple[str, ...]


@dataclass(frozen=True)
class Cavity:
    """ВЧ-резонатор установки (`rf_cavities`): плоскость (точка + нормаль) и доля напряжения."""
    cavity_id: str
    point: Tuple[float, float, float]       # м
    normal: Tuple[float, float, float]
    volt_share: float                       # доля
    note: str = ""


@dataclass(frozen=True)
class PulsedKicker:
    """Импульсный кикер установки (`pulsed_kickers`); удар и момент — в опыте (`kicker_timing`)."""
    kicker_id: str
    point: Tuple[float, float, float]
    normal: Tuple[float, float, float]
    shape: str
    duration: float                         # с
    rise: float                             # с
    note: str = ""


@dataclass(frozen=True)
class BpmSpec:
    """Несовершенства датчика положения (`bpm_config`)."""
    sensor_id: str
    dx: float
    dy: float
    tilt: float
    gain_x: float
    gain_y: float
    setpoint_x: float
    setpoint_y: float


@dataclass(frozen=True)
class Screen:
    """Фотонный экран (`screens`): центр, нормаль и оси плоскости; размеры и углы — в СИ."""
    screen_id: str
    source_ref: str
    source_kind: str
    center: Tuple[float, float, float]      # м
    normal: Tuple[float, float, float]
    u: Tuple[float, float, float]
    v: Tuple[float, float, float]
    width: float                            # м (в таблице width_mm)
    height: float                           # м
    mode: str
    dist: float                             # м
    extract: float                          # рад (в таблице *_mrad)
    los_h: float
    los_v: float
    tilt_h: float
    tilt_v: float
    roll: float


@dataclass(frozen=True)
class ErrorSigma:
    """Допуск по типу магнита (`element_errors`): σ сдвигов (м), поворотов (рад) и поля (доля)."""
    kind: str
    dx: float
    dy: float
    dz: float
    dyaw: float
    dpitch: float
    droll: float
    field_rel: float                        # доля (в таблице field_pct — проценты)


@dataclass(frozen=True)
class MachineMeta:
    """Паспортные данные установки. Сейчас живут в коде (`facility_report.LATTICE_INFO`), не в БД."""
    label: str
    design_energy_mev: float
    rf_volt_mv: float
    rf_harmonic: int
    superperiods: int


# ---------------------------------------------------------------------------------------------- установка
@dataclass(frozen=True, eq=False)
class Machine:
    """Установка: решётка + типы устройств + устройства на кольце. Неизменяема; см. `fingerprint`."""
    lattice_id: int
    instances: InstanceTable
    elements: ElementTable
    devices: Dict[str, DeviceType]                  # только типы, встречающиеся в решётке
    groups: Tuple[ElementGroup, ...] = ()
    cavities: Tuple[Cavity, ...] = ()
    kickers: Tuple[PulsedKicker, ...] = ()
    bpms: Tuple[BpmSpec, ...] = ()
    screens: Tuple[Screen, ...] = ()
    meta: Optional[MachineMeta] = None
    source: str = field(default="", metadata={"fingerprint": False})   # откуда загружено (папка таблиц)
    fingerprint: str = field(default="", metadata={"fingerprint": False})

    def __post_init__(self):
        object.__setattr__(self, "devices", dict(self.devices))
        if not self.fingerprint:
            object.__setattr__(self, "fingerprint", self.compute_fingerprint())

    def compute_fingerprint(self) -> str:
        return fingerprint_of("Machine/1", self.lattice_id, self.instances,
                              {k: self.devices[k] for k in sorted(self.devices)}, self.groups, self.cavities,
                              self.kickers, self.bpms, self.screens, self.meta)

    @property
    def n_elements(self) -> int:
        return len(self.elements)

    def kinds(self) -> Dict[str, int]:
        """Число магнитов по видам."""
        k, n = np.unique(self.elements.kind, return_counts=True)
        return {str(a): int(b) for a, b in zip(k, n)}

    def __repr__(self):
        label = self.meta.label if self.meta else "—"
        return (f"Machine(lattice={self.lattice_id}, «{label}», {len(self.elements)} магнитов из "
                f"{len(self.instances)} строк, {len(self.devices)} типов, групп {len(self.groups)}, "
                f"экранов {len(self.screens)}, fingerprint={self.fingerprint[:12]}…)")
