# -*- coding: utf-8 -*-
"""Документы SynchrotronSim — формат обмена библиотеки с хранилищем (docs/DATA_PROTOCOL.md, раздел 4).

Документ = КОНВЕРТ (JSON: формат, вид, id, ревизия, родитель, отпечаток, дата, заметка, содержимое) + числовые
СТОЛБЦЫ (бинарные массивы). Упаковка для диска и сети — zip: `document.json` + `arrays.npz` (без pickle).

    from synchrotronsim.documents import machine_to_bytes, machine_from_bytes
    b = machine_to_bytes(m)                       # Machine → байты (≈ несколько МБ у решётки 26)
    m2 = machine_from_bytes(b)                    # отпечаток сверяется: m2.fingerprint == m.fingerprint

Правила кодирования (стандартный JSON — читается не только Python):
  * float: конечное — число; NaN — null; ±inf — строка "inf"/"-inf"; при чтении тип поля берётся из dataclass;
  * кортежи — списки (обратно — кортежи по аннотации поля);
  * строковые столбцы — юникодные массивы в npz; dtype исходника (object) восстанавливается по `dtypes`.
"""
from __future__ import annotations

import io
import json
import math
import typing
import zipfile
from dataclasses import fields, is_dataclass
from datetime import datetime, timezone
from typing import Any, Dict, Optional, Tuple

import numpy as np

from .model import (Machine, InstanceTable, ElementTable, DeviceType, FieldMapRef, ElementGroup, Cavity,
                    PulsedKicker, BpmSpec, Screen, MachineMeta)

MACHINE_FORMAT = "synchrotronsim.machine/1"
_FORMAT_MAJOR = {"synchrotronsim.machine": 1}


class DocumentError(ValueError):
    """Документ не читается: незнакомый формат/старшая версия, нарушена целостность (отпечаток)."""


# ---------------------------------------------------------------------------------------------- dataclass ↔ JSON
def _enc(v):
    if is_dataclass(v) and not isinstance(v, type):
        return {f.name: _enc(getattr(v, f.name)) for f in fields(v)}
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer)):
        return int(v)
    if isinstance(v, (float, np.floating)):
        x = float(v)
        if math.isnan(x):
            return None
        if math.isinf(x):
            return "inf" if x > 0 else "-inf"
        return x
    if isinstance(v, (tuple, list)):
        return [_enc(x) for x in v]
    if isinstance(v, dict):
        return {str(k): _enc(x) for k, x in v.items()}
    if v is None or isinstance(v, str):
        return v
    raise TypeError(f"документ: неподдержанный тип {type(v).__name__}")


def _dec(tp, v):
    """Значение JSON → значение типа tp (аннотация поля dataclass)."""
    origin = typing.get_origin(tp)
    args = typing.get_args(tp)
    if origin is typing.Union:                                   # Optional[X]
        if v is None and type(None) in args:
            return None
        tp = next(a for a in args if a is not type(None))
        return _dec(tp, v)
    if tp is float:
        if v is None:
            return math.nan
        return float(v)
    if tp is int:
        return int(v)
    if tp is bool:
        return bool(v)
    if tp is str:
        return "" if v is None else str(v)
    if origin in (tuple, typing.Tuple):
        if len(args) == 2 and args[1] is Ellipsis:
            return tuple(_dec(args[0], x) for x in v)
        return tuple(_dec(a, x) for a, x in zip(args, v))
    if origin in (dict, typing.Dict):
        return {k: _dec(args[1], x) for k, x in v.items()}
    if is_dataclass(tp):
        return from_json_dc(tp, v)
    return v


def from_json_dc(cls, d: Dict[str, Any]):
    hints = typing.get_type_hints(cls)
    kw = {f.name: _dec(hints[f.name], d[f.name]) for f in fields(cls) if f.init and f.name in d}
    return cls(**kw)


# ---------------------------------------------------------------------------------------------- упаковка
def pack(envelope: Dict[str, Any], arrays: Dict[str, np.ndarray]) -> bytes:
    """Конверт + столбцы → байты (zip: document.json, arrays.npz)."""
    buf = io.BytesIO()
    np.savez(buf, **arrays)
    out = io.BytesIO()
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("document.json", json.dumps(envelope, ensure_ascii=False, allow_nan=False, indent=1))
        z.writestr("arrays.npz", buf.getvalue())
    return out.getvalue()


def unpack(data: bytes) -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        env = json.loads(z.read("document.json").decode("utf-8"))
        with np.load(io.BytesIO(z.read("arrays.npz")), allow_pickle=False) as npz:
            arrays = {k: npz[k] for k in npz.files}
    return env, arrays


def check_format(env: Dict[str, Any], family: str) -> None:
    fmt = str(env.get("format", ""))
    name, _, ver = fmt.partition("/")
    if name != family:
        raise DocumentError(f"ожидался документ {family}, получен {fmt or 'без формата'}")
    try:
        major = int(ver.split(".")[0])
    except ValueError:
        raise DocumentError(f"неразборчивая версия формата {fmt}")
    if major > _FORMAT_MAJOR[family]:
        raise DocumentError(f"формат {fmt} новее, чем понимает эта версия SynchrotronSim "
                            f"(до {family}/{_FORMAT_MAJOR[family]})")


def now_iso() -> str:
    return datetime.now(timezone.utc).astimezone().isoformat(timespec="seconds")


# ---------------------------------------------------------------------------------------------- машина
def _table_arrays(prefix: str, t, arrays: dict, dtypes: dict):
    for c in t.column_names():
        a = np.asarray(getattr(t, c))
        key = f"{prefix}.{c}"
        if a.dtype.kind == "O":                                   # строки: юникодный массив без pickle
            arrays[key] = np.asarray(["" if s is None else str(s) for s in a], dtype=str)
            dtypes[key] = "str"
        else:
            arrays[key] = a
            dtypes[key] = a.dtype.str


def _table_from(cls, prefix: str, arrays: dict, dtypes: dict):
    kw = {}
    for c in cls._COLUMNS:
        key = f"{prefix}.{c}"
        if key not in arrays:
            raise DocumentError(f"в документе нет столбца {key}")
        a = arrays[key]
        if dtypes.get(key) == "str":
            a = np.array(a.tolist(), dtype=object)
        kw[c] = a
    return cls(**kw)


def machine_document(m: Machine, machine_id: Optional[str] = None, parent: Optional[str] = None,
                     note: str = "", author: str = "") -> Tuple[Dict[str, Any], Dict[str, np.ndarray]]:
    """Machine → (конверт, столбцы). Ревизия = отпечаток содержимого."""
    arrays: Dict[str, np.ndarray] = {}
    dtypes: Dict[str, str] = {}
    _table_arrays("instances", m.instances, arrays, dtypes)
    _table_arrays("elements", m.elements, arrays, dtypes)
    content = {
        "lattice_id": int(m.lattice_id),
        "devices": {k: _enc(m.devices[k]) for k in sorted(m.devices)},
        "groups": _enc(m.groups), "cavities": _enc(m.cavities), "kickers": _enc(m.kickers),
        "bpms": _enc(m.bpms), "screens": _enc(m.screens), "meta": _enc(m.meta),
        "dtypes": dtypes,
    }
    env = {"format": MACHINE_FORMAT, "kind": "machine", "id": str(machine_id or m.lattice_id),
           "rev": m.fingerprint, "parent": parent, "fingerprint": m.fingerprint, "created": now_iso(),
           "author": author, "note": note, "source": m.source, "content": content}
    return env, arrays


def machine_from_document(env: Dict[str, Any], arrays: Dict[str, np.ndarray], source: str = "",
                          verify: bool = True) -> Machine:
    """(конверт, столбцы) → Machine. verify — пересчитать отпечаток и сверить с записанным."""
    check_format(env, "synchrotronsim.machine")
    c = env["content"]
    dt = c.get("dtypes") or {}
    m = Machine(
        lattice_id=int(c["lattice_id"]),
        instances=_table_from(InstanceTable, "instances", arrays, dt),
        elements=_table_from(ElementTable, "elements", arrays, dt),
        devices={k: from_json_dc(DeviceType, v) for k, v in (c.get("devices") or {}).items()},
        groups=tuple(from_json_dc(ElementGroup, g) for g in c.get("groups") or ()),
        cavities=tuple(from_json_dc(Cavity, g) for g in c.get("cavities") or ()),
        kickers=tuple(from_json_dc(PulsedKicker, g) for g in c.get("kickers") or ()),
        bpms=tuple(from_json_dc(BpmSpec, g) for g in c.get("bpms") or ()),
        screens=tuple(from_json_dc(Screen, g) for g in c.get("screens") or ()),
        meta=from_json_dc(MachineMeta, c["meta"]) if c.get("meta") else None,
        source=source or f"document:{env.get('id')}@{str(env.get('rev', ''))[:12]}")
    want = env.get("fingerprint")
    if verify and want and m.fingerprint != want:
        raise DocumentError(f"документ машины {env.get('id')}: отпечаток содержимого {m.fingerprint[:12]}… "
                            f"не совпадает с записанным {str(want)[:12]}… — документ повреждён или изменён")
    return m


def machine_to_bytes(m: Machine, **kw) -> bytes:
    return pack(*machine_document(m, **kw))


def machine_from_bytes(data: bytes, **kw) -> Machine:
    return machine_from_document(*unpack(data), **kw)


def save_machine(m: Machine, path: str, **kw) -> str:
    with open(path, "wb") as f:
        f.write(machine_to_bytes(m, **kw))
    return path


def load_machine(path: str, verify: bool = True) -> Machine:
    with open(path, "rb") as f:
        return machine_from_bytes(f.read(), source=f"document:{path}", verify=verify)
