# -*- coding: utf-8 -*-
"""Сетевой формат SynchrotronSim (docs/DATA_PROTOCOL.md, этап 3) — без pickle и без сторонних библиотек.

Сообщение = магия `SSW1` + длины + JSON-заголовок (UTF-8) + бинарные вложения. В заголовке — обычный JSON;
то, чего в JSON нет, кодируется метками: массив numpy {"$nd": [dtype, форма, № вложения]}, байты {"$b": №},
кортеж {"$t": [...]}, NaN/±inf {"$f": "nan"|"inf"|"-inf"}, словарь с нестроковыми ключами {"$m": [[k, v], …]}.
Поток (вычислитель, живой кокпит) — кадры: 8 байт длины + сообщение.
"""
from __future__ import annotations

import json
import math
import struct
from typing import Any, List

import numpy as np

MAGIC = b"SSW1"
PROTOCOL = "synchrotronsim.wire/1"
MAX_FRAME = 1 << 31                         # 2 ГБ — предохранитель от мусора в потоке


def encode(obj: Any) -> bytes:
    att: List[bytes] = []

    def enc(o):
        if o is None or isinstance(o, (bool, str)):
            return o
        if isinstance(o, (int, np.integer)) and not isinstance(o, np.bool_):
            return int(o)
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, (float, np.floating)):
            x = float(o)
            if math.isnan(x):
                return {"$f": "nan"}
            if math.isinf(x):
                return {"$f": "inf" if x > 0 else "-inf"}
            return x
        if isinstance(o, np.ndarray):
            if o.dtype.kind in "OUS":
                return {"$t": [enc(x) for x in o.tolist()]} if o.ndim else enc(o.item())
            att.append(np.ascontiguousarray(o).tobytes())
            return {"$nd": [o.dtype.str, list(o.shape), len(att) - 1]}
        if isinstance(o, (bytes, bytearray, memoryview)):
            att.append(bytes(o))
            return {"$b": len(att) - 1}
        if isinstance(o, tuple):
            return {"$t": [enc(x) for x in o]}
        if isinstance(o, (list, set, frozenset)):
            return [enc(x) for x in o]
        if isinstance(o, dict):
            if all(isinstance(k, str) for k in o):
                return {k: enc(v) for k, v in o.items()}
            return {"$m": [[enc(k), enc(v)] for k, v in o.items()]}
        if hasattr(o, "__dataclass_fields__"):
            from dataclasses import fields
            return {f.name: enc(getattr(o, f.name)) for f in fields(o)}
        raise TypeError(f"wire: неподдержанный тип {type(o).__name__}")

    head = json.dumps(enc(obj), ensure_ascii=False, allow_nan=False, separators=(",", ":")).encode("utf-8")
    parts = [MAGIC, struct.pack("<II", len(head), len(att)), head]
    for a in att:
        parts.append(struct.pack("<Q", len(a)))
        parts.append(a)
    return b"".join(parts)


def decode(data: bytes) -> Any:
    mv = memoryview(data)
    if bytes(mv[:4]) != MAGIC:
        raise ValueError("wire: не сообщение SynchrotronSim (нет метки SSW1)")
    nh, na = struct.unpack_from("<II", mv, 4)
    pos = 12
    head = json.loads(bytes(mv[pos:pos + nh]).decode("utf-8"))
    pos += nh
    att = []
    for _ in range(na):
        (n,) = struct.unpack_from("<Q", mv, pos)
        pos += 8
        att.append(bytes(mv[pos:pos + n]))
        pos += n

    def dec(o):
        if isinstance(o, list):
            return [dec(x) for x in o]
        if isinstance(o, dict):
            if len(o) == 1:
                (k, v), = o.items()
                if k == "$nd":
                    dt, shape, i = v
                    return np.frombuffer(att[i], dtype=np.dtype(dt)).reshape(shape).copy()
                if k == "$b":
                    return att[v]
                if k == "$t":
                    return tuple(dec(x) for x in v)
                if k == "$f":
                    return float(v)
                if k == "$m":
                    return {_hashable(dec(a)): dec(b) for a, b in v}
            return {k: dec(v) for k, v in o.items()}
        return o
    return dec(head)


def _hashable(x):
    return tuple(_hashable(i) for i in x) if isinstance(x, list) else x


# ---------------------------------------------------------------------------------------------- кадры потока
def send_frame(sock, obj) -> None:
    data = encode(obj)
    sock.sendall(struct.pack("<Q", len(data)) + data)


def _recv_exact(sock, n: int) -> bytes:
    buf = bytearray()
    while len(buf) < n:
        chunk = sock.recv(min(n - len(buf), 1 << 20))
        if not chunk:
            raise ConnectionError("соединение закрыто")
        buf += chunk
    return bytes(buf)


def recv_frame(sock) -> Any:
    (n,) = struct.unpack("<Q", _recv_exact(sock, 8))
    if n > MAX_FRAME:
        raise ConnectionError(f"wire: кадр {n} байт — больше предела")
    return decode(_recv_exact(sock, n))
