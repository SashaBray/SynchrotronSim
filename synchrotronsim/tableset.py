# -*- coding: utf-8 -*-
"""Таблицы как содержимое ревизий хранилища (docs/DATA_PROTOCOL.md, этап 2).

Ревизия машины хранит СРЕЗ её таблиц — строки с её lattice_id, ячейки как текст, без потерь (пусто и «-» —
разные вещи, единицы таблиц, столбцы вне модели — всё сохраняется). Документ `Machine` для счёта строится из
среза тем же разбором, что и из CSV (`MemoryRepository`). Справочник устройств (hardware_library,
element_transforms, field_maps) — отдельный вид документа «library», общий для машин; ревизия машины
ссылается на ревизию справочника.

Рабочая копия — папка `workspace/tables` (как в git: скрипты и Excel работают с ней; история — в хранилище).
`write_*` вливают срез обратно: строки других решёток не трогаются, таблица без изменений не переписывается,
чтение → запись воспроизводит файл байт-в-байт (проверено на всех таблицах проекта).
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

MACHINE_TABLES = ("lattice_configs", "element_connectivity", "element_groups", "rf_cavities", "pulsed_kickers",
                  "bpm_config", "screens")
LIBRARY_TABLES = ("hardware_library", "element_transforms", "field_maps")
# ключ строки для разницы (по нему строки сопоставляются между ревизиями)
ROW_KEY = {"lattice_configs": ("instance_id",), "element_connectivity": ("instance_id",),
           "element_groups": ("group_id",), "rf_cavities": ("cavity_id",), "pulsed_kickers": ("kicker_id",),
           "bpm_config": ("id_sensor",), "screens": ("screen_id",), "hardware_library": ("device_type_id",),
           "element_transforms": ("device_type_id",), "field_maps": ("map_id",)}
SLICE_FORMAT = "synchrotronsim.tables/1"


# ---------------------------------------------------------------------------------------------- лист CSV
@dataclass
class Sheet:
    """Файл таблицы: заголовок + строки (текст ячеек как есть) + свойства файла (BOM, перевод строки)."""
    header: List[str]
    rows: List[List[str]]
    bom: bool = False
    newline: str = "\r\n"

    def col(self, name: str) -> int:
        return self.header.index(name)


def read_sheet(path: str) -> Sheet:
    with open(path, "rb") as f:
        raw = f.read()
    bom = raw.startswith(b"\xef\xbb\xbf")
    txt = raw.decode("utf-8-sig")
    nl = "\r\n" if "\r\n" in txt else "\n"
    rows = list(csv.reader(io.StringIO(txt, newline=""), delimiter=";"))
    return Sheet(rows[0] if rows else [], rows[1:], bom, nl)


def sheet_bytes(s: Sheet) -> bytes:
    out = io.StringIO(newline="")
    w = csv.writer(out, delimiter=";", lineterminator=s.newline)
    w.writerow(s.header)
    w.writerows(s.rows)
    return (b"\xef\xbb\xbf" if s.bom else b"") + out.getvalue().encode("utf-8")


def write_sheet(path: str, s: Sheet) -> None:
    data = sheet_bytes(s)
    tmp = path + ".part"
    with open(tmp, "wb") as f:
        f.write(data)
    os.replace(tmp, path)


def _nonempty(r: List[str]) -> bool:
    return bool(r) and any(c.strip() for c in r)


# ---------------------------------------------------------------------------------------------- срезы
def tables_dir(folder: str) -> str:
    """Папка таблиц: сама папка, если в ней есть tables/ — то tables/."""
    t = os.path.join(folder, "tables")
    if os.path.isdir(t):
        return t
    w = os.path.join(folder, "workspace", "tables")
    return w if os.path.isdir(w) else folder


def lattices_in(folder: str) -> List[int]:
    p = os.path.join(tables_dir(folder), "lattice_configs.csv")
    s = read_sheet(p)
    i = s.col("lattice_id")
    return sorted({int(r[i]) for r in s.rows if _nonempty(r) and r[i].strip()})


def machine_slice(folder: str, lattice_id) -> Dict[str, Any]:
    """Срез машины из рабочей копии: строки её решётки во всех таблицах машины (+ паспорт из machines.py)."""
    from .machines import LATTICE_INFO
    tb = tables_dir(folder)
    lat = str(int(lattice_id))
    tables = {}
    for name in MACHINE_TABLES:
        p = os.path.join(tb, name + ".csv")
        if not os.path.exists(p):
            continue
        s = read_sheet(p)
        if "lattice_id" not in s.header:
            continue
        i = s.col("lattice_id")
        tables[name] = {"header": list(s.header), "rows": [r for r in s.rows if _nonempty(r) and len(r) > i
                                                            and r[i] == lat]}
    if not tables.get("lattice_configs", {}).get("rows"):
        raise KeyError(f"решётки {lat} нет в {tb}")
    info = LATTICE_INFO.get(int(lat))
    return {"format": SLICE_FORMAT, "kind": "machine", "id": lat, "tables": tables,
            "passport": {k: (list(v) if isinstance(v, tuple) else v) for k, v in info.items()} if info else None}


def library_slice(folder: str) -> Dict[str, Any]:
    """Справочник устройств рабочей копии — таблицы целиком (с BOM и переводом строки — для точной записи)."""
    tb = tables_dir(folder)
    tables = {}
    for name in LIBRARY_TABLES:
        p = os.path.join(tb, name + ".csv")
        if os.path.exists(p):
            s = read_sheet(p)
            tables[name] = {"header": list(s.header), "rows": s.rows, "bom": s.bom, "newline": s.newline}
    return {"format": SLICE_FORMAT, "kind": "library", "id": "devices", "tables": tables}


def canonical(content: Dict[str, Any]) -> bytes:
    return json.dumps(content, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def content_rev(content: Dict[str, Any]) -> str:
    """Ревизия = sha256 канонической кодировки содержимого (без даты/автора/заметки)."""
    c = {k: v for k, v in content.items() if k not in ("created", "author", "note", "parent")}
    return hashlib.sha256(canonical(c)).hexdigest()


# ---------------------------------------------------------------------------------------------- запись в рабочую копию
def write_machine(folder: str, sl: Dict[str, Any]) -> List[str]:
    """Влить срез машины в рабочую копию: строки этой решётки заменяются (на месте первой прежней строки),
    строки других решёток не трогаются; таблица без изменений не переписывается. → список изменённых таблиц."""
    tb = tables_dir(folder)
    lat = str(sl["id"])
    changed = []
    for name, t in sl["tables"].items():
        p = os.path.join(tb, name + ".csv")
        if os.path.exists(p):
            s = read_sheet(p)
        else:
            s = Sheet(list(t["header"]), [], False, "\r\n")
        rows_new = _remap(t["rows"], t["header"], s.header) if t["header"] != s.header else [list(r) for r in t["rows"]]
        if t["header"] != s.header:                               # столбцы ревизии, которых нет в файле — добавить
            extra = [c for c in t["header"] if c not in s.header]
            if extra:
                s.header = s.header + extra
                s.rows = [r + [""] * len(extra) for r in s.rows]
                rows_new = _remap(t["rows"], t["header"], s.header)
        i = s.col("lattice_id")
        old_idx = [k for k, r in enumerate(s.rows) if _nonempty(r) and len(r) > i and r[i] == lat]
        if [s.rows[k] for k in old_idx] == rows_new:
            continue
        old = set(old_idx)
        keep = [r for k, r in enumerate(s.rows) if k not in old]
        at = old_idx[0] if old_idx else len(keep)              # до первой прежней строки — только чужие строки
        s.rows = keep[:at] + rows_new + keep[at:]
        write_sheet(p, s)
        changed.append(name)
    return changed


def write_library(folder: str, sl: Dict[str, Any]) -> List[str]:
    tb = tables_dir(folder)
    changed = []
    for name, t in sl["tables"].items():
        p = os.path.join(tb, name + ".csv")
        s = Sheet(list(t["header"]), [list(r) for r in t["rows"]], bool(t.get("bom")), t.get("newline", "\r\n"))
        if os.path.exists(p) and read_sheet(p) == s:
            continue
        write_sheet(p, s)
        changed.append(name)
    return changed


def _remap(rows, src_header, dst_header):
    ix = {c: k for k, c in enumerate(src_header)}
    return [[(r[ix[c]] if c in ix and ix[c] < len(r) else "") for c in dst_header] for r in rows]


# ---------------------------------------------------------------------------------------------- разница
@dataclass
class TableDiff:
    table: str
    added: List[str] = field(default_factory=list)            # ключи строк
    removed: List[str] = field(default_factory=list)
    cells: List[Tuple[str, str, str, str]] = field(default_factory=list)   # (ключ, столбец, было, стало)
    header: Optional[Tuple[List[str], List[str]]] = None
    reordered: bool = False

    @property
    def empty(self) -> bool:
        return not (self.added or self.removed or self.cells or self.header or self.reordered)

    def n(self) -> int:
        return len(self.added) + len(self.removed) + len(self.cells) + (1 if self.header else 0)


def _keyed(t: Dict[str, Any], name: str):
    h = t["header"]
    cols = [h.index(c) for c in ROW_KEY.get(name, ()) if c in h]
    out, order = {}, []
    for k, r in enumerate(t["rows"]):
        if not _nonempty(r):
            continue
        key = "|".join(r[c] if c < len(r) else "" for c in cols) if cols else f"#{k}"
        if key in out:                                          # неуникальный ключ — с номером повтора
            j = 2
            while f"{key}~{j}" in out:
                j += 1
            key = f"{key}~{j}"
        out[key] = {c: (r[i] if i < len(r) else "") for i, c in enumerate(h)}
        order.append(key)
    return out, order


def diff_slices(a: Optional[Dict[str, Any]], b: Optional[Dict[str, Any]]) -> List[TableDiff]:
    """Разница двух срезов по ячейкам (a — было, b — стало; None — пусто). Только непустые таблицы-разницы."""
    ta = (a or {}).get("tables", {})
    tb = (b or {}).get("tables", {})
    out = []
    for name in list(dict.fromkeys(list(ta) + list(tb))):
        A = ta.get(name, {"header": [], "rows": []})
        B = tb.get(name, {"header": [], "rows": []})
        d = TableDiff(name)
        if A["header"] != B["header"] and A["header"] and B["header"]:
            d.header = (list(A["header"]), list(B["header"]))
        ka, oa = _keyed(A, name)
        kb, ob = _keyed(B, name)
        d.added = [k for k in ob if k not in ka]
        d.removed = [k for k in oa if k not in kb]
        for k in ob:
            if k in ka:
                ra, rb = ka[k], kb[k]
                for c in dict.fromkeys(list(B["header"]) + list(A["header"])):
                    va, vb = ra.get(c, ""), rb.get(c, "")
                    if va != vb:
                        d.cells.append((k, c, va, vb))
        common_a = [k for k in oa if k in kb]
        common_b = [k for k in ob if k in ka]
        d.reordered = common_a != common_b
        if not d.empty:
            out.append(d)
    return out


def format_diff(diffs: List[TableDiff], limit: int = 40) -> str:
    """Человекочитаемая разница: по таблицам, первые `limit` ячеек каждой."""
    if not diffs:
        return "без изменений"
    lines = []
    for d in diffs:
        parts = []
        if d.cells:
            parts.append(f"ячеек {len(d.cells)}")
        if d.added:
            parts.append(f"строк добавлено {len(d.added)}")
        if d.removed:
            parts.append(f"удалено {len(d.removed)}")
        if d.header:
            parts.append("заголовок изменён")
        if d.reordered:
            parts.append("порядок строк изменён")
        lines.append(f"{d.table}: " + ", ".join(parts))
        for k, c, va, vb in d.cells[:limit]:
            lines.append(f"    {k}  {c}: {va!r} → {vb!r}")
        if len(d.cells) > limit:
            lines.append(f"    … ещё {len(d.cells) - limit}")
        for k in d.added[:limit]:
            lines.append(f"    + {k}")
        for k in d.removed[:limit]:
            lines.append(f"    − {k}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------------------------- обмен с Excel
# Выгрузка для правки в Excel: числа с десятичной ЗАПЯТОЙ и «;» (русская локаль), UTF-8 с BOM (иначе Excel не узнает
# кодировку кириллицы). Рядом — `.base.json` (срез, с которого выгружено) и `export.json` (что и откуда).
# Загрузка сравнивает файлы со срезом-базой по ячейкам и даёт ПРАВКИ (было → стало); правки накладываются на текущую
# рабочую копию (независимые правки других не теряются, одна и та же ячейка, изменённая с двух сторон, — конфликт).
# Excel хранит числа с 15 значащими цифрами: ячейка, совпавшая с базой на 15 значащих цифрах (и «001» = «1»,
# «1,00E-05» = «1e-05»), правкой не считается — остаётся текст базы.
import re as _re

_NUM_DOT = _re.compile(r"^[+-]?(\d+\.\d*|\.\d+|\d+)([eE][+-]?\d+)?$")
_NUM_COMMA = _re.compile(r"^[+-]?(\d+,\d*|,\d+|\d+)([eE][+-]?\d+)?$")
EXPORT_FORMAT = "synchrotronsim.export/1"


def to_excel(v: str) -> str:
    return v.replace(".", ",") if "." in v and _NUM_DOT.match(v) else v


def from_excel(v: str) -> str:
    return v.replace(",", ".") if "," in v and _NUM_COMMA.match(v.strip()) else v


def _num(v: str):
    v = v.strip()
    if not (_NUM_DOT.match(v)):
        return None
    try:
        return float(v)
    except ValueError:
        return None


def same_in_excel(base: str, new: str) -> bool:
    """Ячейка не правилась: тот же текст или то же число на 15 значащих цифрах (точность Excel)."""
    if base == new:
        return True
    a, b = _num(base), _num(new)
    return a is not None and b is not None and f"{a:.15g}" == f"{b:.15g}"


def export_slice(sl: Dict[str, Any], folder: str, excel: bool = True, origin: str = "") -> List[str]:
    """Срез → папка CSV (по файлу на таблицу) + `.base.json` + `export.json`. → список файлов."""
    os.makedirs(folder, exist_ok=True)
    files = []
    for name, t in sl["tables"].items():
        rows = [[to_excel(c) for c in r] for r in t["rows"]] if excel else [list(r) for r in t["rows"]]
        p = os.path.join(folder, name + ".csv")
        write_sheet(p, Sheet(list(t["header"]), rows, bom=excel, newline="\r\n"))
        files.append(p)
    with open(os.path.join(folder, ".base.json"), "wb") as f:
        f.write(canonical(sl))
    with open(os.path.join(folder, "export.json"), "w", encoding="utf-8") as f:
        json.dump({"format": EXPORT_FORMAT, "kind": sl["kind"], "id": sl["id"], "rev": content_rev(sl),
                   "origin": origin, "excel": bool(excel), "tables": list(sl["tables"])}, f, ensure_ascii=False,
                  indent=1)
    return files


def read_export(folder: str) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """Папка выгрузки → (сведения export.json, срез-база, срез с правками из файлов)."""
    with open(os.path.join(folder, "export.json"), encoding="utf-8") as f:
        info = json.load(f)
    with open(os.path.join(folder, ".base.json"), encoding="utf-8") as f:
        base = json.load(f)
    excel = bool(info.get("excel"))
    edited = json.loads(json.dumps(base))
    for name, t in base["tables"].items():
        p = os.path.join(folder, name + ".csv")
        if not os.path.exists(p):
            continue                                              # файл убрали — таблица не правилась
        s = read_sheet(p)
        hdr = [h.strip() for h in s.header]
        if hdr != list(t["header"]):
            raise ValueError(f"{name}.csv: заголовок изменён ({hdr[:4]}… ≠ {t['header'][:4]}…) — столбцы в выгрузке "
                             f"не правятся")
        n = len(hdr)
        rows = [(r + [""] * (n - len(r)))[:n] for r in s.rows if _nonempty(r)]
        if excel:
            rows = [[from_excel(c) for c in r] for r in rows]
        kb, ob = _keyed(t, name)
        ke, oe = _keyed({"header": hdr, "rows": rows}, name)
        out_rows = []
        for k in oe:                                              # порядок — как в файле
            re_ = ke[k]
            if k in kb:
                rb = kb[k]
                out_rows.append([rb[c] if (excel and same_in_excel(rb[c], re_[c])) or rb[c] == re_[c] else re_[c]
                                 for c in hdr])
            else:
                out_rows.append([re_[c] for c in hdr])
        edited["tables"][name]["rows"] = out_rows
    return info, base, edited


@dataclass
class Conflict:
    table: str
    key: str
    column: str
    base: str
    ours: str      # в рабочей копии
    theirs: str    # в выгрузке


def apply_edits(current: Dict[str, Any], base: Dict[str, Any], edited: Dict[str, Any]):
    """Наложить правки «база → выгрузка» на текущий срез (рабочую копию). → (новый срез, конфликты).
    Ячейка: текущая = базе → берётся правка; текущая = правке → ничего; иначе — конфликт (текущая остаётся).
    Строка добавлена в выгрузке — добавляется (если ключа нет); удалена — удаляется, если в рабочей копии не правилась."""
    out = json.loads(json.dumps(current))
    conflicts: List[Conflict] = []
    for d in diff_slices(base, edited):
        t = out["tables"].get(d.table)
        if t is None:
            continue
        h = t["header"]
        ki = [h.index(c) for c in ROW_KEY.get(d.table, ()) if c in h]
        index = {"|".join(r[c] for c in ki): r for r in t["rows"] if _nonempty(r)}
        braw = {k: v for k, v in _keyed(base["tables"][d.table], d.table)[0].items()}
        eraw = _keyed(edited["tables"][d.table], d.table)[0]
        for key, col, va, vb in d.cells:
            r = index.get(key)
            if r is None or col not in h:
                conflicts.append(Conflict(d.table, key, col, va, "(строки нет)", vb))
                continue
            cur = r[h.index(col)]
            if cur == va:
                r[h.index(col)] = vb
            elif cur != vb:
                conflicts.append(Conflict(d.table, key, col, va, cur, vb))
        for key in d.added:
            if key not in index:
                t["rows"].append([eraw[key].get(c, "") for c in h])
        for key in d.removed:
            r = index.get(key)
            if r is not None and [r[h.index(c)] for c in h] == [braw[key].get(c, "") for c in h]:
                t["rows"].remove(r)
            elif r is not None:
                conflicts.append(Conflict(d.table, key, "(строка)", "удалена в выгрузке", "изменена в рабочей копии", ""))
    return out, conflicts


# ---------------------------------------------------------------------------------------------- машина из среза
def machine_from_slices(machine_sl: Dict[str, Any], library_sl: Dict[str, Any], source: str = ""):
    """Документ Machine из среза машины и справочника — тем же разбором, что и из CSV (Repository).
    source — метка происхождения (по умолчанию store:machine/<id>@<ревизия среза>+<ревизия справочника>)."""
    from .db import MemoryRepository
    label = source or f"store:machine/{machine_sl['id']}@{content_rev(machine_sl)[:16]}+{content_rev(library_sl)[:16]}"
    repo = MemoryRepository({**library_sl.get("tables", {}), **machine_sl.get("tables", {})},
                            passport={int(machine_sl["id"]): machine_sl.get("passport")}, label=label)
    return repo.machine(int(machine_sl["id"]))
