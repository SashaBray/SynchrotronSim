# -*- coding: utf-8 -*-
"""Журнал правок проекта (docs/DATA_PROTOCOL.md, этап 4, раздел 8).

Правка = изменение ячеек СРЕЗА машины (synchrotronsim.tableset): таблица, строки (instance_id), столбец, было →
стало (текст ячеек). Адрес правки — `design.<instance_id первой обмотки>.<столбец>`: параметр ПРОЕКТА (проектная
сила, положение, длина — сохраняется в машину), в отличие от живого состояния (ток, зазор, усиление контура).

Цена правки:
  * hot  — сила существующего элемента (mat_G, mat_Gv, mat_S, mat_O): модель применяет её сразу ручкой ядра
           (Δk = QP·Δполе: DK/DKV/DK2/DK3 — базовым возмущением, переживает шаги и рестарт);
  * cold — всё остальное (поле диполя, положение и углы, длина, срезы, грани): применяется пересборкой модели
           с изменённой машиной (кэш откликов — по отпечатку, пересчитается сам).

Журнал (`Journal`) живёт у клиента: срез-база + правки → текущий срез → документ `Machine` (тем же разбором, что из
CSV). Отмена/повтор; «сохранить» — правки накладываются на рабочую копию таблиц и фиксируются ревизией хранилища.
Модель получает каждую горячую правку и пересчитывает свою машину (`apply_machine`) — паспорт прогона пишет
машину, которая реально считается, и список правок.
"""
from __future__ import annotations

import copy
import math
from dataclasses import dataclass, field, replace
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# столбец lattice_configs → (какие строки элемента: first — первая обмотка, all — все; цена)
EDITABLE: Dict[str, Tuple[str, str]] = {
    "mat_G": ("first", "hot"), "mat_Gv": ("first", "hot"), "mat_S": ("first", "hot"), "mat_O": ("first", "hot"),
    "mat_B": ("first", "cold"), "eff_length": ("first", "cold"), "n_slices": ("first", "cold"),
    "edge1": ("first", "cold"), "edge2": ("first", "cold"),
    "x": ("all", "cold"), "y": ("all", "cold"), "z": ("all", "cold"),          # обмотки склеены по x,y,z — все
    "yaw": ("all", "cold"), "pitch": ("all", "cold"), "roll": ("all", "cold"),
}
UNITS = {"mat_B": "Тл", "mat_G": "Тл/м", "mat_Gv": "Тл/м", "mat_S": "Тл/м²", "mat_O": "Тл/м³", "eff_length": "м",
         "n_slices": "", "edge1": "рад", "edge2": "рад", "x": "м", "y": "м", "z": "м", "yaw": "рад", "pitch": "рад",
         "roll": "рад"}
TITLES = {"mat_B": "поле B", "mat_G": "градиент G", "mat_Gv": "градиент G верт.", "mat_S": "секступоль S",
          "mat_O": "октуполь O", "eff_length": "длина L", "n_slices": "срезов", "edge1": "грань входа",
          "edge2": "грань выхода", "x": "x", "y": "y", "z": "z", "yaw": "рыскание", "pitch": "тангаж", "roll": "крен"}
# что показывать в окне элемента (по виду): сила вида + геометрия
FIELDS_BY_KIND = {"quad": ("mat_G", "mat_Gv"), "sext": ("mat_S",), "oct": ("mat_O",), "dipole": ("mat_B", "mat_G"),
                  "dq": ("mat_B", "mat_G", "mat_S"), "corrector": (), "kicker": (), "monitor": ()}
GEOMETRY = ("x", "y", "z", "yaw", "pitch", "roll", "eff_length", "n_slices")


class EditError(ValueError):
    """Правка не применяется: столбец не правится, ячейка уже не та (конфликт), структура машины изменилась."""


@dataclass
class Edit:
    """Правка ячеек: строки `keys` таблицы `table`, столбец `column`, было `old` → стало `new` (текст ячеек)."""
    column: str
    keys: Tuple[str, ...]
    old: Tuple[str, ...]
    new: Tuple[str, ...]
    element: int = -1                  # индекс магнита (порядок пучка) — для ручек ядра
    name: str = ""                     # имя элемента для людей
    table: str = "lattice_configs"
    op: str = "set"

    @property
    def addr(self) -> str:
        return f"design.{self.keys[0]}.{self.column}"

    @property
    def cost(self) -> str:
        return EDITABLE.get(self.column, ("first", "cold"))[1]

    def inverse(self) -> "Edit":
        return replace(self, old=self.new, new=self.old)

    def to_dict(self) -> Dict[str, Any]:
        return {"op": self.op, "table": self.table, "column": self.column, "keys": list(self.keys),
                "old": list(self.old), "new": list(self.new), "element": int(self.element), "name": self.name}

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Edit":
        return cls(column=d["column"], keys=tuple(d["keys"]), old=tuple(d["old"]), new=tuple(d["new"]),
                   element=int(d.get("element", -1)), name=d.get("name", ""), table=d.get("table", "lattice_configs"),
                   op=d.get("op", "set"))

    def describe(self) -> str:
        u = UNITS.get(self.column, "")
        return (f"{self.name or self.keys[0]} · {TITLES.get(self.column, self.column)}: {self.old[0] or '—'} → "
                f"{self.new[0] or '—'}{(' ' + u) if u else ''}" + (f" (обмоток/строк {len(self.keys)})" if len(self.keys) > 1 else ""))


def _text(value) -> str:
    """Значение → текст ячейки: строка как есть (обрезка пробелов), число — repr (точно возвращается при чтении)."""
    if isinstance(value, str):
        return value.strip()
    if isinstance(value, (int, np.integer)) and not isinstance(value, bool):
        return str(int(value))
    x = float(value)
    if not math.isfinite(x):
        raise EditError(f"значение {value!r}: нужно конечное число")
    return repr(x)


def element_rows(machine, j: int) -> List[str]:
    """instance_id строк (обмоток) магнита j."""
    E = machine.elements
    s, n = int(E.member_start[j]), int(E.member_count[j])
    return [str(v) for v in machine.instances.instance_id[s:s + n]]


def _lc(sl) -> Dict[str, Any]:
    t = sl["tables"].get("lattice_configs")
    if t is None:
        raise EditError("в срезе машины нет lattice_configs")
    return t


def slice_cells(sl, keys, column) -> List[str]:
    t = _lc(sl)
    h = t["header"]
    if column not in h:
        raise EditError(f"в lattice_configs нет столбца {column}")
    ki, ci = h.index("instance_id"), h.index(column)
    rows = {r[ki]: r for r in t["rows"] if len(r) > ki}
    out = []
    for k in keys:
        r = rows.get(str(k))
        if r is None:
            raise EditError(f"нет строки instance_id={k}")
        out.append(r[ci] if ci < len(r) else "")
    return out


def make_edit(machine, sl, j: int, column: str, value, name: str = "") -> Edit:
    """Правка столбца `column` магнита j на `value` (число или текст ячейки): строки — по EDITABLE."""
    if column not in EDITABLE:
        raise EditError(f"столбец {column} не правится из клиента (правятся: {', '.join(EDITABLE)})")
    rows = element_rows(machine, j)
    keys = tuple(rows[:1] if EDITABLE[column][0] == "first" else rows)
    old = tuple(slice_cells(sl, keys, column))
    txt = _text(value)
    if column == "n_slices" and txt:
        try:
            if int(float(txt)) != float(txt) or int(float(txt)) < 0:
                raise ValueError
            txt = str(int(float(txt)))
        except ValueError:
            raise EditError(f"срезов — целое ≥ 0, не {value!r}") from None
    elif txt and column != "n_slices":
        try:
            float(txt)
        except ValueError:
            raise EditError(f"{column}: нужно число, не {value!r}") from None
    if column in ("x", "y", "z", "yaw") and not txt:
        raise EditError(f"{column} не может быть пустым")
    return Edit(column=column, keys=keys, old=old, new=tuple(txt for _ in keys), element=int(j), name=name)


def apply_slice(sl, edit: Edit, check: bool = True):
    """Срез с правкой (копия). check — ячейки должны быть равны `old` (иначе EditError: конфликт)."""
    out = copy.deepcopy(sl)
    t = out["tables"][edit.table]
    h = t["header"]
    ki, ci = h.index("instance_id"), h.index(edit.column)
    idx = {r[ki]: n for n, r in enumerate(t["rows"]) if len(r) > ki}
    for k, o, nv in zip(edit.keys, edit.old, edit.new):
        n = idx.get(str(k))
        if n is None:
            raise EditError(f"нет строки instance_id={k}")
        r = t["rows"][n]
        while len(r) <= ci:
            r.append("")
        if check and r[ci] != o:
            raise EditError(f"{edit.table} {k} {edit.column}: в ячейке {r[ci]!r}, а правка сделана от {o!r}")
        r[ci] = nv
    return out


# ---------------------------------------------------------------------------------------------- Machine
def _f(s: str, default=math.nan) -> float:
    s = (s or "").strip()
    return default if s in ("", "-") else float(s)


def _n0(x: float) -> float:
    return 0.0 if x != x else float(x)


def apply_machine(m, edit: Edit):
    """Документ Machine с правкой — теми же операциями, что разбор таблиц (synchrotronsim.db: строки → магниты),
    без таблиц: модель на вычислителе применяет горячие правки к своей машине. Структура (склейка обмоток) не
    меняется — правки геометрии переносят все обмотки элемента вместе."""
    from .model import Machine, InstanceTable, ElementTable
    I, E = m.instances, m.elements
    pos = {str(v): n for n, v in enumerate(I.instance_id)}
    cols = {c: np.array(getattr(I, c), copy=True) for c in I.column_names()}
    rows = []
    for k, nv in zip(edit.keys, edit.new):
        r = pos.get(str(k))
        if r is None:
            raise EditError(f"нет строки instance_id={k} в машине")
        rows.append(r)
        if edit.column == "n_slices":
            cols["n_slices"][r] = int(float(nv)) if nv.strip() else 0
        else:
            cols[edit.column][r] = _f(nv)
    I2 = InstanceTable(**cols)
    j = int(edit.element)
    if j < 0:                                                    # элемент по первой строке правки
        j = int(np.searchsorted(E.member_start, rows[0], side="right") - 1)
    s, n = int(E.member_start[j]), int(E.member_count[j])
    if not all(s <= r < s + n for r in rows):
        raise EditError("строки правки не из одного магнита")
    ecol = {c: np.array(getattr(E, c), copy=True) for c in E.column_names()}
    B = G = Gv = S = O = 0.0
    for r in range(s, s + n):                                    # суммы — как Repository._elements
        B += _n0(I2.mat_B[r]); G += _n0(I2.mat_G[r])
        gv = I2.mat_Gv[r]
        Gv += gv if gv == gv else _n0(I2.mat_G[r])
        S += _n0(I2.mat_S[r]); O += _n0(I2.mat_O[r])
    dev = str(E.device_type[j]); dt = m.devices.get(dev)
    ecol["mat_B"][j], ecol["mat_G"][j], ecol["mat_Gv"][j], ecol["mat_S"][j], ecol["mat_O"][j] = B, G, Gv, S, O
    r0 = s
    for k in ("x", "y", "z", "yaw"):
        ecol[k][j] = float(I2.__getattribute__(k)[r0])
    ecol["pitch"][j] = _n0(I2.pitch[r0]); ecol["roll"][j] = _n0(I2.roll[r0])
    el = I2.eff_length[r0]                                       # пусто ⇒ длина устройства; «0» — это ноль
    ecol["eff_length"][j] = float(el) if el == el else float((dt.eff_length if dt else 0.0) or 0.0)
    ecol["n_slices"][j] = int(I2.n_slices[r0]) or (dt.n_slices if dt else 0) or 0
    for k in ("edge1", "edge2", "lf1", "lf2"):
        ecol[k][j] = _n0(getattr(I2, k)[r0])
    E2 = ElementTable(**ecol)
    return Machine(lattice_id=m.lattice_id, instances=I2, elements=E2, devices=m.devices, groups=m.groups,
                   cavities=m.cavities, kickers=m.kickers, bpms=m.bpms, screens=m.screens, meta=m.meta,
                   source=f"{m.source}+правки" if "+правки" not in str(m.source) else m.source)


def strengths(m, j: int) -> List[float]:
    """Силы магнита j так, как их считает ядро при сборке: k = QP·поле (k1, k1v, k2, k3; QP побитово как в C++)."""
    import matrix_tracker as mt
    qp = float(mt.QP)
    E = m.elements
    j = int(j)
    return [qp * float(E.mat_G[j]), qp * float(E.mat_Gv[j]), qp * float(E.mat_S[j]), qp * float(E.mat_O[j])]


def hot_commands(m_old, m_new, j: int, cur: Optional[List[float]] = None):
    """Горячая правка → ручки ядра [(cmd, магнит, Δ)] и силы ядра после них. Ядро ПРИБАВЛЯЕТ Δ к текущей силе
    (DK: k1 и k1v; DKV: k1v; DK2: k2; DK3: k3), поэтому Δ = цель − ТЕКУЩЕЕ (cur — силы ядра сейчас; None — как
    собрано из m_old), а текущее обновляется теми же операциями, что в ядре. Цель — QP·поле новой машины, ровно
    как ядро посчитало бы при пересборке: при изменении силы меньше чем вдвое текущее + Δ = цель точно
    (лемма Стербенца) — кольцо с правкой побитово как пересобранное."""
    j = int(j)
    k1, k1v, k2, k3 = list(cur) if cur is not None else strengths(m_old, j)
    t1, t1v, t2, t3 = strengths(m_new, j)
    out = []
    if t1 != k1:
        d = t1 - k1
        out.append(("DK", j, d))                                  # DK — обе плоскости
        k1 = k1 + d
        k1v = k1v + d
    if t1v != k1v:
        d = t1v - k1v
        out.append(("DKV", j, d))                                 # остаток вертикали
        k1v = k1v + d
    if t2 != k2:
        d = t2 - k2
        out.append(("DK2", j, d))
        k2 = k2 + d
    if t3 != k3:
        d = t3 - k3
        out.append(("DK3", j, d))
        k3 = k3 + d
    return out, [k1, k1v, k2, k3]


# ---------------------------------------------------------------------------------------------- журнал
@dataclass
class Entry:
    edit: Edit
    applied: bool                        # уже в модели (горячая — сразу; холодная — после пересборки)


class Journal:
    """Правки проекта поверх среза-базы (клиент). Текущий срез = база + правки; машина — из него (тем же разбором,
    что из CSV)."""

    def __init__(self, base_slice: Dict[str, Any], library_slice: Dict[str, Any], machine_id: str,
                 base_rev: str = "", source: str = ""):
        from . import tableset as T
        self.base = copy.deepcopy(base_slice)
        self.lib = library_slice
        self.machine_id = str(machine_id)
        self.base_rev = base_rev or T.content_rev(self.base)
        self.source = source
        self.entries: List[Entry] = []
        self.redo_stack: List[Edit] = []
        self._behind = False            # модель отстаёт не только на неприменённые холодные правки (отменена применённая)
        self._cur = copy.deepcopy(self.base)
        self._m = None

    @property
    def rebuild_needed(self) -> bool:
        """Модель отличается от журнала тем, что применяется только пересборкой."""
        return self._behind or any(e.edit.cost != "hot" and not e.applied for e in self.entries)

    # ---- состояние
    @property
    def current(self) -> Dict[str, Any]:
        return self._cur

    def machine(self):
        if self._m is None:
            from . import tableset as T
            self._m = T.machine_from_slices(self._cur, self.lib,
                                            source=f"journal:{self.machine_id}@{self.base_rev[:12]}+{len(self.entries)}")
        return self._m

    def base_machine(self):
        from . import tableset as T
        return T.machine_from_slices(self.base, self.lib, source=f"journal:{self.machine_id}@{self.base_rev[:12]}")

    @property
    def edits(self) -> List[Edit]:
        return [e.edit for e in self.entries]

    def cell(self, keys, column) -> str:
        return slice_cells(self._cur, keys, column)[0]

    def _set(self, sl):
        self._cur = sl
        self._m = None

    # ---- правки
    def make(self, j: int, column: str, value, name: str = "") -> Edit:
        return make_edit(self.machine(), self._cur, j, column, value, name)

    def push(self, edit: Edit) -> Edit:
        """Применить правку к журналу. Горячая — считается применённой в модели (её шлёт кокпит), холодная —
        ждёт пересборки. Структура машины (число магнитов) меняться не должна."""
        n0 = len(self.machine().elements)
        new = apply_slice(self._cur, edit)
        old_cur, old_m = self._cur, self._m
        self._set(new)
        if len(self.machine().elements) != n0:
            self._cur, self._m = old_cur, old_m
            raise EditError("правка изменила склейку обмоток (число магнитов) — так править нельзя; это структурная "
                            "правка (редактор решётки)")
        self.entries.append(Entry(edit, edit.cost == "hot"))
        self.redo_stack.clear()
        return edit

    def undo(self) -> Optional[Tuple[Edit, bool]]:
        """Отменить последнюю правку → (обратная правка, слать ли её модели ручкой) или None."""
        if not self.entries:
            return None
        e = self.entries.pop()
        self._set(apply_slice(self._cur, e.edit.inverse()))
        self.redo_stack.append(e.edit)
        if e.edit.cost != "hot" and e.applied:
            self._behind = True                                   # модель уже с ней — снимет только пересборка
        return e.edit.inverse(), (e.edit.cost == "hot")

    def redo(self) -> Optional[Edit]:
        if not self.redo_stack:
            return None
        return self.push_keep_redo(self.redo_stack.pop())

    def push_keep_redo(self, edit: Edit) -> Edit:
        stack = list(self.redo_stack)
        self.push(edit)
        self.redo_stack = stack
        return edit

    def mark_rebuilt(self):
        """Модель пересобрана по machine(): все правки в ней."""
        for e in self.entries:
            e.applied = True
        self._behind = False

    def to_dict(self) -> Dict[str, Any]:
        return {"machine": self.machine_id, "base": self.base_rev, "edits": [e.edit.to_dict() for e in self.entries]}

    # ---- сохранение
    def save(self, store, note: str = ""):
        """Правки → рабочая копия таблиц (наложением: независимые правки рабочей копии сохраняются, одна ячейка,
        изменённая и там, и в журнале, — конфликт) → ревизия хранилища. Журнал переходит на сохранённое."""
        from . import tableset as T
        if not getattr(store, "working", None):
            raise EditError("у хранилища нет рабочей копии — правки сохраняются через рабочую копию таблиц")
        if not self.entries:
            return []
        cur_w = store.working_slices([int(self.machine_id)])[1][self.machine_id]
        new, conflicts = T.apply_edits(cur_w, self.base, self._cur)
        if conflicts:
            c = conflicts[0]
            raise EditError(f"рабочая копия изменилась: {c.table} {c.key} {c.column} — там {c.ours!r}, в журнале "
                            f"{c.theirs!r} (база {c.base!r}); конфликтов {len(conflicts)}")
        T.write_machine(store.working, new)
        refs = store.commit_working(note=note or f"правки из кокпита: {len(self.entries)}",
                                    lattices=[int(self.machine_id)])
        self._behind = self.rebuild_needed                       # неприменённое в модели остаётся неприменённым
        self.base = copy.deepcopy(self._cur)
        self.base_rev = T.content_rev(store.working_slices([int(self.machine_id)])[1][self.machine_id])
        self.entries = []
        self.redo_stack = []
        return refs
