# -*- coding: utf-8 -*-
"""Единый слой доступа к БД SynchrotronSim (`workspace/tables`, `;`-CSV) — шаг 1 docs/LIBRARY_ARCHITECTURE.md.

Формат таблиц НЕ меняется: Repository читает те же файлы, что `matrix_tracker`, `beam_injection`, `fb_config`
и C++, и собирает из них объекты `synchrotronsim.model`. Отличия от прежних загрузчиков:
  * пути — от корня проекта, а не от текущей папки (скрипт работает из любой папки);
  * `overlay` — папка возмущённой машины (ансамбль ошибок `errors.py`): ЛЮБАЯ таблица, которая в ней есть,
    перекрывает основную (как у C++; прежний Python подменял только lattice_configs и связность);
  * чтение по ИМЕНИ столбца, заголовок с BOM допускается;
  * таблица кэшируется по (путь, время изменения, размер) — повторное чтение без разбора.

    db = Repository()                       # корень проекта по умолчанию
    m = db.machine(26)                      # Machine: строки решётки, магниты, устройства, экраны, …
    m.elements.kind, m.elements.mat_G       # столбцы numpy (только чтение)
    m.fingerprint                           # sha256 содержимого — ключ кэшей и паспорт прогона
"""
from __future__ import annotations

import csv
import math
import os
from typing import Dict, List, Optional, Sequence

import numpy as np

from .model import (BpmSpec, Cavity, DeviceType, ElementGroup, ElementTable, ErrorSigma, FieldMapRef,
                    InstanceTable, Machine, MachineMeta, PulsedKicker, Screen)

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


class Table:
    """Сырой `;`-CSV: заголовок + строки-строки. Доступ по имени столбца; отсутствующий столбец — пустые ячейки."""

    def __init__(self, path: str, header: List[str], rows: List[List[str]]):
        self.path, self.header, self.rows = path, header, rows
        self.index = {c: i for i, c in enumerate(header)}

    def __len__(self):
        return len(self.rows)

    def has(self, col: str) -> bool:
        return col in self.index

    def cell(self, row: Sequence[str], col: str) -> str:
        i = self.index.get(col)
        return row[i] if i is not None and i < len(row) else ""

    def column(self, col: str) -> List[str]:
        i = self.index.get(col)
        if i is None:
            return [""] * len(self.rows)
        return [r[i] if i < len(r) else "" for r in self.rows]

    def where(self, col: str, value) -> List[List[str]]:
        i = self.index.get(col)
        if i is None:
            return []
        v = str(value)
        return [r for r in self.rows if i < len(r) and r[i] == v]

    def dicts(self, rows=None) -> List[Dict[str, str]]:
        return [{c: (r[i] if i < len(r) else "") for c, i in self.index.items()} for r in (self.rows if rows is None else rows)]


_MISSING = ("", "-")                                   # «не задано» в таблицах (arg_val у элементов без карты — «-»)


def _f(s: str, default=math.nan) -> float:
    """Число из ячейки: пусто или «-» ⇒ default (NaN — «не задано»)."""
    s = (s or "").strip()
    return default if s in _MISSING else float(s)


def _i(s: str, default=0) -> int:
    s = (s or "").strip()
    return int(float(s)) if s else default


def _b(s: str) -> bool:
    return (s or "").strip().lower() in ("1", "true", "yes", "да")


def _vec(t: Table, r, cols) -> tuple:
    return tuple(_f(t.cell(r, c), 0.0) for c in cols)


class Repository:
    """Доступ к таблицам БД и сборка объектов модели.

    root    — корень проекта (с `workspace/tables`) или сама рабочая папка (с `tables/`); None — корень проекта.
    overlay — рабочая папка возмущённой машины (с `tables/`): её таблицы перекрывают основные.
    """
    _cache: Dict[str, tuple] = {}                      # путь → (mtime_ns, size, Table) — общий на процесс

    def __init__(self, root: Optional[str] = None, overlay: Optional[str] = None):
        root = os.path.abspath(root or PROJECT_ROOT)
        ws = root if os.path.isdir(os.path.join(root, "tables")) else os.path.join(root, "workspace")
        self.workspace = ws
        self.tables_dir = os.path.join(ws, "tables")
        if not os.path.isdir(self.tables_dir):
            raise FileNotFoundError(f"нет папки таблиц {self.tables_dir}")
        self.overlay_dir = None
        if overlay:
            ov = os.path.abspath(overlay)
            ov_tb = os.path.join(ov, "tables") if os.path.isdir(os.path.join(ov, "tables")) else ov
            if not os.path.isdir(ov_tb):
                raise FileNotFoundError(f"нет папки таблиц overlay {ov_tb}")
            self.overlay_dir = ov_tb

    # ------------------------------------------------------------------------------------------ таблицы
    def path(self, name: str) -> str:
        """Путь к таблице с учётом overlay (`name` — без или с `.csv`)."""
        fn = name if name.endswith(".csv") else name + ".csv"
        if self.overlay_dir and os.path.exists(os.path.join(self.overlay_dir, fn)):
            return os.path.join(self.overlay_dir, fn)
        return os.path.join(self.tables_dir, fn)

    def exists(self, name: str) -> bool:
        return os.path.exists(self.path(name))

    def table(self, name: str, required: bool = True) -> Table:
        p = self.path(name)
        if not os.path.exists(p):
            if required:
                raise FileNotFoundError(f"нет таблицы {p}")
            return Table(p, [], [])
        st = os.stat(p)
        hit = self._cache.get(p)
        if hit and hit[0] == st.st_mtime_ns and hit[1] == st.st_size:
            return hit[2]
        with open(p, encoding="utf-8-sig", newline="") as fh:
            rd = csv.reader(fh, delimiter=";")
            header = next(rd, [])
            rows = [r for r in rd if r and any(c.strip() for c in r)]
        t = Table(p, [h.strip() for h in header], rows)
        self._cache[p] = (st.st_mtime_ns, st.st_size, t)
        return t

    def lattices(self) -> List[int]:
        """Номера решёток в lattice_configs (по возрастанию)."""
        return sorted({int(v) for v in self.table("lattice_configs").column("lattice_id") if v.strip()})

    # ------------------------------------------------------------------------------------------ устройства
    def device_library(self) -> Dict[str, DeviceType]:
        hw = self.table("hardware_library")
        tr = self.table("element_transforms", required=False)
        fm = self.table("field_maps", required=False)
        trd = {tr.cell(r, "device_type_id"): r for r in tr.rows}
        maps: Dict[str, list] = {}
        for r in fm.rows:
            maps.setdefault(fm.cell(r, "device_type_id"), []).append(FieldMapRef(
                map_id=fm.cell(r, "map_id"), element=fm.cell(r, "element"), winding=fm.cell(r, "winding"),
                amp_turns=_f(fm.cell(r, "amp_turns")), current_a=fm.cell(r, "current_a"),
                excitation=fm.cell(r, "excitation"), split=fm.cell(r, "split"), grid=fm.cell(r, "grid"),
                material_curve=fm.cell(r, "material_curve"), raw_path=fm.cell(r, "raw_path"),
                bin_path=fm.cell(r, "bin_path"), source_fld=fm.cell(r, "source_fld"), note=fm.cell(r, "note")))
        out = {}
        for r in hw.rows:
            dev = hw.cell(r, "device_type_id")
            t = trd.get(dev)
            c = hw.cell
            out[dev] = DeviceType(
                device_type_id=dev, active=_b(c(r, "Active")), field_type=c(r, "field_type"),
                kind=tr.cell(t, "kind") if t else "", active_coeffs=tr.cell(t, "active_coeffs") if t else "",
                n_slices=_i(tr.cell(t, "n_slices")) if t else 0,
                eff_length=_f(c(r, "eff_length"), 0.0), nominal_arg=_f(c(r, "nominal_arg")),
                field_direction=_f(c(r, "field_direction"), 1.0), correct_coeff=_f(c(r, "correct_coeff"), 1.0),
                edge_ff=_f(c(r, "edge_ff"), 0.0), edge_delta=_f(c(r, "edge_delta")),
                use_symmetry=_b(c(r, "use_symmetry")), synthetic=_b(c(r, "synthetic")),
                scaled_from=c(r, "scaled_from"), map_raw_path=c(r, "map_raw_path"),
                map_bin_path=c(r, "map_bin_path"), stl_path=c(r, "stl_path"),
                turns_per_coil=_f(c(r, "turns_per_coil")), num_coils=_f(c(r, "num_coils")),
                total_turns=_f(c(r, "total_turns")), nominal_current=_f(c(r, "nominal_current_a")),
                rot_x=math.radians(_f(c(r, "rot_x_deg"), 0.0)), rot_y=math.radians(_f(c(r, "rot_y_deg"), 0.0)),
                rot_z=math.radians(_f(c(r, "rot_z_deg"), 0.0)),
                notes=tr.cell(t, "notes") if t else "", maps=tuple(maps.get(dev, ())))
        return out

    def error_model(self) -> Dict[str, ErrorSigma]:
        """Допуски по типу магнита (`element_errors`); field_pct (проценты) → доля."""
        t = self.table("element_errors", required=False)
        out = {}
        for r in t.rows:
            k = t.cell(r, "kind").strip()
            if k:
                out[k] = ErrorSigma(kind=k, **{c: _f(t.cell(r, c), 0.0) for c in
                                               ("dx", "dy", "dz", "dyaw", "dpitch", "droll")},
                                    field_rel=_f(t.cell(r, "field_pct"), 0.0) / 100.0)
        return out

    # ------------------------------------------------------------------------------------------ установка
    _MACHINE_TABLES = ("lattice_configs", "element_connectivity", "hardware_library", "element_transforms",
                       "field_maps", "element_groups", "rf_cavities", "pulsed_kickers", "bpm_config", "screens")
    _machines: Dict[tuple, tuple] = {}                 # (папки, решётка) → (подпись таблиц, Machine) — на процесс

    def _signature(self, lat: int) -> tuple:
        from .machines import LATTICE_INFO
        sig = []
        for name in self._MACHINE_TABLES:
            p = self.path(name)
            try:
                st = os.stat(p); sig.append((p, st.st_mtime_ns, st.st_size))
            except FileNotFoundError:
                sig.append((p, None, None))
        return tuple(sig) + (repr(LATTICE_INFO.get(lat)),)

    def machine(self, lattice_id: int) -> Machine:
        """Установка по номеру решётки. Machine неизменяема, поэтому кэшируется на процесс: повторный вызов
        отдаёт тот же объект, пока не изменились файлы таблиц (время изменения, размер)."""
        lat = int(lattice_id)
        key = (self.tables_dir, self.overlay_dir, lat)
        sig = self._signature(lat)
        hit = self._machines.get(key)
        if hit and hit[0] == sig:
            return hit[1]
        m = self._build_machine(lat)
        self._machines[key] = (sig, m)
        return m

    def _build_machine(self, lattice_id: int) -> Machine:
        lat = str(int(lattice_id))
        lc = self.table("lattice_configs")
        cn = self.table("element_connectivity")
        conf = {lc.cell(r, "instance_id"): r for r in lc.where("lattice_id", lat)}
        if not conf:
            raise KeyError(f"решётки {lat} нет в {lc.path}")
        order = [cn.cell(r, "instance_id") for r in cn.where("lattice_id", lat)]
        order.sort(key=int)                                       # порядок обхода — instance_id как целое
        missing = [i for i in order if i not in conf]
        if missing:
            raise ValueError(f"решётка {lat}: в element_connectivity есть instance_id без строки lattice_configs: "
                             f"{missing[:5]}{'…' if len(missing) > 5 else ''}")
        rows = [conf[i] for i in order]
        inst = self._instances(lc, rows)
        devices_all = self.device_library()
        elements = self._elements(lc, rows, inst, devices_all)
        used = {d: devices_all[d] for d in dict.fromkeys(inst.device_type.tolist()) if d in devices_all}
        return Machine(
            lattice_id=int(lat), instances=inst, elements=elements, devices=used,
            groups=self._groups(lat), cavities=self._cavities(lat), kickers=self._kickers(lat),
            bpms=self._bpms(lat), screens=self._screens(lat), meta=self._meta(int(lat)),
            source=self.overlay_dir or self.tables_dir)

    @staticmethod
    def _instances(lc: Table, rows) -> InstanceTable:
        c = lc.cell
        S = lambda col: np.array([c(r, col) for r in rows], dtype=object)
        F = lambda col: np.array([_f(c(r, col)) for r in rows], dtype=float)
        return InstanceTable(
            instance_id=S("instance_id"), device_type=S("device_type_id"), kind=S("kind"),
            x=F("x"), y=F("y"), z=F("z"), yaw=F("yaw"), pitch=F("pitch"), roll=F("roll"),
            arg_val=F("arg_val"), correct_coeff=F("correct_coeff"),
            mat_B=F("mat_B"), mat_G=F("mat_G"), mat_Gv=F("mat_Gv"), mat_S=F("mat_S"), mat_O=F("mat_O"),
            eff_length=F("eff_length"), n_slices=np.array([_i(c(r, "n_slices")) for r in rows], dtype=np.int64),
            edge1=F("edge1"), edge2=F("edge2"), lf1=F("lf1"), lf2=F("lf2"))

    @staticmethod
    def _elements(lc: Table, rows, inst: InstanceTable, devices: Dict[str, DeviceType]) -> ElementTable:
        """Склейка обмоток: подряд идущие строки с ОДИНАКОВЫМИ СТРОКАМИ x, y, z (как build_line и C++).
        Числа — теми же операциями и в том же порядке, что matrix_tracker.Magnet (бит-в-бит)."""
        c = lc.cell
        cols = {k: [] for k in ElementTable._COLUMNS}
        n = len(rows); i = 0
        while i < n:
            r0 = rows[i]; key = (c(r0, "x"), c(r0, "y"), c(r0, "z")); j = i + 1
            while j < n and (c(rows[j], "x"), c(rows[j], "y"), c(rows[j], "z")) == key:
                j += 1
            dev = c(r0, "device_type_id"); dt = devices.get(dev)
            B = G = Gv = S = O = 0.0
            for r in rows[i:j]:                                   # суммы — как Magnet.__init__
                B += float(c(r, "mat_B") or 0); G += float(c(r, "mat_G") or 0)
                Gv += float(c(r, "mat_Gv") or c(r, "mat_G") or 0)
                S += float(c(r, "mat_S") or 0); O += float(c(r, "mat_O") or 0)
            cols["instance_id"].append(c(r0, "instance_id")); cols["device_type"].append(dev)
            cols["kind"].append(c(r0, "kind") or (dt.kind if dt else "") or "")
            for k in ("x", "y", "z"):
                cols[k].append(float(c(r0, k)))
            cols["yaw"].append(float(c(r0, "yaw")))
            cols["pitch"].append(float(c(r0, "pitch") or 0.0)); cols["roll"].append(float(c(r0, "roll") or 0.0))
            cols["eff_length"].append(float(c(r0, "eff_length") or (dt.eff_length if dt else 0.0) or 0.0))
            cols["n_slices"].append(_i(c(r0, "n_slices")) or (dt.n_slices if dt else 0) or 0)
            cols["mat_B"].append(B); cols["mat_G"].append(G); cols["mat_Gv"].append(Gv)
            cols["mat_S"].append(S); cols["mat_O"].append(O)
            for k in ("edge1", "edge2", "lf1", "lf2"):
                cols[k].append(float(c(r0, k) or 0.0))
            cols["member_start"].append(i); cols["member_count"].append(j - i)
            i = j
        obj = ("instance_id", "device_type", "kind")
        ints = ("n_slices", "member_start", "member_count")
        return ElementTable(**{k: np.array(v, dtype=object if k in obj else np.int64 if k in ints else float)
                               for k, v in cols.items()})

    def _groups(self, lat: str):
        t = self.table("element_groups", required=False)
        return tuple(ElementGroup(group_id=t.cell(r, "group_id"), kind=t.cell(r, "kind"),
                                  pivot=_vec(t, r, ("pivot_x", "pivot_y", "pivot_z")),
                                  yaw=_f(t.cell(r, "yaw"), 0.0), pitch=_f(t.cell(r, "pitch"), 0.0),
                                  roll=_f(t.cell(r, "roll"), 0.0),
                                  instance_ids=tuple(s for s in t.cell(r, "instance_ids").split("|") if s))
                     for r in t.where("lattice_id", lat))

    def _cavities(self, lat: str):
        """Пустая доля напряжения ⇒ 1/N (как C++ readCavities и beam_injection.load_cavities)."""
        t = self.table("rf_cavities", required=False)
        rows = t.where("lattice_id", lat)
        return tuple(Cavity(cavity_id=t.cell(r, "cavity_id"), point=_vec(t, r, ("x", "y", "z")),
                            normal=_vec(t, r, ("nx", "ny", "nz")),
                            volt_share=_f(t.cell(r, "volt_share"), 1.0 / len(rows)),
                            note=t.cell(r, "note")) for r in rows)

    def _kickers(self, lat: str):
        t = self.table("pulsed_kickers", required=False)
        return tuple(PulsedKicker(kicker_id=t.cell(r, "kicker_id"), point=_vec(t, r, ("x", "y", "z")),
                                  normal=_vec(t, r, ("nx", "ny", "nz")), shape=t.cell(r, "shape"),
                                  duration=_f(t.cell(r, "duration_s")), rise=_f(t.cell(r, "rise_s")),
                                  note=t.cell(r, "note")) for r in t.where("lattice_id", lat))

    def _bpms(self, lat: str):
        t = self.table("bpm_config", required=False)
        return tuple(BpmSpec(sensor_id=t.cell(r, "id_sensor"),
                             **{k: _f(t.cell(r, k), d) for k, d in (("dx", 0.0), ("dy", 0.0), ("tilt", 0.0),
                                                                    ("gain_x", 1.0), ("gain_y", 1.0),
                                                                    ("setpoint_x", 0.0), ("setpoint_y", 0.0))})
                     for r in t.where("lattice_id", lat))

    def _screens(self, lat: str):
        t = self.table("screens", required=False)
        mm, mrad = 1e-3, 1e-3
        out = []
        for r in t.where("lattice_id", lat):
            g = lambda col, d=0.0: _f(t.cell(r, col), d)
            out.append(Screen(
                screen_id=t.cell(r, "screen_id"), source_ref=t.cell(r, "source_ref"),
                source_kind=t.cell(r, "source_kind"), center=_vec(t, r, ("cx", "cy", "cz")),
                normal=_vec(t, r, ("nx", "ny", "nz")), u=_vec(t, r, ("ux", "uy", "uz")),
                v=_vec(t, r, ("vx", "vy", "vz")), width=g("width_mm") * mm, height=g("height_mm") * mm,
                mode=t.cell(r, "mode"), dist=g("dist_m"), extract=g("extract_mrad") * mrad,
                los_h=g("los_h_mrad") * mrad, los_v=g("los_v_mrad") * mrad, tilt_h=g("tilt_h_mrad") * mrad,
                tilt_v=g("tilt_v_mrad") * mrad, roll=g("roll_mrad") * mrad))
        return tuple(out)

    @staticmethod
    def _meta(lat: int) -> Optional[MachineMeta]:
        from .machines import LATTICE_INFO
        inf = LATTICE_INFO.get(lat)
        if not inf:
            return None
        return _meta_from_info(inf)


def _meta_from_info(inf) -> Optional[MachineMeta]:
    if not inf:
        return None
    volt, harm = inf["rf"]
    return MachineMeta(label=inf["label"], design_energy_mev=float(inf["energy_mev"]),
                       rf_volt_mv=float(volt), rf_harmonic=int(harm), superperiods=int(inf["superperiods"]))


class MemoryRepository(Repository):
    """Repository над таблицами В ПАМЯТИ (срезы ревизий хранилища, synchrotronsim.tableset): тот же разбор, что из
    CSV, без файлов. tables — {имя: {"header": [...], "rows": [[...]]}}; passport — {lattice_id: паспорт} (как в
    LATTICE_INFO; паспорт хранится в ревизии — не зависит от кода machines.py); label — уникальная метка
    содержимого (ключ кэша машин, источник Machine.source)."""

    def __init__(self, tables: Dict[str, dict], passport: Optional[Dict[int, dict]] = None, label: str = "memory"):
        self.workspace = None
        self.tables_dir = label
        self.overlay_dir = None
        self._mem = {}
        for name, t in tables.items():
            rows = [list(r) for r in t.get("rows", []) if r and any(c.strip() for c in r)]
            self._mem[name] = Table(f"{label}/{name}.csv", [h.strip() for h in t.get("header", [])], rows)
        self._passport = dict(passport or {})

    def path(self, name: str) -> str:
        return f"{self.tables_dir}/{name if name.endswith('.csv') else name + '.csv'}"

    def exists(self, name: str) -> bool:
        return (name[:-4] if name.endswith(".csv") else name) in self._mem

    def table(self, name: str, required: bool = True) -> Table:
        key = name[:-4] if name.endswith(".csv") else name
        if key in self._mem:
            return self._mem[key]
        if required:
            raise FileNotFoundError(f"нет таблицы {key} в {self.tables_dir}")
        return Table(self.path(name), [], [])

    def _signature(self, lat: int) -> tuple:
        return ("memory", self.tables_dir, repr(self._passport.get(lat)))

    def _meta(self, lat: int) -> Optional[MachineMeta]:
        return _meta_from_info(self._passport.get(int(lat)))


def open_db(root: Optional[str] = None, overlay: Optional[str] = None) -> Repository:
    """Короткая форма `Repository(root, overlay)`."""
    return Repository(root, overlay)
