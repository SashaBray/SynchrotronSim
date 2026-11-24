# -*- coding: utf-8 -*-
"""Прогон (`Run`) — результат с паспортом происхождения, и подписка на запись (`Watch`).
Слой 4 модели данных (docs/LIBRARY_ARCHITECTURE.md, разделы 3.6–3.7; шаг 4).

Папка прогона `workspace/runs/<run_id>/`:
    manifest.json   паспорт: вид прогона, ПОЛНАЯ копия конфигурации, отпечаток машины, seed, git sha (+ признак
                    незакоммиченных правок), сборка движка (версия компилятора, хэши модуля и сервера), версии
                    Python/numpy, время начала/конца, статус (in_progress | done | stopped | failed), состав рядов
    signals.npz     временные ряды по адресам: для ряда i — s<i>_step, s<i>_t, s<i>_value; имена — `names`
    tables/signals.csv  те же скалярные ряды для человека (десятичная ЗАПЯТАЯ, «;» — открывается в Excel)
    artifacts/      тяжёлые массивы (.npz) и картинки
    log.txt         журнал прогона
Индекс прогонов — `workspace/tables/runs.csv` (машинная таблица: десятичная точка).

    run = Run.create("fb_mode4", config=cfg_dict, machine=m, seed=1, source="fb_exp_id=1")
    run.record("obs.ebpm.x:rms", step, value, t=model_time)
    run.artifact("screen", {"I": arr}, meta="…")
    run.finish("done")
    Run.open(run.run_id).series("obs.ebpm.x:rms")  →  (step, t, value)

`Watch` — что записывать: адрес (единое адресное пространство, `synchrotronsim.addresses`), каденс и хранение.
Один и тот же объект создают блокнот (`exp.record.watch(...)`) и галочка «записывать» в клиенте (шаг 5).
"""
from __future__ import annotations

import csv
import datetime as _dt
import hashlib
import json
import os
import platform
import secrets
import subprocess
import sys
from dataclasses import asdict, dataclass, is_dataclass
from typing import Any, Dict, List, Optional

import numpy as np

from . import addresses as AD

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
INDEX_COLS = ["run_id", "kind", "source", "status", "started", "finished", "path", "note"]
STATUSES = ("in_progress", "done", "stopped", "failed")
CADENCES = ("per_step", "per_window", "final")


# ---------------------------------------------------------------------------------------------- Watch
@dataclass(frozen=True)
class Watch:
    """Подписка на запись одной величины.

    addr    — адрес (каноническая форма; псевдонимы приводятся) или имя рекордера (`recorder` ≠ "");
    cadence — per_step (каждый такт) | per_window (каждое окно усреднения) | final (один раз в конце);
    every   — прореживание: записывать каждый every-й раз;
    store   — series (временной ряд в signals.npz) | artifact (массив в artifacts/).
    recorder/params — рекордер-обработчик (control/recorders.py) для производных и дорогих величин."""
    addr: str
    cadence: str = "per_window"
    every: int = 1
    store: str = "series"
    recorder: str = ""
    params: str = ""

    def __post_init__(self):
        if self.cadence not in CADENCES and not self.recorder:   # рекордер сам ведёт свой каденс
            raise ValueError(f"Watch: каденс {self.cadence!r}, допустимы {CADENCES}")
        if int(self.every) < 1:
            raise ValueError("Watch: every ≥ 1")
        if not self.recorder:
            if not AD.is_address(self.addr):
                raise ValueError(f"Watch: {self.addr!r} — не адрес (области: {sorted(AD.AREAS)})")
            object.__setattr__(self, "addr", AD.canonical(self.addr))

    @property
    def meta(self) -> AD.AddrMeta:
        return AD.meta(self.addr)


# ---------------------------------------------------------------------------------------------- паспорт
def _sha16(path: str) -> Optional[str]:
    try:
        h = hashlib.sha256()
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        return h.hexdigest()[:16]
    except OSError:
        return None


def _git(root: str) -> Dict[str, Any]:
    def run(*args):
        try:
            return subprocess.run(["git", *args], cwd=root, capture_output=True, text=True, timeout=10,
                                  encoding="utf-8", errors="replace").stdout.strip()
        except Exception:
            return ""
    sha = run("rev-parse", "HEAD")
    if not sha:
        return {"sha": None}
    dirty = run("status", "--porcelain", "--untracked-files=no")
    return {"sha": sha, "branch": run("rev-parse", "--abbrev-ref", "HEAD"), "dirty": bool(dirty)}


def provenance(root: str = PROJECT_ROOT) -> Dict[str, Any]:
    """Паспорт кода и сборки: git, движок (модуль и сервер — хэши файлов), версии Python и numpy, платформа.
    Имя пользователя и машины не пишется."""
    eng = {"extension": None, "server": None}
    try:
        import synchrotronsim as ss
        if ss.ENGINE_AVAILABLE:
            eng["extension"] = {"build_info": ss._engine.build_info, "sha16": _sha16(ss._engine.__file__),
                                "file": os.path.basename(ss._engine.__file__)}
    except Exception:
        pass
    exe = os.path.join(root, "build", "MatrixTracker.exe" if os.name == "nt" else "MatrixTracker")
    if os.path.exists(exe):
        eng["server"] = {"file": os.path.basename(exe), "sha16": _sha16(exe)}
    return {"git": _git(root), "engine": eng, "python": sys.version.split()[0], "numpy": np.__version__,
            "platform": platform.platform(terse=True)}


def _jsonable(o):
    if is_dataclass(o) and not isinstance(o, type):
        return _jsonable(asdict(o))
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (str, int, float, bool)) or o is None:
        return o
    return str(o)


def _now() -> str:
    return _dt.datetime.now().isoformat(timespec="seconds")


def _human(v) -> str:
    """Число для человеко-читаемого CSV: десятичная запятая (Excel), 10 значащих."""
    return f"{float(v):.10g}".replace(".", ",")


# ---------------------------------------------------------------------------------------------- Run
class Run:
    """Прогон: папка с паспортом, рядами и артефактами + строка в индексе runs.csv."""

    def __init__(self, path: str, manifest: Dict[str, Any], writable: bool, root: str):
        self.path, self.manifest, self.writable, self.root = path, manifest, writable, root
        self._series: Dict[str, Dict[str, list]] = {}
        self._loaded = None

    # ---- создание ---------------------------------------------------------------------------------
    @classmethod
    def create(cls, kind: str, config: Any = None, machine=None, seed=None, source: str = "", note: str = "",
               root: Optional[str] = None, run_id: Optional[str] = None) -> "Run":
        """Новый прогон: папка, manifest (статус in_progress — прерванный прогон так и останется) и строка индекса.
        machine — synchrotronsim.model.Machine (берутся номер решётки, подпись, отпечаток) или None."""
        root = os.path.abspath(root or PROJECT_ROOT)
        ws = root if os.path.isdir(os.path.join(root, "tables")) else os.path.join(root, "workspace")
        stamp = _dt.datetime.now().strftime("%Y%m%d-%H%M%S")
        run_id = run_id or f"{stamp}-{kind}-{secrets.token_hex(2)}"
        path = os.path.join(ws, "runs", run_id)
        if os.path.exists(path):
            raise FileExistsError(f"прогон {run_id} уже есть: {path}")
        os.makedirs(os.path.join(path, "artifacts"))
        os.makedirs(os.path.join(path, "tables"))
        mach = None
        if machine is not None:
            mach = {"lattice_id": machine.lattice_id, "fingerprint": machine.fingerprint,
                    "label": machine.meta.label if machine.meta else None, "source": machine.source}
        man = {"format": "synchrotronsim.run/1", "run_id": run_id, "kind": kind, "source": source, "note": note,
               "status": "in_progress", "started": _now(), "finished": None, "seed": _jsonable(seed),
               "machine": mach, "config": _jsonable(config), "provenance": provenance(root),
               "signals": {}, "artifacts": {}}
        run = cls(path, man, True, ws)
        run._write_manifest()
        run._index_upsert()
        run.log(f"создан прогон {run_id} ({kind}{', ' + source if source else ''})")
        return run

    @classmethod
    def open(cls, run_id_or_path: str, root: Optional[str] = None) -> "Run":
        """Прогон на чтение: по run_id (в workspace/runs) или по пути к папке."""
        p = run_id_or_path
        if not os.path.isdir(p):
            r = os.path.abspath(root or PROJECT_ROOT)
            ws = r if os.path.isdir(os.path.join(r, "tables")) else os.path.join(r, "workspace")
            p = os.path.join(ws, "runs", run_id_or_path)
        with open(os.path.join(p, "manifest.json"), encoding="utf-8") as f:
            man = json.load(f)
        return cls(p, man, False, os.path.dirname(os.path.dirname(p)))

    @property
    def run_id(self) -> str:
        return self.manifest["run_id"]

    @property
    def status(self) -> str:
        return self.manifest["status"]

    # ---- запись -----------------------------------------------------------------------------------
    def _check(self):
        if not self.writable:
            raise RuntimeError(f"прогон {self.run_id} открыт на чтение")

    def record(self, name: str, step: int, value, t: Optional[float] = None) -> None:
        """Точка ряда: имя — адрес (`obs.ebpm.x:rms`, `rec.orbit_rms.x`, …; псевдонимы приводятся)."""
        self._check()
        key = AD.canonical(name)
        s = self._series.setdefault(key, {"step": [], "t": [], "value": []})
        s["step"].append(int(step)); s["t"].append(np.nan if t is None else float(t))
        s["value"].append(np.asarray(value, float) if not np.isscalar(value) else float(value))

    def record_rows(self, rows, prefix: str = "rec.", t_of_step: Optional[Dict[int, float]] = None) -> int:
        """Строки рекордеров (signal, step, value) → ряды `rec.<signal>`; t — по номеру окна, если дан."""
        n = 0
        for sig, step, val in rows:
            self.record(prefix + str(sig), step, val, (t_of_step or {}).get(int(step)))
            n += 1
        return n

    def artifact(self, name: str, arrays: Dict[str, Any], meta: str = "") -> str:
        """Тяжёлые массивы → artifacts/<name>.npz; запись в manifest. Возврат пути относительно папки прогона."""
        self._check()
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
        rel = f"artifacts/{safe}.npz"
        clean = {}
        for k, v in arrays.items():
            if isinstance(v, dict):                          # вложенный словарь (напр. спектр) — плоско
                for k2, v2 in v.items():
                    clean[f"{k}.{k2}"] = np.asarray(v2)
            elif v is not None:
                clean[k] = np.asarray(v)
        np.savez(os.path.join(self.path, rel), **clean)
        self.manifest["artifacts"][name] = {"path": rel, "meta": str(meta),
                                            "arrays": {k: list(np.shape(v)) for k, v in clean.items()}}
        return rel

    def log(self, msg: str) -> None:
        if not self.writable:
            return
        with open(os.path.join(self.path, "log.txt"), "a", encoding="utf-8") as f:
            f.write(f"{_now()}  {msg}\n")

    def flush(self) -> None:
        """Сохранить ряды и паспорт сейчас (контрольная точка; прогон остаётся in_progress)."""
        self._check()
        names = sorted(self._series)
        out = {"names": np.array(names, dtype=str)}
        info = {}
        for i, nm in enumerate(names):
            s = self._series[nm]
            try:
                val = np.asarray(s["value"], float)
            except ValueError:                               # векторы разной длины — по одному
                val = np.array([np.asarray(v, float).ravel() for v in s["value"]], dtype=object)
            if val.dtype == object:
                self.artifact(f"series_{i}", {f"v{j}": v for j, v in enumerate(val)},
                              meta=f"ряд {nm}: значения разной формы")
                val = np.full(len(s["step"]), np.nan)
            out[f"s{i}_step"] = np.asarray(s["step"], np.int64)
            out[f"s{i}_t"] = np.asarray(s["t"], float)
            out[f"s{i}_value"] = val
            m = AD.meta(nm)
            info[nm] = {"index": i, "n": len(s["step"]), "shape": list(val.shape[1:]), "units": m.units,
                        "name": m.name}
        np.savez(os.path.join(self.path, "signals.npz"), **out)
        self._write_human(names)
        self.manifest["signals"] = info
        self._write_manifest()

    def finish(self, status: str = "done", note: Optional[str] = None) -> "Run":
        """Завершить: ряды → signals.npz и tables/signals.csv, статус и время конца → manifest и runs.csv."""
        if status not in STATUSES:
            raise ValueError(f"статус {status!r}, допустимы {STATUSES}")
        self._check()
        self.manifest["status"] = status
        self.manifest["finished"] = _now()
        if note is not None:
            self.manifest["note"] = note
        self.flush()
        self._index_upsert()
        self.log(f"завершён: {status}")
        self.writable = False
        return self

    # ---- чтение -----------------------------------------------------------------------------------
    def names(self) -> List[str]:
        if self.writable:
            return sorted(self._series)
        return sorted(self.manifest.get("signals", {}))

    def series(self, name: str):
        """(step, t, value) ряда по адресу (псевдонимы приводятся)."""
        key = AD.canonical(name)
        if self.writable:
            s = self._series[key]
            return np.asarray(s["step"]), np.asarray(s["t"]), np.asarray(s["value"], float)
        if self._loaded is None:
            self._loaded = dict(np.load(os.path.join(self.path, "signals.npz"), allow_pickle=False))
        i = self.manifest["signals"][key]["index"]
        z = self._loaded
        return z[f"s{i}_step"], z[f"s{i}_t"], z[f"s{i}_value"]

    def load_artifact(self, name: str) -> Dict[str, np.ndarray]:
        a = self.manifest["artifacts"][name]
        return dict(np.load(os.path.join(self.path, a["path"]), allow_pickle=False))

    # ---- служебное --------------------------------------------------------------------------------
    def _write_manifest(self):
        tmp = os.path.join(self.path, "manifest.json.tmp")
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(self.manifest, f, ensure_ascii=False, indent=1)
        os.replace(tmp, os.path.join(self.path, "manifest.json"))

    def _write_human(self, names):
        """Скалярные ряды — в tables/signals.csv для человека (десятичная запятая)."""
        with open(os.path.join(self.path, "tables", "signals.csv"), "w", encoding="utf-8-sig", newline="") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(["адрес", "величина", "единицы", "шаг", "время_с", "значение"])
            for nm in names:
                s = self._series[nm]
                if s["value"] and not np.isscalar(s["value"][0]):
                    continue                                 # векторные ряды — только в signals.npz
                m = AD.meta(nm)
                for st, t, v in zip(s["step"], s["t"], s["value"]):
                    w.writerow([nm, m.name, m.units, st, "" if t != t else _human(t), _human(v)])

    def _index_upsert(self):
        tb = os.path.join(self.root, "tables")
        os.makedirs(tb, exist_ok=True)
        path = os.path.join(tb, "runs.csv")
        rows = []
        if os.path.exists(path):
            with open(path, encoding="utf-8-sig", newline="") as f:
                rows = [r for r in csv.DictReader(f, delimiter=";") if r.get("run_id") != self.run_id]
        m = self.manifest
        rows.append({"run_id": m["run_id"], "kind": m["kind"], "source": m.get("source", ""), "status": m["status"],
                     "started": m["started"], "finished": m.get("finished") or "",
                     "path": os.path.relpath(self.path, os.path.dirname(self.root)).replace(os.sep, "/"),
                     "note": str(m.get("note", "")).replace(";", ",").replace("\n", " ")})
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=INDEX_COLS, delimiter=";")
            w.writeheader(); w.writerows(rows)

    def __repr__(self):
        return f"Run({self.run_id}, {self.status}, рядов {len(self.names())})"


def list_runs(root: Optional[str] = None) -> List[Dict[str, str]]:
    """Строки индекса runs.csv (новые в конце)."""
    r = os.path.abspath(root or PROJECT_ROOT)
    ws = r if os.path.isdir(os.path.join(r, "tables")) else os.path.join(r, "workspace")
    path = os.path.join(ws, "tables", "runs.csv")
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter=";"))
