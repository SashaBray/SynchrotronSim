# -*- coding: utf-8 -*-
"""Хранилище с ревизиями на SQLite (docs/DATA_PROTOCOL.md, этап 2) — источник правды проекта.

Один файл `workspace/store/store.sqlite`: ревизии (машины, справочник устройств, опыты), указатели последних
ревизий, блобы (содержимое ревизий, документы машин). Ревизии неизменяемы; `commit` поверх устаревшей ревизии —
`ConflictError`. Содержимое ревизии машины — СРЕЗ её таблиц (synchrotronsim.tableset: строки её решётки, текст
ячеек без потерь, + паспорт и ссылка на ревизию справочника); документ `Machine` строится из среза.

РАБОЧАЯ КОПИЯ (как в git) — папка `workspace/tables`: скрипты и Excel работают с ней. `status` — что в ней
изменено относительно базовой ревизии (с которой она выгружена), `commit` — изменённое → новые ревизии
(база = родитель: чужая фиксация между выгрузкой и фиксацией — конфликт), `checkout` — ревизию в рабочую копию.
Модель SynchrotronSim на этом компьютере по умолчанию считает рабочую копию (`get("machine", id)` без ревизии);
прогон на незафиксированной рабочей копии фиксирует её ЧЕРНОВОЙ ревизией (`pin`) — воспроизводимость сохраняется.
"""
from __future__ import annotations

import contextlib
import json
import os
import re
import sqlite3
import threading
import zlib
from typing import Any, Dict, List, Optional

from .store import Store, StoreError, ConflictError, Ref, experiment_fingerprint
from . import tableset as T

SCHEMA_VERSION = 1
_SCHEMA = """
CREATE TABLE IF NOT EXISTS meta(key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS blobs(hash TEXT PRIMARY KEY, size INTEGER, data BLOB NOT NULL);
CREATE TABLE IF NOT EXISTS revisions(
    seq INTEGER PRIMARY KEY AUTOINCREMENT,
    kind TEXT NOT NULL, id TEXT NOT NULL, rev TEXT NOT NULL, parent TEXT, created TEXT, author TEXT, note TEXT,
    fingerprint TEXT, draft INTEGER NOT NULL DEFAULT 0, blob TEXT NOT NULL, refs TEXT,
    UNIQUE(kind, id, rev));
CREATE INDEX IF NOT EXISTS rev_fp ON revisions(kind, id, fingerprint);
CREATE TABLE IF NOT EXISTS heads(kind TEXT NOT NULL, id TEXT NOT NULL, rev TEXT NOT NULL, PRIMARY KEY(kind, id));
"""
LIB_ID = "devices"


def _now():
    from .documents import now_iso
    return now_iso()


class SqliteStore(Store):
    """Хранилище ревизий SQLite + (необязательная) рабочая копия таблиц.

    path    — файл хранилища (по умолчанию <workspace>/store/store.sqlite);
    working — корень рабочей копии (проект или папка с tables/); None — без рабочей копии (сервер хранилища)."""

    def __init__(self, path: Optional[str] = None, working: Optional[str] = None, author: str = ""):
        from .db import PROJECT_ROOT
        self.working = os.path.abspath(working) if working else None
        ws = None
        if self.working:
            ws = self.working if os.path.isdir(os.path.join(self.working, "tables")) else \
                os.path.join(self.working, "workspace")
        if path is None:
            if ws is None:
                ws = os.path.join(PROJECT_ROOT, "workspace")
            path = os.path.join(ws, "store", "store.sqlite")
        self.path = os.path.abspath(path)
        self.dir = os.path.dirname(self.path)
        self.workspace = ws
        self.author = author or os.environ.get("USERNAME") or os.environ.get("USER") or ""
        self.site = "local"
        self._lock = threading.RLock()
        self._mcache: Dict[str, Any] = {}                         # ревизия → Machine
        os.makedirs(self.dir, exist_ok=True)
        with self._db() as c:
            c.executescript(_SCHEMA)
            c.execute("INSERT OR IGNORE INTO meta(key, value) VALUES ('schema', ?)", (str(SCHEMA_VERSION),))
            v = int(c.execute("SELECT value FROM meta WHERE key='schema'").fetchone()[0])
        if v > SCHEMA_VERSION:
            raise StoreError(f"хранилище {self.path}: схема {v} новее, чем понимает эта версия SynchrotronSim")

    def __repr__(self):
        return f"SqliteStore({self.path}, рабочая копия {self.working or 'нет'})"

    @contextlib.contextmanager
    def _db(self):
        with self._lock:
            c = sqlite3.connect(self.path, timeout=60, isolation_level=None)
            try:
                c.execute("PRAGMA journal_mode=WAL")
                yield c
            finally:
                c.close()

    # ---------------------------------------------------------------------------------- блобы
    def blob_put(self, data: bytes, c=None) -> str:
        import hashlib
        h = hashlib.sha256(data).hexdigest()
        z = zlib.compress(data, 6)
        if c is not None:
            c.execute("INSERT OR IGNORE INTO blobs(hash, size, data) VALUES (?, ?, ?)", (h, len(data), z))
        else:
            with self._db() as cc:
                cc.execute("INSERT OR IGNORE INTO blobs(hash, size, data) VALUES (?, ?, ?)", (h, len(data), z))
        return h

    def blob_get(self, h: str) -> bytes:
        with self._db() as c:
            r = c.execute("SELECT data FROM blobs WHERE hash=?", (h,)).fetchone()
        if r is None:
            raise StoreError(f"нет блоба {h[:12]}")
        return zlib.decompress(r[0])

    def _content(self, row) -> Any:
        return json.loads(self.blob_get(row["blob"]).decode("utf-8"))

    # ---------------------------------------------------------------------------------- ревизии
    def head(self, kind: str, id) -> Optional[str]:
        with self._db() as c:
            r = c.execute("SELECT rev FROM heads WHERE kind=? AND id=?", (kind, str(id))).fetchone()
        return r[0] if r else None

    def _rows(self, sql, args=()):
        with self._db() as c:
            c.row_factory = sqlite3.Row
            return [dict(r) for r in c.execute(sql, args).fetchall()]

    def revision(self, kind: str, id, rev: str) -> Dict[str, Any]:
        """Строка ревизии по точной ревизии, её префиксу или (машина) отпечатку документа Machine."""
        id = str(id)
        hits = self._rows("SELECT * FROM revisions WHERE kind=? AND id=? AND (rev=? OR rev LIKE ?) ORDER BY seq",
                          (kind, id, rev, rev + "%"))
        if not hits and kind == "machine":
            hits = self._rows("SELECT * FROM revisions WHERE kind IN ('machine', 'machine_doc') AND id=? AND "
                              "(fingerprint=? OR fingerprint LIKE ?) ORDER BY seq", (id, rev, rev + "%"))
            if len({h["fingerprint"] for h in hits}) == 1:
                hits = hits[-1:]                                  # один документ — любая его ревизия
        if not hits:
            raise StoreError(f"нет ревизии {rev[:12]} у {kind} {id} в {self.path}")
        if len(hits) > 1:
            raise StoreError(f"ревизия {rev} у {kind} {id} неоднозначна: {[h['rev'][:12] for h in hits[:3]]}")
        return hits[0]

    def history(self, kind: str, id) -> List[Dict[str, Any]]:
        rows = self._rows("SELECT seq, kind, id, rev, parent, created, author, note, fingerprint, draft, refs "
                          "FROM revisions WHERE kind=? AND id=? ORDER BY seq", (kind, str(id)))
        head = self.head(kind, id)
        for r in rows:
            r["head"] = r["rev"] == head
        return rows

    def list(self, kind: str) -> List[Dict[str, Any]]:
        if kind == "run":
            from .run import list_runs
            return list_runs(self.working)
        out = self._rows("SELECT h.id AS id, h.rev AS rev, r.created AS created, r.note AS note, "
                         "r.fingerprint AS fingerprint FROM heads h JOIN revisions r ON r.kind=h.kind AND "
                         "r.id=h.id AND r.rev=h.rev WHERE h.kind=? ORDER BY h.id", (kind,))
        if kind == "machine":
            out.sort(key=lambda r: int(r["id"]) if r["id"].isdigit() else 10 ** 9)
        return out

    def _insert(self, c, kind, id, rev, content_bytes, parent, note, fingerprint=None, draft=False, refs=None,
                head=True, expect=None, check_parent=True):
        """В одной транзакции: проверка родителя, блоб, ревизия, указатель. → (rev, создана_ли)."""
        c.execute("BEGIN IMMEDIATE")
        try:
            cur = c.execute("SELECT rev FROM heads WHERE kind=? AND id=?", (kind, id)).fetchone()
            cur = cur[0] if cur else None
            if head and check_parent and cur is not None and cur != (expect if expect is not None else parent):
                raise ConflictError(f"{kind} {id}: последняя ревизия {cur[:12]}, а правка сделана поверх "
                                    f"{(expect or parent or 'ничего')[:12]} — перенесите правку на последнюю ревизию")
            h = self.blob_put(content_bytes, c)
            new = c.execute("SELECT 1 FROM revisions WHERE kind=? AND id=? AND rev=?", (kind, id, rev)).fetchone() is None
            if new:
                c.execute("INSERT INTO revisions(kind, id, rev, parent, created, author, note, fingerprint, draft, "
                          "blob, refs) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                          (kind, id, rev, parent, _now(), self.author, note, fingerprint, int(bool(draft)), h,
                           json.dumps(refs, ensure_ascii=False) if refs else None))
            elif head:                                            # черновик становится зафиксированной ревизией
                c.execute("UPDATE revisions SET draft=0, note=CASE WHEN ?<>'' THEN ? ELSE note END, author=? "
                          "WHERE kind=? AND id=? AND rev=?", (note or "", note or "", self.author, kind, id, rev))
            if head:
                c.execute("INSERT INTO heads(kind, id, rev) VALUES (?,?,?) ON CONFLICT(kind, id) DO UPDATE SET "
                          "rev=excluded.rev", (kind, id, rev))
            c.execute("COMMIT")
            return rev, new
        except BaseException:
            c.execute("ROLLBACK")
            raise

    # ---------------------------------------------------------------------------------- машины и справочник
    def _library_at(self, rev: str) -> Dict[str, Any]:
        return self._content(self.revision("library", LIB_ID, rev))

    def machine_slice_at(self, id, rev: str) -> Dict[str, Any]:
        row = self.revision("machine", id, rev)
        if row["kind"] != "machine":
            raise StoreError(f"ревизия {rev[:12]} машины {id} — документ без среза таблиц (из пакета/черновика)")
        return self._content(row)

    def machine_at(self, id, rev: str):
        """Machine ревизии (из среза таблиц + ревизии справочника; документ из пакета — как есть). Кэш по ревизии."""
        row = self.revision("machine", id, rev)
        key = row["kind"] + ":" + row["rev"]
        m = self._mcache.get(key)
        if m is None:
            if row["kind"] == "machine_doc":
                from .documents import machine_from_bytes
                m = machine_from_bytes(self.blob_get(row["blob"]), source=f"store:machine/{id}@{row['rev'][:12]}")
            else:
                sl = self._content(row)
                lib = self._library_at(sl["library"])
                m = T.machine_from_slices(sl, lib, source=f"store:machine/{id}@{row['rev']}")
                if row["fingerprint"] and m.fingerprint != row["fingerprint"]:
                    raise StoreError(f"машина {id}@{row['rev'][:12]}: отпечаток собранного документа "
                                     f"{m.fingerprint[:12]} ≠ записанного {row['fingerprint'][:12]}")
            self._mcache[key] = m
        return m

    def _working_machine(self, id):
        from .db import Repository
        return Repository(self.working).machine(int(id))

    # ---- состояние рабочей копии: базовые ревизии (с какой ревизии выгружена каждая машина/справочник)
    def _wfile(self) -> str:
        return os.path.join(self.dir, "working.json")

    def working_base(self) -> Dict[str, str]:
        try:
            with open(self._wfile(), encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _set_base(self, **kv):
        b = self.working_base()
        b.update({k: v for k, v in kv.items() if v})
        tmp = self._wfile() + ".part"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump(b, f, ensure_ascii=False, indent=1, sort_keys=True)
        os.replace(tmp, self._wfile())

    def _need_working(self):
        if not self.working:
            raise StoreError("у этого хранилища нет рабочей копии таблиц (сервер)")

    def working_slices(self, lattices=None):
        """Срезы рабочей копии: (справочник, {id: срез машины с ссылкой на ревизию справочника})."""
        self._need_working()
        lib = T.library_slice(self.working)
        lrev = T.content_rev(lib)
        out = {}
        for lat in (lattices or T.lattices_in(self.working)):
            sl = T.machine_slice(self.working, lat)
            sl["library"] = lrev
            out[str(int(lat))] = sl
        return lib, out

    # ---------------------------------------------------------------------------------- чтение (протокол)
    def get(self, kind: str, id, rev: Optional[str] = None):
        id = str(id)
        if kind == "machine":
            if not rev and self.working:
                return self._working_machine(id)                  # рабочая копия (правки видны сразу)
            if not rev or rev == "head":
                h = self.head("machine", id)
                if h is None:
                    raise StoreError(f"у машины {id} нет зафиксированных ревизий")
                rev = h
            try:
                return self.machine_at(id, rev)
            except StoreError:
                legacy = self._legacy_snapshot(id, rev)           # снимки этапа 1 (workspace/store/machine/…)
                if legacy is not None:
                    return legacy
                raise
        if kind == "library":
            return self._library_at(rev or self.head("library", LIB_ID) or "")
        if kind == "experiment":
            return self._experiment(id, rev)
        if kind == "run":
            from .run import Run
            return Run.open(id, root=self.working)
        raise StoreError(f"неизвестный вид документа {kind!r}")

    def _legacy_snapshot(self, id, rev):
        d = os.path.join(self.dir, "machine", str(id))
        if not os.path.isdir(d):
            return None
        hits = [f for f in os.listdir(d) if f.endswith(".ssd") and f.startswith(rev)]
        if len(hits) != 1:
            return None
        from .documents import load_machine
        return load_machine(os.path.join(d, hits[0]))

    def _experiment(self, id: str, rev: Optional[str]):
        from .experiment import Experiment
        m = re.fullmatch(r"(?:fb:|fb_exp_)(\d+)", id)
        if m and not rev:
            return Experiment.from_fb(int(m.group(1)))
        if not rev:
            rev = self.head("experiment", id)
            if rev is None:
                if os.path.exists(id) and id.endswith(".json"):
                    return Experiment.from_json(id)
                raise StoreError(f"нет опыта {id!r} (ни fb:N, ни ревизий в хранилище, ни файла)")
        return Experiment.from_dict(self._content(self.revision("experiment", id, rev)))

    # ---------------------------------------------------------------------------------- запись (протокол)
    def commit(self, kind: str, id, doc=None, parent: Optional[str] = None, note: str = "") -> Ref:
        id = str(id)
        if kind == "experiment":
            rev = experiment_fingerprint(doc)
            data = json.dumps(doc.to_dict(), ensure_ascii=False, sort_keys=True).encode("utf-8")
            with self._db() as c:
                self._insert(c, "experiment", id, rev, data, parent, note, fingerprint=rev)
            return Ref("experiment", id, rev)
        if kind == "machine":
            if doc is None:                                       # из рабочей копии
                refs = self.commit_working(note=note, lattices=[int(id)], parents={f"machine/{id}": parent}
                                           if parent is not None else None)
                return next((r for r in refs if r.kind == "machine"), Ref("machine", id, self.head("machine", id)))
            if isinstance(doc, dict):
                return self._commit_machine_slice(id, doc, parent, note)
            raise StoreError("commit машины: срез таблиц (tableset.machine_slice) или рабочая копия; документ Machine "
                             "без таблиц фиксируется как черновик — pin(machine)")
        if kind == "library":
            return self._commit_library(doc, parent, note)
        raise StoreError(f"commit: машина, справочник или опыт, не {kind!r}")

    def _commit_library(self, lib, parent, note, c=None, draft=False):
        rev = T.content_rev(lib)
        data = T.canonical(lib)
        with self._db() as cc:
            self._insert(cc, "library", LIB_ID, rev, data, parent, note, head=not draft)
        return Ref("library", LIB_ID, rev)

    def _commit_machine_slice(self, id, sl, parent, note, draft=False):
        if "library" not in sl:
            raise StoreError("срез машины без ссылки на ревизию справочника (поле library)")
        lib = self._library_at(sl["library"])
        m = T.machine_from_slices(sl, lib)
        rev = T.content_rev(sl)
        with self._db() as c:
            self._insert(c, "machine", id, rev, T.canonical(sl), parent, note, fingerprint=m.fingerprint,
                         refs={"library": sl["library"]}, head=not draft)
        self._mcache["machine:" + rev] = m
        return Ref("machine", id, rev)

    def commit_working(self, note: str = "", lattices=None, parents: Optional[Dict[str, str]] = None) -> List[Ref]:
        """Изменённое в рабочей копии → новые ревизии (справочник, затем машины). Родитель — базовая ревизия
        рабочей копии; чужая фиксация с тех пор — ConflictError (ничего не фиксируется по этой машине)."""
        lib, slices = self.working_slices(lattices)
        base = dict(self.working_base())
        base.update(parents or {})
        out = []
        lrev = T.content_rev(lib)
        if lrev not in (self.head("library", LIB_ID), base.get("library")):   # не изменён локально — не трогаем
            out.append(self._commit_library(lib, base.get("library"), note))
            self._set_base(library=lrev)
        for id, sl in slices.items():
            rev = T.content_rev(sl)
            key = f"machine/{id}"
            if rev in (self.head("machine", id), base.get(key)):
                continue
            out.append(self._commit_machine_slice(id, sl, base.get(key), note))
            self._set_base(**{key: rev})
        for id, sl in slices.items():                            # первые ревизии / совпавшие с последней — база
            if self.head("machine", id) == T.content_rev(sl):
                self._set_base(**{f"machine/{id}": T.content_rev(sl)})
        if self.head("library", LIB_ID) == lrev:
            self._set_base(library=lrev)
        return out

    def pin(self, doc, id: Optional[str] = None, note: str = "") -> Ref:
        """Зафиксировать ревизию документа для прогона: машина — уже зафиксированная ревизия с тем же отпечатком,
        иначе черновая ревизия из рабочей копии (если документ из неё), иначе документ целиком (пакет).
        Опыт — черновая ревизия. Указатель последней ревизии не трогается."""
        from .model import Machine
        from .experiment import Experiment
        if isinstance(doc, Machine):
            mid = str(id or doc.lattice_id)
            hit = self._rows("SELECT rev FROM revisions WHERE kind IN ('machine', 'machine_doc') AND id=? AND "
                             "fingerprint=? ORDER BY draft, seq LIMIT 1", (mid, doc.fingerprint))
            if hit:
                return Ref("machine", mid, hit[0]["rev"])
            if self.working:
                try:
                    lib, slices = self.working_slices([int(mid)])
                    sl = slices[mid]
                    if T.machine_from_slices(sl, lib).fingerprint == doc.fingerprint:
                        lrev = T.content_rev(lib)
                        if lrev not in {r["rev"] for r in self.history("library", LIB_ID)}:
                            self._commit_library(lib, None, note or "черновик (прогон)", draft=True)
                        return self._commit_machine_slice(mid, sl, self.head("machine", mid),
                                                          note or "черновик (прогон на рабочей копии)", draft=True)
                except (KeyError, StoreError):
                    pass
            from .documents import machine_to_bytes
            with self._db() as c:
                self._insert(c, "machine_doc", mid, doc.fingerprint, machine_to_bytes(doc, machine_id=mid), None,
                             note or "документ машины (черновик)", fingerprint=doc.fingerprint, draft=True, head=False)
            return Ref("machine", mid, doc.fingerprint)
        if isinstance(doc, Experiment):
            eid = str(id or doc.name)
            rev = experiment_fingerprint(doc)
            data = json.dumps(doc.to_dict(), ensure_ascii=False, sort_keys=True).encode("utf-8")
            with self._db() as c:
                self._insert(c, "experiment", eid, rev, data, None, note or "черновик", fingerprint=rev, draft=True,
                             head=False)
            return Ref("experiment", eid, rev)
        raise StoreError(f"pin: документы машины и опыта, не {type(doc).__name__}")

    # ---------------------------------------------------------------------------------- рабочая копия: статус/разница/выгрузка
    def status(self, lattices=None) -> List[Dict[str, Any]]:
        """Состояние рабочей копии: по справочнику и каждой машине — база, последняя, рабочая; clean / modified /
        new (нет в хранилище) / behind (рабочая = база, но в хранилище есть новее) / missing (нет в таблицах)."""
        lib, slices = self.working_slices(lattices)
        base = self.working_base()
        out = []

        def st(kind, id, wrev, key):
            head = self.head(kind, id)
            b = base.get(key)
            if head is None:
                s = "new"
            elif wrev == head:
                s = "clean"
            elif b and wrev == b:
                s = "behind"
            else:
                s = "modified"
            return dict(kind=kind, id=id, state=s, working=wrev, head=head, base=b)
        out.append(st("library", LIB_ID, T.content_rev(lib), "library"))
        for id, sl in slices.items():
            out.append(st("machine", id, T.content_rev(sl), f"machine/{id}"))
        if lattices is None:
            for r in self.list("machine"):
                if r["id"] not in slices:
                    out.append(dict(kind="machine", id=r["id"], state="missing", working=None, head=r["rev"],
                                    base=base.get(f"machine/{r['id']}")))
        return out

    def diff(self, kind: str, id, a: Optional[str] = None, b: Optional[str] = None):
        """Разница по ячейкам: a — ревизия (по умолчанию последняя), b — ревизия или None (рабочая копия)."""
        id = str(id)

        def load(rev):
            if rev is None:
                lib, sl = self.working_slices([int(id)] if kind == "machine" else [])
                return lib if kind == "library" else sl[id]
            if kind == "library":
                return self._library_at(rev)
            return self.machine_slice_at(id, rev)
        a = a or self.head(kind, LIB_ID if kind == "library" else id)
        A = load(a) if a else None
        return T.diff_slices(A, load(b))

    def checkout(self, kind: str, id, rev: Optional[str] = None, force: bool = False,
                 with_library: bool = False) -> List[str]:
        """Ревизию — в рабочую копию (по умолчанию последнюю). Изменённая рабочая копия — отказ (force — затереть).
        Машина: строки её решётки; справочник ревизии — только with_library (он общий для всех машин)."""
        self._need_working()
        id = str(id)
        key = "library" if kind == "library" else f"machine/{id}"
        rid = LIB_ID if kind == "library" else id
        rev = self.revision(kind, rid, rev).get("rev") if rev else self.head(kind, rid)
        if rev is None:
            raise StoreError(f"у {kind} {id} нет ревизий")
        if not force:
            s = next((r for r in self.status([int(id)] if kind == "machine" else []) if r["kind"] == kind and
                      (kind == "library" or r["id"] == id)), None)
            if s and s["state"] in ("modified", "new"):
                raise StoreError(f"рабочая копия {kind} {id} изменена и не зафиксирована — commit или force")
        changed = []
        if kind == "library":
            changed += T.write_library(self.working, self._library_at(rev))
            self._set_base(library=rev)
            return changed
        sl = self.machine_slice_at(id, rev)
        changed += T.write_machine(self.working, sl)
        if with_library:
            changed += T.write_library(self.working, self._library_at(sl["library"]))
            self._set_base(library=sl["library"])
        self._set_base(**{key: rev})
        return changed

    def init_from_working(self, note: str = "импорт таблиц рабочей копии") -> List[Ref]:
        """Первые ревизии из рабочей копии: справочник + все решётки (повторный вызов фиксирует только изменения)."""
        return self.commit_working(note=note)


def find_store_file(root: Optional[str] = None) -> Optional[str]:
    from .db import PROJECT_ROOT
    r = os.path.abspath(root or PROJECT_ROOT)
    ws = r if os.path.isdir(os.path.join(r, "tables")) else os.path.join(r, "workspace")
    p = os.path.join(ws, "store", "store.sqlite")
    return p if os.path.exists(p) else None
