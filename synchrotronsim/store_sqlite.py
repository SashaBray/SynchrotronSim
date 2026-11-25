# -*- coding: utf-8 -*-
"""Хранилище с ревизиями на SQLite (docs/DATA_PROTOCOL.md, этап 2) — источник правды проекта.

Один файл `workspace/store/store.sqlite`: ревизии (машины, справочник устройств, опыты, документы машин из
пакетов), указатели последних ревизий, сжатые блобы. Ревизии неизменяемы; запись поверх устаревшей ревизии —
`ConflictError`. Содержимое ревизии машины — СРЕЗ её таблиц (synchrotronsim.tableset); документ `Machine`
строится из среза. Рабочая копия, статус, разница, фиксация и черновики — synchrotronsim.working (общая логика
для локального файла и сервера хранилища); здесь — только простые операции над файлом.
"""
from __future__ import annotations

import contextlib
import hashlib
import json
import os
import sqlite3
import threading
import zlib
from typing import Any, Dict, List, Optional

from .store import StoreError, ConflictError, experiment_fingerprint
from .working import WorkingCopyStore, LIB_ID
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


def _now():
    from .documents import now_iso
    return now_iso()


class SqliteStore(WorkingCopyStore):
    """Хранилище ревизий в файле SQLite + (необязательная) рабочая копия таблиц.

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
        self.dir = self.state_dir = os.path.dirname(self.path)
        self.state_name = "working.json"
        self.workspace = ws
        self.author = author or os.environ.get("USERNAME") or os.environ.get("USER") or ""
        self.site = "local"
        self._lock = threading.RLock()
        self._init_cache()
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

    def _rows(self, sql, args=()):
        with self._db() as c:
            c.row_factory = sqlite3.Row
            return [dict(r) for r in c.execute(sql, args).fetchall()]

    # ---------------------------------------------------------------------------------- блобы
    def blob_put(self, data: bytes, c=None) -> str:
        h = hashlib.sha256(data).hexdigest()
        z = zlib.compress(data, 6)
        sql = "INSERT OR IGNORE INTO blobs(hash, size, data) VALUES (?, ?, ?)"
        if c is not None:
            c.execute(sql, (h, len(data), z))
        else:
            with self._db() as cc:
                cc.execute(sql, (h, len(data), z))
        return h

    def blob_get(self, h: str) -> bytes:
        with self._db() as c:
            r = c.execute("SELECT data FROM blobs WHERE hash=?", (h,)).fetchone()
        if r is None:
            raise StoreError(f"нет блоба {h[:12]}")
        return zlib.decompress(r[0])

    def _content(self, row) -> Any:
        return json.loads(self.blob_get(row["blob"]).decode("utf-8"))

    # ---------------------------------------------------------------------------------- простые операции: чтение
    def head(self, kind: str, id) -> Optional[str]:
        with self._db() as c:
            r = c.execute("SELECT rev FROM heads WHERE kind=? AND id=?", (kind, str(id))).fetchone()
        return r[0] if r else None

    def revision(self, kind: str, id, rev: str) -> Dict[str, Any]:
        """Строка ревизии по точной ревизии, её префиксу или (машина) отпечатку документа Machine."""
        id = str(id)
        kinds = ("machine", "machine_doc") if kind == "machine" else (kind,)
        ph = ",".join("?" * len(kinds))
        hits = self._rows(f"SELECT * FROM revisions WHERE kind IN ({ph}) AND id=? AND (rev=? OR rev LIKE ?) ORDER BY seq",
                          (*kinds, id, rev, rev + "%"))
        if not hits and kind == "machine":
            hits = self._rows("SELECT * FROM revisions WHERE kind IN ('machine', 'machine_doc') AND id=? AND "
                              "(fingerprint=? OR fingerprint LIKE ?) ORDER BY seq", (id, rev, rev + "%"))
            if len({h["fingerprint"] for h in hits}) == 1:
                hits = [next((h for h in reversed(hits) if h["kind"] == "machine"), hits[-1])]
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
            return self.list_runs()
        out = self._rows("SELECT h.id AS id, h.rev AS rev, r.created AS created, r.note AS note, "
                         "r.fingerprint AS fingerprint FROM heads h JOIN revisions r ON r.kind=h.kind AND "
                         "r.id=h.id AND r.rev=h.rev WHERE h.kind=? ORDER BY h.id", (kind,))
        if kind == "machine":
            out.sort(key=lambda r: int(r["id"]) if r["id"].isdigit() else 10 ** 9)
        return out

    def library_content(self, rev: str) -> Dict[str, Any]:
        return self._content(self.revision("library", LIB_ID, rev))

    def machine_payload(self, id, rev: str) -> Dict[str, Any]:
        row = self.revision("machine", str(id), rev)
        if row["kind"] == "machine_doc":
            return {"kind": "machine_doc", "rev": row["rev"], "fingerprint": row["fingerprint"],
                    "doc": self.blob_get(row["blob"])}
        return {"kind": "machine", "rev": row["rev"], "fingerprint": row["fingerprint"], "slice": self._content(row)}

    def experiment_content(self, id, rev: str) -> Dict[str, Any]:
        return self._content(self.revision("experiment", str(id), rev))

    def find_machine_fp(self, id, fp: str) -> Optional[str]:
        hit = self._rows("SELECT rev FROM revisions WHERE kind IN ('machine', 'machine_doc') AND id=? AND "
                         "fingerprint=? ORDER BY CASE kind WHEN 'machine' THEN 0 ELSE 1 END, draft, seq LIMIT 1",
                         (str(id), fp))
        return hit[0]["rev"] if hit else None

    # ---------------------------------------------------------------------------------- простые операции: запись
    def _insert(self, c, kind, id, rev, content_bytes, parent, note, fingerprint=None, draft=False, refs=None,
                head=True):
        """В одной транзакции: проверка родителя (если двигаем указатель), блоб, ревизия, указатель. → создана ли."""
        c.execute("BEGIN IMMEDIATE")
        try:
            cur = c.execute("SELECT rev FROM heads WHERE kind=? AND id=?", (kind, id)).fetchone()
            cur = cur[0] if cur else None
            if head and cur is not None and cur != parent and cur != rev:
                raise ConflictError(f"{kind} {id}: последняя ревизия {cur[:12]}, а правка сделана поверх "
                                    f"{(parent or 'ничего')[:12]} — перенесите правку на последнюю ревизию")
            h = self.blob_put(content_bytes, c)
            new = c.execute("SELECT 1 FROM revisions WHERE kind=? AND id=? AND rev=?", (kind, id, rev)).fetchone() is None
            if new:
                c.execute("INSERT INTO revisions(kind, id, rev, parent, created, author, note, fingerprint, draft, "
                          "blob, refs) VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                          (kind, id, rev, parent if head else (parent or cur), _now(), self.author, note, fingerprint,
                           int(bool(draft)), h, json.dumps(refs, ensure_ascii=False) if refs else None))
            elif head:                                            # черновик становится зафиксированной ревизией
                c.execute("UPDATE revisions SET draft=0, note=CASE WHEN ?<>'' THEN ? ELSE note END, author=? "
                          "WHERE kind=? AND id=? AND rev=?", (note or "", note or "", self.author, kind, id, rev))
            if head:
                c.execute("INSERT INTO heads(kind, id, rev) VALUES (?,?,?) ON CONFLICT(kind, id) DO UPDATE SET "
                          "rev=excluded.rev", (kind, id, rev))
            c.execute("COMMIT")
            return new
        except BaseException:
            c.execute("ROLLBACK")
            raise

    def put_library(self, lib, parent, note, draft=False) -> str:
        rev = T.content_rev(lib)
        with self._db() as c:
            self._insert(c, "library", LIB_ID, rev, T.canonical(lib), parent, note, head=not draft)
        return rev

    def put_machine_slice(self, id, sl, parent, note, draft=False) -> str:
        if "library" not in sl:
            raise StoreError("срез машины без ссылки на ревизию справочника (поле library)")
        lib = self.library_content(sl["library"])
        m = T.machine_from_slices(sl, lib)                         # отпечаток считает хранилище (и проверяет срез)
        rev = T.content_rev(sl)
        with self._db() as c:
            self._insert(c, "machine", str(id), rev, T.canonical(sl), parent, note, fingerprint=m.fingerprint,
                         refs={"library": sl["library"]}, head=not draft)
        self._mcache["machine:" + rev] = m
        return rev

    def put_experiment(self, id, d, parent, note, draft=False) -> str:
        from .experiment import Experiment
        exp = Experiment.from_dict(d)
        rev = experiment_fingerprint(exp)
        data = json.dumps(exp.to_dict(), ensure_ascii=False, sort_keys=True).encode("utf-8")
        with self._db() as c:
            self._insert(c, "experiment", str(id), rev, data, parent, note, fingerprint=rev, draft=draft,
                         head=not draft)
        return rev

    def put_machine_doc(self, id, doc: bytes, note) -> str:
        from .documents import machine_from_bytes
        m = machine_from_bytes(doc)                                # документ проверяется (отпечаток) при приёме
        with self._db() as c:
            self._insert(c, "machine_doc", str(id), m.fingerprint, doc, None, note, fingerprint=m.fingerprint,
                         draft=True, head=False)
        return m.fingerprint

    # ---- простые JSON-документы (профили площадок и т.п.; секретов в них нет)
    JSON_KINDS = ("site",)

    def put_json(self, kind: str, id, d: Dict[str, Any], parent: Optional[str] = None, note: str = "") -> str:
        if kind not in self.JSON_KINDS:
            raise StoreError(f"put_json: вид {kind!r} не из {self.JSON_KINDS}")
        data = json.dumps(d, ensure_ascii=False, sort_keys=True).encode("utf-8")
        rev = hashlib.sha256(data).hexdigest()
        with self._db() as c:
            self._insert(c, kind, str(id), rev, data, parent, note)
        return rev

    def json_content(self, kind: str, id, rev: Optional[str] = None) -> Dict[str, Any]:
        rev = rev or self.head(kind, id)
        if rev is None:
            raise StoreError(f"нет документа {kind} {id}")
        return self._content(self.revision(kind, str(id), rev))

    # ---------------------------------------------------------------------------------- прогоны (этап 3: с вычислителя)
    def put_run(self, run_id: str, bundle: bytes, manifest: Dict[str, Any]) -> str:
        """Прогон целиком (zip папки прогона) — блобом; ревизия = sha256 содержимого."""
        rev = hashlib.sha256(bundle).hexdigest()
        with self._db() as c:
            self._insert(c, "run", str(run_id), rev, bundle, None, manifest.get("note", ""),
                         fingerprint=(manifest.get("machine") or {}).get("fingerprint"),
                         refs={"status": manifest.get("status"), "kind": manifest.get("kind"),
                               "machine": (manifest.get("machine") or {}).get("ref"),
                               "started": manifest.get("started"), "finished": manifest.get("finished")})
        return rev

    def run_bundle(self, run_id: str) -> bytes:
        h = self.head("run", run_id)
        if h is None:
            raise StoreError(f"нет прогона {run_id} в хранилище")
        return self.blob_get(self.revision("run", run_id, h)["blob"])

    def list_runs(self) -> List[Dict[str, Any]]:
        rows = self._rows("SELECT r.id AS run_id, r.created, r.note, r.refs FROM heads h JOIN revisions r ON "
                          "r.kind='run' AND r.id=h.id AND r.rev=h.rev WHERE h.kind='run' ORDER BY r.seq")
        for r in rows:
            r.update(json.loads(r.pop("refs") or "{}"))
        if self.working:
            from .run import list_runs
            have = {r["run_id"] for r in rows}
            rows += [dict(r, where="local") for r in list_runs(self.working) if r.get("run_id") not in have]
        return rows

    def get_run(self, id):
        """Прогон: из хранилища (распаковывается в кэш workspace/runs) или локальная папка прогона."""
        from .run import Run
        if self.head("run", id) is not None:
            return unpack_run(self.run_bundle(id), self.working or os.path.dirname(self.dir), id)
        return Run.open(id, root=self.working)


def pack_run(run) -> bytes:
    """Папка прогона → zip (манифест, ряды, таблицы, артефакты, журнал)."""
    import io
    import zipfile
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        for dp, _, fn in os.walk(run.path):
            for f in fn:
                p = os.path.join(dp, f)
                z.write(p, os.path.relpath(p, run.path))
    return buf.getvalue()


def unpack_run(bundle: bytes, root: str, run_id: str):
    """zip прогона → папка <workspace>/runs/<run_id> (если её нет) → Run."""
    import io
    import zipfile
    from .run import Run
    r = os.path.abspath(root)
    ws = r if os.path.isdir(os.path.join(r, "tables")) else os.path.join(r, "workspace")
    dest = os.path.join(ws, "runs", run_id)
    if not os.path.isdir(dest):
        tmp = dest + ".part"
        with zipfile.ZipFile(io.BytesIO(bundle)) as z:
            for n in z.namelist():
                if n.startswith(("/", "\\")) or ".." in n.replace("\\", "/").split("/"):
                    raise StoreError(f"прогон {run_id}: недопустимый путь в архиве {n!r}")
            z.extractall(tmp)
        os.replace(tmp, dest)
    return Run.open(dest)


def find_store_file(root: Optional[str] = None) -> Optional[str]:
    from .db import PROJECT_ROOT
    r = os.path.abspath(root or PROJECT_ROOT)
    ws = r if os.path.isdir(os.path.join(r, "tables")) else os.path.join(r, "workspace")
    p = os.path.join(ws, "store", "store.sqlite")
    return p if os.path.exists(p) else None
