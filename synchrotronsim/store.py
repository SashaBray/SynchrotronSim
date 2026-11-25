# -*- coding: utf-8 -*-
"""Хранилище SynchrotronSim — протокол обмена библиотеки с БД документами (docs/DATA_PROTOCOL.md, раздел 5).

    from synchrotronsim.store import open_store
    st = open_store()                                  # по умолчанию — локальное хранилище над нынешними CSV
    m = st.get("machine", "26")                        # последняя ревизия машины (Machine)
    ref = st.pin(m)                                    # зафиксировать ревизию: снимок документа → Ref с rev
    st.get("machine", "26", ref.rev)                   # та же машина, даже если таблицы потом поправят
    e = st.get("experiment", "fb:3")                   # опыт ОС из таблиц fb_* (или имя JSON-опыта)
    st.commit("experiment", "my_exp", e, parent=None)  # новая ревизия опыта (конфликт, если parent устарел)

Этап 1 (`LocalCsvStore`): источник ПОСЛЕДНЕЙ ревизии машины — нынешние таблицы (`synchrotronsim.db.Repository`),
прежние ревизии — снимки документов `workspace/store/machine/<id>/<rev>.ssd` (кладутся `pin`, в т.ч. при каждом
прогоне — прогон воспроизводим, даже если таблицы потом изменятся). Запись машин ревизиями — этап 2 (SQLite);
сейчас `commit("machine", …)` отказывает: источник правды машин — CSV до импорта в хранилище.
Опыты — ревизии JSON `workspace/store/experiment/<id>/<rev>.json` + указатель последней; `fb:N` — из таблиц fb_*.
Прогоны — `synchrotronsim.run` (`workspace/runs`, индекс `tables/runs.csv`); блобы — `workspace/store/blobs/<sha256>`.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
from dataclasses import dataclass
from typing import Any, Dict, List, Optional

KINDS = ("machine", "experiment", "run")


class StoreError(RuntimeError):
    """Отказ хранилища: нет документа/ревизии, операция не поддержана на этом этапе."""


class ConflictError(StoreError):
    """commit поверх устаревшей ревизии: последняя ревизия ≠ parent (перезаписи молча не бывает)."""


@dataclass(frozen=True)
class Ref:
    """Ссылка на документ: вид + id + ревизия ("" — последняя). Текстовая форма `machine:26@6e176fb5f6fd`."""
    kind: str
    id: str
    rev: str = ""

    def __str__(self):
        return f"{self.kind}:{self.id}" + (f"@{self.rev}" if self.rev else "")

    @property
    def short(self) -> str:
        return f"{self.kind}:{self.id}" + (f"@{self.rev[:12]}" if self.rev else "")

    @classmethod
    def parse(cls, s: str) -> "Ref":
        m = re.fullmatch(r"(machine|experiment|run):([^@]+)(?:@(.+))?", str(s).strip())
        if not m:
            raise ValueError(f"ссылка вида «вид:id[@ревизия]», получено {s!r}")
        return cls(m.group(1), m.group(2), m.group(3) or "")


def experiment_fingerprint(exp) -> str:
    """Отпечаток содержимого опыта (каноническая кодировка to_dict)."""
    from .model import fingerprint_of
    return fingerprint_of("Experiment/1", exp.to_dict())


class Store:
    """Протокол хранилища (раздел 5 docs/DATA_PROTOCOL.md). Реализации: LocalCsvStore (этап 1); далее — SQLite
    (этап 2) и сетевой клиент (этап 3) с теми же операциями."""

    def get(self, kind: str, id: str, rev: Optional[str] = None):
        raise NotImplementedError

    def list(self, kind: str) -> List[Dict[str, Any]]:
        raise NotImplementedError

    def history(self, kind: str, id: str) -> List[Dict[str, Any]]:
        raise NotImplementedError

    def commit(self, kind: str, id: str, doc, parent: Optional[str] = None, note: str = "") -> Ref:
        raise NotImplementedError

    def pin(self, doc, id: Optional[str] = None, note: str = "") -> Ref:
        raise NotImplementedError

    def blob_put(self, data: bytes) -> str:
        raise NotImplementedError

    def blob_get(self, h: str) -> bytes:
        raise NotImplementedError

    def resolve(self, ref):
        ref = Ref.parse(ref) if isinstance(ref, str) else ref
        return self.get(ref.kind, ref.id, ref.rev or None)

    def machine(self, id, rev: Optional[str] = None):
        return self.get("machine", str(id), rev)


class LocalCsvStore(Store):
    """Локальное хранилище этапа 1: последняя ревизия машины — нынешние таблицы, прежние — снимки документов."""

    def __init__(self, root: Optional[str] = None, overlay: Optional[str] = None, store_dir: Optional[str] = None):
        from .db import Repository
        self.repo = Repository(root, overlay)
        self.root = root
        self.dir = os.path.abspath(store_dir) if store_dir else os.path.join(self.repo.workspace, "store")
        self.site = "local"

    def __repr__(self):
        return f"LocalCsvStore({self.repo.workspace}{', overlay ' + self.repo.overlay_dir if self.repo.overlay_dir else ''})"

    # ---- пути ----
    def _mdir(self, id: str) -> str:
        return os.path.join(self.dir, "machine", str(id))

    def _edir(self, id: str) -> str:
        return os.path.join(self.dir, "experiment", _safe(id))

    # ---- чтение ----
    def get(self, kind: str, id, rev: Optional[str] = None):
        id = str(id)
        if kind == "machine":
            return self._machine(id, rev)
        if kind == "experiment":
            return self._experiment(id, rev)
        if kind == "run":
            from .run import Run
            return Run.open(id, root=self.root)
        raise StoreError(f"неизвестный вид документа {kind!r} (есть {', '.join(KINDS)})")

    def _machine(self, id: str, rev: Optional[str]):
        try:
            head = self.repo.machine(int(id))
        except (KeyError, ValueError) as e:
            head, err = None, e
        if not rev or (head is not None and head.fingerprint.startswith(rev)):
            if head is None:
                raise StoreError(f"нет машины {id}: {err}")
            return head
        path = self._find_rev(self._mdir(id), rev, ".ssd")
        if path is None:
            raise StoreError(f"нет ревизии {rev[:12]} машины {id}: снимка нет в {self._mdir(id)}, а последняя "
                             f"ревизия — {head.fingerprint[:12] if head is not None else 'нет'}")
        from .documents import load_machine
        return load_machine(path)

    @staticmethod
    def _find_rev(d: str, rev: str, ext: str) -> Optional[str]:
        if not os.path.isdir(d):
            return None
        hits = [f for f in os.listdir(d) if f.endswith(ext) and f.startswith(rev)]
        if len(hits) > 1:
            raise StoreError(f"ревизия {rev} неоднозначна в {d}: {hits[:3]}")
        return os.path.join(d, hits[0]) if hits else None

    def _experiment(self, id: str, rev: Optional[str]):
        from .experiment import Experiment
        m = re.fullmatch(r"(?:fb:|fb_exp_)(\d+)", id)
        if m and not rev:
            return Experiment.from_fb(int(m.group(1)))
        d = self._edir(id)
        if not rev:
            head = os.path.join(d, "HEAD")
            if os.path.exists(head):
                rev = _read_text(head)
            elif os.path.exists(id) and id.endswith(".json"):          # путь к JSON-опыту
                return Experiment.from_json(id)
            else:
                raise StoreError(f"нет опыта {id!r} (ни fb:N, ни ревизий в {d}, ни файла)")
        path = self._find_rev(d, rev, ".json")
        if path is None:
            raise StoreError(f"нет ревизии {rev[:12]} опыта {id}")
        return Experiment.from_json(path)

    def list(self, kind: str) -> List[Dict[str, Any]]:
        if kind == "machine":
            out = []
            for lat in self.repo.lattices():
                from .machines import LATTICE_INFO
                info = LATTICE_INFO.get(lat) or {}
                out.append({"id": str(lat), "label": info.get("label", ""), "revisions": len(self.history("machine", lat)) })
            return out
        if kind == "experiment":
            out = [{"id": f"fb:{r['fb_exp_id']}", "source": "fb_experiments", "note": r.get("note", "")}
                   for r in self.repo.table("fb_experiments", required=False).dicts()]
            d = os.path.join(self.dir, "experiment")
            if os.path.isdir(d):
                for name in sorted(os.listdir(d)):
                    head = os.path.join(d, name, "HEAD")
                    if os.path.exists(head):
                        out.append({"id": name, "source": "store", "rev": _read_text(head)})
            return out
        if kind == "run":
            from .run import list_runs
            return list_runs(self.root)
        raise StoreError(f"неизвестный вид документа {kind!r}")

    def history(self, kind: str, id) -> List[Dict[str, Any]]:
        id = str(id)
        if kind == "machine":
            d, ext = self._mdir(id), ".ssd"
        elif kind == "experiment":
            d, ext = self._edir(id), ".json"
        else:
            raise StoreError(f"история есть у машин и опытов, не у {kind!r}")
        out = []
        idx = os.path.join(d, "index.jsonl")
        if os.path.exists(idx):
            with open(idx, encoding="utf-8") as f:
                out = [json.loads(line) for line in f if line.strip()]
        return out

    # ---- запись ----
    def _append_index(self, d: str, rec: Dict[str, Any]):
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "index.jsonl"), "a", encoding="utf-8") as f:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")

    def pin(self, doc, id: Optional[str] = None, note: str = "") -> Ref:
        """Зафиксировать ревизию документа (снимок по отпечатку; повтор — без записи). → Ref с ревизией."""
        from .model import Machine
        from .experiment import Experiment
        from .documents import machine_to_bytes, now_iso
        if isinstance(doc, Machine):
            mid = str(id or doc.lattice_id)
            d = self._mdir(mid)
            path = os.path.join(d, f"{doc.fingerprint}.ssd")
            if not os.path.exists(path):
                os.makedirs(d, exist_ok=True)
                tmp = path + ".part"
                with open(tmp, "wb") as f:
                    f.write(machine_to_bytes(doc, machine_id=mid, note=note))
                os.replace(tmp, path)
                self._append_index(d, {"rev": doc.fingerprint, "created": now_iso(), "note": note,
                                       "source": doc.source})
            return Ref("machine", mid, doc.fingerprint)
        if isinstance(doc, Experiment):
            eid = str(id or doc.name)
            rev = experiment_fingerprint(doc)
            d = self._edir(eid)
            path = os.path.join(d, f"{rev}.json")
            if not os.path.exists(path):
                os.makedirs(d, exist_ok=True)
                doc.to_json(path)
                self._append_index(d, {"rev": rev, "created": now_iso(), "note": note, "parent": None})
            return Ref("experiment", eid, rev)
        raise StoreError(f"pin: документы машины и опыта, не {type(doc).__name__}")

    def commit(self, kind: str, id, doc, parent: Optional[str] = None, note: str = "") -> Ref:
        id = str(id)
        if kind == "machine":
            raise StoreError("этап 1: машины пишутся в таблицы (CSV — источник до импорта в хранилище, этап 2); "
                             "зафиксировать ревизию для прогона — pin(machine)")
        if kind != "experiment":
            raise StoreError(f"commit: машина или опыт, не {kind!r}")
        from .documents import now_iso
        d = self._edir(id)
        head_p = os.path.join(d, "HEAD")
        head = _read_text(head_p) if os.path.exists(head_p) else None
        if head is not None and (not parent or not head.startswith(parent)):
            raise ConflictError(f"опыт {id}: последняя ревизия {head[:12]}, а правка сделана поверх "
                                f"{(parent or 'ничего')[:12]} — перенесите правку на последнюю ревизию")
        rev = experiment_fingerprint(doc)
        path = os.path.join(d, f"{rev}.json")
        os.makedirs(d, exist_ok=True)
        if not os.path.exists(path):
            doc.to_json(path)
        self._append_index(d, {"rev": rev, "created": now_iso(), "note": note, "parent": head})
        with open(head_p + ".part", "w", encoding="utf-8") as f:
            f.write(rev)
        os.replace(head_p + ".part", head_p)
        return Ref("experiment", id, rev)

    # ---- блобы ----
    def blob_put(self, data: bytes) -> str:
        h = hashlib.sha256(data).hexdigest()
        d = os.path.join(self.dir, "blobs", h[:2])
        p = os.path.join(d, h)
        if not os.path.exists(p):
            os.makedirs(d, exist_ok=True)
            with open(p + ".part", "wb") as f:
                f.write(data)
            os.replace(p + ".part", p)
        return h

    def blob_get(self, h: str) -> bytes:
        p = os.path.join(self.dir, "blobs", h[:2], h)
        if not os.path.exists(p):
            raise StoreError(f"нет блоба {h[:12]}")
        with open(p, "rb") as f:
            return f.read()


def _read_text(path: str) -> str:
    with open(path, encoding="utf-8") as f:
        return f.read().strip()


def _safe(id: str) -> str:
    s = re.sub(r"[^\w.\-]+", "_", str(id), flags=re.UNICODE).strip("._")
    return s or "_"


_STORES: Dict[tuple, Store] = {}


def open_store(site=None, root: Optional[str] = None, overlay: Optional[str] = None) -> Store:
    """Хранилище площадки: site=None / "local" — локальное (LocalCsvStore над таблицами проекта, или root /
    overlay); сетевые площадки — этап 3 (`synchrotronsim.Site`)."""
    if site not in (None, "local"):
        store = getattr(site, "store", None)
        if store is not None:
            return store
        raise StoreError(f"площадка {site!r}: сетевое хранилище — этап 3 протокола")
    key = (os.path.abspath(root) if root else None, os.path.abspath(overlay) if overlay else None)
    st = _STORES.get(key)
    if st is None:
        st = _STORES[key] = LocalCsvStore(root, overlay)
    return st
