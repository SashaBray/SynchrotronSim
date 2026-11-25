# -*- coding: utf-8 -*-
"""Хранилище ревизий с рабочей копией таблиц — общая логика поверх ПРОСТЫХ операций хранилища
(docs/DATA_PROTOCOL.md, этапы 2–3). Простые операции реализуют `SqliteStore` (локальный файл) и `RemoteStore`
(сервер хранилища по сети) — одинаково; рабочая копия (`workspace/tables`), её статус, разница, фиксация,
выгрузка, черновики прогонов и сборка документа `Machine` из среза — здесь, на компьютере пользователя.

Простые операции (их и только их передаёт сеть):
    head(kind, id)                         → ревизия | None
    revision(kind, id, rev)                → строка ревизии (по номеру, префиксу или отпечатку Machine)
    history(kind, id), list(kind)
    library_content(rev)                   → срез справочника
    machine_payload(id, rev)               → {"kind": "machine", "slice", "fingerprint", "rev"} | {"kind": "machine_doc", "doc"}
    experiment_content(id, rev)            → dict опыта
    put_library(lib, parent, note, draft)  → ревизия
    put_machine_slice(id, sl, parent, note, draft) → ревизия (отпечаток Machine считает само хранилище)
    put_experiment(id, d, parent, note, draft)     → ревизия
    put_machine_doc(id, doc, note)         → отпечаток (документ машины без среза — пакет опыта, черновик)
    find_machine_fp(id, fp)                → ревизия | None
    blob_put(data) / blob_get(hash)
"""
from __future__ import annotations

import json
import os
import re
from typing import Any, Dict, List, Optional

from .store import Store, StoreError, ConflictError, Ref, experiment_fingerprint
from . import tableset as T

LIB_ID = "devices"


class WorkingCopyStore(Store):
    """Логика ревизий + рабочая копия. Подкласс задаёт простые операции, `self.working` (корень рабочей копии
    или None), `self.state_dir` (где хранить базы рабочей копии) и `self.state_name`."""

    working: Optional[str] = None
    state_dir: str = ""
    state_name: str = "working.json"

    def _init_cache(self):
        self._mcache: Dict[str, Any] = {}

    # ---------------------------------------------------------------------------------- простые операции
    def head(self, kind, id):                       # pragma: no cover — в подклассе
        raise NotImplementedError

    def revision(self, kind, id, rev):
        raise NotImplementedError

    def library_content(self, rev):
        raise NotImplementedError

    def machine_payload(self, id, rev):
        raise NotImplementedError

    def experiment_content(self, id, rev):
        raise NotImplementedError

    def put_library(self, lib, parent, note, draft=False):
        raise NotImplementedError

    def put_machine_slice(self, id, sl, parent, note, draft=False):
        raise NotImplementedError

    def put_experiment(self, id, d, parent, note, draft=False):
        raise NotImplementedError

    def put_machine_doc(self, id, doc: bytes, note):
        raise NotImplementedError

    def find_machine_fp(self, id, fp):
        raise NotImplementedError

    # ---------------------------------------------------------------------------------- машины
    def _library_at(self, rev: str) -> Dict[str, Any]:
        return self.library_content(rev)

    def machine_slice_at(self, id, rev: str) -> Dict[str, Any]:
        p = self.machine_payload(str(id), rev)
        if p["kind"] != "machine":
            raise StoreError(f"ревизия {rev[:12]} машины {id} — документ без среза таблиц (из пакета/черновика)")
        return p["slice"]

    def machine_at(self, id, rev: str):
        """Machine ревизии (из среза таблиц + ревизии справочника; документ из пакета — как есть). Кэш по ревизии."""
        id = str(id)
        p = self.machine_payload(id, rev)
        key = f"{p['kind']}:{p['rev']}"
        m = self._mcache.get(key)
        if m is None:
            if p["kind"] == "machine_doc":
                from .documents import machine_from_bytes
                m = machine_from_bytes(p["doc"], source=f"store:machine/{id}@{p['rev'][:12]}")
            else:
                sl = p["slice"]
                lib = self._library_at(sl["library"])
                m = T.machine_from_slices(sl, lib, source=f"store:machine/{id}@{p['rev']}")
                if p.get("fingerprint") and m.fingerprint != p["fingerprint"]:
                    raise StoreError(f"машина {id}@{p['rev'][:12]}: отпечаток собранного документа "
                                     f"{m.fingerprint[:12]} ≠ записанного {p['fingerprint'][:12]}")
            self._mcache[key] = m
        return m

    def _working_machine(self, id):
        from .db import Repository
        return Repository(self.working).machine(int(id))

    def _legacy_snapshot(self, id, rev):
        d = os.path.join(self.state_dir, "machine", str(id))
        if not os.path.isdir(d):
            return None
        hits = [f for f in os.listdir(d) if f.endswith(".ssd") and f.startswith(rev)]
        if len(hits) != 1:
            return None
        from .documents import load_machine
        return load_machine(os.path.join(d, hits[0]))

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
            return self.get_run(id)
        raise StoreError(f"неизвестный вид документа {kind!r}")

    def get_run(self, id):
        from .run import Run
        return Run.open(id, root=self.working)

    def _experiment(self, id: str, rev: Optional[str]):
        from .experiment import Experiment
        m = re.fullmatch(r"(?:fb:|fb_exp_)(\d+)", id)
        if m and not rev:
            if not self.working:
                raise StoreError("опыт fb:N читается из таблиц fb_* рабочей копии — у этого хранилища её нет; "
                                 "передайте опыт документом (commit) или ссылкой на его ревизию")
            return Experiment.from_fb(int(m.group(1)))
        if not rev:
            rev = self.head("experiment", id)
            if rev is None:
                if os.path.exists(id) and id.endswith(".json"):
                    return Experiment.from_json(id)
                raise StoreError(f"нет опыта {id!r} (ни fb:N, ни ревизий в хранилище, ни файла)")
        return Experiment.from_dict(self.experiment_content(id, rev))

    # ---------------------------------------------------------------------------------- запись (протокол)
    def commit(self, kind: str, id, doc=None, parent: Optional[str] = None, note: str = "") -> Ref:
        id = str(id)
        if kind == "experiment":
            return Ref("experiment", id, self.put_experiment(id, doc.to_dict(), parent, note))
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

    def _commit_library(self, lib, parent, note, draft=False) -> Ref:
        return Ref("library", LIB_ID, self.put_library(lib, parent, note, draft))

    def _commit_machine_slice(self, id, sl, parent, note, draft=False) -> Ref:
        if "library" not in sl:
            raise StoreError("срез машины без ссылки на ревизию справочника (поле library)")
        return Ref("machine", str(id), self.put_machine_slice(str(id), sl, parent, note, draft))

    def pin(self, doc, id: Optional[str] = None, note: str = "") -> Ref:
        """Зафиксировать ревизию документа для прогона: машина — уже зафиксированная ревизия с тем же отпечатком,
        иначе черновая ревизия из рабочей копии (если документ из неё), иначе документ целиком (пакет).
        Опыт — черновая ревизия. Указатель последней ревизии не трогается."""
        from .model import Machine
        from .experiment import Experiment
        if isinstance(doc, Machine):
            mid = str(id or doc.lattice_id)
            hit = self.find_machine_fp(mid, doc.fingerprint)
            if hit:
                return Ref("machine", mid, hit)
            if self.working:
                try:
                    lib, slices = self.working_slices([int(mid)])
                    sl = slices[mid]
                    if T.machine_from_slices(sl, lib).fingerprint == doc.fingerprint:
                        lrev = T.content_rev(lib)
                        try:
                            self.revision("library", LIB_ID, lrev)
                        except StoreError:
                            self._commit_library(lib, None, note or "черновик (прогон)", draft=True)
                        return self._commit_machine_slice(mid, sl, self.head("machine", mid),
                                                          note or "черновик (прогон на рабочей копии)", draft=True)
                except (KeyError, StoreError):
                    pass
            from .documents import machine_to_bytes
            fp = self.put_machine_doc(mid, machine_to_bytes(doc, machine_id=mid), note or "документ машины (черновик)")
            return Ref("machine", mid, fp)
        if isinstance(doc, Experiment):
            eid = str(id or doc.name)
            return Ref("experiment", eid, self.put_experiment(eid, doc.to_dict(), None, note or "черновик", draft=True))
        raise StoreError(f"pin: документы машины и опыта, не {type(doc).__name__}")

    # ---------------------------------------------------------------------------------- рабочая копия
    def _wfile(self) -> str:
        return os.path.join(self.state_dir, self.state_name)

    def working_base(self) -> Dict[str, str]:
        try:
            with open(self._wfile(), encoding="utf-8") as f:
                return json.load(f)
        except (OSError, ValueError):
            return {}

    def _set_base(self, **kv):
        b = self.working_base()
        b.update({k: v for k, v in kv.items() if v})
        os.makedirs(self.state_dir, exist_ok=True)
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
        for lat in (lattices if lattices is not None else T.lattices_in(self.working)):
            sl = T.machine_slice(self.working, lat)
            sl["library"] = lrev
            out[str(int(lat))] = sl
        return lib, out

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
            head = self.head("machine", id)
            if rev in (head, base.get(key)):
                continue
            if head is not None and not base.get(key):
                raise ConflictError(f"машина {id}: рабочая копия этого компьютера не сверялась с хранилищем "
                                    f"{self!r} (нет базовой ревизии), а в хранилище уже есть {head[:12]} — "
                                    f"неизвестно, поверх чего сделана правка. Сохраните правку (export), выгрузите "
                                    f"последнюю ревизию (checkout --force) и наложите правку (import)")
            out.append(self._commit_machine_slice(id, sl, base.get(key), note))
            self._set_base(**{key: rev})
        for id, sl in slices.items():                            # первые ревизии / совпавшие с последней — база
            if self.head("machine", id) == T.content_rev(sl):
                self._set_base(**{f"machine/{id}": T.content_rev(sl)})
        if self.head("library", LIB_ID) == lrev:
            self._set_base(library=lrev)
        return out

    def init_from_working(self, note: str = "импорт таблиц рабочей копии") -> List[Ref]:
        """Первые ревизии из рабочей копии: справочник + все решётки (повторный вызов фиксирует только изменения)."""
        return self.commit_working(note=note)

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
        clean = {("library" if r["kind"] == "library" else f"machine/{r['id']}"): r["head"] for r in out
                 if r["state"] == "clean" and r["base"] != r["head"]}
        if clean:                                                # совпадает с последней ревизией — это и есть база
            self._set_base(**clean)
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
