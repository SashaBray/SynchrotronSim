# -*- coding: utf-8 -*-
"""Хранилище ревизий SynchrotronSim (docs/DATA_PROTOCOL.md, этап 2) — команды как у git.

Рабочая копия — workspace/tables (скрипты и Excel работают с ней), история — workspace/store/store.sqlite.

    py -3.14 ss_store.py init                         # первые ревизии: справочник устройств + все решётки
    py -3.14 ss_store.py status                       # что в рабочей копии изменено и не зафиксировано
    py -3.14 ss_store.py diff 26                      # разница по ячейкам: последняя ревизия → рабочая копия
    py -3.14 ss_store.py commit -m "заметка"          # изменённое → новые ревизии (конфликт, если кто-то успел)
    py -3.14 ss_store.py log 26                       # история машины (library — справочника)
    py -3.14 ss_store.py checkout 26@<ревизия>        # ревизию — в рабочую копию (--force — затереть правки)
    py -3.14 ss_store.py export 26 out/sila26         # для Excel: запятая, «;», UTF-8 с BOM; 26@<рев> | 26@working
    py -3.14 ss_store.py import out/sila26            # правки из Excel → рабочая копия (--commit -m … — сразу ревизия)
    py -3.14 ss_store.py --site lab status            # то же против сервера хранилища площадки (ss_server.py)
"""
import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8")

from synchrotronsim import tableset as T                        # noqa: E402
from synchrotronsim.store import StoreError, ConflictError     # noqa: E402
from synchrotronsim.store_sqlite import SqliteStore, find_store_file, LIB_ID   # noqa: E402

STATE = {"clean": "без изменений", "modified": "ИЗМЕНЕНА", "new": "НОВАЯ (нет в хранилище)",
         "behind": "отстала (в хранилище есть новее)", "missing": "нет в таблицах"}


def _store(a, create=False):
    root = os.path.abspath(a.root)
    if a.site:                                                # хранилище площадки (сервер) + эта рабочая копия
        from synchrotronsim import Site
        from synchrotronsim.remote import RemoteStore
        s = Site.load(a.site)
        if not s.store_url:
            raise SystemExit(f"у площадки «{a.site}» нет сервера хранилища")
        return RemoteStore(s.store_url, s.token(), working=root, site_name=s.name)
    f = find_store_file(root)
    if f is None and not create:
        raise SystemExit("хранилища нет — сначала: py -3.14 ss_store.py init")
    return SqliteStore(f, working=root) if f else SqliteStore(working=root)


def _split(spec):
    id, _, rev = str(spec).partition("@")
    return id, (rev or None)


def cmd_init(a):
    st = _store(a, create=True)
    refs = st.init_from_working(note=a.message or "импорт таблиц рабочей копии")
    print(f"хранилище: {st.path}")
    print(f"зафиксировано ревизий: {len(refs)}" + ("" if refs else " (всё уже было)"))
    for r in refs:
        print("  ", r.short)


def cmd_status(a):
    st = _store(a)
    rows = st.status([int(x) for x in a.lattice] if a.lattice else None)
    dirty = [r for r in rows if r["state"] != "clean"]
    for r in rows:
        if a.all or r["state"] != "clean":
            name = "справочник устройств" if r["kind"] == "library" else f"машина {r['id']}"
            print(f"  {name:24s} {STATE.get(r['state'], r['state'])}   последняя {str(r['head'] or '—')[:12]}")
    print("рабочая копия совпадает с хранилищем" if not dirty else f"не зафиксировано: {len(dirty)}")


def cmd_diff(a):
    st = _store(a)
    kind = "library" if a.id in ("library", LIB_ID) else "machine"
    id, rev = _split(a.id)
    to = None if a.to in (None, "working") else a.to
    print(T.format_diff(st.diff(kind, LIB_ID if kind == "library" else id, rev, to), limit=a.limit))


def cmd_commit(a):
    st = _store(a)
    for r in st.status([int(x) for x in a.lattice] if a.lattice else None):
        if r["state"] in ("modified", "new"):
            kind = r["kind"]; id = r["id"]
            print(f"— {'справочник' if kind == 'library' else 'машина ' + id}:")
            print("  " + T.format_diff(st.diff(kind, id), limit=a.limit).replace("\n", "\n  "))
    try:
        refs = st.commit_working(note=a.message, lattices=[int(x) for x in a.lattice] if a.lattice else None)
    except ConflictError as e:
        raise SystemExit(f"КОНФЛИКТ: {e}")
    print(f"зафиксировано: {len(refs)}" + ("" if refs else " — изменений нет"))
    for r in refs:
        print("  ", r.short)


def cmd_log(a):
    st = _store(a)
    kind, id = ("library", LIB_ID) if a.id in ("library", LIB_ID) else ("machine", a.id)
    for h in st.history(kind, id)[-a.n:]:
        mark = "*" if h["head"] else ("ч" if h["draft"] else " ")
        print(f" {mark} {h['rev'][:12]}  {h['created']}  {h['author'] or '—':10s} {h['note'] or ''}")
    print("(* — последняя, ч — черновик прогона)")


def cmd_checkout(a):
    st = _store(a)
    id, rev = _split(a.id)
    kind = "library" if id in ("library", LIB_ID) else "machine"
    try:
        ch = st.checkout(kind, LIB_ID if kind == "library" else id, rev, force=a.force, with_library=a.with_library)
    except StoreError as e:
        raise SystemExit(str(e))
    print("переписаны таблицы:", ", ".join(ch) if ch else "ничего (рабочая копия уже такая)")


def cmd_export(a):
    st = _store(a)
    id, rev = _split(a.id)
    if id in ("library", LIB_ID):
        sl = st.get("library", LIB_ID, None if rev in (None, "working") else rev) if rev != "working" else \
            st.working_slices([])[0]
        origin = f"library@{rev or 'head'}"
    elif rev == "working":
        sl = st.working_slices([int(id)])[1][id]
        origin = f"machine:{id}@working"
    else:
        r = rev or st.head("machine", id)
        if r is None:
            raise SystemExit(f"у машины {id} нет ревизий — сначала init/commit")
        sl = st.machine_slice_at(id, r)
        origin = f"machine:{id}@{st.revision('machine', id, r)['rev']}"
    files = T.export_slice(sl, a.folder, excel=not a.no_excel, origin=origin)
    print(f"выгружено {len(files)} таблиц в {a.folder} ({'Excel: запятая, «;», UTF-8 с BOM' if not a.no_excel else 'точка'})")
    print(f"источник: {origin}. Поправьте и загрузите: py -3.14 ss_store.py import {a.folder}")


def cmd_import(a):
    st = _store(a)
    info, base, edited = T.read_export(a.folder)
    diffs = T.diff_slices(base, edited)
    print(f"выгрузка {info.get('origin')}: правок")
    print("  " + T.format_diff(diffs, limit=a.limit).replace("\n", "\n  "))
    if not diffs:
        return
    if info["kind"] == "library":
        cur = st.working_slices([])[0]
    else:
        cur = st.working_slices([int(info["id"])])[1][str(info["id"])]
    new, conflicts = T.apply_edits(cur, base, edited)
    if conflicts:
        print(f"КОНФЛИКТЫ ({len(conflicts)}) — ячейка изменена и в рабочей копии, и в выгрузке; не применено:")
        for c in conflicts[:a.limit]:
            print(f"   {c.table} {c.key} {c.column}: база {c.base!r}, рабочая {c.ours!r}, выгрузка {c.theirs!r}")
        if not a.partial:
            raise SystemExit("ничего не записано (--partial — применить бесконфликтные правки)")
    ch = (T.write_library if info["kind"] == "library" else T.write_machine)(st.working, new)
    print("рабочая копия: переписаны", ", ".join(ch) if ch else "ничего")
    if a.commit:
        refs = st.commit_working(note=a.message or f"правки из {a.folder}",
                                 lattices=[] if info["kind"] == "library" else [int(info["id"])])
        print("зафиксировано:", ", ".join(r.short for r in refs) or "—")


def main():
    ap = argparse.ArgumentParser(description="хранилище ревизий SynchrotronSim (рабочая копия — workspace/tables)")
    ap.add_argument("--root", default=os.path.dirname(os.path.abspath(__file__)), help="корень проекта/рабочей копии")
    ap.add_argument("--site", default="", help="площадка (ss_server.py site add …): хранилище — её сервер, рабочая "
                                               "копия — этого компьютера; без ключа — локальное хранилище")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("init", help="первые ревизии из рабочей копии"); p.add_argument("-m", "--message", default="")
    p.set_defaults(f=cmd_init)
    p = sub.add_parser("status", help="изменения рабочей копии"); p.add_argument("--lattice", nargs="*")
    p.add_argument("--all", action="store_true", help="показать и неизменённые"); p.set_defaults(f=cmd_status)
    p = sub.add_parser("diff", help="разница по ячейкам"); p.add_argument("id", help="решётка[@ревизия] | library")
    p.add_argument("--to", default=None, help="ревизия (по умолчанию рабочая копия)")
    p.add_argument("--limit", type=int, default=40); p.set_defaults(f=cmd_diff)
    p = sub.add_parser("commit", help="зафиксировать изменённое"); p.add_argument("-m", "--message", required=True)
    p.add_argument("--lattice", nargs="*"); p.add_argument("--limit", type=int, default=20); p.set_defaults(f=cmd_commit)
    p = sub.add_parser("log", help="история"); p.add_argument("id"); p.add_argument("-n", type=int, default=30)
    p.set_defaults(f=cmd_log)
    p = sub.add_parser("checkout", help="ревизию — в рабочую копию"); p.add_argument("id", help="решётка[@ревизия]")
    p.add_argument("--force", action="store_true"); p.add_argument("--with-library", action="store_true")
    p.set_defaults(f=cmd_checkout)
    p = sub.add_parser("export", help="выгрузка для Excel"); p.add_argument("id", help="решётка[@ревизия|@working] | library")
    p.add_argument("folder"); p.add_argument("--no-excel", action="store_true", help="точка вместо запятой, без BOM")
    p.set_defaults(f=cmd_export)
    p = sub.add_parser("import", help="правки из выгрузки → рабочая копия"); p.add_argument("folder")
    p.add_argument("--commit", action="store_true"); p.add_argument("-m", "--message", default="")
    p.add_argument("--partial", action="store_true"); p.add_argument("--limit", type=int, default=40)
    p.set_defaults(f=cmd_import)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
