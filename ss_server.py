# -*- coding: utf-8 -*-
"""Серверы SynchrotronSim и площадки (docs/DATA_PROTOCOL.md, этап 3) — только стандартная библиотека Python.

На компьютере с базой (хранилище):
    py -3.14 ss_server.py token --file lab.token                       # токен площадки (один раз)
    py -3.14 ss_server.py store --host 0.0.0.0 --port 8701 --token-file lab.token
На вычислительном сервере:
    py -3.14 ss_server.py compute --host 0.0.0.0 --port 8702 --token-file lab.token --store http://db-pc:8701
На своём компьютере (библиотека, клиент):
    py -3.14 ss_server.py site add lab --store http://db-pc:8701 --compute tcp://calc-pc:8702 --token-file lab.token
    py -3.14 ss_server.py site check lab                              # связь, задержка, версии, сборка ядра
    py -3.14 ss_store.py --site lab status                            # рабочая копия против хранилища площадки

Профили — ~/.synchrotronsim/sites.json (адреса), токены — ~/.synchrotronsim/credentials.json (только локально;
в хранилище и в git не попадают). Без токена сервер слушает только 127.0.0.1.
"""
import argparse
import os
import secrets
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.stdout.reconfigure(encoding="utf-8")

ROOT = os.path.dirname(os.path.abspath(__file__))


def _read_token(a):
    if getattr(a, "token_file", None):
        with open(os.path.expanduser(a.token_file), encoding="ascii") as f:
            return f.read().strip()
    return os.environ.get(getattr(a, "token_env", None) or "SS_TOKEN") or None


def cmd_token(a):
    tok = secrets.token_urlsafe(32)
    if a.file:
        p = os.path.expanduser(a.file)
        with open(p, "w", encoding="ascii") as f:
            f.write(tok + "\n")
        try:
            os.chmod(p, 0o600)
        except OSError:
            pass
        print(f"токен записан в {p} — передайте файл на серверы и клиенты площадки защищённым способом")
    if a.save:
        from synchrotronsim.site import Site
        print("токен сохранён для площадки", a.save, "в", Site.set_token(a.save, tok))
    if not a.file and not a.save:
        print(tok)


def cmd_store(a):
    from synchrotronsim.store_sqlite import SqliteStore, find_store_file
    from synchrotronsim.remote import StoreServer
    path = a.path or find_store_file(ROOT)
    if path is None:
        raise SystemExit("нет файла хранилища — создайте: py -3.14 ss_store.py init (или --path)")
    st = SqliteStore(path, working=None)                      # сервер — без рабочей копии
    tok = _read_token(a)
    srv = StoreServer(st, a.host, a.port, token=tok, log_fn=print if a.verbose else None)
    print(f"сервер хранилища {srv.url} → {path}; токен {'есть' if tok else 'нет (только этот компьютер)'}")
    print("остановка — Ctrl+C")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.server_close()


def cmd_compute(a):
    from synchrotronsim.compute import ComputeServer
    tok = _read_token(a)
    srv = ComputeServer(a.host, a.port, token=tok, store_url=a.store, store_token=tok, root=a.root or ROOT,
                        log_fn=print)
    print(f"вычислитель tcp://{srv.host}:{srv.port}; хранилище {a.store or 'локальное ' + (a.root or ROOT)}; "
          f"токен {'есть' if tok else 'нет (только этот компьютер)'}")
    print("остановка — Ctrl+C")
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        srv.close()


def cmd_site(a):
    from synchrotronsim.site import Site
    if a.action == "add":
        s = Site(a.name, store_url=a.store or "", compute_url=a.compute or "", note=a.note or "")
        print("профиль сохранён в", s.save())
        tok = _read_token(a)
        if tok:
            print("токен сохранён в", Site.set_token(a.name, tok))
    elif a.action == "list":
        sites = Site.list_sites()
        if not sites:
            print("площадок нет (ss_server.py site add …)")
        for n, d in sites.items():
            print(f"  {n:12s} хранилище {d.get('store_url') or 'локальное'} · вычислитель {d.get('compute_url') or 'в процессе'}"
                  + (f" · {d['note']}" if d.get("note") else ""))
    elif a.action == "check":
        t0 = time.time()
        rep = Site.load(a.name).check()
        print(rep)
        print(f"(проверка {time.time() - t0:.1f} с)")
        if not rep.ok:
            raise SystemExit(1)
    elif a.action == "remove":
        Site.remove(a.name)
        print("площадка удалена из профилей:", a.name)
    elif a.action == "publish":
        s = Site.load(a.name)
        print("профиль (без токена) в хранилище площадки:", s.publish(s.store)[:12])


def main():
    ap = argparse.ArgumentParser(description="серверы SynchrotronSim и площадки")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("token", help="новый токен площадки")
    p.add_argument("--file", help="записать в файл (для серверов)"); p.add_argument("--save", help="сохранить для площадки")
    p.set_defaults(f=cmd_token)
    for name, port, fn, hlp in (("store", 8701, cmd_store, "сервер хранилища"), ("compute", 8702, cmd_compute, "вычислитель")):
        p = sub.add_parser(name, help=hlp)
        p.add_argument("--host", default="127.0.0.1", help="адрес (0.0.0.0 — все сети; нужен токен)")
        p.add_argument("--port", type=int, default=port)
        p.add_argument("--token-file"); p.add_argument("--token-env", default="SS_TOKEN")
        p.add_argument("--verbose", action="store_true")
        if name == "store":
            p.add_argument("--path", help="файл хранилища (по умолчанию workspace/store/store.sqlite)")
        else:
            p.add_argument("--store", default="", help="адрес сервера хранилища (http://…); пусто — локальное")
            p.add_argument("--root", default="", help="корень проекта на вычислителе (кэш, прогоны)")
        p.set_defaults(f=fn)
    p = sub.add_parser("site", help="площадки: add | list | check | remove | publish")
    p.add_argument("action", choices=["add", "list", "check", "remove", "publish"])
    p.add_argument("name", nargs="?", default="")
    p.add_argument("--store", help="http://хост:8701"); p.add_argument("--compute", help="tcp://хост:8702")
    p.add_argument("--note"); p.add_argument("--token-file"); p.add_argument("--token-env", default="SS_TOKEN_NEW")
    p.set_defaults(f=cmd_site)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
