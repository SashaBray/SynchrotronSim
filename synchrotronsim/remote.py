# -*- coding: utf-8 -*-
"""Сервер хранилища и клиент к нему (docs/DATA_PROTOCOL.md, этап 3) — HTTP из стандартной библиотеки.

    py -3.14 ss_server.py store --host 0.0.0.0 --port 8701 --token-file ~/.synchrotronsim/lab.token

Сервер отдаёт по сети ПРОСТЫЕ операции хранилища (synchrotronsim.working): ревизии, их содержимое, запись с
проверкой родителя, блобы, прогоны. Рабочая копия, статус, разница, сборка `Machine` из среза — у клиента
(`RemoteStore` — тот же `WorkingCopyStore`, что и локальный `SqliteStore`).

HTTP: GET /v1/ping (без токена: протокол, служба), POST /v1/rpc — тело в формате synchrotronsim.wire
{"op", "args", "kw"} → {"ok": true, "result"} | {"ok": false, "error": имя, "message"}. Токен — заголовок
`Authorization: Bearer …` (не в адресе). Слушать не только 127.0.0.1 без токена сервер отказывается.
"""
from __future__ import annotations

import hmac
import json
import os
import platform
import socket
import threading
import time
import urllib.error
import urllib.request
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from typing import Any, Dict, List, Optional

from . import wire
from .store import StoreError, ConflictError
from .working import WorkingCopyStore

STORE_OPS = ("head", "revision", "history", "list", "library_content", "machine_payload", "experiment_content",
             "put_library", "put_machine_slice", "put_experiment", "put_machine_doc", "find_machine_fp", "blob_put",
             "blob_get", "put_run", "run_bundle", "list_runs", "put_json", "json_content")
_ERRORS = {"ConflictError": ConflictError, "StoreError": StoreError}


def _loopback(host: str) -> bool:
    return host in ("127.0.0.1", "localhost", "::1")


def server_info(service: str, extra: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """Сведения о сервере для Site.check(): протокол, ОС, Python, сборка ядра, коммит кода, ядра, время."""
    info = {"service": service, "protocol": wire.PROTOCOL, "host": socket.gethostname(),
            "os": f"{platform.system()} {platform.release()}", "machine": platform.machine(),
            "python": platform.python_version(), "cpus": os.cpu_count(), "time": time.time()}
    try:
        from .run import provenance
        from .db import PROJECT_ROOT
        pv = provenance(PROJECT_ROOT) or {}
        info["engine"] = (pv.get("engine") or {}).get("extension")   # build_info + sha16 файла расширения
        info["code"] = pv.get("git")
        info["numpy"] = pv.get("numpy")
    except Exception as e:                                     # noqa: BLE001
        info["engine"], info["code"] = None, f"нет: {type(e).__name__}"
    info.update(extra or {})
    return info


# ---------------------------------------------------------------------------------------------- сервер
class _StoreHandler(BaseHTTPRequestHandler):
    server_version = "SynchrotronSimStore/1"

    def log_message(self, fmt, *args):                       # тихо; журнал — в server.log_fn
        fn = getattr(self.server, "log_fn", None)
        if fn:
            fn(f"{self.address_string()} {fmt % args}")

    def _send(self, code: int, body: bytes, ctype: str = "application/x-synchrotronsim"):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        if self.path.rstrip("/") == "/v1/ping":
            body = json.dumps({"service": "store", "protocol": wire.PROTOCOL}).encode()
            return self._send(200, body, "application/json")
        self._send(404, b"not found", "text/plain")

    def _authorized(self) -> bool:
        tok = self.server.token
        if not tok:
            return True
        got = self.headers.get("Authorization", "")
        return got.startswith("Bearer ") and hmac.compare_digest(got[7:].strip(), tok)

    def do_POST(self):
        if self.path.rstrip("/") != "/v1/rpc":
            return self._send(404, b"not found", "text/plain")
        if not self._authorized():
            return self._send(401, wire.encode({"ok": False, "error": "AuthError", "message": "нет доступа (токен)"}))
        try:
            n = int(self.headers.get("Content-Length", "0"))
            msg = wire.decode(self.rfile.read(n))
            op = msg.get("op")
            if op == "info":
                res = server_info("store", self.server.store_info())
            elif op in STORE_OPS:
                res = getattr(self.server.store, op)(*msg.get("args", []), **msg.get("kw", {}))
            else:
                raise StoreError(f"неизвестная операция {op!r}")
            out = {"ok": True, "result": res}
        except Exception as e:                                  # noqa: BLE001 — ошибку отдаём клиенту
            out = {"ok": False, "error": type(e).__name__, "message": str(e)}
        self._send(200, wire.encode(out))


class StoreServer(ThreadingHTTPServer):
    """HTTP-сервер над хранилищем (обычно SqliteStore без рабочей копии)."""
    daemon_threads = True

    def __init__(self, store, host: str = "127.0.0.1", port: int = 8701, token: Optional[str] = None, log_fn=None):
        if not token and not _loopback(host):
            raise ValueError(f"сервер хранилища на {host} без токена — откажусь: в сети нужен токен (--token-file)")
        super().__init__((host, int(port)), _StoreHandler)
        self.store, self.token, self.log_fn = store, token or "", log_fn

    def store_info(self) -> Dict[str, Any]:
        st = self.store
        out = {"store": getattr(st, "path", repr(st))}
        try:
            out["machines"] = len(st.list("machine"))
            out["schema"] = 1
        except Exception:
            pass
        return out

    @property
    def url(self) -> str:
        h, p = self.server_address[:2]
        return f"http://{h}:{p}"

    def start(self) -> "StoreServer":
        threading.Thread(target=self.serve_forever, name="store-server", daemon=True).start()
        return self


# ---------------------------------------------------------------------------------------------- клиент
class RemoteStore(WorkingCopyStore):
    """Хранилище на сервере (те же операции, что у SqliteStore) + своя рабочая копия таблиц на этом компьютере.
    Базы рабочей копии — `workspace/store/working-<площадка>.json` (у каждой площадки свои)."""

    def __init__(self, url: str, token: Optional[str] = None, working: Optional[str] = None, site_name: str = "remote",
                 timeout: float = 120.0):
        self.url = url.rstrip("/")
        self._token = token or ""
        if self._token and not self._token.isascii():
            raise StoreError("токен площадки — только латинские буквы, цифры и знаки (ss_server.py token)")
        self.timeout = float(timeout)
        self.working = os.path.abspath(working) if working else None
        self.site = site_name
        ws = None
        if self.working:
            ws = self.working if os.path.isdir(os.path.join(self.working, "tables")) else \
                os.path.join(self.working, "workspace")
        self.workspace = ws
        if ws:
            self.state_dir = os.path.join(ws, "store")
        else:
            import tempfile
            self.state_dir = os.path.join(tempfile.gettempdir(), "synchrotronsim_remote")
        self.state_name = f"working-{site_name}.json"
        self._init_cache()
        self._lib: Dict[str, Any] = {}
        self._payload: Dict[tuple, Any] = {}

    def __repr__(self):
        return f"RemoteStore({self.url}, площадка {self.site}, рабочая копия {self.working or 'нет'})"

    def _call(self, op: str, *args, **kw):
        req = urllib.request.Request(self.url + "/v1/rpc", data=wire.encode({"op": op, "args": list(args), "kw": kw}),
                                     method="POST", headers={"Content-Type": "application/x-synchrotronsim"})
        if self._token:
            req.add_header("Authorization", f"Bearer {self._token}")
        try:
            with urllib.request.urlopen(req, timeout=self.timeout) as r:
                out = wire.decode(r.read())
        except urllib.error.HTTPError as e:
            body = e.read()
            try:
                out = wire.decode(body)
            except Exception:
                raise StoreError(f"сервер хранилища {self.url}: HTTP {e.code}") from None
        except (urllib.error.URLError, OSError) as e:
            raise StoreError(f"сервер хранилища {self.url} недоступен: {e}") from None
        if not out.get("ok"):
            cls = _ERRORS.get(out.get("error"), StoreError)
            raise cls(f"{out.get('message')}" + ("" if cls is not StoreError else f" [{out.get('error')}]"))
        return out.get("result")

    def info(self) -> Dict[str, Any]:
        return self._call("info")

    # ---- простые операции (кэш — только у неизменяемого по ревизии содержимого)
    def head(self, kind, id):
        return self._call("head", kind, str(id))

    def revision(self, kind, id, rev):
        return self._call("revision", kind, str(id), rev)

    def history(self, kind, id):
        return self._call("history", kind, str(id))

    def list(self, kind):
        return self._call("list", kind)

    def library_content(self, rev):
        if rev not in self._lib:
            self._lib[rev] = self._call("library_content", rev)
        return self._lib[rev]

    def machine_payload(self, id, rev):
        key = (str(id), rev)
        if key not in self._payload:
            p = self._call("machine_payload", str(id), rev)
            self._payload[(str(id), p["rev"])] = p                 # по точной ревизии — неизменяемо
            if p["rev"] != rev:
                return p
        return self._payload[key]

    def experiment_content(self, id, rev):
        return self._call("experiment_content", str(id), rev)

    def put_library(self, lib, parent, note, draft=False):
        return self._call("put_library", lib, parent, note, draft)

    def put_machine_slice(self, id, sl, parent, note, draft=False):
        return self._call("put_machine_slice", str(id), sl, parent, note, draft)

    def put_experiment(self, id, d, parent, note, draft=False):
        return self._call("put_experiment", str(id), d, parent, note, draft)

    def put_machine_doc(self, id, doc, note):
        return self._call("put_machine_doc", str(id), doc, note)

    def find_machine_fp(self, id, fp):
        return self._call("find_machine_fp", str(id), fp)

    def blob_put(self, data: bytes) -> str:
        return self._call("blob_put", data)

    def blob_get(self, h: str) -> bytes:
        return self._call("blob_get", h)

    def put_json(self, kind, id, d, parent=None, note=""):
        return self._call("put_json", kind, str(id), d, parent, note)

    def json_content(self, kind, id, rev=None):
        return self._call("json_content", kind, str(id), rev)

    # ---- прогоны
    def put_run(self, run_id, bundle, manifest):
        return self._call("put_run", run_id, bundle, manifest)

    def list_runs(self) -> List[Dict[str, Any]]:
        return self._call("list_runs")

    def get_run(self, id):
        from .store_sqlite import unpack_run
        import tempfile
        root = self.working or os.path.join(tempfile.gettempdir(), "synchrotronsim_remote")
        if not self.working:
            os.makedirs(os.path.join(root, "tables"), exist_ok=True)
        return unpack_run(self._call("run_bundle", id), root, id)


def ping(url: str, timeout: float = 5.0) -> Dict[str, Any]:
    """GET /v1/ping (без токена): {"service", "protocol"} — проверка, что по адресу наш сервер."""
    with urllib.request.urlopen(url.rstrip("/") + "/v1/ping", timeout=timeout) as r:
        return json.loads(r.read().decode("utf-8"))
