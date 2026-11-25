# -*- coding: utf-8 -*-
"""Вычислитель площадки (docs/DATA_PROTOCOL.md, этап 3, раздел 6) — TCP из стандартной библиотеки.

    py -3.14 ss_server.py compute --host 0.0.0.0 --port 8702 --token-file lab.token --store http://db-pc:8701

Клиент присылает опыт и ССЫЛКУ на ревизию машины (`machine:26@…`); вычислитель сам берёт машину из хранилища
площадки (кэш по ревизии), считает (`Simulation`) и кладёт прогон в то же хранилище; клиент получает прогон по
номеру. Нет сервера хранилища у площадки — клиент присылает документ машины целиком.

Поток — кадры synchrotronsim.wire. Начало соединения: {"op": "hello", "token", "protocol"} → {"ok"}. Дальше
запрос → ответ: ping, info, submit, job, jobs, cancel. Задания выполняются по одному (ядро и так занимает все
потоки), остальные ждут в очереди. Живой сеанс кокпита (этап 3c) — отдельное соединение «live».
"""
from __future__ import annotations

import hmac
import itertools
import os
import queue
import socket
import socketserver
import threading
import time
import traceback
from typing import Any, Callable, Dict, List, Optional

from . import wire

STATES = ("queued", "running", "done", "failed", "cancelled")


class ComputeError(RuntimeError):
    pass


def _loopback(host: str) -> bool:
    return host in ("127.0.0.1", "localhost", "::1")


# ---------------------------------------------------------------------------------------------- сервер
class _Job:
    def __init__(self, jid: str, exp: dict, machine_ref: Optional[str], machine_doc: Optional[bytes],
                 steps: Optional[int], note: str):
        self.id, self.exp, self.machine_ref, self.machine_doc, self.steps, self.note = \
            jid, exp, machine_ref, machine_doc, steps, note
        self.state = "queued"
        self.progress: List[str] = []
        self.run_id: Optional[str] = None
        self.error: Optional[str] = None
        self.created = time.time()
        self.started: Optional[float] = None
        self.finished: Optional[float] = None

    def public(self) -> Dict[str, Any]:
        return {"job": self.id, "state": self.state, "progress": list(self.progress[-200:]), "run_id": self.run_id,
                "error": self.error, "created": self.created, "started": self.started, "finished": self.finished,
                "experiment": self.exp.get("name"), "machine": self.machine_ref}


class _Handler(socketserver.BaseRequestHandler):
    def handle(self):
        srv: "ComputeServer" = self.server.owner
        sock = self.request
        try:
            hello = wire.recv_frame(sock)
        except (ConnectionError, ValueError):
            return
        if hello.get("op") != "hello" or not srv._auth(hello.get("token")):
            wire.send_frame(sock, {"ok": False, "error": "AuthError", "message": "нет доступа (токен)"})
            return
        if hello.get("protocol") != wire.PROTOCOL:
            wire.send_frame(sock, {"ok": False, "error": "ProtocolError",
                                   "message": f"протокол {hello.get('protocol')} ≠ {wire.PROTOCOL}"})
            return
        wire.send_frame(sock, {"ok": True, "protocol": wire.PROTOCOL})
        if hello.get("mode") == "live":                         # этап 3c: живой сеанс кокпита
            from .live import serve_live
            return serve_live(srv, sock, hello)
        while True:
            try:
                msg = wire.recv_frame(sock)
            except (ConnectionError, OSError, ValueError):
                return
            try:
                res = srv.dispatch(msg)
                out = {"ok": True, "result": res}
            except Exception as e:                              # noqa: BLE001 — ошибку отдаём клиенту
                out = {"ok": False, "error": type(e).__name__, "message": str(e)}
            try:
                wire.send_frame(sock, out)
            except OSError:
                return


class _TCP(socketserver.ThreadingTCPServer):
    daemon_threads = True
    allow_reuse_address = True


class ComputeServer:
    """Сервер вычислений: очередь заданий Simulation, хранилище площадки (сервер или локальное), прогоны — туда."""

    def __init__(self, host: str = "127.0.0.1", port: int = 8702, token: Optional[str] = None,
                 store_url: str = "", store_token: Optional[str] = None, root: Optional[str] = None,
                 log_fn: Optional[Callable[[str], None]] = None, store=None):
        if not token and not _loopback(host):
            raise ValueError(f"вычислитель на {host} без токена — откажусь: в сети нужен токен (--token-file)")
        from .db import PROJECT_ROOT
        self.token = token or ""
        self.root = os.path.abspath(root or PROJECT_ROOT)
        self.store_url = store_url or ""
        self.log = log_fn or (lambda s: None)
        if store is not None:
            self.store = store
        elif self.store_url:
            from .remote import RemoteStore
            self.store = RemoteStore(self.store_url, store_token or token, working=None, site_name="compute")
        else:
            from .store import open_store
            self.store = open_store(root=self.root)
        self._jobs: Dict[str, _Job] = {}
        self._q: "queue.Queue[_Job]" = queue.Queue()
        self._ids = itertools.count(1)
        self._lock = threading.Lock()
        self._tcp = _TCP((host, int(port)), _Handler)
        self._tcp.owner = self
        self.host, self.port = self._tcp.server_address[:2]
        self._worker = threading.Thread(target=self._work, name="compute-worker", daemon=True)
        self._worker.start()
        self.live_sessions = 0

    @property
    def url(self) -> str:
        return f"tcp://{self.host}:{self.port}"

    def _auth(self, tok) -> bool:
        return not self.token or (isinstance(tok, str) and hmac.compare_digest(tok, self.token))

    def serve_forever(self):
        self._tcp.serve_forever()

    def start(self) -> "ComputeServer":
        threading.Thread(target=self.serve_forever, name="compute-server", daemon=True).start()
        return self

    def close(self):
        self._tcp.shutdown()
        self._tcp.server_close()

    # ---- запросы
    def dispatch(self, msg: Dict[str, Any]):
        op = msg.get("op")
        if op == "ping":
            return {"t": time.time()}
        if op == "info":
            from .remote import server_info
            return server_info("compute", {"store_url": self.store_url, "queue": self._q.qsize(),
                                           "jobs": len(self._jobs), "root": self.root, "live": self.live_sessions})
        if op == "submit":
            with self._lock:
                jid = f"job-{time.strftime('%Y%m%d-%H%M%S')}-{next(self._ids)}"
                job = _Job(jid, msg["experiment"], msg.get("machine_ref"), msg.get("machine_doc"), msg.get("steps"),
                           msg.get("note", ""))
                self._jobs[jid] = job
            self._q.put(job)
            self.log(f"задание {jid}: опыт «{job.exp.get('name')}», машина {job.machine_ref or 'документ'}")
            return {"job": jid}
        if op == "job":
            j = self._jobs.get(msg.get("job"))
            if j is None:
                raise ComputeError(f"нет задания {msg.get('job')}")
            return j.public()
        if op == "jobs":
            return [j.public() for j in self._jobs.values()]
        if op == "cancel":
            j = self._jobs.get(msg.get("job"))
            if j is not None and j.state == "queued":
                j.state = "cancelled"
            return {"state": j.state if j else None}
        if op == "run_bundle":                                  # прогон с диска вычислителя (нет сервера хранилища)
            rid = str(msg.get("run_id"))
            if hasattr(self.store, "run_bundle"):
                try:
                    return self.store.run_bundle(rid)
                except Exception:                               # noqa: BLE001 — нет в хранилище — с диска
                    pass
            from .run import Run
            from .store_sqlite import pack_run
            return pack_run(Run.open(rid, root=self.root))
        raise ComputeError(f"неизвестная операция {op!r}")

    # ---- исполнение
    def live_session(self, exp, machine):
        """Модель живого сеанса кокпита: (plant, state, sch). Подменяется в тестах (лёгкая модель)."""
        from .simulation import build_session
        S = build_session(exp, machine=machine, store=self.store)
        return S.plant, S.state, S.sch

    def machine_for(self, exp, machine_ref: Optional[str], machine_doc: Optional[bytes]):
        from .store import Ref
        if machine_doc is not None:
            from .documents import machine_from_bytes
            return machine_from_bytes(machine_doc, source="document:прислан клиентом")
        if machine_ref:
            return self.store.resolve(Ref.parse(machine_ref))
        return self.store.get("machine", str(exp.lattice_id), exp.machine_rev or "head")

    def _work(self):
        while True:
            job = self._q.get()
            if job.state == "cancelled":
                continue
            job.state, job.started = "running", time.time()
            try:
                self._execute(job)
                job.state = "done"
            except Exception as e:                              # noqa: BLE001
                job.state, job.error = "failed", f"{type(e).__name__}: {e}"
                self.log(f"задание {job.id}: ОШИБКА {job.error}\n{traceback.format_exc()}")
            job.finished = time.time()

    def _execute(self, job: _Job):
        from .experiment import Experiment
        from .simulation import Simulation
        from .store_sqlite import pack_run
        exp = Experiment.from_dict(job.exp)
        m = self.machine_for(exp, job.machine_ref, job.machine_doc)
        job.progress.append(f"машина {m.lattice_id}: {m.fingerprint[:12]} ({m.source})")
        sim = Simulation(exp, root=self.root, legacy_fb=False, machine=m, store=self.store,
                         progress=job.progress.append)
        run = sim.run(steps=job.steps)
        job.run_id = run.run_id
        if hasattr(self.store, "put_run"):
            bundle = pack_run(run)
            job.progress.append(f"прогон {run.run_id}: {len(bundle) // 1024} КБ → хранилище")
            self.store.put_run(run.run_id, bundle, run.manifest)
        else:
            job.progress.append(f"прогон {run.run_id} — на диске вычислителя (у хранилища нет приёма прогонов)")
        self.log(f"задание {job.id}: готово, прогон {run.run_id}")


# ---------------------------------------------------------------------------------------------- клиент
def _parse_tcp(url: str):
    u = url.strip()
    if u.startswith("tcp://"):
        u = u[6:]
    host, _, port = u.rpartition(":")
    return host or "127.0.0.1", int(port or 8702)


def connect(url: str, token: Optional[str], mode: str = "rpc", timeout: float = 30.0, **hello):
    """Соединение с вычислителем после рукопожатия (токен, версия протокола)."""
    host, port = _parse_tcp(url)
    try:
        sock = socket.create_connection((host, port), timeout=timeout)
    except OSError as e:
        raise ComputeError(f"вычислитель {url} недоступен: {e}") from None
    wire.send_frame(sock, dict(hello, op="hello", token=token or "", protocol=wire.PROTOCOL, mode=mode))
    ans = wire.recv_frame(sock)
    if not ans.get("ok"):
        sock.close()
        raise ComputeError(ans.get("message") or "отказ вычислителя")
    return sock


class RemoteCompute:
    """Клиент вычислителя площадки: задания Simulation, ожидание, прогон — из хранилища площадки."""

    def __init__(self, url: str, token: Optional[str] = None, site=None, timeout: float = 60.0):
        self.url, self._token, self.site, self.timeout = url, token or "", site, float(timeout)
        self._sock = None
        self._lock = threading.Lock()

    def __repr__(self):
        return f"RemoteCompute({self.url})"

    def _call(self, op: str, **kw):
        with self._lock:
            for attempt in (1, 2):                              # соединение могло закрыться — одна повторная попытка
                if self._sock is None:
                    self._sock = connect(self.url, self._token, timeout=self.timeout)
                try:
                    wire.send_frame(self._sock, dict(kw, op=op))
                    out = wire.recv_frame(self._sock)
                    break
                except (ConnectionError, OSError):
                    try:
                        self._sock.close()
                    except OSError:
                        pass
                    self._sock = None
                    if attempt == 2:
                        raise ComputeError(f"вычислитель {self.url}: соединение потеряно") from None
        if not out.get("ok"):
            raise ComputeError(f"{out.get('message')} [{out.get('error')}]")
        return out.get("result")

    def close(self):
        if self._sock is not None:
            try:
                self._sock.close()
            finally:
                self._sock = None

    def ping(self):
        return self._call("ping")

    def info(self):
        return self._call("info")

    def jobs(self):
        return self._call("jobs")

    def job(self, jid: str):
        return self._call("job", job=jid)

    def submit(self, exp, steps: Optional[int] = None, machine=None, note: str = "") -> str:
        """Задание: опыт + ссылка на ревизию машины. Машина берётся из хранилища площадки: рабочая копия этого ПК
        фиксируется там (черновик, если не зафиксирована), вычислителю уходит только ссылка. Нет сервера
        хранилища у площадки — уходит документ машины целиком."""
        store = getattr(self.site, "store", None) if self.site is not None and getattr(self.site, "store_url", "") else None
        ref, doc = None, None
        if store is not None:
            m = machine if machine is not None else store.get("machine", str(exp.lattice_id), exp.machine_rev or None)
            ref = str(store.pin(m))
        else:
            from .simulation import resolve_machine
            from .documents import machine_to_bytes
            m = machine if machine is not None else resolve_machine(exp)
            doc = machine_to_bytes(m)
        return self._call("submit", experiment=exp.to_dict(), machine_ref=ref, machine_doc=doc, steps=steps,
                          note=note)["job"]

    def wait(self, jid: str, progress: Optional[Callable[[str], None]] = print, poll: float = 1.0,
             timeout: Optional[float] = None) -> Dict[str, Any]:
        seen, t0 = 0, time.time()
        while True:
            j = self.job(jid)
            lines = j.get("progress") or []
            if progress:
                for s in lines[seen:]:
                    progress(s)
            seen = len(lines)
            if j["state"] in ("done", "failed", "cancelled"):
                return j
            if timeout and time.time() - t0 > timeout:
                raise ComputeError(f"задание {jid}: не завершилось за {timeout} с")
            time.sleep(poll)

    def run(self, exp, steps: Optional[int] = None, machine=None, progress: Optional[Callable[[str], None]] = print):
        """Задание → ожидание → прогон (synchrotronsim.run.Run, локальная копия из хранилища площадки)."""
        jid = self.submit(exp, steps=steps, machine=machine)
        j = self.wait(jid, progress=progress)
        if j["state"] != "done":
            raise ComputeError(f"задание {jid}: {j['state']} — {j.get('error')}")
        if self.site is not None and getattr(self.site, "store_url", ""):
            return self.site.store.get("run", j["run_id"])       # из хранилища площадки
        from .store_sqlite import unpack_run
        from .db import PROJECT_ROOT
        return unpack_run(self._call("run_bundle", run_id=j["run_id"]), PROJECT_ROOT, j["run_id"])
