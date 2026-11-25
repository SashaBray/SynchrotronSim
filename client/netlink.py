# -*- coding: utf-8 -*-
"""NetLink — шов Link по сети: кокпит здесь, модель на вычислителе площадки (docs/DATA_PROTOCOL.md, этап 3c).

Клиентская половина `Link`: те же вызовы, что у InProcLink (subscribe, command, control, set_screen, watch,
record, poll_snapshot), только сообщения уходят по сокету (synchrotronsim.wire), а снимки приходят из потока
чтения. NetLink ПОМНИТ состояние интерфейса — подписки, отметки «● запись», настройки экранов и записи, пуск — и
повторяет его при каждом новом сеансе (`open`): пересборка модели на сервере для окон незаметна.
Движковая половина — на сервере (synchrotronsim.live: InProcLink + EngineThread).
"""
from __future__ import annotations

import threading
from typing import Any, Dict, List, Optional, Set

from .link import Link, Snapshot, Message, RUN, STOP


class NetLink(Link):
    def __init__(self, url: str, token: Optional[str] = None):
        self.url, self.token = url, token or ""
        self._sock = None
        self._send_lock = threading.Lock()
        self._snap: Optional[Snapshot] = None
        self._snap_lock = threading.Lock()
        self._subs: Set[str] = set()
        self._watches: Set[str] = set()
        self._screens: Dict[str, dict] = {}
        self._record: Optional[dict] = None
        self._running = False
        self._reader: Optional[threading.Thread] = None
        self.error: Optional[str] = None
        self.desc: Optional[dict] = None

    def __repr__(self):
        return f"NetLink({self.url}, {'на связи' if self.connected else 'нет связи'})"

    @property
    def connected(self) -> bool:
        return self._sock is not None

    # ---- сеанс
    def open(self, hello: Dict[str, Any], build_timeout: float = 1800.0, replay_run: bool = True) -> dict:
        """Новый живой сеанс на вычислителе: опыт + машина → (модель строится там) → описание сеанса. Блокирует
        до готовности модели (зовите из фонового потока). Прежний сеанс закрывается. Состояние интерфейса
        (подписки, отметки, экраны, запись, пуск) повторяется в новом сеансе; replay_run=False — пуск не
        повторять сразу (кокпит повторит его `resume()`, когда интерфейс готов — как в локальном режиме)."""
        from synchrotronsim.compute import connect, ComputeError
        from synchrotronsim import wire
        self._close_socket()
        sock = connect(self.url, self.token, mode="live", timeout=30.0, **hello)
        sock.settimeout(build_timeout)
        try:
            first = wire.recv_frame(sock)
        except Exception:
            sock.close()
            raise
        if first.get("op") != "session":
            sock.close()
            raise ComputeError(first.get("message") or f"вычислитель не собрал модель: {first}")
        sock.settimeout(None)
        self.desc = first["desc"]
        self.error = None
        with self._snap_lock:
            self._snap = None
        self._sock = sock
        self._reader = threading.Thread(target=self._read, args=(sock,), name="netlink-read", daemon=True)
        self._reader.start()
        self._replay(run=replay_run)
        return self.desc

    def resume(self) -> None:
        """Повторить пуск, если до пересборки модель бежала."""
        if self._running:
            self._send({"op": "control", "value": RUN})

    def _replay(self, run: bool = True):
        for a in sorted(self._subs):
            self._send({"op": "subscribe", "addr": a})
        for a in sorted(self._watches):
            self._send({"op": "watch", "addr": a})
        for a, v in self._screens.items():
            self._send(dict(v, op="set_screen", addr=a))
        if self._record is not None:
            self._send(dict(self._record, op="record"))
        if run and self._running:
            self._send({"op": "control", "value": RUN})

    def _read(self, sock):
        from synchrotronsim import wire
        try:
            while True:
                msg = wire.recv_frame(sock)
                op = msg.get("op")
                if op == "snapshot":
                    snap = Snapshot(values=msg.get("values") or {}, meta=msg.get("meta") or {})
                    with self._snap_lock:
                        self._snap = snap
                elif op == "error":
                    self.error = msg.get("message")
        except (ConnectionError, OSError, ValueError) as e:
            if self._sock is sock:
                self.error = self.error or f"связь с вычислителем потеряна: {e}"
                self._sock = None

    def _send(self, msg: Dict[str, Any]) -> None:
        sock = self._sock
        if sock is None:
            return                                       # нет сеанса — состояние запомнено, повторится в open()
        from synchrotronsim import wire
        try:
            with self._send_lock:
                wire.send_frame(sock, msg)
        except OSError as e:
            self.error = f"связь с вычислителем потеряна: {e}"
            self._sock = None

    def _close_socket(self):
        sock, self._sock = self._sock, None
        if sock is not None:
            try:
                from synchrotronsim import wire
                with self._send_lock:
                    wire.send_frame(sock, {"op": "close"})
            except OSError:
                pass
            try:
                sock.close()
            except OSError:
                pass

    def close(self) -> None:
        self._close_socket()

    # ---- клиентская сторона Link
    def subscribe(self, addr: str) -> None:
        self._subs.add(addr)
        self._send({"op": "subscribe", "addr": addr})

    def unsubscribe(self, addr: str) -> None:
        self._subs.discard(addr)
        self._send({"op": "unsubscribe", "addr": addr})

    def command(self, addr: str, value: Any) -> None:
        self._send({"op": "command", "addr": addr, "value": value})

    def control(self, action: str) -> None:
        if action == RUN:
            self._running = True
        elif action == STOP:
            self._running = False
        self._send({"op": "control", "value": action})

    def set_screen(self, addr: str, area: Any, res: Any, point: Any = None, quantity: Any = None) -> None:
        v = {"area": area, "res": res, "point": point, "quantity": quantity}
        self._screens[addr] = v
        self._send(dict(v, op="set_screen", addr=addr))

    def watch(self, addr: str) -> None:
        self._watches.add(addr)
        self._send({"op": "watch", "addr": addr})

    def unwatch(self, addr: str) -> None:
        self._watches.discard(addr)
        self._send({"op": "unwatch", "addr": addr})

    def record(self, on: bool, max_turns: Optional[int] = None, repeat: bool = False) -> None:
        self._record = {"on": bool(on), "max_turns": int(max_turns) if max_turns else None, "repeat": bool(repeat)}
        self._send(dict(self._record, op="record"))

    def poll_snapshot(self) -> Optional[Snapshot]:
        with self._snap_lock:
            return self._snap

    def remember_running(self, on: bool) -> None:
        """Пуск/пауза для следующего сеанса (пересборка): как было в прежнем."""
        self._running = bool(on)

    # ---- движковая сторона — на сервере
    def pull(self) -> List[Message]:
        raise NotImplementedError("NetLink — клиентская половина; движок на вычислителе")

    def subscriptions(self) -> Set[str]:
        return set(self._subs)

    def publish(self, snapshot: Snapshot) -> None:
        raise NotImplementedError("NetLink — клиентская половина; движок на вычислителе")
