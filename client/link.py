# -*- coding: utf-8 -*-
"""Шов Link — граница между потоком-движком и UI (она же будущая клиент-серверная граница).

Движок крутит модель с полной скоростью и владеет доской; UI в главном потоке Qt читает снимок и
рисует. Между ними — Link. Через него ходят ТОЛЬКО сериализуемые данные (числа, строки, их списки/словари):
никаких живых объектов доски. Это и позволит позже подменить InProcLink на NetLink (сокет, этап 7) не
трогая ни UI, ни движок — обе стороны видят лишь контракт Link.

Направления (эскиз из docs/CLIENT_PLAN.md §1):
  клиент→движок:  subscribe(addr) / unsubscribe(addr) / command(addr, value)
                  control("run"|"stop"|"restart"|"shutdown") / set_screen(addr, area, res, point, quantity)
  движок→клиент:  poll_snapshot() -> Snapshot{ values:{addr:val}, meta:{model_time, wall_time,
                                               speed, cycle, seq} }   # только подписанное + мета

Сообщения клиент→движок едут одной упорядоченной очередью как типизированные Message (так их легко
сериализовать в NetLink). Снимок движок→клиент — один слот «последний побеждает» (двойной буфер):
UI всегда берёт самый свежий кадр, промежуточные не копятся.
"""
from __future__ import annotations

import queue
import threading
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Set

# действия канала control (строки — чтобы ехали по проводу как есть)
RUN = "run"
STOP = "stop"
RESTART = "restart"
SHUTDOWN = "shutdown"

# типы сообщений клиент→движок
SUBSCRIBE = "subscribe"
UNSUBSCRIBE = "unsubscribe"
COMMAND = "command"
CONTROL = "control"
SET_SCREEN = "set_screen"


@dataclass(frozen=True)
class Message:
    """Типизированное сообщение клиент→движок (сериализуемо: строки + сериализуемое value).
      subscribe/unsubscribe: addr — адрес сигнала.
      command:   addr — адрес ручки, value — новое значение (число/строка/bool).
      control:   value — действие ("run"/"stop"/"restart"/"shutdown"), addr пуст.
      set_screen: addr — адрес экрана, value — {"area", "res", "point", "quantity"} (медленный такт);
                 quantity — что рисовать: "flux" (поток на E₁, интерференция) | "power" (мощность, Вт/мм²)."""
    kind: str
    addr: str = ""
    value: Any = None


@dataclass
class Snapshot:
    """Кадр движок→клиент: значения ТОЛЬКО подписанных адресов + мета прогона. Всё сериализуемо."""
    values: Dict[str, Any] = field(default_factory=dict)
    meta: Dict[str, Any] = field(default_factory=dict)


class Link(ABC):
    """Абстрактный транспорт. Клиентская сторона (UI) и движковая сторона зовут разные методы;
    в InProcLink это один объект, в NetLink — две половины по разные стороны сокета."""

    # ---- клиентская сторона (UI) ----
    @abstractmethod
    def subscribe(self, addr: str) -> None: ...

    @abstractmethod
    def unsubscribe(self, addr: str) -> None: ...

    @abstractmethod
    def command(self, addr: str, value: Any) -> None: ...

    @abstractmethod
    def control(self, action: str) -> None: ...

    @abstractmethod
    def set_screen(self, addr: str, area: Any, res: Any, point: Any = None, quantity: Any = None) -> None: ...

    @abstractmethod
    def poll_snapshot(self) -> Optional[Snapshot]: ...

    # ---- движковая сторона ----
    @abstractmethod
    def pull(self) -> List[Message]: ...

    @abstractmethod
    def subscriptions(self) -> Set[str]: ...

    @abstractmethod
    def publish(self, snapshot: Snapshot) -> None: ...

    def close(self) -> None:  # освобождение ресурсов транспорта (для NetLink — закрыть сокет)
        pass


class InProcLink(Link):
    """Внутрипроцессный Link на потоках: очередь сообщений (клиент→движок) + однослотовый буфер снимка
    (движок→клиент). Множество подписок ведётся транспортом: движок зовёт pull() (осушает очередь,
    применяя subscribe/unsubscribe к набору) и читает subscriptions(). Потокобезопасно."""

    def __init__(self):
        self._q: "queue.Queue[Message]" = queue.Queue()
        self._subs: Set[str] = set()
        self._subs_lock = threading.Lock()
        self._snap: Optional[Snapshot] = None
        self._snap_lock = threading.Lock()

    # ---- клиентская сторона ----
    def subscribe(self, addr: str) -> None:
        self._q.put(Message(SUBSCRIBE, addr))

    def unsubscribe(self, addr: str) -> None:
        self._q.put(Message(UNSUBSCRIBE, addr))

    def command(self, addr: str, value: Any) -> None:
        self._q.put(Message(COMMAND, addr, value))

    def control(self, action: str) -> None:
        self._q.put(Message(CONTROL, "", action))

    def set_screen(self, addr: str, area: Any, res: Any, point: Any = None, quantity: Any = None) -> None:
        self._q.put(Message(SET_SCREEN, addr, {"area": area, "res": res, "point": point, "quantity": quantity}))

    def poll_snapshot(self) -> Optional[Snapshot]:
        with self._snap_lock:
            return self._snap

    # ---- движковая сторона ----
    def pull(self) -> List[Message]:
        """Осушить очередь сообщений. subscribe/unsubscribe применяются к набору подписок транспорта
        прямо здесь (чтобы после pull() движок видел актуальный subscriptions()); остальные — движку."""
        msgs: List[Message] = []
        while True:
            try:
                m = self._q.get_nowait()
            except queue.Empty:
                break
            if m.kind == SUBSCRIBE:
                with self._subs_lock:
                    self._subs.add(m.addr)
            elif m.kind == UNSUBSCRIBE:
                with self._subs_lock:
                    self._subs.discard(m.addr)
            msgs.append(m)
        return msgs

    def subscriptions(self) -> Set[str]:
        with self._subs_lock:
            return set(self._subs)

    def publish(self, snapshot: Snapshot) -> None:
        with self._snap_lock:
            self._snap = snapshot
