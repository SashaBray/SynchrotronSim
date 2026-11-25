# -*- coding: utf-8 -*-
"""Живой сеанс кокпита на вычислителе (docs/DATA_PROTOCOL.md, этап 3c): движковая половина шва Link по сети.

Соединение `mode="live"` с вычислителем (synchrotronsim.compute): в рукопожатии — опыт, ссылка на ревизию
машины (или документ), настройки кокпита (скорость привода зазора, перенесённые зазоры, темп публикации).
Сервер собирает модель (`build_session`), отвечает описанием сеанса (`client.engine.describe_session`: геометрия
сцены, контуры, оптика, экраны…), запускает EngineThread за внутренним InProcLink и дальше пересылает:
клиент → движок — сообщения Link; движок → клиент — снимки (при каждом новом кадре движка, не чаще 60 Гц).
Закрытие соединения — остановка движка и закрытие модели. Прогоны, записанные в живом сеансе, уходят в
хранилище площадки по окончании прохода.
"""
from __future__ import annotations

import threading
import traceback

from . import wire


def serve_live(srv, sock, hello) -> None:
    from client.link import InProcLink
    from client.engine import EngineThread, init_undulators, apply_gaps, describe_session
    from .experiment import Experiment
    send_lock = threading.Lock()

    def send(msg):
        with send_lock:
            wire.send_frame(sock, msg)

    try:
        exp = Experiment.from_dict(hello["experiment"])
        m = srv.machine_for(exp, hello.get("machine_ref"), hello.get("machine_doc"))
        srv.log(f"живой сеанс: опыт «{exp.name}», машина {m.lattice_id}@{m.fingerprint[:12]} — сборка…")
        plant, state, sch = srv.live_session(exp, m)
        sch.experiment = exp
        init_undulators(state, float(hello.get("gap_slew_mm_s", 500.0)))
        apply_gaps(state, hello.get("gaps"))
        desc = describe_session(plant, state, sch, exp)
    except Exception as e:                                   # noqa: BLE001 — ошибку сборки — клиенту
        srv.log(f"живой сеанс: ошибка сборки {type(e).__name__}: {e}\n{traceback.format_exc()}")
        try:
            send({"op": "error", "message": f"сборка модели на вычислителе: {type(e).__name__}: {e}"})
        except OSError:
            pass
        return
    link = InProcLink()
    for w in exp.record or []:                                # подписка опыта (адреса) — отмеченные величины
        if not w.recorder:
            link.watch(w.addr)
    if exp.run.max_turns:
        link.record(False, exp.run.max_turns, exp.run.repeat)
    eng = EngineThread(link, sch, steps_per_cycle=int(exp.run.steps_per_cycle or 9),
                       publish_every=int(hello.get("publish_every", 3)), run_root=srv.root)
    eng.runner.store = srv.store                               # ревизия машины в паспорте — в хранилище площадки
    if hasattr(srv.store, "put_run"):
        from .store_sqlite import pack_run

        def upload(run):
            srv.store.put_run(run.run_id, pack_run(run), run.manifest)
        eng.runner.on_finish = upload
    eng.start()
    srv.live_sessions += 1
    stop = threading.Event()

    def reader():
        try:
            while not stop.is_set():
                msg = wire.recv_frame(sock)
                op = msg.get("op")
                if op == "close":
                    break
                elif op == "subscribe":
                    link.subscribe(msg["addr"])
                elif op == "unsubscribe":
                    link.unsubscribe(msg["addr"])
                elif op == "command":
                    link.command(msg["addr"], msg.get("value"))
                elif op == "control":
                    link.control(msg.get("value"))
                elif op == "set_screen":
                    link.set_screen(msg["addr"], msg.get("area"), msg.get("res"), msg.get("point"), msg.get("quantity"))
                elif op == "watch":
                    link.watch(msg["addr"])
                elif op == "unwatch":
                    link.unwatch(msg["addr"])
                elif op == "record":
                    link.record(bool(msg.get("on")), msg.get("max_turns"), bool(msg.get("repeat")))
        except (ConnectionError, OSError, ValueError):
            pass
        finally:
            stop.set()

    try:
        send({"op": "session", "desc": desc})
        threading.Thread(target=reader, name="live-read", daemon=True).start()
        last = None
        while not stop.is_set():
            snap = link.poll_snapshot()
            if snap is not None and snap.meta.get("seq") != last:
                last = snap.meta.get("seq")
                send({"op": "snapshot", "values": snap.values, "meta": snap.meta})
            stop.wait(1.0 / 60.0)
    except OSError:
        pass
    finally:
        stop.set()
        eng.stop()
        eng.join(timeout=20)
        try:
            plant.close()
        except Exception:
            pass
        srv.live_sessions -= 1
        srv.log(f"живой сеанс «{exp.name}» закрыт")
        try:
            sock.close()
        except OSError:
            pass
