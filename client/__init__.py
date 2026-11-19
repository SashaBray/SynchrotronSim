# -*- coding: utf-8 -*-
"""client — живой интерактивный кокпит бесконечного режима (режим 4) SynchrotronSim.

Хребет (см. docs/CLIENT_PLAN.md):
  link.py       — шов Link (ABC) между движком и UI: сериализуемые данные через границу
                  (она же клиент-серверная). InProcLink — внутрипроцессный транспорт на потоках.
  addresses.py  — КОНТРАКТ адресов сигналов: канонический, СТАБИЛЬНЫЙ между запусками адрес →
                  перечисление/резолв значения/метаданные по StateView.
  engine.py     — EngineThread: тело mode5 (TbtPlant+Scheduler) в потоке, очередь команд,
                  модельное время Σdt vs стенные часы → оценка скорости, публикация снимка.

Ядро (диффур, карты поля, токи) ЗАМОРОЖЕНО: клиент только читает StateView и шлёт ручки штатным
путём (команда → Limits/привод → plant). Через шов ходят только числа + строковые адреса.
"""
from .link import Link, InProcLink, Message, Snapshot
from . import addresses

__all__ = ["Link", "InProcLink", "Message", "Snapshot", "addresses"]
