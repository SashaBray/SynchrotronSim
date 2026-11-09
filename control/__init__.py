# -*- coding: utf-8 -*-
"""control — модульная платформа управления орбитой (скелет рефакторинга).

Слои (см. память [[control-framework-modular-design]]):
  state.py      — ExperimentState (blackboard) + Limits (ограничители slew/ОДЗ)
  controller.py — Controller (контракт плагина) + Resources (провайдер матриц отклика)
  scheduler.py  — Scheduler (многотемповый такт, сумма вкладов по слотам)
  plant.py      — Plant (абстракция установки) + DemoPlant (игрушечная, чтобы скелет бегал)
  registry.py   — обнаружение плагинов в control/controllers/*.py

Диффурное ядро и инструменты полей/токов НЕ трогаются — это отдельная ветка.
"""
from .state import (Schema, StateVector, Limits, History, ExperimentState, StateView,
                    SensorConfig, NoiseModel)
from .controller import Controller, Resources, ResponseMatrix
from .scheduler import Scheduler, LoopSpec
from .plant import Plant, DemoPlant
from .registry import discover

__all__ = ["Schema", "StateVector", "Limits", "History", "ExperimentState", "StateView",
           "SensorConfig", "NoiseModel",
           "Controller", "Resources", "ResponseMatrix", "Scheduler", "LoopSpec",
           "Plant", "DemoPlant", "discover"]
