# -*- coding: utf-8 -*-
"""Планировщик — многотемповый такт над общим состоянием.

Пользователь объявляет список регуляторов и их частоты (add). Планировщик на каждом БАЗОВОМ такте:
  1. даёт должным (по их темпу) контурам read-only view → каждый пишет свой слот target;
  2. эффективная цель = config + Σ слотов (сумма по слотам), клип Limits применяет plant-конвейер;
  3. plant делает шаг (впрыск dynamics, интегрирование), обновляет observables + history;
  4. лог выбранных сигналов.
Причинность естественная: контур на такте k читает измерение ПРОШЛОГО шага, его вклад влияет на
измерение ТЕКУЩЕГО шага (задержка в один такт, без алгебраической петли).

Частота/вкл-выкл — здесь, а не в алгоритме: это ручки эксперимента (меняем темп, не трогая код).
"""
from dataclasses import dataclass

import numpy as np

from .state import StateView
from .controller import Resources, Controller, SlotWriter
from .arbiter import Arbiter


@dataclass
class LoopSpec:
    controller: Controller
    period: int                     # период в базовых тактах
    enabled: bool = True
    mode: str = "release"           # при выключении: release (обнулить слот) | freeze (держать вклад)
    rate_hz: float = 0.0


class Scheduler:
    def __init__(self, plant, state, dt, resources=None, log=None, arbiter=None):
        self.plant = plant
        self.state = state
        self.dt = dt
        self.resources = resources if resources is not None else Resources(plant)
        self.arbiter = arbiter if arbiter is not None else Arbiter()   # сборка команды (дефолт = Σ слотов)
        self.loops = {}
        self.tick = 0
        self._log_ids = list(log or [])
        self.log = {k: [] for k in self._log_ids}
        self.plant.reset(self.state)                         # начальные observables + history

    def add(self, controller, rate_hz=None, enabled=True, mode="release"):
        rate = rate_hz or controller.default_rate_hz
        if not rate:
            raise ValueError(f"не задан темп контура {controller.name}")
        controller.setup(self.resources, self.state.schema)
        controller.reset()
        period = max(1, round((1.0 / rate) / self.dt))
        self.loops[controller.name] = LoopSpec(controller, period, enabled, mode, rate)
        self.state.slot(controller.name)                     # создать слот
        return self.loops[controller.name]

    def set_enabled(self, name, on, mode=None):
        spec = self.loops[name]
        spec.enabled = on
        if mode:
            spec.mode = mode

    def step(self):
        view = StateView(self.state, self.arbiter)
        for name, spec in self.loops.items():
            slot = self.state.slot(name)
            if spec.enabled and (self.tick % spec.period == 0):
                slot.v[:] = 0.0                                             # обнулить СВОЙ слот ПЕРЕД апдейтом
                if getattr(spec.controller, "wants_ctl", False):
                    spec.controller.update(self.state.t, view, SlotWriter(self.state, slot))  # пишет в слот/refs сам
                else:
                    for pid, val in (spec.controller.update(self.state.t, view) or {}).items():
                        if isinstance(pid, str) and pid.startswith("ref:"):
                            self.state.refs[pid[4:]] = np.asarray(val, float)   # уставка другому контуру
                        else:
                            slot[pid] = val                                     # вклад в target (сумма по слотам)
            elif not spec.enabled and spec.mode == "release":
                slot.v[:] = 0.0
            # freeze: слот остаётся как был
        req = self.arbiter.combine(self.state)                         # слоты → единая команда (дефолт = Σ)
        self.state.actual.v = self.state.limits.apply(self.state.actual.v, req, self.dt)
        self.plant.step(self.state, self.dt)
        self.arbiter.publish(self.state)                               # выгрузка наружу (дефолт no-op)
        for k in self._log_ids:
            val = self.state.observables.get(k)
            self.log[k].append(float(np.sqrt(np.mean(np.asarray(val) ** 2))) if val is not None else np.nan)
        self.state.t += self.dt
        self.tick += 1

    def run(self, n_ticks):
        for _ in range(n_ticks):
            self.step()
        return self.log

    def table(self):
        """Строки (имя, частота, период_тактов, вкл) — для отчёта/интерфейса."""
        return [(n, s.rate_hz, s.period, s.enabled) for n, s in self.loops.items()]
