# -*- coding: utf-8 -*-
"""Пример плагина-регулятора: интегральная коррекция по SVD матрицы отклика.

Закон: u ← u − gain·R⁺·m, вклад контура = u по его актуаторам. Это ядро нашего FOFB/SOFB
(ср. fofb_sim.py / regulator.py) в форме сменного плагина. Пользователь копирует этот файл и
пишет свой алгоритм — контракт тот же.
"""
import numpy as np

from control.controller import Controller


class IntegralSVD(Controller):
    reads = ["bpm"]

    def __init__(self, name="fofb", sensor="bpm", actuator_key="corr.", gain=0.5, rcond=1e-2, ref=None):
        super().__init__(name, gain=gain, rcond=rcond)
        self.sensor = sensor
        self.actuator_key = actuator_key        # подстрока id актуаторов (напр. "hkick"/"vkick")
        self.ref = ref                          # имя уставки (golden orbit) от внешнего контура, или None

    def setup(self, resources, schema):
        self.act = [p for p in schema.ids if self.actuator_key in p]
        self.writes = self.act
        rm = resources.response(self.act, [self.sensor])
        self.P = rm.pinv(self.params["rcond"])
        self.u = np.zeros(len(self.act))

    def update(self, t, view):
        m = view.read(self.sensor)
        if m is None:
            return {}
        if self.ref is not None:                        # следим за golden orbit, а не за нулём
            r = view.ref(self.ref)
            if r is not None:
                m = m - r
        self.u = self.u - self.params["gain"] * (self.P @ m)
        return {pid: float(self.u[i]) for i, pid in enumerate(self.act)}

    def reset(self):
        if hasattr(self, "u"):
            self.u[:] = 0.0
