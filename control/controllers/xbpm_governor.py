# -*- coding: utf-8 -*-
"""Внешний МЕДЛЕННЫЙ контур каскада B: XBPM-губернатор.

Читает положение фотонного пучка (XBPM), медленно подправляет корректоры, чтобы свести фотон на цель.
Чтобы быстрый FOFB не затёр это, губернатор пишет НЕ только свой вклад в корректоры, но и УСТАВКУ
(golden orbit) внутреннему контуру: ref = Hc·θ_slow — та e-орбита, которую создаёт медленный вклад.
FOFB затем держит e-орбиту НА этой уставке и медленную коррекцию не отменяет.

Закон: θ_slow ← θ_slow − g·R_xbpm⁺·(фотон). Возврат: вклад в корректоры + ref:<golden>.
"""
import numpy as np

from control.controller import Controller


class XBPMGovernor(Controller):
    reads = ["photon.x"]

    def __init__(self, name="xbpm", sensor="photon.x", actuator_key="hkick",
                 ref_sensor="ebpm.x", ref_name="golden.x", gain=0.5, rcond=1e-4):
        super().__init__(name, gain=gain, rcond=rcond)
        self.sensor = sensor
        self.actuator_key = actuator_key
        self.ref_sensor = ref_sensor
        self.ref_name = ref_name

    def setup(self, resources, schema):
        self.act = [p for p in schema.ids if self.actuator_key in p]
        self.writes = self.act
        self.P = resources.response(self.act, [self.sensor]).pinv(self.params["rcond"])   # фотон→корректор
        self.Hc = resources.response(self.act, [self.ref_sensor]).R                        # корректор→e-орбита (для golden)
        self.ts = np.zeros(len(self.act))

    def update(self, t, view):
        ph = view.read(self.sensor)
        if ph is None:
            return {}
        self.ts = self.ts - self.params["gain"] * (self.P @ ph)
        out = {pid: float(self.ts[i]) for i, pid in enumerate(self.act)}
        out["ref:" + self.ref_name] = self.Hc @ self.ts      # golden orbit = e-орбита от медленного вклада
        return out

    def reset(self):
        if hasattr(self, "ts"):
            self.ts[:] = 0.0
