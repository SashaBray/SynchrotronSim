# -*- coding: utf-8 -*-
"""Контур ЧАСТОТ (квадруполи → νx,νy к проектным). Самодостаточный файл-плагин: закон регулирования
И ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ своей наблюдаемой `tune` (как она считается из кольца) — требование модульности.
Платформа собирает реестр Observer'а из `provides()` активных контуров (см. control/obs_assembly.py)."""
import numpy as np

from control.controller import ObsContribution
from control.controllers.integral_svd import IntegralSVD


# --- ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ (живёт в файле контура) ---
def _b_chroma_raw(ctx):
    """Промежуточное: nux/nuy + xi_x/xi_y за ОДИН проход карты оборота (делится с контуром chroma)."""
    return ctx.get("ringmap").chromaticity(ctx.get("one_turn"))

def _tune(ctx):
    ch = ctx.get("chroma_raw"); return np.array([float(ch["nux"]), float(ch["nuy"])])


class TuneLoop(IntegralSVD):
    """Квадруполи (`quad.`) → [νx,νy] к уставке `tune.target`. Закон — интегральный SVD (ядро общее с FOFB)."""

    def __init__(self, name="tune", gain=0.7, rcond=1e-2, ref="tune.target", **kw):
        super().__init__(name, sensor="tune", actuator_key="quad.", gain=gain, rcond=rcond, ref=ref, **kw)

    def provides(self):
        return ObsContribution(observables={self.sensor: _tune},
                               intermediates={"chroma_raw": _b_chroma_raw})
