# -*- coding: utf-8 -*-
"""Контур ХРОМАТИЗМА (секступоли → ξx,ξy к проектным). Самодостаточный файл-плагин: закон регулирования
И ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ своей наблюдаемой `chroma`. Промежуточное `chroma_raw` то же, что у контура tune
(дедуп по ключу при слиянии), поэтому дублируем его одной строкой — файл остаётся самодостаточным."""
import numpy as np

from control.controller import ObsContribution
from control.controllers.integral_svd import IntegralSVD


# --- ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ (живёт в файле контура) ---
def _b_chroma_raw(ctx):
    """nux/nuy + xi_x/xi_y за один проход карты оборота (тот же ключ `chroma_raw`, что у tune → дедуп)."""
    return ctx.get("ringmap").chromaticity(ctx.get("one_turn"))

def _chroma(ctx):
    ch = ctx.get("chroma_raw"); return np.array([float(ch["xi_x"][1]), float(ch["xi_y"][1])])


class ChromaLoop(IntegralSVD):
    """Секступоли (`sext.`) → [ξx,ξy] к уставке `chroma.target`."""

    def __init__(self, name="chroma", gain=0.7, rcond=1e-2, ref="chroma.target", **kw):
        super().__init__(name, sensor="chroma", actuator_key="sext.", gain=gain, rcond=rcond, ref=ref, **kw)

    def provides(self):
        return ObsContribution(observables={self.sensor: _chroma},
                               intermediates={"chroma_raw": _b_chroma_raw})
