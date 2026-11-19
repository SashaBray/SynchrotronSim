# -*- coding: utf-8 -*-
"""Контур ADTS (октуполи → ∂νx/∂Jx, ∂νy/∂Jy, ∂νx/∂Jy к проектным). Самодостаточный файл-плагин: закон
регулирования И ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ `adts`. ТЯЖЁЛЫЙ проход нормальной формы (order-4) ЖИВЁТ ЗДЕСЬ —
именно он дорогой; платформа решает, КОГДА его запускать (в линеаризованном режиме — один раз на старте,
чтобы построить матрицы отклика). Октуполи СИЛА-26 ранг-1 → контур по умолчанию ВЫКЛ (авторитет ограничен)."""
import numpy as np

import ring_map as rmap

from control.controller import ObsContribution
from control.controllers.integral_svd import IntegralSVD


# --- ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ (живёт в файле контура): тяжёлая нормальная форма 4-го порядка ---
def _b_normal_form(ctx):
    return rmap.normal_form(ctx.ring, ctx.order, rmap=ctx.get("ringmap"))

def _adts(ctx):
    nf = ctx.get("normal_form")
    return np.array([float(nf["dnux/dJx"]), float(nf["dnuy/dJy"]), float(nf["dnux/dJy"])])


class AdtsLoop(IntegralSVD):
    """Октуполи (`oct.`) → [∂νx/∂Jx,∂νy/∂Jy,∂νx/∂Jy] к уставке `adts.target`."""

    def __init__(self, name="adts", gain=0.7, rcond=1e-3, ref="adts.target", **kw):
        super().__init__(name, sensor="adts", actuator_key="oct.", gain=gain, rcond=rcond, ref=ref, **kw)

    def provides(self):
        return ObsContribution(observables={self.sensor: _adts},
                               intermediates={"normal_form": _b_normal_form})
