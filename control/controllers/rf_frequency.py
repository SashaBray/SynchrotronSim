# -*- coding: utf-8 -*-
"""Контур ВЧ-ЧАСТОТЫ (энергия пучка). Самодостаточный файл-плагин: закон регулирования И вход в схему
вычисления своей наблюдаемой `rf.delta` (требование модульности, см. control/obs_assembly.py).

ФИЗИКА. Пучок обязан набирать за оборот ровно h периодов ВЧ: путь/β = C_RF = h·c/f_RF. Если кольцо
«поплыло» (приливы, температура: ΔC) или частота не та, равновесная энергия сдвигается
δ ≈ −(ΔC/C + Δf/f)/η, и орбита в мониторах получает ДИСПЕРСИОННУЮ компоненту x = η_x·δ. Дипольные
корректоры её почти не берут (дисперсия вне бетатронного подпространства) — снимает её ЧАСТОТА ВЧ.

ДАТЧИК `rf.delta` — оценка энергии по орбите, как на реальной машине: проекция показаний BPM на модельную
дисперсию в мониторах, δ̂ = (η·x)/(η·η). Модельная η — с доски (state.optics, линейная оптика plant'а),
иначе из ra.twiss кольца (один раз). Калибровку усиления даёт матрица отклика dδ̂/dΔf, которую plant
снимает через ЭТОТ ЖЕ датчик (6D-равновесие, ±Δf) — неточность модельной η в отклик входит сама.

ЗАКОН. Интегральный (как FOFB/SOFB, ядро IntegralSVD): Δf ← Δf − g·R⁺·(δ̂ − δ_ref). Пучок отвечает на
шаг частоты синхротронным колебанием (f_s ≈ 1 кГц, затухание τ_δ ≈ 17 мс), поэтому показание
усредняется по n_avg окнам ≈ синхротронному периоду (273 оборота ≈ 2 окна по 135) — бокс-фильтр гасит
f_s, а полоса контура держится много ниже f_s. Актуатор `rf.freq` — отклонение от номинала h·c/C, Гц.
"""
from collections import deque

import numpy as np

from control.controller import ObsContribution
from control.controllers.integral_svd import IntegralSVD


# --- ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ (живёт в файле контура) ---
def _b_eta_bpm(ctx):
    """Промежуточное: модельная дисперсия η_x в мониторах (с доски, иначе twiss кольца — кэш на plant)."""
    mon = np.asarray(ctx.plant.mon, int)
    opt = getattr(ctx.state, "optics", None) if ctx.state is not None else None
    if opt is not None and np.asarray(opt.etax).size > int(mon.max()):
        return np.asarray(opt.etax, float)[mon]
    cached = getattr(ctx.plant, "_rf_eta_bpm", None)
    if cached is None:
        import ring_analysis as ra
        cached = np.asarray(ra.twiss(ctx.ring)["etax"], float)[mon]
        ctx.plant._rf_eta_bpm = cached
    return cached


def delta_estimate(orbit_x, eta):
    """Оценка энергии δ̂ = (η·x)/(η·η) — проекция горизонтальной орбиты на дисперсию в мониторах."""
    eta = np.asarray(eta, float); x = np.asarray(orbit_x, float)
    return float(eta @ x / (eta @ eta))


def _rf_delta(ctx):
    return np.array([delta_estimate(ctx.get("orbit_x"), ctx.get("rf_eta_bpm"))])


class RfFrequencyLoop(IntegralSVD):
    """Частота ВЧ (`rf.freq`) → энергия δ̂ (`rf.delta`) к уставке (ref, по умолчанию 0 — проектная энергия,
    т.е. нулевая дисперсионная компонента орбиты). Закон — интегральный SVD с усреднением по n_avg окнам."""

    def __init__(self, name="rf", gain=0.3, rcond=1e-2, ref=None, sensor="rf.delta", actuator_key="rf.freq",
                 n_avg=2, **kw):
        super().__init__(name, sensor=sensor, actuator_key=actuator_key, gain=gain, rcond=rcond, ref=ref, **kw)
        self.n_avg = max(1, int(n_avg))
        self._buf = deque(maxlen=self.n_avg)

    def update(self, t, view):
        m = view.read(self.sensor)
        if m is None:
            return {}
        self._buf.append(np.asarray(m, float))
        m = np.mean(self._buf, axis=0)                     # бокс-фильтр ≈ синхротронный период
        if self.ref is not None:
            r = view.ref(self.ref)
            if r is not None:
                m = m - r
        self.u = self.u - self.params["gain"] * (self.P @ m)
        return {pid: float(self.u[i]) for i, pid in enumerate(self.act)}

    def reset(self):
        super().reset()
        self._buf.clear()

    def provides(self):
        return ObsContribution(observables={self.sensor: _rf_delta},
                               intermediates={"rf_eta_bpm": _b_eta_bpm})
