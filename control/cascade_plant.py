# -*- coding: utf-8 -*-
"""Plant для каскада B: орбита e⁻ + положение ФОТОННОГО пучка, обе от корректоров.

Матрицы сняты с ТЕКУЩЕГО кольца единым построителем (RingModel): Hc — корректор→орбита в мониторах,
Csrc — корректор→центроид фотона в источниках (x + d·x'). Два возмущения:
  • быстрое — в корректорном пространстве (вибрация/квадрупольный feed-down), ВИДНО и e-BPM, и фотону;
  • медленное — увод НАВЕДЕНИЯ ИД (дрейф зазора/фазы), добавляется ТОЛЬКО к фотону и e-BPM не виден.
Смысл каскада: быстрый контур (e-BPM) давит первое, медленный (XBPM) — второе, которого e-BPM не видит.
Только горизонталь (механика вертикали идентична).
"""
import numpy as np

from .state import Schema, History
from .controller import ResponseMatrix


class CascadePlant:
    def __init__(self, n_photon=8, photon_d=15.0, meas_noise=1e-7, photon_noise=2e-8,
                 sig_fast=2e-8, rho_fast=0.9, sig_slow=3e-6, rho_slow=0.999, seed=1, rebuild=False):
        from .ring_model import RingModel
        rm = RingModel(lattice=26, nsl=8, n_photon=n_photon, photon_d=photon_d)
        out = rm.response(rm.hk_ids, ["ebpm.x", "photon.x"], rebuild=rebuild)   # один проход, 120 зондирований
        self.Hc = out["ebpm.x"]                              # [n_mon × n_corr]
        self.Csrc = out["photon.x"]                          # [n_photon × n_corr]
        self.hk_ids = list(rm.hk_ids)
        rm.close()
        self.n_corr = self.Hc.shape[1]; self.n_mon = self.Hc.shape[0]; self.n_xbpm = self.Csrc.shape[0]
        self.mon_ids = [f"ebpm.{j:03d}" for j in range(self.n_mon)]
        self.xbpm_ids = [f"xbpm.{i:03d}" for i in range(self.n_xbpm)]
        self._rng = np.random.default_rng(seed)
        self.meas_noise = meas_noise; self.photon_noise = photon_noise
        self.sig_fast = sig_fast; self.rho_fast = rho_fast
        self.sig_slow = sig_slow; self.rho_slow = rho_slow
        self.dfast = np.zeros(self.n_corr)                   # быстрое возмущение (корректорное простр.)
        self.pslow = np.zeros(self.n_xbpm)                   # медленный увод наведения ИД (только фотон)

    def schema(self):
        return Schema(self.hk_ids)

    def history(self):
        return History(self.mon_ids, moments=("x",))

    def _advance(self):
        af = np.sqrt(1 - self.rho_fast ** 2); asl = np.sqrt(1 - self.rho_slow ** 2)
        self.dfast = self.rho_fast * self.dfast + af * self.sig_fast * self._rng.standard_normal(self.n_corr)
        self.pslow = self.rho_slow * self.pslow + asl * self.sig_slow * self._rng.standard_normal(self.n_xbpm)

    def _measure(self, state):
        tot = state.actual.select(self.hk_ids) + self.dfast
        ox = self.Hc @ tot
        ph = self.Csrc @ tot + self.pslow                    # увод ИД виден ТОЛЬКО фотону
        state.history.data[:, 0] = ox
        state.observables["ebpm.x"] = ox + self.meas_noise * self._rng.standard_normal(self.n_mon)
        state.observables["photon.x"] = ph + self.photon_noise * self._rng.standard_normal(self.n_xbpm)

    def reset(self, state):
        self.dfast[:] = 0.0; self.pslow[:] = 0.0
        self._advance(); self._measure(state)

    def step(self, state, dt):
        self._advance(); self._measure(state)

    def build_response(self, actuators, observables):
        cols = [self.hk_ids.index(a) for a in actuators]
        obs = observables[0]
        if obs == "ebpm.x":
            return ResponseMatrix(self.Hc[:, cols], self.mon_ids, actuators)
        if obs == "photon.x":
            return ResponseMatrix(self.Csrc[:, cols], self.xbpm_ids, actuators)
        raise ValueError(f"CascadePlant: наблюдаемая {observables} не поддержана")

    def fingerprint(self):
        return ("cascade26", self.Hc.shape, self.Csrc.shape)
