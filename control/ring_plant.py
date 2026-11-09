# -*- coding: utf-8 -*-
"""Настоящий plant: линейная модель замкнутой орбиты РЕАЛЬНОЙ решётки СИЛА-26.

Матрицы отклика снимаются с ТЕКУЩЕГО кольца единым построителем (RingModel) — свежие, а не из
протухающего npz-кэша ([[matrix-core-refactor-single-source]]):
  Hc / Hcy  — корректор→орбита в мониторах (гор. hkick→ebpm.x, верт. vkick→ebpm.y), [n_mon × n_corr];
  Gx / Gy   — сдвиг квадруполя→орбита через feed-down: dx→ebpm.x, dy→ebpm.y. Интеграл силы K1L и знак
              уже внутри (актуаторы quad.NNN.dx/dy RingModel вносят кик K1L·сдвиг), поэтому отдельный
              массив K1L здесь больше не нужен.
Орбита на µm-уровне линейна по кикам ([[dynamics-time-domain-foundation]]), поэтому это штатный plant
для ОПЫТОВ С ОБРАТНОЙ СВЯЗЬЮ. Нелинейный track_tbt — отдельный, более тяжёлый вариант того же Plant.

schema = кики корректоров hkick/vkick (n_corr×2). Наблюдаемые = орбита в мониторах ebpm.x/ebpm.y.
Возмущение — дрейф позиций квадруполей (AR(1)) через feed-down; опц. пользовательский впрыск
(state.dynamics, моде 5); опц. заданная извне последовательность dist_seq (для точной сверки/parity).

model= можно передать общий RingModel (строится один раз, кольцо дорого) — он НЕ закрывается здесь;
без него plant строит свой и закрывает кольцо сразу после снятия матриц (как CascadePlant).
"""
import numpy as np

from .state import Schema, History
from .controller import ResponseMatrix
from .ring_model import RingModel


class RingPlant:
    def __init__(self, model=None, lattice=26, nsl=8, meas_noise=0.0,
                 sig_drift=2e-8, rho=0.99, seed=1, dist_seq=None, rebuild=False):
        own = model is None
        rm = model if model is not None else RingModel(lattice=lattice, nsl=nsl)
        self.Hc = rm.response(rm.hk_ids, ["ebpm.x"], rebuild=rebuild)["ebpm.x"]    # hkick → ebpm.x
        self.Hcy = rm.response(rm.vk_ids, ["ebpm.y"], rebuild=rebuild)["ebpm.y"]   # vkick → ebpm.y
        self.Gx = rm.response(rm.qx_ids, ["ebpm.x"], rebuild=rebuild)["ebpm.x"]    # dx → ebpm.x (feed-down, K1L внутри)
        self.Gy = rm.response(rm.qy_ids, ["ebpm.y"], rebuild=rebuild)["ebpm.y"]    # dy → ebpm.y (знак −K1L внутри)
        self.hk_ids = list(rm.hk_ids); self.vk_ids = list(rm.vk_ids)
        if own:
            rm.close()                                        # матрицы — обычный numpy, живое кольцо больше не нужно
        self.n_mon, self.n_corr = self.Hc.shape
        self.n_quad = self.Gx.shape[1]
        self.mon_ids = [f"ebpm.{j:03d}" for j in range(self.n_mon)]
        self._col = {p: i for i, p in enumerate(self.hk_ids)}
        self._col.update({p: i for i, p in enumerate(self.vk_ids)})
        self.meas_noise = meas_noise; self.sig_drift = sig_drift; self.rho = rho
        self._rng = np.random.default_rng(seed)
        self.dist_seq = dist_seq                              # (dx, dy) [n_tick+1, n_quad] или None
        self._p = 0
        self.dx = np.zeros(self.n_quad); self.dy = np.zeros(self.n_quad)

    def schema(self):
        return Schema(self.hk_ids + self.vk_ids)

    def history(self):
        return History(self.mon_ids, moments=("x", "y"))

    def _advance_drift(self):
        if self.dist_seq is not None:
            self.dx = self.dist_seq[0][self._p]; self.dy = self.dist_seq[1][self._p]
        else:
            a = np.sqrt(1.0 - self.rho ** 2)
            self.dx = self.rho * self.dx + a * self.sig_drift * self._rng.standard_normal(self.n_quad)
            self.dy = self.rho * self.dy + a * self.sig_drift * self._rng.standard_normal(self.n_quad)

    def _measure(self, state):
        th = state.actual.select(self.hk_ids); tv = state.actual.select(self.vk_ids)
        dyn_h = state.dynamics.select(self.hk_ids); dyn_v = state.dynamics.select(self.vk_ids)  # впрыск (моде 5)
        ox = self.Hc @ (th + dyn_h) + self.Gx @ self.dx
        oy = self.Hcy @ (tv + dyn_v) + self.Gy @ self.dy
        if self.meas_noise:
            ox = ox + self.meas_noise * self._rng.standard_normal(self.n_mon)
            oy = oy + self.meas_noise * self._rng.standard_normal(self.n_mon)
        state.history.data[:, 0] = ox; state.history.data[:, 1] = oy
        state.observables["ebpm.x"] = ox; state.observables["ebpm.y"] = oy

    def reset(self, state):
        self._p = 0
        self.dx = np.zeros(self.n_quad); self.dy = np.zeros(self.n_quad)
        self._advance_drift()                                # возмущение на p=0
        self._measure(state)

    def step(self, state, dt):
        self._p += 1
        self._advance_drift()
        self._measure(state)

    def build_response(self, actuators, observables):
        cols = [self._col[a] for a in actuators]
        if tuple(observables) == ("ebpm.x",):
            return ResponseMatrix(self.Hc[:, cols], self.mon_ids, actuators)
        if tuple(observables) == ("ebpm.y",):
            return ResponseMatrix(self.Hcy[:, cols], self.mon_ids, actuators)
        raise ValueError(f"RingPlant: наблюдаемая {observables} не поддержана")

    def fingerprint(self):
        return ("ring26", self.Hc.shape)
