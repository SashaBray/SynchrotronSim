# -*- coding: utf-8 -*-
"""Установка (plant): абстракция + игрушечная реализация для проверки скелета.

Plant отвечает за: применить реализованное состояние, впрыснуть возмущения dynamics, проинтегрировать
пучок на один такт, обновить observables (сенсоры) и history (снимок по плоскостям), и построить
матрицу отклика (актуаторы→наблюдаемые). Настоящий plant обернёт ring_cpp.CppRing/track_* — это
следующий шаг рефактора; здесь DemoPlant моделирует линейную орбиту, чтобы потоки данных бегали.
"""
import numpy as np

from .state import Schema, History
from .controller import ResponseMatrix


class Plant:
    """Абстракция установки — контракт для планировщика."""

    def schema(self):
        raise NotImplementedError

    def history(self):
        raise NotImplementedError

    def reset(self, state):
        raise NotImplementedError

    def step(self, state, dt):
        raise NotImplementedError

    def build_response(self, actuators, observables):
        raise NotImplementedError

    def fingerprint(self):
        raise NotImplementedError


class DemoPlant(Plant):
    """Игрушечная линейная установка: корректоры → орбита в BPM, фикс. матрица отклика R.
    Возмущение — дрейф AR(1) в пространстве BPM (+ опц. пользовательский впрыск через dynamics).
    Сенсор bpm = истинная орбита (history) + шум измерения. Зеркалит нашу реальную постановку FOFB
    в миниатюре (ср. fofb_sim.py), поэтому интегральный SVD-регулятор давит дрейф как на настоящей машине."""

    def __init__(self, n_corr=12, n_bpm=40, seed=0, sig_drift=1e-7, meas_noise=1e-8, rho=0.999):
        rng = np.random.default_rng(seed)
        self.corr_ids = [f"corr.{i:02d}" for i in range(n_corr)]
        self.bpm_ids = [f"bpm.{i:02d}" for i in range(n_bpm)]
        self.R = rng.normal(0.0, 1.0, (n_bpm, n_corr))       # отклик корректор→орбита (усл. ед.)
        self._rng = np.random.default_rng(seed + 1)
        # возмущение живёт в КОРРЕКТИРУЕМОМ подпространстве (как feed-down квадруполей) — иначе
        # 12 корректоров не достанут изотропный шум в 40 BPM. На реальной машине есть и
        # некорректируемая часть; здесь для наглядности петли её опускаем.
        self.drift = np.zeros(n_corr)
        self.sig_drift = sig_drift
        self.meas_noise = meas_noise
        self.rho = rho

    def schema(self):
        return Schema(self.corr_ids)                         # машинные параметры = корректоры

    def history(self):
        return History(self.bpm_ids, moments=("x",))

    def _measure(self, state):
        corr = state.actual.select(self.corr_ids)            # реализованные кики (после Limits)
        orbit = self.R @ (corr + self.drift)                 # истинная орбита: кики + возмущение
        state.history.data[:, 0] = orbit                     # снимок по плоскостям
        state.observables["bpm"] = orbit + self.meas_noise * self._rng.standard_normal(len(orbit))

    def reset(self, state):
        self.drift[:] = 0.0
        self._measure(state)

    def step(self, state, dt):
        a = np.sqrt(1.0 - self.rho ** 2)
        self.drift = self.rho * self.drift + a * self.sig_drift * self._rng.standard_normal(len(self.drift))
        # dynamics: пользовательский впрыск (моде 5) — возмущение в пространстве корректоров
        dkick = state.dynamics.select(self.corr_ids)
        if np.any(dkick):
            self.drift = self.drift + dkick
        self._measure(state)

    def build_response(self, actuators, observables):
        if tuple(observables) != ("bpm",):
            raise ValueError(f"DemoPlant знает только наблюдаемую 'bpm', запрошено {observables}")
        cols = [self.corr_ids.index(a) for a in actuators]
        return ResponseMatrix(self.R[:, cols], self.bpm_ids, actuators)

    def fingerprint(self):
        return ("demo", self.R.shape)
