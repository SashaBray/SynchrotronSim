# -*- coding: utf-8 -*-
"""Квазистатический plant ОПТИКИ: контуры квадруполей/секступолей/октуполей за интерфейсом Plant.

Отдельно от turn-by-turn орбитного FOFB (`tbt_plant`) — оптика медленная. Наблюдаемые считаются
АНАЛИТИЧЕСКИ (нормальная форма / дифф. алгебра `ring_map`) на ВОЗМУЩЁННОМ аналитическом кольце:
  • `tune`   = [νx, νy] (дробные)         — контур КВАДРУПОЛЕЙ (quad.<fam>.dk → Δk1);
  • `chroma` = [ξx, ξy]                   — контур СЕКСТУПОЛЕЙ (sext.<fam>.dk2 → Δk2);
  • `adts`   = [∂νx/∂Jx, ∂νy/∂Jy, ∂νx/∂Jy] — контур ОКТУПОЛЕЙ (oct.<fam>.dk3 → Δk3).

ПОЧЕМУ аналитически, а не через apply_perturb: apply_perturb меняет C++-ТРЕКЕР (частоты через трекинг
это видят), но аналитическое кольцо `_python_ring` НЕ отражает DK/DK2/DK3 (только YAW/θ₀) — значит
RingMap/normal_form, читающие k1/k2/k3 магнитов, возмущения не увидят. Поэтому Δk вносим ПРЯМО в силы
магнитов аналитической копии и считаем RingMap на ней: один проход даёт частоты+хроматизм+ADTS, все
реагируют на градиентные актуаторы, без шума трекинга. Один RingMap-проход на оценку.

Возмущение — фракционный дрейф градиентов (Δk = frac·k_номин по всем квад/секступолям/октуполям):
k1→частоты, k2→хроматизм, k3→ADTS. Матрицы отклика — аналитическая конечная разность (build_response).
Контуры — плагины IntegralSVD со своей уставкой (ref). model= — общий RingModel, НЕ закрывается здесь.
"""
import copy

import numpy as np

import ring_map as rmap

from .state import Schema, History
from .controller import ResponseMatrix
from .ring_model import RingModel
from .observer import Observer, ObsContext

FD_EPS = 1e-3                                                  # шаг конечной разности матриц отклика (абс. Δk)


class OpticsPlant:
    def __init__(self, model=None, lattice=26, nsl=8, observe=("tune", "chroma", "adts"),
                 drift_k1=0.0, drift_k2=0.0, drift_k3=0.0, rho=1.0, seed=1, order=4, observer=None):
        self.own = model is None
        self.rm = model if model is not None else RingModel(lattice=lattice, nsl=nsl)
        self.observe = tuple(observe); self.order = order
        self.base = self.rm.ring._python_ring()                # аналитическое кольцо (номинал)
        self.ids = list(self.rm.dk_ids) + list(self.rm.dk2_ids) + list(self.rm.dk3_ids)
        self._amap = {}                                        # актуатор → (поле k, [индексы семейства])
        for f, idx in self.rm.quad_fams.items():
            self._amap[f"quad.{f}.dk"] = ("k1", idx)
        for f, idx in self.rm.sext_fams.items():
            self._amap[f"sext.{f}.dk2"] = ("k2", idx)
        for f, idx in self.rm.oct_fams.items():
            self._amap[f"oct.{f}.dk3"] = ("k3", idx)
        mags = self.base.mags
        self.q_idx = [i for f in self.rm.quad_fams.values() for i in f]
        self.s_idx = [i for f in self.rm.sext_fams.values() for i in f]
        self.o_idx = [i for f in self.rm.oct_fams.values() for i in f]
        self.k1_nom = {i: float(getattr(mags[i], "k1", 0.0) or 0.0) for i in self.q_idx}
        self.k2_nom = {i: float(getattr(mags[i], "k2", 0.0) or 0.0) for i in self.s_idx}
        self.k3_nom = {i: float(getattr(mags[i], "k3", 0.0) or 0.0) for i in self.o_idx}
        self.drift = np.array([drift_k1, drift_k2, drift_k3]); self.rho = rho
        self._rng = np.random.default_rng(seed); self._a = np.zeros(3)
        # Observer собирает наблюдаемые из возмущённого кольца; дефолт — подписка = observe
        self.observer = observer or Observer(subscribe=list(self.observe))

    # ---- аналитическая оценка оптики на возмущённом кольце ----
    def _perturbed(self, pk1, pk2, pk3):
        pr = copy.copy(self.base); mags = list(self.base.mags)
        touched = {}                                          # копируем магнит один раз, даже если его трогают несколько каналов

        def _edit(i):
            if i not in touched:
                touched[i] = copy.copy(mags[i])
            return touched[i]

        # k1 = изменение ВОЗБУЖДЕНИЯ квадруполя → ОБЕ плоскости: k1 (νx) И k1v (νy) на ту же Δ.
        # Без k1v νy оставался инвариантным (у ядра раздельные k1/k1v; реальный квадруполь двигает оба).
        for i, dv in pk1.items():
            m = _edit(i)
            m.k1 = float(getattr(m, "k1", 0.0) or 0.0) + dv
            m.k1v = float(getattr(m, "k1v", 0.0) or 0.0) + dv
        for fld, pk in (("k2", pk2), ("k3", pk3)):             # секступоль/октуполь — одна компонента, без v-пары
            for i, dv in pk.items():
                m = _edit(i); setattr(m, fld, float(getattr(m, fld, 0.0) or 0.0) + dv)
        for i, m in touched.items():
            mags[i] = m
        pr.mags = mags
        return pr

    def _ctx(self, pk1, pk2, pk3, state=None):
        """Контекст наблюдателя на ВОЗМУЩЁННОМ кольце (общие промежуточные: карта оборота считается раз)."""
        return ObsContext(state, self, ring=self._perturbed(pk1, pk2, pk3), order=self.order)

    def _optics(self, pk1, pk2, pk3, only=None):
        """Посчитать наблюдаемые через реестр Observer'а (для build_response/nominal)."""
        names = list(only) if only is not None else list(self.observe)
        return self.observer.compute(self._ctx(pk1, pk2, pk3), names)

    def _contrib_maps(self, state):
        pk1, pk2, pk3 = {}, {}, {}
        tgt = {"k1": pk1, "k2": pk2, "k3": pk3}
        for aid in self.ids:                                   # вклад контуров (Δk по семействам)
            v = float(state.actual[aid])
            if v == 0.0:
                continue
            fld, idx = self._amap[aid]; d = tgt[fld]
            for i in idx:
                d[i] = d.get(i, 0.0) + v
        a1, a2, a3 = self._a                                   # возмущение (Δk = frac·k_номин)
        for a, nom, d in ((a1, self.k1_nom, pk1), (a2, self.k2_nom, pk2), (a3, self.k3_nom, pk3)):
            if a:
                for i, kn in nom.items():
                    d[i] = d.get(i, 0.0) + a * kn
        return pk1, pk2, pk3

    def nominal(self):
        """Оптика чистого кольца (уставки контуров)."""
        return self._optics({}, {}, {})

    def _advance(self):
        if self.rho >= 1.0:
            self._a = self.drift.copy()
        else:
            self._a = self.rho * self._a + np.sqrt(1.0 - self.rho ** 2) * self.drift * self._rng.standard_normal(3)

    def _measure(self, state):
        pk1, pk2, pk3 = self._contrib_maps(state)
        self.observer.update(state, self._ctx(pk1, pk2, pk3, state))   # Observer собирает наблюдаемые

    def schema(self):
        return Schema(self.ids)

    def history(self):
        return History(["opt"], moments=("x",))               # оптике history не нужен; формальная заглушка

    def reset(self, state):
        self._a = np.zeros(3); self._advance(); self._measure(state)

    def step(self, state, dt):
        self._advance(); self._measure(state)

    def build_response(self, actuators, observables):
        if len(observables) != 1:
            raise ValueError("OpticsPlant.build_response: один наблюдаемый на матрицу")
        obs = observables[0]
        base = self._optics({}, {}, {}, only=(obs,))[obs]
        cols = []
        for a in actuators:
            fld, idx = self._amap[a]
            pk = {i: FD_EPS for i in idx}
            mp = {"k1": {}, "k2": {}, "k3": {}}; mp[fld] = pk
            o = self._optics(mp["k1"], mp["k2"], mp["k3"], only=(obs,))[obs]
            cols.append((o - base) / FD_EPS)
        R = np.column_stack(cols)
        return ResponseMatrix(R, [f"{obs}.{i}" for i in range(R.shape[0])], list(actuators))

    def fingerprint(self):
        return ("optics", self.rm.lattice, self.rm.nsl, tuple(self.observe), self.order)

    def close(self):
        if self.own:
            self.rm.close()
