# -*- coding: utf-8 -*-
"""Observer — настраиваемый слой СБОРКИ наблюдаемых между симуляцией и контурами.

После шага plant отдаёт СЫРЬЁ (кольцо / усреднённую орбиту / снимок) в `ObsContext`; `Observer` по РЕЕСТРУ
`имя→ObsSpec(needs, compute)` и ПОДПИСКЕ заполняет `state.observables`. Тяжёлые ПРОМЕЖУТОЧНЫЕ (карта оборота,
нормальная форма) считаются ОДИН раз (кэш в `ObsContext`) и делятся между наблюдаемыми. Пользователь
кастомизирует реестр/подписку — это «обёртка, где пишешь логику данных», не трогая plant и контуры.

Слои: plant (сырьё + сенсорная модель) → Observer (сборка наблюдаемых) → контур (управление).
"""
from dataclasses import dataclass

import numpy as np

import ring_map as rmap


class ObsContext:
    """Контекст одного такта: СЫРЬЁ plant'а (`set`) + лениво-кэшируемые ПРОМЕЖУТОЧНЫЕ (`get`).

    `get(key)` отдаёт: сырьё (если plant его положил) → кэш → построенное билдером `_INTERMEDIATES[key]`
    (считается РАЗ, кэшируется). Так tune и chroma делят одну карту оборота, а adts — ту же RingMap."""

    def __init__(self, state, plant, ring=None, order=4):
        self.state = state
        self.plant = plant
        self.ring = ring
        self.order = order
        self._raw = {}
        self._cache = {}

    def set(self, key, value):
        self._raw[key] = value
        return self

    def get(self, key):
        if key in self._raw:
            return self._raw[key]
        if key in self._cache:
            return self._cache[key]
        if key not in _INTERMEDIATES:
            raise KeyError(f"ObsContext: нет сырья/билдера для {key!r}")
        self._cache[key] = _INTERMEDIATES[key](self)
        return self._cache[key]


# --- билдеры ПРОМЕЖУТОЧНЫХ (один раз на такт, кэш в ObsContext) --------------------------------------
def _b_ringmap(ctx):
    return rmap.RingMap(ctx.ring, ctx.order)

def _b_one_turn(ctx):
    return ctx.get("ringmap").one_turn()

def _b_chroma_raw(ctx):
    return ctx.get("ringmap").chromaticity(ctx.get("one_turn"))   # даёт nux/nuy и xi_x/xi_y сразу

def _b_normal_form(ctx):
    return rmap.normal_form(ctx.ring, ctx.order, rmap=ctx.get("ringmap"))

_INTERMEDIATES = {
    "ringmap": _b_ringmap,
    "one_turn": _b_one_turn,
    "chroma_raw": _b_chroma_raw,
    "normal_form": _b_normal_form,
}


@dataclass
class ObsSpec:
    needs: tuple           # какие промежуточные нужны (документация; реальный шеринг — через кэш ObsContext)
    compute: object        # fn(ctx) -> ndarray


# --- экстракторы НАБЛЮДАЕМЫХ ------------------------------------------------------------------------
def _tune(ctx):
    ch = ctx.get("chroma_raw"); return np.array([float(ch["nux"]), float(ch["nuy"])])

def _chroma(ctx):
    ch = ctx.get("chroma_raw"); return np.array([float(ch["xi_x"][1]), float(ch["xi_y"][1])])

def _adts(ctx):
    nf = ctx.get("normal_form")
    return np.array([float(nf["dnux/dJx"]), float(nf["dnuy/dJy"]), float(nf["dnux/dJy"])])

OBSERVABLES = {
    # оптика (OpticsPlant): tune/chroma делят промежуточный chroma_raw (одна карта оборота)
    "tune":     ObsSpec(("chroma_raw",), _tune),
    "chroma":   ObsSpec(("chroma_raw",), _chroma),
    "adts":     ObsSpec(("normal_form",), _adts),
    # орбита (TbtPlant): пас-сквозь усреднённого сырья (сенсорная модель — в plant'е)
    "ebpm.x":   ObsSpec(("orbit_x",), lambda ctx: ctx.get("orbit_x")),
    "ebpm.y":   ObsSpec(("orbit_y",), lambda ctx: ctx.get("orbit_y")),
    "photon.x": ObsSpec(("photon_x",), lambda ctx: ctx.get("photon_x")),
    "photon.y": ObsSpec(("photon_y",), lambda ctx: ctx.get("photon_y")),
}


class Observer:
    """Сборщик наблюдаемых: по подписке заполняет `state.observables` из `ObsContext`. Кастомизируется
    своим `registry` (добавить ObsSpec) и/или `subscribe` (что считать). Промежуточные — через кэш ctx."""

    def __init__(self, registry=None, subscribe=None):
        self.registry = dict(registry if registry is not None else OBSERVABLES)
        self.subscribe = None if subscribe is None else list(subscribe)

    def names(self):
        return list(self.subscribe) if self.subscribe is not None else list(self.registry)

    def compute(self, ctx, names=None):
        """Посчитать названные наблюдаемые из ctx → dict (для build_response/nominal). Промежуточные общие."""
        out = {}
        for n in (names if names is not None else self.names()):
            spec = self.registry.get(n)
            if spec is not None:
                out[n] = spec.compute(ctx)
        return out

    def update(self, state, ctx):
        """Заполнить state.observables подписанными наблюдаемыми. Нет сырья для сигнала → пропускаем."""
        for n in self.names():
            spec = self.registry.get(n)
            if spec is None:
                continue
            try:
                state.observables[n] = spec.compute(ctx)
            except KeyError:                                   # нет нужного сырья на этом plant'е — сигнал недоступен
                continue
