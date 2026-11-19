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

    def __init__(self, state, plant, ring=None, order=4, intermediates=None):
        self.state = state
        self.plant = plant
        self.ring = ring
        self.order = order
        self._raw = {}
        self._cache = {}
        # карта билдеров ПРОМЕЖУТОЧНЫХ — инъекция (её собирают из contributions контуров); дефолт — базис
        self._intermediates = intermediates if intermediates is not None else _INTERMEDIATES

    def set(self, key, value):
        self._raw[key] = value
        return self

    def get(self, key):
        if key in self._raw:
            return self._raw[key]
        if key in self._cache:
            return self._cache[key]
        if key not in self._intermediates:
            raise KeyError(f"ObsContext: нет сырья/билдера для {key!r}")
        self._cache[key] = self._intermediates[key](self)
        return self._cache[key]


# --- билдеры ПРОМЕЖУТОЧНЫХ (один раз на такт, кэш в ObsContext) --------------------------------------
# ЗДЕСЬ — только ОБЩИЕ примитивы кольца (не привязаны к наблюдаемой). Оптические промежуточные
# (chroma_raw, normal_form) ПЕРЕЕХАЛИ в файлы-контуры (optics_tune.py / optics_adts.py): их вводит тот
# контур, которому они нужны, через Controller.provides() — модульность (observer.py НЕ знает про tune/adts).
def _b_ringmap(ctx):
    return rmap.RingMap(ctx.ring, ctx.order)

def _b_one_turn(ctx):
    return ctx.get("ringmap").one_turn()

_INTERMEDIATES = {
    "ringmap": _b_ringmap,
    "one_turn": _b_one_turn,
}


@dataclass
class ObsSpec:
    needs: tuple           # какие промежуточные нужны (документация; реальный шеринг — через кэш ObsContext)
    compute: object        # fn(ctx) -> ndarray


# --- наблюдаемые БАЗИСА ------------------------------------------------------------------------------
# Оптику (tune/chroma/adts) здесь НЕТ — её несут файлы-контуры. Орбита — пас-сквозь усреднённого сырья
# (сенсорная модель живёт в plant'е; это НЕ «вычисление контура», а сырьё plant'а — остаётся базисом).
OBSERVABLES = {
    "ebpm.x":   ObsSpec(("orbit_x",), lambda ctx: ctx.get("orbit_x")),
    "ebpm.y":   ObsSpec(("orbit_y",), lambda ctx: ctx.get("orbit_y")),
    "photon.x": ObsSpec(("photon_x",), lambda ctx: ctx.get("photon_x")),
    "photon.y": ObsSpec(("photon_y",), lambda ctx: ctx.get("photon_y")),
}


class Observer:
    """Сборщик наблюдаемых: по подписке заполняет `state.observables` из `ObsContext`. Кастомизируется
    своим `registry` (добавить ObsSpec) и/или `subscribe` (что считать). Промежуточные — через кэш ctx."""

    def __init__(self, registry=None, subscribe=None, intermediates=None):
        self.registry = dict(registry if registry is not None else OBSERVABLES)
        self.subscribe = None if subscribe is None else list(subscribe)
        # билдеры ПРОМЕЖУТОЧНЫХ, которые plant прокинет в ObsContext (собираются из contributions контуров)
        self.intermediates = dict(intermediates) if intermediates is not None else dict(_INTERMEDIATES)

    def names(self):
        return list(self.subscribe) if self.subscribe is not None else list(self.registry)

    @staticmethod
    def _compute(spec, ctx):
        """Реестр терпит и ObsSpec(.compute), и голый callable fn(ctx) — плагины дают простые функции."""
        fn = getattr(spec, "compute", None) or spec
        return fn(ctx)

    def compute(self, ctx, names=None):
        """Посчитать названные наблюдаемые из ctx → dict (для build_response/nominal). Промежуточные общие."""
        out = {}
        for n in (names if names is not None else self.names()):
            spec = self.registry.get(n)
            if spec is not None:
                out[n] = self._compute(spec, ctx)
        return out

    def update(self, state, ctx):
        """Заполнить state.observables подписанными наблюдаемыми. Нет сырья для сигнала → пропускаем."""
        for n in self.names():
            spec = self.registry.get(n)
            if spec is None:
                continue
            try:
                state.observables[n] = self._compute(spec, ctx)
            except KeyError:                                   # нет нужного сырья на этом plant'е — сигнал недоступен
                continue
