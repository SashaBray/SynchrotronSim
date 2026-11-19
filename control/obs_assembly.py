# -*- coding: utf-8 -*-
"""Сборка схемы вычисления наблюдаемых ИЗ АКТИВНЫХ КОНТУРОВ (модульность).

Каждый контур несёт свой вход в вычисление через `Controller.provides()` (ObsContribution). Здесь мы
собираем из набора контуров единый реестр наблюдаемых + билдеры промежуточных и строим `Observer`,
подписанный на те сенсоры, для которых КОНКРЕТНЫЙ plant даёт сырьё (`serves`). Так ни plant, ни observer.py
не хардкодят, что такое tune/chroma/adts — это приносят файлы-контуры."""
from .observer import Observer, OBSERVABLES, _INTERMEDIATES


def collect_obs(loops):
    """Слить вклады `provides()` всех контуров → (реестр наблюдаемых, билдеры промежуточных). Дедуп по ключу."""
    reg, interm = {}, {}
    for lp in loops:
        c = lp.provides()
        reg.update(c.observables)
        interm.update(c.intermediates)
    return reg, interm


def observer_for(loops, serves):
    """Observer для plant'а: БАЗИС (ebpm/photon + ringmap/one_turn) + вклады контуров; подписан на `serves`
    (сенсоры, для которых этот plant даёт сырьё). Наблюдаемая без сырья на этом plant'е просто пропускается."""
    reg = dict(OBSERVABLES)                     # базис: орбитные пас-сквозь
    interm = dict(_INTERMEDIATES)               # базис: общие примитивы кольца
    loop_reg, loop_interm = collect_obs(loops)
    reg.update(loop_reg)                        # наблюдаемые из контуров (tune/chroma/adts)
    interm.update(loop_interm)                  # их промежуточные (chroma_raw, normal_form)
    serves = set(serves)
    return Observer(registry=reg, intermediates=interm,
                    subscribe=[s for s in reg if s in serves])
