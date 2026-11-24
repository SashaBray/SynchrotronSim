# -*- coding: utf-8 -*-
"""Контракт регулятора (то, что пишет пользователь) и провайдер ресурсов (матрицы отклика).

Регулятор — ЧИСТЫЙ алгоритм: читает показания/состояние, возвращает свой ВКЛАД в target.
Он не видит кольцо и не трогает plant; темп задаёт планировщик (не сам регулятор).
Матрица отклика — не хардкод, а РЕСУРС: регулятор просит её у Resources в setup().

МОДУЛЬНОСТЬ: `provides()` даёт ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ СВОЕЙ наблюдаемой — как она считается из `ObsContext`
(и какие тяжёлые промежуточные вводит). Так контур — ОДИН самодостаточный файл (закон + вычисление сенсора);
платформа собирает реестр Observer'а из активных контуров, а не хардкодит его в plant'е/observer.py.
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class ObsContribution:
    """Вклад контура в схему вычисления наблюдаемых (см. control/observer.py `Observer`/`ObsContext`).
      observables   — {obs_id: fn(ctx)->ndarray}  как посчитать наблюдаемую этого контура из контекста;
      intermediates — {key: fn(ctx)->value}       тяжёлые ОБЩИЕ шаги, которые контур вводит (кэш в ctx)."""
    observables: dict = field(default_factory=dict)
    intermediates: dict = field(default_factory=dict)


class ResponseMatrix:
    """Матрица отклика: строки — наблюдаемые, столбцы — актуаторы."""

    def __init__(self, R, rows, cols):
        self.R = np.asarray(R, float)
        self.rows = list(rows)
        self.cols = list(cols)

    def pinv(self, rcond=1e-2):
        return np.linalg.pinv(self.R, rcond=rcond)


class Resources:
    """Провайдер ресурсов. Строит/кэширует матрицы отклика по отпечатку конфигурации.
    Матрица строится ОДИН раз и делится между контурами; в такте — никогда.
    На реальной машине тот же интерфейс может отдать ИЗМЕРЕННУЮ матрицу (LOCO) — sim/real parity."""

    def __init__(self, plant):
        self.plant = plant
        self._cache = {}

    def response(self, actuators, observables):
        key = (tuple(actuators), tuple(observables), self.plant.fingerprint())
        if key not in self._cache:
            self._cache[key] = self.plant.build_response(list(actuators), list(observables))
        return self._cache[key]


class SlotWriter:
    """Императивная запись контура в СВОЙ слот target (+ уставки refs) — для контуров, которые пишут
    В КОДЕ, а не возвращают словарь. Арбитраж сохранён: планировщик обнулил слот ПЕРЕД update и суммирует
    слоты ПОСЛЕ (с Limits/slew). Писать можно ТОЛЬКО свой слот — чужие/actual напрямую не трогаем."""

    def __init__(self, state, slot):
        self._state = state
        self._slot = slot

    def set(self, pid, value):
        """Задать вклад контура по актуатору pid (в своём слоте)."""
        self._slot[pid] = float(value)

    def add(self, pid, value):
        """Добавить к текущему вкладу по pid."""
        self._slot[pid] = self._slot[pid] + float(value)

    def ref(self, name, value):
        """Задать уставку другому контуру (напр. golden orbit)."""
        self._state.refs[name] = np.asarray(value, float)


class Controller:
    """Базовый контракт плагина-регулятора. Пользователь наследует и реализует update().

      reads            — id сенсоров, которые читает (над observables);
      writes           — id параметров (актуаторов), в которые пишет вклад; можно задать в setup();
      default_rate_hz  — желаемый темп по умолчанию (планировщик может переопределить);
      wants_ctl        — True → update(t, view, ctl) пишет ИМПЕРАТИВНО в свой слот через ctl (SlotWriter);
                         False (деф.) → update(t, view) ВОЗВРАЩАЕТ словарь {id: значение}.
    """

    reads = []
    writes = []
    default_rate_hz = None
    wants_ctl = False

    def __init__(self, name=None, **params):
        self.name = name or type(self).__name__
        self.params = params

    def setup(self, resources, schema):
        """Разовая подготовка: запросить матрицы отклика, собрать псевдообратные, выделить состояние."""
        pass

    def update(self, t, view):
        """Главный метод. view — read-only состояние. Возврат: {id_параметра: значение} — ВКЛАД контура.
        Если wants_ctl=True — сигнатура update(t, view, ctl), пишем в ctl (возврат не нужен)."""
        return {}

    def reset(self):
        """Сброс внутреннего состояния (интеграторы и т.п.)."""
        pass

    def track(self, contribution):
        """БЕЗУДАРНОЕ включение: перед первым апдейтом после включения планировщик отдаёт ТЕКУЩИЙ вклад контура
        (его слот: {id: значение}; после «отпустить» — нули, после «заморозить» — удержанный вклад). Регулятор с
        памятью (интегратор) подхватывает его как своё состояние, иначе продолжит с накопленного до выключения
        и даст бросок. Дефолт — ничего (регулятор без памяти)."""
        pass

    def provides(self):
        """ВХОД В СХЕМУ ВЫЧИСЛЕНИЯ своей наблюдаемой (модульность): как считать свой сенсор из `ObsContext`
        и какие тяжёлые промежуточные контур вводит. Дефолт — ничего (контур читает уже готовую наблюдаемую,
        её даёт plant/другой контур). Платформа собирает реестр Observer'а из `provides()` активных контуров."""
        return ObsContribution()
