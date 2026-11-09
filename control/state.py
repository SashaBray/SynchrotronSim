# -*- coding: utf-8 -*-
"""Состояние эксперимента (blackboard) и ограничители.

Одна структура передаётся планировщиком каждому регулятору. Регионы:
  config    — статика машины (позиции/углы магнитов, поля, геометрия); база, меняется редко;
  dynamics  — мгновенное на такте t (шум поля, вибрация, дрейф) — впрыск возмущений в plant;
  target    — ВЕКТОР ТРЕБУЕМОГО СОСТОЯНИЯ, та же схема, что config/actual: регулятор может
              командовать ЧЕМ УГОДНО; хранится как СУММА ВКЛАДОВ ПО СЛОТАМ (по одному на контур);
  actual    — текущее реализованное состояние машинных параметров (пишет plant через Limits);
  observables — показания сенсоров (выборки над history + шум); это читают регуляторы;
  history   — усреднённый за итерацию снимок пучка во ВСЕХ плоскостях (вход/выход магнитов).

Дисциплина записи: регулятор читает всё, но пишет ТОЛЬКО свой слот target (через возврат из update).
actual/history/observables/config для регулятора — только чтение.
"""
import numpy as np


class Schema:
    """Упорядоченный список id параметров машины (один для config/actual/target)."""

    def __init__(self, ids):
        self.ids = list(ids)
        self._idx = {p: i for i, p in enumerate(self.ids)}

    def __len__(self):
        return len(self.ids)

    def index(self, pid):
        return self._idx[pid]

    def indices(self, ids):
        return np.array([self._idx[p] for p in ids], dtype=int)

    def __contains__(self, pid):
        return pid in self._idx


class StateVector:
    """Значения над Schema: доступ по id и как numpy-массив (.v)."""

    def __init__(self, schema, values=None):
        self.schema = schema
        self.v = np.zeros(len(schema)) if values is None else np.asarray(values, float).copy()

    def __getitem__(self, pid):
        return float(self.v[self.schema.index(pid)])

    def __setitem__(self, pid, val):
        self.v[self.schema.index(pid)] = val

    def select(self, ids):
        return self.v[self.schema.indices(ids)]

    def copy(self):
        return StateVector(self.schema, self.v.copy())


class Limits:
    """Ограничители по каждому параметру — модель реализуемости актуатора.
      odz  — область допустимых значений (жёсткий клип);
      slew — макс |Δ| в СЕКУНДУ (задаёт темп актуатора; inf = мгновенно);
      act  — управляем ли параметр вообще (False = константа, не двигается);
      res  — шаг квантования БП/ЦАП (0 = нет).
    Применяется plant ОДИН раз за такт, ПОСЛЕ суммирования вкладов контуров."""

    def __init__(self, schema):
        n = len(schema)
        self.schema = schema
        self.odz_lo = np.full(n, -np.inf)
        self.odz_hi = np.full(n, np.inf)
        self.slew = np.full(n, np.inf)
        self.act = np.ones(n, bool)
        self.res = np.zeros(n)

    def set(self, pid, odz=None, slew=None, act=None, res=None):
        i = self.schema.index(pid)
        if odz is not None:
            self.odz_lo[i], self.odz_hi[i] = odz
        if slew is not None:
            self.slew[i] = slew
        if act is not None:
            self.act[i] = act
        if res is not None:
            self.res[i] = res

    def apply(self, actual_v, requested_v, dt):
        """actual → requested с ограничением скорости (slew·dt) и клипом в ОДЗ. Возврат нового actual_v."""
        step = requested_v - actual_v
        m = self.slew * dt                                   # inf·dt = inf → без ограничения
        step = np.clip(step, -m, m)
        x = actual_v + step
        x = np.clip(x, self.odz_lo, self.odz_hi)
        rr = self.res > 0
        if np.any(rr):
            x[rr] = np.round(x[rr] / self.res[rr]) * self.res[rr]
        return np.where(self.act, x, actual_v)               # неуправляемые не двигаем


class History:
    """Усреднённый за последнюю итерацию снимок пучка в плоскостях.
    planes — id плоскостей (вход/выход магнитов); moments — имена (x, xp, y, yp, …).
    Сенсоры — это ВЫБОРКИ отсюда: e-BPM = строки мониторов, XBPM = строка источника + плечо."""

    def __init__(self, planes, moments=("x", "xp", "y", "yp")):
        self.planes = list(planes)
        self.moments = list(moments)
        self._pidx = {p: i for i, p in enumerate(self.planes)}
        self._midx = {m: i for i, m in enumerate(self.moments)}
        self.data = np.zeros((len(self.planes), len(self.moments)))

    def at(self, plane):
        return self.data[self._pidx[plane]]

    def column(self, moment):
        return self.data[:, self._midx[moment]]


class ExperimentState:
    """Blackboard эксперимента. См. модуль-docstring."""

    def __init__(self, schema, limits=None, history=None):
        self.schema = schema
        self.config = StateVector(schema)            # база (из .lte, дословно)
        self.actual = StateVector(schema)            # реализованное машинное состояние
        self.dynamics = StateVector(schema)          # мгновенное возмущение параметров на такте t
        self.limits = limits if limits is not None else Limits(schema)
        self.target_slots = {}                       # имя контура -> StateVector вкладов (дельты от config)
        self.observables = {}                        # id сенсора -> ndarray (показания)
        self.refs = {}                               # уставки между контурами (напр. golden orbit): имя -> ndarray
        self.history = history                       # History или None
        self.t = 0.0

    def slot(self, name):
        if name not in self.target_slots:
            self.target_slots[name] = StateVector(self.schema)
        return self.target_slots[name]

    def effective_target(self):
        """config_base + Σ вкладов всех контуров (сумма по слотам)."""
        v = self.config.v.copy()
        for s in self.target_slots.values():
            v = v + s.v
        return v


class StateView:
    """Read-only представление состояния для регулятора (пишет он только через возврат update)."""

    def __init__(self, state):
        self._s = state

    @property
    def t(self):
        return self._s.t

    @property
    def config(self):
        return self._s.config

    @property
    def actual(self):
        return self._s.actual

    @property
    def dynamics(self):
        return self._s.dynamics

    @property
    def history(self):
        return self._s.history

    def read(self, observable_id):
        """Показание сенсора (ndarray) или None, если ещё не измерено."""
        return self._s.observables.get(observable_id)

    def ref(self, name):
        """Уставка, заданная другим контуром (напр. golden orbit от внешнего контура), или None."""
        return self._s.refs.get(name)

    def target(self):
        return StateVector(self._s.schema, self._s.effective_target())
