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


class SensorConfig:
    """Живые настройки сенсора/группы (напр. e-BPM): кладутся в `ExperimentState.sensors[<id>]`.
    plant читает их НА ИЗМЕРЕНИИ, поэтому правка в бесконечном режиме подхватывается. Поля — то, что
    описывает сам датчик (не его показание): шум, окно усреднения, вкл/выкл, калибровка, уставка."""

    def __init__(self, meas_noise=None, avg_chunks=None, enabled=True, gain=1.0, setpoint=0.0, **extra):
        self.meas_noise = meas_noise          # СКО шума измерения (None = дефолт plant'а)
        self.avg_chunks = avg_chunks          # окно усреднения, тактов (None = дефолт plant'а)
        self.enabled = enabled
        self.gain = gain
        self.setpoint = setpoint
        self.extra = dict(extra)              # место под произвольные настройки конкретного датчика


class NoiseModel:
    """ЗАКОН генерации шума как ДАННЫЕ: список компонент-спеков (формат `disturbance_spec`) + ревизия `rev`.
    Кладётся в `ExperimentState.noise`; plant перечитывает его и ПЕРЕСОБИРАЕТ возмущение при смене `rev`.
    В бесконечном режиме пользователь правит `spec` и зовёт `touch()` (или `set(new_spec)`) — изменение
    закона подхватывается на следующем такте без пересоздания plant'а."""

    def __init__(self, spec=None):
        self.spec = list(spec or [])          # компоненты (name/target/plane/form/correlation/quantity/amp/freq/phase)
        self.rev = 0                          # счётчик ревизий: plant следит за ним

    def set(self, spec):
        """Заменить весь закон и пометить изменение."""
        self.spec = list(spec); self.rev += 1
        return self

    def touch(self):
        """Пометить, что spec правили на месте (plant пересоберёт возмущение)."""
        self.rev += 1
        return self


_SR_MOD = None


def _undulator_sr():
    """Ленивая загрузка `undulator_sr` (лежит в undulators/) — чтобы ядро control/ не зависело от него жёстко."""
    global _SR_MOD
    if _SR_MOD is None:
        try:
            import os
            import sys
            sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "undulators"))
            import undulator_sr as _m
            _SR_MOD = _m
        except Exception:
            _SR_MOD = False
    return _SR_MOD or None


class UndulatorConfig:
    """ЖИВАЯ настройка ондулятора на доске: ЗАЗОР (мм) + производные B0/K (lazy, через `undulator_sr`).

    Кладётся в `ExperimentState.undulators[<имя>]` (имя = group_id ондулятора, напр. 'und_01'). plant и
    рекордеры читают её НА КАЖДОМ такте/окне, поэтому правка зазора на лету (оператором, расписанием,
    возмущением или будущим контуром зазора) подхватывается без пересоздания plant'а. Контур видит её
    через `StateView.undulators`: читать `gap_mm` и ЗАДАВАТЬ уставку `command(gap)` — её контур/оператор.

    ПРИВОД С ОГРАНИЧЕНИЕМ СКОРОСТИ (как `Limits.slew` у корректоров): реализованный `gap_mm` едет к уставке
    `gap_target` со скоростью ≤ `slew_mm_s` (клип в ОДЗ `gap_min..gap_max`), plant двигает его каждый такт
    (`advance`). Зазор НЕ телепортируется; привод челюстей медленный (секунды) против динамики пучка (мкс).
    K/B0 считаются ОТ реализованного `gap_mm` на запрос, поэтому отражают реальное положение привода."""

    def __init__(self, gap_mm, name=None, slew_mm_s=0.5, gap_min=5.0, gap_max=20.0, **extra):
        self.gap_mm = float(gap_mm)           # РЕАЛИЗОВАННЫЙ зазор — его видит модель (B0/K/картина/фокус)
        self.gap_target = float(gap_mm)       # УСТАВКА привода — её задаёт оператор/контур (command)
        self.slew_mm_s = float(slew_mm_s)     # макс |Δзазор|/с: скорость привода челюстей; inf = мгновенно
        self.gap_min = float(gap_min)         # ОДЗ (ход привода): минимальный зазор
        self.gap_max = float(gap_max)         # ОДЗ: максимальный зазор
        self.name = name
        self.extra = dict(extra)              # место под пер-прибор настройки (фаза, качество, вкл/выкл)

    def command(self, gap_mm):
        """Задать УСТАВКУ зазора (клип в ОДЗ). Привод поедет к ней со скоростью slew_mm_s — см. advance."""
        self.gap_target = min(max(float(gap_mm), self.gap_min), self.gap_max)
        return self

    def advance(self, dt):
        """Продвинуть РЕАЛИЗОВАННЫЙ зазор к уставке за время dt: шаг ≤ slew·dt, клип в ОДЗ. Зовёт plant каждый такт.
        Уставка==реализованный (дефолт) → no-op, поведение не меняется. inf-скорость → мгновенно."""
        tgt = min(max(self.gap_target, self.gap_min), self.gap_max)
        step = tgt - self.gap_mm
        lim = self.slew_mm_s * float(dt)
        if lim < abs(step):
            step = lim if step > 0 else -lim
        self.gap_mm = min(max(self.gap_mm + step, self.gap_min), self.gap_max)
        return self.gap_mm

    def B0(self):
        """Пиковое поле при текущем зазоре, Т (ролл-офф из undulator_sr); None, если модуль недоступен."""
        sr = _undulator_sr()
        return float(sr.B0_of_gap(self.gap_mm)) if sr is not None else None

    def K(self):
        """Параметр отклонения K при текущем зазоре; None, если модуль недоступен."""
        sr = _undulator_sr()
        return float(sr.K_of_B0(sr.B0_of_gap(self.gap_mm))) if sr is not None else None


class BeamSnapshot:
    """Пучок по ВСЕЙ линии: центроид (x,px,y,py,z,δ) на КАЖDОЙ границе элемента + бегущее среднее по окну.
    Хранит ОКНО снимков (ссылки на массивы plant'а, не копии); `snap` = последний, `avg` считается lazy.
    Запрос по строке или по имени магнита (вход = строка index, выход = index+1) через карту magnets.
    Принцип «текущее + агрегат, без истории»: окно ограничено (deque maxlen)."""

    _COMP = {"x": 0, "px": 1, "y": 2, "py": 3, "z": 4, "delta": 5}

    def __init__(self, window, magnets=None):
        self._win = window            # deque снимков (ndarray (n_mag+1, 6)); обновляется plant'ом на месте
        self.magnets = magnets        # карта имён магнитов (name→index) или None

    @property
    def snap(self):
        return self._win[-1]

    @property
    def avg(self):
        return np.mean(np.asarray(self._win), axis=0)         # бегущее среднее по окну (lazy)

    def at(self, row, avg=False):
        return (self.avg if avg else self.snap)[int(row)]

    def plane(self, row, comp, avg=False):
        return float(self.at(row, avg)[self._COMP[comp]])

    def at_magnet(self, magnet, where="entry", avg=False):
        """Центроид на входе/выходе магнита (по имени/индексу через карту magnets)."""
        i = self.magnets.index(magnet)
        return self.at(i if where == "entry" else i + 1, avg)


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
        self.sensors = {}                            # id сенсора/группы -> SensorConfig (ЖИВЫЕ настройки сенсоров)
        self.undulators = {}                         # имя ондулятора -> UndulatorConfig (ЖИВОЙ зазор, видят контуры)
        self.noise = None                            # NoiseModel — ЗАКОН генерации шума (живой), или None
        self.beam = None                             # BeamSnapshot — пучок по всей линии (снимок+среднее), ставит plant
        self.machine = None                          # lazy-доступ к магнитам (позиция/поле), ставит plant
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

    def __init__(self, state, arbiter=None):
        self._s = state
        self._arb = arbiter               # тот же сборщик команды, что у планировщика (None → сумма слотов)

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

    @property
    def beam(self):
        """BeamSnapshot: пучок по всей линии (снимок/среднее), запрос по плоскости/имени магнита. Или None."""
        return self._s.beam

    @property
    def machine(self):
        """Lazy-доступ к магнитам: .position(magnet)/.field(magnet) по имени/индексу. Или None."""
        return self._s.machine

    @property
    def sensors(self):
        """Живые настройки сенсоров (id → SensorConfig)."""
        return self._s.sensors

    @property
    def undulators(self):
        """Живые настройки ондуляторов (имя → UndulatorConfig): зазор + K/B0. Контур читает и может двигать зазор."""
        return self._s.undulators

    def read(self, observable_id):
        """Показание сенсора (ndarray) или None, если ещё не измерено."""
        return self._s.observables.get(observable_id)

    def ref(self, name):
        """Уставка, заданная другим контуром (напр. golden orbit от внешнего контура), или None."""
        return self._s.refs.get(name)

    def target(self):
        v = self._arb.combine(self._s) if self._arb is not None else self._s.effective_target()
        return StateVector(self._s.schema, v)
