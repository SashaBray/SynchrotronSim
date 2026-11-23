# -*- coding: utf-8 -*-
"""CompositePlant — орбита (TbtPlant, turn-by-turn, инлайн) + ЛИНЕАРИЗОВАННАЯ оптика (tune/chroma/adts).

Зачем: в живом кокпите нужны ВСЕ контуры — орбитный FOFB И оптические (квадруполи/секступоли/октуполи).
Оптика считается дорогой нормальной формой (order-4, секунды–минуты) и не может идти каждый такт движка.
Поэтому: тяжёлую математику (из файлов-контуров через их `provides()`) гоняем ОДИН раз на старте — строим
проектную оптику (уставки), сдвиг от статического дрейфа и матрицы отклика (кэш на диск). В такте оптику
обновляем ЛИНЕЙНО: obs = nominal + drift_offset + R·u (микросекунды). Кокпит не замирает. У тюна и
хроматизма ещё живой член от зазоров ондуляторов: Δobs = R_gap·(fac − 1), fac = (K/K_ref)²
(OpticsPlant.gap_response, кэш gapmodel_*.pkl).

Контур Plant (утиный, как у Tbt/OpticsPlant): schema/history/reset/step/build_response/fingerprint/close.
Схема — ОБЪЕДИНЕНИЕ актуаторов (корректоры + quad./sext./oct.); StateVector id-индексный, каждый суб-plant
читает своё. Орбитные наблюдаемые (ebpm.x/y) даёт TbtPlant как раньше; tune/chroma/adts — линейная модель.
"""
import hashlib
import os
import pickle

import numpy as np

from .state import Schema, ExperimentState
from .controller import ResponseMatrix
from .tbt_plant import TbtPlant
from .optics_plant import OpticsPlant, FD_EPS
from .response import ring_fingerprint, CACHE_DIR

OPTICS_OBS = ("tune", "chroma", "adts")


def _optics_model(optics, order, drift, rebuild=False):
    """Разово (кэш на диск) построить линейную модель оптики: {nominal, drift_offset, R} — каждая наблюдаемая
    против СВОЕГО семейства (tune←quad, chroma←sext, adts←oct). Вычисления идут ЧЕРЕЗ `provides()` контуров
    (optics.observer собран из них). Ключ кэша — контентный отпечаток кольца + order + дрейф → правка решётки
    инвалидирует кэш сама. Это ЕДИНСТВЕННОЕ дорогое место (десятки проходов нормальной формы), один раз."""
    fam = {"tune": list(optics.rm.dk_ids), "chroma": list(optics.rm.dk2_ids), "adts": list(optics.rm.dk3_ids)}
    tag = hashlib.md5(repr((ring_fingerprint(optics.rm.ring), "optmodel", order, tuple(float(d) for d in drift),
                            fam, FD_EPS)).encode()).hexdigest()[:12]
    path = os.path.join(CACHE_DIR, f"optmodel_{tag}.pkl")
    if not rebuild and os.path.exists(path):
        try:
            with open(path, "rb") as f:
                nominal, drift_off, R = pickle.load(f)
            return nominal, drift_off, R, fam
        except Exception:
            pass
    nominal = {o: np.asarray(v, float) for o, v in optics.nominal().items() if o in fam}   # уставки (чистое кольцо)
    tmp = ExperimentState(optics.schema())                       # оптика с ДРЕЙФОМ (u=0) — что корректируют контуры
    optics._a = np.zeros(3); optics._advance()                   # дрейф градиентов — как в OpticsPlant.reset, но БЕЗ
    drift = optics._optics(*optics._contrib_maps(tmp))           # зазоров ИД: их вклад — живой линейный член (_gap_tune)
    drift_off = {o: np.asarray(drift[o], float) - nominal[o] for o in fam}
    R = {o: np.asarray(optics.build_response(fam[o], [o]).R, float) for o in fam}            # отклик obs←своё семейство
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump((nominal, drift_off, R), f)
    return nominal, drift_off, R, fam


class CompositePlant:
    """Орбита (TbtPlant) + линеаризованная оптика. Один scheduler/движок гоняет ВСЕ контуры."""

    def __init__(self, lattice=26, nsl=8, settle_turns=250, meas_noise=1e-7, n_photon=1, gap_mm=5.0,
                 stride=15, disturbance_spec=None, relative_orbit=True, optics_nsl=2,
                 drift=(1e-3, 1e-2, 3e-2), order=4, tbt_observer=None, optics_observer=None, **tbt_kw):
        self.tbt = TbtPlant(lattice=lattice, nsl=nsl, settle_turns=settle_turns, meas_noise=meas_noise,
                            n_photon=n_photon, gap_mm=gap_mm, stride=stride, disturbance_spec=disturbance_spec,
                            relative_orbit=relative_orbit, observer=tbt_observer, **tbt_kw)
        # Оптика (tune/chroma/adts) — ГЛОБАЛЬНЫЕ величины кольца, тонкая нарезка не нужна: считаем на
        # КРУПНОМ кольце (optics_nsl). ПРОВЕРЕНО: nsl=2 даёт ТУ ЖЕ оптику, что nsl=4 (сходимость контуров
        # к машинной точности). NB стоимость order-4 нормальной формы от nsl НЕ зависит (мономиальная
        # алгебра 210 мономов) — грубое кольцо экономит лишь сборку кольца, а матрицы всё равно ~15 мин.
        # ЭТО ОДИН раз: результат кэшируется на диск по отпечатку кольца (resp_/optmodel_*.pkl). Орбита — своим nsl.
        self.optics = OpticsPlant(lattice=lattice, nsl=optics_nsl, observe=OPTICS_OBS, drift_k1=drift[0],
                                  drift_k2=drift[1], drift_k3=drift[2], rho=1.0, order=order, gap_mm=gap_mm,
                                  observer=optics_observer)                 # своё аналитич. кольцо (нет гонки за C++)
        self.stride, self.TREV = self.tbt.stride, self.tbt.TREV             # движок читает их для dt
        self.order = order
        self.drift = tuple(drift)
        self._opt_ids = list(self.optics.ids)                               # dk/dk2/dk3 id в объединённой схеме
        self._nominal = {}; self._drift_off = {}; self._R = {}; self._fam = {}
        self._gap_ids = []; self._gap_R = {}                               # тюн/хроматизм ← живые зазоры ИД (линейно)
        self._gap_key = None; self._gap_d = {}

    def __getattr__(self, name):
        """Делегировать НЕопределённые атрибуты суб-plant'у орбиты: `ring`/`correctors`/`corr`/`magnets`/
        `_und_anchor`/`lattice`/`gap_mm`/`turn`/`alive` и т.п. → `self.tbt`. Так `build_geometry` рисует сцену
        и `EngineThread._meta` (turn/alive) работают при optics=True (CompositePlant не хранит их сам).
        Контрактные методы (schema/step/reset/build_response/…) определены на классе — сюда НЕ попадают."""
        try:
            tbt = self.__dict__["tbt"]
        except KeyError:
            raise AttributeError(name)
        return getattr(tbt, name)

    # ---- контур Plant ---------------------------------------------------------------------------
    def schema(self):
        return Schema(self.tbt.schema().ids + self._opt_ids)                # ОБЪЕДИНЕНИЕ (ключи не коллидируют)

    def history(self):
        return self.tbt.history()                                          # фазовый портрет — орбитный

    def reset(self, state):
        self.tbt.reset(state)                                              # ebpm/beam/machine/optics/undulators
        self._nominal, self._drift_off, self._R, self._fam = _optics_model(self.optics, self.order, self.drift)
        gap_obs = [o for o in self.optics.GAP_OBS if o in self._fam]      # тюн/хроматизм ← зазоры ИД (кэш на диске;
        if gap_obs:                                                        # первый раз ~2 мин)
            self._gap_ids, self._gap_R = self.optics.gap_response(gap_obs)
        self._gap_key = None; self._gap_d = {}
        for o in self._nominal:
            state.refs[f"{o}.target"] = self._nominal[o].copy()            # уставки оптических контуров = проектные
        self._measure_inline(state)                                        # seed observables (u=0 → nominal+drift)

    def _measure_inline(self, state):
        """Линейная оптика в такте: obs = nominal + drift_offset + R·u (u — реализованные Δk семейства);
        у тюна и хроматизма ещё + сдвиг от живых зазоров ондуляторов (_gap_shift)."""
        gap = self._gap_shift(state) if self.__dict__.get("_gap_ids") else {}   # __dict__: не уходить в делегат tbt
        for o, fids in self._fam.items():
            u = state.actual.select(fids)
            obs = self._nominal[o] + self._drift_off[o] + self._R[o] @ u
            if o in gap:
                obs = obs + gap[o]
            state.observables[o] = obs

    def _gap_shift(self, state):
        """Сдвиги наблюдаемых от живых (реализованных) зазоров ондуляторов: {obs: R_gap·(fac − 1)},
        fac = (K/K_ref)². Опора — проектное кольцо, все ИД на опорном зазоре (5 мм). Пересчёт только при
        смене зазоров (мкс)."""
        und = getattr(state, "undulators", None) or {}
        gaps = tuple(float(und[g].gap_mm) if g in und else None for g in self._gap_ids)
        if gaps != self._gap_key:
            ref = self.optics.gap_ref
            fac = np.array([self.optics.gap_fac(g) if g is not None and abs(g - ref) > 1e-9 else 1.0 for g in gaps])
            self._gap_d = {o: R @ (fac - 1.0) for o, R in self._gap_R.items()}
            self._gap_key = gaps
        return self._gap_d

    def step(self, state, dt):
        self.tbt.step(state, dt)                                           # БЫСТРО инлайн (орбита/пучок)
        self._measure_inline(state)                                        # мкс

    def build_response(self, actuators, observables):
        obs = tuple(observables)
        rf_ids = getattr(self.tbt, "rf_ids", None)
        if rf_ids and list(actuators) == rf_ids:                         # частота ВЧ → её наблюдаемые (орбитный суб-plant)
            return self.tbt.build_response(actuators, observables)
        if obs in (("ebpm.x",), ("ebpm.y",)):
            return self.tbt.build_response(actuators, observables)          # орбитный отклик (в tbt, кэширован)
        o = obs[0]
        if o not in self._R:
            raise ValueError(f"CompositePlant: наблюдаемая {observables} не поддержана")
        R = self._R[o]
        rows = [f"{o}.{i}" for i in range(R.shape[0])]
        return ResponseMatrix(R, rows, list(actuators))                    # оптический отклик (кэш на диске)

    def fingerprint(self):
        return ("composite", self.tbt.fingerprint(), self.optics.fingerprint())

    def close(self):
        try:
            self.tbt.close()
        finally:
            self.optics.close()
