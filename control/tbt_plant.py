# -*- coding: utf-8 -*-
"""НЕЛИНЕЙНЫЙ plant: реальный пучок в радиационном кольце СИЛА-26 за тем же интерфейсом Plant.

В отличие от RingPlant (линейная матрица орбиты — для БЫСТРЫХ опытов с ОС), здесь динамика настоящая:
кольцо с СИ+ВЧ (фаза 165°, пучок живёт тысячи оборотов), пучок из NP частиц ведём track_tbt чанками
по STRIDE оборотов (самосогласованный z), а показания сенсоров снимаем СНИМКОМ track_elements —
центроид пучка вдоль линии (готча: выход track_elements НЕ переносим в следующий чанк, его z = путь).
Проверенный рецепт из fofb_sim.py ([[dynamics-time-domain-foundation]]), обёрнутый в контракт Plant,
чтобы ТОТ ЖЕ планировщик и ТЕ ЖЕ плагины гоняли против реальной нелинейной машины (режимы 2/3).

Ключевое разделение:
  • ДИНАМИКА — нелинейная, на 6D-кольце с СИ+ВЧ (track_tbt);
  • МОДЕЛЬ регулятора (build_response) — ЛИНЕЙНАЯ матрица отклика с RingModel (4D), как на реальной
    машине: SVD-петля обращает линейную модель и применяет её к нелинейной реальности. Зазор
    модель↔реальность — ровно предмет опытов режима 2.

Актуаторы (schema): кики корректоров hkick/vkick; opt-in (rf_actuator=True) — ЧАСТОТА ВЧ `rf.freq` (Δf от
номинала h·c/C, Гц: актуатор энергии пучка, контур control/controllers/rf_frequency.py). Возмущение длины
кольца (disturbance quantity=circumference, ΔC) plant вносит как эквивалентную длину ВЧ-оборота C_RF − ΔC.
Сенсоры (observables): ebpm.x/ebpm.y (орбита в
мониторах), photon.x/photon.y (центроид фотонного пучка в прямых промежутках, x+d·x'); history хранит
снимок (x, px, y, py) во всех мониторах → фазовый портрет и датчики плоскости для режима 2. Возмущение —
feed-down квадруполей по фикс. пространственному паттерну (амплитуда статична/медленна); опц. впрыск
пользователя через state.dynamics (режим 5).
"""
import os
from collections import deque

os.environ.setdefault("SS_SECTOR", "1")

import numpy as np

import ring_cpp
import facility_report as fr
import ring_analysis as ra
import disturbance as dist

from .state import Schema, History, BeamSnapshot, UndulatorConfig, OpticsView
from .controller import ResponseMatrix
from .ring_model import RingModel, ID_POLE_DEVS, id_devices
from .observer import Observer, ObsContext

STRIDE = 15                                                   # оборотов в одном базовом такте планировщика
NP = 12                                                       # частиц в пучке
PHASE_DEG = 165.0                                             # фаза ВЧ для 6D-трекинга ([[dynamics-time-domain-foundation]])
# Полюса ондуляторов (короткие секторные диполи-«жёсткий край», θ ≈ 3e-4): им хватает 1 среза — на СИЛА-26 это
# 2280 из 4882 элементов и половина цены трекинга. С 1 срезом против 8: тюн −2e-8, орбита 5e-12, частица через
# 200 оборотов (СИ+ВЧ) 2e-5 от амплитуды; трекинг вдвое быстрее. Включается id_pole_nsl (клиент — по умолчанию 1).


class TbtPlant:
    def __init__(self, lattice=26, nsl=8, n_part=NP, stride=STRIDE, dist_amp=0.3e-6,
                 dist_rho=1.0, meas_noise=0.0, settle_turns=300, seed=1, dist_seed=11, n_photon=8,
                 avg_chunks=9, resp_dims=4, sr=True, disturbance_spec=None, gap_mm=5.0, observer=None,
                 relative_orbit=False, rf_actuator=False, center_bunch=None, id_pole_nsl=None):
        nsl_dev = {d: int(id_pole_nsl) for d in ID_POLE_DEVS} if id_pole_nsl else None   # срезы полюсов ИД
        self.stride = stride
        self.gap_mm = float(gap_mm)                           # стартовый зазор ондуляторов (мм); живёт в state.undulators
        self.sr = bool(sr)                                    # СИ вкл/выкл (режим 2: затухание vs сохранение)
        self.meas_noise = meas_noise
        self.settle_turns = settle_turns
        self.avg_chunks = avg_chunks                          # окно усреднения орбиты ≈ такт (бетатрон усредняется)
        self.lattice = lattice; self.nsl = nsl
        self.resp_dims = resp_dims                            # 4 — модель контура 4D (δ=0); 6 — на своём 6D-кольце (LOCO-подобно)
        self._rng = np.random.default_rng(seed)

        # --- индексы/источники/K1L из RingModel (его 4D-кольцо) + матрицы отклика 4D-модели
        model = RingModel(lattice=lattice, nsl=nsl, n_photon=n_photon)
        self.corr = model.corr; self.mon = model.mon; self.quads = model.quads
        self.K1L = np.asarray(model.K1L, float)
        self.hk_ids = list(model.hk_ids); self.vk_ids = list(model.vk_ids)
        self.SIDX = np.asarray(model.photon_rows, int)                 # источники ИД (как в каскаде B)
        self.DVEC = float(model.photon_d)                              # плечо источник→экран, м
        self.bpm = model.bpm                                           # карта имён BPM ↔ индексов орбиты
        self.correctors = model.correctors                             # карта имён корректоров ↔ индексов актуаторов
        self.magnets = model.magnets                                   # карта имён ВСЕХ магнитов ↔ индексов элементов
        if resp_dims == 4:
            self.Hc = model.response(self.hk_ids, ["ebpm.x"])["ebpm.x"]    # 4D-модель корректор→орбита (δ=0)
            self.Hcy = model.response(self.vk_ids, ["ebpm.y"])["ebpm.y"]
        model.close()
        self._col = {p: i for i, p in enumerate(self.hk_ids)}
        self._col.update({p: i for i, p in enumerate(self.vk_ids)})
        self.n_mon = len(self.mon); self.n_corr = len(self.corr); self.n_quad = len(self.quads)
        self.mon_ids = [f"ebpm.{j:03d}" for j in range(self.n_mon)]

        # --- нелинейное кольцо (6D, СИ+ВЧ): динамика пучка
        inf = dict(fr.info_of(lattice)); inf["_id"] = lattice
        self.ring = ring_cpp.CppRing(lattice=lattice, sr=self.sr,
            rf={"volt_mv": inf["rf"][0], "harmonic": inf["rf"][1], "phase_deg": PHASE_DEG},
            design_mev=inf["energy_mev"], nsl=nsl, nsl_dev=nsl_dev)
        if resp_dims == 6:                                    # модель контура снята на 6D-равновесии своего кольца
            self.Hc = self._build_resp_6d("ebpm.x", "hkick")  # (дисперсия/энергия в горизонтали учтены — как LOCO)
            self.Hcy = self._build_resp_6d("ebpm.y", "vkick")
        self.CO = ra.closed_orbit(self.ring, dims=6)
        # --- ВЧ: номинальная частота h·c/C; opt-in актуатор rf.freq (Δf, Гц) + возмущение длины кольца ΔC
        self._rf0 = dict(self.ring.rf)                        # исходный ВЧ (для точного возврата к номиналу)
        self.F0 = float(self.ring.rf_frequency())             # номинал f_RF = h·c/C, Гц
        self.rf_ids = ["rf.freq"] if rf_actuator else []
        self._rf_key = (0.0, 0.0); self.rf_df = 0.0; self.rf_dC = 0.0   # применённые Δf (Гц) и ΔC (м)
        self.TREV = float(inf.get("t_rev", 0.0)) or self._trev_from_ring()
        # ондуляторы: полюса устройства (s-порядок ↔ group_id) — им живой зазор масштабирует изгиб (BSCALE); якорь —
        # выходной полюс (остаточный орбитный кик от неидеального замыкания I1, привязка сцены)
        self._und_poles = dict(id_devices(self.ring.mags))
        self._und_anchor = {gid: poles[-1] for gid, poles in self._und_poles.items()}
        try:
            import undulator_sr as _sr
            self._und_sr = _sr; self._und_Kref = float(_sr.K_of_B0(_sr.B0_of_gap(_sr.P.gap_ref)))
        except Exception:
            self._und_sr = None; self._und_Kref = None

        # --- многоканальное возмущение из disturbance_spec (режим 4): plant владеет Disturbance,
        #     вливаем её в kick_lists каждый шаг. При spec=None — прежний скалярный путь PATX/PATY.
        self.dist_seed = dist_seed
        self.dist = dist.Disturbance(disturbance_spec, self.ring, seed=dist_seed) if disturbance_spec else None
        self._noise_rev = None                               # ревизия ЖИВОГО закона шума (state.noise), подхват на лету
        self._applied = {}                                   # последнее применённое возмущение (offsets/field) — read-only для state.machine
        self._t = 0.0; self._dt = self.stride * self.TREV    # время/шаг для временных законов возмущения
        rb = np.random.default_rng(seed + 4)
        self._bunch0 = np.array([self.CO + np.array(
            [rb.normal(0, 5e-6), 0, rb.normal(0, 5e-6), 0, 0, rb.normal(0, 1e-3)]) for _ in range(n_part)])

        # --- возмущение: фикс. пространственный паттерн feed-down (СКО = dist_amp), амплитуда во времени
        rp = np.random.default_rng(dist_seed)
        self.PATX = rp.normal(0, 1, self.n_quad); self.PATX *= 1.0 / np.std(self.PATX)
        self.PATY = rp.normal(0, 1, self.n_quad); self.PATY *= 1.0 / np.std(self.PATY)
        self.dist_amp = dist_amp; self.dist_rho = dist_rho
        self._amp = dist_amp

        self.bunch = None; self.alive = True; self.turn = 0
        self._last_kick = None                                # кэш применённого набора (пропуск apply_perturb)
        self.centroid_tbt = []                                # пооборотный центроид (x,px,y,py) — фазовый портрет
        self._buf = deque(maxlen=self.avg_chunks)             # окно орбит для усреднения бетатрона (как такт в fofb_sim)
        self._snap_buf = deque(maxlen=self.avg_chunks)        # окно снимков ВСЕЙ линии → state.beam (снимок+среднее)
        # Observer собирает наблюдаемые из усреднённого сырья; дефолт — орбита (+ свет, если есть источники)
        subs = ["ebpm.x", "ebpm.y"] + (["photon.x", "photon.y"] if self.SIDX.size else [])
        self.observer = observer or Observer(subscribe=subs)

        # опорная (golden) орбита на BPM: невозмущённое кольцо (design), замкнутая орбита в мониторах.
        # При relative_orbit показания BPM отдаём ОТНОСИТЕЛЬНО неё → отклонения (≈0 на невозмущённой машине;
        # возмущение/коррекция видны как уход от нуля). НЕ трогает трекинг/карты/токи — только отчёт о показании.
        self.relative_orbit = bool(relative_orbit)
        self.ref_x = np.zeros(self.n_mon); self.ref_y = np.zeros(self.n_mon)
        # center_bunch: центроид пучка — ТОЧНО на замкнутую орбиту (убираем КОГЕРЕНТНЫЙ стартовый сдвиг центроида
        # по энергии ~σ_δ/√N ≈ 350 ppm при N=8 → через дисперсию горизонт. синхротронный переходник ±10 мкм,
        # затухающий ~τ_δ ≈ 4660 об). Разброс по частицам (размер пучка) сохраняется. None → как раньше: только
        # при relative_orbit (кокпит — старт у нуля, без долгого звона); True — и при абсолютных BPM.
        self.center_bunch = self.relative_orbit if center_bunch is None else bool(center_bunch)
        if self.center_bunch:
            self._bunch0 = self._bunch0 - (self._bunch0.mean(axis=0) - self.CO)
        if self.relative_orbit:
            try:                                                  # ТОЙ ЖЕ измерит. цепью, что _measure, но без возмущения:
                self.ring.apply_perturb()                         # design (корректоры=0, возмущение off)
                self._last_kick = None
                rb = self._bunch0.copy()
                for _ in range(max(0, int(round(self.settle_turns / self.stride)))):   # декогеренция бетатрона
                    rb = np.asarray(self.ring.track_tbt(rb.tolist(), self.stride))[:, self.stride, :]
                snaps = [np.mean(np.asarray(self.ring.track_elements(rb.tolist())), axis=0)
                         for _ in range(max(1, self.avg_chunks))]  # усреднение окна ≈ такт (как BPM-показание)
                ref = np.mean(snaps, axis=0)
                self.ref_x = ref[self.mon, 0].copy(); self.ref_y = ref[self.mon, 2].copy()
            except Exception:
                pass

    def _trev_from_ring(self):
        C = float(getattr(self.ring, "circumference", 1110.089))
        return C / 2.99792458e8

    # --- контракт Plant --------------------------------------------------------------------------
    def schema(self):
        return Schema(self.hk_ids + self.vk_ids + self.rf_ids)

    def history(self):
        return History(self.mon_ids, moments=("x", "px", "y", "py"))

    def build_response(self, actuators, observables):
        if list(actuators) == self.rf_ids and self.rf_ids:     # частота ВЧ → наблюдаемые (через Observer)
            return self._build_rf_resp(list(observables))
        cols = [self._col[a] for a in actuators]
        if tuple(observables) == ("ebpm.x",):
            return ResponseMatrix(self.Hc[:, cols], self.mon_ids, actuators)
        if tuple(observables) == ("ebpm.y",):
            return ResponseMatrix(self.Hcy[:, cols], self.mon_ids, actuators)
        raise ValueError(f"TbtPlant: наблюдаемая {observables} не поддержана (линейная модель знает ebpm.x/y)")

    def fingerprint(self):
        return ("tbt", self.resp_dims, self.Hc.shape)

    # --- ВЧ-частота: применение и отклик --------------------------------------------------------------
    def _set_rf(self, df, dC):
        """Резонатор на f = F0 + Δf при длине кольца, растянутой на ΔC. Растяжение = лишний путь ΔC за оборот
        ⇔ часы ВЧ уходят на C_RF − ΔC: f_eff = h·c/(h·c/f − ΔC). Δf = ΔC = 0 ⇒ ровно исходный ВЧ (бит-в-бит)."""
        if df == 0.0 and dC == 0.0:
            self.ring.rf = dict(self._rf0); self.ring._send_rf()
        else:
            hc = float(self._rf0["harmonic"]) * 2.99792458e8
            self.ring.set_rf(freq_hz=hc / (hc / (self.F0 + df) - dC))
        self.rf_df, self.rf_dC = float(df), float(dC)

    def _apply_rf(self, state):
        """Частота ВЧ (актуатор rf.freq) + возмущение длины кольца → резонатор. ΔC зовём РАЗ за такт (walk)."""
        df = float(state.actual["rf.freq"] + state.dynamics["rf.freq"]) if self.rf_ids else 0.0
        dC = self.dist.circumference(self._t, self._dt) if self.dist is not None else 0.0
        if (df, dC) != self._rf_key:
            self._set_rf(df, dC)
            self._rf_key = (df, dC)

    def _build_rf_resp(self, observables, eps=5.0):
        """Отклик наблюдаемых на Δf (центральная разность ±eps Гц) на 6D-равновесии СОБСТВЕННОГО кольца (СИ+ВЧ):
        замкнутая орбита → снимок → та же цепь, что в такте (relative_orbit, Observer из provides() контуров).
        Кэш на диск по контентному отпечатку кольца. Кольцо на входе чистое (корректоры=0, ВЧ номинал)."""
        import hashlib, pickle
        from .observer import ObsContext
        from .response import ring_fingerprint, CACHE_DIR
        self.ring.apply_perturb(); self._last_kick = None        # чистое кольцо (отпечаток без транзитных)
        self._set_rf(0.0, 0.0); self._rf_key = (0.0, 0.0)
        tag = hashlib.md5(repr((ring_fingerprint(self.ring), "rfresp", self.relative_orbit, list(observables),
                                eps)).encode()).hexdigest()[:12]
        path = os.path.join(CACHE_DIR, f"resp_{tag}.pkl")
        if os.path.exists(path):
            try:
                with open(path, "rb") as f:
                    R = pickle.load(f)
                return ResponseMatrix(R, [f"{o}.{i}" for o in observables for i in range(R.shape[0])], self.rf_ids)
            except Exception:
                pass

        def read(df):
            self._set_rf(df, 0.0)
            co = ra.closed_orbit(self.ring, dims=6)
            snap = np.asarray(self.ring.track_elements([co]))[0]
            ex = snap[self.mon, 0] - (self.ref_x if self.relative_orbit else 0.0)
            ey = snap[self.mon, 2] - (self.ref_y if self.relative_orbit else 0.0)
            ctx = ObsContext(None, self, ring=self.ring, intermediates=self.observer.intermediates)
            ctx.set("orbit_x", ex).set("orbit_y", ey)
            out = self.observer.compute(ctx, observables)
            return np.concatenate([np.ravel(np.asarray(out[o], float)) for o in observables])

        try:
            R = ((read(+eps) - read(-eps)) / (2.0 * eps)).reshape(-1, 1)
        finally:
            self._set_rf(0.0, 0.0); self._rf_key = (0.0, 0.0)
        os.makedirs(CACHE_DIR, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(R, f)
        return ResponseMatrix(R, [f"{o}.{i}" for o in observables for i in range(R.shape[0])], self.rf_ids)

    def _build_resp_6d(self, obs, kind):
        """Матрица отклика корректор→орбита, снятая с СОБСТВЕННОГО 6D-кольца (СИ+ВЧ, dims=6): дисперсия и
        энергия учтены — это модель «как на машине» (аналог измеренной LOCO). Кэш по контентному отпечатку
        кольца (опция 5) — строится раз. Кольцо на входе чистое (вызов в __init__ до возмущений)."""
        from . import response as rsp
        acts = self.hk_ids if kind == "hkick" else self.vk_ids
        idx = {a: int(self.corr[j]) for j, a in enumerate(acts)}
        xcol = 0 if obs.endswith(".x") else 2

        def applier(ring, aid, eps):
            ring.apply_perturb(**{("kicks" if kind == "hkick" else "vkicks"): [(idx[aid], eps)]})

        def reader(_oid, orbit):
            return orbit[self.mon, xcol]

        key = (self.lattice, self.nsl, "tbt6d", kind)
        return rsp.cached_build(key, self.ring, list(acts), [obs], applier, reader, dims=6)[obs]

    # --- живой подхват настроек из доски (бесконечный режим) -------------------------------------
    def _sync_noise(self, state):
        """Подхватить ЖИВОЙ закон шума из state.noise: при смене ревизии ПЕРЕСОБРАТЬ Disturbance.
        state.noise=None → оставить возмущение из конструктора (обратная совместимость)."""
        nm = getattr(state, "noise", None)
        if nm is None or nm.rev == self._noise_rev:
            return
        self.dist = dist.Disturbance(nm.spec, self.ring, seed=self.dist_seed) if nm.spec else None
        self._noise_rev = nm.rev
        self._last_kick = None                               # закон изменился — форсируем повторный apply_perturb

    def _meas_noise(self, state):
        """СКО шума измерения e-BPM: из живых настроек state.sensors['ebpm'], иначе дефолт plant'а."""
        sc = state.sensors.get("ebpm") if getattr(state, "sensors", None) else None
        if sc is not None and sc.meas_noise is not None:
            return sc.meas_noise
        return self.meas_noise

    def set_drive_speed(self, state, corrector, slew, plane="both"):
        """ЖИВО изменить скорость привода (slew, ед/с) КОНКРЕТНОГО корректора — подхват на следующем такте
        (Scheduler зовёт state.limits.apply каждый такт, читая slew заново; кэша нет — rev не нужен).
        corrector: имя SP03.COR05 | sensor_id | индекс (self.correctors). plane: 'x'/'h', 'y'/'v', 'both'.
        slew=inf → мгновенный привод (без ограничения). Возврат — индекс корректора."""
        i = self.correctors.index(corrector)
        if plane in ("x", "h", "both"):
            state.limits.set(f"corr.{i:03d}.hkick", slew=slew)
        if plane in ("y", "v", "both"):
            state.limits.set(f"corr.{i:03d}.vkick", slew=slew)
        return i

    # --- возмущение + актуация -------------------------------------------------------------------
    def _advance_dist(self):
        if self.dist_rho < 1.0:                               # медленный дрейф амплитуды (AR(1))
            a = np.sqrt(1.0 - self.dist_rho ** 2)
            self._amp = self.dist_rho * self._amp + a * self.dist_amp * self._rng.standard_normal()
        # dist_rho == 1 → статичная амплитуда (статическое искажение орбиты)

    def kick_lists(self, state):
        """Полный ТЕКУЩИЙ набор возмущений машины (корректоры + впрыск + feed-down/поле возмущений) как
        аргументы apply_perturb: (kk, vv, extra), где extra — dict прочих каналов (dks/dk2s/dk3s/dsks).
        Это 'конфигурация машины' для анализа (режим 3): тот же набор можно отдать построителю 4D-кольца.

        Возмущение: при spec=None — прежний фикс. паттерн feed-down квадруполей (скалярная амплитуда _amp);
        при заданном disturbance_spec — многоканальная Disturbance (offsets→feed-down квадруполей,
        field→прямые каналы apply_perturb), вычисляется во времени (self._t, self._dt)."""
        th = state.actual.select(self.hk_ids) + state.dynamics.select(self.hk_ids)
        tv = state.actual.select(self.vk_ids) + state.dynamics.select(self.vk_ids)
        kk = [(self.corr[j], float(th[j])) for j in range(self.n_corr)]
        vv = [(self.corr[j], float(tv[j])) for j in range(self.n_corr)]
        extra = {}
        if self.dist is None:                                 # скалярный путь (дефолт)
            dkx = self.K1L * (self.PATX * self._amp); dky = -self.K1L * (self.PATY * self._amp)
            kk += [(self.quads[j], float(dkx[j])) for j in range(self.n_quad)]
            vv += [(self.quads[j], float(dky[j])) for j in range(self.n_quad)]
            self._applied = {}
        else:                                                 # многоканальная Disturbance (режим 4)
            dx, dy = self.dist.offsets(self._t, self._dt)     # сдвиги магнитов [м] → feed-down квадруполей
            fld = self.dist.field(self._t, self._dt)          # рябь поля (вызвать РАЗ: walk продвигает интегратор)
            kk += [(self.quads[j], float(self.K1L[j] * dx[self.quads[j]])) for j in range(self.n_quad)]
            vv += [(self.quads[j], float(-self.K1L[j] * dy[self.quads[j]])) for j in range(self.n_quad)]
            for key, vec in fld.items():                      # рябь поля: прямые каналы apply_perturb
                rows = [(int(i), float(vec[i])) for i in np.nonzero(vec)[0]]
                if key == "kicks":
                    kk += rows
                elif key == "vkicks":
                    vv += rows
                else:
                    extra[key] = extra.get(key, []) + rows
            self._applied = {"dx": dx, "dy": dy, "field": fld}   # read-only снимок возмущения для state.machine
        for gid, idx in self._und_anchor.items():             # остаточный кик орбиты от зазора ондулятора (неид. замыкание)
            u = state.undulators.get(gid)
            if u is not None:
                k = u.orbit_kick_rad()
                if k:
                    vv.append((idx, float(k)))
        bs = self._und_bscale(state)                          # поле полюсов ИД от зазора → орбита-качка, фокусировка, СИ
        if bs:
            extra["bscales"] = bs
        return kk, vv, extra

    UND_F_QUANT = 5e-4     # шаг квантования f = K/K_ref: фокусировка ∝ f² ⇒ Δνy ≤ ~1e-7 на ИД

    def _und_bscale(self, state):
        """Масштаб изгиба полюсов ИД от живых зазоров: [(индекс полюса, f)], f = K/K_ref квантованный UND_F_QUANT.
        Поле полюсов ∝ K ⇒ угол и грани θ/2 ∝ f (рамки заново на сервере): вертикальная фокусировка граней ∝ f²,
        горизонталь — дрейф. На опорном зазоре (сборки, 5 мм) — пусто ⇒ прогоны без смены зазора как раньше."""
        if not self._und_Kref or not getattr(state, "undulators", None):
            return []
        out = []
        q = self.UND_F_QUANT
        for gid, poles in self._und_poles.items():
            u = state.undulators.get(gid)
            if u is None:
                continue
            f = self._und_sr.K_of_B0(self._und_sr.B0_of_gap(float(u.gap_mm))) / self._und_Kref
            f = round(f / q) * q
            if abs(f - 1.0) > 0.5 * q:
                out.extend((j, f) for j in poles)
        return out

    def _apply(self, state):
        """Применить текущий набор возмущений к кольцу. Полный apply_perturb дорог (RESET + набор, ~0,3 с) —
        зовём его ТОЛЬКО когда изменились «прочие» каналы (dks/dk2s/dk3s/dsks — ряби поля, feed-down через
        градиент). Кики KICK/VKICK (корректоры FOFB, раз в цикл) «задают», а не прибавляют, — шлём только
        изменившиеся (set_kicks, без RESET). Масштаб полюсов ИД от зазора (bscales) меняется почти каждый такт,
        пока едут челюсти, — шлём только изменившиеся полюса (set_bend_scale). Оба дешёвых пути — миллисекунды."""
        kk, vv, extra = self.kick_lists(state)
        self._apply_rf(state)                                 # частота ВЧ + длина кольца (дёшево: команда RF)
        bs = dict(extra.pop("bscales", []))                   # масштаб изгиба полюсов ИД (set-семантика, база 1)
        kd, vd = dict(kk), dict(vv)                           # последний выигрывает — как при replay сервера
        xkey = tuple(sorted((k, tuple(v)) for k, v in extra.items()))
        last = self._last_kick
        if last is None or last[2] != xkey:
            self.ring.apply_perturb(kicks=kk, vkicks=vv, bscales=sorted(bs.items()), **extra)
            self._last_kick = (kd, vd, xkey); self._bs_applied = bs
            return
        okd, ovd, _ = last
        ch_h = [(j, kd.get(j, 0.0)) for j in sorted(set(kd) | set(okd)) if kd.get(j, 0.0) != okd.get(j, 0.0)]
        ch_v = [(j, vd.get(j, 0.0)) for j in sorted(set(vd) | set(ovd)) if vd.get(j, 0.0) != ovd.get(j, 0.0)]
        if ch_h or ch_v:
            self.ring.set_kicks(kicks=ch_h, vkicks=ch_v)
            self._last_kick = (kd, vd, xkey)
        old = getattr(self, "_bs_applied", {})
        ch = [(j, bs.get(j, 1.0)) for j in sorted(set(bs) | set(old)) if bs.get(j, 1.0) != old.get(j, 1.0)]
        if ch:
            self.ring.set_bend_scale(ch)
            self._bs_applied = bs

    def _snapshot(self):
        """Центроид пучка вдоль линии (track_elements, read-only — выход НЕ переносим)."""
        return np.mean(np.asarray(self.ring.track_elements(self.bunch.tolist())), axis=0)

    def _measure(self, state, snap):
        self._snap_buf.append(snap)                           # окно снимков ВСЕЙ линии → state.beam (снимок+среднее)
        # history — МГНОВЕННЫЙ истинный снимок пучка (бетатрон виден): фазовый портрет и плоскостные датчики
        state.history.data[:, 0] = snap[self.mon, 0]; state.history.data[:, 1] = snap[self.mon, 1]
        state.history.data[:, 2] = snap[self.mon, 2]; state.history.data[:, 3] = snap[self.mon, 3]
        # BPM-показание — орбита (+ электронный шум), УСРЕДНЁННАЯ по окну ≈ такт: бетатрон и шум сбиваются
        ex = snap[self.mon, 0].copy(); ey = snap[self.mon, 2].copy()
        if self.relative_orbit:                               # показания ОТНОСИТЕЛЬНО опорной (golden) орбиты
            ex -= self.ref_x; ey -= self.ref_y
        mn = self._meas_noise(state)                          # живой шум датчика из state.sensors, иначе дефолт
        if mn:
            ex += mn * self._rng.standard_normal(self.n_mon)
            ey += mn * self._rng.standard_normal(self.n_mon)
        phx = snap[self.SIDX, 0] + snap[self.SIDX, 1] * self.DVEC if self.SIDX.size else None
        phy = snap[self.SIDX, 2] + snap[self.SIDX, 3] * self.DVEC if self.SIDX.size else None
        self._buf.append((ex, ey, phx, phy))
        # усреднённое сырьё (сенсорная модель выше) → в ObsContext; наблюдаемые собирает Observer
        ctx = ObsContext(state, self, ring=self.ring, intermediates=self.observer.intermediates)
        ctx.set("orbit_x", np.mean([e[0] for e in self._buf], axis=0))
        ctx.set("orbit_y", np.mean([e[1] for e in self._buf], axis=0))
        if self.SIDX.size:
            ctx.set("photon_x", np.mean([e[2] for e in self._buf], axis=0))
            ctx.set("photon_y", np.mean([e[3] for e in self._buf], axis=0))
        self.observer.update(state, ctx)

    def _track(self):
        """Провести пучок stride оборотов; накопить пооборотный центроид (фазовый портрет)."""
        adv = np.asarray(self.ring.track_tbt(self.bunch.tolist(), self.stride))   # (NP, stride+1, 6)
        self.bunch = adv[:, self.stride, :]
        cen = np.mean(adv[:, 1:, :], axis=0)                  # центроид по частицам, обороты 1..stride
        self.centroid_tbt.extend(cen[:, [0, 1, 2, 3]].tolist())
        self.turn += self.stride
        if not np.all(np.isfinite(self.bunch)):
            self.alive = False

    # --- такты -----------------------------------------------------------------------------------
    def reset(self, state):
        self.bunch = self._bunch0.copy(); self.alive = True; self.turn = 0
        self._last_kick = None; self.centroid_tbt = []; self._amp = self.dist_amp
        self._buf.clear(); self._snap_buf.clear(); self._t = 0.0   # возмущение держим в t=0 на время декогеренции
        self._sync_noise(state)                               # подхватить закон шума из доски (state.noise)
        self._apply(state)                                    # корректоры=0 + возмущение
        n_settle = max(0, int(round(self.settle_turns / self.stride)))
        for _ in range(n_settle):                             # декогеренция бетатрона к искажённой орбите
            self._track()
        self._measure(state, self._snapshot())
        state.beam = BeamSnapshot(self._snap_buf, self.magnets)    # пучок по всей линии (снимок+среднее), окно обновляется на месте
        state.machine = MachineView(self)                          # lazy-доступ к положению/полю магнитов
        self._seed_undulators(state)                               # живые зазоры ондуляторов на доску (видят контуры)
        try:                                                       # линейная оптика на доску (β/η/тюн) — статический снимок
            state.optics = OpticsView(ra.twiss(self.ring), s_out=np.asarray(self.ring.s_out, float), magnets=self.magnets)
        except Exception:
            pass

    def _seed_undulators(self, state):
        """Заполнить state.undulators живыми зазорами по group_id ондуляторов (дефолт self.gap_mm).
        Контуры видят их через StateView.undulators; зазор правится на лету и читается рекордерами."""
        try:
            import errors as er
            for g in er.load_groups(self.lattice):
                gid = g.get("group_id")
                if gid and gid not in state.undulators:
                    state.undulators[gid] = UndulatorConfig(self.gap_mm, name=gid)
        except Exception:
            pass

    def step(self, state, dt):
        self._t = float(state.t); self._dt = float(dt)        # время/шаг для временных законов Disturbance
        self._sync_noise(state)                               # живой подхват закона шума (бесконечный режим)
        for c in state.undulators.values():                   # привод зазора: реализованный зазор едет к уставке ≤ slew
            c.advance(dt)
        self._advance_dist()
        self._apply(state)                                    # перезовёт apply_perturb лишь при смене киков
        self._measure(state, self._snapshot())                # сенсоры ПЕРЕД продвижением (задержка в такт)
        self._track()

    def close(self):
        try:
            self.ring.apply_perturb(); self.ring.close()
        except Exception:
            pass


class MachineView:
    """Lazy-доступ к состоянию МАГНИТОВ для контура: позиция и поле на ТЕКУЩИЙ такт, по имени/индексу.
    Считает ТОЛЬКО запрошенное из живых данных plant'а (кольцо + последнее применённое возмущение) —
    ничего не материализуется заранее. Адресация через карту plant.magnets (SP<nn>.<dev>)."""

    _FKEY = {"theta0": "kicks", "k1": "dks", "k2": "dk2s", "k3": "dk3s"}   # поле → канал возмущения

    def __init__(self, plant):
        self._p = plant

    def index(self, magnet):
        return self._p.magnets.index(magnet)

    def position(self, magnet):
        """Положение магнита (x, y, z), м: номинальное выравнивание + динамический сдвиг возмущения."""
        i = self.index(magnet)
        c = [float(v) for v in getattr(self._p.ring.mags[i], "center", [0.0, 0.0, 0.0])]
        ap = self._p._applied
        if ap:
            c[0] += float(ap["dx"][i]); c[1] += float(ap["dy"][i])
        return tuple(c)

    def field(self, magnet):
        """Силы магнита на текущий такт: {theta0,k1,k2,k3} = номинал + применённая Δ (рябь поля возмущения)."""
        i = self.index(magnet); m = self._p.ring.mags[i]
        out = {k: float(getattr(m, k, 0.0) or 0.0) for k in self._FKEY}
        fld = (self._p._applied or {}).get("field", {})
        for k, chan in self._FKEY.items():
            if chan in fld:
                out[k] += float(fld[chan][i])
        return out

    def orientation(self, magnet):
        """Углы ориентации прибора (yaw, pitch, roll), рад, из решётки: yaw — азимут (куда смотрит элемент
        в глобальной раме, накапливается по кольцу), pitch/tilt — наклоны (0 у идеальной выставки, ненулевые
        при разъюстировке). Базовое возмущение правит только dx/dy/поле, не ориентацию, — поэтому статично."""
        i = self.index(magnet); m = self._p.ring.mags[i]
        return {k: float(getattr(m, k, 0.0) or 0.0) for k in ("yaw", "pitch", "roll")}
