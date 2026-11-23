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
import ring_analysis as ra

from .state import Schema, History, UndulatorConfig, OpticsView, _undulator_sr
from .controller import ResponseMatrix
from .ring_model import RingModel
from .observer import Observer, ObsContext

FD_EPS = 1e-3                                                  # шаг конечной разности матриц отклика (абс. Δk)


class OpticsPlant:
    def __init__(self, model=None, lattice=26, nsl=8, observe=("tune", "chroma", "adts"),
                 drift_k1=0.0, drift_k2=0.0, drift_k3=0.0, rho=1.0, seed=1, order=4, gap_mm=5.0, observer=None):
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
        # Observer собирает наблюдаемые из возмущённого кольца. Схему вычисления tune/chroma/adts несут
        # ФАЙЛЫ-КОНТУРЫ (модульность): дефолтный observer собираем из их `provides()`, а не из observer.py.
        if observer is None:
            from .obs_assembly import observer_for
            from .controllers.optics_tune import TuneLoop
            from .controllers.optics_chroma import ChromaLoop
            from .controllers.optics_adts import AdtsLoop
            observer = observer_for([TuneLoop(), ChromaLoop(), AdtsLoop()], serves=set(self.observe))
        self.observer = observer
        # зазор ондулятора → вертикальный фокус UFC (k1v ∝ K²) → отклик ТЮНА; UFC↔ондулятор по s-порядку
        self.gap_mm = float(gap_mm); self.gap_ref = 5.0
        _sr = _undulator_sr()
        self._K_ref = float(_sr.K_of_B0(_sr.B0_of_gap(self.gap_ref))) if _sr is not None else None
        _ufc_idx = sorted(i for i, m in enumerate(mags) if str(getattr(m, "dev", "") or "") == "und_UFC")
        self._ufc = [(i, float(getattr(mags[i], "k1", 0.0) or 0.0), float(getattr(mags[i], "k1v", 0.0) or 0.0),
                      f"und_{k + 1:02d}") for k, i in enumerate(_ufc_idx)]

    # ---- аналитическая оценка оптики на возмущённом кольце ----
    def _perturbed(self, pk1, pk2, pk3, gaps=None, ufc_fac=None):
        """Копия аналитического кольца с Δk1/Δk2/Δk3 (pk*), фокусом UFC по живым зазорам (gaps: group_id → мм)
        и/или прямым множителем фокуса ufc_fac (group_id → (K/K_ref)²; имеет приоритет над gaps)."""
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
        # зазор ондулятора → ВЕРТИКАЛЬНЫЙ фокус UFC: k1v ∝ K² (естественный фокус планарного ИД k_y = ½(B0/Bρ)²).
        # Горизонтальный k1 UFC НЕ трогаем: он гасит слабую фокусировку 1/ρ² секторных полюсов, а полюса в
        # кольце от живого зазора не зависят (их углы заданы при сборке). Масштаб k1 ломал бы гашение → ложный
        # Δνx (+0,013 при раскрытии всех 38 ИД 5→20 мм), тогда как у планарного ИД горизонталь — дрейф при любом зазоре.
        if self._ufc and (ufc_fac or (gaps and self._K_ref)):
            for i, k1n, k1vn, gid in self._ufc:
                if ufc_fac and gid in ufc_fac:
                    fac = float(ufc_fac[gid])
                else:
                    g = gaps.get(gid) if gaps else None
                    if g is None or not self._K_ref or abs(g - self.gap_ref) < 1e-9:
                        continue
                    fac = self.gap_fac(g)
                m = _edit(i); m.k1v = k1vn * fac
        for i, m in touched.items():
            mags[i] = m
        pr.mags = mags
        return pr

    def gap_fac(self, gap_mm):
        """Множитель вертикального фокуса UFC при зазоре gap_mm относительно опорного (K/K_ref)². 1 — без модели СИ."""
        _sr = _undulator_sr()
        if _sr is None or not self._K_ref:
            return 1.0
        return float(_sr.K_of_B0(_sr.B0_of_gap(float(gap_mm))) / self._K_ref) ** 2

    GAP_OBS = ("tune", "chroma")

    def gap_response(self, observables=GAP_OBS, order=2, only=None, rebuild=False):
        """Линейный отклик наблюдаемых на фокус ондуляторов: (group_id по s-порядку, {obs: R k×n_ИД}),
        Δobs = R·(fac − 1), fac = gap_fac(зазор) = (K/K_ref)². Столбец ИД — конечная разность «вертикальный
        фокус ИД снят» (fac 1 → 0): отклик линеен по fac (на единицу (1 − fac) одинаков при 10 и 50 мм).
        Хроматизм идёт в основном через бета-биение в секступолях (у СИЛА-26 ×9 больше прямого вклада UFC),
        поэтому формулой −β·k·L/4π его не взять — только разностью. Порядок 2 точен для ν и ξ₁.
        only — подмножество group_id (без кэша, для проверок). Полный набор кэшируется на диск по отпечатку
        кольца: первый раз n_ИД + 1 проход ring_map (СИЛА-26: 39 × ~3 с ≈ 2 мин), дальше — мгновенно."""
        import hashlib, os, pickle
        from .response import ring_fingerprint, CACHE_DIR
        obs = [o for o in observables if o in self.observe]
        gids = [gid for *_, gid in self._ufc if only is None or gid in only]
        if not gids or not obs:
            return [], {}
        path = None
        if only is None:
            key = (ring_fingerprint(self.rm.ring), "gapmodel", int(order), tuple(obs),
                   tuple((int(i), float(k1vn), gid) for i, _, k1vn, gid in self._ufc))
            path = os.path.join(CACHE_DIR, f"gapmodel_{hashlib.md5(repr(key).encode()).hexdigest()[:12]}.pkl")
            if not rebuild and os.path.exists(path):
                try:
                    with open(path, "rb") as f:
                        return pickle.load(f)
                except Exception:
                    pass

        def run(fac):
            ctx = ObsContext(None, self, ring=self._perturbed({}, {}, {}, ufc_fac=fac), order=int(order),
                             intermediates=self.observer.intermediates)
            return self.observer.compute(ctx, obs)

        base = run(None)
        cols = {o: [] for o in obs}
        for gid in gids:
            d = run({gid: 0.0})                                   # фокус ИД снят: fac − 1 = −1
            for o in obs:
                cols[o].append(-(np.asarray(d[o], float) - np.asarray(base[o], float)))
        out = (gids, {o: np.column_stack(cols[o]) for o in obs})
        if path:
            os.makedirs(CACHE_DIR, exist_ok=True)
            with open(path, "wb") as f:
                pickle.dump(out, f)
        return out

    def _ctx(self, pk1, pk2, pk3, state=None):
        """Контекст наблюдателя на ВОЗМУЩЁННОМ кольце (общие промежуточные: карта оборота считается раз)."""
        gaps = ({gid: c.gap_mm for gid, c in state.undulators.items()}
                if state is not None and state.undulators else None)      # живые зазоры → фокус UFC → тюн
        return ObsContext(state, self, ring=self._perturbed(pk1, pk2, pk3, gaps), order=self.order,
                          intermediates=self.observer.intermediates)       # билдеры optics-промежуточных — из контуров

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
        self._seed_undulators(state)
        try:                                                  # линейная оптика на доску (β/η/тюн) — статический снимок
            state.optics = OpticsView(ra.twiss(self.rm.ring), s_out=np.asarray(self.rm.ring.s_out, float),
                                      magnets=self.rm.magnets)
        except Exception:
            pass
        self._a = np.zeros(3); self._advance(); self._measure(state)

    def _seed_undulators(self, state):
        """Живые зазоры ондуляторов на доску (по group_id), дефолт self.gap_mm — видят контуры."""
        try:
            import errors as er
            for g in er.load_groups(self.rm.lattice):
                gid = g.get("group_id")
                if gid and gid not in state.undulators:
                    state.undulators[gid] = UndulatorConfig(self.gap_mm, name=gid)
        except Exception:
            pass

    def step(self, state, dt):
        for c in state.undulators.values():                   # привод зазора ондулятора (ограничение скорости челюстей)
            c.advance(dt)
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
