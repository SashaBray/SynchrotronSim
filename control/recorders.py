# -*- coding: utf-8 -*-
"""Рекордеры режима 4 — РЕАЛИЗАЦИЯ ПОДПИСКИ на показания (`fb_record.csv`).

Строка подписки (signal, params, cadence, enabled) → объект-рекордер. Рантайм на каждом окне зовёт
`sample(state, widx, plant)` (дешёвые сигналы читаются из `state.observables`), в конце — `finalize(...)`
(тяжёлые/производные диагностики, напр. интерференционная картина — считается ТОЛЬКО если подписан).
Не подписался (enabled=0 или неизвестный signal) — рекордер не инстанцируется, диагностика не считается.

Выход рекордера: `result_rows()` → [(signal, step, value)] в `fb_results.csv`; `artifact()` →
(signal, {имя: ndarray}, meta) в `fb_artifacts.csv`+`.npz`. Дешёвые: orbit_rms, light, bpm (из
observables). Тяжёлые: phase (фазовый портрет), screen (интерференция — клон `interference_compare`).
tune/chroma/sf в v1 не реализованы (их сигнал просто пропускается — см. REGISTRY)."""
import csv
import math
import os

import numpy as np

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

# параметры ондулятора/оптики для screen-рекордера (как в interference_compare, ПРЕДПОЛОЖИТЕЛЬНЫЕ)
GAP_MM = 5.0
BETA_X, BETA_Y = 4.0, 2.0
D_DEFAULT = 15.0                                              # плечо источник→экран, м (если нет в screens.csv)


def _rms(v):
    return float(np.sqrt(np.mean(np.asarray(v, float) ** 2))) if v is not None and len(np.ravel(v)) else 0.0


class Recorder:
    """База. По умолчанию ничего не делает; наследники переопределяют sample/finalize."""

    def __init__(self, params="", cadence="per_window", **kw):
        self.params = params
        self.cadence = cadence
        self.rows = []                                       # (signal, step, value)
        self._art = None                                     # (signal, {name: ndarray}, meta) или None

    def sample(self, state, widx, plant):
        pass

    def finalize(self, state, plant):
        pass

    def result_rows(self):
        return self.rows

    def artifact(self):
        return self._art


class OrbitRmsRecorder(Recorder):
    """СКО орбиты в мониторах по окнам: orbit_rms.x / orbit_rms.y (мкм)."""

    def sample(self, state, widx, plant):
        self.rows.append(("orbit_rms.x", widx, _rms(state.observables.get("ebpm.x")) * 1e6))
        self.rows.append(("orbit_rms.y", widx, _rms(state.observables.get("ebpm.y")) * 1e6))


class LightRecorder(Recorder):
    """СКО центроида фотонного пучка в источниках: light.x / light.y (мкм). Нет photon → 0."""

    def sample(self, state, widx, plant):
        self.rows.append(("light.x", widx, _rms(state.observables.get("photon.x")) * 1e6))
        self.rows.append(("light.y", widx, _rms(state.observables.get("photon.y")) * 1e6))


class BpmRecorder(Recorder):
    """Орбита по мониторам: финальный снимок x/y как артефакт (per-sensor массивы)."""

    def finalize(self, state, plant):
        ex = np.asarray(state.observables.get("ebpm.x", []), float)
        ey = np.asarray(state.observables.get("ebpm.y", []), float)
        self._art = ("bpm", {"x": ex, "y": ey}, f"мониторов={ex.size}")


class PhaseRecorder(Recorder):
    """Фазовый портрет центроида (x, x′, y, y′) по оборотам → артефакт."""

    def finalize(self, state, plant):
        cen = np.asarray(getattr(plant, "centroid_tbt", []), float)
        if cen.size:
            self._art = ("phase", {"centroid": cen}, f"оборотов={cen.shape[0]}")


class ScreenRecorder(Recorder):
    """Интерференционная картина на экране за ондулятором — ТЯЖЁЛАЯ опциональная диагностика.

    Клон физики `control/interference_compare.py`, но зазор теперь ЖИВОЙ: каждое окно читаем зазор своего
    ондулятора из `state.undulators` (источник экрана — колонка `source` в `screens.csv`), поэтому смена
    зазора в прогоне СДВИГАЕТ картину (K→E₁→конус/поток). Угловой разброс пучка берём из живого прогона
    (блуждание пятна photon.x/y по окнам, σ/D = смаз угла), складываем в квадратуре к расходимости и
    размываем `screen_map`. В артефакт пишем картины при НАЧАЛЬНОМ и КОНЕЧНОМ зазоре + измеренную, в
    результаты — ряд E₁(окно) и резкость. Нет живого зазора на доске → fallback на зашитый GAP_MM.
    УПРОЩЕНИЕ v1: картина на оси, источник = photon[0]; геометрия экрана из `screens.csv` даёт плечо D."""

    def __init__(self, params="", cadence="final", mesh=121, **kw):
        super().__init__(params, cadence)
        self.mesh = int(mesh)
        self._spot = []
        self._gapser = []                                    # (widx, gap_mm) по окнам — как зазор менялся
        self._srcname = None                                 # кэш имени ондулятора-источника экрана

    def sample(self, state, widx, plant):
        px = state.observables.get("photon.x"); py = state.observables.get("photon.y")
        if px is not None and py is not None and len(np.ravel(px)):
            self._spot.append((float(np.ravel(px)[0]), float(np.ravel(py)[0])))
        self._gapser.append((int(widx), self._live_gap(state)))

    def _screen_row(self, screen_id):
        path = os.path.join(_ROOT, "workspace", "tables", "screens.csv")
        if screen_id and os.path.exists(path):
            with open(path, encoding="utf-8-sig") as f:
                for r in csv.DictReader(f, delimiter=";"):
                    if (r.get("screen_id") or "").strip() == screen_id:
                        return r
        return {}

    def _screen_dist(self, screen_id):
        try:
            return float(self._screen_row(screen_id).get("dist_m") or D_DEFAULT)
        except (ValueError, TypeError):
            return D_DEFAULT

    def _screen_source(self, screen_id):
        if self._srcname is None:                            # замороженная screens.csv: source_ref; манифест: source
            r = self._screen_row(screen_id)
            self._srcname = (r.get("source_ref") or r.get("source") or "").strip()
        return self._srcname

    def _live_gap(self, state):
        """Зазор ондулятора-источника экрана из живой доски; нет записи → зашитый GAP_MM."""
        src = self._screen_source(self.params)
        u = (getattr(state, "undulators", None) or {}).get(src)
        return float(u.gap_mm) if u is not None else GAP_MM

    def finalize(self, state, plant):
        import sys
        sys.path.insert(0, os.path.join(_ROOT, "undulators"))
        import undulator_sr as sr
        from undulator_field import trajectory, calibrate
        from screen_spectrum import beam_params, screen_map

        D = self._screen_dist(self.params); g = sr.gamma()
        gap_start = self._gapser[0][1] if self._gapser else self._live_gap(state)
        gap_end = self._gapser[-1][1] if self._gapser else gap_start

        spot = np.asarray(self._spot, float) if self._spot else np.zeros((1, 2))
        jx = float(np.std(spot[:, 0])) / D; jy = float(np.std(spot[:, 1])) / D   # смаз угла из блуждания пятна
        axis = np.array([0.0, 0.0, 1.0]); uhat = np.array([1.0, 0.0, 0.0]); vhat = np.array([0.0, 1.0, 0.0])
        P = np.zeros(3)
        beam0 = beam_params(sr.P, BETA_X, BETA_Y, sr.P.kappa, sr.P.sig_E)

        def e1_of(gap_mm):
            return float(sr.harmonic_energy_keV(1, sr.K_of_B0(sr.B0_of_gap(gap_mm)), 0.0))

        def pattern_at(gap_mm, ax_, ay_):
            K = sr.K_of_B0(sr.B0_of_gap(gap_mm))
            tr = trajectory(K, sr.P.lam_u, sr.P.N, g, 64); scale, _ = calibrate(tr)
            E0 = float(sr.harmonic_energy_keV(1, K, 0.0)); srp = float(sr.sigma_photon(E0)[1])
            window = max(8.0 * srp, 4.0 * math.hypot(ax_, ay_) + 1e-12)
            b = dict(beam0); b["sxp"] = math.hypot(beam0["sxp"], ax_); b["syp"] = math.hypot(beam0["syp"], ay_)
            mp = screen_map(tr, scale, P, (axis, uhat, vhat), axis, D, np.array([E0]), self.mesh, window, b)
            I = mp["I"][..., 0] if np.asarray(mp["I"]).ndim == 3 else mp["I"]
            return np.asarray(I, float), float(mp["peak"]), np.asarray(mp["extent"], float), E0, window

        Istart, pk_s, ext, E1_s, _ = pattern_at(gap_start, 0.0, 0.0)     # идеал при НАЧАЛЬНОМ зазоре
        Iend, pk_e, ext_e, E1_e, win_e = pattern_at(gap_end, 0.0, 0.0)   # идеал при КОНЕЧНОМ зазоре
        Imeas, pk_m, _, _, _ = pattern_at(gap_end, jx, jy)               # измеренный (джиттер) при конечном зазоре

        series = np.asarray([[w, gp, e1_of(gp)] for w, gp in self._gapser], float) if self._gapser \
            else np.zeros((0, 3))
        self._art = ("screen", {"ideal_start": Istart, "ideal_end": Iend, "measured": Imeas,
                                "extent": ext_e, "jitter_urad": np.array([jx * 1e6, jy * 1e6]),
                                "D_m": np.array([D]), "gap_start_mm": np.array([gap_start]),
                                "gap_end_mm": np.array([gap_end]), "energy_start_keV": np.array([E1_s]),
                                "energy_end_keV": np.array([E1_e]), "gap_series": series,
                                "window_mm": np.array([win_e * D * 1e3])},
                     f"screen={self.params} зазор {gap_start:.1f}->{gap_end:.1f}мм E1 {E1_s:.2f}->{E1_e:.2f}кэВ")
        self.rows += [("screen.gap_start_mm", 0, gap_start), ("screen.gap_end_mm", 0, gap_end),
                      ("screen.E1_start_keV", 0, E1_s), ("screen.E1_end_keV", 0, E1_e),
                      ("screen.peak_meas", 0, pk_m), ("screen.jitter_x_urad", 0, jx * 1e6),
                      ("screen.jitter_y_urad", 0, jy * 1e6),
                      ("screen.sharpness_pct", 0, (pk_m / pk_e * 100.0) if pk_e else 0.0)]
        for w, gp in self._gapser:                                       # ряд: как E₁ и зазор едут по окнам
            self.rows.append(("screen.gap_mm", int(w), float(gp)))
            self.rows.append(("screen.E1_keV", int(w), e1_of(gp)))


class VectorRecorder(Recorder):
    """Компоненты ВЕКТОРНОЙ наблюдаемой (tune/chroma/adts с OpticsPlant) по окнам: строки <signal>.<label>."""

    SIGNAL = None
    LABELS = None

    def sample(self, state, widx, plant):
        v = state.observables.get(self.SIGNAL)
        if v is None:
            return
        v = np.asarray(v, float).ravel()
        labels = self.LABELS or [str(i) for i in range(v.size)]
        for lab, val in zip(labels, v):
            self.rows.append((f"{self.SIGNAL}.{lab}", widx, float(val)))


class TuneRecorder(VectorRecorder):
    SIGNAL = "tune"; LABELS = ["x", "y"]                  # νx, νy (дробные)


class ChromaRecorder(VectorRecorder):
    SIGNAL = "chroma"; LABELS = ["x", "y"]               # ξx, ξy


class AdtsRecorder(VectorRecorder):
    SIGNAL = "adts"; LABELS = ["xx", "yy", "xy"]         # ∂νx/∂Jx, ∂νy/∂Jy, ∂νx/∂Jy


REGISTRY = {
    "orbit_rms": OrbitRmsRecorder,
    "light": LightRecorder,
    "bpm": BpmRecorder,
    "phase": PhaseRecorder,
    "screen": ScreenRecorder,
    "tune": TuneRecorder,
    "chroma": ChromaRecorder,
    "adts": AdtsRecorder,
    # sf — зарезервирован (пост-анализ S(f)), в v1 не реализован
}


def build(records, mesh=121):
    """records = list[fb_config.Record] → list рекордеров (только enabled и известные signal)."""
    out = []
    for r in records:
        if not r.enabled:
            continue
        cls = REGISTRY.get(r.signal)
        if cls is None:                                      # неизвестный/нереализованный сигнал — пропускаем
            continue
        out.append(cls(params=r.params, cadence=r.cadence, mesh=mesh))
    return out
