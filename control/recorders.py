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

    Клон физики `control/interference_compare.py`: чистую когерентную картину 1-й гармоники считаем раз,
    а угловой разброс пучка берём ИЗ ЖИВОГО ПРОГОНА — блуждание пятна источника photon.x/y по окнам
    (σ пятна / D = угловой смаз), складываем в квадратуре к собственной расходимости и штатно размываем
    `screen_map`. Пишем идеальную и измеренную карты в `.npz` + пики/джиттер в результаты.
    УПРОЩЕНИЕ v1: картина на оси (как interference_compare), источник = photon[0]; геометрия экрана из
    `screens.csv` даёт только плечо D (params = screen_id)."""

    def __init__(self, params="", cadence="final", mesh=121, **kw):
        super().__init__(params, cadence)
        self.mesh = int(mesh)
        self._spot = []

    def sample(self, state, widx, plant):
        px = state.observables.get("photon.x"); py = state.observables.get("photon.y")
        if px is not None and py is not None and len(np.ravel(px)):
            self._spot.append((float(np.ravel(px)[0]), float(np.ravel(py)[0])))

    def _screen_dist(self, screen_id):
        path = os.path.join(_ROOT, "workspace", "tables", "screens.csv")
        if screen_id and os.path.exists(path):
            with open(path, encoding="utf-8-sig") as f:
                for r in csv.DictReader(f, delimiter=";"):
                    if (r.get("screen_id") or "").strip() == screen_id:
                        try:
                            return float(r.get("dist_m") or D_DEFAULT)
                        except ValueError:
                            break
        return D_DEFAULT

    def finalize(self, state, plant):
        import sys
        sys.path.insert(0, os.path.join(_ROOT, "undulators"))
        import undulator_sr as sr
        from undulator_field import trajectory, calibrate
        from screen_spectrum import beam_params, screen_map

        D = self._screen_dist(self.params)
        B0 = sr.B0_of_gap(GAP_MM); K = sr.K_of_B0(B0); g = sr.gamma()
        tr = trajectory(K, sr.P.lam_u, sr.P.N, g, 64); scale, _ = calibrate(tr)
        E0 = float(sr.harmonic_energy_keV(1, K, 0.0)); srp = float(sr.sigma_photon(E0)[1])
        beam0 = beam_params(sr.P, BETA_X, BETA_Y, sr.P.kappa, sr.P.sig_E)

        spot = np.asarray(self._spot, float) if self._spot else np.zeros((1, 2))
        jx = float(np.std(spot[:, 0])) / D; jy = float(np.std(spot[:, 1])) / D   # смаз угла из блуждания пятна
        window = max(8.0 * srp, 4.0 * math.hypot(jx, jy) + 1e-12)
        axis = np.array([0.0, 0.0, 1.0]); uhat = np.array([1.0, 0.0, 0.0]); vhat = np.array([0.0, 1.0, 0.0])
        P = np.zeros(3)

        def pattern(ax_, ay_):
            b = dict(beam0); b["sxp"] = math.hypot(beam0["sxp"], ax_); b["syp"] = math.hypot(beam0["syp"], ay_)
            mp = screen_map(tr, scale, P, (axis, uhat, vhat), axis, D, np.array([E0]), self.mesh, window, b)
            I = mp["I"][..., 0] if np.asarray(mp["I"]).ndim == 3 else mp["I"]
            return np.asarray(I, float), float(mp["peak"]), np.asarray(mp["extent"], float)

        Ii, pk_i, ext = pattern(0.0, 0.0)
        Im, pk_m, _ = pattern(jx, jy)
        self._art = ("screen", {"ideal": Ii, "measured": Im, "extent": ext,
                                "jitter_urad": np.array([jx * 1e6, jy * 1e6]), "D_m": np.array([D]),
                                "energy_keV": np.array([E0]), "window_mm": np.array([window * D * 1e3])},
                     f"screen={self.params} E={E0:.2f}keV")
        self.rows += [("screen.peak_ideal", 0, pk_i), ("screen.peak_meas", 0, pk_m),
                      ("screen.jitter_x_urad", 0, jx * 1e6), ("screen.jitter_y_urad", 0, jy * 1e6),
                      ("screen.sharpness_pct", 0, (pk_m / pk_i * 100.0) if pk_i else 0.0)]


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
