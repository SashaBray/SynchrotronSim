# -*- coding: utf-8 -*-
"""Единый построитель матриц отклика — один механизм вместо россыпи
(orbit_correct.response_matrix / threading_response, regulator.build_response, xbpm_svd.angle_response,
orbit_infer.shift_basis, build_loco_jac).

Матрица отклика = (Δнаблюдаемая) / (единичное возмущение актуатора), снятая с кольца; строки —
наблюдаемые, столбцы — актуаторы. Два расширяемых интерфейса:
  applier(ring, actuator_id, eps) — как внести единичное возмущение в кольцо (через apply_perturb,
                                    с RESET: каждый актуатор зондируется независимо);
  reader(observable_id, orbit)    — как считать вектор наблюдаемой из снимка track_elements.
Добавить новый тип актуатора/наблюдаемой = одна функция, а не новый построитель. Кэш на диске.
"""
import hashlib
import os
import pickle

import numpy as np

import ring_analysis as ra

EPS_DEFAULT = 1e-6
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workspace", "results")


def _orbit(ring, dims=4):
    """Снимок орбиты: замкнутая орбита (dims=4 поперечная δ=0; dims=6 — 6D равновесие с ВЧ/СИ,
    захватывает дисперсию и энергию) → track_elements, [n_mag+1 × 6]."""
    co = ra.closed_orbit(ring, dims=dims)
    return ring.track_elements([co])[0]


def build_response(ring, actuators, observables, applier, reader, eps=EPS_DEFAULT, dims=4):
    """Снять {observable_id: R[rows × actuators]} с кольца единичными возмущениями актуаторов.
    dims=6 снимает отклик на 6D-равновесии (дисперсия/энергия в горизонтали) — модель «как на машине»."""
    ring.apply_perturb()                                     # чистое кольцо (RESET киков)
    base = _orbit(ring, dims)
    base_obs = {o: np.asarray(reader(o, base), float) for o in observables}
    out = {o: np.zeros((len(base_obs[o]), len(actuators))) for o in observables}
    for j, a in enumerate(actuators):
        applier(ring, a, eps)                                # apply_perturb с RESET внутри
        orb = _orbit(ring, dims)
        for o in observables:
            out[o][:, j] = (np.asarray(reader(o, orb), float) - base_obs[o]) / eps
    ring.apply_perturb()                                     # вернуть кольцо в исходное
    return out


def ring_fingerprint(ring):
    """КОНТЕНТНЫЙ отпечаток кольца: конфиг (решётка/срезы/СИ/ВЧ/энергия/ws) + геометрия и силы ВСЕХ кусков
    (piece_meta: Leff, θ0, k1, k2, k3). Берётся прямо с C++-сервера, детерминирован. Меняется при любой
    правке таблиц решётки, поэтому входит в ключ кэша — правка решётки инвалидирует resp_*.pkl САМА
    (закрывает готчу: раньше ключ зависел от НОМЕРА решётки, а не содержимого). См. [[matrix-core-refactor-single-source]]."""
    cfg = (getattr(ring, "lattice", None), getattr(ring, "nsl", None), getattr(ring, "order", None),
           bool(getattr(ring, "sr", False)), bool(getattr(ring, "isr", False)),
           round(float(getattr(ring, "design_mev", 0.0)), 6), getattr(ring, "ws", None),
           repr(getattr(ring, "rf", None)))
    geom = (round(float(getattr(ring, "circumference", 0.0)), 9),
            int(getattr(ring, "n_mag", 0)), int(getattr(ring, "n_pieces", 0)))
    try:
        meta = ring.piece_meta()                             # (магнит, стык, изгиб, Leff, θ0, k1, k2, k3) на кусок
    except Exception:
        meta = None                                          # не-CppRing: отпечаток по конфигу+геометрии
    return hashlib.md5(repr((cfg, geom, meta)).encode()).hexdigest()[:16]


def cached_build(key, ring, actuators, observables, applier, reader, eps=EPS_DEFAULT, rebuild=False, dims=4):
    """build_response с кэшем на диске по отпечатку (контент кольца, key, dims, актуаторы, наблюдаемые, eps).
    Контентный отпечаток кольца в теге → правка таблиц решётки инвалидирует кэш автоматически."""
    tag = hashlib.md5(repr((ring_fingerprint(ring), key, dims, list(actuators), list(observables), eps))
                      .encode()).hexdigest()[:12]
    path = os.path.join(CACHE_DIR, f"resp_{tag}.pkl")
    if not rebuild and os.path.exists(path):
        try:
            with open(path, "rb") as f:
                return pickle.load(f)
        except Exception:
            pass
    out = build_response(ring, actuators, observables, applier, reader, eps, dims=dims)
    os.makedirs(CACHE_DIR, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(out, f)
    return out
