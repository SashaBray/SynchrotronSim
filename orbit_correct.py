# -*- coding: utf-8 -*-
"""Коррекция замкнутой орбиты корректорами (глобальная, по матрице отклика + SVD).

Погрешности установки (рассогласование магнитов, ошибки поля — см. errors.py) искажают замкнутую орбиту
всего кольца. Корректоры — это магниты-кикеры (в решётках СИЛА роль играют слоты SH1A/SH2B/SH3E, заведённые
как kind=kicker), которые задают локальные углы (гор. и верт.), чтобы выправить орбиту глобально.

Схема (стандартная для накопителей):
  1. Матрица отклика R (BPM × корректоры) — на ИДЕАЛЬНОЙ решётке, ОДИН раз: единичный кик каждого
     корректора → смещение орбиты во всех точках наблюдения. Отдельно для x (KICK) и y (VKICK).
  2. На каждую машину с ошибками: измерить искажённую орбиту o в точках наблюдения, решить θ = −R⁺·o
     (псевдообратная через SVD с обрезкой малых сингулярных чисел — гасит неэффективные комбинации
     корректоров, иначе они дают гигантские углы), применить кики, повторить 2–4 раза (отклик линейный
     лишь приближённо на искажённой орбите). Настройка СТАЦИОНАРНА для данной машины.

Точки наблюдения (BPM) — выходы всех магнитов (track_elements): least-squares минимизирует орбиту по
всему кольцу, а не в отдельных мониторах. Один сервер ядра переиспользуется через CppRing.apply_perturb
(RESET сохраняет ошибки в ws, сбрасывает лишь кики), так что ни матрица отклика, ни итерации не плодят
процессы.
"""
import os
import pickle

import numpy as np

import ring_analysis as ra

RESP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workspace", "results")


def corrector_indices(ring):
    """Индексы корректоров (kind=kicker) среди магнитов кольца."""
    return [i for i, k in enumerate(ring.kinds) if k == "kicker"]


def _orbit(ring, bpm):
    """Замкнутая орбита в точках наблюдения: (ox, oy) или None, если машина неустойчива."""
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return None, None
    orb = ring.track_elements([co])[0]
    return co, orb


def response_matrix(ideal, corr=None, bpm=None, theta=1e-5):
    """Матрицы отклика на ИДЕАЛЬНОЙ решётке: Rx, Ry [м/рад] размера (BPM × корректоры).
    ideal мутируется (apply_perturb: RESET + единичный кик), в конце сбрасывается. bpm — индексы строк
    в выводе track_elements (по умолчанию выходы всех магнитов, строки 1..n_mag)."""
    corr = corrector_indices(ideal) if corr is None else corr
    bpm = np.arange(1, ideal.n_mag + 1) if bpm is None else np.asarray(bpm)
    ideal.apply_perturb()                                       # чистая решётка, без киков
    co0, base = _orbit(ideal, bpm)
    if co0 is None:
        raise RuntimeError("идеальная решётка неустойчива — матрица отклика не считается")
    bx0, by0 = base[bpm, 0], base[bpm, 2]
    Rx = np.zeros((len(bpm), len(corr)))
    Ry = np.zeros((len(bpm), len(corr)))
    for j, c in enumerate(corr):
        ideal.apply_perturb(kicks=[(c, theta)])
        co, orb = _orbit(ideal, bpm)
        if co is not None:
            Rx[:, j] = (orb[bpm, 0] - bx0) / theta
        ideal.apply_perturb(vkicks=[(c, theta)])
        co, orb = _orbit(ideal, bpm)
        if co is not None:
            Ry[:, j] = (orb[bpm, 2] - by0) / theta
    ideal.apply_perturb()                                       # вернуть в исходное состояние
    return dict(corr=list(corr), bpm=bpm, Rx=Rx, Ry=Ry)


def cached_response_matrix(ideal, lattice, theta=1e-5, rebuild=False):
    """Матрица отклика с кэшем на диске (workspace/results/orbit_resp_exp{lattice}.pkl). Инвалидируется,
    если сменились геометрия (n_mag/периметр) или число корректоров — считается заново (≈2-3 мин)."""
    path = os.path.join(RESP_DIR, f"orbit_resp_exp{int(lattice)}.pkl")
    corr = corrector_indices(ideal)
    if not rebuild and os.path.exists(path):
        try:
            with open(path, "rb") as f:
                r = pickle.load(f)
            if (r.get("n_mag") == ideal.n_mag and abs(r.get("circ", -1) - ideal.circumference) < 1e-6
                    and r.get("corr") == list(corr)):
                return r
        except Exception:
            pass
    r = response_matrix(ideal, corr=corr, theta=theta)
    r["n_mag"] = ideal.n_mag; r["circ"] = ideal.circumference
    os.makedirs(RESP_DIR, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(r, f)
    return r


def correct(ring, resp, iters=4, rcond=1e-2):
    """Скорректировать орбиту кольца ring (с ошибками в ws) корректорами по матрице отклика resp.
    ring мутируется: в конце на нём применены итоговые кики (готов для трекинга ДА/диффузии). Возврат
    dict(co, thx, thy, orb_before, orb_after) или None, если машина неустойчива даже в ходе коррекции.
      thx, thy — углы корректоров [рад]; orb_before/after — (ox, oy) в точках наблюдения [м]."""
    corr, bpm = resp["corr"], resp["bpm"]
    Px = np.linalg.pinv(resp["Rx"], rcond=rcond)
    Py = np.linalg.pinv(resp["Ry"], rcond=rcond)
    thx = np.zeros(len(corr))
    thy = np.zeros(len(corr))
    orb_before = None
    for _ in range(iters):
        ring.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
        co, orb = _orbit(ring, bpm)
        if co is None:
            return None
        ox, oy = orb[bpm, 0], orb[bpm, 2]
        if orb_before is None:                                  # первый проход: кики нулевые → орбита ДО
            orb_before = (ox.copy(), oy.copy())
        thx = thx - Px @ ox
        thy = thy - Py @ oy
    ring.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
    co, orb = _orbit(ring, bpm)
    if co is None:
        return None
    return dict(co=co, thx=thx, thy=thy, orb_before=orb_before, orb_after=(orb[bpm, 0], orb[bpm, 2]))


def rms(a):
    a = np.asarray(a, float)
    return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")
