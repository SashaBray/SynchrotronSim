# -*- coding: utf-8 -*-
"""Слепой инференс рассогласования по орбите: предсказать сдвиги магнитов из показаний BPM.

Задача (валидация метода наладки реальной установки): задаём случайный шум (errors.py), меряем орбиту
в BPM, и ВСЛЕПУЮ (не заглядывая в истинные ошибки) восстанавливаем, какие сдвиги были применены. Истина
(errors_truth.csv из perturbed_workspace) нужна лишь чтобы ОЦЕНИТЬ качество предсказания.

Физика: искажённая замкнутая орбита рождается дипольными киками ошибок. Сдвиг квадруполя на dx даёт
дипольный кик θx = +K1·L·dx (feed-down), сдвиг на dy — θy = −K1·L·dy. Значит
    o = B · s,   B[:,j] = ±K1ⱼ·Lⱼ · (отклик орбиты на единичный кик в квадруполе j),
и сдвиги восстанавливаются регуляризованным инверсом s = B⁺·o (SVD с обрезкой).

ЧТО ОРБИТА ВИДИТ, А ЧТО НЕТ (честно):
  • сдвиги квадруполей/DQ, крен и ошибка поля ДИПОЛЕЙ — кикают орбиту → восстановимы;
  • ошибка ГРАДИЕНТА квадруполей (поле/ток) — НЕ кикает орбиту (даёт бета-биение) → только из матрицы
    отклика (LOCO), не отсюда;
  • сдвиги секступолей/октуполей — дают градиент, не кик → орбите не видны.
Версия 1: только пуре-квадруполи (несекторные, кикаются напрямую). DQ/диполи (секторные) — отдельным
базисом через сдвиг в ws (следующий шаг).
"""
import os
import pickle

import numpy as np

import ring_analysis as ra

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workspace", "results")


def quad_indices(ring, kinds=("quad",)):
    """Индексы несекторных квадруполей с ненулевым K1 (их сдвиг кикуем через feed-down)."""
    out = []
    for i, k in enumerate(ring.kinds):
        if k in kinds and abs(float(getattr(ring.mags[i], "k1", 0.0) or 0.0)) > 1e-9:
            out.append(i)
    return out


def shift_basis(ideal, elems=None, bpm=None, theta=1e-5, lattice=None, nsl=None, rebuild=False):
    """Базис B (орбита в BPM от единичного СДВИГА каждого элемента), [n_bpm × n_elem], для x и y.
    Строится на ИДЕАЛЬНОЙ решётке через кик (apply_perturb на одном сервере) с масштабом K1·L. Кэш на диске.
    Возврат dict(elems, bpm, Bx, By, k1l)."""
    nsl = int(getattr(ideal, "nsl", 40)) if nsl is None else nsl
    path = os.path.join(CACHE, f"shift_basis_exp{int(lattice)}_nsl{nsl}.pkl") if lattice is not None else None
    elems = quad_indices(ideal) if elems is None else list(elems)
    bpm = np.arange(1, ideal.n_mag + 1) if bpm is None else np.asarray(bpm)
    if path and not rebuild and os.path.exists(path):
        try:
            with open(path, "rb") as f:
                b = pickle.load(f)
            if b.get("elems") == elems and len(b.get("bpm", [])) == len(bpm) and b.get("n_mag") == ideal.n_mag:
                return b
        except Exception:
            pass
    k1l = np.array([float(ideal.mags[j].k1) * float(ideal.mags[j].Leff) for j in elems])
    ideal.apply_perturb()
    co0 = ra.closed_orbit(ideal, dims=4)
    base = ideal.track_elements([co0])[0]
    bx0, by0 = base[bpm, 0], base[bpm, 2]
    Bx = np.zeros((len(bpm), len(elems)))
    By = np.zeros((len(bpm), len(elems)))
    for k, j in enumerate(elems):
        ideal.apply_perturb(kicks=[(j, theta)])
        co = ra.closed_orbit(ideal, dims=4); orb = ideal.track_elements([co])[0]
        Bx[:, k] = k1l[k] * (orb[bpm, 0] - bx0) / theta          # сдвиг dx → θx=+K1L·dx
        ideal.apply_perturb(vkicks=[(j, theta)])
        co = ra.closed_orbit(ideal, dims=4); orb = ideal.track_elements([co])[0]
        By[:, k] = -k1l[k] * (orb[bpm, 2] - by0) / theta         # сдвиг dy → θy=−K1L·dy
    ideal.apply_perturb()
    b = dict(elems=elems, bpm=bpm, Bx=Bx, By=By, k1l=k1l, n_mag=ideal.n_mag, nsl=nsl)
    if path:
        os.makedirs(CACHE, exist_ok=True)
        with open(path, "wb") as f:
            pickle.dump(b, f)
    return b


def bpm_rows_at_quads(ring):
    """Строки track_elements (выходы квадруполей) как реалистичный набор точек наблюдения ~322, по BPM
    на квадруполь (естественная пара к BBA). Выход магнита i — строка i+1."""
    return np.array([i + 1 for i in quad_indices(ring)])


def draw_bpm_offsets(n_bpm, seed, sigma, cutoff=2.5):
    """Статичные поперечные смещения датчиков (x,y) на машину: BPM меряет положение относительно СВОЕГО
    (смещённого) центра. Возврат (ox_off, oy_off) [м]. Это то, что на реальной машине разводит BBA."""
    rng = np.random.default_rng(seed)
    def col():
        a = rng.normal(0.0, sigma, n_bpm) if sigma else np.zeros(n_bpm)
        if sigma:
            m = np.abs(a) > cutoff * sigma
            while np.any(m):
                a[m] = rng.normal(0.0, sigma, int(m.sum())); m = np.abs(a) > cutoff * sigma
        return a
    return col(), col()


def measure_orbit(ring, bpm, bpm_noise=0.0, offset=None, rng=None):
    """Показания BPM возмущённой машины: (орбита в точке) − (смещение датчика) + гауссов шум датчика [м].
    offset=(ox_off,oy_off) — статичные смещения BPM на машину (см. draw_bpm_offsets). None, если орбиты нет
    (на больших ошибках — тогда инференс по орбите неприменим без протяжки/траектории)."""
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return None
    orb = ring.track_elements([co])[0]
    ox, oy = orb[bpm, 0].copy(), orb[bpm, 2].copy()
    if offset is not None:                                      # датчик меряет относительно своего центра
        ox = ox - np.asarray(offset[0]); oy = oy - np.asarray(offset[1])
    if bpm_noise:
        rng = rng or np.random.default_rng()
        ox += rng.normal(0.0, bpm_noise, ox.shape)
        oy += rng.normal(0.0, bpm_noise, oy.shape)
    return ox, oy


def infer_shifts(basis, ox, oy, rcond=1e-2):
    """Слепое восстановление сдвигов: s = B⁺·o (SVD-обрезка rcond). Возврат (dx_est, dy_est) по элементам."""
    dx = np.linalg.pinv(basis["Bx"], rcond=rcond) @ ox
    dy = np.linalg.pinv(basis["By"], rcond=rcond) @ oy
    return dx, dy


def load_truth(ws, elems):
    """Истинные сдвиги (dx,dy) по элементам базиса из errors_truth.csv (для ОЦЕНКИ, не для инференса)."""
    import csv
    path = os.path.join(ws, "tables", "errors_truth.csv")
    dxs = {}; dys = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter=";"):
            i = int(r["idx"]); dxs[i] = float(r["dx"]); dys[i] = float(r["dy"])
    return (np.array([dxs.get(j, 0.0) for j in elems]),
            np.array([dys.get(j, 0.0) for j in elems]))


def score(est, true):
    """Качество предсказания: корреляция, СКО ошибки, СКО истины, доля объяснённой дисперсии R²."""
    est = np.asarray(est); true = np.asarray(true)
    err = est - true
    ss_res = float(np.sum(err ** 2)); ss_tot = float(np.sum((true - true.mean()) ** 2))
    corr = float(np.corrcoef(est, true)[0, 1]) if est.std() > 0 and true.std() > 0 else float("nan")
    return dict(corr=corr, rms_err=float(np.sqrt(np.mean(err ** 2))),
                rms_true=float(np.sqrt(np.mean(true ** 2))),
                r2=(1.0 - ss_res / ss_tot) if ss_tot > 0 else float("nan"))
