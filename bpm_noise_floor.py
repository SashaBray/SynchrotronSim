# -*- coding: utf-8 -*-
"""ПУНКТ 2 статики медленной системы: ПОЛ СТАБИЛЬНОСТИ ОТ ШУМА BPM.

Коррекция орбиты зануляет показания датчиков; шум чтения σ_bpm заставляет корректоры гасить ШУМ,
впрыскивая реальное искажение орбиты: δθ = R⁺·шум. Остаточное движение в точках источника (центры
ИД-промежутков) от этого — минимально достижимая «тишина» наведения, ограниченная лишь разрешением BPM.
Независимо от юстировочных ошибок (линейно) → считаем на ИДЕАЛЬНОЙ решётке = чистый пол.

Модель: реальные мониторы (kind=monitor) с шумом → δθ=R_mon⁺·шум на корректорах → замкнутая орбита →
положение/угол в центрах ИД. Monte-Carlo по реализациям, свип σ_bpm. Запуск: py -3.14 bpm_noise_floor.py
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

NREAL, RCOND = 200, 1e-2
inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
resp = oc.cached_response_matrix(r, 26)
corr = resp["corr"]
mon = oi.monitor_rows(r)                                   # реальные датчики (kind=monitor)
Rmx = resp["Rx"][mon - 1]; Rmy = resp["Ry"][mon - 1]       # отклик монитор×корректор (гор/вер)
Pmx = np.linalg.pinv(Rmx, rcond=RCOND); Pmy = np.linalg.pinv(Rmy, rcond=RCOND)
# ИД-промежутки (центры длинных дрейфов) — как в combo_pointing
s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
print(f"ПОЛ ШУМА BPM (реш.26): мониторов {len(mon)}, корректоров {len(corr)}, ИД-промежутков {len(SIDX)}, реализаций {NREAL}", flush=True)
rng = np.random.default_rng(12345)
# БАЗОВАЯ орбита идеальной решётки (нулевой шум) — вычитаем, чтобы выделить ЧИСТЫЙ вклад шума
r.apply_perturb(); _co0 = ra.closed_orbit(r, dims=4); _o0 = r.track_elements([_co0])[0]
BX = _o0[SIDX, 0] + _o0[SIDX, 1] * DVEC; BXP = _o0[SIDX, 1]
BY = _o0[SIDX, 2] + _o0[SIDX, 3] * DVEC; BYP = _o0[SIDX, 3]
print(f"базовая орбита идеала в ИД: σx={np.sqrt(np.mean(BX**2))*1e6:.3f} σy={np.sqrt(np.mean(BY**2))*1e6:.3f} мкм (вычитается)", flush=True)


def floor(sigma):
    X = []; XP = []; Y = []; YP = []
    for _ in range(NREAL):
        nx = rng.normal(0, sigma, len(mon)); ny = rng.normal(0, sigma, len(mon))
        dthx = Pmx @ nx; dthy = Pmy @ ny                   # паразитные кики от шума
        r.apply_perturb(kicks=list(zip(corr, dthx.tolist())), vkicks=list(zip(corr, dthy.tolist())))
        co = ra.closed_orbit(r, dims=4)
        if not np.all(np.isfinite(co)):
            continue
        o = r.track_elements([co])[0]
        X.append(o[SIDX, 0] + o[SIDX, 1] * DVEC - BX); XP.append(o[SIDX, 1] - BXP)   # ОТКЛОНЕНИЕ от базы
        Y.append(o[SIDX, 2] + o[SIDX, 3] * DVEC - BY); YP.append(o[SIDX, 3] - BYP)
    rms = lambda a: float(np.sqrt(np.mean(np.concatenate(a) ** 2)))
    return rms(X) * 1e6, rms(Y) * 1e6, rms(XP) * 1e6, rms(YP) * 1e6   # мкм, мкрад


print(f"\n{'σ_bpm':>8} | {'σx':>7} {'σy':>7} мкм | {'σx′':>7} {'σy′':>7} мкрад  (в центрах ИД)", flush=True)
for sig in (0.1e-6, 0.3e-6, 1.0e-6, 3.0e-6):
    sx, sy, sxp, syp = floor(sig)
    print(f"{sig*1e6:6.1f}мкм | {sx:7.3f} {sy:7.3f}     | {sxp:7.3f} {syp:7.3f}", flush=True)
r.apply_perturb(); r.close()
print("\n[готово] (линейно по σ_bpm; усреднение по N чтениям в ОС делит на √N)", flush=True)
