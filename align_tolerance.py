# -*- coding: utf-8 -*-
"""ПУНКТ 2: ПОЭЛЕМЕНТНАЯ юстировка по слепому алгоритму — границы допусков установки.
Машина со случайными сдвигами квадруполей σ → слепой инференс (orbit_infer.infer_shifts по орбите) →
ДВИГАЕМ магниты на предсказанное (итеративно) → остаток. Развёртка по σ: порог, с которого статнастройка
поднимает машину (замкнутая орбита существует И инференс сходится к малому остатку).
Теоретическая граница: идеальные BPM (малый шум, без смещений датчиков). Запуск: py -3.14 align_tolerance.py
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_infer as oi, ring_analysis as ra

inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
bpm = oi.monitor_rows(r)                                          # реальные датчики (kind=monitor)
print("строю базис сдвигов (кэш) ...", flush=True)
basis = oi.shift_basis(r, bpm=bpm, lattice=26, nsl=8)
elems = basis["elems"]; k1l = basis["k1l"]; ne = len(elems)
print(f"элементов (квадруполей) {ne}, датчиков {len(bpm)}", flush=True)

BPM_NOISE = 0.1e-6                                                # шум датчика (теор. граница: без смещений BPM)
NITER = 6


def apply_shifts(dx, dy):
    kk = [(int(elems[i]), float(k1l[i] * dx[i])) for i in range(ne)]
    vv = [(int(elems[i]), float(-k1l[i] * dy[i])) for i in range(ne)]
    r.apply_perturb(kicks=kk, vkicks=vv)


def orbit_rms():
    co = ra.closed_orbit(r, dims=4)
    if not np.all(np.isfinite(co)):
        return None
    o = r.track_elements([co])[0]
    return float(np.sqrt(np.mean(o[bpm, 0] ** 2 + o[bpm, 2] ** 2)))


def realign(dx0, dy0, seed):
    rng = np.random.default_rng(seed)
    dx, dy = dx0.copy(), dy0.copy()                              # текущий ОСТАТОЧНЫЙ рассинхрон
    for _ in range(NITER):
        apply_shifts(dx, dy)
        meas = oi.measure_orbit(r, bpm, BPM_NOISE, None, rng)
        if meas is None:
            return None                                         # орбита не замкнулась → не поднять по орбите
        ex, ey = oi.infer_shifts(basis, meas[0], meas[1])
        dx = dx - ex; dy = dy - ey                              # ДВИГАЕМ магниты на предсказанное
    apply_shifts(dx, dy)
    return dx, dy, orbit_rms()


SIGS = [5e-6, 10e-6, 20e-6, 50e-6, 100e-6, 200e-6]; NSEED = 6
print(f"\n{'σ_нач, мкм':>9} | {'орбита закрылась':>16} | {'остаток σ, мкм':>14} | {'×сниж':>6} | {'ост. орбита, мкм':>16}")
for sig in SIGS:
    res_sig = []; res_orb = []; closed = 0
    for s in range(NSEED):
        rng = np.random.default_rng(100 + s)
        dx0 = rng.normal(0, sig, ne); dy0 = rng.normal(0, sig, ne)
        out = realign(dx0, dy0, seed=s)
        if out is None:
            continue
        dx, dy, orb = out
        closed += 1
        res_sig.append(np.sqrt(np.mean(dx ** 2 + dy ** 2)))
        if orb is not None:
            res_orb.append(orb)
    if res_sig:
        rs = float(np.mean(res_sig)) * 1e6; ro = (float(np.mean(res_orb)) * 1e6 if res_orb else float("nan"))
        red = (sig * 1e6) / max(rs, 1e-9)
        print(f"{sig*1e6:9.0f} | {f'{closed}/{NSEED}':>16} | {rs:14.2f} | ×{red:4.0f} | {ro:16.2f}", flush=True)
    else:
        print(f"{sig*1e6:9.0f} | {f'0/{NSEED}':>16} | {'—':>14} | {'—':>6} | {'не закрылась':>16}", flush=True)
r.apply_perturb(); r.close()
print("\n[готово] порог допуска = макс σ, где орбита закрывается и остаток мал (инференс сошёлся)", flush=True)
