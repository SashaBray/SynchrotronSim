# -*- coding: utf-8 -*-
"""ПУНКТ 3 (часть 2): УПРАВЛЯЕМОСТЬ ТОЧЕК ИСТОЧНИКА корректорами.

Могут ли 120 корректоров независимо управлять положением пучка в центрах ИД-промежутков?
Строим отклик «корректор → положение источника» S (38 ИД × 120 корр) трекингом единичных киков,
SVD → ранг (сколько источников независимо управляемы) и обусловленность. 38 источников < 120 корр ⇒
ждём полный ранг 38 без вырождения. Запуск: py -3.14 source_controllability.py
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, ring_analysis as ra

TH = 1e-5
inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8); resp = oc.cached_response_matrix(r, 26); corr = resp["corr"]
s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
nID = len(SIDX)
print(f"УПРАВЛЯЕМОСТЬ ИСТОЧНИКА: {nID} ИД-промежутков, {len(corr)} корректоров", flush=True)


def src(co):
    o = r.track_elements([co])[0]
    return o[SIDX, 0] + o[SIDX, 1] * DVEC, o[SIDX, 2] + o[SIDX, 3] * DVEC


r.apply_perturb(); co0 = ra.closed_orbit(r, dims=4); x0, y0 = src(co0)
Sx = np.zeros((nID, len(corr))); Sy = np.zeros((nID, len(corr)))
for j, c in enumerate(corr):
    r.apply_perturb(kicks=[(c, TH)]); co = ra.closed_orbit(r, dims=4)
    if np.all(np.isfinite(co)): Sx[:, j] = (src(co)[0] - x0) / TH
    r.apply_perturb(vkicks=[(c, TH)]); co = ra.closed_orbit(r, dims=4)
    if np.all(np.isfinite(co)): Sy[:, j] = (src(co)[1] - y0) / TH
r.apply_perturb(); r.close()

for name, S in (("ИСТОЧНИК X (гор-кики)", Sx), ("ИСТОЧНИК Y (вер-кики)", Sy)):
    s = np.linalg.svd(S, compute_uv=False)
    smax = s[0]; neff = int(np.sum(s > 1e-2 * smax))
    print(f"\n=== {name} === S {S.shape}")
    print(f"  σ_max={smax:.3e} σ_min={s[-1]:.3e} обусл.={smax/s[-1]:.1e}")
    print(f"  ранг (σ>1e-2·σmax)={neff}/{nID}  → {'ВСЕ источники независимо управляемы' if neff>=nID else 'есть слабоуправляемые'}")
    print(f"  слабейшие 3 σ/σmax: {[f'{v/smax:.2e}' for v in s[-3:]]}")
print("\n[готово]", flush=True)
