# -*- coding: utf-8 -*-
"""ПУНКТ 2 (геодезия): требование к разрешению съёмки R. Геодезия двигает магниты по ВНЕШНЕЙ метрологии
(не по пучку) → сводит монтажный разброс к остатку σ=R на всех юстируемых типах. Затем beam-based
настройка (орбита+LOCO, tuner.tune) → поднялась ли машина. Ансамбль: доля ok (vs бутстрап-стена), ДА+x,
остаточная орбита, запас корректоров. R=20 мкм. Запуск: py -3.14 survey_tolerance.py
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import ring_cpp, facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra
import errors as er, tuner as tn

R = 20e-6; NSEED = 6; APER = 0.012
KINDS = ("quad", "dipole", "sext", "oct", "dq")                  # геодезия юстирует все магниты до R
inf = dict(fr.info_of(26)); inf["_id"] = 26
ideal = fr._ring(inf, nsl=8)
resp = oc.cached_response_matrix(ideal, 26)
loco = tn.load_loco(26, 8); octx = tn.octupole_context(ideal)
bpm = oi.monitor_rows(ideal)
corr = np.array(resp["corr"])
SCR = os.path.join("workspace", "results", "_ws")


def da_x(ring, co, nturns=384, amax=0.012):
    base = np.array([co[0], co[1], co[2], co[3], 0.0, 0.0])
    def alive(x):
        st = base.copy(); st[0] += x
        tb = ring.track_tbt([st.tolist()], nturns, amax)[0]
        return np.all(np.isfinite(tb[-1]))
    lo, hi = 0.0, amax
    if not alive(0.2e-3):
        return 0.0
    for _ in range(14):                                         # бинпоиск границы
        mid = 0.5 * (lo + hi)
        if alive(mid):
            lo = mid
        else:
            hi = mid
    return lo


def resid_orbit(ring, co):
    o = ring.track_elements([co])[0]
    return float(np.sqrt(np.mean(o[bpm, 0] ** 2 + o[bpm, 2] ** 2)))


# идеал (R=0) для опоры
ide = da_x(ideal, ra.closed_orbit(ideal, dims=4))
print(f"R={R*1e6:.0f} мкм на {KINDS}; идеальная ДА+x {ide*1e3:.2f} мм, ансамбль {NSEED}\n", flush=True)
print(f"{'seed':>4} | {'ok':>3} | {'ДА+x, мм':>8} | {'ост.орбита, мкм':>15} | {'корр СКО, мкрад':>15} | стадии", flush=True)
res_ok = 0; das = []; orbs = []; cks = []
for s in range(NSEED):
    tol = {k: {"dx": R, "dy": R} for k in KINDS}
    wsd = os.path.join(SCR, f"ws_surv_{s}")
    er.perturbed_workspace(26, 1000 + s, tol, wsd)
    ring = ring_cpp.CppRing(lattice=26, ws=wsd, design_mev=inf["energy_mev"], nsl=8)
    out = tn.tune(ring, resp, APER, loco=loco, oct_ctx=octx, octupole_g=0.0)
    if not out["ok"]:
        print(f"{s:4d} | {'НЕТ':>3} | {'—':>8} | {'бутстрап-стена':>15} | {'':>15} | {'+'.join(out['stages'])}", flush=True)
        ring.close(); continue
    co = out["co"]; da = da_x(ring, co); orb = resid_orbit(ring, co)
    ck = float(np.sqrt(np.mean(np.array([v for _, v in out["bk"]]) ** 2 + np.array([v for _, v in out["bv"]]) ** 2)))
    res_ok += 1; das.append(da); orbs.append(orb); cks.append(ck)
    print(f"{s:4d} | {'да':>3} | {da*1e3:8.2f} | {orb*1e6:15.1f} | {ck*1e6:15.1f} | {'+'.join(out['stages'])}", flush=True)
    ring.close()
ideal.close()
print(f"\nИТОГ R=20 мкм: поднялось {res_ok}/{NSEED}; ДА+x медиана {np.median(das)*1e3 if das else float('nan'):.2f} мм "
      f"(идеал {ide*1e3:.2f}); ост.орбита мед {np.median(orbs)*1e6 if orbs else float('nan'):.1f} мкм; "
      f"корр СКО мед {np.median(cks)*1e6 if cks else float('nan'):.1f} мкрад", flush=True)
print("[готово] требование к геодезии: доля поднявшихся + ДА при остатке R", flush=True)
