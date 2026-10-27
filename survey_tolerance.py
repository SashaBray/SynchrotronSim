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

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
RS = [20e-6] if SMOKE else [5e-6, 10e-6, 15e-6, 20e-6]
NSEED = 6; APER = 0.012; GOCT = tn.OCTUPOLE_PARETO_G             # g=+0.25 (бесплатная ДА, Парето)
KINDS = ("quad", "dipole", "sext", "oct", "dq")                  # геодезия юстирует все магниты до R
inf = dict(fr.info_of(26)); inf["_id"] = 26
ideal = fr._ring(inf, nsl=8)
resp = oc.cached_response_matrix(ideal, 26)
loco = tn.load_loco(26, 8); octx = tn.octupole_context(ideal)
bpm = oi.monitor_rows(ideal)
corr = np.array(resp["corr"])
Px = np.linalg.pinv(resp["Rx"][bpm - 1], rcond=1e-2); Py = np.linalg.pinv(resp["Ry"][bpm - 1], rcond=1e-2)
SCR = os.path.join("workspace", "results", "_ws")


def orbit_recorrect(ring, out):
    """Финальный орбитный шаг ПОСЛЕ LOCO: LOCO-подстройки градиентов заново искажают орбиту (feed-down),
    добавляем корректорные кики, чтобы занулить остаточную орбиту в BPM (точность)."""
    o = ring.track_elements([out["co"]])[0]
    dthx = -(Px @ o[bpm, 0]); dthy = -(Py @ o[bpm, 2])
    kk = dict(out["bk"]); vv = dict(out["bv"])
    for j, c in enumerate(corr):
        kk[c] = kk.get(c, 0.0) + float(dthx[j]); vv[c] = vv.get(c, 0.0) + float(dthy[j])
    ring.apply_perturb(kicks=list(kk.items()), vkicks=list(vv.items()), dks=out["bd"], dk3s=out["dk3"])
    co = ra.closed_orbit(ring, dims=4)
    return co if np.all(np.isfinite(co)) else None


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


ide = da_x(ideal, ra.closed_orbit(ideal, dims=4))
print(f"идеальная ДА+x {ide*1e3:.2f} мм; октуполь g={GOCT:+.2f}; юстируемые {KINDS}", flush=True)


def one(sigma, seed):
    tol = {k: {"dx": sigma, "dy": sigma} for k in KINDS}
    wsd = os.path.join(SCR, f"ws_surv_{int(sigma*1e9)}_{seed}")
    er.perturbed_workspace(26, 1000 + seed, tol, wsd)
    ring = ring_cpp.CppRing(lattice=26, ws=wsd, design_mev=inf["energy_mev"], nsl=8)
    out = tn.tune(ring, resp, APER, loco=loco, oct_ctx=octx, octupole_g=GOCT)
    if not out["ok"]:
        ring.close(); return None
    co2 = orbit_recorrect(ring, out)                            # финальный орбитный шаг после LOCO
    co = co2 if co2 is not None else out["co"]
    r = (da_x(ring, co), resid_orbit(ring, co), resid_orbit(ring, out["co"]))
    ring.close(); return r


if SMOKE:                                                        # сравнить с/без финальной перекоррекции (R=20, seed0)
    tol = {k: {"dx": 20e-6, "dy": 20e-6} for k in KINDS}
    wsd = os.path.join(SCR, "ws_surv_sm")
    er.perturbed_workspace(26, 1000, tol, wsd)
    ring = ring_cpp.CppRing(lattice=26, ws=wsd, design_mev=inf["energy_mev"], nsl=8)
    out = tn.tune(ring, resp, APER, loco=loco, oct_ctx=octx, octupole_g=GOCT)
    orb0 = resid_orbit(ring, out["co"]); da0 = da_x(ring, out["co"])
    co2 = orbit_recorrect(ring, out)
    orb1 = resid_orbit(ring, co2); da1 = da_x(ring, co2)
    print(f"[smoke R=20 seed0 g={GOCT:+.2f}] БЕЗ финала: орбита {orb0*1e6:.1f} мкм, ДА {da0*1e3:.2f} мм; "
          f"С финалом: орбита {orb1*1e6:.1f} мкм, ДА {da1*1e3:.2f} мм", flush=True)
    ring.close(); ideal.close(); sys.exit(0)

print(f"\n{'R, мкм':>6} | {'поднялось':>9} | {'ДА+x мед, мм':>12} | {'ост.орбита мед, мкм':>19}", flush=True)
for sig in RS:
    rr = [one(sig, s) for s in range(NSEED)]
    ok = [x for x in rr if x is not None]
    if not ok:
        print(f"{sig*1e6:6.0f} | {f'0/{NSEED}':>9} | {'—':>12} | {'—':>19}", flush=True); continue
    da = np.median([x[0] for x in ok]) * 1e3; orb = np.median([x[1] for x in ok]) * 1e6
    print(f"{sig*1e6:6.0f} | {f'{len(ok)}/{NSEED}':>9} | {da:12.2f} | {orb:19.1f}", flush=True)
ideal.close()
print(f"\n[готово] требование к геодезии: ДА+x (идеал {ide*1e3:.2f} мм) и орбита как функция R", flush=True)
