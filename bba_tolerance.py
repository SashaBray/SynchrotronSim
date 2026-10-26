# -*- coding: utf-8 -*-
"""ПУНКТ 2 (BBA): поэлементная юстировка через beam-based alignment — границы допусков.
Возмущённое кольцо (сдвиги квадруполей σ встроены в ws через errors.perturbed_workspace) → BBA меряет
положение пучка в каждом квадруполе xb → реконструкция сдвига dx_est = (орбита в центре квадруполя) − xb →
остаток dx_true − dx_est. Развёртка по σ: насколько BBA снижает рассинхрон и где потолок замыкания орбиты.
Запуск: py -3.14 bba_tolerance.py [smoke]
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import ring_cpp, facility_report as fr, orbit_infer as oi, ring_analysis as ra, errors as er

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
SCR = os.environ.get("CLAUDE_SCRATCH", os.path.join("workspace", "results", "_ws"))
inf = dict(fr.info_of(26)); inf["_id"] = 26
ideal = fr._ring(inf, nsl=8)
bpm = oi.monitor_rows(ideal)
print("базис сдвигов (кэш) ...", flush=True)
basis = oi.shift_basis(ideal, bpm=bpm, lattice=26, nsl=8)
elems = np.array(basis["elems"]); k1l = basis["k1l"]; ne = len(elems)
print(f"квадруполей {ne}, датчиков {len(bpm)}", flush=True)


def cx_cy(ring):
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return None
    o = ring.track_elements([co])[0]
    cx = 0.5 * (o[elems, 0] + o[elems + 1, 0]); cy = 0.5 * (o[elems, 2] + o[elems + 1, 2])
    return cx, cy, float(np.sqrt(np.mean(o[bpm, 0] ** 2 + o[bpm, 2] ** 2)))


def one_machine(sigma, seed):
    tol = {"quad": {"dx": sigma, "dy": sigma}}                   # только квадруполи, σ по dx/dy
    wsd = os.path.join(SCR, f"ws_a{int(sigma*1e9)}_{seed}")
    er.perturbed_workspace(26, seed, tol, wsd)
    ring = ring_cpp.CppRing(lattice=26, ws=wsd, design_mev=inf["energy_mev"], nsl=8)
    dxt, dyt = oi.load_truth(wsd, list(elems))
    cc = cx_cy(ring)
    if cc is None:
        ring.close(); return None
    cx, cy, orb0 = cc
    xb, yb = oi.bba(ring, basis)                                # beam-in-quad (322 решения)
    ring.close()
    ok = np.isfinite(xb) & np.isfinite(yb)
    dxe = cx - xb; dye = cy - yb                                # реконструкция сдвига
    rx = dxt - dxe; ry = dyt - dye                              # остаток
    init = np.sqrt(np.mean(dxt[ok] ** 2 + dyt[ok] ** 2))
    resid = np.sqrt(np.mean(rx[ok] ** 2 + ry[ok] ** 2))
    cor = float(np.corrcoef(np.r_[dxe[ok], dye[ok]], np.r_[dxt[ok], dyt[ok]])[0, 1])
    return init, resid, cor, orb0, int(ok.sum())


if SMOKE:
    for sig in (0.5e-6, 1e-6, 2e-6):
        tol = {"quad": {"dx": sig, "dy": sig}}
        wsd = os.path.join(SCR, f"ws_sm_{int(sig*1e9)}")
        er.perturbed_workspace(26, 0, tol, wsd)
        ring = ring_cpp.CppRing(lattice=26, ws=wsd, design_mev=inf["energy_mev"], nsl=8)
        dxt, dyt = oi.load_truth(wsd, list(elems))
        cc = cx_cy(ring)
        if cc is None:
            print(f"σ={sig*1e6:.1f}: орбита не замкнулась", flush=True); ring.close(); continue
        cx, cy, orb0 = cc
        btx, bty = oi.beam_in_quads_true(ring, basis)           # истинный beam-in-quad
        xb, yb = oi.bba(ring, basis); ring.close()
        ok = np.isfinite(xb) & np.isfinite(yb)
        cbba = float(np.corrcoef(np.r_[xb[ok], yb[ok]], np.r_[btx[ok], bty[ok]])[0, 1])  # BBA vs истина
        dxe = cx - xb; dye = cy - yb
        cest = float(np.corrcoef(np.r_[dxe[ok], dye[ok]], np.r_[dxt[ok], dyt[ok]])[0, 1])
        init = np.sqrt(np.mean(dxt[ok] ** 2 + dyt[ok] ** 2)); resid = np.sqrt(np.mean((dxt - dxe)[ok] ** 2 + (dyt - dye)[ok] ** 2))
        print(f"σ={sig*1e6:.1f} мкм: орбита0 {orb0*1e6:6.1f} мкм; BBA↔истина corr {cbba:.3f}; "
              f"сдвиг corr {cest:.3f}; init {init*1e6:.2f}→ост {resid*1e6:.2f} (×{init/max(resid,1e-12):.1f})", flush=True)
    ideal.close(); sys.exit(0)

SIGS = [5e-6, 10e-6, 20e-6, 30e-6]; NSEED = 3
print(f"\n{'σ, мкм':>7} | {'закрыл':>6} | {'init σ':>7} | {'остаток σ':>9} | {'×сниж':>6} | {'корр BBA':>8} | {'орбита0':>8}")
for sig in SIGS:
    res = [one_machine(sig, s) for s in range(NSEED)]
    ok = [x for x in res if x is not None]
    if not ok:
        print(f"{sig*1e6:7.0f} | {f'0/{NSEED}':>6} | {'—':>7} | {'не закрылась':>9}", flush=True); continue
    init = np.mean([x[0] for x in ok]) * 1e6; resid = np.mean([x[1] for x in ok]) * 1e6
    cor = np.mean([x[2] for x in ok]); orb = np.mean([x[3] for x in ok]) * 1e6
    print(f"{sig*1e6:7.0f} | {f'{len(ok)}/{NSEED}':>6} | {init:7.2f} | {resid:9.2f} | "
          f"×{init/max(resid,1e-9):5.1f} | {cor:8.3f} | {orb:8.1f}", flush=True)
ideal.close()
print("\n[готово] BBA: поэлементный остаток и потолок замыкания", flush=True)
