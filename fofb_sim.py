# -*- coding: utf-8 -*-
"""ДИНАМИКА v3: ПЕР-ОБОРОТНЫЙ FOFB (быстрая петля, betatron-режим).

Отличие от квазистатики (SOFB v1/v2): такт < декогеренции → пучок НЕ установился, видим betatron-движение.
Центроид трекается оборот-за-оборотом (track_elements на РАДИАЦИОННОМ кольце sr+ВЧ, phase 165° → затухание и
выживание пучка). BPM = орбита на мониторах, усреднённая по такту (фильтр betatron). Регулятор раз в такт
(N оборотов) с задержкой; ступенчатое возмущение (сдвиг квадруполей) в t=0.

Наблюдаем СКО/знаковый источник оборот-за-оборотом: betatron-вздрог + коррекция FOFB + радиац. затухание.
Один центроид (консервативно для устойчивости: когерентная мода без Ландау-затухания пучка). Пишет CSV.
Запуск: py -3.14 fofb_sim.py
"""
import os, sys, csv
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import ring_cpp, facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

NCYC, NTURN = 60, 50        # тактов × оборотов/такт (такт=50 об ≈ 0.19 мс ≈ 5.4 кГц; < декогеренции ~140)
GAINS, LAT = [0.1, 0.2, 0.3, 0.5], 1    # скан усиления (граница устойчивости), задержка (тактов)
STEP = 3e-6                 # м, ступень сдвига квадруполей в t=0
inf = dict(fr.info_of(26)); inf["_id"] = 26
# 4D-кольцо: матрица отклика регулятора, feed-down, геометрия источника
r4 = fr._ring(inf, nsl=8); resp = oc.cached_response_matrix(r4, 26); corr = resp["corr"]
mon = oi.monitor_rows(r4)
Pmx = np.linalg.pinv(resp["Rx"][mon - 1], rcond=1e-2); Pmy = np.linalg.pinv(resp["Ry"][mon - 1], rcond=1e-2)
quads = oi.quad_indices(r4, kinds=("quad",))
K1L = np.array([float(r4.mags[q].k1) * float(r4.mags[q].Leff) for q in quads])
s_out = np.array(r4.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r4.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
r4.close()
rng = np.random.default_rng(11)
DXP = rng.normal(0, STEP, len(quads)); DYP = rng.normal(0, STEP, len(quads))   # паттерн ступени
# радиационное кольцо (выживающее)
ring = ring_cpp.CppRing(lattice=26, sr=True, rf={"volt_mv": inf["rf"][0], "harmonic": inf["rf"][1], "phase_deg": 165.0},
                        design_mev=inf["energy_mev"], nsl=8)
CO = ra.closed_orbit(ring, dims=6)
NP = 48                                                     # частиц в пучке (Ландау-затухание = декогеренция)
rb = np.random.default_rng(5)
BUNCH0 = np.array([CO + np.array([rb.normal(0, 5e-6), 0.0, rb.normal(0, 5e-6), 0.0, 0.0, rb.normal(0, 1e-3)])
                   for _ in range(NP)])                     # энергоразброс σδ=0.1% → хром. декогеренция
print(f"FOFB turn-by-turn: {NCYC} тактов × {NTURN} об (такт≈{NTURN*3.7:.0f} мкс≈{1e3/(NTURN*3.7):.1f} кГц), "
      f"задержка {LAT}, ступень {STEP*1e6:.0f} мкм/квадруполь; скан усиления {GAINS}", flush=True)


def src_x(o):
    return float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC))


def run(feedback, gain):
    bunch = BUNCH0.copy()
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))
    buf = [(np.zeros(len(corr)), np.zeros(len(corr)))] * (LAT + 1)
    hist = []; stable = True
    for c in range(NCYC):
        dkx = K1L * DXP; dky = -K1L * DYP                   # ступень включена с t=0 (постоянна)
        tgtx, tgty = buf.pop(0); thx, thy = tgtx.copy(), tgty.copy()
        kk = [(corr[j], thx[j]) for j in range(len(corr))] + [(quads[j], dkx[j]) for j in range(len(quads))]
        vv = [(corr[j], thy[j]) for j in range(len(corr))] + [(quads[j], dky[j]) for j in range(len(quads))]
        ring.apply_perturb(kicks=kk, vkicks=vv)             # поля постоянны внутри такта
        mon_x = np.zeros(len(mon)); mon_y = np.zeros(len(mon)); nok = 0
        for t in range(NTURN):
            res = ring.track_elements(bunch.tolist())       # пучок параллельно (OpenMP)
            finals = np.array([o[-1] for o in res])
            if not np.all(np.isfinite(finals)):
                stable = False; break
            bunch = finals
            mo = np.mean(np.array(res), axis=0)             # средняя орбита по пучку = центроид
            sx = float(np.mean(mo[SIDX, 0] + mo[SIDX, 1] * DVEC))
            if abs(sx) > 2e-3:
                stable = False; break
            hist.append((c * NTURN + t, sx))
            mon_x += mo[mon, 0]; mon_y += mo[mon, 2]; nok += 1
        if not stable:
            break
        if feedback and nok:
            ox = mon_x / nok; oy = mon_y / nok              # BPM = среднее по пучку и по такту
            buf.append((thx - gain * (Pmx @ ox), thy - gain * (Pmy @ oy)))
        else:
            buf.append((thx, thy))
    return hist, stable


def rms_tail(h):
    v = [x for _, x in h[len(h)//2:]]
    return float(np.sqrt(np.mean(np.square(v)))) * 1e6 if v else float("nan")


print("прогон: FOFB ВЫКЛ ...", flush=True)
off, _ = run(False, 0.0); off_rms = rms_tail(off)
print(f"  ВЫКЛ: СКО источника x = {off_rms:.1f} мкм\n", flush=True)
runs = {"off": off}
print(f"{'усиление':>8} | {'устойчив':>9} | {'СКО x, мкм':>11} | подавление")
for g in GAINS:
    h, st = run(True, g); runs[f"g{g}"] = h
    rr = rms_tail(h)
    tag = f"×{off_rms/max(rr,1e-9):.1f}" if st else "РАСХОДИТСЯ"
    print(f"{g:8.2f} | {'да' if st else 'НЕТ':>9} | {rr:11.1f} | {tag}", flush=True)
ring.apply_perturb(); ring.close()
out = os.path.join("workspace", "results", "sensitivity_study", "fofb_sim_exp26.csv")
cols = ["off"] + [f"g{g}" for g in GAINS]
L = max(len(v) for v in runs.values())
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter=";"); w.writerow(["turn"] + cols)
    for i in range(L):
        row = [i]
        for c in cols:
            v = runs[c]; row.append(f"{v[i][1]*1e6:.3f}" if i < len(v) else "")
        w.writerow(row)
print(f"\nCSV: {out}\n[готово] неустойчивость на высоком усилении = то, чего квазистатика не показывала", flush=True)
