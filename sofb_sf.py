# -*- coding: utf-8 -*-
"""ДИНАМИКА v2: частотный паспорт медленной ОС — функция подавления S(f) + задержка петли.

Усложнение v1: (а) ПЕРИОДИЧЕСКОЕ возмущение (фикс. пространственный паттерн × sin 2πft); (б) ЗАДЕРЖКА
регулятора L тактов (замер на такте n применяется на n+L — главный драйвер устойчивости); (в) замер
амплитуды колебаний источника на частоте f при петле ВКЛ/ВЫКЛ → S(f)=A_вкл/A_выкл.

Квазистатика (такт ≫ время затухания) → замкнутая орбита на такт; вся динамика в дискретной петле.
Разрешимо до частоты Найквиста 1/(2·dt)=5 Гц (медленная петля; выше — забота FOFB). Пишет CSV.
Запуск: py -3.14 sofb_sf.py
"""
import os, sys, csv
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

DT = 0.1                    # с/такт (10 Гц)
GAIN = 0.3
BPM_NOISE = 0.0            # для чистого S(f) шум выкл (линейный отклик)
AMP = 1e-6                 # м, амплитуда дрейф-паттерна (линейно: S(f) от амплитуды не зависит)
FREQS = [0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2, 2.0, 3.0]
LATENCIES = [0, 1, 3]      # тактов задержки регулятора
inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
resp = oc.cached_response_matrix(r, 26); corr = resp["corr"]
mon = oi.monitor_rows(r)
Pmx = np.linalg.pinv(resp["Rx"][mon - 1], rcond=1e-2); Pmy = np.linalg.pinv(resp["Ry"][mon - 1], rcond=1e-2)
quads = oi.quad_indices(r, kinds=("quad",))
K1L = np.array([float(r.mags[q].k1) * float(r.mags[q].Leff) for q in quads])
s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
# фикс. пространственный паттерн дрейфа (один и тот же во всех опытах)
rng = np.random.default_rng(7)
DXP = rng.normal(0, AMP, len(quads)); DYP = rng.normal(0, AMP, len(quads))
print(f"S(f) медленной ОС: {len(FREQS)} частот, задержки {LATENCIES} тактов, усиление {GAIN}, dt={DT*1e3:.0f} мс (Найквист {1/(2*DT):.0f} Гц)", flush=True)


def source_mean_x(dxq, dyq, thx, thy):
    kk = [(corr[j], thx[j]) for j in range(len(corr))] + [(quads[j], K1L[j] * dxq[j]) for j in range(len(quads))]
    vv = [(corr[j], thy[j]) for j in range(len(corr))] + [(quads[j], -K1L[j] * dyq[j]) for j in range(len(quads))]
    r.apply_perturb(kicks=kk, vkicks=vv)
    co = ra.closed_orbit(r, dims=4)
    if not np.all(np.isfinite(co)):
        return None, None
    o = r.track_elements([co])[0]
    return float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC)), o   # знаковый скаляр источника


def run_tone(f, latency, feedback):
    per = 1.0 / (f * DT)                                   # тактов в периоде
    ncyc = int(min(600, max(60, 3.0 * per + 40)))
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))
    buf = [(np.zeros(len(corr)), np.zeros(len(corr)))] * (latency + 1)   # очередь целей (задержка)
    s_hist = []
    for n in range(ncyc):
        ph = np.sin(2 * np.pi * f * n * DT)
        dxq = DXP * ph; dyq = DYP * ph
        tgtx, tgty = buf.pop(0)                            # цель, посчитанная latency тактов назад
        thx, thy = tgtx.copy(), tgty.copy()                # slew не связываем (чистый S(f) от gain/задержки)
        s, o = source_mean_x(dxq, dyq, thx, thy)
        s_hist.append(s if s is not None else 0.0)
        if feedback and o is not None:
            ox = o[mon, 0]; oy = o[mon, 2]
            if BPM_NOISE:
                ox = ox + rng.normal(0, BPM_NOISE, len(mon)); oy = oy + rng.normal(0, BPM_NOISE, len(mon))
            buf.append((thx - GAIN * (Pmx @ ox), thy - GAIN * (Pmy @ oy)))
        else:
            buf.append((thx, thy))
    tail = s_hist[-int(max(per, 10)):]                     # последний период
    return (max(tail) - min(tail)) / 2.0                   # амплитуда колебаний источника


A_off = run_tone(0.1, 0, False)                            # разомкнутая амплитуда (квазистатика ⇒ от f не зависит)
print(f"разомкнутая амплитуда источника A_выкл = {A_off*1e6:.2f} мкм\n", flush=True)
print(f"{'f, Гц':>7} |" + "".join(f" S(L={L})" for L in LATENCIES), flush=True)
rows = []
for f in FREQS:
    line = f"{f:7.2f} |"; rec = [f]
    for L in LATENCIES:
        A_on = run_tone(f, L, True); S = A_on / A_off
        line += f"  {S:6.2f}"; rec.append(S)
    print(line, flush=True); rows.append(rec)
r.apply_perturb(); r.close()
out = os.path.join("workspace", "results", "sensitivity_study", "sofb_sf_exp26.csv")
with open(out, "w", newline="", encoding="utf-8") as fcsv:
    w = csv.writer(fcsv, delimiter=";"); w.writerow(["f_Hz"] + [f"S_L{L}" for L in LATENCIES]); w.writerows(rows)
print(f"\nCSV: {out}\n[готово] S<1 — подавление, S≈1 — петля не действует, S>1 — усиление (горб неустойчивости)", flush=True)
