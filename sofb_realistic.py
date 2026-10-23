# -*- coding: utf-8 -*-
"""ДИНАМИКА (пункт b): SOFB с РЕАЛИСТИЧНЫМ возмущением (корреляция + спектр из БД disturbance_spec.csv).

Многомасштабное время: мелкий шаг dt_sim (разрешает вибрации) для возмущения+наблюдения; регулятор
SOFB реже (такт N_sim шагов = 10 Гц). Возмущение — модель disturbance.Disturbance (тон/ломаная ×
independent/girder/ground). Позиции квадруполей → орбита feed-down киком. Смотрим СКО источника и
СПЕКТР остатка (off vs on): медленный дрейф петля давит, вибрации выше полосы — проходят (нужен FOFB).
Запуск: py -3.14 sofb_realistic.py
"""
import os, sys, csv
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra
import disturbance as dist

DT_SIM, T_TOTAL = 2e-3, 3.0         # шаг 2 мс (Найквист 250 Гц), опыт 3 с
LOOP_EVERY = 50                     # такт регулятора = 50 шагов = 100 мс = 10 Гц
GAIN, BPM_NOISE = 0.3, 0.1e-6
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
spec = dist.load_spec()
gid, ng = dist.girder_groups(r)
nstep = int(T_TOTAL / DT_SIM)
print(f"SOFB реалистичный: шаг {DT_SIM*1e3:.0f} мс × {nstep} ({T_TOTAL} с), такт 10 Гц; опор-гирдеров {ng}", flush=True)
print("спека возмущений:", ", ".join(f"{c['name']}({c['form']},{c['correlation']},{c['freq_hz']}Гц)" for c in spec), flush=True)


def src_x(thx, thy, dqx, dqy):
    kk = [(corr[j], thx[j]) for j in range(len(corr))] + [(quads[j], K1L[j] * dqx[j]) for j in range(len(quads))]
    vv = [(corr[j], thy[j]) for j in range(len(corr))] + [(quads[j], -K1L[j] * dqy[j]) for j in range(len(quads))]
    r.apply_perturb(kicks=kk, vkicks=vv)
    co = ra.closed_orbit(r, dims=4)
    if not np.all(np.isfinite(co)):
        return None, None
    o = r.track_elements([co])[0]
    return float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC)), o


def run(feedback):
    D = dist.Disturbance(spec, r, seed=3)                   # одна реализация на оба прогона
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr)); tgtx = thx.copy(); tgty = thy.copy()
    sig = []; accx = np.zeros(len(mon)); accy = np.zeros(len(mon)); nacc = 0
    for n in range(nstep):
        t = n * DT_SIM
        dx, dy = D.offsets(t, DT_SIM)                       # смещения всех магнитов
        thx, thy = tgtx, tgty                               # привод (slew не связываем здесь)
        sx, o = src_x(thx, thy, dx[quads], dy[quads])
        sig.append(sx if sx is not None else 0.0)
        if o is not None:
            accx += o[mon, 0]; accy += o[mon, 2]; nacc += 1
        if feedback and (n + 1) % LOOP_EVERY == 0 and nacc:
            ox = accx / nacc + np.random.default_rng(n).normal(0, BPM_NOISE, len(mon))
            oy = accy / nacc + np.random.default_rng(n + 1).normal(0, BPM_NOISE, len(mon))
            tgtx = thx - GAIN * (Pmx @ ox); tgty = thy - GAIN * (Pmy @ oy)
            accx[:] = 0; accy[:] = 0; nacc = 0
    return np.array(sig)


print("прогон: петля ВЫКЛ ...", flush=True); off = run(False)
print("прогон: петля ВКЛ ...", flush=True); on = run(True)
r.apply_perturb(); r.close()
fs = 1.0 / DT_SIM
def band(sig, f0, f1):
    F = np.fft.rfft(sig - np.mean(sig)); fr_ = np.fft.rfftfreq(len(sig), DT_SIM)
    m = (fr_ >= f0) & (fr_ < f1)
    return float(np.sqrt(np.sum(np.abs(F[m]) ** 2)) / len(sig)) * 1e6 * 2   # ~амплитуда в полосе, мкм
rms = lambda s: float(np.std(s)) * 1e6
out = os.path.join("workspace", "results", "sensitivity_study", "sofb_realistic_exp26.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter=";"); w.writerow(["t_s", "off_um", "on_um"])
    for i in range(len(off)):
        w.writerow([f"{i*DT_SIM:.4f}", f"{off[i]*1e6:.4f}", f"{on[i]*1e6:.4f}"])
print(f"\nСКО источника x: ВЫКЛ {rms(off):.2f} → ВКЛ {rms(on):.2f} мкм  (×{rms(off)/max(rms(on),1e-9):.1f})")
print(f"\nпо полосам (амплитуда, мкм):   {'ВЫКЛ':>8} {'ВКЛ':>8}  трактовка")
for f0, f1, nm in [(0, 0.5, "дрейф <0.5 Гц"), (10, 14, "вибрация ~12 Гц"), (38, 42, "тон ~40 Гц")]:
    print(f"  {nm:16} {band(off,f0,f1):8.3f} {band(on,f0,f1):8.3f}", flush=True)
print(f"\nCSV: {out}\n[готово] ожидаем: дрейф петля давит, вибрации/тон выше полосы проходят (нужен FOFB)", flush=True)
