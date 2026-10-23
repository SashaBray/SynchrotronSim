# -*- coding: utf-8 -*-
"""ДИНАМИКА v1: медленная орбитная обратная связь (SOFB), КВАЗИСТАТИКА во времени.

Первый срез модельного времени. Петля ~Гц: такт (27 000 об @10 Гц) ≫ время затухания (5800 об) ⇒ пучок
на каждом такте в своей замкнутой орбите → трекать пооборотно НЕ нужно, считаем замкнутую орбиту на такт.

Эффективность: всё кик-вектором на ОДНОМ идеальном кольце, без пересборок.
  • ВОЗМУЩЕНИЕ (дрейф позиций квадруполей dx,dy(t)) ≡ feed-down кики: θx=K1·L·dx, θy=−K1·L·dy
    (ядро применяет кик и на полевых элементах — фикс из orbit_infer). Модель дрейфа — случайная ломаная
    (кусочно-постоянная скорость, рывок через случайный интервал), т.е. медленное блуждание (тепло/грунт).
  • ПРИВОД — корректоры: интегральный регулятор θ←θ−g·R⁺·(орбита_BPM+шум), с рампом к цели на slew-лимите.
  • НАБЛЮДАЕМ — положение источника (центры ИД) σ(t): петля ВЫКЛ (дрейф уносит) vs ВКЛ (держит). Пишет CSV.

Параметры — аргументами/в начале. Запуск: py -3.14 sofb_sim.py
"""
import os, sys, csv
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

# --- параметры опыта ---
NCYC = 120                 # тактов регулятора
DT = 0.1                   # с/такт (10 Гц SOFB)
GAIN = 0.3                 # интегральное усиление петли
BPM_NOISE = 0.3e-6         # м, шум датчика (на чтение)
SLEW = 2e-3                # рад/с, предел скорости корректора (dθ/dt); для SOFB заведомо не связывает
DRIFT_RATE = 2e-6          # м/с, СКО скорости дрейфа позиции квадруполя (тепло/грунт)
JERK = 2.0                 # с, средний интервал смены скорости дрейфа (рывок)
SEED = 20261023
inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
resp = oc.cached_response_matrix(r, 26); corr = resp["corr"]
mon = oi.monitor_rows(r)
Rmx = resp["Rx"][mon - 1]; Rmy = resp["Ry"][mon - 1]
Pmx = np.linalg.pinv(Rmx, rcond=1e-2); Pmy = np.linalg.pinv(Rmy, rcond=1e-2)
# квадруполи-источники дрейфа: feed-down коэф K1·L (с идеала)
quads = oi.quad_indices(r, kinds=("quad",))
K1L = np.array([float(r.mags[q].k1) * float(r.mags[q].Leff) for q in quads])
# точки источника
s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
print(f"SOFB v1: {NCYC} тактов × {DT*1e3:.0f} мс, усиление {GAIN}, шум BPM {BPM_NOISE*1e6:.1f} мкм; "
      f"квадруполей-источников {len(quads)}, корректоров {len(corr)}, ИД {len(SIDX)}", flush=True)


def source_rms(thx, thy, dkx, dky):
    """СКО положения источника при киках корректоров (thx,thy) и дрейф-киках квадруполей (dkx,dky)."""
    kk = [(corr[j], thx[j]) for j in range(len(corr))] + [(quads[j], dkx[j]) for j in range(len(quads))]
    vv = [(corr[j], thy[j]) for j in range(len(corr))] + [(quads[j], dky[j]) for j in range(len(quads))]
    r.apply_perturb(kicks=kk, vkicks=vv)
    co = ra.closed_orbit(r, dims=4)
    if not np.all(np.isfinite(co)):
        return None, None, None
    o = r.track_elements([co])[0]
    x = o[SIDX, 0] + o[SIDX, 1] * DVEC; y = o[SIDX, 2] + o[SIDX, 3] * DVEC
    return float(np.sqrt(np.mean(x ** 2))), float(np.sqrt(np.mean(y ** 2))), o


def run(feedback):
    rng = np.random.default_rng(SEED)                 # ОДНА реализация дрейфа для обоих прогонов
    dx = np.zeros(len(quads)); dy = np.zeros(len(quads))        # позиции квадруполей (дрейф)
    vx = rng.normal(0, DRIFT_RATE, len(quads)); vy = rng.normal(0, DRIFT_RATE, len(quads))
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))        # уставки корректоров (привод)
    tgtx = np.zeros(len(corr)); tgty = np.zeros(len(corr))      # цель регулятора
    hist = []
    for n in range(NCYC):
        # 1. эволюция дрейфа (случайная ломаная)
        dx += vx * DT; dy += vy * DT
        reroll = rng.random(len(quads)) < (DT / JERK)           # рывок: смена скорости
        vx[reroll] = rng.normal(0, DRIFT_RATE, int(reroll.sum()))
        vy[reroll] = rng.normal(0, DRIFT_RATE, int(reroll.sum()))
        # 2. привод: рамп к цели на slew-лимите
        thx += np.clip(tgtx - thx, -SLEW * DT, SLEW * DT); thy += np.clip(tgty - thy, -SLEW * DT, SLEW * DT)
        # 3. замкнутая орбита (квазистатика): дрейф-кики + уставки корректоров
        dkx = K1L * dx; dky = -K1L * dy
        sx, sy, o = source_rms(thx, thy, dkx, dky)
        if o is None:
            hist.append((n * DT, None, None)); continue
        hist.append((n * DT, sx, sy))
        # 4. регулятор (интегральный, с шумом BPM): цель ← текущие кики − g·R⁺·(орбита+шум)
        if feedback:
            ox = o[mon, 0] + rng.normal(0, BPM_NOISE, len(mon))
            oy = o[mon, 2] + rng.normal(0, BPM_NOISE, len(mon))
            tgtx = thx - GAIN * (Pmx @ ox); tgty = thy - GAIN * (Pmy @ oy)
    return hist


print("прогон: петля ВЫКЛ ...", flush=True); off = run(False)
print("прогон: петля ВКЛ ...", flush=True); on = run(True)
r.apply_perturb(); r.close()

out = os.path.join("workspace", "results", "sensitivity_study", "sofb_sim_exp26.csv")
with open(out, "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f, delimiter=";"); w.writerow(["t_s", "off_x_um", "off_y_um", "on_x_um", "on_y_um"])
    for (t, ox, oy), (_, nx, ny) in zip(off, on):
        w.writerow([f"{t:.3f}", "" if ox is None else f"{ox*1e6:.4f}", "" if oy is None else f"{oy*1e6:.4f}",
                    "" if nx is None else f"{nx*1e6:.4f}", "" if ny is None else f"{ny*1e6:.4f}"])

def last_rms(h, i):
    v = [row[i] for row in h[len(h)//2:] if row[i] is not None]
    return float(np.sqrt(np.mean(np.square(v)))) * 1e6 if v else float("nan")
print(f"\nУстановившееся СКО источника (вторая половина опыта):")
print(f"  петля ВЫКЛ: σx={last_rms(off,1):.2f}  σy={last_rms(off,2):.2f} мкм")
print(f"  петля ВКЛ : σx={last_rms(on,1):.2f}  σy={last_rms(on,2):.2f} мкм")
print(f"  подавление: ×{last_rms(off,1)/max(last_rms(on,1),1e-9):.0f} (x), ×{last_rms(off,2)/max(last_rms(on,2),1e-9):.0f} (y)")
print(f"CSV: {out}\n[готово]", flush=True)
