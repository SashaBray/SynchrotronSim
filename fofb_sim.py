# -*- coding: utf-8 -*-
"""ДИНАМИКА v3: ПЕР-ОБОРОТНЫЙ FOFB (быстрая петля) — развёртка по частоте, граница полосы.

Постановка (исправлена после диагностики живучести + бага переноса z):
  * возмущение — ТОН частоты f на позициях квадруполей (A·sin(2πft), фикс. пространственный паттерн):
    реалистичная вибрация/сетевая наводка; орбита колеблется на f через feed-down. НЕ ступень
    (ступень 3мкм давала 200мкм орбиты + когерентный вздрог 220мкм — нефизичная для FOFB постановка).
  * радиационное кольцо (sr+ВЧ, фаза 165°) — пучок живёт тысячи оборотов (проверено).
  * пучок ВЕДЁМ track_tbt чанками по STRIDE об (самосогласован, корректный z); орбиту на мониторах/
    источнике читаем СНИМКОМ track_elements раз в чанк (его выход НЕ переносим — у track_elements z =
    накопленный путь, +периметр/об, перенос раздувал ВЧ-фазу и ронял пучок).
  * такт FOFB 135 об ≈ 2 кГц (= декогеренция ~140 об, предел квазистатики): BPM = орбита мониторов,
    усреднённая по всему такту (betatron усреднён); корректоры интегрально, с задержкой LAT тактов.
  * амплитуду тона на f меряем ЛОК-ИНОМ (корреляция с sin/cos на f) по устоявшемуся окну (~1 период).

Смысл: честная S(f) петли из ПОЛНОГО трекинга. 50 Гц (5400 об/период ≫ задержки петли) — давится;
гармоники 150/250 Гц попадают в горб водяного матраса (S>1, усиление) — их медленной петлёй не взять.
Запуск: py -3.14 fofb_sim.py   |   дымовой тест: py -3.14 fofb_sim.py smoke
"""
import os, sys, csv, math, time
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import ring_cpp, facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
_argf = [float(a) for a in sys.argv[1:] if a.replace(".", "").isdigit()]   # частоты из argv (точечный прогон)
FREQS = [150.0] if SMOKE else (_argf or [50.0, 150.0, 250.0, 500.0, 1000.0])  # 50 — сеть; 150/250 — её гармоники
TACT, GAIN, LAT = 135, 0.8, 0          # такт=135 об≈2 кГц (= декогеренция, предел квазистатики), усиление, задержка
STRIDE = 15                            # тон+корректоры обновляем раз в STRIDE об (apply_perturb дорог ~113мс)
A_POS, NP = 0.3e-6, 12                 # амплитуда тона позиции квадруполей (м); частиц в пучке
inf = dict(fr.info_of(26)); inf["_id"] = 26

# 4D-кольцо: индексы, матрица отклика, feed-down, геометрия источника
r4 = fr._ring(inf, nsl=8); resp = oc.cached_response_matrix(r4, 26); corr = resp["corr"]
mon = oi.monitor_rows(r4)
Pmx = np.linalg.pinv(resp["Rx"][mon - 1], rcond=1e-2); Pmy = np.linalg.pinv(resp["Ry"][mon - 1], rcond=1e-2)
# rcond=1e-2 — компромисс авторитет/устойчивость: 1e-3 включает слабые моды SVD → огромные кики → петля расходится
quads = oi.quad_indices(r4, kinds=("quad",))
K1L = np.array([float(r4.mags[q].k1) * float(r4.mags[q].Leff) for q in quads])
s_out = np.array(r4.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r4.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0
FREV = float(getattr(r4, "circumference", 1110.0)); FREV = 2.99792458e8 / FREV; TREV = 1.0 / FREV
r4.close()

# пространственный паттерн тона (фикс. случайные знаки, СКО = A_POS), колеблется во времени как sin(2πft)
rp = np.random.default_rng(11)
PATX = rp.normal(0, 1, len(quads)); PATX *= A_POS / np.std(PATX)
PATY = rp.normal(0, 1, len(quads)); PATY *= A_POS / np.std(PATY)

ring = ring_cpp.CppRing(lattice=26, sr=True,
    rf={"volt_mv": inf["rf"][0], "harmonic": inf["rf"][1], "phase_deg": 165.0},
    design_mev=inf["energy_mev"], nsl=8)
CO = ra.closed_orbit(ring, dims=6)
rb = np.random.default_rng(5)
BUNCH0 = np.array([CO + np.array([rb.normal(0, 5e-6), 0, rb.normal(0, 5e-6), 0, 0, rb.normal(0, 1e-3)])
                   for _ in range(NP)])
print(f"FOFB turn-by-turn: f_об={FREV/1e3:.1f} кГц, такт {TACT} об≈{FREV/TACT/1e3:.2f} кГц, "
      f"усиление {GAIN}, задержка {LAT}; тон A={A_POS*1e6:.2f} мкм, пучок {NP}ч; частоты {FREQS}", flush=True)


DTCH = STRIDE * TREV                                     # шаг дискретизации сигнала (один чанк), с


def lockin(sig, f, i0):
    """Амплитуда компоненты частоты f в sig (дискретизован по чанкам) на окне [i0:] (лок-ин)."""
    w = sig[i0:]; fin = np.isfinite(w)
    if fin.sum() < 8:
        return float("nan")
    idx = np.arange(i0, i0 + len(w))[fin]; w = w[fin] - np.mean(w[fin]); t = idx * DTCH
    c = np.sum(w * np.cos(2 * np.pi * f * t)); s = np.sum(w * np.sin(2 * np.pi * f * t))
    return 2.0 * math.hypot(c, s) / len(w)


def plan(freq):
    per = FREV / freq                                   # оборотов в периоде тона
    settle = 1620                                       # 6 тактов — хватает для установления петли и декогеренции
    meas = int(max(per, 3 * TACT))                      # не меньше одного периода тона
    return settle, settle + meas


def run(freq, feedback):
    """Пучок ведём track_tbt чанками по STRIDE об (самосогласован, корректный z); орбиту на мониторах/
    источнике читаем СНИМКОМ track_elements раз в чанк (его выход НЕ переносим — так нет бага накопления z).
    Корректоры держим в такте, обновляем раз в TACT по BPM устоявшегося хвоста. Возврат: сигнал источника
    (по чанкам), индекс конца установления (в чанках), флаг устойчивости."""
    settle, nturn = (540, 1620) if SMOKE else plan(freq)
    nch = nturn // STRIDE                                            # число чанков
    bunch = BUNCH0.copy()
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))
    tgtx, tgty = thx.copy(), thy.copy()
    buf = [(thx.copy(), thy.copy()) for _ in range(LAT + 1)]
    sig = np.full(nch, np.nan); monx = np.full((nch, len(mon)), np.nan)   # источник; орбита по мониторам (гор.)
    accx = np.zeros(len(mon)); accy = np.zeros(len(mon)); nacc = 0
    tpc = TACT // STRIDE                                             # чанков в такте (135/15=9)
    for c in range(nch):
        n0 = c * STRIDE
        osc = math.sin(2 * math.pi * freq * (n0 * TREV))
        dkx = K1L * (PATX * osc); dky = -K1L * (PATY * osc)          # feed-down тона (пост. внутри чанка)
        kk = [(corr[j], thx[j]) for j in range(len(corr))] + [(quads[j], dkx[j]) for j in range(len(quads))]
        vv = [(corr[j], thy[j]) for j in range(len(corr))] + [(quads[j], dky[j]) for j in range(len(quads))]
        ring.apply_perturb(kicks=kk, vkicks=vv)
        snap = np.mean(np.array(ring.track_elements(bunch.tolist())), axis=0)   # центроид по линии (read-only)
        sig[c] = float(np.mean(snap[SIDX, 0] + snap[SIDX, 1] * DVEC))
        monx[c] = snap[mon, 0]                                       # орбита в каждом мониторе (для по-помониторного лок-ина)
        adv = ring.track_tbt(bunch.tolist(), STRIDE)                 # (NP, STRIDE+1, 6) — ведём пучок
        bunch = adv[:, STRIDE, :]
        if not np.all(np.isfinite(bunch)):
            return sig, monx, settle // STRIDE, False
        accx += snap[mon, 0]; accy += snap[mon, 2]; nacc += 1        # BPM по всему такту (betatron ~ усреднён)
        if (c + 1) % tpc == 0:
            if feedback and nacc:
                tgtx = thx - GAIN * (Pmx @ (accx / nacc)); tgty = thy - GAIN * (Pmy @ (accy / nacc))
            buf.append((tgtx.copy(), tgty.copy()))
            thx, thy = buf.pop(0)                                    # применяем с задержкой LAT тактов
            accx[:] = 0; accy[:] = 0; nacc = 0
    return sig, monx, settle // STRIDE, True


def bpm_tone_rms(monx, f, i0):
    """СКО (по мониторам) тоновой амплитуды орбиты на частоте f — лок-ин каждого монитора."""
    amps = np.array([lockin(monx[:, j], f, i0) for j in range(monx.shape[1])])
    amps = amps[np.isfinite(amps)]
    return float(np.sqrt(np.mean(amps ** 2))) if len(amps) else float("nan")


if SMOKE:
    t0 = time.time(); off, ob, st, _ = run(150.0, False); on, nb, st2, ok = run(150.0, True); dt = time.time() - t0
    print(f"[smoke] 150 Гц, 2×1620 об (STRIDE={STRIDE}): {dt:.0f} с ({dt/3240*1e3:.0f} мс/об), устойчив={ok}; "
          f"источник выкл={lockin(off,150.0,st)*1e6:.3f} вкл={lockin(on,150.0,st)*1e6:.3f} мкм; "
          f"BPM-тон(СКО) выкл={bpm_tone_rms(ob,150.0,st)*1e6:.3f} вкл={bpm_tone_rms(nb,150.0,st2)*1e6:.3f} мкм", flush=True)
    ring.apply_perturb(); ring.close(); sys.exit(0)

rows = []
print(f"\n{'f,Гц':>6} | {'период':>6} | источник A, мкм (выкл→вкл)  подавл | BPM-тон СКО, мкм (выкл→вкл)  подавл | уст")
for f in FREQS:
    t0 = time.time()
    off, ob, st, ok0 = run(f, False); on, nb, _, ok1 = run(f, True)
    ao = lockin(off, f, st); an = lockin(on, f, st); supp = ao / an if (an and an > 0) else float("nan")
    bo = bpm_tone_rms(ob, f, st); bn = bpm_tone_rms(nb, f, st); bsupp = bo / bn if (bn and bn > 0) else float("nan")
    rows.append((f, FREV / f, ao, an, supp, bo, bn, bsupp, ok0 and ok1))
    print(f"{f:6.0f} | {FREV/f:6.0f} | {ao*1e6:7.3f} → {an*1e6:6.3f}   ×{supp:4.1f} | "
          f"{bo*1e6:7.3f} → {bn*1e6:6.3f}   ×{bsupp:4.1f} | {'да' if ok0 and ok1 else 'НЕТ'}"
          f"  [{time.time()-t0:.0f}с]", flush=True)
ring.apply_perturb(); ring.close()
out = os.path.join("workspace", "results", "sensitivity_study", "fofb_sweep_exp26.csv")
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter=";")
    w.writerow(["freq_hz", "period_turns", "src_off_um", "src_on_um", "src_supp", "bpm_off_um", "bpm_on_um", "bpm_supp"])
    for f, per, ao, an, supp, bo, bn, bsupp, _ in rows:
        w.writerow([f"{f:.0f}", f"{per:.0f}", f"{ao*1e6:.4f}", f"{an*1e6:.4f}", f"{supp:.3f}",
                    f"{bo*1e6:.4f}", f"{bn*1e6:.4f}", f"{bsupp:.3f}"])
print(f"\nCSV: {out}\n[готово] источник упирается в ПОЛ утечки BPM→источник, BPM-орбита давится полосой петли", flush=True)
