# -*- coding: utf-8 -*-
"""Режим 3 через платформу: охарактеризовать машину, которую оставил опыт с ОС.

Один слой анализа (control/analysis) считает динамическую апертуру/частоты для ТРЁХ конфигураций той же
решётки СИЛА-26: проектная, с искажённой орбитой (feed-down квадруполей), и скорректированная (орбита
сведена корректорами). Смысл: орбитное искажение НЕ двигает линейную частоту (кики не меняют фокусировку),
но через feed-down сквозь секступоли СЪЕДАЕТ динамическую апертуру — а коррекция орбиты её ВОЗВРАЩАЕТ.
Так режим 3 (анализ) смыкается с режимами ОС: анализируем ровно ту машину, что получилась после коррекции.

Кики искажения/коррекции берём линейной моделью RingModel (быстро, детерминированно) — тот же набор,
что отдал бы plant.kick_lists после опыта с FOFB. ЗАПУСК: py -3.14 control/mode3_demo.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import facility_report as fr
import ring_analysis as ra

from control import analysis as an
from control.ring_model import RingModel

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
TURNS = 128 if SMOKE else 256
TARGET_ORBIT = 200e-6                                         # целевое СКО искажённой орбиты, м (видимый feed-down)
RCOND = 1e-2
OUT = os.path.join("workspace", "results", "screens", "mode3_da.png")


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def orbit_rms(ring, mon):
    co = ra.closed_orbit(ring, dims=4)
    snap = np.asarray(ring.track_elements([co.tolist()]))[0]
    return rms(snap[mon, 0]), rms(snap[mon, 2])


def main():
    print("=== РЕЖИМ 3 через платформу: ДА проектной / искажённой / скорректированной машины ===")
    print(f"решётка 26, {TURNS} оборотов; искажение орбиты ~{TARGET_ORBIT*1e6:.0f} мкм (feed-down квадруполей)\n")

    # --- кики искажения и коррекции из линейной модели (тот же набор, что дал бы plant.kick_lists)
    rm = RingModel(lattice=26, nsl=8)
    corr, quads, K1L, mon = rm.corr, rm.quads, np.asarray(rm.K1L), rm.mon
    Hc = rm.response(rm.hk_ids, ["ebpm.x"])["ebpm.x"]         # корректор→орбита
    Gx = rm.response(rm.qx_ids, ["ebpm.x"])["ebpm.x"]         # сдвиг квадруполя→орбита (feed-down, K1L внутри)
    rm.close()

    rp = np.random.default_rng(11); pat = rp.normal(0, 1, len(quads)); pat /= np.std(pat)
    amp = 1e-6
    amp *= TARGET_ORBIT / max(rms(Gx @ (amp * pat)), 1e-30)   # масштаб под целевую орбиту (линейно по amp)
    dxpat = amp * pat
    orbit_mon = Gx @ dxpat                                    # орбита на мониторах от искажения (линейная модель)
    u = -np.linalg.pinv(Hc, rcond=RCOND) @ orbit_mon          # одношаговая SVD-коррекция (как делает контур)

    dist_k = [(int(quads[j]), float(K1L[j] * dxpat[j])) for j in range(len(quads))]   # feed-down кики
    corr_k = [(int(corr[j]), float(u[j])) for j in range(len(corr))]                   # кики корректоров
    configs = {"проектная": [], "искажённая": dist_k, "скорректированная": dist_k + corr_k}

    # --- один анализ-слой на одном кольце, apply_perturb между конфигурациями
    inf = dict(fr.info_of(26)); inf["_id"] = 26
    ring = fr._ring(inf, nsl=8)
    rows = {}; da_by = {}
    try:
        for name, kk in configs.items():
            ring.apply_perturb(kicks=kk)
            ox, oy = orbit_rms(ring, mon)
            r = an.analyze(ring, inf, turns=TURNS, do=("optics", "da"))
            m = an.da_metrics(r["da"])
            rows[name] = (ox, oy, r["tunes"], m)
            da_by[name] = r["da"]
    finally:
        ring.apply_perturb(); ring.close()

    print(f"{'конфигурация':18} | {'орбита X/Y, мкм':>16} | {'частоты νx/νy':>16} | {'ДА x±, мм':>12} | {'ДА площадь, мм²':>15}")
    for name, (ox, oy, (nux, nuy), m) in rows.items():
        print(f"{name:18} | {ox*1e6:7.1f}/{oy*1e6:6.1f} | {nux:8.4f}/{nuy:7.4f} | "
              f"{m['x_pos']:5.2f}/{m['x_neg']:5.2f} | {m['area_mm2']:15.1f}")

    ad, ai, ac = (rows['проектная'][3]['area_mm2'], rows['искажённая'][3]['area_mm2'],
                  rows['скорректированная'][3]['area_mm2'])
    print(f"\nплощадь ДА: проект {ad:.1f} → искажённая {ai:.1f} ({ai/ad*100:.0f}% проекта) "
          f"→ скорректированная {ac:.1f} ({ac/ad*100:.0f}% проекта)")
    print("частоты не двигаются (кики не меняют фокусировку) — эффект искажения чисто НЕЛИНЕЙНЫЙ, в ДА")

    if an.save_da_png(da_by, OUT):
        print(f"\nконтуры ДА: {OUT}")
    print("\n[готово] слой анализа (режим 3) характеризует машину после опыта с ОС: тот же движок, та же решётка")


if __name__ == "__main__":
    main()
