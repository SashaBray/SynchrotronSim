# -*- coding: utf-8 -*-
"""ПУНКТ 1: рябь ПОЛЯ магнитов (имитация тока) → движение источника. Проверяет канал field:* модели
disturbance.py и показывает, как рябь РАЗНЫХ типов/компонент по-разному доходит до источника.
Для каждой (тип × компонента × корреляция) задаём относительную рябь 1e-4 (50 Гц, пик) и меряем
смещение источника квазистатически (замкнутая орбита). Запуск: py -3.14 field_noise_demo.py
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, ring_analysis as ra
import disturbance as dist

inf = dict(fr.info_of(26)); inf["_id"] = 26
r = fr._ring(inf, nsl=8)
s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0


def source(o):
    return float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC))


r.apply_perturb()
base = source(r.track_elements([ra.closed_orbit(r, dims=4)])[0])
REL = 1e-4; F = 50.0; TPK = 1.0 / (4 * F)                     # пик тона (sin=1)
print(f"рябь поля 1e-4 (50 Гц, пик); базовый источник {base*1e6:.3f} мкм\n", flush=True)
print(f"{'тип×компонента':>18} | {'корреляция':>11} | {'|Δисточник|, мкм':>16}")
cases = [("dipole", "b1", "common"), ("dipole", "b1", "independent"),
         ("quad", "k1", "common"), ("quad", "k1", "independent"),
         ("sext", "k2", "independent"), ("oct", "k3", "independent")]
for target, comp, corr in cases:
    vals = []
    for seed in range(8 if corr == "independent" else 1):
        row = dict(name="t", target=target, plane="", form="tone", correlation=corr,
                   quantity=f"field:{comp}", amplitude=REL, freq_hz=F, phase_deg=0.0)
        D = dist.Disturbance([row], r, seed=seed)
        fd = D.field(TPK, 0.0)                                # {арг apply_perturb: ΔK по магнитам}
        kw = {key: [(i, float(v)) for i, v in enumerate(arr) if v != 0.0] for key, arr in fd.items()}
        r.apply_perturb(**kw)
        o = r.track_elements([ra.closed_orbit(r, dims=4)])[0]
        vals.append(abs(source(o) - base))
    val = float(np.sqrt(np.mean(np.square(vals))))
    print(f"{target+'·'+comp:>18} | {corr:>11} | {val*1e6:16.4f}", flush=True)
r.apply_perturb(); r.close()
print("\n[готово] ожидаем: диполь b1 (independent) — главный орбитный эффект; common диполь мал у "
      "бездисперсного источника; квадруполь k1 — через feed-down; sext/oct — почти ноль на орбите", flush=True)
