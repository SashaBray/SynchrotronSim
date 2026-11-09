# -*- coding: utf-8 -*-
"""Каскад B через платформу: быстрый FOFB (e-BPM, кГц) + медленный XBPM-губернатор (фотон, суб-кГц).

Демонстрирует многотемповость, сумму по слотам и уставку между контурами (golden orbit) на реальной
задаче. Ключевой смысл: медленный увод НАВЕДЕНИЯ ИД виден фотону (XBPM), но НЕ e-BPM — его берёт только
внешний контур, а быстрый при этом держит орбиту и ему не мешает.

Три режима: без контуров / только FOFB / FOFB+XBPM. Смотрим СКО орбиты e⁻ и положения фотона.
ЗАПУСК: py -3.14 control/cascade_demo.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control import ExperimentState, Limits, Scheduler
from control.cascade_plant import CascadePlant
from control.controllers.integral_svd import IntegralSVD
from control.controllers.xbpm_governor import XBPMGovernor

DT = 1.0 / 10000.0                                          # базовый такт 10 кГц
N = 20000


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def tail(series):
    return rms(series[len(series) // 4:])


def run(mode):
    """Пошаговый прогон с ЧЕСТНЫМИ метриками: ошибка слежения FOFB (ebpm−golden), смещение golden
    (намеренная цена наведения фотона) и СКО фотона."""
    plant = CascadePlant(n_photon=8, photon_d=15.0, seed=5)
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, DT)
    sch.add(IntegralSVD("fofb", sensor="ebpm.x", actuator_key="hkick", gain=0.5, ref="golden.x"),
            rate_hz=10000, enabled=(mode in ("fofb", "cascade")))
    if mode == "cascade":
        sch.add(XBPMGovernor("xbpm", sensor="photon.x", actuator_key="hkick",
                             ref_sensor="ebpm.x", ref_name="golden.x", gain=0.5, rcond=1e-4), rate_hz=200)
    track, gold, ph = [], [], []
    for _ in range(N):
        sch.step()
        o = state.observables["ebpm.x"]
        r = state.refs.get("golden.x")
        track.append(rms(o if r is None else o - r))         # ошибка слежения FOFB
        gold.append(0.0 if r is None else rms(r))            # смещение golden orbit (цена)
        ph.append(rms(state.observables["photon.x"]))        # фотон (цель)
    return tail(track), tail(gold), tail(ph), sch


def main():
    print("=== КАСКАД B через платформу: FOFB (10 кГц) + XBPM-губернатор (200 Гц) ===")
    print(f"{N} тактов; быстрое возмущение (видит e-BPM) + медленный увод наведения ИД (видит ТОЛЬКО фотон)\n")
    print(f"{'режим':14} | {'слежение FOFB':>13} | {'golden-смещ.':>12} | {'фотон':>8}  (мкм СКО)")
    rows = {}
    sch_last = None
    for mode, label in [("off", "без контуров"), ("fofb", "только FOFB"), ("cascade", "FOFB + XBPM")]:
        tr, gd, ph, sch_last = run(mode)
        rows[mode] = (tr, gd, ph)
        print(f"{label:14} | {tr*1e6:13.3f} | {gd*1e6:12.3f} | {ph*1e6:8.3f}")

    ph_f = rows["fofb"][2]; ph_c = rows["cascade"][2]
    tr_f = rows["fofb"][0]; tr_c = rows["cascade"][0]; gd_c = rows["cascade"][1]
    print(f"\nфотон: только FOFB {ph_f*1e6:.3f} → +XBPM {ph_c*1e6:.3f} мкм  (×{ph_f/max(ph_c,1e-30):.1f}) — "
          f"медленный контур берёт увод ИД, которого e-BPM не видит")
    print(f"слежение FOFB: {tr_f*1e6:.3f} → {tr_c*1e6:.3f} мкм — быстрый контур держит орбиту НА golden одинаково")
    print(f"цена: golden-смещение вырастает до {gd_c*1e6:.3f} мкм — намеренный сдвиг орбиты ради наведения фотона")
    print("\nтаблица контуров (имя, частота Гц, период тактов, вкл):")
    for name, rate, period, en in sch_last.table():
        print(f"   {name:6} | {rate:7.0f} | {period:6d} | {'да' if en else 'нет'}")
    print("\n[готово] многотемповый каскад: сумма по слотам + уставка golden orbit между контурами")


if __name__ == "__main__":
    main()
