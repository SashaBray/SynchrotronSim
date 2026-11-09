# -*- coding: utf-8 -*-
"""Программа-планировщик (пример): объявляем регуляторы и частоты, гоняем на DemoPlant.

Показывает весь конвейер скелета: ExperimentState + Limits + Scheduler + Controller + Resources.
Сравнение «петля ВЫКЛ → ВКЛ» на одинаковом дрейфе (как наши реальные прогоны), таблица контуров
(многотемповость), обнаружение плагинов. Настоящий plant подключается вместо DemoPlant позже.

ЗАПУСК: py -3.14 control/demo.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control import ExperimentState, Limits, Scheduler, DemoPlant, discover
from control.controllers.integral_svd import IntegralSVD

DT = 1.0 / 10000.0          # базовый такт 10 кГц
N = 4000


def build():
    plant = DemoPlant(n_corr=12, n_bpm=40, seed=3)
    schema = plant.schema()
    limits = Limits(schema)
    for pid in schema.ids:                               # реалистичные рамки корректора
        limits.set(pid, odz=(-1e-3, 1e-3), slew=1.0)     # ±1 мрад, slew 1 рад/с
    state = ExperimentState(schema, limits=limits, history=plant.history())
    sch = Scheduler(plant, state, DT, log=["bpm"])
    return plant, sch


def rms_tail(series):
    a = np.asarray(series[len(series) // 4:])
    return float(np.sqrt(np.nanmean(a ** 2)))


def main():
    # 1. петля ВЫКЛ (только дрейф)
    _, sch_off = build()
    sch_off.add(IntegralSVD("fofb", gain=0.5), rate_hz=10000, enabled=False)
    off = sch_off.run(N)
    r_off = rms_tail(off["bpm"])

    # 2. петля ВКЛ
    _, sch_on = build()
    sch_on.add(IntegralSVD("fofb", gain=0.5), rate_hz=10000, enabled=True)
    sch_on.add(IntegralSVD("sofb", gain=0.1), rate_hz=1000, enabled=False)   # второй контур, др. темп
    on = sch_on.run(N)
    r_on = rms_tail(on["bpm"])

    print("=== СКЕЛЕТ ПЛАТФОРМЫ УПРАВЛЕНИЯ — демо на DemoPlant ===")
    print(f"базовый такт {1/DT/1e3:.0f} кГц, {N} тактов; корректоров 12, BPM 40")
    print(f"\nорбита BPM (СКО хвоста): петля ВЫКЛ {r_off:.3e} → ВКЛ {r_on:.3e}  ×{r_off/max(r_on,1e-30):.1f}")
    print("\nтаблица контуров (имя, частота Гц, период тактов, вкл):")
    for name, rate, period, en in sch_on.table():
        print(f"   {name:6} | {rate:7.0f} | {period:5d} | {'да' if en else 'нет'}")

    plugins = discover(os.path.join(os.path.dirname(__file__), "controllers"))
    print(f"\nнайдено плагинов-регуляторов в control/controllers: {list(plugins)}")
    print("\n[готово] скелет бегает: сумма по слотам, Limits, многотемповый планировщик, плагины")


if __name__ == "__main__":
    main()
