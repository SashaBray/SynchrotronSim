# -*- coding: utf-8 -*-
"""Демо трёх новых контуров через платформу: КВАДРУПОЛИ (частоты), СЕКСТУПОЛИ (хроматизм), ОКТУПОЛИ (ADTS).

Вносим фракционный дрейф градиентов (k1/k2/k3 всех соотв. магнитов) → уводит νx/νy, ξx/ξy, ∂ν/∂J от
проектных. Три контура-плагина IntegralSVD (каждый со своей матрицей отклика и уставкой=проектное
значение) возвращают оптику к уставке. Показываем «цель → до (увод) → после (коррекция)».

Контуры — те же IntegralSVD, что и орбитный FOFB, только с векторными сенсорами tune/chroma/adts и
актуаторами quad./sext./oct. — подтверждает, что платформа расширяется без нового кода регулятора.

ЗАПУСК: py -3.14 control/optics_loops_demo.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control.ring_model import RingModel
from control.optics_plant import OpticsPlant
from control.state import ExperimentState
from control.scheduler import Scheduler
from control.controllers.integral_svd import IntegralSVD

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
N_TICKS = 8 if SMOKE else 14
DRIFT_K1 = 1e-3           # +0.1% дрейф градиента квадруполей → увод (дробных) частот
DRIFT_K2 = 1e-2           # +1% дрейф секступолей → увод хроматизма (меньше — слабее паразитный увод ADTS)
DRIFT_K3 = 3e-2           # +3% дрейф октуполей → увод ADTS


def fmt(v):
    return "[" + ", ".join(f"{x:+.4g}" for x in np.asarray(v)) + "]"


def main():
    print("=== ДЕМО КОНТУРОВ ОПТИКИ: квадруполи(частоты) / секступоли(хроматизм) / октуполи(ADTS) ===")
    rm = RingModel(lattice=26, nsl=8)
    plant = OpticsPlant(model=rm, observe=("tune", "chroma", "adts"),
                        drift_k1=DRIFT_K1, drift_k2=DRIFT_K2, drift_k3=DRIFT_K3, rho=1.0)

    # уставки = оптика ЧИСТОГО кольца (до возмущения); tune — ДРОБНЫЕ частоты
    nom = plant.nominal()
    tune0, chroma0, adts0 = nom["tune"], nom["chroma"], nom["adts"]
    print(f"проектные (уставки): tune(дробн.)={fmt(tune0)}  chroma={fmt(chroma0)}  adts={fmt(adts0)}")
    print(f"возмущение: Δk1={DRIFT_K1*100:.1f}% квадруполей, Δk2={DRIFT_K2*100:.1f}% секступолей, "
          f"Δk3={DRIFT_K3*100:.1f}% октуполей\n", flush=True)

    state = ExperimentState(plant.schema(), history=plant.history())
    state.refs["tune.target"] = np.asarray(tune0)
    state.refs["chroma.target"] = np.asarray(chroma0)
    state.refs["adts.target"] = np.asarray(adts0)

    sch = Scheduler(plant, state, dt=1.0)                      # reset() уже применил возмущение
    before = {k: state.observables[k].copy() for k in ("tune", "chroma", "adts")}

    # СТАДИЙНОСТЬ (как при наладке машины): сперва частоты+хроматизм, затем ADTS октуполями.
    # Секступоли сильно двигают ADTS (2-й порядок) — пока хроматизм не установлен, октуполь гонится за
    # движущейся целью. Октуполи НЕ трогают ξ1/частоты, поэтому поздний запуск ADTS устойчив.
    sch.add(IntegralSVD(name="tune", sensor="tune", actuator_key="quad.", gain=0.7, rcond=1e-2,
                        ref="tune.target"), rate_hz=1.0)
    sch.add(IntegralSVD(name="chroma", sensor="chroma", actuator_key="sext.", gain=0.7, rcond=1e-2,
                        ref="chroma.target"), rate_hz=1.0)
    sch.add(IntegralSVD(name="adts", sensor="adts", actuator_key="oct.", gain=0.7, rcond=1e-3,
                        ref="adts.target"), rate_hz=1.0, enabled=False)                 # включим во 2-й фазе

    def err(k):
        return float(np.linalg.norm(state.observables[k] - state.refs[k + ".target"]))

    # диагностика авторитета октуполей: ранг/сингулярные матрицы ADTS-отклика (3 цели × n_oct семейств)
    Ra = plant.build_response(rm.dk3_ids, ["adts"]).R
    U, sv, _ = np.linalg.svd(Ra, full_matrices=False)
    dof = int(np.sum(sv > 1e-3 * sv[0]))                       # эффективных управляемых направлений ADTS
    print(f"ADTS-отклик {Ra.shape[0]}×{Ra.shape[1]} (цели ∂νx/∂Jx,∂νy/∂Jy,∂νx/∂Jy × {len(rm.dk3_ids)} окт. семейств):")
    print(f"  сингулярные {np.array2string(sv, precision=3)}  → эффективных управляемых DOF: {dof} из 3 "
          f"(октуполи СИЛА-26 почти коллинеарны — это DA-ручка, не полный контроль ∂ν/∂J)\n", flush=True)
    Usig = U[:, :max(dof, 1)]                                  # достижимое подпространство ADTS

    n1 = N_TICKS; n2 = (6 if SMOKE else 10)
    print(f"ФАЗА 1 — частоты + хроматизм ({n1} тиков):")
    print(f"{'тик':>3} | {'‖Δtune‖':>10} | {'‖Δchroma‖':>10} | {'‖Δadts‖':>12}")
    for t in range(n1):
        sch.step()
        print(f"{t+1:3d} | {err('tune'):10.3e} | {err('chroma'):10.3e} | {err('adts'):12.3e}", flush=True)
    sch.set_enabled("adts", True)                              # ФАЗА 2: добавляем контур октуполей
    print(f"ФАЗА 2 — + ADTS октуполями ({n2} тиков):")
    for t in range(n2):
        sch.step()
        print(f"{n1+t+1:3d} | {err('tune'):10.3e} | {err('chroma'):10.3e} | {err('adts'):12.3e}", flush=True)
    after = {k: state.observables[k].copy() for k in ("tune", "chroma", "adts")}
    plant.close()

    print("\nСВОДКА (цель → до коррекции → после):")
    for k, tgt in (("tune", tune0), ("chroma", chroma0), ("adts", adts0)):
        print(f"  {k:7}: цель {fmt(tgt)}")
        print(f"          до   {fmt(before[k])}")
        print(f"          после{fmt(after[k])}")

    def red(k, basis=None):                                    # доля остатка «после/до» (в подпространстве basis)
        eb = before[k] - state.refs[k + ".target"]; ea = after[k] - state.refs[k + ".target"]
        if basis is not None:
            eb = basis.T @ eb; ea = basis.T @ ea
        return np.linalg.norm(ea) / max(np.linalg.norm(eb), 1e-30)
    tq, tc = red("tune"), red("chroma")
    ta_full, ta_reach = red("adts"), red("adts", Usig)
    print(f"\nостаток после/до: tune {tq:.1e}, chroma {tc:.1e}  (машинная точность — контуры квад/секступолей)")
    print(f"ADTS: достижимая компонента {ta_reach:.2f} (её октуполи давят), полный вектор {ta_full:.2f} "
          f"(остаток — в НЕуправляемых {3-dof} направлениях: нехватка независимых октуполей)")
    ok = tq < 0.05 and tc < 0.05 and ta_reach < 0.5
    print(f"\n[{'ok' if ok else 'ВНИМАНИЕ'}] квадрупольный и секступольный контуры свели оптику к уставке; "
          f"октупольный контур давит ADTS в пределах авторитета решётки ({dof}/3 DOF)")


if __name__ == "__main__":
    main()
