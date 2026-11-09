# -*- coding: utf-8 -*-
"""Нелинейный plant через ПЛАТФОРМУ: тот же FOFB-плагин давит реальную орбиту реального пучка.

Показывает, что ExperimentState/Scheduler/Controller (та же платформа, что гоняла линейный RingPlant и
каскад B) работают и против НЕЛИНЕЙНОЙ машины — пучок из частиц в кольце с СИ+ВЧ, ведомый track_tbt.
Статическое искажение орбиты (feed-down квадруполей) давится интегральным SVD-контуром, который внутри
обращает ЛИНЕЙНУЮ модель (Hc с RingModel), а применяется к нелинейной реальности — как на машине.

Сенсоры режима 2: орбита e⁻ (ebpm.x/y), фотон (photon.x/y), и пооборотный центроид (x,px) — фазовый
портрет. ЗАПУСК: py -3.14 control/tbt_demo.py   |   быстро: py -3.14 control/tbt_demo.py smoke
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control import ExperimentState, Limits, Scheduler
from control.tbt_plant import TbtPlant
from control.controllers.integral_svd import IntegralSVD

SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
N = 60 if SMOKE else 200                                      # базовых тактов (чанков по 15 об)
TACT_CHUNKS = 9                                               # контур обновляется раз в 9 чанков = 135 об ≈ 2 кГц
OUT = os.path.join("workspace", "results", "screens", "tbt_fofb_phase.png")


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def run(enabled, capture=False):
    plant = TbtPlant(lattice=26, nsl=8, settle_turns=(150 if SMOKE else 300), dist_amp=0.3e-6, meas_noise=1e-7)
    dt = plant.stride * plant.TREV
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, dt)
    sch.add(IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5, rcond=1e-2),
            rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    sch.add(IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=0.5, rcond=1e-2),
            rate_hz=1.0 / (dt * TACT_CHUNKS), enabled=enabled)
    ox, oy, px = [], [], []
    for _ in range(N):
        sch.step()
        ox.append(np.asarray(state.observables["ebpm.x"]))
        oy.append(np.asarray(state.observables["ebpm.y"]))
        px.append(np.asarray(state.observables["photon.x"]))
    tail = slice(len(ox) // 2, None)                          # хвост (контур сошёлся)
    # метрика DC-искажения: СКО СРЕДНЕЙ по хвосту орбиты (остаточный бетатрон/шум усредняется)
    out = (rms(np.mean(ox[tail], axis=0)), rms(np.mean(oy[tail], axis=0)), rms(np.mean(px[tail], axis=0)),
           plant.alive, sch, np.array(plant.centroid_tbt) if capture else None)
    if not capture:
        plant.close()
    else:
        out = out + (plant,)
    return out


def save_phase(centroid, path):
    """Фазовый портрет: горизонтальный центроид x–px, цвет = оборот (видно декогеренцию+затухание к орбите)."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(matplotlib недоступен: {e} — фазовый портрет пропущен)")
        return False
    x = centroid[:, 0] * 1e6; pxv = centroid[:, 1] * 1e6; t = np.arange(len(x))
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    sc = ax.scatter(x, pxv, c=t, cmap="viridis", s=6)
    ax.plot(x, pxv, lw=0.3, color="0.6", alpha=0.5)
    ax.set_xlabel("центроид x, мкм"); ax.set_ylabel("центроид px, мкрад")
    ax.set_title("Фазовый портрет центроида (FOFB вкл): декогеренция → затухание к орбите")
    fig.colorbar(sc, ax=ax, label="оборот")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    return True


def main():
    print("=== НЕЛИНЕЙНЫЙ plant через платформу: FOFB-плагин на реальном пучке (СИ+ВЧ, track_tbt) ===")
    print(f"{N} чанков×{15} об; статическое искажение feed-down 0.3 мкм; контур раз в {TACT_CHUNKS} чанков ≈ 2 кГц\n")

    oxf, oyf, phf, alf, _, _ = run(enabled=False)
    oxn, oyn, phn, aln, sch, cen, plant = run(enabled=True, capture=True)

    print(f"{'режим':16} | {'орбита X':>9} | {'орбита Y':>9} | {'фотон X':>9} | пучок")
    print(f"{'FOFB выкл':16} | {oxf*1e6:8.3f} | {oyf*1e6:8.3f} | {phf*1e6:8.3f} | {'жив' if alf else 'ПОТЕРЯН'}")
    print(f"{'FOFB вкл':16} | {oxn*1e6:8.3f} | {oyn*1e6:8.3f} | {phn*1e6:8.3f} | {'жив' if aln else 'ПОТЕРЯН'}")
    print(f"\nподавление орбиты на РЕАЛЬНОМ пучке: X ×{oxf/max(oxn,1e-30):.1f}, Y ×{oyf/max(oyn,1e-30):.1f}; "
          f"фотон ×{phf/max(phn,1e-30):.1f} (мкм СКО)")
    print("\nтаблица контуров (имя, частота Гц, период тактов, вкл):")
    for name, rate, period, en in sch.table():
        print(f"   {name:7} | {rate:8.1f} | {period:6d} | {'да' if en else 'нет'}")

    ok = save_phase(cen, OUT) if cen is not None and len(cen) else False
    plant.close()
    if ok:
        print(f"\nфазовый портрет (режим 2): {OUT}")
    print("\n[готово] тот же каркас/плагин гоняет НЕЛИНЕЙНУЮ машину — динамика настоящая, модель контура линейная")


if __name__ == "__main__":
    main()
