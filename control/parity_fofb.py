# -*- coding: utf-8 -*-
"""Сверка (parity): FOFB через ПЛАТФОРМУ против независимой реализации закона regulator.py.

На ОДНИХ матрицах реальной решётки СИЛА-26 и ОДНОЙ последовательности возмущения, при нулевом шуме
измерения, трассы орбиты обязаны совпасть до машинной точности — это доказывает, что каркас
(ExperimentState/Scheduler/Controller) воспроизводит установленную физику петли, ничего не ломая.
Плюс реалистичный прогон с шумом: коэффициент подавления орбиты (петля ВЫКЛ→ВКЛ).

Матрицы снимаются с ТЕКУЩЕГО кольца единым построителем (RingModel) — свежие, не из протухающего npz.
И эталон, и платформа берут ОДНИ и те же свежие Hc/Hcy/Gx/Gy (общий RingModel), поэтому сверка остаётся
бит-в-бит. ЗАПУСК: py -3.14 control/parity_fofb.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

from control import ExperimentState, Limits, Scheduler
from control.ring_model import RingModel
from control.ring_plant import RingPlant
from control.controllers.integral_svd import IntegralSVD

N = 3000
GAIN = 0.5
RCOND = 1e-2
DT = 1.0 / 2000.0                                            # 2 кГц (такт 135 об, как в regulator)


def make_dist(n_quad, n, seed=7, rho=0.99, sig=2e-8):
    """Фикс. последовательность дрейфа позиций квадруполей (AR(1)), общая для обоих прогонов."""
    rng = np.random.default_rng(seed)
    dx = np.zeros((n + 1, n_quad)); dy = np.zeros((n + 1, n_quad)); a = np.sqrt(1.0 - rho ** 2)
    for k in range(1, n + 1):
        dx[k] = rho * dx[k - 1] + a * sig * rng.standard_normal(n_quad)
        dy[k] = rho * dy[k - 1] + a * sig * rng.standard_normal(n_quad)
    return dx, dy


def reference(Hc, Hcy, Gx, Gy, dx, dy, gain, rcond):
    """Независимая реализация закона regulator.py: интегральный SVD с задержкой в один такт.
    Возмущение орбиты = Gx·dx (Gx — сдвиг квадруполя→орбита, K1L уже внутри), симметрично по Y."""
    Px = np.linalg.pinv(Hc, rcond=rcond); Py = np.linalg.pinv(Hcy, rcond=rcond)
    uh = np.zeros(Hc.shape[1]); uv = np.zeros(Hcy.shape[1])
    ox = Gx @ dx[0]; oy = Gy @ dy[0]
    rx = np.empty(N); ry = np.empty(N)
    for k in range(N):
        uh = uh - gain * (Px @ ox); uv = uv - gain * (Py @ oy)
        p = k + 1
        ox = Hc @ uh + Gx @ dx[p]; oy = Hcy @ uv + Gy @ dy[p]
        rx[k] = np.sqrt(np.mean(ox ** 2)); ry[k] = np.sqrt(np.mean(oy ** 2))
    return rx, ry


def framework(model, dx, dy, gain, rcond, meas_noise=0.0, enabled=True):
    plant = RingPlant(model=model, meas_noise=meas_noise, dist_seq=(dx, dy))
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, DT, log=["ebpm.x", "ebpm.y"])
    sch.add(IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=gain, rcond=rcond),
            rate_hz=1.0 / DT, enabled=enabled)
    sch.add(IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=gain, rcond=rcond),
            rate_hz=1.0 / DT, enabled=enabled)
    log = sch.run(N)
    return np.array(log["ebpm.x"]), np.array(log["ebpm.y"])


def rms_tail(a):
    return float(np.sqrt(np.mean(np.asarray(a[len(a) // 4:]) ** 2)))


def main():
    print("=== PARITY: платформа vs regulator.py на решётке 26 (свежие матрицы с кольца) ===")
    rm = RingModel(lattice=26, nsl=8)
    print("снимаю Hc/Hcy/Gx/Gy единым построителем с текущего кольца…", flush=True)
    Hc = rm.response(rm.hk_ids, ["ebpm.x"])["ebpm.x"]         # корректор hkick → ebpm.x
    Hcy = rm.response(rm.vk_ids, ["ebpm.y"])["ebpm.y"]        # корректор vkick → ebpm.y
    Gx = rm.response(rm.qx_ids, ["ebpm.x"])["ebpm.x"]         # сдвиг квадруполя dx → ebpm.x (K1L внутри)
    Gy = rm.response(rm.qy_ids, ["ebpm.y"])["ebpm.y"]         # сдвиг квадруполя dy → ebpm.y (−K1L внутри)
    n_quad = Gx.shape[1]
    print(f"мониторов {Hc.shape[0]}, корректоров {Hc.shape[1]}, квадруполей {n_quad}; "
          f"{N} тактов @ {1/DT/1e3:.0f} кГц, усиление {GAIN}, rcond {RCOND:g}")

    dx, dy = make_dist(n_quad, N)

    # 1) точная сверка (шум измерения = 0): эталон и платформа на ОДНИХ свежих матрицах
    rx, ry = reference(Hc, Hcy, Gx, Gy, dx, dy, GAIN, RCOND)
    fx, fy = framework(rm, dx, dy, GAIN, RCOND, meas_noise=0.0, enabled=True)
    ex = float(np.max(np.abs(rx - fx))); ey = float(np.max(np.abs(ry - fy)))
    sx = max(rx.max(), 1e-30); sy = max(ry.max(), 1e-30)
    print(f"\nточная сверка (шум=0): max|Δ| X {ex:.2e} ({ex/sx:.1e} отн.), "
          f"Y {ey:.2e} ({ey/sy:.1e} отн.)")
    print("  → " + ("СОВПАЛО до машинной точности — каркас воспроизводит петлю"
                     if max(ex / sx, ey / sy) < 1e-9 else "РАСХОЖДЕНИЕ, разбираемся"))

    # 2) реалистичный прогон с шумом измерения: ВЫКЛ → ВКЛ
    ox_off, oy_off = framework(rm, dx, dy, GAIN, RCOND, meas_noise=1e-7, enabled=False)
    ox_on, oy_on = framework(rm, dx, dy, GAIN, RCOND, meas_noise=1e-7, enabled=True)
    rox, roy = rms_tail(ox_off), rms_tail(oy_off)
    rnx, rny = rms_tail(ox_on), rms_tail(oy_on)
    print(f"\nреалистично (шум BPM 0.1 мкм): орбита СКО")
    print(f"  X: ВЫКЛ {rox*1e6:.3f} → ВКЛ {rnx*1e6:.3f} мкм  ×{rox/max(rnx,1e-30):.1f}")
    print(f"  Y: ВЫКЛ {roy*1e6:.3f} → ВКЛ {rny*1e6:.3f} мкм  ×{roy/max(rny,1e-30):.1f}")
    rm.close()
    print("\n[готово] настоящий plant на СВЕЖИХ матрицах решётки 26 + FOFB-плагин через платформу")


if __name__ == "__main__":
    main()
