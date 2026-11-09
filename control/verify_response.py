# -*- coding: utf-8 -*-
"""Проверка единого построителя: матрицы, снятые RingModel, сверяем с НЕЗАВИСИМЫМ построителем
orbit_correct.response_matrix на ТОМ ЖЕ текущем кольце (это тест КОРРЕКТНОСТИ, не затронутый тем,
что кэш regulator_resp_exp26.npz мог протухнуть после правок таблиц решётки). Дополнительно печатаем
расхождение со старым npz — как меру его протухания.

ЗАПУСК: py -3.14 control/verify_response.py
"""
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import orbit_correct as oc
from control.ring_model import RingModel

# Легаси-кэш regulator.py (мера протухания; RingPlant на него больше не опирается — снимает свежее с кольца)
LEGACY_NPZ = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                          "workspace", "results", "regulator_resp_exp26.npz")


def reldiff(A, B):
    s = max(np.max(np.abs(B)), 1e-30)
    return float(np.max(np.abs(A - B)) / s)


def main():
    print("=== VERIFY: единый построитель RingModel vs orbit_correct на ТОМ ЖЕ кольце (решётка 26) ===")
    rm = RingModel(lattice=26, nsl=8)
    out = rm.response(rm.hk_ids + rm.vk_ids, ["ebpm.x", "ebpm.y"])      # из кэша (снято ранее на этом кольце)
    nh = len(rm.hk_ids)
    Hc, Hcy = out["ebpm.x"][:, :nh], out["ebpm.y"][:, nh:]

    print(f"независимый эталон orbit_correct.response_matrix (тот же шаг 1e-6) на том же кольце…", flush=True)
    r = oc.response_matrix(rm.ring, corr=rm.corr, theta=1e-6)           # тот же пробный шаг, что у построителя
    Rx, Ry = r["Rx"][rm.mon - 1], r["Ry"][rm.mon - 1]

    dx, dy = reldiff(Hc, Rx), reldiff(Hcy, Ry)
    print(f"\nHc  vs orbit_correct Rx[mon]: отн. расхождение {dx:.2e}")
    print(f"Hcy vs orbit_correct Ry[mon]: отн. расхождение {dy:.2e}")
    ok = max(dx, dy) < 1e-6
    print("→ " + ("СОВПАЛО бит-в-бит — единый построитель корректен (россыпь сводится к нему)"
                  if ok else "РАСХОЖДЕНИЕ — есть баг в построителе"))
    print("(при РАЗНОМ шаге H-расхождение ~3e-3 — нелинейность отклика, не баг)")

    ref = np.load(LEGACY_NPZ)
    print(f"\n(контекст) расхождение с кэшем regulator_resp_exp26.npz (26.10): "
          f"Hc {reldiff(Hc, ref['Hc']):.2e}, Hcy {reldiff(Hcy, ref['Hcy']):.2e} — "
          f"мера протухания кэша после правок таблиц решётки")
    rm.close()


if __name__ == "__main__":
    main()
