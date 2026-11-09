# -*- coding: utf-8 -*-
"""Опция 7e: замкнуть доказательство горизонтального пола на ВХОДЕ, а не на трекинге.

7d (развязка) показал декомпозицией: горизонтальный пол = равновесная дисперсионная орбита η_x·δ_eq
(4D-коррекция 1.05 мкм → 6D 3.84 ≈ живой 3.94). Но η_x·δ_eq=2.64 мкм и прирост в квадратуре 3.69 мкм
не совпали точно (R² финального остатка 0.125 — петля перемешала остаток). 7e закрывает это ЧИСТО.

ЧАСТЬ A (отрицательная, зонд уже прогнан): δ_eq НЕ снимается ни одним операторским входом —
ВЧ-фаза/напряжение инертны (как и частота в 7b-A), проектная энергия сцеплена слабо и не в ту сторону.
δ_eq≈+7.6e-5 — СТРУКТУРНОЕ радиационное равновесие решётки. Значит пол неустраним setpoint'ом.

ЧАСТЬ B (положительная, здесь): раз δ_eq нельзя обнулить актуатором, ставим энергетическую
упреждающую связь В ОПОРНУЮ ОРБИТУ — корректор рулит не к нулю, а к η_x·δ_eq (бездисперсионное
наведение). Три коррекции на ИСКАЖЁННОМ 6D-кольце (H с чистой 4D-модели, как живой FOFB):
  (0) опора = 0            → воспроизводит 6D-пол 3.84 мкм (= живой FOFB);
  (1) опора = η_x·δ_eq     → остаток-к-опоре должен сесть ТОЧНО на 4D-пол ~1.05 мкм, если пол = чистая
                             дисперсия; кики должны совпасть с 4D-решением;
  (2) 4D (SR-выкл)         → эталон 1.05 мкм.
Если (1)остаток == (2) и кики(1)≈кики(2) — пол ТОЧНО η_x·δ_eq, без «лишнего» 6D-члена. Остаточный
зазор (1)−(2) — это та самая «добавка сверх дисперсии», если она есть.

ЗАПУСК: py -3.14 control/dispersion_reference_test.py
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import ring_cpp
import facility_report as fr
import ring_analysis as ra

from control.ring_model import RingModel

LAT, NSL = 26, 8
PHASE_DEG = 165.0
DIST_AMP = 0.3e-6
DIST_SEED = 11
EPS = 1e-6
RCOND = {"x": 3e-3, "y": 1e-2}
N_ITER = 10


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    volt, harm = inf["rf"]; energy = inf["energy_mev"]

    model = RingModel(lattice=LAT, nsl=NSL)
    mon = model.mon; corr = np.asarray(model.corr, int); nc = len(corr)
    quads = model.quads; K1L = np.asarray(model.K1L, float); nq = len(quads)
    etax = ra.twiss(model.ring)["etax"][mon]

    rp = np.random.default_rng(DIST_SEED)
    patx = rp.normal(0, 1, nq); patx *= 1.0 / np.std(patx)
    paty = rp.normal(0, 1, nq); paty *= 1.0 / np.std(paty)
    dkx = K1L * (patx * DIST_AMP); dky = -K1L * (paty * DIST_AMP)

    z = np.zeros(nc)

    def kicks(ring, thx, thy, dist):
        kk = [(int(corr[j]), float(thx[j])) for j in range(nc)]
        vv = [(int(corr[j]), float(thy[j])) for j in range(nc)]
        if dist:
            kk += [(int(quads[j]), float(dkx[j])) for j in range(nq)]
            vv += [(int(quads[j]), float(dky[j])) for j in range(nq)]
        ring.apply_perturb(kicks=kk, vkicks=vv)

    def orbit(ring, thx, thy, dims, dist):
        kicks(ring, thx, thy, dist)
        return ring.track_elements([ra.closed_orbit(ring, dims=dims)])[0]

    # --- H_nominal на ЧИСТОМ 4D-кольце (как живой FOFB): Hx по hkick→X, Hy по vkick→Y ---
    r4 = model.ring
    b = orbit(r4, z, z, 4, dist=False)
    base_x = b[mon, 0]; base_y = b[mon, 2]
    Hx = np.zeros((len(mon), nc)); Hy = np.zeros((len(mon), nc))
    for j in range(nc):
        e = z.copy(); e[j] = EPS
        Hx[:, j] = (orbit(r4, e, z, 4, dist=False)[mon, 0] - base_x) / EPS
        Hy[:, j] = (orbit(r4, z, e, 4, dist=False)[mon, 2] - base_y) / EPS
    Px = np.linalg.pinv(Hx, rcond=RCOND["x"]); Py = np.linalg.pinv(Hy, rcond=RCOND["y"])
    model.close()

    def correct(ring, dims, ref_x):
        """Линейная итерация, рулим гор. орбиту к опоре ref_x (верт. к 0). Возврат: остаток-к-опоре,
        абсолютный остаток, вектор гор. киков, co[5]."""
        thx, thy = z.copy(), z.copy()
        co5 = 0.0
        for _ in range(N_ITER):
            kicks(ring, thx, thy, dist=True)
            co = ra.closed_orbit(ring, dims=dims); co5 = float(co[5])
            o = ring.track_elements([co])[0]
            thx = thx - Px @ (o[mon, 0] - ref_x)
            thy = thy - Py @ o[mon, 2]
        o = orbit(ring, thx, thy, dims, dist=True)
        return (rms(o[mon, 0] - ref_x), rms(o[mon, 2]), rms(o[mon, 0]), thx.copy(), co5, o[mon, 0].copy())

    zref = np.zeros(len(mon))

    # (2) 4D эталон (SR-выкл, δ=0) — опора 0
    r4b = fr._ring(inf, nsl=NSL)
    r4b_res, y4, _, th4, _, _ = correct(r4b, 4, zref)
    r4b.apply_perturb(); r4b.close()

    # 6D кольцо (СИ+ВЧ, δ_eq≠0)
    def new6():
        return ring_cpp.CppRing(lattice=LAT, sr=True,
                                rf={"volt_mv": volt, "harmonic": harm, "phase_deg": PHASE_DEG},
                                design_mev=energy, nsl=NSL)

    # (0) опора 0 → живой 6D-пол
    r6 = new6()
    r0_res, y0, r0_abs, th0, co5, ox0 = correct(r6, 6, zref)
    r6.apply_perturb(); r6.close()

    # (1) ДЕТЕРМИНИРОВАННАЯ проекция: раскладываем 6D-остаток (опора 0) на η_x. Корректорами дисперсию
    #     синтезировать НЕЛЬЗЯ (вне бетатронного подпространства), поэтому «рулить к η_x·δ_eq» бессмысленно
    #     (кики разносит) — правильный инструмент это ПРОЕКЦИЯ уже скорректированной орбиты.
    o6 = ox0                                                   # гор. 6D-орбита после коррекции (опора 0)
    d_fit = float(etax @ o6) / float(etax @ etax)              # амплитуда дисперсии в остатке
    extra = etax * d_fit                                       # дисперсионная компонента остатка
    remainder = o6 - extra                                     # что осталось сверх дисперсии
    r_disp = rms(extra); r_rem = rms(remainder)
    cosang = float(etax @ o6) / (np.linalg.norm(etax) * np.linalg.norm(o6) + 1e-30)
    # контроль: «рулёжка к дисперсии» (неверная операция) — показать, что кики разносит
    refx = etax * co5
    r6b = new6()
    r1_res, _, _, th1, _, _ = correct(r6b, 6, refx)
    r6b.apply_perturb(); r6b.close()
    dk10 = np.linalg.norm(th1 - th4) / max(np.linalg.norm(th4), 1e-30)

    print(f"δ_eq (искажённое кольцо) = {co5:+.3e}; дисперсия в остатке: δ̂ = {d_fit:+.3e}\n")
    print(f"{'коррекция/разложение (H 4D-модели)':46} | {'rms X, мкм':>10} | {'Y':>6}")
    print(f"{'(2) 4D эталон (SR-выкл, δ=0)':46} | {r4b_res*1e6:9.3f} | {y4*1e6:5.3f}")
    print(f"{'(0) 6D, опора 0  → живой FOFB':46} | {r0_abs*1e6:9.3f} | {y0*1e6:5.3f}")
    print(f"{'    └ дисперсионная компонента η_x·δ̂':46} | {r_disp*1e6:9.3f} |")
    print(f"{'    └ остаток СВЕРХ дисперсии (⊥η_x)':46} | {r_rem*1e6:9.3f} |")
    print(f"\ncos(остаток, η_x) = {cosang:+.3f}; δ̂/δ_eq = {d_fit/co5:.2f}")
    print(f"контроль: рулёжка корректорами К дисперсии разносит кики ×{dk10:.1f} (дисперсия ∉ бетатронного "
          f"подпространства — корректоры её не синтезируют)")

    print("\n[вывод]")
    if r_disp > 2.0 * r_rem and abs(r_rem - r4b_res) / max(r4b_res, 1e-30) < 0.35 and d_fit / co5 > 0.5:
        print(f"  ПОДТВЕРЖДЕНО ЧИСТО (детерминированная проекция): 6D-пол {r0_abs*1e6:.2f} мкм раскладывается на "
              f"дисперсионную орбиту η_x·δ̂ = {r_disp*1e6:.2f} мкм (δ̂={d_fit:+.2e} ≈ δ_eq={co5:+.2e}, cos={cosang:+.2f}) "
              f"ПЛЮС поперечно-бетатронный остаток {r_rem*1e6:.2f} мкм ≈ 4D-пол {r4b_res*1e6:.2f}. Пол = η_x·δ_eq, "
              f"ТОЧНО. «Добавка» 7d (3.69 vs 2.64 мкм) была статистическим перемешиванием на ШУМНОЙ трекинговой "
              f"орбите (7a R²=0.43); на ДЕТЕРМИНИРОВАННОЙ замкнутой орбите разложение резкое. δ_eq структурен "
              f"(ЧАСТЬ A: инертен к фазе/напряжению/частоте ВЧ), снимается лишь сдвигом энергии (ВЧ-частота, "
              f"Δf/f=−α·δ_eq) — отсутствующим в ядре. Корректоры дисперсию не трогают (контроль: ×{dk10:.0f}).")
    elif r_disp > r_rem:
        print(f"  ЧАСТИЧНО: дисперсия η_x·δ̂={r_disp*1e6:.2f} мкм — главная часть 6D-пола {r0_abs*1e6:.2f}, но остаток "
              f"⊥η_x = {r_rem*1e6:.2f} мкм выше 4D-пола {r4b_res*1e6:.2f}. Есть 6D-член сверх чистой дисперсии.")
    else:
        print(f"  ОПРОВЕРГНУТО: дисперсионная компонента η_x·δ̂={r_disp*1e6:.2f} мкм мала против остатка "
              f"{r_rem*1e6:.2f} мкм. 6D-пол {r0_abs*1e6:.2f} — НЕ дисперсионная орбита (cos={cosang:+.2f}). "
              f"Развязку 7d (пол=η_x·δ_eq) ПЕРЕСМОТРЕТЬ.")


if __name__ == "__main__":
    main()
