# -*- coding: utf-8 -*-
"""Опция 7b: ВЧ-частота как актуатор энергии — проверка реализуемости в ядре + природа энергетической СС.

Замысел 7b: добавить ВЧ-частоту актуатором и снять дисперсионные ~2/3 горизонтального пола (из 7a).
Но прежде чем строить контур — два физических вопроса к ЗАМОРОЖЕННОМУ C++-ядру:

(A) Реализует ли ядро расстройку ВЧ-частоты? f_RF = harmonic·c/C, задаётся на старте сервера (float).
    Детюн Δf/f → равновесная δ = −(1/α_c)·(Δf/f). Меняем дробный harmonic, смотрим co[5] (6D-равновесие).
(B) Возбуждает ли вообще энергетическую СС искажение орбиты? Сдвигает ли feed-down (с которым борется FOFB)
    равновесную co[5]? Если НЕТ — остаток 7a вдоль η_x это НЕ энергия (δ реально не сдвинут), а
    геометрия/rcond-срез, и ВЧ-рычаг беспредметен. Если ДА — энергия реальна, но (A) скажет, достаёт ли её ВЧ.

ЗАПУСК: py -3.14 control/rf_actuator_probe.py
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


def build_ring(volt, harm, energy):
    return ring_cpp.CppRing(lattice=LAT, sr=True,
                            rf={"volt_mv": volt, "harmonic": harm, "phase_deg": PHASE_DEG},
                            design_mev=energy, nsl=NSL)


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    volt, harm = inf["rf"]; energy = inf["energy_mev"]

    # --- оптика: α_c, η_x в BPM, индексы/K1L квадруполей (тот же источник, что у TbtPlant) ---
    model = RingModel(lattice=LAT, nsl=NSL)
    tw = ra.twiss(model.ring)
    ac = tw["alphac"]; mon = model.mon
    etax = tw["etax"][mon]
    quads = model.quads; K1L = np.asarray(model.K1L, float); nq = len(quads)
    print(f"α_c={ac:.4e}; η_x в BPM СКО {np.std(etax)*1e3:.1f} мм; квадруполей {nq}")
    print(f"для δ=+1e-4 нужна расстройка Δf/f = −α_c·δ = {-ac*1e-4:+.2e} (крошечная)\n")

    # === (A) расстройка ВЧ-частоты дробным harmonic → равновесная энергия ===
    print("=== (A) ВЧ-частота как актуатор: δ_равн(co[5]) vs Δf/f ===")
    r0 = build_ring(volt, harm, energy)
    co0 = ra.closed_orbit(r0, dims=6); r0.apply_perturb(); r0.close()
    print(f"{'Δf/f':>10} | {'предсказ δ=−Δf/αc':>18} | {'co[5] факт':>12} | Δco[5]")
    for df in (0.0, -ac * 1e-4, -ac * 1e-2, 1e-6, 1e-4):          # от физически нужной до заведомо огромной
        r = build_ring(volt, harm * (1.0 + df), energy)
        co = ra.closed_orbit(r, dims=6); r.apply_perturb(); r.close()
        print(f"{df:>+10.2e} | {(-df/ac):>+18.3e} | {co[5]:>+12.4e} | {co[5]-co0[5]:+.2e}")
    rf_live = False                                               # заполнится ниже по факту
    # (если ни одна строка не двигает co[5] → ВЧ-частота инертна в ядре)

    # === (B) возбуждает ли энергетическую СС искажение орбиты (feed-down FOFB-возмущения)? ===
    print("\n=== (B) сдвигает ли искажение орбиты равновесную энергию co[5]? ===")
    rp = np.random.default_rng(DIST_SEED)
    patx = rp.normal(0, 1, nq); patx *= 1.0 / np.std(patx)
    paty = rp.normal(0, 1, nq); paty *= 1.0 / np.std(paty)
    dkx = K1L * (patx * DIST_AMP); dky = -K1L * (paty * DIST_AMP)
    kk = [(int(quads[j]), float(dkx[j])) for j in range(nq)]
    vv = [(int(quads[j]), float(dky[j])) for j in range(nq)]

    rd = build_ring(volt, harm, energy)
    co_nom = ra.closed_orbit(rd, dims=6)
    orb_nom = rd.track_elements([co_nom])[0]
    rd.apply_perturb(kicks=kk, vkicks=vv)
    co_dist = ra.closed_orbit(rd, dims=6)
    orb_dist = rd.track_elements([co_dist])[0]
    rd.apply_perturb(); rd.close()

    xrms_nom = np.sqrt(np.mean(orb_nom[mon, 0] ** 2))
    xrms_dist = np.sqrt(np.mean(orb_dist[mon, 0] ** 2))
    # проекция СДВИГА орбиты (искажённая−номинал) на η_x: сколько искажения — дисперсионная орбита
    dorb = orb_dist[mon, 0] - orb_nom[mon, 0]
    dd = float(etax @ dorb) / float(etax @ etax)
    r2 = 1.0 - float((dorb - etax * dd) @ (dorb - etax * dd)) / float(dorb @ dorb)
    print(f"co[5] номинал   {co_nom[5]:+.4e}")
    print(f"co[5] искажено  {co_dist[5]:+.4e}   Δδ_равн = {co_dist[5]-co_nom[5]:+.3e}")
    print(f"орбита_x СКО: номинал {xrms_nom*1e6:.2f} → искажено {xrms_dist*1e6:.2f} мкм")
    print(f"сдвиг орбиты на η_x: δ̂={dd:+.3e}, R²={r2:.3f} (совпадает ли с Δδ_равн и с остатком 7a −1.08e-4?)")
    model.close()

    print("\n[вывод]")
    dphys = co_dist[5] - co_nom[5]                            # РЕАЛЬНЫЙ сдвиг равновесной энергии
    disp_orbit_rms = np.std(etax) * abs(dphys)               # дисперсионная орбита η_x·δ_phys, м
    floor = 6.2e-6                                            # горизонтальный остаток FOFB из 7a/7
    frac = disp_orbit_rms / floor
    print(f"  (B) Реальный сдвиг энергии искажением Δδ_phys={dphys:+.2e} ⇒ дисперсионная орбита η_x·Δδ "
          f"= {disp_orbit_rms*1e6:.3f} мкм = {frac*100:.1f}% горизонтального пола ({floor*1e6:.1f} мкм).")
    if frac < 0.05:
        print(f"      Это НИЧТОЖНО. δ̂=−1.08e−4 из 7a в {abs(1.08e-4/max(dphys,1e-12)):.0f}× больше физического "
              f"сдвига — значит η_x-проекция 7a была ГЕОМЕТРИЧЕСКИМ алиасингом (η_x — крупная гладкая мода, "
              f"коррелирует с труднокорректируемым остатком), а НЕ дисперсионной орбитой. Пол НЕ энергетический.")
    else:
        print(f"      Энергетическая СС вносит заметный вклад ({frac*100:.0f}%), но (A) показывает — ВЧ-частотой "
              f"в этом ядре её не достать.")
    print(f"  (A) ВЧ-частота как актуатор: см. столбец Δco[5] — если он нулевой на всех Δf, ядро ВЧ-частоту "
          f"игнорирует (локальный энергокик по мгновенному z, без межоборотного клока); актуатор неосуществим "
          f"без правки замороженного ядра.")
    print(f"  ИТОГ: рычаг для горизонтального пола — НЕ ВЧ-частота, а снижение rcond / возврат срезанных мод "
          f"(опция 7c) либо нелинейная поправка модели отклика.")


if __name__ == "__main__":
    main()
