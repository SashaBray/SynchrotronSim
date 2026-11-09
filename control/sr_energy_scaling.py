# -*- coding: utf-8 -*-
"""Опция 7f: радиационный член горизонтального пола — реальное 6-ГэВ число или артефакт суррогатной γ?

7e показал: 6D-пол 3.84 мкм = off-momentum (~1.9) + СИНХРОТРОННОЕ ИЗЛУЧЕНИЕ (~3.1 мкм). Оговорка: SR ∝ γ⁴,
а в проекте есть суррогатная γ≈7 ([[energy-convention-surrogate-gamma]], V_BEAM в element_compare/bip_matrix,
путь СРАВНЕНИЯ диффур↔кик, SR=0). Контрольный путь строит CppRing(design_mev=6000, sr=True) напрямую,
ring_cpp.py:169 pc=design_mev·(1+δ) → реальные 6 ГэВ. ПРОВЕРКА по существу: радиационный член пола обязан
масштабироваться как U0 ∝ E⁴. Скан design_mev:
  • радиационный член ∝ E⁴ (показатель ≈4) → SR считается на РЕАЛЬНОЙ энергии кольца → 3.1 мкм — число 6 ГэВ;
  • член плоский по E → энергия SR пришпилена к суррогату → оговорка в силе.
H_4D (модель живого FOFB) фиксирована (строится раз на 6000, как на реальной машине). Искажение фикс.

ЗАПУСК: py -3.14 control/sr_energy_scaling.py
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
ENERGIES = [3000.0, 4242.6, 6000.0, 8485.3, 12000.0]   # фактор √2, декада вокруг 6 ГэВ


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    volt, harm = inf["rf"]

    model = RingModel(lattice=LAT, nsl=NSL)
    mon = model.mon; corr = np.asarray(model.corr, int); nc = len(corr)
    quads = model.quads; K1L = np.asarray(model.K1L, float); nq = len(quads)

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

    # H_4D фикс (как у живого FOFB) — на чистом 4D-кольце 6 ГэВ
    r4 = model.ring
    b = orbit(r4, z, z, 4, dist=False)
    bx, by = b[mon, 0], b[mon, 2]
    Hx = np.zeros((len(mon), nc)); Hy = np.zeros((len(mon), nc))
    for j in range(nc):
        e = z.copy(); e[j] = EPS
        Hx[:, j] = (orbit(r4, e, z, 4, dist=False)[mon, 0] - bx) / EPS
        Hy[:, j] = (orbit(r4, z, e, 4, dist=False)[mon, 2] - by) / EPS
    Px = np.linalg.pinv(Hx, rcond=RCOND["x"]); Py = np.linalg.pinv(Hy, rcond=RCOND["y"])
    model.close()

    def correct(ring, dims):
        thx, thy = z.copy(), z.copy()
        for _ in range(N_ITER):
            kicks(ring, thx, thy, dist=True)
            o = ring.track_elements([ra.closed_orbit(ring, dims=dims)])[0]
            thx = thx - Px @ o[mon, 0]
            thy = thy - Py @ o[mon, 2]
        o = orbit(ring, thx, thy, dims, dist=True)
        return rms(o[mon, 0]), rms(o[mon, 2])

    def ring6(E, sr):
        return ring_cpp.CppRing(lattice=LAT, sr=sr,
                                rf={"volt_mv": volt, "harmonic": harm, "phase_deg": PHASE_DEG},
                                design_mev=E, nsl=NSL)

    print(f"{'E, МэВ':>8} | {'U0, МэВ':>9} | {'SR-выкл X':>9} | {'SR-вкл X':>9} | {'рад.член √(on²−off²)':>20} | {'Y(on)':>6}")
    Es, rads, u0s = [], [], []
    for E in ENERGIES:
        inf_e = dict(inf); inf_e["energy_mev"] = E
        try:
            rr = ring6(E, True); u0 = float(fr.radiation(rr, inf_e).get("U0", float("nan"))); rr.close()
        except Exception:
            u0 = float("nan")
        r = ring6(E, False); xoff, _ = correct(r, 6); r.apply_perturb(); r.close()
        r = ring6(E, True);  xon, yon = correct(r, 6); r.apply_perturb(); r.close()
        rad = np.sqrt(max(xon ** 2 - xoff ** 2, 0.0))
        print(f"{E:8.0f} | {u0:9.4f} | {xoff*1e6:8.3f} | {xon*1e6:8.3f} | {rad*1e6:19.3f} | {yon*1e6:5.3f}")
        Es.append(E); rads.append(rad); u0s.append(u0)

    Es = np.array(Es); rads = np.array(rads); u0s = np.array(u0s)
    ok = rads > 1e-10
    if ok.sum() >= 2:
        p = np.polyfit(np.log(Es[ok]), np.log(rads[ok]), 1)[0]            # рад.член ∝ E^p
        pu = np.polyfit(np.log(Es[ok]), np.log(u0s[ok]), 1)[0] if np.all(u0s[ok] > 0) else float("nan")
        print(f"\nпоказатель степени: рад.член ∝ E^{p:.2f}   (U0 ∝ E^{pu:.2f}; теория U0 ∝ E⁴)")
        print("[вывод]")
        if 3.0 < p < 5.0:
            print(f"  ПОДТВЕРЖДЕНО: радиационный член пола масштабируется ∝ E^{p:.1f} ≈ E⁴ — SR считается на "
                  f"РЕАЛЬНОЙ энергии кольца (design_mev), не на суррогате. 3.1 мкм при 6 ГэВ — ФИЗИЧЕСКОЕ число "
                  f"реальной машины, оговорка снята. Контрольный путь (CppRing sr=True) живёт на 6 ГэВ; "
                  f"суррогат γ≈7 — только в element_compare/bip_matrix (путь извлечения коэффициентов, SR=0).")
        elif abs(p) < 1.0:
            print(f"  ОГОВОРКА В СИЛЕ: радиационный член плоский по E (∝ E^{p:.1f}) — энергия SR пришпилена, не "
                  f"следует design_mev. 3.1 мкм нельзя переносить на реальную машину как есть.")
        else:
            print(f"  ЧАСТИЧНО: показатель E^{p:.1f} не E⁴ и не плоско — SR зависит от E нетривиально, разобрать.")


if __name__ == "__main__":
    main()
