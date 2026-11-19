# -*- coding: utf-8 -*-
"""Зонд: снимает ли контур ВЧ-ЧАСТОТЫ горизонтальный 6D-радиационный пол дипольной коррекции (~3.8 мкм)?

Продолжение control/sr_floor_origin.py (то же искажение: feed-down квадруполей 0.3 мкм, seed 11; та же
6D+СИ+ВЧ модель, фаза 165°; 10→20 итераций коррекции на 6D-замкнутой орбите — детерминированно, без шума).
Там показано: пол физичен (идеальная 6D-матрица оставляет 3.82), проекция остатка на η_x — лишь ~1.4 мкм.
Здесь к дипольной коррекции добавляем ЧАСТОТУ ВЧ (Δf, Гц) — актуатор rf.freq:

  (A) FOFB, H_4D, rcond_x 3e-3            — эталон пола 3.84 (как sr_floor_origin)
  (B) FOFB, H_6D                          — идеальная дипольная матрица
  (C) FOFB H_4D + контур ВЧ (закон RfFrequencyLoop: δ̂=(η·x)/(η·η) → 0, Δf ← Δf − g·δ̂/(dδ̂/dΔf))
  (D) FOFB H_6D + контур ВЧ
  (E) совместный МНК [H_6D | ∂x/∂Δf] — лучшее, что дают диполи + частота вместе
  (A'), (C') — то же, что A, C, но с rcond_x = 1e-2, как у живого FOFB в fb_loops (опыты 1/5/6)

Отдельно печатаем СИ-орбиту НЕИСКАЖЁННОГО кольца (без коррекции): какой пол есть у идеальной машины.

ЗАПУСК: py -3.14 control/sr_floor_rf.py
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
from control.controllers.rf_frequency import delta_estimate

LAT, NSL = 26, 8
PHASE_DEG = 165.0
DIST_AMP = 0.3e-6
DIST_SEED = 11
EPS = 1e-6
EPS_F = 5.0                       # Гц — шаг для отклика на частоту (как TbtPlant._build_rf_resp)
RCOND = {"x": 3e-3, "y": 1e-2}
RCOND_LIVE_X = 1e-2               # rcond горизонтального FOFB в fb_loops (живые опыты режима 4)
N_ITER = 20
G_RF = 0.5                        # усиление контура ВЧ на итерацию (FOFB — 1, как в sr_floor_origin)


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def main():
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    volt, harm = inf["rf"]; energy = inf["energy_mev"]

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

    def set_df(ring, df):
        if hasattr(ring, "set_rf") and ring.rf:
            ring.set_rf(freq_hz=(F0 + df) if df else 0.0)     # 0 ⇒ номинал h·c/C (бит-в-бит прежнее)

    def full_orbit(ring, thx, thy, dims, dist, df=0.0):
        kicks(ring, thx, thy, dist); set_df(ring, df)
        co = ra.closed_orbit(ring, dims=dims)
        return ring.track_elements([co])[0], co

    def build_H(ring, dims):
        b, _ = full_orbit(ring, z, z, dims, dist=False)
        bx, by = b[mon, 0], b[mon, 2]
        Hx = np.zeros((len(mon), nc)); Hy = np.zeros((len(mon), nc))
        for j in range(nc):
            e = z.copy(); e[j] = EPS
            Hx[:, j] = (full_orbit(ring, e, z, dims, dist=False)[0][mon, 0] - bx) / EPS
            Hy[:, j] = (full_orbit(ring, z, e, dims, dist=False)[0][mon, 2] - by) / EPS
        return Hx, Hy

    # --- матрицы: H_4D (чистое 4D, как живой FOFB) и H_6D (чистое 6D+СИ) ---
    Hx4, Hy4 = build_H(model.ring, dims=4)
    model.close()
    Px4, Py4 = np.linalg.pinv(Hx4, rcond=RCOND["x"]), np.linalg.pinv(Hy4, rcond=RCOND["y"])
    Px4L = np.linalg.pinv(Hx4, rcond=RCOND_LIVE_X)

    def new6(sr=True):
        return ring_cpp.CppRing(lattice=LAT, sr=sr,
                                rf={"volt_mv": volt, "harmonic": harm, "phase_deg": PHASE_DEG},
                                design_mev=energy, nsl=NSL)

    r = new6()
    F0 = float(r.rf_frequency())
    eta = np.asarray(ra.twiss(r)["etax"], float)[mon]          # модельная η_x в BPM (как датчик rf.delta)
    Hx6, Hy6 = build_H(r, dims=6)
    Px6, Py6 = np.linalg.pinv(Hx6, rcond=RCOND["x"]), np.linalg.pinv(Hy6, rcond=RCOND["y"])
    # отклик на частоту: орбита x в BPM и датчик δ̂ (центральная разность ±EPS_F)
    op, cop = full_orbit(r, z, z, 6, dist=False, df=+EPS_F)
    om, com = full_orbit(r, z, z, 6, dist=False, df=-EPS_F)
    dx_df = (op[mon, 0] - om[mon, 0]) / (2 * EPS_F)
    g_rf = (delta_estimate(op[mon, 0], eta) - delta_estimate(om[mon, 0], eta)) / (2 * EPS_F)
    ddel_df = (cop[5] - com[5]) / (2 * EPS_F)
    # СИ-орбита неискажённого кольца (без коррекции, Δf=0)
    o0, co0 = full_orbit(r, z, z, 6, dist=False)
    x0 = o0[mon, 0]
    r.apply_perturb(); set_df(r, 0.0); r.close()

    print(f"СИЛА-{LAT}, nsl {NSL}: f0={F0/1e6:.6f} МГц; dδ_CO/dΔf={ddel_df*1e6:.3f} ppm/Гц; "
          f"dδ̂/dΔf={g_rf*1e6:.3f} ppm/Гц (датчик/истина {g_rf/ddel_df:.3f})")
    print(f"неискажённое кольцо, без коррекции: x rms {rms(x0)*1e6:.3f} мкм, δ̂={delta_estimate(x0, eta)*1e6:.2f} ppm, "
          f"δ_CO={co0[5]*1e6:.2f} ppm; η-часть {abs(delta_estimate(x0, eta))*rms(eta)*1e6:.3f} мкм\n")

    def correct(dims, Px, Py, rf_loop=False, joint=None):
        rr = new6()
        thx, thy, df = z.copy(), z.copy(), 0.0
        for _ in range(N_ITER):
            o, _ = full_orbit(rr, thx, thy, dims, dist=True, df=df)
            ox, oy = o[mon, 0], o[mon, 2]
            if joint is not None:                              # совместный МНК: [кики | Δf]
                u = joint @ ox
                thx = thx - u[:nc]; df = df - u[nc]
            else:
                thx = thx - Px @ ox
                if rf_loop:
                    df = df - G_RF * delta_estimate(ox, eta) / g_rf
            thy = thy - Py @ oy
        o, co = full_orbit(rr, thx, thy, dims, dist=True, df=df)
        ox = o[mon, 0]
        res = dict(x=rms(ox), y=rms(o[mon, 2]), df=df, dhat=delta_estimate(ox, eta), d_co=co[5],
                   eta_part=abs(delta_estimate(ox, eta)) * rms(eta))
        rr.apply_perturb(); set_df(rr, 0.0); rr.close()
        return res

    # совместная псевдообратная: столбцы нормированы (кики ~м/рад, частота ~м/Гц)
    A = np.hstack([Hx6, dx_df[:, None]])
    sc = np.linalg.norm(A, axis=0); sc[sc == 0] = 1.0
    Pj = (np.linalg.pinv(A / sc, rcond=RCOND["x"]).T / sc).T

    cases = [("(A) FOFB H_4D, rcond_x 3e-3", dict(Px=Px4, Py=Py4)),
             ("(B) FOFB H_6D (идеальн. дипольн.)", dict(Px=Px6, Py=Py6)),
             ("(C) FOFB H_4D + контур ВЧ", dict(Px=Px4, Py=Py4, rf_loop=True)),
             ("(D) FOFB H_6D + контур ВЧ", dict(Px=Px6, Py=Py6, rf_loop=True)),
             ("(E) совместный МНК [H_6D | Δf]", dict(Px=Px6, Py=Py6, joint=Pj)),
             ("(A') FOFB H_4D, rcond_x 1e-2", dict(Px=Px4L, Py=Py4)),
             ("(C') FOFB H_4D, rcond_x 1e-2 + ВЧ", dict(Px=Px4L, Py=Py4, rf_loop=True))]
    print(f"{'конфигурация':36} | {'rms X':>6} | {'rms Y':>6} | {'Δf, Гц':>7} | {'δ̂, ppm':>8} | {'δ_CO,ppm':>8} | {'η-часть':>7}")
    out = {}
    for lab, kw in cases:
        q = correct(6, **kw); out[lab.split()[0]] = q
        print(f"{lab:36} | {q['x']*1e6:6.3f} | {q['y']*1e6:6.3f} | {q['df']:7.3f} | {q['dhat']*1e6:8.3f} | "
              f"{q['d_co']*1e6:8.3f} | {q['eta_part']*1e6:7.3f}", flush=True)

    a, c, e = out["(A)"], out["(C)"], out["(E)"]
    print("\n[разбор]")
    print(f"  • контур ВЧ к FOFB (rcond_x 3e-3, как sr_floor_origin): X {a['x']*1e6:.2f} → {c['x']*1e6:.2f} мкм "
          f"(−{(1 - c['x']/a['x'])*100:.0f}%), Δf = {c['df']:.2f} Гц; η-часть {a['eta_part']*1e6:.2f} → {c['eta_part']*1e6:.2f}.")
    print(f"  • предел «диполи + частота» (совместный МНК): {e['x']*1e6:.2f} мкм; остаток ⊥ (корректоры ⊕ η) — "
          f"пилообразная орбита энергопотерь между резонаторами, ни диполям, ни частоте не доступна.")
    a2, c2 = out["(A')"], out["(C')"]
    print(f"  • с rcond_x живого FOFB (1e-2): X {a2['x']*1e6:.2f} → {c2['x']*1e6:.2f} мкм с контуром ВЧ "
          f"(Δf = {c2['df']:.2f} Гц).")


if __name__ == "__main__":
    main()
