# -*- coding: utf-8 -*-
"""Анализ ПОВОРОТНЫХ (bend) экранов — волноотводы: веер поворотного СИ на экране.

Для bend-экрана (замороженного) и его диполя-источника:
  1) АВТОПОИСК: какие диполи реально светят на экран (касательная дуги смотрит на экран в пределах веера);
  2) МОЩНОСТЬ: проекция веера (гориз. = развёртка дуги, верт. ~1/γ) → карта Вт/мм² (полоса);
  3) СПЕКТР в центре экрана: СПЛОШНОЙ спектр поворотного СИ (ф/с/мрад²/0.1%BW vs энергия) + ε_c.
Геометрия как у ондуляторных экранов: экран неподвижен, источник берётся по ФАКТ. геометрии (--ws) ⇒
перекос диполя двигает/наклоняет веер по экрану. Несколько диполей — по ИНТЕНСИВНОСТИ (некогерентны).

ЗАПУСК
  py -3.14 undulators/bend_screen.py --lattice 27 --screen SCR_BM08
  py -3.14 undulators/bend_screen.py --lattice 27 --screen SCR_BM08 --ws <папка>   # увод при перекосе
"""
import os, sys, csv, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt
import bend_sr as bsr
from build_screens import _rot, _unit, bend_source, WORLD_UP
from screen_analysis import load_screens, RESDIR

TB = os.path.join("workspace", "tables")


def load_lattice(lattice, ws):
    path = os.path.join(ws, "tables", "lattice_configs.csv") if ws else os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.reader(open(path, encoding="utf-8-sig"), delimiter=";"))
    ci = {c: i for i, c in enumerate(rows[0])}
    byid = {r[ci["instance_id"]]: r for r in rows[1:] if r and r[ci["lattice_id"]] == str(lattice)}
    return ci, byid


def dipole_arc(r, ci, E):
    """Геометрия диполя-строки: центр, рамка, θ0, ρ, B, ε_c, P_tot. None если не изгибает."""
    if r[ci["kind"]].strip() not in ("dipole", "dq"):
        return None
    L = float(r[ci["eff_length"]] or 0); matB = float(r[ci["mat_B"]] or 0)
    th = mt.QP * matB * L
    if abs(th) < 1e-12:
        return None
    C = np.array([float(r[ci[a]] or 0) for a in ("x", "y", "z")])
    axis, uhat, vhat = (np.asarray(v, float) for v in
                        mt.frame(*(float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))))
    rho, B, ec = bsr.bend_geom(th, L, E)
    return dict(C=C, axis=axis, uhat=uhat, vhat=vhat, theta0=th, rho=rho, B=B, ec=ec,
                Ptot=bsr.total_power_kW(B, L, E, 0.2) * 1e3, L=L, dev=r[ci["device_type_id"]],
                iid=r[ci["instance_id"]])


def illuminates(d, Cs, gamma, nfan=8.0):
    """Светит ли диполь d на точку экрана Cs: касательная дуги достаёт в пределах веера (верт.~nfan/γ)."""
    dC = Cs - d["C"]; D = np.linalg.norm(dC)
    if D == 0:
        return False
    dh = dC / D
    if float(np.dot(dh, d["axis"])) <= 0:                      # экран позади
        return False
    psi = math.asin(max(-1.0, min(1.0, float(np.dot(dh, d["vhat"])))))   # верт. угол к плоскости орбиты
    if abs(psi) > nfan / gamma:
        return False
    dh_h = _unit(dh - np.dot(dh, d["vhat"]) * d["vhat"])        # горизонт. проекция направления
    Dh = math.atan2(float(np.dot(np.cross(d["axis"], dh_h), d["vhat"])), float(np.dot(d["axis"], dh_h)))
    return -abs(d["theta0"]) / 2 - 1e-3 <= Dh <= abs(d["theta0"]) / 2 + 1e-3


def footprint(scr, emitters, gamma, I, mesh):
    """Карта плотности мощности [Вт/м²] на экране: веера всех осветителей (некогерентно)."""
    hw, hh = scr["w"] / 2e3, scr["h"] / 2e3
    su = np.linspace(-hw, hw, mesh); sv = np.linspace(-hh, hh, mesh)
    UU, VV = np.meshgrid(su, sv, indexing="ij")
    cell = (2 * hw / (mesh - 1)) * (2 * hh / (mesh - 1))
    pts = (scr["C"][None, None, :] + UU[..., None] * scr["u"][None, None, :]
           + VV[..., None] * scr["v"][None, None, :])
    vnorm = bsr.vprofile_norm(gamma)
    S = np.zeros((mesh, mesh))
    for em in emitters:
        P, t, vn = em["P"], em["t"], em["vn"]                  # vn — нормаль плоскости орбиты (наклоняется с диполем)
        dh_lo, dh_hi = em["dh_lo"], em["dh_hi"]
        Ptot = em["Ptot"] * (I / 0.2)                          # масштаб тока
        d = pts - P[None, None, :]
        D = np.linalg.norm(d, axis=2); dh = d / D[..., None]
        psi = np.arcsin(np.clip((dh * vn).sum(2), -1, 1))      # верт. угол ОТ ПЛОСКОСТИ ОРБИТЫ источника
        dh_h = dh - (dh * vn).sum(2)[..., None] * vn
        dh_h /= np.linalg.norm(dh_h, axis=2, keepdims=True) + 1e-30
        Dh = np.arctan2((np.cross(t, dh_h) * vn).sum(2), (dh_h * t).sum(2))
        g = bsr.power_vprofile(psi, gamma) / vnorm             # ∫dψ=1
        inside = (Dh >= dh_lo) & (Dh <= dh_hi) & ((dh * t).sum(2) > 0)
        dPdO = np.where(inside, (Ptot / abs(em["theta0"])) * g, 0.0)   # Вт/рад²
        cos_inc = np.abs((dh * scr["n"][None, None, :]).sum(2))
        S += dPdO * cos_inc / D ** 2
    tot = float(S.sum() * cell); wsum = float(S.sum()) or 1.0
    cv = float((S * VV).sum() / wsum) * 1e3                    # верт. центроид полосы, мм
    cu = float((S * UU).sum() / wsum) * 1e3
    return dict(S=S * 1e-6, extent=(-hw * 1e3, hw * 1e3, -hh * 1e3, hh * 1e3),
                peak=float(S.max() * 1e-6), total=tot, cv=cv, cu=cu)


def emitter_geom(d, Cs, gamma):
    """Из диполя d и центра экрана: точка отбора (там, где касательная смотрит на экран) + гор. развёртка веера."""
    dC = _unit(Cs - d["C"])
    dh_h = _unit(dC - np.dot(dC, d["vhat"]) * d["vhat"])
    alpha_ex = math.atan2(float(np.dot(np.cross(d["axis"], dh_h), d["vhat"])), float(np.dot(d["axis"], dh_h)))
    alpha_ex = max(-abs(d["theta0"]) / 2, min(abs(d["theta0"]) / 2, alpha_ex))
    O = d["C"] + d["rho"] * d["uhat"]
    P = O + _rot(d["C"] - O, d["vhat"], alpha_ex)
    t = _unit(_rot(d["axis"], d["vhat"], alpha_ex))
    dh1 = (-d["theta0"] / 2) - alpha_ex; dh2 = (d["theta0"] / 2) - alpha_ex
    return dict(P=P, t=t, vn=d["vhat"], theta0=d["theta0"], Ptot=d["Ptot"], dev=d["dev"], iid=d["iid"],
                ec=d["ec"], dh_lo=min(dh1, dh2), dh_hi=max(dh1, dh2))


def plot_footprint(scr, fp, emitters, path, ec):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5.6, 3.2))
    im = ax.imshow(fp["S"].T, origin="lower", extent=fp["extent"], aspect="auto", cmap="inferno")
    ax.set_xlabel("u (горизонт.), мм"); ax.set_ylabel("v (верт.), мм")
    dev = emitters[0]["dev"] if emitters else "?"
    ax.set_title(f"{scr['sid']} ← диполь {dev} (ε_c={ec:.2f} кэВ): веер поворотного СИ\n"
                 f"пик {fp['peak']:.2f} Вт/мм², перехвачено {fp['total']:.1f} Вт", fontsize=8.5)
    cb = fig.colorbar(im, ax=ax, shrink=0.9); cb.set_label("Вт/мм²", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def plot_spectrum(Es, Fa, ec, scr, dev, path, psi_c):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.loglog(Es, np.maximum(Fa, 1), color="#1b7f5a", lw=1.6)
    ax.axvline(ec, color="#c0392b", ls="--", lw=1); ax.annotate(f"ε_c={ec:.2f} кэВ", (ec, Fa.max()),
               fontsize=8, color="#c0392b", ha="left", va="top")
    ax.set_xlabel("энергия фотона, кэВ"); ax.set_ylabel("ф/с/мрад²/0.1%BW")
    ax.set_title(f"{scr['sid']} ← диполь {dev}: спектр поворотного СИ в центре экрана "
                 f"(ψ={psi_c*1e6:+.0f} мкрад)", fontsize=8.5)
    ax.grid(alpha=0.25, which="both")
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="анализ bend-экранов (волноотводы): веер + спектр поворотного СИ")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--screen", required=True)
    ap.add_argument("--ws", default=None)
    ap.add_argument("--E", type=float, default=6.0); ap.add_argument("--I", type=float, default=0.2)
    ap.add_argument("--mesh", type=int, default=161)
    a = ap.parse_args()

    scr = next((s for s in load_screens(a.lattice) if s["sid"] == a.screen), None)
    if scr is None:
        raise SystemExit(f"экран {a.screen} не найден в screens.csv реш.{a.lattice}")
    ci, byid = load_lattice(a.lattice, a.ws)
    gamma = bsr.gamma(a.E)

    # привязанный источник (для ε_c и подписи)
    bound = bend_source(scr["src"], None, ci, byid, a.E)

    # АВТОПОИСК: какие диполи светят на центр экрана
    emitters = []
    for iid, r in byid.items():
        d = dipole_arc(r, ci, a.E)
        if d is None or not illuminates(d, scr["C"], gamma):
            continue
        emitters.append(emitter_geom(d, scr["C"], gamma))
    tag = "номинал" if not a.ws else os.path.basename(a.ws.rstrip("/\\"))
    print(f"решётка {a.lattice} [{tag}] bend-экран {a.screen} ← диполь {bound['dev']} "
          f"(θ={bound['theta0']*1e3:+.2f} мрад, ε_c={bound['ec']:.2f} кэВ, отбор {bound['phi']*1e3:.2f} мрад)")
    print(f"  осветителей найдено: {len(emitters)} ({','.join(e['dev']+'#'+e['iid'] for e in emitters[:5])}"
          f"{'…' if len(emitters)>5 else ''})")
    if not emitters:
        raise SystemExit("нет осветителей — проверьте геометрию экрана")

    os.makedirs(RESDIR, exist_ok=True); sfx = "" if not a.ws else "_pert"
    fp = footprint(scr, emitters, gamma, a.I, a.mesh)
    fpng = os.path.join(RESDIR, f"L{a.lattice}_{a.screen}_bendfoot{sfx}.png")
    plot_footprint(scr, fp, emitters, fpng, bound["ec"])

    # спектр в центре экрана (верт. угол ψ_c относительно плоскости орбиты источника)
    dC = _unit(scr["C"] - bound["P"]); psi_c = math.asin(max(-1.0, min(1.0, float(np.dot(dC, bound["vhat"])))))
    Es = np.geomspace(0.03 * bound["ec"], 12 * bound["ec"], 400)
    Fa = bsr.flux_angular(Es, psi_c, bound["ec"], a.E, a.I, gamma)
    spng = os.path.join(RESDIR, f"L{a.lattice}_{a.screen}_bendspec{sfx}.png")
    plot_spectrum(Es, Fa, bound["ec"], scr, bound["dev"], spng, psi_c)
    Fh = bsr.flux_horiz(np.array([bound["ec"]]), bound["ec"], a.E, a.I)[0]

    print(f"  веер: пик {fp['peak']:.2f} Вт/мм², перехвачено {fp['total']:.1f} Вт, "
          f"центроид (u={fp['cu']:+.3f}, v={fp['cv']:+.3f}) мм")
    print(f"  спектр в центре (ψ={psi_c*1e6:+.0f} мкрад): пик {Fa.max():.3e} ф/с/мрад²/0.1%; "
          f"поток/гор при ε_c {Fh:.3e} ф/с/мрад/0.1%")
    print(f"\n[OK] {fpng}\n[OK] {spng}")


if __name__ == "__main__":
    main()
