# -*- coding: utf-8 -*-
"""ЭТАП #2 — МУЛЬТИИСТОЧНИКОВЫЙ спектр на экране: когерентная группировка.

Один экран может ловить свет НЕСКОЛЬКИХ источников сразу (напр. ондулятор + веер соседнего диполя: вход
диполя за прямолинейкой смотрит вдоль оси ИД, его СИ попадает на тот же экран). Правило сложения:
  • ВНУТРИ устройства — амплитуды (интерференция → гармоники): это уже делает интеграл поля;
  • МЕЖДУ устройствами (разные ондуляторы, диполи) — ВЗАИМНО НЕКОГЕРЕНТНЫ ⇒ складываем ИНТЕНСИВНОСТИ.
Автопоиск находит ВСЕ источники (ондуляторы по конусу, диполи по вееру); спектр в центре экрана =
Σ|A_ИД|² (каждый — свой когерентный интеграл) + Σ F_диполь (сплошной). Единицы общие: ф/с/мрад²/0.1%BW.

ЗАПУСК: py -3.14 undulators/screen_multi.py --lattice 27 --screen SCR_U01 --gap 5
"""
import os, sys, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import undulator_sr as sr
import bend_sr as bsr
from undulator_field import trajectory, calibrate, apply_phase_errors
from screen_spectrum import center_spectrum, beam_params
from screen_analysis import load_screens, load_groups, emitter, power_params, _unit, RESDIR
from bend_screen import load_lattice, dipole_arc, illuminates
from build_screens import bend_source


def find_undulators(screen, groups, ci, byid, sig, cut=6.0):
    """Ондуляторы, чей конус достаёт центр экрана."""
    found = []
    for gid, insts in groups.items():
        P, axis = emitter(insts, ci, byid)
        if P is None:
            continue
        vec = screen["C"] - P; D = float(np.linalg.norm(vec)); dh = vec / D
        if float(np.dot(dh, axis)) <= 0:
            continue
        theta = math.acos(max(-1.0, min(1.0, float(np.dot(dh, axis)))))
        half = 0.5 * math.hypot(screen["w"], screen["h"]) / 1e3 / D
        if theta <= cut * sig + half:
            found.append((gid, P, axis))
    return found


def find_dipoles(screen, ci, byid, gamma, E, exclude):
    """Диполи, чей веер достаёт центр экрана. exclude — инстансы полюсов ондуляторов (они КОГЕРЕНТНЫ,
    их нельзя считать отдельными некогерентными источниками)."""
    found = []
    for iid, r in byid.items():
        if iid in exclude:
            continue
        d = dipole_arc(r, ci, E)
        if d is None or not illuminates(d, screen["C"], gamma):
            continue
        found.append((iid, d))
    return found


def main():
    ap = argparse.ArgumentParser(description="мультиисточниковый спектр на экране (#2: когерентная группировка)")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--screen", required=True)
    ap.add_argument("--ws", default=None)
    ap.add_argument("--gap", type=float, default=5.0)
    ap.add_argument("--E", type=float, default=6.0); ap.add_argument("--I", type=float, default=0.2)
    ap.add_argument("--beta-x", type=float, default=None); ap.add_argument("--beta-y", type=float, default=None)
    ap.add_argument("--kappa", type=float, default=sr.P.kappa); ap.add_argument("--sigE", type=float, default=sr.P.sig_E)
    ap.add_argument("--phase-error-deg", type=float, default=0.0)
    ap.add_argument("--ppp", type=int, default=64)
    ap.add_argument("--no-beam", action="store_true")
    a = ap.parse_args()

    scr = next((s for s in load_screens(a.lattice) if s["sid"] == a.screen), None)
    if scr is None:
        raise SystemExit(f"экран {a.screen} не найден")
    groups = load_groups(a.lattice); ci, byid = load_lattice(a.lattice, a.ws)
    g = bsr.gamma(a.E)
    B0 = sr.B0_of_gap(a.gap); K = sr.K_of_B0(B0)
    _, _, _, _, sig = power_params(a.gap)

    unds = find_undulators(scr, groups, ci, byid, sig)
    pole_insts = set().union(*groups.values()) if groups else set()   # полюса ИД — когерентны, не считать диполями
    dips = find_dipoles(scr, ci, byid, g, a.E, pole_insts)
    tag = "номинал" if not a.ws else os.path.basename(a.ws.rstrip("/\\"))
    print(f"решётка {a.lattice} [{tag}] экран {a.screen}: источников — ондуляторов {len(unds)}, диполей {len(dips)}")
    for gid, _, _ in unds:
        print(f"    ондулятор {gid} (когерентно внутри)")
    for iid, d in dips:
        print(f"    диполь {d['dev']}#{iid} θ={d['theta0']*1e3:+.2f} мрад ε_c={d['ec']:.2f} кэВ (некогерентно)")

    # общая энергосетка
    E1 = sr.harmonic_energy_keV(1, K, 0.0)
    Es = np.linspace(0.4 * E1, sr.harmonic_energy_keV(9, K, 0.0) * 1.05, 3200)

    # пучок (реальные β)
    etaxp = 0.0
    if (a.beta_x is None or a.beta_y is None) and not a.no_beam and unds:
        from lattice_optics import undulator_optics
        o = undulator_optics(a.lattice, unds[0][0])
        if a.beta_x is None: a.beta_x = o["betax"]
        if a.beta_y is None: a.beta_y = o["betay"]
        etaxp = o["etaxp"]
    beam = None if a.no_beam else beam_params(sr.P, a.beta_x or 5.0, a.beta_y or 2.5, a.kappa, a.sigE, etaxp)

    tr = trajectory(K, sr.P.lam_u, sr.P.N, g, a.ppp); scale, _ = calibrate(tr)
    if a.phase_error_deg:
        tr = apply_phase_errors(tr, a.phase_error_deg)

    contribs = []                                              # (label, spectrum)
    total = np.zeros_like(Es)
    # ── ондуляторы: |A|² в направлении на центр экрана (внутри — когерентно) ──
    for gid, P, axis in unds:
        uhat = _unit(np.cross([0, 0, 1.0], axis)); vhat = _unit(np.cross(axis, uhat))
        dc = _unit(scr["C"] - P); thx0 = float(np.dot(dc, uhat)); thy0 = float(np.dot(dc, vhat))
        Ise, Ibeam = center_spectrum(tr, scale, thx0, thy0, Es, beam)
        S = Ibeam if Ibeam is not None else Ise
        contribs.append((f"ондулятор {gid}", S, "#185fa5")); total = total + S      # ИНТЕНСИВНОСТИ
    # ── диполи: сплошной спектр при верт. угле ψ центра экрана (между источниками — некогерентно) ──
    for iid, d in dips:
        dc = _unit(scr["C"] - d["C"])
        psi = math.asin(max(-1.0, min(1.0, float(np.dot(dc, d["vhat"])))))
        Fa = bsr.flux_angular(Es, psi, d["ec"], a.E, a.I, g)
        contribs.append((f"диполь {d['dev']}#{iid} (фон)", Fa, "#1b7f5a")); total = total + Fa

    # график
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.6, 4.0))
    for lbl, S, c in contribs:
        ax.semilogy(Es, np.maximum(S, 1), color=c, lw=1.0, alpha=0.8, label=lbl)
    ax.semilogy(Es, np.maximum(total, 1), color="#111", lw=1.6, label="СУММА (интенсивности)")
    ax.set_xlabel("энергия фотона, кэВ"); ax.set_ylabel("ф/с/мрад²/0.1%BW")
    ax.set_title(f"{a.screen}: мультиисточниковый спектр (ИД когерентно внутри, между источниками — мощности)",
                 fontsize=8.5)
    ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both")
    os.makedirs(RESDIR, exist_ok=True)
    png = os.path.join(RESDIR, f"L{a.lattice}_{a.screen}_multi{'' if not a.ws else '_pert'}.png")
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)

    # сводка: вклад на энергии 1-й гармоники
    j = int(np.argmin(np.abs(Es - E1)))
    print(f"\n  на E₁={E1:.2f} кэВ: " + ", ".join(f"{lbl}={S[j]:.2e}" for lbl, S, _ in contribs)
          + f"; сумма={total[j]:.2e} ф/с/мрад²/0.1%")
    print(f"[OK] {png}")


if __name__ == "__main__":
    main()
