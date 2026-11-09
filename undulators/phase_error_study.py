# -*- coding: utf-8 -*-
"""ЭТАП #1 — деградация спектра ондулятора от ФАЗОВЫХ ОШИБОК (качество реального ИД).

Идеальный ондулятор даёт острые гармоники; реальный — несовершенные полюса ⇒ электрон приходит к периодам
не в такт ⇒ ФАЗОВАЯ ОШИБКА (RMS σ_φ, градусы) ⇒ высокие гармоники валятся как exp(−(n·σ_φ)²). Низкие почти
не страдают. σ_φ — паспортный показатель качества ИД (хороший: <5°).

Строит: (A) гребёнку спектра идеал vs σ_φ (на оси, один электрон — чтобы изолировать эффект фаз);
        (B) пик гармоники / идеал vs n для набора σ_φ (среднее по реализациям) + теоретические кривые.
Выход: workspace/results/screens/phase_error_{spectrum,harmonics}.png + phase_error_table.csv

ЗАПУСК: py -3.14 undulators/phase_error_study.py --gap 5 --seeds 24
"""
import os, sys, math, argparse, csv
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import undulator_sr as sr
from undulator_field import trajectory, on_axis_spectrum, calibrate, apply_phase_errors
from screen_analysis import RESDIR


def harm_peak(tr, scale, K, n, p=sr.P):
    En = sr.harmonic_energy_keV(n, K, 0.0, p)
    Es = np.linspace(0.96 * En, 1.04 * En, 220)
    s = on_axis_spectrum(tr, Es) * scale
    j = int(np.argmax(s))
    return float(Es[j]), float(s[j])


def main():
    ap = argparse.ArgumentParser(description="деградация гармоник от фазовых ошибок ИД")
    ap.add_argument("--gap", type=float, default=5.0)
    ap.add_argument("--ppp", type=int, default=64)
    ap.add_argument("--seeds", type=int, default=24, help="реализаций для усреднения")
    ap.add_argument("--sigmas", type=float, nargs="*", default=[3, 5, 8, 10], help="σ_φ, градусы")
    ap.add_argument("--demo-sigma", type=float, default=8.0, help="σ_φ для гребёнки-примера")
    ap.add_argument("--mode", default="iid", choices=["iid", "rw"])
    a = ap.parse_args()
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

    B0 = sr.B0_of_gap(a.gap); K = sr.K_of_B0(B0); g = sr.gamma()
    tr0 = trajectory(K, sr.P.lam_u, sr.P.N, g, a.ppp); scale, E1 = calibrate(tr0)
    ns = [1, 3, 5, 7, 9]
    base = {n: harm_peak(tr0, scale, K, n)[1] for n in ns}
    os.makedirs(RESDIR, exist_ok=True)
    print(f"ондулятор: зазор {a.gap:g} мм, K={K:.3f}, N={sr.P.N}; усреднение по {a.seeds} реализациям ({a.mode})")
    print(f"  деградация пика гармоники I_n/I_n0 (изм ≈ теор exp(−(n·σ_φ)²)):")

    # (B) таблица/кривые: усреднённая деградация vs n
    tbl = {}
    for sphi in a.sigmas:
        acc = {n: 0.0 for n in ns}
        for seed in range(a.seeds):
            trp = apply_phase_errors(tr0, sphi, seed=seed, mode=a.mode)
            for n in ns:
                acc[n] += harm_peak(trp, scale, K, n)[1] / base[n]
        tbl[sphi] = {n: acc[n] / a.seeds for n in ns}
        print(f"  σ_φ={sphi:4.1f}°: " + "  ".join(
            f"n{n}={tbl[sphi][n]:.3f}(т{math.exp(-(n*math.radians(sphi))**2):.3f})" for n in ns))

    fig, ax = plt.subplots(figsize=(5.8, 3.8))
    cols = plt.cm.viridis(np.linspace(0.1, 0.85, len(a.sigmas)))
    for sphi, c in zip(a.sigmas, cols):
        ax.plot(ns, [tbl[sphi][n] for n in ns], "o-", color=c, label=f"σ_φ={sphi:g}° (изм)")
        ax.plot(ns, [math.exp(-(n * math.radians(sphi)) ** 2) for n in ns], "--", color=c, alpha=0.5)
    ax.set_xlabel("номер гармоники n"); ax.set_ylabel("пик гармоники / идеал")
    ax.set_title(f"Деградация гармоник от фазовой ошибки (пунктир — exp(−(n·σ_φ)²))", fontsize=9)
    ax.set_xticks(ns); ax.set_ylim(0, 1.05); ax.grid(alpha=0.25); ax.legend(fontsize=8)
    hpng = os.path.join(RESDIR, "phase_error_harmonics.png")
    fig.tight_layout(); fig.savefig(hpng, dpi=120); plt.close(fig)

    # (A) гребёнка: идеал vs demo-σ_φ (лёгкое инстр. сглаживание ~0.4%, чтобы гармоники читались пиками)
    from screen_spectrum import energy_smear
    Emax = sr.harmonic_energy_keV(9, K, 0.0) * 1.05
    Es = np.linspace(0.5 * E1, Emax, 4000)
    S0 = energy_smear(Es, on_axis_spectrum(tr0, Es) * scale, 0.004)
    Sd = energy_smear(Es, on_axis_spectrum(apply_phase_errors(tr0, a.demo_sigma, seed=1, mode=a.mode), Es) * scale, 0.004)
    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.semilogy(Es, np.maximum(S0, 1), color="#9fb4c8", lw=1.0, label="идеал (σ_φ=0)")
    ax.semilogy(Es, np.maximum(Sd, 1), color="#c0392b", lw=1.3, label=f"σ_φ={a.demo_sigma:g}°")
    ax.set_xlabel("энергия фотона, кэВ"); ax.set_ylabel("ф/с/мрад²/0.1%BW")
    ax.set_title("Спектр на оси: идеальный vs с фазовой ошибкой (высокие гармоники валятся)", fontsize=9)
    ax.legend(fontsize=8); ax.grid(alpha=0.25, which="both")
    spng = os.path.join(RESDIR, "phase_error_spectrum.png")
    fig.tight_layout(); fig.savefig(spng, dpi=120); plt.close(fig)

    with open(os.path.join(RESDIR, "phase_error_table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(["sigma_phi_deg"] + [f"n{n}" for n in ns] +
                                                      [f"n{n}_theor" for n in ns])
        for sphi in a.sigmas:
            w.writerow([sphi] + [f"{tbl[sphi][n]:.4f}" for n in ns] +
                       [f"{math.exp(-(n*math.radians(sphi))**2):.4f}" for n in ns])
    print(f"\n[OK] {hpng}\n[OK] {spng}\n[OK] {os.path.join(RESDIR,'phase_error_table.csv')}")


if __name__ == "__main__":
    main()
