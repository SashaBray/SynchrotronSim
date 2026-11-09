# -*- coding: utf-8 -*-
"""ЭТАП B — ядро излучения ондулятора: КОГЕРЕНТНЫЙ интеграл поля Лиенара–Вихерта вдоль траектории.

Это та самая «интерференционная решётка»: гармоники рождаются фазовым множителем exp[iω(t−n·r/c)],
проинтегрированным по траектории сквозь ВСЕ периоды. Траектория — ИДЕАЛЬНАЯ синусоида (из K, λ_u, N, γ;
решение по выбору пользователя), планарный ондулятор ⇒ развёртка в горизонтали, поляризация линейная H.
Перекос устройства применяется ЖЁСТКО к рамке траектории снаружи (в screen_spectrum.py) — здесь физика
в СОБСТВЕННОЙ рамке устройства: z — вдоль оси, x — горизонталь (развёртка), y — вертикаль.

ФОРМУЛА (велосити-форма far-field, Jackson/Kim):
  A(ω,n) = ∫ [n×(n×β)] e^{iω(t−n·r/c)} dt,   интенсивность ∝ ω²|A|².
Для конечных N это автоматически даёт множитель решётки (sin Nπδ/sin πδ)². Абсолютная шкала калибруется
на аналитический поток на оси flux_density_onaxis(1,K) из undulator_sr (единый источник нормировки);
ПОЛОЖЕНИЯ и ОТНОШЕНИЯ гармоник — из самого интеграла (независимая проверка).

Самотест: py -3.14 undulators/undulator_field.py        # сверка энергий и отношений F_n с undulator_sr
"""
import os, sys, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
import numpy as np
import undulator_sr as sr

E_CH = 1.602176634e-19          # Кл
HBAR = 1.054571817e-34          # Дж·с
C = 299792458.0                 # м/с


def omega_of_keV(E_keV):
    return np.asarray(E_keV, float) * 1e3 * E_CH / HBAR         # рад/с


def trajectory(K, lam_u, N, gamma, ppp=64):
    """Идеальная планарная синусоида. Возврат словаря с z, x, t, beta (Ns,3), веса интегрирования w=dz/(cβz)."""
    k_u = 2 * math.pi / lam_u
    Ns = int(ppp * N) + 1
    z = np.linspace(0.0, N * lam_u, Ns)
    beta_tot2 = 1.0 - 1.0 / gamma ** 2
    bx = -(K / gamma) * np.cos(k_u * z)                          # поперечная скорость (развёртка)
    bz = np.sqrt(np.clip(beta_tot2 - bx ** 2, 1e-30, None))
    # x(z)=∫(βx/βz)dz, t(z)=∫dz/(cβz) — кумулятивная трапеция
    dz = np.diff(z)
    def cum(integrand):
        mid = 0.5 * (integrand[1:] + integrand[:-1]) * dz
        return np.concatenate([[0.0], np.cumsum(mid)])
    x = cum(bx / bz)
    t = cum(1.0 / (C * bz))
    beta = np.stack([bx, np.zeros_like(bx), bz], axis=1)         # (Ns,3)
    r = np.stack([x, np.zeros_like(x), z], axis=1)               # (Ns,3)
    # веса трапеции для ∫ ... dz/(cβz)
    wz = np.empty(Ns); wz[0] = 0.5 * dz[0]; wz[-1] = 0.5 * dz[-1]
    wz[1:-1] = 0.5 * (dz[1:] + dz[:-1])
    w = wz / (C * bz)
    return dict(z=z, r=r, beta=beta, t=t, w=w, K=K, gamma=gamma, N=N, lam_u=lam_u, L=N * lam_u)


def apply_phase_errors(tr, sigma_phi_deg, seed=0, mode="iid", p=sr.P):
    """Наложить ФАЗОВЫЕ ОШИБКИ по периодам на идеальную траекторию (качество реального ИД).

    Каждому периоду k — случайная фазовая ошибка φ_k (RMS σ_φ, градусы); в траекторию вносится сдвиг
    тайминга δt=φ_k/ω₁ (ω₁ — фундамент). Интеграл поля сам даёт расфазировку ~exp(−(n·σ_φ)²): низкие
    гармоники почти не страдают, высокие (n=5,7,9) валятся. mode: iid (независимые) | rw (случ. блуждание).
    Среднее и линейный тренд убираются (это сдвиг энергии/резонанса, не «ошибка»)."""
    if not sigma_phi_deg:
        return tr
    rng = np.random.default_rng(seed)
    sig = math.radians(sigma_phi_deg); N = tr["N"]
    phis = np.cumsum(rng.normal(0, sig, N)) if mode == "rw" else rng.normal(0, sig, N)
    k = np.arange(N); A = np.polyfit(k, phis, 1); phis = phis - (A[0] * k + A[1])   # детренд
    if mode == "rw" and phis.std() > 0:                       # нормировать RMS к заданному σ_φ
        phis *= sig / phis.std()
    om1 = omega_of_keV(sr.harmonic_energy_keV(1, tr["K"], 0.0, p))
    kper = np.clip((tr["z"] / tr["lam_u"]).astype(int), 0, N - 1)
    tr2 = dict(tr); tr2["t"] = tr["t"] + phis[kper] / om1
    tr2["sigma_phi_deg"] = sigma_phi_deg
    return tr2


def _dirs(thx, thy):
    """Единичные направления n из малых углов (thx гор., thy верт.). Формы broadcast-совместимы."""
    thx = np.asarray(thx, float); thy = np.asarray(thy, float)
    nz = np.sqrt(np.clip(1.0 - thx ** 2 - thy ** 2, 0.0, None))
    return np.stack([thx, thy, nz], axis=-1)                     # (...,3)


def amplitude(tr, E_keV, dirs, chunk=256):
    """|A|²·ω² по сетке направлений dirs(...,3) и энергий E_keV(nE). Возврат shape (*dirs_shape, nE).

    A(ω,n)=Σ_s [n(n·β)−β] e^{iω(t−n·r/c)} w_s.  Векторизовано по s, чанк по направлениям."""
    r = tr["r"]; beta = tr["beta"]; t = tr["t"]; w = tr["w"]     # (Ns,3),(Ns,3),(Ns,),(Ns,)
    om = omega_of_keV(np.atleast_1d(E_keV))                      # (nE,)
    flat = dirs.reshape(-1, 3)
    out = np.empty((flat.shape[0], om.size))
    for i0 in range(0, flat.shape[0], chunk):
        n = flat[i0:i0 + chunk]                                  # (g,3)
        ndotb = n @ beta.T                                       # (g,Ns)
        ndotr = n @ r.T                                          # (g,Ns)
        Vperp = n[:, None, :] * ndotb[:, :, None] - beta[None, :, :]   # (g,Ns,3)
        arg = t[None, :] - ndotr / C                             # (g,Ns)
        # A(g,nE,3) = Σ_s Vperp(g,Ns,3) * exp(i om[nE] arg(g,Ns)) * w[Ns]
        phase = np.exp(1j * om[None, None, :] * arg[:, :, None]) * w[None, :, None]  # (g,Ns,nE)
        A = np.einsum("gsc,gse->gec", Vperp, phase)             # (g,nE,3)
        out[i0:i0 + chunk] = (om[None, :, None] ** 2 * np.abs(A) ** 2).sum(axis=2)
    return out.reshape(*dirs.shape[:-1], om.size)


def on_axis_spectrum(tr, E_keV):
    """Спектр на оси (одно направление n=ẑ): shape (nE,)."""
    return amplitude(tr, E_keV, _dirs(0.0, 0.0)[None, :])[0]


def calibrate(tr, p=sr.P):
    """Множитель Iraw→(ф/с/мрад²/0.1%BW): привязка пика 1-й гармоники на оси к flux_density_onaxis(1,K)."""
    E1 = sr.harmonic_energy_keV(1, tr["K"], 0.0, p)
    Es = np.linspace(0.90 * E1, 1.05 * E1, 400)
    s = on_axis_spectrum(tr, Es)
    Iraw1 = float(s.max())
    target = sr.flux_density_onaxis(1, tr["K"], p)
    return target / Iraw1 if Iraw1 > 0 else 0.0, E1


def selftest(gap, ppp):
    B0 = sr.B0_of_gap(gap); K = sr.K_of_B0(B0); g = sr.gamma()
    tr = trajectory(K, sr.P.lam_u, sr.P.N, g, ppp)
    scale, E1 = calibrate(tr)
    print(f"ядро поля: зазор {gap:g} мм → B0={B0:.3f} Т, K={K:.3f}, γ={g:.0f}, N={sr.P.N}, "
          f"точек/период={ppp} (Ns={tr['r'].shape[0]})")
    print(f"калибровка: пик 1-й гарм. на оси ≡ flux_density_onaxis(1,K)={sr.flux_density_onaxis(1,K):.3e}\n")
    print(f"  {'гарм n':>6} | {'E (sr), кэВ':>11} | {'E (интеграл)':>12} | {'поток/ось (sr)':>14} | "
          f"{'поток (интеграл)':>15} | {'откл':>6}")
    for n in (1, 2, 3, 5, 7):
        En = sr.harmonic_energy_keV(n, K, 0.0)
        Es = np.linspace(0.97 * En, 1.03 * En, 240)
        s = on_axis_spectrum(tr, Es) * scale
        j = int(np.argmax(s)); Epk = Es[j]; Fpk = float(s[j])
        Fsr = sr.flux_density_onaxis(n, K)                       # 0 для чётных
        dev = (Fpk - Fsr) / Fsr * 100 if Fsr > 0 else float("nan")
        tagE = (Epk - En) / En * 100
        print(f"  {n:>6} | {En:>11.3f} | {Epk:>10.3f} ({tagE:+.1f}%) | {Fsr:>14.3e} | {Fpk:>15.3e} | "
              f"{dev:>5.1f}%" if Fsr > 0 else
              f"  {n:>6} | {En:>11.3f} | {Epk:>10.3f} ({tagE:+.1f}%) | {'0 (чётная)':>14} | {Fpk:>15.3e} | "
              f"{'подавл':>6}")
    print("\nпроверка: нечётные — энергии и потоки совпадают с undulator_sr (интеграл независим, "
          "нормирован только n=1); чётные на оси подавлены.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap", type=float, default=5.0)
    ap.add_argument("--ppp", type=int, default=64, help="точек траектории на период")
    a = ap.parse_args()
    selftest(a.gap, a.ppp)
