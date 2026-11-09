# -*- coding: utf-8 -*-
"""Излучение ПОВОРОТНОГО МАГНИТА (диполя) — для волноотводов (экранов на поворотном СИ).

В отличие от ондулятора: НЕТ гармоник (сплошной спектр), свет — вертикально узкий (~1/γ) ВЕЕР, широкий по
горизонтали (= угол изгиба θ_bend, электрон «разворачивает» касательную вдоль дуги). Пучок отбирают под
конкретным УГЛОМ ОТБОРА φ (касательная в точке дуги смотрит на линию). Формулы — стандартная теория СИ
поворота (K-функции Бесселя 2/3, 5/3; Wiedemann / X-ray Data Booklet).

ε_c = 0.665·E²[ГэВ]·B[Т] кэВ;  поток(гор.) ∝ G1(y)=y∫_y^∞K_{5/3};  угл.плотность ∝ K_{2/3},K_{1/3};
полная мощность P[кВт]=1.266·E²·B²·L·I.

Самотест: py -3.14 undulators/bend_sr.py --theta -9.573e-3 --L 0.4178
"""
import sys, math, argparse, functools
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy.special import kv
from scipy.integrate import quad

MC2 = 0.510998950e-3            # ГэВ


def brho(E):
    return E / 0.299792458      # Т·м


def gamma(E):
    return E / MC2


def bend_geom(theta, L, E):
    """ρ[м], B[Т], ε_c[кэВ] из геометрического угла θ и длины L при энергии E."""
    rho = L / theta if theta else float("inf")
    B = brho(E) / abs(rho) if rho and math.isfinite(rho) else 0.0
    ec = 0.665 * E ** 2 * B
    return rho, B, ec


def total_power_kW(B, L, E, I):
    return 1.266 * E ** 2 * B ** 2 * L * I


@functools.lru_cache(maxsize=4096)
def _G1(y):
    """G1(y)=y·∫_y^∞ K_{5/3}(t)dt — спектральная функция поворотного СИ (гор.-интегрированная)."""
    if y <= 0 or y > 60:
        return 0.0
    val, _ = quad(lambda t: kv(5.0 / 3.0, t), y, np.inf, limit=200)
    return y * val


def flux_horiz(E_keV, Ec, E_beam, I):
    """Поток на 1 мрад ГОРИЗОНТА: ф/с/мрад/0.1%BW = 2.457e13·E·I·G1(E/Ec)."""
    y = np.atleast_1d(E_keV) / Ec
    g = np.array([_G1(float(v)) for v in y])
    return 2.457e13 * E_beam * I * g


def flux_angular(E_keV, psi, Ec, E_beam, I, g):
    """Угловая плотность потока на верт. угле ψ[рад]: ф/с/мрад²/0.1%BW (σ+π-компоненты)."""
    y = np.atleast_1d(E_keV) / Ec
    X = g * psi
    xi = (y / 2.0) * (1 + X ** 2) ** 1.5
    K23 = kv(2.0 / 3.0, xi); K13 = kv(1.0 / 3.0, xi)
    return 1.327e13 * E_beam ** 2 * I * y ** 2 * (1 + X ** 2) ** 2 * (K23 ** 2 + (X ** 2 / (1 + X ** 2)) * K13 ** 2)


def power_vprofile(psi, g):
    """НОРМИРОВАННЫЙ вертикальный профиль плотности мощности g(ψ), ∫g dψ=1 (пик на плоскости орбиты)."""
    X = g * np.asarray(psi, float)
    f = (1 + X ** 2) ** -2.5 * (1 + (5.0 / 7.0) * X ** 2 / (1 + X ** 2))
    return f


def vprofile_norm(g, half=8.0, n=2001):
    """Константа нормировки ∫ power_vprofile dψ по ψ∈[−half/γ, half/γ]."""
    psi = np.linspace(-half / g, half / g, n)
    f = power_vprofile(psi, g)
    return np.trapezoid(f, psi)


def selftest(theta, L, E, I):
    g = gamma(E); rho, B, ec = bend_geom(theta, L, E)
    P = total_power_kW(B, L, E, I)
    print(f"диполь: θ={theta*1e3:+.3f} мрад, L={L:.4f} м, E={E} ГэВ, I={I} А, γ={g:.0f}")
    print(f"  ρ={rho:+.3f} м, B={B:.4f} Т, ε_c={ec:.3f} кэВ, полная мощность {P*1e3:.1f} Вт "
          f"(плотн. {P*1e3/(abs(theta)*1e3):.2f} Вт/мрад гор.)")
    print(f"  верт. раскрыв фотонов ~1/γ={1e6/g:.1f} мкрад")
    for frac in (0.1, 0.3, 1.0, 3.0):
        Ek = frac * ec
        fh = flux_horiz(Ek, ec, E, I)[0]
        fa = flux_angular(Ek, 0.0, ec, E, I, g)[0]
        print(f"  E={Ek:7.3f} кэВ (y={frac:.1f}): поток/гор {fh:.3e} ф/с/мрад/0.1%, "
              f"на оси {fa:.3e} ф/с/мрад²/0.1%")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", type=float, default=-9.573e-3, help="угол изгиба, рад")
    ap.add_argument("--L", type=float, default=0.4178)
    ap.add_argument("--E", type=float, default=6.0)
    ap.add_argument("--I", type=float, default=0.2)
    a = ap.parse_args()
    selftest(a.theta, a.L, a.E, a.I)
