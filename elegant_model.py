# -*- coding: utf-8 -*-
"""elegant_model.py — порт канонического мультипольного трекинга elegant (KQUAD/KSEXT/KOCT)
из ELEGANT/elegant-master/src/multipole.{c,h}. Симплектический интегратор Форест-Рута
4-го порядка (drift-kick) + стандартное мультипольное ядро Re/Im[(x+iy)ⁿ].

Аналог bip_matrix.py, но для elegant: позволяет сравнить модель elegant с нашим Boris
теми же средствами. elegant принимает физические силы K_n (а не токи обмоток).

Состояние частицы [x, x', y, y']; x — плоскость с положит. K1-фокусировкой. δ=0 (монохром.).
"""
import math

# --- интегратор Форест-Рута 4-го порядка (multipole.c) ---
BETA = 2.0 ** (1.0 / 3.0)                       # = 1.259921... (#define BETA в elegant)
_D = 2.0 - BETA
DRIFT4 = [0.5 / _D, (1 - BETA) / _D / 2, (1 - BETA) / _D / 2, 0.5 / _D]
KICK4 = [1.0 / _D, -BETA / _D, 1.0 / _D, 0.0]

# --- expansion_coefficients(n): coef[i] = (-1)^floor(i/2) / (i!·(n-i)!) (multipole.h) ---
COEF = {
    1: [1.0, 1.0],                              # квадруполь
    2: [1.0 / 2, 1.0, -1.0 / 2],                # секступоль:  By∝x²−y²
    3: [1.0 / 6, 1.0 / 2, -1.0 / 2, -1.0 / 6],  # октуполь:    By∝x³−3xy² (СТАНДАРТ)
}
ORDER = {"quad": 1, "sext": 2, "oct": 3}


def field_sums(x, y, order):
    """Стандартное мультипольное ядро: ΣFx (нечётные i), ΣFy (чётные i) из Re/Im[(x+iy)ⁿ].
    Общее для elegant/MAD/BDSIM (у всех ядро стандартное). Кик: Δqx=−KnL·ΣFy, Δqy=+KnL·ΣFx."""
    coef = COEF[order]
    sFx = sFy = 0.0
    for i in range(order + 1):
        f = coef[i] * (x ** (order - i)) * (y ** i)
        if i & 1:
            sFx += f
        else:
            sFy += f
    return sFx, sFy


def _kick(qx, qy, x, y, order, KnL, skew=False):
    """apply_canonical_multipole_kicks_ret: qx-=KnL·ΣFy, qy+=KnL·ΣFx (multipole.h)."""
    sFx, sFy = field_sums(x, y, order)
    if skew:
        sFx, sFy = -sFy, sFx
    return qx - KnL * sFy, qy + KnL * sFx


def _slopes_to_momenta(xp, yp):
    """convertSlopesToMomenta (δ=0): qx=xp/√(1+x'²+y'²)."""
    d = math.sqrt(1.0 + xp * xp + yp * yp)
    return xp / d, yp / d


def _momenta_to_slopes(qx, qy):
    """convertMomentaToSlopes (δ=0): x'=qx/√(1−qx²−qy²)."""
    d = math.sqrt(max(1e-18, 1.0 - qx * qx - qy * qy))
    return qx / d, qy / d


def track_multipole(particles, kind, KnL, L, d1, d2, n_slices=20):
    """Прогон [x,x',y,y'] через drift(d1) + мультиполь(KnL,L) + drift(d2) методом elegant.
    KnL — интегральная сила K_n·L (для kick-ядра); n_slices ломтей × 4 подшага Форест-Рута."""
    order = ORDER[kind]
    ds = L / n_slices
    knl_s = KnL / n_slices
    out = []
    for (x, xp, y, yp) in particles:
        x += xp * d1; y += yp * d1                       # дрейф до
        qx, qy = _slopes_to_momenta(xp, yp)
        for _ in range(n_slices):
            for st in range(4):
                xp, yp = _momenta_to_slopes(qx, qy)
                dsh = ds * DRIFT4[st]
                x += xp * dsh; y += yp * dsh
                if KICK4[st] == 0.0:
                    break
                qx, qy = _kick(qx, qy, x, y, order, knl_s * KICK4[st])
        xp, yp = _momenta_to_slopes(qx, qy)
        x += xp * d2; y += yp * d2                       # дрейф после
        out.append((x, xp, y, yp))
    return out


def trajectory(particle, kind, KnL, L, d1, d2, n_slices=30, n_draw=30):
    """Полилиния (X,Y,Z) для рисунка: дрейф d1 + элемент (по ломтям) + дрейф d2."""
    order = ORDER[kind]
    x, xp, y, yp = particle
    pts = [(0.0, x, y)]
    x += xp * d1; y += yp * d1
    pts.append((d1, x, y))
    qx, qy = _slopes_to_momenta(xp, yp)
    ds = L / n_draw; knl_s = KnL / n_draw
    for j in range(n_draw):
        for st in range(4):
            xp, yp = _momenta_to_slopes(qx, qy)
            dsh = ds * DRIFT4[st]
            x += xp * dsh; y += yp * dsh
            if KICK4[st] == 0.0:
                break
            qx, qy = _kick(qx, qy, x, y, order, knl_s * KICK4[st])
        pts.append((d1 + ds * (j + 1), x, y))
    xp, yp = _momenta_to_slopes(qx, qy)
    x += xp * d2; y += yp * d2
    pts.append((d1 + L + d2, x, y))
    return pts
