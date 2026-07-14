#!/usr/bin/env python3
"""single_element_boris_fld.py — DROP-IN замена single_element_boris для НОВЫХ карт KD (.fld).
Тот же API (run/axis_field_integral/map_length, WS/TOTAL/DRIFT/N_SAT, пучок 37 частиц), но
трекинг НАШИМ Python-Boris прямо по .fld (без C++ и без .bin). Нелинейные (SP/OP) — по мелкой
карте v3 (2мм), чтобы избежать артефакта интерполяции грубой сетки; квадруполи/DQ — по полной
v2 (линейное поле интерполируется точно); диполи/SP/OP — полукарты X≥0 (зеркало X<0).

Геометрия как в оригинале: старт пучка x=0, ЦЕНТР элемента x=DRIFT=1м, финиш x=TOTAL=5м.
Ради скорости Boris-интегрируем ТОЛЬКО зону поля (элемент±margin), дрейфы — аналитически
(прямые): до элемента пучок параллелен (y'=0), после — сохраняем угол."""
import os, math
from pathlib import Path
import numpy as np
from field_quality import load
from _ct_asym_track import make_field, C, ME, QE, V, GAMMA0

WS = Path(__file__).resolve().parent / "workspace"
KD = Path(r"C:/Sasha/Work/SILA/ansis approximations/KD")
TOTAL = 5.0
DRIFT = 1.0
RADII = (0.0005, 0.001, 0.002)
N_PER = 12
N_SAT = len(RADII) * N_PER            # 36 спутников + центр = 37
V_BEAM = V
DX = 1e-3                             # шаг Boris в зоне поля

# device -> (относит. путь .fld, halfmap(зеркало X<0), nominal_A|None(=коэф.))
MAPS = {
    "QP_L_1A":  ("QP_L/QP_L_100A_full_v2.fld",  False, 100.0),
    "QP_sh_1A": ("QP_sh/QP_sh_100A_full_v2.fld", False, 100.0),
    "SP_I1_1A": ("SP/SP_I1_100A_v3_symsm.bin",   True,  100.0),   # симметризов.+сглажено
    "SP_I2_1A": ("SP/SP_I2_100A_v4_symsm.bin",   True,  100.0),   # v4 симметризов.+сглажено (рабочая)
    "OP_1A":    ("OP/OP_100A_v3_sym.bin",        True,  100.0),   # симметризованная (сглаж. пере-глаживает r³)
    "DQ_I1_1A": ("DQ/DQ_I1_100A_v2_sym.bin",     False, 100.0),   # симметризов. (средняя плоскость z=0)
    "DQ_I2_1A": ("DQ/DQ_I2_100A_v2_sym.bin",     False, 100.0),
    "DQ_I3_1A": ("DQ/DQ_I3_100A_v2_sym.bin",     False, 100.0),
    "DP_L_PM":  ("DP_L/dp_L_PM_1A.fld",          True,  None),
    "DP_sh_PM": ("DP_sh/dp_sh_PM_1mm.fld",       True,  None),
}
_CACHE = {}


def _get(device):
    if device not in _CACHE:
        rel, half, nom = MAPS[device]
        h, B = load(str(KD / rel))
        fn = make_field(h, B, "mirror" if half else "full")
        _CACHE[device] = (h, fn, h["mxx"], nom)
    return _CACHE[device]


def map_length(device):
    return _get(device)[2]                                    # maxX карты, м


def _satellites():
    sats = []
    for r in RADII:
        for k in range(N_PER):
            th = 2 * math.pi * k / N_PER
            sats.append((r * math.cos(th), r * math.sin(th)))
    return sats


def _scale(device, value, is_coeff):
    _, _, _, nom = _get(device)
    return value if (is_coeff or nom is None) else value / nom


def _field(device, value, is_coeff):
    """Глобальная функция поля: карта центрирована в x=DRIFT, масштаб по току/коэф., 0 вне."""
    _, fn, _, _ = _get(device)
    s = _scale(device, value, is_coeff)
    def f(pos):
        return fn(np.array([pos[0] - DRIFT, pos[1], pos[2]])) * s
    return f


def _track(y0, z0, ffield, xlo, xhi, x_fin, dx=DX):
    """Boris (x,y,z) от xlo до xhi через ffield, затем аналитич. дрейф до x_fin.
    Возврат (Y,Y',Z,Z') на финиш-плоскости."""
    pos = np.array([xlo, y0, z0]); mom = GAMMA0 * ME * np.array([V, 0.0, 0.0])
    dt = dx / V; nstep = int(round((xhi - xlo) / dx))
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = ffield(pos)
        gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
    yp = mom[1] / mom[0]; zp = mom[2] / mom[0]
    xe = pos[0]
    return (pos[1] + yp * (x_fin - xe), yp, pos[2] + zp * (x_fin - xe), zp)


def _track_maxr(y0, z0, ffield, xlo, xhi, dx=DX):
    """Макс. поперечный радиус √(y²+z²), достигнутый В ЗОНЕ ПОЛЯ (для авто-масштаба свипа)."""
    pos = np.array([xlo, y0, z0]); mom = GAMMA0 * ME * np.array([V, 0.0, 0.0])
    dt = dx / V; nstep = int(round((xhi - xlo) / dx)); mr = math.hypot(y0, z0)
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = ffield(pos)
        gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        r = math.hypot(pos[1], pos[2])
        if r > mr: mr = r
    return mr


def probe_excursion(device, value, is_coeff):
    """Макс. поперечный вынос внешнего кольца (r=2мм) + центра в зоне поля при данном value.
    Ловит и изгиб (центр), и фокус/дефокус (спутники) — критерий безопасного свипа."""
    _, _, maxX, _ = _get(device); ffield = _field(device, value, is_coeff)
    xlo = DRIFT - maxX - 0.02; xhi = DRIFT + maxX + 0.02
    pts = [(0.0, 0.0)] + [(0.002 * math.cos(a), 0.002 * math.sin(a))
                          for a in (0, math.pi / 2, math.pi, 3 * math.pi / 2)]
    return max(_track_maxr(y0, z0, ffield, xlo, xhi) for y0, z0 in pts)


def run(device, value, is_coeff=False, dt=None):
    """Прогон 37 частиц. Возврат dict: L(=maxX), x_fin, value, init[(Y0,Z0)], final[(Y,Y',Z,Z')]."""
    _, _, maxX, _ = _get(device)
    ffield = _field(device, value, is_coeff)
    xlo = DRIFT - maxX - 0.02; xhi = DRIFT + maxX + 0.02          # зона поля (± margin)
    init = [(0.0, 0.0)] + _satellites()
    final = [_track(y0, z0, ffield, xlo, xhi, TOTAL) for (y0, z0) in init]
    return dict(device=device, value=value, is_coeff=is_coeff, L=maxX, x_fin=TOTAL,
                init=init, final=final, drift=DRIFT)


def _field_multi(devices, currents):
    """Суперпозиция полей нескольких обмоток в одной позиции (центр x=DRIFT), масштаб по току."""
    fns = []
    for dev, cur in zip(devices, currents):
        _, fn, _, nom = _get(dev); fns.append((fn, cur / nom))
    def f(pos):
        lp = np.array([pos[0] - DRIFT, pos[1], pos[2]])
        out = np.zeros(3)
        for fn, s in fns: out = out + fn(lp) * s
        return out
    return f


def run_multi(devices, currents):
    """Прогон 37 частиц через СУММУ полей обмоток (для диполь-квадруполя: I1/I2/I3 сразу)."""
    _, _, maxX, _ = _get(devices[0])
    ffield = _field_multi(devices, currents)
    xlo = DRIFT - maxX - 0.02; xhi = DRIFT + maxX + 0.02
    init = [(0.0, 0.0)] + _satellites()
    final = [_track(y0, z0, ffield, xlo, xhi, TOTAL) for (y0, z0) in init]
    return dict(devices=devices, currents=currents, L=maxX, x_fin=TOTAL,
                init=init, final=final, drift=DRIFT)


def traj(device, value, y0, z0, is_coeff=False, dx=DX):
    """Полная траектория [(x,y,z)] частицы 0..TOTAL: прямой дрейф-вход + Boris(зона поля) +
    прямой дрейф-выход. Для рисунков сравнения с матричной моделью."""
    _, _, maxX, _ = _get(device); ffield = _field(device, value, is_coeff)
    xlo = DRIFT - maxX - 0.02; xhi = DRIFT + maxX + 0.02
    pts = [(0.0, y0, z0), (xlo, y0, z0)]                          # дрейф-вход (пучок параллелен)
    pos = np.array([xlo, y0, z0]); mom = GAMMA0 * ME * np.array([V, 0.0, 0.0])
    dt = dx / V; nstep = int(round((xhi - xlo) / dx))
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = ffield(pos)
        gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        pts.append((pos[0], pos[1], pos[2]))
    yp = mom[1] / mom[0]; zp = mom[2] / mom[0]
    pts.append((TOTAL, pos[1] + yp * (TOTAL - pos[0]), pos[2] + zp * (TOTAL - pos[0])))  # дрейф-выход
    return np.array(pts)


def traj_multi(devices, currents, y0, z0, dx=DX):
    """Полная траектория частицы через СУММУ полей обмоток (для DQ combined-возбуждения)."""
    _, _, maxX, _ = _get(devices[0]); ffield = _field_multi(devices, currents)
    xlo = DRIFT - maxX - 0.02; xhi = DRIFT + maxX + 0.02
    pts = [(0.0, y0, z0), (xlo, y0, z0)]
    pos = np.array([xlo, y0, z0]); mom = GAMMA0 * ME * np.array([V, 0.0, 0.0])
    dt = dx / V; nstep = int(round((xhi - xlo) / dx))
    for _ in range(nstep):
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        Bv = ffield(pos)
        gm = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2)
        t = (QE * dt / (2 * ME * gm)) * Bv; s = 2 * t / (1 + t.dot(t))
        pp = mom + np.cross(mom, t); mom = mom + np.cross(pp, s)
        g = math.sqrt(1 + mom.dot(mom) / (ME * C) ** 2); pos = pos + mom * (dt / (2 * g * ME))
        pts.append((pos[0], pos[1], pos[2]))
    yp = mom[1] / mom[0]; zp = mom[2] / mom[0]
    pts.append((TOTAL, pos[1] + yp * (TOTAL - pos[0]), pos[2] + zp * (TOTAL - pos[0])))
    return np.array(pts)


def axis_field_integral(device, value, is_coeff):
    """∫By·dl, ∫Bz·dl (Тл·м) приложенного поля вдоль оси по ПОЛНОЙ длине элемента."""
    _, fn, maxX, _ = _get(device)
    s = _scale(device, value, is_coeff)
    xs = np.linspace(-maxX + 1e-6, maxX - 1e-6, 801)
    By = np.array([fn(np.array([x, 0.0, 0.0]))[1] for x in xs])
    Bz = np.array([fn(np.array([x, 0.0, 0.0]))[2] for x in xs])
    return float(np.trapezoid(By, xs) * s), float(np.trapezoid(Bz, xs) * s)


if __name__ == "__main__":
    import sys
    dev = sys.argv[1] if len(sys.argv) > 1 else "QP_L_1A"
    val = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    isc = dev.endswith("_PM")
    r = run(dev, val, is_coeff=isc)
    cy, cyp, cz, czp = r["final"][0]
    print(f"{dev}@{val}: L_map={r['L']:.3f}м  центр Y={cy*1e3:+.4f}мм Y'={cyp*1e3:+.4f}мрад "
          f"Z={cz*1e3:+.4f}мм Z'={czp*1e3:+.4f}мрад")
    for i in (1, 13, 25):
        y0, z0 = r["init"][i]; y, yp, z, zp = r["final"][i]
        print(f"  #{i} r={math.hypot(y0,z0)*1e3:.1f}мм: dY={(y-cy)*1e3:+.4f} dY'={(yp-cyp)*1e3:+.4f} "
              f"dZ={(z-cz)*1e3:+.4f} dZ'={(zp-czp)*1e3:+.4f} (отн.центра)")
