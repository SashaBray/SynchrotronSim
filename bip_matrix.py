#!/usr/bin/env python3
"""bip_matrix.py — Python-копия матричной модели BIP (transfer-matrix оптика) для
сравнения с нашим Boris-трекером на полевых картах.

Формулы СКОПИРОВАНЫ из BIP (BIP не трогаем):
  - drift          : [[1,L],[0,1]]                       (driftspace.h)
  - quad/DQ (лин.) : positive_k cos/sin·√κ, negative_k cosh/sinh   (quadrupole.h, dipolequadrupole.h)
                     κx = k + θ²/L²,  κy = −k,  k=(q/p)·G            (omegax²/omegay², см. ниже)
  - sextupole/oct. : нормированный RK4 ровно как в BIP (sextrupole.h, octopole.h),
                     lx=5e-4, lalfa=lx·1e-2, ls=L, 10 шагов, RK4
  - kicker         : x'+=dx, y'+=dy                       (kicker.h)

ВАЖНО (энергия): импульс p берём из НАШЕГО пучка v=2.967e8 м/с (γ≈6.98), НЕ из примера
BIP (3 МэВ). k=(q/p)G; из .lte G=−Kn·p/q ⇒ k=−Kn (Kn — нормированная сила из .lte).

Связь BIP-κ (focus-convention, κ>0 ⇒ фокусировка cos/sin) и BIP-матриц:
  x-плоскость: omegax²=θ²/L²+Gθ/(BL); при самосогласованном изгибе (B=θp/(qL)) ⇒ θ²/L²+k.
               omegax²>0 → positive_k (фокус) ⇒ κx=omegax²=k+θ²/L².
  y-плоскость: omegay²=Gθ/(BL)=k; BIP при omegay²>0 берёт negative_k (ДЕфокус) ⇒ κy=−k.
Для чистого квадруполя θ=0 ⇒ κx=k, κy=−k (классика). Диполь k=0 ⇒ κx=θ²/L², κy=0.
"""
import math

# ----------------------------- пучок / энергия -----------------------------
C_LIGHT = 299792458.0
ME      = 9.1093837e-31          # масса электрона, кг
E_CH    = 1.602176634e-19        # элементарный заряд, Кл
Q       = -E_CH                  # заряд электрона
V_BEAM  = 2.967e8                # м/с — из particle_groups.csv (наш пучок)

BETA  = V_BEAM / C_LIGHT
GAMMA = 1.0 / math.sqrt(1.0 - BETA * BETA)
P_MOM = GAMMA * ME * V_BEAM      # импульс, кг·м/с
BRHO  = P_MOM / E_CH             # магнитная жёсткость, Тл·м
QP    = Q / P_MOM                # q/p, 1/(Тл·м)  (= −1/Brho)

# характерные константы нормировки BIP (bunchimitationblock.h)
LX    = 5e-4
LALFA = LX * 1e-2

def G_from_K(Kn):
    """Нормированная сила .lte Kn → физический градиент G (СИ). G = −Kn·p/q."""
    return -Kn * P_MOM / Q

def K_from_G(G):
    """Физический градиент G → нормированная сила Kn (= −(q/p)G = −k)."""
    return -QP * G

def k_from_G(G):
    """k = (q/p)·G  [1/m^(n+1)], внутренняя сила BIP."""
    return QP * G


# --------- матрица связи токи→(B,G) диполь-квадруполя (калибровка по Boris) ---------
# Аддитивная модель BIP (B=BI(a1·I1), G=GI(a2·I2+a3·I3)) занижает совпадение с Boris:
# она считает, что гнёт только I1, а фокусирует только I2,I3. Реально ВСЕ три обмотки
# дают и изгиб, и градиент (combined-function). Единицы: на ампер.
# ПЕРЕКАЛИБРОВАНО на НОВЫХ картах KD (_dq_coupling_refit.py, L_eff=0.434 м):
# combined-возбуждения → Boris; поз. RMS 116→33 мм (×3.6), R² 0.53→0.96 vs аддитивная.
DQ_COUPLING_B = (5.692757e-05, 5.093815e-05, 4.923006e-05)     # B[Тл]   от (I1,I2,I3)
DQ_COUPLING_G = (-5.137492e-03, 4.750480e-03, 4.807284e-03)    # G[Тл/м] от (I1,I2,I3)
# прежняя (старые карты, L=1.0): B=(6.27e-6,4.49e-6,4.11e-6) G=(-4.11e-4,3.29e-4,3.24e-4)
# СЕКСТУПОЛЬНАЯ строка: у DQ есть реальный секступоль — при бенде пучок уходит off-axis, эфф.
# градиент меняется (feed-down G(орбиты)=G(0)+S·y). С распределённым секступ. киком в трекере DQ
# (_dq_coupling_sext.py, nsl срезов drift-kick-drift) и траекторно-калиброванной S матрица связи
# даёт R²=0.9999 (RMS 1мм) vs 0.998 без S — улучшение ×4.6. I2≈I3, у I1 слабее. ∂²Bz[Тл/м²] на ампер.
# (поле-фит давал -0.057,-0.192,-0.188; траект.-фит крупнее — впитывает feed-down/высш. поправки.)
DQ_COUPLING_S = (-1.887e-01, -2.819e-01, -2.709e-01)           # ∂²Bz[Тл/м²] от (I1,I2,I3), траект.-калибр.


def dq_BG_from_currents(I1, I2, I3, mB=DQ_COUPLING_B, mG=DQ_COUPLING_G):
    """Физическое (B, G) диполь-квадруполя через матрицу связи (а не аддитивную модель)."""
    B = mB[0] * I1 + mB[1] * I2 + mB[2] * I3
    G = mG[0] * I1 + mG[1] * I2 + mG[2] * I3
    return B, G


def dq_kt_from_currents(I1, I2, I3, L, mB=DQ_COUPLING_B, mG=DQ_COUPLING_G):
    """(k, θ) для DQ-матрицы из токов обмоток: θ=(q/p)·B·L, k=(q/p)·G."""
    B, G = dq_BG_from_currents(I1, I2, I3, mB, mG)
    return QP * G, QP * B * L


# ----------------------------- линейные элементы -----------------------------
def _plane(u, up, kappa, L):
    """2×2 перенос (u,u') на длину L при силе kappa (focus-convention: kappa>0 ⇒ фокус)."""
    if kappa > 1e-12:
        w = math.sqrt(kappa); wl = w * L
        c, s = math.cos(wl), math.sin(wl)
        return c * u + (s / w) * up, -w * s * u + c * up
    if kappa < -1e-12:
        w = math.sqrt(-kappa); wl = w * L
        ch, sh = math.cosh(wl), math.sinh(wl)
        return ch * u + (sh / w) * up, w * sh * u + ch * up
    return u + L * up, up               # дрейф


def drift(state, L):
    x, xp, y, yp = state
    return (x + L * xp, xp, y + L * yp, yp)


def linear_element(state, k, theta, L):
    """Квадруполь/диполь-квадруполь/диполь в криволинейных координатах (δ=0).
    κx = k + θ²/L²,  κy = −k.  k=(q/p)G; θ — геометрический угол изгиба (0 для квадруполя)."""
    x, xp, y, yp = state
    kx = k + ((theta * theta) / (L * L) if L > 0 else 0.0)
    ky = -k
    x, xp = _plane(x, xp, kx, L)
    y, yp = _plane(y, yp, ky, L)
    return (x, xp, y, yp)


def kicker(state, dx, dy):
    x, xp, y, yp = state
    return (x, xp + dx, y, yp + dy)


def track_planes(particles, ky, kz, L, d1, d2):
    """Прогон с НЕЗАВИСИМОЙ фокусировкой по плоскостям (κy в Y, κz в Z) + дрейфы d1,d2.
    Для разбора диполя: слабая фокусировка в плоскости изгиба vs краевая в перпенд."""
    out = []
    for (x, xp, y, yp) in particles:
        x, xp = _plane(x, xp, 0.0, d1)          # дрейф d1
        y, yp = _plane(y, yp, 0.0, d1)
        x, xp = _plane(x, xp, ky, L)            # элемент: своя κ в каждой плоскости
        y, yp = _plane(y, yp, kz, L)
        x, xp = _plane(x, xp, 0.0, d2)          # дрейф d2
        y, yp = _plane(y, yp, 0.0, d2)
        out.append((x, xp, y, yp))
    return out


# ----------------------- нелинейные элементы (RK4 как в BIP) -----------------------
def _rhs_sext(X, k, A, B, C):
    x0, x1, x2, x3 = X[0], X[1], X[2], X[3]
    dzds = math.sqrt(1.0 / C + x1 * x1 + x3 * x3)
    r0 = A * x1
    r1 = B * k * dzds * (-(1.0 / C + x1 * x1) * (x0 * x0 - x2 * x2) / 2.0 + x0 * x1 * x2 * x3)
    r2 = A * x3
    r3 = -B * k * dzds * (-(1.0 / C + x3 * x3) * x0 * x2 + x1 * x3 * (x0 * x0 - x2 * x2) / 2.0)
    return (r0, r1, r2, r3, LALFA * dzds, 0.0)


def _rhs_oct(X, k, A, B, C):
    x0, x1, x2, x3 = X[0], X[1], X[2], X[3]
    dzds = math.sqrt(1.0 / C + x1 * x1 + x3 * x3)
    By = x0 ** 3 - 3.0 * x0 * x0 * x2          # ВЕРБАТИМ из BIP octopole.h (нестандартная форма)
    Bx = 3.0 * x0 * x2 * x2 - x2 ** 3
    r0 = A * x1
    r1 = B * k * dzds * (-(1.0 / C + x1 * x1) * By + x1 * x3 * Bx) / 6.0
    r2 = A * x3
    r3 = -B * k * dzds * (-(1.0 / C + x3 * x3) * Bx + x1 * x3 * By) / 6.0
    return (r0, r1, r2, r3, LALFA * dzds, 0.0)


def _rhs_oct_std(X, k, A, B, C):
    """СТАНДАРТНЫЙ октуполь из потенциала Re[(x+iy)⁴]: By∝x³−3xy², Bx∝3x²y−y³
    (для проверки гипотезы, что ядро BIP `x³−3x²y` — опечатка)."""
    x0, x1, x2, x3 = X[0], X[1], X[2], X[3]
    dzds = math.sqrt(1.0 / C + x1 * x1 + x3 * x3)
    By = x0 ** 3 - 3.0 * x0 * x2 * x2           # x³ − 3xy²  (стандарт)
    Bx = 3.0 * x0 * x0 * x2 - x2 ** 3           # 3x²y − y³
    r0 = A * x1
    r1 = B * k * dzds * (-(1.0 / C + x1 * x1) * By + x1 * x3 * Bx) / 6.0
    r2 = A * x3
    r3 = -B * k * dzds * (-(1.0 / C + x3 * x3) * Bx + x1 * x3 * By) / 6.0
    return (r0, r1, r2, r3, LALFA * dzds, 0.0)


def rk_multipole(state, kind, k, L, n_iter=10):
    """Секступоль (kind='S') / октуполь BIP (kind='O') / октуполь СТАНДАРТНЫЙ (kind='Ostd').
    Нормированные координаты, RK4, n_iter шагов по s∈[0,1]. k=(q/p)G (G в Тл/м² или Тл/м³)."""
    ls = L
    A = ls * LALFA / LX
    C = LALFA * LALFA
    if kind == 'S':
        B = ls * LX * LX * LALFA * LALFA
        rhs = _rhs_sext
    elif kind == 'Ostd':
        B = ls * LX * LX * LX * LALFA * LALFA
        rhs = _rhs_oct_std
    else:                                        # 'O' — ядро BIP (нестандартное)
        B = ls * LX * LX * LX * LALFA * LALFA
        rhs = _rhs_oct
    # нормировка координат: X0=x/lx, X1=x'/lalfa, X2=y/lx, X3=y'/lalfa
    x, xp, y, yp = state
    X = [x / LX, xp / LALFA, y / LX, yp / LALFA, 0.0, 0.0]
    h = 1.0 / n_iter
    for _ in range(n_iter):
        k1 = rhs(X, k, A, B, C)
        k2 = rhs([X[i] + 0.5 * h * k1[i] for i in range(6)], k, A, B, C)
        k3 = rhs([X[i] + 0.5 * h * k2[i] for i in range(6)], k, A, B, C)
        k4 = rhs([X[i] + h * k3[i] for i in range(6)], k, A, B, C)
        X = [X[i] + (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6)]
    return (X[0] * LX, X[1] * LALFA, X[2] * LX, X[3] * LALFA)


# ----------------------------- сборка элемента -----------------------------
def track_element(particles, spec):
    """Прогон списка частиц [x,x',y,y'] через drift(d1)+элемент+drift(d2).
    spec: {'kind': 'quad'|'dq'|'dipole'|'sext'|'oct'|'kicker', 'L','d1','d2', и коэф.}
      quad/dq/dipole: 'k' (=(q/p)G), 'theta' (рад, 0 для quad)
      sext/oct       : 'k'
      kicker         : 'dx','dy' (рад), L игнор (тонкий)
    Возвращает список финальных [x,x',y,y']."""
    kind = spec['kind']; d1 = spec.get('d1', 1.0); d2 = spec.get('d2', 1.0)
    L = spec.get('L', 0.0)
    out = []
    for p in particles:
        s = drift(p, d1)
        if kind in ('quad', 'dq', 'dipole'):
            s = linear_element(s, spec.get('k', 0.0), spec.get('theta', 0.0), L)
        elif kind == 'sext':
            s = rk_multipole(s, 'S', spec.get('k', 0.0), L, spec.get('n_iter', 10))
        elif kind == 'oct':
            s = rk_multipole(s, 'O', spec.get('k', 0.0), L, spec.get('n_iter', 10))
        elif kind == 'oct_std':
            s = rk_multipole(s, 'Ostd', spec.get('k', 0.0), L, spec.get('n_iter', 10))
        elif kind == 'kicker':
            s = kicker(s, spec.get('dx', 0.0), spec.get('dy', 0.0))
        else:
            raise ValueError(f"unknown kind {kind}")
        s = drift(s, d2)
        out.append(s)
    return out


def matrix_trajectory(particle, spec, n_elem=24):
    """Полилиния (X, Y, Z) частицы по матричной модели через drift(d1)+элемент+drift(d2),
    в ЛАБ-системе: фокусировка (κ) + общий ИЗГИБ (bend_y/bend_z, распределён по элементу).
    X — продольная координата 0..d1+L+d2. Для рисунка сравнения с Boris."""
    kind = spec["kind"]; d1 = spec.get("d1", 1.0); d2 = spec.get("d2", 1.0); L = spec.get("L", 0.0)
    by, bz = spec.get("bend_y", 0.0), spec.get("bend_z", 0.0)   # полные углы изгиба, рад
    st = tuple(particle)
    pts = [(0.0, st[0], st[2])]
    st = drift(st, d1); pts.append((d1, st[0], st[2]))          # вход элемента
    yen, ypen, zen, zpen = st
    if kind in ("quad", "dq", "dipole", "kicker"):
        if kind == "kicker":
            kx = ky = 0.0
        else:
            kx = spec.get("k", 0.0) + ((spec.get("theta", 0.0) ** 2) / (L * L) if L > 0 else 0.0)
            ky = -spec.get("k", 0.0)
        for j in range(1, n_elem + 1):
            s = L * j / n_elem
            yf, _ = _plane(yen, ypen, kx, s); zf, _ = _plane(zen, zpen, ky, s)
            yb = by * s * s / (2 * L) if L > 0 else 0.0         # параболический изгиб
            zb = bz * s * s / (2 * L) if L > 0 else 0.0
            pts.append((d1 + s, yf + yb, zf + zb))
        yex, ypex = _plane(yen, ypen, kx, L); zex, zpex = _plane(zen, zpen, ky, L)
        st = (yex + by * L / 2, ypex + by, zex + bz * L / 2, zpex + bz)   # +изгиб на выходе
    elif kind in ("sext", "oct"):
        kk = "S" if kind == "sext" else "O"
        for j in range(1, n_elem + 1):
            st = rk_multipole(st, kk, spec.get("k", 0.0), L / n_elem, n_iter=4)
            pts.append((d1 + L * j / n_elem, st[0], st[2]))
    st = drift(st, d2); pts.append((d1 + L + d2, st[0], st[2]))
    return pts


if __name__ == "__main__":
    print(f"Пучок: v={V_BEAM:.4e} м/с  β={BETA:.6f}  γ={GAMMA:.4f}")
    print(f"p={P_MOM:.4e} кг·м/с  Bρ={BRHO:.6f} Тл·м  q/p={QP:.4f} 1/(Тл·м)")
    # Валидация: квадруполь L=0.295, K1=2.462 (.lte QF1A) → ожидаемая фокусировка
    K1 = 2.462177415336321
    G = G_from_K(K1); k = k_from_G(G)
    print(f"\nQF1A: K1={K1:.4f} → G={G:.6e} Тл/м → k={k:.4f} 1/м²  (ожид. k=−K1={-K1:.4f})")
    L = 0.295
    parts = [(1e-3, 0.0, 0.0, 0.0), (0.0, 0.0, 1e-3, 0.0)]   # x=1мм и y=1мм пробы
    res = track_element(parts, {'kind': 'quad', 'L': L, 'd1': 0.0, 'd2': 0.0, 'k': k, 'theta': 0.0})
    print(f"  чистый квадруполь (без дрейфов):")
    print(f"   x=1мм проба: x_out={res[0][0]*1e3:+.4f}мм x'={res[0][1]*1e3:+.4f}мрад  (фокус если x' противоположен x)")
    print(f"   y=1мм проба: y_out={res[1][2]*1e3:+.4f}мм y'={res[1][3]*1e3:+.4f}мрад  (другая плоскость)")
    # дрейф 1м не меняет угол
    d = track_element([(0.0, 1e-3, 0.0, 0.0)], {'kind': 'quad', 'L': 0.0, 'd1': 1.0, 'd2': 0.0, 'k': 0.0})
    print(f"  дрейф 1м, x'=1мрад: x_out={d[0][0]*1e3:.4f}мм (ожид 1.0), x'={d[0][1]*1e3:.4f}мрад (ожид 1.0)")
