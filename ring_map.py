# -*- coding: utf-8 -*-
"""Отображения высших порядков — аналитический механизм SynchrotronSim сверх линейных матриц ring_matrix.

Отображение оборота раскладывается в ряд Тейлора по фазовым координатам до порядка N (по умолчанию 4)
методом дифференциальной алгебры: каждая величина — усечённый степенной ряд от шести переменных
z = (x, px, y, py, t, δ) ядра SynchrotronSim (px, py — канонические импульсы, t — путь). Ряды складываются,
перемножаются и подставляются друг в друга точно до порядка N, поэтому коэффициенты — это точные производные
отображения, а не разности по траекториям.

Составляющие оборота:
  СТЫК (модифицированный дрейф) — точная формула пролёта по прямой от выходной грани одного элемента до
      входной грани следующего (как modified_drift ядра), вычисленная над рядами. Грани могут быть
      непараллельны: у отображения есть свободный член, и он входит во все порядки.
  ДРЕЙФ — точная формула.
  ТЕЛО ЭЛЕМЕНТА — точное решение уравнений движения внутри элемента: z(L) = exp(L·D)·z, где D = Σ fᵢ·∂/∂zᵢ —
      производная Ли по векторному полю элемента. Поле то же, что интегрирует ядро (секторный изгиб — порт
      CSBEND: точный криволинейный дрейф + поле с множителем 1 + x/ρ; прямой мультиполь — точный дрейф +
      мультипольная сила), но решение точное, без срезов интегратора. Экспонента берётся от матрицы D
      в базисе мономов (scipy.linalg.expm).
  КРАЯ ИЗГИБА — тонкие линзы px += k·x, py −= k·y, как в ядре.

Из отображения оборота:
  closed_orbit_series — замкнутая орбита как ряд по δ;
  normal_form         — сдвиг частот с действием ∂ν/∂J (нормальная форма отображения);
  chromaticity        — ν(δ) = ν0 + ξ1·δ + ξ2·δ²/2 + ξ3·δ³/6: по линеаризации отображения вокруг орбиты
                        с данным δ; точны производные до порядка N − 1;
  momentum_compaction — ΔC/C = αc·δ + αc2·δ² + …

Пример:
  py -3.14 ring_map.py --lattice 2 --lead 1.4335 --mev 3000
"""
import math, argparse, itertools
import numpy as np
from scipy.linalg import expm


# ============================================================================ алгебра усечённых рядов
class DA:
    """Усечённые степенные ряды nv переменных до порядка N. Ряд — вектор коэффициентов по мономам."""

    def __init__(self, nv=6, order=4):
        self.nv, self.N = nv, order
        monos = [e for deg in range(order + 1)
                 for e in sorted((e for e in itertools.product(range(deg + 1), repeat=nv) if sum(e) == deg),
                                 reverse=True)]
        self.monos = monos
        self.n = len(monos)
        self.index = {e: i for i, e in enumerate(monos)}
        self.deg = np.array([sum(e) for e in monos])
        I, J, K = [], [], []
        for i, a in enumerate(monos):
            for j, b in enumerate(monos):
                if self.deg[i] + self.deg[j] <= order:
                    I.append(i); J.append(j); K.append(self.index[tuple(x + y for x, y in zip(a, b))])
        self.I, self.J, self.K = np.array(I), np.array(J), np.array(K)
        self._dmat = []                                   # ∂/∂z_v матрицей
        for v in range(nv):
            M = np.zeros((self.n, self.n))
            for i, e in enumerate(monos):
                if e[v]:
                    t = list(e); t[v] -= 1
                    M[self.index[tuple(t)], i] = e[v]
            self._dmat.append(M)
        # родитель монома: моном = родитель · z_var (для подстановки рядов)
        self.parent = np.zeros(self.n, int); self.pvar = np.zeros(self.n, int)
        for i, e in enumerate(monos[1:], 1):
            v = next(k for k in range(nv) if e[k])
            t = list(e); t[v] -= 1
            self.parent[i] = self.index[tuple(t)]; self.pvar[i] = v

    # ---- создание
    def const(self, c):
        a = np.zeros(self.n); a[0] = c
        return a

    def var(self, v, c=0.0):
        a = np.zeros(self.n); a[0] = c; a[1 + v] = 1.0   # мономы первой степени идут после константы
        assert self.monos[1 + v][v] == 1
        return a

    def identity(self):
        return [self.var(v) for v in range(self.nv)]

    # ---- арифметика
    def mul(self, a, b):
        return np.bincount(self.K, weights=a[self.I] * b[self.J], minlength=self.n)

    def mulmat(self, f):
        """Матрица умножения на ряд f: mulmat(f) @ b = f·b."""
        M = np.zeros((self.n, self.n))
        np.add.at(M, (self.K, self.J), f[self.I])
        return M

    def deriv(self, a, v):
        return self._dmat[v] @ a

    def dmat(self, v):
        return self._dmat[v]

    def _series(self, a, derivs):
        """f(a) по ряду Тейлора вокруг свободного члена: derivs[k] = f⁽ᵏ⁾(a₀)."""
        eps = a.copy(); eps[0] = 0.0
        out = self.const(derivs[0]); p = self.const(1.0)
        for k in range(1, self.N + 1):
            p = self.mul(p, eps)
            out = out + derivs[k] / math.factorial(k) * p
        return out

    def inv(self, a):
        a0 = a[0]
        return self._series(a, [(-1) ** k * math.factorial(k) / a0 ** (k + 1) for k in range(self.N + 1)])

    def div(self, a, b):
        return self.mul(a, self.inv(b))

    def sqrt(self, a):
        a0 = a[0]
        if a0 <= 0:
            raise ValueError("корень из неположительного ряда")
        d = [math.sqrt(a0)]
        c = 0.5
        for k in range(1, self.N + 1):                  # (x^½)⁽ᵏ⁾ = ½·(½−1)·…·(½−k+1)·x^(½−k)
            d.append(math.prod(0.5 - i for i in range(k)) * a0 ** (0.5 - k))
        return self._series(a, d)

    # ---- подстановка
    def monomials(self, state):
        """Все мономы от рядов state (подставленных переменных)."""
        out = np.empty((self.n, self.n))
        out[0] = self.const(1.0)
        for i in range(1, self.n):
            out[i] = self.mul(out[self.parent[i]], state[self.pvar[i]])
        return out

    def compose(self, coefs, state):
        """Ряды coefs (список векторов коэффициентов) с переменными, замещёнными рядами state."""
        mono = self.monomials(state)
        return [c @ mono for c in coefs]

    def coef(self, a, exps):
        return a[self.index[tuple(exps)]]


# ============================================================================ составляющие оборота
def drift(da, st, L):
    """Точный дрейф ядра (δ-корректный): x += L·px/ps, t += L·(1+δ)/ps."""
    x, px, y, py, t, d = st
    one = da.const(1.0) + d
    ips = da.inv(da.sqrt(da.mul(one, one) - da.mul(px, px) - da.mul(py, py)))
    return [x + L * da.mul(px, ips), px, y + L * da.mul(py, ips), py, t + L * da.mul(one, ips), d]


def junction(da, st, A, B, same_plane=False, para_tol=1e-4):
    """Модифицированный дрейф ядра над рядами: выход рамки A → вход рамки B (рамка — (p, axis, u, v))."""
    aP, aAx, aU, aV = (np.asarray(v, float) for v in A)
    bP, bAx, bU, bV = (np.asarray(v, float) for v in B)
    x, px, y, py, t, d = st
    one = da.const(1.0) + d
    ps = da.sqrt(da.mul(one, one) - da.mul(px, px) - da.mul(py, py))
    iop = da.inv(one)
    P = [da.const(aP[k]) + aU[k] * x + aV[k] * y for k in range(3)]
    dirv = [da.mul(aU[k] * px + aV[k] * py + aAx[k] * ps, iop) for k in range(3)]
    dot = lambda vec, w: sum(w[k] * vec[k] for k in range(3))
    if same_plane:
        Q = P
    elif np.linalg.norm(np.cross(aAx, bAx)) < para_tol:          # оси параллельны: плоскость ⊥ aAx через bP
        num = dot([da.const(bP[k]) - P[k] for k in range(3)], aAx)
        tt = da.div(num, dot(dirv, aAx))
        Q = [P[k] + da.mul(dirv[k], tt) for k in range(3)]
    else:                                                        # клин: пролёт до плоскости грани B
        num = dot([da.const(bP[k]) - P[k] for k in range(3)], bAx)
        tt = da.div(num, dot(dirv, bAx))
        Q = [P[k] + da.mul(dirv[k], tt) for k in range(3)]
    path = sum(da.mul(Q[k] - P[k], dirv[k]) for k in range(3))
    rel = [Q[k] - da.const(bP[k]) for k in range(3)]
    return [dot(rel, bU), da.mul(one, dot(dirv, bU)), dot(rel, bV), da.mul(one, dot(dirv, bV)), t + path, d]


def edge(st, k):
    """Тонкая краевая линза грани изгиба (как в ядре): px += k·x, py −= k·y."""
    x, px, y, py, t, d = st
    return [x, px + k * x, y, py - k * y, t, d]


def _csbend_fields(da, X, Y, b1, b2, b3, h):
    """Нормированное поле комбинированной функции (Fx, Fy) над рядами — csbend_fields ядра при b4 = 0."""
    h2, h3 = h * h, h ** 3
    x2, y2 = da.mul(X, X), da.mul(Y, Y)
    x3, xy2 = da.mul(x2, X), da.mul(X, y2)
    Fy = (da.const(1.0) + b1 * X + 0.5 * b2 * x2 + (b3 / 6.0) * x3
          + (-(h * b1) / 2.0 - b2 / 2.0) * y2
          + ((h2 * b1) / 2.0 - (h * b2) / 2.0 - b3 / 2.0) * xy2
          + (-(h3 * b1) / 2.0 + (h2 * b2) / 2.0 - (h * b3) / 4.0) * da.mul(x2, y2)
          + ((h3 * b1) / 24.0 - (h2 * b2) / 24.0 + (h * b3) / 12.0) * da.mul(y2, y2))
    Fx = (b1 * Y + b2 * da.mul(X, Y) + (b3 / 2.0) * da.mul(x2, Y)
          + ((h2 * b1) / 6.0 - (h * b2) / 6.0 - b3 / 6.0) * da.mul(y2, Y)
          + (-(h3 * b1) / 3.0 + (h2 * b2) / 3.0 - (h * b3) / 6.0) * da.mul(xy2, Y))
    return Fx, Fy


def field_sector(da, h, b1, b2, b3):
    """Векторное поле секторного тела ядра (порт CSBEND): точный криволинейный дрейф + кик поля."""
    X, QX, Y, QY, T, D = da.identity()
    one = da.const(1.0) + D
    ps = da.sqrt(da.mul(one, one) - da.mul(QX, QX) - da.mul(QY, QY))
    ips = da.inv(ps)
    hx = da.const(1.0) + h * X
    Fx, Fy = _csbend_fields(da, X, Y, b1, b2, b3, h)
    return [da.mul(hx, da.mul(QX, ips)), h * ps - h * da.mul(hx, Fy), da.mul(hx, da.mul(QY, ips)),
            h * da.mul(hx, Fx), da.mul(hx, da.mul(one, ips)), da.const(0.0)]


def field_multipole(da, k1, k1v, k2, k3):
    """Векторное поле прямого мультиполя ядра (KQUAD/KSEXT/KOCT): точный дрейф + сила."""
    X, QX, Y, QY, T, D = da.identity()
    one = da.const(1.0) + D
    ips = da.inv(da.sqrt(da.mul(one, one) - da.mul(QX, QX) - da.mul(QY, QY)))
    x2, y2 = da.mul(X, X), da.mul(Y, Y)
    fpx = -k1 * X - 0.5 * k2 * (x2 - y2) - (k3 / 6.0) * (da.mul(x2, X) - 3.0 * da.mul(X, y2))
    fpy = k1v * Y + k2 * da.mul(X, Y) + (k3 / 6.0) * (3.0 * da.mul(x2, Y) - da.mul(y2, Y))
    return [da.mul(QX, ips), fpx, da.mul(QY, ips), fpy, da.mul(one, ips), da.const(0.0)]


def flow(da, f, L):
    """Точное решение z′ = f(z) на длине L: коэффициенты рядов z(L) = exp(L·D)·z, D = Σ fᵢ·∂ᵢ."""
    D = sum(da.mulmat(f[v]) @ da.dmat(v) for v in range(da.nv))
    E = expm(L * D)
    return [E @ da.var(v) for v in range(da.nv)]


# ============================================================================ оборот
class RingMap:
    """Отображение оборота кольца до порядка N по описанию магнитов (атрибуты mags/plane кольца)."""

    def __init__(self, ring, order=4):
        self.ring, self.da = ring, DA(6, order)
        self._body = {}

    def _body_map(self, m):
        """Ряды тела магнита (кэш по параметрам: одинаковые магниты считаются один раз)."""
        L = float(m.Leff)
        if getattr(m, "sector", False):
            rho0 = -L / m.theta0; h = 1.0 / rho0
            key = ("s", L, h, m.k1, m.k2, m.k3)
            if key not in self._body:
                self._body[key] = flow(self.da, field_sector(self.da, h, m.k1 * rho0, m.k2 * rho0, m.k3 * rho0), L)
            return self._body[key]
        key = ("k", L, m.k1, m.k1v, m.k2, m.k3)
        if key not in self._body:
            self._body[key] = flow(self.da, field_multipole(self.da, m.k1, m.k1v, m.k2, m.k3), L)
        return self._body[key]

    def element(self, st, m):
        da = self.da
        if getattr(m, "sector", False):
            st = edge(st, m.kEdge1)
            st = da.compose(self._body_map(m), st)
            return edge(st, m.kEdge2)
        if getattr(m, "kick_sector", False):
            return da.compose(self._body_map(m), st)
        return drift(da, st, float(m.Leff))

    def one_turn(self, state=None):
        """Ряды отображения оборота: плоскость → стык → магнит → … → стык → плоскость."""
        r, da = self.ring, self.da
        st = state or da.identity()
        prev = r.plane
        for m in r.mags:
            st = junction(da, st, prev.fr_out, m.fr_in, getattr(m, "same_plane_as_prev", False))
            st = self.element(st, m)
            prev = m
        return junction(da, st, prev.fr_out, r.plane.fr_in, False)

    # ---- производные показатели
    def closed_orbit_series(self, M):
        """Замкнутая орбита (x, px, y, py) как ряд по δ: z = M(z, δ) при t = 0."""
        da = self.da
        R = np.array([[da.coef(M[i], np.eye(6, dtype=int)[j]) for j in range(4)] for i in range(4)])
        inv = np.linalg.inv(np.eye(4) - R)
        dvar = da.var(5)
        z = [da.const(0.0) for _ in range(4)]
        for _ in range(da.N + 2):
            full = z + [da.const(0.0), dvar]
            Mz = da.compose(M[:4], full)
            rhs = [Mz[i] - sum(R[i, j] * z[j] for j in range(4)) for i in range(4)]
            z = [sum(inv[i, j] * rhs[j] for j in range(4)) for i in range(4)]
        return z + [da.const(0.0), dvar]

    def _delta_coefs(self, a):
        return np.array([self.da.coef(a, (0, 0, 0, 0, 0, k)) for k in range(self.da.N + 1)])

    def chromaticity(self, M=None, nu_int=(0, 0)):
        """ν(δ) из линеаризации отображения вокруг орбиты с δ. Возврат dict: nux, nuy, xi_x, xi_y
        (производные dᵏν/dδᵏ, k = 0…N−1; точны до N−1), orbit (ряды замкнутой орбиты)."""
        da = self.da
        M = M or self.one_turn()
        co = self.closed_orbit_series(M)
        A = [[da.compose([da.deriv(M[i], j)], co)[0] for j in range(4)] for i in range(4)]
        out = dict(orbit=co)
        for pl, (i, k) in (("x", (0, 1)), ("y", (2, 3))):
            c = 0.5 * (A[i][i] + A[k][k])                     # cos μ(δ) — ряд по δ
            sgn = 1.0 if A[i][k][0] > 0 else -1.0
            s = sgn * da.sqrt(da.const(1.0) - da.mul(c, c))    # sin μ(δ)
            mu0 = math.atan2(s[0], c[0]) % (2 * math.pi)
            dmu = -da.mul(da.deriv(c, 5), da.inv(s))          # dμ/dδ
            cd = self._delta_coefs(dmu)
            mu = np.zeros(da.N + 1); mu[0] = mu0
            for n in range(da.N):
                mu[n + 1] = cd[n] / (n + 1)
            nu = mu / (2 * math.pi)
            nu[0] += nu_int[0 if pl == "x" else 1]
            fac = np.array([math.factorial(n) for n in range(da.N + 1)])
            out["nu" + pl] = nu[0]
            out["xi_" + pl] = (nu * fac)[: da.N]              # ν, ξ1, ξ2, ξ3 (при N = 4)
        return out

    def momentum_compaction(self, M=None, co=None):
        """Путь по замкнутой орбите: C(δ) = C·(1 + αc·δ + αc2·δ² + …). Возврат (C, [αc, αc2, …])."""
        da = self.da
        M = M or self.one_turn()
        co = co or self.closed_orbit_series(M)
        P = self._delta_coefs(da.compose([M[4]], co)[0])
        return P[0], P[1:] / P[0]


# ============================================================================ нормальная форма
S4 = np.array([[0.0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]])   # симплектическая форма 4×4


class CDA:
    """Комплексные ряды над той же алгеброй мономов (умножение — через вещественные части)."""

    def __init__(self, da):
        self.da = da
        self.n = da.n

    def zero(self):
        return np.zeros(self.n, complex)

    def var(self, v):
        return self.da.var(v).astype(complex)

    def identity(self):
        return [self.var(v) for v in range(self.da.nv)]

    def mul(self, a, b):
        m = self.da.mul
        ar, ai, br, bi = a.real, a.imag, b.real, b.imag
        return (m(ar, br) - m(ai, bi)) + 1j * (m(ar, bi) + m(ai, br))

    def monomials(self, state):
        out = np.empty((self.n, self.n), complex)
        out[0] = 0.0; out[0][0] = 1.0
        for i in range(1, self.n):
            out[i] = self.mul(out[self.da.parent[i]], state[self.da.pvar[i]])
        return out

    def compose(self, coefs, state):
        mono = self.monomials(state)
        return [c @ mono for c in coefs]

    def coef(self, a, exps):
        return a[self.da.index[tuple(exps)]]


def _map4(rmap, M=None):
    """Поперечное отображение 4×4 вокруг замкнутой орбиты при δ = 0 (t и δ выброшены)."""
    da6 = rmap.da
    M = M or rmap.one_turn()
    co = rmap.closed_orbit_series(M)
    z0 = [float(c[0]) for c in co[:4]]
    st = [da6.var(i, z0[i]) for i in range(4)] + [da6.const(0.0), da6.const(0.0)]
    Ms = da6.compose(M[:4], st)
    da = DA(4, da6.N)
    idx = {}                                         # моном 6 переменных без t и δ → моном 4 переменных
    for i, e in enumerate(da6.monos):
        if e[4] == 0 and e[5] == 0:
            idx[i] = da.index[e[:4]]
    out = []
    for c in Ms:
        v = np.zeros(da.n)
        for i, j in idx.items():
            v[j] += c[i]
        v[0] = 0.0                                   # свободный член — сама орбита
        out.append(v)
    return da, out


def _eigen_basis(R):
    """Комплексный базис, в котором линейная часть — поворот: A (столбцы v₁, v̄₁, v₂, v̄₂) и λ."""
    lam, V = np.linalg.eig(R)
    picked = []
    used = np.zeros(len(lam), bool)
    for i in range(len(lam)):
        if used[i] or lam[i].imag <= 0:
            continue
        v = V[:, i]
        plane = 0 if (abs(v[0]) ** 2 + abs(v[1]) ** 2) >= (abs(v[2]) ** 2 + abs(v[3]) ** 2) else 1
        picked.append((plane, lam[i], v)); used[i] = True
    if len(picked) != 2:
        raise RuntimeError("не удалось разделить плоскости: движение связано или неустойчиво")
    picked.sort(key=lambda t: t[0])
    cols, lams = [], []
    for _, l, v in picked:
        n = np.conj(v) @ S4 @ v                      # чисто мнимое; нормируем на −2i
        v = v / math.sqrt(abs(n) / 2.0)
        if (np.conj(v) @ S4 @ v).imag > 0:           # знак: v̄ᵀ S v = −2i
            v = np.conj(v); l = np.conj(l)
        cols += [v, np.conj(v)]; lams += [l, np.conj(l)]
    return np.array(cols).T, np.array(lams)


def _invert_near_identity(cda, t, order):
    """Обратное к T = I + t (t степени ≥ 2): U = I − t∘U, сходится за число шагов порядка."""
    u = cda.identity()
    for _ in range(order):
        tu = cda.compose(t, u)
        u = [cda.var(i) - tu[i] for i in range(4)]
    return u


def normal_form(ring, order=4, rmap=None):
    """Нормальная форма отображения оборота: сдвиг частот с действием ∂ν/∂J.

    Отображение вокруг замкнутой орбиты переводится в комплексный базис, где линейная часть — поворот
    e^{iμ}, после чего нерезонансные члены убираются заменой переменных степень за степенью. Остаются
    только резонансные: ζ₁ → e^{iμx}·ζ₁·(1 + i·axx·|ζ₁|² + i·axy·|ζ₂|²). Действие J = |ζ|², то есть
    J = (q² + (αq + βp)²)/(2β) — как у elegant (ε = 2J).

    Возврат dict: nux, nuy (дробные), dnux/dJx, dnux/dJy, dnuy/dJx, dnuy/dJy (1/м)."""
    r = rmap or RingMap(ring, order)
    da, M4 = _map4(r)
    return normal_form_of_map(da, M4)


def normal_form_of_map(da, M4):
    """Нормальная форма готового 4-мерного отображения (см. normal_form)."""
    cda = CDA(da)
    unit = np.eye(4, dtype=int)
    R = np.array([[da.coef(M4[i], unit[j]) for j in range(4)] for i in range(4)])
    A, lam = _eigen_basis(R)
    Ainv = np.linalg.inv(A)
    zeta = cda.identity()
    st = [sum(A[i, j] * zeta[j] for j in range(4)) for i in range(4)]        # z = A·ζ
    N = cda.compose([c.astype(complex) for c in M4], st)
    N = [sum(Ainv[i, j] * N[j] for j in range(4)) for i in range(4)]         # ζ' = A⁻¹·M(A·ζ)
    for k in range(2, da.N + 1):                     # убрать нерезонансные члены степени k
        t = [cda.zero() for _ in range(4)]
        hit = False
        for i in range(4):
            for m, e in enumerate(da.monos):
                if da.deg[m] != k or N[i][m] == 0:
                    continue
                lp = np.prod(lam ** np.array(e))
                if abs(lp - lam[i]) < 1e-6:          # резонансный член: остаётся в нормальной форме
                    continue
                t[i][m] = N[i][m] / (lp - lam[i]); hit = True
        if not hit:
            continue
        T = [cda.var(i) + t[i] for i in range(4)]
        Tinv = _invert_near_identity(cda, t, da.N)
        N = cda.compose(cda.compose(Tinv, N), T)     # T⁻¹ ∘ N ∘ T
    # ветвь частоты: у собственного вектора с нормировкой v̄ᵀSv = −2i фаза может идти в обратную сторону,
    # поэтому знак выбирается по частоте Куранта–Снайдера того же блока 2×2
    sgn, nu = {}, {}
    for pl, i in (("x", 0), ("y", 2)):
        m = R[i:i + 2, i:i + 2]; cmu = 0.5 * float(np.trace(m))
        smu = math.copysign(math.sqrt(max(0.0, 1.0 - cmu * cmu)), m[0, 1])
        nu_cs = math.atan2(smu, cmu) % (2 * math.pi) / (2 * math.pi)
        ph = math.atan2(lam[i].imag, lam[i].real) / (2 * math.pi)
        sgn[pl] = 1.0 if abs((ph % 1.0) - nu_cs) < 1e-6 else -1.0
        nu[pl] = (sgn[pl] * ph) % 1.0
    out = dict(nux=nu["x"], nuy=nu["y"])
    # действие Куранта–Снайдера: при нормировке v̄ᵀSv = −2i получается |ζ|² = J/2 (ε = 2J, как у elegant)
    terms = {("x", "x"): (0, (2, 1, 0, 0)), ("x", "y"): (0, (1, 0, 1, 1)),
             ("y", "x"): (2, (1, 1, 1, 0)), ("y", "y"): (2, (0, 0, 2, 1))}
    for (pl, act), (i, e) in terms.items():
        c = cda.coef(N[i], e) / lam[i]
        out[f"dnu{pl}/dJ{act}"] = sgn[pl] * float(c.imag) / (4 * math.pi)
        out[f"|re|dnu{pl}/dJ{act}"] = float(abs(c.real)) / (4 * math.pi)   # контроль: должен быть мал
    return out


def main():
    import os
    os.environ.setdefault("SS_SECTOR", "1")
    import ring_analysis as ra
    ap = argparse.ArgumentParser(description="хроматичность и уплотнение орбит по отображению высших порядков")
    ap.add_argument("--lattice", type=int, default=1)
    ap.add_argument("--lead", type=float, default=None)
    ap.add_argument("--mev", type=float, default=6000.0)
    ap.add_argument("--order", type=int, default=4)
    a = ap.parse_args()
    ring = ra.Ring(a.lattice, start_before_first_m=a.lead, design_mev=a.mev)
    rm = RingMap(ring, a.order)
    M = rm.one_turn()
    ch = rm.chromaticity(M)
    C, alph = rm.momentum_compaction(M, ch["orbit"])
    print(f"C = {C:.9f} м; αc = {alph[0]:.6e}, αc2 = {alph[1]:.6e}")
    for pl in ("x", "y"):
        print(f"ν{pl} (дробная) = {ch['nu' + pl]:.8f};  ξ1..ξ{a.order - 1} = "
              + ", ".join(f"{v:+.6g}" for v in ch["xi_" + pl][1:]))


if __name__ == "__main__":
    main()
