# -*- coding: utf-8 -*-
"""Диагностика кольца трекингом — одни и те же алгоритмы поверх двух трекеров: SynchrotronSim и elegant.

Все величины получаются ИЗ ТРЕКИНГА, без аналитических формул оптики: линейные функции — из матриц, снятых
центральными разностями по траекториям, нелинейные — из пооборотных координат. Трекер подставляется
объектом с общим интерфейсом, поэтому различие результатов означает различие самих трекеров:

  ring_cpp.CppRing (SynchrotronSim, продакшн-ядро C++ через MatrixTracker --ring-serve) — рабочий трекер
  Ring (SynchrotronSim, Python-эталон ядра matrix_tracker, SS_SECTOR=1) — только для сверки C++ с эталоном
  elegant_diag.ElegantRing (elegant, та же решётка)

Интерфейс трекера:
  track(states, n_turns=1)      → (n, 6) после n_turns, NaN у потерянных
  track_elements(states)        → (n, n_mag+1, 6) на выходе каждого магнита за один оборот + у плоскости
  s_out, circumference, names   — опорная s выхода магнитов, периметр, имена

Координаты — канонические SynchrotronSim: x, px=(1+δ)·x'/√(1+x'²+y'²), y, py, t, δ.
t в track/track_tbt — запаздывание относительно света за целые обороты, в метрах: Σ(путь/β(δ)) − n·C
(β — по проектной энергии, как время пролёта в elegant); в track_elements — путь от плоскости, м.

ВЧ-резонатор (rf=dict(volt_mv, harmonic, phase_deg)) стоит у плоскости наблюдения в конце оборота и
повторяет RFCA elegant с FIDUCIAL="light", FREQ = h·c/C: ΔE = e·V·sin(φ + 2πh·t/C). Излучение (sr, isr) —
распределённое в изгибах ядра (SS_SR / SS_ISR), по проектной энергии design_mev.

Показатели:
  closed_orbit(ring, δ)          замкнутая орбита — Ньютон по отображению за оборот
  transfer_matrices(ring, δ)     матрицы 6×6 от точки наблюдения до выхода каждого магнита и за оборот
  twiss(ring, δ)                 β, α, фазовый набег, η, η', частоты, αc
  chromaticity(ring, deltas)     ν(δ) и ряд Тейлора ν0 + ξ1·δ + ξ2·δ²/2 + ξ3·δ³/6
  element_matrices(ring)         матрица 6×6 каждого магнита (с дрейфом перед ним)

Пример:
  py -3.14 ring_analysis.py twiss
"""
import os, io, math, argparse, contextlib
import numpy as np

os.environ.setdefault("SS_SECTOR", "1")                 # секторная модель обязательна (сверка с elegant)
import matrix_tracker as mt
from beam_injection import PlaneMagnet


# ============================================================================ трекер SynchrotronSim
class Ring:
    """Замкнутое кольцо решётки БД с плоскостью наблюдения перед первым магнитом линии.

    start_before_first_m — путь от плоскости до входа первого магнита по опорной орбите (None ⇒ половина
    замыкающего промежутка; у СИЛА промежуток симметричен, и это совпадает с s=0 elegant)."""

    label = "SynchrotronSim"
    MC2_MEV = 0.51099895069

    def __init__(self, lattice=1, start_before_first_m=None, settings=None, sr=False, isr=False, rf=None,
                 design_mev=6000.0):
        eff, phys, coup, trans = mt.load_db()
        env = {"SS_SR": "1" if (sr or isr) else "", "SS_ISR": "1" if isr else "", "SS_DESIGN_MEV": str(design_mev)}
        saved = {k: os.environ.get(k) for k in env}
        os.environ.update(env)                               # коэффициенты излучения считаются при сборке
        try:
            with contextlib.redirect_stdout(io.StringIO()):  # сборка печатает предупреждения стыков
                self.mags = mt.build_line(lattice, eff, phys, coup, trans, settings=settings)
        finally:
            for k, v in saved.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
        self.sr, self.isr, self.rf, self.design_mev = bool(sr or isr), bool(isr), rf, float(design_mev)
        m0, m1 = self.mags[0], self.mags[-1]
        gap = float((m0.fr_in[0] - m1.fr_out[0]) @ m1.fr_out[1])
        d0 = gap / 2.0 if start_before_first_m is None else float(start_before_first_m)
        self.plane = PlaneMagnet(m0.fr_in[0] - d0 * m0.fr_in[1], m0.fr_in[1])
        self.lattice = lattice
        s = d0; s_in = []; s_out = []
        for j, m in enumerate(self.mags):
            if j:
                s += float((m.fr_in[0] - self.mags[j - 1].fr_out[0]) @ self.mags[j - 1].fr_out[1])
            s_in.append(s); s += m.L; s_out.append(s)
        self.s_in = np.array(s_in); self.s_out = np.array(s_out)
        self.circumference = s + float((self.plane.fr_in[0] - m1.fr_out[0]) @ m1.fr_out[1])
        self.names = [m.dev for m in self.mags]
        self.kinds = [m.kind for m in self.mags]

    def beta(self, delta):
        """Скорость частицы (v/c) при импульсе (1+δ)·p0 проектной энергии."""
        pc = self.design_mev * (1.0 + delta)
        return pc / math.sqrt(pc * pc + self.MC2_MEV ** 2)

    def turn(self, st):
        """Один оборот одной частицы: плоскость → все магниты → плоскость (+ ВЧ). Бросает mt.ParticleLost."""
        mags = self.mags
        t0 = float(st[4])
        st = mt.modified_drift(tuple(st), self.plane, mags[0]); st = mags[0].body(st)
        for j in range(1, len(mags)):
            st = mt.modified_drift(st, mags[j - 1], mags[j]); st = mags[j].body(st)
        x, px, y, py, t, d = mt.modified_drift(st, mags[-1], self.plane)
        t = t0 + (t - t0) / self.beta(d) - self.circumference   # запаздывание относительно света, м
        if self.rf:
            d = self.rf_kick(t, d)
        return (x, px, y, py, t, d)

    def rf_kick(self, t, d):
        """RFCA: ΔE = e·V·sin(φ + 2πh·t/C); канонические поперечные импульсы не меняются."""
        rf = self.rf
        dE = rf["volt_mv"] * math.sin(math.radians(rf["phase_deg"]) + 2 * math.pi * rf["harmonic"] * t
                                      / self.circumference)
        pc = self.design_mev * (1.0 + d)
        E = math.sqrt(pc * pc + self.MC2_MEV ** 2) + dE
        return math.sqrt(E * E - self.MC2_MEV ** 2) / self.design_mev - 1.0

    def track(self, states, n_turns=1):
        out = np.full((len(states), 6), np.nan)
        for p, st in enumerate(states):
            cur = tuple(map(float, st))
            try:
                for _ in range(n_turns):
                    cur = self.turn(cur)
                out[p] = cur
            except (mt.ParticleLost, ValueError, OverflowError):
                pass
        return out

    def track_elements(self, states):
        mags = self.mags; n = len(mags)
        out = np.full((len(states), n + 1, 6), np.nan)
        for p, st in enumerate(states):
            try:
                st = mt.modified_drift(tuple(map(float, st)), self.plane, mags[0]); st = mags[0].body(st)
                out[p, 0] = st
                for j in range(1, n):
                    st = mt.modified_drift(st, mags[j - 1], mags[j]); st = mags[j].body(st); out[p, j] = st
                out[p, n] = mt.modified_drift(st, mags[-1], self.plane)
            except (mt.ParticleLost, ValueError, OverflowError):
                pass
        return out

    def track_tbt(self, states, n_turns, limit=None):
        """Пооборотные координаты у плоскости: (n, n_turns+1, 6), NaN после потери.
        limit — порог |x| или |y| в точке наблюдения, м: превысившая частица считается потерянной."""
        out = np.full((len(states), n_turns + 1, 6), np.nan)
        for p, st in enumerate(states):
            cur = tuple(map(float, st)); out[p, 0] = cur
            try:
                for k in range(n_turns):
                    cur = self.turn(cur)
                    if not all(map(math.isfinite, cur)):
                        break
                    if limit is not None and (abs(cur[0]) > limit or abs(cur[2]) > limit):
                        break
                    out[p, k + 1] = cur
            except (mt.ParticleLost, ValueError, OverflowError):
                pass
        return out


def yaw_element(m, alpha):
    """Развернуть ПРЯМОЙ элемент на угол alpha вокруг его вертикали (вход остаётся на месте) — имитация
    рассогласования установки. Грани элемента наклоняются к опорной орбите, и стыки становятся клиньями."""
    if getattr(m, "sector", False):
        raise ValueError("поворот рамок описан для прямых элементов; у изгиба грани лежат на дуге")
    P, ax, u, v = m.fr_in
    ax2 = mt._rot(ax, v, alpha); u2 = mt._rot(u, v, alpha)
    m.fr_in = (P, ax2, u2, v)
    m.fr_out = (P + m.L * ax2, ax2, u2, v)
    return m


# ============================================================================ алгоритмы (общие для трекеров)
def closed_orbit(ring, delta=0.0, tol=1e-13, eps=1e-8, max_iter=10, dims=4):
    """Замкнутая орбита: Ньютон по T(z)−z, якобиан — разностями.
    dims=4 — поперечная (x,px,y,py) при заданном δ (без ВЧ); dims=6 — неподвижная точка 6D (с ВЧ и излучением:
    синхронная частица с её δ_s и фазой)."""
    z = np.zeros(dims)
    if dims == 6:
        z[5] = delta
    full = (lambda v: np.r_[v, 0.0, delta]) if dims == 4 else (lambda v: np.asarray(v, float))
    for _ in range(max_iter):
        pts = [full(z)]
        for k in range(dims):
            q = z.copy(); q[k] += eps; pts.append(full(q))
        fin = ring.track(pts)[:, [0, 1, 2, 3] if dims == 4 else slice(0, 6)]
        F = fin[0] - z
        if np.max(np.abs(F)) < tol:
            break
        J = ((fin[1:] - fin[0]) / eps).T
        z = z - np.linalg.solve(J - np.eye(dims), F)
    return full(z)


def _probes(orbit, eps):
    pts = [orbit.copy()]
    for k in range(6):
        for sgn in (1, -1):
            p = orbit.copy(); p[k] += sgn * eps; pts.append(p)
    return pts


# Разностная схема. Орбита за оборот набирает ошибку округления ~1e-12 м (тысячи элементов), поэтому шаг
# 1e-6 даёт шум ~1e-6 в матрице — столько же, сколько радиационное затухание за оборот. Схема 4-го порядка
# с шагом 1e-5: шум округления ~1e-7, погрешность от нелинейности ~h⁴ — пренебрежимо.
_STENCIL = ((2, -1.0 / 12), (1, 8.0 / 12), (-1, -8.0 / 12), (-2, 1.0 / 12))
FD_EPS = 1e-5


def _fd_columns(orbit, cols, eps):
    """Пробные частицы 4-точечной схемы: орбита, затем по 4 на столбец."""
    pts = [np.asarray(orbit, float).copy()]
    for k in cols:
        for m, _ in _STENCIL:
            p = np.asarray(orbit, float).copy(); p[k] += m * eps[k]; pts.append(p)
    return pts


def _fd_collect(tr, cols, eps, axis_shape):
    D = np.zeros(axis_shape)
    for i, k in enumerate(cols):
        acc = 0.0
        for j, (_, w) in enumerate(_STENCIL):
            acc = acc + w * tr[1 + 4 * i + j]
        D[..., k] = acc / eps[k]
    return D


def transfer_matrices(ring, delta=0.0, orbit=None, eps=FD_EPS):
    """(orbit, R, orbit_s): R формы (n_mag+1, 6, 6) — от плоскости до выхода магнита j; последняя — оборот.
    Столбец t (путь) не влияет на поперечное движение без ВЧ ⇒ его производная берётся как единичная."""
    if orbit is None:
        orbit = closed_orbit(ring, delta)
    cols = [0, 1, 2, 3, 5]
    eps = np.broadcast_to(np.asarray(eps, float), (6,))
    tr = ring.track_elements(_fd_columns(orbit, cols, eps))
    R = _fd_collect(tr, cols, eps, (tr.shape[1], 6, 6))
    R[:, 4, 4] = 1.0
    return orbit, R, tr[0]


def twiss_from_matrices(R, s, circumference):
    """β, α, ψ, η по накопленным матрицам R (последняя — оборот). Плоскости не связаны (скошенных нет)."""
    M = R[-1]
    out = dict(s=np.asarray(s), M=M)
    for pl, i in (("x", 0), ("y", 2)):
        m = M[i:i + 2, i:i + 2]
        cmu = 0.5 * np.trace(m)
        if abs(cmu) >= 1:
            raise RuntimeError(f"плоскость {pl}: движение неустойчиво (cos μ = {cmu:.6f})")
        smu = math.copysign(math.sqrt(1 - cmu * cmu), m[0, 1])
        b0 = m[0, 1] / smu; a0 = (m[0, 0] - m[1, 1]) / (2 * smu)
        r11 = R[:, i, i]; r12 = R[:, i, i + 1]; r21 = R[:, i + 1, i]; r22 = R[:, i + 1, i + 1]
        u = r11 * b0 - r12 * a0; v = r21 * b0 - r22 * a0
        out[f"beta{pl}"] = (u * u + r12 * r12) / b0
        out[f"alpha{pl}"] = -(u * v + r12 * r22) / b0
        ph = np.unwrap(np.arctan2(r12, u))              # набег за магнит ≪ π ⇒ развёртка однозначна
        out[f"psi{pl}"] = ph
        out[f"nu{pl}"] = ph[-1] / (2 * np.pi)
    D0 = np.linalg.solve(np.eye(4) - M[:4, :4], M[:4, 5])
    D = np.einsum("nij,j->ni", R[:, :4, :4], D0) + R[:, :4, 5]
    out["etax"], out["etaxp"], out["etay"], out["etayp"] = D[:, 0], D[:, 1], D[:, 2], D[:, 3]
    # уплотнение орбит: удлинение пути за оборот на единицу δ по дисперсионной орбите
    out["alphac"] = float(M[4, :4] @ D0 + M[4, 5]) / circumference
    return out


def twiss(ring, delta=0.0, orbit=None, eps=FD_EPS):
    orbit, R, orb_s = transfer_matrices(ring, delta, orbit, eps=eps)
    tw = twiss_from_matrices(R, np.append(ring.s_out, ring.circumference), ring.circumference)
    tw.update(orbit=orb_s, closed_orbit=orbit, delta=delta, R=R)
    return tw


def one_turn_matrix(ring, delta=0.0, orbit=None, eps=FD_EPS):
    """Матрица оборота 6×6 вокруг орбиты (4-точечная схема)."""
    if orbit is None:
        orbit = closed_orbit(ring, delta)
    cols = [0, 1, 2, 3, 4, 5]
    eps = np.broadcast_to(np.asarray(eps, float), (6,))
    fin = ring.track(_fd_columns(orbit, cols, eps))
    return orbit, _fd_collect(fin, cols, eps, (6, 6))


def damping_from_matrix(M, circumference):
    """Собственные числа матрицы оборота с излучением: |λ| = exp(−T0/τ) ⇒ τ (оборотов) и частоты arg λ/2π.
    Моды сопоставляются плоскостям по наибольшей доле собственного вектора. Возврат dict по x, y, s."""
    lam, vec = np.linalg.eig(M)
    out = {}
    used = set()
    for name, idx in (("x", (0, 1)), ("y", (2, 3)), ("s", (4, 5))):
        best = None
        for k in range(6):
            if k in used or lam[k].imag < 0:
                continue
            w = np.sum(np.abs(vec[idx, k]) ** 2) / np.sum(np.abs(vec[:, k]) ** 2)
            if best is None or w > best[0]:
                best = (w, k)
        k = best[1]; used.add(k)
        conj = int(np.argmin(np.abs(lam - np.conj(lam[k])) + (np.arange(6) == k) * 1e9)); used.add(conj)
        out[name] = dict(lam=lam[k], tau_turns=-1.0 / math.log(abs(lam[k])) if abs(lam[k]) < 1 else np.inf,
                         tune=abs(np.angle(lam[k])) / (2 * np.pi))
    return out


def tunes_from_matrix(M, strict=True):
    """Дробные частоты по матрице оборота (знак sin μ — по R12). Неустойчивая плоскость (|cos μ| ≥ 1)
    даёт NaN: обрезать косинус до ±1 нельзя, иначе неустойчивость выглядит как частота 0 или 1/2."""
    out = []
    for i in (0, 2):
        m = M[i:i + 2, i:i + 2]
        cmu = 0.5 * np.trace(m)
        if strict and abs(cmu) >= 1:
            out.append(float("nan")); continue
        fr = math.acos(max(-1.0, min(1.0, cmu))) / (2 * np.pi)
        out.append(1 - fr if m[0, 1] < 0 else fr)
    return tuple(out)


def chromaticity(ring, deltas=None, order=3):
    """Дробные ν(δ) по замкнутым орбитам вне импульса и ряд Тейлора (ν0, ξ1, ξ2, ξ3)."""
    deltas = np.array(deltas if deltas is not None else np.linspace(-0.01, 0.01, 9))
    nx = []; ny = []
    for d in deltas:
        _, M = one_turn_matrix(ring, d)
        qx, qy = tunes_from_matrix(M); nx.append(qx); ny.append(qy)
    nx = np.unwrap(np.array(nx) * 2 * np.pi) / (2 * np.pi)
    ny = np.unwrap(np.array(ny) * 2 * np.pi) / (2 * np.pi)
    fac = np.array([math.factorial(k) for k in range(order + 1)])
    return dict(delta=deltas, nux=nx, nuy=ny,
                xi_x=np.polyfit(deltas, nx, order)[::-1] * fac, xi_y=np.polyfit(deltas, ny, order)[::-1] * fac)


def element_matrices(ring, delta=0.0, orbit=None):
    """Матрица 6×6 каждого магнита вместе с дрейфом перед ним: R_j·R_{j−1}⁻¹."""
    _, R, _ = transfer_matrices(ring, delta, orbit)
    out = [R[0]] + [R[j] @ np.linalg.inv(R[j - 1]) for j in range(1, len(R))]
    return np.array(out)


# ============================================================================ равновесные моменты (метод Чао)
def _bend_slices(m):
    """Изгиб, разрезанный на свои N срезов интегратора: копии магнита с L/N, θ/N и краями только на
    крайних срезах. Последовательный проход копий повторяет тело целиком (те же операции в том же порядке)."""
    import copy
    N = m._nsl_sector
    out = []
    for i in range(N):
        c = copy.copy(m)
        c.Leff = m.Leff / N; c.L = m.L / N; c.theta0 = m.theta0 / N; c._nsl_sector = 1
        c.kEdge1 = m.kEdge1 if i == 0 else 0.0
        c.kEdge2 = m.kEdge2 if i == N - 1 else 0.0
        out.append(c)
    return out


def isr_diffusion_matrix(ring, orbit, eps=FD_EPS):
    """Матрица диффузии за оборот у плоскости наблюдения от квантового излучения, B = Σ σ²·r·rᵀ.

    σ² среза — дисперсия скачка δ, которую даёт ядро SynchrotronSim в изгибе (ISR: isr_coef²·F^1.5·ds по
    подшагам среза), r — столбец δ матрицы переноса от середины среза до плоскости после ВЧ. Матрицы до
    срезов снимаются трекингом разностями, как и прочие матрицы; излучение при этом только классическое."""
    gam = mt._design_gamma(getattr(ring, "design_mev", None))
    if hasattr(ring, "track_pieces"):
        return _isr_diffusion_cpp(ring, orbit, gam, eps)
    pieces = []                                              # (магнит-кусок, стык перед ним?, isr σ²/ds)
    for j, m in enumerate(ring.mags):
        if m.sector:
            rho0 = -m.Leff / m.theta0
            isr = mt._ISR_K * gam ** 2.5 / abs(rho0) ** 1.5
            for i, c in enumerate(_bend_slices(m)):
                pieces.append((c, i == 0, isr, j))
        else:
            pieces.append((m, True, 0.0, j))
    cols = [0, 1, 2, 3, 5]
    epsv = np.full(6, float(eps))
    pts = _fd_columns(orbit, cols, epsv)
    n = len(pieces)
    rec = np.full((len(pts), n + 1, 6), np.nan)
    for p, st in enumerate(pts):
        cur = tuple(map(float, st)); rec[p, 0] = cur
        prev_mag = ring.plane
        for q, (piece, junction, _, j) in enumerate(pieces):
            if junction:
                cur = mt.modified_drift(cur, prev_mag, ring.mags[j])
                prev_mag = ring.mags[j]
            cur = piece.body(cur); rec[p, q + 1] = cur
    R = _fd_collect(rec, cols, epsv, (n + 1, 6, 6))
    R[:, 4, 4] = 1.0
    _, M = one_turn_matrix(ring, orbit=orbit, eps=eps)
    B = np.zeros((6, 6))
    for q, (piece, _, isr, j) in enumerate(pieces):
        if not isr:
            continue
        mid = 0.5 * (np.linalg.inv(R[q]) + np.linalg.inv(R[q + 1]))          # вход и выход среза
        r = M @ mid[:, 5]
        x, y, d = 0.5 * (rec[0, q, 0] + rec[0, q + 1, 0]), 0.5 * (rec[0, q, 2] + rec[0, q + 1, 2]), rec[0, q, 5]
        rho0 = -piece.Leff / piece.theta0
        Fx, Fy = mt._csbend_fields(x, y, piece.k1 * rho0, piece.k2 * rho0, piece.k3 * rho0, 0.0, 1.0 / rho0)
        sig2 = isr ** 2 * (Fx * Fx + Fy * Fy) ** 1.5 * piece.Leff * (1 + d) ** 4
        B += sig2 * np.outer(r, r)
    return M, B


def _isr_diffusion_cpp(ring, orbit, gam, eps):
    """isr_diffusion_matrix на C++-ядре: срезы изгибов и траектории по ним считает сервер кольца."""
    meta = ring.piece_meta()                                 # (магнит, стык, изгиб, Leff, θ0, k1, k2, k3)
    cols = [0, 1, 2, 3, 5]
    epsv = np.full(6, float(eps))
    pts = _fd_columns(orbit, cols, epsv)
    rec = ring.track_pieces(pts)
    n = len(meta)
    R = _fd_collect(rec, cols, epsv, (n + 1, 6, 6))
    R[:, 4, 4] = 1.0
    _, M = one_turn_matrix(ring, orbit=orbit, eps=eps)
    B = np.zeros((6, 6))
    for q, (_, _, sector, Leff, th, k1, k2, k3) in enumerate(meta):
        if not sector:
            continue
        rho0 = -Leff / th
        isr = mt._ISR_K * gam ** 2.5 / abs(rho0) ** 1.5
        mid = 0.5 * (np.linalg.inv(R[q]) + np.linalg.inv(R[q + 1]))
        r = M @ mid[:, 5]
        x, y, d = 0.5 * (rec[0, q, 0] + rec[0, q + 1, 0]), 0.5 * (rec[0, q, 2] + rec[0, q + 1, 2]), rec[0, q, 5]
        Fx, Fy = mt._csbend_fields(x, y, k1 * rho0, k2 * rho0, k3 * rho0, 0.0, 1.0 / rho0)
        sig2 = isr ** 2 * (Fx * Fx + Fy * Fy) ** 1.5 * Leff * (1 + d) ** 4
        B += sig2 * np.outer(r, r)
    return M, B


def equilibrium_moments(M, B):
    """Равновесная матрица моментов Σ = M·Σ·Mᵀ + B и собственные эмиттансы (|Im λ(Σ·J)|)."""
    from scipy.linalg import solve_discrete_lyapunov
    S = solve_discrete_lyapunov(M, B)
    J = np.zeros((6, 6))
    for i in (0, 2, 4):
        J[i, i + 1] = 1.0; J[i + 1, i] = -1.0
    lam, vec = np.linalg.eig(S @ J)
    em = np.zeros(3)                                         # (εx, εy, εz) — мода по наибольшей доле вектора
    for k in range(6):
        if lam[k].imag <= 0:
            continue
        w = [np.sum(np.abs(vec[i:i + 2, k]) ** 2) for i in (0, 2, 4)]
        em[int(np.argmax(w))] = abs(lam[k].imag)
    return S, em


# ============================================================================ спектр: NAFF
def naff(z, n_freq=1, window=1, search=None):
    """Частоты сигнала за оборот методом NAFF (Ласкар): окно Ханнинга^window, пик БПФ, уточнение
    максимума |⟨z, e^{2πiνk}⟩| золотым сечением, вычитание найденной гармоники. Возврат (ν[], a[] комплексн.).
    Точность ~1/N⁴ для изолированной линии. search=(lo, hi) ограничивает поиск первой частоты."""
    z = np.asarray(z, complex); N = len(z)
    k = np.arange(N)
    w = (1.0 - np.cos(2 * np.pi * k / N)) ** window
    w /= w.sum()
    zz = z - np.sum(w * z)                                    # без постоянной составляющей
    freqs = []; amps = []

    def proj(nu, sig):
        return np.sum(w * sig * np.exp(-2j * np.pi * nu * k))

    for _ in range(n_freq):
        pad = 8 * N
        spec = np.abs(np.fft.fft(w * zz, pad))
        grid = np.fft.fftfreq(pad)
        if search is not None and not freqs:
            lo, hi = search
            g = np.mod(grid, 1.0); mask = (g >= lo) & (g <= hi)
            spec = np.where(mask, spec, 0.0)
        i = int(np.argmax(spec))
        a, b = grid[i] - 1.0 / pad, grid[i] + 1.0 / pad
        gr = (math.sqrt(5) - 1) / 2
        c, d = b - gr * (b - a), a + gr * (b - a)
        fc, fd = abs(proj(c, zz)), abs(proj(d, zz))
        for _ in range(60):
            if fc > fd:
                b, d, fd = d, c, fc; c = b - gr * (b - a); fc = abs(proj(c, zz))
            else:
                a, c, fc = c, d, fd; d = a + gr * (b - a); fd = abs(proj(d, zz))
        nu = 0.5 * (a + b)
        amp = proj(nu, zz)
        freqs.append(nu % 1.0); amps.append(amp)
        zz = zz - amp * np.exp(2j * np.pi * nu * k)
    return np.array(freqs), np.array(amps)


def damping_by_envelope(tbt, plane, beta=1.0, alpha=0.0, orbit=None, window=512, step=128, search=None):
    """Время затухания (оборотов) по спаду амплитуды пооборотных колебаний: амплитуда NAFF в скользящих окнах,
    наклон ln|a| по центру окна. plane: 0 — x, 2 — y, 4 — продольная (δ, t). Возврат (τ, центры, амплитуды)."""
    tbt = np.asarray(tbt, float)
    if orbit is not None:
        tbt = tbt - np.asarray(orbit)[None, :]
    if plane in (0, 2):
        z = (tbt[:, plane] - 1j * (alpha * tbt[:, plane] + beta * tbt[:, plane + 1])) / math.sqrt(beta)
    else:
        dd = tbt[:, 5]; tt = tbt[:, 4]
        z = dd + 1j * tt * (np.std(dd) / max(np.std(tt), 1e-30))
    centers = []; amps = []
    for a in range(0, len(z) - window + 1, step):
        seg = z[a:a + window]
        f, amp = naff(seg, 1, search=search)
        centers.append(a + window / 2.0); amps.append(abs(amp[0]))
    centers = np.array(centers); amps = np.array(amps)
    slope = np.polyfit(centers, np.log(amps), 1)[0]
    return -1.0 / slope, centers, amps


def tbt_tunes(tbt, beta=(1.0, 1.0), alpha=(0.0, 0.0), orbit=None, n_lines=4, tol=2e-3):
    """Дробные частоты (νx, νy) по пооборотным координатам (n_turns, 6) у плоскости наблюдения.

    Сигнал нормирован: z = (x − i(α·x + β·px))/√β ⇒ модуль ~ √(2J), частота бетатронных колебаний положительна.
    Сначала берётся плоскость с большей амплитудой (её сильнейшая линия), затем в другой плоскости —
    сильнейшая линия, НЕ совпадающая с ±ν, ±2ν первой: при малой амплитуде вынужденные линии связи
    (напр. 1−2νy в x) бывают сильнее собственной."""
    tbt = np.asarray(tbt, float)
    if orbit is not None:
        tbt = tbt - np.asarray(orbit)[None, :]
    lines = []
    for (i, b, a) in ((0, beta[0], alpha[0]), (2, beta[1], alpha[1])):
        x = tbt[:, i]; px = tbt[:, i + 1]
        if not np.all(np.isfinite(x)) or np.max(np.abs(x)) == 0:
            return (np.nan, np.nan)
        f, amp = naff((x - 1j * (a * x + b * px)) / math.sqrt(b), n_lines)
        lines.append((f, np.abs(amp)))
    first = 0 if lines[0][1][0] >= lines[1][1][0] else 1
    nu = [np.nan, np.nan]
    nu[first] = float(lines[first][0][0])
    other = 1 - first
    ref = nu[first]
    bad = [(k * ref) % 1.0 for k in (1, -1, 2, -2)]
    for f in lines[other][0]:
        if all(min(abs(f - b), 1 - abs(f - b)) > tol for b in bad):
            nu[other] = float(f); break
    else:
        nu[other] = float(lines[other][0][0])
    return tuple(nu)


def survived_turns(tbt, limit=None):
    """Число полных оборотов до потери (NaN или |x|,|y| > limit) по пооборотным данным (n, T+1, 6)."""
    ok = np.all(np.isfinite(tbt), axis=2)
    if limit is not None:
        ok &= (np.abs(np.nan_to_num(tbt[:, :, 0], nan=1e9)) <= limit) & \
              (np.abs(np.nan_to_num(tbt[:, :, 2], nan=1e9)) <= limit)
    bad = ~ok
    first = np.where(bad.any(axis=1), bad.argmax(axis=1), tbt.shape[1])
    return first - 1


def frequency_map(tbt, beta, alpha, orbit, n_half):
    """Карта частот: NAFF по первой и второй половинам (по n_half оборотов); диффузия
    D = log10 √(Δνx² + Δνy²). Возврат (νx, νy, D) для выживших частиц, NaN для прочих."""
    n = tbt.shape[0]
    nux = np.full(n, np.nan); nuy = np.full(n, np.nan); diff = np.full(n, np.nan)
    for i in range(n):
        tr = tbt[i]
        if tr.shape[0] < 2 * n_half + 1 or not np.all(np.isfinite(tr[:2 * n_half + 1])):
            continue
        q1 = tbt_tunes(tr[:n_half], beta, alpha, orbit)
        q2 = tbt_tunes(tr[n_half:2 * n_half], beta, alpha, orbit)
        nux[i], nuy[i] = q1
        dx = (q2[0] - q1[0] + 0.5) % 1 - 0.5; dy = (q2[1] - q1[1] + 0.5) % 1 - 0.5
        diff[i] = math.log10(max(math.hypot(dx, dy), 1e-15))
    return nux, nuy, diff


def resonance_lines(order_max, box):
    """Линии m·νx + n·νy = p (|m|+|n| ≤ order_max) в прямоугольнике box=(x0,x1,y0,y1) дробных частот.
    Возврат [(m, n, p, [(x,y),(x,y)])] — отрезки, обрезанные прямоугольником."""
    x0, x1, y0, y1 = box
    out = []
    seen = set()
    for m in range(-order_max, order_max + 1):
        for n in range(0, order_max + 1):
            if abs(m) + n == 0 or abs(m) + n > order_max or (n == 0 and m < 0):
                continue
            g = math.gcd(abs(m), n)
            vals = [m * x + n * y for x in (x0, x1) for y in (y0, y1)]
            for pp in range(math.floor(min(vals)), math.ceil(max(vals)) + 1):
                key = (m, n, pp)
                if key in seen or math.gcd(g, abs(pp)) != 1:       # кратные резонансы низшего порядка не рисуем
                    continue
                pts = []
                if n != 0:
                    for x in (x0, x1):
                        y = (pp - m * x) / n
                        if y0 - 1e-12 <= y <= y1 + 1e-12:
                            pts.append((x, y))
                if m != 0:
                    for y in (y0, y1):
                        x = (pp - n * y) / m
                        if x0 - 1e-12 <= x <= x1 + 1e-12:
                            pts.append((x, y))
                pts = sorted(set((round(a, 12), round(b, 12)) for a, b in pts))
                if len(pts) >= 2:
                    seen.add(key); out.append((m, n, pp, [pts[0], pts[-1]]))
    return out


# ============================================================================ параллельный трекинг SynchrotronSim
_W_RINGS = {}
_W_LATTICE = 1


def _w_init(lattice):
    global _W_LATTICE
    _W_LATTICE = lattice


def _w_ring(env, ring_kwargs=None):
    """Кольцо процесса под набор переменных ядра (SS_SECTOR_NSL, SS_SECTOR_ORDER…) и параметры Ring (rf, sr…)
    — собирается один раз."""
    key = (tuple(sorted((env or {}).items())), repr(sorted((ring_kwargs or {}).items())))
    if key not in _W_RINGS:
        base = {k: os.environ.get(k) for k in dict(key[0])}
        os.environ.update(dict(key[0]))
        try:
            _W_RINGS[key] = Ring(_W_LATTICE, **(ring_kwargs or {}))
        finally:
            for k, v in base.items():
                if v is None:
                    os.environ.pop(k, None)
                else:
                    os.environ[k] = v
    return _W_RINGS[key]


def _w_task(task):
    name, args, kwargs, env, ring_kwargs = (tuple(task) + (None, None))[:5]
    ring = _w_ring(env, ring_kwargs)
    if name.startswith("ring."):
        return getattr(ring, name[5:])(*args, **kwargs)
    return globals()[name](ring, *args, **kwargs)


def pool_map(tasks, lattice=1, processes=None):
    """Задачи [(имя, args, kwargs[, env[, ring_kwargs]])] по процессам. Имя 'ring.метод' — метод Ring, иначе функция этого
    модуля с кольцом первым аргументом. env — переменные ядра на время сборки кольца задачи (кольца
    кэшируются в процессе). Порядок результатов = порядок задач."""
    from concurrent.futures import ProcessPoolExecutor
    processes = processes or max(1, min(len(tasks), (os.cpu_count() or 2) - 2))
    with ProcessPoolExecutor(max_workers=processes, initializer=_w_init, initargs=(lattice,)) as ex:
        return list(ex.map(_w_task, tasks))


def tbt_parallel(states, n_turns, lattice=1, env=None, processes=None, ring_kwargs=None, limit=None):
    """Пооборотный трекинг частиц по процессам (частица — задача). (n, n_turns+1, 6)."""
    res = pool_map([("ring.track_tbt", ([s], n_turns), {"limit": limit}, env, ring_kwargs) for s in states],
                   lattice, processes)
    return np.concatenate(res, axis=0)


def analytic_orbit_response(ring, alpha, index=None):
    """Отклик замкнутой орбиты на поворот одного магнита: матрицами (ring_matrix) и трекингом.

    Считается в процессе пула, потому что требует своего кольца (сошедшийся интегратор) и двух поисков
    замкнутой орбиты. Возврат dict с орбитами идеальной решётки и смещениями Δx после поворота."""
    import ring_matrix as rm
    import copy
    if index is None:
        index = next(i for i, k in enumerate(ring.kinds) if k == "quad" and i > 5)
    lm0 = rm.line_matrices(ring)
    co_m = rm.closed_orbit(lm0["one_turn_C"], lm0["one_turn_R"])
    co_t = closed_orbit(ring)
    if hasattr(ring, "with_yaw"):                        # C++: поворот и в ядре, и в описании магнитов
        pert = ring.with_yaw(index, alpha)
    else:
        pert = copy.copy(ring)
        pert.mags = list(ring.mags)
        pert.mags[index] = yaw_element(copy.copy(ring.mags[index]), alpha)
    lm = rm.line_matrices(pert)
    com = rm.closed_orbit(lm["one_turn_C"], lm["one_turn_R"])
    cot = closed_orbit(pert)
    return dict(index=index, name=ring.names[index], alpha=alpha, co_matrix=co_m, co_track=co_t,
                dx_matrix=com[0] - co_m[0], dx_track=cot[0] - co_t[0],
                C_transverse=float(np.max(np.abs(lm["one_turn_C"][:4]))))


def kicked(ring, j, theta):
    """Кольцо с однородным горизонтальным киком θ в прямом неизгибающем магните j (пробный корректор)."""
    import copy
    if hasattr(ring, "with_kick"):
        return ring.with_kick(j, theta)
    pert = copy.copy(ring); pert.mags = list(ring.mags)
    m = copy.copy(ring.mags[j]); m.theta0 = theta; m.nonzero = True; m.xref = None
    pert.mags[j] = m
    return pert


def orm_tracking_check(ring, picks=(0, 41), theta=1e-6, kinds=("kicker",)):
    """Сверка матрицы отклика с трекингом: кик θ в корректоре — замкнутая орбита матрицами и трекингом.
    Считается в процессе пула с сошедшимся интегратором (иначе сравнение упирается в недосходимость)."""
    import copy
    import ring_matrix as rm
    import ring_response as rr
    lm = rm.line_matrices(ring)
    corr = rr.correctors(ring, lm, kinds)
    M = lm["one_turn_R"]; inv = np.linalg.inv(np.eye(4) - M[:4, :4])
    co0 = closed_orbit(ring)
    out = []
    for pick in picks:
        j, s_k, Rc = corr[pick]
        e = np.zeros(6); e[1] = theta
        z0 = inv @ (M @ np.linalg.inv(Rc) @ e)[:4]
        pert = kicked(ring, j, theta)
        d = (closed_orbit(pert) - co0)[:4]
        out.append(dict(pick=pick, s=s_k, x_matrix=z0[0], x_track=d[0], px_matrix=z0[1], px_track=d[1]))
    return out


def tunes_point(ring, delta=0.0):
    """(νx, νy дробные, замкнутая орбита) при импульсе δ — по матрице оборота; NaN, если орбиты нет."""
    try:
        orbit, M = one_turn_matrix(ring, delta)
        if not np.all(np.isfinite(M)):
            raise RuntimeError
        return (*tunes_from_matrix(M), orbit)
    except (RuntimeError, np.linalg.LinAlgError, ValueError):
        return (np.nan, np.nan, np.full(6, np.nan))


def main():
    ap = argparse.ArgumentParser(description="линейная оптика кольца трекингом SynchrotronSim")
    ap.add_argument("what", choices=["twiss", "chrom", "orbit"])
    ap.add_argument("--lattice", type=int, default=1)
    ap.add_argument("--delta", type=float, default=0.0)
    a = ap.parse_args()
    ring = Ring(a.lattice)
    if a.what == "orbit":
        print("замкнутая орбита:", closed_orbit(ring, a.delta))
    elif a.what == "twiss":
        tw = twiss(ring, a.delta)
        print(f"C = {ring.circumference:.6f} м, νx = {tw['nux']:.6f}, νy = {tw['nuy']:.6f}, αc = {tw['alphac']:.6e}")
        print(f"βx max {tw['betax'].max():.4f} м, βy max {tw['betay'].max():.4f} м, ηx max {tw['etax'].max():.5f} м")
    else:
        c = chromaticity(ring)
        print("ξx (ξ1, ξ2, ξ3):", c["xi_x"][1:], "ξy:", c["xi_y"][1:])


if __name__ == "__main__":
    main()
