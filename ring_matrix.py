# -*- coding: utf-8 -*-
"""Аналитические матрицы элементов и их перемножение — механизм анализа установки, как в elegant.

Отличие от `ring_analysis`: там матрицы снимаются разностями по ТРЕКИНГУ, здесь — считаются ФОРМУЛАМИ
для каждого элемента и перемножаются вдоль линии. Это даёт мгновенную линейную оптику, производные по
параметрам (матрицы отклика) и разложение «кто во что вносит», а заодно независимую проверку интегратора:
матрица — точное решение, трекинг — его приближение с конечным числом срезов.

Соглашения — elegant (проверено по matrix_output): z = (x, x′, y, y′, s, δ), пятая координата — ПУТЬ
(у дрейфа R56 = 0), K1 > 0 фокусирует по x, края изгиба — тонкие линзы tan(E)/ρ (edge_order = 1),
секступоли и октуполи в первом порядке — дрейфы. Матрицы симплектичны с формой, у которой продольная
сопряжённая пара — (−s, δ) (следствие того, что пятая координата — путь, а не время). На опорной орбите канонический импульс ядра совпадает
с наклоном, поэтому матрицы сравнимы с `ring_analysis` напрямую.

НОВОЕ ПРОТИВ elegant — СТЫКИ. В SynchrotronSim элементы соединены модифицированными дрейфами: частица
летит по прямой от выходной грани одного элемента до входной грани другого, а грани могут быть НЕ
параллельны (боксовая модель, повороты установки, неидеальная геометрия). Такое отображение не линейное,
а АФФИННОЕ: z → R·z + C. Поэтому линия собирается из пар (C, R), а замкнутая орбита получается из C:
x_замк = (I − R)⁻¹·C — то есть тот же механизм сразу даёт искажение орбиты от рассогласований.

Пример:
  py -3.14 ring_matrix.py twiss           # оптика кольца по матрицам
  py -3.14 ring_matrix.py check           # сверка с матрицами из трекинга
"""
import math, argparse
import numpy as np

import ring_analysis as ra
import matrix_tracker as mt


# ============================================================================ матрицы элементов
def _foc_block(k2, L):
    """Блок 2×2 фокусирующего участка: k2 = k² (>0 фокус, <0 дефокус, 0 дрейф). Возврат (c, s, c′, s′),
    где c = cos(kL) либо ch(kL), s = sin(kL)/k либо sh(kL)/k, и производные по длине."""
    if k2 > 1e-14:
        k = math.sqrt(k2)
        return math.cos(k * L), math.sin(k * L) / k, -k2 * math.sin(k * L) / k, math.cos(k * L)
    if k2 < -1e-14:
        k = math.sqrt(-k2)                              # дефокус: R21 = +k·sh(kL) = −k²·sh/k
        return math.cosh(k * L), math.sinh(k * L) / k, -k2 * math.sinh(k * L) / k, math.cosh(k * L)
    return 1.0, L, 0.0, 1.0


def drift(L):
    R = np.eye(6); R[0, 1] = L; R[2, 3] = L
    return R


def quad(L, K1, K1v=None):
    """Толстый квадруполь: K1 > 0 — фокус по x. K1v (вертикальный градиент) по умолчанию = K1."""
    K1v = K1 if K1v is None else K1v
    R = np.eye(6)
    R[0, 0], R[0, 1], R[1, 0], R[1, 1] = _foc_block(K1, L)
    R[2, 2], R[2, 3], R[3, 2], R[3, 3] = _foc_block(-K1v, L)
    return R


def sector_body(L, h, K1):
    """Тело секторного изгиба с градиентом: кривизна h = 1/ρ, K1 — нормальный градиент.
    Горизонтально работает kx² = h² + K1 (слабая фокусировка изгиба), вертикально ky² = −K1."""
    R = np.eye(6)
    kx2 = h * h + K1
    c, s, cp, sp = _foc_block(kx2, L)
    R[0, 0], R[0, 1], R[1, 0], R[1, 1] = c, s, cp, sp
    R[2, 2], R[2, 3], R[3, 2], R[3, 3] = _foc_block(-K1, L)
    if h:
        # дисперсия: частный интеграл по изогнутой орбите; путь — симплектически сопряжённые члены
        if abs(kx2) > 1e-14:
            d16 = h * (1 - c) / kx2
            d56 = h * h * (L - s) / kx2
        else:
            d16 = h * L * L / 2
            d56 = h * h * L ** 3 / 6
        R[0, 5] = d16; R[1, 5] = h * s
        R[4, 0] = h * s; R[4, 1] = d16; R[4, 5] = d56
    return R


def edge(kedge):
    """Тонкая краевая линза грани изгиба: x′ += k·x, y′ −= k·y (k = tan(E)/ρ, edge_order = 1)."""
    R = np.eye(6); R[1, 0] = kedge; R[3, 2] = -kedge
    return R


def element_matrix(m):
    """(C, R) элемента решётки SynchrotronSim по его физкоэффициентам.

    ТОЛЬКО секторный режим ядра (SS_SECTOR=1). В боксовой модели грани изгиба ⊥ оси ящика, опорная орбита
    входит в них под углом θ/2, и телу соответствует другая матрица — здесь она не реализована. Ядро в секторном режиме:
    изгиб — край·тело·край, прямой мультиполь — толстая линза (секступоль/октуполь линейно = дрейф)."""
    if not mt.SECTOR:
        raise RuntimeError("аналитические матрицы описаны для секторного режима ядра: нужен SS_SECTOR=1 "
                           "ДО импорта matrix_tracker")
    L = float(getattr(m, "Leff", 0.0) or 0.0)
    C = np.array([0.0, 0.0, 0.0, 0.0, L, 0.0])          # путь опорной частицы (пятая координата — абсолютная)
    if L <= 0 or not getattr(m, "nonzero", False):
        return C, drift(L)
    if getattr(m, "theta0", 0.0):                      # изгиб: h > 0 в рамке ядра (rho0 = −Leff/θ0)
        h = -m.theta0 / L
        R = edge(m.kEdge2) @ sector_body(L, h, m.k1) @ edge(m.kEdge1)
        return C, R
    if m.kind == "quad":
        return C, quad(L, m.k1, m.k1v)
    return C, drift(L)                                  # секступоль, октуполь, корректор


# ============================================================================ стыки (модифицированные дрейфы)
def junction_matrix(A, B):
    """(C, R) модифицированного дрейфа: прямой пролёт от выходной грани A до входной грани B.

    Грани заданы рамками (точка, ось, орт u, орт v). Частица летит из P = a + x·u_a + y·v_a по направлению
    n = (px·u_a + py·v_a + p_s·a_ax)/(1+δ) до плоскости ⊥ b_ax через b, затем координаты берутся в базисе B.
    На опорной орбите (x = px = 0, δ = 0) направление n не зависит от δ, поэтому столбец δ у стыка нулевой;
    зато при непараллельных гранях появляется НЕНУЛЕВОЙ свободный член C — наклон грани к опорной орбите.
    """
    aP, aAx, aU, aV = A.fr_out
    bP, bAx, bU, bV = B.fr_in
    if getattr(B, "same_plane_as_prev", False):
        t0 = 0.0
    else:
        den = float(bAx @ aAx)
        if abs(den) < 1e-12:
            raise ValueError("грани стыка почти параллельны оси пролёта")
        t0 = float(bAx @ (bP - aP)) / den
    Q0 = aP + t0 * aAx
    rel = Q0 - bP
    C = np.array([float(rel @ bU), float(aAx @ bU), float(rel @ bV), float(aAx @ bV), t0, 0.0])
    # производные пролёта: ∂t/∂x = −(b_ax·u_a)/(b_ax·a_ax), ∂t/∂px = −t0·(b_ax·u_a)/(b_ax·a_ax)
    den = float(bAx @ aAx) if t0 or not getattr(B, "same_plane_as_prev", False) else 1.0
    tx = -float(bAx @ aU) / den; ty = -float(bAx @ aV) / den
    R = np.zeros((6, 6))
    dP = {0: aU, 1: aU * t0, 2: aV, 3: aV * t0}          # ∂Q/∂z без члена вдоль оси
    dt = {0: tx, 1: t0 * tx, 2: ty, 3: t0 * ty}
    for j in range(4):
        dQ = dP[j] + dt[j] * aAx
        R[0, j] = float(dQ @ bU); R[2, j] = float(dQ @ bV)
        R[4, j] = dt[j]
    for j, e in ((1, aU), (3, aV)):                      # направление меняют только px, py
        R[1, j] = float(e @ bU); R[3, j] = float(e @ bV)
    # px₂ = (1+δ)·(n·u_b), а n от δ не зависит ⇒ у наклонной грани появляется столбец δ (он же C)
    R[1, 5] = float(aAx @ bU); R[3, 5] = float(aAx @ bV)
    R[4, 4] = 1.0; R[5, 5] = 1.0
    return C, R


# ============================================================================ сборка линии
def compose(first, second):
    """Композиция аффинных отображений: сначала first, затем second."""
    C1, R1 = first; C2, R2 = second
    return C2 + R2 @ C1, R2 @ R1


def line_matrices(ring):
    """Матрицы от точки наблюдения до выхода каждого магнита и за оборот.
    Возврат dict: C, R (n_mag+1, …), one_turn_C, one_turn_R, s."""
    mags = ring.mags
    accC = np.zeros(6); accR = np.eye(6)
    Cs = np.empty((len(mags) + 1, 6)); Rs = np.empty((len(mags) + 1, 6, 6))
    prev = ring.plane
    for j, m in enumerate(mags):
        accC, accR = compose((accC, accR), junction_matrix(prev, m))
        accC, accR = compose((accC, accR), element_matrix(m))
        Cs[j] = accC; Rs[j] = accR
        prev = m
    accC, accR = compose((accC, accR), junction_matrix(prev, ring.plane))
    Cs[-1] = accC; Rs[-1] = accR
    return dict(C=Cs, R=Rs, one_turn_C=accC, one_turn_R=accR,
                s=np.append(ring.s_out, ring.circumference))


# ============================================================================ интегралы излучения
C_LIGHT = 299792458.0
MC2_MEV = 0.51099895069
C_GAMMA = 8.846273435551005e-5      # 4π·r_e/3/(m c²)³, м/ГэВ³ — потеря энергии
C_Q = 3.8319402812e-13              # 55/(32√3)·ħ/(m c), м — квантовое возбуждение


def radiation_integrals(ring, nodes=8):
    """Интегралы излучения I1…I5 по АНАЛИТИЧЕСКОЙ оптике (без трекинга) и всё, что из них следует.

    Подынтегральные функции нужны ВНУТРИ изгибов, поэтому каждый изгиб проходится узлами Гаусса–Лежандра:
    для узла строится матрица от точки наблюдения (край + часть тела), по набору матриц обычным разбором
    получаются β, α, η, η′, и интегралы берутся квадратурой. У прямых элементов h = 0, и они не дают вклада.

        I1 = ∮ η·h ds            I2 = ∮ h² ds            I3 = ∮ |h|³ ds
        I4 = ∮ η·h·(h² + 2K₁) ds − Σ η·h²·tg(E) по граням       (краевой член)
        I5 = ∮ |h|³·H ds,  H = γη² + 2αηη′ + βη′²

    Возврат dict: I1…I5, U0 (МэВ), Jx/Jy/Js, taux/tauy/taudelta (с), ex0 (м·рад), Sdelta0, alphac;
    при rf=dict(volt_mv, harmonic) — ещё Qs, phi_s (°) и Sz0 (м)."""
    mags = ring.mags
    u_n, w_n = np.polynomial.legendre.leggauss(nodes)
    accC = np.zeros(6); accR = np.eye(6)
    prev = ring.plane
    Rn, wn, hn, kn = [], [], [], []
    Re, fe = [], []                                     # грани изгибов: краевой член I4
    for m in mags:
        accC, accR = compose((accC, accR), junction_matrix(prev, m))
        L = float(getattr(m, "Leff", 0.0) or 0.0)
        th = float(getattr(m, "theta0", 0.0) or 0.0)
        if L > 0 and getattr(m, "nonzero", False) and th:
            h = -th / L
            e1 = edge(m.kEdge1)
            for uu, ww in zip(u_n, w_n):
                t = 0.5 * (uu + 1.0) * L
                Rn.append(sector_body(t, h, m.k1) @ e1 @ accR)
                wn.append(0.5 * ww * L); hn.append(h); kn.append(m.k1)
            if m.kEdge1:                                # −η·h·k_края на входной грани (k = tg(E)·h)
                Re.append(accR); fe.append(h * m.kEdge1)
        accC, accR = compose((accC, accR), element_matrix(m))
        if L > 0 and getattr(m, "nonzero", False) and th and m.kEdge2:
            Re.append(accR); fe.append((-th / L) * m.kEdge2)
        prev = m
    accC, accR = compose((accC, accR), junction_matrix(prev, ring.plane))
    if not Rn:
        raise RuntimeError("в кольце нет изгибающих элементов")
    n_n, n_e = len(Rn), len(Re)
    tw = ra.twiss_from_matrices(np.array(Rn + Re + [accR]), np.arange(n_n + n_e + 1), ring.circumference)
    w = np.array(wn); h = np.array(hn); k1 = np.array(kn)
    b, a = tw["betax"][:n_n], tw["alphax"][:n_n]
    eta, etap = tw["etax"][:n_n], tw["etaxp"][:n_n]
    i4_edge = float(np.sum(tw["etax"][n_n:n_n + n_e] * np.array(fe))) if n_e else 0.0
    H = (1.0 + a * a) / b * eta ** 2 + 2.0 * a * eta * etap + b * etap ** 2
    I = [float(np.sum(w * eta * h)), float(np.sum(w * h * h)), float(np.sum(w * np.abs(h) ** 3)),
         float(np.sum(w * eta * h * (h * h + 2.0 * k1))) - i4_edge,
         float(np.sum(w * np.abs(h) ** 3 * H))]
    E_gev = ring.design_mev / 1e3
    gamma = ring.design_mev / MC2_MEV
    u0 = C_GAMMA * E_gev ** 4 * I[1] / (2 * math.pi) * 1e3          # МэВ
    jx = 1.0 - I[3] / I[1]; jy = 1.0; js = 2.0 + I[3] / I[1]
    t0 = ring.circumference / C_LIGHT
    out = dict(zip(("I1", "I2", "I3", "I4", "I5"), I))
    out.update(U0=u0, Jx=jx, Jy=jy, Jdelta=js,
               taux=2 * ring.design_mev * t0 / (jx * u0), tauy=2 * ring.design_mev * t0 / (jy * u0),
               taudelta=2 * ring.design_mev * t0 / (js * u0),
               ex0=C_Q * gamma ** 2 * I[4] / (jx * I[1]),
               Sdelta0=math.sqrt(C_Q * gamma ** 2 * I[2] / (js * I[1])),
               alphac=I[0] / ring.circumference)
    rf = getattr(ring, "rf", None)
    if rf:
        volt = float(rf["volt_mv"]); harm = float(rf["harmonic"])
        phi_s = math.pi - math.asin(u0 / volt)                      # выше перехода: синхронная фаза тупая
        qs = math.sqrt(harm * out["alphac"] * volt * abs(math.cos(phi_s)) / (2 * math.pi * ring.design_mev))
        out.update(Qs=qs, phi_s=math.degrees(phi_s),
                   Sz0=out["alphac"] * ring.circumference * out["Sdelta0"] / (2 * math.pi * qs))
    return out


def closed_orbit(one_turn_C, one_turn_R):
    """Замкнутая орбита из аффинной матрицы оборота: z = (I − R)⁻¹·C (поперечные 4 координаты)."""
    M = one_turn_R[:4, :4]; c = one_turn_C[:4]
    z = np.linalg.solve(np.eye(4) - M, c)
    return np.array([z[0], z[1], z[2], z[3], 0.0, 0.0])


def twiss(ring, lm=None):
    """Функции Твисса и частоты по перемноженным матрицам (тот же разбор, что у матриц из трекинга)."""
    lm = lm or line_matrices(ring)
    tw = ra.twiss_from_matrices(lm["R"], lm["s"], ring.circumference)
    tw["closed_orbit"] = closed_orbit(lm["one_turn_C"], lm["one_turn_R"])
    tw["C"] = lm["C"]
    return tw


# ============================================================================ CLI
def main():
    ap = argparse.ArgumentParser(description="оптика кольца по аналитическим матрицам элементов")
    ap.add_argument("what", choices=["twiss", "check", "radiation"])
    ap.add_argument("--lattice", type=int, default=1)
    a = ap.parse_args()
    ring = ra.Ring(a.lattice)
    if a.what == "twiss":
        tw = twiss(ring)
        print(f"C = {ring.circumference:.6f} м, νx = {tw['nux']:.6f}, νy = {tw['nuy']:.6f}, "
              f"αc = {tw['alphac']:.6e}")
        print(f"βx max {tw['betax'].max():.4f} м, βy max {tw['betay'].max():.4f} м, "
              f"ηx max {tw['etax'].max():.5f} м, замкнутая орбита x = {tw['closed_orbit'][0] * 1e6:.3f} мкм")
    elif a.what == "radiation":
        r = radiation_integrals(ring)
        print("  ".join(f"{k} = {r[k]:.6g}" for k in ("I1", "I2", "I3", "I4", "I5")))
        print(f"U0 = {r['U0']:.6g} МэВ, Jx = {r['Jx']:.6g}, τx = {r['taux'] * 1e3:.5g} мс, "
              f"εx = {r['ex0'] * 1e9:.5g} нм, σδ = {r['Sdelta0']:.6g}, αc = {r['alphac']:.6g}")
    else:
        lm = line_matrices(ring)
        el = ra.element_matrices(ring)                  # из трекинга (конечное число срезов)
        an = [compose(junction_matrix(ring.plane if j == 0 else ring.mags[j - 1], m), element_matrix(m))[1]
              for j, m in enumerate(ring.mags)]
        d = np.array([np.max(np.abs(a_ - e[:6, :6])) for a_, e in zip(an, el)])
        print(f"матрицы элементов: медиана |Δ| {np.median(d):.2e}, максимум {d.max():.2e} "
              f"(элемент {ring.names[int(np.argmax(d))]})")


if __name__ == "__main__":
    main()
