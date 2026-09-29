# -*- coding: utf-8 -*-
"""Матрицы отклика установки: орбита на корректоры, частоты и бета-биения на ошибки градиента.

Строится поверх аналитического механизма `ring_matrix` (матрицы элементов и аффинные стыки). Отклик
считается ОДНИМ перемножением матриц, поэтому вся матрица 400×120 получается за секунды, а не трекингом.

  orbit_response(ring)     отклик орбиты в точках наблюдения на единичный кик каждого корректора, м/рад
  svd_solve(M, orbit)      токи/кики корректоров по орбите, обращение SVD с отсечкой малых чисел
  tune_response(ring)      dν/d(K1·L) на каждом элементе (из β)
  beta_beat(ring, j, dKL)  бета-биение вдоль кольца от ошибки градиента одного элемента
  tolerances(ring)         сводка по типам устройств: сдвиг частот и бета-биение на относительную ошибку

Кик корректора связан с физкоэффициентом и током библиотеки: θ = (q/p)·B·L_eff, а токи обмоток даёт
линейная матрица корректора (`saturation_solve.inverse`), так что отклик замыкается на реальные токи.

Пример:
  py -3.14 ring_response.py orbit        # матрица отклика орбиты и её спектр SVD
  py -3.14 ring_response.py tolerances   # допуска на градиенты по типам устройств
"""
import math, argparse
import numpy as np

import ring_analysis as ra
import ring_matrix as rm
import matrix_tracker as mt


# ============================================================================ точки наблюдения и корректоры
def bpm_stations(config=None, source=None, beamline=None, kind="MONI"):
    """Продольные координаты датчиков положения (MONI) по описанию линии .lte.
    config/source/beamline — описание установки (elegant_diag.machine_lattice)."""
    import elegant_diag as ed
    par = ed.element_params(config, source)
    names = ed.line_names(config, source, beamline)
    s = 0.0; out = []
    for n in names:
        q = par[n.upper()]
        if q["type"].upper() == kind:
            out.append((n, s))
        s += q.get("L", 0.0)
    return out


def station_affine(ring, lm, stations):
    """Аффинные отображения от точки наблюдения до каждой станции: (C, R). Станция лежит в дрейфе после
    последнего магнита перед ней, поэтому к накопленной паре дописывается дрейф на остаток длины."""
    s_out = ring.s_out
    Cs = []; Rs = []
    for name, s in stations:
        j = int(np.searchsorted(s_out, s, side="right")) - 1        # последний магнит перед станцией
        C0 = lm["C"][j] if j >= 0 else np.zeros(6)
        R0 = lm["R"][j] if j >= 0 else np.eye(6)
        d = s - (s_out[j] if j >= 0 else 0.0)
        Cd = np.array([0.0, 0.0, 0.0, 0.0, d, 0.0])
        C, R = rm.compose((C0, R0), (Cd, rm.drift(d)))
        Cs.append(C); Rs.append(R)
    return np.array(Cs), np.array(Rs)


def station_matrices(ring, lm, stations):
    """Накопленные матрицы от точки наблюдения до каждой станции (магнит + остаток дрейфа)."""
    return station_affine(ring, lm, stations)[1]


def station_orbit(ring, lm, stations, plane="x"):
    """Замкнутая орбита в точках наблюдения по аффинной матрице оборота (рассогласования сидят в C)."""
    C, R = station_affine(ring, lm, stations)
    z0 = rm.closed_orbit(lm["one_turn_C"], lm["one_turn_R"])
    row = 0 if plane == "x" else 2
    return (R[:, row, :4] @ z0[:4]) + C[:, row]


def correctors(ring, lm, kinds=("kicker",)):
    """Корректоры кольца: (индекс магнита, s центра, матрица от точки наблюдения до центра).
    kinds — типы устройств, работающих корректорами (у решёток без штатных корректоров это, например,
    секступоли нулевой силы)."""
    idx = [j for j, k in enumerate(ring.kinds) if k in kinds]
    out = []
    for j in idx:
        m = ring.mags[j]
        Rprev = lm["R"][j - 1] if j > 0 else np.eye(6)
        _, Rj = rm.junction_matrix(ring.mags[j - 1] if j > 0 else ring.plane, m)
        Rc = rm.drift(m.Leff / 2.0) @ Rj @ Rprev                     # кик — в середине элемента
        out.append((j, ring.s_out[j] - m.Leff / 2.0, Rc))
    return out


# ============================================================================ отклик орбиты
def orbit_response(ring, lm=None, stations=None, plane="x", kinds=("kicker",)):
    """Матрица отклика орбиты (станции × корректоры), м на радиан кика.

    Замкнутая орбита от кика θ в точке k: z₀ = (I − M)⁻¹·M·R_k⁻¹·e·θ у точки наблюдения, далее переносится
    до станции; станциям ПОСЛЕ корректора добавляется прямой вклад R_b·R_k⁻¹·e·θ."""
    lm = lm or rm.line_matrices(ring)
    stations = stations or bpm_stations()
    Rb = station_matrices(ring, lm, stations)
    M = lm["one_turn_R"]
    inv = np.linalg.inv(np.eye(4) - M[:4, :4])
    e = np.zeros(6); e[1 if plane == "x" else 3] = 1.0
    row = 0 if plane == "x" else 2
    s_b = np.array([s for _, s in stations])
    cs = correctors(ring, lm, kinds)
    out = np.zeros((len(stations), len(cs)))
    for c, (j, s_k, Rc) in enumerate(cs):
        Rk_inv = np.linalg.inv(Rc)
        v = (M @ Rk_inv @ e)[:4]                                     # кик, донесённый до конца оборота
        z0 = inv @ v                                                 # замкнутая орбита в точке наблюдения
        direct = (Rb @ (Rk_inv @ e))[:, row]                         # прямой вклад для станций после кика
        out[:, c] = (Rb[:, row, :4] @ z0) + np.where(s_b > s_k, direct, 0.0)
    return out


def svd_solve(R, orbit, n_sv=None, rcond=1e-6):
    """Кики корректоров, компенсирующие орбиту: минимальная норма по SVD с отсечкой малых чисел.
    n_sv — сколько сингулярных чисел оставить (None — все выше rcond). Возврат (кики, сингулярные числа)."""
    U, S, Vt = np.linalg.svd(R, full_matrices=False)
    keep = len(S) if n_sv is None else int(n_sv)
    inv_s = np.where(np.arange(len(S)) < keep, np.where(S > rcond * S[0], 1.0 / np.maximum(S, 1e-300), 0.0), 0.0)
    return -(Vt.T * inv_s) @ (U.T @ np.asarray(orbit, float)), S


def rigidity(design_mev=6000.0):
    """Магнитная жёсткость Bρ = p/q проектного пучка, Тл·м (p в МэВ/c)."""
    return design_mev * 1e6 / 299792458.0


def kick_to_currents(ring, j, theta, design_mev=6000.0):
    """Кик корректора (рад) → поле B (Тл) и токи обмоток через линейную матрицу корректора.

    ВАЖНО: пересчёт идёт по ПРОЕКТНОЙ жёсткости (6 ГэВ), а не по суррогатной энергии трекинга
    ([[energy-convention-surrogate-gamma]]): B = θ·Bρ/L_eff. Иначе поле выйдет в ~1700 раз меньше."""
    m = ring.mags[j]
    B = theta * rigidity(design_mev) / m.Leff
    try:
        import saturation_solve as ss
        cur = ss.inverse("corrector", {"B": B})
    except Exception as exc:                                        # таблица матрицы может быть не собрана
        cur = {"ошибка": str(exc)}
    return B, cur


# ============================================================================ отклик частот и бета-биения
def tune_response(ring, lm=None):
    """dνx/d(K1·L) и dνy/d(K1·L) на каждом элементе: ±β/(4π) по функциям Твисса из матриц."""
    lm = lm or rm.line_matrices(ring)
    tw = rm.twiss(ring, lm)
    bx = np.array([0.5 * (tw["betax"][j - 1] + tw["betax"][j]) if j else tw["betax"][0]
                   for j in range(len(ring.mags))])
    by = np.array([0.5 * (tw["betay"][j - 1] + tw["betay"][j]) if j else tw["betay"][0]
                   for j in range(len(ring.mags))])
    return bx / (4 * np.pi), -by / (4 * np.pi), tw


def beta_beat(ring, j, dKL, lm=None, plane="x", tw=None):
    """Бета-биение вдоль кольца от ошибки интеграла градиента ΔK·L на элементе j.

    Δβ/β(s) = −(ΔK·L)·β_k/(2·sin 2πν)·cos(2|ψ(s) − ψ_k| − 2πν) — стандартный отклик на тонкую линзу.
    tw — готовые функции Твисса (чтобы не пересчитывать их на каждый элемент)."""
    if tw is None:
        tw = rm.twiss(ring, lm or rm.line_matrices(ring))
    b = tw["beta" + plane]; ps = tw["psi" + plane]; nu = tw["nu" + plane]
    bk = 0.5 * (b[j - 1] + b[j]) if j else b[0]
    pk = 0.5 * (ps[j - 1] + ps[j]) if j else ps[0]
    sign = 1.0 if plane == "x" else -1.0
    return -sign * dKL * bk / (2 * math.sin(2 * math.pi * nu)) * np.cos(2 * np.abs(ps - pk) - 2 * math.pi * nu), tw


def tolerances(ring, rel_error=1e-3, lm=None):
    """Чувствительность к относительной ошибке градиента по типам устройств.

    Для каждого типа: суммарный интеграл |K1·L|, сдвиг частот при СОГЛАСОВАННОЙ ошибке (все элементы типа
    ошиблись одинаково) и СКО бета-биения при случайной независимой ошибке того же размера."""
    lm = lm or rm.line_matrices(ring)
    dnx, dny, tw = tune_response(ring, lm)
    out = {}
    for j, m in enumerate(ring.mags):
        KL = m.k1 * m.Leff
        if not KL:
            continue
        d = out.setdefault(m.dev, dict(n=0, KL=0.0, dnux=0.0, dnuy=0.0, beat2x=0.0, beat2y=0.0))
        d["n"] += 1; d["KL"] += abs(KL)
        d["dnux"] += dnx[j] * KL * rel_error
        d["dnuy"] += dny[j] * KL * rel_error
        for pl, key in (("x", "beat2x"), ("y", "beat2y")):
            bb, _ = beta_beat(ring, j, KL * rel_error, lm, pl, tw)
            d[key] += float(np.mean(bb ** 2))                        # независимые ошибки складываются по мощности
    for d in out.values():
        d["beat_rms_x"] = math.sqrt(d.pop("beat2x")); d["beat_rms_y"] = math.sqrt(d.pop("beat2y"))
    return out, tw


# ============================================================================ CLI
def main():
    ap = argparse.ArgumentParser(description="матрицы отклика установки (аналитические)")
    ap.add_argument("what", choices=["orbit", "tolerances"])
    ap.add_argument("--lattice", type=int, default=1)
    ap.add_argument("--rel-error", type=float, default=1e-3)
    a = ap.parse_args()
    ring = ra.Ring(a.lattice); lm = rm.line_matrices(ring)
    if a.what == "orbit":
        R = orbit_response(ring, lm)
        S = np.linalg.svd(R, compute_uv=False)
        print(f"матрица отклика: {R.shape[0]} станций × {R.shape[1]} корректоров, "
              f"макс {np.max(np.abs(R)):.4f} м/рад")
        print(f"сингулярные числа: max {S[0]:.4f}, min {S[-1]:.4g}, число обусловленности {S[0] / S[-1]:.3g}")
    else:
        tol, _ = tolerances(ring, a.rel_error)
        print(f"относительная ошибка градиента {a.rel_error:g}:")
        for dev, d in sorted(tol.items(), key=lambda kv: -abs(kv[1]["dnux"])):
            print(f"  {dev:22s} n={d['n']:4d} Σ|K1L|={d['KL']:8.3f}  Δνx={d['dnux']:+.5f} Δνy={d['dnuy']:+.5f}  "
                  f"СКО Δβ/β: x {d['beat_rms_x'] * 100:.3f}% y {d['beat_rms_y'] * 100:.3f}%")


if __name__ == "__main__":
    main()
