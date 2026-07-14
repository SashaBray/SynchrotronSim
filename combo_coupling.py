# -*- coding: utf-8 -*-
"""Калибровка матрицы связи МНОГООБМОТОЧНЫХ элементов по КОМБИНАЦИЯМ токов (через диффур-ядро).

Проблема одиночной калибровки: match_coefficients свипирует каждую обмотку ПООТДЕЛЬНОСТИ при малых
токах (пучок у оси) и собирает матрицу в предположении суперпозиции. Для combined-function DQ этого
мало: при рабочих токах пучок идёт по ИЗОГНУТОЙ дуге (смещён от оси карты), сильный градиент даёт
feed-down — эффект, зависящий от КОМБИНАЦИИ и положения, которого одиночные приосевые свипы не видят.

Здесь (по ТЗ пользователя, обязательно для всех элементов с несколькими входными токами):
  1. РАБОЧИЙ ДИАПАЗОН каждого тока — по выносу пучка за апертуру (probe_excursion).
  2. МНОЖЕСТВО КОМБИНАЦИЙ токов (латинский гиперкуб в этих диапазонах, вне-апертурные отбрасываются).
  3. Для каждой комбинации — диффур-прогон 37 частиц (seb.run_multi, суперпозиция полей обмоток) и
     подгонка 8 сил (θy,θz,k1,k2,k3 + скью) universal_element.track под финальные углы (Gauss-Newton).
  4. По ВСЕМ комбинациям — линейная подгонка матрицы связи «ток обмотки → сила» (lstsq), т.е. силы
     калибруются в РАБОЧЕМ режиме, а не экстраполяцией со слабых одиночных.
Выход — матрица связи (обмотка → 8 физкоэф.) в формате coupling_matrix / element_coupling.

Запуск:  python combo_coupling.py DQ         # калибровать один элемент, печать сравнения
         python combo_coupling.py --all --write   # все многообмоточные, записать в БД
"""
import os
import sys
import csv
import numpy as np
import single_element_boris_fld as seb
import universal_element as ue
import bip_matrix as bm

QP = bm.QP
TB = os.path.join("workspace", "tables")

# многообмоточные элементы (несколько входных токов) — обязательная комбо-калибровка
GROUPS = {
    "SP": ["SP_I1_1A", "SP_I2_1A"],
    "DQ": ["DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A"],
}
STR_NAMES = ("θy", "θz", "k1", "k2", "k3", "ks1", "ks2", "ks3")
# сила → (физ.столбец coupling, множитель из phys_from_strength). θy→per_B, k1→per_G и т.д.
PHYS = {"θy": ("per_B_Tl", None), "θz": ("per_Bs_Tl", None), "k1": ("per_G_Tlm", None),
        "k2": ("per_S_Tlm2", None), "k3": ("per_O_Tlm3", None),
        "ks1": ("per_Gs_Tlm", None), "ks2": ("per_Ss_Tlm2", None), "ks3": ("per_Os_Tlm3", None)}


_EFF = None


def _eff_len(dev):
    """Эфф. длина устройства из hardware_library (та же, что кольцо использует для k·L_eff)."""
    global _EFF
    if _EFF is None:
        _EFF = {r["device_type_id"]: float(r["eff_length"]) for r in csv.DictReader(
            open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";")}
    return _EFF[dev]


def elem_geom(dev):
    """L = L_eff из библиотеки (согласованность per_* с кольцом: θ=QP·per_B·L_eff, k1·L_eff);
    aper = ПОПЕРЕЧНАЯ полу-апертура (45мм, как APER в match_coefficients) — для капа токов и
    отключения поля за краем в форвард-модели."""
    return _eff_len(dev), 0.045


def excursion_multi(devs, currents, aper):
    """Макс. поперечный вынос комбо-пучка (центр + кольцо 2мм) в зоне поля — критерий апертуры."""
    _, _, maxX, _ = seb._get(devs[0]); ff = seb._field_multi(devs, currents)
    xlo = seb.DRIFT - maxX - 0.02; xhi = seb.DRIFT + maxX + 0.02
    import math
    pts = [(0.0, 0.0)] + [(2e-3 * math.cos(a), 2e-3 * math.sin(a)) for a in (0, math.pi/2, math.pi, 3*math.pi/2)]
    return max(seb._track_maxr(y0, z0, ff, xlo, xhi) for y0, z0 in pts)


def working_range(dev, aper, cap=60.0):
    """Рабочий диапазон тока обмотки: наибольший |I|, при котором одиночная обмотка держит пучок
    в 0.5·апертуры (как sweep_for). Гарантирует, что комбо у этого предела ещё в карте."""
    v = cap
    for _ in range(40):
        if excursion_multi([dev], [v], aper) <= 0.5 * aper or v < 1e-3:
            break
        v *= 0.75
    return v


def sample_combos(devs, ranges, aper, n=24, seed=0):
    """Латинский гиперкуб в [-r,r] по обмоткам + одиночные якоря ±r; вне-апертурные отброшены."""
    rng = np.random.default_rng(seed); N = len(devs)
    grid = np.linspace(-1, 1, n)                                  # простой LHS: перестановка сетки на обмотку
    cols = np.stack([rng.permutation(grid) for _ in range(N)], axis=1)  # (n, N) в [-1,1]
    combos = [cols[i] * np.array(ranges) for i in range(n)]
    for j in range(N):                                            # якоря: одиночная обмотка ±r
        for s in (-1, 1):
            c = np.zeros(N); c[j] = s * ranges[j]; combos.append(c)
    valid = [c for c in combos if excursion_multi(devs, c, aper) <= 0.9 * aper and np.any(np.abs(c) > 1e-9)]
    return valid


def fit_strengths(init, diff_final, L, aper, nsl=40, iters=12):
    """Подгонка 8 сил (θy,θz,k1,k2,k3,ks1,ks2,ks3) под финальные УГЛЫ диффур-прогона (Gauss-Newton).
    Углы дрейф-инвариантны ⇒ d1=d2=0. init — [(Y0,Z0)], diff_final — [(Y,Yp,Z,Zp)]."""
    target = np.array([[f[1], f[3]] for f in diff_final]).ravel()   # (2·npart,) целевые Yp,Zp
    def model(s):
        out = ue.track(init, s[0], s[1], s[2], s[3], s[4], L, 0.0, 0.0, nsl=nsl,
                       ks1=s[5], ks2=s[6], ks3=s[7], aper=aper)
        return np.array([[o[1], o[3]] for o in out]).ravel()
    s = np.zeros(8)
    for _ in range(iters):
        r = model(s) - target
        J = np.zeros((len(r), 8))
        for k in range(8):
            ds = max(1e-6, abs(s[k]) * 1e-3); sp = s.copy(); sp[k] += ds
            J[:, k] = (model(sp) - (r + target)) / ds
        try:
            delta, *_ = np.linalg.lstsq(J, -r, rcond=None)
        except np.linalg.LinAlgError:
            break
        s = s + delta
        if np.linalg.norm(delta) < 1e-12:
            break
    return s


def calibrate(name, n=24, verbose=True):
    devs = GROUPS[name]; N = len(devs)
    L, aper = elem_geom(devs[0])
    ranges = [working_range(d, aper) for d in devs]
    if verbose:
        print(f"[{name}] обмоток {N}, L={L*1e3:.0f}мм, апертура {aper*1e3:.0f}мм; "
              f"рабочие диапазоны токов: {[round(r,1) for r in ranges]} А")
    combos = sample_combos(devs, ranges, aper, n=n)
    if verbose:
        print(f"   комбинаций в апертуре: {len(combos)}")
    Imat = []; Smat = []
    for c in combos:
        ref = seb.run_multi(devs, list(c))
        s = fit_strengths(ref["init"], ref["final"], L, aper)
        Imat.append(c); Smat.append(s)
    Imat = np.array(Imat); Smat = np.array(Smat)               # (ncombo, N), (ncombo, 8)
    # линейная матрица связи: силы = токи · M  (M: N×8), по всем комбинациям сразу
    M, res, rank, sv = np.linalg.lstsq(Imat, Smat, rcond=None)  # M[w, strength]
    if verbose:                                                # ЛИНЕЙНОСТЬ: R² подгонки каждой силы
        pred = Imat @ M
        print("   R² линейной матрицы связи (силы = токи·M; <1 ⇒ нелинейно, матрица недостаточна):")
        for j, nm in enumerate(STR_NAMES):
            sr = ((Smat[:, j] - pred[:, j]) ** 2).sum(); st = ((Smat[:, j] - Smat[:, j].mean()) ** 2).sum()
            r2 = 1 - sr / st if st > 1e-30 else 1.0
            if abs(Smat[:, j]).max() > 1e-9:
                print(f"      {nm:4}: R²={r2:+.4f}  размах силы {Smat[:,j].min():+.3g}..{Smat[:,j].max():+.3g}")
    # физкоэф. на 1А обмотки: θ→per_B=−θ/(QP·L) ... но per_* нужны в конвенции coupling.
    # phys_from_strength даёт на ЕДИНИЧНУЮ силу; здесь сила уже на 1А ⇒ применяем ту же формулу.
    # силы → per_* В КОНВЕНЦИИ МАТРИЦЫ СВЯЗИ (θ0=QP·per_B·L_eff, k1=QP·per_G, …):
    #   per_B=θy/(QP·L), per_Bs=θz/(QP·L), per_G=k1/QP, per_S=k2/QP, per_O=k3/QP, скью аналогично.
    # (НЕ phys_from_strength.Bz0=−θy/(QP·L) — там др. знак «физ. поля», не совпадает с coupling.)
    rows = {}
    for w, dev in enumerate(devs):
        s = M[w]
        rows[dev] = {"per_B_Tl": s[0] / (QP * L), "per_Bs_Tl": s[1] / (QP * L),
                     "per_G_Tlm": s[2] / QP, "per_Gs_Tlm": s[5] / QP,
                     "per_S_Tlm2": s[3] / QP, "per_Ss_Tlm2": s[6] / QP,
                     "per_O_Tlm3": s[4] / QP, "per_Os_Tlm3": s[7] / QP}
    return rows, Imat, Smat, M, L


def _compare(name, rows):
    """Сравнить комбо-калибровку с текущей coupling_matrix (главные мультиполи)."""
    cur = {r["device_type_id"]: r for r in csv.DictReader(
        open(os.path.join(TB, "coupling_matrix.csv"), encoding="utf-8"), delimiter=";")}
    print(f"\n[{name}] комбо-калибровка vs текущая (per-амперные коэф., комбо/текущ.):")
    for col in ("per_B_Tl", "per_G_Tlm", "per_S_Tlm2"):
        print(f"   {col}:")
        for dev, r in rows.items():
            c = float(cur.get(dev, {}).get(col, 0))
            rel = (abs(r[col] - c) / abs(c) * 100) if abs(c) > 1e-12 else 0
            print(f"      {dev:14} {r[col]:+13.4e} / {c:+13.4e}   Δ={rel:5.0f}%")


def main():
    args = sys.argv[1:]
    write = "--write" in args
    names = list(GROUPS) if "--all" in args else [a for a in args if a in GROUPS]
    if not names:
        names = ["DQ"]
    allrows = {}
    for nm in names:
        rows, *_ = calibrate(nm)
        _compare(nm, rows)
        allrows.update(rows)
    if write:
        _write(allrows)


def _write(allrows):
    # распространить на СИНТЕТИЧЕСКИЕ элементы (per-амперные коэф. длино-независимы): scaled_from ∈ allrows
    hw = {r["device_type_id"]: r for r in csv.DictReader(
        open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";")}
    for dev, r in hw.items():
        if r.get("synthetic") == "True" and r.get("scaled_from") in allrows and dev not in allrows:
            allrows[dev] = allrows[r["scaled_from"]]
    path = os.path.join(TB, "coupling_matrix.csv")
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";")); head = rows[0]
    idx = {c: head.index(c) for c in head}
    for r in rows[1:]:
        d = r[0]
        if d in allrows:
            for col, val in allrows[d].items():
                r[idx[col]] = f"{val:.6e}"
            r[idx["notes"]] = "комбо-калибровка (диффур, множ. комбинаций)"
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    print(f"\n[OK] coupling_matrix.csv обновлён для {len(allrows)} обмоток.")


if __name__ == "__main__":
    main()
