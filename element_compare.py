# -*- coding: utf-8 -*-
"""Сравнение ДИФФУРНОЙ (C++ SynchrotronTracker, Boris по картам поля) и КИКОВОЙ (Python matrix_tracker в
режиме SS_SECTOR — симплектический drift-kick-drift, НЕ передаточные матрицы) моделей по одиночным элементам
БД. Датасет траекторий при разных токах + подгонка кика (L_eff+физкоэф).

ФАЗА 1 (этот файл пока): фундамент — собрать ВРЕМЕННЫЙ одно-элементный конфиг в основной БД
(пучок → элемент → финиш-плоскость; дрейфы = свободный пролёт, длины НАШИ, из диффура не берём),
прогнать C++-диффур, прочитать траекторию, вычистить временные строки. Дальше — свип токов, фит кика,
Word-отчёт. Диффур eff_length НЕ использует (трек по карте поля) ⇒ спека «без диффурных длин» соблюдена.
"""
import os, csv, math, struct, subprocess, shutil, sys, itertools
import numpy as np
import matrix_tracker as mt
import bip_matrix as bm
from scipy.optimize import minimize, minimize_scalar

ROOT = os.path.dirname(os.path.abspath(__file__))
TB = os.path.join(ROOT, "workspace", "tables")
RESULTS = os.path.join(ROOT, "workspace", "results")
EXE = os.path.join(ROOT, "build", "SynchrotronTracker.exe")
V_BEAM = 2.967e8            # скорость пучка (particle_groups)
ME, QE = 9.109e-31, -1.602e-19

# временные id (не пересекаются с реальными 1/2/3); всё под ними чистится
TMP_LAT, TMP_GRP, TMP_EXP = 900, 900, 900
CENTER = 2.0               # положение центра элемента по оси x, м
ENTRY, EXIT = 1.5, 1.5     # БОЛЬШОЙ фиксированный дрейф до/после, м (за полем любого элемента; НЕ из диффура)

CSV4 = ["lattice_configs.csv", "element_connectivity.csv", "particle_groups.csv", "experiments.csv"]


def _rows(name):
    p = os.path.join(TB, name)
    with open(p, encoding="utf-8") as f:
        r = list(csv.reader(f, delimiter=";"))
    return r[0], r[1:]


def _write(name, header, rows):
    with open(os.path.join(TB, name), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(header); w.writerows(rows)


def make_beam(x0, radius=1e-3, n_circle=12, grp_id=TMP_GRP):
    """13 частиц у стартовой плоскости x=x0: 12 по кругу Ø2·radius в поперечной (y,z) + 1 в центре.
    Все летят вдоль +x со скоростью пучка. Возврат список строк particle_groups."""
    parts = [(0.0, 0.0)]
    for k in range(n_circle):
        a = 2 * math.pi * k / n_circle
        parts.append((radius * math.cos(a), radius * math.sin(a)))
    rows = []
    for pid, (y, z) in enumerate(parts, start=1):
        rows.append([str(grp_id), str(pid), repr(ME), repr(QE),
                     repr(x0), repr(y), repr(z), repr(V_BEAM), "0.0", "0.0"])
    return rows


def build_config(dev, cur, yaw=0.0):
    """Одно-элементный конфиг (dev при токе cur). Обёртка над build_windings (одна обмотка)."""
    return build_windings([(dev, cur)], yaw)


def build_windings(windings, yaw=0.0, lat_id=TMP_LAT, grp_id=TMP_GRP, exp_id=TMP_EXP):
    """Записать конфиг из N СО-ЛОКОВАННЫХ обмоток (все в CENTER, ось +x, поворот yaw) в 4 CSV под id
    lat/grp/exp. windings = [(dev, cur|'-'), …]. Диффур суперпозирует поля со-локованных элементов
    (getTotalFields). id по умолчанию TMP (900, чистится cleanup); свои id (напр. 950) — ПОСТОЯННЫЙ конфиг
    (для visualize.py), cleanup не трогает. Возврат dict с геометрией."""
    lat_h, lat_r = _rows("lattice_configs.csv")
    con_h, con_r = _rows("element_connectivity.csv")
    grp_h, grp_r = _rows("particle_groups.csv")
    exp_h, exp_r = _rows("experiments.csv")
    ci = {c: i for i, c in enumerate(lat_h)}
    coni = {c: i for i, c in enumerate(con_h)}
    lat_rows, con_rows = [], []
    for k, (dev, cur) in enumerate(windings, start=1):     # instance_id 1..N, все в одной точке CENTER
        row = [""] * len(lat_h)
        row[ci["lattice_id"]] = str(lat_id); row[ci["instance_id"]] = str(k); row[ci["device_type_id"]] = dev
        row[ci["x"]] = repr(CENTER); row[ci["y"]] = "0.0"; row[ci["z"]] = "0.0"
        row[ci["yaw"]] = repr(yaw); row[ci["pitch"]] = "0.0"; row[ci["roll"]] = "0.0"
        row[ci["arg_val"]] = cur if cur == "-" else repr(float(cur)); row[ci["correct_coeff"]] = "1.0"  # float(): np.float64 repr ломает C++ stod
        for c in ("mat_B", "mat_G", "mat_S", "mat_O"):
            if c in ci: row[ci[c]] = "0"
        lat_rows.append(row)
        crow = [""] * len(con_h)
        crow[coni["lattice_id"]] = str(lat_id); crow[coni["instance_id"]] = str(k); crow[coni["device_type_id"]] = dev
        con_rows.append(crow)
    lat_new = [r for r in lat_r if r[0] != str(lat_id)] + lat_rows
    con_new = [r for r in con_r if r[0] != str(lat_id)] + con_rows
    # пучок у стартовой плоскости
    grp_new = [r for r in grp_r if r[0] != str(grp_id)] + make_beam(CENTER - ENTRY, grp_id=grp_id)
    # опыт: финиш-плоскость на x=CENTER+EXIT, нормаль +x (3 точки + направление)
    xf = CENTER + EXIT
    p1 = (xf, 1.0, 1.0); p2 = (xf, -1.0, 1.0); p3 = (xf, 0.0, -1.0)
    erow = [str(exp_id), str(lat_id), str(grp_id), "1.00E-11", "0", "1",
            f"{p1[0]:.6f}", f"{p1[1]:.6f}", f"{p1[2]:.6f}",
            f"{p2[0]:.6f}", f"{p2[1]:.6f}", f"{p2[2]:.6f}",
            f"{p3[0]:.6f}", f"{p3[1]:.6f}", f"{p3[2]:.6f}", "1.0", "0.0", "0.0"]
    erow += [""] * (len(exp_h) - len(erow))
    exp_new = [r for r in exp_r if r and r[0] != str(exp_id)] + [erow]
    _write("lattice_configs.csv", lat_h, lat_new)
    _write("element_connectivity.csv", con_h, con_new)
    _write("particle_groups.csv", grp_h, grp_new)
    _write("experiments.csv", exp_h, exp_new)
    return dict(center=CENTER, entry=ENTRY, exit=EXIT, xstart=CENTER - ENTRY, xfinish=xf)


def read_bin(exp_id):
    """results/trajectories_exp{id}.bin → список массивов (N,3) X,Y,Z на частицу (формат C++)."""
    p = os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin")
    out = []
    with open(p, "rb") as f:
        data = f.read()
    off = 0
    while off < len(data):
        (n,) = struct.unpack_from("I", data, off); off += 4
        a = np.frombuffer(data, np.float64, count=3 * n, offset=off).reshape(n, 3); off += 3 * n * 8
        out.append(a.copy())
    return out


def run_diffur(dev, cur, yaw=0.0):
    """Собрать конфиг, прогнать C++-диффур, вернуть траектории (список (N,3)). Конфиг остаётся в БД —
    чистит вызывающий (cleanup). Возврат (trajs, geom)."""
    geom = build_config(dev, cur, yaw)
    r = subprocess.run([EXE, str(TMP_EXP)], capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        raise RuntimeError(f"диффур упал: {r.stderr[-400:]}")
    return read_bin(TMP_EXP), geom


def cleanup():
    """Убрать все временные строки (TMP_*) из 4 CSV."""
    for name, key in (("lattice_configs.csv", TMP_LAT), ("element_connectivity.csv", TMP_LAT),
                      ("particle_groups.csv", TMP_GRP), ("experiments.csv", TMP_EXP)):
        h, rows = _rows(name)
        _write(name, h, [r for r in rows if not (r and r[0] == str(key))])
    for f in (f"trajectories_exp{TMP_EXP}.bin", f"debug_exp{TMP_EXP}.csv"):
        p = os.path.join(RESULTS, f)
        if os.path.exists(p): os.remove(p)


# ============================ МАТРИЧНОЕ ЯДРО (Python) + ФИТ ==================================
os.environ.setdefault("SS_SECTOR", "1")                # секторный/кик режим ядра (drift-kick-drift, НЕ матрицы)
EFF, PHYS, COUP, TRANS = mt.load_db()
HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
LCROW = {r["device_type_id"]: r for r in mt._read("lattice_configs.csv") if r.get("lattice_id") == "1"}
DELTA = (bm.GAMMA * 9.109e-31 * V_BEAM - bm.P_MOM) / bm.P_MOM   # δ пучка (усечённая масса, как диффур)


def is_pm(dev):
    """Постоянный магнит (не токовый): nominal_arg='-' ⇒ диффур ставит currentScale=1.0 (поле как снято)."""
    return (HW.get(dev, {}).get("nominal_arg") or "").strip() == "-"
KIND_COEFF = {"quad": ["G"], "sext": ["S"], "oct": ["O"], "dipole": ["B", "G"], "dq": ["B", "G"]}
SGRID = np.linspace(0.0, ENTRY + EXIT, 240)            # общая сетка пути s для сравнения


def build_matrix_magnet(dev, coeffs, eff_length, kind):
    """Magnet matrix_tracker (source=physical) с заданными физкоэф. и ПОДБИРАЕМОЙ eff_length."""
    g = lambda k: repr(float(coeffs.get(k, 0.0)))      # Python float (repr(np.float64) ломает matrix_tracker)
    inst = {"device_type_id": dev, "x": "0", "y": "0", "z": "0", "yaw": "0", "pitch": "0", "roll": "0",
            "instance_id": "1", "arg_val": "0", "correct_coeff": "1",
            "mat_B": g("B"), "mat_G": g("G"), "mat_Gv": g("G"), "mat_S": g("S"), "mat_O": g("O"),
            "eff_length": repr(float(eff_length)), "kind": kind, "n_slices": "5"}
    return mt.Magnet([inst], EFF, PHYS, COUP, TRANS, source="physical")


def track_matrix(dev, coeffs, eff_length, kind, geom, probe, M=24):
    """Трек одной частицы (u,v) через дрейф→элемент(M подэлементов для гладкой траектории)→дрейф.
    ПРЯМОЙ элемент (ось +x): локальные x=u, y=v ≡ поперечные глоб. Возврат (K,3): (s, u, v)."""
    u0, v0 = probe
    st = (u0, 0.0, v0, 0.0, 0.0, DELTA)
    di = geom["entry"] - eff_length / 2.0; do = geom["exit"] - eff_length / 2.0
    pts = [(0.0, u0, v0)]; s = 0.0

    def drift_rec(st, s, L, n=8):
        for _ in range(n):
            st = mt._drift(st, L / n); s += L / n; pts.append((s, st[0], st[2]))
        return st, s
    st, s = drift_rec(st, s, di)
    sub = build_matrix_magnet(dev, coeffs, eff_length / M, kind)
    for _ in range(M):
        st = sub.body(st); s += eff_length / M; pts.append((s, st[0], st[2]))
    st, s = drift_rec(st, s, do)
    return np.array(pts)


def _resample(traj):
    """(K,3)=(s,u,v) → (u,v) на общей сетке SGRID (интерполяция)."""
    return np.interp(SGRID, traj[:, 0], traj[:, 1]), np.interp(SGRID, traj[:, 0], traj[:, 2])


def track_set(dev, coeffs, L, kind, geom, probes):
    """Матрица для ВСЕХ проб (K,3)=(s,u,v). Для ИЗГИБА возвращает отклонение от МАТРИЧНОЙ опорной (probe 0)
    на сетке SGRID — тот же смысл, что Френе-проекция диффура от ч0. Иначе опорная уходит от прямой оси на
    θ·L_дрейф (~8.7мм у DP_L), т.к. track_matrix меряет x от прямой оси, а не от изогнутой дуги. Для прямых
    (опорная ≡ 0) — как есть."""
    mats = [track_matrix(dev, coeffs, L, kind, geom, p) for p in probes]
    if kind not in ("dipole", "dq"):
        return mats
    ur, vr = _resample(mats[0])                            # киковая опорная (ч0) на SGRID
    out = []
    for m in mats:
        um, vm = _resample(m)
        out.append(np.column_stack([SGRID, um - ur, vm - vr]))
    return out


def diff_local(trajs, xstart, bending=False):
    """Диффур .bin (X,Y,Z глоб.) → на частицу (s, u, v). ПРЯМОЙ: u=Y, v=Z, s=X−xstart (опорная на оси).
    ИЗГИБ: проекция ФРЕНЕ на опорную дугу ч0 — s=ДЛИНА ДУГИ вдоль ч0, u=⊥смещение В ПЛОСКОСТИ ИЗГИБА
    (кольцо X-Y), v=вертикаль (Z−Z_ч0). Именно так работает киковый сектор (отклонение от изогнутой
    опорной по длине дуги). Прежний вычет «при том же X» для изогнутой орбиты давал спурьёзный член
    (две параллельные дуги одного радиуса при проекции на X расходятся) — артефакт ~мм на чистом изгибе."""
    if not bending:
        return [np.column_stack([a[:, 0] - xstart, a[:, 1], a[:, 2]]) for a in trajs]
    ref = trajs[0]; Xr, Yr, Zr = ref[:, 0], ref[:, 1], ref[:, 2]
    Sref = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(Xr), np.diff(Yr)))])   # длина дуги ч0 (изгиб в X-Y)
    # единичная касательная ч0 в каждой вершине (центр.разность) и ⊥нормаль в плоскости X-Y
    tx = np.gradient(Xr); ty = np.gradient(Yr); tn = np.hypot(tx, ty); tn[tn == 0] = 1.0
    tx /= tn; ty /= tn; nx, ny = -ty, tx
    out = []
    for a in trajs:
        X, Y, Z = a[:, 0], a[:, 1], a[:, 2]
        # ближайшая вершина ч0 к каждой точке частицы (в плоскости X-Y)
        i = np.argmin((Xr[None, :] - X[:, None]) ** 2 + (Yr[None, :] - Y[:, None]) ** 2, axis=1)
        u = (X - Xr[i]) * nx[i] + (Y - Yr[i]) * ny[i]     # ⊥смещение в плоскости изгиба
        v = Z - Zr[i]                                      # вертикаль
        out.append(np.column_stack([Sref[i], u, v]))
    return out


diff_to_su = diff_local                                              # обратная совместимость (прямые элементы)


def traj_mse(diff_su, mat_trajs, smax=None):
    """MSE по частицам и точкам сетки: (u_d−u_m)²+(v_d−v_m)², м². smax — ограничить сравнение окном s≤smax
    (для ФИТА: элемент + небольшой хвост, без дальнего дрейфа, где угл.ошибка сильно усиливается)."""
    mask = (SGRID <= smax) if smax is not None else np.ones_like(SGRID, bool)
    tot = 0.0; n = 0
    for d, m in zip(diff_su, mat_trajs):
        ud, vd = _resample(d); um, vm = _resample(m)
        tot += np.sum((ud[mask] - um[mask]) ** 2 + (vd[mask] - vm[mask]) ** 2); n += 2 * int(mask.sum())
    return tot / n


def _exit_slope(traj, comp):
    """Наклон компоненты (1=u,2=v) в ВЫХОДНОМ дрейфе (последние 0.3м — прямая)."""
    s = traj[:, 0]; m = s > (s[-1] - 0.3)
    return np.polyfit(s[m], traj[m, comp], 1)[0] if m.sum() >= 2 else 0.0


def _estimate_init(dev, diff_su, kind, geom, probes, L0, theta0=0.0):
    """Физический начальный guess из диффура. Прямой: доминантный коэф. из наклона отклонения частицы 1.
    Изгиб: B из полного изгиба опорной (theta0=2·yaw; отклонение ч0 от себя ≈0, потому берём угол), у dq
    + G из наклона отклонения частицы. Коэф = kL/(QP·L0)."""
    QP = mt.QP; u0 = probes[1][0]                       # частица 1 = (u=1мм, v=0)
    if kind in ("dipole", "dq"):
        est = {"B": theta0 / (QP * L0)} if QP else {}
        if kind == "dq" and QP:                        # + градиент из наклона отклонения частицы
            est["G"] = (-_exit_slope(diff_su[1], 1) / u0) / (QP * L0)
        return est
    du = _exit_slope(diff_su[1], 1)
    kL = {"quad": -du / u0, "sext": -2 * du / u0 ** 2, "oct": -6 * du / u0 ** 3}.get(kind, 0.0)
    key = {"quad": "G", "sext": "S", "oct": "O"}.get(kind)
    return {key: kL / (QP * L0)} if (key and QP) else {}


# физический диапазон КОАРС-СКАНА доминантного коэф. по типу (робастный старт вместо аналит. оценки)
SCAN_RANGE = {"quad": ("G", 3.0), "sext": ("S", 400.0), "oct": ("O", 2e4),
              "dipole": ("G", 2.0), "dq": ("G", 5.0)}   # изгибающие: B из угла, скан по G (градиент)


def fit_element(dev, cur, diff_su, geom, kind, probes, theta0=0.0):
    """Последовательное приближение: подобрать физкоэф (по kind) + L_eff, минимизируя MSE траекторий.
    L_eff — ОДИН на элемент; физкоэф — на этот ток. theta0 — изгиб опорной (для B у изгибающих).
    Старт доминантного коэф. — КОАРС-СКАН по физ.диапазону (робастно, симметрично; аналит.оценка из наклона
    ненадёжна). Возврат (coeffs, Leff, mse)."""
    keys = KIND_COEFF[kind]
    L0 = float(HW[dev].get("eff_length") or 0.3)       # старт L_eff — .lte длина модели (НЕ диффурная)
    base = {"B": theta0 / (mt.QP * L0)} if (kind in ("dipole", "dq") and mt.QP) else {}   # B изгиба фиксирован из угла
    key, span = SCAN_RANGE[kind]
    best_v, best_e = 0.0, 1e18                          # коарс-скан доминантного коэф. (61 точка ±span)
    for val in np.linspace(-span, span, 61):
        cf = dict(base); cf[key] = float(val)
        try:
            mat = track_set(dev, cf, L0, kind, geom, probes); e = traj_mse(diff_su, mat)
        except mt.ParticleLost:
            continue
        if e < best_e:
            best_e, best_v = e, float(val)
    init = dict(base); init[key] = best_v
    x0 = [init.get(k, 0.0) for k in keys] + [L0]

    def unpack(x):
        return {k: x[i] for i, k in enumerate(keys)}, x[-1]

    def obj(x):
        coeffs, Leff = unpack(x)
        if not (0.6 * L0 <= Leff <= 1.4 * L0):          # L_eff в узких границах вокруг .lte (не схлопывать/раздувать)
            return 1e3
        try:
            mat = track_set(dev, coeffs, Leff, kind, geom, probes)
        except mt.ParticleLost:
            return 1e3
        return traj_mse(diff_su, mat)                   # полный путь: дальний дрейф констрейнит выходной угол

    res = minimize(obj, x0, method="Nelder-Mead",
                   options=dict(xatol=1e-8, fatol=1e-20, maxiter=12000, maxfev=12000))
    coeffs, Leff = unpack(res.x)
    try:                                                # ПОЛНАЯ невязка (весь путь) с фитированными коэф. — для отчёта
        mat = track_set(dev, coeffs, Leff, kind, geom, probes)
        full_mse = traj_mse(diff_su, mat)
    except mt.ParticleLost:
        full_mse = res.fun
    return coeffs, Leff, full_mse


def _coeff_at_L(dev, kind, geom, probes, diff_su, L, theta0=0.0):
    """1D-фит ДОМИНАНТНОГО физкоэф при ФИКСИРОВАННОМ L_eff (коарс-скан 41 точка + brent-refine). B изгиба —
    из угла (theta0), не фитится. Возврат (coeffs, mse). Ядро общего-L_eff фита: L задаётся снаружи."""
    base = {"B": theta0 / (mt.QP * L)} if (kind in ("dipole", "dq") and mt.QP) else {}
    key, span = SCAN_RANGE[kind]

    def m(val):
        cf = dict(base); cf[key] = float(val)
        try:
            return traj_mse(diff_su, track_set(dev, cf, L, kind, geom, probes))
        except mt.ParticleLost:
            return 1e3

    grid = np.linspace(-span, span, 41)
    ev = [m(v) for v in grid]
    v0 = float(grid[int(np.argmin(ev))]); e0 = min(ev)
    step = span / 20.0                                  # локальный refine вокруг лучшего узла скана
    res = minimize_scalar(m, bounds=(v0 - step, v0 + step), method="bounded", options=dict(xatol=step * 1e-4))
    val, err = (float(res.x), float(res.fun)) if res.fun < e0 else (v0, e0)
    coeffs = dict(base); coeffs[key] = val
    return coeffs, err


def fit_shared(dev, kind, probes, runs):
    """ОДИН общий L_eff (геом. константа элемента, НЕ зависит от тока) + физкоэф на КАЖДЫЙ ток. Внешняя
    1D-оптимизация L_eff по СУММАРНОЙ невязке всех токов; внутри — доминантный коэф. при фикс. L_eff. runs —
    список dict(cur, diff_su, geom, theta0). Заполняет r['coeffs','Leff','rms_mm']; возврат Lshared."""
    L0 = float(HW[dev].get("eff_length") or 0.3)
    if is_pm(dev):                                      # PM/чистый изгиб: L_eff ВЫРОЖДЕН (важно только G·L_eff,
        Lshared = L0                                    # частицы едут параллельно) — не фитим, берём .lte
    else:
        strong = [r for r in runs if abs(r["cur"]) >= 30] or runs   # сильные токи констрейнят L_eff надёжнее

        def total(L):                                   # суммарная невязка по сильным токам при общем L_eff
            return sum(_coeff_at_L(dev, kind, r["geom"], probes, r["diff_su"], L, r["theta0"])[1] for r in strong)

        res = minimize_scalar(total, bounds=(0.6 * L0, 1.4 * L0), method="bounded", options=dict(xatol=L0 * 1e-4))
        Lshared = float(res.x)
    for r in runs:                                      # физкоэф каждого тока при ОБЩЕМ L_eff
        coeffs, mse = _coeff_at_L(dev, kind, r["geom"], probes, r["diff_su"], Lshared, r["theta0"])
        r["coeffs"] = coeffs; r["Leff"] = Lshared; r["rms_mm"] = math.sqrt(mse) * 1e3
    return Lshared


# ============================ ИЗГИБАЮЩИЕ: ПОДБОР ПОВОРОТА ===================================
def tune_bend_yaw(dev, cur, iters=6, tol=1e-6):
    """Итеративно подобрать yaw элемента, чтобы ОПОРНАЯ орбита (ч0, на оси) входила и выходила под
    симметричным углом (угол входа отн. оси элемента = −угол выхода). Пучок летит +x; при повороте yaw
    вход отн. элемента = −yaw, выход = φ−yaw (φ — угол выхода в лаб.). Симметрия ⇒ yaw = φ/2. Сходится
    за 1-2 итерации. Каждый шаг — прогон диффура (только ч0 важна), очистка. Возврат yaw (рад)."""
    yaw = 0.0
    for _ in range(iters):
        trajs, _ = run_diffur(dev, cur, yaw)
        ref = trajs[0]                                     # ч0 = опорная орбита
        s = ref[:, 0]; m = s > (s[-1] - 0.3)              # выходной дрейф — прямая
        phi = math.atan(np.polyfit(ref[m, 0], ref[m, 1], 1)[0])   # угол выхода в лаб. (dY/dX)
        cleanup()
        new = phi / 2.0
        if abs(new - yaw) < tol:
            yaw = new; break
        yaw = new
    return yaw


# ============================ СВИП ТОКОВ + ДАТАСЕТ ==========================================
CUR_RANGE = list(range(-50, 51, 10))                   # −50..+50 шаг 10 (прямые токовые)


def dq_operating_current(dev):
    """Рабочий ток DQ: где изгиб КАРТЫ = проектному (θ0=QP·mat_B·L). Наклон изгиб/ток — из одного пробного
    прогона. У DQ карта согласована с mat_B (в отличие от PM-диполей), рабочая точка ~несколько А."""
    r = LCROW[dev]; B = float(r.get("mat_B") or 0); L = float(r.get("eff_length") or 1.0)
    th_design = mt.QP * B * L
    Itest = 5.0; yaw = tune_bend_yaw(dev, Itest); cleanup()
    slope = (2.0 * yaw) / Itest                          # изгиб на ампер (рад/А)
    return abs(th_design / slope) if abs(slope) > 1e-12 else Itest


def element_currents(dev, kind):
    """Токи свипа по типу. PM-диполь → одна конфигурация (поле фиксировано). DQ → вокруг РАБОЧЕГО тока
    (±50А = перегруз ×17, пучок раздувается). Прямой токовый → ±50А."""
    if is_pm(dev):
        return [0.0]                                     # одна конфигурация (arg_val игнор — элемент неактивен)
    if kind == "dq":
        Iop = dq_operating_current(dev)                  # ~несколько ампер
        return [round(f * Iop, 3) for f in (0.4, 0.6, 0.8, 1.0, 1.2, 1.4)]
    return [float(c) for c in CUR_RANGE]


def sweep_currents(dev, kind, probes):
    """Свип токов: диффур на каждый ток (сбор датасета), затем ОДИН общий L_eff + физкоэф на каждый ток
    (fit_shared — L_eff геом.константа, одна на все токи). Изгибающие — с подбором yaw (симметрия) и
    сравнением по ОТКЛОНЕНИЮ от опорной. Токи — element_currents (PM=1 конфиг, DQ=вокруг рабочей точки).
    Возврат записей (Leff у всех = общий)."""
    bending = kind in ("dipole", "dq")
    runs = []
    for cur in element_currents(dev, kind):
        yaw = tune_bend_yaw(dev, cur) if bending else 0.0     # изгиб: симметричный поворот
        try:
            trajs, geom = run_diffur(dev, cur, yaw)
        except Exception as e:
            print(f"  ток {cur:+d}: диффур упал ({e})"); cleanup(); continue
        diff_su = diff_local(trajs, geom["xstart"], bending)
        runs.append(dict(cur=cur, diff_su=diff_su, geom=geom, yaw=yaw, theta0=2.0 * yaw))
        cleanup()
    if not runs:
        return []
    Lshared = fit_shared(dev, kind, probes, runs)             # ОДИН L_eff на элемент (общий для всех токов)
    for r in runs:
        r["Leff_common"] = Lshared
        cs = " ".join(f"{k}={v:+.3e}" for k, v in r["coeffs"].items())
        print(f"  ток {r['cur']:+7.2f}А: {cs}  L_eff(общий)={Lshared:.4f}  RMS={r['rms_mm']:.4f}мм")
    return runs


def save_dataset(dev, kind, recs, outdir):
    """Датасет current↔физкоэф+MSE в CSV (запятая-десятичная для Excel, ; поле)."""
    keys = KIND_COEFF[kind]
    Lc = recs[0].get("Leff_common", recs[0]["Leff"])
    path = os.path.join(outdir, f"dataset_{dev}.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["ток_А"] + [f"{k}" for k in keys] + ["RMS_мм"])   # физкоэф на ток; L_eff — общий (ниже)
        for r in recs:
            row = [r["cur"]] + [r["coeffs"][k] for k in keys] + [r["rms_mm"]]
            w.writerow([str(x).replace(".", ",") for x in row])
        w.writerow([]); w.writerow(["L_eff_общий_м", str(Lc).replace(".", ",")])
    return path


# ============================ ГРАФИКИ + WORD-ОТЧЁТ ==========================================
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = os.path.join(RESULTS, "element_compare")
DATASETS = os.path.join(OUTDIR, "datasets")        # свип- и комбо-датасеты матриц связи
HTML_OUT = os.path.join(RESULTS, "html")           # HTML-отчёты элементов (viz_html)


def _show_currents(recs):
    """Показательные токи для графиков: 2 шт БЕЗ нулевого (ток 0 = чистый дрейф, неинформативен). Знаковый
    свип (прямые) → ток ближайший к ±50% макс. Положительный свип (DQ) → 50% и 100% макс. PM/малый — как есть."""
    curs = sorted({r["cur"] for r in recs})
    if len(curs) <= 2:
        return curs
    amax = max(abs(c) for c in curs); half = 0.5 * amax
    nearest = lambda t: min(curs, key=lambda c: (abs(c - t), -abs(c)))   # ближайший к t, тай → больший модуль
    return [nearest(-half), nearest(half)] if any(c < 0 for c in curs) else [nearest(half), nearest(amax)]


def _mat_for_plot(dev, coeffs, L, kind, geom, probes):
    """Матрица по пробам для графика: (s,u,v) на пробу, для ИЗГИБА — отклонение от МАТРИЧНОЙ опорной (probe 0).
    None для потерянных частиц (частица за пределом параксиальности). Возврат список len(probes)."""
    bending = kind in ("dipole", "dq")
    ref_uv = None
    if bending:
        try:
            ref_uv = _resample(track_matrix(dev, coeffs, L, kind, geom, probes[0]))   # опорная обязана выжить
        except mt.ParticleLost:
            return [None] * len(probes)
    out = []
    for p in probes:
        try:
            m = track_matrix(dev, coeffs, L, kind, geom, p)
        except mt.ParticleLost:
            out.append(None); continue
        if bending:
            um, vm = _resample(m); out.append(np.column_stack([SGRID, um - ref_uv[0], vm - ref_uv[1]]))
        else:
            out.append(m)
    return out


def plot_element(dev, kind, recs, probes, outdir):
    """Графики траекторий x/y (диффур сплошная синяя vs кик пунктир красная, 13 частиц) при показательных
    токах. Кик — с фитированными физкоэф и ОБЩИМ L_eff. Изгиб: отклонение от опорной. Возврат [(cur, png)]."""
    Lc = recs[0].get("Leff_common", recs[0]["Leff"]); figs = []
    for cur in _show_currents(recs):
        r = next((x for x in recs if x["cur"] == cur), None)
        if r is None: continue
        fig, (axx, axy) = plt.subplots(1, 2, figsize=(11, 4.2))
        mats = _mat_for_plot(dev, r["coeffs"], Lc, kind, r["geom"], probes)
        for d, mat in zip(r["diff_su"], mats):
            axx.plot(d[:, 0], d[:, 1] * 1e3, color="#2b7bba", lw=1.1)   # диффур всегда
            axy.plot(d[:, 0], d[:, 2] * 1e3, color="#2b7bba", lw=1.1)
            if mat is None: continue                                    # кик разошёлся — только диффур
            axx.plot(mat[:, 0], mat[:, 1] * 1e3, color="#d1495b", lw=1.0, ls="--")
            axy.plot(mat[:, 0], mat[:, 2] * 1e3, color="#d1495b", lw=1.0, ls="--")
        ttl2 = " (откл. от опорной)" if kind in ("dipole", "dq") else ""
        for ax in (axx, axy):
            ax.set_xlabel("s, м"); ax.grid(alpha=0.3)
        axx.set_ylabel("x (гориз.), мм"); axy.set_ylabel("y (верт.), мм")
        axx.set_title(f"x-плоскость — ток {cur:+.2f}А (— диффур, -- кик)")
        axy.set_title(f"y-плоскость{ttl2}")
        fig.suptitle(f"{dev} ({kind})  L_eff={Lc:.4f}м  RMS={r['rms_mm']:.4f}мм", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        tok = f"{cur:+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")   # только токен тока (не .png!)
        png = os.path.join(outdir, f"traj_{dev}_{tok}.png")
        fig.savefig(png, dpi=110); plt.close(fig); figs.append((cur, png))
    return figs


def plot_xy(dev, kind, recs, probes, outdir):
    """Проекция пучка в плоскости xy ДО элемента (стартовый круг Ø2мм) и ЗА элементом (диффур ○ vs кик ✕)
    при максимальном по модулю токе. Показывает, как элемент преобразует сечение пучка."""
    r = max(recs, key=lambda x: abs(x["cur"])); Lc = recs[0].get("Leff_common", r["Leff"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(9.5, 4.6))
    a0.scatter([p[0] * 1e3 for p in probes], [p[1] * 1e3 for p in probes], c="#2b7bba", s=32)
    a0.set_title("до элемента (старт)")
    mats = _mat_for_plot(dev, r["coeffs"], Lc, kind, r["geom"], probes)
    for d, m in zip(r["diff_su"], mats):
        a1.scatter(d[-1, 1] * 1e3, d[-1, 2] * 1e3, c="#2b7bba", s=34, marker="o")
        if m is not None:
            a1.scatter(m[-1, 1] * 1e3, m[-1, 2] * 1e3, c="#d1495b", s=44, marker="x")
    a1.set_title(f"за элементом, {r['cur']:+.2f} А  (o диффур, x кик)")
    for a in (a0, a1):
        a.set_xlabel("x (гориз.), мм"); a.set_ylabel("y (верт.), мм"); a.grid(alpha=0.3); a.set_aspect("equal", "datalim")
    fig.suptitle(f"{dev} ({kind}) — сечение пучка в плоскости xy", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    png = os.path.join(outdir, f"xy_{dev}.png"); fig.savefig(png, dpi=110); plt.close(fig)
    return png


def process_element(dev, kind, probes, outdir, bend_yaw=0.0):
    """Полный опыт по одному элементу: свип+фит → датасет + графики. Возврат dict для отчёта."""
    print(f"\n=== элемент {dev} ({kind}) ===")
    recs = sweep_currents(dev, kind, probes)
    if not recs:
        return None
    ds = save_dataset(dev, kind, recs, outdir)
    figs = plot_element(dev, kind, recs, probes, outdir)
    xy = plot_xy(dev, kind, recs, probes, outdir)
    return dict(dev=dev, kind=kind, recs=recs, figs=figs, xy=xy, dataset=ds,
                Leff_common=recs[0].get("Leff_common"))


def build_report(elements, path):
    """Word-отчёт: глава на элемент (таблица метрик current↔физкоэф+L_eff+RMS + графики x/y)."""
    from docx import Document
    from docx.shared import Inches, Pt
    doc = Document()
    doc.add_heading("SynchrotronSim — сопоставление киковой и диффурной моделей по элементам", level=0)
    doc.add_paragraph("Диффур (C++ SynchrotronTracker, Boris по картам поля) — эталон. Киковая модель (Python "
                      "matrix_tracker в режиме SS_SECTOR — симплектический drift-kick-drift, порт elegant "
                      "KQUAD/KSEXT/KOCT/CSBEND, а НЕ передаточные матрицы) подогнана последовательными "
                      "приближениями: физкоэф на каждый ток + один L_eff на элемент. Пучок 13 частиц (12 по "
                      "кругу Ø2мм + центр). Прямые токовые — свип −50..+50 А; изгибающие DQ — вокруг рабочего "
                      "тока; PM-диполи — одна конфигурация (поле фиксировано). Изгибающие сравниваются в "
                      "проекции Френе (отклонение от опорной дуги).")
    for el in elements:
        if el is None: continue
        keys = KIND_COEFF[el["kind"]]
        doc.add_heading(f"{el['dev']}  ({el['kind']})", level=1)
        doc.add_paragraph(f"Общий L_eff = {el['Leff_common']:.4f} м (одна геом. константа на все токи; "
                          f"подобрана по суммарной невязке всех токов). Физкоэф — на каждый ток. "
                          f"Макс. невязка RMS по токам = {max(r['rms_mm'] for r in el['recs']):.4f} мм.")
        t = doc.add_table(rows=1, cols=len(keys) + 2); t.style = "Light Grid Accent 1"
        hdr = t.rows[0].cells; hdr[0].text = "ток, А"
        for j, k in enumerate(keys): hdr[j + 1].text = k
        hdr[len(keys) + 1].text = "RMS, мм"
        pm = is_pm(el["dev"])
        for r in el["recs"]:
            c = t.add_row().cells; c[0].text = "PM" if pm else f"{r['cur']:+.2f}"
            for j, k in enumerate(keys): c[j + 1].text = f"{r['coeffs'][k]:+.4e}"
            c[len(keys) + 1].text = f"{r['rms_mm']:.4f}"
        if el.get("xy") and os.path.exists(el["xy"]):
            doc.add_paragraph("Сечение пучка в плоскости xy (до / за элементом):").runs[0].font.size = Pt(10)
            doc.add_picture(el["xy"], width=Inches(5.6))
        for cur, png in el["figs"]:
            if not os.path.exists(png): continue          # защита: не ронять отчёт из-за одной картинки
            cap = "Траектории (PM):" if pm else f"Траектории при токе {cur:+.2f} А:"
            doc.add_paragraph(cap).runs[0].font.size = Pt(10)
            doc.add_picture(png, width=Inches(6.3))
    doc.save(path)
    return path


# ============================ ПЕРЕЧИСЛЕНИЕ + ПОЛНЫЙ ЦИКЛ ====================================
def candidates(kinds=None):
    """Устройства с картой поля И mat_* И ТОКОВЫЕ (nominal_arg — число) → [(dev, kind)]. Фильтр по kinds.
    nominal_arg='-' ⇒ постоянный/PM (Сибирь s2_, диполи DP) — свип тока бессмыслен, исключаем."""
    def current_driven(d):
        na = (HW[d].get("nominal_arg") or "").strip()
        try: return float(na) > 0
        except ValueError: return False
    maps = {d for d, r in HW.items() if (r.get("map_bin_path") or "").strip() and current_driven(d)}
    kd = {}
    for r in mt._read("lattice_configs.csv"):
        d = r["device_type_id"]
        if d in maps and d not in kd and any((r.get(k) or "").strip() not in ("", "0")
                                             for k in ("mat_B", "mat_G", "mat_S", "mat_O")):
            kd[d] = r.get("kind") or ""
    out = [(d, k) for d, k in sorted(kd.items()) if k]
    return [(d, k) for d, k in out if (kinds is None or k in kinds)]


def candidates_bending(kinds=None):
    """Изгибающие для сравнения: PM-диполи (nominal_arg='-') + ГЛАВНЫЕ обмотки DQ (I1 — несут mat_B/mat_G).
    I2/I3 (вспомог., mat_*=0) и комбо — отдельным этапом (нет проектных mat_* для сравнения)."""
    out = []
    for d, r in sorted(LCROW.items()):
        k = r.get("kind") or ""
        has_map = (HW.get(d, {}).get("map_bin_path") or "").strip()
        if not has_map: continue
        if k == "dipole" and is_pm(d):
            out.append((d, "dipole"))
        elif k == "dq" and "_I1_" in d:
            out.append((d, "dq"))
    return [(d, k) for d, k in out if (kinds is None or k in kinds)]


PROBES = [(0.0, 0.0)] + [(1e-3 * math.cos(2 * math.pi * k / 12), 1e-3 * math.sin(2 * math.pi * k / 12)) for k in range(12)]


def run_all(kinds=("quad", "sext", "oct")):
    """Полный цикл: обработать все кандидаты (по kinds) → общий Word-отчёт + датасеты. Прямые токовые —
    candidates; изгибающие (dipole/dq) — candidates_bending (PM-диполи + главные обмотки DQ)."""
    os.makedirs(OUTDIR, exist_ok=True)
    straight = [k for k in kinds if k in ("quad", "sext", "oct")]
    bend = [k for k in kinds if k in ("dipole", "dq")]
    devs = candidates(straight) + (candidates_bending(bend) if bend else [])
    print(f"кандидатов ({','.join(kinds)}): {len(devs)}")
    bak = os.path.join(TB, "_cmp_backup"); os.makedirs(bak, exist_ok=True)
    for f in CSV4: shutil.copy2(os.path.join(TB, f), os.path.join(bak, f))
    elements = []
    try:
        for dev, kind in devs:
            try:
                el = process_element(dev, kind, PROBES, OUTDIR)
                if el: elements.append(el)
            except Exception as e:
                print(f"  ! {dev}: {type(e).__name__}: {e}"); cleanup()
    finally:
        cleanup()
        for f in CSV4:                                # восстановить при сбое очистки
            b = os.path.join(bak, f)
            if os.path.exists(b) and os.path.getsize(os.path.join(TB, f)) == 0:
                shutil.copy2(b, os.path.join(TB, f))
    rep = build_report(elements, os.path.join(OUTDIR, "element_compare_report.docx"))
    print(f"\nГОТОВО: {len(elements)} элементов → {rep}")
    return elements


# ============================ КОМБО-ДАТАСЕТ DQ (все комбинации токов обмоток) =================
def run_diffur_windings(windings, yaw=0.0):
    """Диффур для N со-локованных обмоток (суперпозиция полей). Возврат (trajs, geom). Чистит вызывающий."""
    geom = build_windings(windings, yaw)
    r = subprocess.run([EXE, str(TMP_EXP)], capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        raise RuntimeError(f"диффур упал: {r.stderr[-400:]}")
    return read_bin(TMP_EXP), geom


def tune_bend_yaw_windings(windings, iters=6, tol=1e-6):
    """yaw для симметричного изгиба ОПОРНОЙ (ч0) комбинированного поля обмоток (аналог tune_bend_yaw)."""
    yaw = 0.0
    for _ in range(iters):
        trajs, _ = run_diffur_windings(windings, yaw)
        ref = trajs[0]; s = ref[:, 0]; m = s > (s[-1] - 0.3)
        phi = math.atan(np.polyfit(ref[m, 0], ref[m, 1], 1)[0]); cleanup()
        new = phi / 2.0
        if abs(new - yaw) < tol:
            yaw = new; break
        yaw = new
    return yaw


def _fit_gs_combo(main, geom, diff_su, L, theta0):
    """Комбо-фит DQ при ФИКС. L: B из угла; (G,S) — старт G из 1D-скана (_coeff_at_L, S=0) → Nelder-Mead по
    (G,S) с масштабным симплексом (шаг G~0.01, S~1). Секступоль S ловит feed-down комбинированной функции
    (память dq-gk-diagnosis). Возврат (coeffs{B,G,S}, mse)."""
    c0, _ = _coeff_at_L(main, "dq", geom, PROBES, diff_su, L, theta0)    # G0 при S=0
    B = c0["B"]; x0 = np.array([c0["G"], 0.0])
    simplex = np.array([x0, x0 + [0.01, 0.0], x0 + [0.0, 1.0]])          # масштабы G и S

    def obj(x):
        try:
            return traj_mse(diff_su, track_set(main, {"B": B, "G": x[0], "S": x[1]}, L, "dq", geom, PROBES))
        except mt.ParticleLost:
            return 1e3

    res = minimize(obj, x0, method="Nelder-Mead",
                   options=dict(initial_simplex=simplex, xatol=1e-6, fatol=1e-16, maxiter=4000, maxfev=4000))
    return {"B": B, "G": float(res.x[0]), "S": float(res.x[1])}, float(res.fun)


def save_combo_dataset(base, Leff, recs, outdir):
    """CSV комбо: I1;I2;I3;B;G[;S];RMS (запятая-десятичная для Excel, ; поле). Столбец S — только если фит
    его давал (sext=True). L_eff общий в шапке."""
    os.makedirs(outdir, exist_ok=True)
    has_s = bool(recs) and "S" in recs[0]["coeffs"]
    keys = ["B", "G", "S"] if has_s else ["B", "G"]
    path = os.path.join(outdir, f"combo_dataset_{base}.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["L_eff_общий_м", str(Leff).replace(".", ",")]); w.writerow([])
        w.writerow(["I1_А", "I2_А", "I3_А"] + keys + ["RMS_мм"])
        for r in recs:
            row = [r["i1"], r["i2"], r["i3"]] + [r["coeffs"][k] for k in keys] + [r["rms_mm"]]
            w.writerow([str(x).replace(".", ",") for x in row])
    print(f"  датасет: {path} ({len(recs)} комбо, {'B,G,S' if has_s else 'B,G'})")
    return path


def run_dq_combo(base="L400", ngrid=4, sext=False):
    """Комбо-датасет DQ: перебор ВСЕХ комбинаций токов обмоток I1/I2/I3. Шаг подобран под лимит ≤120 примеров:
    ngrid=4 ⇒ 4³=64 (5³=125 превысило бы). Диапазон [−Iмакс,+Iмакс] на обмотку, Iмакс = рабочий ток I1
    (dq_operating_current). L_eff — общий (из свипа I1). На комбо: yaw(симметрия) → диффур(суперпозиция
    полей обмоток) → фит B(из угла)+G [+S если sext] при фикс. L_eff. sext=False (по умолч.) — DQ как
    чистый диполь+квадруполь (как elegant/.lte: CSBEND с K1, K2=0; секступоль — физ.феномен карт, не elegant,
    см. dq_gk_diagnosis). Датасет (I1,I2,I3)→(B,G[,S],RMS)."""
    wdev = [f"DQ_I{k}_1A_{base}" for k in (1, 2, 3)]
    main = wdev[0]; N3 = ngrid ** 3
    Imax = dq_operating_current(main)
    print(f"DQ {base}: рабочий ток I1≈{Imax:.2f}А; сетка ±{Imax:.2f}А × {ngrid}/обмотка ⇒ {N3} комбо (лимит ≤120)")
    Leff = sweep_currents(main, "dq", PROBES)[0]["Leff_common"]   # общий L_eff из свипа главной обмотки
    print(f"  общий L_eff (из свипа I1) = {Leff:.4f}м; перебор комбинаций…")
    grid = np.linspace(-Imax, Imax, ngrid)
    bak = os.path.join(TB, "_cmp_backup"); os.makedirs(bak, exist_ok=True)
    for f in CSV4: shutil.copy2(os.path.join(TB, f), os.path.join(bak, f))
    recs = []
    try:
        for n, (i1, i2, i3) in enumerate(itertools.product(grid, grid, grid), 1):
            wnd = [(wdev[0], round(float(i1), 4)), (wdev[1], round(float(i2), 4)), (wdev[2], round(float(i3), 4))]
            try:
                yaw = tune_bend_yaw_windings(wnd)
                trajs, geom = run_diffur_windings(wnd, yaw)
                diff_su = diff_local(trajs, geom["xstart"], bending=True)
                if sext:
                    coeffs, mse = _fit_gs_combo(main, geom, diff_su, Leff, 2.0 * yaw)          # B+G+S
                else:
                    coeffs, mse = _coeff_at_L(main, "dq", geom, PROBES, diff_su, Leff, theta0=2.0 * yaw)  # B+G
                recs.append(dict(i1=float(i1), i2=float(i2), i3=float(i3), coeffs=coeffs, rms_mm=math.sqrt(mse) * 1e3))
            except Exception as e:
                print(f"  [{n}/{N3}] I=({i1:.2f},{i2:.2f},{i3:.2f}): {type(e).__name__}: {e}")
            cleanup()
            if n % 8 == 0:
                print(f"  {n}/{N3}  (посл. RMS={recs[-1]['rms_mm']:.4f}мм)" if recs else f"  {n}/{N3}")
    finally:
        cleanup()
        for f in CSV4:
            b = os.path.join(bak, f)
            if os.path.exists(b) and os.path.getsize(os.path.join(TB, f)) == 0:
                shutil.copy2(b, os.path.join(TB, f))
    rmss = [r["rms_mm"] for r in recs]
    print(f"\nГОТОВО комбо {base}: {len(recs)}/{N3}, RMS {min(rmss):.4f}..{max(rmss):.4f}мм")
    save_combo_dataset(base, Leff, recs, OUTDIR)
    return recs


# ============================ ПОСТОЯННЫЙ КОНФИГ ДЛЯ VISUALIZE.PY (диффур + кик в 3D) =========
VIZ_START = 4   # id постоянных конфигов визуализации С 4 (1/2/3 — реальные опыты, TMP=900 чистится)


def kick_to_lab(kick_trajs, ref_lab):
    """Киковые ЛОКАЛЬНЫЕ (s,u,v) → лаб XYZ наложением на опорную орбиту ref_lab (диффур ч0, лаб XYZ) в рамке
    Френе: точка = ч0(s) + u·n(s) + v·ẑ, n — ⊥ в плоскости изгиба X-Y. Обратно к diff_local(bending). Для
    прямых ч0 идёт вдоль +x ⇒ n=(0,1) ⇒ (X=xstart+s, Y=u, Z=v)."""
    Xr, Yr, Zr = ref_lab[:, 0], ref_lab[:, 1], ref_lab[:, 2]
    Sref = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(Xr), np.diff(Yr)))])
    tx = np.gradient(Xr); ty = np.gradient(Yr); tn = np.hypot(tx, ty); tn[tn == 0] = 1.0
    tx /= tn; ty /= tn; nx, ny = -ty, tx
    out = []
    for kt in kick_trajs:
        s, u, v = kt[:, 0], kt[:, 1], kt[:, 2]
        X = np.interp(s, Sref, Xr) + u * np.interp(s, Sref, nx)
        Y = np.interp(s, Sref, Yr) + u * np.interp(s, Sref, ny)
        Z = np.interp(s, Sref, Zr) + v
        out.append(np.column_stack([X, Y, Z]))
    return out


def write_bin(path, trajs):
    """Траектории (список (N,3) лаб XYZ) → .bin формата C++ (uint32 count + N×3 float64 на частицу)."""
    with open(path, "wb") as f:
        for a in trajs:
            a = np.ascontiguousarray(a, dtype=np.float64)
            f.write(struct.pack("I", a.shape[0])); f.write(a.tobytes())


def design_current(dev, kind, recs):
    """Ток, при котором фитируемый ДОМИНАНТНЫЙ физкоэф = ПРОЕКТНОМУ mat_* (.lte СИЛЫ, лежит в БД). Линейная
    экстраполяция по свипу (coeff ≈ slope·I). DQ — по изгибу mat_B (рабочая точка). PM → '-' (без тока).
    Знак тока даёт полярность (у части квадов проектный градиент противоположного знака)."""
    if is_pm(dev):
        return "-"
    dom = {"quad": "G", "sext": "S", "oct": "O", "dq": "B"}[kind]
    design = float(LCROW[dev].get({"G": "mat_G", "S": "mat_S", "O": "mat_O", "B": "mat_B"}[dom]) or 0.0)
    Is = np.array([r["cur"] for r in recs], float); Cs = np.array([r["coeffs"][dom] for r in recs], float)
    slope = np.polyfit(Is, Cs, 1)[0]                    # coeff ≈ slope·I
    return round(design / slope, 4) if abs(slope) > 1e-30 else 0.0


def _finalize_demo_config(lat_id, dev, kind, coeffs, Leff):
    """Записать в демо-конфиг ФИТИРОВАННЫЕ физкоэф (mat_*) + геометрию (eff_length=Leff, kind, edges из
    lattice-1). build_windings ОБНУЛЯЕТ mat_* (диффур берёт карту), но КИКОВАЯ модель (viz_html) читает
    mat_* — без этого показывала бы НУЛЕВОЙ кик (прямые линии). Кик = фит ⇒ совпадает с диффуром."""
    src = LCROW.get(dev, {})
    matmap = {"B": "mat_B", "G": "mat_G", "S": "mat_S", "O": "mat_O", "Gv": "mat_Gv"}
    h, rows = _rows("lattice_configs.csv"); ci = {c: i for i, c in enumerate(h)}
    for r in rows:
        if not (r and r[ci["lattice_id"]] == str(lat_id) and r[ci["device_type_id"]] == dev):
            continue
        for k, mc in matmap.items():
            if mc in ci: r[ci[mc]] = repr(float(coeffs.get(k, 0.0)))
        if "eff_length" in ci: r[ci["eff_length"]] = repr(float(Leff))
        if "kind" in ci: r[ci["kind"]] = kind
        for c in ("n_slices", "edge1", "edge2"):                     # геометрия изгибов (края) из lattice-1
            if c in ci and src.get(c) not in (None, ""): r[ci[c]] = src[c]
    _write("lattice_configs.csv", h, rows)


def setup_viz(dev, cur, viz_id, Leff=None):
    """ПОСТОЯННЫЙ конфиг id=viz_id для visualize.py: диффур (trajectories_exp<id>.bin) + кик
    (trajectories_exp<id>_kick.bin) в лаб-XYZ при ТОЧНОМ токе cur ('-' ⇒ PM). Leff — общий L_eff (если None,
    берётся из свипа). Конфиг НЕ чистится. mat_* = фит (для viz_html кик-модели). Возврат dict со сводкой."""
    kind = LCROW[dev].get("kind") or "quad"
    bending = kind in ("dipole", "dq"); pm = (cur == "-")
    if Leff is None:
        rc = sweep_currents(dev, kind, PROBES); Leff = rc[0]["Leff_common"]
    cur_val = 0.0 if pm else float(cur)
    yaw = tune_bend_yaw(dev, cur_val) if bending else 0.0
    geom = build_windings([(dev, "-" if pm else cur_val)], yaw, lat_id=viz_id, grp_id=viz_id, exp_id=viz_id)
    rr = subprocess.run([EXE, str(viz_id)], capture_output=True, text=True, cwd=ROOT)   # диффур → .bin
    if rr.returncode != 0:
        raise RuntimeError(f"диффур упал: {rr.stderr[-400:]}")
    trajs = read_bin(viz_id); ref_lab = trajs[0]
    diff_su = diff_local(trajs, geom["xstart"], bending)             # фит коэф при ТОЧНОМ токе (общий L_eff)
    if kind == "dq":
        coeffs, mse = _fit_gs_combo(dev, geom, diff_su, Leff, 2.0 * yaw)
    else:
        coeffs, mse = _coeff_at_L(dev, kind, geom, PROBES, diff_su, Leff, theta0=2.0 * yaw)
    kick_lab = kick_to_lab(track_set(dev, coeffs, Leff, kind, geom, PROBES), ref_lab)   # кик → лаб XYZ по ч0
    write_bin(os.path.join(RESULTS, f"trajectories_exp{viz_id}_kick.bin"), kick_lab)
    _finalize_demo_config(viz_id, dev, kind, coeffs, Leff)          # mat_* + геометрия в конфиг → viz_html кик ≠ 0
    rms = math.sqrt(mse) * 1e3; cs = " ".join(f"{k}={v:+.3e}" for k, v in coeffs.items())
    print(f"  id={viz_id} {dev} ({kind}) ток={'PM' if pm else f'{cur_val:+.3f}А'}  L_eff={Leff:.4f}  {cs}  RMS={rms:.4f}мм")
    return dict(id=viz_id, dev=dev, kind=kind, cur=("PM" if pm else cur_val), Leff=Leff, coeffs=coeffs, rms=rms)


def _drop_config(*ids):
    """Убрать конфиги с данными id из 4 CSV (при переезде id/пересборке)."""
    for key in ids:
        for name in CSV4:
            h, rows = _rows(name); _write(name, h, [r for r in rows if not (r and r[0] == str(key))])
        for suf in ("", "_kick"):
            p = os.path.join(RESULTS, f"trajectories_exp{key}{suf}.bin")
            if os.path.exists(p): os.remove(p)


def setup_viz_all(start_id=VIZ_START):
    """Постоянные конфиги визуализации для ВСЕХ элементов, id подряд с start_id. Ток каждого = ПРОЕКТНЫЙ
    (design_current: фит.коэф = mat_* из .lte СИЛЫ). Диполи — PM. Печатает таблицу id↔элемент↔ток."""
    cands = candidates(("quad", "sext", "oct")) + candidates_bending()
    print(f"=== визуализация ВСЕХ элементов ({len(cands)}), id с {start_id}, токи по ПРОЕКТНЫМ коэф ===")
    _drop_config(950, *range(start_id, start_id + len(cands) + 5))   # снять старый 950 и прежние viz-конфиги
    summary = []
    try:
        for i, (dev, kind) in enumerate(cands):
            vid = start_id + i
            try:
                if is_pm(dev):
                    summary.append(setup_viz(dev, "-", vid, Leff=None))
                else:
                    recs = sweep_currents(dev, kind, PROBES)
                    summary.append(setup_viz(dev, design_current(dev, kind, recs), vid, Leff=recs[0]["Leff_common"]))
            except Exception as e:
                print(f"  ! id={vid} {dev}: {type(e).__name__}: {e}"); cleanup()
    finally:
        cleanup()
    print("\n=== СВОДКА (визуализировать: py -3.14 visualize.py <id>) ===")
    for s in summary:
        c = s["cur"] if s["cur"] == "PM" else f"{s['cur']:+.3f}А"
        print(f"  id {s['id']:2d}  {s['dev']:18s} ({s['kind']:6s})  ток {c:>11s}  RMS {s['rms']:.4f}мм")
    return summary


# ============================ МАТРИЦЫ СВЯЗИ ТОК→ФИЗКОЭФ + РАЗДЕЛ ОТЧЁТА =======================
# Мультипольные моменты, извлекаемые из КИКА (изменения угла): гориз/верт диполь, норм/скью квад, секступ, окт.
MULTIPOLES = ["B", "Bv", "G", "Gv", "S", "O"]


def extract_coeffs(windings, L):
    """Диффур комбинации СО-ЛОКОВАННЫХ обмоток (прямой элемент, yaw=0) → мультипольные mat_* из КИКА:
    выходной наклон (sy=dY/dX гориз, sz=dZ/dX верт) 13 частиц раскладывается по (Y,Z) до октуполя.
    Даёт ВСЕ моменты, вкл. ВЕРТИКАЛЬНЫЙ диполь (Bv), которого нет в столбцах БД. Возврат dict MULTIPOLES→коэф."""
    trajs, _ = run_diffur_windings(windings, 0.0)
    Y = np.array([p[0] for p in PROBES]); Z = np.array([p[1] for p in PROBES])
    sy, sz = [], []
    for a in trajs:
        s = a[:, 0]; m = s > (s[-1] - 0.3)                           # выходной дрейф — прямая
        sy.append(np.polyfit(a[m, 0], a[m, 1], 1)[0]); sz.append(np.polyfit(a[m, 0], a[m, 2], 1)[0])
    cleanup()
    A = np.column_stack([np.ones(len(Y)), Y, Z, 0.5 * (Y * Y - Z * Z), Y * Z,
                         (Y ** 3 - 3 * Y * Z * Z) / 6.0, (3 * Y * Y * Z - Z ** 3) / 6.0])
    ay = np.linalg.lstsq(A, np.array(sy), rcond=None)[0]
    az = np.linalg.lstsq(A, np.array(sz), rcond=None)[0]
    QPL = mt.QP * L                                                  # физкоэф = -кик/(QP·L) (диполь без минуса)
    return {"B": ay[0] / QPL, "Bv": az[0] / QPL, "G": -ay[1] / QPL, "Gv": -ay[2] / QPL,
            "S": -ay[3] / QPL, "O": -ay[5] / QPL}


DOM = {"quad": "G", "sext": "S", "oct": "O"}              # доминантный (проектный) момент чистого элемента
LIN_FRACS = 3                                            # первые N точек лестницы (около рабочей точки) — линейны, по ним НАКЛОН


def _ladder_fracs(kind):
    """Доли рабочего тока для свипа. Чистый элемент — вверх за рабочую точку (ищем начало нелинейности);
    корректор (kicker) — доли номинала (не гонять за паспорт)."""
    if kind in DOM:
        return np.array([0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0])
    return np.array([0.25, 0.5, 0.75, 1.0, 1.25])


def element_operating_current(kind, devs, L):
    """Ток масштаба РАБОЧЕЙ точки (где designed-момент = проектному mat_* из .lte). Корректор → номинал из
    имени (…_100A_…). Нужен лишь для МАСШТАБА лестницы (наклон инвариантен к масштабу, связь линейна через 0)."""
    import re
    if kind in DOM:
        mp = DOM[kind]; col = {"G": "mat_G", "S": "mat_S", "O": "mat_O"}[mp]
        design, dwd = 0.0, devs[0]
        for d in devs:                                       # доминантная обмотка = с ненулевым проектным mat_*
            v = float(LCROW[d].get(col) or 0.0)
            if abs(v) > abs(design): design, dwd = v, d
        sl = extract_coeffs([(dwd, 1.0)], L)[mp]             # наклон коэф/ток около нуля (1 А заведомо линейно)
        if abs(design) > 1e-30 and abs(sl) > 1e-30:
            return abs(design / sl)
    m = re.search(r"_(\d+(?:\.\d+)?)A", devs[0])             # номинал из имени устройства
    return float(m.group(1)) if m else 45.0


def save_sweep_dataset(base, labels, curs, Vals):
    """Свип-датасет матрицы: для КАЖДОЙ обмотки — лестница токов → извлечённые мультиполи (по которым построена
    M). Запятая-десятичная (Excel), ';'. → datasets/sweep_<base>.csv."""
    os.makedirs(DATASETS, exist_ok=True)
    path = os.path.join(DATASETS, f"sweep_{base}.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["обмотка", "ток_А"] + MULTIPOLES)
        for j, lbl in enumerate(labels):
            for t, cur in enumerate(curs):
                row = [lbl, cur] + [Vals[i, j, t] for i in range(len(MULTIPOLES))]
                w.writerow([row[0]] + [str(round(float(x), 8)).replace(".", ",") for x in row[1:]])
    return path


def winding_matrix(wdevs, L, i_op, fracs):
    """Матрица связи [мультиполи × обмотки] + сырьё для оценки релевантного диапазона. На КАЖДУЮ обмотку —
    свип по лестнице i_op·fracs; НАКЛОН mat_*/ток берётся в ЛИНЕЙНОЙ области (первые LIN_FRACS точек, около
    рабочей точки), т.к. связь линейна через 0 (суперпозиция полей). Возврат (M, curs, Vals[момент,обмотка,ток])."""
    curs = i_op * np.asarray(fracs, float)
    nlin = min(LIN_FRACS, len(curs))
    M = np.zeros((len(MULTIPOLES), len(wdevs)))
    Vals = np.zeros((len(MULTIPOLES), len(wdevs), len(curs)))
    for j, wd in enumerate(wdevs):
        for t, cur in enumerate(curs):
            c = extract_coeffs([(wd, cur)], L)
            for i, mp in enumerate(MULTIPOLES): Vals[i, j, t] = c[mp]
        for i in range(len(MULTIPOLES)):
            M[i, j] = np.polyfit(curs[:nlin], Vals[i, j, :nlin], 1)[0]
    return M, curs, Vals


def _prefix_r2(curs, y, t):
    """R² линейного фита coef=slope·ток по префиксу из t точек (t≥2; 2 точки → 1.0)."""
    cx = curs[:t]; yy = y[:t]; sl, bb = np.polyfit(cx, yy, 1)
    sst = ((yy - yy.mean()) ** 2).sum()
    return (1.0 - ((yy - (sl * cx + bb)) ** 2).sum() / sst) if sst > 1e-30 else 1.0


def _relevant_imax(curs, vals_rows, thr=0.95):
    """Макс. ток лестницы, при котором линейный фит по префиксу [0..t] держит R²≥thr для ВСЕХ значимых
    (row,обмотка). Возврат (Imax, R²_на_Imax). vals_rows — список массивов коэф(ток)."""
    n = len(curs); best = 2                                   # 2 точки → линия точна (R²=1)
    for t in range(3, n + 1):
        if all(_prefix_r2(curs, y, t) >= thr for y in vals_rows): best = t
        else: break
    r2 = min((_prefix_r2(curs, y, best) for y in vals_rows), default=1.0)
    return float(curs[best - 1]), float(r2)


def _combo_grid(nw, lvl):
    """Узлы комбо-сетки токов обмоток в кубе [−lvl,lvl]^nw (без нулевого узла). Малая размерность (3^nw≤30) —
    полная 3-уровневая сетка; средняя (2^nw≤40) — углы куба; большая — случайная выборка (фикс. seed)."""
    import itertools
    if 3 ** nw <= 30:
        pts = list(itertools.product([-lvl, 0.0, lvl], repeat=nw))
    elif 2 ** nw <= 40:
        pts = list(itertools.product([-lvl, lvl], repeat=nw))
    else:
        rng = np.random.default_rng(0); pts = [tuple(rng.uniform(-lvl, lvl, nw)) for _ in range(32)]
    return [p for p in pts if any(abs(x) > 1e-9 for x in p)]


def save_combo_verify(base, devs, L, mps, recs):
    """CSV комбо-проверки суперпозиции: токи обмоток; извлечённые моменты; предсказание матрицы; |невязка|.
    Запятая-десятичная (Excel, память csv-comma-decimal-for-excel), ';' поле."""
    import re
    lbls = [(re.search(r"_(I\d+|c\d+|main)", d).group(1) if re.search(r"_(I\d+|c\d+|main)", d) else "I") for d in devs]
    os.makedirs(DATASETS, exist_ok=True)
    path = os.path.join(DATASETS, f"combo_{base}.csv")
    with open(path, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["L_eff_м", str(round(L, 6)).replace(".", ",")]); w.writerow([])
        w.writerow([f"I_{l}_А" for l in lbls] + [f"{m}_извл" for m in mps]
                   + [f"{m}_матрица" for m in mps] + [f"|Δ{m}|" for m in mps])
        for I, act, pred in recs:
            row = list(I) + list(act) + list(pred) + list(np.abs(np.asarray(act) - np.asarray(pred)))
            w.writerow([str(round(float(x), 8)).replace(".", ",") for x in row])
    return path


def combo_verify(base, devs, L, M, rows, imax, save=True):
    """Проверка СУПЕРПОЗИЦИИ многообмоточного элемента: комбо-сетка токов обмоток в ±0.5·Imax → диффур
    (со-локованные обмотки суперпозируют поля, getTotalFields) → извлечённые designed-моменты vs предсказание
    ЛИНЕЙНОЙ матрицы M·I. Малая невязка = матрица корректно описывает комбинированное возбуждение. Датасет →
    combo_verify_<base>.csv. Возврат dict(n, max_abs, rel, moments, path)."""
    grid = _combo_grid(len(devs), 0.5 * imax)
    Mk = M[rows]; mps = [MULTIPOLES[i] for i in rows]
    recs = []; maxa = np.zeros(len(rows)); num = np.zeros(len(rows)); den = np.zeros(len(rows))
    for I in grid:
        c = extract_coeffs([(devs[j], float(I[j])) for j in range(len(devs))], L)
        actual = np.array([c[MULTIPOLES[i]] for i in rows]); pred = Mk @ np.asarray(I, float)
        d = np.abs(actual - pred); maxa = np.maximum(maxa, d)
        num += d ** 2; den += actual ** 2
        recs.append((I, actual, pred))
    per = {mps[k]: float(math.sqrt(num[k] / max(den[k], 1e-30))) for k in range(len(rows))}  # отн. невязка по каждому моменту
    path = save_combo_verify(base, devs, L, mps, recs) if save else None
    return dict(n=len(grid), per=per, max_abs=float(maxa.max()), moments=mps, path=path)


def element_windings(kinds=("quad", "sext", "oct", "kicker"), include_pm=False):
    """Токовые элементы реш.1 с картами → (base, kind, [обмотки], L, [метки обмоток]). Группировка по базовому
    имени (суффикс обмотки _I<k>/_c<k>/_main). `kinds` — какие виды включать; `include_pm=True` добавляет
    PM-диполи (фикс. поле, без свипа — ветка fixed). DQ всегда отдельно (combo, compute_coupling)."""
    import re
    grp = {}
    for d, r in LCROW.items():
        if not (HW.get(d, {}).get("map_bin_path") or "").strip(): continue
        pm = (HW.get(d, {}).get("nominal_arg") or "").strip() == "-"
        if pm and not include_pm: continue
        base = re.sub(r"_I\d+|_c\d+|_main", "", d)
        grp.setdefault((base, r.get("kind") or ""), []).append(d)
    out = []
    for (base, kind), devs in sorted(grp.items()):
        if kind not in kinds: continue                                           # DQ отдельно (combo)
        devs = sorted(devs)
        labels = [(re.search(r"_(I\d+|c\d+|main)", d).group(1) if re.search(r"_(I\d+|c\d+|main)", d) else "I") for d in devs]
        out.append((base, kind, devs, float(LCROW[devs[0]].get("eff_length") or 0.1), labels))
    return out


def library_windings():
    """Полный охват матриц: все элементы реш.1 с картами И проектными lte-данными — прямые (quad/sext/oct/
    kicker) + PM-диполи (dipole, fixed). DQ добавляется отдельно (combo) на уровне таблицы/отчёта."""
    return element_windings(kinds=("quad", "sext", "oct", "kicker", "dipole"), include_pm=True)


# Чистый элемент несёт ТОЛЬКО свой мультипольный момент (как elegant: KQUAD=K1, KSEXT=K2, KOCT=K3 — без
# посторонних членов). Корректор (kicker) — многофункциональный, показываем все значимые моменты.
# dipole (PM) — фикс. поле: диполь B + слабый краевой градиент G (свипа тока нет, ветка fixed).
DESIGNED_MOMENT = {"quad": ["G"], "sext": ["S"], "oct": ["O"], "dipole": ["B", "G"]}


def compute_element_matrices(elements=None):
    """Матрицы связи [мультиполи × обмотки] токовых элементов (кроме DQ — combo). Для чистых элементов —
    ТОЛЬКО designed-момент (как elegant); корректор — все значимые по кику; PM-диполи — ФИКС. состав (свипа
    нет, ветка pm). + R² линейности, комбо-проверка. `elements` — список из element_windings/library_windings
    (None → element_windings()). Возврат список dict (pm=True → есть fixed{момент:коэф}, нет M/imax)."""
    res = []
    r = 1e-3                                                          # опорный радиус пучка (1мм) для значимости
    scale = np.array([1.0, 1.0, r, r, 0.5 * r * r, r ** 3 / 6.0])    # вклад в КИК на единицу mat_* при 1мм (n=0..3)
    for base, kind, devs, L, labels in (element_windings() if elements is None else elements):
        if is_pm(devs[0]):                                           # PM: фикс. поле, без свипа → ФИКС. состав
            c = extract_coeffs([(devs[0], "-")], L)
            rows = [MULTIPOLES.index(mp) for mp in DESIGNED_MOMENT.get(kind, ["B", "G"])]
            res.append(dict(base=base, kind=kind, labels=labels, L=L, pm=True, rows=rows,
                            fixed={MULTIPOLES[i]: float(c[MULTIPOLES[i]]) for i in rows}))
            continue
        i_op = element_operating_current(kind, devs, L)              # масштаб лестницы = рабочая точка / номинал
        M, curs, Vals = winding_matrix(devs, L, i_op, _ladder_fracs(kind))
        save_sweep_dataset(base, labels, curs, Vals)                 # датасет свипа (токи обмоток → мультиполи)
        if kind in DESIGNED_MOMENT:                                  # чистый элемент — ТОЛЬКО свой момент
            rows = [MULTIPOLES.index(mp) for mp in DESIGNED_MOMENT[kind]]
            pairs = [(i, j) for i in rows for j in range(len(labels))]
        else:                                                        # корректор — все значимые по физ.вкладу в кик
            kick = np.abs(M).max(axis=1) * abs(mt.QP * L) * scale
            kmax = kick.max() or 1.0
            rows = [i for i in range(len(MULTIPOLES)) if kick[i] > 1e-4 * kmax and i != 5] or [int(kick.argmax())]
            kij = np.abs(M) * (abs(mt.QP * L) * scale)[:, None]      # линейность по ДОМИНАНТНОМУ (по кику) выходу обмотки
            pairs = [(int(kij[:, j].argmax()), j) for j in range(len(labels))]
        imax, r2 = _relevant_imax(curs, [Vals[i, j] for i, j in pairs])
        combo = combo_verify(base, devs, L, M, rows, imax) if len(devs) > 1 else None   # суперпозиция многообмоточных
        res.append(dict(base=base, kind=kind, labels=labels, M=M, rows=rows, L=L, r2=r2,
                        imax=imax, i_op=i_op, combo=combo, pm=False))
    return res


def _load_combo(base):
    """Возврат (X:токи Nx3, Y:выходы NxK, syms:имена выходов). Выходы — по ЗАГОЛОВКУ (B,G или B,G,S)."""
    fv = lambda s: float(s.replace(",", "."))
    rows = list(csv.reader(open(os.path.join(OUTDIR, f"combo_dataset_{base}.csv"), encoding="utf-8-sig"), delimiter=";"))
    hdr = next(r for r in rows if r and r[0] == "I1_А")
    syms = [c for c in hdr[3:] if c in ("B", "G", "S", "O")]           # выходные коэф. (без RMS)
    d = [r for r in rows if len(r) >= 4 + len(syms) and r[0] and r[0][0] in "+-0123456789" and "," in r[0]]
    X = np.array([[fv(r[0]), fv(r[1]), fv(r[2])] for r in d])
    Y = np.array([[fv(r[3 + j]) for j in range(len(syms))] for r in d])
    return X, Y, syms


def _load_single(dev):
    fv = lambda s: float(s.replace(",", "."))
    rows = list(csv.reader(open(os.path.join(OUTDIR, f"dataset_{dev}.csv"), encoding="utf-8-sig"), delimiter=";"))
    d = [r for r in rows if r and r[0] and r[0][0] in "+-0123456789" and "," in r[0]]
    return np.array([fv(r[0]) for r in d]), np.array([fv(r[1]) for r in d])


def _r2(y, yhat):
    return 1.0 - ((y - yhat) ** 2).sum() / max(1e-30, ((y - y.mean()) ** 2).sum())


def compute_coupling():
    """Матрицы связи ток→физкоэф (МНК по датасетам element_compare). DQ: (I1,I2,I3)→(B,G[,S]), матрица Kx3
    (K=число выходов в датасете). Одиночные: coef = k·I (скаляр). Возврат dict(dq=[…], single=[…])."""
    out = {"dq": [], "single": []}
    for base in ("L400", "L1000"):
        if not os.path.exists(os.path.join(OUTDIR, f"combo_dataset_{base}.csv")):
            continue
        X, Y, syms = _load_combo(base)
        MT = np.linalg.lstsq(X, Y, rcond=None)[0]                    # Y = X·MTᵀ, MT: 3×K
        R2 = [_r2(Y[:, k], X @ MT[:, k]) for k in range(Y.shape[1])]
        out["dq"].append(dict(base=base, M=MT.T, R2=R2, syms=syms, Imax=float(np.abs(X).max())))
    dom = {"quad": "G", "sext": "S", "oct": "O"}
    for dev, kind in candidates(("quad", "sext", "oct")):
        if not os.path.exists(os.path.join(OUTDIR, f"dataset_{dev}.csv")):
            continue
        I, C = _load_single(dev)
        k, b = np.polyfit(I, C, 1)
        out["single"].append(dict(dev=dev, kind=kind, sym=dom[kind], k=float(k),
                                  R2=_r2(C, k * I + b), Imax=float(np.abs(I).max())))
    return out


def _fmt_sci(x, sig=3):
    """Число → LaTeX научной записи: 5.746e-5 → '5.75\\cdot10^{-5}', ~0 → '0'."""
    if abs(x) < 1e-30:
        return "0"
    e = int(math.floor(math.log10(abs(x)))); m = x / 10 ** e
    ms = f"{m:.{sig-1}f}"
    return ms if e == 0 else f"{ms}\\cdot10^{{{e}}}"


def _matrix_latex(M, label="M"):
    """ОДНА числовая матрица (форма по M.shape) с меткой. Office MML2OMML.XSL надёжно рендерит лишь ОДИН
    \\begin{pmatrix} на формулу — вектор-уравнение даём текстом, не формулой."""
    nr, nc = M.shape
    body = r"\\".join(" & ".join(_fmt_sci(M[i, j]) for j in range(nc)) for i in range(nr))
    return label + r"=\begin{pmatrix}" + body + r"\end{pmatrix}"


def _sym_matrix_latex(nr, nc=3, label="M"):
    """Символьная матрица nr×nc: M=(m_{ij})."""
    body = r"\\".join(" & ".join(f"m_{{{i+1}{j+1}}}" for j in range(nc)) for i in range(nr))
    return label + r"=\begin{pmatrix}" + body + r"\end{pmatrix}"


def _drop_section(doc, sentinel):
    """Удалить абзацы/таблицы от заголовка sentinel до конца (кроме финального sectPr) — идемпотентность."""
    from docx.oxml.ns import qn
    body = doc.element.body; kids = list(body); start = None
    for i, el in enumerate(kids):
        if el.tag == qn("w:p") and sentinel in "".join(n.text or "" for n in el.iter(qn("w:t"))):
            start = i; break
    if start is None:
        return False
    for el in kids[start:]:
        if el.tag != qn("w:sectPr"):
            body.remove(el)
    return True


# Заголовки, которыми начинаются наши приложения (сводка/корректор/секции сводки). Всё от ПЕРВОГО такого
# заголовка до конца документа — приложения матриц (идут ПОСЛЕ реальных глав элементов из run_all).
MATRIX_SENT = "Сводка: матрицы связи и диапазоны валидности"
MATRIX_CORR = "Матрицы связи ток→физкоэф (элементы без отдельной главы)"
_ARTIFACT_HEADINGS = {MATRIX_SENT, MATRIX_CORR, "Матрицы связи токов и физкоэффициентов",
                      "Обозначения", "Диапазоны тестов и валидности", "Все матрицы", "Замечания"}


_REAL_CHAPTER = __import__("re").compile(r"\((oct|dipole|dq|quad|sext|kicker)\)\s*$")


def _purge_matrix_artifacts(doc):
    """Надёжная идемпотентность: удалить ВСЕ прежние приложения матриц — от первого артефактного заголовка до
    конца документа. Артефакт = заголовок из фикс. набора / с «, обмоток:» ЛИБО любая H1 ПОСЛЕ реальных глав
    элементов, не являющаяся главой элемента (не оканчивается на «(oct|dipole|dq|quad|sext|kicker)»). Ловит
    приложения любого прошлого названия (напр. переименованную сборную главу) и осиротевшие секции."""
    from docx.oxml.ns import qn
    seen_real = False; target = None
    for p in doc.paragraphs:
        st = p.style.name; txt = p.text.strip()
        is_h1 = (st == "Heading 1"); is_head = st.startswith(("Heading", "Заголовок"))
        if is_h1 and _REAL_CHAPTER.search(txt):
            seen_real = True; continue
        artifact = (is_head and (txt in _ARTIFACT_HEADINGS or ", обмоток:" in p.text)) \
            or (is_h1 and seen_real and txt and not _REAL_CHAPTER.search(txt))
        if artifact:
            target = p._p; break
    if target is None:
        return
    body = doc.element.body; kids = list(body)
    try:
        start = kids.index(target)
    except ValueError:
        return
    for el in kids[start:]:
        if el.tag != qn("w:sectPr"):
            body.remove(el)


def _ltag(name):
    import re
    m = re.search(r"L\d+", name); return m.group(0) if m else name


def _curname(lbl):
    """Имя тока в вектор-уравнении: I1/I2/I3 как есть, c1/main → I_c1/I_main, одиночная → I."""
    return lbl if lbl.startswith("I") else ("I_" + lbl)


def _match_ltag(text, ltag):
    """ltag в тексте с границей: за ним НЕ цифра (иначе L100 матчит L1000)."""
    import re
    return re.search(re.escape(ltag) + r"(?!\d)", text) is not None


def _has_chapter(doc, ltag):
    return any(p.style.name == "Heading 1" and _match_ltag(p.text, ltag) for p in doc.paragraphs)


def _chapter_next_p(doc, ltag):
    """_p следующего H1 после главы с ltag (точка вставки в конец главы). None если главы/след. H1 нет."""
    ps = doc.paragraphs; idx = None
    for i, p in enumerate(ps):
        if p.style.name == "Heading 1" and _match_ltag(p.text, ltag): idx = i; break
    if idx is None: return None
    for j in range(idx + 1, len(ps)):
        if ps[j].style.name == "Heading 1": return ps[j]._p
    return None


def _drop_inserted_matrices(doc):
    """Убрать ранее вставленные в главы блоки матриц (прямая+M, обратная+M⁺, диапазоны = 5 абзацев) —
    идемпотентность. Совместимо со старым 3-абзацным блоком (лишние берутся до след. не-нашего абзаца)."""
    from docx.oxml.ns import qn
    body = doc.element.body; kids = list(body); rm = []
    for i, el in enumerate(kids):
        if el.tag != qn("w:p"): continue
        txt = "".join(n.text or "" for n in el.iter(qn("w:t")))
        if txt.startswith("Матрица связи (ток→коэф)"):
            rm += [k for k in kids[i:i + 5] if k.tag == qn("w:p")]
        elif txt.startswith("PM-диполь — ФИКСИРОВАННОЕ"):             # PM-блок = 1 абзац
            rm.append(el)
    for el in rm:
        if el.getparent() is not None: body.remove(el)


def _matrix_block(doc, m, before, lo):
    """Блок матрицы (5 абзацев: прямое вектор-уравнение + M, обратное вектор-уравнение + M⁺/M⁻¹, диапазоны).
    before=None → в конец. Прямая: коэф = M·ток; обратная: ток = M⁺·коэф (M⁻¹ для квадратной невырожденной).
    PM-диполь (fixed): один абзац с фикс. составом (матрицы связи нет — свипа тока нет)."""
    from docx.shared import Pt
    if m.get("pm"):                                                  # PM: фикс. поле, матрицы нет
        fx = ";  ".join(f"{k}={v:+.3e}" for k, v in m["fixed"].items())
        p = doc.add_paragraph(f"PM-диполь — ФИКСИРОВАННОЕ поле (свипа тока нет ⇒ матрицы связи ток→коэф нет). "
                              f"Мультипольный состав из карты (диффур): {fx}.  {m.get('note', '')}")
        p.runs[0].font.size = Pt(9)
        if before is not None: before.addprevious(p._p)
        return
    ov = ", ".join(m["outputs"]); iv = ", ".join(_curname(l) for l in m["labels"])
    M = np.asarray(m["M"], float); nr, nc = M.shape
    imax = float(m.get("imax", 0.0)); r2 = m.get("r2", 1.0)
    Minv = np.linalg.pinv(M)                                          # псевдообратная (= M⁻¹ для квадратной невырожденной)
    square = (nr == nc and abs(np.linalg.det(M)) > 1e-30)
    ilab = "M^{-1}" if square else "M^{+}"
    ihint = ("" if square else
             (" (M⁺ — псевдообратная минимальной нормы: обмоток больше моментов)" if nc > nr
              else " (M⁺ — псевдообратная МНК: моментов больше обмоток)"))
    crng = ";  ".join(f"|{m['outputs'][i]}|≤{np.abs(M[i]).sum() * imax:.2e}" for i in range(nr))
    valid = (f"весь тест-диапазон (R²={r2:.3f}≥0.95)" if r2 >= 0.95
             else f"УЖЕ — R²={r2:.3f}<0.95 (связь нелинейна выше)")
    cb = m.get("combo")
    combo_txt = (f"  Проверка суперпозиции ({cb['n']} комбо обмоток в ±0.5·Imax; отн. невязка матрица↔диффур по "
                 f"моментам): " + ", ".join(f"{k} {v * 100:.1f}%" for k, v in cb["per"].items()) + "." if cb else "")
    p1 = doc.add_paragraph(f"Матрица связи (ток→коэф):  ({ov})ᵀ = M · ({iv})ᵀ  — строки M: {ov};  столбцы: {iv}.")
    p2 = doc.add_paragraph(); lo.add_formula(p2, _matrix_latex(M, "M"))
    p3 = doc.add_paragraph(f"Обратный расчёт (коэф→ток):  ({iv})ᵀ = {'M⁻¹' if square else 'M⁺'} · ({ov})ᵀ.{ihint}")
    p4 = doc.add_paragraph(); lo.add_formula(p4, _matrix_latex(Minv, ilab))
    p5 = doc.add_paragraph(f"Релевантный ток (критерий R²≥0.95): |ток| ≤ {imax:.1f} А — {valid}.  "
                           f"Диапазон коэф.: {crng}.  {m.get('note', '')}{combo_txt}")
    p5.runs[0].font.size = Pt(9)
    if before is not None:
        for p in (p1, p2, p3, p4, p5): before.addprevious(p._p)


def _coupling_summary(doc, mats, lo):
    """Сводная глава: обозначения, таблица тест-диапазонов (ток→коэф = валидность), все матрицы, замечания."""
    from docx.shared import Pt
    doc.add_paragraph("Все матрицы связи ток→физкоэф в одном месте. Указан диапазон ТОКОВ, в котором проведены "
                      "тесты, и соответствующий диапазон ФИЗКОЭФ (в этих пределах линейная связь справедлива). "
                      "Прямой расчёт: (моменты)ᵀ = M · (токи обмоток)ᵀ, элемент M_ij = ∂(коэф_i)/∂(ток_j).  "
                      "Обратный расчёт (требуемые коэф → токи): (токи)ᵀ = M⁺ · (моменты)ᵀ, где M⁺ — псевдообратная "
                      "Мура–Пенроуза (= M⁻¹ для квадратной невырожденной M; при избытке обмоток — решение "
                      "минимальной нормы, при недостатке — МНК). Числовая M⁺/M⁻¹ приведена в разделе каждого элемента.")
    doc.add_heading("Обозначения", level=2)
    leg = [("I, I₁,I₂,I₃, I_c1..main", "токи обмоток, А (одиночная / секступоль-DQ / корректор)"),
           ("B, Bv", "гориз./верт. дипольный момент (mat_B; Bv извлечён из вертикального кика)"),
           ("G, Gv", "норм./скью квадрупольный градиент (mat_G / mat_Gv)"),
           ("S, O", "секступольный / октупольный коэф. (mat_S / mat_O)"),
           ("M", "прямая матрица: строки — моменты, столбцы — токи обмоток; [коэф]/А"),
           ("M⁻¹, M⁺", "обратная матрица (коэф→ток): (токи)ᵀ = M⁺·(моменты)ᵀ; M⁺ = M⁻¹ для квадратной невырожденной"),
           ("QP", "q/p₀ — опорная жёсткость (≈ −84.9); θ=QP·B·L, k₁=QP·G, k₂=QP·S")]
    t = doc.add_table(rows=1, cols=2); t.style = "Light Grid Accent 1"
    t.rows[0].cells[0].text = "символ"; t.rows[0].cells[1].text = "смысл"
    for s, mn in leg:
        c = t.add_row().cells; c[0].text = s; c[1].text = mn
    doc.add_heading("Диапазоны тестов и валидности", level=2)
    doc.add_paragraph("Релевантный диапазон тока — тот, в котором линейная связь справедлива. КРИТЕРИЙ: R² ≥ 0.95 "
                      "(R² линейного фита коэф.=slope·ток по тест-точкам от рабочей точки вверх). Ток свипа привязан "
                      "к РАБОЧЕЙ точке каждого элемента (не фиксирован): у слабых квадов рабочий ток ~2 А, гнать их "
                      "до десятков А физически бессмысленно (перегруз, насыщение). Релевантный |ток| — максимальный, "
                      "при котором R² ещё ≥ 0.95; в этих пределах матрица (и обратная к ней) справедлива.")
    st = doc.add_table(rows=1, cols=6); st.style = "Light Grid Accent 1"
    for j, h in enumerate(("элемент", "обмотки", "моменты", "релевантный ток, А (критерий R²≥0.95)", "R²", "диапазон |коэф|")):
        st.rows[0].cells[j].text = h
    for m in mats:
        c = st.add_row().cells
        c[0].text = m["name"]; c[2].text = ", ".join(m["outputs"])
        if m.get("pm"):                                              # PM: фикс. поле, диапазона нет
            c[1].text = "— (PM)"; c[3].text = "фикс. поле"; c[4].text = "—"
            c[5].text = ";  ".join(f"{k}={v:+.1e}" for k, v in m["fixed"].items())
            continue
        imax = float(m.get("imax", 0.0)); r2 = m.get("r2", 1.0)
        c[1].text = ", ".join(m["labels"])
        c[3].text = f"|I|≤{imax:.1f}" + ("" if r2 >= 0.95 else "  (R²<0.95 — уже)")
        c[4].text = f"{r2:.3f}"
        c[5].text = ";  ".join(f"|{m['outputs'][i]}|≤{np.abs(m['M'][i]).sum() * imax:.1e}" for i in range(len(m["outputs"])))
    doc.add_heading("Все матрицы", level=2)
    for m in mats:
        doc.add_heading(m["name"], level=3)
        _matrix_block(doc, m, None, lo)
    doc.add_heading("Замечания", level=2)
    for txt in (
        "Связь линейна (суперпозиция полей обмоток), поэтому столбцы матрицы получены НЕЗАВИСИМЫМ свипом каждой "
        "обмотки. Для многообмоточных (секступоли, корректор) суперпозиция ПОДТВЕРЖДЕНА комбо-сеткой токов: "
        "предсказание матрицы M·I сверено с диффуром на комбинациях возбуждения (строка «Проверка суперпозиции» "
        "в разделе элемента, датасеты combo_verify_<элемент>.csv). Секступоли — невязка <0.1%. Корректор: сильные "
        "моменты (гориз./верт. диполь) сверяются на 1–4%, а слабые квадрупольный и особенно секступольный моменты "
        "у корректора на ШУМОВОМ ПОЛУ карты (невязка не убывает с амплитудой) — они мало влияют на коррекцию.",
        "Вертикальный диполь Bv извлечён из вертикального кика (столбца mat_Bv в БД нет). Октуполь O показан "
        "только у октуполя (у прочих — шум кубического члена при 1 мм).",
        "Вырождения нормальны (обе обмотки секступоля дают один момент; у DQ I₂≡I₃): прямая задача ток→коэф "
        "корректна, ОБРАТНАЯ коэф→ток недоопределена — псевдообратная M⁺ даёт решение минимальной нормы "
        "(ток распределяется между обмотками по их чувствительности); при недостатке обмоток M⁺ — МНК-решение.",
        "DQ смоделирован как B+G (без секступоля = elegant/.lte K2=0); у L1000 на рабочем токе линейность "
        "ограничена (R²≈0.7, при |I|≳2.5 А нужна интерполяция/квадратичные члены)."):
        doc.add_paragraph(txt, style="List Bullet")


# ---- Унифицированный список матриц + таблица БД coupling_matrices (прямые + обратные, upsert) ----
COUPLING_TABLE = "coupling_matrices.csv"
COUPLING_HEADER = ["element", "kind", "n_windings", "n_moments", "direction",
                   "out_label", "in_label", "value", "units", "eff_length", "i_op", "imax", "r2"]


def _num(x):
    """Число → строка с ТОЧКОЙ (машинная таблица, память csv-comma-decimal); пусто/None → '-'."""
    return "-" if x in ("", None) else repr(float(x))


def build_all_mats(elements=None, dq_bases=None):
    """Унифицированный список матриц связи: прямые элементы (мультиполи из кика) + DQ (combo B+G) + PM (fixed,
    состав карты). elements=None → library_windings() (весь охват); dq_bases=None → все DQ, [] → без DQ,
    [<base>] → только эти. Каждый dict: name/kind/ltag/outputs/labels + либо (M,imax,i_op,r2,combo,L)
    [токовый], либо (fixed,L,pm=True) [PM]."""
    els = library_windings() if elements is None else elements
    data_el = compute_element_matrices(els)
    mats = []
    for el in data_el:
        if el.get("pm"):
            mats.append(dict(name=el["base"], kind=el["kind"], ltag=_ltag(el["base"]), pm=True, L=el["L"],
                             outputs=[MULTIPOLES[i] for i in el["rows"]], labels=el["labels"],
                             fixed=el["fixed"], note="PM: фикс. поле (свипа тока нет); состав из карты."))
        else:
            rows = el["rows"]
            mats.append(dict(name=el["base"], kind=el["kind"], ltag=_ltag(el["base"]), pm=False, L=el["L"],
                             outputs=[MULTIPOLES[i] for i in rows], labels=el["labels"],
                             M=el["M"][rows], imax=el["imax"], i_op=el["i_op"], r2=el["r2"],
                             combo=el.get("combo"), note=""))
    data = compute_coupling()                                        # DQ (combo-датасеты, дёшево)
    for dq in data["dq"]:
        if dq_bases is not None and dq["base"] not in dq_bases: continue
        r2 = dq["R2"]; sm = dq["syms"]
        L = LCROW.get(f"DQ_I1_1A_{dq['base']}", {}).get("eff_length", "")
        mats.append(dict(name=f"DQ_I1_1A_{dq['base']}", kind="dq", ltag=dq["base"], pm=False, L=L,
                         outputs=sm, labels=["I1", "I2", "I3"], M=dq["M"], imax=dq["Imax"], i_op=dq["Imax"],
                         r2=float(min(r2)), combo=None,
                         note="R²: " + ", ".join(f"{sm[k]}={r2[k]:.3f}" for k in range(len(sm))) +
                              f" (рабочий ток ~{dq['Imax']:.1f}А)."))
    return mats


def _matrix_table_rows(m):
    """Строки таблицы coupling_matrices для одного элемента: forward (коэф←ток) + inverse (ток←коэф);
    PM → fixed. value с ТОЧКОЙ. Восстановление матрицы: группировка по (element,direction)."""
    el, kind = m["name"], m["kind"]
    outs, labs = m["outputs"], m["labels"]
    nw, nmo = len(labs), len(outs); L = _num(m.get("L", ""))
    if m.get("pm"):
        return [[el, kind, "0", str(nmo), "fixed", mp, "-", repr(float(m["fixed"][mp])),
                 "coef", L, "-", "-", "-"] for mp in outs]
    M = np.asarray(m["M"], float); Minv = np.linalg.pinv(M)          # forward (nmo×nw), inverse (nw×nmo)
    iop, imax, r2 = _num(m.get("i_op", "")), _num(m.get("imax", "")), _num(m.get("r2", ""))
    rows = []
    for i, mp in enumerate(outs):                                    # прямая: момент ← токи обмоток
        for j, lb in enumerate(labs):
            rows.append([el, kind, str(nw), str(nmo), "forward", mp, _curname(lb),
                         repr(float(M[i, j])), "coef/A", L, iop, imax, r2])
    for j, lb in enumerate(labs):                                    # обратная: ток обмотки ← моменты
        for i, mp in enumerate(outs):
            rows.append([el, kind, str(nw), str(nmo), "inverse", _curname(lb), mp,
                         repr(float(Minv[j, i])), "A/coef", L, iop, imax, r2])
    return rows


def upsert_coupling_table(rows_by_element):
    """UPSERT в workspace/tables/coupling_matrices.csv: нет файла → создать; есть → удалить строки ТОЛЬКО
    обновляемых элементов (колонка element) и дописать новые. rows_by_element: dict element→list[row]."""
    path = os.path.join(TB, COUPLING_TABLE)
    existing = []
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            rd = list(csv.reader(f, delimiter=";"))
        existing = rd[1:] if rd else []
    upd = set(rows_by_element)
    keep = [r for r in existing if r and r[0] not in upd]            # чужие элементы не трогаем
    new = keep + [row for el in rows_by_element for row in rows_by_element[el]]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(COUPLING_HEADER); w.writerows(new)
    return path, len(upd), len(new)


def write_coupling_table(mats):
    """Записать/обновить таблицу БД coupling_matrices из списка mats (прямые+обратные / fixed для PM)."""
    path, nupd, ntot = upsert_coupling_table({m["name"]: _matrix_table_rows(m) for m in mats})
    print(f"  таблица БД: обновлено элементов {nupd}, всего строк {ntot} → workspace/tables/{COUPLING_TABLE}")
    return path


def build_coupling_report(path=None, mats=None):
    """Матрицы связи ток→физкоэф в Word-отчёте: (1) блок матрицы (прямая M + обратная M⁺/M⁻¹ + диапазоны) в
    главу КАЖДОГО элемента; элементы без своей главы — сборная глава; (2) сводная глава в конце (все матрицы +
    тест-диапазоны токов/коэф = валидность, критерий R²≥0.95). Нативные формулы Word. Идемпотентно.
    mats=None → строит по всему охвату (library) через build_all_mats."""
    from docx import Document
    from docx.shared import Pt
    import latex2omml as lo
    path = path or os.path.join(OUTDIR, "element_compare_report.docx")
    if mats is None:
        mats = build_all_mats()
    doc = Document(path) if os.path.exists(path) else Document()
    SENT = MATRIX_SENT
    CORR = MATRIX_CORR
    _drop_inserted_matrices(doc)                                 # снять вставленные в главы блоки матриц (маркер)
    _purge_matrix_artifacts(doc)                                 # снять все прежние приложения (сводка/корректор/сироты)
    # элементы БЕЗ своей главы в отчёте (квады, секступоли, корректор, DQ_L400) → сборная глава с матрицами
    orphan = [m for m in mats if not _has_chapter(doc, m["ltag"])]
    if orphan:
        doc.add_heading(CORR, level=1)
        doc.add_paragraph("Прямая (ток→коэф) и обратная (коэф→ток) матрицы для элементов, у которых нет "
                          "собственной главы в отчёте: квадруполи (→G), секступоли (2 обмотки →S), корректор, DQ, "
                          "PM-диполи. Корректор (kicker) многофункционален: обмотки дают горизонтальный и вертикальный "
                          "дипольные моменты, нормальный и скью квадруполь, секступоль. PM-диполи — фикс. поле "
                          "(матрицы связи нет, приведён извлечённый из карты состав).")
        for m in orphan:
            tag = "PM, фикс. поле" if m.get("pm") else f"обмоток: {len(m['labels'])}"
            doc.add_heading(f"{m['name']}  ({m['kind']}, {tag})", level=2)
            _matrix_block(doc, m, None, lo)
    # сводная глава В КОНЦЕ (все матрицы + диапазоны)
    doc.add_heading(SENT, level=1)
    _coupling_summary(doc, mats, lo)
    # разнести матрицы по главам своих элементов (у кого глава в отчёте есть) — вектор-уравнение + M + диапазоны
    for m in mats:
        nxt = _chapter_next_p(doc, m["ltag"])
        if nxt is not None:
            _matrix_block(doc, m, nxt, lo)
    try:
        doc.save(path)
    except PermissionError:                              # отчёт открыт в Word → писать в фолбэк
        path = os.path.splitext(path)[0] + "_new.docx"
        doc.save(path)
        print(f"  (основной .docx заблокирован — записан {os.path.basename(path)}; закрой Word и переименуй)")
    print(f"матрицы разнесены по главам + сводка «{SENT}» → {path}")
    return path


# ---- Единый конвейер: матрицы для отдельного элемента / всей библиотеки + таблица БД + demo + Word + HTML ----
def _resolve_single(name):
    """<element> → (elements для compute_element_matrices, dq_bases). Prime/PM элемент по base ИЛИ DQ (L400/L1000)."""
    for tup in library_windings():
        if tup[0] == name:
            return [tup], []
    for b in ("L400", "L1000"):
        if name == b or name.endswith(b):
            return [], [b]
    raise SystemExit(f"элемент '{name}' не в охвате. Ожидается имя из library (напр. QP_L_1A_L295, "
                     f"SP_1A_L166, corrector_100A_L100, DP_sh_PM_L49) или DQ (L400 / L1000).")


def print_matrices(mats):
    """Печать прямой M и обратной M⁺/M⁻¹ каждого элемента (single-режим)."""
    for m in mats:
        if m.get("pm"):
            print(f"  {m['name']} (PM, фикс. поле): " + ", ".join(f"{k}={v:+.4e}" for k, v in m["fixed"].items()))
            continue
        M = np.asarray(m["M"], float); Minv = np.linalg.pinv(M)
        print(f"  {m['name']} ({m['kind']}): моменты={m['outputs']} обмотки={m['labels']} "
              f"Imax={m.get('imax', 0):.1f}А R²={m.get('r2', 1):.3f}")
        print("    M (коэф←ток):   " + " | ".join("  ".join(f"{x:+.3e}" for x in row) for row in M))
        print("    M⁺ (ток←коэф):  " + " | ".join("  ".join(f"{x:+.3e}" for x in row) for row in Minv))


def _config_device(eid):
    """device_type_id первого инстанса конфига lattice_id=eid (для имени HTML-файла)."""
    _, rows = _rows("lattice_configs.csv")
    for r in rows:
        if r and r[0] == str(eid):
            return r[2]
    return None


def generate_html_reports(ids=None):
    """HTML-отчёты (viz_html: детальный 3D-магнит + траектории в относит. координатах внизу) для демо-конфигов.
    ids=None → все id≥VIZ_START из experiments.csv. ПОЛНЫЙ STL (glyph=False — один магнит, показываем в
    качестве) + РЕАЛЬНЫЙ пучок конфига (config_beam=True — те же 13 частиц, что диффур). Устойчиво к сбоям."""
    os.makedirs(HTML_OUT, exist_ok=True)
    try:
        import viz_html                                            # тяжёлые pyvista/plotly — импорт лениво
    except Exception as e:
        print(f"  ! HTML пропущен (viz_html не импортируется: {type(e).__name__}: {e})"); return []
    if ids is None:
        _, ex = _rows("experiments.csv")
        ids = sorted({int(r[0]) for r in ex if r and r[0].isdigit() and int(r[0]) >= VIZ_START})
    done = []
    for eid in ids:
        dev = _config_device(eid) or str(eid)
        out = os.path.join(HTML_OUT, f"element_{dev}.html")
        try:
            viz_html.main(eid, glyph=False, out=out, config_beam=True); done.append(out)
            print(f"  HTML id={eid} {dev} → {os.path.relpath(out, ROOT)}")
        except Exception as e:
            print(f"  ! HTML id={eid} {dev}: {type(e).__name__}: {e}")
    return done


def build_library(single=None, do_report=None, do_demo=None, do_html=None):
    """Единый конвейер матриц связи. single=None → ВСЯ библиотека (матрицы+датасеты+таблица БД+Word+demo+HTML);
    single=<element> → ОДИН элемент (матрица+датасет+UPSERT его записи+печать M/M⁺; Word/demo/HTML — по флагам,
    и они пересобирают ОБЩИЕ артефакты). do_report/demo/html: None → all=True, single=False."""
    allmode = (single is None)
    if do_report is None: do_report = allmode
    if do_demo is None: do_demo = allmode
    if do_html is None: do_html = allmode
    els, dq_bases = (None, None) if allmode else _resolve_single(single)
    print(f"=== матрицы связи ток→физкоэф: {'ВСЯ БИБЛИОТЕКА' if allmode else single} ===")
    mats = build_all_mats(els, dq_bases)                          # матрицы (+ свип/комбо датасеты внутри)
    write_coupling_table(mats)                                    # таблица БД (upsert)
    if not allmode:
        print_matrices(mats)
    if do_report:                                                # Word — всегда полный охват (общий артефакт)
        build_coupling_report(mats=(mats if allmode else build_all_mats()))
    ids = None
    if do_demo:                                                  # демо-конфиги (setup_viz_all — весь набор)
        try:
            ids = [s["id"] for s in setup_viz_all()]
        except Exception as e:
            print(f"  ! демо-конфиги: {type(e).__name__}: {e}")
    if do_html:
        generate_html_reports(ids)
    print(f"=== готово: элементов {len(mats)} ===")
    return mats


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "library":
        a = sys.argv[2:]; flags = {x for x in a if x.startswith("--")}; pos = [x for x in a if not x.startswith("--")]
        kw = {}
        for name, key in (("report", "do_report"), ("demo", "do_demo"), ("html", "do_html")):
            if f"--{name}" in flags: kw[key] = True
            if f"--no-{name}" in flags: kw[key] = False
        build_library(pos[0] if pos else None, **kw); sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "matrices":
        build_coupling_report(); sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "viz-all":
        setup_viz_all(); sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "viz":
        dev = sys.argv[2] if len(sys.argv) > 2 else "DQ_I1_1A_L400"
        vid = int(sys.argv[4]) if len(sys.argv) > 4 else VIZ_START
        kind = LCROW[dev].get("kind") or "quad"
        if len(sys.argv) > 3:                              # явный ток: viz <dev> <ток> [id]
            setup_viz(dev, float(sys.argv[3]), vid)
        elif is_pm(dev):                                   # проектный ток: viz <dev> [_ id]
            setup_viz(dev, "-", vid)
        else:
            recs = sweep_currents(dev, kind, PROBES)
            setup_viz(dev, design_current(dev, kind, recs), vid, Leff=recs[0]["Leff_common"])
        print(f"  визуализировать:  py -3.14 visualize.py {vid}")
        sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "combo":
        base = sys.argv[2] if len(sys.argv) > 2 else "L400"
        ng = int(sys.argv[3]) if len(sys.argv) > 3 else 4
        run_dq_combo(base, ng); sys.exit(0)
    if len(sys.argv) > 1 and sys.argv[1] == "all":
        kinds = tuple(sys.argv[2].split(",")) if len(sys.argv) > 2 else ("quad", "sext", "oct")
        run_all(kinds); sys.exit(0)
    dev = sys.argv[1] if len(sys.argv) > 1 else "QP_L_1A_L295"
    cur = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    kind = {r["device_type_id"]: r.get("kind") for r in mt._read("lattice_configs.csv") if r.get("kind")}.get(dev, "quad")
    probes = [(0.0, 0.0)] + [(1e-3 * math.cos(2 * math.pi * k / 12), 1e-3 * math.sin(2 * math.pi * k / 12)) for k in range(12)]
    bak = os.path.join(TB, "_cmp_backup"); os.makedirs(bak, exist_ok=True)   # страховка
    for f in CSV4: shutil.copy2(os.path.join(TB, f), os.path.join(bak, f))
    try:
        trajs, geom = run_diffur(dev, cur)
        diff_su = diff_to_su(trajs, geom["xstart"])
        print(f"диффур {dev} ({kind}) @ {cur}А: {len(trajs)} частиц, точек/частицу {len(trajs[0])}")
        coeffs, Leff, mse = fit_element(dev, cur, diff_su, geom, kind, probes)
        cs = "  ".join(f"{k}={v:+.4e}" for k, v in coeffs.items())
        print(f"ФИТ кика: {cs}  L_eff={Leff:.4f}м  (старт .lte {HW[dev].get('eff_length')}м)")
        print(f"невязка траекторий: RMS={math.sqrt(mse)*1e3:.4f} мм  (MSE={mse:.3e} м²)")
        # финальные позиции: диффур vs кик (проверка)
        mat = track_set(dev, coeffs, Leff, kind, geom, probes)
        dfin = np.array([_resample(d) for d in diff_su])[:, :, -1]   # (u,v) на финише
        mfin = np.array([_resample(m) for m in mat])[:, :, -1]
        worst = np.max(np.abs(dfin - mfin)) * 1e3
        print(f"макс. расхождение финальных позиций (13 частиц): {worst:.4f} мм")
    finally:
        cleanup()
        for f in CSV4:                                    # восстановить на случай сбоя очистки
            b = os.path.join(bak, f)
            if os.path.exists(b) and os.path.getsize(os.path.join(TB, f)) == 0:
                shutil.copy2(b, os.path.join(TB, f))
