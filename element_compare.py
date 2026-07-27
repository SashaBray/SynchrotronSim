# -*- coding: utf-8 -*-
"""Сравнение ДИФФУРНОЙ (C++ SynchrotronTracker, Boris по картам поля) и МАТРИЧНОЙ (Python matrix_tracker)
моделей по одиночным элементам БД. Датасет траекторий при разных токах + подгонка матрицы (L_eff+физкоэф).

ФАЗА 1 (этот файл пока): фундамент — собрать ВРЕМЕННЫЙ одно-элементный конфиг в основной БД
(пучок → элемент → финиш-плоскость; дрейфы = свободный пролёт, длины НАШИ, из диффура не берём),
прогнать C++-диффур, прочитать траекторию, вычистить временные строки. Дальше — свип токов, фит матрицы,
Word-отчёт. Диффур eff_length НЕ использует (трек по карте поля) ⇒ спека «без диффурных длин» соблюдена.
"""
import os, csv, math, struct, subprocess, shutil, sys
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


def make_beam(x0, radius=1e-3, n_circle=12):
    """13 частиц у стартовой плоскости x=x0: 12 по кругу Ø2·radius в поперечной (y,z) + 1 в центре.
    Все летят вдоль +x со скоростью пучка. Возврат список строк particle_groups."""
    parts = [(0.0, 0.0)]
    for k in range(n_circle):
        a = 2 * math.pi * k / n_circle
        parts.append((radius * math.cos(a), radius * math.sin(a)))
    rows = []
    for pid, (y, z) in enumerate(parts, start=1):
        rows.append([str(TMP_GRP), str(pid), repr(ME), repr(QE),
                     repr(x0), repr(y), repr(z), repr(V_BEAM), "0.0", "0.0"])
    return rows


def build_config(dev, cur, yaw=0.0):
    """Записать ВРЕМЕННЫЙ одно-элементный конфиг (dev при токе cur) в 4 CSV. Возврат dict с геометрией."""
    lat_h, lat_r = _rows("lattice_configs.csv")
    con_h, con_r = _rows("element_connectivity.csv")
    grp_h, grp_r = _rows("particle_groups.csv")
    exp_h, exp_r = _rows("experiments.csv")
    ci = {c: i for i, c in enumerate(lat_h)}
    # строка элемента (позиция — центр CENTER по x, ось +x; поворот yaw для изгибов)
    row = [""] * len(lat_h)
    row[ci["lattice_id"]] = str(TMP_LAT); row[ci["instance_id"]] = "1"; row[ci["device_type_id"]] = dev
    row[ci["x"]] = repr(CENTER); row[ci["y"]] = "0.0"; row[ci["z"]] = "0.0"
    row[ci["yaw"]] = repr(yaw); row[ci["pitch"]] = "0.0"; row[ci["roll"]] = "0.0"
    row[ci["arg_val"]] = repr(cur); row[ci["correct_coeff"]] = "1.0"
    for c in ("mat_B", "mat_G", "mat_S", "mat_O"):
        if c in ci: row[ci[c]] = "0"
    lat_new = [r for r in lat_r if r[0] != str(TMP_LAT)] + [row]
    # связность: один элемент
    coni = {c: i for i, c in enumerate(con_h)}
    crow = [""] * len(con_h)
    crow[coni["lattice_id"]] = str(TMP_LAT); crow[coni["instance_id"]] = "1"
    crow[coni["device_type_id"]] = dev
    con_new = [r for r in con_r if r[0] != str(TMP_LAT)] + [crow]
    # пучок у стартовой плоскости
    grp_new = [r for r in grp_r if r[0] != str(TMP_GRP)] + make_beam(CENTER - ENTRY)
    # опыт: финиш-плоскость на x=CENTER+EXIT, нормаль +x (3 точки + направление)
    xf = CENTER + EXIT
    p1 = (xf, 1.0, 1.0); p2 = (xf, -1.0, 1.0); p3 = (xf, 0.0, -1.0)
    erow = [str(TMP_EXP), str(TMP_LAT), str(TMP_GRP), "1.00E-11", "0", "1",
            f"{p1[0]:.6f}", f"{p1[1]:.6f}", f"{p1[2]:.6f}",
            f"{p2[0]:.6f}", f"{p2[1]:.6f}", f"{p2[2]:.6f}",
            f"{p3[0]:.6f}", f"{p3[1]:.6f}", f"{p3[2]:.6f}", "1.0", "0.0", "0.0"]
    erow += [""] * (len(exp_h) - len(erow))
    exp_new = [r for r in exp_r if r and r[0] != str(TMP_EXP)] + [erow]
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
os.environ.setdefault("SS_SECTOR", "1")                # секторный/кик режим матрицы
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
    ur, vr = _resample(mats[0])                            # матричная опорная (ч0) на SGRID
    out = []
    for m in mats:
        um, vm = _resample(m)
        out.append(np.column_stack([SGRID, um - ur, vm - vr]))
    return out


def diff_local(trajs, xstart, bending=False):
    """Диффур .bin (X,Y,Z глоб.) → на частицу (s, u, v). ПРЯМОЙ: u=Y, v=Z, s=X−xstart (опорная на оси).
    ИЗГИБ: проекция ФРЕНЕ на опорную дугу ч0 — s=ДЛИНА ДУГИ вдоль ч0, u=⊥смещение В ПЛОСКОСТИ ИЗГИБА
    (кольцо X-Y), v=вертикаль (Z−Z_ч0). Именно так работает матричный сектор (отклонение от изогнутой
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
    L0 = float(HW[dev].get("eff_length") or 0.3)       # старт L_eff — .lte/матричная длина (НЕ диффурная)
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


def _show_currents(recs):
    """3 показательных тока: край−середина−край (для изгиба/PM токи не [−50,0,50])."""
    curs = sorted({r["cur"] for r in recs})
    return curs if len(curs) <= 3 else [curs[0], curs[len(curs) // 2], curs[-1]]


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
    """Графики траекторий x/y (диффур сплошная синяя vs матрица пунктир красная, 13 частиц) при показательных
    токах. Матрица — с фитированными физкоэф и ОБЩИМ L_eff. Изгиб: отклонение от опорной. Возврат [(cur, png)]."""
    Lc = recs[0].get("Leff_common", recs[0]["Leff"]); figs = []
    for cur in _show_currents(recs):
        r = next((x for x in recs if x["cur"] == cur), None)
        if r is None: continue
        fig, (axx, axy) = plt.subplots(1, 2, figsize=(11, 4.2))
        mats = _mat_for_plot(dev, r["coeffs"], Lc, kind, r["geom"], probes)
        for d, mat in zip(r["diff_su"], mats):
            axx.plot(d[:, 0], d[:, 1] * 1e3, color="#2b7bba", lw=1.1)   # диффур всегда
            axy.plot(d[:, 0], d[:, 2] * 1e3, color="#2b7bba", lw=1.1)
            if mat is None: continue                                    # матрица разошлась — только диффур
            axx.plot(mat[:, 0], mat[:, 1] * 1e3, color="#d1495b", lw=1.0, ls="--")
            axy.plot(mat[:, 0], mat[:, 2] * 1e3, color="#d1495b", lw=1.0, ls="--")
        ttl2 = " (откл. от опорной)" if kind in ("dipole", "dq") else ""
        for ax in (axx, axy):
            ax.set_xlabel("s, м"); ax.grid(alpha=0.3)
        axx.set_ylabel("x (гориз.), мм"); axy.set_ylabel("y (верт.), мм")
        axx.set_title(f"x-плоскость — ток {cur:+.2f}А (— диффур, -- матрица)")
        axy.set_title(f"y-плоскость{ttl2}")
        fig.suptitle(f"{dev} ({kind})  L_eff={Lc:.4f}м  RMS={r['rms_mm']:.4f}мм", fontsize=11)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        tok = f"{cur:+.1f}".replace("+", "p").replace("-", "m").replace(".", "_")   # только токен тока (не .png!)
        png = os.path.join(outdir, f"traj_{dev}_{tok}.png")
        fig.savefig(png, dpi=110); plt.close(fig); figs.append((cur, png))
    return figs


def plot_xy(dev, kind, recs, probes, outdir):
    """Проекция пучка в плоскости xy ДО элемента (стартовый круг Ø2мм) и ЗА элементом (диффур ○ vs матрица ✕)
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
    a1.set_title(f"за элементом, {r['cur']:+.2f} А  (o диффур, x матрица)")
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
    doc.add_heading("SynchrotronSim — сопоставление матричной и диффурной моделей по элементам", level=0)
    doc.add_paragraph("Диффур (C++ SynchrotronTracker, Boris по картам поля) — эталон. Матрица (Python "
                      "matrix_tracker) подогнана последовательными приближениями: физкоэф на каждый ток + один "
                      "L_eff на элемент. Пучок 13 частиц (12 по кругу Ø2мм + центр). Прямые токовые — свип "
                      "−50..+50 А; изгибающие DQ — вокруг рабочего тока; PM-диполи — одна конфигурация (поле "
                      "фиксировано). Изгибающие сравниваются в проекции Френе (отклонение от опорной дуги).")
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


if __name__ == "__main__":
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
        print(f"ФИТ матрицы: {cs}  L_eff={Leff:.4f}м  (старт .lte {HW[dev].get('eff_length')}м)")
        print(f"невязка траекторий: RMS={math.sqrt(mse)*1e3:.4f} мм  (MSE={mse:.3e} м²)")
        # финальные позиции: диффур vs матрица (проверка)
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
