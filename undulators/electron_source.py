# -*- coding: utf-8 -*-
"""Источник света по ОРБИТЕ электрона (а не по рамке магнита) — для тепловых карт и экранов.

Свет летит туда, куда летит электрон. Поэтому ось конуса ондулятора и веер диполя задаются ЗАМКНУТОЙ
ОРБИТОЙ машины в месте источника, а не установочной рамкой железа. Разница принципиальна при перекосе:
у ондулятора нулевой интеграл поля — наклонённое устройство электрон почти не поворачивает, и пятно на
экране не уходит на D·θ_наклона. Пятно двигает УГОЛ ОРБИТЫ в источнике (сдвиги линз, корректоры), а
коррекция орбиты возвращает его — ровно то, что видят фотонные BPM (`control.ring_model`: photon = x + D·x′).

Как считается:
  1. Python-описание решётки строится из таблиц машины (для --ws — возмущённые lattice_configs/connectivity
     из папки ws, остальное — из БД), отсюда рамки входа/выхода каждого магнита (fr_in/fr_out).
  2. Замкнутая орбита — C++-ядром (продакшн): ring_cpp.CppRing(ws) + ra.closed_orbit; состояния на выходах
     магнитов — track_elements; на входах — тем же геометрическим пролётом, что в ядре (mt.modified_drift).
     Опционально коррекция орбиты (orbit_correct.correct, SVD по корректорам-кикерам).
  3. Глобальные положение и направление электрона: P = O + x·U + y·V, n = (px·U + py·V + ps·Ax)/(1+δ).

ЗАПУСК (сравнение рамки железа и орбиты для групп-ондуляторов):
  py -3.14 undulators/electron_source.py --lattice 27 [--ws <папка>] [--correct]
"""
import os, sys, csv, math, argparse, contextlib
sys.stdout.reconfigure(encoding="utf-8")
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt

WS_TABLES = ("lattice_configs.csv", "element_connectivity.csv", "particle_groups.csv", "rf_cavities.csv", "pulsed_kickers.csv", "kicker_timing.csv")
# Срезы для орбиты: 12 — как ансамбли ошибок (тюн точен до 0,001, отклик коррекции кэшируется по nsl).
# Дефолтные 5 (паритет с elegant) смещают тюн на ~0,03 и η на ~4%.
ORBIT_NSL = 12


@contextlib.contextmanager
def ws_tables(ws):
    """Python-сборка решётки из таблиц машины ws (как C++ buildLine: lattice_configs/connectivity — из ws)."""
    if not ws:
        yield
        return
    orig = mt._read

    def _rd(name):
        p = os.path.join(ws, "tables", name)
        if name in WS_TABLES and os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                return list(csv.DictReader(f, delimiter=";"))
        return orig(name)
    mt._read = _rd
    try:
        yield
    finally:
        mt._read = orig


def _mag_instances(lattice):
    """Списки instance_id по магнитам в порядке трекинга — та же группировка, что в mt.build_line."""
    conf = {(r["lattice_id"], r["instance_id"]): r for r in mt._read("lattice_configs.csv")}
    conn = sorted((r for r in mt._read("element_connectivity.csv") if r["lattice_id"] == str(lattice)),
                  key=lambda r: int(r["instance_id"]))
    out = []; i = 0
    while i < len(conn):
        r = conf[(str(lattice), conn[i]["instance_id"])]; grp = [conn[i]["instance_id"]]
        while i + 1 < len(conn):
            rn = conf[(str(lattice), conn[i + 1]["instance_id"])]
            if (rn["x"], rn["y"], rn["z"]) == (r["x"], r["y"], r["z"]):
                grp.append(conn[i + 1]["instance_id"]); i += 1
            else:
                break
        out.append(grp); i += 1
    return out


def _glob(fr, st):
    """Глобальные положение и единичное направление частицы по её состоянию st в рамке fr=(O, Ax, U, V)."""
    O, Ax, U, V = fr
    x, px, y, py, _, d = st
    ps = math.sqrt((1.0 + d) ** 2 - px * px - py * py)
    return O + x * U + y * V, (px * U + py * V + ps * Ax) / (1.0 + d)


_CACHE = {}


def electron_track(lattice, ws=None, correct=False):
    """Замкнутая орбита машины в глобальных координатах на входах/выходах всех магнитов.
    Возврат dict(inst2mag, pos_in, dir_in, pos_out, dir_out, st_in, st_out, rms_x, rms_y, corrected).
    rms_x/rms_y — по выходам магнитов в ИХ рамках: внутри перекошенного магнита сюда входит и смещение
    прямого пути электрона относительно наклонённых полюсов (это не искажение орбиты, коррекция его не снимает)."""
    key = (int(lattice), os.path.abspath(ws) if ws else None, bool(correct))
    if key in _CACHE:
        return _CACHE[key]
    import ring_analysis as ra
    import ring_cpp as rc
    with ws_tables(ws):
        py_ring = ra.Ring(lattice)                               # рамки элементов машины
        groups = _mag_instances(lattice)
    mags = py_ring.mags
    cpp = rc.CppRing(lattice=lattice, ws=ws, nsl=ORBIT_NSL)
    try:
        if cpp.n_mag != len(mags) or abs(cpp.circumference - py_ring.circumference) > 1e-6:
            raise RuntimeError("electron_track: Python-рамки и C++-ядро разошлись по геометрии машины")
        corrected = False
        if correct:
            import orbit_correct as oc
            ideal = rc.CppRing(lattice=lattice, nsl=ORBIT_NSL)
            try:
                resp = oc.cached_response_matrix(ideal, lattice)
            finally:
                ideal.close()
            res = oc.correct(cpp, resp)
            if res is None:
                raise RuntimeError("electron_track: коррекция орбиты не сошлась")
            co = np.asarray(res["co"], float); corrected = True
        else:
            co = ra.closed_orbit(cpp, dims=4)
        if not np.all(np.isfinite(co)):
            raise RuntimeError("electron_track: замкнутой орбиты нет (машина неустойчива — попробуйте --correct)")
        z0 = np.r_[co[:4], 0.0, 0.0]
        st_out = cpp.track_elements([z0])[0]                     # (n+1, 6): выходы магнитов, затем плоскость
    finally:
        cpp.close()
    n = len(mags)
    st_in = np.empty((n, 6))
    prev_st, prev = tuple(z0), py_ring.plane
    for j, m in enumerate(mags):
        st_in[j] = mt.modified_drift(tuple(map(float, prev_st)), prev, m)
        prev_st, prev = st_out[j], m
    pos_in = np.empty((n, 3)); dir_in = np.empty((n, 3)); pos_out = np.empty((n, 3)); dir_out = np.empty((n, 3))
    for j, m in enumerate(mags):
        pos_in[j], dir_in[j] = _glob(m.fr_in, st_in[j])
        pos_out[j], dir_out[j] = _glob(m.fr_out, st_out[j])
    inst2mag = {iid: j for j, grp in enumerate(groups) for iid in grp}
    out = dict(inst2mag=inst2mag, pos_in=pos_in, dir_in=dir_in, pos_out=pos_out, dir_out=dir_out,
               st_in=st_in, st_out=st_out[:n], corrected=corrected,
               rms_x=float(np.sqrt(np.mean(st_out[:n, 0] ** 2))), rms_y=float(np.sqrt(np.mean(st_out[:n, 2] ** 2))))
    _CACHE[key] = out
    return out


def device_source(et, insts):
    """Источник устройства (ондулятора) по орбите: среднее положение и среднее направление электрона по
    входам и выходам всех его магнитов (знакопеременный изгиб полюсов при этом усредняется). None, если
    магнитов устройства в машине нет."""
    js = sorted({et["inst2mag"][i] for i in insts if i in et["inst2mag"]})
    if not js:
        return None, None
    P = np.mean(np.r_[et["pos_in"][js], et["pos_out"][js]], axis=0)
    d = np.mean(np.r_[et["dir_in"][js], et["dir_out"][js]], axis=0)
    return P, d / np.linalg.norm(d)


def arc_state(et, iid, frac):
    """Состояние орбиты (x, x′, y, y′) в рамке дуги диполя iid на доле дуги frac∈[0,1] (0 — вход, 1 — выход):
    линейная интерполяция между входом и выходом магнита. None, если магнита нет."""
    j = et["inst2mag"].get(iid)
    if j is None:
        return None
    a = et["st_in"][j]; b = et["st_out"][j]
    f = min(1.0, max(0.0, float(frac)))

    def ang(s):
        ps = math.sqrt((1.0 + s[5]) ** 2 - s[1] ** 2 - s[3] ** 2)
        return s[1] / ps, s[3] / ps
    (xa, ya), (xb, yb) = ang(a), ang(b)
    return ((1 - f) * a[0] + f * b[0], (1 - f) * xa + f * xb, (1 - f) * a[2] + f * b[2], (1 - f) * ya + f * yb)


def fan_on_orbit(em, et):
    """Веер диполя (bend_screen.emitter_geom) по орбите: точка отбора смещается на (x, y) орбиты, касательная
    поворачивается на (x′, y′), плоскость веера наклоняется вместе с ней. Возврат нового dict (em не трогается)."""
    th = em["theta0"]
    frac = (em["alpha"] + th / 2) / th if th else 0.5            # 0 — вход дуги, 1 — выход
    s = arc_state(et, em["iid"], frac)
    if s is None:
        return dict(em)
    x, xp, y, yp = s
    U, V, t = em["U"], em["vn"], em["t"]
    d = t + xp * U + yp * V
    d = d / np.linalg.norm(d)
    vn = V - float(np.dot(V, d)) * d
    return dict(em, P=em["P"] + x * U + y * V, t=d, vn=vn / np.linalg.norm(vn))


def main():
    ap = argparse.ArgumentParser(description="источник света по орбите электрона vs по рамке железа")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--ws", default=None)
    ap.add_argument("--correct", action="store_true", help="сначала скорректировать орбиту (SVD)")
    a = ap.parse_args()
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from screen_analysis import load_groups, load_lattice, emitter_frame
    et = electron_track(a.lattice, a.ws, a.correct)
    ci, byid = load_lattice(a.lattice, a.ws)
    print(f"решётка {a.lattice} [{'номинал' if not a.ws else a.ws}]"
          f"{' + коррекция' if et['corrected'] else ''}: орбита в рамках магнитов rms x={et['rms_x']*1e6:.2f} мкм, "
          f"y={et['rms_y']*1e6:.2f} мкм")
    print(f"{'группа':>8s} | угол электрона к оси железа (гор., верт.), мкрад | смещение электрона от центра железа, мкм")
    for gid, insts in sorted(load_groups(a.lattice).items()):
        Pm, ax, u, v = emitter_frame(insts, ci, byid)
        Pe, de = device_source(et, insts)
        if Pm is None or Pe is None:
            continue
        print(f"{gid:>8s} | {float(np.dot(de, u))*1e6:+8.2f} {float(np.dot(de, v))*1e6:+8.2f} | "
              f"{float(np.dot(Pe - Pm, u))*1e6:+8.2f} {float(np.dot(Pe - Pm, v))*1e6:+8.2f}")


if __name__ == "__main__":
    main()
