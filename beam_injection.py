# -*- coding: utf-8 -*-
"""Конфигурация опыта для матричного ядра: ЯВНЫЕ старт/финиш-плоскости, пучок в координатах старта, резонаторы.

ОПЫТ (experiments.csv) задаёт две плоскости — точку и нормаль (направление пролёта):
  start_x..start_nz    старт: начало локальных координат пучка;
  finish_x..finish_nz  финиш; ПУСТО ⇒ финиш = старт (КОЛЬЦО).
Один «оборот» = пролёт старт → финиш. Кольцо — частица возвращается в ту же плоскость. Плоскости разные —
РАЗОМКНУТАЯ система с ПОРТАЛОМ: состояние на финише (в его рамке) становится состоянием на старте (в рамке
старта) следующего пролёта — тот же кортеж (x, px, y, py, t, δ).

Рамка плоскости (PlaneMagnet): ось = нормаль, u — горизонталь (n повёрнута на +90° в XY), v = n×u. Ядро
плоскость НЕ подправляет: modified_drift переводит состояние между любыми рамками точно. Если точка лежит на
опорной орбите и нормаль вдоль неё, локальные x, y — отклонения от опорной.

ПУЧОК (particle_groups.csv: group_id;part_id;m;q;x;xp;y;yp;ct;delta) — в рамке старт-плоскости, СИ:
x, y — смещения (м), xp = dx/ds, yp = dy/ds (как elegant), ct — начальное t (м), δ = (p − P0)/P0.
Каноника ядра: px = (1+δ)·xp/√(1+xp²+yp²), py аналогично.

РЕЗОНАТОРЫ — свойство УСТАНОВКИ (rf_cavities.csv: cavity_id;lattice_id;x;y;z;nx;ny;nz;volt_share;note):
плоскость нулевой длины в промежутке между элементами (пустой volt_share ⇒ доля 1/N). Встают в путь, только
если трекер включает ВЧ (rf=True). Нет резонаторов у решётки ⇒ ВЧ в конце пролёта (как всегда было).
Резонатор в плоскости финиша (у кольца — старта) — тоже в конце пролёта.

ИМПУЛЬСНЫЕ КИКЕРЫ (инжекция/выпуск; не путать со стационарными корректорами kind=kicker) — тонкий дипольный
удар, амплитуда которого зависит от времени прихода частицы. «Железо» — у УСТАНОВКИ (pulsed_kickers.csv:
kicker_id;lattice_id;x;y;z;nx;ny;nz;shape;duration_s;rise_s;note — плоскость, форма half_sine | trapezoid,
длительность основания, фронт = спад трапеции), работа — у ОПЫТА (kicker_timing.csv: exp_id;kicker_id;
theta_x_rad;theta_y_rad;delay_s;period_s;note — пиковый удар, задержка начала импульса, период; пусто — один
импульс). Кикер без строки у опыта в путь не встаёт. Удар (как elegant KICKER при L = 0): px += θx·f(τ),
py += θy·f(τ) в канонике, τ = T/c − delay; T — абсолютное время опыта в c·t (м) с момента, когда опорная
частица пересекла старт на пролёте 0 (начально T = ct частицы; копится Δt_raw/β(δ), портал мгновенный).
Кикер в плоскости финиша (у кольца — старта) бьёт в конце пролёта, в старте разомкнутой линии — в начале.

Положение плоскости в линии (locate_plane): промежуток между элементами (включая шов кольца), где опорная
прямая [выход mags[i] → вход mags[i+1]] пересекает плоскость по направлению нормали, ближайший к её точке;
иначе луч до первого / после последнего элемента (открытая линия). Плоскость внутри элемента — ошибка.
"""
import os, csv, math
import numpy as np
import matrix_tracker as mt
import bip_matrix as bm

ME, C = bm.ME, bm.C_LIGHT
P0 = bm.P_MOM                                    # опорный импульс γ₀·m·v₀
START_COLS = ("start_x", "start_y", "start_z", "start_nx", "start_ny", "start_nz")
FINISH_COLS = ("finish_x", "finish_y", "finish_z", "finish_nx", "finish_ny", "finish_nz")
LOCAL_COLS = ("x", "xp", "y", "yp", "ct", "delta")
CAVITY_COLS = ("cavity_id", "lattice_id", "x", "y", "z", "nx", "ny", "nz", "volt_share", "note")
KICKER_COLS = ("kicker_id", "lattice_id", "x", "y", "z", "nx", "ny", "nz", "shape", "duration_s", "rise_s", "note")
TIMING_COLS = ("exp_id", "kicker_id", "theta_x_rad", "theta_y_rad", "delay_s", "period_s", "note")
PULSE_SHAPES = ("half_sine", "trapezoid")
SEAM_TOL = 1e-3                                  # шов замкнут: отрезок выход→вход вдоль осей (рад) — классификация


def _gamma(v):
    return 1.0 / math.sqrt(1.0 - (v / C) ** 2)


def _rows(name):
    p = os.path.join(mt.TB, name)                # mt.TB читается при вызове — временные ws подменяют его
    if not os.path.exists(p):
        return []
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f, delimiter=";"))


def _num(s):
    s = (s or "").strip().replace(",", ".")
    return float(s) if s and s != "-" else None


def _vec(row, cols):
    """Плоскость из 6 столбцов → (точка, нормаль) или None, если все пусты. Частично заполненная — ошибка."""
    vals = [_num(row.get(c)) for c in cols]
    if all(v is None for v in vals):
        return None
    if any(v is None for v in vals):
        raise ValueError(f"плоскость {cols[0][:-2]}: заполнены не все столбцы {cols}")
    n = np.array(vals[3:], float)
    if not np.linalg.norm(n) > 0:
        raise ValueError(f"плоскость {cols[0][:-2]}: нулевая нормаль")
    return np.array(vals[:3], float), n


# ---------------- пучок ----------------
def load_group(group_id):
    """particle_groups.csv → список dict(local=(x,xp,y,yp,ct,delta), m, q, part_id) группы."""
    rows = _rows("particle_groups.csv")
    if rows and "xp" not in rows[0]:
        raise ValueError("particle_groups.csv в старом формате (глобальные x,y,z,vx..) — запустите "
                         "migrate_explicit_planes.py")
    out = []
    for r in rows:
        if r["group_id"] == str(group_id):
            loc = tuple((_num(r.get(c)) or 0.0) for c in LOCAL_COLS)
            out.append(dict(local=loc, m=float(r["m"]), q=float(r["q"]), part_id=r.get("part_id")))
    if not out:
        raise ValueError(f"группа {group_id} не найдена в particle_groups.csv")
    return out


def local_to_canonical(loc):
    """(x, xp, y, yp, ct, δ) в рамке плоскости → канонические (x, px, y, py, t, δ) ядра."""
    x, xp, y, yp, ct, d = (float(v) for v in loc[:6])
    k = (1.0 + d) / math.sqrt(1.0 + xp * xp + yp * yp)
    return (x, k * xp, y, k * yp, ct, d)


def canonical_to_local(st):
    """Обратное к local_to_canonical: (x, px, y, py, t, δ) → (x, xp, y, yp, ct, δ)."""
    x, px, y, py, t, d = (float(v) for v in st[:6])
    ps = math.sqrt((1.0 + d) ** 2 - px * px - py * py)
    return (x, px / ps, y, py / ps, t, d)


def local_to_global(plane, loc, m=ME):
    """Локальная частица → глобальные (pos, vel) — как диффурное ядро (src/exp_config.h): |p| = (1+δ)·P0."""
    x, xp, y, yp, _ct, d = (float(v) for v in loc[:6])
    pos = plane.center + x * plane.uhat + y * plane.vhat
    dvec = xp * plane.uhat + yp * plane.vhat + plane.axis
    dvec = dvec / np.linalg.norm(dvec)
    gb = (1.0 + d) * P0 / (m * C)                # γβ
    return pos, C * gb / math.sqrt(1.0 + gb * gb) * dvec


# ---------------- псевдо-магниты «плоскость» и «резонатор» ----------------
class PlaneMagnet:
    """Нулевой элемент-плоскость: рамка (axis=нормаль, uhat=гориз., vhat=верт.), entry=exit=center."""
    is_cavity = False
    is_kicker = False

    def __init__(self, center, axis):
        a = np.array(axis, float); a = a / np.linalg.norm(a)
        u = np.array([-a[1], a[0], 0.0])
        nu = np.linalg.norm(u); u = u / nu if nu > 1e-9 else np.array([0.0, 1.0, 0.0])
        self.center = np.array(center, float); self.axis = a; self.uhat = u
        self.vhat = np.cross(a, u)               # вертикаль (для планарного кольца ≈ +z)
        self.entry = self.center; self.exit = self.center
        self.L = self.Leff = 0.0; self.dev = "INJ_PLANE"; self.devs = ["INJ_PLANE"]
        self.yaw = math.atan2(a[1], a[0]); self.nsl = 1
        self.theta0 = self.k1 = self.k2 = self.k3 = 0.0; self.nonzero = False
        self.sector = False                                   # плоскость — не сектор
        self.fr_in = (self.entry, self.axis, self.uhat, self.vhat)   # рамки стыковки (⊥ оси)
        self.fr_out = (self.exit, self.axis, self.uhat, self.vhat)

    def body(self, st):                          # плоскость длины 0 — состояние не меняется
        return st


class CavityNode(PlaneMagnet):
    """ВЧ-резонатор — плоскость нулевой длины на пути. Кик делает трекер (body — тождество, как у плоскости).
    s_c — опорная дуга от старт-плоскости, share — доля напряжения."""
    is_cavity = True

    def __init__(self, center, axis, cavity_id, share):
        super().__init__(center, axis)
        self.cavity_id = str(cavity_id); self.share = float(share); self.s_c = 0.0
        self.dev = f"RF_CAVITY_{cavity_id}"; self.devs = [self.dev]


def pulse_shape(shape, duration, rise, tau):
    """Форма импульса f(τ) ∈ [0, 1]; τ — время от начала импульса (с). half_sine — sin(π·τ/D) на [0, D];
    trapezoid — линейный фронт rise, плато, линейный спад rise (rise = 0 — прямоугольник). Вне [0, D] — 0."""
    if tau < 0.0 or tau > duration:
        return 0.0
    if shape == "half_sine":
        return math.sin(math.pi * tau / duration)
    if rise > 0.0:
        if tau < rise:
            return tau / rise
        if tau > duration - rise:
            return (duration - tau) / rise
    return 1.0


class KickerNode(PlaneMagnet):
    """Импульсный кикер — плоскость нулевой длины на пути. Удар делает трекер по абсолютному времени T (c·t, м):
    τ = T/c − delay (с периодом — по модулю периода), px += θx·f(τ), py += θy·f(τ). s_k — опорная дуга от старта."""
    is_kicker = True

    def __init__(self, center, axis, spec):
        super().__init__(center, axis)
        self.kicker_id = str(spec["kicker_id"]); self.shape = spec["shape"]
        self.duration = float(spec["duration"]); self.rise = float(spec["rise"])
        self.theta_x = float(spec["theta_x"]); self.theta_y = float(spec["theta_y"])
        self.delay = float(spec["delay"]); self.period = spec["period"]
        self.s_k = 0.0
        self.dev = f"PULSE_KICKER_{self.kicker_id}"; self.devs = [self.dev]

    def amplitude(self, T):
        tau = T / C - self.delay
        if self.period is not None:
            if tau < 0.0:
                return 0.0
            tau = math.fmod(tau, self.period)
        return pulse_shape(self.shape, self.duration, self.rise, tau)

    def kick(self, st, T):
        f = self.amplitude(T)
        if f == 0.0:
            return st
        x, px, y, py, t, d = st
        return (x, px + self.theta_x * f, y, py + self.theta_y * f, t, d)


# ---------------- опыт и установка ----------------
def experiment_row(exp_id):
    for r in _rows("experiments.csv"):
        if r.get("exp_id") == str(exp_id):
            return r
    raise ValueError(f"опыт {exp_id} не найден в experiments.csv")


def planes_of(row):
    """Строка experiments.csv → (start, finish): (точка, нормаль); finish None ⇒ кольцо."""
    if "start_x" not in row:
        raise ValueError("experiments.csv в старом формате (p1..p3/dir) — запустите migrate_explicit_planes.py")
    start = _vec(row, START_COLS)
    if start is None:
        raise ValueError(f"опыт {row.get('exp_id')}: не задана старт-плоскость (start_x..start_nz) — "
                         f"py -3.14 exp_planes.py --exp {row.get('exp_id')} --start-s 0")
    return start, _vec(row, FINISH_COLS)


def find_experiment(lattice, group_id):
    """Единственный опыт с этой парой (решётка, группа) или ValueError (нет / неоднозначно — передайте exp_id)."""
    hits = [r["exp_id"] for r in _rows("experiments.csv")
            if r.get("lattice_id") == str(lattice) and (group_id is None or r.get("group_id") == str(group_id))]
    if len(hits) != 1:
        raise ValueError(f"{'нет опыта' if not hits else 'неоднозначно: опыты ' + ', '.join(hits)} с решёткой "
                         f"{lattice} и группой {group_id} — передайте exp_id")
    return int(hits[0])


def load_cavities(lattice):
    """Резонаторы установки → [dict(cavity_id, center, normal, share)] (пустой volt_share ⇒ 1/N)."""
    mine = [r for r in _rows("rf_cavities.csv") if r.get("lattice_id") == str(lattice)]
    out = []
    for r in mine:
        pl = _vec(r, ("x", "y", "z", "nx", "ny", "nz"))
        if pl is None:
            raise ValueError(f"rf_cavities: у резонатора {r.get('cavity_id')} не задана плоскость")
        sh = _num(r.get("volt_share"))
        out.append(dict(cavity_id=r.get("cavity_id"), center=pl[0], normal=pl[1],
                        share=(1.0 / len(mine)) if sh is None else sh))
    return out


def load_kickers(lattice, exp_id):
    """Импульсы опыта: строки kicker_timing.csv[exp_id] + «железо» pulsed_kickers.csv[lattice] → [dict(kicker_id,
    center, normal, shape, duration, rise, theta_x, theta_y, delay, period)]. Несколько строк одного кикера —
    несколько импульсов (удары складываются). Нет строк у опыта ⇒ []."""
    timing = [r for r in _rows("kicker_timing.csv") if r.get("exp_id") == str(exp_id)]
    if not timing:
        return []
    hw = {r.get("kicker_id"): r for r in _rows("pulsed_kickers.csv") if r.get("lattice_id") == str(lattice)}
    out = []
    for r in timing:
        kid = r.get("kicker_id")
        h = hw.get(kid)
        if h is None:
            raise ValueError(f"kicker_timing: опыт {exp_id} включает кикер {kid}, но у решётки {lattice} его нет "
                             f"в pulsed_kickers.csv")
        pl = _vec(h, ("x", "y", "z", "nx", "ny", "nz"))
        if pl is None:
            raise ValueError(f"pulsed_kickers: у кикера {kid} не задана плоскость")
        shape = (h.get("shape") or "half_sine").strip().lower()
        if shape not in PULSE_SHAPES:
            raise ValueError(f"pulsed_kickers: кикер {kid} — неизвестная форма {shape!r} (допустимо: "
                             f"{', '.join(PULSE_SHAPES)})")
        D = _num(h.get("duration_s"))
        if D is None or not D > 0:
            raise ValueError(f"pulsed_kickers: кикер {kid} — длительность duration_s должна быть > 0")
        rise = _num(h.get("rise_s")) or 0.0
        if shape == "trapezoid" and (rise < 0 or 2.0 * rise > D):
            raise ValueError(f"pulsed_kickers: кикер {kid} — фронт rise_s должен быть в [0; duration_s/2]")
        period = _num(r.get("period_s"))
        if period is not None and not period >= D:
            raise ValueError(f"kicker_timing: опыт {exp_id}, кикер {kid} — период {period} с короче импульса {D} с")
        out.append(dict(kicker_id=kid, center=pl[0], normal=pl[1], shape=shape, duration=D, rise=rise,
                        theta_x=_num(r.get("theta_x_rad")) or 0.0, theta_y=_num(r.get("theta_y_rad")) or 0.0,
                        delay=_num(r.get("delay_s")) or 0.0, period=period))
    return out


# ---------------- положение плоскости в линии и путь старт → финиш ----------------
def locate_plane(mags, center, normal):
    """(slot, along): slot k — промежуток ПОСЛЕ mags[k] (k = N−1 — шов кольца или луч после последнего
    элемента), slot −1 — луч ДО первого элемента; along — путь от выхода mags[k] (у −1 — от входа mags[0],
    ≤ 0) до точки пересечения опорной прямой с плоскостью."""
    n = np.asarray(normal, float); n = n / np.linalg.norm(n); c = np.asarray(center, float)
    N = len(mags); best = None
    for i in range(N):                                      # промежутки (включая шов) — как прежний поиск
        A = mags[i].exit; B = mags[(i + 1) % N].entry
        AB = B - A; L = float(np.linalg.norm(AB))
        if L < 1e-9:
            continue
        axis = AB / L
        denom = float(n @ AB)
        if abs(denom) < 1e-12 or float(n @ axis) <= 0:
            continue
        t = float(n @ (c - A)) / denom
        if -1e-9 <= t <= 1.0 + 1e-9:                         # точка пересечения — внутри промежутка
            dist = float(np.linalg.norm(A + t * AB - c))
            if best is None or dist < best[2]:
                best = (i, t * L, dist)
    if best is None:                                        # открытая линия: лучи до первого / после последнего
        for slot, P, ax, sign in ((-1, mags[0].entry, mags[0].axis, -1), (N - 1, mags[-1].exit, mags[-1].axis, +1)):
            d = float(n @ ax)
            if d <= 1e-12:
                continue
            s = float(n @ (c - P)) / d
            if sign * s >= 0:
                dist = float(np.linalg.norm(P + s * ax - c))
                if best is None or dist < best[2]:
                    best = (slot, s, dist)
    if best is None:
        raise ValueError(f"плоскость {np.round(c, 6).tolist()} n={np.round(n, 6).tolist()} не пересекает опорную "
                         f"орбиту по направлению нормали ни в одном промежутке (внутри элемента или против хода?)")
    return best[0], best[1]


def seam_is_drift(mags):
    """Линия замкнута в кольцо: шов (выход последнего → вход первого) — прямой промежуток вдоль их осей."""
    A = mags[-1].exit; B = mags[0].entry; AB = B - A; L = float(np.linalg.norm(AB))
    if L < 1e-9:
        return True
    ax = AB / L
    return (float(ax @ mags[-1].axis) > 0 and float(np.linalg.norm(np.cross(ax, mags[-1].axis))) < SEAM_TOL
            and float(np.linalg.norm(np.cross(ax, mags[0].axis))) < SEAM_TOL)


def same_plane(p, q, settings=None):
    """Две плоскости (точка, нормаль) совпадают — по допускам опыта (как стыки элементов)."""
    st = settings or mt.Settings()
    a = PlaneMagnet(*p); b = PlaneMagnet(*q)
    return mt._plane_relation(a, b, st)[0] == "same" and float(a.axis @ b.axis) > 0


def _gap(prev, node):
    return float((node.fr_in[0] - prev.fr_out[0]) @ prev.fr_out[1])


def build_path(mags, start, finish=None, cavities=(), settings=None, kickers=()):
    """Путь опыта: dict(plane, ordered, finish, closed, line_len, cavities, end_share, skipped, kickers,
    kick_begin, kick_end, kick_skipped, after, before).
    ordered — элементы от старта до финиша (с переходом через шов кольца) + узлы-резонаторы и узлы-кикеры
    (load_kickers) в своих промежутках; end_share — доля резонаторов в плоскости финиша (кик в конце пролёта),
    None — таких нет; kick_end / kick_begin — кикеры в плоскости финиша (конец пролёта) / старта разомкнутой
    линии (начало пролёта)."""
    st = settings or mt.Settings()
    N = len(mags)
    plane = PlaneMagnet(*start)
    ks, a_s = locate_plane(mags, *start)
    closed = finish is None or same_plane(start, finish, st)
    fin = plane if closed else PlaneMagnet(*finish)
    if closed:
        if not seam_is_drift(mags):
            raise ValueError("финиш не задан (кольцо), но линия решётки не замкнута — задайте финиш-плоскость")
        order = [(ks + 1 + i) % N for i in range(N)]
        kf, a_f = ks, a_s
    else:
        kf, a_f = locate_plane(mags, *finish)
        if kf > ks or (kf == ks and a_f > a_s):
            order = list(range(ks + 1, kf + 1))
        else:
            if not seam_is_drift(mags):
                raise ValueError("финиш-плоскость стоит раньше старта, а линия не замкнута — путь невозможен")
            order = list(range(ks + 1, N)) + list(range(0, kf + 1))
    elements = [mags[i] for i in order]
    slot_of_gap = [ks] + [order[g] for g in range(len(order))]          # промежуток g — перед элементом g
    fin_plane = start if closed else finish

    def place(cpl):
        """Где плоскость события на пути: 'end' (финиш), 'begin' (старт разомкнутой линии), (g, along) —
        промежуток g (перед элементом g), None — вне пути."""
        if same_plane(cpl, fin_plane, st):
            return "end"
        if not closed and same_plane(cpl, start, st):
            return "begin"
        kc, a_c = locate_plane(mags, *cpl)
        for g, slot in enumerate(slot_of_gap):
            if slot != kc:
                continue
            lo = a_s if g == 0 else -math.inf
            hi = a_f if g == len(elements) else math.inf
            if lo < a_c < hi:
                return (g, a_c)
        return None

    ins = {}; end_share = None; skipped = []
    for cv in cavities:
        where = place((cv["center"], cv["normal"]))
        if where == "end":
            end_share = (end_share or 0.0) + cv["share"]; continue
        if where is None or where == "begin":                # на старте разомкнутой линии = на финише прошлого
            skipped.append(cv["cavity_id"]); continue
        ins.setdefault(where[0], []).append(
            (where[1], CavityNode(cv["center"], cv["normal"], cv["cavity_id"], cv["share"])))
    kick_begin = []; kick_end = []; kick_skipped = []
    for kk in kickers:
        where = place((kk["center"], kk["normal"]))
        node = KickerNode(kk["center"], kk["normal"], kk)
        if where == "end":
            kick_end.append(node); continue
        if where == "begin":
            kick_begin.append(node); continue
        if where is None:
            kick_skipped.append(kk["kicker_id"]); continue
        ins.setdefault(where[0], []).append((where[1], node))
    ordered = []
    for g in range(len(elements) + 1):
        ordered += [node for _a, node in sorted(ins.get(g, []), key=lambda t: t[0])]
        if g < len(elements):
            ordered.append(elements[g])
    s = 0.0; prev = plane; cav_on = []; kick_on = []
    for node in ordered:                                     # опорная дуга — как C++ entry_arc / ring_len
        s = s + _gap(prev, node)
        if node.is_cavity:
            node.s_c = s; cav_on.append(node)
        elif node.is_kicker:
            node.s_k = s; kick_on.append(node)
        s = s + node.L; prev = node
    line_len = s + _gap(prev, fin)
    for node in kick_end:
        node.s_k = line_len
    after = mags[ks].dev if ks >= 0 else "НАЧАЛО ЛИНИИ"
    before = mags[(ks + 1) % N].dev if (ks + 1 < N or closed) else "КОНЕЦ ЛИНИИ"
    return dict(plane=plane, ordered=ordered, finish=fin, closed=closed, line_len=line_len, cavities=cav_on,
                end_share=end_share, skipped=skipped, kickers=kick_on, kick_begin=kick_begin, kick_end=kick_end,
                kick_skipped=kick_skipped, after=after, before=before)


def is_node(m):
    """Узел-событие пути (резонатор или импульсный кикер), а не элемент решётки."""
    return getattr(m, "is_cavity", False) or getattr(m, "is_kicker", False)


# ---------------- сборка опыта ----------------
def inject(lattice, group_id=None, source="physical", settings=None, exp_id=None, start=None, finish=None, rf=False,
           kick=False):
    """Опыт → (plane, ordered, states, info). Плоскости — из experiments.csv[exp_id] (exp_id не задан ⇒
    единственный опыт с парой (lattice, group_id)) либо явные start/finish ((точка, нормаль); finish None ⇒
    кольцо). rf=True ⇒ резонаторы установки встают в путь (info['cavities'], info['end_share']); у решётки без
    резонаторов — ВЧ в конце пролёта (end_share = 1). kick=True ⇒ импульсные кикеры опыта (kicker_timing.csv)
    встают в путь (info['kickers'], info['kick_begin'], info['kick_end']); явные плоскости без опыта — без
    кикеров. settings — matrix_tracker.Settings прогона."""
    if start is None:
        if exp_id is None:
            exp_id = find_experiment(lattice, group_id)
        row = experiment_row(exp_id)
        if int(row["lattice_id"]) != int(lattice):
            raise ValueError(f"опыт {exp_id} — решётка {row['lattice_id']}, запрошена {lattice}")
        start, finish = planes_of(row)
        if group_id is None:
            group_id = int(row["group_id"])
    if group_id is None:
        raise ValueError("не задана группа частиц")
    eff, phys, coup, trans = mt.load_db()
    mags = mt.build_line(lattice, eff, phys, coup, trans, source=source, settings=settings)
    cav = []
    if rf:
        cav = load_cavities(lattice)
    kk = load_kickers(lattice, exp_id) if (kick and exp_id is not None) else []
    path = build_path(mags, start, finish, cav, settings, kk)
    if not any(not is_node(m) for m in path["ordered"]):
        raise ValueError(f"на пути старт → финиш опыта {exp_id} нет ни одного элемента решётки {lattice}")
    if rf:
        if not cav:
            path["end_share"] = 1.0                           # прежнее правило: ВЧ в конце пролёта
        elif not path["cavities"] and not path["end_share"]:
            raise ValueError(f"ни один резонатор решётки {lattice} (rf_cavities.csv) не стоит на пути старт → "
                             f"финиш опыта {exp_id}")
        for cid in path["skipped"]:
            print(f"  [ВЧ] резонатор {cid} решётки {lattice} вне пути старт → финиш — пропущен")
    for kid in path["kick_skipped"]:
        print(f"  [кикер] импульсный кикер {kid} решётки {lattice} вне пути старт → финиш — пропущен")
    parts = load_group(group_id)
    states = [local_to_canonical(p["local"]) for p in parts]
    plane = path["plane"]
    info = dict(nparts=len(parts), after_elem=path["after"], before_elem=path["before"],
                center=plane.center, axis=plane.axis, closed=path["closed"], finish=path["finish"],
                line_len=path["line_len"], cavities=path["cavities"], end_share=path["end_share"],
                kickers=path["kickers"], kick_begin=path["kick_begin"], kick_end=path["kick_end"],
                exp_id=exp_id, group=group_id, masses=[p["m"] for p in parts])
    return plane, path["ordered"], states, info


# ---------------- трек одного пролёта (старт → путь → финиш) ----------------
def track_turn(plane, ordered, states, finish=None):
    """Каждая частица: старт-плоскость → элементы → ФИНИШ-плоскость. Возврат [(final_state, lost)].
    finish=None ⇒ финиш совпадает со стартом (замкнутый полный оборот)."""
    if finish is None:
        finish = plane
    res = []
    for st in states:
        s = tuple(st); lost = None
        try:
            prev = plane
            for m in ordered:
                s = mt.modified_drift(s, prev, m); s = m.body(s); prev = m
            s = mt.modified_drift(s, prev, finish)               # последний → ФИНИШ (замкнут: = старт)
        except mt.ParticleLost as e:
            lost = str(e)
        res.append((s, lost))
    return res


if __name__ == "__main__":
    import sys
    exp = int(sys.argv[1]) if len(sys.argv) > 1 else 1
    row = experiment_row(exp)
    lat, gid = int(row["lattice_id"]), int(row["group_id"])
    plane, ordered, states, info = inject(lat, gid, exp_id=exp, rf=True, kick=True)
    print(f"Опыт {exp}: решётка {lat}, группа {gid}, {info['nparts']} частиц")
    print(f"  старт между '{info['after_elem']}' и '{info['before_elem']}': точка {info['center'].round(6)}, "
          f"нормаль {info['axis'].round(6)}")
    kind = "КОЛЬЦО (финиш = старт)" if info['closed'] else "РАЗОМКНУТ (портал финиш → старт)"
    print(f"  {kind}; элементов на пути {sum(1 for m in ordered if not is_node(m))}, "
          f"длина пролёта {info['line_len']:.6f} м")
    for c in info["cavities"]:
        print(f"  резонатор {c.cavity_id}: s = {c.s_c:.6f} м от старта, доля {c.share:g}")
    if info["end_share"]:
        print(f"  резонатор в конце пролёта: доля {info['end_share']:g}")
    for k in info["kickers"] + info["kick_begin"] + info["kick_end"]:
        print(f"  импульсный кикер {k.kicker_id}: s = {k.s_k:.6f} м от старта, θ = ({k.theta_x:g}, {k.theta_y:g}) "
              f"рад, {k.shape} {k.duration:g} с, задержка {k.delay:g} с")
    # track_turn — статическая машина (без ударов кикеров и ВЧ): только геометрия пролёта
    surv = 0
    for k, (st, lost) in enumerate(track_turn(plane, ordered, states, info['finish'])):
        if lost:
            print(f"  частица {k}: ПОТЕРЯНА — {lost}")
        else:
            print(f"  частица {k}: финиш x={st[0]*1e3:+.4f}мм y={st[2]*1e3:+.4f}мм t={st[4]:.4f}м"); surv += 1
    print(f"ВЫЖИЛО {surv}/{info['nparts']}")
