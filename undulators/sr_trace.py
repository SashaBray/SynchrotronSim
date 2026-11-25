# -*- coding: utf-8 -*-
"""ТОЧНЫЙ режим тепловой карты экранов: прямое интегрирование формулы Льенара по траектории электрона.

Быстрый режим (screen_analysis.footprint, bend_screen.footprint) — «точечный эмиттер»: ондулятор светит из одной
точки табличной угловой мощностью (Льенар, усреднённый по периоду бесконечной синусоиды), диполь — веером из одной
точки отбора с жёстким краем, пучок — один электрон. Здесь вместо этого траектория электрона разбивается на
короткие участки, и каждый участок светит на каждый пиксель экрана из СВОЕГО места, СВОИМ направлением:

    dP/dA = (I/e)·Σ_k Δt_k · e²/(16π²ε₀c) · [|β̇|²/q³ − (n·β̇)²/(γ²q⁵)] · |n·N_экр| / R_k²,   q = 1 − n·β.

Это точно для мощности (наблюдатель ↔ момент излучения взаимно однозначны, перекрёстных членов в ∫|E|²dt нет —
тот же Парсеваль, что запрещает интерференцию в тепловой карте), в том числе в ближней зоне: n и R берутся для
каждой точки отдельно. Отброшено только поле скорости (~ρ/(γ²R) ~ 3·10⁻⁸). q = (1−β) + β·|n−t|²/2 — без
вычитания близких единиц (t — единичная касательная).

ЭФФЕКТЫ (каждый можно выключить, см. parse_effects):
  extended — протяжённый источник: каждая точка дуги/ондулятора светит из своего места (выкл.: все точки
             траектории стянуты в одну — центр ИД / точку отбора дуги, как в быстром режиме);
  poles    — реальная структура ИД из БД: полупериоды синусоиды с амплитудами полюсов (основные ±1, концевые
             UE0/UE1 ≈ 0,49), электрон входит по орбите без начального угла (выкл.: N идеальных полных периодов,
             x′ = (K/γ)·sin — как таблица быстрого режима);
  fringe   — мягкий край диполей: продольный профиль поля по карте устройства (hardware_library.map_bin_path;
             полукарта зеркалится), модули соседних диполей складываются (суперпозиция), профиль нормирован
             на ∫f = L_eff — угол изгиба как в решётке (выкл.: жёсткий край, как в трекинге);
  beam     — пучок: гауссова свёртка карты каждого источника с размером пучка на экране
             σx²(D) = εx·(β−2αD+γD²) + σδ²·(η+η′D)²,  σy²(D) = κ·εx·(β−2αD+γD²)  (Твисс в месте источника).
Не входит (пренебрежимо или отдельная задача): энергоразброс в самом излучении (~σδ²), СИ квадруполей и
корректоров от смещённой орбиты, поглощение/отражение по спектру и материалу, рассеяние фотонов, затенение.
Диполи без карты и DQ — жёсткий край.

ЧИСЛЕННОЕ: шаг ИД — ppp точек на период; шаг дуги — поворот касательной ≤ 1/(8γ); поле в шаге — среднее по
шагу (угол изгиба сохраняется точно и при жёстком крае). Отсев: точка траектории, чья касательная не достаёт
экран с запасом cull/γ, отбрасывается (cull=40 → хвост Льенара < 10⁻⁸ пика; cull=0 — без отсева).

ИСПОЛЬЗОВАНИЕ: ctx = make_context(lattice, ws, et, gap); fp = footprint_accurate(scr, ctx, effects=…).
CLI — `screen_analysis.py --calc accurate`, `bend_screen.py --calc accurate`, `screen_report.py --calc accurate`.
"""
import os, sys, csv, math, time, functools
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.abspath(os.path.join(_HERE, ".."))
for _p in (_ROOT, _HERE):
    if _p not in sys.path:
        sys.path.insert(0, _p)
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt
import undulator_sr as sr

_E_CH = 1.602176634e-19; _EPS0 = 8.8541878128e-12; _C = 299792458.0
_LW = _E_CH ** 2 / (16 * math.pi ** 2 * _EPS0 * _C)          # e²/(16π²ε₀c), Дж·с
EFFECTS = ("extended", "poles", "fringe", "beam")
FIELDS_BIN = os.path.join("workspace", "fields_bin")
HWL = os.path.join("workspace", "tables", "hardware_library.csv")
_UP = np.array([0.0, 0.0, 1.0])


def parse_effects(spec=None):
    """Набор эффектов: None/'all' — все; 'all,-beam' или '-fringe' — все, кроме; 'none' — ни одного;
    'extended,poles' — только перечисленные."""
    if spec is None:
        return frozenset(EFFECTS)
    if isinstance(spec, (set, frozenset, list, tuple)):
        spec = ",".join(spec) if spec else "none"
    items = [x.strip() for x in str(spec).split(",") if x.strip()]
    on = set(EFFECTS) if (not items or items[0] == "all" or all(x.startswith("-") for x in items)) else set()
    for x in items:
        if x == "all":
            on |= set(EFFECTS)
            continue
        if x == "none":
            on.clear()
            continue
        name = x.lstrip("+-")
        if name not in EFFECTS:
            raise ValueError(f"неизвестный эффект {name!r}; есть: {', '.join(EFFECTS)}")
        (on.discard if x.startswith("-") else on.add)(name)
    return frozenset(on)


def _unit(v):
    v = np.asarray(v, float); n = float(np.linalg.norm(v))
    return v / n if n else v


def _rotf(v, n, a):
    """Поворот вектора-кортежа v вокруг единичной оси n на угол a (Родриг, на float — для шагов интегратора)."""
    c, s = math.cos(a), math.sin(a); d = (n[0] * v[0] + n[1] * v[1] + n[2] * v[2]) * (1.0 - c)
    return (v[0] * c + (n[1] * v[2] - n[2] * v[1]) * s + n[0] * d,
            v[1] * c + (n[2] * v[0] - n[0] * v[2]) * s + n[1] * d,
            v[2] * c + (n[0] * v[1] - n[1] * v[0]) * s + n[2] * d)


# ─────────────────────────── траектория ондулятора ───────────────────────────

def wiggle_trace(O, t_in, u, sj, a, kap0, lam, ppp, g, centered=False):
    """Траектория в ондуляторе: прямая O + s·t_in плюс качание в плоскости u от поля B0·Σa_j·cos(k(s−s_j)) на
    полупериодах |s−s_j| ≤ λ/4 (kap0 = B0/(Bρ) — кривизна в пике). Качание интегрируется аналитически по лепесткам
    (x′ и x — от нуля на входе). centered=True — вычесть средний угол (идеальная синусоида x′=(K/γ)·sin).
    Возврат dict(r, t, acc=β̇ [1/с], dt, …) по точкам с полем (середины шагов)."""
    O = np.asarray(O, float); t_in = _unit(t_in)
    up = _unit(np.asarray(u, float) - float(np.dot(u, t_in)) * t_in)
    sj = np.asarray(sj, float); a = np.asarray(a, float)
    k = 2 * math.pi / lam; q4 = lam / 4
    s0, s1 = float(sj.min()) - q4, float(sj.max()) + q4
    nseg = max(1, int(math.ceil((s1 - s0) / (lam / ppp))))
    ds = (s1 - s0) / nseg
    s = s0 + (np.arange(nseg) + 0.5) * ds
    x = s[None, :] - sj[:, None]                                   # (лепестки, точки)
    phi = np.clip(k * x, -math.pi / 2, math.pi / 2)
    inside = np.abs(x) <= q4
    kap = kap0 * (a[:, None] * np.cos(k * x) * inside).sum(0)
    hp = (kap0 / k) * (a[:, None] * (np.sin(phi) + 1.0)).sum(0)
    h = (kap0 / k) * (a[:, None] * (-np.cos(phi) / k + np.clip(np.minimum(x, q4) + q4, 0.0, None)
                                    + 2.0 * np.clip(x - q4, 0.0, None))).sum(0)
    if centered:
        m = float(np.mean(hp))
        hp = hp - m; h = h - m * (s - s0)
    tv = t_in[None, :] + hp[:, None] * up[None, :]
    nr = np.linalg.norm(tv, axis=1)
    t = tv / nr[:, None]
    kv = kap[:, None] * (up[None, :] - (t @ up)[:, None] * t) / nr[:, None]   # dt/ds
    r = O[None, :] + s[:, None] * t_in[None, :] + h[:, None] * up[None, :]
    beta = math.sqrt(1.0 - 1.0 / g ** 2)
    keep = kap != 0.0
    sc = float(np.mean(sj))
    return dict(r=r[keep], t=t[keep], acc=beta ** 2 * _C * kv[keep], dt=np.full(int(keep.sum()), ds / (beta * _C)),
                center=O + sc * t_in, dir=t_in, ex=up, ey=np.cross(t_in, up), ds=ds)


def undulator_trace(gid, insts, ci, byid, gap, et=None, poles=True, ppp=128, p=sr.P):
    """Траектория электрона в ондуляторе группы gid. Полюса — строки kind∈{dipole,dq} группы; их центры
    проецируются на линию электрона. С et (орбита): вход — состояние замкнутой орбиты на входе первого полюса;
    без et — прямая по оси железа через центр устройства. Амплитуда — B0(зазор), как в быстром режиме;
    poles=True — относительные силы и знаки полюсов из БД (концевые ≈ 0,49), иначе N идеальных периодов."""
    import screen_analysis as sa
    rows = [(i, byid[i]) for i in insts if i in byid and byid[i][ci["kind"]].strip() in sa.RADIATE]
    if not rows:
        return None
    Pm, axis, u, _ = sa.emitter_frame(insts, ci, byid)
    C = np.array([[float(r[ci[c]] or 0) for c in ("x", "y", "z")] for _, r in rows])
    th = np.array([mt.QP * float(r[ci["mat_B"]] or 0) * float(r[ci["eff_length"]] or 0) for _, r in rows])
    first = int(np.argmin((C - Pm) @ axis))
    O, t_in = Pm, axis
    if et is not None and rows[first][0] in et["inst2mag"]:
        j = et["inst2mag"][rows[first][0]]
        O, t_in = np.asarray(et["pos_in"][j], float), _unit(et["dir_in"][j])
    s_all = (C - O) @ t_in
    lam = p.lam_u; g = sr.gamma(p)
    kap0 = sr.B0_of_gap(gap, p) / (p.E / 0.299792458)
    if poles:
        order = np.argsort(s_all)
        sj = s_all[order]; a = th[order] / float(np.max(np.abs(th)))
        tr = wiggle_trace(O, t_in, u, sj, a, kap0, lam, ppp, g, centered=False)
    else:
        n = 2 * int(p.N); sc = float(np.mean(s_all))
        sj = sc + (np.arange(n) - (n - 1) / 2.0) * lam / 2.0
        a = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
        tr = wiggle_trace(O, t_in, u, sj, a, kap0, lam, ppp, g, centered=True)
    tr.update(kind="und", name=gid, n_poles=len(sj), sum_a2=float(np.sum(a ** 2)))
    return tr


# ─────────────────────────── диполи: модули, цепочки, мягкий край ───────────────────────────

_DEV_MAPS = {}                 # тип устройства → путь карты (или None) из документов машин (register_machine_maps)


@functools.lru_cache(maxsize=None)
def _map_paths():
    out = {}
    if os.path.exists(HWL):
        with open(HWL, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f, delimiter=";"):
                pth = (r.get("map_bin_path") or "").strip()
                if pth and pth != "-":
                    out[r["device_type_id"]] = os.path.join(FIELDS_BIN, pth)
    return out


@functools.lru_cache(maxsize=None)
def map_profile(device_type):
    """Нормированный продольный профиль поля диполя по его карте: f(0)=1 в центре, ось — продольная X карты,
    точка y=z=0, компонента с наибольшим |B|. Полукарта (minX = 0) зеркалится. Возврат dict(x, f, cum, half,
    integral) или None (нет карты)."""
    path = _DEV_MAPS[device_type] if device_type in _DEV_MAPS else _map_paths().get(device_type)
    if not path or not os.path.exists(path):
        return None
    import scan_field_maps as sfm
    h, F = sfm.read_field(path)
    iy = int(round((0.0 - h["minY"]) / h["stepY"])); iz = int(round((0.0 - h["minZ"]) / h["stepZ"]))
    B = F[:, iy, iz, :]
    x = h["minX"] + np.arange(h["nx"]) * h["stepX"]
    c = int(np.argmax(np.abs(B).max(axis=0)))
    b = B[:, c]
    if abs(h["minX"]) < 1e-12:                                     # полукарта от центра магнита
        x = np.concatenate([-x[:0:-1], x]); b = np.concatenate([b[:0:-1], b])
    else:                                                          # полная карта: центр — середина
        x = x - 0.5 * (h["minX"] + h["maxX"])
    f = b / b[np.argmin(np.abs(x))]
    cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])
    return dict(x=x, f=f, cum=cum, half=float(x[-1]), integral=float(cum[-1]), path=path)


def bend_modules(ci, byid, E, exclude=()):
    """Изгибающие модули машины (dipole/dq с θ≠0, без полюсов ИД): геометрия дуги + точки и касательные входа/выхода."""
    import bend_screen as bs
    from build_screens import _rot
    mods = []
    for iid, r in byid.items():
        if iid in exclude:
            continue
        d = bs.dipole_arc(r, ci, E)
        if d is None:
            continue
        O = d["C"] + d["rho"] * d["uhat"]; th = d["theta0"]
        d.update(kind=r[ci["kind"]].strip(), O=O,
                 P_in=O + _rot(d["C"] - O, d["vhat"], -th / 2), P_out=O + _rot(d["C"] - O, d["vhat"], th / 2),
                 t_in=_unit(_rot(d["axis"], d["vhat"], -th / 2)), t_out=_unit(_rot(d["axis"], d["vhat"], th / 2)))
        mods.append(d)
    return mods


def link_chains(mods, tol=1e-3):
    """Цепочки — модули, стоящие встык (выход одного = вход следующего с точностью tol, касательные совпадают):
    продольно-градиентный диполь из нескольких модулей — одна цепочка, их краевые поля складываются."""
    if not mods:
        return []
    Pin = np.array([m["P_in"] for m in mods]); Pout = np.array([m["P_out"] for m in mods])
    Tin = np.array([m["t_in"] for m in mods]); Tout = np.array([m["t_out"] for m in mods])
    nxt = {}; has_prev = set()
    for i in range(len(mods)):
        dist = np.linalg.norm(Pin - Pout[i], axis=1)
        dist[i] = np.inf
        j = int(np.argmin(dist))
        if dist[j] < tol and float(Tin[j] @ Tout[i]) > 0.0 and j not in has_prev:
            nxt[i] = j; has_prev.add(j)
    chains = []
    for i in range(len(mods)):
        if i in has_prev:
            continue
        ch = [i]
        while ch[-1] in nxt and len(ch) <= len(mods):
            ch.append(nxt[ch[-1]])
        chains.append([mods[k] for k in ch])
    return chains


def _cum_hard(x, L):
    return np.clip(x + L / 2.0, 0.0, L)


def chain_trace(chain, g, et=None, fringe=True, profiles=None, step_frac=8.0):
    """Траектория электрона через цепочку модулей. Кривизна κ(s) = Σ κ_i·f_i(s−S_i) вокруг вертикали КАЖДОГО модуля
    (перекос модулей учитывается), f_i — профиль карты (fringe) или жёсткий край; внутри шага берётся среднее поле
    (точный угол). Старт — состояние орбиты на входе первого модуля (et) или жёсткий вход дуги; до мягкого края
    электрон летит прямо. profiles — явные профили (тест), иначе map_profile по типу устройства."""
    n = len(chain)
    S = [chain[0]["L"] / 2.0]
    for i in range(1, n):
        S.append(S[-1] + float(np.linalg.norm(chain[i]["C"] - chain[i - 1]["C"])))
    prof = []
    for i, m in enumerate(chain):
        pr = None
        if fringe and m["kind"] == "dipole":
            pr = profiles[i] if profiles is not None else map_profile(m["dev"])
        prof.append(pr)
    ext_in = max(0.0, prof[0]["half"] - chain[0]["L"] / 2.0) if prof[0] else 0.0
    ext_out = max(0.0, prof[-1]["half"] - chain[-1]["L"] / 2.0) if prof[-1] else 0.0
    s_lo, s_hi = -ext_in, S[-1] + chain[-1]["L"] / 2.0 + ext_out
    rho_min = min(abs(m["rho"]) for m in chain)
    ds0 = rho_min / (step_frac * g)
    nstep = max(1, int(math.ceil((s_hi - s_lo) / ds0)))
    ds = (s_hi - s_lo) / nstep
    edges = s_lo + np.arange(nstep + 1) * ds
    om = np.zeros((nstep, 3))
    for i, m in enumerate(chain):
        kap = m["theta0"] / m["L"]
        xe = edges - S[i]
        if prof[i] is None:
            cum = _cum_hard(xe, m["L"])
        else:
            pr = prof[i]
            cum = np.interp(xe, pr["x"], pr["cum"], left=0.0, right=pr["cum"][-1]) * (m["L"] / pr["integral"])
        fav = np.diff(cum) / ds                                    # среднее f по шагу
        om += (kap * fav)[:, None] * np.asarray(m["vhat"], float)[None, :]
    m0 = chain[0]
    P0, t0 = np.asarray(m0["P_in"], float), np.asarray(m0["t_in"], float)
    if et is not None and m0["iid"] in et["inst2mag"]:
        j = et["inst2mag"][m0["iid"]]
        P0, t0 = np.asarray(et["pos_in"][j], float), _unit(et["dir_in"][j])
    r = tuple(P0 - ext_in * t0); t = tuple(t0)
    R = np.empty((nstep, 3)); T = np.empty((nstep, 3)); A = np.empty((nstep, 3))
    for k in range(nstep):
        w = om[k]; wn = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
        if wn > 0.0:
            ax = (w[0] / wn, w[1] / wn, w[2] / wn)
            tm = _rotf(t, ax, 0.5 * wn * ds); t1 = _rotf(t, ax, wn * ds)
        else:
            tm = t1 = t
        rm = (r[0] + 0.25 * ds * (t[0] + tm[0]), r[1] + 0.25 * ds * (t[1] + tm[1]), r[2] + 0.25 * ds * (t[2] + tm[2]))
        R[k] = rm; T[k] = tm
        A[k] = (w[1] * tm[2] - w[2] * tm[1], w[2] * tm[0] - w[0] * tm[2], w[0] * tm[1] - w[1] * tm[0])
        r = (r[0] + 0.5 * ds * (t[0] + t1[0]), r[1] + 0.5 * ds * (t[1] + t1[1]), r[2] + 0.5 * ds * (t[2] + t1[2]))
        t = t1
    beta = math.sqrt(1.0 - 1.0 / g ** 2)
    keep = np.linalg.norm(om, axis=1) > 0.0
    vmean = _unit(np.sum([np.asarray(m["vhat"], float) for m in chain], axis=0))
    return dict(kind="bend", name=f"{m0['dev']}#{m0['iid']}" + (f"+{n - 1}" if n > 1 else ""),
                r=R[keep], t=T[keep], acc=beta ** 2 * _C * A[keep], dt=np.full(int(keep.sum()), ds / (beta * _C)),
                vhat=vmean, exit_r=np.array(r), exit_t=np.array(t), ds=ds,
                soft=[pr is not None for pr in prof], iids=[m["iid"] for m in chain],
                int_B2=float(np.sum(np.linalg.norm(om, axis=1) ** 2) * ds))   # ∫κ²ds (→ ∫B²ds·(1/Bρ)²)


# ─────────────────────────── проекция на экран ───────────────────────────

def _corners(scr):
    hw, hh = scr["w"] / 2e3, scr["h"] / 2e3
    return np.array([scr["C"] + a * hw * scr["u"] + b * hh * scr["v"] for a in (-1, 1) for b in (-1, 1)])


def cull_mask(tr, scr, g, cull):
    """Маска точек траектории, чья касательная достаёт экран с запасом cull/γ (охватывающий прямоугольник углов
    к углам экрана в поперечной рамке касательной — заведомо шире самого экрана)."""
    K = tr["r"].shape[0]
    if cull is None or cull <= 0 or K == 0:
        return np.ones(K, bool)
    cor = _corners(scr)
    t = tr["t"]
    ref = np.where(np.abs(t[:, 2:3]) < 0.9, _UP[None, :], np.array([[0.0, 1.0, 0.0]]))
    e1 = np.cross(t, ref); e1 /= np.linalg.norm(e1, axis=1, keepdims=True)
    e2 = np.cross(t, e1)
    d = cor[None, :, :] - tr["r"][:, None, :]
    d /= np.linalg.norm(d, axis=2, keepdims=True)
    a1 = (d * e1[:, None, :]).sum(2); a2 = (d * e2[:, None, :]).sum(2); fw = (d * t[:, None, :]).sum(2)
    g1 = np.maximum(0.0, np.maximum(a1.min(1), -a1.max(1)))
    g2 = np.maximum(0.0, np.maximum(a2.min(1), -a2.max(1)))
    return (np.hypot(g1, g2) <= cull / g) & (fw.max(1) > 0.0)


N_THREADS = max(1, min(16, (os.cpu_count() or 2)))


def _power_chunk(r, t, a, dt, px, py, pz, nrm, beta, omb, g2):
    """Вклад пачки точек траектории во все пиксели (покомпонентно: без массивов (k, M, 3))."""
    dx = px[None, :] - r[:, 0:1]; dy = py[None, :] - r[:, 1:2]; dz = pz[None, :] - r[:, 2:3]
    R2 = dx * dx; R2 += dy * dy; R2 += dz * dz
    iR = 1.0 / np.sqrt(R2)
    dx *= iR; dy *= iR; dz *= iR                                   # → n
    cosi = np.abs(dx * nrm[0] + dy * nrm[1] + dz * nrm[2])
    dx -= t[:, 0:1]; dy -= t[:, 1:2]; dz -= t[:, 2:3]               # → n − t
    na = dx * a[:, 0:1]; na += dy * a[:, 1:2]; na += dz * a[:, 2:3]
    dx *= dx; dy *= dy; dz *= dz
    q = dx; q += dy; q += dz; q *= 0.5 * beta; q += omb              # q = (1−β) + β|n−t|²/2
    iq3 = q ** -3
    br = (a * a).sum(1)[:, None] * iq3
    br -= na * na * iq3 / (q * q * g2)
    br *= cosi; br *= iR * iR
    return dt @ br


def power_on_screen(tr, pts, nrm, g, I, chunk=400_000, threads=None):
    """Плотность мощности [Вт/м²] в точках pts (M,3) плоскости с нормалью nrm от траектории tr (Льенар, ближняя
    зона). Пачки точек траектории считаются параллельно потоками (numpy отпускает GIL в поэлементных операциях)."""
    M = pts.shape[0]; K = tr["r"].shape[0]
    if K == 0:
        return np.zeros(M)
    beta = math.sqrt(1.0 - 1.0 / g ** 2); omb = 1.0 / (g ** 2 * (1.0 + beta)); g2 = g * g
    px, py, pz = (np.ascontiguousarray(pts[:, j]) for j in range(3))
    nrm = np.asarray(nrm, float)
    kc = max(1, int(chunk // M))
    jobs = [(tr["r"][i:i + kc], tr["t"][i:i + kc], tr["acc"][i:i + kc], tr["dt"][i:i + kc])
            for i in range(0, K, kc)]
    nth = N_THREADS if threads is None else max(1, int(threads))
    if nth == 1 or len(jobs) == 1:
        parts = [_power_chunk(*j, px, py, pz, nrm, beta, omb, g2) for j in jobs]
    else:
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(max_workers=nth) as ex:
            parts = list(ex.map(lambda j: _power_chunk(*j, px, py, pz, nrm, beta, omb, g2), jobs))
    return np.sum(parts, axis=0) * (I / _E_CH) * _LW


def trace_power(tr, g, I):
    """Полная излучаемая мощность траектории [Вт] (релятивистский Лармор при β̇ ⟂ β: e²γ⁴|β̇|²/(6πε₀c) на
    электрон, × I/e) — эталон для интеграла карты по экрану, перекрывающему весь свет."""
    a2 = np.einsum("kj,kj->k", tr["acc"], tr["acc"])
    return float((I / _E_CH) * _E_CH ** 2 * g ** 4 / (6 * math.pi * _EPS0 * _C) * np.dot(tr["dt"], a2))


def collapse(tr, point):
    """Протяжённость выкл.: все точки траектории стянуты в одну (направления и ускорения сохраняются)."""
    out = dict(tr); out["r"] = np.repeat(np.asarray(point, float)[None, :], tr["r"].shape[0], axis=0)
    return out


def extraction_point(tr, target):
    """Точка траектории, чья касательная ближе всего смотрит на target (точка отбора), и её касательная."""
    d = np.asarray(target, float)[None, :] - tr["r"]
    d /= np.linalg.norm(d, axis=1, keepdims=True)
    k = int(np.argmax((d * tr["t"]).sum(1)))
    return tr["r"][k], tr["t"][k]


# ─────────────────────────── пучок ───────────────────────────

def beam_cov_screen(point, t, ex, ey, scr, optics, emit_x, kappa, sigE):
    """Ковариация 2×2 пятна пучка на плоскости экрана (u, v) [м²]: размер пучка на расстоянии D вдоль t по Твиссу
    источника, перенесённый на (наклонную) плоскость экрана вдоль лучей."""
    t = _unit(t); D = float(np.dot(scr["C"] - point, t))
    bx, ax = optics["betax"], optics["alphax"]; by, ay = optics["betay"], optics["alphay"]
    gx, gy = (1 + ax ** 2) / bx, (1 + ay ** 2) / by
    sx2 = emit_x * (bx - 2 * ax * D + gx * D ** 2) + (sigE * (optics["etax"] + optics["etaxp"] * D)) ** 2
    sy2 = kappa * emit_x * (by - 2 * ay * D + gy * D ** 2)
    N = scr["n"]; tn = float(np.dot(t, N))

    def proj(e):
        e = np.asarray(e, float); ep = e - t * float(np.dot(e, N)) / tn
        return [float(np.dot(ep, scr["u"])), float(np.dot(ep, scr["v"]))]
    J = np.array([proj(ex), proj(ey)]).T
    return J @ np.diag([sx2, sy2]) @ J.T, math.sqrt(sx2), math.sqrt(sy2), D


def gauss_blur(S, cov, du, dv):
    """Свёртка карты S[iu, iv] с 2D-гауссианой ковариации cov (м²) через FFT — точна и при σ меньше пикселя."""
    sig = math.sqrt(max(cov[0, 0] / du ** 2, cov[1, 1] / dv ** 2, 0.0))
    pad = int(math.ceil(6 * sig)) + 2
    Sp = np.pad(S, pad)
    ku = 2 * np.pi * np.fft.fftfreq(Sp.shape[0], du); kv = 2 * np.pi * np.fft.fftfreq(Sp.shape[1], dv)
    KU, KV = np.meshgrid(ku, kv, indexing="ij")
    H = np.exp(-0.5 * (cov[0, 0] * KU ** 2 + 2 * cov[0, 1] * KU * KV + cov[1, 1] * KV ** 2))
    out = np.real(np.fft.ifft2(np.fft.fft2(Sp) * H))
    return out[pad:pad + S.shape[0], pad:pad + S.shape[1]]


# ─────────────────────────── контекст и итоговая карта ───────────────────────────

def sa_lat_id(lattice):
    import screen_analysis as sa
    return sa.lat_id(lattice)


def register_machine_maps(machine):
    """Карты поля устройств машины — из документа (devices[*].map_bin_path), без чтения hardware_library."""
    for d, dt in machine.devices.items():
        pth = (dt.map_bin_path or "").strip()
        _DEV_MAPS[d] = os.path.join(FIELDS_BIN, pth) if pth and pth != "-" else None


def make_context(lattice, ws=None, et=None, gap=5.0, E=None, I=None, ppp=128, ci=None, byid=None, groups=None):
    """Всё, что не зависит от экрана: таблицы машины, модули и цепочки диполей, кэши траекторий и карт.
    lattice — номер решётки (таблицы) или документ машины (synchrotronsim.model.Machine)."""
    import screen_analysis as sa
    if sa.is_machine(lattice):
        register_machine_maps(lattice)
    if ci is None or byid is None:
        ci, byid = sa.load_lattice(lattice, ws)
    groups = sa.load_groups(lattice) if groups is None else groups
    E = sr.P.E if E is None else E; I = sr.P.I if I is None else I
    poles = set().union(*groups.values()) if groups else set()
    mods = bend_modules(ci, byid, E, poles)
    return dict(lattice=lattice, ws=ws, et=et, gap=gap, E=E, I=I, ppp=ppp, ci=ci, byid=byid, groups=groups,
                g=E / sr.MC2, chains=link_chains(mods), tr_cache={}, optics_cache={})


def _und_trace(ctx, gid, poles):
    key = ("und", gid, poles)
    if key not in ctx["tr_cache"]:
        p = sr.P if ctx["E"] == sr.P.E else type("Pe", (sr.P,), dict(E=ctx["E"]))
        ctx["tr_cache"][key] = undulator_trace(gid, ctx["groups"][gid], ctx["ci"], ctx["byid"], ctx["gap"],
                                               ctx["et"], poles, ctx["ppp"], p)
    return ctx["tr_cache"][key]


def _chain_trace(ctx, ic, fringe):
    key = ("bend", ic, fringe)
    if key not in ctx["tr_cache"]:
        ctx["tr_cache"][key] = chain_trace(ctx["chains"][ic], ctx["g"], ctx["et"], fringe)
    return ctx["tr_cache"][key]


def _optics(ctx, point):
    key = tuple(np.round(np.asarray(point, float), 3))
    if key not in ctx["optics_cache"]:
        import lattice_optics as lo
        ctx["optics_cache"][key] = lo.optics_at_point(sa_lat_id(ctx["lattice"]), point)
    return ctx["optics_cache"][key]


def _chain_candidates(ctx, scr, cull):
    """Цепочки, чей веер (с запасом на мягкий край и хвост Льенара) может достать экран."""
    import bend_screen as bs
    g = ctx["g"]; out = []
    for ic, ch in enumerate(ctx["chains"]):
        for m in ch:
            pr = map_profile(m["dev"]) if m["kind"] == "dipole" else None
            ext = max(0.0, pr["half"] - m["L"] / 2.0) if pr else 0.0
            margin = 2e-4 + max(cull, 8.0) / g + ext / abs(m["rho"])
            if bs.reaches_screen(m, scr, g, nfan=max(cull, 8.0), margin=margin):
                out.append(ic)
                break
    return out


def _stats(S, UU, VV, cell):
    tot = float(S.sum() * cell); wsum = float(S.sum()) or 1.0
    cu = float((S * UU).sum() / wsum) * 1e3; cv = float((S * VV).sum() / wsum) * 1e3
    ru = math.sqrt(max(0.0, float((S * (UU * 1e3 - cu) ** 2).sum() / wsum)))
    rv = math.sqrt(max(0.0, float((S * (VV * 1e3 - cv) ** 2).sum() / wsum)))
    return dict(peak=float(S.max() * 1e-6), total=tot, cu=cu, cv=cv, ru=ru, rv=rv)


def footprint_accurate(scr, ctx, effects=None, mesh=161, cull=40.0, emit_x=None, kappa=None, sigE=None,
                       breakdown=False, bends=True, undulators=True):
    """Тепловая карта экрана в ТОЧНОМ режиме. Возврат — те же ключи, что у screen_analysis.footprint
    (S [Вт/мм²], extent, peak, total, und_W, bend_W, n_bend, cu, cv, ru, rv) + calc, effects, illuminators,
    chains, n_samples, time_s, beam (σ пучка у экрана по источникам); при breakdown — dict эффект → итог БЕЗ него."""
    import screen_analysis as sa
    t_start = time.time()
    eff = parse_effects(effects)
    emit_x = sr.P.emit_x if emit_x is None else emit_x
    kappa = sr.P.kappa if kappa is None else kappa
    sigE = sr.P.sig_E if sigE is None else sigE
    g, I = ctx["g"], ctx["I"]
    hw, hh, UU, VV, cell, pts3 = sa._grid(scr, mesh)
    pts = pts3.reshape(-1, 3)
    du, dv = 2 * hw / (mesh - 1), 2 * hh / (mesh - 1)
    raw = {}                                                       # (источник, вариант) → (карта, траектория)
    srcs = []
    if undulators:
        srcs += [("und", gid) for gid in sorted(ctx["groups"])]
    if bends:
        srcs += [("bend", ic) for ic in _chain_candidates(ctx, scr, cull)]
    n_samp = [0]

    def src_map(kind, key, e):
        var = (("extended" in e), ("poles" in e)) if kind == "und" else (("extended" in e), ("fringe" in e))
        mk = (kind, key, var)
        if mk in raw:
            return raw[mk]
        tr = _und_trace(ctx, key, var[1]) if kind == "und" else _chain_trace(ctx, key, var[1])
        if tr is None:
            raw[mk] = None
            return None
        m = cull_mask(tr, scr, g, cull)
        if not m.any():
            raw[mk] = None
            return None
        sub = dict(tr, **{k: tr[k][m] for k in ("r", "t", "acc", "dt")})
        if kind == "und":
            point, tdir = tr["center"], tr["dir"]
            ex, ey = tr["ex"], tr["ey"]
        else:
            point, tdir = extraction_point(sub, scr["C"])
            ex = _unit(np.cross(tr["vhat"], tdir)); ey = _unit(np.cross(tdir, ex))
        if not var[0]:
            sub = collapse(sub, point)
        n_samp[0] += int(m.sum())
        S = power_on_screen(sub, pts, scr["n"], g, I).reshape(mesh, mesh)
        raw[mk] = (S, dict(point=point, t=tdir, ex=ex, ey=ey, name=tr["name"], n=int(m.sum())))
        return raw[mk]

    def assemble(e):
        Su = np.zeros((mesh, mesh)); Sb = np.zeros((mesh, mesh)); ill, chs, bm = [], [], []
        for kind, key in srcs:
            res = src_map(kind, key, e)
            if res is None:
                continue
            S, info = res
            if "beam" in e:
                cov, sx, sy, D = beam_cov_screen(info["point"], info["t"], info["ex"], info["ey"], scr,
                                                 _optics(ctx, info["point"]), emit_x, kappa, sigE)
                S = gauss_blur(S, cov, du, dv)
                bm.append(dict(src=info["name"], sigma_x_um=sx * 1e6, sigma_y_um=sy * 1e6, D=D))
            if kind == "und":
                Su += S; ill.append(key)
            else:
                Sb += S; chs.append(info["name"])
        S = Su + Sb
        st = _stats(S, UU, VV, cell)
        st.update(S=S * 1e-6, und_W=float(Su.sum() * cell), bend_W=float(Sb.sum() * cell), illuminators=ill,
                  chains=chs, n_bend=len(chs), beam=bm)
        return st

    fp = assemble(eff)
    fp.update(extent=(-hw * 1e3, hw * 1e3, -hh * 1e3, hh * 1e3), calc="accurate", effects=sorted(eff),
              n_samples=n_samp[0], mesh=mesh, cull=cull)
    if breakdown:
        fp["breakdown"] = {}
        for e in sorted(eff):
            r = assemble(eff - {e})
            r.pop("S"); fp["breakdown"][e] = r
    fp["time_s"] = time.time() - t_start
    return fp


def stationary_points(tr, pts, chunk=2_000_000):
    """Для точек pts (M,3): участок траектории, чья касательная в плоскости изгиба смотрит на точку (смена знака
    горизонтального угла между касательной и направлением на точку). Возврат (valid (M,), ψ (M,) — вертикальный
    угол от плоскости орбиты, ec (M,) — локальная ε_c [кэВ], интерполированная между соседними точками).
    Плоскость изгиба — tr["vhat"]; ε_c по локальной кривизне: B = |κ|·Bρ, ε_c = 0,665·E²·B."""
    step = max(1, tr["r"].shape[0] // 600)                     # касательная поворачивается плавно — хватит ~600
    r, t, acc = tr["r"][::step], tr["t"][::step], tr["acc"][::step]   # точек (между ними — линейная интерполяция)
    K, M = r.shape[0], pts.shape[0]
    valid = np.zeros(M, bool); psi = np.zeros(M); ec = np.zeros(M)
    if K < 2:
        return valid, psi, ec
    vh = np.asarray(tr["vhat"], float)
    E = sr.P.E
    ec_k = 0.665 * E ** 2 * (np.linalg.norm(acc, axis=1) / _C) * (E / 0.299792458)   # |β̇|/c ≈ |κ| (β≈1)
    mc = max(1, int(chunk // (3 * K)))
    for i in range(0, M, mc):
        p = pts[i:i + mc]
        d = p[None, :, :] - r[:, None, :]
        dv = d @ vh
        dh = d - dv[..., None] * vh
        cr = np.cross(t[:, None, :], dh) @ vh
        dt_ = np.einsum("kmj,kj->km", dh, t)
        dl = np.arctan2(cr, dt_)                                    # горизонт. угол касательная→точка
        ch = (np.sign(dl[:-1]) != np.sign(dl[1:])) & (dt_[:-1] > 0)
        ok = ch.any(axis=0)
        k = np.argmax(ch, axis=0)
        cols = np.arange(p.shape[0])
        a0, a1 = dl[k, cols], dl[k + 1, cols]
        f = np.where(np.abs(a0 - a1) > 0, a0 / (a0 - a1), 0.5)
        dk = d[k, cols]
        valid[i:i + mc] = ok
        psi[i:i + mc] = np.arcsin(np.clip(dv[k, cols] / np.linalg.norm(dk, axis=1), -1, 1))
        ec[i:i + mc] = ec_k[k] + f * (ec_k[np.minimum(k + 1, K - 1)] - ec_k[k])
    return valid, psi, ec


def trace_flux_map(traces, scr, E_keV, mesh, I=None):
    """Карта ПОТОКА поворотного СИ на энергии E_keV [ф/с/мрад²/0.1%BW] по траекториям (протяжённый источник,
    мягкий край: у края поля ε_c падает → поток на фиксированной энергии плавно гаснет; стыков модулей нет)."""
    import screen_analysis as sa
    import bend_sr as bsr
    I = sr.P.I if I is None else I
    hw, hh, UU, VV, cell, pts3 = sa._grid(scr, mesh)
    pts = pts3.reshape(-1, 3)
    S = np.zeros(pts.shape[0])
    g = sr.gamma()
    for tr in traces:
        ok, psi, ec = stationary_points(tr, pts)
        ok &= ec > 0
        if ok.any():
            S[ok] += bsr.flux_angular(float(E_keV), psi[ok], ec[ok], sr.P.E, I, g)
    S = S.reshape(mesh, mesh)
    return dict(S=S, extent=(-hw * 1e3, hw * 1e3, -hh * 1e3, hh * 1e3), peak=float(S.max()))


def trace_point_spectrum(traces, scr, u_mm, v_mm, n=400, I=None):
    """Сплошной спектр поворотного СИ в точке экрана (u, v мм) по траекториям: на каждую — участок, чья
    касательная смотрит на точку, со своей ε_c. Возврат dict(E_keV, flux [ф/с/мрад²/0.1%BW], ecs)."""
    import bend_sr as bsr
    I = sr.P.I if I is None else I
    p = (scr["C"] + (u_mm * 1e-3) * scr["u"] + (v_mm * 1e-3) * scr["v"])[None, :]
    hits = []
    for tr in traces:
        ok, psi, ec = stationary_points(tr, p)
        if bool(ok[0]) and ec[0] > 0:
            hits.append((float(psi[0]), float(ec[0])))
    ecs = [e for _, e in hits] or [float(scr.get("ec") or 10.0)]
    Es = np.geomspace(0.03 * min(ecs), 12.0 * max(ecs), int(n))
    F = np.zeros_like(Es)
    g = sr.gamma()
    for psi, e in hits:
        F += bsr.flux_angular(Es, psi, e, sr.P.E, I, g)
    return dict(E_keV=Es, flux=F, ecs=ecs)


def compare_line(fast, acc):
    """Короткая строка «точный против быстрого»."""
    dp = (acc["peak"] / fast["peak"] - 1) * 100 if fast.get("peak") else float("nan")
    dt = (acc["total"] / fast["total"] - 1) * 100 if fast.get("total") else float("nan")
    return (f"точный: пик {acc['peak']:7.2f} Вт/мм² ({dp:+.2f}%), перехв. {acc['total']:7.1f} Вт ({dt:+.2f}%), "
            f"центроид ({acc['cu']:+.3f},{acc['cv']:+.3f}) мм [Δ {acc['cu'] - fast.get('cu', 0):+.3f},"
            f"{acc['cv'] - fast.get('cv', 0):+.3f}], {acc['n_samples']} точек траектории, {acc['time_s']:.1f} с")


def breakdown_lines(acc):
    """Таблица вкладов эффектов: итог с эффектом минус итог без него."""
    out = []
    for e, r in acc.get("breakdown", {}).items():
        dp = acc["peak"] - r["peak"]; dt = acc["total"] - r["total"]
        out.append(f"      {e:9s}: Δпик {dp:+8.3f} Вт/мм² ({dp / (r['peak'] or 1) * 100:+7.3f}%), "
                   f"Δмощн. {dt:+8.2f} Вт ({dt / (r['total'] or 1) * 100:+7.3f}%), "
                   f"Δцентроид ({acc['cu'] - r['cu']:+.4f},{acc['cv'] - r['cv']:+.4f}) мм, "
                   f"Δσ ({acc['ru'] - r['ru']:+.4f},{acc['rv'] - r['rv']:+.4f}) мм")
    return out
