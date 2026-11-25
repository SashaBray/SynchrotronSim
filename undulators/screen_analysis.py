# -*- coding: utf-8 -*-
"""Анализ ФОТОННЫХ ЭКРАНОВ — ЭТАП A: геометрия + мощность + АВТОПОИСК осветителей (без интерференции).

Для решётки (НОМИНАЛЬНОЙ или ВОЗМУЩЁННОЙ через --ws) и замороженных экранов (screens.csv):
  1) АВТОПОИСК. Какие ондуляторы РЕАЛЬНО светят на экран — тест видимости по конусу от ФАКТИЧЕСКОЙ
     геометрии источника (из lattice_configs ws). Привязка экрана к источнику — лишь стартовое зерно;
     экран берёт всех, чей конус его достаёт.
  2) МОЩНОСТЬ. Проекция угловой плотности мощности каждого осветителя на плоскость экрана и
     НЕКОГЕРЕНТНАЯ сумма по интенсивности → карта Вт/мм². (Гармоник тут нет — они на этапе B,
     когерентный интеграл поля по траектории.)
  3) ИТОГ. Пик Вт/мм², перехваченная мощность, смещение/размер центроида. Перекос машины (--ws)
     двигает ПЯТНО по неподвижному экрану — видно в центроиде.

ИСТОЧНИК СВЕТА (--source orbit, по умолчанию): свет летит по ЭЛЕКТРОНУ, поэтому точка и ось конуса ИД и
касательная веера диполя берутся по ЗАМКНУТОЙ ОРБИТЕ машины (electron_source.py: C++-ядро, рамки элементов
из таблиц ws; --correct — после SVD-коррекции орбиты). Плоскость колебаний ИД — по рамке железа (её задаёт
поле). Жёсткий наклон ондулятора электрон почти не поворачивает (нулевой интеграл поля) — пятно стоит;
пятно двигает угол орбиты от других ошибок. --source magnet — прежнее «свет по оси железа» (для сравнения).

МОДЕЛЬ ИЗЛУЧЕНИЯ (этап A). На устройство — ОДИН эмиттер (все полюса агрегированы). Угловая плотность
мощности — ТОЧНАЯ (по умолчанию, --model lienard): формула Льенара, усреднённая по периоду синусоиды
(undulator_sr.angular_power_table). Пятно вытянуто в плоскости колебаний (±K/γ) и уже поперёк (~1/γ), хвосты
степенные. Интерференция в мощность не входит (перераспределяет энергию по частотам, не по углам —
Парсеваль), поэтому кольца карт потока одной энергии (этап B) в тепловой карте законно отсутствуют.
--model gauss — прежний круглый гаусс σ_θ=√(P/2π·P') (для сравнения; форма и хвосты неверны).
ВЕЕРА ДИПОЛЕЙ (по умолчанию; --no-bends — выкл.): поворотное СИ диполей, чей веер достаёт ХОТЬ ЧАСТЬ экрана
(диполь перед прямолинейкой смотрит выходом вдоль оси ИД), складывается в ту же карту — модель bend_screen.
Полюса ИД из поиска диполей исключены (они уже в эмиттере ИД).
Зазор всех устройств задаётся --gap (как при вставке insert_undulators); позже — по инстансу.

РЕЖИМЫ СЧЁТА (--calc). fast (по умолчанию) — всё описанное выше: точечный эмиттер ИД с табличной угловой
мощностью и веера из точки отбора; ~1 с на экран. accurate — прямое интегрирование формулы Льенара по траектории
электрона (undulators/sr_trace.py): протяжённый источник, реальная структура полюсов ИД (концевые полюса),
мягкий край диполей по карте поля, свёртка с размером пучка; ~10–60 с на экран, --breakdown — вклад каждого
эффекта (--effects all,-beam … — выключить эффект). В точном режиме строка быстрого печатается для сравнения.

ЗАПУСК
  py -3.14 undulators/screen_analysis.py --lattice 27 --gap 5
  py -3.14 undulators/screen_analysis.py --lattice 27 --calc accurate --breakdown   # точный режим + вклады
  py -3.14 undulators/screen_analysis.py --lattice 27 --ws <папка_возмущённой_машины>   # увод при перекосе
  py -3.14 undulators/screen_analysis.py --lattice 27 --model gauss --no-bends           # прежняя модель
"""
import os, sys, csv, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt
import undulator_sr as sr

TB = os.path.join("workspace", "tables")
EG = os.path.join(TB, "element_groups.csv")
SCR = os.path.join(TB, "screens.csv")
RESDIR = os.path.join("workspace", "results", "screens")
RADIATE = {"dipole", "dq"}                           # что светит (изгибы; полюса ИД — свой эмиттер)


def _unit(v):
    v = np.asarray(v, float); n = np.linalg.norm(v)
    return v / n if n else v


# ---- машина: номер решётки (таблицы, как раньше) ИЛИ документ synchrotronsim.model.Machine (docs/DATA_PROTOCOL.md:
# модель во время счёта берёт данные из документа и таблиц не открывает). Адаптеры отдают ТЕ ЖЕ структуры, что
# чтение таблиц: строки lattice_configs — строками (repr числа точно возвращает значение, «-0» — тоже), пусто — "".
_LC_COLS = ("lattice_id", "instance_id", "device_type_id", "x", "y", "z", "yaw", "pitch", "roll", "arg_val",
            "correct_coeff", "mat_B", "mat_G", "mat_S", "mat_O", "mat_Gv", "eff_length", "kind", "n_slices", "edge1",
            "edge2")
_LC_FIELD = {"device_type_id": "device_type"}
_M_ROWS = {}                   # отпечаток машины → (ci, byid)


def is_machine(x) -> bool:
    return hasattr(x, "instances") and hasattr(x, "fingerprint")


def lat_id(lattice) -> int:
    """Номер решётки: у документа машины — его lattice_id."""
    return int(lattice.lattice_id) if is_machine(lattice) else int(lattice)


def lat_key(lattice):
    """Ключ кэшей: номер решётки, у документа — (номер, отпечаток): разные ревизии не смешиваются."""
    return (int(lattice.lattice_id), lattice.fingerprint) if is_machine(lattice) else int(lattice)


def _cell(v):
    if isinstance(v, float) or isinstance(v, np.floating):
        return "" if v != v else repr(float(v))
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    return "" if v is None else str(v)


def machine_rows(m):
    """Строки lattice_configs машины в прежнем виде (ci, byid) — из документа, без таблиц."""
    hit = _M_ROWS.get(m.fingerprint)
    if hit is None:
        I = m.instances
        cols = [np.asarray(getattr(I, _LC_FIELD.get(c, c))) if c != "lattice_id" else None for c in _LC_COLS]
        lat = str(int(m.lattice_id))
        byid = {}
        for j in range(len(I)):
            row = [lat if a is None else _cell(a[j]) for a in cols]
            byid[row[1]] = row
        hit = _M_ROWS[m.fingerprint] = ({c: i for i, c in enumerate(_LC_COLS)}, byid)
    return hit


def machine_screen_rows(m):
    """Строки screens машины в табличных единицах (мм, мрад) — для прежних потребителей строк."""
    out = []
    for s in m.screens:
        r = {"lattice_id": str(int(m.lattice_id)), "screen_id": s.screen_id, "source_ref": s.source_ref,
             "source_kind": s.source_kind, "mode": s.mode, "dist_m": repr(float(s.dist)),
             "width_mm": repr(float(s.width) / 1e-3), "height_mm": repr(float(s.height) / 1e-3)}
        for k, v in (("c", s.center), ("n", s.normal), ("u", s.u), ("v", s.v)):
            for a, x in zip("xyz", v):
                r[k + a] = repr(float(x))
        for k in ("extract", "los_h", "los_v", "tilt_h", "tilt_v", "roll"):
            r[k + "_mrad"] = repr(float(getattr(s, k)) / 1e-3)
        out.append(r)
    return out


def load_groups(lattice):
    """{group_id: set(instance_ids)} для решётки (или документа машины)."""
    if is_machine(lattice):
        return {g.group_id: set(g.instance_ids) for g in lattice.groups}
    out = {}
    for r in csv.DictReader(open(EG, encoding="utf-8-sig"), delimiter=";"):
        if r.get("lattice_id") == str(lattice):
            out[r["group_id"]] = set(r["instance_ids"].split("|"))
    return out


def load_screens(lattice):
    if is_machine(lattice):
        return [dict(sid=s.screen_id, src=s.source_ref, C=np.array(s.center, float), n=np.array(s.normal, float),
                     u=np.array(s.u, float), v=np.array(s.v, float), w=float(s.width) / 1e-3,
                     h=float(s.height) / 1e-3, mode=s.mode) for s in lattice.screens]
    if not os.path.exists(SCR):
        raise SystemExit(f"нет {SCR} — сначала разместите экраны (build_screens.py)")
    out = []
    for r in csv.DictReader(open(SCR, encoding="utf-8-sig"), delimiter=";"):
        if r.get("lattice_id") != str(lattice):
            continue
        g = lambda k: np.array([float(r[k + a]) for a in ("x", "y", "z")])
        out.append(dict(sid=r["screen_id"], src=r["source_ref"],
                        C=g("c"), n=g("n"), u=g("u"), v=g("v"),
                        w=float(r["width_mm"]), h=float(r["height_mm"]), mode=r["mode"]))
    return out


def load_lattice(lattice, ws):
    if is_machine(lattice):                                  # документ машины (возмущённая машина — свой документ)
        return machine_rows(lattice)
    path = os.path.join(ws, "tables", "lattice_configs.csv") if ws else os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.reader(open(path, encoding="utf-8-sig"), delimiter=";"))
    ci = {c: i for i, c in enumerate(rows[0])}
    byid = {r[ci["instance_id"]]: r for r in rows[1:] if r and r[ci["lattice_id"]] == str(lattice)}
    return ci, byid


def emitter(insts, ci, byid):
    """Точка источника и ось устройства из ФАКТических строк полюсов (kind∈RADIATE)."""
    P, axis, _, _ = emitter_frame(insts, ci, byid)
    return P, axis


def emitter_frame(insts, ci, byid):
    """Точка источника и рамка устройства (ось, u — плоскость колебаний, v = ось×u) по ФАКТ. полюсам.
    u берётся из рамок полюсов (знак выровнен, ортогонализован к оси): крен устройства поворачивает эллипс."""
    cs, axes, us = [], [], []
    for iid in insts:
        r = byid.get(iid)
        if r is None or r[ci["kind"]].strip() not in RADIATE:
            continue
        C = np.array([float(r[ci[a]] or 0) for a in ("x", "y", "z")])
        ypr = (float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))
        ax, uh, _ = mt.frame(*ypr)
        cs.append(C); axes.append(np.asarray(ax, float)); us.append(np.asarray(uh, float))
    if not cs:
        return None, None, None, None
    axis = _unit(np.sum(axes, axis=0))
    u = sum(x if np.dot(x, us[0]) >= 0 else -x for x in us)
    u = _unit(u - np.dot(u, axis) * axis)
    return np.mean(cs, axis=0), axis, u, np.cross(axis, u)


def power_params(gap):
    B0 = sr.B0_of_gap(gap); K = sr.K_of_B0(B0)
    P_tot = sr.total_power_kW(B0) * 1e3                       # Вт
    Pd0 = sr.power_density_onaxis(B0, K) * 1e9               # кВт/мрад² → Вт/рад²
    sig = math.sqrt(P_tot / (2 * math.pi * Pd0))             # рад, эффективный конус (только модель gauss)
    return B0, K, P_tot, Pd0, sig


def make_emitters(groups, ci, byid, gap, model="lienard", et=None):
    """{gid: эмиттер} по фактической геометрии. Эмиттер — dict(P, axis, u, v, tab|Pd0+sig).
    et (electron_source.electron_track) — точка и ось берутся по ОРБИТЕ электрона, плоскость колебаний u —
    по рамке железа (поле задаёт плоскость, в которой качается электрон); без et — по рамке железа."""
    B0, K, P_tot, Pd0, sig = power_params(gap)
    tab = sr.angular_power_cached(K) if model == "lienard" else None
    out = {}
    for gid, insts in groups.items():
        P, axis, u, v = emitter_frame(insts, ci, byid)
        if P is None:
            continue
        if et is not None:
            import electron_source as es
            Pe, de = es.device_source(et, insts)
            if Pe is not None:
                P, axis = Pe, de
                u = _unit(u - np.dot(u, axis) * axis); v = np.cross(axis, u)
        out[gid] = dict(gid=gid, P=P, axis=axis, u=u, v=v, tab=tab, Pd0=Pd0, sig=sig)
    return out


def find_illuminators(scr, emit, cone_cut=6.0):
    """Ондуляторы, чья угловая мощность достаёт ХОТЬ ЧАСТЬ экрана (по таблице Льенара или по cone_cut·σ)."""
    found = []
    for gid, em in emit.items():
        vec = scr["C"] - em["P"]; D = float(np.linalg.norm(vec)); dh = vec / D
        ca = float(np.dot(dh, em["axis"]))
        if ca <= 0:
            continue                                         # экран позади источника
        half = 0.5 * math.hypot(scr["w"], scr["h"]) / 1e3 / D
        if em["tab"] is not None:
            tx = abs(float(np.dot(dh, em["u"]))); ty = abs(float(np.dot(dh, em["v"])))
            ok = tx <= em["tab"]["tx"][-1] + half and ty <= em["tab"]["ty"][-1] + half
        else:
            ok = math.acos(max(-1.0, min(1.0, ca))) <= cone_cut * em["sig"] + half
        if ok:
            found.append(gid)
    return found


def find_bends(scr, ci, byid, exclude, E=None, et=None):
    """Веера диполей, достающие экран (полюса ИД — exclude — не считаются: они в эмиттере ИД).
    et — веер по ОРБИТЕ электрона (смещение и поворот касательной в точке отбора), иначе по дуге железа."""
    import bend_screen as bs                                 # ленивый импорт: bend_screen импортирует этот модуль
    E = sr.P.E if E is None else E
    g = E / sr.MC2
    out = []
    for iid, r in byid.items():
        if iid in exclude:
            continue
        d = bs.dipole_arc(r, ci, E)
        if d is None or not bs.reaches_screen(d, scr, g):
            continue
        em = bs.emitter_geom(d, scr["C"], g)
        if et is not None:
            import electron_source as es
            em = es.fan_on_orbit(em, et)
        out.append(em)
    return out


def _grid(scr, mesh):
    hw, hh = scr["w"] / 2e3, scr["h"] / 2e3                  # полуразмеры, м
    su = np.linspace(-hw, hw, mesh); sv = np.linspace(-hh, hh, mesh)
    UU, VV = np.meshgrid(su, sv, indexing="ij")
    cell = (2 * hw / (mesh - 1)) * (2 * hh / (mesh - 1))     # площадь ячейки, м²
    pts = (scr["C"][None, None, :] + UU[..., None] * scr["u"][None, None, :]
           + VV[..., None] * scr["v"][None, None, :])        # (m,m,3) глобальные точки
    return hw, hh, UU, VV, cell, pts


def footprint(scr, emitters, mesh, bends=None, E=None, I=None):
    """Карта плотности мощности на плоскости экрана [Вт/мм²] + сводка. emitters — эмиттеры ИД (make_emitters),
    bends — веера диполей (find_bends) или None. Ондуляторы и диполи складываются по мощности."""
    hw, hh, UU, VV, cell, pts = _grid(scr, mesh)
    S = np.zeros((mesh, mesh))
    for em in emitters:
        d = pts - em["P"][None, None, :]
        D = np.linalg.norm(d, axis=2)
        dhat = d / D[..., None]
        ca = (dhat * em["axis"][None, None, :]).sum(2)
        if em["tab"] is not None:                            # точная угловая мощность (Льенар)
            tx = (dhat * em["u"][None, None, :]).sum(2); ty = (dhat * em["v"][None, None, :]).sum(2)
            dPdO = sr.angular_power(em["tab"], tx, ty)       # Вт/рад²
        else:                                                # прежний круглый гаусс
            theta = np.arccos(np.clip(ca, -1.0, 1.0))
            dPdO = em["Pd0"] * np.exp(-theta ** 2 / (2 * em["sig"] ** 2))
        cos_inc = np.abs((dhat * scr["n"][None, None, :]).sum(2))
        S += np.where(ca > 0, dPdO * cos_inc / D ** 2, 0.0)
    und_W = float(S.sum() * cell)
    bend_W = 0.0
    if bends:
        import bend_screen as bs
        E = sr.P.E if E is None else E; I = sr.P.I if I is None else I
        Sb = bs.footprint(scr, bends, E / sr.MC2, I, mesh)["S"] * 1e6   # Вт/мм² → Вт/м²
        bend_W = float(Sb.sum() * cell)
        S = S + Sb
    tot = float(S.sum() * cell)                              # Вт (перехвачено)
    wsum = float(S.sum()) or 1.0
    cu = float((S * UU).sum() / wsum) * 1e3                  # центроид, мм
    cv = float((S * VV).sum() / wsum) * 1e3
    ru = math.sqrt(max(0.0, float((S * (UU * 1e3 - cu) ** 2).sum() / wsum)))
    rv = math.sqrt(max(0.0, float((S * (VV * 1e3 - cv) ** 2).sum() / wsum)))
    return dict(S=S * 1e-6, extent=(-hw * 1e3, hw * 1e3, -hh * 1e3, hh * 1e3),  # S в Вт/мм²
                peak=float(S.max() * 1e-6), total=tot, und_W=und_W, bend_W=bend_W, n_bend=len(bends or []),
                cu=cu, cv=cv, ru=ru, rv=rv)


SCREEN_DEFAULT_DIST = 20.0     # м: экраны кокпита (client.scene3d, control.recorders) — у всех ондуляторов и линеек
SCREEN_DEFAULT_SIZE_MM = 50.0
BEND_LINE_MODULES = 5          # линейка продольно-градиентного диполя: модулей встык
BEND_SCREEN_MM = (100.0, 30.0)  # ширина × высота экрана у линейки (веер модуля №3 на 20 м ≈ 84 мм)
_LIVE = {}                     # (решётка, экран, веера) → таблицы, геометрия экрана, веера
_BEND_SCR = {}                 # (решётка, dist) → экраны у центральных модулей линеек


def _screen_rows(lattice=None):
    """Строки экранов: у документа машины — только её экраны; по номеру — вся таблица screens (все решётки)."""
    if lattice is not None and is_machine(lattice):
        return machine_screen_rows(lattice)
    p = SCR if os.path.exists(SCR) else None
    if p is None:
        return []
    with open(p, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=";"))


def bend_line_screens(lattice, dist=SCREEN_DEFAULT_DIST, ci=None, byid=None, groups=None, E=None):
    """Экраны у ЦЕНТРАЛЬНЫХ модулей (№3) линеек диполя из BEND_LINE_MODULES модулей встык (продольно-градиентный
    диполь: цепочки sr_trace.link_chains). Источник — центр дуги модуля №3, ось — касательная там (туда летит его
    СИ), экран поперёк на dist, u — горизонталь (плоскость веера). Имя SCR_DP<instance_id модуля>. Кэш по решётке."""
    key = (lat_key(lattice), float(dist))
    if key in _BEND_SCR:
        return _BEND_SCR[key]
    import sr_trace as st
    from build_screens import screen_frame
    if ci is None or byid is None:
        ci, byid = load_lattice(lattice, None)
    groups = load_groups(lattice) if groups is None else groups
    E = sr.P.E if E is None else E
    poles = set().union(*groups.values()) if groups else set()
    out, lines, rest = [], [], []
    for ch in st.link_chains(st.bend_modules(ci, byid, E, poles)):
        (lines if len(ch) == BEND_LINE_MODULES else rest).append(ch)
    lines = [(ch, ch[BEND_LINE_MODULES // 2]) for ch in lines]
    # линейка, разбитая вставкой (инжекция: крайний модуль через ~0,4 м): модули одного семейства с номерами
    # _1.._5 в типе устройства — центральный «_3» (номер — порядок в линейке)
    fam = {}
    for ch in rest:
        for m in ch:
            if m["kind"] == "dipole" and "_" in m["dev"] and m["dev"].rsplit("_", 1)[1].isdigit():
                fam.setdefault(m["dev"].rsplit("_", 1)[0], []).append(m)
    used = set()
    for pre, ms in fam.items():
        nums = sorted(int(m["dev"].rsplit("_", 1)[1]) for m in ms)
        if nums == list(range(1, BEND_LINE_MODULES + 1)):
            mid = next(m for m in ms if int(m["dev"].rsplit("_", 1)[1]) == BEND_LINE_MODULES // 2 + 1)
            lines.append((sorted(ms, key=lambda m: int(m["dev"].rsplit("_", 1)[1])), mid))
            used |= {m["iid"] for m in ms}
    # без номеров в имени (реш.27): куски одного типа, сложенные в BEND_LINE_MODULES, со вставкой короче модуля
    rest = [ch for ch in rest if not ({m["iid"] for m in ch} & used)]
    for i, a in enumerate(rest):
        for b in rest[i + 1:]:
            if len(a) + len(b) != BEND_LINE_MODULES or a[0]["dev"] != b[0]["dev"] or a[0]["kind"] != "dipole":
                continue
            Lmax = max(m["L"] for m in a + b)
            for first, second in ((a, b), (b, a)):
                if float(np.linalg.norm(second[0]["P_in"] - first[-1]["P_out"])) < Lmax:
                    line = first + second
                    lines.append((line, line[BEND_LINE_MODULES // 2]))
    for ch, m in lines:
        P = np.asarray(m["C"], float); t = _unit(m["axis"])
        n = -t; u, v = screen_frame(n)
        out.append(dict(sid=f"SCR_DP{m['iid']}", src=str(m["iid"]), C=P + float(dist) * t, n=n, u=u, v=v,
                        w=BEND_SCREEN_MM[0], h=BEND_SCREEN_MM[1], mode="bend", kind="bend", P=P, t=t,
                        dev=m["dev"], theta=float(m["theta0"]), ec=float(m["ec"]), D=float(dist),
                        line=[c["iid"] for c in ch]))
    _BEND_SCR[key] = out
    return out


def live_screen(lattice, screen_id, groups, ci, byid):
    """Геометрия экрана для живой карты (кокпит). Есть замороженная строка этой решётки — она; экран линейки
    диполя (SCR_DP…) — bend_line_screens; иначе экран ондулятора строится так же, как его рисует сцена кокпита:
    по оси ондулятора-источника (рамка железа) на dist_m строки с тем же screen_id из любой решётки (нет —
    SCREEN_DEFAULT_DIST), размер — из той же строки (нет — 50×50 мм), плоскость поперёк оси. None — не построить."""
    rows = _screen_rows(lattice)
    for r in rows:
        if r.get("lattice_id") == str(lat_id(lattice)) and r.get("screen_id") == screen_id:
            return next(s for s in load_screens(lattice) if s["sid"] == screen_id)
    if str(screen_id).startswith("SCR_DP"):
        return next((s for s in bend_line_screens(lattice, ci=ci, byid=byid, groups=groups) if s["sid"] == screen_id),
                    None)
    row = next((r for r in rows if r.get("screen_id") == screen_id), None)
    src = (row or {}).get("source_ref", "").strip() or (screen_id[4:] if screen_id.startswith("SCR_") else "")
    if (row or {}).get("mode") == "bend" or src not in groups:
        return None
    P, axis, _, _ = emitter_frame(groups[src], ci, byid)
    if P is None:
        return None
    f = lambda k, d: float((row or {}).get(k) or 0) or d
    from build_screens import screen_frame
    n = -axis; u, v = screen_frame(n)
    return dict(sid=screen_id, src=src, C=P + f("dist_m", SCREEN_DEFAULT_DIST) * axis, n=n, u=u, v=v,
                w=f("width_mm", SCREEN_DEFAULT_SIZE_MM), h=f("height_mm", SCREEN_DEFAULT_SIZE_MM), mode="axis")


def live_context(lattice, screen_id, bends=True):
    """Кэш живого экрана: таблицы машины, геометрия экрана, веера диполей, достающие экран. None — не построить."""
    key = (lat_key(lattice), str(screen_id), bool(bends))
    c = _LIVE.get(key)
    if c is None:
        ci, byid = load_lattice(lattice, None)
        groups = load_groups(lattice)
        scr = live_screen(lattice, screen_id, groups, ci, byid)
        if scr is None:
            return None
        poles = set().union(*groups.values()) if groups else set()
        c = dict(ci=ci, byid=byid, groups=groups, scr=scr, fans=find_bends(scr, ci, byid, poles) if bends else [])
        _LIVE[key] = c
    return c


def live_power_map(lattice, screen_id, gaps, mesh=121, default_gap=5.0, bends=True):
    """Карта ПЛОТНОСТИ МОЩНОСТИ (быстрый режим, Льенар) для живого кокпита: зазор КАЖДОГО ондулятора — из gaps
    {gid: мм} (нет — default_gap), источник — по рамке железа (номинальная машина), веера диполей — как в
    footprint. Возврат — dict footprint (S [Вт/мм²], extent, peak, total, und_W, bend_W, n_bend, центроид) +
    illuminators [(gid, зазор)], src, D, kind; None — экран не построить."""
    c = live_context(lattice, screen_id, bends)
    if c is None:
        return None
    scr = c["scr"]; ems, ill = [], []
    for gid, insts in c["groups"].items():
        gap = float(gaps.get(gid, default_gap))
        em = make_emitters({gid: insts}, c["ci"], c["byid"], gap)
        if gid in em and find_illuminators(scr, em):
            ems.append(em[gid]); ill.append((gid, gap))
    fp = footprint(scr, ems, int(mesh), c["fans"])
    if scr.get("mode") == "bend":
        D = float(scr["D"]) if scr.get("D") else None
    else:
        srcP = emitter_frame(c["groups"][scr["src"]], c["ci"], c["byid"])[0] if scr["src"] in c["groups"] else None
        D = float(np.linalg.norm(scr["C"] - srcP)) if srcP is not None else None
    fp.update(illuminators=ill, src=scr["src"], D=D, kind=scr.get("mode", "axis"))
    return fp


_TRACE_CTX = {}               # решётка → контекст sr_trace (модули, цепочки, кэш траекторий)
_BEND_TRACES = {}             # (решётка, экран) → траектории цепочек, достающих экран (отсеянные)
_BEND_VIEW = {}               # (решётка, экран, сетка, величина) → карта


def _bend_traces(lattice, scr, cull=40.0):
    """Траектории цепочек диполей (мягкий край по картам поля), достающих экран, — отсеянные под экран."""
    import sr_trace as st
    key = (lat_key(lattice), scr["sid"])
    if key not in _BEND_TRACES:
        lk = lat_key(lattice)
        if lk not in _TRACE_CTX:
            _TRACE_CTX[lk] = st.make_context(lattice if is_machine(lattice) else int(lattice), None, None, 5.0)
        ctx = _TRACE_CTX[lk]
        trs = []
        for ic in st._chain_candidates(ctx, scr, cull):
            tr = st._chain_trace(ctx, ic, True)
            m = st.cull_mask(tr, scr, ctx["g"], cull)
            if m.any():
                trs.append(dict(tr, **{k: tr[k][m] for k in ("r", "t", "acc", "dt")}))
        _BEND_TRACES[key] = (ctx, trs)
    return _BEND_TRACES[key]


def bend_view(lattice, screen_id, mesh=121, quantity="power"):
    """Карта экрана у ДИПОЛЯ для кокпита — по траектории электрона (точный режим sr_trace: протяжённый источник,
    мягкий край по картам поля), без стыков вееров быстрого режима. От живого состояния не зависит — считается
    один раз на (экран, сетка, величина). «power» — Вт/мм² (footprint_accurate, эффекты extended+fringe);
    «flux» — ф/с/мрад²/0.1%BW на ε_c центрального модуля. Возврат dict(S, extent, peak, [total, n_bend, cu, cv],
    ec, dev, D) или None."""
    key = (lat_key(lattice), str(screen_id), int(mesh), str(quantity))
    if key in _BEND_VIEW:
        return _BEND_VIEW[key]
    import sr_trace as st
    c = live_context(lattice, screen_id)
    if c is None or c["scr"].get("mode") != "bend":
        return None
    scr = c["scr"]
    ctx, trs = _bend_traces(lattice, scr)
    ec0 = float(scr.get("ec") or (max(f["ec"] for f in c["fans"]) if c["fans"] else 0.0))
    if quantity == "power":
        fp = st.footprint_accurate(scr, ctx, "extended,fringe", int(mesh), undulators=False)
        res = dict(S=fp["S"], extent=fp["extent"], peak=fp["peak"], total=fp["total"], bend_W=fp["bend_W"],
                   und_W=0.0, n_bend=len(fp["chains"]), cu=fp["cu"], cv=fp["cv"])
    else:
        res = st.trace_flux_map(trs, scr, ec0, int(mesh))
    res.update(ec=ec0, dev=scr.get("dev", ""), D=scr.get("D"))
    _BEND_VIEW[key] = res
    return res


def bend_spectrum(lattice, screen_id, u_mm, v_mm, n=400):
    """Сплошной спектр в точке экрана у диполя по траектории (та же модель, что карта bend_view)."""
    import sr_trace as st
    c = live_context(lattice, screen_id)
    if c is None or c["scr"].get("mode") != "bend":
        return None
    ctx, trs = _bend_traces(lattice, c["scr"])
    return st.trace_point_spectrum(trs, c["scr"], u_mm, v_mm, n)


def save_png(scr, fp, path, gap, calc="fast"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    im = ax.imshow(fp["S"].T, origin="lower", extent=fp["extent"], aspect="equal", cmap="inferno")
    ax.plot(fp["cu"], fp["cv"], "+", color="cyan", ms=10, mew=1.5)
    ax.set_xlabel("u (горизонт.), мм"); ax.set_ylabel("v (верт.), мм")
    ax.set_title(f"{scr['sid']} ← {scr['src']}  (зазор {gap:g} мм, {'точный' if calc == 'accurate' else 'быстрый'})\n"
                 f"пик {fp['peak']:.1f} Вт/мм²,  перехвачено {fp['total']:.0f} Вт", fontsize=9)
    cb = fig.colorbar(im, ax=ax, shrink=0.85); cb.set_label("Вт/мм²", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=120); plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description="анализ фотонных экранов: автопоиск осветителей + мощность")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--ws", default=None, help="папка возмущённой машины (tables/lattice_configs.csv); иначе номинал")
    ap.add_argument("--gap", type=float, default=5.0, help="зазор всех ондуляторов, мм")
    ap.add_argument("--mesh", type=int, default=161, help="сетка экрана (точек на сторону)")
    ap.add_argument("--cone-cut", type=float, default=6.0, help="порог видимости в σ_θ конуса (только --model gauss)")
    ap.add_argument("--model", choices=["lienard", "gauss"], default="lienard",
                    help="угловая мощность ИД: точная (Льенар) или прежний круглый гаусс")
    ap.add_argument("--no-bends", action="store_true", help="не добавлять веера диполей")
    ap.add_argument("--source", choices=["orbit", "magnet"], default="orbit",
                    help="направление света: по замкнутой орбите электрона (физично) или по рамке железа (прежнее)")
    ap.add_argument("--correct", action="store_true", help="с --source orbit: сначала скорректировать орбиту (SVD)")
    ap.add_argument("--no-png", action="store_true")
    ap.add_argument("--calc", choices=["fast", "accurate"], default="fast",
                    help="fast — точечный эмиттер (как раньше); accurate — Льенар по траектории (sr_trace), медленно")
    ap.add_argument("--effects", default="all",
                    help="эффекты точного режима: all | none | список extended,poles,fringe,beam | all,-beam")
    ap.add_argument("--breakdown", action="store_true", help="точный режим: вклад каждого эффекта (счёт ×~3)")
    ap.add_argument("--ppp", type=int, default=64, help="точный режим: точек траектории на период ИД")
    ap.add_argument("--cull", type=float, default=40.0,
                    help="точный режим: отсев точек траектории, не достающих экран с запасом cull/γ (0 — без отсева)")
    ap.add_argument("--emit-x", type=float, default=sr.P.emit_x, help="точный режим: εx, м·рад (допущение)")
    ap.add_argument("--kappa", type=float, default=sr.P.kappa, help="точный режим: связь εy/εx (допущение)")
    ap.add_argument("--sigE", type=float, default=sr.P.sig_E, help="точный режим: σδ (допущение)")
    a = ap.parse_args()

    groups = load_groups(a.lattice)
    screens = load_screens(a.lattice)
    ci, byid = load_lattice(a.lattice, a.ws)
    B0, K, P_tot, Pd0, sig = power_params(a.gap)
    et = None
    if a.source == "orbit":
        import electron_source as es
        et = es.electron_track(a.lattice, a.ws, a.correct)
    emit = make_emitters(groups, ci, byid, a.gap, a.model, et)   # эмиттер на устройство
    poles = set().union(*groups.values()) if groups else set()

    ctx = None
    if a.calc == "accurate":
        import sr_trace as st
        st.parse_effects(a.effects)                              # проверка списка до долгого счёта
        ctx = st.make_context(a.lattice, a.ws, et, a.gap, ppp=a.ppp, ci=ci, byid=byid, groups=groups)

    os.makedirs(RESDIR, exist_ok=True)
    tag = "номинал" if not a.ws else os.path.basename(a.ws.rstrip("/\\"))
    print(f"решётка {a.lattice} [{tag}]: зазор {a.gap:g} мм → B0={B0:.3f} Т, K={K:.3f}; "
          f"устройство: P={P_tot:.0f} Вт, пик {Pd0/1e9:.1f} кВт/мрад²; модель ИД: {a.model}"
          + ("" if a.model == "lienard" else f" (σ_θ={sig*1e3:.3f} мрад)")
          + ("; веера диполей: выкл." if a.no_bends else "; веера диполей: вкл."))
    src_txt = ("по рамке железа" if et is None else
               f"по орбите электрона{' после коррекции' if et['corrected'] else ''} "
               f"(в рамках магнитов rms x={et['rms_x']*1e6:.2f}, y={et['rms_y']*1e6:.2f} мкм)")
    print(f"  источник света: {src_txt}")
    if ctx is not None:
        print(f"  режим: ТОЧНЫЙ (Льенар по траектории), эффекты {','.join(sorted(st.parse_effects(a.effects))) or '—'}; "
              f"εx={a.emit_x:.3g} м·рад, κ={a.kappa:g}, σδ={a.sigE:g} (допущения); ppp={a.ppp}, cull={a.cull:g}/γ")
    print(f"  экранов: {len(screens)}; эмиттеров-ондуляторов: {len(emit)}\n")
    hdr = ["lattice_id", "screen_id", "source_ref", "n_illum", "illuminators", "gap_mm",
           "peak_W_mm2", "total_W", "centroid_u_mm", "centroid_v_mm", "rms_u_mm", "rms_v_mm",
           "und_W", "bend_W", "n_bend", "model", "calc"]
    out_rows = []
    for scr in screens:
        found = find_illuminators(scr, emit, a.cone_cut)          # АВТОПОИСК: кто светит на этот экран
        bends = [] if a.no_bends else find_bends(scr, ci, byid, poles, et=et)
        emitters = [emit[g] for g in found]
        if not emitters and not bends and ctx is None:
            print(f"  {scr['sid']:12s} ← {scr['src']}: НЕТ осветителей (конус и веера мимо экрана)")
            out_rows.append([a.lattice, scr["sid"], scr["src"], 0, "", f"{a.gap:g}", 0, 0, 0, 0, 0, 0,
                             0, 0, 0, a.model, a.calc])
            continue
        fp = footprint(scr, emitters, a.mesh, bends)
        model = a.model
        if ctx is not None:                                      # быстрый — для сравнения, в CSV/PNG — точный
            fast = fp
            fp = st.footprint_accurate(scr, ctx, a.effects, a.mesh, a.cull, a.emit_x, a.kappa, a.sigE, a.breakdown)
            found = fp["illuminators"]; model = "lienard-trace"
        png = os.path.join(RESDIR, f"L{a.lattice}_{scr['sid']}{'' if not a.ws else '_pert'}"
                                   f"{'_acc' if ctx is not None else ''}.png")
        if not a.no_png:
            save_png(scr, fp, png, a.gap, a.calc)
        bend_txt = f" + {fp['n_bend']} диполь(ей) {fp['bend_W']:.0f} Вт" if fp["n_bend"] else ""
        print(f"  {scr['sid']:12s} ← {scr['src']}  осветит.: {','.join(found) or '—'}{bend_txt}  "
              f"пик {fp['peak']:7.1f} Вт/мм²  перехвач. {fp['total']:6.0f} Вт  "
              f"центроид ({fp['cu']:+.3f},{fp['cv']:+.3f}) мм  σ=({fp['ru']:.2f},{fp['rv']:.2f}) мм")
        if ctx is not None:
            print(f"      быстрый: пик {fast['peak']:7.2f} Вт/мм², перехв. {fast['total']:7.1f} Вт; " + st.compare_line(fast, fp))
            if fp.get("chains"):
                print(f"      цепочки диполей: {', '.join(fp['chains'])}")
            for line in st.breakdown_lines(fp):
                print(line)
        out_rows.append([a.lattice, scr["sid"], scr["src"], len(found), ",".join(found), f"{a.gap:g}",
                         f"{fp['peak']:.4f}", f"{fp['total']:.3f}", f"{fp['cu']:.5f}", f"{fp['cv']:.5f}",
                         f"{fp['ru']:.5f}", f"{fp['rv']:.5f}", f"{fp['und_W']:.3f}", f"{fp['bend_W']:.3f}",
                         fp["n_bend"], model, a.calc])

    rp = os.path.join(RESDIR, f"screen_readings_L{a.lattice}{'' if not a.ws else '_pert'}"
                              f"{'_acc' if ctx is not None else ''}.csv")
    with open(rp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(out_rows)
    print(f"\n[OK] {rp}")
    if not a.no_png:
        print(f"[OK] карты PNG → {RESDIR}")


if __name__ == "__main__":
    main()
