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

ЗАПУСК
  py -3.14 undulators/screen_analysis.py --lattice 27 --gap 5
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


def load_groups(lattice):
    """{group_id: set(instance_ids)} для решётки."""
    out = {}
    for r in csv.DictReader(open(EG, encoding="utf-8-sig"), delimiter=";"):
        if r.get("lattice_id") == str(lattice):
            out[r["group_id"]] = set(r["instance_ids"].split("|"))
    return out


def load_screens(lattice):
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


def save_png(scr, fp, path, gap):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    im = ax.imshow(fp["S"].T, origin="lower", extent=fp["extent"], aspect="equal", cmap="inferno")
    ax.plot(fp["cu"], fp["cv"], "+", color="cyan", ms=10, mew=1.5)
    ax.set_xlabel("u (горизонт.), мм"); ax.set_ylabel("v (верт.), мм")
    ax.set_title(f"{scr['sid']} ← {scr['src']}  (зазор {gap:g} мм)\n"
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
    print(f"  экранов: {len(screens)}; эмиттеров-ондуляторов: {len(emit)}\n")
    hdr = ["lattice_id", "screen_id", "source_ref", "n_illum", "illuminators", "gap_mm",
           "peak_W_mm2", "total_W", "centroid_u_mm", "centroid_v_mm", "rms_u_mm", "rms_v_mm",
           "und_W", "bend_W", "n_bend", "model"]
    out_rows = []
    for scr in screens:
        found = find_illuminators(scr, emit, a.cone_cut)          # АВТОПОИСК: кто светит на этот экран
        bends = [] if a.no_bends else find_bends(scr, ci, byid, poles, et=et)
        emitters = [emit[g] for g in found]
        if not emitters and not bends:
            print(f"  {scr['sid']:12s} ← {scr['src']}: НЕТ осветителей (конус и веера мимо экрана)")
            out_rows.append([a.lattice, scr["sid"], scr["src"], 0, "", f"{a.gap:g}", 0, 0, 0, 0, 0, 0,
                             0, 0, 0, a.model])
            continue
        fp = footprint(scr, emitters, a.mesh, bends)
        png = os.path.join(RESDIR, f"L{a.lattice}_{scr['sid']}{'' if not a.ws else '_pert'}.png")
        if not a.no_png:
            save_png(scr, fp, png, a.gap)
        bend_txt = f" + {fp['n_bend']} диполь(ей) {fp['bend_W']:.0f} Вт" if fp["n_bend"] else ""
        print(f"  {scr['sid']:12s} ← {scr['src']}  осветит.: {','.join(found) or '—'}{bend_txt}  "
              f"пик {fp['peak']:7.1f} Вт/мм²  перехвач. {fp['total']:6.0f} Вт  "
              f"центроид ({fp['cu']:+.3f},{fp['cv']:+.3f}) мм  σ=({fp['ru']:.2f},{fp['rv']:.2f}) мм")
        out_rows.append([a.lattice, scr["sid"], scr["src"], len(found), ",".join(found), f"{a.gap:g}",
                         f"{fp['peak']:.4f}", f"{fp['total']:.3f}", f"{fp['cu']:.5f}", f"{fp['cv']:.5f}",
                         f"{fp['ru']:.5f}", f"{fp['rv']:.5f}", f"{fp['und_W']:.3f}", f"{fp['bend_W']:.3f}",
                         fp["n_bend"], a.model])

    rp = os.path.join(RESDIR, f"screen_readings_L{a.lattice}{'' if not a.ws else '_pert'}.csv")
    with open(rp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(out_rows)
    print(f"\n[OK] {rp}")
    if not a.no_png:
        print(f"[OK] карты PNG → {RESDIR}")


if __name__ == "__main__":
    main()
