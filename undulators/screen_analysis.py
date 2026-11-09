# -*- coding: utf-8 -*-
"""Анализ ФОТОННЫХ ЭКРАНОВ — ЭТАП A: геометрия + мощность + АВТОПОИСК осветителей (без интерференции).

Для решётки (НОМИНАЛЬНОЙ или ВОЗМУЩЁННОЙ через --ws) и замороженных экранов (screens.csv):
  1) АВТОПОИСК. Какие ондуляторы РЕАЛЬНО светят на экран — тест видимости по конусу от ФАКТИЧЕСКОЙ
     геометрии источника (из lattice_configs ws). Привязка экрана к источнику — лишь стартовое зерно;
     экран берёт всех, чей конус его достаёт.
  2) МОЩНОСТЬ. Проекция угловой плотности мощности каждого осветителя на плоскость экрана и
     НЕКОГЕРЕНТНАЯ сумма по интенсивности → карта Вт/мм². (Гармоник тут нет — они на этапе B,
     когерентный интеграл поля по траектории.)
  3) ИТОГ. Пик Вт/мм², перехваченная мощность, смещение/размер центроида. Перекос источника (--ws)
     двигает ПЯТНО по неподвижному экрану — видно в центроиде.

МОДЕЛЬ ИЗЛУЧЕНИЯ (этап A). На устройство — ОДИН эмиттер (все полюса агрегированы): симметричный гаусс
по углу, НОРМИРОВАННЫЙ на пик (плотность на оси power_density_onaxis) и полную мощность total_power ⇒
σ_θ=√(P/2π·P'). Реальная картина слегка эллиптична (гориз. развёртка ±K/γ) — уточняется этапом B.
Зазор всех устройств задаётся --gap (как при вставке insert_undulators); позже — по инстансу.

ЗАПУСК
  py -3.14 undulators/screen_analysis.py --lattice 27 --gap 5
  py -3.14 undulators/screen_analysis.py --lattice 27 --ws <папка_возмущённой_машины>   # увод при перекосе
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
RADIATE = {"dipole", "dq"}                           # что светит (квадруполь-фокус UFC — нет)


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
    cs, axes = [], []
    for iid in insts:
        r = byid.get(iid)
        if r is None or r[ci["kind"]].strip() not in RADIATE:
            continue
        C = np.array([float(r[ci[a]] or 0) for a in ("x", "y", "z")])
        ypr = (float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))
        ax, _, _ = mt.frame(*ypr)
        cs.append(C); axes.append(np.asarray(ax, float))
    if not cs:
        return None, None
    return np.mean(cs, axis=0), _unit(np.sum(axes, axis=0))


def power_params(gap):
    B0 = sr.B0_of_gap(gap); K = sr.K_of_B0(B0)
    P_tot = sr.total_power_kW(B0) * 1e3                       # Вт
    Pd0 = sr.power_density_onaxis(B0, K) * 1e9               # кВт/мрад² → Вт/рад²
    sig = math.sqrt(P_tot / (2 * math.pi * Pd0))             # рад, эффективный конус
    return B0, K, P_tot, Pd0, sig


def footprint(scr, emitters, mesh):
    """Карта плотности мощности S [Вт/м²] на плоскости экрана + сводка."""
    hw, hh = scr["w"] / 2e3, scr["h"] / 2e3                  # полуразмеры, м
    su = np.linspace(-hw, hw, mesh); sv = np.linspace(-hh, hh, mesh)
    UU, VV = np.meshgrid(su, sv, indexing="ij")
    cell = (2 * hw / (mesh - 1)) * (2 * hh / (mesh - 1))     # площадь ячейки, м²
    pts = (scr["C"][None, None, :] + UU[..., None] * scr["u"][None, None, :]
           + VV[..., None] * scr["v"][None, None, :])        # (m,m,3) глобальные точки
    S = np.zeros((mesh, mesh))
    for (Psrc, axis, Pd0, sig) in emitters:
        d = pts - Psrc[None, None, :]
        D = np.linalg.norm(d, axis=2)
        dhat = d / D[..., None]
        costh = np.clip((dhat * axis[None, None, :]).sum(2), -1.0, 1.0)
        theta = np.arccos(costh)
        dPdO = Pd0 * np.exp(-theta ** 2 / (2 * sig ** 2))    # Вт/рад²
        cos_inc = np.abs((dhat * scr["n"][None, None, :]).sum(2))
        S += np.where(costh > 0, dPdO * cos_inc / D ** 2, 0.0)
    tot = float(S.sum() * cell)                              # Вт (перехвачено)
    wsum = float(S.sum()) or 1.0
    cu = float((S * UU).sum() / wsum) * 1e3                  # центроид, мм
    cv = float((S * VV).sum() / wsum) * 1e3
    ru = math.sqrt(max(0.0, float((S * (UU * 1e3 - cu) ** 2).sum() / wsum)))
    rv = math.sqrt(max(0.0, float((S * (VV * 1e3 - cv) ** 2).sum() / wsum)))
    return dict(S=S * 1e-6, extent=(-hw * 1e3, hw * 1e3, -hh * 1e3, hh * 1e3),  # S в Вт/мм²
                peak=float(S.max() * 1e-6), total=tot, cu=cu, cv=cv, ru=ru, rv=rv)


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
    ap.add_argument("--cone-cut", type=float, default=6.0, help="порог видимости в σ_θ конуса")
    ap.add_argument("--no-png", action="store_true")
    a = ap.parse_args()

    groups = load_groups(a.lattice)
    screens = load_screens(a.lattice)
    ci, byid = load_lattice(a.lattice, a.ws)
    B0, K, P_tot, Pd0, sig = power_params(a.gap)

    # эмиттер на устройство (из ФАКТической геометрии)
    emit = {}
    for gid, insts in groups.items():
        P, axis = emitter(insts, ci, byid)
        if P is not None:
            emit[gid] = (P, axis, Pd0, sig)

    os.makedirs(RESDIR, exist_ok=True)
    tag = "номинал" if not a.ws else os.path.basename(a.ws.rstrip("/\\"))
    print(f"решётка {a.lattice} [{tag}]: зазор {a.gap:g} мм → B0={B0:.3f} Т, K={K:.3f}; "
          f"устройство: P={P_tot:.0f} Вт, пик {Pd0/1e9:.1f} кВт/мрад², конус σ_θ={sig*1e3:.3f} мрад")
    print(f"  экранов: {len(screens)}; эмиттеров-ондуляторов: {len(emit)}\n")
    hdr = ["lattice_id", "screen_id", "source_ref", "n_illum", "illuminators", "gap_mm",
           "peak_W_mm2", "total_W", "centroid_u_mm", "centroid_v_mm", "rms_u_mm", "rms_v_mm"]
    out_rows = []
    for scr in screens:
        # АВТОПОИСК: кто светит на этот экран
        found = []
        for gid, (P, axis, pd0, sg) in emit.items():
            vec = scr["C"] - P; D = float(np.linalg.norm(vec)); dh = vec / D
            if float(np.dot(dh, axis)) <= 0:
                continue                                     # экран позади источника
            theta = math.acos(max(-1.0, min(1.0, float(np.dot(dh, axis)))))
            half_ang = 0.5 * math.hypot(scr["w"], scr["h"]) / 1e3 / D
            if theta <= a.cone_cut * sg + half_ang:
                found.append(gid)
        emitters = [emit[g] for g in found]
        if not emitters:
            print(f"  {scr['sid']:12s} ← {scr['src']}: НЕТ осветителей (конус мимо экрана)")
            out_rows.append([a.lattice, scr["sid"], scr["src"], 0, "", f"{a.gap:g}", 0, 0, 0, 0, 0, 0])
            continue
        fp = footprint(scr, emitters, a.mesh)
        png = os.path.join(RESDIR, f"L{a.lattice}_{scr['sid']}{'' if not a.ws else '_pert'}.png")
        if not a.no_png:
            save_png(scr, fp, png, a.gap)
        print(f"  {scr['sid']:12s} ← {scr['src']}  осветит.: {','.join(found)}  "
              f"пик {fp['peak']:7.1f} Вт/мм²  перехвач. {fp['total']:6.0f} Вт  "
              f"центроид ({fp['cu']:+.3f},{fp['cv']:+.3f}) мм  σ=({fp['ru']:.2f},{fp['rv']:.2f}) мм")
        out_rows.append([a.lattice, scr["sid"], scr["src"], len(found), ",".join(found), f"{a.gap:g}",
                         f"{fp['peak']:.4f}", f"{fp['total']:.3f}", f"{fp['cu']:.5f}", f"{fp['cv']:.5f}",
                         f"{fp['ru']:.5f}", f"{fp['rv']:.5f}"])

    rp = os.path.join(RESDIR, f"screen_readings_L{a.lattice}{'' if not a.ws else '_pert'}.csv")
    with open(rp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(out_rows)
    print(f"\n[OK] {rp}")
    if not a.no_png:
        print(f"[OK] карты PNG → {RESDIR}")


if __name__ == "__main__":
    main()
