# -*- coding: utf-8 -*-
"""Режим 3 — HTML-отчёт с картой частот/ДИФФУЗИЕЙ для машины, которую оставил опыт с ОС.

mode3_demo.py показывал только КОНТУР ДА (PNG) трёх конфигураций (проектная / искажённая орбитой /
скорректированная). Здесь — полноценный интерактивный HTML поверх того же слоя анализа (control/analysis),
с РАЗДЕЛАМИ FMA/ДИФФУЗИИ: сплошная карта диффузии D=log10|Δν| в плоскости (x,y), диаграмма рабочей точки с
точками карты частот, ДА-тепловая карта и импульсная апертура — для КАЖДОЙ из трёх машин, плюс сводная
таблица (орбита, частоты, площадь ДА, медиана/хаос диффузии, импульсный аксептанс).

Физический смысл (как в mode3_demo): орбитное искажение НЕ двигает линейную частоту, но через feed-down
сквозь секступоли РАЗГОНЯЕТ ДИФФУЗИЮ и съедает ДА — а коррекция орбиты возвращает регулярность движения.
Карта диффузии показывает это нагляднее контура ДА: видно, КАК хаос заполняет апертуру и отступает.

Движок — продакшн facility_report (freq_map/dynamic_aperture/momentum_aperture/linear_optics на C++-ядре)
через control.analysis; фигуры — штатные fr._fig_* ; каркас — traj_analyzer._assemble. Новой физики нет.

ЗАПУСК: py -3.14 control/mode3_report.py [smoke] [--turns N] [--no-open]
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np

import facility_report as fr
import ring_analysis as ra
import traj_analyzer as ta

from control import analysis as an
from control.ring_model import RingModel

LAT, NSL = 26, 8
TARGET_ORBIT = 200e-6
RCOND = 1e-2
OUT = os.path.join("workspace", "results", "mode3_report.html")


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def orbit_rms(ring, mon):
    co = ra.closed_orbit(ring, dims=4)
    snap = np.asarray(ring.track_elements([co.tolist()]))[0]
    return rms(snap[mon, 0]), rms(snap[mon, 2])


def diffusion_stats(fma):
    """Медиана log10 D по выжившим, доля хаоса (D>−3) и доля потерянных (нефинитных) точек сетки."""
    d = np.asarray(fma["diff"], float)
    fin = np.isfinite(d)
    n = d.size
    if not np.any(fin):
        return dict(median=float("nan"), chaos=1.0, lost=1.0)
    return dict(median=float(np.median(d[fin])),
                chaos=float(np.mean(d[fin] > -3.0)),
                lost=float(np.mean(~fin)))


def titled(fig, prefix):
    """Префикс имени конфигурации к заголовку штатной фигуры facility_report."""
    try:
        cur = fig.layout.title.text or ""
    except Exception:
        cur = ""
    fig.update_layout(title=dict(text=f"<b>{prefix}</b> — {cur}", x=0.5, y=0.97, font=dict(size=14)))
    return fig


def build_configs():
    """Три конфигурации машины киками линейной модели (как plant.kick_lists): проектная, искажённая
    орбитой (feed-down квадруполей), скорректированная (SVD-коррекция орбиты)."""
    rm = RingModel(lattice=LAT, nsl=NSL)
    corr, quads, K1L, mon = rm.corr, rm.quads, np.asarray(rm.K1L), rm.mon
    Hc = rm.response(rm.hk_ids, ["ebpm.x"])["ebpm.x"]
    Gx = rm.response(rm.qx_ids, ["ebpm.x"])["ebpm.x"]
    rm.close()

    rp = np.random.default_rng(11); pat = rp.normal(0, 1, len(quads)); pat /= np.std(pat)
    amp = 1e-6
    amp *= TARGET_ORBIT / max(rms(Gx @ (amp * pat)), 1e-30)
    dxpat = amp * pat
    u = -np.linalg.pinv(Hc, rcond=RCOND) @ (Gx @ dxpat)
    dist_k = [(int(quads[j]), float(K1L[j] * dxpat[j])) for j in range(len(quads))]
    corr_k = [(int(corr[j]), float(u[j])) for j in range(len(corr))]
    return mon, {"проектная": [], "искажённая": dist_k, "скорректированная": dist_k + corr_k}


def main():
    ap = argparse.ArgumentParser(description="режим 3 — HTML-отчёт FMA/диффузии для 3 машин")
    ap.add_argument("mode", nargs="?", default="", help="smoke — быстрый прогон")
    ap.add_argument("--turns", type=int, default=None, help="оборотов в сканах ДА/карты частот")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()
    smoke = a.mode == "smoke"
    turns = a.turns if a.turns else (128 if smoke else 256)

    print("=== РЕЖИМ 3 (HTML): карта частот/диффузия проектной / искажённой / скорректированной машины ===")
    print(f"решётка {LAT}, {turns} оборотов; искажение орбиты ~{TARGET_ORBIT*1e6:.0f} мкм (feed-down)\n", flush=True)

    mon, configs = build_configs()
    inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
    ring = fr._ring(inf, nsl=NSL)

    results = {}
    try:
        for name, kk in configs.items():
            print(f"  анализ: {name} …", flush=True)
            ring.apply_perturb(kicks=kk)
            ox, oy = orbit_rms(ring, mon)
            r = an.analyze(ring, inf, turns=turns, do=("optics", "da", "fma", "mom"))
            results[name] = dict(res=r, ox=ox, oy=oy, dm=an.da_metrics(r["da"]), ds=diffusion_stats(r["fma"]))
    finally:
        ring.apply_perturb(); ring.close()

    # --- сводная таблица
    hdr = ["Конфигурация", "Орбита X/Y, мкм", "νx/νy", "ДА площадь, мм²", "ДА x±, мм",
           "медиана log₁₀D", "хаос (D>−3)", "импульс ±δ, %"]
    rows = []
    for name, R in results.items():
        r = R["res"]; dm = R["dm"]; ds = R["ds"]; mom = r["mom"]
        rows.append([name, f"{R['ox']*1e6:.1f} / {R['oy']*1e6:.1f}",
                     f"{r['tunes'][0]:.4f} / {r['tunes'][1]:.4f}",
                     fr._num(dm["area_mm2"], 4), f"{dm['x_pos']:.2f} / {dm['x_neg']:.2f}",
                     fr._num(ds["median"], 3), f"{ds['chaos']*100:.0f}%",
                     f"+{mom['acc_pos']*100:.2f} / {mom['acc_neg']*100:.2f}"])
    figs = [(fr._table(hdr, rows, [16, 13, 12, 12, 11, 11, 8, 12]), ta._table_h(len(rows) + 1))]

    # --- по разделам: диффузия, рабочая точка, ДА, импульс для каждой конфигурации
    order = list(results.keys())
    for name in order:
        figs.append((titled(fr._fig_fma_space(results[name]["res"]["fma"]), name), 500))
    for name in order:
        figs.append((titled(fr._fig_tunes(results[name]["res"]["optics"], None,
                                          results[name]["res"]["fma"]), name), 560))
    for name in order:
        figs.append((titled(fr._fig_da(results[name]["res"]["da"]), name), 480))
    for name in order:
        figs.append((titled(fr._fig_momentum(results[name]["res"]["mom"],
                                             inf.get("rf") and None), name), 380))

    ad = results["проектная"]["dm"]["area_mm2"]
    ai = results["искажённая"]["dm"]["area_mm2"]
    ac = results["скорректированная"]["dm"]["area_mm2"]
    ci = results["искажённая"]["ds"]["chaos"]; cc = results["скорректированная"]["ds"]["chaos"]
    title = f"SynchrotronSim — режим 3: карта частот / диффузия (решётка {LAT})"
    subtitle = (f"3 машины · {turns} об · искажение орбиты ~{TARGET_ORBIT*1e6:.0f} мкм (feed-down квадруполей) · "
                f"ДА проект {ad:.1f} → искажённая {ai:.1f} ({ai/ad*100:.0f}%) → скорр. {ac:.1f} ({ac/ad*100:.0f}%) мм² · "
                f"хаос искажённая {ci*100:.0f}% → скорр. {cc*100:.0f}%")
    html = ta._assemble(LAT, 0, title, subtitle, figs)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"\nплощадь ДА: проект {ad:.1f} → искажённая {ai:.1f} ({ai/ad*100:.0f}%) → скорр. {ac:.1f} ({ac/ad*100:.0f}%) мм²")
    print(f"хаос (D>−3): искажённая {ci*100:.0f}% → скорректированная {cc*100:.0f}%")
    print(f"[готово] HTML: {a.out}")
    if not a.no_open:
        try:
            import webbrowser
            webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception:
            pass


if __name__ == "__main__":
    main()
