# -*- coding: utf-8 -*-
"""Слой анализа машины (режим 3): динамическая апертура, карта частот/диффузия, линейная оптика —
поверх ЛЮБОГО сконфигурированного кольца, единым стандартным результатом.

Платформенная идея: 'охарактеризуй вот ЭТУ машину' — ту же, что оставил опыт с ОС (искажённую или
скорректированную), или любую конфигурацию. Машина задаётся набором киков apply_perturb (корректоры +
feed-down возмущения) — ровно то, что отдаёт plant.kick_lists(state). Никакой новой физики: внутри —
продакшн-расчёты facility_report (dynamic_aperture/freq_map/momentum_aperture/linear_optics на ядре
C++), сведённые к одному вызову и стандартному словарю результата + метрикам и фигуре.

ЗАПУСК как библиотека; демо — control/mode3_demo.py.
"""
import os

os.environ.setdefault("SS_SECTOR", "1")

import numpy as np

import facility_report as fr

SECTIONS = ("optics", "da", "fma", "mom")


def analyze(ring, inf, turns=512, do=("optics", "da"), fast=False):
    """Стандартный результат анализа сконфигурированного кольца `ring`. do — какие разделы считать:
    optics (твисс/частоты), da (динамическая апертура), fma (карта частот/диффузия), mom (импульсная
    апертура). fast=True — только аналитика без трекинга (optics)."""
    res = {"turns": turns, "label": inf.get("label", ""), "sections": tuple(do)}
    if "optics" in do:
        lin = fr.linear_optics(ring, inf)
        res["optics"] = lin
        res["tunes"] = (lin["nux"], lin["nuy"])
    track = set(do) & {"da", "fma", "mom"}
    if fast or not track:
        return res
    co, xe, ye, dpos, dneg = fr.auto_ranges(ring, inf, turns)
    res["ranges"] = dict(co=co, xe=xe, ye=ye, dpos=dpos, dneg=dneg)
    if "da" in do:
        res["da"] = fr.dynamic_aperture(ring, inf, turns, co, xe, ye)
    if "fma" in do:
        lin = res.get("optics") or fr.linear_optics(ring, inf)
        res["fma"] = fr.freq_map(ring, inf, lin, turns, co, xe, ye)
    if "mom" in do:
        res["mom"] = fr.momentum_aperture(ring, inf, turns, co, xe, dpos, dneg)
    return res


def analyze_config(lattice=26, nsl=8, kicks=(), vkicks=(), turns=512, do=("optics", "da"), fast=False):
    """Удобство: построить 4D-кольцо решётки, применить набор киков (конфигурацию машины) и анализировать.
    Для нескольких конфигураций подряд дешевле держать одно кольцо и звать analyze() с apply_perturb между."""
    inf = dict(fr.info_of(lattice)); inf["_id"] = lattice
    ring = fr._ring(inf, nsl=nsl)
    try:
        ring.apply_perturb(kicks=list(kicks), vkicks=list(vkicks))
        return analyze(ring, inf, turns=turns, do=do, fast=fast)
    finally:
        ring.apply_perturb(); ring.close()


def da_metrics(da):
    """Скаляры ДА: гор. границы (+/−, мм), макс. верт. (мм), площадь контура выживания (мм²)."""
    ys = da["ys"] * 1e3; ep = da["edge_pos"] * 1e3; en = da["edge_neg"] * 1e3
    live = (ep - en) > 1e-9
    area = float(np.trapezoid(ep - en, ys)) if len(ys) > 1 else 0.0   # ∫(x⁺−x⁻) dy, мм²
    return dict(x_pos=float(np.max(ep)) if ep.size else 0.0,
                x_neg=float(np.min(en)) if en.size else 0.0,
                y_max=float(np.max(np.abs(ys[live]))) if np.any(live) else 0.0,
                area_mm2=area)


def save_da_png(da_by_name, path, title="Динамическая апертура (контур выживания)"):
    """Наложить контуры ДА нескольких машин (имя→da) на одну (x, y)-плоскость."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"(matplotlib недоступен: {e} — фигура ДА пропущена)")
        return False
    fig, ax = plt.subplots(figsize=(6.8, 5.4))
    for name, da in da_by_name.items():
        ys = da["ys"] * 1e3; ep = da["edge_pos"] * 1e3; en = da["edge_neg"] * 1e3
        xc = np.concatenate([ep, en[::-1], ep[:1]]); yc = np.concatenate([ys, ys[::-1], ys[:1]])
        ax.plot(xc, yc, lw=1.9, label=name)
    ax.axhline(0, color="0.8", lw=0.6); ax.axvline(0, color="0.8", lw=0.6)
    ax.set_xlabel("x, мм"); ax.set_ylabel("y, мм"); ax.grid(alpha=0.3); ax.legend()
    ax.set_title(title)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)
    return True
