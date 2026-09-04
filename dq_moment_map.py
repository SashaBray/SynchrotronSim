# -*- coding: utf-8 -*-
"""Интерактивная карта моментов DQ: диполь и квадруполь как поверхности над плоскостью токов (I1, I2).

Считает моменты сетевой моделью насыщения (`branch_network`) на сетке токов при I3 = 0, переводит
гармоники в физические коэффициенты `mat_B` и `mat_G` по калибровке на снятых картах и собирает
самодостаточный HTML (plotly офлайн) с двумя поверхностями. На поверхностях отмечены точки, где
реально сняты карты поля, и проектные точки, посчитанные ЭТОЙ ЖЕ моделью (а не линейной из базы).

    py -3.14 dq_moment_map.py
"""
from __future__ import annotations

import argparse
import os

import numpy as np

import branch_network as bn
import saturation_model as sm

OUT = os.path.join("workspace", "results", "dq_moment_map.html")
TURNS = {"I1": 53, "I2": 13, "I3": 90}
N_GRID = 55


def calibrate(order, key):
    """Множитель перехода «гармоника → физкоэффициент»: mat_* = k · |C_n|, по снятым картам."""
    d = sm.load_dataset("DQ")
    harm = {t[3]: t[1] for t in bn.read_maps("DQ", order)}
    ks = []
    for r in d["records"]:
        if not r["coeffs"] or key not in r["coeffs"] or r["map_id"] not in harm:
            continue
        c = abs(harm[r["map_id"]])
        if c > 0:
            ks.append(abs(r["coeffs"][key]) / c)     # ЗНАК несёт не отношение, а обмотка:
    ks = np.array(ks)                                # у I1 градиент одного знака, у I2/I3 другого
    return float(ks.mean()), float(ks.std() / ks.mean() * 100)


def fit_order(order):
    """Веса модели под гармонику данного порядка."""
    cir, p, w = bn.fit("DQ", order, verbose=False)
    return cir, p, w


def signed(cir, p, w, ref, cur):
    """Момент СО ЗНАКОМ: проекция на направление, которое задаёт обмотка I1.

    Модель даёт комплексную гармонику; её фаза и несёт знак момента. Проекция на опорное направление
    возвращает знак естественно — у I1 градиент одного знака, у I2/I3 противоположного."""
    v = bn.predict(cir, cur, p, w)
    return float(np.real(v * np.conj(ref))) / abs(ref)


def surface(cir, p, w, k, ref, i1_a, i2_a):
    """Момент на сетке токов (амперы → ампер-витки через число витков)."""
    Z = np.zeros((len(i2_a), len(i1_a)))
    for a, i2 in enumerate(i2_a):
        for b, i1 in enumerate(i1_a):
            cur = {"I1": i1 * TURNS["I1"], "I2": i2 * TURNS["I2"]}
            Z[a, b] = k * signed(cir, p, w, ref, cur)
    return Z


def measured_points(order, k, ref):
    """Точки снятых карт: (I1, I2 в амперах, измеренный момент со знаком)."""
    pts = []
    for cur, meas, sp, mid in bn.read_maps("DQ", order):
        if cur.get("I3"):                       # на карте нужен I3 = 0
            continue
        i1 = cur.get("I1", 0.0) / TURNS["I1"]
        i2 = cur.get("I2", 0.0) / TURNS["I2"]
        pts.append((i1, i2, k * float(np.real(meas * np.conj(ref))) / abs(ref), mid))
    return pts


def design_points(models, ks):
    """Проектные точки — решение СИСТЕМЫ из двух уравнений (диполь И квадруполь) по I1 и I2.

    Раньше на каждой поверхности отмечалась своя независимая точка, и они не совпадали: одной
    обмоткой обе величины не получить. Проектный градиент от I1 при токе, дающем проектный изгиб,
    выходит выше нужного, и вторая обмотка его гасит — это и есть смысл совмещённой функции."""
    from scipy.optimize import fsolve
    import element_compare as ec
    (c1, p1, w1, r1), (c2, p2, w2, r2) = models
    k1, k2 = ks
    lc = {r["device_type_id"]: r for r in ec.mt._read("lattice_configs.csv")
          if str(r.get("lattice_id")) == "1"}
    out = []
    for dev, guess in (("DQ_I1_1A_L1000", (150.0, 30.0)), ("DQ_I1_1A_L400", (95.0, -60.0))):
        row = lc.get(dev, {})
        B, G = (row.get("mat_B") or "").strip(), (row.get("mat_G") or "").strip()
        if not B or not G:
            continue
        B, G = float(B), float(G)

        def eqs(x):
            cur = {"I1": x[0] * TURNS["I1"], "I2": x[1] * TURNS["I2"]}
            return [k1 * signed(c1, p1, w1, r1, cur) - B,
                    k2 * signed(c2, p2, w2, r2, cur) - G]

        sol, _, ier, _ = fsolve(eqs, guess, full_output=True)
        if ier != 1:
            continue
        cur = {"I1": sol[0] * TURNS["I1"], "I2": sol[1] * TURNS["I2"]}
        out.append((sol[0], sol[1],
                    k1 * signed(c1, p1, w1, r1, cur),
                    k2 * signed(c2, p2, w2, r2, cur),
                    dev.replace("DQ_I1_1A_", "")))
    return out


def build(path=OUT):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    k1, s1 = calibrate(1, "B")
    k2, s2 = calibrate(2, "G")
    print(f"калибровка: mat_B = {k1:.5g}·|C1| (разброс {s1:.1f}%), mat_G = {k2:.5g}·|C2| ({s2:.1f}%)")

    i1_a = np.linspace(0.0, 300.0, N_GRID)          # ампер, с учётом витков
    i2_a = np.linspace(-140.0, 140.0, N_GRID)

    models = []
    for order in (1, 2):
        cir, p, w = fit_order(order)
        models.append((cir, p, w, bn.predict(cir, {"I1": 5000.0}, p, w)))
    dp = design_points(models, (k1, k2))                # ОДНА точка на оба момента

    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "surface"}, {"type": "surface"}]],
                        subplot_titles=("Дипольный момент mat_B, Тл",
                                        "Квадрупольный момент mat_G, Тл/м"))
    for col, (order, key, k, name) in enumerate(
            ((1, "B", k1, "диполь"), (2, "G", k2, "квадруполь")), start=1):
        cir, p, w, ref = models[order - 1]
        Z = surface(cir, p, w, k, ref, i1_a, i2_a)
        fig.add_trace(go.Surface(x=i1_a, y=i2_a, z=Z, colorscale="Viridis", showscale=False,
                                 opacity=0.92, name=name,
                                 hovertemplate="I1 %{x:.1f} А<br>I2 %{y:.1f} А<br>%{z:.5g}<extra></extra>"),
                      row=1, col=col)
        mp = measured_points(order, k, ref)
        if mp:
            fig.add_trace(go.Scatter3d(
                x=[q[0] for q in mp], y=[q[1] for q in mp], z=[q[2] for q in mp],
                mode="markers", marker=dict(size=5, color="#c0392b", symbol="circle"),
                name="снятые карты", showlegend=(col == 1),
                text=[q[3] for q in mp], hovertemplate="%{text}<br>I1 %{x:.1f} А, I2 %{y:.1f} А"
                                                       "<br>измерено %{z:.5g}<extra></extra>"),
                row=1, col=col)
        if dp:
            zi = 2 if key == "B" else 3
            fig.add_trace(go.Scatter3d(
                x=[q[0] for q in dp], y=[q[1] for q in dp], z=[q[zi] for q in dp],
                mode="markers+text", marker=dict(size=9, color="#e8a33d", symbol="diamond"),
                text=[q[4] for q in dp], textposition="top center",
                name="проектная точка (оба момента сразу)", showlegend=(col == 1),
                hovertemplate="%{text}<br>I1 %{x:.1f} А, I2 %{y:.1f} А"
                              "<br>%{z:.5g}<extra></extra>"),
                row=1, col=col)
        fig.update_scenes(dict(xaxis_title="I1, А", yaxis_title="I2, А",
                               zaxis_title=f"mat_{key}"), row=1, col=col)
    fig.update_layout(
        title="DQ: моменты как функция токов обмоток (I3 = 0). Поверхность — сетевая модель "
              "насыщения; красные точки — снятые карты, жёлтые ромбы — проектные точки",
        height=760, margin=dict(l=0, r=0, t=90, b=0))
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fig.write_html(path, include_plotlyjs=True, full_html=True)
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args(argv)
    print("карта ->", build(a.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
