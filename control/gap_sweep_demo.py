# -*- coding: utf-8 -*-
"""Демо ПРИВОДА ЗАЗОРА с ограничением скорости: оператор командует уставку, зазор едет к ней ≤ slew.

Показывает две вещи сразу:
  • ПРИВОД — реализованный зазор НЕ телепортируется: `command(15 мм)`, привод едет со скоростью 0.5 мм/с
    (как Limits.slew у корректоров), клип в ОДЗ; зазор(t) — рампа, а не скачок;
  • ОТКЛИК — `ScreenRecorder` читает РЕАЛИЗОВАННЫЙ зазор своего ондулятора → картина и E₁ едут за ним.

Квазистатически (без кольца): привод челюстей медленный (секунды) против динамики пучка (мкс), поэтому
зазор разумно гонять на операторском масштабе. В реальном прогоне plant зовёт `advance(dt)` каждый такт —
здесь advance по сетке времени делает то же. На выходе — HTML: рампа зазора(t), E₁(t) и картина при
начальном/конечном зазоре.

ЗАПУСК: py -3.14 control/gap_sweep_demo.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

from control import ExperimentState, Schema, UndulatorConfig
from control.recorders import ScreenRecorder

SCREEN = "SCR_U19"            # ондуляторный экран (источник und_19)
SOURCE = "und_19"


def main():
    smoke = "smoke" in sys.argv
    slew = 0.5                                   # мм/с — скорость привода челюстей (ПРЕДПОЛОЖЕНИЕ)
    dt_s = 2.0                                   # шаг времени, с
    n_win = 7 if smoke else 13                   # окон (0..24 с) — хватит на ход 5→15 мм при 0.5 мм/с

    st = ExperimentState(Schema(["a"]))
    st.undulators[SOURCE] = UndulatorConfig(5.0, name=SOURCE, slew_mm_s=slew)  # старт 5 мм
    st.undulators[SOURCE].command(15.0)          # ← оператор КОМАНДУЕТ уставку 15 мм (привод поедет)
    rng = np.random.default_rng(0)
    rec = ScreenRecorder(params=SCREEN, cadence="final")

    print(f"=== ДЕМО ПРИВОДА ЗАЗОРА: {SCREEN} ({SOURCE}), уставка 15 мм при {slew} мм/с ===")
    for w in range(n_win):
        st.undulators[SOURCE].advance(dt_s)      # привод продвигает реализованный зазор к уставке (как plant.step)
        st.observables["photon.x"] = np.array([rng.normal(0, 1e-7)])   # блуждание пятна (джиттер) для размытия
        st.observables["photon.y"] = np.array([rng.normal(0, 1e-7)])
        rec.sample(st, w, None)
        u = st.undulators[SOURCE]
        print(f"  t={w * dt_s:4.0f} с: зазор={u.gap_mm:5.2f} мм (уставка {u.gap_target:.0f})  K={u.K():.3f}", flush=True)
    rec.finalize(st, None)

    art = rec.artifact()[1]
    e1s, e1e = float(art["energy_start_keV"][0]), float(art["energy_end_keV"][0])
    print(f"\nзазор {art['gap_start_mm'][0]:.1f}→{art['gap_end_mm'][0]:.1f} мм  →  E₁ {e1s:.2f}→{e1e:.2f} кэВ")
    print(f"[OK] {_render(art, dt_s)}")


def _render(art, dt_s):
    import plotly.graph_objects as go
    from plotly.io import to_html
    from plotly.subplots import make_subplots

    ser = art["gap_series"]                              # (окно, зазор, E₁)
    t = ser[:, 0] * dt_s if ser.size else np.array([])
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        f"картина @ {art['gap_start_mm'][0]:.1f} мм (E₁={art['energy_start_keV'][0]:.1f} кэВ)",
        f"картина @ {art['gap_end_mm'][0]:.1f} мм (E₁={art['energy_end_keV'][0]:.1f} кэВ)",
        "привод: зазор(t) с ограничением скорости"))
    for col, key in ((1, "ideal_start"), (2, "ideal_end")):
        fig.add_trace(go.Heatmap(z=art[key], colorscale="Inferno", showscale=False), row=1, col=col)
        fig.update_xaxes(title_text="мм", row=1, col=col); fig.update_yaxes(title_text="мм", row=1, col=col)
    if ser.size:
        fig.add_trace(go.Scatter(x=t, y=ser[:, 1], mode="lines+markers", name="зазор",
                                 line=dict(color="#185fa5", width=2)), row=1, col=3)
        fig.add_hline(y=15.0, line=dict(color="#c0392b", dash="dash", width=1),
                      annotation_text="уставка 15 мм", annotation_font_size=10, row=1, col=3)
        fig.update_xaxes(title_text="время, с", row=1, col=3); fig.update_yaxes(title_text="зазор, мм", row=1, col=3)
    fig.update_layout(template="plotly_white", height=360, width=1000, showlegend=False,
                      title_text=f"Привод зазора (≤0.5 мм/с) → картина едет за реализованным зазором ({SCREEN})")
    out = os.path.join("workspace", "results", "screens", "gap_sweep_demo.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").write(to_html(fig, include_plotlyjs=True, full_html=True))
    return out


if __name__ == "__main__":
    main()
