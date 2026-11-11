# -*- coding: utf-8 -*-
"""Демо ЖИВОГО ЗАЗОРА: во время прогона TbtPlant двигаем зазор ондулятора → интерференционная картина едет.

Показывает связь, которую мы завели: `state.undulators[<имя>].gap_mm` правится на лету (здесь — рампа
5→15 мм по окнам), `ScreenRecorder` читает зазор своего ондулятора-источника каждое окно и пересчитывает
картину (K→E₁→конус/поток). На выходе — самодостаточный HTML: картина при НАЧАЛЬНОМ и КОНЕЧНОМ зазоре +
кривая E₁(зазор). Контур увидел бы тот же зазор через StateView.undulators.

ЗАПУСК: py -3.14 control/gap_sweep_demo.py [smoke]
"""
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")

import numpy as np

from control import ExperimentState, Scheduler
from control.tbt_plant import TbtPlant
from control.recorders import ScreenRecorder

SCREEN = "SCR_U19"            # ондуляторный экран (источник und_19), 15 м
SOURCE = "und_19"


def main():
    smoke = "smoke" in sys.argv
    n_win = 5 if smoke else 9
    settle = 60 if smoke else 120
    gaps = np.linspace(5.0, 15.0, n_win)         # рампа зазора по окнам, мм

    plant = TbtPlant(lattice=26, n_photon=8, settle_turns=settle, gap_mm=float(gaps[0]))
    state = ExperimentState(plant.schema(), history=plant.history())
    sch = Scheduler(plant, state, dt=plant.stride * plant.TREV)   # reset() сеет state.undulators
    rec = ScreenRecorder(params=SCREEN, cadence="final")

    print(f"=== ДЕМО ЖИВОГО ЗАЗОРА: {SCREEN} (источник {SOURCE}), рампа {gaps[0]:.0f}→{gaps[-1]:.0f} мм ===")
    for w, gmm in enumerate(gaps):
        state.undulators[SOURCE].gap_mm = float(gmm)   # ← двигаем зазор на лету (как оператор/контур)
        sch.step()
        rec.sample(state, w, plant)
        k = state.undulators[SOURCE].K()
        print(f"  окно {w}: зазор={gmm:5.2f} мм  K={k:.3f}", flush=True)
    rec.finalize(state, plant)

    art = rec.artifact()[1]
    rows = {s: v for s, _, v in rec.result_rows() if _ == 0}
    e1s, e1e = float(art["energy_start_keV"][0]), float(art["energy_end_keV"][0])
    print(f"\nE₁: {e1s:.2f} кэВ (зазор {art['gap_start_mm'][0]:.1f}) → {e1e:.2f} кэВ (зазор {art['gap_end_mm'][0]:.1f})")
    print(f"резкость (джиттер пучка): {rows.get('screen.sharpness_pct', 0):.1f}%")

    out = _render(art)
    print(f"[OK] {out}")


def _render(art):
    import plotly.graph_objects as go
    from plotly.io import to_html
    from plotly.subplots import make_subplots

    ext = art["extent"]; w_mm = float(art["window_mm"][0])
    ser = art["gap_series"]                              # (окно, зазор, E₁)
    fig = make_subplots(rows=1, cols=3, subplot_titles=(
        f"картина @ зазор {art['gap_start_mm'][0]:.1f} мм (E₁={art['energy_start_keV'][0]:.1f} кэВ)",
        f"картина @ зазор {art['gap_end_mm'][0]:.1f} мм (E₁={art['energy_end_keV'][0]:.1f} кэВ)",
        "E₁ vs зазор (по окнам)"))
    for col, key in ((1, "ideal_start"), (2, "ideal_end")):
        fig.add_trace(go.Heatmap(z=art[key], colorscale="Inferno", showscale=False), row=1, col=col)
        fig.update_xaxes(title_text="мм", row=1, col=col); fig.update_yaxes(title_text="мм", row=1, col=col)
    if ser.size:
        fig.add_trace(go.Scatter(x=ser[:, 1], y=ser[:, 2], mode="lines+markers",
                                 line=dict(color="#185fa5", width=2)), row=1, col=3)
        fig.update_xaxes(title_text="зазор, мм", row=1, col=3); fig.update_yaxes(title_text="E₁, кэВ", row=1, col=3)
    fig.update_layout(template="plotly_white", height=360, width=1000, showlegend=False,
                      title_text=f"Живой зазор → интерференционная картина ({SCREEN}); окно экрана ≈ {w_mm:.1f} мм")
    out = os.path.join("workspace", "results", "screens", "gap_sweep_demo.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").write(to_html(fig, include_plotlyjs=True, full_html=True))
    return out


if __name__ == "__main__":
    main()
