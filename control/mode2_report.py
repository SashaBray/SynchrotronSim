# -*- coding: utf-8 -*-
"""Режим 2 — фазовый портрет + плоскостные сенсоры, синхротронное излучение ВКЛ/ВЫКЛ.

Один пучок (одна частица — чистый портрет без филаментации) инжектируется со смещением и ведётся много
оборотов на 6D-кольце СИЛА-26 (ВЧ фаза 165°) через движок TbtPlant. СИ ВЫКЛ → консервативное движение
(замкнутый эллипс в (x, x′)); СИ ВКЛ → радиационное ЗАТУХАНИЕ (спираль внутрь к равновесию). Это и есть
смысл режима 2 — увидеть динамику и затухание. Сенсоры плоскости старта берём из пооборотного центроида.

Движок — тот же TbtPlant (track_tbt), без обратной связи (свободная динамика); СИ переключается флагом sr.
Рендер — интерактивный HTML (plotly) через traj_analyzer._assemble. ЗАПУСК: py -3.14 control/mode2_report.py [smoke]
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import plotly.graph_objects as go

import traj_analyzer as ta
from control import ExperimentState
from control.tbt_plant import TbtPlant


def run(sr, turns, x0=1e-4, y0=3e-5):
    """Инжекция одной частицы со смещением, свободная динамика `turns` оборотов; возврат центроида (x,px,y,py)."""
    plant = TbtPlant(lattice=26, nsl=8, n_part=1, sr=sr, dist_amp=0.0, settle_turns=turns,
                     resp_dims=4, n_photon=0, meas_noise=0.0)
    off = plant.CO.copy(); off[0] += x0; off[2] += y0        # бетатронное смещение при инжекции
    plant._bunch0 = np.array([off])
    state = ExperimentState(plant.schema(), history=plant.history())
    plant.reset(state)                                       # трекинг settle_turns накапливает центроид (транзиент)
    cen = np.asarray(plant.centroid_tbt)                     # (turns, 4) [x, px, y, py]
    plant.close()
    return cen


def fig_phase(cen, title):
    x = cen[:, 0] * 1e6; px = cen[:, 1] * 1e6; t = np.arange(len(x))
    fig = go.Figure(go.Scattergl(x=x, y=px, mode="markers",
                                 marker=dict(size=3, color=t, colorscale="Viridis",
                                             colorbar=dict(title="оборот", thickness=12, len=0.9))))
    fig.update_layout(title=dict(text=title, x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=58, r=18, t=40, b=48), plot_bgcolor="white")
    fig.update_xaxes(title_text="центроид x, мкм"); fig.update_yaxes(title_text="центроид x′, мкрад")
    return ta._grid(fig, legend=False)


def fig_amp(cen_off, cen_on):
    def env(c):
        x = c[:, 0] * 1e6; px = c[:, 1] * 1e6
        return np.sqrt(x ** 2 + px ** 2)                     # огибающая (норм. бетатронная амплитуда, у.е.)
    fig = go.Figure()
    for c, name, col in ((cen_off, "СИ выкл (сохранение)", "#e15759"), (cen_on, "СИ вкл (затухание)", "#2b7bba")):
        fig.add_trace(go.Scattergl(x=np.arange(len(c)), y=env(c), mode="lines", name=name, line=dict(color=col, width=1.4)))
    fig.update_layout(title=dict(text="Огибающая амплитуды центроида (√(x²+x′²)) — радиационное затухание",
                                 x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=58, r=18, t=40, b=60), plot_bgcolor="white",
                      legend=dict(orientation="h", y=-0.18, x=0.5, xanchor="center"))
    fig.update_xaxes(title_text="оборот"); fig.update_yaxes(title_text="амплитуда, у.е.", rangemode="tozero")
    return ta._grid(fig)


def main():
    ap = argparse.ArgumentParser(description="режим 2 — фазовый портрет, СИ вкл/выкл")
    ap.add_argument("mode", nargs="?", default="")
    ap.add_argument("--turns", type=int, default=None)
    ap.add_argument("--out", default=os.path.join("workspace", "results", "mode2_phase.html"))
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()
    smoke = a.mode == "smoke"
    turns = a.turns if a.turns else (1500 if smoke else 5000)

    print("=== РЕЖИМ 2: фазовый портрет центроида, СИ вкл/выкл ===")
    print(f"инжекция 1 частицы со смещением x=100 мкм; {turns} оборотов на 6D-кольце (ВЧ 165°)\n", flush=True)
    print("СИ ВЫКЛ (консервативно) …", flush=True)
    cen_off = run(sr=False, turns=turns)
    print("СИ ВКЛ (радиационное затухание) …", flush=True)
    cen_on = run(sr=True, turns=turns)

    def amp0_end(c):
        e = np.sqrt((c[:, 0] * 1e6) ** 2 + (c[:, 1] * 1e6) ** 2)
        n = max(len(e) // 20, 1)
        return float(np.mean(e[:n])), float(np.mean(e[-n:]))
    a0_off, ae_off = amp0_end(cen_off); a0_on, ae_on = amp0_end(cen_on)
    print(f"амплитуда центроида начало→конец: СИ выкл {a0_off:.2f}→{ae_off:.2f} у.е. (×{ae_off/max(a0_off,1e-9):.2f}); "
          f"СИ вкл {a0_on:.2f}→{ae_on:.2f} (×{ae_on/max(a0_on,1e-9):.2f} — затухание)")

    figs = [(fig_phase(cen_off, "СИ ВЫКЛ: (x, x′) — замкнутый эллипс (сохранение)"), 480),
            (fig_phase(cen_on, "СИ ВКЛ: (x, x′) — спираль внутрь (радиационное затухание)"), 480),
            (fig_amp(cen_off, cen_on), 420)]
    title = "SynchrotronSim — режим 2: фазовый портрет центроида (СИ вкл/выкл)"
    sub = (f"решётка 26 · инжекция x=100 мкм · {turns} оборотов 6D (ВЧ 165°) · СИ выкл {ae_off/max(a0_off,1e-9):.2f}× "
           f"vs СИ вкл {ae_on/max(a0_on,1e-9):.2f}× (затухание) · графики интерактивны")
    html = ta._assemble(26, 0, title, sub, figs)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w", encoding="utf-8").write(html)
    print(f"\n[OK] {a.out}")
    if not a.no_open:
        try:
            import webbrowser
            webbrowser.open("file://" + os.path.abspath(a.out))
        except Exception:
            pass


if __name__ == "__main__":
    main()
