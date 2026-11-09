# -*- coding: utf-8 -*-
"""Режим 5 — бесконечный real-time контур с ручками, безголовый авто-HTML-дашборд (без GUI).

Непрерывный цикл движка (TbtPlant + Scheduler, орбитный FOFB) гоняет реальную машину; раз в цикл пишет
АВТО-ОБНОВЛЯЕМЫЙ HTML-дашборд (meta-refresh): сенсоры (орбита X/Y rms), свет (photon rms), фазовый портрет,
состояние контуров. РУЧКИ — в конфиг-файле control/mode5_knobs.json (FOFB вкл/выкл, амплитуда возмущения):
цикл перечитывает его каждый шаг по mtime (как live_visualizer опрашивает файлы), так что правка на лету
подхватывается — открытый в браузере дашборд сам обновляется и показывает реакцию.

В headless-демо главный поток САМ флипает ручку FOFB в середине прогона (имитация оператора), чтобы в
дашборде было видно, как орбита падает при включении СУ. В реальном использовании — просто правьте JSON.
ЗАПУСК: py -3.14 control/mode5_realtime.py [smoke]   (открыть workspace/results/mode5_dashboard.html в браузере)
"""
import argparse
import json
import os
import sys
import time

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import numpy as np
import plotly.graph_objects as go

from control import ExperimentState, Limits, Scheduler
from control.tbt_plant import TbtPlant
from control.controllers.integral_svd import IntegralSVD

KNOBS = os.path.join(os.path.dirname(__file__), "mode5_knobs.json")
OUT = os.path.join("workspace", "results", "mode5_dashboard.html")
DEFAULTS = {"fofb_on": False, "disturbance_amp_um": 0.3, "_note": "ручки режима 5: правьте на лету — цикл подхватит"}


def rms(v):
    return float(np.sqrt(np.mean(np.asarray(v) ** 2)))


def read_knobs():
    try:
        with open(KNOBS, encoding="utf-8") as f:
            k = json.load(f)
        return bool(k.get("fofb_on", False)), float(k.get("disturbance_amp_um", 0.3))
    except Exception:
        return DEFAULTS["fofb_on"], DEFAULTS["disturbance_amp_um"]


def write_knobs(fofb_on, amp_um):
    with open(KNOBS, "w", encoding="utf-8") as f:
        json.dump({"fofb_on": fofb_on, "disturbance_amp_um": amp_um, "_note": DEFAULTS["_note"]}, f, ensure_ascii=False, indent=2)


def dashboard(cyc, hist, cen, fofb_on, amp_um, refresh):
    """Собрать авто-обновляемый HTML: KPI + орбита X/Y rms во времени + фазовый портрет."""
    t = np.arange(len(hist["x"]))
    fo = go.Figure()
    for key, name, col in (("x", "орбита X rms", "#2b7bba"), ("y", "орбита Y rms", "#e15759"),
                           ("ph", "свет (photon) rms", "#59a14f")):
        fo.add_trace(go.Scattergl(x=t, y=np.array(hist[key]) * 1e6, mode="lines", name=name, line=dict(width=1.6, color=col)))
    fo.update_layout(height=320, margin=dict(l=56, r=12, t=10, b=40), plot_bgcolor="white",
                     legend=dict(orientation="h", y=1.02, yanchor="bottom"))
    fo.update_xaxes(title_text="цикл"); fo.update_yaxes(title_text="мкм", rangemode="tozero")
    fp = go.Figure(go.Scattergl(x=cen[:, 0] * 1e6, y=cen[:, 1] * 1e6, mode="markers",
                                marker=dict(size=3, color=np.arange(len(cen)), colorscale="Viridis")))
    fp.update_layout(height=320, margin=dict(l=56, r=12, t=10, b=40), plot_bgcolor="white")
    fp.update_xaxes(title_text="центроид x, мкм"); fp.update_yaxes(title_text="центроид x′, мкрад")

    def kpi(l, v, u):
        return f'<div class="kpi"><div class="l">{l}</div><div class="v">{v}</div><div class="u">{u}</div></div>'
    kpis = "".join([
        kpi("цикл", str(cyc), "такт"),
        kpi("СУ (FOFB)", "ВКЛ" if fofb_on else "ВЫКЛ", "ручка"),
        kpi("возмущение", f"{amp_um:.2f}", "мкм"),
        kpi("орбита X", f"{hist['x'][-1]*1e6:.2f}", "мкм rms"),
        kpi("орбита Y", f"{hist['y'][-1]*1e6:.2f}", "мкм rms"),
        kpi("свет", f"{hist['ph'][-1]*1e6:.2f}", "мкм rms"),
    ])
    o1 = fo.to_html(include_plotlyjs="inline", full_html=False)
    o2 = fp.to_html(include_plotlyjs=False, full_html=False)
    return f"""<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta http-equiv="refresh" content="{refresh}"><title>Режим 5 — real-time</title><style>
body{{margin:0;background:#fcfcfb;color:#0b0b0b;font-family:system-ui,sans-serif;padding:18px}}
h1{{font-size:19px;margin:0 0 2px}}.sub{{color:#52514e;font-size:12px;margin:0 0 12px}}
.kpis{{display:grid;grid-template-columns:repeat(auto-fit,minmax(120px,1fr));gap:8px;margin:10px 0}}
.kpi{{background:#fff;border:1px solid rgba(0,0,0,.14);border-radius:9px;padding:8px 11px}}
.kpi .l{{font-size:11px;color:#52514e}}.kpi .v{{font-size:18px;font-weight:600;color:#185fa5}}.kpi .u{{font-size:10px;color:#898781}}
.pic{{background:#fff;border:1px solid rgba(0,0,0,.14);border-radius:10px;padding:10px;margin:8px 0}}
h2{{font-size:14px;margin:14px 0 4px}}</style></head><body>
<h1>Режим 5 — real-time контур орбиты (безголовый дашборд)</h1>
<p class="sub">авто-обновление каждые {refresh} с · ручки: control/mode5_knobs.json (FOFB вкл/выкл, амплитуда) — правьте на лету</p>
<div class="kpis">{kpis}</div>
<h2>Сенсоры во времени: орбита X/Y и свет (photon) rms</h2><div class="pic">{o1}</div>
<h2>Фазовый портрет центроида (последние обороты)</h2><div class="pic">{o2}</div>
</body></html>"""


def main():
    ap = argparse.ArgumentParser(description="режим 5 — real-time дашборд")
    ap.add_argument("mode", nargs="?", default="")
    ap.add_argument("--cycles", type=int, default=None, help="циклов (headless-демо; в реале — до Ctrl+C)")
    ap.add_argument("--refresh", type=int, default=2, help="период авто-обновления HTML, с")
    ap.add_argument("--out", default=OUT)
    a = ap.parse_args()
    smoke = a.mode == "smoke"
    cycles = a.cycles if a.cycles else (16 if smoke else 26)
    steps_per_cycle = 9                                       # = окно усреднения (один шаг FOFB за цикл)

    write_knobs(False, 0.3)                                   # старт: СУ выкл (демо флипнет в середине)
    print("=== РЕЖИМ 5: real-time контур + авто-HTML-дашборд ===")
    print(f"ручки: {KNOBS} · дашборд: {a.out} (meta-refresh {a.refresh} с)")
    print(f"{cycles} циклов демо; в середине ИМИТИРУЕМ оператора — флип FOFB ВЫКЛ→ВКЛ\n", flush=True)

    plant = TbtPlant(lattice=26, nsl=8, dist_amp=0.3e-6, dist_rho=1.0, settle_turns=(120 if smoke else 250),
                     meas_noise=1e-7, n_photon=1)
    dt = plant.stride * plant.TREV
    schema = plant.schema()
    state = ExperimentState(schema, limits=Limits(schema), history=plant.history())
    sch = Scheduler(plant, state, dt)
    rate = 1.0 / (dt * 9)
    sch.add(IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5, rcond=1e-2), rate_hz=rate, enabled=False)
    sch.add(IntegralSVD("fofb_y", sensor="ebpm.y", actuator_key="vkick", gain=0.5, rcond=1e-2), rate_hz=rate, enabled=False)

    hist = {"x": [], "y": [], "ph": []}
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    for cyc in range(cycles):
        if cyc == cycles // 2:                               # ИМИТАЦИЯ ОПЕРАТОРА: включить FOFB (как правка JSON)
            write_knobs(True, 0.3); print(f"  [цикл {cyc}] оператор включил FOFB (ручка mode5_knobs.json)", flush=True)
        fofb_on, amp_um = read_knobs()                       # перечитываем ручки каждый цикл (на лету)
        sch.set_enabled("fofb_x", fofb_on); sch.set_enabled("fofb_y", fofb_on)
        plant.dist_amp = amp_um * 1e-6; plant._amp = amp_um * 1e-6
        for _ in range(steps_per_cycle):
            sch.step()
        hist["x"].append(rms(state.observables["ebpm.x"])); hist["y"].append(rms(state.observables["ebpm.y"]))
        hist["ph"].append(rms(state.observables.get("photon.x", [0.0])))
        cen = np.asarray(plant.centroid_tbt[-400:]) if plant.centroid_tbt else np.zeros((1, 4))
        open(a.out, "w", encoding="utf-8").write(dashboard(cyc + 1, hist, cen, fofb_on, amp_um, a.refresh))
        print(f"  цикл {cyc+1:3d}: FOFB {'ВКЛ' if fofb_on else 'ВЫКЛ'} | орбита X {hist['x'][-1]*1e6:6.2f} "
              f"Y {hist['y'][-1]*1e6:6.2f} | свет {hist['ph'][-1]*1e6:6.2f} мкм", flush=True)
    plant.close()
    xo = np.mean(hist["x"][:cycles//2]) * 1e6; xn = np.mean(hist["x"][cycles//2:]) * 1e6
    print(f"\nорбита X: до включения {xo:.2f} → после {xn:.2f} мкм (×{xo/max(xn,1e-9):.1f})")
    print(f"[OK] {a.out} — открыть в браузере (авто-обновляется); ручки — {KNOBS}")


if __name__ == "__main__":
    main()
