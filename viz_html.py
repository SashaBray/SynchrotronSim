# -*- coding: utf-8 -*-
"""Экспорт 3D-сцены (магниты кольца + траектории пучка) в САМОДОСТАТОЧНЫЙ интерактивный HTML
(plotly, встроенный JS) — файл можно передать другому человеку, он крутит/зумит в браузере БЕЗ
установок и интернета. Переиспользует геометрию visualize.py (STL-модели, чтение траекторий).

Магниты — plotly Mesh3d из тех же STL (окрашены по типу), пучок — Scatter3d линии по частицам.
Две степени детализации:
  • ДЕТАЛЬНАЯ — настоящий меш STL (при необходимости прорежен, но форма сохраняется);
  • ЛЁГКАЯ (--glyph) — магнит = блок по габаритам STL (12 треуг.), для полного кольца в разы легче.
--both пишет ОБЕ версии (_light + _detailed). Полный STL кольца (101) — миллионы треуг., поэтому
детальная версия авто-прореживается под цель ~300k треуг. (--decimate -1, по умолчанию); --decimate 0
форсит полный STL, --decimate 0.9 — вручную.
Запуск: py -3.14 viz_html.py <exp_id> [--both | --glyph] [--decimate 0..1] [--out путь.html]"""
import sys, os, csv, argparse
os.environ.setdefault("SS_SECTOR", "1")     # график отклонения считаем СЕКТОРНЫМ телом (бит-в-бит; иначе
#                                             боксовое тело + криволин. рамки разбегаются). 3D из bin не зависит.
import numpy as np
import pyvista as pv
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from collections import defaultdict, Counter
import visualize as vz
import matrix_traj_export as mte

WS = "workspace"
HTML_DIR = os.path.join(WS, "results", "html")               # отдельная папка для передаваемых HTML-экспортов
COLORS = vz.VIZ_CFG.get("particle_colors", ["magenta", "cyan", "lime", "yellow", "red"])
LABELS = vz.VIZ_CFG.get("particle_labels", ["P1", "P2", "P3", "P4", "P5"])
KINDCOL = {"quad": "#4e79a7", "dipole": "#e15759", "sext": "#59a14f", "oct": "#b07aa1",
           "dq": "#f28e2b", "kicker": "#9c755f", "monitor": "#bab0ac"}
TRAJ_MAX = 8000                                              # макс. точек на частицу (длинную траекторию прореживаем)
TARGET_TRI = 300_000                                        # цель по треуг. для авто-прореживания детальной версии


def _read(name):
    return list(csv.DictReader(open(os.path.join(WS, "tables", name), encoding="utf-8"), delimiter=";"))


def _total_tris(lc, hw):                                     # оценка полного числа треуг. STL (бинарный STL: (size-84)/50)
    tot = 0
    for dev, n in Counter(r["device_type_id"] for r in lc).items():
        stl = hw.get(dev, {}).get("stl_path", "-")
        if stl in ("-", "", "nan"):
            continue
        p = os.path.join(WS, "geometry", stl)
        if os.path.exists(p):
            tot += max(0, (os.path.getsize(p) - 84) // 50) * n
    return tot


def _deviation_fig(dev, labels):
    """2D-график отклонения пробного пучка от опорной орбиты: x(s) сверху, y(s) снизу, линия на частицу."""
    col = ["#555555", "#e15759", "#f28e2b", "#4e79a7", "#59a14f"]     # опорная серая, ±x тёплые, ±y холодные
    fig2 = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.09,
                         subplot_titles=("Горизонтальное отклонение x(s)", "Вертикальное отклонение y(s)"))
    for i, d in enumerate(dev):
        if len(d) < 2:
            continue
        s = d[:, 0]; x = d[:, 1] * 1e3; y = d[:, 2] * 1e3            # мм
        nm = labels[i] if i < len(labels) else f"P{i+1}"
        c = col[i % len(col)]; w = 2.0 if i == 0 else 1.1
        fig2.add_trace(go.Scattergl(x=s, y=x, mode="lines", line=dict(width=w, color=c),
                                    name=nm, legendgroup=nm), row=1, col=1)
        fig2.add_trace(go.Scattergl(x=s, y=y, mode="lines", line=dict(width=w, color=c),
                                    name=nm, legendgroup=nm, showlegend=False), row=2, col=1)
    fig2.update_xaxes(title_text="s, м", row=2, col=1)
    fig2.update_yaxes(title_text="x, мм", row=1, col=1); fig2.update_yaxes(title_text="y, мм", row=2, col=1)
    fig2.update_xaxes(gridcolor="#eee", zeroline=False)
    fig2.update_yaxes(gridcolor="#eee", zeroline=True, zerolinecolor="#bbb", zerolinewidth=1)
    fig2.update_layout(title="Отклонение пробного пучка от опорной орбиты (±1мм по x и y)",
                       margin=dict(l=60, r=20, t=70, b=45), plot_bgcolor="white",
                       legend=dict(orientation="h", y=1.12, x=0.5, xanchor="center", title="частица"))
    return fig2


def build(exp_id, lat, hw, kind, lc, decimate, glyph, out):
    """Собрать сцену и записать самодостаточный HTML. Возвращает (путь, число_треуг)."""
    stl_cache, grouped = {}, defaultdict(list)               # STL-кэш + группировка инстансов по устройству
    for r in lc:
        dev = r["device_type_id"]; stl = hw.get(dev, {}).get("stl_path", "-")
        if stl in ("-", "", "nan"):
            continue
        if stl not in stl_cache:
            p = os.path.join(WS, "geometry", stl)
            m = pv.read(p).scale(0.001, inplace=False) if os.path.exists(p) else None
            if m is not None and glyph:                      # лёгкая версия: магнит = блок по габаритам STL
                m = pv.Box(bounds=m.bounds).triangulate()
            elif m is not None and decimate > 0:             # детальная: прореженный настоящий меш (форма сохраняется)
                m = m.triangulate().decimate(decimate)
            stl_cache[stl] = m
        if stl_cache[stl] is None:
            continue
        grouped[dev].append((float(r["yaw"]), float(r["pitch"]), float(r["roll"]),
                             float(r["x"]), float(r["y"]), float(r["z"])))

    fig = go.Figure()
    ntri = 0
    for dev, tlist in grouped.items():                       # магниты → Mesh3d, цвет по типу
        big = vz.build_combined_mesh(stl_cache[hw[dev]["stl_path"]], np.asarray(tlist, float))
        V = np.asarray(big.points); F = np.asarray(big.faces).reshape(-1, 4)[:, 1:]
        ntri += len(F)
        fig.add_trace(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                                color=KINDCOL.get(kind.get(dev, ""), "silver"), opacity=0.55,
                                name=dev, hovertext=dev, hoverinfo="text", showlegend=False))

    # опорный путь — ломаная через центры магнитов в порядке кольца (проектная орбита, приближённо)
    lc_ord = sorted(lc, key=lambda r: int(r["instance_id"]))
    ref = np.array([(float(r["x"]), float(r["y"]), float(r["z"])) for r in lc_ord])
    if len(ref) >= 2:
        fig.add_trace(go.Scatter3d(x=ref[:, 0], y=ref[:, 1], z=ref[:, 2], mode="lines",
                                   line=dict(width=2, color="dimgray"), name="опорный путь", opacity=0.6))
    # BPM — маркеры на датчиках положения (kind=monitor)
    bpm = np.array([(float(r["x"]), float(r["y"]), float(r["z"])) for r in lc if kind.get(r["device_type_id"]) == "monitor"])
    if len(bpm):
        fig.add_trace(go.Scatter3d(x=bpm[:, 0], y=bpm[:, 1], z=bpm[:, 2], mode="markers",
                                   marker=dict(size=3, color="black", symbol="diamond"), name=f"BPM ({len(bpm)})"))

    parts = vz.load_trajectories(os.path.join(WS, "results", f"trajectories_exp{exp_id}.bin"))
    for i, p in enumerate(parts):                            # пучок → линии по частицам
        if len(p) < 2:
            continue
        if len(p) > TRAJ_MAX:                                # прореживание длинной траектории (форма кольца сохраняется)
            step = len(p) // TRAJ_MAX + 1
            p = np.vstack([p[::step], p[-1]])
        fig.add_trace(go.Scatter3d(x=p[:, 0], y=p[:, 1], z=p[:, 2], mode="lines",
                                   line=dict(width=5, color=COLORS[i % len(COLORS)]),
                                   name=LABELS[i] if i < len(LABELS) else f"P{i+1}"))

    mode = "блоки" if glyph else ("STL, прорежено " + f"{decimate:g}" if decimate > 0 else "полный STL")
    fig.update_layout(title=f"SynchrotronSim — эксперимент {exp_id} (решётка {lat}, {len(grouped)} типов магнитов; {mode})",
                      scene=dict(aspectmode="data", xaxis_title="X, м", yaxis_title="Y, м", zaxis_title="Z, м",
                                 bgcolor="white", domain=dict(x=[0, 1], y=[0, 1])),  # сцена на всю область
                      autosize=True,                          # тянется на всю ширину окна браузера (без пустого поля справа)
                      legend=dict(title="частицы", x=0.99, xanchor="right", y=0.99, yanchor="top",  # оверлей, не резервирует колонку
                                  bgcolor="rgba(255,255,255,0.6)", bordercolor="lightgray", borderwidth=1),
                      margin=dict(l=0, r=0, t=40, b=0))
    # 2D-график отклонения ПРОБНОГО пучка от опорной орбиты (перетрек x,y секторным телом) — под 3D-сценой.
    dev_html = ""
    try:
        dev, dlabels = mte.deviation_paths(int(lat))
        dev_html = _deviation_fig(dev, dlabels).to_html(
            include_plotlyjs=False, full_html=False, default_width="100%", default_height="640",
            config={"responsive": True})
    except Exception as e:                                        # не рушим 3D, если отклонение не собралось
        print(f"[!] график отклонения пропущен: {type(e).__name__}: {e}")

    s3d = fig.to_html(include_plotlyjs="inline", full_html=False, default_width="100%",
                      default_height="78vh", config={"responsive": True})
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:                   # 3D сверху, график отклонения снизу — один файл
        f.write("<!doctype html><html><head><meta charset='utf-8'>"
                f"<title>SynchrotronSim — эксперимент {exp_id}</title></head>"
                f"<body style='margin:0;font-family:sans-serif'>{s3d}{dev_html}</body></html>")
    print(f"[OK] {out}")
    print(f"     {os.path.getsize(out)/1e6:.1f} МБ ({mode}); {ntri:,} треуг. магнитов, {len(parts)} частиц"
          f"{' + график отклонения' if dev_html else ''}. Самодостаточный — крутить/зумить, передавать файлом.")
    return out, ntri


def _resolve_decimate(decimate, glyph, lc, hw):
    """decimate<0 (авто): для детальной версии прорежаем полное кольцо под TARGET_TRI, малое — полный STL."""
    if glyph or decimate >= 0:
        return max(0.0, decimate)
    tot = _total_tris(lc, hw)
    if tot > TARGET_TRI:
        dd = round(1 - TARGET_TRI / tot, 3)
        print(f"[авто] детальная версия прорежена --decimate {dd} (полный STL {tot:,} треуг. → цель {TARGET_TRI:,}); "
              f"--decimate 0 форсит полный STL")
        return dd
    return 0.0


def main(exp_id, decimate=-1.0, out=None, glyph=False, both=False):
    exp = {r["exp_id"]: r for r in _read("experiments.csv")}[str(exp_id)]
    lat = exp["lattice_id"]
    hw = {r["device_type_id"]: r for r in _read("hardware_library.csv")}
    kind = {r["device_type_id"]: r["kind"] for r in _read("element_transforms.csv")}
    lc = [r for r in _read("lattice_configs.csv") if r["lattice_id"] == lat]

    base = out or os.path.join(HTML_DIR, f"exp{exp_id}_3d.html")
    stem = base[:-5] if base.lower().endswith(".html") else base
    if both:                                                 # обе версии: лёгкая (блоки) + детальная (STL)
        build(exp_id, lat, hw, kind, lc, 0.0, True, stem + "_light.html")
        build(exp_id, lat, hw, kind, lc, _resolve_decimate(decimate, False, lc, hw), False, stem + "_detailed.html")
    else:
        build(exp_id, lat, hw, kind, lc, _resolve_decimate(decimate, glyph, lc, hw), glyph, base)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("exp", type=int)
    ap.add_argument("--decimate", type=float, default=-1.0,
                    help="доля удаляемых треуг. STL (0..1); -1=авто (детальную прорежает под ~300k треуг.), 0=полный STL")
    ap.add_argument("--glyph", action="store_true", help="лёгкая версия: магниты блоками по габаритам STL")
    ap.add_argument("--both", action="store_true", help="сгенерировать ОБЕ версии: _light (блоки) и _detailed (STL)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    main(a.exp, a.decimate, a.out, a.glyph, a.both)
