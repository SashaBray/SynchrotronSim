# -*- coding: utf-8 -*-
"""АНАЛИЗАТОР ТРАЕКТОРИЙ SynchrotronSim в ОТНОСИТЕЛЬНЫХ (локальных) координатах матричного ядра.

Читает results/trajectories_local_exp{N}.npz (пишется `track3.py --mode kick-phys --exp N --local`):
на точку — S (опорная арк-длина, м), x, x'≡px, y, y'≡py, ct, δ, индекс элемента. Работает В ЛОКАЛЬНЫХ
координатах: x,y = отклонение от ОПОРНОЙ ОРБИТЫ (секторная модель m.body, бит-в-бит elegant).

Считает по каждой частице СРЕДНЕЕ/МАКС/МИН отклонение от опорной орбиты и делает ЧИСЛЕННОЕ ПРЕОБРАЗОВАНИЕ
ФУРЬЕ x(S), y(S) — доминирующие частоты колебаний (число колебаний за оборот ≈ бетатронный тюн).
Результат — ИНТЕРАКТИВНОЕ ОКНО (самодостаточный HTML, авто-открытие в браузере): сверху поле с общей
информацией, ниже визуализация координат x(S)/y(S), ниже спектр колебаний (Фурье).

Запуск:  py -3.14 traj_analyzer.py <exp_id> [--no-open] [--out файл.html] [--csv]
Сначала выгрузить локальные координаты:  py -3.14 track3.py --mode kick-phys --exp <N> --local
"""
import os
import sys
import argparse
import webbrowser
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import matrix_traj_export as mte

RESULTS = os.path.join("workspace", "results")
COLORS = ["#555555", "#e15759", "#f28e2b", "#4e79a7", "#59a14f", "#b07aa1", "#9c755f", "#edc948"]


def fft_spectrum(S, v):
    """Амплитудный спектр v(S): ресемплинг на равномерную сетку, вычет среднего, |rFFT|·2/N.
    Возврат (Q, amp): Q — число колебаний на пройденном пути (гармоника ≈ тюн за оборот), amp — амплитуда (те же ед., что v)."""
    n = len(S)
    if n < 8 or S[-1] <= S[0]:
        return np.array([0.0]), np.array([0.0])
    Su = np.linspace(float(S[0]), float(S[-1]), n)
    vu = np.interp(Su, S, v)
    vu = vu - np.mean(vu)
    amp = np.abs(np.fft.rfft(vu)) * 2.0 / n
    return np.arange(len(amp), dtype=float), amp        # k = число полных колебаний за пройденный путь


def _dominant(k, amp):
    """Доминирующая частота (гармоника с макс. амплитудой, кроме DC=0)."""
    if len(amp) < 2:
        return 0, 0.0
    j = int(np.argmax(amp[1:])) + 1
    return int(k[j]), float(amp[j])


def metrics(p, closed):
    """Метрики частицы: отклонение от опорной орбиты r=√(x²+y²), спектры x/y, доминирующие частоты."""
    x, y, S = p["x"], p["y"], p["S"]
    r_mm = np.hypot(x, y) * 1e3
    kx, ax = fft_spectrum(S, x); ky, ay = fft_spectrum(S, y)
    dqx, dax = _dominant(kx, ax); dqy, day = _dominant(ky, ay)
    clo = float(np.hypot(x[-1] - x[0], y[-1] - y[0]) * 1e3) if (closed and p["survived"] and len(x) >= 2) else None
    return dict(survived=p["survived"], n=len(S), S_end=float(S[-1]) if len(S) else 0.0,
                mean_mm=float(np.mean(r_mm)), max_mm=float(np.max(r_mm)), min_mm=float(np.min(r_mm)),
                dqx=dqx, dqy=dqy, clo=clo, kx=kx, ax=ax * 1e3, ky=ky, ay=ay * 1e3)   # амплитуды в мм


def build_figure(d, res):
    """3-панельная интерактивная фигура: таблица общей информации / координаты x(S),y(S) / спектр Фурье."""
    parts, closed, ring = d["parts"], d["closed"], d["ring_len"]
    P = len(parts)
    fig = make_subplots(
        rows=3, cols=1, row_heights=[0.22, 0.40, 0.38], vertical_spacing=0.09,
        specs=[[{"type": "table"}], [{"type": "xy"}], [{"type": "xy"}]],
        subplot_titles=("", "Координаты: отклонение от опорной орбиты  x(S), y(S)",
                        "Спектр колебаний (Фурье): амплитуда по числу колебаний за оборот"))

    # --- панель 1: таблица общей информации по частицам ---
    hdr = ["№", "дошла", "ср.|откл|, мм", "макс, мм", "мин, мм", "дом. Qx", "дом. Qy", "замык. |Δ|, мм"]
    col = [[str(i) for i in range(P)],
           ["да" if r["survived"] else "НЕТ" for r in res],
           [f"{r['mean_mm']:.3f}" for r in res], [f"{r['max_mm']:.3f}" for r in res],
           [f"{r['min_mm']:.3f}" for r in res], [str(r["dqx"]) for r in res], [str(r["dqy"]) for r in res],
           [("—" if r["clo"] is None else f"{r['clo']:.3f}") for r in res]]
    fig.add_trace(go.Table(
        columnwidth=[26, 40, 70, 55, 55, 50, 50, 75],
        header=dict(values=hdr, fill_color="#4e79a7", font=dict(color="white", size=12), align="center"),
        cells=dict(values=col, align="center", height=22,
                   fill_color=[["#f4f6f8" if i % 2 else "white" for i in range(P)]])), row=1, col=1)

    # --- панель 2: координаты x(S) (сплошная), y(S) (пунктир) по частицам ---
    for i, p in enumerate(parts):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=p["S"], y=p["x"] * 1e3, mode="lines", legendgroup=f"ч{i}",
                                   name=f"ч{i} · x", line=dict(color=c, width=1.6)), row=2, col=1)
        fig.add_trace(go.Scattergl(x=p["S"], y=p["y"] * 1e3, mode="lines", legendgroup=f"ч{i}",
                                   name=f"ч{i} · y", line=dict(color=c, width=1.2, dash="dot"),
                                   showlegend=False), row=2, col=1)
    # --- панель 3: спектр Фурье x (сплошная), y (пунктир) ---
    qmax = 0
    for i, r in enumerate(res):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=r["kx"], y=r["ax"], mode="lines", legendgroup=f"ч{i}",
                                   name=f"ч{i} · x", line=dict(color=c, width=1.4), showlegend=False), row=3, col=1)
        fig.add_trace(go.Scattergl(x=r["ky"], y=r["ay"], mode="lines", legendgroup=f"ч{i}",
                                   name=f"ч{i} · y", line=dict(color=c, width=1.1, dash="dot"), showlegend=False), row=3, col=1)
        qmax = max(qmax, r["dqx"], r["dqy"])
    fig.update_xaxes(title_text="S, м", row=2, col=1)
    fig.update_yaxes(title_text="отклонение, мм", row=2, col=1, zeroline=True, zerolinecolor="#bbb")
    fig.update_xaxes(title_text="число колебаний за оборот (Q)", row=3, col=1,
                     range=[0, max(40, qmax * 3)])            # видимая зона у пиков; интерактивно зумится дальше
    fig.update_yaxes(title_text="амплитуда, мм", row=3, col=1)

    surv = sum(r["survived"] for r in res)
    mean_all = np.mean([r["mean_mm"] for r in res]); max_all = max(r["max_mm"] for r in res)
    kind = "замкнутое кольцо" if closed else "открытый сегмент"
    fig.update_layout(
        title=dict(text=f"SynchrotronSim — анализ траектории, эксперимент {d['lattice']}·{d['group']} "
                        f"({kind})<br><sub>оборот {ring:.1f} м · элементов {len(d['elem_names'])-1} · частиц {P} · "
                        f"выжило {surv}/{P} · ср.откл. {mean_all:.3f} мм · макс {max_all:.3f} мм · "
                        f"локальные координаты (отклонение от опорной орбиты)</sub>", x=0.5, xanchor="center"),
        margin=dict(l=70, r=30, t=90, b=50), plot_bgcolor="white", hovermode="closest",
        legend=dict(orientation="h", y=-0.04, x=0.5, xanchor="center"))
    return fig


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp", type=int, help="id эксперимента (читает results/trajectories_local_exp{N}.npz)")
    ap.add_argument("--out", default=None, help="путь HTML (по умолч. results/traj_analysis_exp{N}.html)")
    ap.add_argument("--no-open", action="store_true", help="не открывать окно в браузере автоматически")
    ap.add_argument("--csv", action="store_true", help="доп. сохранить сводку в results/traj_analysis_exp{N}.csv")
    a = ap.parse_args()

    path = os.path.join(RESULTS, f"trajectories_local_exp{a.exp}.npz")
    if not os.path.exists(path):
        sys.exit(f"[!] нет {path} — сначала: py -3.14 track3.py --mode kick-phys --exp {a.exp} --local")
    d = mte.read_local_npz(path)
    res = [metrics(p, d["closed"]) for p in d["parts"]]

    # --- краткая сводка в консоль ---
    surv = sum(r["survived"] for r in res)
    print(f"=== Анализ траектории — эксперимент {a.exp} (решётка {d['lattice']}, пучок {d['group']}, "
          f"{'замкнуто' if d['closed'] else 'открыто'}) ===")
    print(f"оборот {d['ring_len']:.2f}м, элементов {len(d['elem_names'])-1}, частиц {len(res)}, выжило {surv}")
    print("\n№  дошла  ср.|откл|,мм  макс,мм  мин,мм   дом.Qx  дом.Qy  замык.,мм")
    for i, r in enumerate(res):
        clo = "—" if r["clo"] is None else f"{r['clo']:.3f}"
        print(f"{i:<2} {'да ' if r['survived'] else 'НЕТ'}   {r['mean_mm']:>9.3f}  {r['max_mm']:>7.3f}  "
              f"{r['min_mm']:>6.3f}   {r['dqx']:>5}  {r['dqy']:>6}   {clo:>8}")

    # --- интерактивное окно (HTML) ---
    out = a.out or os.path.join(RESULTS, f"traj_analysis_exp{a.exp}.html")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    fig = build_figure(d, res)
    fig.write_html(out, include_plotlyjs="inline", full_html=True, config={"responsive": True})
    print(f"\n[окно] {out}")
    if not a.no_open:
        webbrowser.open("file:///" + os.path.abspath(out).replace(os.sep, "/"))

    if a.csv:                                                # человеко-читаемая сводка (запятая-десятичная, ; — рус. Excel)
        cpath = os.path.join(RESULTS, f"traj_analysis_exp{a.exp}.csv")
        with open(cpath, "w", encoding="utf-8") as f:
            f.write("частица;дошла;ср_откл_мм;макс_мм;мин_мм;дом_Qx;дом_Qy;замыкание_мм\n")
            for i, r in enumerate(res):
                row = [i, "да" if r["survived"] else "нет", r["mean_mm"], r["max_mm"], r["min_mm"],
                       r["dqx"], r["dqy"], "" if r["clo"] is None else r["clo"]]
                f.write(";".join(str(c).replace(".", ",") if isinstance(c, float) else str(c) for c in row) + "\n")
        print(f"[CSV] {cpath}")


if __name__ == "__main__":
    main()
