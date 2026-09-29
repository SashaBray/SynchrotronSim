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
import matrix_traj_export as mte

RESULTS = os.path.join("workspace", "results")
COLORS = ["#2b7bba", "#e15759", "#f28e2b", "#59a14f", "#b07aa1", "#9c755f", "#edc948", "#76b7b2"]
GRID = "#e9e9e9"                                             # цвет фоновой сетки графиков


def _grid(fig, legend=True):
    """Единый стиль осей: светлая сетка на фон + (по умолч.) принудительная легенда — plotly скрывает её
    при единственной кривой, а нам нужна всегда. legend=False для панелей с цветовой шкалой (фазовый портрет)."""
    fig.update_xaxes(showgrid=True, gridcolor=GRID, gridwidth=1)
    fig.update_yaxes(showgrid=True, gridcolor=GRID, gridwidth=1)
    if legend:
        fig.update_layout(showlegend=True)
    return fig


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
        return 0.0, 0.0
    j = int(np.argmax(amp[1:])) + 1
    return float(k[j]), float(amp[j])


def metrics(p, closed, n_turns=1):
    """Метрики частицы: отклонение от опорной орбиты r=√(x²+y²), спектр x/y. За n_turns оборотов гармоника FFT
    делится на n_turns ⇒ ось = колебаний за ОБОРОТ (тюн), точнее при большем числе оборотов."""
    x, y, S = p["x"], p["y"], p["S"]
    r_mm = np.hypot(x, y) * 1e3
    kx, ax = fft_spectrum(S, x); ky, ay = fft_spectrum(S, y)
    kx = kx / n_turns; ky = ky / n_turns                # гармоника за весь прогон → колебаний за ОБОРОТ (тюн)
    dqx, _ = _dominant(kx, ax); dqy, _ = _dominant(ky, ay)
    clo = (float(np.hypot(x[-1] - x[0], y[-1] - y[0]) * 1e3)
           if (closed and n_turns == 1 and p["survived"] and len(x) >= 2) else None)
    return dict(survived=p["survived"], surv_turns=p.get("surv_turns", 1), n=len(S),
                mean_mm=float(np.mean(r_mm)), max_mm=float(np.max(r_mm)), min_mm=float(np.min(r_mm)),
                x0=float(x[0]) * 1e3 if len(x) else 0.0, y0=float(y[0]) * 1e3 if len(y) else 0.0,
                dqx=dqx, dqy=dqy, clo=clo, kx=kx, ax=ax * 1e3, ky=ky, ay=ay * 1e3)


def _fig_info(res, n_turns=1):
    """Панель 1: таблица общей информации по частицам. При n_turns>1 — выживание за N оборотов + дробный тюн."""
    P = len(res)
    # начальные x и y — СО ЗНАКОМ: у пучка «центр + спутники ±» модуль делает пары неразличимыми
    x0 = [f"{r.get('x0', 0.0):+.3f}" for r in res]
    y0 = [f"{r.get('y0', 0.0):+.3f}" for r in res]
    if n_turns > 1:
        hdr = ["№", "оборотов", "нач. x, мм", "нач. y, мм", "ср.|откл|, мм", "макс, мм", "тюн Qx", "тюн Qy"]
        col = [[str(i) for i in range(P)],
               [f"{r['surv_turns']}/{n_turns}" for r in res], x0, y0,
               [f"{r['mean_mm']:.3f}" for r in res], [f"{r['max_mm']:.3f}" for r in res],
               [f"{r['dqx']:.2f}" for r in res], [f"{r['dqy']:.2f}" for r in res]]
        cw = [26, 60, 62, 62, 70, 55, 60, 60]
    else:
        hdr = ["№", "дошла", "нач. x, мм", "нач. y, мм", "ср.|откл|, мм", "макс, мм", "мин, мм",
               "дом. Qx", "дом. Qy", "замык. |Δ|, мм"]
        col = [[str(i) for i in range(P)],
               ["да" if r["survived"] else "НЕТ" for r in res], x0, y0,
               [f"{r['mean_mm']:.3f}" for r in res], [f"{r['max_mm']:.3f}" for r in res],
               [f"{r['min_mm']:.3f}" for r in res], [f"{r['dqx']:.0f}" for r in res], [f"{r['dqy']:.0f}" for r in res],
               [("—" if r["clo"] is None else f"{r['clo']:.3f}") for r in res]]
        cw = [26, 40, 62, 62, 70, 55, 55, 50, 50, 75]
    fig = go.Figure(go.Table(
        columnwidth=cw,
        header=dict(values=hdr, fill_color="#4e79a7", font=dict(color="white", size=12), align="center", height=26),
        cells=dict(values=col, align="center", height=22,
                   fill_color=[["#f4f6f8" if i % 2 else "white" for i in range(P)]])))
    fig.update_layout(margin=dict(l=6, r=6, t=6, b=6))
    return fig


def _fig_dev(parts, comp, title):
    """График отклонения от опорной орбиты по одной компоненте (comp='x' или 'y'): кривая на частицу."""
    fig = go.Figure()
    for i, p in enumerate(parts):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=p["S"].astype(np.float32), y=(p[comp] * 1e3).astype(np.float32),  # float32 — вдвое легче
                                   mode="lines", name=f"ч{i}", line=dict(color=c, width=1.5),
                                   hovertemplate="S=%{x:.1f} м<br>" + comp + "=%{y:.3f} мм<extra>ч" + str(i) + "</extra>"))
    fig.update_layout(title=dict(text=title, x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=40, b=72), plot_bgcolor="white", hovermode="closest",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="S, м")
    fig.update_yaxes(title_text=f"отклонение {comp}, мм", zeroline=True, zerolinecolor="#bbb")
    return _grid(fig)


LIGHT_C = 299792458.0                                            # м/с — пучок ультрарелятивистский (v≈c)


def _hz(f):
    if f >= 1e6:
        return f"{f/1e6:.2f} МГц"
    if f >= 1e3:
        return f"{f/1e3:.1f} кГц"
    return f"{f:.0f} Гц"


def _spectrum_peaks(res, n=6, rel=0.10, sep=2):
    """Самые явные пики спектра — ОТДЕЛЬНО по x и по y (огибающая max по частицам на каждой гармонике).
    x и y ищутся независимо: близкие Qx и Qy (напр. 8 и 7) — РАЗНЫЕ тюны, оба помечаются. Раздвижка sep
    убирает соседние бины одного пика ВНУТРИ компоненты. Возврат [(Q, amp, comp)] по убыванию амплитуды."""
    if not res or min(len(r["ax"]) for r in res) < 3:
        return []
    L = min(min(len(r["ax"]) for r in res), min(len(r["ay"]) for r in res),
            len(res[0]["kx"]), len(res[0]["ky"]))
    qx = np.asarray(res[0]["kx"])[:L]; qy = np.asarray(res[0]["ky"])[:L]   # ось Q ЗА ОБОРОТ (гармоника/n_turns)

    def env(key):
        e = np.zeros(L)
        for r in res:
            e = np.maximum(e, np.asarray(r[key])[:L])
        return e

    ex, ey = env("ax"), env("ay")
    gmax = max(float(ex.max()), float(ey.max()), 1e-30)

    def local(e, qgrid, comp):                                   # позиция пика — Q ЗА ОБОРОТ, не сырой индекс
        return [(float(qgrid[j]), float(e[j]), comp) for j in range(1, len(e) - 1)
                if e[j] > e[j - 1] and e[j] >= e[j + 1] and e[j] >= rel * gmax]

    cand = sorted(local(ex, qx, "x") + local(ey, qy, "y"), key=lambda t: -t[1])
    sel = []
    for q, a, comp in cand:
        if all(not (comp == s[2] and abs(q - s[0]) < sep) for s in sel):   # раздвижка только внутри x или y
            sel.append((q, a, comp))
        if len(sel) >= n:
            break
    return sel


def _fig_spectrum(res, ring_len):
    """Спектр Фурье: амплитуда по ЧАСТОТЕ (МГц), x — сплошная, y — пунктир; легенда по частицам. Частота =
    Q·f_оборота, f_оборота=c/кольцо (пучок v≈c). У явных пиков — вертикальный пунктир, ранг, тюн Q и частота."""
    frev = LIGHT_C / ring_len if ring_len > 0 else 0.0           # частота оборота, Гц

    def fMHz(q):                                                 # Q (колебаний за оборот) → частота, МГц
        return q * frev / 1e6

    fig = go.Figure()
    for i, r in enumerate(res):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=(np.asarray(r["kx"]) * frev / 1e6).astype(np.float32), y=r["ax"].astype(np.float32),
                                   mode="lines", name=f"ч{i}", legendgroup=f"ч{i}", line=dict(color=c, width=2.6)))
        fig.add_trace(go.Scattergl(x=(np.asarray(r["ky"]) * frev / 1e6).astype(np.float32), y=r["ay"].astype(np.float32),
                                   mode="lines", name=f"ч{i} · y", legendgroup=f"ч{i}",
                                   line=dict(color=c, width=2.0, dash="dot"), showlegend=False))
    tunemax = max((max(r["dqx"], r["dqy"]) for r in res), default=40.0)   # видимая зона — по ТЮНАМ, не по гармоникам
    xmax = fMHz(max(40.0, tunemax * 1.6))                        # верхняя граница по частоте, МГц
    peaks = [p for p in _spectrum_peaks(res) if fMHz(p[0]) <= xmax]      # подписываем только пики в видимой зоне
    lab = sorted([(fMHz(q), a, comp, q, rk) for rk, (q, a, comp) in enumerate(peaks, 1)], key=lambda t: t[0])
    prev_f = -1e9; row = 0                                       # близкие по частоте подписи разносим по высоте
    for f, a, comp, q, rk in lab:
        row = (row + 1) if (f - prev_f) < 0.05 * xmax else 0
        extra = (row % 3) * 28
        fig.add_vline(x=f, line=dict(color="#6b6b6b", width=1, dash="dash"))
        fig.add_annotation(x=f, y=a, yanchor="bottom", yshift=6 + extra, showarrow=False, align="center",
                           text=f"<b>#{rk}</b> · Q{comp}={round(q, 2):g}<br>{_hz(q * frev)}",
                           font=dict(size=10, color="#333"), bgcolor="rgba(255,255,255,0.8)",
                           bordercolor="#ccc", borderwidth=1, borderpad=1)
        prev_f = f
    fig.update_layout(title=dict(text="Спектр колебаний (Фурье)  (сплошная — x, пунктир — y)",
                                 x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=48, b=72), plot_bgcolor="white",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="частота, МГц", range=[0, xmax])
    ytop = peaks[0][1] * 1.3 if peaks else None                  # запас сверху под подписи (пики заполняют панель полнее)
    fig.update_yaxes(title_text="амплитуда, мм", range=([0, ytop] if ytop else None))
    return _grid(fig)


# --- ДЛИННЫЙ режим (много оборотов): огибающая амплитуды / фазовый портрет / затухание ------------
def _running_mean(a, W):
    """Скользящее среднее окном W (медленная РАВНОВЕСНАЯ орбита). Центрированное, с краевой нормировкой
    (convolve 'same' занижает концы — делим на число реально усреднённых точек)."""
    if W < 2 or len(a) <= W:
        return np.full_like(a, float(np.mean(a)) if len(a) else 0.0)
    k = np.ones(W) / W
    return np.convolve(a, k, mode="same") / np.convolve(np.ones_like(a), k, mode="same")


def metrics_long(p):
    """Метрики частицы в ДЛИННОМ режиме (turn-by-turn у плоскости): положение по оборотам + огибающая
    БЕТАТРОННОЙ амплитуды Ax=√(xb²+(xb'·σ/σ')²) (радиус в фазовом пр-ве), где xb = x − равновесная орбита
    (вычет скользящего среднего). Спад Ax = затухание. Вычет орбиты делает метрику устойчивой к сдвигу
    равновесной орбиты от потерь энергии на излучение (иначе офсет маскирует затухание колебания)."""
    x = np.asarray(p["tt_x"]) * 1e3; xp = np.asarray(p["tt_xp"]) * 1e3      # мм, мрад
    y = np.asarray(p["tt_y"]) * 1e3; yp = np.asarray(p["tt_yp"]) * 1e3
    turns = np.arange(len(x), dtype=float)
    W = int(np.clip(len(x) // 40, 8, 200))                                  # окно ~ несколько бетатронных периодов
    xb = x - _running_mean(x, W); xpb = xp - _running_mean(xp, W)           # бетатрон = вычет равновесной орбиты
    yb = y - _running_mean(y, W); ypb = yp - _running_mean(yp, W)
    sx = np.std(xb) or 1.0; sxp = np.std(xpb) or 1.0
    sy = np.std(yb) or 1.0; syp = np.std(ypb) or 1.0
    Ax = np.sqrt(xb ** 2 + (xpb * sx / sxp) ** 2); Ay = np.sqrt(yb ** 2 + (ypb * sy / syp) ** 2)
    awin = lambda A: (float(np.mean(A[:W])), float(np.mean(A[-W:]))) if len(A) >= W else (
        float(A[0]) if len(A) else 0.0, float(A[-1]) if len(A) else 0.0)   # оконная амплитуда старт/конец (не 1 точка)
    Ax0, Axend = awin(Ax)
    return dict(survived=p["survived"], surv_turns=p["surv_turns"], turns=turns, x=x, xp=xp, y=y, yp=yp,
                Ax=Ax, Ay=Ay, x0=float(x[0]) if len(x) else 0.0, y0=float(y[0]) if len(y) else 0.0,
                xmax=float(np.max(np.abs(x))) if len(x) else 0.0, Ax0=Ax0, Axend=Axend)


def _fig_info_long(res, n_turns):
    P = len(res)
    # начальные x и y — СО ЗНАКОМ: у пучка «центр + спутники ±» модуль делает пары неразличимыми
    hdr = ["№", "оборотов", "нач. x, мм", "нач. y, мм", "макс |x|, мм", "ампл Ax нач", "ампл Ax кон", "Δампл, %"]
    col = [[str(i) for i in range(P)], [f"{r['surv_turns']}/{n_turns}" for r in res],
           [f"{r['x0']:+.3f}" for r in res], [f"{r.get('y0', 0.0):+.3f}" for r in res],
           [f"{r['xmax']:.3f}" for r in res],
           [f"{r['Ax0']:.3f}" for r in res], [f"{r['Axend']:.3f}" for r in res],
           [f"{(r['Axend']/r['Ax0']-1)*100:+.1f}" if r['Ax0'] > 0 else "—" for r in res]]
    fig = go.Figure(go.Table(columnwidth=[26, 60, 62, 62, 65, 65, 65, 60],
        header=dict(values=hdr, fill_color="#4e79a7", font=dict(color="white", size=12), align="center", height=26),
        cells=dict(values=col, align="center", height=22,
                   fill_color=[["#f4f6f8" if i % 2 else "white" for i in range(P)]])))
    fig.update_layout(margin=dict(l=6, r=6, t=6, b=6))
    return fig


def _fig_tbt(res, comp, title, unit):
    """Turn-by-turn: положение (x или y) у плоскости старта по номеру оборота."""
    fig = go.Figure()
    for i, r in enumerate(res):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=r["turns"].astype(np.float32), y=r[comp].astype(np.float32), mode="lines",
                                   name=f"ч{i}", line=dict(color=c, width=1.0)))
    fig.update_layout(title=dict(text=title, x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=40, b=72), plot_bgcolor="white", hovermode="closest",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="номер оборота")
    fig.update_yaxes(title_text=unit, zeroline=True, zerolinecolor="#bbb")
    return _grid(fig)


def _fig_envelope(res):
    """Огибающая амплитуды Ax/Ay по оборотам — прямой индикатор затухания (спад) или роста."""
    fig = go.Figure()
    for i, r in enumerate(res):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=r["turns"].astype(np.float32), y=r["Ax"].astype(np.float32), mode="lines",
                                   name=f"ч{i} · Ax", line=dict(color=c, width=1.6)))
        fig.add_trace(go.Scattergl(x=r["turns"].astype(np.float32), y=r["Ay"].astype(np.float32), mode="lines",
                                   name=f"ч{i} · Ay", line=dict(color=c, width=1.1, dash="dot")))
    fig.update_layout(title=dict(text="Огибающая амплитуды по оборотам  (сплошная — Ax, пунктир — Ay)",
                                 x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=40, b=72), plot_bgcolor="white", hovermode="closest",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="номер оборота")
    fig.update_yaxes(title_text="амплитуда (норм.), мм", rangemode="tozero")
    return _grid(fig)


def _fig_phase(res):
    """Фазовый портрет (x, x') у плоскости. Одна частица — цвет по обороту (эллипс=сохранение, спираль=затухание)."""
    fig = go.Figure()
    single = len(res) == 1
    for i, r in enumerate(res):
        if single:
            fig.add_trace(go.Scattergl(x=r["x"].astype(np.float32), y=r["xp"].astype(np.float32), mode="markers",
                                       marker=dict(size=2.6, color=r["turns"], colorscale="Viridis", showscale=True,
                                                   colorbar=dict(title="оборот", thickness=12, len=0.9)),
                                       name="(x, x')", hoverinfo="skip"))
        else:
            c = COLORS[i % len(COLORS)]
            fig.add_trace(go.Scattergl(x=r["x"].astype(np.float32), y=r["xp"].astype(np.float32), mode="markers",
                                       marker=dict(size=2.2, color=c), name=f"ч{i}", hoverinfo="skip"))
    fig.update_layout(title=dict(text="Фазовый портрет (x, x') в плоскости старта  (эллипс — сохранение, спираль — затухание/рост)",
                                 x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=40, b=72), plot_bgcolor="white", hovermode="closest",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="x, мм", zeroline=True, zerolinecolor="#bbb")
    fig.update_yaxes(title_text="x', мрад", zeroline=True, zerolinecolor="#bbb")
    return _grid(fig, legend=False)                          # фазовый портрет: легенда-эквивалент = цветовая шкала


# JS: у каждой панели СВОЯ высота; ручка под панелью тянет только её (страница прокручивается).
_RESIZE_JS = r"""
(function(){
 var IDS=__IDS__, DEF=__DEF__, MIN=140;
 function pdiv(id){return document.querySelector('#'+id+' .plotly-graph-div');}
 function rs(id){var d=pdiv(id); if(d&&window.Plotly){try{Plotly.Plots.resize(d);}catch(e){}}}
 function rsAll(){IDS.forEach(rs);}
 function setH(id,h){var e=document.getElementById(id); if(!e)return;
  e.style.height=Math.max(MIN,Math.round(h))+'px'; rs(id);}
 function reset(){IDS.forEach(function(id,i){setH(id,DEF[i]);});}
 function init(){reset(); setTimeout(rsAll,250); setTimeout(rsAll,800);}   // добить, если plotly ещё грузился
 var hs=document.querySelectorAll('.handle');
 for(var k=0;k<hs.length;k++){(function(h){
  var id=h.getAttribute('data-a'), idx=IDS.indexOf(id);
  h.addEventListener('mousedown',function(e){
   e.preventDefault();
   var a=document.getElementById(id), y0=e.clientY, h0=a.getBoundingClientRect().height;
   document.body.style.cursor='row-resize';
   function mv(ev){setH(id,h0+(ev.clientY-y0));}
   function up(){document.body.style.cursor='';
    document.removeEventListener('mousemove',mv);document.removeEventListener('mouseup',up);rs(id);}
   document.addEventListener('mousemove',mv);document.addEventListener('mouseup',up);
  });
  h.addEventListener('dblclick',function(){setH(id,DEF[idx]);});          // двойной клик — исходная высота
 })(hs[k]);}
 var btn=document.getElementById('reset'); if(btn)btn.addEventListener('click',reset);
 var rt=null;
 window.addEventListener('resize',function(){                             // окно шире/уже — только перерисовка
  clearTimeout(rt); rt=setTimeout(rsAll,120);                             // высоты НЕ трогаем: их задал человек
 });
 if(document.readyState!=='loading')init(); else document.addEventListener('DOMContentLoaded',init);
})();
"""

_PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>__TITLE_ATTR__</title>
<style>
html{scroll-behavior:smooth;}
body{margin:0;font-family:sans-serif;background:#fff;}
#hdr{position:sticky;top:0;z-index:20;background:#fff;padding:8px 16px 7px;text-align:center;
 border-bottom:1px solid #e3e6ea;box-shadow:0 1px 4px rgba(0,0,0,.06);}
#hdr .t{font-size:18px;color:#2a3f5f;font-weight:600;}
#hdr .s{font-size:12px;color:#666;margin-top:2px;}
#hdr .hint{font-size:11px;color:#8a96a3;margin-top:3px;}
#reset{font:inherit;font-size:11px;color:#41597d;background:#eef1f5;border:1px solid #d5dae1;
 border-radius:4px;padding:1px 8px;margin-left:8px;cursor:pointer;}
#reset:hover{background:#e2e7ee;}
#wrap{padding-bottom:48px;}
.panel{position:relative;overflow:hidden;}
.panel>.plotly-graph-div{width:100%;height:100%;}
.handle{height:14px;background:#eef0f2;border-top:1px solid #dce0e5;border-bottom:1px solid #dce0e5;
 cursor:row-resize;user-select:none;}
.handle:hover{background:#ccd3da;}
.handle::after{content:"\2261";display:block;text-align:center;color:#8a96a3;font-size:11px;line-height:13px;}
</style></head><body>
<div id="hdr"><div class="t">__TITLE__</div><div class="s">__SUBTITLE__</div></div>
<div id="wrap">
__WRAP__
</div>
<script>__JS__</script>
</body></html>"""


def _table_h(n_rows):
    """Высота панели-таблицы под её содержимое: шапка 26 + строки по 22 + поля."""
    return int(26 + 22 * max(1, n_rows) + 34)


def _rad_str(d):
    """Поле «фотонное излучение» для подзаголовка: учтён ли SR/ISR + режим ВЧ + проектная энергия SR."""
    if not d.get("sr"):
        return "фотонное излучение: не учитывается"
    tags = ["SR"]
    if d.get("isr"):
        tags.append("ISR")
    rf = d.get("rf") or "none"
    if rf != "none":
        tags.append(f"ВЧ {rf}")
    e = d.get("design_mev")
    etxt = f", {e / 1000:.1f} ГэВ" if e else ""
    return f"фотонное излучение: учтено ({' + '.join(tags)}{etxt})"


def build_html(d, res):
    """Самодостаточный HTML: заголовок + 4 панели (таблица · x(S) · y(S) · спектр) с перетаскиваемыми
    разделителями высоты. Каждая панель — отдельный plotly-график, растягивается под контейнер. Для
    многооборотного опыта x(S)/y(S) — ПОЛНАЯ траектория за все обороты (склейка), спектр уточняет тюн."""
    import json
    parts, closed, ring, nt = d["parts"], d["closed"], d["ring_len"], d["n_turns"]
    P = len(res); surv = sum(r["survived"] for r in res)
    mean_all = float(np.mean([r["mean_mm"] for r in res])); max_all = max(r["max_mm"] for r in res)
    kind = "замкнутое кольцо" if closed else "открытый сегмент"
    title = f"SynchrotronSim — анализ траектории, эксперимент {d['lattice']}·{d['group']} ({kind})"
    if nt > 1:
        subtitle = (f"{nt} оборотов · выжило {surv}/{P} за все {nt} об · кольцо {ring:.1f} м · частиц {P} · "
                    f"ср.откл. {mean_all:.3f} мм · макс {max_all:.3f} мм · ПОЛНАЯ траектория (склейка {nt} оборотов) · {_rad_str(d)}")
    else:
        subtitle = (f"оборот {ring:.1f} м · элементов {len(d['elem_names'])-1} · частиц {P} · выжило {surv}/{P} · "
                    f"ср.откл. {mean_all:.3f} мм · макс {max_all:.3f} мм · отклонение от опорной орбиты · {_rad_str(d)}")
    xt = f" ({nt} оборотов)" if nt > 1 else ""
    figs = [(_fig_info(res, nt), _table_h(P)),
            (_fig_dev(parts, "x", f"Отклонение по X от опорной орбиты: x(S) по частицам{xt}"), 400),
            (_fig_dev(parts, "y", f"Отклонение по Y от опорной орбиты: y(S) по частицам{xt}"), 400),
            (_fig_spectrum(res, ring), 460)]
    return _assemble(d["lattice"], d["group"], title, subtitle, figs)


def _assemble(lattice, group, title, subtitle, figs):
    """Собрать самодостаточный HTML из списка (fig, высота_px): панели идут ОДНА ПОД ДРУГОЙ, страница
    прокручивается, под каждой — ручка, меняющая высоту ТОЛЬКО этой панели. plotly.js встроен ОДИН раз."""
    import json
    heights = [int(h) for _, h in figs]
    to = dict(full_html=False, default_width="100%", default_height="100%", config={"responsive": True})
    panels = [fig.to_html(include_plotlyjs=("inline" if i == 0 else False), **to) for i, (fig, _) in enumerate(figs)]
    wrap = []
    for i, ph in enumerate(panels):                              # панель + ручка высоты под ней
        wrap.append(f'<div class="panel" id="p{i}" style="height:{heights[i]}px">{ph}</div>')
        wrap.append(f'<div class="handle" data-a="p{i}" title="тяните — изменить высоту этого графика; '
                    f'двойной клик — вернуть исходную"></div>')
    js = (_RESIZE_JS.replace("__IDS__", json.dumps([f"p{i}" for i in range(len(panels))]))
          .replace("__DEF__", json.dumps(heights)))
    return (_PAGE.replace("__JS__", js).replace("__TITLE_ATTR__", f"Анализ траектории {lattice}·{group}")
            .replace("__TITLE__", title).replace("__SUBTITLE__", subtitle).replace("__WRAP__", "\n".join(wrap)))


def build_long_html(d, res):
    """ДЛИННЫЙ отчёт (много оборотов, turn-by-turn): таблица · x(оборот) · огибающая амплитуды · фазовый портрет.
    Отвечает на вопрос «угасают ли колебания»: спад огибающей / спираль фазового портрета = затухание."""
    nt = d["n_turns"]; P = len(res); surv = sum(r["survived"] for r in res)
    dA = [(r["Axend"] / r["Ax0"] - 1) * 100 for r in res if r["Ax0"] > 0]
    dA_txt = f"Δампл. {min(dA):+.1f}…{max(dA):+.1f}%" if dA else "Δампл. —"
    kind = "замкнутое кольцо" if d["closed"] else "открытый сегмент"
    title = f"SynchrotronSim — анализ траектории, эксперимент {d['lattice']}·{d['group']} ({kind})"
    subtitle = (f"{nt} оборотов · выжило {surv}/{P} за все {nt} об · кольцо {d['ring_len']:.1f} м · "
                f"частиц {P} · {dA_txt} (спад = затухание) · {_rad_str(d)}")
    figs = [(_fig_info_long(res, nt), _table_h(P)),
            (_fig_tbt(res, "x", "Положение в плоскости старта по оборотам", "x, мм"), 380),
            (_fig_tbt(res, "y", "Положение в плоскости старта по оборотам", "y, мм"), 380),
            (_fig_envelope(res), 440),
            (_fig_phase(res), 560)]
    return _assemble(d["lattice"], d["group"], title, subtitle, figs)


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
    nt = d["n_turns"]; long_mode = d.get("mode") == "long"
    out = a.out or os.path.join(RESULTS, f"traj_analysis_exp{a.exp}.html")
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)

    if long_mode:                                           # === ДЛИННЫЙ режим: огибающая / фазовый портрет ===
        res = [metrics_long(p) for p in d["parts"]]
        surv = sum(r["survived"] for r in res)
        print(f"=== анализ траектории — эксперимент {a.exp} (решётка {d['lattice']}, {nt} оборотов) ===")
        print(f"кольцо {d['ring_len']:.2f}м, частиц {len(res)}, выжило {surv}/{len(res)} · {_rad_str(d)}")
        print("\n№  оборотов    нач.x,мм  нач.y,мм  макс|x|,мм  ампл_нач  ампл_кон  Δампл,%")
        for i, r in enumerate(res):
            dA = f"{(r['Axend']/r['Ax0']-1)*100:+.1f}" if r["Ax0"] > 0 else "—"
            print(f"{i:<2} {r['surv_turns']:>5}/{nt:<6} {r['x0']:>+8.3f}  {r['y0']:>+8.3f}   {r['xmax']:>8.3f}   "
                  f"{r['Ax0']:>7.3f}  {r['Axend']:>7.3f}   {dA:>6}")
        with open(out, "w", encoding="utf-8") as f:
            f.write(build_long_html(d, res))
        print(f"\n[окно] {out}  (страница прокручивается, высота графиков тянется; фазовый портрет — эллипс=сохранение, спираль=затухание)")
        if not a.no_open:
            webbrowser.open("file:///" + os.path.abspath(out).replace(os.sep, "/"))
        return

    res = [metrics(p, d["closed"], nt) for p in d["parts"]]
    # --- краткая сводка в консоль ---
    surv = sum(r["survived"] for r in res)
    print(f"=== Анализ траектории — эксперимент {a.exp} (решётка {d['lattice']}, пучок {d['group']}, "
          f"{'замкнуто' if d['closed'] else 'открыто'}, оборотов {nt}) ===")
    print(f"кольцо {d['ring_len']:.2f}м, элементов {len(d['elem_names'])-1}, частиц {len(res)}, выжило {surv}/{len(res)}")
    if nt > 1:
        print("\n№  оборотов   нач.x,мм  нач.y,мм  ср.|откл|,мм  макс,мм   тюн Qx   тюн Qy")
        for i, r in enumerate(res):
            print(f"{i:<2} {r['surv_turns']:>4}/{nt:<5} {r['x0']:>+8.3f}  {r['y0']:>+8.3f}  "
                  f"{r['mean_mm']:>9.3f}  {r['max_mm']:>7.3f}   {r['dqx']:>6.2f}   {r['dqy']:>6.2f}")
    else:
        print("\n№  дошла  нач.x,мм  нач.y,мм  ср.|откл|,мм  макс,мм  мин,мм   дом.Qx  дом.Qy  замык.,мм")
        for i, r in enumerate(res):
            clo = "—" if r["clo"] is None else f"{r['clo']:.3f}"
            print(f"{i:<2} {'да ' if r['survived'] else 'НЕТ'}  {r['x0']:>+8.3f}  {r['y0']:>+8.3f}  "
                  f"{r['mean_mm']:>9.3f}  {r['max_mm']:>7.3f}  "
                  f"{r['min_mm']:>6.3f}   {r['dqx']:>5.0f}  {r['dqy']:>6.0f}   {clo:>8}")

    # --- интерактивное окно (HTML) ---
    with open(out, "w", encoding="utf-8") as f:
        f.write(build_html(d, res))
    print(f"\n[окно] {out}  (страница прокручивается; полоска под графиком тянет его высоту)")
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
