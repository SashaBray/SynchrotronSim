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


def _fig_info(res):
    """Панель 1: таблица общей информации по частицам (fill'ит свой контейнер)."""
    P = len(res)
    hdr = ["№", "дошла", "ср.|откл|, мм", "макс, мм", "мин, мм", "дом. Qx", "дом. Qy", "замык. |Δ|, мм"]
    col = [[str(i) for i in range(P)],
           ["да" if r["survived"] else "НЕТ" for r in res],
           [f"{r['mean_mm']:.3f}" for r in res], [f"{r['max_mm']:.3f}" for r in res],
           [f"{r['min_mm']:.3f}" for r in res], [str(r["dqx"]) for r in res], [str(r["dqy"]) for r in res],
           [("—" if r["clo"] is None else f"{r['clo']:.3f}") for r in res]]
    fig = go.Figure(go.Table(
        columnwidth=[26, 40, 70, 55, 55, 50, 50, 75],
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
        fig.add_trace(go.Scattergl(x=p["S"], y=p[comp] * 1e3, mode="lines", name=f"ч{i}",
                                   line=dict(color=c, width=1.5),
                                   hovertemplate="S=%{x:.1f} м<br>" + comp + "=%{y:.3f} мм<extra>ч" + str(i) + "</extra>"))
    fig.update_layout(title=dict(text=title, x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=34, b=42), plot_bgcolor="white", hovermode="closest",
                      legend=dict(orientation="h", y=1.0, x=1.0, xanchor="right", yanchor="bottom"))
    fig.update_xaxes(title_text="S, м")
    fig.update_yaxes(title_text=f"отклонение {comp}, мм", zeroline=True, zerolinecolor="#bbb")
    return fig


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

    def env(key):
        L = min(len(r[key]) for r in res); e = np.zeros(L)
        for r in res:
            e = np.maximum(e, np.asarray(r[key])[:L])
        return e

    ex, ey = env("ax"), env("ay")
    gmax = max(float(ex.max()), float(ey.max()), 1e-30)

    def local(e, comp):
        return [(j, float(e[j]), comp) for j in range(1, len(e) - 1)
                if e[j] > e[j - 1] and e[j] >= e[j + 1] and e[j] >= rel * gmax]

    cand = sorted(local(ex, "x") + local(ey, "y"), key=lambda t: -t[1])
    sel = []
    for q, a, comp in cand:
        if all(not (comp == s[2] and abs(q - s[0]) < sep) for s in sel):   # раздвижка только внутри x или y
            sel.append((q, a, comp))
        if len(sel) >= n:
            break
    return sel


def _fig_spectrum(res, ring_len):
    """Спектр Фурье: амплитуда по частоте (число колебаний за оборот), x — сплошная, y — пунктир; легенда по
    частицам. У самых явных пиков — вертикальный пунктир, ранг (#1..) и частота в Гц (Q·c/длина кольца)."""
    fig = go.Figure(); qmax = 0
    for i, r in enumerate(res):
        c = COLORS[i % len(COLORS)]
        fig.add_trace(go.Scattergl(x=r["kx"], y=r["ax"], mode="lines", name=f"ч{i}", legendgroup=f"ч{i}",
                                   line=dict(color=c, width=1.4)))                        # сплошная = x (в легенде)
        fig.add_trace(go.Scattergl(x=r["ky"], y=r["ay"], mode="lines", name=f"ч{i} · y", legendgroup=f"ч{i}",
                                   line=dict(color=c, width=1.1, dash="dot"), showlegend=False))   # пунктир = y
        qmax = max(qmax, r["dqx"], r["dqy"])
    peaks = _spectrum_peaks(res)                                 # разметка самых явных частот
    lab = sorted([(q, a, comp, rk) for rk, (q, a, comp) in enumerate(peaks, 1)], key=lambda t: t[0])
    prev_q = -1e9; row = 0                                       # близкие по Q подписи разносим по высоте (цикл 0/28/56)
    for q, a, comp, rk in lab:
        row = (row + 1) if (q - prev_q) < 2.0 else 0
        extra = (row % 3) * 28
        fig.add_vline(x=q, line=dict(color="#6b6b6b", width=1, dash="dash"))
        hz = _hz(q * LIGHT_C / ring_len) if ring_len > 0 else "—"
        fig.add_annotation(x=q, y=a, yanchor="bottom", yshift=6 + extra, showarrow=False, align="center",
                           text=f"<b>#{rk}</b> · Q{comp}={q}<br>{hz}",
                           font=dict(size=10, color="#333"), bgcolor="rgba(255,255,255,0.8)",
                           bordercolor="#ccc", borderwidth=1, borderpad=1)
        prev_q = q; qmax = max(qmax, q)
    fig.update_layout(title=dict(text="Спектр колебаний (Фурье): амплитуда по Q  (сплошная — x, пунктир — y)",
                                 x=0.5, xanchor="center", y=0.97, yanchor="top", font=dict(size=14)),
                      margin=dict(l=62, r=18, t=48, b=44), plot_bgcolor="white",
                      legend=dict(orientation="h", y=1.0, x=1.0, xanchor="right", yanchor="bottom"))
    fig.update_xaxes(title_text="Q — число колебаний за оборот", range=[0, max(40, qmax * 1.5)])
    ytop = peaks[0][1] * 1.6 if peaks else None                  # запас сверху под подписи (в т.ч. разнесённые)
    fig.update_yaxes(title_text="амплитуда, мм", range=([0, ytop] if ytop else None))
    return fig


# JS: перетаскиваемые разделители между панелями (row-resize) + reflow графиков (Plotly.Plots.resize).
_RESIZE_JS = r"""
(function(){
 var SP=10, PROP=[0.15,0.26,0.26,0.33], IDS=['p0','p1','p2','p3'];
 function pdiv(id){return document.querySelector('#'+id+' .plotly-graph-div');}
 function rs(id){var d=pdiv(id); if(d&&window.Plotly){try{Plotly.Plots.resize(d);}catch(e){}}}
 function rsAll(){IDS.forEach(rs);}
 function avail(){var w=document.getElementById('wrap'); return w?w.clientHeight-(IDS.length-1)*SP:0;}
 function init(){
  var w=document.getElementById('wrap'); if(!w||!w.clientHeight){setTimeout(init,60);return;}
  var H=avail();
  IDS.forEach(function(id,i){document.getElementById(id).style.height=Math.max(60,H*PROP[i])+'px';});
  rsAll(); setTimeout(rsAll,250); setTimeout(rsAll,700);        // добить, если plotly ещё грузился
 }
 var sps=document.querySelectorAll('.splitter');
 for(var k=0;k<sps.length;k++){(function(s){
  s.addEventListener('mousedown',function(e){
   e.preventDefault();
   var a=document.getElementById(s.getAttribute('data-a')), b=document.getElementById(s.getAttribute('data-b'));
   var y0=e.clientY, ha=a.getBoundingClientRect().height, hb=b.getBoundingClientRect().height;
   function mv(ev){var dy=ev.clientY-y0, na=ha+dy, nb=hb-dy; if(na<60||nb<60)return;
    a.style.height=na+'px'; b.style.height=nb+'px'; rs(a.id); rs(b.id);}
   function up(){document.removeEventListener('mousemove',mv);document.removeEventListener('mouseup',up);rsAll();}
   document.addEventListener('mousemove',mv);document.addEventListener('mouseup',up);
  });
 })(sps[k]);}
 var rt=null;
 window.addEventListener('resize',function(){                   // окно меняется → масштабируем панели пропорционально
  clearTimeout(rt); rt=setTimeout(function(){
   var w=document.getElementById('wrap'); if(!w)return;
   var cur=IDS.map(function(id){return document.getElementById(id).getBoundingClientRect().height;});
   var sum=cur.reduce(function(x,y){return x+y;},0); if(sum<=0){init();return;}
   var H=avail();
   IDS.forEach(function(id,i){document.getElementById(id).style.height=Math.max(60,H*cur[i]/sum)+'px';});
   rsAll();
  },120);
 });
 if(document.readyState!=='loading')init(); else document.addEventListener('DOMContentLoaded',init);
})();
"""

_PAGE = r"""<!doctype html><html><head><meta charset="utf-8"><title>__TITLE_ATTR__</title>
<style>
html,body{margin:0;height:100%;}
body{display:flex;flex-direction:column;font-family:sans-serif;}
#hdr{padding:6px 16px 5px;text-align:center;border-bottom:1px solid #eee;flex:0 0 auto;}
#hdr .t{font-size:18px;color:#2a3f5f;font-weight:600;}
#hdr .s{font-size:12px;color:#666;margin-top:2px;}
#wrap{flex:1 1 auto;display:flex;flex-direction:column;min-height:0;}
.panel{overflow:hidden;position:relative;}
.panel>.plotly-graph-div{width:100%;height:100%;}
.splitter{height:10px;background:#eef0f2;border-top:1px solid #dce0e5;border-bottom:1px solid #dce0e5;
 cursor:row-resize;flex:0 0 auto;}
.splitter:hover{background:#ccd3da;}
.splitter::after{content:"\2261";display:block;text-align:center;color:#8a96a3;font-size:11px;line-height:9px;}
</style></head><body>
<div id="hdr"><div class="t">__TITLE__</div><div class="s">__SUBTITLE__</div></div>
<div id="wrap">
<div class="panel" id="p0">__INFO__</div>
<div class="splitter" data-a="p0" data-b="p1" title="потяните, чтобы изменить высоту"></div>
<div class="panel" id="p1">__DEVX__</div>
<div class="splitter" data-a="p1" data-b="p2" title="потяните, чтобы изменить высоту"></div>
<div class="panel" id="p2">__DEVY__</div>
<div class="splitter" data-a="p2" data-b="p3" title="потяните, чтобы изменить высоту"></div>
<div class="panel" id="p3">__SPEC__</div>
</div>
<script>__JS__</script>
</body></html>"""


def build_html(d, res):
    """Самодостаточный HTML: заголовок + 3 панели (таблица/координаты/спектр) с перетаскиваемыми
    разделителями высоты. Каждая панель — отдельный plotly-график, растягивается под свой контейнер."""
    parts, closed, ring = d["parts"], d["closed"], d["ring_len"]
    P = len(res); surv = sum(r["survived"] for r in res)
    mean_all = float(np.mean([r["mean_mm"] for r in res])); max_all = max(r["max_mm"] for r in res)
    kind = "замкнутое кольцо" if closed else "открытый сегмент"
    title = f"SynchrotronSim — анализ траектории, эксперимент {d['lattice']}·{d['group']} ({kind})"
    subtitle = (f"оборот {ring:.1f} м · элементов {len(d['elem_names'])-1} · частиц {P} · выжило {surv}/{P} · "
                f"ср.откл. {mean_all:.3f} мм · макс {max_all:.3f} мм · локальные координаты (отклонение от опорной орбиты)")
    to = dict(full_html=False, default_width="100%", default_height="100%", config={"responsive": True})
    info = _fig_info(res).to_html(include_plotlyjs="inline", **to)      # plotly.js встроен ОДИН раз (первая панель)
    devx = _fig_dev(parts, "x", "Отклонение по X от опорной орбиты: x(S) по частицам").to_html(include_plotlyjs=False, **to)
    devy = _fig_dev(parts, "y", "Отклонение по Y от опорной орбиты: y(S) по частицам").to_html(include_plotlyjs=False, **to)
    spec = _fig_spectrum(res, ring).to_html(include_plotlyjs=False, **to)
    return (_PAGE.replace("__JS__", _RESIZE_JS)
            .replace("__TITLE_ATTR__", f"Анализ траектории {d['lattice']}·{d['group']}")
            .replace("__TITLE__", title).replace("__SUBTITLE__", subtitle)
            .replace("__INFO__", info).replace("__DEVX__", devx).replace("__DEVY__", devy).replace("__SPEC__", spec))


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
    with open(out, "w", encoding="utf-8") as f:
        f.write(build_html(d, res))
    print(f"\n[окно] {out}  (панели тянутся мышкой за разделители)")
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
