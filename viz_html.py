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
DEV_COL = ["#555555", "#e15759", "#f28e2b", "#4e79a7", "#59a14f"]   # цвета пробного пучка (общие 3D↔2D):
#                                                                     опорная серая, ±x тёплые, ±y холодные
KINDCOL = {"quad": "#4e79a7", "dipole": "#e15759", "sext": "#59a14f", "oct": "#b07aa1",
           "dq": "#f28e2b", "kicker": "#9c755f", "monitor": "#bab0ac"}
TRAJ_MAX = 8000                                              # макс. точек на частицу (длинную траекторию прореживаем)
TARGET_TRI = 650_000                                        # цель по треуг. авто-прореживания (float32 ⇒ вдвое больше
#                                                             деталей при прежнем весе; μm-точность на масштабе метров)

# JS-связка 2D↔3D: наведение на график отклонения → временная подсветка точки на 3D-траектории (снимается
# при уходе мыши); клик → стойкая подсветка (повторный клик по той же точке снимает, клик по другой — переносит).
# Точка 2D несёт глоб. XYZ в customdata[2..4]; d[0]=3D-сцена, d[1]=2D-график (порядок в теле HTML).
_LINK_JS = """<script>
(function(){
 function wire(){
  var d=document.querySelectorAll('.plotly-graph-div');
  if(d.length<2||!window.Plotly){setTimeout(wire,150);return;}
  var g3=d[0],g2=d[1];
  if(!g3.data||!g2.on){setTimeout(wire,150);return;}
  var hi=g3.data.length,ci=g3.data.length+1;
  Plotly.addTraces(g3,[
   {type:'scatter3d',x:[],y:[],z:[],mode:'markers',hoverinfo:'skip',showlegend:false,marker:{size:7,color:'#111',line:{color:'#fff',width:2}}},
   {type:'scatter3d',x:[],y:[],z:[],mode:'markers',hoverinfo:'skip',showlegend:false,marker:{size:10,color:'#ffd400',line:{color:'#111',width:2}}}
  ]);
  var ck=null;
  function xyz(pt){var c=pt.customdata;return c?[[c[2]],[c[3]],[c[4]]]:null;}
  g2.on('plotly_hover',function(e){var p=xyz(e.points[0]);if(p)Plotly.restyle(g3,{x:[p[0]],y:[p[1]],z:[p[2]]},[hi]);});
  g2.on('plotly_unhover',function(){Plotly.restyle(g3,{x:[[]],y:[[]],z:[[]]},[hi]);});
  g2.on('plotly_click',function(e){
   var pt=e.points[0],k=pt.curveNumber+':'+pt.pointIndex,p=xyz(pt);
   if(k===ck){Plotly.restyle(g3,{x:[[]],y:[[]],z:[[]]},[ci]);ck=null;}
   else if(p){Plotly.restyle(g3,{x:[p[0]],y:[p[1]],z:[p[2]]},[ci]);ck=k;}
  });
 }
 if(document.readyState!=='loading')wire();else document.addEventListener('DOMContentLoaded',wire);
})();
</script>"""


def _anim_html(traj_idx, traj_col):
    """Кометная анимация 3D: по каждой трассе траектории бежит частица с УГАСАЮЩИМ следом (маркеры,
    цвет→белый + размер→0 к хвосту). Кнопка «animate» ВКЛ/ВЫКЛ: по умолчанию статичный вид (полные треки),
    по клику — линии приглушаются (ghost) и пускается комета; повторный клик возвращает статику. Панель:
    animate, скорость. Работает на УЖЕ встроенных массивах траектории (g.data[i].x/y/z) — файл не растёт.
    traj_idx — индексы трасс, traj_col — их цвета."""
    import json
    return (r"""<div id="ssac" style="position:fixed;left:12px;bottom:12px;z-index:9999;font:13px sans-serif;
background:rgba(255,255,255,.9);border:1px solid #ccc;border-radius:8px;padding:8px 11px;box-shadow:0 1px 5px rgba(0,0,0,.15)">
<button id="ssap" style="cursor:pointer;border:1px solid #888;border-radius:5px;background:#fff;padding:3px 10px;font:13px sans-serif">&#9654; animate</button>
&nbsp;скорость <input id="ssas" type="range" min="1" max="40" value="8" style="vertical-align:middle;width:90px">
</div>
<script>
(function(){
 var TRAJ=__IDX__, COL=__COL__;
 function wire(){
  var d=document.querySelectorAll('.plotly-graph-div');
  if(!d.length||!window.Plotly){setTimeout(wire,160);return;}
  var g=d[0];
  if(!g.data||!g.on){setTimeout(wire,160);return;}
  function arr(a){                                             // массив / типизир. / base64-спека plotly typed array
   if(a&&a.length!==undefined)return a;                        // уже готов (plotly декодировал на newPlot)
   if(a&&a.bdata){var b=atob(a.bdata),L=b.length,u=new Uint8Array(L);for(var q=0;q<L;q++)u[q]=b.charCodeAt(q);
    var D=a.dtype,B=u.buffer;                                  // plotly пишет numpy-коды: f4/f8/i4/…
    if(D==='f4')return new Float32Array(B); if(D==='f8')return new Float64Array(B);
    if(D==='i4')return new Int32Array(B);   if(D==='u4')return new Uint32Array(B);
    if(D==='i2')return new Int16Array(B);   if(D==='u2')return new Uint16Array(B);
    if(D==='i1')return new Int8Array(B);    if(D==='u1'||D==='u1c')return new Uint8Array(B);}
   return a;
  }
  var paths=TRAJ.map(function(ix){var t=g.data[ix],x=arr(t.x);return{x:x,y:arr(t.y),z:arr(t.z),n:x.length};});
  var n=paths.length?paths[0].n:0;
  if(n<2)return;
  var TAIL=Math.min(300,Math.max(20,Math.round(n*0.1)));       // длина следа ~10% пути (голова→хвост)
  var base=g.data.length;
  Plotly.addTraces(g,TRAJ.map(function(){return{type:'scatter3d',mode:'markers',x:[],y:[],z:[],
    hoverinfo:'skip',showlegend:false,marker:{size:[],color:[]}};}));
  var comet=TRAJ.map(function(_,i){return base+i;});
  var rgb=COL.map(function(h){h=h.replace('#','');return[parseInt(h.substr(0,2),16),parseInt(h.substr(2,2),16),parseInt(h.substr(4,2),16)];});
  var head=TRAJ.map(function(){return 0;}), RUN=false, last=0;  // анимация ВЫКЛ по умолчанию (статичный вид)
  var sp=document.getElementById('ssas'), speed=Math.max(1,Math.round(n/480)); sp.value=speed;
  sp.oninput=function(){speed=+sp.value;};
  function empty(){return comet.map(function(){return[];});}
  var pb=document.getElementById('ssap');
  pb.onclick=function(){                                        // кнопка animate: вкл/выкл анимацию
   RUN=!RUN;
   if(RUN){pb.innerHTML='&#9209; стоп'; Plotly.restyle(g,{opacity:0.16},TRAJ);}  // приглушить линии, пустить комет
   else{pb.innerHTML='&#9654; animate';                                          // вернуть статичный вид (полные треки)
    for(var i=0;i<comet.length;i++)head[i]=0;
    Plotly.restyle(g,{x:empty(),y:empty(),z:empty()},comet);
    Plotly.restyle(g,{opacity:1},TRAJ);}
  };
  function step(ts){
   requestAnimationFrame(step);
   if(!RUN||ts-last<33)return;                                 // ~30 к/с
   last=ts;
   var xs=[],ys=[],zs=[],cs=[],ss=[];
   for(var i=0;i<comet.length;i++){
    var p=paths[i],h=head[i],lo=Math.max(0,h-TAIL+1),r=rgb[i];
    var wx=[],wy=[],wz=[],wc=[],ws=[];
    for(var k=lo;k<=h;k++){
     var a=1-(h-k)/(TAIL-1), t=1-a;                            // a: 1 голова .. 0 хвост; t=1-a
     wx.push(p.x[k]);wy.push(p.y[k]);wz.push(p.z[k]);
     wc.push('rgba('+Math.round(r[0]+(255-r[0])*t)+','+Math.round(r[1]+(255-r[1])*t)+','+Math.round(r[2]+(255-r[2])*t)+','+(0.2+0.8*a).toFixed(3)+')');
     ws.push(2.5+7*a);                                         // размер 9.5 голова .. 2.5 хвост
    }
    xs.push(wx);ys.push(wy);zs.push(wz);cs.push(wc);ss.push(ws);
    head[i]=h+speed>=p.n?0:h+speed;                            // на конце — заново от плоскости инжекции
   }
   Plotly.restyle(g,{x:xs,y:ys,z:zs,'marker.color':cs,'marker.size':ss},comet);
  }
  requestAnimationFrame(step);
 }
 if(document.readyState!=='loading')wire();else document.addEventListener('DOMContentLoaded',wire);
})();
</script>""").replace("__IDX__", json.dumps(traj_idx)).replace("__COL__", json.dumps(traj_col))


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


def _deviation_fig(dev, labels, meta):
    """2D-график отклонения пробного пучка от опорной орбиты: x(s) сверху, y(s) снизу, линия на частицу.
    Ховер: координаты + № элемента по порядку + имя + суперпериод (из meta, общего для всех частиц)."""
    name, midx, msp = meta["name"], meta["idx"], meta["sp"]
    fig2 = make_subplots(rows=2, cols=1, shared_xaxes=True, vertical_spacing=0.1,
                         subplot_titles=("Горизонтальное отклонение x(s)", "Вертикальное отклонение y(s)"))
    for i, d in enumerate(dev):
        if len(d) < 2:
            continue
        L = len(d)
        s = d[:, 0].astype(np.float32); x = (d[:, 1] * 1e3).astype(np.float32); y = (d[:, 2] * 1e3).astype(np.float32)
        nm = labels[i] if i < len(labels) else f"P{i+1}"
        c = DEV_COL[i % len(DEV_COL)]; w = 2.0 if i == 0 else 1.1
        # %{text}=имя, customdata[0]=№, [1]=СП; [2..4]=глоб. XYZ (для JS-связки клик/наведение → 3D). float32:
        # −43% веса; №/СП стали float ⇒ формат :.0f. XYZ в μm-точности хватает для подсветки.
        text = name[:L]; cust = np.column_stack([midx[:L], msp[:L], d[:, 3], d[:, 4], d[:, 5]]).astype(np.float32)
        ht = ("s = %{x:.2f} м · {ax} = %{y:+.3f} мм<br>элемент №%{customdata[0]:.0f}: %{text}"
              "<br>суперпериод %{customdata[1]:.0f}<extra>%{fullData.name}</extra>")
        fig2.add_trace(go.Scattergl(x=s, y=x, mode="lines", line=dict(width=w, color=c), name=nm, legendgroup=nm,
                                    text=text, customdata=cust, hovertemplate=ht.replace("{ax}", "x")), row=1, col=1)
        fig2.add_trace(go.Scattergl(x=s, y=y, mode="lines", line=dict(width=w, color=c), name=nm, legendgroup=nm,
                                    showlegend=False, text=text, customdata=cust,
                                    hovertemplate=ht.replace("{ax}", "y")), row=2, col=1)
    fig2.update_xaxes(title_text="s, м", row=2, col=1)
    fig2.update_yaxes(title_text="x, мм", row=1, col=1); fig2.update_yaxes(title_text="y, мм", row=2, col=1)
    fig2.update_xaxes(gridcolor="#eee", zeroline=False)
    fig2.update_yaxes(gridcolor="#eee", zeroline=True, zerolinecolor="#bbb", zerolinewidth=1)
    fig2.update_layout(title=dict(text="Отклонение пробного пучка от опорной орбиты (±1мм по x и y)",
                                  x=0.5, xanchor="center", y=0.99, yanchor="top"),
                       margin=dict(l=65, r=135, t=64, b=55), plot_bgcolor="white", hovermode="closest",
                       legend=dict(orientation="v", x=1.005, y=1.0, xanchor="left", yanchor="top",
                                   title="частица", bgcolor="rgba(255,255,255,0.7)", bordercolor="#ddd", borderwidth=1))
    return fig2


def build(exp_id, lat, group, hw, kind, lc, decimate, glyph, out):
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

    # ПРОБНЫЙ пучок (сектор): та же перетрекованная 3D-траектория и 2D-отклонение (один пучок ⇒ связка
    # клик↔3D работает) + метаданные точек (имя/№/СП, глоб. XYZ) и магнитов.
    dev = dlabels = dmeta = dmags = None
    try:
        dev, dlabels, dmeta, dmags = mte.deviation_paths(int(lat), int(group))
    except Exception as e:
        print(f"[!] пробный пучок не собрался: {type(e).__name__}: {e}")

    fig = go.Figure()
    ntri = 0
    for dv, tlist in grouped.items():                        # магниты → Mesh3d, цвет по типу
        big = vz.build_combined_mesh(stl_cache[hw[dv]["stl_path"]], np.asarray(tlist, float))
        V = np.asarray(big.points, dtype=np.float32)         # float32: −43% веса, μm-точность (масштаб метров) — с запасом
        F = np.asarray(big.faces).reshape(-1, 4)[:, 1:].astype(np.int32)
        ntri += len(F)
        fig.add_trace(go.Mesh3d(x=V[:, 0], y=V[:, 1], z=V[:, 2], i=F[:, 0], j=F[:, 1], k=F[:, 2],
                                color=KINDCOL.get(kind.get(dv, ""), "silver"), opacity=0.5,
                                name=dv, hovertext=dv, hoverinfo="text", showlegend=False))

    # МАРКЕРЫ МАГНИТОВ в центрах — ховер «имя · № · суперпериод» (полная инфа по магниту в 3D)
    if dmags:
        mc = np.array(dmags["center"], np.float32)
        mh = [f"{dmags['name'][k]} · №{dmags['idx'][k]} · суперпериод {dmags['sp'][k]}" for k in range(len(mc))]
        fig.add_trace(go.Scatter3d(x=mc[:, 0], y=mc[:, 1], z=mc[:, 2], mode="markers",
                                   marker=dict(size=2.5, color="rgba(70,70,70,0.35)"), name="магниты",
                                   text=mh, hovertemplate="%{text}<extra>магнит</extra>", showlegend=False))
    # BPM — маркеры на датчиках положения (kind=monitor)
    bpm = np.array([(float(r["x"]), float(r["y"]), float(r["z"])) for r in lc if kind.get(r["device_type_id"]) == "monitor"], np.float32)
    if len(bpm):
        fig.add_trace(go.Scatter3d(x=bpm[:, 0], y=bpm[:, 1], z=bpm[:, 2], mode="markers",
                                   marker=dict(size=3, color="black", symbol="diamond"), name=f"BPM ({len(bpm)})"))

    # ТРАЕКТОРИЯ = пробный пучок (глоб. X,Y,Z из dev[:,3:6]); ховер «имя · № · суперпериод · s»
    traj_idx, traj_col = [], []                              # индексы/цвета трасс — для кометной анимации
    if dev:
        for i, d in enumerate(dev):
            if len(d) < 2:
                continue
            ht = [f"{dmeta['name'][k]} · №{dmeta['idx'][k]} · суперпериод {dmeta['sp'][k]}<br>s = {d[k, 0]:.1f} м"
                  for k in range(len(d))]
            col = DEV_COL[i % len(DEV_COL)]
            fig.add_trace(go.Scatter3d(x=d[:, 3].astype(np.float32), y=d[:, 4].astype(np.float32),
                                       z=d[:, 5].astype(np.float32), mode="lines",
                                       line=dict(width=(6 if i == 0 else 4), color=col),
                                       name=dlabels[i] if i < len(dlabels) else f"P{i+1}",
                                       text=ht, hovertemplate="%{text}<extra>%{fullData.name}</extra>"))
            traj_idx.append(len(fig.data) - 1); traj_col.append(col)

    # ПОДЗАГОЛОВОК: сколько оборотов, замкнута ли орбита, невязка (секторным телом — как график отклонения).
    sub = ""
    try:
        c = mte.experiment_closure(int(lat), int(group))
        st = f"Оборотов: {c['turns']} · "
        if c["closed"]:
            st += "орбита замкнута" + (f" (невязка {c['residual_mm']:.2f} мм)" if c["residual_mm"] is not None else "")
        else:
            st += "открытый сегмент (орбита не замкнута)"
        st += f" · выжило {c['survived']}/{c['n']}"
        sub = f"<br><sub>{st}</sub>"
    except Exception as e:                                        # подпись не критична — не рушим виз
        print(f"[!] сводка замыкания пропущена: {type(e).__name__}: {e}")

    mode = "блоки" if glyph else ("STL, прорежено " + f"{decimate:g}" if decimate > 0 else "полный STL")
    fig.update_layout(title=dict(text=f"SynchrotronSim — эксперимент {exp_id} "
                                      f"(решётка {lat}, {len(grouped)} типов магнитов; {mode}){sub}", y=0.96, yanchor="top"),
                      scene=dict(aspectmode="data", xaxis_title="X, м", yaxis_title="Y, м", zaxis_title="Z, м",
                                 bgcolor="white", domain=dict(x=[0, 1], y=[0, 1])),  # сцена на всю область
                      autosize=True,                          # тянется на всю ширину окна браузера (без пустого поля справа)
                      legend=dict(title="частицы", x=0.99, xanchor="right", y=0.99, yanchor="top",  # оверлей, не резервирует колонку
                                  bgcolor="rgba(255,255,255,0.6)", bordercolor="lightgray", borderwidth=1),
                      margin=dict(l=0, r=0, t=66, b=0))          # t↑ под двухстрочный заголовок
    # 2D-график отклонения (ТОТ ЖЕ пробный пучок) — под 3D; связь с 3D через customdata (глоб. XYZ точки)
    dev_html = ""
    if dev:
        try:
            dev_html = _deviation_fig(dev, dlabels, dmeta).to_html(
                include_plotlyjs=False, full_html=False, default_width="100%", default_height="640",
                config={"responsive": True})
        except Exception as e:                                    # не рушим 3D, если график не собрался
            print(f"[!] график отклонения пропущен: {type(e).__name__}: {e}")

    s3d = fig.to_html(include_plotlyjs="inline", full_html=False, default_width="100%",
                      default_height="78vh", config={"responsive": True})
    anim = _anim_html(traj_idx, traj_col) if traj_idx else ""    # кометная анимация частиц с угасающим следом
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    with open(out, "w", encoding="utf-8") as f:                   # 3D сверху, график отклонения снизу — один файл
        f.write("<!doctype html><html><head><meta charset='utf-8'>"
                f"<title>SynchrotronSim — эксперимент {exp_id}</title></head>"
                f"<body style='margin:0;padding-top:16px;font-family:sans-serif'>{s3d}{dev_html}"
                f"{_LINK_JS if dev_html else ''}{anim}</body></html>")
    print(f"[OK] {out}")
    print(f"     {os.path.getsize(out)/1e6:.1f} МБ ({mode}); {ntri:,} треуг. магнитов, "
          f"{len(dev) if dev else 0} частиц (пробный пучок){' + график отклонения' if dev_html else ''}. "
          f"Самодостаточный — крутить/зумить, передавать файлом.")
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
    lat = exp["lattice_id"]; group = exp["group_id"]
    hw = {r["device_type_id"]: r for r in _read("hardware_library.csv")}
    kind = {r["device_type_id"]: r["kind"] for r in _read("element_transforms.csv")}
    lc = [r for r in _read("lattice_configs.csv") if r["lattice_id"] == lat]

    base = out or os.path.join(HTML_DIR, f"exp{exp_id}_3d.html")
    stem = base[:-5] if base.lower().endswith(".html") else base
    if both:                                                 # обе версии: лёгкая (блоки) + детальная (STL)
        build(exp_id, lat, group, hw, kind, lc, 0.0, True, stem + "_light.html")
        build(exp_id, lat, group, hw, kind, lc, _resolve_decimate(decimate, False, lc, hw), False, stem + "_detailed.html")
    else:
        build(exp_id, lat, group, hw, kind, lc, _resolve_decimate(decimate, glyph, lc, hw), glyph, base)


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
