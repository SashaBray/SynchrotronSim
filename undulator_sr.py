# -*- coding: utf-8 -*-
"""Паспорт синхротронного излучения ОНДУЛЯТОРА — СИЛА.

Инструмент анализа СИ одного ондулятора: энергии гармоник E_n(K), кривые перестройки по зазору,
спектральный поток (на оси и в центральном конусе), яркость (свёртка с эмиттансом/бетой пучка),
полная мощность и плотность мощности, критическая энергия (вигглерный предел).

Формулы — стандартная теория ондуляторного излучения (Kim; X-ray Data Booklet). Параметры ПУЧКА
(ток, εy, σ_E, β в месте установки) — ТИПОВЫЕ ПРЕДПОЛОЖЕНИЯ для СИЛА-26, помечены ниже; подставляются
реальные из решётки/паспорта, когда появятся. Поле B0(зазор) — предполагаемое (якорь 1,1 Т @ 5 мм),
см. [[undulator-parameters]].

Запуск: py -3.14 undulator_sr.py                      # паспорт при рабочем зазоре 5 мм (K≈1.85)
        py -3.14 undulator_sr.py --gap 10             # при другом зазоре
        py -3.14 undulator_sr.py --b0 1.0 --emit-y 1.0  # переопределить допущения
"""
import os, sys, argparse, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy.special import jv

MC2 = 0.510998950e-3            # ГэВ, масса электрона

# ───── параметры (ПРЕДПОЛОЖЕНИЯ для СИЛА, помечены) ─────
class P:
    E = 6.0                     # ГэВ — энергия пучка
    I = 0.2                     # А   — ток (ПРЕДПОЛОЖЕНИЕ)
    emit_x = 42e-12            # м·рад — εx (СИЛА-26, реш.26)
    kappa = 0.01               # связь → εy = kappa·εx (ПРЕДПОЛОЖЕНИЕ)
    sig_E = 1.0e-3             # σ_E/E — энергетический разброс (ПРЕДПОЛОЖЕНИЕ)
    beta_x = 5.0               # м — β_x в прямолинейке ИД (ПРЕДПОЛОЖЕНИЕ)
    beta_y = 2.5               # м — β_y в прямолинейке ИД (ПРЕДПОЛОЖЕНИЕ)
    # ондулятор (in-vacuum CPMU, паспорт ЭСК 2180)
    lam_u = 0.018              # м — период
    N = 30                     # число периодов
    L = 0.55                   # м — длина структуры
    B0_ref = 1.1               # Т — поле при зазоре gap_ref (ПРЕДПОЛОЖЕНИЕ, якорь)
    gap_ref = 5.0              # мм
    gap_min, gap_max = 5.0, 20.0


def B0_of_gap(gap_mm, p=P):
    """Пиковое поле vs зазор: ролл-офф первой гармоники ~exp(-π·g/λ_u), якорь B0_ref при gap_ref."""
    C = p.B0_ref / math.exp(-math.pi * p.gap_ref / (p.lam_u * 1e3))
    return C * math.exp(-math.pi * gap_mm / (p.lam_u * 1e3))


def K_of_B0(B0, p=P):
    return 0.0934 * B0 * (p.lam_u * 1e3)       # K = 0.0934·B0[Т]·λ_u[мм]


def gamma(p=P):
    return p.E / MC2


def harmonic_energy_keV(n, K, theta=0.0, p=P):
    """Энергия n-й гармоники под углом theta (рад). На оси (theta=0) — нечётные n."""
    g = gamma(p)
    lam_n = (p.lam_u / (2 * n * g ** 2)) * (1 + K ** 2 / 2 + (g * theta) ** 2)   # м
    return 1.239841984e-6 / lam_n / 1e3        # кэВ (hc=1.2398 эВ·мкм)


def _Fn(n, K):
    """Функция ондулятора F_n(K) для потока на оси (нечётные n)."""
    if K == 0:
        return 0.0
    x = n * K ** 2 / (4 * (1 + K ** 2 / 2))
    a = (n - 1) // 2
    b = (n + 1) // 2
    jj = jv(a, x) - jv(b, x) if n % 2 == 1 else 0.0
    return (n ** 2 * K ** 2 / (1 + K ** 2 / 2) ** 2) * jj ** 2


def _Qn(n, K):
    """Функция для углоинтегрированного потока в центральном конусе."""
    if K == 0:
        return 0.0
    return (1 + K ** 2 / 2) / n * _Fn(n, K)


def _Gpow(K):
    """Функция плотности мощности на оси G(K) (Kim)."""
    if K == 0:
        return 0.0
    return K * (K ** 6 + 24.0 / 7 * K ** 4 + 4 * K ** 2 + 16.0 / 7) / (1 + K ** 2) ** 3.5


def flux_density_onaxis(n, K, p=P):
    """Спектральная плотность потока на оси, ф/с/мрад²/0.1%BW."""
    return 1.744e14 * p.N ** 2 * p.E ** 2 * p.I * _Fn(n, K)


def central_flux(n, K, p=P):
    """Поток в центральном конусе n-й гармоники, ф/с/0.1%BW."""
    return 1.431e14 * p.N * p.I * _Qn(n, K)


def sigma_photon(E_keV, p=P):
    """Дифракционные размер и расходимость фотонного пучка (м, рад)."""
    lam = 1.239841984e-6 / (E_keV * 1e3)       # м
    sr = math.sqrt(2 * lam * p.L) / (4 * math.pi)
    srp = math.sqrt(lam / (2 * p.L))
    return sr, srp


def brightness(n, K, p=P):
    """Спектральная яркость n-й гармоники, ф/с/мрад²/мм²/0.1%BW."""
    En = harmonic_energy_keV(n, K, 0.0, p)
    sr, srp = sigma_photon(En, p)
    ey = p.kappa * p.emit_x
    Sx = math.sqrt(p.emit_x * p.beta_x + sr ** 2)
    Sy = math.sqrt(ey * p.beta_y + sr ** 2)
    Sxp = math.sqrt(p.emit_x / p.beta_x + srp ** 2)
    Syp = math.sqrt(ey / p.beta_y + srp ** 2)
    flux = central_flux(n, K, p)               # ф/с/0.1%BW
    # ф/с/0.1%BW / (2π Σx Σx' · 2π ... ) → приводим к /мрад²/мм²
    denom = (2 * math.pi) ** 2 * Sx * Sxp * Sy * Syp
    B_SI = flux / denom                        # ф/с/0.1%BW/м²/рад²
    return B_SI * 1e-6 * 1e-6                   # → /мм²/мрад²


def total_power_kW(B0, p=P):
    return 0.633 * p.E ** 2 * B0 ** 2 * p.L * p.I


def power_density_onaxis(B0, K, p=P):
    """Плотность мощности на оси, кВт/мрад² (формула Кима 10.84·B0·E⁴·I·N·G(K) даёт Вт/мрад²)."""
    return 10.84 * B0 * p.E ** 4 * p.I * p.N * _Gpow(K) / 1e3


def critical_energy_keV(B0, p=P):
    return 0.665 * p.E ** 2 * B0           # вигглерный предел (ε_c)


def passport(gap_mm, p=P, nmax=9):
    B0 = B0_of_gap(gap_mm, p); K = K_of_B0(B0, p)
    harm = []
    for n in range(1, nmax + 1, 2):            # нечётные на оси
        harm.append(dict(n=n, E=harmonic_energy_keV(n, K, 0, p),
                         flux_dens=flux_density_onaxis(n, K, p),
                         flux=central_flux(n, K, p), bright=brightness(n, K, p)))
    return dict(gap=gap_mm, B0=B0, K=K, gamma=gamma(p),
                P_tot=total_power_kW(B0, p), P_dens=power_density_onaxis(B0, K, p),
                eps_c=critical_energy_keV(B0, p), harm=harm)


def _fmt(x):
    return f"{x:.3e}" if (x and (abs(x) >= 1e4 or abs(x) < 1e-2)) else f"{x:.3f}"


def summary(gap, p=P):
    d = passport(gap, p)
    print(f"ОНДУЛЯТОР СИЛА (in-vacuum CPMU): λ_u={p.lam_u*1e3:.0f} мм, N={p.N}, L={p.L} м; "
          f"пучок {p.E} ГэВ, I={p.I} А, εx={p.emit_x*1e12:.0f} пм, εy={p.kappa*p.emit_x*1e12:.2f} пм "
          f"(связь {p.kappa*100:.0f}%), β=({p.beta_x},{p.beta_y}) м  [ПУЧОК/ПОЛЕ — ПРЕДПОЛОЖЕНИЯ]", flush=True)
    print(f"\nЗазор {d['gap']:.1f} мм:  B0={d['B0']:.3f} Т,  K={d['K']:.3f},  γ={d['gamma']:.0f}")
    print(f"  Полная мощность СИ: {d['P_tot']:.3f} кВт;  плотность на оси: {d['P_dens']:.2f} кВт/мрад²;  "
          f"ε_c (вигглер) {d['eps_c']:.2f} кэВ")
    print(f"\n  {'гарм n':>6} | {'E_n, кэВ':>9} | {'поток на оси':>13} | {'поток (конус)':>13} | {'яркость':>11}")
    print(f"  {'':>6} | {'':>9} | {'ф/с/мрад²/0.1%':>13} | {'ф/с/0.1%BW':>13} | {'ф/с/мрад²/мм²':>11}")
    for h in d["harm"]:
        print(f"  {h['n']:>6} | {h['E']:>9.2f} | {h['flux_dens']:>13.3e} | {h['flux']:>13.3e} | {h['bright']:>11.3e}")
    # кривая перестройки 1-й гармоники по зазору
    print(f"\n  перестройка E1 по зазору:", end="")
    for g in [5, 7, 10, 13, 16, 20]:
        B = B0_of_gap(g, p); Kk = K_of_B0(B, p)
        print(f"  {g}мм→{harmonic_energy_keV(1,Kk,0,p):.1f}кэВ(K={Kk:.2f})", end="")
    print(flush=True)


def curves(p=P, ng=60):
    """Сетки по зазору для кривых перестройки/мощности."""
    gaps = np.linspace(p.gap_min, p.gap_max, ng)
    ns = [1, 3, 5, 7]
    En = {n: [] for n in ns}; Bn = {n: [] for n in ns}
    Pt = []; Pd = []
    for g in gaps:
        B0 = B0_of_gap(g, p); K = K_of_B0(B0, p)
        for n in ns:
            En[n].append(harmonic_energy_keV(n, K, 0, p)); Bn[n].append(brightness(n, K, p))
        Pt.append(total_power_kW(B0, p)); Pd.append(power_density_onaxis(B0, K, p))
    return dict(gaps=[round(float(x), 3) for x in gaps], ns=ns,
                En={n: [round(v, 3) for v in En[n]] for n in ns},
                Bn={n: [float(f"{v:.4g}") for v in Bn[n]] for n in ns},
                Pt=[round(v, 4) for v in Pt], Pd=[round(v, 2) for v in Pd])


def write_html(gap, path, p=P):
    import json
    d = passport(gap, p); cv = curves(p)
    spec = [dict(E=round(h["E"], 2), B=float(f"{h['bright']:.4g}"), F=float(f"{h['flux']:.4g}"), n=h["n"]) for h in d["harm"]]
    D = dict(gap=gap, B0=round(d["B0"], 3), K=round(d["K"], 3), Ptot=round(d["P_tot"], 3),
             Pdens=round(d["P_dens"], 1), epsc=round(d["eps_c"], 2), spec=spec, cv=cv,
             par=dict(E=p.E, I=p.I, ex=p.emit_x * 1e12, ey=p.kappa * p.emit_x * 1e12, kap=p.kappa * 100,
                      sE=p.sig_E * 100, bx=p.beta_x, by=p.beta_y, lam=p.lam_u * 1e3, N=p.N, L=p.L,
                      gmin=p.gap_min, gmax=p.gap_max))
    html = _TPL.replace("/*DATA*/", json.dumps(D, ensure_ascii=False))
    with open(path, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"HTML-паспорт: {path}", flush=True)


_TPL = r"""<!doctype html>
<html lang="ru"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Паспорт СИ ондулятора СИЛА</title>
<style>
:root{--bg:#fcfcfb;--s1:#fff;--tp:#0b0b0b;--ts:#52514e;--tm:#898781;--bd:rgba(11,11,11,.14);--ac:#185fa5;--c1:#2a78d6;--c3:#eb6834;--c5:#1b7f5a;--c7:#8a6d00}
@media (prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#1a1a19;--s1:#201f1e;--tp:#f0efec;--ts:#c3c2b7;--tm:#898781;--bd:rgba(255,255,255,.14);--ac:#378add;--c1:#5a9fe6;--c3:#f0905a;--c5:#3fae82;--c7:#c7a53a}}
:root[data-theme=dark]{--bg:#1a1a19;--s1:#201f1e;--tp:#f0efec;--ts:#c3c2b7;--tm:#898781;--bd:rgba(255,255,255,.14);--ac:#378add;--c1:#5a9fe6;--c3:#f0905a;--c5:#3fae82;--c7:#c7a53a}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--tp);font-family:system-ui,-apple-system,"Segoe UI",sans-serif;line-height:1.5;padding:22px 16px}
.wrap{max-width:900px;margin:0 auto}h1{font-size:21px;font-weight:500;margin:0 0 2px}h2{font-size:15px;font-weight:500;margin:22px 0 4px}
.sub{color:var(--ts);font-size:12px;margin:0 0 10px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(135px,1fr));gap:8px;margin:10px 0}
.kpi{background:var(--s1);border:1px solid var(--bd);border-radius:9px;padding:9px 11px}.kpi .l{font-size:11px;color:var(--ts)}.kpi .v{font-size:18px;font-weight:500;color:var(--ac)}.kpi .u{font-size:11px;color:var(--tm)}
.pic{background:var(--s1);border:1px solid var(--bd);border-radius:10px;padding:10px;margin:8px 0}
.lg{font-size:11px;color:var(--ts);margin:0 0 4px;display:flex;gap:10px;flex-wrap:wrap}.sw{width:11px;height:3px;display:inline-block;vertical-align:middle;margin-right:3px}
.cg{display:grid;grid-template-columns:1fr 1fr;gap:10px}@media(max-width:620px){.cg{grid-template-columns:1fr}}
table{border-collapse:collapse;font-size:12px;margin:4px 0}td,th{border:1px solid var(--bd);padding:3px 8px;text-align:left}
.note{font-size:11px;color:var(--ts);margin:8px 0}
</style></head><body><div class="wrap">
<h1>Паспорт синхротронного излучения ондулятора — СИЛА</h1>
<p class="sub">In-vacuum CPMU: λ_u=<b id="lam"></b> мм, N=<b id="N"></b>, L=<b id="L"></b> м. Пучок <b id="E"></b> ГэВ, I=<b id="I"></b> А, ε_x=<b id="ex"></b> пм. <b style="color:var(--c3)">Поле и параметры пучка — предположения</b> (см. память undulator-parameters).</p>
<div class="kpis" id="kpis"></div>

<h2>Спектр при зазоре <span id="g1"></span> мм (гармоники на оси)</h2>
<div class="pic"><div class="lg"><span>вертикальные линии — нечётные гармоники E_n, высота = яркость (лог)</span></div>
  <canvas id="spec" width="880" height="230"></canvas></div>

<div class="cg">
  <div class="pic"><div class="lg"><b>Перестройка E_n по зазору</b> &nbsp;<span><span class="sw" style="background:var(--c1)"></span>n=1</span><span><span class="sw" style="background:var(--c3)"></span>3</span><span><span class="sw" style="background:var(--c5)"></span>5</span><span><span class="sw" style="background:var(--c7)"></span>7</span></div>
    <canvas id="tune" width="430" height="240"></canvas></div>
  <div class="pic"><div class="lg"><b>Яркость — энергия (перестройка)</b> &nbsp;<span>лог-лог, кривые по гармоникам</span></div>
    <canvas id="bright" width="430" height="240"></canvas></div>
</div>

<h2>Мощность vs зазор</h2>
<div class="pic"><div class="lg"><span><span class="sw" style="background:var(--c1)"></span>полная P, кВт</span><span><span class="sw" style="background:var(--c3)"></span>плотность на оси, кВт/мрад²</span></div>
  <canvas id="pow" width="880" height="220"></canvas></div>

<h2>Параметры</h2>
<table id="partab"></table>
<p class="note">Формулы — стандартная теория ондуляторного излучения (Kim; X-ray Data Booklet). Яркость — свёртка дифракционного фотонного пучка с эмиттансом/бетой электронного. Параметры пучка (ток, ε_y, σ_E, β) и поле B_0(зазор) — типовые предположения 4-го поколения; подставляются реальные из решётки и магнитных измерений.</p>
<script>
const D=/*DATA*/;
const cssv=n=>getComputedStyle(document.documentElement).getPropertyValue(n).trim();
function axes(g,ml,mt,pw,ph){g.strokeStyle=cssv('--tm');g.lineWidth=1;g.beginPath();g.moveTo(ml,mt);g.lineTo(ml,mt+ph);g.lineTo(ml+pw,mt+ph);g.stroke();}
function grid(g,xs,ys,ml,mt,pw,ph){g.strokeStyle=cssv('--bd');g.lineWidth=1;g.beginPath();xs.forEach(x=>{g.moveTo(x,mt);g.lineTo(x,mt+ph);});ys.forEach(y=>{g.moveTo(ml,y);g.lineTo(ml+pw,y);});g.stroke();}
function ytit(g,t,x,yc){g.save();g.translate(x,yc);g.rotate(-Math.PI/2);g.textAlign='center';g.fillStyle=cssv('--tm');g.font='10px sans-serif';g.fillText(t,0,0);g.restore();}
const L10=Math.log10;
// универсальный линейный/лог график линий
function plot(id,series,o){const cv=document.getElementById(id),g=cv.getContext('2d'),W=cv.width,H=cv.height;g.clearRect(0,0,W,H);const tm=cssv('--tm');
  const ml=o.ml||54,mb=28,mt=10,mr=14,pw=W-ml-mr,ph=H-mt-mb;
  let xs=[],ys=[];series.forEach(s=>s.x.forEach((x,i)=>{xs.push(x);ys.push(s.y[i]);}));
  let x0=Math.min(...xs),x1=Math.max(...xs),y0,y1;
  const yl=o.ylog, xl=o.xlog;
  ys=ys.filter(v=>v>0||!yl);
  if(yl){y0=Math.min(...ys),y1=Math.max(...ys);y0=Math.pow(10,Math.floor(L10(y0)));y1=Math.pow(10,Math.ceil(L10(y1)));}
  else{y0=Math.min(0,...ys);y1=Math.max(...ys)*1.08;}
  if(xl){x0=Math.pow(10,Math.floor(L10(x0)));x1=Math.pow(10,Math.ceil(L10(x1)));}
  const X=v=>ml+(xl?(L10(v)-L10(x0))/(L10(x1)-L10(x0)):(v-x0)/(x1-x0))*pw;
  const Y=v=>mt+(1-(yl?(L10(v)-L10(y0))/(L10(y1)-L10(y0)):(v-y0)/(y1-y0)))*ph;
  // сетка + подписи
  const gx=[],gy=[];
  if(xl){for(let e=Math.floor(L10(x0));e<=L10(x1);e++){const v=Math.pow(10,e);gx.push(X(v));}}
  else{for(let i=0;i<=4;i++)gx.push(ml+i/4*pw);}
  if(yl){for(let e=Math.floor(L10(y0));e<=L10(y1);e++)gy.push(Y(Math.pow(10,e)));}
  else{for(let i=0;i<=4;i++)gy.push(mt+i/4*ph);}
  grid(g,gx,gy,ml,mt,pw,ph);axes(g,ml,mt,pw,ph);
  g.fillStyle=tm;g.font='9px sans-serif';g.textAlign='center';
  if(xl){for(let e=Math.floor(L10(x0));e<=L10(x1);e++){const v=Math.pow(10,e);g.fillText(v>=1?v:v.toFixed(2),X(v),H-15);}}
  else{for(let i=0;i<=4;i++){const v=x0+i/4*(x1-x0);g.fillText(v.toFixed(v<10?1:0),ml+i/4*pw,H-15);}}
  g.textAlign='right';
  if(yl){for(let e=Math.floor(L10(y0));e<=L10(y1);e++)g.fillText('1e'+e,ml-3,Y(Math.pow(10,e))+3);}
  else{for(let i=0;i<=4;i++){const v=y0+i/4*(y1-y0);g.fillText(v.toFixed(v<10?1:0),ml-3,mt+(1-i/4)*ph+3);}}
  g.fillStyle=tm;g.textAlign='center';g.fillText(o.xt,ml+pw/2,H-2);ytit(g,o.yt,11,mt+ph/2);
  series.forEach(s=>{g.strokeStyle=s.c;g.lineWidth=1.5;g.beginPath();let st=false;s.x.forEach((x,i)=>{const v=s.y[i];if(yl&&v<=0)return;const px=X(x),py=Y(v);st?g.lineTo(px,py):g.moveTo(px,py);st=true;});g.stroke();
    if(s.dot)s.x.forEach((x,i)=>{const v=s.y[i];if(yl&&v<=0)return;g.fillStyle=s.c;g.beginPath();g.arc(X(x),Y(v),2,0,6.29);g.fill();});});}
// стем-спектр
function stem(id,pts,o){const cv=document.getElementById(id),g=cv.getContext('2d'),W=cv.width,H=cv.height;g.clearRect(0,0,W,H);const tm=cssv('--tm');
  const ml=58,mb=28,mt=12,mr=14,pw=W-ml-mr,ph=H-mt-mb;
  const x1=Math.max(...pts.map(p=>p.E))*1.1;const ys=pts.map(p=>p.B);
  let y0=Math.pow(10,Math.floor(L10(Math.min(...ys)))),y1=Math.pow(10,Math.ceil(L10(Math.max(...ys))));
  const X=v=>ml+v/x1*pw,Y=v=>mt+(1-(L10(v)-L10(y0))/(L10(y1)-L10(y0)))*ph;
  const gy=[];for(let e=Math.floor(L10(y0));e<=L10(y1);e++)gy.push(Y(Math.pow(10,e)));
  const gx=[];for(let i=0;i<=5;i++)gx.push(ml+i/5*pw);
  grid(g,gx,gy,ml,mt,pw,ph);axes(g,ml,mt,pw,ph);
  pts.forEach(p=>{g.strokeStyle=cssv('--c1');g.lineWidth=2;g.beginPath();g.moveTo(X(p.E),mt+ph);g.lineTo(X(p.E),Y(p.B));g.stroke();
    g.fillStyle=cssv('--c1');g.beginPath();g.arc(X(p.E),Y(p.B),3,0,6.29);g.fill();
    g.fillStyle=tm;g.font='9px sans-serif';g.textAlign='center';g.fillText('n='+p.n,X(p.E),Y(p.B)-6);});
  g.fillStyle=tm;g.font='9px sans-serif';g.textAlign='center';for(let i=0;i<=5;i++)g.fillText((i/5*x1).toFixed(0),ml+i/5*pw,H-15);g.fillText('энергия фотона, кэВ',ml+pw/2,H-2);
  g.textAlign='right';for(let e=Math.floor(L10(y0));e<=L10(y1);e++)g.fillText('1e'+e,ml-3,Y(Math.pow(10,e))+3);
  ytit(g,'яркость, ф/с/мрад²/мм²',11,mt+ph/2);}
function run(){
  const p=D.par;
  document.getElementById('lam').textContent=p.lam;document.getElementById('N').textContent=p.N;document.getElementById('L').textContent=p.L;
  document.getElementById('E').textContent=p.E;document.getElementById('I').textContent=p.I;document.getElementById('ex').textContent=p.ex;
  document.getElementById('g1').textContent=D.gap;
  document.getElementById('kpis').innerHTML=[
    ['B₀',D.B0,'Т'],['K',D.K,''],['E₁ (1-я гарм.)',D.spec[0].E,'кэВ'],
    ['Полная мощность',D.Ptot,'кВт'],['Плотн. на оси',D.Pdens,'кВт/мрад²'],['ε_c (вигглер)',D.epsc,'кэВ']
  ].map(k=>`<div class="kpi"><div class="l">${k[0]}</div><div class="v">${k[1]}</div><div class="u">${k[2]}</div></div>`).join('');
  stem('spec',D.spec,{});
  const col={1:cssv('--c1'),3:cssv('--c3'),5:cssv('--c5'),7:cssv('--c7')};
  plot('tune',D.cv.ns.map(n=>({x:D.cv.gaps,y:D.cv.En[n],c:col[n]})),{xt:'зазор, мм',yt:'E_n, кэВ'});
  plot('bright',D.cv.ns.map(n=>({x:D.cv.En[n],y:D.cv.Bn[n],c:col[n],dot:false})),{xlog:true,ylog:true,xt:'энергия фотона, кэВ',yt:'яркость'});
  plot('pow',[{x:D.cv.gaps,y:D.cv.Pt,c:cssv('--c1')},{x:D.cv.gaps,y:D.cv.Pd,c:cssv('--c3')}],{ylog:true,xt:'зазор, мм',yt:'мощность (лог)'});
  const rows=[['Период λ_u',p.lam+' мм'],['Число периодов N',p.N],['Длина L',p.L+' м'],['Зазор',p.gmin+'–'+p.gmax+' мм'],
    ['Энергия пучка',p.E+' ГэВ'],['Ток I',p.I+' А (предпол.)'],['ε_x',p.ex+' пм'],['ε_y',p.ey.toFixed(2)+' пм (связь '+p.kap+'%, предпол.)'],
    ['σ_E/E',p.sE+'% (предпол.)'],['β_x, β_y',p.bx+', '+p.by+' м (предпол.)']];
  document.getElementById('partab').innerHTML=rows.map(r=>`<tr><td>${r[0]}</td><td>${r[1]}</td></tr>`).join('');
}
run();
let _t;new MutationObserver(()=>{clearTimeout(_t);_t=setTimeout(run,60);}).observe(document.documentElement,{attributes:true,attributeFilter:['data-theme']});
</script>
</div></body></html>"""


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--gap", type=float, default=5.0, help="зазор, мм")
    ap.add_argument("--b0", type=float, default=None, help="переопределить B0 при якорном зазоре, Т")
    ap.add_argument("--current", type=float, default=None, help="ток, А")
    ap.add_argument("--emit-y", type=float, default=None, help="εy, пм (переопределяет связь)")
    ap.add_argument("--beta-x", type=float, default=None, help="β_x, м (напр. реальная из twiss)")
    ap.add_argument("--beta-y", type=float, default=None, help="β_y, м")
    ap.add_argument("--from-lattice", nargs=2, metavar=("LAT", "GROUP"), default=None,
                    help="взять реальные β из twiss решётки: напр. --from-lattice 27 und_19")
    ap.add_argument("--html", default="workspace/results/undulator_passport.html", help="путь HTML-паспорта (пусто — не писать)")
    a = ap.parse_args()
    if a.b0 is not None: P.B0_ref = a.b0
    if a.current is not None: P.I = a.current
    if a.emit_y is not None: P.kappa = (a.emit_y * 1e-12) / P.emit_x
    if a.from_lattice is not None:
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "undulators"))
        from lattice_optics import undulator_optics
        o = undulator_optics(int(a.from_lattice[0]), a.from_lattice[1])
        P.beta_x, P.beta_y = o["betax"], o["betay"]
        print(f"[реальные β из реш.{a.from_lattice[0]} {a.from_lattice[1]}]: "
              f"βx={P.beta_x:.3f}, βy={P.beta_y:.3f} м", flush=True)
    if a.beta_x is not None: P.beta_x = a.beta_x
    if a.beta_y is not None: P.beta_y = a.beta_y
    summary(a.gap)
    if a.html:
        write_html(a.gap, a.html)
