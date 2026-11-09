# -*- coding: utf-8 -*-
"""Сводный HTML-ПАСПОРТ ЛИНИИ для выбранного экрана: тепло + спектр + мультиисточник в одном файле.

Объединяет этапы A/B/#1/#2/#3 для одного экрана:
  ОНДУЛЯТОР: тепловой след (Вт/мм²), спектр на оси (гармоники, качество ИД, реальные β), мультиисточник
             (ИД + фон соседних диполей);
  BEND (волноотвод): веер поворотного СИ + сплошной спектр.
Картинки встраиваются base64 (самодостаточный файл). Стиль — как у undulator_sr (темы, KPI-карточки).

ЗАПУСК: py -3.14 undulators/screen_report.py --lattice 27 --screen SCR_U01 --gap 5 --phase-error-deg 4
"""
import os, sys, io, csv, math, base64, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

import undulator_sr as sr
import bend_sr as bsr
from screen_analysis import load_screens, load_groups, emitter, power_params, footprint as power_footprint, _unit, RESDIR
from screen_spectrum import trajectory, calibrate, center_spectrum, find_harmonics, beam_params, energy_smear
from undulator_field import apply_phase_errors
from bend_screen import load_lattice, dipole_arc, illuminates
from bend_screen import footprint as bend_footprint, emitter_geom
from build_screens import bend_source
from screen_multi import find_undulators, find_dipoles


_CFG = dict(scrollZoom=True, displaylogo=False, responsive=True,
            modeBarButtonsToRemove=["lasso2d", "select2d"],
            toImageButtonOptions=dict(format="png", scale=2))
_LAYOUT = dict(template="plotly_white", paper_bgcolor="rgba(0,0,0,0)",
               margin=dict(l=62, r=12, t=12, b=46), font=dict(size=11))


def fig_div(fig):
    """Интерактивный plotly-div (зум/пан/ховер); plotly.js подключается один раз в шапке."""
    return fig.to_html(full_html=False, include_plotlyjs=False, config=_CFG)


def heatmap(S, extent, zlabel, equal=True):
    """2D-карта: data-backed heatmap, зум перерисовывает чётко (а не пикселизует картинку)."""
    x0, x1, y0, y1 = extent
    nu, nv = S.shape
    x = np.linspace(x0, x1, nu); y = np.linspace(y0, y1, nv)
    fig = go.Figure(go.Heatmap(z=S.T, x=x, y=y, colorscale="Inferno",
                               colorbar=dict(title=dict(text=zlabel, side="right")),
                               hovertemplate="u=%{x:.3f} мм<br>v=%{y:.3f} мм<br>%{z:.4g}<extra></extra>"))
    fig.update_xaxes(title_text="u (горизонт.), мм", constrain="domain")
    fig.update_yaxes(title_text="v (верт.), мм", **({"scaleanchor": "x", "scaleratio": 1} if equal else {}))
    fig.update_layout(height=430, **_LAYOUT)
    return fig


def spectrum(traces, ylog=True, xlog=False, height=360, ylabel="ф/с/мрад²/0.1%BW"):
    fig = go.Figure()
    for x, y, name, color, width in traces:
        fig.add_trace(go.Scatter(x=x, y=np.maximum(y, 1), mode="lines", name=name,
                                 line=dict(color=color, width=width),
                                 hovertemplate="%{x:.2f} кэВ<br>%{y:.3g}<extra>" + name + "</extra>"))
    fig.update_xaxes(title_text="энергия фотона, кэВ", type="log" if xlog else "linear")
    fig.update_yaxes(title_text=ylabel, type="log" if ylog else "linear")
    fig.update_layout(height=height, legend=dict(orientation="h", y=1.02, yanchor="bottom",
                                                 font=dict(size=10)), **_LAYOUT)
    return fig


def kpi(label, val, unit=""):
    return f'<div class="kpi"><div class="l">{label}</div><div class="v">{val}</div><div class="u">{unit}</div></div>'


def section(title, div, caption, note=""):
    n = f'<p class="note">{note}</p>' if note else ""
    return (f'<h2>{title}</h2><div class="pic">{div}'
            f'<div class="cap">{caption}</div></div>{n}')


# Пояснение единицы оси Y — вставляется ОДИН раз, перед первым спектром (тепловые карты до него — в Вт/мм²)
UNIT_NOTE = (
    '<div class="unitbox"><b>Единица оси Y — ф/с/мрад²/0.1%BW</b> '
    '(спектральная угловая плотность потока): число фотонов в секунду, летящих в телесный угол 1 мрад², '
    'в полосе энергий 0.1% (ΔE/E) вокруг данной энергии. Это <b>плотность</b> — не зависит от апертуры и '
    'монохроматора; для реальных отсчётов умножьте на свой телесный угол и полосу. Дополнительное деление '
    'на площадь источника (/мм²) даёт <b>яркость</b> (см. KPI). Ондулятор и диполи — в одной единице, '
    'поэтому их законно складывать на одной оси.</div>')


def undulator_report(scr, lattice, gap, args):
    ci, byid = load_lattice(lattice, None)
    groups = load_groups(lattice)
    B0 = sr.B0_of_gap(gap); K = sr.K_of_B0(B0); g = sr.gamma()
    _, _, P_tot, Pd0, sig = power_params(gap)

    # реальные β
    from lattice_optics import undulator_optics
    o = undulator_optics(lattice, scr["src"])
    sr.P.beta_x, sr.P.beta_y = o["betax"], o["betay"]
    beam = beam_params(sr.P, o["betax"], o["betay"], args.kappa, args.sigE, o["etaxp"])

    # ── тепловой след (A) ──
    emits = {}
    for gid, insts in groups.items():
        P, axis = emitter(insts, ci, byid)
        if P is not None:
            emits[gid] = (P, axis, Pd0, sig)
    found = []
    for gid, (P, axis, pd0, sg) in emits.items():
        vec = scr["C"] - P; D = float(np.linalg.norm(vec)); dh = vec / D
        if float(np.dot(dh, axis)) <= 0:
            continue
        th = math.acos(max(-1.0, min(1.0, float(np.dot(dh, axis)))))
        half = 0.5 * math.hypot(scr["w"], scr["h"]) / 1e3 / D
        if th <= 6 * sg + half:
            found.append(gid)
    emitters = [emits[gid] for gid in found]
    fp = power_footprint(scr, emitters, 301)                    # мельче сетка → зум в конус остаётся чётким
    heat_div = fig_div(heatmap(fp["S"], fp["extent"], "Вт/мм²"))

    # ── спектр на оси (B) с качеством ИД (#1) ──
    tr = trajectory(K, sr.P.lam_u, sr.P.N, g, args.ppp); scale, E1 = calibrate(tr)
    if args.phase_error_deg:
        tr = apply_phase_errors(tr, args.phase_error_deg)
    P0, axis0 = emits[scr["src"]][0], emits[scr["src"]][1]
    uhat = _unit(np.cross([0, 0, 1.0], axis0)); vhat = _unit(np.cross(axis0, uhat))
    dc = _unit(scr["C"] - P0); thx0 = float(np.dot(dc, uhat)); thy0 = float(np.dot(dc, vhat))
    Es = np.linspace(0.5 * E1, sr.harmonic_energy_keV(9, K, 0.0) * 1.05, 3200)
    Ise, Ibeam = center_spectrum(tr, scale, thx0, thy0, Es, beam)
    harm = find_harmonics(tr, scale, K)
    spec_div = fig_div(spectrum([(Es, Ise, "один электрон", "#9fb4c8", 0.8),
                                 (Es, Ibeam, "пучок", "#185fa5", 1.6)]))

    # ── мультиисточник (#2) ──
    pole_insts = set().union(*groups.values()) if groups else set()
    unds = find_undulators(scr, groups, ci, byid, sig)
    dips = find_dipoles(scr, ci, byid, g, args.E, pole_insts)
    total = np.zeros_like(Es); parts = []
    for gid, Pu, au in unds:
        uh = _unit(np.cross([0, 0, 1.0], au)); vh = _unit(np.cross(au, uh))
        dcu = _unit(scr["C"] - Pu)
        _, Ib = center_spectrum(tr, scale, float(np.dot(dcu, uh)), float(np.dot(dcu, vh)), Es, beam)
        total = total + Ib; parts.append((f"ИД {gid}", Ib, "#185fa5"))
    for iid, d in dips:
        dcu = _unit(scr["C"] - d["C"]); psi = math.asin(max(-1, min(1, float(np.dot(dcu, d["vhat"])))))
        Fa = bsr.flux_angular(Es, psi, d["ec"], args.E, args.I, g); total = total + Fa
        parts.append((f"диполь {d['dev']}#{iid}", Fa, "#1b7f5a"))
    mtr = [(Es, S, lbl, c, 1.0) for lbl, S, c in parts] + [(Es, total, "сумма", "#111111", 1.8)]
    multi_div = fig_div(spectrum(mtr))

    # KPI + таблицы
    br1 = sr.brightness(1, K)
    kpis = "".join([
        kpi("K", f"{K:.3f}"), kpi("E₁", f"{harm[0][1]:.2f}", "кэВ"),
        kpi("Пик тепла", f"{fp['peak']:.1f}", "Вт/мм²"), kpi("Перехвачено", f"{fp['total']:.0f}", "Вт"),
        kpi("Яркость n=1", f"{br1:.2e}", "ф/с/мрад²/мм²/0.1%"), kpi("Осветителей", f"{len(unds)}+{len(dips)}", "ИД+бенд")])
    htab = "".join(f"<tr><td>n={n}</td><td>{En:.2f} кэВ</td><td>{Fn:.2e}</td></tr>" for n, En, Fn in harm)
    ptab = (f"<tr><td>Источник</td><td>{scr['src']} (ондулятор)</td></tr>"
            f"<tr><td>Расстояние D</td><td>{np.linalg.norm(scr['C']-P0):.2f} м</td></tr>"
            f"<tr><td>Экран</td><td>{scr['w']:.0f}×{scr['h']:.0f} мм</td></tr>"
            f"<tr><td>β в месте ИД</td><td>({o['betax']:.2f}, {o['betay']:.2f}) м [twiss реш.{lattice}]</td></tr>"
            f"<tr><td>Зазор</td><td>{gap:g} мм (B₀={B0:.3f} Т)</td></tr>"
            f"<tr><td>Качество ИД σ_φ</td><td>{args.phase_error_deg:g}°</td></tr>")

    body = (f'<div class="kpis">{kpis}</div>'
            + section("Тепловая нагрузка (веер мощности)", heat_div,
                      f"Плотность мощности на экране; пик {fp['peak']:.1f} Вт/мм², перехвачено {fp['total']:.0f} Вт. "
                      f"Конус ~0.9 мм на 50-мм экране — зумом (колесо/рамка) рассмотрите центр.",
                      "Некогерентная проекция угловой мощности устройства на плоскость экрана.")
            + UNIT_NOTE
            + section("Спектр на оси", spec_div,
                      "Гармоники (острые — один электрон; размытые — пучок с реальными β, σ_E и качеством ИД).",
                      "Гармоники из когерентного интеграла поля; фазовые ошибки бьют по высоким.")
            + section("Мультиисточник (ИД + фон диполей)", multi_div,
                      "Чёрное — сумма интенсивностей: гармоники ИД на гладком фоне поворотного СИ соседних диполей.",
                      "Внутри ИД — когерентно; между источниками — некогерентно (мощности).")
            + '<div class="tcols"><div><h3>Гармоники на оси</h3><table><tr><th>n</th><th>E</th>'
              f'<th>поток, ф/с/мрад²/0.1%</th></tr>{htab}</table></div>'
              f'<div><h3>Параметры</h3><table>{ptab}</table></div></div>')
    return body, f"{scr['sid']} — ондуляторная линия"


def bend_report(scr, lattice, gap, args):
    ci, byid = load_lattice(lattice, None)
    g = bsr.gamma(args.E)
    bound = bend_source(scr["src"], None, ci, byid, args.E)
    emitters = []
    for iid, r in byid.items():
        d = dipole_arc(r, ci, args.E)
        if d is None or not illuminates(d, scr["C"], g):
            continue
        emitters.append(emitter_geom(d, scr["C"], g))
    fp = bend_footprint(scr, emitters, g, args.I, 241)
    foot_div = fig_div(heatmap(fp["S"], fp["extent"], "Вт/мм²", equal=False))

    dc = _unit(scr["C"] - bound["P"]); psi = math.asin(max(-1, min(1, float(np.dot(dc, bound["vhat"])))))
    Es = np.geomspace(0.03 * bound["ec"], 12 * bound["ec"], 400)
    Fa = bsr.flux_angular(Es, psi, bound["ec"], args.E, args.I, g)
    sfig = spectrum([(Es, Fa, "поворотное СИ", "#1b7f5a", 1.6)], xlog=True)
    sfig.add_vline(x=bound["ec"], line=dict(color="#c0392b", dash="dash", width=1),
                   annotation_text=f"ε_c={bound['ec']:.2f} кэВ", annotation_font_size=10)
    spec_div = fig_div(sfig)

    kpis = "".join([
        kpi("ε_c", f"{bound['ec']:.2f}", "кэВ"), kpi("B", f"{bound['B']:.3f}", "Т"),
        kpi("θ изгиба", f"{bound['theta0']*1e3:+.2f}", "мрад"), kpi("Угол отбора", f"{bound['phi']*1e3:.2f}", "мрад"),
        kpi("Пик тепла", f"{fp['peak']:.2f}", "Вт/мм²"), kpi("Перехвачено", f"{fp['total']:.0f}", "Вт")])
    ptab = (f"<tr><td>Источник</td><td>диполь {bound['dev']} (волноотвод)</td></tr>"
            f"<tr><td>Расстояние D</td><td>{np.linalg.norm(scr['C']-bound['P']):.2f} м</td></tr>"
            f"<tr><td>Экран</td><td>{scr['w']:.0f}×{scr['h']:.0f} мм</td></tr>"
            f"<tr><td>Осветителей</td><td>{len(emitters)} диполь(ей)</td></tr>")
    body = (f'<div class="kpis">{kpis}</div>'
            + section("Веер поворотного СИ", foot_div,
                      f"Горизонтальная полоса (развёртка дуги), вертикаль ~1/γ; перехвачено {fp['total']:.0f} Вт. "
                      f"Зум (колесо/рамка) — рассмотреть вертикальный профиль.")
            + UNIT_NOTE
            + section("Спектр поворотного СИ в центре экрана", spec_div,
                      f"Сплошной спектр, пик у ε_c={bound['ec']:.2f} кэВ (пунктир), завал выше.")
            + f'<h3>Параметры</h3><table>{ptab}</table>')
    return body, f"{scr['sid']} — поворотная линия (волноотвод)"


def main():
    ap = argparse.ArgumentParser(description="сводный HTML-паспорт линии для экрана")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--screen", required=True)
    ap.add_argument("--gap", type=float, default=5.0)
    ap.add_argument("--E", type=float, default=6.0); ap.add_argument("--I", type=float, default=0.2)
    ap.add_argument("--kappa", type=float, default=sr.P.kappa); ap.add_argument("--sigE", type=float, default=sr.P.sig_E)
    ap.add_argument("--phase-error-deg", type=float, default=0.0)
    ap.add_argument("--ppp", type=int, default=64)
    ap.add_argument("--out", default=None)
    a = ap.parse_args()

    scr = next((s for s in load_screens(a.lattice) if s["sid"] == a.screen), None)
    if scr is None:
        raise SystemExit(f"экран {a.screen} не найден в screens.csv реш.{a.lattice}")
    kind = scr.get("source_kind") or ("bend" if scr["mode"] == "bend" else "undulator")
    print(f"паспорт линии: экран {a.screen}, источник-{kind}, зазор {a.gap:g} мм …")
    body, title = (bend_report if kind == "bend" else undulator_report)(scr, a.lattice, a.gap, a)

    html = (_TPL.replace("/*TITLE*/", title).replace("/*BODY*/", body)
            .replace("/*SUB*/", f"Решётка {a.lattice} · ток {a.I*1e3:.0f} мА · E={a.E} ГэВ · "
                                f"графики интерактивны (зум колесом/рамкой, пан, ховер)")
            .replace("/*PLOTLYJS*/", get_plotlyjs()))
    out = a.out or os.path.join(RESDIR, f"report_{a.screen}.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").write(html)
    print(f"[OK] {out}")


_TPL = r"""<!doctype html><html lang="ru"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1"><title>/*TITLE*/</title>
<script>/*PLOTLYJS*/</script><style>
:root{--bg:#fcfcfb;--s1:#fff;--tp:#0b0b0b;--ts:#52514e;--tm:#898781;--bd:rgba(11,11,11,.14);--ac:#185fa5}
@media(prefers-color-scheme:dark){:root:not([data-theme=light]){--bg:#1a1a19;--s1:#201f1e;--tp:#f0efec;--ts:#c3c2b7;--tm:#898781;--bd:rgba(255,255,255,.14);--ac:#5a9fe6}}
:root[data-theme=dark]{--bg:#1a1a19;--s1:#201f1e;--tp:#f0efec;--ts:#c3c2b7;--tm:#898781;--bd:rgba(255,255,255,.14);--ac:#5a9fe6}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--tp);font-family:system-ui,"Segoe UI",sans-serif;line-height:1.5;padding:22px 16px}
.wrap{max-width:900px;margin:0 auto}h1{font-size:21px;font-weight:600;margin:0 0 2px}h2{font-size:15px;font-weight:600;margin:22px 0 6px}
h3{font-size:13px;font-weight:600;margin:10px 0 4px;color:var(--ts)}.sub{color:var(--ts);font-size:12px;margin:0 0 12px}
.kpis{display:grid;grid-template-columns:repeat(auto-fit,minmax(130px,1fr));gap:8px;margin:12px 0}
.kpi{background:var(--s1);border:1px solid var(--bd);border-radius:9px;padding:9px 11px}
.kpi .l{font-size:11px;color:var(--ts)}.kpi .v{font-size:17px;font-weight:600;color:var(--ac)}.kpi .u{font-size:10px;color:var(--tm)}
.pic{background:var(--s1);border:1px solid var(--bd);border-radius:10px;padding:10px;margin:6px 0}
.pic .plotly-graph-div{width:100%}.cap{font-size:11px;color:var(--ts);margin-top:6px}
.note{font-size:11px;color:var(--tm);margin:2px 0 4px;font-style:italic}
.unitbox{background:color-mix(in srgb,var(--ac) 9%,var(--s1));border:1px solid var(--bd);border-left:3px solid var(--ac);border-radius:8px;padding:9px 12px;margin:10px 0;font-size:12px;color:var(--ts);line-height:1.5}
.unitbox b{color:var(--tp)}
.tcols{display:grid;grid-template-columns:1fr 1fr;gap:14px}@media(max-width:620px){.tcols{grid-template-columns:1fr}}
table{border-collapse:collapse;font-size:12px;width:100%}td,th{border:1px solid var(--bd);padding:3px 8px;text-align:left}
th{color:var(--ts);font-weight:600}
</style></head><body><div class="wrap">
<h1>/*TITLE*/</h1><p class="sub">/*SUB*/</p>/*BODY*/
<p class="note" style="margin-top:20px">SynchrotronSim · единая модель источника и экранов. Пучковые допущения (ε_x, связь, σ_E) помечены; β — реальные из twiss. Экран заморожен в глобале, источник — по фактической геометрии.</p>
</div></body></html>"""


if __name__ == "__main__":
    main()
