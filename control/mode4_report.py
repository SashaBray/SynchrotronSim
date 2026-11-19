# -*- coding: utf-8 -*-
"""Режим 4 — HTML-отчёт эксперимента ОС ПО ТРЕБОВАНИЮ (отдельно от прогона).

Читает результаты из БД по ключу (fb_results.csv + fb_artifacts.csv/.npz) и собирает самодостаточный
интерактивный HTML: KPI, временной ряд орбиты/света по окнам, интерференционная карта (идеал vs
измеренная) и фазовый портрет — если эти сигналы были в подписке. Это ПЕРВЫЙ минимальный генератор;
он изолирован от физики (Stage 1–4) и будет улучшаться отдельно. Переиспользует plotly-тулкит
screen_report (_TPL/fig_div/heatmap/get_plotlyjs).

ЗАПУСК: py -3.14 control/mode4_report.py --exp 1 [--no-open]
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, _ROOT)
sys.path.insert(0, os.path.join(_ROOT, "undulators"))

import numpy as np
import plotly.graph_objects as go

import screen_report as rep
from control import fb_config, fb_store


def _kpi(label, val, unit=""):
    return f'<div class="kpi"><div class="l">{label}</div><div class="v">{val}</div><div class="u">{unit}</div></div>'


def _last(res, sig):
    if sig in res and res[sig][1]:
        order = np.argsort(res[sig][0])
        return float(np.asarray(res[sig][1])[order][-1])
    return None


def main():
    ap = argparse.ArgumentParser(description="режим 4 — HTML-отчёт эксперимента ОС")
    ap.add_argument("--exp", type=int, required=True)
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()

    cfg = fb_config.load(a.exp)
    res = fb_store.read_results(a.exp)
    art = fb_store.read_artifacts(a.exp)
    if not res and not art:
        raise SystemExit(f"mode4_report: нет результатов для fb_exp_id={a.exp} — сперва запустите mode4_run")

    body = []

    # --- KPI ---
    kpis = []
    for lab, sig, unit in (("орбита X", "orbit_rms.x", "мкм rms"), ("орбита Y", "orbit_rms.y", "мкм rms"),
                           ("орбита X, среднее по времени", "orbit_dc.x", "мкм rms"),
                           ("орбита Y, среднее по времени", "orbit_dc.y", "мкм rms"),
                           ("свет X", "light.x", "мкм rms"), ("свет Y", "light.y", "мкм rms"),
                           ("резкость конуса", "screen.sharpness_pct", "% идеала"),
                           ("джиттер угла X", "screen.jitter_x_urad", "мкрад"),
                           ("удлинение кольца", "rf.circ_um", "мкм"), ("энергия δ", "rf.delta_true_ppm", "ppm"),
                           ("частота ВЧ Δf", "rf.df_hz", "Гц")):
        v = _last(res, sig)
        if v is not None:
            kpis.append(_kpi(lab, f"{v:.2f}", unit))
    if kpis:
        body.append(f'<div class="kpis">{"".join(kpis)}</div>')
    body.append(f'<div class="unitbox">Эксперимент ОС fb_exp_id={cfg.fb_exp_id}, решётка {cfg.lattice_id}. '
                f'{cfg.note or ""} Возмущения: {", ".join(c["name"] for c in cfg.disturbance_spec) or "—"}. '
                f'Контуры: {", ".join(l.name for l in cfg.enabled_loops()) or "—"}. '
                f'Конфиг и результаты — в БД (таблицы fb_*). Генератор отчёта минимальный (v1).</div>')

    # --- временной ряд орбиты/света по окнам ---
    series = [("orbit_rms.x", "орбита X", "#2b7bba"), ("orbit_rms.y", "орбита Y", "#e15759"),
              ("light.x", "свет X", "#59a14f"), ("light.y", "свет Y", "#b07aa1")]
    present = [(s, n, c) for s, n, c in series if s in res and res[s][1]]
    if present:
        fo = go.Figure()
        for sig, name, col in present:
            order = np.argsort(res[sig][0])
            x = np.asarray(res[sig][0])[order]; y = np.asarray(res[sig][1])[order]
            fo.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=name, line=dict(width=1.8, color=col)))
        fo.update_layout(height=340, margin=dict(l=60, r=16, t=10, b=42), plot_bgcolor="white",
                         legend=dict(orientation="h", y=1.02, yanchor="bottom"))
        fo.update_xaxes(title_text="окно усреднения"); fo.update_yaxes(title_text="мкм rms", rangemode="tozero")
        body.append(f'<h2>Сенсоры во времени (контур ОС по окнам)</h2><div class="pic">{rep.fig_div(fo)}</div>')

    # --- резонатор и энергия (рекордер rf): ΔC, δ (истинная и оценка по BPM), Δf ---
    rfs = [("rf.delta_true_ppm", "δ истинная, ppm", "#2b7bba", "y"), ("rf.delta_ppm", "δ̂ по BPM, ppm", "#76b7b2", "y"),
           ("rf.circ_um", "ΔC, мкм", "#9c755f", "y2"), ("rf.df_hz", "Δf ВЧ, Гц", "#e15759", "y2")]
    prf = [(s, n, c, ax) for s, n, c, ax in rfs if s in res and res[s][1]]
    if prf:
        frf = go.Figure()
        for sig, name, col, ax in prf:
            order = np.argsort(res[sig][0])
            x = np.asarray(res[sig][0])[order]; y = np.asarray(res[sig][1])[order]
            frf.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=name, yaxis=ax,
                                     line=dict(width=1.8, color=col, dash=("dot" if sig == "rf.delta_ppm" else None))))
        frf.update_layout(height=360, margin=dict(l=60, r=60, t=10, b=42), plot_bgcolor="white",
                          legend=dict(orientation="h", y=1.02, yanchor="bottom"),
                          yaxis=dict(title="энергия δ, ppm"),
                          yaxis2=dict(title="ΔC, мкм  ·  Δf, Гц", overlaying="y", side="right"))
        frf.update_xaxes(title_text="окно усреднения")
        body.append('<h2>Резонатор и энергия пучка</h2><div class="unitbox">Длина кольца ΔC (возмущение '
                    'circumference) сдвигает равновесную энергию δ ≈ −ΔC/(αc·C); контур ВЧ-частоты (если включён) '
                    'держит δ, двигая f_RF на Δf ≈ −f·ΔC/C. δ̂ — оценка по дисперсионной компоненте орбиты BPM '
                    '(датчик контура), δ — средняя энергия частиц пучка относительно проектного равновесия.</div>'
                    f'<div class="pic">{rep.fig_div(frf)}</div>')

    # --- оптика во времени (tune/chroma/adts) — для plant=optics ---
    optics = [("tune.x", "νx", "#2b7bba"), ("tune.y", "νy", "#e15759"),
              ("chroma.x", "ξx", "#59a14f"), ("chroma.y", "ξy", "#b07aa1"),
              ("adts.xx", "∂νx/∂Jx", "#edc948"), ("adts.yy", "∂νy/∂Jy", "#ff9da7"), ("adts.xy", "∂νx/∂Jy", "#9c755f")]
    popt = [(s, n, c) for s, n, c in optics if s in res and res[s][1]]
    if popt:
        fop = go.Figure()
        for sig, name, col in popt:
            order = np.argsort(res[sig][0])
            x = np.asarray(res[sig][0])[order]; y = np.asarray(res[sig][1])[order]
            fop.add_trace(go.Scatter(x=x, y=y, mode="lines+markers", name=name, line=dict(width=1.8, color=col)))
        fop.update_layout(height=340, margin=dict(l=60, r=16, t=10, b=42), plot_bgcolor="white",
                          legend=dict(orientation="h", y=1.02, yanchor="bottom"))
        fop.update_xaxes(title_text="окно"); fop.update_yaxes(title_text="частоты / хроматизм / ADTS")
        body.append(f'<h2>Оптика во времени (контуры квад/секступоль/октуполь)</h2><div class="pic">{rep.fig_div(fop)}</div>')

    # --- интерференция: идеал vs измеренная ---
    if "screen" in art:
        s = art["screen"]; ext = s["extent"].tolist() if "extent" in s else None
        zl = "ф/с/мрад²/0.1%BW"
        e0 = float(s["energy_keV"][0]) if "energy_keV" in s else 0.0
        jx, jy = (s["jitter_urad"].tolist() if "jitter_urad" in s else [0.0, 0.0])
        for key, title in (("ideal", "Идеал (стабильный пучок)"), ("measured", "С измеренным джиттером пучка")):
            if key in s:
                body.append(f'<h2>Интерференция — {title}</h2><div class="pic">'
                            f'{rep.fig_div(rep.heatmap(np.asarray(s[key]), ext, zl))}'
                            f'<div class="cap">1-я гармоника, {e0:.2f} кэВ; измеренный джиттер угла '
                            f'σx′={jx:.2f}, σy′={jy:.2f} мкрад</div></div>')

    # --- фазовый портрет ---
    if "phase" in art and "centroid" in art["phase"]:
        cen = np.asarray(art["phase"]["centroid"], float)
        if cen.size:
            fp = go.Figure(go.Scattergl(x=cen[:, 0] * 1e6, y=cen[:, 1] * 1e6, mode="markers",
                                        marker=dict(size=3, color=np.arange(len(cen)), colorscale="Viridis")))
            fp.update_layout(height=340, margin=dict(l=60, r=16, t=10, b=42), plot_bgcolor="white")
            fp.update_xaxes(title_text="центроид x, мкм"); fp.update_yaxes(title_text="центроид x′, мкрад")
            body.append(f'<h2>Фазовый портрет центроида</h2><div class="pic">{rep.fig_div(fp)}</div>')

    html = (rep._TPL.replace("/*TITLE*/", f"Эксперимент ОС fb_exp_id={cfg.fb_exp_id}")
            .replace("/*SUB*/", f"решётка {cfg.lattice_id} · СИ {'вкл' if cfg.sr else 'выкл'} · "
                                f"результаты из БД · графики интерактивны")
            .replace("/*BODY*/", "".join(body))
            .replace("/*PLOTLYJS*/", rep.get_plotlyjs()))
    out = os.path.join(_ROOT, "workspace", "results", "fb", str(a.exp), "report.html")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").write(html)
    rel = os.path.relpath(out, _ROOT).replace(os.sep, "/")
    print(f"[OK] {rel}")
    if not a.no_open:
        try:
            import webbrowser
            webbrowser.open("file://" + os.path.abspath(out))
        except Exception:
            pass


if __name__ == "__main__":
    main()
