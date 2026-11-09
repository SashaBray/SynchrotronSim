# -*- coding: utf-8 -*-
"""ПАНЕЛЬ ЭФФЕКТИВНОСТИ — качество превращения энергии пучка в полезное излучение.

Собирает в одном месте показатели КПД установки СИЛА и каждого ВКЛЮЧЁННОГО ондулятора (= присутствующего
в element_groups выбранной решётки). «КПД» у накопителя — не энергетический в буквальном смысле (электроны
не расходуются; ВЧ восполняет потери U₀·I каждый оборот, и бОльшая часть СИ — широкополосное дипольное
излучение в поглотители). Поэтому меряем КАЧЕСТВО превращения РАДИАЦИОННОЙ мощности:

УСТАНОВКА (в целом):
  • P_SR  = Σ мощность диполей + Σ мощность ондуляторов  (= полная излучаемая мощность = I·U₀);
  • U₀    = P_SR/I — потери энергии на оборот (их восполняет ВЧ);
  • η_ИД  = P_ИД / P_SR — доля излучаемой энергии, попадающая во вставные устройства (добываемая линиями),
            остальное — дипольный фон (волноотводы + тепло поглотителей);
  • η_полезн = ΣP_конус / P_SR — доля ПОЛНОЙ излучаемой мощности, попадающая в полезные центральные конусы
            ондуляторов (острое, используемое линиями излучение). Мала по сути — это честный предел.

ОНДУЛЯТОР (каждый включённый, при своём зазоре):
  • P_tot (кВт) — полная излучаемая мощность (тепловая нагрузка front-end, «цена» в энергии пучка);
  • плотность мощности на оси (кВт/мрад²) — пиковое тепло;
  • центральный конус n=1: полуугол θ_cen, мощность P_конус, доля η_конус=P_конус/P_tot (угловое качество);
  • поток в конусе и ЯРКОСТЬ n=1/3/5 на РЕАЛЬНЫХ β (из twiss) — главная фигура достоинства источника.
  Зазор общий (--gap) ИЛИ индивидуальный (--gaps und_01=8,und_19=12 или файл gid;gap). При общем зазоре и
  симметричной решётке все прямолинейки эквивалентны (β одинаковы); различие появляется при РАЗНЫХ зазорах
  (линии на разные энергии) — тогда у каждого ИД свои K, E₁, мощность, яркость. Качество полюсов задаётся
  --phase-error-deg (RMS σ_φ): высокие гармоники валятся ~exp(−(n·σ_φ)²), мощность из них уходит в фон.

ЗАПУСК
  py -3.14 undulators/efficiency.py --lattice 27 --gap 5
  py -3.14 undulators/efficiency.py --lattice 27 --gaps und_01=6,und_10=9,und_19=14 --phase-error-deg 4
"""
import os, sys, csv, math, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

import undulator_sr as sr
import bend_sr as bsr
from screen_analysis import load_groups, RESDIR
from bend_screen import load_lattice, dipole_arc
from screen_report import fig_div, kpi, section, _TPL, _LAYOUT, _CFG


# ─────────────────────────── физика эффективности ───────────────────────────
def cone_halfangle(K, n, p=sr.P):
    """Полуугол центрального конуса n-й гармоники, рад:  (1/γ)·√((1+K²/2)/(2nN))."""
    return (1.0 / sr.gamma(p)) * math.sqrt((1 + K ** 2 / 2) / (2 * n * p.N))


def cone_power_kW(B0, K, p=sr.P):
    """Мощность в центральном конусе n=1 ≈ плотность_на_оси · телесный угол конуса. (кВт, мрад, мрад²)."""
    th = cone_halfangle(K, 1, p) * 1e3                       # мрад
    omega = math.pi * th ** 2                                 # мрад² (конус малого полуугла)
    Pd0 = sr.power_density_onaxis(B0, K, p)                   # кВт/мрад²
    return Pd0 * omega, th, omega


def qphase(n, sigma_phi_deg):
    """Доля сохранённого пика n-й гармоники при фазовой ошибке σ_φ (°): exp(−(n·σ_φ)²)."""
    return math.exp(-(n * math.radians(sigma_phi_deg)) ** 2) if sigma_phi_deg else 1.0


def undulator_metrics(gid, lattice, gap, sig_phi, p=sr.P):
    """Показатели одного ондулятора при его зазоре (β и яркость — по реальной оптике места установки)."""
    from lattice_optics import undulator_optics
    o = undulator_optics(lattice, gid)
    p.beta_x, p.beta_y = o["betax"], o["betay"]
    B0 = sr.B0_of_gap(gap, p); K = sr.K_of_B0(B0, p)
    Pcone, th_cen, omega = cone_power_kW(B0, K, p)
    P_tot = sr.total_power_kW(B0, p)
    rec = dict(gid=gid, s=o["s"], betax=o["betax"], betay=o["betay"], gap=gap, B0=B0, K=K,
               E1=sr.harmonic_energy_keV(1, K, 0.0, p), P_tot=P_tot,
               Pd0=sr.power_density_onaxis(B0, K, p), th_cen_urad=th_cen * 1e3,
               Pcone=Pcone, eta_cone=Pcone / P_tot if P_tot else 0.0,
               flux1=sr.central_flux(1, K, p) * qphase(1, sig_phi))
    for n in (1, 3, 5):
        rec[f"bright{n}"] = sr.brightness(n, K, p) * qphase(n, sig_phi)
    return rec


def facility_metrics(lattice, default_gap, gap_over, E, I, sig_phi):
    sr.P.E, sr.P.I = E, I
    groups = load_groups(lattice)                            # {und_xx: {pole instance_ids}}
    if not groups:
        raise SystemExit(f"нет ондуляторов (element_groups) для решётки {lattice}")
    pole_insts = set().union(*groups.values())              # полюса ИД — НЕ считать поворотными магнитами
    ci, byid = load_lattice(lattice, None)

    # мощность поворотных магнитов (диполи + DQ; корректоры/квадруполи/… не светят или пренебрежимо)
    P_dip = 0.0; ndip = 0
    for iid, r in byid.items():
        if iid in pole_insts:
            continue
        d = dipole_arc(r, ci, E)
        if d is None:
            continue
        P_dip += bsr.total_power_kW(d["B"], d["L"], E, I)    # кВт
        ndip += 1

    unds = sorted(groups)
    unknown = set(gap_over) - set(unds)
    if unknown:
        raise SystemExit(f"--gaps: нет таких ондуляторов в реш.{lattice}: {', '.join(sorted(unknown))}")
    rows = [undulator_metrics(g, lattice, gap_over.get(g, default_gap), sig_phi) for g in unds]

    P_unds = sum(r["P_tot"] for r in rows)
    P_useful = sum(r["Pcone"] for r in rows)
    flux_total = sum(r["flux1"] for r in rows)
    P_SR = P_unds + P_dip
    U0_MeV = P_SR * 1e3 / I / 1e6 if I else 0.0             # P[Вт]/I = U₀[эВ]
    uniform = len({r["gap"] for r in rows}) == 1
    fac = dict(lattice=lattice, E=E, I=I, gamma=sr.gamma(), N_und=len(unds), n_dipoles=ndip,
               P_unds=P_unds, P_dip=P_dip, P_SR=P_SR, U0_MeV=U0_MeV, P_useful=P_useful,
               eta_id=P_unds / P_SR if P_SR else 0.0, eta_useful=P_useful / P_SR if P_SR else 0.0,
               flux_total=flux_total, sig_phi=sig_phi, uniform=uniform,
               gap=rows[0]["gap"] if uniform else None, K=rows[0]["K"] if uniform else None,
               E1=rows[0]["E1"] if uniform else None)
    return fac, rows


# ─────────────────────────── графики (plotly) ───────────────────────────
def funnel_fig(fac):
    f = go.Figure(go.Funnel(
        y=["Излучено всего (СИ кольца)", "В ондуляторах (ИД)", "В центральных конусах (полезное)"],
        x=[fac["P_SR"], fac["P_unds"], fac["P_useful"]],
        textposition="inside", textinfo="value+percent initial",
        marker=dict(color=["#6b7280", "#185fa5", "#eb6834"]),
        hovertemplate="%{y}<br>%{x:.3g} кВт<extra></extra>"))
    f.update_layout(height=300, **_LAYOUT)
    return f


def donut_fig(fac):
    f = go.Figure(go.Pie(labels=["Ондуляторы (ИД)", "Поворотные магниты"],
                         values=[fac["P_unds"], fac["P_dip"]], hole=0.52, sort=False,
                         marker=dict(colors=["#185fa5", "#1b7f5a"]),
                         texttemplate="%{label}<br>%{percent}",
                         hovertemplate="%{label}: %{value:.3g} кВт (%{percent})<extra></extra>"))
    f.update_layout(height=300, showlegend=False,
                    annotations=[dict(text=f"P_SR<br>{fac['P_SR']:.0f} кВт", showarrow=False, font=dict(size=13))],
                    **_LAYOUT)
    return f


def und_bars(rows, key, ylabel, color="#185fa5"):
    x = [r["gid"].replace("und_", "") for r in rows]
    f = go.Figure(go.Bar(x=x, y=[r[key] for r in rows], marker_color=color,
                         customdata=[[r["betax"], r["betay"], r["s"], r["gap"], r["E1"]] for r in rows],
                         hovertemplate="ИД %{x}<br>%{y:.3g}<br>зазор %{customdata[3]:.1f} мм, "
                                       "E₁=%{customdata[4]:.2f} кэВ<br>β=(%{customdata[0]:.2f}, "
                                       "%{customdata[1]:.2f}) м, s=%{customdata[2]:.1f} м<extra></extra>"))
    f.update_xaxes(title_text="ондулятор №")
    f.update_yaxes(title_text=ylabel)
    f.update_layout(height=320, **_LAYOUT)
    return f


# ─────────────────────────── вывод ───────────────────────────
def console(fac, rows):
    print(f"\nУСТАНОВКА СИЛА (реш.{fac['lattice']}): E={fac['E']} ГэВ, I={fac['I']*1e3:.0f} мА, "
          f"γ={fac['gamma']:.0f}; ондуляторов включено {fac['N_und']}, поворотных магнитов {fac['n_dipoles']}")
    if fac["uniform"]:
        print(f"  зазор ИД {fac['gap']:g} мм → B₀={rows[0]['B0']:.3f} Т, K={fac['K']:.3f}, E₁={fac['E1']:.2f} кэВ"
              + (f", качество σ_φ={fac['sig_phi']:g}°" if fac['sig_phi'] else ""))
    else:
        gs = sorted({r["gap"] for r in rows})
        print(f"  зазоры РАЗНЫЕ: {', '.join(f'{g:g}' for g in gs)} мм (см. таблицу по ондуляторам)"
              + (f"; качество σ_φ={fac['sig_phi']:g}°" if fac['sig_phi'] else ""))
    print(f"  полная излучаемая мощность P_SR = {fac['P_SR']:.1f} кВт  (U₀ = {fac['U0_MeV']:.3f} МэВ/оборот)")
    print(f"    из них в ондуляторах  {fac['P_unds']:.1f} кВт ({fac['eta_id']*100:.1f}%)  — η_ИД")
    print(f"    в поворотных магнитах {fac['P_dip']:.1f} кВт ({(1-fac['eta_id'])*100:.1f}%)")
    print(f"  в полезных центральных конусах ΣP = {fac['P_useful']:.2f} кВт  → η_полезн = {fac['eta_useful']*100:.3f}%")
    print(f"  суммарный поток n=1 (все ИД): {fac['flux_total']:.3e} ф/с/0.1%BW")
    print(f"\n  по ондуляторам:")
    print(f"  {'ИД':>6} | {'s,м':>7} | {'зазор':>5} | {'K':>5} | {'E₁,кэВ':>7} | {'P_tot':>6} | "
          f"{'η_кон':>6} | {'ярк n=1':>10} | {'ярк n=5':>10}")
    for r in rows:
        print(f"  {r['gid']:>6} | {r['s']:>7.1f} | {r['gap']:>5.1f} | {r['K']:>5.2f} | {r['E1']:>7.2f} | "
              f"{r['P_tot']:>6.3f} | {r['eta_cone']*100:>5.2f}% | {r['bright1']:>10.2e} | {r['bright5']:>10.2e}")


def write_csv(fac, rows, lattice):
    os.makedirs(RESDIR, exist_ok=True)
    fp = os.path.join(RESDIR, f"efficiency_facility_L{lattice}.csv")
    with open(fp, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["metric", "value", "unit"])
        for k, v, u in [
            ("E_beam", fac["E"], "GeV"), ("current", fac["I"], "A"), ("gamma", f"{fac['gamma']:.1f}", ""),
            ("gap_uniform_mm", fac["gap"] if fac["uniform"] else "mixed", "mm"),
            ("n_undulators", fac["N_und"], ""), ("n_dipoles", fac["n_dipoles"], ""),
            ("P_SR_total", f"{fac['P_SR']:.4f}", "kW"), ("U0_per_turn", f"{fac['U0_MeV']:.5f}", "MeV"),
            ("P_undulators", f"{fac['P_unds']:.4f}", "kW"), ("P_dipoles", f"{fac['P_dip']:.4f}", "kW"),
            ("P_useful_cones", f"{fac['P_useful']:.5f}", "kW"), ("eta_ID", f"{fac['eta_id']:.5f}", "fraction"),
            ("eta_useful", f"{fac['eta_useful']:.6f}", "fraction"),
            ("flux_total_n1", f"{fac['flux_total']:.5e}", "ph/s/0.1%BW"), ("sigma_phi", fac["sig_phi"], "deg")]:
            w.writerow([k, v, u])

    up = os.path.join(RESDIR, f"efficiency_undulators_L{lattice}.csv")
    with open(up, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["undulator", "s_m", "gap_mm", "B0_T", "K", "E1_keV", "betax_m", "betay_m", "P_tot_kW",
                    "Pdens_kW_mrad2", "theta_cone_urad", "P_cone_W", "eta_cone", "flux_n1",
                    "bright_n1", "bright_n3", "bright_n5"])
        for r in rows:
            w.writerow([r["gid"], f"{r['s']:.4f}", f"{r['gap']:.3f}", f"{r['B0']:.4f}", f"{r['K']:.4f}",
                        f"{r['E1']:.4f}", f"{r['betax']:.4f}", f"{r['betay']:.4f}", f"{r['P_tot']:.4f}",
                        f"{r['Pd0']:.4f}", f"{r['th_cen_urad']:.4f}", f"{r['Pcone']*1e3:.4f}",
                        f"{r['eta_cone']:.5f}", f"{r['flux1']:.4e}", f"{r['bright1']:.4e}",
                        f"{r['bright3']:.4e}", f"{r['bright5']:.4e}"])
    return fp, up


_EFF_NOTE = (
    '<div class="unitbox"><b>Что здесь «КПД».</b> У накопителя электроны не расходуются — ВЧ каждый оборот '
    'возвращает им потерянную на СИ энергию U₀·I. Поэтому энергетического КПД в буквальном смысле нет; меряем '
    '<b>качество превращения излучаемой мощности</b>. Воронка ниже: из всей излучаемой мощности P_SR часть '
    'рождается в ондуляторах (η_ИД), а из неё лишь малая доля попадает в полезные <b>центральные конусы</b> '
    '(острое используемое линиями излучение) — это η_полезн. Остальная мощность — широкополосный дипольный фон '
    '(волноотводы + тепло поглотителей) и внеосевые/межгармонические части ондуляторов. Низкое η_полезн — не '
    'дефект: настоящая фигура достоинства источника 4-го поколения — <b>яркость</b>, а не энергетический КПД.</div>')


def _rng(rows, key, scale=1.0, fmt="{:.3g}"):
    vals = [r[key] * scale for r in rows]
    lo, hi = min(vals), max(vals)
    return fmt.format(lo) if abs(hi - lo) < 1e-12 else f"{fmt.format(lo)} … {fmt.format(hi)}"


def per_undulator_section(rows, uniform):
    """Секция «по ондуляторам»: бар яркости (всегда) + бар E₁ при разных зазорах."""
    br = sorted(rows, key=lambda r: r["bright1"])
    best, worst = br[-1], br[0]
    bar_bright = fig_div(und_bars(rows, "bright1", "яркость n=1, ф/с/мрад²/мм²/0.1%"))
    if uniform:
        cap = (f"Общий зазор ⇒ мощность и поток у всех ИД одинаковы; различие — только β места установки. "
               f"В симметричной решётке прямолинейки эквивалентны (разброс яркости "
               f"{best['bright1']/worst['bright1']:.3f}×). Различие появится при разных зазорах или перекосах.")
        return section("Яркость по ондуляторам (n=1)", bar_bright, cap,
                       "Яркость — на РЕАЛЬНЫХ β из twiss; ε_x, связь, σ_E — помеченные допущения.")
    bar_e1 = fig_div(und_bars(rows, "E1", "E₁ (1-я гарм.), кэВ", color="#eb6834"))
    cap = (f"Разные зазоры ⇒ каждый ИД отдаёт свою энергию и мощность. Самый яркий {best['gid']} "
           f"(зазор {best['gap']:g} мм, E₁={best['E1']:.2f} кэВ) — {best['bright1']/worst['bright1']:.2f}× ярче "
           f"слабейшего {worst['gid']} (зазор {worst['gap']:g} мм).")
    return (section("Энергия 1-й гармоники по ондуляторам", bar_e1,
                    "Перестройка по зазору: каждая линия настроена на свою рабочую энергию.")
            + section("Яркость по ондуляторам (n=1)", bar_bright, cap,
                      "Яркость — на РЕАЛЬНЫХ β из twiss; ε_x, связь, σ_E — помеченные допущения."))


def build_html(fac, rows, lattice, out):
    kpis = "".join([
        kpi("Ондуляторов", f"{fac['N_und']}", "включено"),
        kpi("P_SR (всё СИ)", f"{fac['P_SR']:.0f}", "кВт"),
        kpi("U₀ / оборот", f"{fac['U0_MeV']:.2f}", "МэВ"),
        kpi("η_ИД", f"{fac['eta_id']*100:.1f}", "% мощности в ИД"),
        kpi("ΣP ондуляторов", f"{fac['P_unds']:.0f}", "кВт"),
        kpi("η_полезн", f"{fac['eta_useful']*100:.2f}", "% в конусах")])

    # правый столбец: общие значения (если зазор единый) или диапазоны по ИД
    if fac["uniform"]:
        thr = rows[0]
        und_tab = [("Полная мощность P_tot", f"{thr['P_tot']:.3f} кВт"),
                   ("Плотность на оси", f"{thr['Pd0']:.1f} кВт/мрад²"),
                   ("Полуугол центр. конуса", f"{thr['th_cen_urad']:.1f} мкрад (n=1)"),
                   ("Мощность в конусе", f"{thr['Pcone']*1e3:.1f} Вт"),
                   ("η_конус = P_конус/P_tot", f"{thr['eta_cone']*100:.2f} %"),
                   ("Поток в конусе n=1", f"{thr['flux1']:.2e} ф/с/0.1%BW"),
                   ("Яркость n=1", f"{thr['bright1']:.2e}")]
        und_head = "На один ондулятор (общий зазор)"
    else:
        und_tab = [("P_tot (диапазон)", f"{_rng(rows,'P_tot')} кВт"),
                   ("E₁ (диапазон)", f"{_rng(rows,'E1')} кэВ"),
                   ("η_конус (диапазон)", f"{_rng(rows,'eta_cone',100,'{:.2f}')} %"),
                   ("Мощн. в конусе (диапазон)", f"{_rng(rows,'Pcone',1e3,'{:.1f}')} Вт"),
                   ("Поток n=1 (диапазон)", f"{_rng(rows,'flux1')} ф/с/0.1%BW"),
                   ("Яркость n=1 (диапазон)", f"{_rng(rows,'bright1')}"),
                   ("Полезная мощн. ΣP_конус", f"{fac['P_useful']:.2f} кВт")]
        und_head = "По ондуляторам (разные зазоры)"

    gap_txt = f"{fac['gap']:g} мм (K={fac['K']:.3f}, E₁={fac['E1']:.2f} кэВ)" if fac["uniform"] else "разные (см. бары/CSV)"
    body = (
        f'<div class="kpis">{kpis}</div>'
        + _EFF_NOTE
        + section("Качество превращения: воронка мощности", fig_div(funnel_fig(fac)),
                  f"P_SR={fac['P_SR']:.0f} кВт → в ондуляторах {fac['P_unds']:.0f} кВт ({fac['eta_id']*100:.1f}%) "
                  f"→ в центральных конусах {fac['P_useful']:.2f} кВт ({fac['eta_useful']*100:.2f}% от всего). "
                  "Проценты — доля от полной излучаемой мощности.",
                  "Вложенные подмножества: полезный конус ⊂ излучение ИД ⊂ всё СИ кольца.")
        + section("Куда уходит излучаемая мощность", fig_div(donut_fig(fac)),
                  f"Ондуляторы {fac['P_unds']:.0f} кВт против поворотных магнитов {fac['P_dip']:.0f} кВт "
                  f"({fac['n_dipoles']} диполей/DQ). Дипольный фон частично используется волноотводами, "
                  "в основном — тепло поглотителей.")
        + per_undulator_section(rows, fac["uniform"])
        + f'<div class="tcols"><div><h3>Показатели установки</h3><table>'
        + "".join(f"<tr><td>{a}</td><td>{b}</td></tr>" for a, b in [
            ("Решётка / энергия / ток", f"{lattice} · {fac['E']} ГэВ · {fac['I']*1e3:.0f} мА"),
            ("Ондуляторов включено", f"{fac['N_und']}"), ("Зазор ИД", gap_txt),
            ("Поворотных магнитов", f"{fac['n_dipoles']} (диполи + DQ)"),
            ("P_SR (полная излучаемая)", f"{fac['P_SR']:.1f} кВт"),
            ("U₀ (потери/оборот)", f"{fac['U0_MeV']:.3f} МэВ"),
            ("η_ИД (доля в ИД)", f"{fac['eta_id']*100:.2f} %"),
            ("η_полезн (в центр. конусах)", f"{fac['eta_useful']*100:.3f} %"),
            ("Σ поток n=1 (все ИД)", f"{fac['flux_total']:.2e} ф/с/0.1%BW"),
            ("Качество полюсов σ_φ", f"{fac['sig_phi']:g}°")])
        + f'</table></div><div><h3>{und_head}</h3><table>'
        + "".join(f"<tr><td>{a}</td><td>{b}</td></tr>" for a, b in und_tab)
        + '</table></div></div>')

    html = (_TPL.replace("/*TITLE*/", f"Эффективность установки — СИЛА (реш.{lattice})")
            .replace("/*SUB*/", f"Качество превращения энергии пучка в полезное излучение · ток {fac['I']*1e3:.0f} мА "
                                f"· графики интерактивны (зум, пан, ховер)")
            .replace("/*BODY*/", body)
            .replace("/*PLOTLYJS*/", get_plotlyjs()))
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out, "w", encoding="utf-8").write(html)
    return out


def parse_gaps(spec):
    """'und_01=8,und_19=12' или путь к файлу 'gid;gap' → {gid: gap_mm}."""
    out = {}
    if not spec:
        return out
    if os.path.exists(spec):
        for r in csv.reader(open(spec, encoding="utf-8-sig"), delimiter=";"):
            if len(r) >= 2 and r[0].strip() and not r[0].lstrip().startswith("#"):
                try:
                    out[r[0].strip()] = float(str(r[1]).replace(",", "."))
                except ValueError:
                    pass
    else:
        for tok in spec.split(","):
            if "=" in tok:
                k, v = tok.split("=", 1)
                out[k.strip()] = float(v.strip().replace(",", "."))
    return out


def main():
    ap = argparse.ArgumentParser(description="панель эффективности: КПД установки и каждого ондулятора")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--gap", type=float, default=5.0, help="общий зазор ондуляторов, мм (если нет --gaps)")
    ap.add_argument("--gaps", default=None, help="индивидуальные зазоры: 'und_01=8,und_19=12' или файл gid;gap")
    ap.add_argument("--E", type=float, default=6.0); ap.add_argument("--I", type=float, default=0.2)
    ap.add_argument("--kappa", type=float, default=sr.P.kappa, help="связь εy/εx (допущение)")
    ap.add_argument("--sigE", type=float, default=sr.P.sig_E, help="σ_E/E (допущение)")
    ap.add_argument("--phase-error-deg", type=float, default=0.0, help="RMS фазовая ошибка полюсов σ_φ, °")
    ap.add_argument("--no-html", action="store_true")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    sr.P.kappa, sr.P.sig_E = a.kappa, a.sigE
    gap_over = parse_gaps(a.gaps)

    print(f"панель эффективности: решётка {a.lattice}, ток {a.I*1e3:.0f} мА"
          + (f", зазоры индивидуальные ({len(gap_over)})" if gap_over else f", зазор {a.gap:g} мм") + " …", flush=True)
    fac, rows = facility_metrics(a.lattice, a.gap, gap_over, a.E, a.I, a.phase_error_deg)
    console(fac, rows)
    fp, up = write_csv(fac, rows, a.lattice)
    print(f"\n[OK] {fp}\n[OK] {up}")
    if not a.no_html:
        out = a.out or os.path.join(RESDIR, f"efficiency_L{a.lattice}.html")
        build_html(fac, rows, a.lattice, out)
        print(f"[OK] {out}")


if __name__ == "__main__":
    main()
