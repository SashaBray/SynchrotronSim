# -*- coding: utf-8 -*-
"""element_bench.py — ЕДИНЫЙ ПАЙПЛАЙН сравнения траекторий через элемент в разных МОДЕЛЯХ (ядрах).
Развиваемый инструмент, а не одноразовые скрипты. Задача — понять, как получить проектные физкоэф.
через токи, сравнивая ядра на одинаковых опытных конфигурациях.

Опытная конфигурация = старт → входной дрейф → ЭЛЕМЕНТ(токи обмоток) → выходной дрейф (напр. 10 м) +
probe-пучок. Хранится в БД (bench_experiments.csv). Прогон в РАЗНЫХ ядрах, траектории сохраняются
с ПОМЕТКОЙ ЯДРА (bench/traj_<exp>.csv: core;particle;S;u;v;…). Визуализация — накладывает все ядра
на один график (PNG + вращаемый HTML) ОДНОВРЕМЕННО.

Реестр моделей MODELS (расширяемый одной функцией):
  diff    — Boris по 3D-картам поля (суперпозиция обмоток × токи). ЭТАЛОН.
  kick    — matrix_tracker (токи → матрица связи → физкоэф → кики). Наше киковое ядро.
  elegant — elegant CSBEND (Форест-Рут + станд. ядро) на тех же физкоэф из матрицы связи.
Физкоэф. для kick/elegant берутся ИЗ МАТРИЦЫ СВЯЗИ по токам (по умолчанию) — сравнение с diff и
показывает, верна ли связь. Можно задать coeffs вручную в конфиге (проверить другие значения).

CLI:
  py element_bench.py create <name> --dev DQ_I1_1A=50,DQ_I2_1A=0,DQ_I3_1A=0 --entry 0.5 --exit 10
  py element_bench.py run <exp> [--cores diff,kick,elegant]
  py element_bench.py viz <exp>                     # 2D PNG + вращаемый HTML (с подписью: элемент/токи/физкоэф.)
  py element_bench.py viz3d <exp> [--view persp|top|side|full] [--shot out.png]  # в общем 3D-визуализаторе
  py element_bench.py list
"""
import os, sys, csv, json, math, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
sys.path.insert(0, "C:/Sasha/Work/SILA/SynchrotronSim_archive/legacy/scripts")
import elegant_model as elm
import matrix_tracker as mt
import dq_sector_profile as dsp
from bip_matrix import QP

TB = os.path.join("workspace", "tables")
BENCH = os.path.join("workspace", "research", "bench"); os.makedirs(BENCH, exist_ok=True)
CFG_CSV = os.path.join(TB, "bench_experiments.csv")
HW = {r["device_type_id"]: r for r in mt._read("hardware_library.csv")}
CM = {r["device_type_id"]: r for r in mt._read("coupling_matrix.csv")}
EFF, PHYS, COUP, TRANS = mt.load_db()
CFG_COLS = ["exp_id", "name", "windings", "entry_m", "exit_m", "probes", "coeffs", "notes"]
DEF_PROBES = [[0, 0, 0, 0], [2e-3, 0, 0, 0], [-2e-3, 0, 0, 0], [4e-3, 0, 0, 0], [-4e-3, 0, 0, 0]]

# ------------------------------------------------------------------ конфиги (БД) --
def load_cfgs():
    if not os.path.exists(CFG_CSV):
        return {}
    return {r["exp_id"]: r for r in csv.DictReader(open(CFG_CSV, encoding="utf-8-sig"), delimiter=";")}


def save_cfg(cfg):
    cfgs = load_cfgs(); cfgs[cfg["exp_id"]] = cfg
    with open(CFG_CSV, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, CFG_COLS, delimiter=";"); w.writeheader()
        for c in cfgs.values():
            w.writerow(c)


def create(name, windings, entry=0.5, exit=10.0, probes=None, coeffs=None, notes=""):
    cfgs = load_cfgs(); eid = f"e{len(cfgs)+1}"
    while eid in cfgs:
        eid = f"e{int(eid[1:])+1}"
    cfg = dict(exp_id=eid, name=name, windings=json.dumps(windings), entry_m=str(entry), exit_m=str(exit),
               probes=json.dumps(probes or DEF_PROBES), coeffs=json.dumps(coeffs or {}), notes=notes)
    save_cfg(cfg)
    print(f"[создан] {eid}: {name}  обмотки={windings}  вход={entry}м выход={exit}м")
    return eid


def _parse(cfg):
    return (json.loads(cfg["windings"]), float(cfg["entry_m"]), float(cfg["exit_m"]),
            json.loads(cfg["probes"]), json.loads(cfg.get("coeffs") or "{}"))

# ------------------------------------------------------- поле/магнит из обмоток --
def build_field(windings):
    """Суперпозиция полей обмоток × токи → (field_fn, half)."""
    fns = []; half = None
    for dev, cur in windings.items():
        h, B, nom = dsp._load_full(dev, HW[dev]); fns.append((dsp._field_fn(h, B), cur / nom)); half = float(h["mxx"])
    return (lambda p: sum(w * f(p) for f, w in fns)), half


def build_magnet(windings, coeffs):
    """Magnet matrix_tracker. coeffs пуст → source=current (токи через матрицу связи);
    иначе source=physical с заданными mat_*."""
    dev0 = next(iter(windings))
    if coeffs:
        inst = {"device_type_id": dev0, "x": "0", "y": "0", "z": "0", "yaw": "0", "pitch": "0", "roll": "0",
                "instance_id": "1", "arg_val": "0", "correct_coeff": "1",
                "mat_B": str(coeffs.get("B", 0)), "mat_G": str(coeffs.get("G", 0)),
                "mat_S": str(coeffs.get("S", 0)), "mat_O": str(coeffs.get("O", 0))}
        return mt.Magnet([inst], EFF, PHYS, COUP, TRANS, source="physical")
    insts = [{"device_type_id": d, "x": "0", "y": "0", "z": "0", "yaw": "0", "pitch": "0", "roll": "0",
              "instance_id": str(i + 1), "arg_val": str(c), "correct_coeff": "1"} for i, (d, c) in enumerate(windings.items())]
    return mt.Magnet(insts, EFF, PHYS, COUP, TRANS, source="current")

# ------------------------------------------------------------------- трекеры --
def _drift(pts, s, u, v, up, vp, L, n=2):
    for t in np.linspace(0, L, n)[1:]:
        pts.append((s + t, u + up * t, v + vp * t))
    return s + L, u + up * L, v + vp * L


def track_diff(windings, entry, exit, probe):
    field, half = build_field(windings); u, v, up, vp = probe
    pts = [(0.0, u, v)]; s = 0.0
    s, u, v = _drift(pts, s, u, v, up, vp, entry)             # входной дрейф
    x = -half; ds = 2e-4; n = int(round(2 * half / ds))       # зона поля — Boris
    for _ in range(n):
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        B = field(np.array([x, u, v])); up += -QP * B[2] * ds; vp += QP * B[1] * ds
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        s += ds; pts.append((s, u, v))
    _drift(pts, s, u, v, up, vp, exit)                        # выходной дрейф
    return np.array(pts)


def track_kick(windings, coeffs, entry, exit, probe):
    m = build_magnet(windings, coeffs); half = _half(windings); u0, v0, up0, vp0 = probe
    N = m.nsl
    dsg = m.L / N; dsf = m.Leff / N                          # РИНГ-ядро: дрейф по L_phys, кик по Leff
    th, k1, k1v, k2, k3 = m.theta0, m.k1, m.k1v, m.k2, m.k3
    lead = max(0.0, (2 * half - m.L) / 2)                    # центрируем физ.элемент (L_phys) в зоне поля 2·half
    pts = [(0.0, u0, v0)]; s = 0.0
    s, uu, vv = _drift(pts, s, u0, v0, up0, vp0, entry + lead); st = (uu, up0, vv, vp0, 0.0, 0.0)
    for _ in range(N):
        st = mt._drift(st, dsg / 2.0); s += dsg / 2.0
        x, px, y, py, t, d = st
        px += th / N - k1 * x * dsf - 0.5 * k2 * (x * x - y * y) * dsf - (k3 / 6.0) * (x**3 - 3 * x * y * y) * dsf
        py += k1v * y * dsf + k2 * x * y * dsf + (k3 / 6.0) * (3 * x * x * y - y**3) * dsf
        st = mt._drift((x, px, y, py, t, d), dsg / 2.0); s += dsg / 2.0
        pts.append((s, st[0], st[2]))
    _drift(pts, s, st[0], st[2], st[1], st[3], lead + exit)
    return np.array(pts)


def track_elegant(windings, coeffs, entry, exit, probe):
    field, half = build_field(windings)
    m = build_magnet(windings, coeffs); Leff = m.Leff
    th, k1, k2, k3 = m.theta0, m.k1, m.k2, m.k3
    Xo, Uo = _orbit(field, half)                              # опорная дуга — из diff на оси
    u0, v0, up0, vp0 = probe; N = 40; ds = Leff / N           # elegant CSBEND: дрейф=кик по Leff
    lead = max(0.0, (2 * half - Leff) / 2)                    # центрируем элемент (Leff) в зоне поля
    pts = [(0.0, u0, v0)]; s = 0.0
    s, u, v = _drift(pts, s, u0, v0, up0, vp0, entry + lead)
    x = -Leff / 2.0; qx, qy = elm._slopes_to_momenta(up0, vp0)   # X от −Leff/2 (центр зоны)
    for _ in range(N):
        for stp in range(4):
            up, vp = elm._momenta_to_slopes(qx, qy); dg = ds * elm.DRIFT4[stp]
            x += dg; u += up * dg; v += vp * dg
            if elm.KICK4[stp] == 0.0:
                break
            du = u - float(np.interp(x, Xo, Uo)); duc = min(0.03, max(-0.03, du)); vc = min(0.03, max(-0.03, v))
            kf = ds * elm.KICK4[stp]                          # конвенция matrix_tracker: +θ0/L, −k1·du, …
            qx += (th / Leff - k1 * duc - 0.5 * k2 * (duc * duc - vc * vc) - (k3 / 6.0) * (duc**3 - 3 * duc * vc * vc)) * kf
            qy += (k1 * vc + k2 * duc * vc + (k3 / 6.0) * (3 * duc * duc * vc - vc**3)) * kf
        s += ds; pts.append((s, u, v))
    up, vp = elm._momenta_to_slopes(qx, qy)
    _drift(pts, s, u, v, up, vp, lead + exit)
    return np.array(pts)


def _half(windings):
    h, _, _ = dsp._load_full(next(iter(windings)), HW[next(iter(windings))]); return float(h["mxx"])


def coeff_summary(windings, coeffs):
    """Сводка: элемент, токи обмоток и ФИЗКОЭФ. + кики, которыми отработали kick/elegant.
    Физкоэф. берутся из собранного Magnet (source=current → из матрицы связи; coeffs → заданные)."""
    m = build_magnet(windings, coeffs); Leff, L = m.Leff, m.L; QPl = QP * Leff
    phys = {"B": (m.theta0 / QPl if QPl else 0.0), "G": m.k1 / QP, "S": m.k2 / QP, "O": m.k3 / QP}
    kick = {"θ0": m.theta0, "k1": m.k1, "k2": m.k2, "k3": m.k3}
    dev0 = next(iter(windings)); kind = TRANS.get(dev0, {}).get("kind", "?")
    return dict(devices=list(windings), kind=kind, currents=windings, Leff=Leff, L=L,
                phys=phys, kick=kick, src=("вручную (coeffs)" if coeffs else "матрица связи ← токи"))


def summary_lines(windings, coeffs):
    """Список строк сводки для подписи (PNG/HTML). Скью-строка — только если скью ≠ 0."""
    s = coeff_summary(windings, coeffs); p, k = s["phys"], s["kick"]
    cur = ", ".join(f"{d}={c:g} А" for d, c in s["currents"].items())
    lines = [f"Элемент: {'+'.join(s['devices'])}  ({s['kind']}) | Leff={s['Leff']:.3f} м, L_phys={s['L']:.3f} м",
             f"Токи обмоток: {cur}",
             f"Физкоэф. ({s['src']}):  B={p['B']:+.3e} Тл   G={p['G']:+.3e} Тл/м   "
             f"S={p['S']:+.3e} Тл/м²   O={p['O']:+.3e} Тл/м³"]
    if any(abs(p[c]) > 1e-12 for c in ("Bs", "Gs", "Ss", "Os")):
        lines.append(f"   скью:  Bs={p['Bs']:+.3e} Тл   Gs={p['Gs']:+.3e} Тл/м   "
                     f"Ss={p['Ss']:+.3e} Тл/м²   Os={p['Os']:+.3e} Тл/м³")
    lines.append(f"Кики (интеграл): θ0={k['θ0'] * 1e3:+.4g} мрад   θy={k['θy'] * 1e3:+.4g} мрад   "
                 f"k1={k['k1']:+.4g} м⁻²   k2={k['k2']:+.4g} м⁻³   k3={k['k3']:+.4g} м⁻⁴")
    return lines


def _orbit(field, half, ds=1e-3):
    x = -half; u = v = up = vp = 0.0; Xs = [x]; Us = [u]
    for _ in range(int(round(2 * half / ds))):
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps)
        B = field(np.array([x, u, 0.0])); up += -QP * B[2] * ds
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps)
        Xs.append(x); Us.append(u)
    return np.array(Xs), np.array(Us)


MODELS = {
    "diff": lambda w, c, en, ex, p: track_diff(w, en, ex, p),
    "kick": lambda w, c, en, ex, p: track_kick(w, c, en, ex, p),
    "elegant": lambda w, c, en, ex, p: track_elegant(w, c, en, ex, p),
}

# ------------------------------------------------------------------- прогон --
def run(exp_id, cores=("diff", "kick", "elegant")):
    cfg = load_cfgs()[exp_id]; windings, entry, exit, probes, coeffs = _parse(cfg)
    print(f"[прогон] {exp_id} ({cfg['name']}): ядра {list(cores)}, probe {len(probes)}")
    rows = []
    for core in cores:
        for pi, probe in enumerate(probes):
            tr = MODELS[core](windings, coeffs, entry, exit, probe)
            for (S, u, v) in tr:
                rows.append([core, pi, f"{S:.5f}", f"{u:.6e}", f"{v:.6e}"])
    p = os.path.join(BENCH, f"traj_{exp_id}.csv")
    with open(p, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(["core", "particle", "S_m", "u_m", "v_m"]); w.writerows(rows)
    print(f"[OK] траектории с пометкой ядра → {p} ({len(rows)} точек)")
    return p

# ------------------------------------------------------------- визуализация --
CORE_STYLE = {"diff": ("-", 2.4, 0.5), "kick": ("--", 1.1, 1.0), "elegant": (":", 1.4, 1.0)}
CORE_HTML = {"diff": "solid", "kick": "dash", "elegant": "dot"}


def visualize(exp_id):
    cfg = load_cfgs()[exp_id]; windings, entry, exit, _, coeffs = _parse(cfg)
    info = summary_lines(windings, coeffs)                    # элемент, токи, физкоэф., кики
    p = os.path.join(BENCH, f"traj_{exp_id}.csv")
    R = list(csv.DictReader(open(p, encoding="utf-8-sig"), delimiter=";"))
    cores = sorted(set(r["core"] for r in R), key=lambda c: list(CORE_STYLE).index(c) if c in CORE_STYLE else 9)
    parts = sorted(set(int(r["particle"]) for r in R)); cmap = plt.cm.viridis(np.linspace(0, 0.9, len(parts)))
    half = _half(windings)
    fig, ax = plt.subplots(1, 2, figsize=(15, 6))
    for ci, core in enumerate(cores):
        ls, lw, al = CORE_STYLE.get(core, ("-", 1, 1))
        for pi in parts:
            d = [(float(r["S_m"]), float(r["u_m"]), float(r["v_m"])) for r in R if r["core"] == core and int(r["particle"]) == pi]
            d = np.array(d)
            ax[0].plot(d[:, 0], d[:, 1] * 1e3, ls, lw=lw, alpha=al, color=cmap[pi])
            ax[1].plot(d[:, 0], d[:, 2] * 1e3, ls, lw=lw, alpha=al, color=cmap[pi])
    for a, lab in zip(ax, ["u (гор.)", "v (верт.)"]):
        a.axvspan(entry, entry + 2 * half, color="0.92", zorder=0); a.grid(alpha=0.3)
        a.set_xlabel("S вдоль оси, м"); a.set_ylabel(f"{lab}, мм")
    from matplotlib.lines import Line2D
    leg = [Line2D([0], [0], color="k", ls=CORE_STYLE.get(c, ("-",))[0], lw=2, label=c) for c in cores]
    ax[0].legend(handles=leg, title="ядро", loc="best")
    ax[0].set_title(f"{cfg['name']}  (обмотки {json.loads(cfg['windings'])})"); ax[1].set_title("серая зона = элемент")
    fig.suptitle(f"element_bench {exp_id}: сравнение ядер (одновременно)", weight="bold")
    fig.text(0.5, 0.015, "\n".join(info), ha="center", va="bottom", fontsize=8.5, family="monospace",
             bbox=dict(boxstyle="round", fc="#f7f7f7", ec="0.7"))          # подпись: элемент/токи/физкоэф.
    png = os.path.join(BENCH, f"bench_{exp_id}.png"); fig.tight_layout(rect=[0, 0.14, 1, 0.97])
    fig.savefig(png, dpi=130); plt.close(fig)
    # интерактивный HTML — с той же сводкой (заголовок = элемент; блок снизу = токи/физкоэф./кики)
    figh = go.Figure()
    for core in cores:
        for pi in parts:
            d = np.array([(float(r["S_m"]), float(r["u_m"]) * 1e3, float(r["v_m"]) * 1e3)
                          for r in R if r["core"] == core and int(r["particle"]) == pi])
            figh.add_trace(go.Scatter(x=d[:, 0], y=d[:, 1], mode="lines", name=f"{core} p{pi}",
                           legendgroup=core, line=dict(dash=CORE_HTML.get(core, "solid")),
                           customdata=d[:, 2], hovertemplate="S=%{x:.3f} м<br>u=%{y:.4f} мм<br>v=%{customdata:.4f} мм<extra></extra>"))
    figh.update_layout(title=f"element_bench {exp_id}: {cfg['name']} — u(S) по ядрам<br><sub>{info[0]}</sub>",
                       xaxis_title="S, м", yaxis_title="u (гор.), мм", width=1150, height=760,
                       margin=dict(b=170), legend=dict(title="ядро · частица"))
    figh.add_annotation(text="<br>".join(info[1:]), xref="paper", yref="paper", x=0.0, y=-0.20,
                        showarrow=False, align="left", xanchor="left", yanchor="top",
                        font=dict(family="monospace", size=12), bordercolor="#999", borderwidth=1,
                        borderpad=6, bgcolor="#f7f7f7")
    html = os.path.join(BENCH, f"bench_{exp_id}.html"); figh.write_html(html, include_plotlyjs="inline")
    print(f"[OK] {png}\n[OK] {html} (вращаемый/зумируемый)")
    return png


CORE_COL3D = {"diff": "black", "kick": "crimson", "elegant": "deepskyblue"}
CORE_LW3D = {"diff": 5, "kick": 3, "elegant": 3}


def viz3d(exp_id, shot=None, view="persp"):
    """Опыт в ОБЩЕМ 3D-визуализаторе синхротрона: STL-геометрия элемента (как в кольце) + пучок
    по всем ядрам (diff/kick/elegant) в глобальных координатах X=S, Y=u(гор.), Z=v(верт.).
    Переиспользует visualize.build_combined_mesh и стиль сцены. shot=путь → off-screen PNG."""
    import pyvista as pv
    import visualize as vz
    cfg = load_cfgs()[exp_id]; windings, entry, exit, _, coeffs = _parse(cfg)
    info = summary_lines(windings, coeffs)
    half = _half(windings); Xc = entry + half                # центр зоны поля по оси X
    R = list(csv.DictReader(open(os.path.join(BENCH, f"traj_{exp_id}.csv"), encoding="utf-8-sig"), delimiter=";"))
    if not R:
        print(f"[!] нет траекторий — сначала: py element_bench.py run {exp_id}"); return
    cores = sorted(set(r["core"] for r in R), key=lambda c: list(CORE_STYLE).index(c) if c in CORE_STYLE else 9)
    parts = sorted(set(int(r["particle"]) for r in R))
    geo = os.path.join("workspace", "geometry")

    plotter = pv.Plotter(off_screen=bool(shot), title=f"element_bench {exp_id} — {cfg['name']} (3D)")
    plotter.set_background(vz.VIZ_CFG.get("background", "white"))
    # --- геометрия элемента: STL (silver, как в кольце) если есть, + полупрозрачная зона поля ---
    placed = []
    for dev in windings:
        stl = str(HW[dev].get("stl_path", "-")).strip()
        path = os.path.join(geo, stl)
        if stl not in ("-", "", "nan") and os.path.exists(path) and dev not in placed:
            m = pv.read(path); m.scale(0.001, inplace=True)   # мм→м
            big = vz.build_combined_mesh(m, np.array([[0.0, 0.0, 0.0, Xc, 0.0, 0.0]]))  # yaw=0, центр Xc
            plotter.add_mesh(big, color="silver", opacity=0.30, show_edges=True, label=f"{dev} (STL)")
            placed.append(dev)
    box = pv.Box(bounds=(entry, entry + 2 * half, -0.03, 0.03, -0.03, 0.03))
    # без STL — подписываем коробку именем прибора (латиница, чтобы легенда PyVista не искажала кириллицу)
    plotter.add_mesh(box, color="#8c8c8c", opacity=0.12, label=(None if placed else f"{next(iter(windings))} (field zone)"))
    # --- пучок по ядрам: X=S, Y=u, Z=v; цвет по ядру, легенда один раз на ядро ---
    for core in cores:
        first = True
        for pi in parts:
            d = np.array([(float(r["S_m"]), float(r["u_m"]), float(r["v_m"]))
                          for r in R if r["core"] == core and int(r["particle"]) == pi])
            if len(d) < 2:
                continue
            poly = pv.PolyData(d)
            cells = np.full((1, len(d) + 1), len(d), dtype=np.int_); cells[0, 1:] = np.arange(len(d)); poly.lines = cells
            plotter.add_mesh(poly, color=CORE_COL3D.get(core, "gray"), line_width=CORE_LW3D.get(core, 3),
                             label=f"ядро: {core}" if first else None)
            first = False
    # --- подпись (элемент/токи/физкоэф.) + оси/камера ---
    plotter.add_text("\n".join(info), position="lower_left", font_size=8, color="black", font="courier")
    plotter.show_bounds(grid="back", location="outer", xtitle="X=S (м)", ytitle="Y=u гор. (м)",
                        ztitle="Z=v верт. (м)", font_size=10, color="black")
    plotter.add_axes()
    try:
        plotter.add_legend(bcolor=None)
    except Exception:
        pass
    allx = np.array([float(r["S_m"]) for r in R]); cx = float((allx.min() + allx.max()) / 2); span = float(allx.max() - allx.min()) or 2.0
    win = max(2 * half * 4, 2.5)                              # окно вокруг элемента (element-виды): неск. длин магнита
    # ПЕРСПЕКТИВНАЯ проекция: (направление focus→камера, фокус, up, размер кадра). Дистанция считается
    # из угла обзора, чтобы кадр вмещал ~framesize. persp/top/side наведены на ЭЛЕМЕНТ; full — весь веер.
    VIEWS = {"persp": ((-0.6, -0.7, 0.55), (Xc, 0, 0), (0, 0, 1), win),
             "top":   ((0.0, 0.0, 1.0),    (Xc, 0, 0), (0, 1, 0), win),
             "side":  ((0.0, -1.0, 0.05),  (Xc, 0, 0), (0, 0, 1), win),
             "full":  ((-0.35, -0.6, 0.4), (cx, 0, 0), (0, 0, 1), span * 1.1)}
    dv, focus, up, fs = VIEWS.get(view, VIEWS["persp"])
    dv = np.array(dv, float); dv /= np.linalg.norm(dv); va = 30.0
    dist = (fs * 0.5) / np.tan(np.radians(va * 0.5))         # дистанция под перспективу (кадр ≈ fs)
    plotter.camera_position = [tuple(np.array(focus, float) + dv * dist), tuple(focus), up]
    plotter.camera.view_angle = va                           # ПЕРСПЕКТИВА (parallel projection НЕ включаем)
    if shot:
        plotter.screenshot(shot, window_size=(1600, 900)); plotter.close()
        print(f"[OK] 3D-снимок → {shot}")
    else:
        print(f"[3D] {cfg['name']}: элемент {'+'.join(windings)} + пучок по ядрам {cores}. Окно PyVista (крутить мышью).")
        plotter.show()


def _cli_windings(s):
    return {kv.split("=")[0]: float(kv.split("=")[1]) for kv in s.split(",")}


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd")
    c = sub.add_parser("create"); c.add_argument("name"); c.add_argument("--dev", required=True)
    c.add_argument("--entry", type=float, default=0.5); c.add_argument("--exit", type=float, default=10.0); c.add_argument("--notes", default="")
    r = sub.add_parser("run"); r.add_argument("exp"); r.add_argument("--cores", default="diff,kick,elegant")
    v = sub.add_parser("viz"); v.add_argument("exp")
    v3 = sub.add_parser("viz3d"); v3.add_argument("exp")
    v3.add_argument("--shot", default=None, help="путь PNG (off-screen); без него — интерактивное окно")
    v3.add_argument("--view", default="persp", choices=["persp", "top", "side", "full"])
    sub.add_parser("list")
    a = ap.parse_args()
    if a.cmd == "create":
        create(a.name, _cli_windings(a.dev), a.entry, a.exit, notes=a.notes)
    elif a.cmd == "run":
        run(a.exp, tuple(a.cores.split(",")))
    elif a.cmd == "viz":
        visualize(a.exp)
    elif a.cmd == "viz3d":
        viz3d(a.exp, shot=a.shot, view=a.view)
    elif a.cmd == "list":
        for e, c in load_cfgs().items():
            print(f"  {e}: {c['name']}  обмотки={c['windings']} вход={c['entry_m']} выход={c['exit_m']}")
    else:
        ap.print_help()
