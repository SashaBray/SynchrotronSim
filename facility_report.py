# -*- coding: utf-8 -*-
"""Паспорт установки — самодостаточный HTML-отчёт с полным анализом решётки.

Собирает в одну прокручиваемую страницу всё, что известно об установке из ОДНОГО описания решётки
(без сверки с elegant): геометрию, линейную оптику, хроматичность, излучение и равновесие, нелинейную
динамику (динамическая апертура, карта частот, сдвиг частот с амплитудой, импульсная апертура) и
таблицу сил по семействам магнитов. Переиспользует продакшн-расчёты:
  ring_matrix     — линейные матрицы, функции Твисса, интегралы излучения I₁…I₅;
  ring_map        — хроматизм высших порядков ξ₁…ξ₃, αc2, нормальная форма (∂ν/∂J);
  ring_cpp        — трекинг на ядре C++ (динамическая апертура, карта частот, импульсная апертура);
  ring_analysis   — частоты по матрице/NAFF, диффузия, линии резонансов.
HTML-каркас (прокрутка, тянущиеся по высоте панели) — из traj_analyzer.

Запуск:
  py -3.14 facility_report.py 26                 → workspace/results/html/facility_exp26.html
  py -3.14 facility_report.py 26 --fast          только аналитика (без трекинга: быстро)
  py -3.14 facility_report.py 2 --turns 1024     длиннее прогон динамической апертуры/карты частот
"""
import argparse
import math
import os
import sys
import webbrowser

import numpy as np
import plotly.graph_objects as go

for _s in (sys.stdout, sys.stderr):                           # греческие буквы в лог не должны падать на cp1251
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass

import errors as errmod                                       # ошибки установки: допуски + возмущённая решётка
import ring_analysis as ra
import ring_cpp
import ring_map as rmap
import ring_matrix as rm
import traj_analyzer as ta                                    # HTML-каркас: _assemble, _grid, COLORS, GRID

RESULTS = os.path.join("workspace", "results")
HTML_DIR = os.path.join(RESULTS, "html")
C_LIGHT = 299792458.0
MC2_MEV = 0.51099895069

# Параметры установок по id решётки: энергия, ВЧ, число суперпериодов и диапазоны сканов нелинейной
# динамики. Диапазоны — в мм по амплитуде и % по импульсу; порог потери — в метрах.
LATTICE_INFO = {
    25: dict(label="СИЛА 2025.08.29", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
    26: dict(label="СИЛА 2026.06.11", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
    2:  dict(label="Сибирь", energy_mev=3000.0, rf=(2.5, 75), superperiods=6,
             # секступоли нулевые: устойчивость теряется только от кинематики, на ~0.4–0.5 м —
             # поэтому предел пробы амплитуды сотни мм, иначе карта выживания однотонная (всё выжило)
             x_probe_mm=600.0, y_probe_mm=800.0, d_probe_pct=15.0, loss_m=1.0),
    1:  dict(label="СИЛА (реш. 1)", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
}
KIND_RU = {"dipole": "диполь", "dq": "диполь-квадруполь", "quad": "квадруполь", "sext": "секступоль",
           "oct": "октуполь", "kicker": "корректор", "monitor": "датчик положения"}


def info_of(lattice):
    if lattice not in LATTICE_INFO:
        raise SystemExit(f"нет описания решётки {lattice} в LATTICE_INFO — добавьте энергию/ВЧ/суперпериоды")
    return LATTICE_INFO[lattice]


def _ring(inf, nsl=40, order=4, rf=False, ws=None):
    kw = dict(lattice=inf["_id"], design_mev=inf["energy_mev"], nsl=nsl, order=order)
    if rf:
        kw["rf"] = dict(volt_mv=inf["rf"][0], harmonic=inf["rf"][1])
    if ws:                                                    # возмущённая рабочая папка (ансамбль ошибок)
        kw["ws"] = ws
    return ring_cpp.CppRing(**kw)


# ------------------------------------------------------------------ РАСЧЁТ ПАРАМЕТРОВ
def geometry(ring, inf):
    """Паспортные величины машины из геометрии и энергии."""
    C = ring.circumference
    E_gev = inf["energy_mev"] / 1000.0
    gamma = inf["energy_mev"] / MC2_MEV
    brho = math.sqrt(inf["energy_mev"] ** 2 - MC2_MEV ** 2) / 1000.0 / C_LIGHT * 1e9  # p[GeV]/c → Тл·м
    T0 = C / C_LIGHT
    f0 = 1.0 / T0
    volt, harm = inf["rf"]
    kinds = {}
    theta = 0.0
    for m in ring.mags:
        kinds[m.kind] = kinds.get(m.kind, 0) + 1
        theta += abs(float(getattr(m, "theta0", 0.0) or 0.0))
    return dict(circumference=C, energy_mev=inf["energy_mev"], gamma=gamma, brho=brho, T0=T0, f0=f0,
                superperiods=inf["superperiods"], rf_volt_mv=volt, rf_harmonic=harm, rf_freq=harm * f0,
                kinds=kinds, n_mag=sum(v for k, v in kinds.items() if k != "monitor"),
                sum_angle=theta, sum_angle_2pi=theta / (2 * math.pi))


def linear_optics(ring, inf):
    """Функции Твисса вдоль кольца + сводные величины."""
    lm = rm.line_matrices(ring)
    tw = rm.twiss(ring, lm)
    n = inf["superperiods"]
    return dict(s=np.asarray(tw["s"]), betax=np.asarray(tw["betax"]), betay=np.asarray(tw["betay"]),
                alphax=np.asarray(tw["alphax"]), alphay=np.asarray(tw["alphay"]),
                etax=np.asarray(tw["etax"]),
                nux=float(tw["nux"]), nuy=float(tw["nuy"]),
                mux_cell=float(tw["nux"]) / n, muy_cell=float(tw["nuy"]) / n,
                betax0=float(tw["betax"][0]), betay0=float(tw["betay"][0]),
                alphax0=float(tw["alphax"][0]), alphay0=float(tw["alphay"][0]),
                etax0=float(tw["etax"][0]), betax_max=float(np.max(tw["betax"])),
                betay_max=float(np.max(tw["betay"])), etax_max=float(np.max(tw["etax"])),
                alphac=float(tw["alphac"]), closed_orbit=np.asarray(tw["closed_orbit"]))


def chroma(ring):
    """Хроматизм ξ₁…ξ₃ и второй порядок αc2 (ring_map, дифференциальная алгебра)."""
    try:
        rr = rmap.RingMap(ring, order=4)
        M = rr.one_turn()
        ch = rr.chromaticity(M)
        _, ac = rr.momentum_compaction(M)                     # (C, [αc, αc2, …])
        out = {}
        for pl in ("x", "y"):
            xi = ch["xi_" + pl]
            for k in range(1, min(4, len(xi))):
                out[f"xi{k}_{pl}"] = float(xi[k])
        out["alphac2"] = float(ac[1]) if len(ac) > 1 else float("nan")
        return out
    except Exception as e:
        return dict(error=str(e))


def radiation(ring, inf):
    """Интегралы излучения и всё равновесие; при ВЧ — синхротрон и энергетический аксептанс."""
    ring.rf = dict(volt_mv=inf["rf"][0], harmonic=inf["rf"][1])
    try:
        d = dict(rm.radiation_integrals(ring))
    finally:
        ring.rf = None
    volt, harm = inf["rf"]
    q = volt / d["U0"] if d["U0"] > 0 else 0.0                 # перенапряжение ВЧ
    if q > 1:                                                  # энергетический аксептанс ВЧ-сепаратрисы
        F = math.sqrt(q * q - 1) - math.acos(1.0 / q)
        d["delta_rf"] = math.sqrt(2 * d["U0"] / (math.pi * d["alphac"] * harm * inf["energy_mev"]) * F)
    d["overvoltage"] = q
    return d


def detuning(ring):
    """Сдвиг частот с действием ∂ν/∂J — нормальная форма отображения оборота."""
    try:
        nf = rmap.normal_form(ring, 4)
        return {k: float(v) for k, v in nf.items() if k.startswith("dnu")}
    except Exception as e:
        return dict(error=str(e))


def _grid_states(orbit, xs, ys):
    """Сетка начальных частиц на замкнутой орбите: смещения по x и y (м)."""
    out = []
    for y in ys:
        for x in xs:
            st = np.array(orbit, float); st[0] += x
            st[2] += (y if abs(y) > 1e-9 else 1e-6)
            out.append(st)
    return np.asarray(out)


def _mirror_y(ys, grid):
    """Дополнить верхнюю половину (y ≥ 0) зеркальной нижней; движение симметрично относит. медианы."""
    ys = np.asarray(ys); pos = np.where(ys > 1e-12)[0]
    ys_sym = np.concatenate([-ys[pos][::-1], ys])
    grid_sym = np.concatenate([grid[pos[::-1]], grid], axis=0)
    return ys_sym, grid_sym


# ---------- аналитика нелинейной динамики: tune footprint из нормальной формы (без трекинга) ----------
def _foot_amp(lin, dt, budget=0.4, cap_mm=20.0):
    """Амплитуды (ax, ay), при которых линейный по действию сдвиг частоты достигает `budget` — дальше
    линеаризация нормальной формы теряет смысл, footprint строится только в этом валидном диапазоне."""
    bx, by = lin["betax0"], lin["betay0"]
    ax = math.sqrt(2 * bx * budget / max(abs(dt["dnux/dJx"]), 1e-30))
    ay = math.sqrt(2 * by * budget / max(abs(dt["dnuy/dJy"]), 1e-30))
    return min(ax, cap_mm * 1e-3), min(ay, cap_mm * 1e-3)


def footprint(lin, dt, xe=None, ye=None, n=17):
    """Аналитический tune footprint: сетка амплитуд (ax, ay) → частоты по нормальной форме
    νx = νx0 + (∂νx/∂Jx)·Jx + (∂νx/∂Jy)·Jy,  J = a²/(2β). ЛИНЕАРИЗАЦИЯ по действию — верна на малых
    амплитудах; диапазон сам ограничивается бюджетом |Δν| ≈ 0.4 (или реальной апертурой из трекинга)."""
    if dt.get("error"):
        return None
    bx, by = lin["betax0"], lin["betay0"]
    axmax, aymax = _foot_amp(lin, dt)
    if xe:                                                    # если есть апертура из трекинга — не выходить за неё
        axmax = min(axmax, xe); aymax = min(aymax, ye)
    AX, AY = np.meshgrid(np.linspace(0.0, axmax, n), np.linspace(0.0, aymax, n))
    Jx = AX ** 2 / (2 * bx); Jy = AY ** 2 / (2 * by)
    nux = lin["nux"] + dt["dnux/dJx"] * Jx + dt["dnux/dJy"] * Jy
    nuy = lin["nuy"] + dt["dnuy/dJx"] * Jx + dt["dnuy/dJy"] * Jy
    return dict(ax=AX.ravel() * 1e3, ay=AY.ravel() * 1e3, nux=nux.ravel(), nuy=nuy.ravel(),
                amp=np.hypot(AX, AY).ravel() * 1e3, axmax=axmax * 1e3, aymax=aymax * 1e3,
                nux0=lin["nux"], nuy0=lin["nuy"])


# ---------- трекинг: автоподбор масштаба сканов ----------
def _edge(ring, co, axis, turns, limit, amp_signed, seed_axis=None, seed=0.0, n=48):
    """Граница устойчивости вдоль оси: наименьшая по модулю амплитуда гибели; если до предела пробы
    никто не гибнет — |amp_signed| (предел)."""
    amps = np.linspace(amp_signed / n, amp_signed, n)
    st = np.tile(co, (n, 1)).astype(float); st[:, axis] += amps
    if seed_axis is not None:
        st[:, seed_axis] += seed
    tbt = ring.track_tbt(st, turns, limit=limit)
    ok = ra.survived_turns(tbt, limit=limit) >= turns
    lost = np.where(~ok)[0]
    return abs(float(amps[lost[0]])) if lost.size else abs(float(amp_signed))


def auto_ranges(ring, inf, turns):
    """Автоподбор границ сканов пробегом по +x/−x, +y и ±δ. Возврат (co, x_edge, y_edge, d_pos, d_neg)."""
    co = ra.closed_orbit(ring, dims=4)
    Xp = inf["x_probe_mm"] * 1e-3; Yp = inf["y_probe_mm"] * 1e-3; Dp = inf["d_probe_pct"] / 100.0
    L = inf["loss_m"]
    xe = max(_edge(ring, co, 0, turns, L, +Xp), _edge(ring, co, 0, turns, L, -Xp))
    ye = _edge(ring, co, 2, turns, L, +Yp)
    seed = 0.3 * xe
    dpos = _edge(ring, co, 5, turns, L, +Dp, seed_axis=0, seed=seed)
    dneg = _edge(ring, co, 5, turns, L, -Dp, seed_axis=0, seed=seed)
    return co, xe, ye, dpos, dneg


def dynamic_aperture(ring, inf, turns, co, xe, ye):
    """Скан выживания в (x, y): диапазон 1.5× автограницы (полуторный запас), y симметризован."""
    xs = np.linspace(-1.5 * xe, 1.5 * xe, 41)
    ys = np.linspace(0.0, 1.5 * ye, 16)
    tbt = ring.track_tbt(_grid_states(co, xs, ys), turns, limit=inf["loss_m"])
    surv = ra.survived_turns(tbt, limit=inf["loss_m"]).reshape(len(ys), len(xs)).astype(float)
    ys, surv = _mirror_y(ys, surv)
    alive = surv >= turns; x0 = int(np.argmin(np.abs(xs)))
    ep = []; en = []
    for r in range(len(ys)):
        row = alive[r]; rp = xs[x0:][row[x0:]]; rn = xs[:x0 + 1][row[:x0 + 1]]
        ep.append(rp.max() if rp.size else 0.0); en.append(rn.min() if rn.size else 0.0)
    return dict(xs=xs, ys=ys, surv=surv, turns=turns, edge_pos=np.array(ep), edge_neg=np.array(en))


def freq_map(ring, inf, lin, turns, co, xe, ye):
    """Карта частот: ПЛОТНАЯ регулярная (x, y)-сетка до 1.5× апертуры, NAFF по двум половинам,
    диффузия D = log10|Δν|. Возвращает и 2D-массив diff2d (для сплошной тепловой карты), и плоские
    nux/nuy/diff (для точек на диаграмме частот). y симметризован."""
    n_half = turns // 2
    xs = np.linspace(-1.5 * xe, 1.5 * xe, 61); ys = np.linspace(0.0, 1.5 * ye, 16)
    tbt = ring.track_tbt(_grid_states(co, xs, ys), 2 * n_half, limit=inf["loss_m"])
    beta = (lin["betax0"], lin["betay0"]); alpha = (lin["alphax0"], lin["alphay0"])
    nux, nuy, diff = ra.frequency_map(tbt, beta, alpha, co, n_half)
    sh = (len(ys), len(xs))
    ys_sym, nux = _mirror_y(ys, nux.reshape(sh))
    _, nuy = _mirror_y(ys, nuy.reshape(sh))
    _, diff = _mirror_y(ys, diff.reshape(sh))
    xg = np.tile(xs, len(ys_sym)); yg = np.repeat(ys_sym, len(xs))
    return dict(xs=xs, ys=ys_sym, diff2d=diff, x=xg, y=yg, nux=nux.ravel(), nuy=nuy.ravel(),
                diff=diff.ravel(), nux0=lin["nux"], nuy0=lin["nuy"])


def momentum_aperture(ring, inf, turns, co, xe, dpos, dneg):
    """Импульсная апертура: скан по δ с бетатронной затравкой, диапазон 1.2× автограницы."""
    D = 1.2 * max(dpos, dneg, 1e-3)
    ds = np.linspace(-D, D, 61)
    st = np.tile(co, (len(ds), 1)).astype(float); st[:, 5] += ds; st[:, 0] += 0.3 * xe
    tbt = ring.track_tbt(st, turns, limit=inf["loss_m"])
    surv = ra.survived_turns(tbt, limit=inf["loss_m"]); alive = surv >= turns
    pos = ds[(ds > 0) & alive]; neg = ds[(ds < 0) & alive]
    return dict(ds=ds, surv=surv, turns=turns,
                acc_pos=pos.max() if pos.size else 0.0, acc_neg=neg.min() if neg.size else 0.0)


# ------------------------------------------------------------------ АНСАМБЛЬ МАШИН С ОШИБКАМИ
def _edge_prob(xs, ys, prob, level=0.5):
    """Граница ДА по уровню вероятности выживания level: дальний по |x| в каждой строке y."""
    x0 = int(np.argmin(np.abs(xs))); ep = []; en = []
    for r in range(len(ys)):
        row = prob[r] >= level
        rp = xs[x0:][row[x0:]]; rn = xs[:x0 + 1][row[:x0 + 1]]
        ep.append(rp.max() if rp.size else 0.0); en.append(rn.min() if rn.size else 0.0)
    return np.array(ep), np.array(en)


def ensemble(lattice, inf, lin, turns, n_seeds, tol, cutoff, seed0=1000, correct=False, nsl=40):
    """Ансамбль из n_seeds случайных машин с ошибками. Возврат агрегатов: ДА как ВЕРОЯТНОСТЬ выживания
    по seed'ам, средняя диффузия, импульсная как доля выживших, и статистика влияния ошибок (искажение
    орбиты, разброс частот, доля неустойчивых машин). Сетки фиксируются по ИДЕАЛЬНОЙ машине (1.5× границ),
    чтобы все seed'ы считались на одной сетке.
    correct=True — на каждой машине выправить замкнутую орбиту корректорами (матрица отклика + SVD,
    orbit_correct.py); отчёт тогда отражает УЖЕ скорректированную установку и добавляет статистику
    коррекции (орбита до/после, СКО и максимум углов корректоров)."""
    import shutil
    import tempfile
    ideal = _ring(inf, nsl=nsl)
    co0, xe, ye, dpos, dneg = auto_ranges(ideal, inf, turns)
    xs = np.linspace(-1.5 * xe, 1.5 * xe, 41); ys_h = np.linspace(0.0, 1.5 * ye, 16)
    fxs = np.linspace(-1.5 * xe, 1.5 * xe, 61); fys_h = np.linspace(0.0, 1.5 * ye, 16)
    D = 1.2 * max(dpos, dneg, 1e-3); ds = np.linspace(-D, D, 61)
    n_half = turns // 2
    beta = (lin["betax0"], lin["betay0"]); alpha = (lin["alphax0"], lin["alphay0"])
    resp = None; n_corr = 0
    if correct:                                               # матрица отклика — один раз, на идеале (кэш)
        import orbit_correct as oc
        resp = oc.cached_response_matrix(ideal, lattice)
        n_corr = len(resp["corr"])
    bpm_rows = np.arange(1, ideal.n_mag + 1)                   # выходы всех магнитов — точки наблюдения орбиты
    da_sum = np.zeros((len(ys_h), len(xs)))
    fdiff_sum = np.zeros((len(fys_h), len(fxs))); fdiff_cnt = np.zeros_like(fdiff_sum)
    mom_sum = np.zeros(len(ds))
    co_x = []; co_y = []; nux = []; nuy = []; n_unstable = 0
    co_x0 = []; co_y0 = []; kx = []; ky = []; kx_max = 0.0; ky_max = 0.0; n_threaded = 0   # статистика коррекции
    root = tempfile.mkdtemp(prefix="ens_")
    try:
        for s in range(n_seeds):
            wsdir = os.path.join(root, f"s{s}")
            errmod.perturbed_workspace(lattice, seed0 + s, tol, wsdir, cutoff)
            r = _ring(inf, ws=wsdir, nsl=nsl)
            if correct:                                       # выправить орбиту корректорами (стационарно)
                import orbit_correct as oc
                try:
                    res = oc.correct(r, resp, aperture=inf["loss_m"])
                except Exception:
                    res = None
                if res is None:
                    n_unstable += 1; r.close(); continue
                co = res["co"]
                if res.get("threaded"):
                    n_threaded += 1
                co_x0.append(oc.rms(res["orb_before"][0])); co_y0.append(oc.rms(res["orb_before"][1]))
                co_x.append(oc.rms(res["orb_after"][0])); co_y.append(oc.rms(res["orb_after"][1]))
                kx.append(oc.rms(res["thx"])); ky.append(oc.rms(res["thy"]))
                kx_max = max(kx_max, float(np.max(np.abs(res["thx"]))) if len(res["thx"]) else 0.0)
                ky_max = max(ky_max, float(np.max(np.abs(res["thy"]))) if len(res["thy"]) else 0.0)
            else:
                try:
                    co = ra.closed_orbit(r, dims=4)
                    stable = bool(np.all(np.isfinite(co)))
                except Exception:
                    co, stable = None, False
                if not stable:                                # машина неустойчива → всюду 0 выживших
                    n_unstable += 1; r.close(); continue
                orb = r.track_elements([co])[0]               # орбита по всему кольцу (СКО искажения)
                co_x.append(float(np.sqrt(np.mean(np.square(orb[bpm_rows, 0])))))
                co_y.append(float(np.sqrt(np.mean(np.square(orb[bpm_rows, 2])))))
            try:                                              # частоты этой машины (для разброса тюна)
                _, M = ra.one_turn_matrix(r, orbit=co)
                qx, qy = ra.tunes_from_matrix(M)
                if math.isfinite(qx):
                    nux.append(qx); nuy.append(qy)
            except Exception:
                pass
            da_sum += _da_alive(r, inf, turns, co, xs, ys_h).astype(float)
            st = _grid_states(co, fxs, fys_h)
            tbt = r.track_tbt(st, 2 * n_half, limit=inf["loss_m"])
            _, _, diff = ra.frequency_map(tbt, beta, alpha, co, n_half)
            diff = diff.reshape(len(fys_h), len(fxs)); ok = np.isfinite(diff)
            fdiff_sum[ok] += diff[ok]; fdiff_cnt[ok] += 1
            mom_sum += _mom_alive(r, inf, turns, co, xe, ds).astype(float)
            r.close()
    finally:
        ideal.close()
        shutil.rmtree(root, ignore_errors=True)
    prob = da_sum / max(n_seeds, 1)
    ys, prob = _mirror_y(ys_h, prob)
    ep, en = _edge_prob(xs, ys, prob, 0.5)
    fdiff = np.where(fdiff_cnt > 0, fdiff_sum / np.maximum(fdiff_cnt, 1), np.nan)
    fys, fdiff = _mirror_y(fys_h, fdiff)
    da = dict(xs=xs, ys=ys, prob=prob, turns=turns, edge_pos=ep, edge_neg=en, ensemble=n_seeds)
    fma = dict(xs=fxs, ys=fys, diff2d=fdiff, ensemble=n_seeds)
    momp = mom_sum / max(n_seeds - n_unstable, 1)
    mom = dict(ds=ds, prob=momp, turns=turns, ensemble=n_seeds,
               acc_pos=float(ds[(ds > 0) & (momp >= 0.5)].max()) if np.any((ds > 0) & (momp >= 0.5)) else 0.0,
               acc_neg=float(ds[(ds < 0) & (momp >= 0.5)].min()) if np.any((ds < 0) & (momp >= 0.5)) else 0.0)
    def _rms(a):                                              # СКО по машинам, nan-безопасно (одна машина
        a = np.asarray(a, float)                              # с потерянной частицей не должна портить агрегат)
        a = a[np.isfinite(a)]
        return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")
    err = dict(n_seeds=n_seeds, n_unstable=n_unstable, n_stable=n_seeds - n_unstable,
               co_x_rms=_rms(co_x), co_y_rms=_rms(co_y),
               nux_mean=float(np.mean(nux)) if nux else float("nan"),
               nux_std=float(np.std(nux)) if nux else float("nan"),
               nuy_mean=float(np.mean(nuy)) if nuy else float("nan"),
               nuy_std=float(np.std(nuy)) if nuy else float("nan"),
               da_ideal_x=1.5 * xe, corrected=correct, n_corr=n_corr)
    if correct:                                               # орбита ДО коррекции и силы корректоров
        err.update(co_x_rms0=_rms(co_x0), co_y_rms0=_rms(co_y0),
                   corr_x_rms=_rms(kx), corr_y_rms=_rms(ky), corr_x_max=kx_max, corr_y_max=ky_max,
                   n_threaded=n_threaded)
    return dict(da=da, fma=fma, mom=mom, err=err)


def _da_alive(ring, inf, turns, co, xs, ys):
    tbt = ring.track_tbt(_grid_states(co, xs, ys), turns, limit=inf["loss_m"])
    return ra.survived_turns(tbt, limit=inf["loss_m"]).reshape(len(ys), len(xs)) >= turns


def _mom_alive(ring, inf, turns, co, xe, ds):
    st = np.tile(co, (len(ds), 1)).astype(float); st[:, 5] += ds; st[:, 0] += 0.3 * xe
    return ra.survived_turns(ring.track_tbt(st, turns, limit=inf["loss_m"]), limit=inf["loss_m"]) >= turns


def families(ring):
    """Таблица магнитов по семействам: тип, число, длина, сила (K1/K2/K3)."""
    fam = {}
    for m in ring.mags:
        base = _family(m.dev)
        r = fam.setdefault(base, dict(kind=m.kind, n=0, L=float(getattr(m, "Leff", 0.0) or 0.0),
                                      k1=0.0, k2=0.0, k3=0.0, theta=0.0))
        r["n"] += 1
        for a, src in (("k1", "k1"), ("k2", "k2"), ("k3", "k3"), ("theta", "theta0")):
            v = float(getattr(m, src, 0.0) or 0.0)
            if abs(v) > abs(r[a]):                            # характерное значение семейства — по модулю
                r[a] = v
    return fam


import collections


def _elem_sig(m):
    """Сигнатура элемента для распознавания структуры: тип, длина и силы (округлённые)."""
    return (m.kind, round(float(getattr(m, "Leff", 0.0) or 0.0), 3),
            round(float(getattr(m, "k1", 0.0) or 0.0), 4), round(float(getattr(m, "k2", 0.0) or 0.0), 3),
            round(float(getattr(m, "k3", 0.0) or 0.0), 1), round(float(getattr(m, "theta0", 0.0) or 0.0), 5))


def _cyc(seq, i0, i1):
    """Циклический срез [i0, i1) — с переходом через шов кольца, если i1 <= i0."""
    return seq[i0:i1] if i0 < i1 else (list(seq[i0:]) + list(seq[:i1]))


def detect_structures(ring, lin, n_sp):
    """Разбить кольцо на суперпериоды и сгруппировать уникальные структуры. Возврат списка
    {label, count, s, betax, betay, etax, length} — по одной записи на уникальную структуру
    (стандартный суперпериод + область(и) инжекции); None, если регулярная структура не распознаётся."""
    mags = ring.mags
    fam = collections.Counter(_family(m.dev) for m in mags)
    cand = [f for f, c in fam.items() if c == n_sp]                # семейство раз в суперпериод — якорь
    if not cand:
        return None
    def gaps(f):
        idx = [i for i, m in enumerate(mags) if _family(m.dev) == f]
        return idx, np.diff(idx)
    anchor = min(cand, key=lambda f: (len(set(gaps(f)[1].tolist())), f))   # самый ровный по шагу
    idx, _ = gaps(anchor)
    nc = len(idx)
    cells = [(idx[k], idx[(k + 1) % nc]) for k in range(nc)]       # nc суперпериодов (последний — через шов)
    sig = [tuple(_elem_sig(m) for m in _cyc(mags, i0, i1)) for (i0, i1) in cells]
    modal = collections.Counter(sig).most_common(1)[0][0]
    # повернуть порядок так, чтобы начать со стандартной ячейки — тогда особые идут подряд
    st = sig.index(modal)
    rot = list(range(st, nc)) + list(range(0, st))
    special_rot = [i for i, k in enumerate(rot) if sig[k] != modal]
    s_arr = np.asarray(lin["s"]); C = float(s_arr[-1])

    def extract(i0, i1):                                          # твисс на циклическом диапазоне [i0, i1]
        if i0 < i1:
            sl = slice(i0, i1 + 1)
            return (s_arr[sl] - s_arr[i0], lin["betax"][sl], lin["betay"][sl], lin["etax"][sl])
        s1 = s_arr[i0:] - s_arr[i0]; s2 = (C - s_arr[i0]) + s_arr[:i1 + 1]
        cat = lambda a, b: np.concatenate([a, b])
        return (cat(s1, s2), cat(lin["betax"][i0:], lin["betax"][:i1 + 1]),
                cat(lin["betay"][i0:], lin["betay"][:i1 + 1]), cat(lin["etax"][i0:], lin["etax"][:i1 + 1]))

    out = []
    # стандартный суперпериод — представитель (первая стандартная ячейка)
    i0, i1 = cells[st]
    sx, bx, by, ex = extract(i0, i1)
    out.append(dict(label="Стандартный суперпериод", count=sig.count(modal),
                    s=sx, betax=bx, betay=by, etax=ex, length=float(sx[-1])))
    # особые области — непрерывные (в повёрнутом порядке) прогоны
    if special_rot:
        runs = []; run = [special_rot[0]]
        for a in special_rot[1:]:
            if a == run[-1] + 1:
                run.append(a)
            else:
                runs.append(run); run = [a]
        runs.append(run)
        for r in runs:
            c0 = cells[rot[r[0]]][0]; c1 = cells[rot[r[-1]]][1]
            sx, bx, by, ex = extract(c0, c1)
            out.append(dict(label="Область инжекции" if len(runs) == 1 else "Особый суперпериод",
                            count=len(r), s=sx, betax=bx, betay=by, etax=ex, length=float(sx[-1])))
    return dict(anchor=anchor, structures=out)


def _family(dev):
    """Имя семейства: срезаем префикс sila26_/s2_ и хвостовой индекс среза (_1.._5, _1A…)."""
    base = dev.split("_", 1)[1] if dev.split("_", 1)[0].endswith(("26", "25", "s2", "2")) and "_" in dev else dev
    import re
    return re.sub(r"_(\d+[A-Z]?|[A-Z])$", "", base)


def _lin(t):
    return (t[0], t[1], int(t[2]))


def analyze(lattice, turns=512, fast=False, track_cache=None, errors_cfg=None, nsl=40):
    """Полный расчёт установки. Аналитика (геометрия/оптика/структура/footprint) — всегда мгновенно.
    Трекинг (ДА/карта частот/импульс) считается, если не fast; track_cache подставляет готовый трекинг
    (--rebuild). errors_cfg={n_seeds, tol, cutoff} — ансамбль машин с ошибками вместо идеального трекинга.
    nsl — число срезов магнита в трекинге (40 для «геройских» паспортов; для ансамбля хватает ~12,
    граница ДА не сдвигается, но в ~3× быстрее — см. бенчмарк nsl)."""
    inf = dict(info_of(lattice)); inf["_id"] = lattice
    ring = _ring(inf, nsl=nsl)
    out = dict(lattice=lattice, label=inf["label"], turns=turns,
               fast=fast and track_cache is None and errors_cfg is None)
    out["geom"] = geometry(ring, inf)
    out["lin"] = linear_optics(ring, inf)
    out["struct"] = detect_structures(ring, out["lin"], inf["superperiods"])   # уникальные суперпериоды
    out["chroma"] = chroma(ring)
    out["rad"] = radiation(ring, inf)
    out["detune"] = detuning(ring)
    out["fam"] = families(ring)
    # аналитический tune footprint из нормальной формы — мгновенно, без трекинга (диапазон по бюджету |Δν|)
    out["foot"] = footprint(out["lin"], out["detune"])
    if errors_cfg is not None:                                # ансамбль машин с ошибками (± коррекция орбиты)
        out.update(ensemble(lattice, inf, out["lin"], turns, errors_cfg["n_seeds"],
                            errors_cfg["tol"], errors_cfg["cutoff"],
                            correct=errors_cfg.get("correct", False), nsl=nsl))
        out["errors"] = errors_cfg
    elif track_cache is not None:                             # готовый трекинг из кэша (--rebuild)
        out.update(track_cache); out["turns"] = int(track_cache.get("turns", turns))
    elif not fast:
        co, xe, ye, dpos, dneg = auto_ranges(ring, inf, turns)
        out["ranges"] = dict(xe=xe, ye=ye, dpos=dpos, dneg=dneg)
        out["da"] = dynamic_aperture(ring, inf, turns, co, xe, ye)
        out["fma"] = freq_map(ring, inf, out["lin"], turns, co, xe, ye)
        out["mom"] = momentum_aperture(ring, inf, turns, co, xe, dpos, dneg)
    return out


# ------------------------------------------------------------------ ГРАФИКИ
def _fig_twiss(lin, sp_len):
    fig = go.Figure()
    s = lin["s"]
    fig.add_trace(go.Scattergl(x=s, y=lin["betax"], mode="lines", name="βx", line=dict(color="#2b7bba", width=1.6)))
    fig.add_trace(go.Scattergl(x=s, y=lin["betay"], mode="lines", name="βy", line=dict(color="#e15759", width=1.6)))
    fig.add_trace(go.Scattergl(x=s, y=lin["etax"] * 100, mode="lines", name="ηx × 100",
                               line=dict(color="#59a14f", width=1.4, dash="dot"), yaxis="y"))
    fig.update_layout(title=dict(text="Функции Твисса вдоль кольца (βx, βy, ηx)", x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=62, r=18, t=34, b=42), plot_bgcolor="white", hovermode="x unified",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="s, м", range=[0, min(s[-1], sp_len * 1.02)])
    fig.update_yaxes(title_text="β, м  ·  ηx × 100, м")
    return ta._grid(fig)


def _fig_twiss_struct(st):
    """Функции Твисса на ОДНОЙ уникальной структуре (стандартный суперпериод или область инжекции)."""
    s = st["s"]
    fig = go.Figure()
    fig.add_trace(go.Scattergl(x=s, y=st["betax"], mode="lines", name="βx", line=dict(color="#2b7bba", width=1.7)))
    fig.add_trace(go.Scattergl(x=s, y=st["betay"], mode="lines", name="βy", line=dict(color="#e15759", width=1.7)))
    fig.add_trace(go.Scattergl(x=s, y=st["etax"] * 100, mode="lines", name="ηx × 100",
                               line=dict(color="#59a14f", width=1.5, dash="dot")))
    fig.update_layout(title=dict(text=f"Функции Твисса — {st['label'].lower()} (×{st['count']}, {st['length']:.2f} м)",
                                 x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=62, r=18, t=40, b=70), plot_bgcolor="white", hovermode="x unified",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="s, м (от начала структуры)", range=[0, float(s[-1])])
    fig.update_yaxes(title_text="β, м  ·  ηx × 100, м")
    return ta._grid(fig)


def _fig_tunes(lin, foot=None, fma=None):
    """Диаграмма рабочей точки: резонансы + рабочая точка + аналитический tune footprint (∂ν/∂J);
    при наличии трекинга поверх — точки карты частот, окрашенные диффузией."""
    fx, fy = lin["nux"] % 1, lin["nuy"] % 1
    fig = go.Figure()
    for m, n, p, seg in ra.resonance_lines(5, (0.0, 1.0, 0.0, 1.0)):
        order = abs(m) + abs(n)
        (x0, y0), (x1, y1) = seg
        fig.add_trace(go.Scatter(x=[x0, x1], y=[y0, y1], mode="lines", hoverinfo="skip", showlegend=False,
                                 line=dict(color="rgba(120,120,120,%.2f)" % (0.55 - 0.07 * order),
                                           width=max(0.4, 2.2 - 0.35 * order))))
    if foot is not None:                                       # аналитический footprint: цвет = амплитуда
        fig.add_trace(go.Scattergl(x=foot["nux"] % 1, y=foot["nuy"] % 1, mode="markers",
                                   marker=dict(size=5, color=foot["amp"], colorscale="Viridis", showscale=True,
                                               colorbar=dict(title="ампл., мм", thickness=12, len=0.9)),
                                   name="footprint (аналитика)", hoverinfo="skip"))
    if fma is not None and "nux" in fma:                       # трекинг: точки, окрашенные диффузией
        ok = np.isfinite(fma["nux"])
        fig.add_trace(go.Scattergl(x=fma["nux"][ok] % 1, y=fma["nuy"][ok] % 1, mode="markers",
                                   marker=dict(size=3, color="rgba(20,20,20,0.35)"),
                                   name="карта частот (трекинг)", hoverinfo="skip"))
    fig.add_trace(go.Scatter(x=[fx], y=[fy], mode="markers",
                             marker=dict(size=13, color="#e15759", symbol="star"), name="рабочая точка"))
    fig.update_layout(title=dict(text="Диаграмма рабочей точки: резонансы (до 5-го) + tune footprint",
                                 x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=54, r=18, t=40, b=70), plot_bgcolor="white",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="νx (дробная)", range=[0, 1], constrain="domain")
    fig.update_yaxes(title_text="νy (дробная)", range=[0, 1], scaleanchor="x", scaleratio=1)
    return ta._grid(fig, legend=False)


def _fig_da_estimate(est, lin):
    """Аналитическая оценка ДА: радиус устойчивости по направлениям (footprint до резонанса)."""
    ang = np.concatenate([est["ang"], math.pi - est["ang"][::-1]])   # отразить в первый+второй квадрант
    r = np.concatenate([est["r"], est["r"][::-1]])
    x = r * np.cos(ang); y = r * np.sin(ang)
    x = np.concatenate([x, x[::-1]]); y = np.concatenate([y, -y[::-1]])  # и нижняя половина
    fig = go.Figure(go.Scatter(x=x, y=y, mode="lines", fill="toself",
                               line=dict(color="#59a14f", width=2), fillcolor="rgba(89,161,79,0.15)",
                               name="оценка ДА"))
    fig.update_layout(title=dict(text=f"Аналитическая оценка ДА (footprint до резонанса ≤{est['order_max']}-го), "
                                      f"мин. {est['r_min']:.2f} мм", x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=54, r=18, t=34, b=42), plot_bgcolor="white")
    fig.update_xaxes(title_text="x, мм")
    fig.update_yaxes(title_text="y, мм", scaleanchor="x", scaleratio=1)
    return ta._grid(fig, legend=False)


def _fig_da(da):
    """Карта выживания в (x, y): число оборотов (идеал) или ВЕРОЯТНОСТЬ выживания по ансамблю + граница."""
    ens = da.get("ensemble")
    if ens:                                                   # ансамбль: доля выживших машин 0..1
        z = da["prob"]; cb = dict(title="доля машин", thickness=12, len=0.9)
        hm = go.Heatmap(x=da["xs"] * 1e3, y=da["ys"] * 1e3, z=z, colorscale="Viridis",
                        zmin=0, zmax=1, zsmooth="best", colorbar=cb)
        ttl = f"Динамическая апертура: доля выживших из {ens} машин с ошибками ({da['turns']} об.)"
    else:
        hm = go.Heatmap(x=da["xs"] * 1e3, y=da["ys"] * 1e3, z=da["surv"],
                        colorscale="Viridis", colorbar=dict(title="оборотов", thickness=12, len=0.9))
        ttl = f"Динамическая апертура: выживание за {da['turns']} оборотов"
    fig = go.Figure(hm)
    lab = "граница 50% машин" if ens else "граница +x"
    fig.add_trace(go.Scatter(x=da["edge_pos"] * 1e3, y=da["ys"] * 1e3, mode="lines+markers",
                             line=dict(color="#e15759", width=2), name=lab, marker=dict(size=4)))
    fig.add_trace(go.Scatter(x=da["edge_neg"] * 1e3, y=da["ys"] * 1e3, mode="lines+markers",
                             line=dict(color="#e15759", width=2, dash="dot"), name="", marker=dict(size=4)))
    fig.update_layout(title=dict(text=ttl,
                                 x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=54, r=18, t=40, b=70), plot_bgcolor="white",
                      legend=dict(orientation="h", y=-0.16, x=0.5, xanchor="center", yanchor="top"))
    fig.update_xaxes(title_text="x, мм")
    fig.update_yaxes(title_text="y, мм")
    return ta._grid(fig, legend=False)


def _fig_fma_space(fma):
    """Диффузия в конфигурационном пространстве (x, y) — сплошная сглаженная тепловая карта; где
    движение регулярно (синий), где хаотично/потеряно (красный/пусто)."""
    fig = go.Figure(go.Heatmap(x=fma["xs"] * 1e3, y=fma["ys"] * 1e3, z=fma["diff2d"],
                               colorscale="Turbo", zmin=-10, zmax=-2, zsmooth="best", connectgaps=False,
                               colorbar=dict(title="log₁₀ D", thickness=12, len=0.9)))
    fig.update_layout(title=dict(text="Карта диффузии в плоскости (x, y): регулярность движения",
                                 x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=54, r=18, t=34, b=42), plot_bgcolor="#f2f2f2")
    fig.update_xaxes(title_text="x, мм")
    fig.update_yaxes(title_text="y, мм")
    return ta._grid(fig, legend=False)


def _fig_momentum(mom, delta_rf=None):
    ens = mom.get("ensemble")
    if ens:                                                   # ансамбль: доля выживших машин 0..1
        y = mom["prob"]; ytitle = "доля выживших машин"; ttl = f"Импульсная апертура: доля из {ens} машин с ошибками"
    else:
        y = np.minimum(mom["surv"], mom["turns"]); ytitle = "оборотов до потери"
        ttl = f"Импульсная апертура: выживание за {mom['turns']} оборотов (красное — трекинг, зелёное — ВЧ)"
    fig = go.Figure(go.Scatter(x=mom["ds"] * 100, y=y, mode="lines+markers",
                               line=dict(color="#2b7bba", width=2), marker=dict(size=5), name="выживание"))
    fig.add_vline(x=mom["acc_pos"] * 100, line=dict(color="#e15759", dash="dash"), annotation_text="ДА +δ")
    fig.add_vline(x=mom["acc_neg"] * 100, line=dict(color="#e15759", dash="dash"), annotation_text="ДА −δ")
    if delta_rf:
        fig.add_vline(x=delta_rf * 100, line=dict(color="#59a14f", dash="dot"), annotation_text="ВЧ")
        fig.add_vline(x=-delta_rf * 100, line=dict(color="#59a14f", dash="dot"))
    fig.update_layout(title=dict(text=ttl, x=0.5, y=0.97, font=dict(size=14)),
                      margin=dict(l=62, r=18, t=34, b=42), plot_bgcolor="white")
    fig.update_xaxes(title_text="δ = Δp/p, %")
    fig.update_yaxes(title_text=ytitle, rangemode="tozero")
    return ta._grid(fig, legend=False)


def _tbl_errors(err, tol):
    """Сводка влияния ошибок: заданные допуски, устойчивость ансамбля, искажение орбиты, разброс частот,
    и (если была коррекция) орбита до/после + силы корректоров."""
    corr = err.get("corrected")
    rows = [
        ("Машин в ансамбле", str(err["n_seeds"])),
        ("Неустойчивых машин", f"{err['n_unstable']} из {err['n_seeds']}"),
    ]
    if corr:
        rows += [
            ("Корректоров орбиты (гор.+верт.)", str(err["n_corr"])),
            ("Машин через протяжку (threading)", str(err.get("n_threaded", 0))),
            ("СКО орбиты x: до → после коррекции",
             f"{_num(err['co_x_rms0'] * 1e3, 3)} → {_num(err['co_x_rms'] * 1e3, 3)} мм"),
            ("СКО орбиты y: до → после коррекции",
             f"{_num(err['co_y_rms0'] * 1e3, 3)} → {_num(err['co_y_rms'] * 1e3, 3)} мм"),
            ("Углы корректоров x: СКО / макс.",
             f"{_num(err['corr_x_rms'] * 1e6, 3)} / {_num(err['corr_x_max'] * 1e6, 3)} мкрад"),
            ("Углы корректоров y: СКО / макс.",
             f"{_num(err['corr_y_rms'] * 1e6, 3)} / {_num(err['corr_y_max'] * 1e6, 3)} мкрад"),
        ]
    else:
        rows += [
            ("СКО искажения орбиты x (по кольцу)", f"{_num(err['co_x_rms'] * 1e3, 3)} мм"),
            ("СКО искажения орбиты y (по кольцу)", f"{_num(err['co_y_rms'] * 1e3, 3)} мм"),
        ]
    rows += [
        ("νx: среднее ± разброс", f"{_num(err['nux_mean'], 6)} ± {_num(err['nux_std'], 2)}"),
        ("νy: среднее ± разброс", f"{_num(err['nuy_mean'], 6)} ± {_num(err['nuy_std'], 2)}"),
    ]
    for k, t in tol.items():
        if any(t.get(c) for c in errmod.TOL_COLS):
            rows.append((f"Допуск {KIND_RU.get(k, k)} (СКО)",
                         f"сдвиг {t.get('dx', 0)*1e6:.0f} мкм · крен {t.get('droll', 0)*1e3:.2f} мрад · "
                         f"поле {t.get('field_pct', 0):.2f} %"))
    head = "Ошибки установки (ансамбль" + (", после коррекции орбиты" if corr else ", без коррекции орбиты") + ")"
    return _table([head, "Значение"], rows, [55, 45])


def _table(header, rows, widths=None):
    fig = go.Figure(go.Table(columnwidth=widths,
        header=dict(values=header, fill_color="#4e79a7", font=dict(color="white", size=12), align="left", height=26),
        cells=dict(values=list(map(list, zip(*rows))) if rows else [[] for _ in header],
                   align="left", height=22, font=dict(size=12),
                   fill_color=[["#f4f6f8" if i % 2 else "white" for i in range(len(rows))]])))
    fig.update_layout(margin=dict(l=6, r=6, t=6, b=6))
    return fig


def _num(v, d=4):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "—"
    return f"{v:.{d}g}".replace(".", ",")


# ------------------------------------------------------------------ СБОРКА ТАБЛИЦ
def _tbl_geometry(g):
    rows = [
        ("Энергия пучка", f"{_num(g['energy_mev']/1000, 4)} ГэВ"),
        ("Релятивистский фактор γ", _num(g["gamma"], 6)),
        ("Периметр", f"{_num(g['circumference'], 9)} м"),
        ("Магнитная жёсткость Bρ", f"{_num(g['brho'], 6)} Тл·м"),
        ("Число суперпериодов", str(g["superperiods"])),
        ("Период обращения T₀", f"{_num(g['T0']*1e9, 5)} нс"),
        ("Частота обращения f₀", f"{_num(g['f0']/1e3, 5)} кГц"),
        ("ВЧ: напряжение", f"{_num(g['rf_volt_mv'], 4)} МВ"),
        ("ВЧ: кратность", str(g["rf_harmonic"])),
        ("ВЧ: частота", f"{_num(g['rf_freq']/1e6, 5)} МГц"),
        ("Суммарный угол изгиба", f"{_num(g['sum_angle_2pi'], 6)} × 2π"),
    ]
    for k, n in sorted(g["kinds"].items()):
        rows.append((f"Магниты: {KIND_RU.get(k, k)}", str(n)))
    return _table(["Параметр", "Значение"], rows, [50, 50])


def _tbl_linear(lin):
    rows = [
        ("Бетатронные частоты νx / νy", f"{_num(lin['nux'], 7)} / {_num(lin['nuy'], 7)}"),
        ("Набег фазы на суперпериод μx / μy", f"{_num(lin['mux_cell'], 5)} / {_num(lin['muy_cell'], 5)}"),
        ("βx / βy в точке наблюдения", f"{_num(lin['betax0'], 5)} / {_num(lin['betay0'], 5)} м"),
        ("αx / αy в точке наблюдения", f"{_num(lin['alphax0'], 4)} / {_num(lin['alphay0'], 4)}"),
        ("max βx / max βy", f"{_num(lin['betax_max'], 5)} / {_num(lin['betay_max'], 5)} м"),
        ("ηx в точке наблюдения / max", f"{_num(lin['etax0']*1e3, 4)} / {_num(lin['etax_max']*1e3, 4)} мм"),
        ("Коэффициент уплотнения орбит αc", _num(lin["alphac"], 6)),
        ("Замкнутая орбита x", f"{_num(lin['closed_orbit'][0]*1e6, 3)} мкм"),
    ]
    return _table(["Линейная оптика", "Значение"], rows, [55, 45])


def _tbl_chroma(ch):
    if ch.get("error"):
        return _table(["Хроматичность", "Значение"], [("ошибка расчёта", ch["error"])], [40, 60])
    rows = [
        ("Натуральная ξ₁ (x / y)", f"{_num(ch.get('xi1_x'), 5)} / {_num(ch.get('xi1_y'), 5)}"),
        ("2-й порядок ξ₂ (x / y)", f"{_num(ch.get('xi2_x'), 5)} / {_num(ch.get('xi2_y'), 5)}"),
        ("3-й порядок ξ₃ (x / y)", f"{_num(ch.get('xi3_x'), 5)} / {_num(ch.get('xi3_y'), 5)}"),
        ("αc второго порядка αc2", _num(ch.get("alphac2"), 5)),
    ]
    return _table(["Хроматичность (ν(δ))", "Значение"], rows, [55, 45])


def _tbl_radiation(r):
    rows = [
        ("Интеграл I₁", f"{_num(r['I1'], 5)} м"),
        ("Интеграл I₂", f"{_num(r['I2'], 5)} м⁻¹"),
        ("Интеграл I₃", f"{_num(r['I3'], 5)} м⁻²"),
        ("Интеграл I₄", f"{_num(r['I4'], 5)} м⁻¹"),
        ("Интеграл I₅", f"{_num(r['I5'], 5)} м⁻¹"),
        ("Потеря энергии за оборот U₀", f"{_num(r['U0'], 5)} МэВ"),
        ("Времена затухания τx / τy / τs", f"{_num(r['taux']*1e3, 4)} / {_num(r['tauy']*1e3, 4)} / {_num(r['taudelta']*1e3, 4)} мс"),
        ("Числа Робинсона Jx / Jy / Jδ", f"{_num(r['Jx'], 4)} / {_num(r['Jy'], 4)} / {_num(r['Jdelta'], 4)}"),
        ("Равновесный эмиттанс εx", f"{_num(r['ex0']*1e12, 4)} пм·рад"),
        ("Энергетический разброс σδ", _num(r["Sdelta0"], 4)),
    ]
    if "Qs" in r:
        rows += [("Синхротронная частота Qs", _num(r["Qs"], 5)),
                 ("Синхронная фаза φs", f"{_num(r['phi_s'], 5)}°"),
                 ("Длина сгустка σz", f"{_num(r.get('Sz0', float('nan'))*1e3, 4)} мм"),
                 ("Перенапряжение ВЧ eV/U₀", _num(r.get("overvoltage"), 4))]
    if "delta_rf" in r:
        rows.append(("Аксептанс ВЧ δ_rf", f"±{_num(r['delta_rf']*100, 4)} %"))
    return _table(["Излучение и равновесие", "Значение"], rows, [55, 45])


def _tbl_detune(dt, mom=None):
    if dt.get("error"):
        rows = [("ошибка расчёта нормальной формы", dt["error"])]
    else:
        rows = [(k.replace("dnu", "∂ν").replace("/dJ", "/∂J"), _num(v, 5) + " м⁻¹·рад⁻¹")
                for k, v in dt.items() if not k.startswith("|re|")]
    if mom is not None:
        rows.append(("Импульсная апертура (трекинг)",
                     f"+{_num(mom['acc_pos']*100, 3)} % / {_num(mom['acc_neg']*100, 3)} %"))
    return _table(["Нелинейная динамика", "Значение"], rows, [55, 45])


def _tbl_families(fam):
    """Таблица семейств: в графе «Сила» перечислены ВСЕ ненулевые компоненты магнита — угол изгиба θ
    и/или нормированные градиенты K1/K2/K3. У комбинированных (диполь-квадруполь) их несколько."""
    rows = []
    for name in sorted(fam, key=lambda k: (fam[k]["kind"], k)):
        r = fam[name]
        parts = []
        if abs(r["theta"]) > 1e-12:
            parts.append(f"θ = {_num(r['theta'], 5)} рад")
        if abs(r["k1"]) > 1e-12:
            parts.append(f"K1 = {_num(r['k1'], 5)} м⁻²")
        if abs(r["k2"]) > 1e-12:
            parts.append(f"K2 = {_num(r['k2'], 5)} м⁻³")
        if abs(r["k3"]) > 1e-12:
            parts.append(f"K3 = {_num(r['k3'], 5)} м⁻⁴")
        rows.append((name, KIND_RU.get(r["kind"], r["kind"]), str(r["n"]), _num(r["L"], 4),
                     " · ".join(parts) or "—"))
    return _table(["Семейство", "Тип", "Число", "L, м", "Сила (θ, K1, K2, K3)"], rows, [20, 20, 10, 12, 38])


# ------------------------------------------------------------------ HTML
def build_html(d):
    g, lin = d["geom"], d["lin"]
    sp_len = g["circumference"] / g["superperiods"]
    title = f"SynchrotronSim — паспорт установки: {d['label']} (решётка {d['lattice']})"
    parts = [f"периметр {_num(g['circumference'], 7)} м", f"{g['superperiods']} суперпериодов",
             f"νx/νy {_num(lin['nux'], 6)}/{_num(lin['nuy'], 6)}",
             f"εx {_num(d['rad']['ex0']*1e12, 4)} пм", f"U₀ {_num(d['rad']['U0'], 4)} МэВ"]
    err = d.get("errors")
    if err:
        parts.append(f"ансамбль {d['err']['n_seeds']} машин с ошибками ({d['turns']} об.)")
    elif not d["fast"]:
        parts.append(f"трекинг {d['turns']} об.")
    subtitle = " · ".join(parts)

    figs = [(_tbl_geometry(g), ta._table_h(len(g["kinds"]) + 11)),
            (_tbl_linear(lin), ta._table_h(8))]
    if d.get("struct"):                                       # отдельный график Твисса на каждую уник. структуру
        for st in d["struct"]["structures"]:
            figs.append((_fig_twiss_struct(st), 420))
    else:                                                     # структура не распознана — один график по суперпериоду
        figs.append((_fig_twiss(lin, sp_len), 420))
    figs += [
            (_tbl_chroma(d["chroma"]), ta._table_h(4)),
            (_tbl_radiation(d["rad"]), ta._table_h(15)),
            (_tbl_detune(d["detune"], d.get("mom")), ta._table_h(5)),
            (_fig_tunes(lin, d.get("foot"), d.get("fma")), 620)]
    if err:                                                    # ансамбль машин с ошибками
        n_erows = (10 if d["err"].get("corrected") else 6) + len(err["tol"])
        figs.append((_tbl_errors(d["err"], err["tol"]), ta._table_h(n_erows)))
    if err or not d["fast"]:                                   # ДА / диффузия / импульс (идеал или ансамбль)
        figs += [(_fig_da(d["da"]), 520),
                 (_fig_fma_space(d["fma"]), 520),
                 (_fig_momentum(d["mom"], d["rad"].get("delta_rf")), 400)]
    figs.append((_tbl_families(d["fam"]), ta._table_h(len(d["fam"]) + 1)))
    return ta._assemble(d["lattice"], 0, title, subtitle, figs)


def _cache_path(lattice):
    return os.path.join(RESULTS, f"facility_cache_exp{lattice}.pkl")


def main():
    import pickle
    ap = argparse.ArgumentParser(description="паспорт установки — HTML-отчёт с анализом решётки")
    ap.add_argument("lattice", type=int, help="id решётки (25, 26, 2, 1)")
    ap.add_argument("--turns", type=int, default=512, help="оборотов в сканах динамической апертуры/карты частот")
    ap.add_argument("--fast", action="store_true", help="только аналитика (без трекинга)")
    ap.add_argument("--rebuild", action="store_true",
                    help="пересобрать HTML на КЭШЕ трекинга (аналитика свежая, обороты не гоняются)")
    ap.add_argument("--errors", action="store_true",
                    help="учесть ошибки установки: ансамбль случайных машин по допускам из element_errors.csv")
    ap.add_argument("--correct", action="store_true",
                    help="выправить замкнутую орбиту корректорами на каждой машине (только с --errors)")
    ap.add_argument("--seeds", type=int, default=32, help="машин в ансамбле ошибок (--errors)")
    ap.add_argument("--nsl", type=int, default=40,
                    help="срезов магнита в трекинге (40 — паспорт идеала; ~12 для ансамбля: ×3 быстрее, "
                         "граница ДА не сдвигается)")
    ap.add_argument("--cutoff", type=float, default=2.5, help="обрезка гауссовых ошибок, ед. σ")
    ap.add_argument("--errors-table", default=None, help="путь к таблице допусков (по умолч. element_errors.csv)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--no-open", action="store_true")
    a = ap.parse_args()
    tag = (' + ошибки' if a.errors else '') + (' + коррекция' if a.errors and a.correct else '')
    print(f"=== паспорт установки: решётка {a.lattice}{tag} ===", flush=True)
    if a.errors:                                             # ансамбль машин с ошибками
        tol = errmod.load_tolerances(a.errors_table)
        print(f"  допуски: {errmod.error_summary(tol)}", flush=True)
        print(f"  срезов магнита nsl={a.nsl}", flush=True)
        d = analyze(a.lattice, turns=a.turns, nsl=a.nsl,
                    errors_cfg=dict(n_seeds=a.seeds, tol=tol, cutoff=a.cutoff, correct=a.correct))
    elif a.rebuild:                                          # быстрая пересборка: трекинг из кэша
        cp = _cache_path(a.lattice)
        if not os.path.exists(cp):
            sys.exit(f"[!] нет кэша {cp} — сначала полный прогон без --rebuild")
        with open(cp, "rb") as f:
            d = analyze(a.lattice, track_cache=pickle.load(f))
    else:
        d = analyze(a.lattice, turns=a.turns, fast=a.fast, nsl=a.nsl)
        if not a.fast:                                        # сохранить трекинг для будущих --rebuild
            with open(_cache_path(a.lattice), "wb") as f:
                pickle.dump({k: d[k] for k in ("da", "fma", "mom", "ranges", "turns")}, f)
    g, lin, r = d["geom"], d["lin"], d["rad"]
    print(f"  {d['label']}: периметр {g['circumference']:.3f} м, νx/νy {lin['nux']:.4f}/{lin['nuy']:.4f}, "
          f"εx {r['ex0']*1e12:.2f} пм, U₀ {r['U0']:.4f} МэВ, τx {r['taux']*1e3:.3f} мс", flush=True)
    if a.errors:
        e = d["err"]
        if e.get("corrected"):
            print(f"  ансамбль {e['n_seeds']} машин (коррекция, {e['n_corr']} корректоров, протяжка "
                  f"{e.get('n_threaded', 0)}): неустойчивых {e['n_unstable']}, орбита СКО x "
                  f"{e['co_x_rms0']*1e3:.3f}→{e['co_x_rms']*1e3:.3f} / y "
                  f"{e['co_y_rms0']*1e3:.3f}→{e['co_y_rms']*1e3:.3f} мм; углы корр. СКО x "
                  f"{e['corr_x_rms']*1e6:.1f}/y {e['corr_y_rms']*1e6:.1f} мкрад; ДА +x "
                  f"{d['da']['edge_pos'].max()*1e3:.2f} мм (50% машин)", flush=True)
        else:
            print(f"  ансамбль {e['n_seeds']} машин: неустойчивых {e['n_unstable']}, орбита СКО "
                  f"x {e['co_x_rms']*1e3:.3f}/y {e['co_y_rms']*1e3:.3f} мм; ДА +x "
                  f"{d['da']['edge_pos'].max()*1e3:.2f} мм (50% машин)", flush=True)
    elif not a.fast:
        print(f"  динамическая апертура +x {d['da']['edge_pos'].max()*1e3:.2f} мм, "
              f"импульсная +{d['mom']['acc_pos']*100:.2f}/{d['mom']['acc_neg']*100:.2f} %", flush=True)
    os.makedirs(HTML_DIR, exist_ok=True)
    suffix = ("_errors_corr" if a.correct else "_errors") if a.errors else ""
    out = a.out or os.path.join(HTML_DIR, f"facility_exp{a.lattice}{suffix}.html")
    with open(out, "w", encoding="utf-8") as f:
        f.write(build_html(d))
    print(f"  → {out}", flush=True)
    if not a.no_open:
        webbrowser.open("file://" + os.path.abspath(out))


if __name__ == "__main__":
    main()
