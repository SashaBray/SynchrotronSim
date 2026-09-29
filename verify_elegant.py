# -*- coding: utf-8 -*-
"""Опыты сверки SynchrotronSim с elegant на одной установке (ключ --machine, см. elegant_diag.MACHINES).

Каждый опыт считает показатель ОДНИМ алгоритмом из ring_analysis поверх двух трекеров — SynchrotronSim
(продакшн-ядро C++: ring_cpp.CppRing → MatrixTracker --ring-serve) и elegant (elegant_diag.ElegantRing) — и
дополнительно встроенным инструментом elegant, если он есть. Python-ядро (ring_analysis.Ring) трекинг опытов
не ведёт: с ним C++ сверяется опытом kernel.
Результаты: workspace/results/verify_elegant/<опыт>/ — картинки PNG, таблицы CSV (запятая, «;»), кэш .npz.

  py -3.14 verify_elegant.py kernel      C++-ядро против Python-эталона на решётке установки
  py -3.14 verify_elegant.py linear      β, α, фаза, η вдоль кольца; частоты; αc; матрица оборота
  py -3.14 verify_elegant.py matrices    матрицы 6×6 отдельных элементов
  py -3.14 verify_elegant.py analytic    аналитические матрицы элементов и аффинные матрицы стыков
  py -3.14 verify_elegant.py response    матрицы отклика: орбита на корректоры, коррекция SVD, допуска
  py -3.14 verify_elegant.py chrom       ν(δ): хроматичность 1–3 порядка, αc2 (+ отображение ring_map)
  py -3.14 verify_elegant.py integrator  частоты при разном числе срезов и порядке интегратора
  py -3.14 verify_elegant.py tbt         пооборотный трекинг: фазовые портреты, спектр, частоты от амплитуды
  py -3.14 verify_elegant.py long        ВЧ и излучение: U0, синхротронная частота, времена затухания
  py -3.14 verify_elegant.py emit        равновесный эмиттанс, энергоразброс, длина сгустка
  py -3.14 verify_elegant.py fma         динамическая апертура, карта частот и сетка резонансов

Установка задаётся ключом --machine: sila (описание + линия, решётка 1 БД, 6 ГэВ) или siberia (чужой
s2lattice.lte целиком, решётка 2 БД, 3 ГэВ). Результаты чужой установки — в подпапке её имени.
"""
import os, sys, json, math, argparse, datetime

try:                                          # сводки печатают ν, ε, τ, Δ — консоль Windows их не кодирует
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import ring_analysis as ra
import ring_cpp as rc
import elegant_diag as ed

ROOT = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(ROOT, "workspace", "results", "verify_elegant")
C_AN = "#2f9e44"                                            # аналитика SynchrotronSim
C_SS, C_EL, C_TW = "#1f5fbf", "#d9480f", "#7a7a7a"          # SynchrotronSim / elegant трекинг / elegant twiss


# ---------------------------------------------------------------- установка опыта
# Все опыты параметризованы установкой из elegant_diag.MACHINES: решётка БД и энергия для SynchrotronSim,
# файл .lte и имя линии для elegant, ВЧ, корректоры и амплитуды. Переключается ключом --machine.
MACHINE = "sila"
M = ed.machine(MACHINE)


def use_machine(name):
    global MACHINE, M
    MACHINE, M = name, ed.machine(name)
    ed.E_MEV = M["e_mev"]                       # write_particles по умолчанию берёт её


def _dir(name):
    d = os.path.join(OUT, name) if MACHINE == "sila" else os.path.join(OUT, MACHINE, name)
    os.makedirs(d, exist_ok=True); return d


def _stem(name):
    """Имя файла elegant, своё у каждой установки (иначе опыты затирают файлы друг друга)."""
    return name if MACHINE == "sila" else f"{MACHINE}_{name}"


def _ring_kwargs(**kw):
    """Аргументы кольца SynchrotronSim для текущей установки."""
    kw.setdefault("start_before_first_m", ed.lead_drift(*ed.machine_lattice(M)))
    kw.setdefault("design_mev", M["e_mev"])
    return kw


def _env(extra=None):
    """Переменные ядра SynchrotronSim под интегратор установки (+ точечные добавки опыта)."""
    n, o = M["slices"]
    e = {"SS_SECTOR_NSL": str(n), "SS_SECTOR_ORDER": str(o)}
    if extra:
        e.update(extra)
    return e


_RINGS = {}


def _ring(default_slices=False, nsl=None, order=None, **kw):
    """Кольцо SynchrotronSim на ПРОДАКШН-ЯДРЕ C++ (ring_cpp.CppRing поверх MatrixTracker --ring-serve).
    Интегратор — установки (default_slices — умолчания ядра: 5 срезов, порядок 4; nsl/order — явно).
    Кольца кэшируются: одно C++-кольцо на набор параметров."""
    if not default_slices:
        n0, o0 = M["slices"]
        nsl, order = nsl or n0, order or o0
    kw = _ring_kwargs(**kw)
    key = (MACHINE, nsl, order, repr(sorted(kw.items())))
    if key not in _RINGS:
        _RINGS[key] = rc.CppRing(M["lattice"], nsl=nsl, order=order, **kw)
    return _RINGS[key]


def _tbt(states, n_turns, limit=None, **kw):
    """Пооборотный трекинг SynchrotronSim (C++, частицы параллельно): (n, n_turns+1, 6)."""
    return _ring(**kw).track_tbt(np.asarray(states, float), n_turns, limit=limit)


def _er(tag, default_slices=False, **kw):
    """Трекер elegant на том же интеграторе (число срезов ставится в описании устройств)."""
    lk = dict(kw.pop("lattice_kwargs", None) or {})
    if not default_slices:
        lk.setdefault("integrator", M["slices"])
    return ed.ElegantRing(tag, machine_name=MACHINE, lattice_kwargs=lk, **kw)


def _pool(tasks, **kw):
    """Задачи [(имя, args, kwargs, env, ring_kwargs)] на C++-кольцах: имя 'ring.метод' — метод кольца, иначе
    функция ring_analysis с кольцом первым аргументом; env задаёт интегратор (SS_SECTOR_NSL/ORDER; None —
    умолчания ядра). Частицы C++-ядро считает параллельно само, поэтому задачи идут подряд."""
    out = []
    for task in tasks:
        name, args, kwargs, env, ring_kwargs = (tuple(task) + (None, None))[:5]
        rk = dict(ring_kwargs or {})
        if env is None:
            ring = _ring(default_slices=True, **rk)
        else:
            ring = _ring(nsl=int(env["SS_SECTOR_NSL"]), order=int(env["SS_SECTOR_ORDER"]), **rk)
        fn = getattr(ring, name[5:]) if name.startswith("ring.") else (lambda *a, **k: getattr(ra, name)(ring, *a, **k))
        out.append(fn(*args, **(kwargs or {})))
    return out


def _lattice(**kw):
    """Собрать .lte установки для встроенных инструментов elegant; вернуть имя файла."""
    kw.setdefault("integrator", M["slices"])
    if M["source"]:
        kw.update(source=M["source"], beamline=M["beamline"], use_beamline=M["use_beamline"])
    else:
        kw.setdefault("config", M["config"])
    return os.path.basename(ed.write_lattice(_stem(M["lte_name"]), **kw))


def _setup(order, lte=None, extra=""):
    """Строка &run_setup для установки."""
    return (f'&run_setup lattice="{lte or _stem(M["lte_name"])}", use_beamline="{M["use_beamline"]}", '
            f'p_central_mev={M["e_mev"]:g}, default_order={order}, {extra} &end\n')


def _period():
    """Длина одного супериода, м — по ней режется график функций Твисса."""
    return ed.line_length(*ed.machine_lattice(M)) / M["superperiods"]


def _csv(path, header, rows):
    fmt = lambda v: (f"{v:.10g}".replace(".", ",") if isinstance(v, (float, np.floating)) else str(v))
    with open(path, "w", encoding="utf-8-sig", newline="") as f:
        f.write(";".join(header) + "\n")
        for r in rows:
            f.write(";".join(fmt(v) for v in r) + "\n")


def _cached(path, fn, refresh):
    if os.path.exists(path) and not refresh:
        return dict(np.load(path, allow_pickle=True))
    res = fn()
    np.savez(path, **{k: np.asarray(v) for k, v in res.items()})
    return res


# ============================================================================ линейная оптика
def exp_linear(refresh=False):
    d = _dir("linear")

    def run_ss():
        ring = _ring()
        tw = ra.twiss(ring)
        tw["names"] = np.array(ring.names + ["конец оборота"])
        return tw

    def run_el():
        er = _er("lin")
        tw = ra.twiss(er)
        tw["names"] = np.array(er.names)
        return tw

    def run_twi():
        _lattice()
        ed.run(_stem("twiss_closed") + ".ele", _setup(3) +
               '&twiss_output filename="%s.twi", matched=1, radiation_integrals=1, higher_order_chromaticity=1, &end\n')
        twi = _stem("twiss_closed") + ".twi"
        c = ed.columns(twi, ["s", "betax", "alphax", "psix", "etax", "etaxp", "betay", "alphay",
                             "psiy", "ElementType"])
        mask = np.array([t not in ed._MAGNET_TYPES_SKIP for t in c["ElementType"]]); mask[-1] = True
        out = {k: v[mask] for k, v in c.items() if k != "ElementType"}
        out.update(ed.parameters(twi, ["nux", "nuy", "alphac", "dnux/dp", "dnuy/dp"]))
        return out

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el, refresh)
    tw = _cached(os.path.join(d, "elegant_twiss.npz"), run_twi, refresh)

    s = el["s"]
    if len(ss["s"]) != len(s):
        raise SystemExit(f"число магнитов не совпало: SynchrotronSim {len(ss['s'])}, elegant {len(s)}")
    ds = float(np.max(np.abs(ss["s"] - s)))

    # --- рисунок 1: функции на первом супериоде
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    sel = s <= _period()
    for a, key, lab in ((ax[0], "betax", "βx, м"), (ax[1], "betay", "βy, м"), (ax[2], "etax", "ηx, м")):
        a.plot(tw["s"][sel], tw[key][sel], "--", color=C_TW, lw=1.2, label="elegant twiss_output (аналитика)")
        a.plot(s[sel], el[key][sel], "o", color=C_EL, ms=4, label="elegant, трекинг")
        a.plot(ss["s"][sel], ss[key][sel], "-", color=C_SS, lw=1.4, label="SynchrotronSim, трекинг")
        a.set_ylabel(lab); a.grid(alpha=0.3)
    ax[0].legend(loc="upper right", fontsize=9)
    ax[2].set_xlabel("s, м (значения на выходе магнитов)")
    fig.suptitle(f"Функции Твисса на первом супериоде ({M['label']})")
    fig.tight_layout(); fig.savefig(os.path.join(d, "twiss_superperiod.png"), dpi=130); plt.close(fig)

    # --- рисунок 2: расхождения вдоль всего кольца
    fig, ax = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    for key, lab, a in (("betax", "Δβx/βx", ax[0]), ("betay", "Δβy/βy", ax[1])):
        a.plot(s, (tw[key] - el[key]) / el[key] * 100, color=C_TW, lw=0.8, label="twiss_output − трекинг elegant")
        a.plot(s, (ss[key] - el[key]) / el[key] * 100, color=C_SS, lw=1.0, label="SynchrotronSim − трекинг elegant")
        a.set_ylabel(lab + ", %"); a.grid(alpha=0.3)
    ax[0].legend(loc="upper right", fontsize=9)
    ax[2].plot(s, (tw["etax"] - el["etax"]) * 1e3, color=C_TW, lw=0.8)
    ax[2].plot(s, (ss["etax"] - el["etax"]) * 1e3, color=C_SS, lw=1.0)
    ax[2].set_ylabel("Δηx, мм"); ax[2].grid(alpha=0.3)
    for key, ls in (("psix", "-"), ("psiy", ":")):
        ax[3].plot(s, (tw[key] - el[key]) / (2 * np.pi), ls, color=C_TW, lw=0.8)
        ax[3].plot(s, (ss[key] - el[key]) / (2 * np.pi), ls, color=C_SS, lw=1.0)
    ax[3].set_ylabel("Δ фазы, ед. 2π\n(сплошная x, пунктир y)"); ax[3].grid(alpha=0.3)
    ax[3].set_xlabel("s, м")
    fig.suptitle("Расхождение с трекингом elegant вдоль кольца")
    fig.tight_layout(); fig.savefig(os.path.join(d, "twiss_difference.png"), dpi=130); plt.close(fig)

    rel = lambda a, b: float(np.max(np.abs((a - b) / b)))
    rows = [
        ("νx", float(ss["nux"]), float(el["nux"]), float(tw["nux"])),
        ("νy", float(ss["nuy"]), float(el["nuy"]), float(tw["nuy"])),
        ("αc", float(ss["alphac"]), float(el["alphac"]), float(tw["alphac"])),
        ("βx на плоскости наблюдения, м", float(ss["betax"][-1]), float(el["betax"][-1]), float(tw["betax"][-1])),
        ("βy на плоскости наблюдения, м", float(ss["betay"][-1]), float(el["betay"][-1]), float(tw["betay"][-1])),
        ("max βx, м", float(ss["betax"].max()), float(el["betax"].max()), float(tw["betax"].max())),
        ("max βy, м", float(ss["betay"].max()), float(el["betay"].max()), float(tw["betay"].max())),
        ("max ηx, м", float(ss["etax"].max()), float(el["etax"].max()), float(tw["etax"].max())),
        ("замкнутая орбита x, мкм", float(ss["closed_orbit"][0]) * 1e6, float(el["closed_orbit"][0]) * 1e6, 0.0),
    ]
    _csv(os.path.join(d, "summary.csv"), ["показатель", "SynchrotronSim", "elegant трекинг", "elegant twiss_output"], rows)
    summ = dict(
        max_rel_beta_x_ss=rel(ss["betax"], el["betax"]), max_rel_beta_y_ss=rel(ss["betay"], el["betay"]),
        max_rel_beta_x_twiss=rel(tw["betax"], el["betax"]), max_rel_beta_y_twiss=rel(tw["betay"], el["betay"]),
        max_deta_x_mm_ss=float(np.max(np.abs(ss["etax"] - el["etax"]))) * 1e3,
        max_s_mismatch_m=ds)
    json.dump(summ, open(os.path.join(d, "summary.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)
    for r in rows:
        print(f"{r[0]:32s} SS {r[1]:.8g}   elegant-трек {r[2]:.8g}   elegant-twiss {r[3]:.8g}")
    print(json.dumps(summ, ensure_ascii=False, indent=1))


def _el_parallel(fn, items, workers=8):
    """elegant — внешние процессы, поэтому хватает потоков (ожидание subprocess не держит GIL)."""
    from concurrent.futures import ThreadPoolExecutor
    with ThreadPoolExecutor(max_workers=workers) as ex:
        return list(ex.map(fn, items))


def _fdist(a, b):
    """Расстояние дробных частот по модулю 1."""
    return np.abs((np.asarray(a, float) - np.asarray(b, float) + 0.5) % 1.0 - 0.5)


def _unwrap_frac(v, i0):
    """Дробные частоты без скачков через 0/1: развёртка по непрерывности от точки i0 в обе стороны."""
    v = np.asarray(v, float).copy()
    for rng in (range(i0 + 1, len(v)), range(i0 - 1, -1, -1)):
        prev = v[i0]
        for i in rng:
            if np.isfinite(v[i]):
                v[i] -= np.round(v[i] - prev); prev = v[i]
    return v


# ============================================================================ матрицы элементов
def exp_matrices(refresh=False):
    """Матрица 6×6 каждого магнита вместе с дрейфом перед ним: трекинг SynchrotronSim, трекинг elegant,
    аналитика elegant matrix_output. Сравнение по поперечно-дисперсионному блоку (x, px, y, py, δ)."""
    d = _dir("matrices")
    lin = _dir("linear")
    ss = np.load(os.path.join(lin, "ss.npz"), allow_pickle=True)
    el = np.load(os.path.join(lin, "elegant_track.npz"), allow_pickle=True)

    def run_mat():
        _lattice()
        ed.run(_stem("matrix_closed") + ".ele", _setup(1) +
               '&matrix_output SDDS_output="%s.mat", SDDS_output_order=1, individual_matrices=1, &end\n')
        cols = [f"R{i}{j}" for i in range(1, 7) for j in range(1, 7)]
        c = ed.columns(_stem("matrix_closed") + ".mat", ["ElementName", "ElementType"] + cols)
        R = np.stack([c[k] for k in cols], 1).reshape(-1, 6, 6)
        return dict(R=R, names=np.array(c["ElementName"]), types=np.array(c["ElementType"]))

    mo = _cached(os.path.join(d, "elegant_matrix_output.npz"), run_mat, refresh)
    # композиция аналитики elegant: дрейфы перед магнитом × магнит
    comp = []; names = []; types = []
    acc = np.eye(6)
    for R, nm, tp in zip(mo["R"], mo["names"], mo["types"]):
        if tp in ("MARK", "WATCH"):
            continue
        acc = R @ acc
        if tp not in ed._MAGNET_TYPES_SKIP:
            comp.append(acc); names.append(str(nm)); types.append(str(tp)); acc = np.eye(6)
    comp = np.array(comp)

    def elem(Rcum):
        out = [Rcum[0]] + [Rcum[j] @ np.linalg.inv(Rcum[j - 1]) for j in range(1, len(Rcum) - 1)]
        return np.array(out)

    import ring_matrix as rm
    ring = _ring()
    lmr = rm.line_matrices(ring)
    E_an = np.array([rm.compose(rm.junction_matrix(ring.plane if j == 0 else ring.mags[j - 1], m),
                                rm.element_matrix(m))[1] for j, m in enumerate(ring.mags)])
    E_ss = elem(ss["R"]); E_el = elem(el["R"])
    n = min(len(E_ss), len(E_el), len(comp))
    E_ss, E_el, comp = E_ss[:n], E_el[:n], comp[:n]
    # elegant: столбец 5 — путь s; у трекинга t — путь (track_elements), знаки совпадают
    idx = [0, 1, 2, 3, 5]
    sub = lambda A: A[:, idx][:, :, idx]
    n2 = min(n, len(E_an))
    d_ss = np.abs(sub(E_ss) - sub(E_el)); d_an = np.abs(sub(comp) - sub(E_el))
    d_pair = np.abs(sub(E_an[:n2]) - sub(comp[:n2]))          # аналитика против аналитики

    # по типам элементов
    # первый магнит включает дрейф от точки наблюдения, а его длина у SynchrotronSim больше на половину
    # невязки замыкания геометрии БД — показывается отдельной строкой, чтобы не заслонять тела элементов
    kinds = sorted(set(types))
    rows = []
    rest = np.arange(n) > 0
    for kd in kinds:
        m = np.array([t == kd for t in types]) & rest
        rows.append((kd, int(m.sum()), float(d_ss[m].max()), float(np.median(d_ss[m].max(axis=(1, 2)))),
                     float(d_pair[m[:n2]].max()), float(d_an[m].max())))
    rows.append((f"первый: {names[0]} с дрейфом от s = 0", 1, float(d_ss[0].max()), float(d_ss[0].max()),
                 float(d_pair[0].max()), float(d_an[0].max())))
    _csv(os.path.join(d, "by_type.csv"),
         ["тип elegant", "число", "трекинг: max |ΔR| SynchrotronSim − elegant",
          "трекинг: медиана max |ΔR|", "аналитика: max |ΔR| ring_matrix − matrix_output",
          "внутри elegant: max |ΔR| matrix_output − трекинг"], rows)

    # представительные элементы
    reps = ["QF8B", "QD2A", "DQ1B", "DQ2C_1", "DL1A_5", "SF2A", "OD1B"]     # СИЛА: типичные устройства
    if not any(nm in names for nm in reps):          # чужая решётка: первое вхождение каждого устройства
        reps = []
        for nm in names[1:]:
            if nm not in reps:
                reps.append(nm)
    terms = [("R11", 0, 0), ("R12", 0, 1), ("R21", 1, 0), ("R22", 1, 1), ("R33", 2, 2), ("R34", 2, 3),
             ("R43", 3, 2), ("R16", 0, 5), ("R26", 1, 5), ("R51", 4, 0), ("R56", 4, 5)]
    rep_rows = []
    for nm in reps:
        if nm not in names:
            continue
        j = names.index(nm)
        for tn, a, b in terms:
            rep_rows.append((nm, tn, float(E_ss[j, a, b]), float(E_el[j, a, b]), float(comp[j, a, b])))
    _csv(os.path.join(d, "representative.csv"), ["элемент", "член", "SynchrotronSim трекинг", "elegant трекинг",
                                                 "elegant matrix_output"], rep_rows)

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    lab = ["x", "px", "y", "py", "δ"]
    for a, D, title in ((axs[0], d_ss[1:], "SynchrotronSim − elegant (оба трекинг)"),
                        (axs[1], d_an[1:], "matrix_output (аналитика) − трекинг elegant")):
        im = a.imshow(np.log10(np.maximum(D.max(axis=0), 1e-16)), cmap="magma", vmin=-12, vmax=-2)
        a.set_xticks(range(5)); a.set_xticklabels(lab); a.set_yticks(range(5)); a.set_yticklabels(lab)
        for i in range(5):
            for k in range(5):
                a.text(k, i, f"{np.log10(max(D.max(axis=0)[i, k], 1e-16)):.0f}", ha="center", va="center",
                       color="w", fontsize=9)
        a.set_title(title + "\nlog10 max|ΔR| по магнитам (кроме первого)"); fig.colorbar(im, ax=a)
    fig.suptitle("Матрицы элементов (магнит + дрейф перед ним)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "element_matrices.png"), dpi=130); plt.close(fig)
    for r in rows:
        print(f"{r[0]:8s} {r[1]:5d}  трекинг {r[2]:.2e} (медиана {r[3]:.2e})  аналитика {r[4]:.2e}  "
              f"внутри elegant {r[5]:.2e}")


# ============================================================================ матричный механизм (аналитика)
MISALIGN_RAD = 1e-4                                  # поворот магнита для отклика орбиты
CONVERGED = {"SS_SECTOR_NSL": "40", "SS_SECTOR_ORDER": "6"}   # сошедшийся интегратор для честной сверки


def _lte_matrices():
    """Матрицы, собранные ПРЯМО по .lte (без геометрии БД) — проверка самих формул против elegant."""
    import ring_matrix as rm
    lat = ed.machine_lattice(M)
    par = ed.element_params(lat[0], lat[1]); names = ed.line_names(*lat)
    R = np.eye(6); Rs = []
    for n in names:
        q = par[n.upper()]; t = q["type"].upper(); L = q.get("L", 0.0)
        if t == "CSBEND":
            h = q.get("ANGLE", 0.0) / L
            Mel = (rm.edge(math.tan(q.get("E2", 0.0)) * h) @ rm.sector_body(L, h, q.get("K1", 0.0))
                   @ rm.edge(math.tan(q.get("E1", 0.0)) * h))
        elif t == "KQUAD":
            Mel = rm.quad(L, q.get("K1", 0.0))
        elif t in ("MONI", "MARK"):
            Mel = np.eye(6)
        else:
            Mel = rm.drift(L)                             # секступоль и октуполь линейно — дрейфы
        R = Mel @ R
        if t not in ("EDRIFT", "DRIF", "MONI", "MARK"):
            Rs.append(R.copy())
    Rs.append(R.copy())
    return ra.twiss_from_matrices(np.array(Rs), np.arange(len(Rs), dtype=float), ed.line_length(*lat))


def exp_analytic(refresh=False):
    """Аналитические матрицы элементов и аффинные матрицы стыков против elegant и против трекинга."""
    import ring_matrix as rm
    d = _dir("analytic")
    tw_el = np.load(os.path.join(_dir("linear"), "elegant_twiss.npz"), allow_pickle=True)
    tr = np.load(os.path.join(_dir("linear"), "ss.npz"), allow_pickle=True)

    def run_mat():
        ring = _ring()
        tw = rm.twiss(ring); lm = rm.line_matrices(ring)
        return dict(s=tw["s"], betax=tw["betax"], betay=tw["betay"], etax=tw["etax"],
                    alphax=tw["alphax"], alphay=tw["alphay"], etaxp=tw["etaxp"],
                    nux=tw["nux"], nuy=tw["nuy"], alphac=tw["alphac"], orbit=tw["closed_orbit"],
                    C_transverse=np.max(np.abs(lm["one_turn_C"][:4])), C_path=lm["one_turn_C"][4],
                    circumference=ring.circumference)

    def run_resp():
        """Отклик орбиты на поворот магнита: матрицы против трекинга на сошедшемся интеграторе."""
        res = _pool([("analytic_orbit_response", (MISALIGN_RAD,), {}, CONVERGED, _ring_kwargs())],
                    processes=1)
        return res[0]

    mat = _cached(os.path.join(d, "ss_matrix.npz"), run_mat, refresh)
    resp = _cached(os.path.join(d, "response.npz"), run_resp, refresh)
    lte = _lte_matrices()

    # --- проверка матриц стыков на четырёх геометриях граней
    from beam_injection import PlaneMagnet
    import matrix_tracker as mt
    junc = []; geoms = []
    ax = np.array([1.0, 0.0, 0.0]); vv = np.array([0.0, 0.0, 1.0])
    for alpha, roll, off in ((0.0, 0.0, 0.0), (5e-3, 0.0, 0.0), (5e-2, 0.0, 1e-3), (2e-2, 3e-2, 5e-4)):
        A = PlaneMagnet(np.zeros(3), ax)
        b_ax = mt._rot(ax, vv, alpha)
        B = PlaneMagnet(np.array([0.8, off, 0.0]), b_ax)
        if roll:
            u = mt._rot(B.uhat, b_ax, roll)
            B.uhat = u; B.vhat = np.cross(b_ax, u); B.fr_in = (B.center, B.axis, u, B.vhat)
        C, R = rm.junction_matrix(A, B)
        z0 = np.zeros(6); h = 1e-6
        pts = np.vstack([z0] + [z0 + h * e for e in np.eye(6)] + [z0 - h * e for e in np.eye(6)])
        out = _ring().junction(pts, A.fr_out, B.fr_in)          # пролёт между плоскостями — C++-ядро
        Cn = out[0]; Rn = ((out[1:7] - out[7:13]) / (2 * h)).T
        junc.append((f"клин {alpha * 1e3:g} мрад, крен {roll * 1e3:g} мрад, смещение {off * 1e3:g} мм".replace(".", ","),
                     float(np.max(np.abs(C - Cn))), float(np.max(np.abs(R - Rn))), float(np.max(np.abs(C[:4])))))
        geoms.append((junc[-1][0], A, B, C, R))

    # --- стык в высших порядках: ряд ring_map (порядок 5) и аффинная матрица против точного пролёта C++
    import ring_map as rmap
    da = rmap.DA(6, 5)
    direction = np.array([1.0, 0.5, -0.7, 0.3, 0.0, 0.8])
    junc_ho = []
    for label, A, B, C, R in geoms:
        J = rmap.junction(da, da.identity(), A.fr_out, B.fr_in)
        mono = lambda z: np.array([np.prod([zi ** e for zi, e in zip(z, ex)]) for ex in da.monos])
        for amp in (1e-3, 3e-3, 1e-2):
            z = amp * direction
            exact = _ring().junction([z], A.fr_out, B.fr_in)[0]
            lin = C + R @ z
            ser = np.array([c @ mono(z) for c in J])
            junc_ho.append((label, amp * 1e3, float(np.max(np.abs(lin[:4] - exact[:4]))),
                            float(np.max(np.abs(ser[:4] - exact[:4])))))

    # --- рисунки
    s = mat["s"]
    fig, axs = plt.subplots(3, 1, figsize=(11, 8.5), sharex=True)
    for a, key, lab in ((axs[0], "betax", "Δβx/βx"), (axs[1], "betay", "Δβy/βy")):
        a.plot(s, (mat[key] - tw_el[key]) / tw_el[key] * 1e6, color=C_SS, lw=0.9,
               label="матрицы SynchrotronSim − twiss_output elegant")
        a.plot(s, (tr[key] - tw_el[key]) / tw_el[key] * 1e6, color=C_EL, lw=0.6, alpha=0.8,
               label="трекинг SynchrotronSim (5 срезов) − twiss_output elegant")
        a.set_ylabel(lab + ", 10⁻⁶"); a.grid(alpha=0.3)
    axs[0].legend(fontsize=9)
    axs[2].plot(s, (mat["etax"] - tw_el["etax"]) * 1e6, color=C_SS, lw=0.9)
    axs[2].plot(s, (tr["etax"] - tw_el["etax"]) * 1e6, color=C_EL, lw=0.6, alpha=0.8)
    axs[2].set_ylabel("Δηx, мкм"); axs[2].set_xlabel("s, м"); axs[2].grid(alpha=0.3)
    fig.suptitle(f"Аналитические матрицы против twiss_output elegant "
                 f"(для сравнения — трекинг, {M['slices'][0]} срезов)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "analytic_vs_elegant.png"), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(3, 1, figsize=(11, 8), sharex=True)
    axs[0].plot(s, (mat["betax"] - tw_el["betax"]) / tw_el["betax"] * 1e6, color=C_SS, lw=0.8)
    axs[0].set_ylabel("Δβx/βx, 10⁻⁶")
    axs[1].plot(s, (mat["betay"] - tw_el["betay"]) / tw_el["betay"] * 1e6, color=C_SS, lw=0.8)
    axs[1].set_ylabel("Δβy/βy, 10⁻⁶")
    axs[2].plot(s, (mat["etax"] - tw_el["etax"]) * 1e9, color=C_SS, lw=0.8)
    axs[2].set_ylabel("Δηx, нм"); axs[2].set_xlabel("s, м")
    for a in axs:
        a.grid(alpha=0.3)
    fig.suptitle(f"Тот же остаток в собственном масштабе: периметр решётки БД отличается от линии .lte "
                 f"на {(float(mat['circumference']) - ed.line_length(*ed.machine_lattice(M))) * 1e6:.2f} мкм")
    fig.tight_layout(); fig.savefig(os.path.join(d, "analytic_residual.png"), dpi=130); plt.close(fig)

    # --- таблицы
    rows = [
        ("νx (матрицы по .lte)", float(lte["nux"]), float(tw_el["nux"]), abs(float(lte["nux"]) - float(tw_el["nux"]))),
        ("νy (матрицы по .lte)", float(lte["nuy"]), float(tw_el["nuy"]), abs(float(lte["nuy"]) - float(tw_el["nuy"]))),
        ("αc (матрицы по .lte)", float(lte["alphac"]), float(tw_el["alphac"]),
         abs(float(lte["alphac"]) / float(tw_el["alphac"]) - 1)),
        ("max |Δβx/βx| (матрицы по .lte)", "", "",
         float(np.max(np.abs((lte["betax"][:len(tw_el["betax"])] - tw_el["betax"]) / tw_el["betax"])))),
        ("νx (матрицы по решётке БД)", float(mat["nux"]), float(tw_el["nux"]),
         abs(float(mat["nux"]) - float(tw_el["nux"]))),
        ("νy (матрицы по решётке БД)", float(mat["nuy"]), float(tw_el["nuy"]),
         abs(float(mat["nuy"]) - float(tw_el["nuy"]))),
        ("max |Δβx/βx| (решётка БД)", "", "",
         float(np.max(np.abs((mat["betax"] - tw_el["betax"]) / tw_el["betax"])))),
        ("max |Δβy/βy| (решётка БД)", "", "",
         float(np.max(np.abs((mat["betay"] - tw_el["betay"]) / tw_el["betay"])))),
        ("max Δηx, нм (решётка БД)", "", "", float(np.max(np.abs(mat["etax"] - tw_el["etax"]))) * 1e9),
        ("периметр, м", float(mat["circumference"]), ed.line_length(*ed.machine_lattice(M)),
         float(mat["circumference"]) - ed.line_length(*ed.machine_lattice(M))),
        ("поперечный C за оборот (идеальная решётка)", float(mat["C_transverse"]), 0.0, float(mat["C_transverse"])),
        ("путь в C за оборот, м", float(mat["C_path"]), float(mat["circumference"]),
         abs(float(mat["C_path"]) - float(mat["circumference"]))),
    ]
    _csv(os.path.join(d, "summary.csv"), ["показатель", "SynchrotronSim, матрицы", "elegant / ожидание",
                                          "расхождение"], rows)
    _csv(os.path.join(d, "junctions.csv"),
         ["геометрия граней стыка", "max |ΔC| к точному отображению", "max |ΔR| к точному отображению",
          "поперечный C"], junc)
    _csv(os.path.join(d, "junctions_high_order.csv"),
         ["геометрия граней стыка", "амплитуда, мм", "аффинная матрица: max |Δz|, м",
          "ряд 5-го порядка (ring_map): max |Δz|, м"], junc_ho)
    r = resp
    _csv(os.path.join(d, "response.csv"),
         ["показатель", "матрицы", "трекинг (40 срезов, порядок 6)", "расхождение"],
         [("замкнутая орбита идеальной решётки, нм", float(r["co_matrix"][0]) * 1e9,
           float(r["co_track"][0]) * 1e9, abs(float(r["co_matrix"][0]) - float(r["co_track"][0])) * 1e9),
          (f"Δx от поворота магнита на {MISALIGN_RAD * 1e3:g} мрад, мкм", float(r["dx_matrix"]) * 1e6,
           float(r["dx_track"]) * 1e6, abs(float(r["dx_matrix"]) / float(r["dx_track"]) - 1))])
    for row in rows:
        print(f"{row[0]:44s} {row[1]!s:>18s} {row[2]!s:>18s} {row[3]:.3e}")
    for row in junc:
        print(f"стык {row[0]:52s} ΔC {row[1]:.1e}  ΔR {row[2]:.1e}")
    print(f"отклик орбиты: матрицы {float(r['dx_matrix']) * 1e6:.4f} мкм, трекинг "
          f"{float(r['dx_track']) * 1e6:.4f} мкм")


# ============================================================================ матрицы отклика
SEED_MISALIGN = 20260916
N_MISALIGNED = 20                                    # сколько квадруполей развернуть в опыте коррекции
SIGMA_YAW = 5e-5                                     # СКО поворота, рад
KEEP_SV = 100                                        # сколько сингулярных чисел оставить при обращении


def _kick_kinds():
    """Виды устройств SynchrotronSim, в которых стоят корректоры установки: по самому кольцу БД в месте
    первого вхождения каждого корректора (тип в .lte не годится: у СИЛА корректоры SH* описаны как KSEXT,
    а в БД это устройства вида kicker)."""
    ring = _ring()
    return tuple(sorted({ring.kinds[_magnet_index(n, 1)] for n in M["correctors"]}))


def _correctors(ring, lm):
    import ring_response as rr
    return rr.correctors(ring, lm, _kick_kinds())


def _kicker_overrides():
    """Подмена устройств-корректоров на KICKER в описании elegant."""
    return {n: "KICKER,L=%g" % M["corrector_len"] for n in M["correctors"]}


def _stations():
    import ring_response as rr
    return rr.bpm_stations(*ed.machine_lattice(M))


def _magnet_index(lte_name, occurrence=1, kinds=None):
    """Номер магнита в решётке БД по имени элемента .lte и номеру его вхождения в линию.

    Сопоставление — по положению: середина вхождения в линии .lte против середин магнитов кольца. По счёту
    магнитов нельзя: в БД СИЛА устройства, стоящие в одной точке, собраны в один магнит, и номера линии и
    кольца расходятся. kinds — ограничить поиск видами устройств (корректоры)."""
    lat = ed.machine_lattice(M)
    par = ed.element_params(lat[0], lat[1])
    s = 0.0; seen = 0; s_mid = None
    for n in ed.line_names(*lat):
        L = float(par[n.upper()].get("L", 0.0))
        if n.upper() == lte_name.upper():
            seen += 1
            if seen == occurrence:
                s_mid = s + L / 2; break
        s += L
    if s_mid is None:
        raise KeyError("%s #%d не найден в линии" % (lte_name, occurrence))
    ring = _ring()
    mid = 0.5 * (ring.s_in + ring.s_out)
    cand = [j for j in range(len(mid)) if kinds is None or ring.kinds[j] in kinds]
    j = min(cand, key=lambda k: abs(mid[k] - s_mid))
    if abs(mid[j] - s_mid) > 0.05:
        raise KeyError("%s #%d: магнит кольца в точке s = %.3f м не найден" % (lte_name, occurrence, s_mid))
    return j


def _elegant_orm(refresh=False):
    """Матрица отклика elegant (датчики MONI × корректоры) для обеих плоскостей."""
    d = _dir("response")

    def run():
        lte = _lattice(overrides=_kicker_overrides())
        stem = _stem("orm")
        ed.run(stem + ".ele", _setup(1, lte) +
               '&steering_element element_type="KICKER", item="HKICK", plane="h", tweek=1e-6, &end\n'
               '&steering_element element_type="KICKER", item="VKICK", plane="v", tweek=1e-6, &end\n'
               '&correct mode="orbit", method="global", n_iterations=0, n_xy_cycles=1, do_correction[0]=0,0, &end\n'
               '&correction_matrix_output response[0] = "%s.hrm","%s.vrm","%s.hvrm","%s.vhrm", &end\n')
        import subprocess
        cols = subprocess.run([ed._exe("sddsquery"), stem + ".hrm", "-col"], cwd=ed.DIAG, capture_output=True,
                              text=True).stdout.split()
        corr = [c for c in cols if c not in ("BPMName", "s")]
        h = ed.columns(stem + ".hrm", ["s"] + corr); v = ed.columns(stem + ".vrm", corr)
        return dict(s=h["s"], names=np.array(corr),
                    H=np.stack([h[k] for k in corr], 1), V=np.stack([v[k] for k in corr], 1))

    return _cached(os.path.join(d, "elegant_orm.npz"), run, refresh)


def _orm_tracking_pairs(refresh=False):
    """Столбец матрицы отклика ТРЕКИНГОМ в обоих кодах: кик в одном корректоре, замкнутая орбита.
    В elegant кик ставится в одно вхождение элемента через &alter_elements."""
    import ring_matrix as rm
    import ring_response as rr
    import copy
    d = _dir("response")

    def run():
        theta = 1e-6
        ov = _kicker_overrides()
        ring = _ring(); lm = rm.line_matrices(ring); corr = _correctors(ring, lm)
        Mt = lm["one_turn_R"]; inv = np.linalg.inv(np.eye(4) - Mt[:4, :4])
        co0 = ra.closed_orbit(ring)
        co0_el = ra.closed_orbit(_er("ormb", lattice_kwargs=dict(overrides=ov)))
        by_index = {j: (j, s_k, Rc) for j, s_k, Rc in corr}
        rows = []
        for k, name in enumerate(M["correctors"][:2]):       # по первому вхождению двух корректоров
            j, s_k, Rc = by_index[_magnet_index(name, 1, _kick_kinds())]
            e = np.zeros(6); e[1] = theta
            z_an = (inv @ (Mt @ np.linalg.inv(Rc) @ e)[:4])[0]
            pert = ra.kicked(ring, j, theta)
            z_ss = (ra.closed_orbit(pert) - co0)[0]
            er = _er("ormt%d" % k, lattice_kwargs=dict(overrides=ov),
                     pre_commands="&alter_elements name=%s, item=HKICK, value=%g, "
                                  "start_occurence=1, end_occurence=1, &end\n" % (name, theta))
            z_el = (ra.closed_orbit(er) - co0_el)[0]
            rows.append((s_k, z_an, z_ss, z_el))
        return dict(rows=np.array(rows), theta=theta)

    return _cached(os.path.join(d, "tracking_pairs.npz"), run, refresh)


def exp_response(refresh=False):
    """Матрицы отклика: орбита на корректоры (против elegant и против трекинга), коррекция по SVD, допуска."""
    import ring_matrix as rm
    import ring_response as rr
    d = _dir("response")
    el = _elegant_orm(refresh)
    ring = _ring(); lm = rm.line_matrices(ring)
    st = _stations()
    kinds = _kick_kinds()
    ours = {pl: rr.orbit_response(ring, lm, st, pl, kinds) for pl in ("x", "y")}
    picks = tuple(sorted({0, len(_correctors(ring, lm)) // 2}))

    def run_track():
        return dict(rows=np.array(_pool([("orm_tracking_check", (picks,), {"kinds": kinds}, CONVERGED,
                                          _ring_kwargs())], processes=1)[0], dtype=object))

    trk = _cached(os.path.join(d, "tracking_check.npz"), run_track, refresh)["rows"]
    pairs = _orm_tracking_pairs(refresh)["rows"]

    # --- сверка с elegant
    rows = []
    for pl, E in (("x", el["H"]), ("y", el["V"])):
        R = ours[pl]; big = np.abs(E) > 0.3 * np.max(np.abs(E))
        rows.append((f"{pl}: max |отклик| elegant, м/рад", float(np.max(np.abs(E))), float(np.max(np.abs(R))),
                     float(np.max(np.abs(R - E)))))
        rows.append((f"{pl}: медианное относительное расхождение (крупные члены)", "", "",
                     float(np.median(np.abs((R - E)[big] / E[big])))))
    nsl, order = M["slices"]
    for r in trk:
        rows.append((f"аналитика ↔ трекинг (40 срезов, порядок 6), кик на s = {float(r['s']):.1f} м: орбита, мкм",
                     float(r["x_matrix"]) * 1e6, float(r["x_track"]) * 1e6,
                     abs(float(r["x_matrix"]) / float(r["x_track"]) - 1)))
    for s_k, z_an, z_ss, z_el in pairs:
        rows.append((f"трекинг ↔ трекинг ({nsl} срезов), кик на s = {s_k:.1f} м: орбита, мкм",
                     z_ss * 1e6, z_el * 1e6, abs(z_ss / z_el - 1)))
        rows.append((f"внутри кода: аналитика против трекинга ({nsl} срезов) на s = {s_k:.1f} м, мкм",
                     z_an * 1e6, z_ss * 1e6, abs(z_an / z_ss - 1)))
    _csv(os.path.join(d, "summary.csv"),
         ["показатель", "матрицы SynchrotronSim", "elegant / трекинг", "расхождение"], rows)

    # --- коррекция орбиты: разворачиваем часть квадруполей и правим орбиту корректорами
    rng = np.random.default_rng(SEED_MISALIGN)
    quads = [j for j, k in enumerate(ring.kinds) if k == "quad"]
    picked = rng.choice(quads, N_MISALIGNED, replace=False)
    yaws = rng.normal(0.0, SIGMA_YAW, N_MISALIGNED)
    pert = _ring().perturbed(yaws=[(int(j), float(a)) for j, a in zip(picked, yaws)])
    lmp = rm.line_matrices(pert)
    orbit0 = rr.station_orbit(pert, lmp, st, "x")
    Rx = rr.orbit_response(pert, lmp, st, "x", kinds)
    kicks, S = rr.svd_solve(Rx, orbit0, n_sv=min(KEEP_SV, Rx.shape[1]))
    orbit1 = orbit0 + Rx @ kicks
    corr = _correctors(pert, lmp)
    B, cur = rr.kick_to_currents(pert, corr[int(np.argmax(np.abs(kicks)))][0], float(np.max(np.abs(kicks))),
                                 design_mev=M["e_mev"])
    crows = [("СКО орбиты до коррекции, мкм", float(np.std(orbit0)) * 1e6, "", ""),
             ("СКО орбиты после коррекции, мкм", float(np.std(orbit1)) * 1e6, "", ""),
             ("максимум орбиты до / после, мкм", float(np.max(np.abs(orbit0))) * 1e6,
              float(np.max(np.abs(orbit1))) * 1e6, ""),
             ("СКО кика корректоров, мкрад", float(np.std(kicks)) * 1e6, "", ""),
             ("максимальный кик, мкрад", float(np.max(np.abs(kicks))) * 1e6, "", ""),
             ("физкоэффициент B при максимальном кике, Тл", B, "", ""),
             ("число обусловленности матрицы отклика", float(S[0] / S[-1]), "", ""),
             ("оставлено сингулярных чисел", min(KEEP_SV, len(S)), len(S), "")]
    for w, v in (cur.items() if isinstance(cur, dict) else []):
        crows.append((f"ток обмотки {w} при максимальном кике, А", v, "", ""))
    _csv(os.path.join(d, "correction.csv"), ["показатель", "значение", "после", ""], crows)

    # --- допуска
    tol, tw = rr.tolerances(ring, 1e-3, lm)
    trows = [(dev, t["n"], t["KL"], t["dnux"], t["dnuy"], t["beat_rms_x"] * 100, t["beat_rms_y"] * 100)
             for dev, t in sorted(tol.items(), key=lambda kv: -abs(kv[1]["dnux"]))]
    _csv(os.path.join(d, "tolerances.csv"),
         ["устройство", "число в кольце", "Σ|K1·L|, 1/м", "Δνx (согласованная ошибка 0,1 %)",
          "Δνy (согласованная)", "СКО Δβx/βx, % (случайная 0,1 %)", "СКО Δβy/βy, %"], trows)

    # --- рисунки
    s_b = np.array([x for _, x in st])
    corr0_s = float(_correctors(ring, lm)[0][1])
    fig, axs = plt.subplots(2, 2, figsize=(14, 8), sharex="col")
    for k, (pl, E) in enumerate((("x", el["H"]), ("y", el["V"]))):
        col = 0
        axs[k, col].plot(s_b, E[:, 0], "o", color=C_EL, ms=4, mfc="none", label="elegant")
        axs[k, col].plot(s_b, ours[pl][:, 0], "-", color=C_SS, lw=1.0, label="SynchrotronSim, матрицы")
        axs[k, col].set_ylabel(f"отклик {pl}, м/рад"); axs[k, col].grid(alpha=0.3)
        axs[k, col].set_title(f"первый корректор (s = {corr0_s:.1f} м), плоскость {pl}")
        axs[k, 1].plot(s_b, (ours[pl][:, 0] - E[:, 0]) * 1e6, color=C_SS, lw=0.9)
        axs[k, 1].set_ylabel("разность, мкм/рад"); axs[k, 1].grid(alpha=0.3)
        axs[k, 1].set_title("SynchrotronSim − elegant")
    axs[0, 0].legend(fontsize=9)
    for a in (axs[1, 0], axs[1, 1]):
        a.set_xlabel("s датчика, м")
    fig.suptitle(f"Отклик орбиты на кик корректора: {len(st)} датчиков")
    fig.tight_layout(); fig.savefig(os.path.join(d, "orbit_response.png"), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    axs[0].semilogy(np.arange(1, len(S) + 1), S, "o-", color=C_SS, ms=3)
    axs[0].axvline(KEEP_SV, color=C_EL, ls="--", label=f"оставлено {KEEP_SV}")
    axs[0].set_xlabel("номер сингулярного числа"); axs[0].set_ylabel("сингулярное число, м/рад")
    axs[0].grid(alpha=0.3, which="both"); axs[0].legend(); axs[0].set_title("спектр матрицы отклика")
    axs[1].plot(s_b, orbit0 * 1e6, color=C_EL, lw=0.9,
                label=f"до коррекции, СКО {np.std(orbit0) * 1e6:.2f} мкм")
    axs[1].plot(s_b, orbit1 * 1e6, color=C_SS, lw=0.9,
                label=f"после коррекции, СКО {np.std(orbit1) * 1e6:.4f} мкм")
    axs[1].set_xlabel("s датчика, м"); axs[1].set_ylabel("x, мкм"); axs[1].grid(alpha=0.3); axs[1].legend()
    axs[1].set_title(f"{N_MISALIGNED} квадруполей развёрнуты, СКО поворота {SIGMA_YAW * 1e6:.0f} мкрад")
    fig.suptitle("Обращение матрицы отклика по SVD и коррекция орбиты")
    fig.tight_layout(); fig.savefig(os.path.join(d, "svd_correction.png"), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    devs = [r[0] for r in trows]; ind = np.arange(len(devs))
    axs[0].barh(ind - 0.2, [abs(r[3]) for r in trows], 0.4, color=C_SS, label="Δνx")
    axs[0].barh(ind + 0.2, [abs(r[4]) for r in trows], 0.4, color=C_EL, label="Δνy")
    axs[0].set_yticks(ind); axs[0].set_yticklabels(devs, fontsize=8); axs[0].invert_yaxis()
    axs[0].set_xlabel("сдвиг частоты при согласованной ошибке 0,1 %"); axs[0].grid(alpha=0.3, axis="x")
    axs[0].legend()
    axs[1].barh(ind - 0.2, [r[5] for r in trows], 0.4, color=C_SS, label="Δβx/βx")
    axs[1].barh(ind + 0.2, [r[6] for r in trows], 0.4, color=C_EL, label="Δβy/βy")
    axs[1].set_yticks(ind); axs[1].set_yticklabels([]); axs[1].invert_yaxis()
    axs[1].set_xlabel("СКО бета-биения, % (случайная ошибка 0,1 %)"); axs[1].grid(alpha=0.3, axis="x")
    axs[1].legend()
    fig.suptitle("Чувствительность кольца к ошибке градиента по типам устройств")
    fig.tight_layout(); fig.savefig(os.path.join(d, "tolerances.png"), dpi=130); plt.close(fig)

    for r in rows:
        print(f"{r[0]:58s} {r[1]!s:>16s} {r[2]!s:>16s} {r[3]!s:>12s}")
    for r in crows:
        print(f"{r[0]:52s} {r[1]!s:>16s} {r[2]!s:>16s}")


# ============================================================================ частота от энергии
MAP_ORDER = 5                     # порядок аналитического отображения: точны ξ1…ξ4
FIT_HALF, FIT_DEG = 0.01, 7       # подгонка ν(δ) трекинга: |δ| ≤ 1 %, многочлен 7-й степени


def _map_chrom():
    """Аналитика SynchrotronSim высших порядков: хроматичность и уплотнение орбит по отображению оборота
    (ring_map — дифференциальная алгебра, точные стыки и тела элементов)."""
    import ring_map as rmap
    r = rmap.RingMap(_ring(), MAP_ORDER)
    Mp = r.one_turn(); ch = r.chromaticity(Mp)
    C, al = r.momentum_compaction(Mp, ch["orbit"])
    return dict(xi_x=ch["xi_x"], xi_y=ch["xi_y"], C=C, alphac=al, order=MAP_ORDER)


def exp_chrom(refresh=False):
    d = _dir("chromaticity")
    deltas = np.unique(np.round(np.concatenate([np.arange(-0.03, 0.03001, 0.0025),
                                                np.arange(-FIT_HALF, FIT_HALF + 1e-9, 0.001)]), 5))

    def run_ss():
        res = _pool([("tunes_point", (float(x),), {}, _env(), _ring_kwargs()) for x in deltas])
        return dict(delta=deltas, qx=np.array([r[0] for r in res]), qy=np.array([r[1] for r in res]),
                    orbit=np.array([r[2] for r in res]))

    def run_el():
        def one(i):
            return ra.tunes_point(_er(f"chr{i}_"), float(deltas[i]))
        res = _el_parallel(one, range(len(deltas)))
        return dict(delta=deltas, qx=np.array([r[0] for r in res]), qy=np.array([r[1] for r in res]),
                    orbit=np.array([r[2] for r in res]))

    def run_path():
        """Путь по замкнутой орбите вне импульса (трекинг SynchrotronSim) — для αc и αc2."""
        ring = _ring()
        out = []
        for x in deltas[np.abs(deltas) <= FIT_HALF + 1e-9]:
            co = ra.closed_orbit(ring, float(x))
            out.append((float(x), float(ring.track_elements([co])[0, -1, 4] - co[4])))
        return dict(path=np.array(out))

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el, refresh)
    mp = _cached(os.path.join(d, "ss_map.npz"), _map_chrom, refresh)
    pth = _cached(os.path.join(d, "ss_path.npz"), run_path, refresh)["path"]
    lin = np.load(os.path.join(_dir("linear"), "elegant_twiss.npz"), allow_pickle=True)
    par = ed.parameters(_stem("twiss_closed") + ".twi",
                        ["nux", "nuy", "dnux/dp", "dnuy/dp", "dnux/dp2", "dnuy/dp2", "dnux/dp3", "dnuy/dp3",
                         "alphac", "alphac2"])
    i0 = int(np.argmin(np.abs(deltas)))
    res = {}
    for code, dat in (("ss", ss), ("el", el)):
        for pl in ("x", "y"):
            q = _unwrap_frac(dat["q" + pl], i0)
            res[(code, pl)] = q + math.floor(float(lin["nu" + pl]))
    # аналитика elegant: ряд Тейлора ν0 + ξ1δ + ξ2δ²/2 + ξ3δ³/6
    dd = np.linspace(deltas[0], deltas[-1], 300)
    tay = {pl: par["nu" + pl] + par[f"dnu{pl}/dp"] * dd + par[f"dnu{pl}/dp2"] * dd ** 2 / 2
           + par[f"dnu{pl}/dp3"] * dd ** 3 / 6 for pl in ("x", "y")}
    # аналитика SynchrotronSim: ряд Тейлора по отображению оборота (до ξ4)
    tay_map = {pl: math.floor(float(lin["nu" + pl])) + sum(float(c) * dd ** k / math.factorial(k)
                                                          for k, c in enumerate(mp["xi_" + pl]))
               for pl in ("x", "y")}

    fig, ax = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for col, pl in enumerate(("x", "y")):
        a = ax[0, col]
        a.plot(dd * 100, tay[pl], "--", color=C_TW, lw=1.2, label="elegant twiss_output, default_order=3")
        a.plot(dd * 100, tay_map[pl], ":", color="k", lw=1.4,
               label=f"SynchrotronSim, отображение {int(mp['order'])}-го порядка")
        a.plot(deltas * 100, res[("el", pl)], "o", color=C_EL, ms=6, label="elegant, трекинг")
        a.plot(deltas * 100, res[("ss", pl)], "-", color=C_SS, lw=1.5, label="SynchrotronSim, трекинг")
        a.set_ylabel(f"ν{pl}"); a.grid(alpha=0.3); a.set_title(f"ν{pl}(δ)")
        b = ax[1, col]
        b.plot(deltas * 100, (res[("ss", pl)] - res[("el", pl)]) * 1e6, "o-", color=C_SS)
        b.set_ylabel(f"ν{pl}: SynchrotronSim − elegant, ×10⁻⁶"); b.grid(alpha=0.3)
        b.set_xlabel("δ = Δp/p, %")
    ax[0, 0].legend(fontsize=9)
    fig.suptitle("Бетатронная частота от энергии (замкнутая орбита вне импульса, матрица оборота)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "tune_vs_delta.png"), dpi=130); plt.close(fig)

    # трекинг: производные в нуле — подгонка многочленом высокой степени на |δ| ≤ 1 % (кубическая
    # подгонка смещена членами δ⁴ и выше и портит ξ2, ξ3 на проценты)
    fac = np.array([math.factorial(k) for k in range(4)])
    fits = {}
    for code in ("ss", "el"):
        for pl in ("x", "y"):
            q = res[(code, pl)]; ok = np.isfinite(q) & (np.abs(deltas) <= FIT_HALF + 1e-9)
            fits[(code, pl)] = np.polyfit(deltas[ok], q[ok], FIT_DEG)[::-1][:4] * fac
    rows = []
    for k, nm in enumerate(("ν0", "ξ1 = dν/dδ", "ξ2 = d²ν/dδ²", "ξ3 = d³ν/dδ³")):
        for pl in ("x", "y"):
            key = {0: "nu" + pl, 1: f"dnu{pl}/dp", 2: f"dnu{pl}/dp2", 3: f"dnu{pl}/dp3"}[k]
            an = float(mp["xi_" + pl][k]) + (math.floor(float(lin["nu" + pl])) if k == 0 else 0.0)
            rows.append((f"{nm} ({pl})", an, float(par[key]), float(fits[("ss", pl)][k]), float(fits[("el", pl)][k])))
    # уплотнение орбит: путь по замкнутой орбите C(δ) = C·(1 + αc·δ + αc2·δ² + …)
    pc = np.polyfit(pth[:, 0], pth[:, 1], FIT_DEG)[::-1]
    for k, nm in ((1, "αc"), (2, "αc2")):
        rows.append((nm, float(mp["alphac"][k - 1]), float(par["alphac" if k == 1 else "alphac2"]),
                     float(pc[k] / pc[0]), ""))
    _csv(os.path.join(d, "summary.csv"), ["коэффициент", "SynchrotronSim аналитика (отображение)",
                                          "elegant twiss_output", "SynchrotronSim трекинг",
                                          "elegant трекинг"], rows)
    tab = [(float(x), float(a), float(b), float(c), float(e)) for x, a, b, c, e in
           zip(deltas, res[("ss", "x")], res[("el", "x")], res[("ss", "y")], res[("el", "y")])]
    _csv(os.path.join(d, "tunes.csv"), ["δ", "νx SynchrotronSim", "νx elegant", "νy SynchrotronSim", "νy elegant"], tab)
    for r in rows:
        print(f"{r[0]:24s} SS-аналит {r[1]:+.7g}   elegant-twiss {r[2]:+.7g}   SS-трек {r[3]:+.7g}   "
              f"elegant-трек {r[4] if r[4] == '' else format(r[4], '+.7g')}")


# ============================================================================ ядро: C++ против Python-эталона
def exp_kernel(refresh=False):
    """Трекинг всех опытов идёт на продакшн-ядре C++ (ring_cpp.CppRing). Здесь оно сверяется с
    Python-эталоном ядра (ring_analysis.Ring поверх matrix_tracker.py) на решётке установки: оборот,
    выходы магнитов, пооборотный трекинг, ВЧ с излучением, кик корректора, разворот линзы, срезы изгибов."""
    import time, copy
    import matrix_tracker as mt
    d = _dir("kernel")
    n, o = M["slices"]
    kw = _ring_kwargs()
    os.environ.update({"SS_SECTOR_NSL": str(n), "SS_SECTOR_ORDER": str(o)})
    try:
        P = ra.Ring(M["lattice"], **kw)
        volt, harm = _rf()
        rf = dict(volt_mv=volt, harmonic=harm, phase_deg=150.0)
        P2 = ra.Ring(M["lattice"], sr=True, rf=rf, **kw)
    finally:
        os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
    C = _ring(); C2 = _ring(sr=True, rf=rf)
    co = ra.closed_orbit(C)
    rng = np.random.default_rng(1)
    st = co + np.vstack([np.zeros(6), rng.normal(0, [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 1e-3], (6, 6))])

    def dmax(a, b):
        return float(np.nanmax(np.abs(np.asarray(a) - np.asarray(b))))

    def timed(fn):
        t = time.perf_counter(); v = fn(); return v, time.perf_counter() - t

    rows = []
    a, tp = timed(lambda: P.track(st, 1)); b, tc = timed(lambda: C.track(st, 1))
    rows.append(("1 оборот, 7 частиц: max |Δz|", dmax(a, b), tp, tc))
    a, tp = timed(lambda: P.track_elements(st)); b, tc = timed(lambda: C.track_elements(st))
    rows.append((f"выходы всех {len(P.mags)} магнитов за оборот: max |Δz|", dmax(a, b), tp, tc))
    a, tp = timed(lambda: P.track_tbt(st[:3], 30)); b, tc = timed(lambda: C.track_tbt(st[:3], 30))
    rows.append(("30 оборотов, 3 частицы: max |Δz|", dmax(a, b), tp, tc))
    a, tp = timed(lambda: P2.track(st, 5)); b, tc = timed(lambda: C2.track(st, 5))
    rows.append(("ВЧ и классическое излучение, 5 оборотов: max |Δz|", dmax(a, b), tp, tc))
    corr = [j for j, k in enumerate(P.kinds) if k in _kick_kinds()]
    if corr:
        j = corr[0]
        a = ra.kicked(P, j, 1e-4).track(st, 1); ck = C.with_kick(j, 1e-4); b = ck.track(st, 1); ck.close()
        rows.append((f"кик 0,1 мрад в пробном корректоре {P.names[j]}: max |Δz|", dmax(a, b), "", ""))
    q = next(j for j, k in enumerate(P.kinds) if k == "quad" and j > 5)
    pert = copy.copy(P); pert.mags = list(P.mags)
    pert.mags[q] = ra.yaw_element(copy.copy(P.mags[q]), 2e-4)
    cy = C.with_yaw(q, 2e-4); b = cy.track(st, 1); cy.close()
    rows.append((f"разворот линзы {P.names[q]} на 0,2 мрад: max |Δz|", dmax(pert.track(st, 1), b), "", ""))
    rec = np.full((len(st), C.n_pieces + 1, 6), np.nan)
    pieces = []
    for j, m in enumerate(P.mags):
        pieces += ([(c, i == 0, j) for i, c in enumerate(ra._bend_slices(m))] if m.sector else [(m, True, j)])
    for p_, s0 in enumerate(st):
        cur = tuple(map(float, s0)); rec[p_, 0] = cur; prev = P.plane
        for qq, (pc, junction, j) in enumerate(pieces):
            if junction:
                cur = mt.modified_drift(cur, prev, P.mags[j]); prev = P.mags[j]
            cur = pc.body(cur); rec[p_, qq + 1] = cur
    rows.append((f"после каждого из {C.n_pieces} срезов изгибов (для диффузии излучения): max |Δz|",
                 dmax(rec, C.track_pieces(st)), "", ""))
    _csv(os.path.join(d, "summary.csv"), ["проверка", "расхождение C++ − Python", "время Python, с",
                                          "время C++, с"], rows)
    for r in rows:
        print(f"{r[0]:70s} {r[1]:.3e}  {r[2] if r[2] == '' else f'{r[2]:.3f}'}  {r[3] if r[3] == '' else f'{r[3]:.4f}'}")


# ============================================================================ порядок кикового интегратора
def exp_integrator(refresh=False):
    d = _dir("integrator")
    combos = [(n, o) for o in (2, 4, 6) for n in (1, 2, 3, 5, 10, 20, 40)]

    def run_ss():
        tasks = [("tunes_point", (0.0,), {}, {"SS_SECTOR_NSL": str(n), "SS_SECTOR_ORDER": str(o)},
                  _ring_kwargs()) for n, o in combos]
        res = _pool(tasks)
        return dict(combos=np.array(combos), qx=np.array([r[0] for r in res]), qy=np.array([r[1] for r in res]))

    def run_el():
        def one(i):
            n, o = combos[i]
            return ra.tunes_point(_er(f"int{i}_", lattice_kwargs=dict(integrator=(n, o))), 0.0)
        res = _el_parallel(one, range(len(combos)))
        return dict(combos=np.array(combos), qx=np.array([r[0] for r in res]), qy=np.array([r[1] for r in res]))

    def run_def():
        """Число срезов ПО УМОЛЧАНИЮ: elegant даёт линзе один шаг интегратора, ядро — пять."""
        q_el = ra.tunes_point(_er("intdef", default_slices=True), 0.0)
        q_ss = _pool([("tunes_point", (0.0,), {}, None, _ring_kwargs())], processes=1)[0]
        return dict(el=np.array(q_el[:2], float), ss=np.array(q_ss[:2], float))

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el, refresh)
    dflt = _cached(os.path.join(d, "default.npz"), run_def, refresh)
    par = ed.parameters(_stem("twiss_closed") + ".twi", ["nux", "nuy"])
    ref = {"x": par["nux"] % 1, "y": par["nuy"] % 1}

    fig, ax = plt.subplots(1, 3, figsize=(16, 5))
    marks = {2: "s", 4: "o", 6: "^"}
    for col, pl in enumerate(("x", "y")):
        a = ax[col]
        for o in (2, 4, 6):
            sel = [i for i, (n, oo) in enumerate(combos) if oo == o]
            ns = np.array([combos[i][0] for i in sel])
            a.loglog(ns, _fdist(el["q" + pl][sel], ref[pl]), marks[o] + "--", color=C_EL, ms=9, mfc="none",
                     label=f"elegant, порядок {o}")
            a.loglog(ns, _fdist(ss["q" + pl][sel], ref[pl]), marks[o] + "-", color=C_SS, ms=4,
                     label=f"SynchrotronSim, порядок {o}")
        a.set_xlabel("срезов на элемент N"); a.set_ylabel(f"|ν{pl}(трекинг) − ν{pl}(twiss_output)|")
        a.grid(alpha=0.3, which="both"); a.set_title(f"плоскость {pl}")
    ax[0].legend(fontsize=8)
    lab = [f"N={n}, п.{o}" for n, o in combos]
    ax[2].semilogy(range(len(combos)), _fdist(ss["qx"], el["qx"]) + 1e-12, "o-", color=C_SS, label="νx")
    ax[2].semilogy(range(len(combos)), _fdist(ss["qy"], el["qy"]) + 1e-12, "s--", color=C_SS, label="νy")
    ax[2].set_xticks(range(len(combos))); ax[2].set_xticklabels(lab, fontsize=7, rotation=60)
    ax[2].set_title("|SynchrotronSim − elegant|, те же N и порядок"); ax[2].grid(alpha=0.3, which="both")
    ax[2].legend()
    fig.suptitle("Сходимость частот по числу срезов и порядку симплектического интегратора")
    fig.tight_layout(); fig.savefig(os.path.join(d, "integrator_convergence.png"), dpi=130); plt.close(fig)
    rows = [("по умолчанию", float(dflt["ss"][0]), float(dflt["el"][0]),
             float(dflt["ss"][1]), float(dflt["el"][1]))]
    rows += [(f"N={n}, порядок {o}", float(ss["qx"][i]), float(el["qx"][i]), float(ss["qy"][i]),
              float(el["qy"][i])) for i, (n, o) in enumerate(combos)]
    rows.append(("twiss_output (аналитика)", ref["x"], ref["x"], ref["y"], ref["y"]))
    _csv(os.path.join(d, "tunes.csv"), ["интегратор", "νx SynchrotronSim", "νx elegant", "νy SynchrotronSim",
                                        "νy elegant"], rows)
    for r in rows:
        print(f"{r[0]:26s} νx SS {r[1]:.7f} el {r[2]:.7f}   νy SS {r[3]:.7f} el {r[4]:.7f}")


# ============================================================================ пооборотный трекинг: портреты, спектр, амплитуды
TBT_TURNS = 1024
SMALL = 1e-5                                   # малое смещение другой плоскости, м (чтобы её частота была видна)


def _amps():
    """Амплитуды сканов по x и по y, мм — свои у каждой установки."""
    return tuple(M["amp_x_mm"]), tuple(M["amp_y_mm"])


def _tbt_states(orbit):
    """Начальные частицы: скан по x (y мал), затем скан по y (x мал); смещения от замкнутой орбиты."""
    ax_mm, ay_mm = _amps()
    out = []; labels = []
    for a in ax_mm:
        st = np.array(orbit, float); st[0] += a * 1e-3; st[2] += SMALL; out.append(st); labels.append(("x", a))
    for a in ay_mm:
        st = np.array(orbit, float); st[0] += SMALL; st[2] += a * 1e-3; out.append(st); labels.append(("y", a))
    return out, labels


_DRIVING_KEYS = ("dnux/dJx", "dnux/dJy", "dnuy/dJy", "h21000", "h30000", "h10110", "h10020", "h10200")


MAP_NF_ORDER = 4                  # порядок отображения для нормальной формы (нужны члены 3-й степени)


def _detuning():
    """Сдвиг частот с действием ∂ν/∂J — аналитика SynchrotronSim: нормальная форма отображения оборота
    (ring_map.normal_form). Действие J = (q² + (αq + βp)²)/(2β), как у elegant."""
    import ring_map as rmap
    nf = rmap.normal_form(_ring(), MAP_NF_ORDER)
    return {k: float(v) for k, v in nf.items()}


def _driving_terms():
    """Аналитика elegant для нелинейной части: сдвиг частот с амплитудой ∂ν/∂J и резонансные члены
    (twiss_output, compute_driving_terms). Считается по мультиполям решётки — кинематические члены
    (точные дрейфы, края) в эти формулы не входят, и у решётки без секступолей ∂ν/∂J выходит нулевым."""
    lte = _lattice()
    stem = _stem("driving")
    ed.run(stem + ".ele", _setup(3, lte) +
           f'&twiss_output filename="%s.twi", matched=1, radiation_integrals=1, compute_driving_terms=1, &end\n'
           % stem)
    return dict(ed.parameters(stem + ".twi", list(_DRIVING_KEYS)))


def exp_tbt(refresh=False):
    d = _dir("tbt")
    lin = _dir("linear")
    ss_lin = np.load(os.path.join(lin, "ss.npz"), allow_pickle=True)
    el_lin = np.load(os.path.join(lin, "elegant_track.npz"), allow_pickle=True)

    def run_ss():
        states, _ = _tbt_states(ss_lin["closed_orbit"])
        return dict(tbt=_tbt(states, TBT_TURNS))

    def run_el():
        states, _ = _tbt_states(el_lin["closed_orbit"])
        chunks = [states[i:i + 2] for i in range(0, len(states), 2)]

        def one(i):
            return _er(f"tbt{i}_").track_tbt(chunks[i], TBT_TURNS)
        return dict(tbt=np.concatenate(_el_parallel(one, range(len(chunks))), axis=0))

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)["tbt"]
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el, refresh)["tbt"]
    _, labels = _tbt_states(np.zeros(6))
    bx, ax_, by, ay = (float(ss_lin["betax"][-1]), float(ss_lin["alphax"][-1]),
                       float(ss_lin["betay"][-1]), float(ss_lin["alphay"][-1]))

    # --- частоты (NAFF) и выживание
    rows = []; tunes = {}
    for i, (pl, a) in enumerate(labels):
        rec = [f"{pl} = {a:g} мм"]
        for code, dat, co in (("ss", ss, ss_lin["closed_orbit"]), ("el", el, el_lin["closed_orbit"])):
            tr = dat[i]; n_ok = int(np.sum(np.all(np.isfinite(tr), axis=1))) - 1
            if n_ok < TBT_TURNS:
                q = (np.nan, np.nan)
            else:
                q = ra.tbt_tunes(tr[:TBT_TURNS], (bx, by), (ax_, ay), orbit=co)
            tunes[(code, i)] = q
            rec += [n_ok, q[0], q[1]]
        rows.append(rec)
    _csv(os.path.join(d, "tunes_vs_amplitude.csv"),
         ["частица", "оборотов SynchrotronSim", "νx SynchrotronSim", "νy SynchrotronSim",
          "оборотов elegant", "νx elegant", "νy elegant"], rows)

    # --- аналитика elegant: сдвиг частот с амплитудой и резонансные члены по формулам мультиполей
    drv = _cached(os.path.join(d, "driving_terms.npz"), _driving_terms, refresh)
    _csv(os.path.join(d, "driving_terms.csv"), ["величина", "elegant twiss_output (аналитика)"],
         [(k, float(drv[k])) for k in _DRIVING_KEYS])

    # --- сдвиг частот с действием: аналитика SynchrotronSim (нормальная форма) против трекинга обеих
    nf = _cached(os.path.join(d, "ss_detuning.npz"), _detuning, refresh)
    act = {"x": lambda a: a * a * (1 + ax_ * ax_) / (2 * bx), "y": lambda a: a * a * (1 + ay * ay) / (2 * by)}
    slope = {}
    for code in ("ss", "el"):
        for pl in ("x", "y"):
            pts = [(act[pl](a_ * 1e-3), *tunes[(code, i)]) for i, (p_, a_) in enumerate(labels)
                   if p_ == pl and np.isfinite(tunes[(code, i)][0])]
            pts.sort()
            if len(pts) < 2:
                continue
            js = np.array([p_[0] for p_ in pts])
            # наклон берётся у САМЫХ МАЛЫХ амплитуд: дальше вступают члены J² и дробная частота может
            # перевалить через период (у СИЛА сдвиг за несколько миллиметров больше половины периода)
            use = js <= 4.0 * js[0] + 1e-30
            for k, out_pl in ((1, "x"), (2, "y")):
                q = np.unwrap(np.array([p_[k] for p_ in pts]) * 2 * np.pi) / (2 * np.pi)
                slope[(code, out_pl, pl)] = float(np.polyfit(js[use], q[use], 1)[0])
    rows_nf = []
    for out_pl, act_pl in (("x", "x"), ("x", "y"), ("y", "y")):
        key = f"dnu{out_pl}/dJ{act_pl}"
        rows_nf.append((f"∂ν{out_pl}/∂J{act_pl}, 1/м", float(nf[key]),
                        slope.get(("ss", out_pl, act_pl), float("nan")),
                        slope.get(("el", out_pl, act_pl), float("nan")), float(drv[key])))
    rows_nf.append(("∂νy/∂Jx, 1/м (должен совпасть с ∂νx/∂Jy)", float(nf["dnuy/dJx"]),
                    slope.get(("ss", "y", "x"), float("nan")), slope.get(("el", "y", "x"), float("nan")), ""))
    _csv(os.path.join(d, "detuning.csv"),
         ["величина", "SynchrotronSim аналитика (нормальная форма)", "SynchrotronSim трекинг (NAFF)",
          "elegant трекинг (NAFF)", "elegant twiss_output"], rows_nf)
    for r in rows_nf:
        print(f"{r[0]:44s} аналитика {r[1]:+.6g}   трекинг SS {r[2]:+.6g}   трекинг elegant {r[3]:+.6g}   "
              f"elegant-twiss {r[4] if r[4] == '' else format(r[4], '+.6g')}")

    # --- рисунок: фазовые портреты
    ax_mm, ay_mm = _amps()
    pick3 = lambda v: (v[len(v) // 4], v[len(v) // 2], v[-2 if len(v) > 2 else -1])
    show = [labels.index((pl, a)) for pl, vs in (("x", ax_mm), ("y", ay_mm)) for a in pick3(vs)]
    fig, axs = plt.subplots(2, 3, figsize=(15, 9))
    for k, i in enumerate(show):
        pl, a = labels[i]; c = 0 if pl == "x" else 2
        axx = axs[k // 3, k % 3]
        for code, dat, col, mk, ms in (("elegant", el, C_EL, "o", 5), ("SynchrotronSim", ss, C_SS, ".", 2)):
            tr = dat[i]
            axx.plot(tr[:, c] * 1e3, tr[:, c + 1] * 1e3, mk, color=col, ms=ms, mfc="none" if mk == "o" else col,
                     alpha=0.6 if mk == "o" else 1.0, label=code, ls="none")
        axx.set_xlabel(f"{pl}, мм"); axx.set_ylabel(f"p{pl}, мрад"); axx.grid(alpha=0.3)
        n_ss = int(np.sum(np.all(np.isfinite(ss[i]), axis=1))) - 1
        n_el = int(np.sum(np.all(np.isfinite(el[i]), axis=1))) - 1
        axx.set_title(f"старт {pl} = {a:g} мм; прошло оборотов: SS {n_ss}, elegant {n_el}", fontsize=10)
    axs[0, 0].legend(markerscale=3)
    fig.suptitle(f"Фазовые портреты в точке наблюдения (s = 0), {TBT_TURNS} оборотов")
    fig.tight_layout(); fig.savefig(os.path.join(d, "phase_portraits.png"), dpi=130); plt.close(fig)

    # --- рисунок: расхождение координат по оборотам
    fig, axs = plt.subplots(1, 2, figsize=(14, 5))
    cmap = plt.get_cmap("viridis")
    for i, (pl, a) in enumerate(labels):
        c = 0 if pl == "x" else 2
        dv = np.abs(ss[i][:, c] - el[i][:, c]) * 1e6
        axx = axs[0 if pl == "x" else 1]
        n = len(ax_mm) if pl == "x" else len(ay_mm)
        j = [aa for p2, aa in labels if p2 == pl].index(a)
        axx.semilogy(np.arange(len(dv)), np.maximum(dv, 1e-6), color=cmap(j / max(1, n - 1)), lw=0.8,
                     label=f"{pl}0 = {a:g} мм")
    for axx, pl in zip(axs, ("x", "y")):
        axx.set_xlabel("оборот"); axx.set_ylabel(f"|{pl}(SynchrotronSim) − {pl}(elegant)|, мкм")
        axx.grid(alpha=0.3, which="both"); axx.legend(fontsize=8)
    fig.suptitle("Расхождение пооборотных координат SynchrotronSim и elegant")
    fig.tight_layout(); fig.savefig(os.path.join(d, "tbt_difference.png"), dpi=130); plt.close(fig)

    # --- рисунок: спектр
    fig, axs = plt.subplots(2, 1, figsize=(12, 8))
    for axx, pl, amps in ((axs[0], "x", pick3(ax_mm)), (axs[1], "y", pick3(ay_mm))):
        c = 0 if pl == "x" else 2
        for j, a in enumerate(amps):
            i = labels.index((pl, a))
            for dat, col, ls in ((el, C_EL, "-"), (ss, C_SS, "--")):
                tr = dat[i][:TBT_TURNS, c]
                if not np.all(np.isfinite(tr)):
                    continue
                w = np.hanning(len(tr)); sp = np.abs(np.fft.rfft((tr - tr.mean()) * w)); f = np.fft.rfftfreq(len(tr))
                axx.semilogy(f, sp / sp.max() * 10 ** (-3 * j), ls, color=col, lw=0.8,
                             label=("elegant" if col == C_EL else "SynchrotronSim") + f", {pl}0={a:g} мм")
        axx.set_xlabel("частота, доли оборотной"); axx.set_ylabel(f"|БПФ {pl}| (кривые сдвинуты)")
        axx.grid(alpha=0.3); axx.legend(fontsize=7, ncol=3)
    fig.suptitle("Спектр пооборотных колебаний")
    fig.tight_layout(); fig.savefig(os.path.join(d, "spectrum.png"), dpi=130); plt.close(fig)

    # --- рисунок: частоты от амплитуды
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    nu_lin = (float(ss_lin["nux"]) % 1, float(ss_lin["nuy"]) % 1)
    for col_i, pl in enumerate(("x", "y")):
        idx = [i for i, (p2, a) in enumerate(labels) if p2 == pl]
        amps = np.array([0.0] + [labels[i][1] for i in idx])
        for row_i, q in enumerate(("νx", "νy")):
            axx = axs[row_i, col_i]
            for code, col, mk in (("el", C_EL, "o"), ("ss", C_SS, "s")):
                # сдвиг частоты от линейной, развёрнутый по непрерывности (дробная часть проходит через 0 и 1)
                v = _unwrap_frac(np.array([nu_lin[row_i]] + [tunes[(code, i)][row_i] for i in idx]), 0) - nu_lin[row_i]
                axx.plot(amps, v, mk + ("--" if code == "el" else "-"), color=col, mfc="none" if code == "el" else col,
                         ms=8 if code == "el" else 4, label="elegant" if code == "el" else "SynchrotronSim")
            axx.set_xlabel(f"начальная амплитуда {pl}, мм"); axx.set_ylabel(f"Δ{q} от линейной частоты")
            axx.grid(alpha=0.3)
    axs[0, 0].legend()
    fig.suptitle(f"Сдвиг частот с амплитудой (NAFF по {TBT_TURNS} оборотам; точка 0 — линейная оптика)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "tunes_vs_amplitude.png"), dpi=130); plt.close(fig)
    for r in rows:
        print(r)


# ============================================================================ продольное движение: ВЧ, излучение, затухание
LONG_TURNS = 2048
DAMP_TURNS = 4096
DELTA_AMPS = (0.5e-3, 1e-3, 2e-3, 4e-3, 6e-3)


def _rf():
    """ВЧ установки: (напряжение, МВ; кратность)."""
    return float(M["rf"]["volt_mv"]), int(M["rf"]["harmonic"])


def _rad_integrals(rf=None):
    """Аналитика излучения SynchrotronSim: интегралы I1…I5 по описанию решётки и всё, что из них следует
    (ring_matrix.radiation_integrals). Трекинг не участвует."""
    import ring_matrix as rm
    ring = _ring()
    if rf:
        ring.rf = dict(rf)
    try:
        return rm.radiation_integrals(ring)
    finally:
        ring.rf = None


def _u0_by_tracking(ring):
    """Потеря энергии за оборот осевой частицей (излучение включено, ВЧ нет), МэВ."""
    return float(-ring.track([np.zeros(6)])[0, 5] * M["e_mev"])


def _long_point(ring):
    """Неподвижная точка 6D, матрица оборота, затухание и частоты — один алгоритм для обоих трекеров."""
    co = ra.closed_orbit(ring, dims=6, eps=1e-7, tol=1e-12)
    _, M = ra.one_turn_matrix(ring, orbit=co)
    return dict(orbit=co, M=M)


def exp_long(refresh=False):
    d = _dir("longitudinal")
    lin = np.load(os.path.join(_dir("linear"), "ss.npz"), allow_pickle=True)
    volt, harm = _rf()
    rf0 = dict(volt_mv=volt, harmonic=harm, phase_deg=180.0)                # без излучения: φs = 180°

    # --- U0 по трекингу (излучение, без ВЧ): у каждого кода свой, фаза ВЧ — от U0 elegant twiss_output
    def run_u0():
        u_ss = _u0_by_tracking(_ring(sr=True))
        u_el = _u0_by_tracking(_er("u0", sr=True))
        tw = ed.parameters(_stem("twiss_closed") + ".twi",
                           ["U0", "taux", "tauy", "taudelta", "Jx", "Jy", "Jdelta", "Sdelta0", "ex0", "alphac"])
        return dict(u_ss=u_ss, u_el=u_el, **{k: v for k, v in tw.items()})

    u = _cached(os.path.join(d, "u0.npz"), run_u0, refresh)
    phase_sr = 180.0 - math.degrees(math.asin(float(u["U0"]) / volt))
    rf_sr = dict(volt_mv=volt, harmonic=harm, phase_deg=phase_sr)

    # --- матрицы оборота: (а) ВЧ без излучения; (б) ВЧ + излучение
    def run_mats():
        out = {}
        for tag, kw in (("rf", dict(rf=rf0)), ("rfsr", dict(rf=rf_sr, sr=True))):
            r = _long_point(_ring(**kw)); out[f"ss_{tag}_M"] = r["M"]; out[f"ss_{tag}_orbit"] = r["orbit"]
            r = _long_point(_er(f"lm{tag}", **kw)); out[f"el_{tag}_M"] = r["M"]
            out[f"el_{tag}_orbit"] = r["orbit"]
        return out

    mats = _cached(os.path.join(d, "matrices.npz"), run_mats, refresh)

    # --- пооборотные синхротронные колебания (ВЧ без излучения)
    def tbt_states(orbit):
        out = []
        for a in DELTA_AMPS:
            st = np.array(orbit, float); st[5] += a; out.append(st)
        return out

    def run_tbt_ss():
        return dict(tbt=_tbt(tbt_states(mats["ss_rf_orbit"]), LONG_TURNS, rf=rf0))

    def run_tbt_el():
        st = tbt_states(mats["el_rf_orbit"])

        def one(i):
            return _er(f"lt{i}_", rf=rf0).track_tbt([st[i]], LONG_TURNS)
        return dict(tbt=np.concatenate(_el_parallel(one, range(len(st))), axis=0))

    tss = _cached(os.path.join(d, "tbt_ss.npz"), run_tbt_ss, refresh)["tbt"]
    tel = _cached(os.path.join(d, "tbt_elegant.npz"), run_tbt_el, refresh)["tbt"]

    # --- затухание по спаду амплитуды (излучение + ВЧ): без разностных производных
    def damp_states(orbit):
        out = []
        for k, a in ((0, 3e-4), (2, 3e-4), (5, 5e-4)):
            st = np.array(orbit, float); st[k] += a; out.append(st)
        return out

    def run_damp_ss():
        return dict(tbt=_tbt(damp_states(mats["ss_rfsr_orbit"]), DAMP_TURNS, rf=rf_sr, sr=True))

    def run_damp_el():
        st = damp_states(mats["el_rfsr_orbit"])

        def one(i):
            return _er(f"dmp{i}_", rf=rf_sr, sr=True).track_tbt([st[i]], DAMP_TURNS)
        return dict(tbt=np.concatenate(_el_parallel(one, range(len(st))), axis=0))

    dss = _cached(os.path.join(d, "damp_ss.npz"), run_damp_ss, refresh)["tbt"]
    delg = _cached(os.path.join(d, "damp_elegant.npz"), run_damp_el, refresh)["tbt"]

    C = float(lin["s"][-1]); T0 = C / ed.C_LIGHT
    rows = []
    res = {}
    for tag in ("rf", "rfsr"):
        for code in ("ss", "el"):
            res[(code, tag)] = ra.damping_from_matrix(mats[f"{code}_{tag}_M"], C)
    rows.append(("U0 по трекингу за оборот, МэВ", float(u["u_ss"]), float(u["u_el"]), float(u["U0"])))
    rows.append(("Qs (ВЧ, без излучения), матрица", res[("ss", "rf")]["s"]["tune"], res[("el", "rf")]["s"]["tune"],
                 math.sqrt(harm * float(u["alphac"]) * volt / (2 * math.pi * M["e_mev"]))))
    rows.append(("Qs (ВЧ + излучение), матрица", res[("ss", "rfsr")]["s"]["tune"], res[("el", "rfsr")]["s"]["tune"],
                 math.sqrt(harm * float(u["alphac"]) * volt * abs(math.cos(math.radians(phase_sr)))
                           / (2 * math.pi * M["e_mev"]))))
    for pl, key in (("x", "taux"), ("y", "tauy"), ("s", "taudelta")):
        rows.append((f"τ{pl}, мс", res[("ss", "rfsr")][pl]["tau_turns"] * T0 * 1e3,
                     res[("el", "rfsr")][pl]["tau_turns"] * T0 * 1e3, float(u[key]) * 1e3))
    tau = {c: [res[(c, "rfsr")][p]["tau_turns"] for p in ("x", "y", "s")] for c in ("ss", "el")}
    for i, (pl, key) in enumerate((("x", "Jx"), ("y", "Jy"), ("s", "Jdelta"))):
        J = {c: (4.0 / sum(1.0 / t for t in tau[c])) / tau[c][i] for c in ("ss", "el")}
        rows.append((f"J{pl} (по τ, ΣJ=4)", J["ss"], J["el"], float(u[key])))
    bx, ax_, by, ay = (float(lin["betax"][-1]), float(lin["alphax"][-1]), float(lin["betay"][-1]),
                       float(lin["alphay"][-1]))
    env = {}
    for code, dat, co in (("ss", dss, mats["ss_rfsr_orbit"]), ("el", delg, mats["el_rfsr_orbit"])):
        env[(code, "x")] = ra.damping_by_envelope(dat[0], 0, bx, ax_, co)
        env[(code, "y")] = ra.damping_by_envelope(dat[1], 2, by, ay, co)
        env[(code, "s")] = ra.damping_by_envelope(dat[2], 4, orbit=co, window=1024, step=256, search=(0.0, 0.05))
    for pl, key in (("x", "taux"), ("y", "tauy"), ("s", "taudelta")):
        rows.append((f"τ{pl} по спаду амплитуды, мс", env[("ss", pl)][0] * T0 * 1e3, env[("el", pl)][0] * T0 * 1e3,
                     float(u[key]) * 1e3))
    rows.append(("δ синхронной частицы у плоскости", float(mats["ss_rfsr_orbit"][5]), float(mats["el_rfsr_orbit"][5]), 0.0))
    rows.append(("t синхронной частицы у плоскости, мм", float(mats["ss_rfsr_orbit"][4]) * 1e3,
                 float(mats["el_rfsr_orbit"][4]) * 1e3, 0.0))

    # частота по пооборотным данным (NAFF δ) и амплитудная зависимость
    qs_t = {}
    for code, dat, co in (("ss", tss, mats["ss_rf_orbit"]), ("el", tel, mats["el_rf_orbit"])):
        for i, a in enumerate(DELTA_AMPS):
            tr = dat[i]
            if not np.all(np.isfinite(tr)):
                qs_t[(code, i)] = np.nan; continue
            dd = tr[:, 5] - co[5]; tt = tr[:, 4] - co[4]
            scale = np.std(dd) / max(np.std(tt), 1e-30)
            qs_t[(code, i)] = float(ra.naff(dd + 1j * tt * scale, search=(0.0, 0.05))[0][0])
            qs_t[(code, i)] = min(qs_t[(code, i)], 1 - qs_t[(code, i)])
    for i, a in enumerate(DELTA_AMPS):
        rows.append((f"Qs по трекингу (NAFF), δ0 = {a:g}", qs_t[("ss", i)], qs_t[("el", i)], np.nan))
    ri = _rad_integrals(rf_sr)
    qs0 = math.sqrt(harm * ri["alphac"] * volt / (2 * math.pi * M["e_mev"]))
    ss_an = {"U0 по трекингу за оборот, МэВ": ri["U0"], "Qs (ВЧ, без излучения), матрица": qs0,
             "Qs (ВЧ + излучение), матрица": ri["Qs"]}
    for pl, key in (("x", "taux"), ("y", "tauy"), ("s", "taudelta")):
        ss_an[f"τ{pl}, мс"] = ri[key] * 1e3
        ss_an[f"τ{pl} по спаду амплитуды, мс"] = ri[key] * 1e3
    for pl, key in (("x", "Jx"), ("y", "Jy"), ("s", "Jdelta")):
        ss_an[f"J{pl} (по τ, ΣJ=4)"] = ri[key]
    rows = [tuple(r) + (ss_an.get(r[0], ""),) for r in rows]
    _csv(os.path.join(d, "summary.csv"),
         ["показатель", "SynchrotronSim трекинг", "elegant трекинг", "elegant аналитика / формула",
          "SynchrotronSim аналитика (интегралы излучения)"], rows)
    for r in rows:
        print(f"{r[0]:42s} SS {r[1]:.7g}   elegant {r[2]:.7g}   аналитика-elegant {r[3]:.7g}   "
              f"аналитика-SS {r[4] if r[4] == '' else format(r[4], '.7g')}")

    # --- рисунки: продольные портреты и Qs от амплитуды
    fig, axs = plt.subplots(1, 3, figsize=(17, 5.5))
    cmap = plt.get_cmap("viridis")
    for i, a in enumerate(DELTA_AMPS):
        for dat, col, mk, ms in ((tel, C_EL, "o", 4), (tss, C_SS, ".", 1.5)):
            axs[0].plot(dat[i][:, 4] * 1e3, dat[i][:, 5] * 1e3, mk, color=col, ms=ms,
                        mfc="none" if mk == "o" else col, alpha=0.5 if mk == "o" else 1)
    axs[0].plot([], [], "o", color=C_EL, mfc="none", label="elegant"); axs[0].plot([], [], ".", color=C_SS, label="SynchrotronSim")
    axs[0].set_xlabel("запаздывание t, мм"); axs[0].set_ylabel("δ, 10⁻³"); axs[0].grid(alpha=0.3); axs[0].legend()
    axs[0].set_title(f"продольные портреты, {LONG_TURNS} оборотов (ВЧ {volt:g} МВ, h = {harm})")
    for i, a in enumerate(DELTA_AMPS):
        dv = np.abs(tss[i][:, 5] - tel[i][:, 5])
        axs[1].semilogy(np.maximum(dv, 1e-12), color=cmap(i / (len(DELTA_AMPS) - 1)), lw=0.8, label=f"δ0 = {a:g}")
    axs[1].set_xlabel("оборот"); axs[1].set_ylabel("|δ(SynchrotronSim) − δ(elegant)|"); axs[1].grid(alpha=0.3, which="both")
    axs[1].legend(fontsize=8); axs[1].set_title("расхождение по оборотам")
    axs[2].plot(np.array(DELTA_AMPS) * 1e3, [qs_t[("el", i)] for i in range(len(DELTA_AMPS))], "o--", color=C_EL,
                mfc="none", ms=9, label="elegant")
    axs[2].plot(np.array(DELTA_AMPS) * 1e3, [qs_t[("ss", i)] for i in range(len(DELTA_AMPS))], "s-", color=C_SS, ms=4,
                label="SynchrotronSim")
    axs[2].set_xlabel("начальное δ0, 10⁻³"); axs[2].set_ylabel("Qs"); axs[2].grid(alpha=0.3); axs[2].legend()
    axs[2].set_title("синхротронная частота от амплитуды")
    fig.tight_layout(); fig.savefig(os.path.join(d, "synchrotron.png"), dpi=130); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.8))
    for a, pl, title in ((axs[0], "x", "x"), (axs[1], "y", "y"), (axs[2], "s", "продольная")):
        for code, col, mk in (("el", C_EL, "o"), ("ss", C_SS, ".")):
            tau, cen, amp = env[(code, pl)]
            a.semilogy(cen, amp / amp[0], mk, color=col, mfc="none" if mk == "o" else col, ms=6 if mk == "o" else 5,
                       label=f"{'elegant' if code == 'el' else 'SynchrotronSim'}: τ = {tau * T0 * 1e3:.4f} мс")
        a.set_xlabel("оборот"); a.set_ylabel("амплитуда / начальная"); a.grid(alpha=0.3, which="both")
        a.set_title(f"плоскость {title}"); a.legend(fontsize=9)
    fig.suptitle(f"Радиационное затухание по спаду амплитуды ({DAMP_TURNS} оборотов, излучение + ВЧ)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "damping.png"), dpi=130); plt.close(fig)


# ============================================================================ портреты с излучением
def _rad_cfg():
    """Длина прогона пробных частиц и ансамбль квантового излучения установки:
    (обороты, (частиц, обороты, выход на равновесие, шаг выборки) или None — без ансамбля)."""
    return int(M.get("rad_turns", 8192)), M.get("rad_ensemble")


def _rad_rf():
    """ВЧ с синхронной фазой под потерю энергии установки (U0 из интегралов излучения elegant)."""
    volt, harm = _rf()
    u0 = float(ed.parameters(_stem("twiss_closed") + ".twi", ["U0"])["U0"])
    return dict(volt_mv=volt, harmonic=harm, phase_deg=180.0 - math.degrees(math.asin(u0 / volt)))


def _rad_states(orbit):
    """Три пробные частицы: отклонение по x, по y и по δ от синхронной точки."""
    ax, ay, ad = M["rad_amp"]
    out = []
    for k, a in ((0, ax * 1e-3), (2, ay * 1e-3), (5, ad)):
        st = np.array(orbit, float); st[k] += a; out.append(st)
    return out, (("x", "px", 0, 1, "x, мм", "px, мрад"),
                 ("y", "py", 2, 3, "y, мм", "py, мрад"),
                 ("t", "δ", 4, 5, "t, мм", "δ, 10⁻³"))


def exp_rad(refresh=False):
    """Фазовые портреты пробных частиц с излучением: классическое затухание и квантовое возбуждение.

    Классическое излучение (SYNCH_RAD) детерминировано — портреты обеих программ должны лечь друг на
    друга. Квантовое (ISR) разыгрывается независимыми генераторами, поэтому траектории сравнивать
    поточечно нельзя: там сверяется размер равновесного облака."""
    d = _dir("radiation")
    rf = _rad_rf()
    RAD_TURNS, ENS = _rad_cfg()
    ENS_N, ENS_TURNS, ENS_SKIP, ENS_STRIDE = ENS or (0, 0, 0, 1)

    def orbit_of(ring):
        return ra.closed_orbit(ring, dims=6, eps=1e-7, tol=1e-12)

    def run_ss(isr):
        co = orbit_of(_ring(rf=rf, sr=True))       # синхронная точка — по классическому излучению
        st, _ = _rad_states(co)
        return dict(orbit=co, tbt=_tbt(st, RAD_TURNS, rf=rf, **({"isr": True} if isr else {"sr": True})))

    def run_el(isr):
        kw = dict(rf=rf, **({"isr": True} if isr else {"sr": True}))
        co = orbit_of(_er("rado" + ("i" if isr else "s"), rf=rf, sr=True))
        st, _ = _rad_states(co)

        def one(i):                                  # у каждого запуска свой генератор квантового шума
            return _er(f"rad{'i' if isr else 's'}{i}_", run_setup_extra=f"random_number_seed={1000003 + 7919 * i},",
                       **kw).track_tbt([st[i]], RAD_TURNS)
        return dict(orbit=co, tbt=np.concatenate(_el_parallel(one, range(len(st))), axis=0))

    def run_ens_ss():
        """Ансамбль с синхронной точки: ENS_N частиц, квантовое излучение, выборка после выхода на равновесие."""
        co = orbit_of(_ring(rf=rf, sr=True))
        tb = _tbt(np.tile(co, (ENS_N, 1)), ENS_TURNS, rf=rf, isr=True)
        return dict(orbit=co, snap=tb[:, ENS_SKIP::ENS_STRIDE, :])

    def run_ens_el():
        co = orbit_of(_er("radoe", rf=rf, sr=True))
        k = 8                                         # частиц на запуск; запуски — параллельно, у каждого свой seed

        def one(i):
            tb = _er(f"rade{i}_", rf=rf, isr=True, run_setup_extra=f"random_number_seed={2000003 + 7919 * i},"
                     ).track_tbt(np.tile(co, (k, 1)), ENS_TURNS)
            return tb[:, ENS_SKIP::ENS_STRIDE, :]
        return dict(orbit=co, snap=np.concatenate(_el_parallel(one, range(ENS_N // k)), axis=0))

    dat = {}
    for tag, isr in (("sr", False), ("isr", True)):
        dat[("ss", tag)] = _cached(os.path.join(d, f"ss_{tag}.npz"), lambda i=isr: run_ss(i), refresh)
        dat[("el", tag)] = _cached(os.path.join(d, f"elegant_{tag}.npz"), lambda i=isr: run_el(i), refresh)
    ens = ({"ss": _cached(os.path.join(d, "ss_ensemble.npz"), run_ens_ss, refresh),
            "el": _cached(os.path.join(d, "elegant_ensemble.npz"), run_ens_el, refresh)} if ENS else None)

    _, planes = _rad_states(np.zeros(6))
    sc = (1e3, 1e3, 1e3, 1e3, 1e3, 1e3)             # x,y,t → мм; px,py → мрад; δ → 10⁻³
    # равновесные размеры из опыта «эмиттанс» — ими проверяется размер квантового облака
    em = {}
    ep = os.path.join(_dir("emittance"), "summary.csv")
    if os.path.exists(ep):
        for line in open(ep, encoding="utf-8-sig").read().splitlines()[1:]:
            c = line.split(";")                                  # SynchrotronSim, elegant-трек, moments, twiss
            em[c[0]] = [float(v.replace(",", ".")) if v.strip() else float("nan") for v in c[1:]]
    def _sig(key, scale, col):
        return em.get(key, [0.0] * 4)[col] * scale
    sigma = {0: _sig("σx у плоскости, мкм", 1e-3, 2),            # мм, moments_output
             4: _sig("σz, мм", 1.0, 2),                          # мм
             5: _sig("σδ, 10⁻⁴", 0.1, 2)}                        # в единицах 10⁻³
    sigma_codes = {0: ("σx у плоскости, мкм", 1e-3), 4: ("σz, мм", 1.0), 5: ("σδ, 10⁻⁴", 0.1)}
    # функции Твисса в точке наблюдения: по ним амплитуда колебания считается инвариантом Куранта–Снайдера
    # (hypot по осям графика зависел бы от фазы, и разница частот программ давала бы ложный разброс)
    lin = np.load(os.path.join(_dir("linear"), "ss.npz"), allow_pickle=True)
    lo = os.path.join(_dir("longitudinal"), "summary.csv")
    qs = 0.0
    if os.path.exists(lo):
        for line in open(lo, encoding="utf-8-sig").read().splitlines():
            if line.startswith("Qs (ВЧ + излучение)"):
                qs = float(line.split(";")[1].replace(",", "."))
    beta_s = (float(lin["alphac"]) * float(lin["s"][-1]) / (2 * np.pi * qs)) if qs else 0.0
    TW = {0: (float(lin["betax"][-1]), float(lin["alphax"][-1])),
          2: (float(lin["betay"][-1]), float(lin["alphay"][-1])),
          4: (beta_s, 0.0)}

    def amplitude(code, tag, col, ia, ib):
        """Амплитуда колебания по оборотам, в единицах координаты (мм): sqrt(q² + (αq + βp)²)."""
        tr = dat[(code, tag)]["tbt"][col]; co = dat[(code, "sr")]["orbit"]
        ok = np.all(np.isfinite(tr), axis=1)
        q = tr[ok, ia] - co[ia]; pp = tr[ok, ib] - co[ib]
        b, al = TW[ia]
        return np.hypot(q, al * q + b * pp) * 1e3

    # --- рисунок: два ряда (классическое / квантовое) × три плоскости
    fig, axs = plt.subplots(2, 3, figsize=(16, 9))
    for row, (tag, title) in enumerate((("sr", "классическое излучение (детерминировано)"),
                                        ("isr", "квантовое излучение (независимые розыгрыши)"))):
        for col, (na, nb, ia, ib, la, lb) in enumerate(planes):
            a = axs[row, col]
            for code, lab, color, ms in (("el", "elegant", C_EL, 1.2), ("ss", "SynchrotronSim", C_SS, 1.2)):
                tr = dat[(code, tag)]["tbt"][col]
                co = dat[(code, "sr")]["orbit"]
                ok = np.all(np.isfinite(tr), axis=1)
                a.plot((tr[ok, ia] - co[ia]) * sc[ia], (tr[ok, ib] - co[ib]) * sc[ib], ".",
                       color=color, ms=ms, alpha=0.6, label=lab)
            if tag == "isr":                       # ожидаемое равновесие по моментам
                for idx, fn in ((ia, a.axvline), (ib, a.axhline)):
                    if sigma.get(idx):
                        for sgn in (-1, 1):
                            fn(sgn * sigma[idx], color="0.35", ls="--", lw=1.0,
                               label="±σ по moments_output" if (sgn > 0 and idx == ia) else None)
            a.set_xlabel(la); a.set_ylabel(lb); a.grid(alpha=0.3)
            a.set_title(f"{na}–{nb}, {title}", fontsize=10)
            if tag == "isr" and col == 0:
                a.legend(markerscale=8, fontsize=8)
    axs[0, 0].legend(markerscale=8, fontsize=9)
    fig.suptitle(f"Фазовые портреты пробных частиц с излучением, {RAD_TURNS} оборотов "
                 f"(смещения от синхронной точки)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "phase_portraits.png"), dpi=130); plt.close(fig)

    # --- рисунок: спад амплитуды (классическое) и выход на равновесие (квантовое)
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    for col, (na, nb, ia, ib, la, lb) in enumerate(planes):
        a = axs[col]
        for tag, ls in (("sr", "-"), ("isr", ":")):
            for code, color in (("el", C_EL), ("ss", C_SS)):
                env = amplitude(code, tag, col, ia, ib)
                w = max(1, len(env) // 512)
                env = env[: len(env) // w * w].reshape(-1, w).mean(axis=1)
                a.semilogy(np.arange(len(env)) * w, env, ls, color=color, lw=1.0,
                           label=f"{'elegant' if code == 'el' else 'SynchrotronSim'}, "
                                 f"{'класс.' if tag == 'sr' else 'квант.'}")
        a.set_xlabel("оборот"); a.set_ylabel(f"амплитуда {na}, мм")
        a.grid(alpha=0.3, which="both")
    axs[0].legend(fontsize=8)
    fig.suptitle("Огибающая амплитуды: классическое излучение уводит к нулю, квантовое — к равновесию")
    fig.tight_layout(); fig.savefig(os.path.join(d, "envelope.png"), dpi=130); plt.close(fig)

    # --- равновесное облако ансамбля: СКО по частицам и оборотам после выхода на равновесие
    ENS_KEYS = (0, 4, 5) if ENS else ()
    par_rad = ed.parameters(_stem("twiss_closed") + ".twi", ["taux", "tauy", "taudelta"])
    T0 = float(lin["s"][-1]) / 299792458.0
    tau_turns = {0: float(par_rad["taux"]) / T0, 2: float(par_rad["tauy"]) / T0,
                 4: float(par_rad["taudelta"]) / T0, 5: float(par_rad["taudelta"]) / T0}
    span = max(1, ENS_TURNS - ENS_SKIP)
    # независимых выборок: частицы × (длина выборки / время жизни действия τ/2); ошибка СКО ≈ 1/√(2n)
    ens_err = {k: 1.0 / math.sqrt(2.0 * ENS_N * max(1.0, span / (tau_turns[k] / 2.0))) for k in ENS_KEYS}

    def ens_cloud(code, i):
        sn = ens[code]["snap"]; co = ens[code]["orbit"]
        v = (sn[..., i] - co[i]).ravel()
        return v[np.isfinite(v)] * sc[i]

    def ens_sigma(code, i):
        v = ens_cloud(code, i)
        return float(np.sqrt(np.mean((v - v.mean()) ** 2)))

    fig, axs = plt.subplots(2, 2, figsize=(13, 10)) if ENS else (None, None)
    for a, (ia, ib, la, lb) in (zip(axs[0], ((0, 1, "x, мм", "px, мрад"), (4, 5, "t, мм", "δ, 10⁻³"))) if ENS
                                else ()):
        for code, lab, color in (("el", "elegant", C_EL), ("ss", "SynchrotronSim", C_SS)):
            a.plot(ens_cloud(code, ia)[::7], ens_cloud(code, ib)[::7], ".", ms=0.8, alpha=0.35, color=color,
                   label=lab)
        a.set_xlabel(la); a.set_ylabel(lb); a.grid(alpha=0.3)
        a.set_title(f"равновесное облако: {ENS_N} частиц, обороты {ENS_SKIP}…{ENS_TURNS}", fontsize=10)
        a.legend(markerscale=12, fontsize=9)
    for a, (i, la) in (zip(axs[1], ((0, "x, мм"), (5, "δ, 10⁻³"))) if ENS else ()):
        rng_ = max(np.max(np.abs(ens_cloud("ss", i))), np.max(np.abs(ens_cloud("el", i))))
        bins = np.linspace(-rng_, rng_, 81)
        for code, lab, color in (("el", "elegant", C_EL), ("ss", "SynchrotronSim", C_SS)):
            a.hist(ens_cloud(code, i), bins=bins, density=True, histtype="step", lw=1.6, color=color,
                   label=f"{lab}: σ = {ens_sigma(code, i):.4g}")
        if sigma.get(i):
            g = np.linspace(-rng_, rng_, 400)
            a.plot(g, np.exp(-g ** 2 / (2 * sigma[i] ** 2)) / (sigma[i] * math.sqrt(2 * math.pi)), "--",
                   color="0.3", lw=1.2, label=f"гаусс, σ = {sigma[i]:.4g} (moments_output elegant)")
        a.set_xlabel(la); a.set_ylabel("плотность"); a.grid(alpha=0.3); a.legend(fontsize=8)
    if ENS:
        fig.suptitle("Квантовое излучение: равновесное облако ансамбля частиц у точки наблюдения")
        fig.tight_layout(); fig.savefig(os.path.join(d, "equilibrium_clouds.png"), dpi=130); plt.close(fig)

    # --- таблица
    rows = []
    for col, (na, nb, ia, ib, la, lb) in enumerate(planes):
        e_ss, e_el = amplitude("ss", "sr", col, ia, ib), amplitude("el", "sr", col, ia, ib)
        n = min(len(e_ss), len(e_el))
        big = e_el[:n] > 0.1 * e_el[0]                # где сигнал ещё велик: ниже 10 % начальной амплитуды
        k1 = min(n - 1, int(round(tau_turns.get(ia, tau_turns[0]))))   # через одно время затухания
        rows.append((f"{na}, классическое: амплитуда через одно время затухания ({k1} об.), начальная "
                     f"{e_ss[0]:.3g} мм".replace(".", ",").replace(",)", ".)"),
                     float(e_ss[k1]), float(e_el[k1]), abs(float(e_ss[k1]) / float(e_el[k1]) - 1)))
        rows.append((f"{na}, классическое: наибольшее расхождение огибающих (амплитуда > 10 % начальной)",
                     "", "", float(np.max(np.abs(e_ss[:n][big] - e_el[:n][big]) / e_el[:n][big]))))
        for i in ((ia, 5) if ia == 4 else (ia,)):
            if i in ENS_KEYS:
                v = [ens_sigma(code, i) for code in ("ss", "el")]
                nm = {0: "σx, мм", 4: "σz, мм", 5: "σδ, 10⁻³"}[i]
                rows.append((f"{nm}: квантовое излучение, ансамбль из {ENS_N} частиц "
                             f"(статистика ±{100 * ens_err[i]:.1f} %)".replace(".", ","), v[0], v[1],
                             abs(v[0] / v[1] - 1)))
            if i in sigma_codes and em:
                key, scale = sigma_codes[i]
                nm = {0: "σx, мм", 4: "σz, мм", 5: "σδ, 10⁻³"}[i]
                rows.append((f"{nm}: расчёт равновесия (опыт «эмиттанс», метод Чао)",
                             em[key][0] * scale, em[key][1] * scale,
                             abs(em[key][0] / em[key][1] - 1) if em[key][1] else ""))
    _csv(os.path.join(d, "summary.csv"),
         ["показатель", "SynchrotronSim", "elegant", "расхождение"], rows)
    for r in rows:
        print(f"{r[0]:62s} {r[1]!s:>16s} {r[2]!s:>16s} {r[3]!s:>12s}")


# ============================================================================ равновесный эмиттанс
EMIT_ISR_PARTICLES = 16000


def exp_emit(refresh=False):
    d = _dir("emittance")
    tw = ed.parameters(_stem("twiss_closed") + ".twi", ["U0", "ex0", "Sdelta0", "Jx", "alphac"])
    volt, harm = _rf()
    phase = 180.0 - math.degrees(math.asin(float(tw["U0"]) / volt))
    rf = dict(volt_mv=volt, harmonic=harm, phase_deg=phase)

    def run_ss():
        ring = _ring(sr=True, rf=rf)
        co = ra.closed_orbit(ring, dims=6, eps=1e-7, tol=1e-12)
        M, B = ra.isr_diffusion_matrix(ring, co)
        return dict(M=M, B=B, orbit=co)

    def run_el_track():
        er = _er("em", sr=True, rf=rf)
        co = ra.closed_orbit(er, dims=6, eps=1e-7, tol=1e-12)
        _, M = ra.one_turn_matrix(er, orbit=co)
        n_chunk = 8; per = EMIT_ISR_PARTICLES // n_chunk

        def one(i):
            return _er(f"emi{i}_", isr=True, rf=rf).track([co] * per)
        fin = np.concatenate(_el_parallel(one, range(n_chunk)), axis=0)
        fin = fin[np.all(np.isfinite(fin), axis=1)]
        return dict(M=M, B=np.cov(fin.T), orbit=co, n=len(fin))

    def run_moments():
        stem = _stem("moments")
        lte = os.path.basename(ed.write_lattice(
            stem + ".lte", csbend_extra="SYNCH_RAD=1,ISR=1", rf=rf,
            **(dict(source=M["source"], beamline=M["beamline"], use_beamline=M["use_beamline"])
               if M["source"] else dict(config=M["config"]))))
        ed.run(stem + ".ele", _setup(2, lte) +
               '&moments_output filename="%s.mom", matched=1, equilibrium=1, radiation=1, &end\n')
        c = ed.columns(stem + ".mom", ["s1", "s5", "s6", "ex", "ez"], page=1)
        return dict(sx=c["s1"][-1], sz=c["s5"][-1], sd=c["s6"][-1], ex=c["ex"][-1], ez=c["ez"][-1])

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el_track, refresh)
    mo = _cached(os.path.join(d, "elegant_moments.npz"), run_moments, refresh)
    S_ss, em_ss = ra.equilibrium_moments(ss["M"], ss["B"])
    S_el, em_el = ra.equilibrium_moments(el["M"], el["B"])
    n = int(el["n"]); stat = math.sqrt(2.0 / n)
    import ring_matrix as rm
    ri = _rad_integrals(rf)                       # аналитика SynchrotronSim: интегралы излучения
    tws = rm.twiss(_ring())                       # β и η у точки наблюдения — для размера пучка
    sx_an = math.sqrt(float(tws["betax"][-1]) * ri["ex0"] + (float(tws["etax"][-1]) * ri["Sdelta0"]) ** 2)
    rows = [
        ("εx, пм", em_ss[0] * 1e12, em_el[0] * 1e12, float(mo["ex"]) * 1e12, float(tw["ex0"]) * 1e12,
         ri["ex0"] * 1e12),
        ("σδ, 10⁻⁴", math.sqrt(S_ss[5, 5]) * 1e4, math.sqrt(S_el[5, 5]) * 1e4, float(mo["sd"]) * 1e4,
         float(tw["Sdelta0"]) * 1e4, ri["Sdelta0"] * 1e4),
        ("σz, мм", math.sqrt(S_ss[4, 4]) * 1e3, math.sqrt(S_el[4, 4]) * 1e3, float(mo["sz"]) * 1e3, np.nan,
         ri["Sz0"] * 1e3),
        ("σx у плоскости, мкм", math.sqrt(S_ss[0, 0]) * 1e6, math.sqrt(S_el[0, 0]) * 1e6, float(mo["sx"]) * 1e6,
         np.nan, sx_an * 1e6),
        ("εz, мкм", em_ss[2] * 1e6, em_el[2] * 1e6, float(mo["ez"]) * 1e6, np.nan,
         ri["Sz0"] * ri["Sdelta0"] * 1e6),
    ]
    _csv(os.path.join(d, "summary.csv"),
         ["показатель", "SynchrotronSim (Чао, диффузия по срезам)",
          f"elegant (Чао, диффузия трекингом {n} частиц, ±{stat * 100:.1f}%)".replace(".", ","),
          "elegant moments_output", "elegant twiss_output (интегралы излучения)",
          "SynchrotronSim (интегралы излучения)"], rows)
    for r in rows:
        print(f"{r[0]:22s} SS {r[1]:.5g}   elegant-трек {r[2]:.5g}   moments_output {r[3]:.5g}   twiss {r[4]:.5g}")
    print(f"статистика диффузии elegant: {n} частиц, ±{stat * 100:.1f}%")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.8))
    labels = ["SynchrotronSim\n(Чао)", "elegant\n(Чао, трекинг)", "elegant\nmoments_output",
              "elegant\ntwiss_output", "SynchrotronSim\n(интегралы)"]
    for a, k, unit in ((ax[0], 0, "εx, пм"), (ax[1], 1, "σδ, 10⁻⁴")):
        vals = rows[k][1:]
        cols = [C_SS, C_EL, C_EL, C_TW, C_AN]
        a.bar(range(len(vals)), vals, color=cols[:len(vals)], alpha=0.85)
        if k == 0:
            a.errorbar([1], [vals[1]], yerr=[vals[1] * stat], color="k", capsize=6)
        for i, v in enumerate(vals):
            if np.isfinite(v):
                a.text(i, v, f"{v:.4g}", ha="center", va="bottom", fontsize=9)
        a.set_xticks(range(len(vals))); a.set_xticklabels(labels[:len(vals)], fontsize=8)
        a.set_ylabel(unit); a.grid(alpha=0.3, axis="y")
        lo = np.nanmin(vals); hi = np.nanmax(vals); a.set_ylim(lo - (hi - lo) * 3 - 0.01 * lo, hi + (hi - lo) * 2 + 0.01 * hi)
    fig.suptitle(f"Равновесный эмиттанс и энергетический разброс (излучение + ВЧ {volt:g} МВ)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "emittance.png"), dpi=130); plt.close(fig)


# ============================================================================ динамическая апертура, карта частот, резонансы
FMA_TURNS_HALF = 256


def _loss_limit():
    """Порог потери в точке наблюдения, м: |x| или |y| больше — частица потеряна (своё у установки)."""
    return float(M.get("loss_limit_m", 0.03))


def _fma_grid():
    """Сетка начальных смещений апертуры установки: (x, мм), (y, мм)."""
    a, b, n = M["fma_x_mm"]
    return np.linspace(a, b, int(n)), np.array(M["fma_y_mm"], float)


def _fma_states(orbit):
    out = []
    FMA_X_MM, FMA_Y_MM = _fma_grid()
    for y in FMA_Y_MM:
        for x in FMA_X_MM:
            st = np.array(orbit, float); st[0] += x * 1e-3; st[2] += y * 1e-3; out.append(st)
    return out


def exp_fma(refresh=False):
    d = _dir("dynamic_aperture")
    lin = np.load(os.path.join(_dir("linear"), "ss.npz"), allow_pickle=True)
    lin_el = np.load(os.path.join(_dir("linear"), "elegant_track.npz"), allow_pickle=True)
    FMA_X_MM, FMA_Y_MM = _fma_grid()
    nt = 2 * FMA_TURNS_HALF

    def run_ss():
        return dict(tbt=_tbt(_fma_states(lin["closed_orbit"]), nt, limit=_loss_limit()))

    def run_el():
        st = _fma_states(lin_el["closed_orbit"])
        chunks = [st[i:i + 24] for i in range(0, len(st), 24)]

        def one(i):
            return _er(f"fma{i}_").track_tbt(chunks[i], nt)
        return dict(tbt=np.concatenate(_el_parallel(one, range(len(chunks))), axis=0))

    def run_builtin():
        _lattice()
        stem = _stem("fma_builtin")
        ed.run(stem + ".ele", _setup(2) +
               f'&run_control n_passes={nt}, &end\n'
               f'&frequency_map output="%s.fma", xmin={FMA_X_MM[0] * 1e-3:g}, xmax={FMA_X_MM[-1] * 1e-3:g}, '
               f'ymin={FMA_Y_MM[0] * 1e-3:g}, ymax={FMA_Y_MM[-1] * 1e-3:g}, nx={len(FMA_X_MM)}, ny=13, '
               'include_changes=1, verbosity=0, &end\n')
        return ed.columns(stem + ".fma", ["x", "y", "nux", "nuy", "diffusionRate"])

    ss = _cached(os.path.join(d, "ss.npz"), run_ss, refresh)["tbt"]
    el = _cached(os.path.join(d, "elegant_track.npz"), run_el, refresh)["tbt"]
    try:
        bi = _cached(os.path.join(d, "elegant_builtin.npz"), run_builtin, refresh)
    except Exception as e:                                  # встроенный инструмент — справочно
        print("frequency_map elegant:", e); bi = None

    nx, ny = len(FMA_X_MM), len(FMA_Y_MM)
    beta = (float(lin["betax"][-1]), float(lin["betay"][-1])); alpha = (float(lin["alphax"][-1]), float(lin["alphay"][-1]))
    res = {}
    for code, dat, co in (("ss", ss, lin["closed_orbit"]), ("el", el, lin_el["closed_orbit"])):
        surv = ra.survived_turns(dat, _loss_limit())
        nux, nuy, dif = ra.frequency_map(dat, beta, alpha, co, FMA_TURNS_HALF)
        res[code] = dict(surv=surv.reshape(ny, nx), nux=nux, nuy=nuy, dif=dif)

    # --- выживание и граница апертуры
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    X, Y = np.meshgrid(FMA_X_MM, FMA_Y_MM)
    for k, (code, title) in enumerate((("ss", "SynchrotronSim"), ("el", "elegant"))):
        im = axs[k].pcolormesh(X, Y, res[code]["surv"], shading="nearest", cmap="viridis", vmin=0, vmax=nt)
        axs[k].set_title(f"{title}: оборотов до потери (из {nt})"); axs[k].set_xlabel("x, мм"); axs[k].set_ylabel("y, мм")
        fig.colorbar(im, ax=axs[k])
    diffm = res["ss"]["surv"] - res["el"]["surv"]
    im = axs[2].pcolormesh(X, Y, diffm, shading="nearest", cmap="coolwarm", vmin=-nt, vmax=nt)
    axs[2].set_title("разность: SynchrotronSim − elegant, оборотов"); axs[2].set_xlabel("x, мм"); axs[2].set_ylabel("y, мм")
    fig.colorbar(im, ax=axs[2])
    fig.suptitle(f"Динамическая апертура в точке наблюдения (δ = 0; потеря — |x| или |y| > {_loss_limit() * 1e3:g} мм)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "dynamic_aperture.png"), dpi=130); plt.close(fig)

    # --- карта частот на сетке резонансов
    # частицы, дожившие до конца: у погибших частоты хаотичны и засыпают всю сетку
    box = (0.0, 1.0, 0.0, 1.0)
    lines = ra.resonance_lines(4, box)
    wp = (float(lin["nux"]) % 1, float(lin["nuy"]) % 1)
    fig, axs = plt.subplots(1, 2, figsize=(15, 7))
    for k, (code, title) in enumerate((("ss", "SynchrotronSim"), ("el", "elegant"))):
        a = axs[k]
        for m, n, pp, pts in lines:
            order = abs(m) + n
            a.plot([pts[0][0], pts[1][0]], [pts[0][1], pts[1][1]], "-" if n % 2 == 0 else "--",
                   color="0.45", lw=2.4 / order, alpha=0.7, zorder=1)
        alive = res[code]["surv"].ravel() == nt
        sc = a.scatter(res[code]["nux"][alive], res[code]["nuy"][alive], c=res[code]["dif"][alive], cmap="jet",
                       vmin=-10, vmax=-2, s=28, zorder=3, edgecolors="k", linewidths=0.2)
        a.plot(*wp, "*", color="k", ms=14, zorder=4, label="рабочая точка (линейная оптика)")
        a.set_xlim(0, 1); a.set_ylim(0, 1); a.set_aspect("equal")
        a.set_xlabel("νx, дробная часть"); a.set_ylabel("νy, дробная часть")
        a.set_title(f"{title}: {int(alive.sum())} частиц дожили до {nt} оборотов")
        a.legend(loc="upper right", fontsize=9)
        fig.colorbar(sc, ax=a, label="диффузия log10|Δν|", shrink=0.8)
        # у слабонелинейного кольца весь след стягивается в точку — показываем его врезкой
        qx, qy = res[code]["nux"][alive], res[code]["nuy"][alive]
        span = max(np.ptp(qx), np.ptp(qy)) if len(qx) else 0.0
        if 0 < span < 0.02:
            ins = a.inset_axes([0.06, 0.06, 0.34, 0.34])
            ins.scatter(qx, qy, c=res[code]["dif"][alive], cmap="jet", vmin=-10, vmax=-2, s=16,
                        edgecolors="k", linewidths=0.2)
            ins.plot(*wp, "*", color="k", ms=10)
            m = 0.15 * span
            ins.set_xlim(min(qx.min(), wp[0]) - m, max(qx.max(), wp[0]) + m)
            ins.set_ylim(min(qy.min(), wp[1]) - m, max(qy.max(), wp[1]) + m)
            ins.ticklabel_format(useOffset=True, style="plain")
            ins.tick_params(labelsize=7); ins.set_title("след в масштабе", fontsize=8)
            a.indicate_inset_zoom(ins, edgecolor="0.3")
    fig.suptitle("Карта частот на сетке резонансов до 4-го порядка (сплошные — нормальные, пунктир — косые)")
    fig.tight_layout(); fig.savefig(os.path.join(d, "frequency_map.png"), dpi=130); plt.close(fig)

    # --- диффузия в координатах и разность частот кодов
    fig, axs = plt.subplots(1, 3, figsize=(17, 5))
    for k, (code, title) in enumerate((("ss", "SynchrotronSim"), ("el", "elegant"))):
        im = axs[k].pcolormesh(X, Y, res[code]["dif"].reshape(ny, nx), shading="nearest", cmap="jet", vmin=-10, vmax=-2)
        axs[k].set_title(f"{title}: диффузия log10|Δν|"); axs[k].set_xlabel("x, мм"); axs[k].set_ylabel("y, мм")
        fig.colorbar(im, ax=axs[k])
    dq = np.hypot((res["ss"]["nux"] - res["el"]["nux"] + 0.5) % 1 - 0.5, (res["ss"]["nuy"] - res["el"]["nuy"] + 0.5) % 1 - 0.5)
    im = axs[2].pcolormesh(X, Y, np.log10(np.maximum(dq, 1e-12)).reshape(ny, nx), shading="nearest", cmap="magma",
                           vmin=-9, vmax=-2)
    axs[2].set_title("log10 |ν(SynchrotronSim) − ν(elegant)|"); axs[2].set_xlabel("x, мм"); axs[2].set_ylabel("y, мм")
    fig.colorbar(im, ax=axs[2])
    fig.tight_layout(); fig.savefig(os.path.join(d, "diffusion_xy.png"), dpi=130); plt.close(fig)

    rows = []
    for j, y in enumerate(FMA_Y_MM):
        for i, x in enumerate(FMA_X_MM):
            q = j * nx + i
            rows.append((float(x), float(y), int(res["ss"]["surv"][j, i]), int(res["el"]["surv"][j, i]),
                         res["ss"]["nux"][q], res["el"]["nux"][q], res["ss"]["nuy"][q], res["el"]["nuy"][q],
                         res["ss"]["dif"][q], res["el"]["dif"][q]))
    _csv(os.path.join(d, "grid.csv"), ["x, мм", "y, мм", "оборотов SS", "оборотов elegant", "νx SS", "νx elegant",
                                       "νy SS", "νy elegant", "диффузия SS", "диффузия elegant"], rows)
    agree = int(np.sum((res["ss"]["surv"] == nt) == (res["el"]["surv"] == nt)))
    print(f"совпадение «выжила/потеряна»: {agree} из {nx * ny}")
    print("медиана |Δν| по выжившим:", float(np.nanmedian(dq)), " максимум:", float(np.nanmax(dq)))


def main():
    ap = argparse.ArgumentParser(description="опыты сверки SynchrotronSim с elegant")
    ap.add_argument("experiment", choices=["kernel", "linear", "matrices", "analytic", "response", "chrom",
                                           "integrator", "tbt", "long", "rad", "emit", "fma"])
    ap.add_argument("--machine", default="sila", choices=sorted(ed.MACHINES),
                    help="установка сверки (описана в elegant_diag.MACHINES)")
    ap.add_argument("--refresh", action="store_true", help="пересчитать, игнорируя кэш")
    a = ap.parse_args()
    use_machine(a.machine)
    print(f"установка: {M['label']} ({a.machine}), решётка {M['lattice']} БД, {M['e_mev']:g} МэВ")
    {"kernel": exp_kernel, "linear": exp_linear, "matrices": exp_matrices, "analytic": exp_analytic, "response": exp_response,
     "chrom": exp_chrom, "integrator": exp_integrator,
     "tbt": exp_tbt, "long": exp_long, "rad": exp_rad, "emit": exp_emit,
     "fma": exp_fma}[a.experiment](a.refresh)


if __name__ == "__main__":
    main()
