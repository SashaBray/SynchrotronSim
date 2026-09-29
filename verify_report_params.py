# -*- coding: utf-8 -*-
"""Отчёт «Сверка SynchrotronSim с elegant по параметрам» — Word, глава на каждый рассчитываемый параметр.

Формат:
  1. Установка и настройки программ, при которых шло сравнение.
  2… По главе на параметр: физический смысл → чем считается (модуль комплекса) → результаты (графики,
     таблицы) → расхождения отдельной таблицей в конце главы.

Данные — те же опыты verify_elegant.py (кэш .npz и CSV в workspace/results/verify_elegant/[<установка>/]).
Свои графики глав (орбита, β, η) строятся здесь же из кэша, в подпапку params_report.

  py -3.14 verify_report_params.py --machine siberia   → workspace/results/verify_elegant_params_siberia.docx
"""
import os, json, math, argparse

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import verify_report as vr
import elegant_diag as ed

C_LIGHT = 299792458.0
C_SS, C_EL, C_TW, C_AN = "#1f5fbf", "#d9480f", "#7a7a7a", "#2f9e44"
DIFF_HEADER = ["пара инструментов", "показатель", "расхождение", "чем объясняется"]


# ============================================================================ данные
def res(*p):
    return os.path.join(vr.RES, *p)


def npz(*p):
    f = res(*p)
    return dict(np.load(f, allow_pickle=True)) if os.path.exists(f) else None


def have(*names):
    return all(vr.exists(n) for n in names)


def fnum(v):
    return float(np.asarray(v))


def figdir():
    d = res("params_report"); os.makedirs(d, exist_ok=True); return d


# ============================================================================ оформление глав
class Doc(vr.Report):
    """Отчёт с нумерованными главами и одинаковыми подзаголовками внутри главы."""

    def __init__(self):
        super().__init__()
        self.n = 0

    def chapter(self, title):
        self.n += 1
        self.h1(f"{self.n}. {title}")

    def meaning(self):
        self.h2("Физический смысл")

    def tools(self, rows, experiment):
        """Чем считается: таблица «уровень × программа» и опыт, который это запускает."""
        self.h2("Чем считается")
        self.table(["уровень", "SynchrotronSim", "elegant"], rows,
                   title=f"Инструменты расчёта (глава {self.n})")
        self.p(f"Запуск: py -3.14 verify_elegant.py {experiment} --machine {vr.MACHINE}.")

    def results(self):
        self.h2("Результаты")

    def diffs(self, rows, note=None):
        """Расхождения — отдельной таблицей в конце главы."""
        self.h2("Расхождения")
        if note:
            self.p(note)
        self.table(DIFF_HEADER, rows, digits=[0, 0, 3, 0],
                   title=f"Расхождения инструментов (глава {self.n})")


# ============================================================================ графики глав
def fig_orbit():
    ss, el = npz("linear", "ss.npz"), npz("linear", "elegant_track.npz")
    s = el["s"]
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    for a, k, lab in ((ax[0], 0, "x"), (ax[1], 2, "y")):
        a.plot(ss["s"], ss["orbit"][:, k] * 1e6, "-", color=C_SS, lw=1.3, label="SynchrotronSim, трекинг")
        a.plot(s, el["orbit"][:, k] * 1e6, "o", color=C_EL, ms=3, label="elegant, трекинг")
        a.set_ylabel(f"{lab} замкнутой орбиты, мкм"); a.grid(alpha=0.3)
    ax[0].legend(fontsize=9); ax[1].set_xlabel("s, м (выходы магнитов)")
    fig.suptitle("Замкнутая орбита идеальной решётки вдоль кольца")
    fig.tight_layout(); p = os.path.join(figdir(), "orbit.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_beta():
    ss, el = npz("linear", "ss.npz"), npz("linear", "elegant_track.npz")
    tw, an = npz("linear", "elegant_twiss.npz"), npz("analytic", "ss_matrix.npz")
    s = tw["s"]; per = s[-1] / vr.M["superperiods"]; k = s <= per
    fig, ax = plt.subplots(3, 1, figsize=(11, 10))
    for key, col, lab in (("betax", C_SS, "βx"), ("betay", C_EL, "βy")):
        ax[0].plot(s[k], an[key][k], "-", color=col, lw=3.0, alpha=0.4, label=f"{lab}: SynchrotronSim, аналитика")
        ax[0].plot(s[k], tw[key][k], "--", color=col, lw=1.1, label=f"{lab}: elegant, twiss_output")
        ax[0].plot(s[k], el[key][k], "o", color=col, ms=3, mfc="none", label=f"{lab}: elegant, трекинг")
        ax[0].plot(s[k], ss[key][k], ".", color="k", ms=2.5)
    ax[0].plot([], [], ".", color="k", label="SynchrotronSim, трекинг")
    ax[0].set_xlabel("s, м"); ax[0].set_ylabel("β, м"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=7, ncol=2)
    ax[0].set_title("один супериод: четыре инструмента")
    for key, ls in (("betax", "-"), ("betay", "--")):
        ax[1].semilogy(s, np.abs(an[key] - tw[key]) / tw[key] + 1e-16, ls, color=C_AN, lw=0.9,
                       label=f"аналитика ↔ аналитика, {key[-1]}")
        ax[1].semilogy(s, np.abs(ss[key] - el[key]) / el[key] + 1e-16, ls, color=C_SS, lw=0.9,
                       label=f"трекинг ↔ трекинг, {key[-1]}")
    ax[1].set_ylabel("|Δβ/β|"); ax[1].grid(alpha=0.3, which="both"); ax[1].legend(fontsize=8, ncol=2)
    ax[1].set_title("относительное расхождение вдоль всего кольца")
    for key, ls in (("psix", "-"), ("psiy", "--")):
        ax[2].plot(s, (ss[key] - el[key]) / (2 * np.pi), ls, color=C_SS, lw=1.0, label=f"трекинг ↔ трекинг, {key[-1]}")
        ax[2].plot(s, (tw[key] - el[key]) / (2 * np.pi), ls, color=C_TW, lw=0.8,
                   label=f"twiss_output − трекинг elegant, {key[-1]}")
    ax[2].set_xlabel("s, м"); ax[2].set_ylabel("Δψ, ед. 2π"); ax[2].grid(alpha=0.3); ax[2].legend(fontsize=8, ncol=2)
    ax[2].set_title("расхождение фазового набега")
    fig.tight_layout(); p = os.path.join(figdir(), "beta.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


def fig_eta():
    ss, el = npz("linear", "ss.npz"), npz("linear", "elegant_track.npz")
    tw, an = npz("linear", "elegant_twiss.npz"), npz("analytic", "ss_matrix.npz")
    s = tw["s"]; per = s[-1] / vr.M["superperiods"]; k = s <= per
    fig, ax = plt.subplots(2, 1, figsize=(11, 7.5))
    ax[0].plot(s[k], an["etax"][k], "-", color=C_AN, lw=3.0, alpha=0.4, label="SynchrotronSim, аналитика")
    ax[0].plot(s[k], tw["etax"][k], "--", color=C_TW, lw=1.1, label="elegant, twiss_output")
    ax[0].plot(s[k], el["etax"][k], "o", color=C_EL, ms=3, mfc="none", label="elegant, трекинг")
    ax[0].plot(s[k], ss["etax"][k], ".", color=C_SS, ms=3, label="SynchrotronSim, трекинг")
    ax[0].set_xlabel("s, м"); ax[0].set_ylabel("ηx, м"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
    ax[0].set_title("один супериод: четыре инструмента")
    ax[1].plot(s, (an["etax"] - tw["etax"]) * 1e6, color=C_AN, lw=1.0, label="аналитика ↔ аналитика")
    ax[1].plot(s, (ss["etax"] - el["etax"]) * 1e6, color=C_SS, lw=1.0, label="трекинг ↔ трекинг")
    ax[1].plot(s, (ss["etax"] - an["etax"]) * 1e6, ":", color="k", lw=1.0, label="внутри SynchrotronSim")
    ax[1].set_xlabel("s, м"); ax[1].set_ylabel("Δηx, мкм"); ax[1].grid(alpha=0.3); ax[1].legend(fontsize=8)
    ax[1].set_title("расхождение вдоль всего кольца")
    fig.tight_layout(); p = os.path.join(figdir(), "eta.png"); fig.savefig(p, dpi=130); plt.close(fig)
    return p


def rad_amplitudes():
    """Амплитуды классического затухания (инвариант Куранта–Снайдера) по кэшу опыта rad: {плоскость: (SS, el)}."""
    lin = npz("linear", "ss.npz")
    lo = vr.csv_dict(res("longitudinal", "summary.csv"))
    qs = lo.get("Qs (ВЧ + излучение), матрица", [0])[0]
    beta_s = fnum(lin["alphac"]) * fnum(lin["s"][-1]) / (2 * np.pi * qs) if qs else 0.0
    tw = {0: (fnum(lin["betax"][-1]), fnum(lin["alphax"][-1])), 2: (fnum(lin["betay"][-1]),
          fnum(lin["alphay"][-1])), 4: (beta_s, 0.0)}
    out = {}
    for col, (pl, ia, ib) in enumerate((("x", 0, 1), ("y", 2, 3), ("t", 4, 5))):
        pair = []
        for code in ("ss", "elegant"):
            d = npz("radiation", f"{code}_sr.npz")
            tr, co = d["tbt"][col], d["orbit"]
            ok = np.all(np.isfinite(tr), axis=1)
            q, pp = tr[ok, ia] - co[ia], tr[ok, ib] - co[ib]
            b, al = tw[ia]
            pair.append(np.hypot(q, al * q + b * pp) * 1e3)
        out[pl] = tuple(pair)
    return out


# ============================================================================ сборка
def build(out=None):
    M = vr.M
    f = vr.lattice_facts()
    nsl, order = M["slices"]
    volt, harm = M["rf"]["volt_mv"], M["rf"]["harmonic"]
    P = vr.pretty
    lin = vr.csv_dict(res("linear", "summary.csv"))
    an = vr.csv_dict(res("analytic", "summary.csv"))
    u0 = npz("longitudinal", "u0.npz")
    lo = vr.csv_dict(res("longitudinal", "summary.csv")) if have("longitudinal") else {}
    em = vr.csv_dict(res("emittance", "summary.csv")) if have("emittance") else {}
    ch = vr.csv_dict(res("chromaticity", "summary.csv")) if have("chromaticity") else {}
    C = f["circumference"]; T0 = C / C_LIGHT
    gamma = M["e_mev"] / ed.ME_MEV
    brho = M["e_mev"] * 1e6 / C_LIGHT
    phi_s = 180.0 - math.degrees(math.asin(fnum(u0["U0"]) / volt)) if u0 else float("nan")
    lat = ed.machine_lattice(M)
    par = ed.element_params(lat[0], lat[1]); names = ed.line_names(*lat)
    count = {}
    for n in names:
        t = par[n.upper()]["type"].upper()
        count.setdefault(t, set()).add(n.upper())
    n_of = lambda t: sum(1 for n in names if par[n.upper()]["type"].upper() == t)

    R = Doc()
    R.doc.add_heading(f"Сверка SynchrotronSim с elegant по параметрам: {M['label']}", level=0)
    R.p("Каждая глава посвящена одному рассчитываемому параметру: что он означает физически, каким "
        "модулем комплекса считается, что получено и насколько программы расходятся. Первая глава "
        "описывает установку и настройки обеих программ.")

    # ================================================================ 1. установка и настройки
    R.chapter("Установка и настройки сравнения")
    R.h2("Установка")
    R.p(f"Кольцо «{M['label']}». Описание — файл "
        f"{os.path.relpath(M['source'], ed.ROOT) if M['source'] else 'configuration_v1_geo'}"
        f"{' (перевод из PyAT)' if M['source'] else ''}. Кольцо из {M['superperiods']} одинаковых "
        f"супериодов, {f['n_line']} {vr.plural(f['n_line'], 'элемент', 'элемента', 'элементов')} линии, из "
        f"них {f['n_mag']} {vr.plural(f['n_mag'], 'магнит', 'магнита', 'магнитов')}.")
    R.p(vr.db_note())
    rows = [
        ("периметр, м", f["circumference"]),
        ("энергия пучка, МэВ", M["e_mev"]),
        ("релятивистский фактор γ", gamma),
        ("магнитная жёсткость Bρ, Тл·м", brho),
        ("период обращения T₀, нс", T0 * 1e9),
        ("число супериодов", M["superperiods"]),
        ("квадруполи (KQUAD)", f"{n_of('KQUAD')}, семейств {len(count.get('KQUAD', ()))}"),
        ("изгибы (CSBEND)", f"{n_of('CSBEND')}, типов {len(count.get('CSBEND', ()))}"),
        ("секступоли (KSEXT)", f"{n_of('KSEXT')}, сила нулевая" if M["source"] else n_of("KSEXT")),
        ("датчики положения (MONI)", n_of("MONI")),
    ]
    tw = npz("linear", "elegant_twiss.npz")
    if tw:
        rows += [("частоты νx / νy", f"{P(tw['nux'], 8)} / {P(tw['nuy'], 8)}"),
                 ("коэффициент уплотнения орбит αc", fnum(tw["alphac"])),
                 ("натуральная хроматичность ξx / ξy",
                  f"{P(tw['dnux/dp'], 6)} / {P(tw['dnuy/dp'], 6)}")]
    if u0:
        rows += [("потеря энергии за оборот U₀, МэВ", fnum(u0["U0"])),
                 ("времена затухания τx / τy / τs, мс",
                  f"{P(fnum(u0['taux']) * 1e3, 5)} / {P(fnum(u0['tauy']) * 1e3, 5)} / "
                  f"{P(fnum(u0['taudelta']) * 1e3, 5)}"),
                 ("равновесный эмиттанс εx, нм", fnum(u0["ex0"]) * 1e9),
                 ("энергетический разброс σδ", fnum(u0["Sdelta0"]))]
    rows += [("ВЧ: напряжение, МВ / кратность", f"{P(volt, 4)} / {harm}"),
             ("ВЧ: частота, МГц", harm * C_LIGHT / C / 1e6)]
    if u0:
        rows.append(("синхронная фаза φs, °", phi_s))
    R.table(["величина", "значение"], rows, digits=[0, 8], title="Параметры установки")
    R.p("Характеристики кольца в таблице — по аналитике elegant (twiss_output с интегралами излучения); "
        "ВЧ задано для сверки: кратность даёт штатную частоту около 181 МГц, напряжение выбрано с запасом "
        "над потерей энергии.")
    if M["source"]:
        R.p("Решётка взята без правок. Секступоли в файле нулевые, поэтому хроматизм натуральный, а "
            "нелинейность кольца — только кинематическая (точные дрейфы, края изгибов).")
        R.p(vr.corrector_note())

    R.h2("Настройки программ")
    R.p("Обе программы считают одну и ту же решётку одним и тем же интегратором; где устройство программ "
        "различается принципиально, это отмечено.")
    R.table(["параметр", "SynchrotronSim", "elegant"], [
        ("версия, ядро", "продакшн-ядро C++ (src/matrix_core.h, src/matrix_ring.h), сервер кольца "
         "build/MatrixTracker.exe --ring-serve; секторный режим", "elegant 2026.3.0 и SDDS ToolKit"),
        ("исполнение трекинга", "C++, частицы параллельно (OpenMP); Python-обёртка ring_cpp.py только "
         "передаёт координаты", "отдельный запуск elegant на пачку частиц"),
        ("описание решётки", f"конфигурация {M['lattice']} базы данных комплекса (импорт того же файла)",
         f"файл целиком, линия «{M['beamline'] or 'RING'}» обёрнута в {M['use_beamline']}"),
        ("энергия", f"{M['e_mev']:g} МэВ — для излучения и времени пролёта",
         f"p_central_mev = {M['e_mev']:g}"),
        ("изгибы", "секторное тело, края — тонкие линзы tan(E)/ρ", "CSBEND"),
        ("прямые мультиполи", "симплектическое тело — порт integrate_kick_multipole_ordn",
         "KQUAD, KSEXT"),
        ("срезов на элемент", f"{nsl} (параметр nsl сервера кольца)",
         f"N_SLICES = {nsl} (по умолчанию 1 у KQUAD и KSEXT, 4 у CSBEND)"),
        ("порядок схемы", f"{order} (параметр order, схема Йошиды)", f"INTEGRATION_ORDER = {order}"),
        ("промежутки между элементами", "модифицированный дрейф: прямой пролёт между гранями, грани "
         "могут быть непараллельны", "DRIF — прямой дрейф"),
        ("излучение", "распределено по срезам изгибов: классическое (флаг sr) и квантовое (флаг isr) "
         "с независимым генератором у каждой частицы",
         "CSBEND: SYNCH_RAD = 1, ISR = 1"),
        ("ВЧ", "в точке наблюдения в конце оборота: ΔE = eV·sin(φ + 2πht/C)", "RFCA, FIDUCIAL = light"),
        ("точка наблюдения", f"s = 0, до входа первого магнита {P(f['lead'], 5)} м",
         "WATCH в начале линии"),
        ("координаты", "x, px = (1+δ)·x′/√(1+x′²+y′²), y, py, t (запаздывание, м), δ",
         "x, x′, y, y′, t (с), p — переводятся в координаты SynchrotronSim"),
        ("аналитика", "ring_matrix.py: линейные матрицы элементов и аффинные стыки; ring_map.py: "
         "отображение оборота до 5-го порядка (дифференциальная алгебра)",
         "twiss_output (default_order = 3), matrix_output"),
    ], title="Настройки обеих программ")
    R.p(f"Число срезов у elegant задано явно. По умолчанию линза проходится одним шагом интегратора, и "
        f"на короткой сильной оптике этой установки собственный трекинг elegant выходит неустойчивым — "
        f"подробно в главе о порядке интегратора. С {nsl} срезами обе программы сходятся к точному "
        f"решению.")

    R.h2("Как строится сравнение")
    R.p("В каждой программе есть два механизма анализа. АНАЛИТИЧЕСКИЙ — матрица каждого элемента "
        "считается формулой, матрицы перемножаются вдоль кольца. ТРЕКИНГОВЫЙ — частицы гоняются по "
        "кольцу, а показатель извлекается из траекторий. Трекинговые показатели обеих программ "
        "обрабатываются ОДНИМ и тем же кодом (ring_analysis.py): он подключается к C++-ядру SynchrotronSim "
        "(ring_cpp.CppRing) и к elegant (elegant_diag.ElegantRing) через общий интерфейс, поэтому "
        "различие результатов — это различие самих трекеров, а не обработки.")
    R.p("Сравниваются однородные пары: аналитика с аналитикой, трекинг с трекингом. Разница между "
        "уровнями внутри одной программы — это ошибка дискретизации интегратора, её надо отличать от "
        "различия программ; где она важна, она указана в таблице расхождений отдельной строкой.")
    R.p("Матрицы переноса из трекинга снимаются центральными разностями 4-го порядка с шагом 10⁻⁵; "
        "частоты по пооборотным координатам — методом NAFF с окном Ханна; замкнутая орбита — методом "
        "Ньютона по отображению за оборот. Расхождение частот указывается абсолютным, остальных величин "
        "— относительным, если не сказано иное.")
    R.formula(r"\frac{\partial f}{\partial z} \approx \frac{-f(z+2h) + 8f(z+h) - 8f(z-h) + f(z-2h)}{12h}",
              [
        (r"f", "координата частицы после оборота (отображение за оборот)", "м или рад"),
        (r"z", "фазовая координата, по которой берётся производная: x, pₓ, y, p_y, t, δ", "м, рад, м, рад, м, —"),
        (r"h", "шаг разности, 10⁻⁵", "в единицах координаты z"),
    ])

    vr.mass_section(R, R.h2)

    # ================================================================ 2. замкнутая орбита
    if have("linear", "analytic"):
        R.chapter("Замкнутая орбита")
        R.meaning()
        R.p("Замкнутая орбита — траектория, которая через оборот возвращается в себя. Относительно неё "
            "отсчитываются все колебания частиц, вокруг неё линеаризуется движение и снимаются матрицы. "
            "В идеальной решётке она совпадает с осью; ошибки положения и поля магнитов её смещают. "
            "Математически это неподвижная точка отображения за оборот; в аффинной модели, где у стыка "
            "есть свободный член C, она находится одной формулой.")
        R.formula(r"T(z_0) = z_0,\qquad z_0 = (I - R)^{-1} C",
                  [
            (r"T", "отображение за оборот (трекинг или матрицы)"),
            (r"z_0", "замкнутая орбита — неподвижная точка отображения", "м, рад"),
            (r"R,\ C", "линейная матрица и свободный член отображения оборота", "—; м, рад"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.closed_orbit: (I − R)⁻¹·C по аффинной матрице оборота",
             "в идеальном описании орбита — ось по построению"),
            ("трекинг", "ring_analysis.closed_orbit: Ньютон по отображению C++-ядра",
             "тот же Ньютон поверх трекинга elegant (elegant_diag.ElegantRing)"),
        ], "linear")
        R.results()
        R.fig(fig_orbit(), "Замкнутая орбита вдоль кольца по трекингу обеих программ.")
        ss_l, el_l = npz("linear", "ss.npz"), npz("linear", "elegant_track.npz")
        resp = vr.csv_dict(res("analytic", "response.csv"))
        co_r = [v for k, v in resp.items() if k.startswith("замкнутая орбита")]
        dx_r = [(k, v) for k, v in resp.items() if k.startswith("Δx от поворота")]
        R.table(["величина", "SynchrotronSim", "elegant"], [
            ("x у точки наблюдения (трекинг), мкм", lin["замкнутая орбита x, мкм"][0],
             lin["замкнутая орбита x, мкм"][1]),
            ("max |x| вдоль кольца (трекинг), мкм", float(np.max(np.abs(ss_l["orbit"][:, 0]))) * 1e6,
             float(np.max(np.abs(el_l["orbit"][:, 0]))) * 1e6),
            ("max |y| вдоль кольца (трекинг), мкм", float(np.max(np.abs(ss_l["orbit"][:, 2]))) * 1e6,
             float(np.max(np.abs(el_l["orbit"][:, 2]))) * 1e6),
            ("x у точки наблюдения (аналитика), мкм",
             co_r[0][0] * 1e-3 if co_r else "", 0.0),
        ], digits=[0, 6, 6], title="Замкнутая орбита идеальной решётки")
        co_ss, co_el = lin["замкнутая орбита x, мкм"][0], lin["замкнутая орбита x, мкм"][1]
        R.p(f"Замкнутая орбита идеальной решётки: у elegant {vr.si_len(abs(co_el) * 1e-6)}, у "
            f"SynchrotronSim {vr.si_len(abs(co_ss) * 1e-6)}" + vr.orbit_verdict(co_ss, co_el) + ". "
            f"Аналитика и трекинг SynchrotronSim дают одно и то же значение — значит, это свойство "
            f"описания решётки, а не счёта.")
        if dx_r:
            k, v = dx_r[0]
            R.p(f"Та же машинерия даёт орбиту от ошибок: при повороте одного квадруполя на 0,1 мрад "
                f"аналитика даёт смещение {P(v[0], 7)} мкм, трекинг ядра на сошедшемся интеграторе "
                f"(40 срезов, порядок 6) — {P(v[1], 7)} мкм.")
        diff = [("трекинг ↔ трекинг", "x у точки наблюдения, мкм (абсолютное)",
                 abs(co_ss - co_el),
                 vr.orbit_label(co_ss, co_el))]
        if co_r:
            diff.append(("внутри SynchrotronSim", "аналитика против трекинга, нм (абсолютное)",
                         abs(co_r[0][0] - co_r[0][1]), "разностная погрешность Ньютона"))
        if dx_r:
            diff.append(("внутри SynchrotronSim", "орбита от поворота квадруполя", dx_r[0][1][2],
                         "аффинная модель стыков против трекинга"))
        R.diffs(diff)

    # ================================================================ 3. бетатронные функции
    if have("linear", "analytic"):
        j = json.load(open(res("linear", "summary.json"), encoding="utf-8"))
        R.chapter("Бетатронные функции")
        R.meaning()
        R.p("Поперечное движение частицы — колебания вокруг замкнутой орбиты с переменной вдоль кольца "
            "амплитудой. Бетатронная функция β(s) задаёт огибающую этих колебаний: размер пучка равен "
            "√(εβ), чувствительность орбиты к ошибке магнита пропорциональна √β в его месте. Функция "
            "α = −β′/2 описывает наклон фазового эллипса, фазовый набег ψ(s) — сколько колебания "
            "набралось от начала кольца. Все три величины однозначно следуют из матрицы оборота.")
        R.formula(r"x(s) = \sqrt{\varepsilon\beta(s)}\cos\left(\psi(s)+\varphi_0\right),\quad "
                  r"\psi(s) = \int_0^s \frac{ds'}{\beta},\quad "
                  r"M = \begin{pmatrix}\cos\mu+\alpha\sin\mu & \beta\sin\mu\\ -\gamma\sin\mu & "
                  r"\cos\mu-\alpha\sin\mu\end{pmatrix}",
                  [
            (r"x(s)", "поперечная координата частицы в сечении s", "м"),
            (r"\varepsilon", "эмиттанс (инвариант Куранта–Снайдера)", "м·рад"),
            (r"\beta(s),\ \alpha(s)", "бетатронная функция и её наклон, α = −β′/2", "м; —"),
            (r"\gamma = (1+\alpha^2)/\beta", "третья функция Твисса", "м⁻¹"),
            (r"\psi(s)", "фазовый набег от начала кольца", "рад"),
            (r"\varphi_0", "начальная фаза частицы", "рад"),
            (r"\mu", "фазовый набег на участке (за оборот — 2πν)", "рад"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.twiss: перемножение матриц элементов и стыков",
             "twiss_output"),
            ("трекинг", "ring_analysis.twiss: матрицы переноса разностями по траекториям C++-ядра",
             "ring_analysis.twiss поверх трекинга elegant"),
        ], "linear и analytic")
        R.results()
        R.fig(fig_beta(), "Бетатронные функции на супериоде, их расхождение вдоль кольца и расхождение "
                          "фазового набега.")
        ssl = npz("linear", "ss.npz"); ell = npz("linear", "elegant_track.npz")
        anm = npz("analytic", "ss_matrix.npz")
        R.table(["величина", "SynchrotronSim, аналитика", "SynchrotronSim, трекинг",
                 "elegant, twiss_output", "elegant, трекинг"], [
            ("βx в точке наблюдения, м", fnum(anm["betax"][-1]), fnum(ssl["betax"][-1]),
             fnum(tw["betax"][-1]), fnum(ell["betax"][-1])),
            ("βy в точке наблюдения, м", fnum(anm["betay"][-1]), fnum(ssl["betay"][-1]),
             fnum(tw["betay"][-1]), fnum(ell["betay"][-1])),
            ("max βx, м", fnum(anm["betax"].max()), fnum(ssl["betax"].max()), fnum(tw["betax"].max()),
             fnum(ell["betax"].max())),
            ("max βy, м", fnum(anm["betay"].max()), fnum(ssl["betay"].max()), fnum(tw["betay"].max()),
             fnum(ell["betay"].max())),
            ("αx в точке наблюдения", fnum(anm["alphax"][-1]) if "alphax" in anm else "",
             fnum(ssl["alphax"][-1]), fnum(tw["alphax"][-1]), fnum(ell["alphax"][-1])),
            ("αy в точке наблюдения", fnum(anm["alphay"][-1]) if "alphay" in anm else "",
             fnum(ssl["alphay"][-1]), fnum(tw["alphay"][-1]), fnum(ell["alphay"][-1])),
        ], digits=[0, 8, 8, 8, 8], title="Бетатронные функции и их экстремумы по кольцу")
        R.p("Кривые всех четырёх инструментов на супериоде неразличимы. Аналитика SynchrotronSim ведёт "
            "функции по формулам матриц — как twiss_output; трекинг ведёт их по матрицам, снятым с "
            "траекторий, — как такой же расчёт поверх elegant.")
        R.p(vr.phase_note())
        dpsi = float(np.max(np.abs(ssl["psix"] - ell["psix"]))) / (2 * np.pi)
        dpsy = float(np.max(np.abs(ssl["psiy"] - ell["psiy"]))) / (2 * np.pi)
        R.diffs([
            ("аналитика ↔ аналитика", "βx по .lte, max", an["max |Δβx/βx| (матрицы по .lte)"][2],
             "формулы совпадают до округления"),
            ("аналитика ↔ аналитика", "βx по решётке БД, max", an["max |Δβx/βx| (решётка БД)"][2],
             vr.geom_label()),
            ("аналитика ↔ аналитика", "βy по решётке БД, max", an["max |Δβy/βy| (решётка БД)"][2],
             "то же"),
            ("трекинг ↔ трекинг", "βx вдоль кольца, max", j["max_rel_beta_x_ss"],
             "предел интегратора, одинаковый у обеих программ"),
            ("трекинг ↔ трекинг", "βy вдоль кольца, max", j["max_rel_beta_y_ss"], "то же"),
            ("трекинг ↔ трекинг", "набег фазы ψx, max, ед. 2π (абсолютное)", dpsi,
             "накопление разницы частот вдоль кольца"),
            ("трекинг ↔ трекинг", "набег фазы ψy, max, ед. 2π (абсолютное)", dpsy, "то же"),
            ("внутри elegant", "twiss_output против трекинга, βx", j["max_rel_beta_x_twiss"],
             f"дискретизация интегратора, {nsl} срезов"),
        ])

    # ================================================================ 4. частоты
    if have("linear", "analytic"):
        R.chapter("Бетатронные частоты")
        R.meaning()
        R.p("Бетатронная частота ν — число поперечных колебаний за оборот, то есть полный фазовый набег "
            "кольца, делённый на 2π. Целая часть задаётся оптикой в целом, дробная определяет, насколько "
            "рабочая точка далека от резонансов: при рациональном ν малые ошибки накапливаются от оборота "
            "к обороту, и движение раскачивается. Поэтому частоты — главная сводная характеристика "
            "линейной оптики и самая чувствительная величина сверки.")
        R.formula(r"\nu = \frac{\mu}{2\pi} = \frac{1}{2\pi}\oint\frac{ds}{\beta(s)},\qquad "
                  r"\cos\mu = \frac{M_{11}+M_{22}}{2}",
                  [
            (r"\nu", "бетатронная частота: число колебаний за оборот", "—"),
            (r"\mu", "полный фазовый набег за оборот", "рад"),
            (r"\beta(s)", "бетатронная функция", "м"),
            (r"M_{11},\ M_{22}", "диагональ поперечного блока матрицы оборота", "—"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.twiss: набег фазы по матрицам", "twiss_output"),
            ("трекинг", "ring_analysis.twiss: набег фазы по матрицам разностями",
             "ring_analysis.twiss поверх трекинга elegant"),
        ], "linear и analytic")
        R.results()
        R.table(["инструмент", "νx", "νy"], [
            ("SynchrotronSim, аналитика по .lte", an["νx (матрицы по .lte)"][0], an["νy (матрицы по .lte)"][0]),
            ("SynchrotronSim, аналитика по решётке БД", an["νx (матрицы по решётке БД)"][0],
             an["νy (матрицы по решётке БД)"][0]),
            ("SynchrotronSim, трекинг", lin["νx"][0], lin["νy"][0]),
            ("elegant, twiss_output", lin["νx"][2], lin["νy"][2]),
            ("elegant, трекинг", lin["νx"][1], lin["νy"][1]),
        ], digits=[0, 11, 11], title="Бетатронные частоты четырьмя инструментами")
        R.p("«Матрицы по .lte» — аналитический механизм SynchrotronSim, собранный прямо по файлу "
            "описания, без базы данных: так проверяются сами формулы. «По решётке БД» — тот же механизм на "
            "той решётке, по которой считает трекинг.")
        R.diffs([
            ("аналитика ↔ аналитика", "νx по .lte", an["νx (матрицы по .lte)"][2], "до округления"),
            ("аналитика ↔ аналитика", "νy по .lte", an["νy (матрицы по .lte)"][2], "до округления"),
            ("аналитика ↔ аналитика", "νx по решётке БД", an["νx (матрицы по решётке БД)"][2],
             vr.geom_label()),
            ("аналитика ↔ аналитика", "νy по решётке БД", an["νy (матрицы по решётке БД)"][2], "то же"),
            ("трекинг ↔ трекинг", "νx", abs(lin["νx"][0] - lin["νx"][1]),
             "предел сходимости интегратора elegant"),
            ("трекинг ↔ трекинг", "νy", abs(lin["νy"][0] - lin["νy"][1]), "то же"),
            ("внутри SynchrotronSim", "νx: аналитика против трекинга",
             abs(an["νx (матрицы по решётке БД)"][0] - lin["νx"][0]), f"дискретизация, {nsl} срезов"),
            ("внутри elegant", "νx: twiss_output против трекинга", abs(lin["νx"][2] - lin["νx"][1]),
             f"дискретизация, {nsl} срезов"),
        ])

    # ================================================================ 5. дисперсия
    if have("linear", "analytic"):
        R.chapter("Дисперсионная функция и уплотнение орбит")
        R.meaning()
        R.p("Частица с отклонением импульса δ по-другому отклоняется в изгибах и идёт по своей замкнутой "
            "орбите, смещённой на η(s)·δ. Дисперсионная функция η возникает только в изгибах; она "
            "добавляет к размеру пучка вклад η·σδ, связывает поперечное движение с энергией и определяет "
            "удлинение пути частицы — коэффициент уплотнения орбит αc. Он управляет всем продольным "
            "движением: синхротронной частотой и длиной сгустка.")
        R.formula(r"x_\delta(s) = \eta(s)\,\delta,\qquad \eta'' + K(s)\,\eta = \frac{1}{\rho(s)},\qquad "
                  r"\alpha_c = \frac{\Delta C / C}{\delta} = \frac{1}{C}\oint\frac{\eta(s)}{\rho(s)}\,ds",
                  [
            (r"x_\delta(s)", "смещение орбиты частицы с импульсом (1+δ)·p₀", "м"),
            (r"\eta(s)", "дисперсионная функция", "м"),
            (r"\delta = \Delta p/p", "относительное отклонение импульса", "—"),
            (r"K(s)", "сила фокусировки: K = 1/ρ² + K₁", "м⁻²"),
            (r"\rho(s)", "радиус изгиба", "м"),
            (r"\alpha_c", "коэффициент уплотнения орбит", "—"),
            (r"C", "периметр кольца", "м"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.twiss: столбец δ матриц элементов и стыков", "twiss_output"),
            ("трекинг", "ring_analysis.twiss: столбец δ матриц разностями, αc — по удлинению пути",
             "ring_analysis.twiss поверх трекинга elegant"),
        ], "linear и analytic")
        R.results()
        R.fig(fig_eta(), "Дисперсия на супериоде и её расхождение вдоль кольца.")
        ssl = npz("linear", "ss.npz"); ell = npz("linear", "elegant_track.npz"); anm = npz("analytic", "ss_matrix.npz")
        R.table(["величина", "SynchrotronSim, аналитика", "SynchrotronSim, трекинг",
                 "elegant, twiss_output", "elegant, трекинг"], [
            ("max ηx, м", fnum(anm["etax"].max()), fnum(ssl["etax"].max()), fnum(tw["etax"].max()),
             fnum(ell["etax"].max())),
            ("ηx в точке наблюдения, м", fnum(anm["etax"][-1]), fnum(ssl["etax"][-1]),
             fnum(tw["etax"][-1]), fnum(ell["etax"][-1])),
            ("αc", fnum(anm["alphac"]), lin["αc"][0], lin["αc"][2], lin["αc"][1]),
        ], digits=[0, 9, 9, 9, 9], title="Дисперсия и коэффициент уплотнения орбит")
        d_tr = float(np.max(np.abs(ssl["etax"] - ell["etax"])))
        d_in = float(np.max(np.abs(ssl["etax"] - anm["etax"])))
        d_el = float(np.max(np.abs(ell["etax"] - tw["etax"])))      # разрыв уровней внутри самого elegant
        R.p(f"Аналитические пары сходятся до {vr.si_len(an['max Δηx, нм (решётка БД)'][2] * 1e-9)}. В "
            f"трекинге остаётся {P(d_tr * 1e3, 2)} мм — это {P(d_tr / fnum(ell['etax'].max()), 2)} от "
            f"максимума дисперсии. Столько же стоит между трекингом elegant и его собственным "
            f"twiss_output ({P(d_el * 1e3, 2)} мм), а внутри SynchrotronSim аналитика и трекинг расходятся "
            f"на {P(d_in * 1e3, 2)} мм. То есть след в дисперсии — это дискретизация интегратора, одинаковая "
            f"у обеих программ, а не различие моделей.")
        R.diffs([
            ("аналитика ↔ аналитика", "ηx по решётке БД, max, мкм (абсолютное)",
             an["max Δηx, нм (решётка БД)"][2] * 1e-3, vr.geom_label()),
            ("аналитика ↔ аналитика", "αc по .lte", an["αc (матрицы по .lte)"][2], "до округления"),
            ("трекинг ↔ трекинг", "ηx вдоль кольца, max, мм (абсолютное)", d_tr * 1e3,
             "дискретизация интегратора: столько же у elegant против его же twiss_output"),
            ("трекинг ↔ трекинг", "αc", vr.rel(lin["αc"][0], lin["αc"][1]), "то же"),
            ("внутри SynchrotronSim", "ηx: аналитика против трекинга, мм", d_in * 1e3, "то же"),
        ])

    # ================================================================ 6. матрицы элементов
    if have("matrices"):
        R.chapter("Матрицы элементов")
        R.meaning()
        R.p("Матрица элемента R — линейное отображение координат частицы через один магнит (вместе с "
            "промежутком перед ним). Верхний левый блок 2×2 — фокусировка по x, следующий — по y, столбец "
            "δ (R16, R26) — рождение дисперсии в изгибах, строка пути (R51, R52, R56) — удлинение "
            "траектории. Вся линейная оптика кольца — произведение этих матриц, поэтому сверка поэлементно "
            "показывает, в каком именно магните рождается расхождение. У стыков SynchrotronSim "
            "отображение аффинное: есть ещё свободный член C.")
        R.formula(r"z_{\text{вых}} = R\,z_{\text{вх}} + C,\qquad "
                  r"(C_2, R_2)\circ(C_1, R_1) = (C_2 + R_2 C_1,\; R_2 R_1)",
                  [
            (r"z_{\text{вх}},\ z_{\text{вых}}", "фазовый вектор на входе и выходе элемента", "м, рад"),
            (r"R", "матрица элемента или стыка 6×6"),
            (r"C", "свободный член стыка (наклон граней)", "м, рад"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.element_matrix и junction_matrix", "matrix_output (individual_matrices)"),
            ("трекинг", "ring_analysis.element_matrices: R_j·R_{j−1}⁻¹ из накопленных разностных матриц",
             "то же поверх трекинга elegant"),
        ], "matrices и analytic")
        R.results()
        R.fig(res("matrices", "element_matrices.png"), "Наибольшее по магнитам расхождение каждого члена "
                                                       "матрицы (log10).")
        hdr, by = vr.read_csv(res("matrices", "by_type.csv"))
        R.table(hdr, by, digits=[0, 6, 2, 2, 2, 2], title="Матрицы элементов по типам устройств")
        rh, rr = vr.read_csv(res("matrices", "representative.csv"))
        kept, per_type = [], {}
        for r in rr:                                   # по одному устройству каждого типа, у изгибов — все виды
            nm = r[0]
            if nm in kept:
                continue
            t = par.get(nm.upper(), {}).get("type", "?").upper()
            if per_type.get(t, 0) < (3 if t == "CSBEND" else 1):
                kept.append(nm); per_type[t] = per_type.get(t, 0) + 1
        R.p(f"Представительные элементы ({', '.join(kept)}), по одному каждого типа: члены фокусировки, "
            f"дисперсии и пути.")
        R.table(rh, [r for r in rr if r[0] in kept], digits=[0, 0, 9, 9, 9],
            title="Матрица представительного элемента")
        R.p("Матрицы стыков (модифицированных дрейфов) проверены против точного отображения ядра на "
            "четырёх геометриях граней:")
        R.csv_table(res("analytic", "junctions.csv"), 3, title="Матрицы стыков против точного отображения ядра")
        if os.path.exists(res("analytic", "junctions_high_order.csv")):
            R.p("Для хроматичности и других нелинейных величин стык нужен и в высших порядках. Модуль "
                "ring_map считает его точной формулой пролёта над степенными рядами; ниже — ошибка ряда "
                "5-го порядка и аффинной матрицы против точного пролёта C++-ядра при смещении частицы на 1, "
                "3 и 10 мм (и соответствующих углах и δ). Ошибка аффинной матрицы растёт как квадрат "
                "амплитуды, ошибка ряда — как шестая степень:")
            R.csv_table(res("analytic", "junctions_high_order.csv"), 3,
                title="Стык в высших порядках: ряд и матрица против точного пролёта")
        body = [r for r in by if not r[0].startswith("первый")]
        R.diffs([
            ("аналитика ↔ аналитика", "max |ΔR| по всем магнитам", max(vr.num(r[4]) for r in body),
             "формулы совпадают, остаток — геометрия базы данных"),
            ("трекинг ↔ трекинг", "медиана max |ΔR|", max(vr.num(r[3]) for r in body),
             "шум разностных производных 10⁻⁷…10⁻⁶"),
            ("трекинг ↔ трекинг", "max |ΔR|", max(vr.num(r[2]) for r in body), "то же"),
            ("внутри elegant", "matrix_output против трекинга", max(vr.num(r[5]) for r in body),
             f"дискретизация, {nsl} срезов"),
            ("стык ↔ точное отображение", "max |ΔR| на четырёх геометриях",
             max(vr.num(r[2]) for r in vr.read_csv(res("analytic", "junctions.csv"))[1]),
             "аффинная матрица стыка точна до округления"),
        ])

    # ================================================================ 7. интегратор
    if have("integrator"):
        R.chapter("Порядок кикового интегратора")
        R.meaning()
        R.p("Уравнение движения в магните решается не точно, а симплектической схемой: элемент режется на "
            "N срезов, внутри среза чередуются точные дрейфы и тонкие кики поля. Порядок схемы p (2, 4 или "
            "6 по Йошиде) задаёт, как быстро ошибка убывает с числом срезов: частоты сходятся к точному "
            "значению как N⁻ᵖ. Симплектичность — сохранение фазового объёма — выполняется при любом N, "
            "поэтому грубая схема даёт правдоподобное, но смещённое движение. Параметр определяет, "
            "насколько трекинг вообще приближает точное решение, и на нём держатся все трекинговые "
            "сравнения.")
        R.formula(r"e^{hH} \approx \prod_i e^{c_i h D}\, e^{d_i h K},\qquad \delta\nu \propto N^{-p}",
                  [
            (r"H", "гамильтониан элемента; D и K — его дрейфовая и киковая части"),
            (r"h", "длина шага интегратора (среза)", "м"),
            (r"c_i,\ d_i", "коэффициенты симплектической схемы Йошиды", "—"),
            (r"N", "число срезов на элемент", "—"),
            (r"p", "порядок схемы (2, 4 или 6)", "—"),
            (r"\delta\nu", "ошибка частоты от дискретизации", "—"),
        ])
        R.tools([
            ("трекинг", "C++-ядро с N срезами и схемой порядка p; частоты — ring_analysis.tunes_point",
             "N_SLICES = N и INTEGRATION_ORDER = p у KQUAD, KSEXT, CSBEND"),
            ("опорное значение", "—", "twiss_output (точные матрицы, от N не зависит)"),
        ], "integrator")
        R.results()
        R.fig(res("integrator", "integrator_convergence.png"),
              "Отклонение частот от аналитики и расхождение программ при одинаковых настройках.")
        R.csv_table(res("integrator", "tunes.csv"), 8, title="Частоты при разном числе срезов и порядке схемы")
        hdr, it = vr.read_csv(res("integrator", "tunes.csv"))
        good = [r for r in it if r[0].startswith("N=") and all(vr.num(v) is not None and math.isfinite(vr.num(v))
                                                               for v in r[1:5])]
        dm = max(max(abs(vr.num(r[1]) - vr.num(r[2])), abs(vr.num(r[3]) - vr.num(r[4]))) for r in good)
        dflt = it[0]
        R.p(f"При одинаковых N и p программы идут вместе во всём диапазоне: и сходимость, и неустойчивость "
            f"при самом грубом интеграторе у них одна и та же. Схема 6-го порядка выходит на аналитику с "
            f"N ≥ 5, схема 4-го — с N ≥ 20.")
        if M["source"]:
            R.p(f"Первая строка — настройки ПО УМОЛЧАНИЮ. У elegant это один шаг интегратора на линзу "
                f"(N_SLICES = 1 у KQUAD и KSEXT): на этой решётке ошибка фокусировки около 0,7 % на линзу, "
                f"и за {f['n_mag']} {vr.plural(f['n_mag'], 'магнит', 'магнита', 'магнитов')} собственный "
                f"трекинг elegant теряет устойчивость по горизонтали "
                f"(частоты νx нет, cos μx = 4,76 вместо 0,139), хотя его же twiss_output даёт нормальные "
                f"функции. Ядро SynchrotronSim по умолчанию берёт 5 срезов 4-го порядка и остаётся "
                f"устойчивым, но ошибается в νx на {P(abs(vr.num(dflt[1]) - lin['νx'][2] % 1), 2)}. Отсюда "
                f"настройка всех глав этого отчёта — {nsl} срезов у обеих программ.")
        R.diffs([
            ("трекинг ↔ трекинг", f"частоты, {len(good)} устойчивых сочетаний N и p, max", dm,
             "предел сходимости интегратора elegant"),
            ("трекинг ↔ аналитика", "частота при N = 40, порядок 6, SynchrotronSim",
             abs(vr.num([r for r in it if r[0] == "N=40, порядок 6"][0][1]) - lin["νx"][2] % 1),
             "сходимость ядра"),
            ("трекинг ↔ аналитика", "частота при N = 40, порядок 6, elegant",
             abs(vr.num([r for r in it if r[0] == "N=40, порядок 6"][0][2]) - lin["νx"][2] % 1),
             "сходимость интегратора elegant"),
        ])

    # ================================================================ 8. хроматичность
    if have("chromaticity"):
        R.chapter("Бетатронная частота от энергии")
        R.meaning()
        R.p("Частица с большим импульсом слабее отклоняется в линзах, поэтому её бетатронная частота "
            "другая. Зависимость ν(δ) раскладывается в ряд; коэффициенты — хроматичности. Первая, ξ1, "
            "задаёт разброс частот в пучке с энергоразбросом и его устойчивость к коллективным эффектам; "
            "натуральная хроматичность всегда отрицательна и компенсируется секступолями в местах с "
            "дисперсией. Старшие ξ2, ξ3 описывают кривизну ν(δ) и важны для энергетической апертуры.")
        R.formula(r"\nu(\delta) = \nu_0 + \xi_1\delta + \frac{\xi_2}{2}\delta^2 + \frac{\xi_3}{6}\delta^3,\qquad "
                  r"\xi_1^{\text{нат}} \approx -\frac{1}{4\pi}\oint\beta(s)\,K(s)\,ds",
                  [
            (r"\nu(\delta)", "бетатронная частота частицы с отклонением импульса", "—"),
            (r"\delta = \Delta p/p", "относительное отклонение импульса", "—"),
            (r"\xi_1,\ \xi_2,\ \xi_3", "хроматичности первого, второго и третьего порядка", "—"),
            (r"\beta(s),\ K(s)", "бета-функция и сила линзы вдоль кольца", "м; м⁻²"),
        ])
        R.tools([
            ("аналитика", "ring_map: отображение оборота 5-го порядка (дифференциальная алгебра) — замкнутая "
             "орбита рядом по δ, ν(δ) по линеаризации вокруг неё; путь по орбите — αc, αc2",
             "twiss_output, higher_order_chromaticity при default_order = 3"),
            ("трекинг", "ring_analysis.tunes_point: замкнутая орбита при δ и матрица оборота (C++-ядро); "
             "производные — подгонка 7-й степени на |δ| ≤ 1 %", "то же поверх трекинга elegant"),
        ], "chrom")
        R.p("Аналитическая ветвь SynchrotronSim — отображение оборота высших порядков (ring_map.py). Каждая "
            "величина там — усечённый степенной ряд по фазовым координатам (x, px, y, py, t, δ) до 5-го "
            "порядка; ряды складываются, перемножаются и подставляются друг в друга точно, поэтому "
            "коэффициенты — точные производные отображения. Составляющие оборота:")
        R.bullet("стык (модифицированный дрейф) — точная формула пролёта по прямой от грани до грани, "
                 "вычисленная над рядами: клин между гранями входит во все порядки;")
        R.bullet("тело элемента — точное решение его уравнений движения z(L) = exp(L·D)·z, где D — "
                 "производная по векторному полю элемента (то же поле, что интегрирует ядро, но без срезов);")
        R.bullet("края изгиба — тонкие линзы, как в ядре.")
        R.p("Замкнутая орбита при отклонении импульса находится как ряд по δ, линеаризация отображения "
            "вокруг неё даёт cos μ(δ) и ν(δ). При порядке отображения 5 точны ξ1…ξ4.")
        R.formula(r"z(L) = e^{L D} z,\qquad D = \sum_i f_i(z)\,\frac{\partial}{\partial z_i},\qquad "
                  r"\cos\mu(\delta) = \tfrac{1}{2}\,\mathrm{tr}\,\frac{\partial M}{\partial z}\Big|_{z_{co}(\delta)}",
                  [
            (r"L", "длина элемента", "м"),
            (r"f_i(z)", "компонента векторного поля элемента (правая часть уравнения движения)", "1/м"),
            (r"D", "производная Ли по этому полю", "м⁻¹"),
            (r"M", "отображение оборота"),
            (r"z_{co}(\delta)", "замкнутая орбита при отклонении импульса δ", "м, рад"),
            (r"\mu(\delta)", "фазовый набег за оборот при этом δ", "рад"),
        ])
        R.results()
        R.fig(res("chromaticity", "tune_vs_delta.png"), "ν(δ): трекинг обеих программ, ряд Тейлора "
                                                        "twiss_output и ряд отображения SynchrotronSim; "
                                                        "внизу — разность трекинга программ.")
        R.csv_table(res("chromaticity", "summary.csv"), 8, title="Хроматичность и коэффициенты уплотнения орбит")
        _, tn = vr.read_csv(res("chromaticity", "tunes.csv"))
        k2x, k3x, k2y, k3y = ("ξ2 = d²ν/dδ² (x)", "ξ3 = d³ν/dδ³ (x)", "ξ2 = d²ν/dδ² (y)", "ξ3 = d³ν/dδ³ (y)")
        R.p(f"Трекинг: {len(tn)} {vr.plural(len(tn), 'точка', 'точки', 'точек')} δ от −3 до +3 % (на |δ| ≤ 1 % "
            f"— с шагом 0,1 %), производные в нуле "
            f"— подгонкой многочленом 7-й степени (кубическая подгонка смещена членами δ⁴ и выше и портит "
            f"ξ2, ξ3 на проценты). "
            + ("Секступоли решётки нулевые, поэтому хроматизм здесь натуральный. " if M["source"] else ""))
        R.p(f"Отображение SynchrotronSim и трекинг обеих программ дают одни и те же старшие хроматичности: "
            f"ξ2x = {P(ch[k2x][0], 6)} по отображению против {P(ch[k2x][2], 6)} и {P(ch[k2x][3], 6)} по "
            f"трекингу, ξ3y = {P(ch[k3y][0], 6)} против {P(ch[k3y][2], 6)} и {P(ch[k3y][3], 6)}. Остаток "
            f"между отображением и трекингом — точность подгонки ν(δ) и сходимость интегратора ({nsl} срезов). "
            f"twiss_output даёт ξ2x = {P(ch[k2x][1], 5)} и ξ3x = {P(ch[k3x][1], 5)}: ξ3 по матрицам "
            f"3-го порядка определить нельзя (нужен 4-й), а ξ2 расходится с трекингом самого elegant — это "
            f"расхождение внутри elegant между его матричным инструментом и его трекингом. Коэффициенты "
            f"уплотнения орбит αc и αc2 отображение и twiss_output дают одинаково.")
        R.diffs([
            ("аналитика ↔ аналитика", "ξ1 (x)", vr.rel(ch["ξ1 = dν/dδ (x)"][0], ch["ξ1 = dν/dδ (x)"][1]),
             "ring_map против twiss_output"),
            ("аналитика ↔ аналитика", "ξ1 (y)", vr.rel(ch["ξ1 = dν/dδ (y)"][0], ch["ξ1 = dν/dδ (y)"][1]), "то же"),
            ("аналитика ↔ аналитика", "ξ2 (x)", vr.rel(ch[k2x][0], ch[k2x][1]),
             "twiss_output расходится и с трекингом elegant"),
            ("аналитика ↔ аналитика", "ξ3 (x)", vr.rel(ch[k3x][0], ch[k3x][1]),
             "для ξ3 у twiss_output мал порядок матриц"),
            ("аналитика ↔ аналитика", "αc", vr.rel(ch["αc"][0], ch["αc"][1]), vr.geom_label()),
            ("аналитика ↔ аналитика", "αc2", vr.rel(ch["αc2"][0], ch["αc2"][1]), "то же"),
            ("трекинг ↔ трекинг", "ξ1 (x)", vr.rel(ch["ξ1 = dν/dδ (x)"][2], ch["ξ1 = dν/dδ (x)"][3]),
             "предел интегратора"),
            ("трекинг ↔ трекинг", "ξ1 (y)", vr.rel(ch["ξ1 = dν/dδ (y)"][2], ch["ξ1 = dν/dδ (y)"][3]), "то же"),
            ("трекинг ↔ трекинг", "ξ2 (x)", vr.rel(ch[k2x][2], ch[k2x][3]), "то же, плюс подгонка"),
            ("трекинг ↔ трекинг", "ξ3 (x)", vr.rel(ch[k3x][2], ch[k3x][3]), "то же, плюс подгонка"),
            ("аналитика ↔ трекинг SynchrotronSim", "ξ2 (x)", vr.rel(ch[k2x][0], ch[k2x][2]),
             "подгонка ν(δ) трекинга"),
            ("аналитика ↔ трекинг SynchrotronSim", "ξ3 (y)", vr.rel(ch[k3y][0], ch[k3y][2]), "то же"),
            ("аналитика ↔ трекинг SynchrotronSim", "αc2", vr.rel(ch["αc2"][0], ch["αc2"][2]), "то же"),
        ])

    # ================================================================ 9. фазовые траектории
    tb_hdr, tb = vr.read_csv(res("tbt", "tunes_vs_amplitude.csv")) if have("tbt") else (None, [])
    if have("tbt"):
        R.chapter("Фазовые траектории в горизонтальной и вертикальной плоскости")
        R.meaning()
        R.p("Фазовая траектория — последовательность точек (x, px) или (y, py) в точке наблюдения, оборот "
            "за оборотом (отображение Пуанкаре). В линейной оптике точки лежат на эллипсе площади πε, "
            "форма которого задана α и β. Искажение эллипса показывает нелинейность, размытие — хаос, "
            "уход в бесконечность — потерю частицы. Сравнение портретов — самая прямая проверка того, "
            "что программы одинаково ведут частицу на больших амплитудах.")
        R.formula(r"\varepsilon = \gamma x^2 + 2\alpha x p_x + \beta p_x^2 = \text{const}",
                  [
            (r"\varepsilon", "инвариант Куранта–Снайдера (площадь эллипса, делённая на π)", "м·рад"),
            (r"x,\ p_x", "координата и канонический импульс в точке наблюдения", "м; рад"),
            (r"\alpha,\ \beta,\ \gamma", "функции Твисса в точке наблюдения", "—; м; м⁻¹"),
        ])
        R.tools([
            ("трекинг", "CppRing.track_tbt (C++-ядро), пооборотно в точке наблюдения",
             "ElegantRing.track_tbt: WATCH в начале линии"),
        ], "tbt")
        R.results()
        ax_mm, ay_mm = M["amp_x_mm"], M["amp_y_mm"]
        R.p(f"Частицы стартуют от замкнутой орбиты со смещением по одной плоскости: {len(ax_mm)} амплитуд "
            f"по x от {P(ax_mm[0], 3)} до {P(ax_mm[-1], 3)} мм и {len(ay_mm)} по y от {P(ay_mm[0], 3)} до "
            f"{P(ay_mm[-1], 3)} мм; другая плоскость смещена на 10 мкм. Прогон — 1024 оборота.")
        R.fig(res("tbt", "phase_portraits.png"), "Фазовые портреты в точке наблюдения; точки elegant и "
                                                 "SynchrotronSim наложены.")
        R.fig(res("tbt", "tbt_difference.png"), "Расхождение пооборотных координат.")
        alive = [r for r in tb if vr.num(r[1]) == vr.num(r[4]) == 1024]
        R.p(f"Все {len(alive)} частиц из {len(tb)} прошли полный прогон в обеих программах. Портреты "
            f"совпадают; расхождение координат растёт линейно по оборотам — это набег фазы от разницы "
            f"частот, а не разъезд траекторий.")
        tr_ss, tr_el = npz("tbt", "ss.npz")["tbt"], npz("tbt", "elegant_track.npz")["tbt"]
        co_ss = npz("linear", "ss.npz")["closed_orbit"]; co_el = npz("linear", "elegant_track.npz")["closed_orbit"]
        dev_ss = tr_ss[:, :, [0, 2]] - co_ss[[0, 2]]; dev_el = tr_el[:, :, [0, 2]] - co_el[[0, 2]]
        amp = np.nanmax(np.abs(dev_el[:, :, :]), axis=(1, 2))                  # амплитуда частицы
        d1 = float(np.nanmax(np.abs(dev_ss[:, 1] - dev_el[:, 1]))) * 1e9
        dN = float(np.nanmax(np.nanmax(np.abs(dev_ss[:, -1] - dev_el[:, -1]), axis=1) / amp))
        dq_tbt = max(abs(vr.num(r[2]) - vr.num(r[5])) for r in alive)
        R.diffs([
            ("трекинг ↔ трекинг", "отклонение от замкнутой орбиты после 1 оборота, max, нм (абсолютное)", d1,
             "различие отображений за оборот; растёт с амплитудой"),
            ("трекинг ↔ трекинг", "то же после 1024 оборотов, в долях амплитуды", dN,
             f"набег фазы 2π·Δν·N при Δν = {P(dq_tbt, 2)}"),
            ("трекинг ↔ трекинг", "частиц, доживших до конца", f"{len(alive)} из {len(tb)}", "—"),
        ])

    # ================================================================ 10. спектр
    if have("tbt"):
        R.chapter("Спектр бетатронных колебаний")
        R.meaning()
        R.p("Спектр — разложение пооборотной координаты по частотам. Главная линия стоит на дробной "
            "бетатронной частоте; гармоники и комбинационные линии m·νx + n·νy появляются от нелинейности "
            "и связи плоскостей, их высота говорит о силе соответствующих резонансных членов. Положение "
            "главной линии методом NAFF определяется много точнее, чем шагом дискретного преобразования "
            "Фурье (ошибка убывает как 1/N⁴ при окне Ханна против 1/N).")
        R.formula(r"\hat{x}(\nu) = \sum_{k=0}^{N-1} w_k\, x_k\, e^{-2\pi i \nu k},\qquad "
                  r"w_k = \sin^2\frac{\pi k}{N}",
                  [
            (r"\hat{x}(\nu)", "амплитуда спектра на частоте ν", "м"),
            (r"x_k", "координата на k-м обороте", "м"),
            (r"w_k", "окно Ханна", "—"),
            (r"N", "число оборотов в выборке", "—"),
            (r"\nu", "частота в единицах обратного оборота", "—"),
        ])
        R.tools([
            ("обработка", "numpy.fft с окном Ханна; частоты — ring_analysis.naff / tbt_tunes",
             "та же обработка данных трекинга elegant"),
        ], "tbt")
        R.results()
        R.fig(res("tbt", "spectrum.png"), "Спектры пооборотных колебаний по x и по y (кривые сдвинуты по "
                                          "вертикали).")
        R.table(["частица", "νx SynchrotronSim", "νx elegant", "νy SynchrotronSim", "νy elegant"],
                [(r[0], r[2], r[5], r[3], r[6]) for r in tb], digits=[0, 10, 10, 10, 10],
                    title="Частоты по спектру для каждой пробной частицы")
        R.p("Линии обеих программ стоят в одних и тех же местах. "
            + ("Секступоли решётки нулевые, поэтому гармоники слабые: кольцо почти линейно, и спектр "
               "определяется главными линиями." if M["source"] else ""))
        dq = [max(abs(vr.num(r[2]) - vr.num(r[5])), abs(vr.num(r[3]) - vr.num(r[6]))) for r in tb
              if all(vr.num(r[i]) is not None and math.isfinite(vr.num(r[i])) for i in (2, 3, 5, 6))]
        R.diffs([("трекинг ↔ трекинг", "положение главных линий (NAFF), max", max(dq),
                  "предел интегратора elegant"),
                 ("трекинг ↔ трекинг", "положение главных линий (NAFF), min", min(dq), "то же")])

    # ================================================================ 11. частоты при разных амплитудах
    if have("tbt"):
        R.chapter("Частоты при разных амплитудах")
        R.meaning()
        R.p("В нелинейной решётке частота зависит от амплитуды колебаний: частица с большой амплитудой "
            "проходит больше нелинейного поля и набирает другую фазу. Этот сдвиг определяет размер следа "
            "пучка на диаграмме частот и то, не заденет ли пучок резонанс; он же — один из способов "
            "проверить, что программы одинаково описывают нелинейную часть поля и кинематику.")
        R.formula(r"\nu_x(J_x, J_y) = \nu_{x0} + \frac{\partial\nu_x}{\partial J_x}J_x + "
                  r"\frac{\partial\nu_x}{\partial J_y}J_y + \dots",
                  [
            (r"\nu_x", "горизонтальная частота частицы с действиями Jₓ и J_y", "—"),
            (r"J_x,\ J_y", "действие: J = (q² + (αq + βp)²)/(2β)", "м·рад"),
            (r"\partial\nu/\partial J", "сдвиг частоты с действием", "м⁻¹"),
            (r"\nu_{x0}", "частота при нулевой амплитуде", "—"),
        ])
        R.tools([
            ("трекинг", "ring_analysis.tbt_tunes (NAFF) по данным CppRing.track_tbt",
             "то же по данным трекинга elegant"),
            ("аналитика", "ring_map.normal_form: ∂ν/∂J нормальной формой отображения оборота",
             "twiss_output, compute_driving_terms: ∂ν/∂J и резонансные члены по мультиполям"),
            ("встроенный инструмент", "—", "tune_shift_with_amplitude (рассчитан на малый сдвиг частот)"),
        ], "tbt")
        R.results()
        R.fig(res("tbt", "tunes_vs_amplitude.png"), "Сдвиг частот с амплитудой; точка 0 — линейная оптика.")
        nx = [r for r in tb if r[0].startswith("x")]
        if nx:
            shift = vr.num(nx[-1][2]) - vr.num(nx[0][2])
            R.p(f"Сдвиг νx от {nx[0][0].split('=')[1].strip()} до {nx[-1][0].split('=')[1].strip()} — "
                f"{P(shift, 3)}"
                + (": секступоли нулевые, и остаётся только кинематическая нелинейность (точные дрейфы, края "
                   "изгибов). Обе программы воспроизводят её одинаково." if M["source"] else "."))
        dn = vr.read_csv(res("tbt", "detuning.csv"))[1]
        if dn:
            R.p("Аналитически сдвиг частот с амплитудой даёт нормальная форма отображения оборота "
                "(ring_map.normal_form): отображение вокруг замкнутой орбиты приводится заменой переменных "
                "к повороту, угол которого зависит от действий, и коэффициенты этой зависимости — ∂ν/∂J. "
                "Трекинговая оценка берётся по самым малым амплитудам сетки: дальше в ν(J) вступают члены "
                "J² и оценка наклона смещается.")
            R.table(["величина", "SynchrotronSim аналитика", "SynchrotronSim трекинг", "elegant трекинг",
                     "elegant twiss_output"],
                    [(r[0], vr.num(r[1]), vr.num(r[2]), vr.num(r[3]), vr.num(r[4]) if r[4] else "")
                     for r in dn], digits=[0, 6, 6, 6, 6], title="Сдвиг частот с действием ∂ν/∂J")
            R.p("Перекрёстные коэффициенты ∂νx/∂Jy и ∂νy/∂Jx обязаны совпадать — у аналитики они совпадают "
                "тождественно, у трекинга сходятся в пределах точности оценки наклона.")
        dt = vr.csv_dict(res("tbt", "driving_terms.csv"))
        if dt:
            R.p("Аналитическая ветвь есть только у elegant: twiss_output с compute_driving_terms считает по "
                "мультиполям решётки сдвиг частот с действием ∂ν/∂J и резонансные члены h_jklm. Кинематика "
                "(точные дрейфы, края изгибов) в эти формулы не входит"
                + (", и при нулевых секступолях ∂ν/∂J выходит нулевым, тогда как трекинг обеих программ "
                   "показывает ненулевой кинематический сдвиг. Резонансные члены при этом ненулевые: у "
                   "изгиба с градиентом есть собственный секступольный член 2·h·K₁."
                   if M["source"] else ". У SynchrotronSim такой ветви пока нет: её можно получить "
                   "нормальной формой того же отображения ring_map."))
            R.table(["величина", "elegant twiss_output (аналитика)"],
                    [(k, v[0]) for k, v in dt.items()], digits=[0, 6],
                        title="Аналитика нелинейности elegant: ∂ν/∂J и резонансные члены")
        rows_a = [("трекинг ↔ трекинг", "частоты на всех амплитудах, max", max(dq),
                   "постоянный сдвиг — предел интегратора, от амплитуды не зависит")]
        for r in dn or []:
            if vr.num(r[1]) and vr.num(r[2]):
                rows_a.append(("аналитика ↔ трекинг SynchrotronSim", r[0].split(",")[0],
                               vr.rel(vr.num(r[1]), vr.num(r[2])), "члены J² в оценке наклона"))
        R.diffs(rows_a)

    # ================================================================ 12. синхротронная частота
    if have("longitudinal") and lo:
        R.chapter("Синхротронная частота")
        R.meaning()
        R.p("ВЧ-резонатор восполняет потерю энергии и одновременно создаёт продольную «яму»: частица с "
            "лишней энергией идёт по длинной орбите, приходит позже, получает меньше энергии — и наоборот. "
            "Возникают продольные колебания энергии и фазы вокруг синхронной частицы; их число за оборот — "
            "синхротронная частота Qs. Она определяет длину сгустка и продольную устойчивость; из-за "
            "несинусоидальности ВЧ-ямы Qs падает с амплитудой.")
        R.formula(r"Q_s = \sqrt{\frac{h\,\alpha_c\,eV\,|\cos\varphi_s|}{2\pi E}},\qquad "
                  r"\Delta E = eV\sin\left(\varphi_s + \frac{2\pi h\,t}{C}\right)",
                  [
            (r"Q_s", "синхротронная частота: колебаний за оборот", "—"),
            (r"h", "кратность ВЧ", "—"),
            (r"\alpha_c", "коэффициент уплотнения орбит", "—"),
            (r"eV", "амплитуда ускоряющего напряжения", "МэВ"),
            (r"\varphi_s", "синхронная фаза", "рад"),
            (r"E", "энергия пучка", "МэВ"),
            (r"t", "запаздывание частицы относительно света", "м"),
            (r"C", "периметр кольца", "м"),
        ])
        R.tools([
            ("трекинг: матрица", "ring_analysis.one_turn_matrix + damping_from_matrix (собственные числа)",
             "то же поверх трекинга elegant с RFCA"),
            ("трекинг: спектр", "ring_analysis.naff по пооборотным (t, δ)", "то же"),
            ("аналитика", "ring_matrix.radiation_integrals: αc по I₁, синхронная фаза по U₀ из I₂",
             "rf_setup, twiss_output"),
        ], "long")
        R.results()
        R.fig(res("longitudinal", "synchrotron.png"),
              "Продольные портреты, расхождение по оборотам и синхротронная частота от амплитуды.")
        qs_rows = [(k, v[0], v[1], v[2]) for k, v in lo.items() if k.startswith("Qs")]
        R.table(["способ", "SynchrotronSim", "elegant", "формула / аналитика"], qs_rows,
                digits=[0, 9, 9, 9], title="Синхротронная частота тремя способами")
        R.p(f"ВЧ {P(volt, 3)} МВ с кратностью {harm}; без излучения синхронная фаза 180°, с излучением — "
            f"{P(phi_s, 5)}°. Qs по матрице оборота и по спектру пооборотных колебаний дают одно и то же "
            f"в обеих программах; формула малых амплитуд совпадает с ними до долей процента.")
        rows_q = [("трекинг ↔ трекинг", k, vr.rel(v[0], v[1]), "предел интегратора")
                  for k, v in lo.items() if k.startswith("Qs")]
        rows_q += [("аналитика ↔ аналитика", k, vr.rel(v[3], v[2]), "интегралы излучения против twiss_output")
                   for k, v in lo.items() if k.startswith("Qs") and len(v) >= 4 and v[3] is not None]
        R.diffs(rows_q)

    # ================================================================ 13. потери и затухание
    if have("longitudinal") and lo:
        R.chapter("Потери энергии на излучение и времена затухания")
        R.meaning()
        R.p("В изгибах частица излучает синхротронный свет и теряет энергию U₀ за оборот. ВЧ возвращает "
            "энергию только вдоль орбиты, поэтому поперечные импульсы, как и отклонение энергии, "
            "постепенно гасятся: колебания затухают. Времена затухания τx, τy, τs задают, как быстро пучок "
            "забывает начальные условия и приходит к равновесию; их сумма декрементов фиксирована "
            "теоремой Робинсона. Обе величины прямо проверяют модель излучения в изгибах.")
        R.formula(r"U_0 = \frac{C_\gamma E^4}{2\pi}\oint\frac{ds}{\rho^2},\qquad "
                  r"\tau_i = \frac{2E\,T_0}{J_i\,U_0},\qquad J_x + J_y + J_s = 4",
                  [
            (r"U_0", "потеря энергии на излучение за оборот", "МэВ"),
            (r"C_\gamma = 8{,}846\cdot10^{-5}", "постоянная излучения для электрона", "м/ГэВ³"),
            (r"E", "энергия пучка", "ГэВ"),
            (r"\rho", "радиус изгиба", "м"),
            (r"\tau_i", "время затухания в плоскости i", "с"),
            (r"T_0", "период обращения", "с"),
            (r"J_x,\ J_y,\ J_s", "числа затухания (теорема Робинсона)", "—"),
        ])
        R.tools([
            ("аналитика", "ring_matrix.radiation_integrals: U₀, Jx/Jy/Js и τ по интегралам I₁…I₅ "
             "(квадратура внутри изгибов, краевой член I₄)", "twiss_output (radiation_integrals)"),
            ("U₀", "трекинг осевой частицы за оборот C++-ядром с классическим излучением",
             "трекинг с SYNCH_RAD"),
            ("τ по матрице", "ring_analysis.damping_from_matrix (|λ| матрицы оборота)", "то же поверх elegant"),
            ("τ по спаду", "подгонка огибающей пооборотных колебаний", "то же"),
        ], "long")
        R.results()
        R.fig(res("longitudinal", "damping.png"), "Затухание амплитуды по оборотам в трёх плоскостях.")
        keys = ["U0 по трекингу за оборот, МэВ", "τx, мс", "τy, мс", "τs, мс", "Jx (по τ, ΣJ=4)",
                "Jy (по τ, ΣJ=4)", "Js (по τ, ΣJ=4)", "τx по спаду амплитуды, мс", "τy по спаду амплитуды, мс",
                "τs по спаду амплитуды, мс"]
        R.table(["величина", "SynchrotronSim трекинг", "elegant трекинг", "аналитика elegant",
                 "аналитика SynchrotronSim"],
                [(k, lo[k][0], lo[k][1], lo[k][2], lo[k][3] if len(lo[k]) > 3 and lo[k][3] is not None else "")
                 for k in keys if k in lo], digits=[0, 8, 8, 8, 8],
                     title="Потери энергии на излучение и времена затухания")
        R.p("Аналитика SynchrotronSim — интегралы излучения I₁…I₅ по описанию решётки (ring_matrix."
            "radiation_integrals): U₀ = CγE⁴I₂/2π, Jx = 1 − I₄/I₂, τ = 2ET₀/(J·U₀). Трекинг в них не "
            "участвует, и это прямая пара к radiation_integrals elegant.")
        R.p("Два независимых способа — по собственным числам матрицы оборота и по спаду амплитуды — дают "
            "одни и те же времена в обеих программах, сумма декрементов равна 4.")
        R.diffs([("аналитика ↔ аналитика", k.replace(" по трекингу за оборот", ""),
                  vr.rel(lo[k][3], lo[k][2]), "интегралы излучения против twiss_output")
                 for k in keys if k in lo and len(lo[k]) > 3 and lo[k][3] is not None and lo[k][2]
                 and "по спаду" not in k]
                + [("трекинг ↔ трекинг", k, vr.rel(lo[k][0], lo[k][1]),
                    "предел интегратора" if "τ" in k or "J" in k else "модель излучения совпадает")
                   for k in keys if k in lo]
                + [("трекинг ↔ аналитика elegant", "U₀", vr.rel(lo[keys[0]][0], lo[keys[0]][2]),
                    "интегралы излучения считаются по проектной орбите")])

    # ================================================================ 14. портреты с излучением
    if have("radiation"):
        R.chapter("Фазовые траектории с излучением")
        R.meaning()
        R.p("С излучением фазовая траектория перестаёт быть замкнутой. Классическое излучение — "
            "непрерывная потеря энергии — сворачивает эллипс в спираль к синхронной точке. Квантовое — "
            "испускание отдельных фотонов со случайной энергией — раскачивает частицу, и спираль "
            "останавливается на равновесном облаке, размер которого и есть эмиттанс. Портреты с "
            "излучением показывают оба механизма разом.")
        R.tools([
            ("трекинг", "CppRing.track_tbt с классическим (sr) или квантовым (isr) излучением, ВЧ в "
             "точке наблюдения", "ElegantRing.track_tbt с SYNCH_RAD или SYNCH_RAD + ISR, RFCA"),
        ], "rad")
        R.results()
        ax, ay, ad = M["rad_amp"]
        rs = npz("radiation", "ss_sr.npz"); ens = npz("radiation", "ss_ensemble.npz")
        n_turns = rs["tbt"].shape[1] - 1
        tx = fnum(u0["taux"]) / T0 if u0 else float("nan")
        ts = fnum(u0["taudelta"]) / T0 if u0 else float("nan")
        R.p(f"Три пробные частицы стартуют от синхронной точки со смещением x = {P(ax, 3)} мм, "
            f"y = {P(ay, 3)} мм и δ = {P(ad, 3)}; прогон — {n_turns} оборотов, это {P(n_turns / tx, 2)} "
            f"времени затухания по x и {P(n_turns / ts, 2)} продольного. Амплитуда колебания считается "
            f"инвариантом Куранта–Снайдера: так она не зависит от фазы, и разница частот программ не даёт "
            f"ложного разброса.")
        R.fig(res("radiation", "phase_portraits.png"),
              "Фазовые портреты: классическое излучение (сверху) и квантовое (снизу).")
        R.fig(res("radiation", "envelope.png"), "Амплитуда по оборотам.")
        tail = slice(int(0.98 * n_turns), None)
        f_el = float(np.nanstd(npz("radiation", "elegant_sr.npz")["tbt"][2][tail, 4]))
        f_ss = float(np.nanstd(rs["tbt"][2][tail, 4]))
        if f_el > 1.5 * f_ss:              # полка видна только когда прогон дошёл до предела точности elegant
            R.p(f"На продольной кривой (третий график) в самом низу программы расходятся: у elegant амплитуда "
                f"выходит на полку около {P(f_el * 1e3, 2)} мм, у SynchrotronSim продолжает экспоненту до "
                f"{P(f_ss * 1e3, 2)} мм. Это не физика, а разная бухгалтерия продольной координаты. elegant "
                f"выдаёт абсолютное время пролёта, и запаздывание получается разностью двух чисел порядка "
                f"{P(n_turns * f['circumference'], 2)} м: шаг двоичного представления такого числа 5·10⁻¹⁰ м, а "
                f"накопленное за {n_turns} оборотов округление даёт случайное блуждание на два порядка больше "
                f"(растёт как корень из числа оборотов). SynchrotronSim вычитает периметр на каждом обороте и "
                f"хранит сразу малую величину, поэтому такого предела не имеет. На сверку это не влияет: "
                f"огибающие сравниваются, пока амплитуда выше десятой части начальной — на четыре порядка выше "
                f"полки.")
        if ens and M.get("rad_ensemble"):
            n_ens, e_turns, e_skip, e_stride = M["rad_ensemble"]
            R.p(f"Равновесное облако квантового излучения проверяется отдельно АНСАМБЛЕМ: {n_ens} "
                f"{vr.plural(n_ens, 'частица стартует', 'частицы стартуют', 'частиц стартуют')} из синхронной точки, у каждой свой генератор квантового шума; первые {e_skip} "
                f"оборотов ({P(e_skip / tx, 2)} τx) частицы выходят на равновесие, затем ещё "
                f"{e_turns - e_skip} оборотов координаты записываются через каждые {e_stride}. СКО по "
                f"частицам и оборотам — размеры облака; статистическая ошибка оценивается по числу "
                f"независимых выборок (частицы × длина записи / время жизни действия τ/2).")
            R.fig(res("radiation", "equilibrium_clouds.png"),
                  "Равновесное облако ансамбля: фазовые плоскости (сверху) и распределения x и δ (снизу); "
                  "штрих — гауссово распределение с σ из moments_output elegant.")
        rd = vr.read_csv(res("radiation", "summary.csv"))[1]
        R.table(["величина", "SynchrotronSim", "elegant", "расхождение"],
                [(r[0], vr.num(r[1]) if r[1] else "", vr.num(r[2]) if r[2] else "", vr.num(r[3]) if r[3] else "")
                 for r in rd], digits=[0, 6, 6, 2], title="Излучение: затухание амплитуд и равновесное облако")
        R.p("Классическое затухание детерминировано, и кривые программ совпадают. Квантовые скачки "
            "разыгрываются независимыми генераторами случайных чисел, поэтому траектории совпадают только "
            "статистически: равновесные размеры ансамбля сравниваются друг с другом и с расчётом равновесия "
            "(глава об эмиттансе). Облако каждой программы совпадает с её же расчётом равновесия в пределах "
            "статистики; облако elegant на 3–4 % меньше облака SynchrotronSim — это тот же систематический "
            "сдвиг квантового шума в трекинге elegant, что и в главе об эмиттансе, где SynchrotronSim "
            "совпадает с moments_output. Вертикаль квантового возбуждения не получает вовсе (нет вертикальной "
            "дисперсии и связи): кривая с ним ложится на классическую.")
        R.diffs([(f"трекинг ↔ трекинг", r[0].split(" (")[0], vr.num(r[3]),
                  "предел интегратора" if "классическое" in r[0] else
                  ("шум в трекинге elegant плюс статистика ансамбля" if "ансамбль" in r[0]
                   else "см. главу об эмиттансе"))
                 for r in rd if r[3]])

    # ================================================================ 15. эмиттанс
    if have("emittance") and em:
        R.chapter("Равновесный эмиттанс и размеры пучка")
        R.meaning()
        R.p("Затухание сжимает пучок, квантовое возбуждение его раскачивает; равновесие между ними задаёт "
            "эмиттанс εx — площадь, занимаемую пучком на фазовой плоскости, — энергоразброс σδ и длину "
            "сгустка σz. Это главные параметры источника излучения: от εx зависит яркость. Эмиттанс "
            "рождается там, где частица излучает в месте с дисперсией, поэтому он прямо проверяет "
            "совместно оптику и модель излучения.")
        R.formula(r"\varepsilon_x = C_q\gamma^2\frac{I_5}{J_x I_2},\qquad "
                  r"\sigma_\delta^2 = C_q\gamma^2\frac{I_3}{J_s I_2},\qquad "
                  r"\Sigma = M\Sigma M^{T} + B",
                  [
            (r"\varepsilon_x", "равновесный горизонтальный эмиттанс", "м·рад"),
            (r"\sigma_\delta", "равновесный разброс по импульсу", "—"),
            (r"C_q = 3{,}832\cdot10^{-13}", "квантовая постоянная", "м"),
            (r"\gamma", "релятивистский фактор", "—"),
            (r"I_2,\ I_3,\ I_5", "интегралы излучения", "м⁻¹; м⁻²; м⁻¹"),
            (r"\Sigma,\ B", "матрица моментов пучка и матрица диффузии за оборот"),
            (r"M", "матрица оборота с излучением и ВЧ"),
        ])
        R.tools([
            ("равновесные моменты", "ring_analysis.isr_diffusion_matrix + equilibrium_moments (метод "
             "Чао): матрица оборота с излучением и ВЧ плюс диффузия по срезам изгибов (траектории по "
             "срезам — C++-ядро, CppRing.track_pieces)",
             "moments_output"),
            ("интегралы излучения", "ring_matrix.radiation_integrals: εx, σδ, σz по I₁…I₅",
             "twiss_output (radiation_integrals)"),
            ("статистика частиц", "—", "16000 частиц с ISR за оборот, ковариация — тот же метод Чао"),
        ], "emit")
        R.results()
        R.fig(res("emittance", "emittance.png"), "Равновесный эмиттанс и энергетический разброс.")
        R.csv_table(res("emittance", "summary.csv"), 6, title="Равновесный эмиттанс и размеры пучка")
        ex = em["εx, пм"]
        R.p(f"Метод Чао поверх C++-ядра даёт εx = {P(ex[0] / 1e3, 5)} нм против "
            f"{P(ex[3] / 1e3, 5)} нм по интегралам излучения elegant; σδ, σz и продольный эмиттанс "
            f"сходятся с moments_output до 0,2 %.")
        if ex[2] and abs(ex[2] / ex[3] - 1) > 0.1:
            R.p(f"Столбец moments_output по εx ({P(ex[2] / 1e3, 5)} нм) стоит особняком не из-за "
                f"расхождения программ: при большой дисперсии в точке наблюдения дисперсионный вклад "
                f"η·σδ в размер пучка сравним с бетатронным, и moments_output относит к «эмиттансу» полный "
                f"размер. Сам σx у точки наблюдения у всех инструментов совпадает.")
        R.p(vr.ISR_NOTE)
        R.p("Последний столбец таблицы — аналитика SynchrotronSim: εx, σδ и σz по интегралам излучения, "
            "без трекинга. Это прямая пара к twiss_output elegant.")
        R.diffs([
            ("аналитика ↔ аналитика", "εx (интегралы излучения)", vr.rel(ex[4], ex[3]),
             "одни и те же формулы по одной решётке"),
            ("аналитика ↔ аналитика", "σδ (интегралы излучения)",
             vr.rel(em["σδ, 10⁻⁴"][4], em["σδ, 10⁻⁴"][3]), "то же"),
            ("SynchrotronSim ↔ интегралы elegant", "εx", vr.rel(ex[0], ex[3]), "оба — точный расчёт"),
            ("SynchrotronSim ↔ moments_output", "σδ", vr.rel(em["σδ, 10⁻⁴"][0], em["σδ, 10⁻⁴"][2]), "то же"),
            ("SynchrotronSim ↔ moments_output", "σz", vr.rel(em["σz, мм"][0], em["σz, мм"][2]), "то же"),
            ("SynchrotronSim ↔ moments_output", "σx у точки наблюдения",
             vr.rel(em["σx у плоскости, мкм"][0], em["σx у плоскости, мкм"][2]), "то же"),
            ("SynchrotronSim ↔ moments_output", "εz", vr.rel(em["εz, мкм"][0], em["εz, мкм"][2]), "то же"),
            ("SynchrotronSim ↔ статистика elegant", "σδ", vr.rel(em["σδ, 10⁻⁴"][0], em["σδ, 10⁻⁴"][1]),
             "квантовый шум в трекинге elegant (см. текст)"),
        ])

    # ================================================================ 16. апертура
    _, grid = vr.read_csv(res("dynamic_aperture", "grid.csv")) if have("dynamic_aperture") else (None, [])
    if grid:
        nt = max(vr.num(r[2]) for r in grid)
        agree = sum((vr.num(r[2]) == nt) == (vr.num(r[3]) == nt) for r in grid)
        lost = sum(1 for r in grid if vr.num(r[2]) < nt)
        R.chapter("Динамическая апертура")
        R.meaning()
        R.p("Динамическая апертура — область начальных смещений, в которой частица остаётся в кольце "
            "заданное число оборотов. Её ограничивает нелинейность: на больших амплитудах движение "
            "становится хаотическим, и частица уходит. Апертура определяет эффективность инжекции и время "
            "жизни пучка; на её границе траектории чувствительны к малейшим различиям, поэтому это строгая "
            "проверка нелинейной части моделей.")
        R.tools([
            ("трекинг", "CppRing.track_tbt (частицы параллельно в C++) + survived_turns по сетке",
             "то же поверх трекинга elegant"),
            ("встроенный инструмент", "—", "frequency_map / find_aperture (справочно)"),
        ], "fma")
        R.results()
        fx, fy = M["fma_x_mm"], M["fma_y_mm"]
        R.p(f"Сетка: x от {P(fx[0], 3)} до {P(fx[1], 3)} мм ({int(fx[2])} "
            f"{vr.plural(fx[2], 'точка', 'точки', 'точек')}), y от {P(fy[0], 3)} до {P(fy[-1], 3)} мм, "
            f"δ = 0; {int(nt)} оборотов; частица потеряна, если |x| или |y| в точке наблюдения превысит "
            f"{P(float(M.get('loss_limit_m', 0.03)) * 1e3, 4)} мм.")
        R.fig(res("dynamic_aperture", "dynamic_aperture.png"), "Число прожитых оборотов по сетке.")
        lost_el = sum(1 for r in grid if vr.num(r[3]) < nt)
        R.p(f"Потеряно частиц: у SynchrotronSim {lost} из {len(grid)}, у elegant {lost_el}.")
        if M["source"]:
            R.p(vr.da_note(grid, nt))
        R.diffs([("трекинг ↔ трекинг", "исход «выжила/потеряна»", f"{agree} из {len(grid)}",
                  "совпадает во всех точках" if agree == len(grid) else "различия на хаотической границе")])

    # ================================================================ 17. карта частот и резонансы
    if grid:
        R.chapter("Карта частот и сетка резонансов")
        R.meaning()
        R.p("Для каждой частицы сетки по первой и второй половине прогона определяются частоты (νx, νy), "
            "и точка ставится на плоскость частот. Резонансные линии m·νx + n·νy = p — места, где "
            "нелинейные члены поля накапливаются от оборота к обороту; чем ниже порядок |m| + |n|, тем "
            "опаснее. Изменение частоты между половинами прогона — скорость диффузии — отмечает хаотические "
            "области. Карта показывает, как далеко след пучка от опасных резонансов.")
        R.formula(r"m\,\nu_x + n\,\nu_y = p,\qquad D = \log_{10}\sqrt{\left(\nu_x^{(2)}-\nu_x^{(1)}\right)^2"
                  r" + \left(\nu_y^{(2)}-\nu_y^{(1)}\right)^2}",
                  [
            (r"D", "показатель диффузии частоты", "—"),
            (r"\nu^{(1)},\ \nu^{(2)}", "частоты по первой и второй половине прогона", "—"),
            (r"m,\ n,\ p", "целые числа резонансного условия", "—"),
        ])
        R.tools([
            ("трекинг", "ring_analysis.frequency_map (NAFF по половинам) и resonance_lines (до 4-го порядка)",
             "то же поверх трекинга elegant; встроенный frequency_map — справочно"),
        ], "fma")
        R.results()
        R.fig(res("dynamic_aperture", "frequency_map.png"), "Карта частот на сетке резонансов до 4-го "
                                                            "порядка; врезка — след в масштабе.")
        R.fig(res("dynamic_aperture", "diffusion_xy.png"), "Скорость диффузии по сетке и разность частот "
                                                          "программ.")
        both = [r for r in grid if vr.num(r[2]) >= nt and vr.num(r[3]) >= nt
                and all(vr.num(r[i]) is not None and math.isfinite(vr.num(r[i])) for i in (4, 5, 8, 9))]
        dqf = lambda r: abs((vr.num(r[4]) - vr.num(r[5]) + 0.5) % 1.0 - 0.5)
        regular = [dqf(r) for r in both if max(vr.num(r[8]), vr.num(r[9])) < -6]
        chaotic = [dqf(r) for r in both if max(vr.num(r[8]), vr.num(r[9])) >= -6]
        dq = [dqf(r) for r in both]
        R.p("Рабочая точка и след частот обеих программ лежат в одном и том же месте диаграммы."
            + (" Нелинейность только кинематическая, поэтому на малых амплитудах след стягивается к рабочей "
               "точке и заметно растягивается лишь на сотнях миллиметров. У границы устойчивости часть "
               "частиц захватывается полуцелым резонансом 2νy = 13 (цепочка точек на νy = 0,5), дальше — "
               "хаотическая область с большой диффузией; картина у обеих программ одна и та же."
               if M["source"] else ""))
        rows_fma = [("трекинг ↔ трекинг", f"νx выживших в обеих программах ({len(dq)}), медиана",
                     float(np.median(dq)), "предел интегратора")]
        if regular:
            rows_fma.append(("трекинг ↔ трекинг", f"νx регулярных частиц (диффузия < 10⁻⁶, {len(regular)}), max",
                             float(np.max(regular)), "то же; растёт с амплитудой"))
        if chaotic:
            rows_fma.append(("трекинг ↔ трекинг", f"νx хаотических частиц ({len(chaotic)}), max",
                             float(np.max(chaotic)), "хаос: малое различие растёт экспоненциально"))
        R.diffs(rows_fma)

    # ================================================================ 18. матрица отклика
    if have("response"):
        rsp = vr.csv_dict(res("response", "summary.csv"))
        crc = vr.csv_dict(res("response", "correction.csv"))
        n_st = len(vr._stations_report())
        R.chapter("Матрица отклика орбиты")
        R.meaning()
        R.p("Матрица отклика — смещение замкнутой орбиты в каждом датчике на единичный кик каждого "
            "корректора. Это основа коррекции орбиты: обращение матрицы даёт кики, возвращающие орбиту на "
            "ось. Её элементы определяются β и фазами в местах датчика и корректора, поэтому матрица "
            "отклика одновременно проверяет всю линейную оптику. Та же машинерия даёт чувствительность "
            "частот и β к ошибкам градиента — допуска на магниты.")
        R.formula(r"R_{bk} = \frac{\sqrt{\beta_b\beta_k}}{2\sin\pi\nu}\cos\left(\left|\psi_b-\psi_k\right|"
                  r"-\pi\nu\right),\qquad \Delta\nu = \pm\frac{\beta}{4\pi}\,\Delta(K_1L)",
                  [
            (r"R_{bk}", "отклик орбиты на датчике b от кика корректора k", "м/рад"),
            (r"\beta_b,\ \beta_k", "бета-функции в местах датчика и корректора", "м"),
            (r"\psi_b,\ \psi_k", "фазовые набеги до датчика и до корректора", "рад"),
            (r"\nu", "бетатронная частота кольца", "—"),
            (r"\Delta(K_1L)", "ошибка интегрированного градиента", "м⁻¹"),
            (r"\beta", "бета-функция в месте магнита", "м"),
        ])
        R.tools([
            ("аналитика", "ring_response.orbit_response: кик тонкий, в середине корректора",
             "correction_matrix_output (функции Твисса усреднены по концам корректора)"),
            ("трекинг", "кик в одном корректоре (CppRing.with_kick) и замкнутая орбита C++-ядра",
             "кик через alter_elements и замкнутая орбита трекинга elegant"),
            ("коррекция и допуска", "ring_response.svd_solve, kick_to_currents, tolerances", "—"),
        ], "response")
        R.results()
        n_cr = int(crc["оставлено сингулярных чисел"][1])
        R.p(f"Матрица {n_st} {vr.plural(n_st, 'датчик', 'датчика', 'датчиков')} × {n_cr} "
            f"{vr.plural(n_cr, 'корректор', 'корректора', 'корректоров')}."
            + (" Корректоры — пробные, в пустых промежутках SX и SZ (см. главу 1)." if M["source"] else ""))
        R.fig(res("response", "orbit_response.png"), "Отклик орбиты на кик первого корректора и разность "
                                                     "с elegant.")
        R.fig(res("response", "svd_correction.png"), "Спектр матрицы отклика и коррекция орбиты по SVD.")
        crow = []
        for k, v in crc.items():
            if k.startswith("максимум орбиты"):
                crow.append(("максимум орбиты до коррекции, мкм", v[0]))
                crow.append(("максимум орбиты после коррекции, мкм", v[1]))
            elif k.startswith("оставлено"):
                crow.append(("сингулярных чисел: оставлено / всего", f"{int(v[0])} / {int(v[1])}"))
            else:
                crow.append((k, v[0]))
        R.p("Опыт коррекции: 20 случайных квадруполей развёрнуты со СКО 50 мкрад, орбита правится по SVD, "
            "кик переводится в поле по проектной жёсткости и в токи обмоток линейной матрицей корректора.")
        R.table(["величина", "значение"], crow, digits=[0, 6],
            title="Коррекция орбиты по SVD после разворота квадруполей")
        R.fig(res("response", "tolerances.png"), "Чувствительность к ошибке градиента 0,1 % по типам "
                                                 "устройств.")
        R.csv_table(res("response", "tolerances.csv"), 5, title="Допуска на ошибку градиента по типам магнитов")
        med = [v[2] for k, v in rsp.items() if "медианное" in k]
        R.diffs(
            [("аналитика ↔ аналитика", k.split(":")[0] + ": медиана по крупным членам", v[2],
              "у elegant функции Твисса корректора усреднены по концам") for k, v in rsp.items()
             if "медианное" in k]
            + [("аналитика ↔ трекинг SynchrotronSim",
                "орбита от кика на " + k.split("кик на ")[-1].replace(": орбита, мкм", "").replace(".", ","),
                v[2], "модель тонкого кика в середине") for k, v in rsp.items()
               if k.startswith("аналитика ↔ трекинг")]
            + [("трекинг ↔ трекинг",
                "орбита от кика на " + k.split("кик на ")[-1].replace(": орбита, мкм", "").replace(".", ","),
                v[2], "предел интегратора") for k, v in rsp.items() if k.startswith("трекинг ↔ трекинг")])

    out = out or os.path.join(ed.ROOT, "workspace", "results", f"verify_elegant_params_{vr.MACHINE}.docx")
    R.save(out)
    return out


def main():
    ap = argparse.ArgumentParser(description="отчёт сверки SynchrotronSim с elegant: глава на параметр")
    ap.add_argument("--machine", default="siberia", choices=sorted(ed.MACHINES))
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    vr.use_machine(a.machine)
    print(build(a.out))


if __name__ == "__main__":
    main()
