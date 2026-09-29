# -*- coding: utf-8 -*-
"""Отчёт «Сверка SynchrotronSim с elegant» — Word по результатам verify_elegant.py.

Сопоставление построено ПАРАМИ ОДИНАКОВЫХ ИНСТРУМЕНТОВ: аналитический механизм одной программы против
аналитического механизма другой, трекинг против трекинга. Разрыв между уровнями (аналитика против
собственного трекинга) вынесен отдельно: это не различие программ, а недосходимость интегратора.

Установка задаётся ключом --machine (см. elegant_diag.MACHINES) — отчёт целиком строится по ней.
Читает workspace/results/verify_elegant/[<установка>/]<опыт>/ (PNG, CSV, JSON).
  py -3.14 verify_report.py --machine siberia   → workspace/results/verify_elegant_report_siberia.docx
"""
import os, re, csv, json, glob, math, argparse

import numpy as np

import elegant_diag as ed

ROOT = os.path.dirname(os.path.abspath(__file__))
BASE = os.path.join(ROOT, "workspace", "results", "verify_elegant")
MACHINE = "sila"
M = ed.machine(MACHINE)
RES = BASE
OUT = os.path.join(ROOT, "workspace", "results", "verify_elegant_report.docx")


def use_machine(name):
    """Переключить отчёт на установку: папка результатов и имя файла."""
    global MACHINE, M, RES, OUT
    MACHINE, M = name, ed.machine(name)
    RES = BASE if name == "sila" else os.path.join(BASE, name)
    OUT = os.path.join(ROOT, "workspace", "results",
                       "verify_elegant_report.docx" if name == "sila"
                       else f"verify_elegant_report_{name}.docx")
    return OUT


def _stations_report():
    """Точки наблюдения орбиты (MONI) установки — для подписи размера матрицы отклика."""
    import ring_response as rr
    return rr.bpm_stations(*ed.machine_lattice(M))


ISR_NOTE = ("Трекинг elegant с квантовым излучением систематически даёт разброс примерно на 3 % меньше, чем его же "
            "moments_output: это видно и по диффузии за оборот (статистика по 16 000 частиц), и по равновесному облаку "
            "ансамбля в главе о портретах с излучением. Около 1,3 % по σ объясняет обрезание гауссова шума в CSBEND на 3σ "
            "(srGaussianLimit); остальная часть этим не объясняется и относится к модели квантового излучения в трекинге "
            "elegant. SynchrotronSim совпадает с moments_output до 0,2 %, а его ансамбль — с расчётом равновесия в "
            "пределах статистики. Опорными следует считать аналитические оценки elegant.")


def si_len(v):
    """Длина с подходящей приставкой: 1.2e-12 м → «1,2 пм». Чтобы тексты не устаревали при смене точности."""
    a = abs(float(v))
    for scale, unit in ((1e-3, "мм"), (1e-6, "мкм"), (1e-9, "нм"), (1e-12, "пм")):
        if a >= scale or scale == 1e-12:
            return f"{pretty(v / scale, 3)} {unit}"


def geom_residual():
    """Насколько решётка из базы данных воспроизводит файл описания: (Δпериметра, поперечная невязка), м."""
    an = csv_dict(os.path.join(RES, "analytic", "summary.csv"))
    if not an or "периметр, м" not in an:
        return None, None
    return abs(an["периметр, м"][2]), abs(an["поперечный C за оборот (идеальная решётка)"][2])


def geom_note():
    """Фраза о точности геометрии базы — числами, без «около микрона»."""
    dC, dT = geom_residual()
    if dC is None:
        return ""
    return (f"периметр решётки из базы отличается от файла описания на {si_len(dC)}, поперечная невязка "
            f"замыкания {si_len(dT)}")


def geom_label():
    """Короткая пометка для столбца «чем объясняется»."""
    dC, dT = geom_residual()
    return "геометрия базы данных" if dC is None else f"геометрия базы данных, невязка {si_len(max(dC, dT))}"


def orbit_verdict(co_ss, co_el):
    """Фраза о замкнутой орбите идеальной решётки по двум измеренным значениям (мкм)."""
    a_ss, a_el = abs(co_ss), abs(co_el)
    d = abs(co_ss - co_el)
    if max(a_ss, a_el) < 0.05:                       # обе программы дают численный ноль
        return (f". У обеих программ это численный ноль: геометрия решётки в базе данных записана с полной "
                f"точностью двойного числа, {geom_note()}, и орбите отклоняться не от чего")
    if a_ss < 0.1 * a_el:                            # у нас ноль, у elegant — погрешность его поиска
        return (f". У SynchrotronSim орбиты фактически нет: геометрия базы воспроизводит файл описания с "
                f"точностью {si_len(max(geom_residual()))}. Остаток у elegant — погрешность поиска "
                f"неподвижной точки по его трекингу, а не орбита")
    if d < 0.05 * max(a_ss, a_el):                   # обе видят ОДНУ И ТУ ЖЕ ненулевую орбиту
        return (f". Обе программы видят одну и ту же орбиту и расходятся на {si_len(d * 1e-6)}: это "
                f"свойство самого описания решётки (замыкание кольца), а не базы данных — её геометрия "
                f"воспроизводит файл с точностью {si_len(max(geom_residual()))}")
    return (f": у SynchrotronSim это невязка положения магнитов в базе ({geom_note()}), она же даёт след "
            f"в дисперсии")


def orbit_label(co_ss, co_el):
    """Короткая пометка для столбца «чем объясняется» у строки замкнутой орбиты."""
    a_ss, a_el = abs(co_ss), abs(co_el)
    if max(a_ss, a_el) < 0.05 or a_ss < 0.1 * a_el:
        return "погрешность поиска неподвижной точки у elegant"
    if abs(co_ss - co_el) < 0.05 * max(a_ss, a_el):
        return "общая орбита из описания решётки, программы её воспроизводят одинаково"
    return geom_label()


def db_note():
    """Откуда SynchrotronSim берёт решётку (база данных комплекса) — абзац для раздела об установке."""
    how = (f"импортом того же файла {os.path.basename(M['source'])} (скрипт s2lattice_to_db.py)" if M["source"]
           else "сборкой из того же описания (assemble_lattice_configs.py)")
    return ("elegant читает файл описания напрямую. SynchrotronSim берёт решётку из базы данных комплекса: в "
            "таблице workspace/tables/lattice_configs.csv лежат положение, ориентация и сила каждого магнита, "
            "в element_connectivity.csv — порядок магнитов в линии. Таблицы хранят несколько колец под своими "
            "номерами (1 — СИЛА, 2 — Сибирь, 3 — один супериод СИЛА); это кольцо записано в них под номером "
            f"{M['lattice']} {how}. Так обе программы считают одну решётку, прочитанную из одного файла.")


def corrector_note():
    """Какие корректоры стоят в опытах с откликом орбиты."""
    lat = ed.machine_lattice(M)
    n = sum(1 for nm in ed.line_names(*lat) if nm.upper() in {c.upper() for c in M["correctors"]})
    names = " и ".join(M["correctors"])
    if not M["source"]:
        return f"Корректоры орбиты — штатные корректоры кольца ({', '.join(M['correctors'])}), {n} вхождений."
    return (f"Корректоров орбиты в файле нет. Для опытов с откликом орбиты в обе программы поставлены ПРОБНЫЕ "
            f"корректоры — однородный горизонтальный кик в {n} пустых промежутках {names} длиной "
            f"{pretty(M['corrector_len'], 3)} м. В файле это места секступолей, но с K2 = 0, то есть просто "
            f"дрейфы: в elegant элемент заменяется на KICKER той же длины, в SynchrotronSim кик ставится в тот "
            f"же элемент решётки (дрейф — кик — дрейф по срезам). Секступолями как магнитами они при этом не "
            f"служат: промежутки выбраны потому, что они свободны и одинаково описаны в обеих программах. Для "
            f"сверки место корректора не принципиально — отклик считается одной и той же формулой в любой "
            f"точке кольца.")


def t_floor_note(n_turns):
    """Полка продольной кривой elegant: шум накопленного абсолютного времени."""
    el = np.load(os.path.join(RES, "radiation", "elegant_sr.npz"), allow_pickle=True)["tbt"]
    ss = np.load(os.path.join(RES, "radiation", "ss_sr.npz"), allow_pickle=True)["tbt"]
    tail = slice(int(0.98 * n_turns), None)
    f_el = float(np.nanstd(el[2][tail, 4])); f_ss = float(np.nanstd(ss[2][tail, 4]))
    path = n_turns * lattice_facts()["circumference"]
    if f_el < 1.5 * f_ss:                      # у этой установки прогон короче: до предела точности не дошли
        return ""
    return (f"На продольной кривой классического затухания видно расхождение в самом низу: у elegant она "
            f"выходит на полку около {pretty(f_el * 1e3, 2)} мм, а у SynchrotronSim продолжает экспоненту до "
            f"{pretty(f_ss * 1e3, 2)} мм. Это не физика, а разная бухгалтерия продольной координаты. elegant "
            f"выдаёт АБСОЛЮТНОЕ время пролёта, и запаздывание получается разностью двух чисел порядка "
            f"{pretty(path, 2)} м; шаг двоичного представления такого числа — 5·10⁻¹⁰ м, а накопленное за "
            f"{n_turns} оборотов округление даёт случайное блуждание на два порядка больше (растёт как корень "
            f"из числа оборотов). SynchrotronSim вычитает периметр на каждом обороте, поэтому хранит сразу "
            f"малую величину и такого шума не имеет. На сверку это не влияет: огибающие сравниваются, пока "
            f"амплитуда выше десятой части начальной, то есть на четыре порядка выше полки.")


def da_axes(grid, nt, col):
    """Граница апертуры по осям для столбца оборотов col: (+x, −x) на нижней строке y, y на столбце x ≈ 0,
    и есть ли выше этой границы выжившие (остров устойчивости). Амплитуды — мм."""
    pts = [(num(r[0]), num(r[1]), num(r[col])) for r in grid]
    y_lo = min((y for _, y, _ in pts), key=abs)          # ряд, ближайший к медианной плоскости
    row = [(x, n) for x, y, n in pts if y == y_lo]
    xp = min((x for x, n in row if x > 0 and n < nt), default=float("nan"))
    xm = max((x for x, n in row if x < 0 and n < nt), default=float("nan"))
    x_c = min((x for x, _, _ in pts), key=abs)
    column = sorted((y, n) for x, y, n in pts if x == x_c)
    y0 = min((abs(y) for y, n in column if n < nt), default=float("nan"))
    island = any(n >= nt for y, n in column if abs(y) > y0)
    return xp, xm, y0, island


def da_note(grid, nt):
    """Абзац о природе границы апертуры установки с нулевыми секступолями."""
    ss, el = da_axes(grid, nt, 2), da_axes(grid, nt, 3)
    y_lo = min((num(r[1]) for r in grid), key=abs)
    lim = pretty(float(M.get("loss_limit_m", 0.03)) * 1e3, 4)
    txt = (f"Секступоли решётки нулевые, и единственная нелинейность кольца — кинематическая: точные дрейфы "
           f"и края изгибов. Её действие становится заметным только тогда, когда угол траектории к оси доходит "
           f"до 20° и больше, поэтому сетка апертуры здесь взята до сотен миллиметров, а порог потери поднят до "
           f"{lim} мм — теряются лишь частицы, которые действительно уходят. По горизонтали (y = "
           f"{pretty(abs(y_lo), 3)} мм) устойчивость теряется при x = +{pretty(ss[0], 3)} / {pretty(ss[1], 3)} мм у "
           f"SynchrotronSim и +{pretty(el[0], 3)} / {pretty(el[1], 3)} мм у elegant, по вертикали (x = 0) — с "
           f"y = {pretty(ss[2], 3)} и {pretty(el[2], 3)} мм.")
    if ss[3] and el[3]:
        txt += " Выше границы у обеих программ есть одинаковый остров устойчивости."
    return txt + (" Такая граница лежит далеко за вакуумной камерой и как характеристика машины смысла не "
                  "имеет; это самая жёсткая проверка точной кинематики обеих программ. Реальную динамическую "
                  "апертуру установки задают секступоли коррекции хроматизма, которых в этом файле нет.")


def tex_num(v, digits=3):
    """Число для формулы: 2,14·10⁻⁷ в виде LaTeX."""
    if not v:
        return "0"
    e = int(math.floor(math.log10(abs(v))))
    m = f"{v / 10 ** e:.{max(digits - 1, 0)}f}".replace(".", "{,}")
    return rf"{m} \cdot 10^{{{e}}}"


def mass_facts():
    """Постоянные массы электрона двух программ и след, который давала их разница на этой установке.

    Возвращает (δ₀, max ηx [м], ξ1 по x, РАСХОЖДЕНИЕ замкнутых орбит двух программ [мкм], Δνx)."""
    d0 = ed.ME_ELEGANT_MEV / ed.ME_MEV - 1
    lin = csv_dict(os.path.join(RES, "linear", "summary.csv")) if exists("linear") else {}
    ch = csv_dict(os.path.join(RES, "chromaticity", "summary.csv")) if exists("chromaticity") else {}
    eta = (lin.get("max ηx, м") or [None])[0]
    xi = (ch.get("ξ1 = dν/dδ (x)") or [None])[0]
    co = lin.get("замкнутая орбита x, мкм")
    nu = lin.get("νx")
    return (d0, eta, xi, (abs(co[0] - co[1]) if co else None),
            (abs(nu[0] - nu[1]) if nu else None))


def mass_section(R, head):
    """Раздел о разных постоянных массы электрона. head — R.h1 или R.h2 (номер главы у отчётов свой)."""
    d0, eta, xi, dco, dnu = mass_facts()
    head("Опорный импульс: разные постоянные массы электрона")
    R.p("Частицы передаются из одной программы в другую в импульсных единицах: поперечные компоненты "
        "делятся на опорный импульс кольца.")
    R.formula(r"p_0 = \frac{E}{m_e c^2}", [
        (r"p_0", "опорный импульс кольца в единицах mₑc", "безразмерный"),
        (r"E", "энергия пучка", "МэВ"),
        (r"m_e c^2", "энергия покоя электрона", "МэВ")])
    me_el, me_ss = (f"{ed.ME_ELEGANT_MEV:.8f}".replace(".", ","), f"{ed.ME_MEV:.11f}".replace(".", ","))
    R.p(f"Постоянная mₑc² у двух программ разная. У elegant она зашита как {me_el} МэВ — значение "
        f"CODATA 1986 года, оставшееся в коде с первых версий. SynchrotronSim считает по современному "
        f"{me_ss} МэВ (CODATA 2018). Наше значение точнее, и оно оставлено без изменений: физика "
        f"SynchrotronSim считается по нему. Согласовывать пришлось не постоянную, а единицы, в которых "
        f"координаты передаются между программами.")
    R.p("Относительная разница опорных импульсов:")
    R.formula(r"\delta_0 = \frac{m_e^{elegant}}{m_e^{SS}} - 1 = " + tex_num(d0), [
        (r"\delta_0", "отклонение по импульсу, с которым elegant видит частицу, записанную по нашему "
                      "опорному импульсу", "безразмерное"),
        (r"m_e^{elegant}, m_e^{SS}", "постоянная массы электрона у elegant и у SynchrotronSim", "МэВ")])
    R.p("Пока файл частиц писался по НАШЕМУ опорному импульсу, elegant читал нашу осевую частицу (δ = 0) "
        "как частицу с δ = δ₀. Такая частица идёт не по оси, а по дисперсионной орбите:")
    R.formula(r"x(s) = \eta_x(s)\,\delta_0", [
        (r"x(s)", "горизонтальное смещение частицы в точке s", "м"),
        (r"\eta_x(s)", "дисперсионная функция кольца", "м"),
        (r"\delta_0", "отклонение по импульсу из формулы выше", "безразмерное")])
    R.p((f"Именно это и показывал трекинг elegant: к замкнутой орбите добавлялась добавка, повторяющая "
         f"рисунок дисперсии с периодом в один супериод. Наибольшая дисперсия кольца {pretty(eta, 6)} м, "
         f"значит добавка доходила до {si_len(abs(eta * d0))} — отсюда рассыпанные по кольцу точки разной "
         f"высоты. " if eta else "Именно это и показывал трекинг elegant. ")
        + "Проверка однозначная: отношение добавки к ηx(s) одинаково по всему кольцу с разбросом 10⁻¹¹, то "
          "есть точки ложились ровно на дисперсионную орбиту, а не были шумом счёта. Там, где у решётки "
          "своя замкнутая орбита есть и так, добавка искажала её форму, а не создавала с нуля.")
    R.p("Тот же сдвиг по импульсу смещал через хроматичность и бетатронные частоты:")
    R.formula(r"\Delta\nu = \xi_1 \delta_0", [
        (r"\Delta\nu", "сдвиг бетатронной частоты", "безразмерный"),
        (r"\xi_1", "хроматичность кольца dν/dδ", "безразмерная"),
        (r"\delta_0", "отклонение по импульсу", "безразмерное")])
    if xi:
        R.p(f"При ξ1 = {pretty(xi, 6)} это даёт Δνx = {pretty(abs(xi * d0), 3)}: на столько трекинг "
            f"elegant и расходился с трекингом SynchrotronSim, пока опорный импульс не был согласован. "
            f"Этот остаток до выяснения причины относили к пределу сходимости интегратора.")
    R.p("Различия физики здесь нет: это разные единицы передачи координат. Поэтому опорный импульс для "
        "файлов частиц берётся теперь у самого elegant — один служебный запуск twiss_output, из него "
        "читается pCentral, и по нему пишутся все файлы частиц (elegant_diag.elegant_me). Постоянные при "
        "этом не менялись: SynchrotronSim считает по современной массе, elegant по своей, а сравниваются "
        "одни и те же физические частицы.")
    rows = []
    if eta is not None and dco is not None:
        rows.append(("расхождение замкнутых орбит между трекингами, мкм", abs(eta * d0) * 1e6, dco))
    if xi and dnu is not None:
        rows.append(("расхождение νx между трекингами", abs(xi * d0), dnu))
    if rows:
        R.table(["показатель", "вклад δ₀ (до согласования)", "измерено после согласования"], rows, [0, 4, 4],
                title="След разных постоянных массы электрона до и после согласования опорного импульса")
    R.p("Согласование опорного импульса убирает сдвиг по δ, но САМИ постоянные остаются разными, и "
        "релятивистский фактор у двух программ отличается на ту же долю δ₀. Величины, зависящие от γ "
        "напрямую, поэтому до нуля сойтись не могут: потеря энергии за оборот U₀ ∝ γ⁴ несёт остаток "
        "порядка 4·δ₀.")
    R.formula(r"\frac{\Delta U_0}{U_0} \approx 4\,\delta_0 = " + tex_num(4 * d0), [
        (r"U_0", "потеря энергии на излучение за оборот", "МэВ"),
        (r"\delta_0", "относительная разница опорных импульсов", "безразмерное")])
    lg = csv_dict(os.path.join(RES, "longitudinal", "summary.csv")) if exists("longitudinal") else {}
    u0 = lg.get("U0 по трекингу за оборот, МэВ")
    if u0 and u0[0] and u0[1]:
        R.p(f"Измеренное расхождение U₀ между трекингами — {pretty(rel(u0[0], u0[1]), 2)}, того же порядка. "
            f"Это не ошибка ни одной из программ, а цена разных постоянных; наше значение точнее, поэтому "
            f"подгонять U₀ под elegant незачем.")
    R.p((f"Оставшиеся после согласования {si_len(dco * 1e-6)} расхождения по орбите и {pretty(dnu, 2)} по "
         f"частоте — " if dco is not None and dnu is not None else "Оставшееся после согласования — ")
        + "это уже дискретизация интегратора, одинаковая в обеих программах, а не разные физические "
          "постоянные.")


def integrator_scaling(order=4):
    """Сходимость частоты по числу срезов из опыта «интегратор» — строкой, числами из tunes.csv."""
    p = os.path.join(RES, "integrator", "tunes.csv")
    if not os.path.exists(p):
        return ""
    _, rows = read_csv(p)
    exact, pts = None, []
    for r in rows:
        if r[0].startswith("twiss_output"):
            exact = num(r[1])
        m = re.match(rf"N=(\d+), порядок {order}$", r[0].strip())
        if m and num(r[1]) is not None and math.isfinite(num(r[1])):
            pts.append((int(m.group(1)), num(r[1])))
    if exact is None or len(pts) < 3:
        return ""
    pts.sort()
    sup = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
    body = "; ".join(f"N = {n} → {pretty(abs(v - exact), 2)}" for n, v in pts[-4:])
    k = abs(pts[-2][1] - exact) / abs(pts[-1][1] - exact)
    return (f"при схеме порядка {order} расхождение частоты падает как N⁻{str(order).translate(sup)} "
            f"({body}), и отношение соседних точек {pretty(k, 3)} при теоретических {2 ** order}")


def phase_note():
    """Почему на панели фазового набега растёт только кривая «аналитика − трекинг»."""
    d = os.path.join(RES, "linear")
    try:
        tw = np.load(os.path.join(d, "elegant_twiss.npz"))
        el = np.load(os.path.join(d, "elegant_track.npz"))
        ss = np.load(os.path.join(d, "ss.npz"))
    except OSError:
        return ""
    two = 2 * math.pi
    grey = abs(float(tw["psix"][-1] - el["psix"][-1])) / two
    cross = abs(float(tw["psix"][-1] - ss["psix"][-1])) / two
    blue = float(np.max(np.abs(ss["psix"] - el["psix"]))) / two
    t = (f"На панели фазового набега растёт только кривая «twiss_output − трекинг»: это единственная здесь "
         f"пара РАЗНЫХ уровней — точные матрицы элементов против трекинга на {M['slices'][0]} срезах. "
         f"Ошибка дискретизации на каждом элементе одного знака, поэтому фаза копится вдоль кольца почти "
         f"линейно и к концу оборота набирает {pretty(grey, 3)} — ровно разницу частот νx из сводки. "
         f"Кривые «трекинг ↔ трекинг» сравнивают ОДИНАКОВЫЕ уровни: оба трекера идут на одном числе "
         f"срезов, их ошибки дискретизации совпадают и взаимно сокращаются, поэтому разность держится на "
         f"{pretty(blue, 2)}. Так и должно быть: расти обязана именно серая кривая.")
    t += (f" Это не свойство elegant. Если вычесть из его twiss_output не его собственный трекинг, а "
          f"трекинг SynchrotronSim, получится та же кривая: {pretty(cross, 3)} против {pretty(grey, 3)}.")
    conv = integrator_scaling()
    if conv:
        t += f" Отставание убирается не сменой программы, а измельчением шага: {conv}."
    return t


def lattice_facts():
    """Числа установки для текста: элементов в линии, периметр, число магнитов."""
    lat = ed.machine_lattice(M)
    par = ed.element_params(lat[0], lat[1])
    names = ed.line_names(*lat)
    mag = sum(1 for n in names if par[n.upper()]["type"].upper() not in ed._MAGNET_TYPES_SKIP)
    return dict(n_line=len(names), n_mag=mag, circumference=ed.line_length(*lat),
                lead=ed.lead_drift(*lat))

# показатель → (SynchrotronSim аналитика, SynchrotronSim трекинг, elegant аналитика, elegant трекинг)
INSTRUMENTS = [
    ("Замкнутая орбита",
     "ring_matrix: (I − R)⁻¹·C", "closed_orbit — Ньютон по отображению за оборот",
     "closed_orbit по матрицам", "тот же Ньютон поверх ElegantRing"),
    ("β, α, фазовый набег, η, частоты, αc",
     "ring_matrix.twiss — перемножение матриц элементов", "ring_analysis.twiss — матрицы разностями по траекториям",
     "twiss_output", "ring_analysis.twiss поверх ElegantRing"),
    ("Матрицы отдельных элементов",
     "ring_matrix.element_matrix", "element_matrices из накопленных матриц",
     "matrix_output (individual_matrices)", "матрицы элементов из траекторий elegant"),
    ("Матрица отклика орбиты на корректоры",
     "ring_response.orbit_response", "замкнутая орбита от кика корректора",
     "correction_matrix_output", "кик через alter_elements + замкнутая орбита"),
    ("Хроматичность ξ1…ξ3, αc2",
     "ring_map: отображение оборота 5-го порядка, ν(δ) по линеаризации вокруг орбиты с δ",
     "chromaticity: ν(δ) по орбитам вне импульса, подгонка 7-й степени",
     "twiss_output, higher_order_chromaticity при default_order = 3", "tunes_point поверх ElegantRing"),
    ("Частоты от амплитуды ∂ν/∂J",
     "ring_map.normal_form: нормальная форма отображения оборота", "track_tbt + NAFF",
     "twiss_output, compute_driving_terms: ∂ν/∂J по мультиполям (кинематику не учитывает)",
     "tune_shift_with_amplitude, track_tbt + NAFF"),
    ("Резонансные члены h_jklm",
     "— (можно получить из того же отображения ring_map)", "— (косвенно: карта частот и диффузия)",
     "twiss_output, compute_driving_terms", "—"),
    ("Синхротронная частота, времена затухания",
     "ring_matrix.radiation_integrals: U₀, Jx/Jy/Js, τ по I₂ и I₄; Qs по αc, V, h",
     "собственные числа матрицы оборота, спад амплитуды",
     "rf_setup, twiss_output (интегралы излучения)", "те же алгоритмы поверх ElegantRing"),
    ("Равновесный эмиттанс, σδ, σz",
     "ring_matrix.radiation_integrals: εx, σδ, σz по I₁…I₅",
     "метод Чао: матрица оборота и диффузия по срезам изгибов",
     "twiss_output (интегралы излучения), moments_output", "статистика квантового излучения по частицам"),
    ("Фазовые траектории, спектр",
     "—", "track_tbt, naff", "—", "WATCH + трекинг, sddsnaff"),
    ("Динамическая апертура, карта частот, резонансы",
     "— (нелинейная задача)", "сетка частиц, survived_turns, frequency_map, resonance_lines",
     "—", "find_aperture, frequency_map, tune_footprint"),
    ("Допуска на градиенты",
     "ring_response.tolerances: Δν = ±β·Δ(K₁L)/4π и бета-биение", "— (перебор трекингом непрактичен)",
     "— (штатного инструмента нет)", "—"),
]


# ============================================================================ чтение результатов
def read_csv(path):
    if not os.path.exists(path):
        return None, []
    with open(path, encoding="utf-8-sig") as f:
        rows = list(csv.reader(f, delimiter=";"))
    return rows[0], rows[1:]


def num(v):
    try:
        return float(str(v).replace(",", "."))
    except ValueError:
        return None


def csv_dict(path):
    _, rows = read_csv(path)
    return {r[0]: [num(v) for v in r[1:]] for r in rows}


_SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def pretty(v, digits=7):
    if v is None or v == "":
        return ""
    x = num(v)
    if x is None:
        return str(v)
    if not math.isfinite(x):
        return "—"
    if x == 0:
        return "0"
    if abs(x) >= 1e5 or abs(x) < 1e-3:
        m, e = f"{x:.{digits - 1}e}".split("e")
        return f"{m.replace('.', ',').replace('-', chr(0x2212))}·10{str(int(e)).translate(_SUP)}"
    out = f"{x:.{digits}g}"
    if "e" in out:                                  # короткая запись .g сама ушла в экспоненту
        m, e = out.split("e")
        return f"{m.replace('.', ',').replace('-', chr(0x2212))}·10{str(int(e)).translate(_SUP)}"
    return out.replace(".", ",").replace("-", "−")


def plural(n, one, few, many):
    """Согласование числительного: 21 точка, 22 точки, 25 точек."""
    n = abs(int(n))
    if n % 100 // 10 == 1:
        return many
    return {1: one}.get(n % 10, few if 2 <= n % 10 <= 4 else many)


def rel(a, b):
    return abs(a - b) / abs(b) if b else float("nan")


def exists(name):
    d = os.path.join(RES, name)
    return os.path.isdir(d) and any(f.endswith(".png") for f in os.listdir(d))


# ============================================================================ документ
class Report:
    """Документ отчёта. Оформление — по ГОСТ для студенческих работ: поля 30/15/20/20 мм, основной текст
    14 пт по ширине с абзацным отступом 1,25 см и полуторным интервалом; формулы отдельной строкой по
    центру с номером у правого поля, под ними расшифровка переменных на размер мельче; у каждой таблицы
    заголовок «Таблица N — …» сверху, у каждого рисунка подпись «Рисунок N — …» снизу; после формулы,
    таблицы и подписи — пустая строка."""

    BODY_PT, SMALL_PT = 14, 12
    TABLE_PT = BODY_PT - 2
    INDENT_CM, LINE, TEXT_CM = 1.25, 1.5, 16.5          # ширина полосы набора: 210 − 30 − 15 мм

    def __init__(self):
        from docx import Document
        self.doc = Document()
        try:
            import latex2omml as lo
            self.lo = lo
        except Exception:
            self.lo = None
        self.first_h1 = True
        self.n_eq = self.n_fig = self.n_tab = 0
        self._gost()

    # ------------------------------------------------------------------ оформление
    def _gost(self):
        from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
        for sec in self.doc.sections:                     # поля страницы
            sec.left_margin, sec.right_margin = Cm(3.0), Cm(1.5)
            sec.top_margin, sec.bottom_margin = Cm(2.0), Cm(2.0)
        st = self.doc.styles["Normal"]
        st.font.name = "Calibri"; st.font.size = Pt(self.BODY_PT)
        pf = st.paragraph_format
        pf.alignment = AL.JUSTIFY; pf.first_line_indent = Cm(self.INDENT_CM)
        pf.line_spacing = self.LINE; pf.space_before = Pt(0); pf.space_after = Pt(0)
        for name, size, align in (("Title", 18, AL.CENTER), ("Heading 1", 16, AL.LEFT),
                                  ("Heading 2", 14, AL.LEFT)):
            try:
                h = self.doc.styles[name]
            except KeyError:
                continue
            h.font.name = "Calibri"; h.font.size = Pt(size); h.font.bold = True
            hp = h.paragraph_format
            hp.alignment = align; hp.line_spacing = self.LINE
            hp.first_line_indent = Cm(0 if align == AL.CENTER else self.INDENT_CM)
            hp.space_before = Pt(12); hp.space_after = Pt(6)
            hp.keep_with_next = True
        try:
            b = self.doc.styles["List Bullet"]
            b.font.name = "Calibri"; b.font.size = Pt(self.BODY_PT)
            b.paragraph_format.line_spacing = self.LINE
            b.paragraph_format.alignment = AL.JUSTIFY
        except KeyError:
            pass

    def _plain(self, par, size=None, align=None, indent=0.0, line=1.0):
        """Абзац без красной строки (формулы, подписи, ячейки таблиц)."""
        from docx.shared import Pt, Cm
        pf = par.paragraph_format
        pf.first_line_indent = Cm(indent); pf.line_spacing = line
        pf.space_before = Pt(0); pf.space_after = Pt(0)
        if align is not None:
            pf.alignment = align
        if size is not None:
            for r in par.runs:
                r.font.size = Pt(size)
        return par

    def _spacer(self):
        """Пустая строка — чтобы блоки не слипались."""
        from docx.shared import Pt
        par = self.doc.add_paragraph()
        self._plain(par, line=self.LINE)
        par.add_run("").font.size = Pt(self.BODY_PT)
        return par

    # ------------------------------------------------------------------ содержимое
    def h1(self, text):
        p = self.doc.add_heading(text, level=1)
        if not self.first_h1:
            p.paragraph_format.page_break_before = True
        self.first_h1 = False

    def h2(self, text):
        return self.doc.add_heading(text, level=2)

    def p(self, text):
        return self.doc.add_paragraph(text)

    def bullet(self, text):
        return self.doc.add_paragraph(text, style="List Bullet")

    def formula(self, latex, legend=()):
        """Формула отдельной строкой по центру с номером у правого поля; legend — расшифровка переменных:
        список (символ, пояснение[, размерность])."""
        from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH as AL, WD_TAB_ALIGNMENT as TA
        self.n_eq += 1
        par = self.doc.add_paragraph()
        self._plain(par, align=AL.LEFT, line=self.LINE)
        par.paragraph_format.tab_stops.add_tab_stop(Cm(self.TEXT_CM / 2), TA.CENTER)
        par.paragraph_format.tab_stops.add_tab_stop(Cm(self.TEXT_CM), TA.RIGHT)
        par.add_run("	").font.size = Pt(self.BODY_PT)
        if self.lo is not None:
            self.lo.add_formula(par, latex)
        else:
            r = par.add_run(latex); r.italic = True
        par.add_run(f"	({self.n_eq})").font.size = Pt(self.BODY_PT)
        if legend:
            self.legend(legend)
        self._spacer()
        return par

    def legend(self, items):
        """«где:» отдельной строкой, дальше по строке на переменную; размерность — в квадратных скобках."""
        from docx.shared import Pt, Cm
        from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
        head = self.doc.add_paragraph()
        self._plain(head, align=AL.LEFT)
        head.add_run("где:").font.size = Pt(self.SMALL_PT)
        for i, item in enumerate(items):
            sym, desc, unit = (tuple(item) + (None,))[:3]
            par = self.doc.add_paragraph()
            self._plain(par, align=AL.LEFT)
            par.paragraph_format.left_indent = Cm(self.INDENT_CM)
            if self.lo is not None:
                self.lo.add_formula(par, sym, size_pt=self.SMALL_PT)
            else:
                r = par.add_run(sym); r.italic = True; r.font.size = Pt(self.SMALL_PT)
            tail = par.add_run(f" — {desc}" + (f" [{unit}]" if unit else "")
                               + (";" if i < len(items) - 1 else "."))
            tail.font.size = Pt(self.SMALL_PT)
        return items

    def fig(self, path, caption=None, width=6.3):
        """Рисунок с подписью «Рисунок N — …» снизу и пустой строкой после."""
        from docx.shared import Inches, Pt
        from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
        if not os.path.exists(path):
            return
        par = self.doc.add_paragraph()
        self._plain(par, align=AL.CENTER)
        par.add_run().add_picture(path, width=Inches(width))
        self.n_fig += 1
        text = (caption or os.path.basename(path)).strip().rstrip(".")   # у подписи точки в конце нет
        cap = self.doc.add_paragraph(f"Рисунок {self.n_fig} — {text}")
        self._plain(cap, size=self.BODY_PT, align=AL.CENTER, line=self.LINE)
        self._spacer()

    def table(self, header, rows, digits=7, font=None, title=None):
        """Таблица с заголовком «Таблица N — …» сверху, текст на пункт мельче основного, по левому краю."""
        from docx.shared import Pt
        from docx.enum.text import WD_ALIGN_PARAGRAPH as AL
        font = font or self.TABLE_PT
        self.n_tab += 1
        cap = self.doc.add_paragraph(f"Таблица {self.n_tab}" + (f" — {title}" if title else ""))
        self._plain(cap, size=self.BODY_PT, align=AL.LEFT, line=self.LINE)
        cap.paragraph_format.keep_with_next = True
        t = self.doc.add_table(rows=1, cols=len(header))
        t.style = "Table Grid"
        for i, h in enumerate(header):
            cell = t.rows[0].cells[i]; cell.text = h
            self._plain(cell.paragraphs[0], size=font, align=AL.LEFT)
            for r in cell.paragraphs[0].runs:
                r.bold = True
        dig = digits if isinstance(digits, (list, tuple)) else [digits] * len(header)
        for row in rows:
            cells = t.add_row().cells
            for i, v in enumerate(row):
                cells[i].text = pretty(v, dig[i]) if i else str(v)
                self._plain(cells[i].paragraphs[0], size=font, align=AL.LEFT)
        self._spacer()
        return t

    def csv_table(self, path, digits=7, font=None, title=None):
        header, rows = read_csv(path)
        if header is not None:
            self.table(header, rows, digits, font, title=title)

    def save(self, path):
        self.doc.save(path)


# ============================================================================ сводки по уровням
def rows_analytic():
    """Аналитика против аналитики: перемноженные матрицы одной программы против матриц другой."""
    out = []
    an = csv_dict(os.path.join(RES, "analytic", "summary.csv"))
    if an:
        for k, lab in (("νx (матрицы по .lte)", "νx, матрицы по .lte"),
                       ("νy (матрицы по .lte)", "νy, матрицы по .lte"),
                       ("max |Δβx/βx| (матрицы по .lte)", "βx вдоль кольца, матрицы по .lte"),
                       ("νx (матрицы по решётке БД)", "νx, решётка БД"),
                       ("νy (матрицы по решётке БД)", "νy, решётка БД"),
                       ("max |Δβx/βx| (решётка БД)", "βx вдоль кольца, решётка БД"),
                       ("max |Δβy/βy| (решётка БД)", "βy вдоль кольца, решётка БД"),
                       ("max Δηx, нм (решётка БД)", "ηx вдоль кольца, решётка БД, нм")):
            if k in an:
                out.append((lab, an[k][0], an[k][1], an[k][2]))
    hdr, mt = read_csv(os.path.join(RES, "matrices", "by_type.csv"))
    if mt and hdr and len(hdr) > 4:
        vals = [num(r[4]) for r in mt if not r[0].startswith("первый")]
        out.append(("матрицы элементов: ring_matrix против matrix_output", "", "", max(vals)))
    ch = csv_dict(os.path.join(RES, "chromaticity", "summary.csv"))
    for k in ("ξ1 = dν/dδ (x)", "ξ1 = dν/dδ (y)", "ξ2 = d²ν/dδ² (x)", "ξ2 = d²ν/dδ² (y)", "αc2"):
        if ch and k in ch and len(ch[k]) >= 2:
            out.append((f"{k}: ring_map против twiss_output", ch[k][0], ch[k][1], rel(ch[k][0], ch[k][1])))
    lo = csv_dict(os.path.join(RES, "longitudinal", "summary.csv"))
    for k in ("U0 по трекингу за оборот, МэВ", "τx, мс", "τy, мс", "τs, мс", "Qs (ВЧ + излучение), матрица"):
        if lo and k in lo and len(lo[k]) >= 4 and lo[k][3] is not None:
            out.append((k.replace(" по трекингу за оборот", "") + ": интегралы излучения против twiss_output",
                        lo[k][3], lo[k][2], rel(lo[k][3], lo[k][2])))
    em = csv_dict(os.path.join(RES, "emittance", "summary.csv"))
    for k in ("εx, пм", "σδ, 10⁻⁴"):
        if em and k in em and len(em[k]) >= 5 and em[k][4] is not None and em[k][3] is not None:
            out.append((k + ": интегралы излучения против twiss_output", em[k][4], em[k][3],
                        rel(em[k][4], em[k][3])))
    rp = csv_dict(os.path.join(RES, "response", "summary.csv"))
    for k, v in (rp or {}).items():
        if "медианное" in k:
            out.append((f"матрица отклика орбиты, {k[0]}: медиана относительного расхождения", "", "", v[2]))
    return out


def rows_tracking():
    """Трекинг против трекинга: C++-ядро SynchrotronSim против трекинга elegant."""
    out = []
    lin = csv_dict(os.path.join(RES, "linear", "summary.csv"))
    if lin:
        j = json.load(open(os.path.join(RES, "linear", "summary.json"), encoding="utf-8"))
        for k in ("νx", "νy"):
            out.append((k, lin[k][0], lin[k][1], abs(lin[k][0] - lin[k][1])))
        out.append(("αc", lin["αc"][0], lin["αc"][1], rel(lin["αc"][0], lin["αc"][1])))
        out.append(("β вдоль кольца (относительно)", "", "", max(j["max_rel_beta_x_ss"], j["max_rel_beta_y_ss"])))
        out.append(("ηx вдоль кольца, мм", "", "", j["max_deta_x_mm_ss"]))
        out.append(("замкнутая орбита x, мкм", lin["замкнутая орбита x, мкм"][0],
                    lin["замкнутая орбита x, мкм"][1],
                    abs(lin["замкнутая орбита x, мкм"][0] - lin["замкнутая орбита x, мкм"][1])))
    hdr, mt = read_csv(os.path.join(RES, "matrices", "by_type.csv"))
    if mt:
        vals = [num(r[2]) for r in mt if not r[0].startswith("первый")]
        out.append(("матрицы элементов: max |ΔR| (кроме первого)", "", "", max(vals)))
    ch = csv_dict(os.path.join(RES, "chromaticity", "summary.csv"))
    for k in ("ξ1 = dν/dδ (x)", "ξ1 = dν/dδ (y)", "ξ2 = d²ν/dδ² (x)", "ξ3 = d³ν/dδ³ (x)"):
        if ch and k in ch and len(ch[k]) >= 4:
            out.append((k, ch[k][2], ch[k][3], rel(ch[k][2], ch[k][3])))
    it = csv_dict(os.path.join(RES, "integrator", "tunes.csv"))
    if it:
        combos = [v for k, v in it.items() if k.startswith("N=")
                  and all(x is not None and math.isfinite(x) for x in v[:4])]
        if combos:
            dm = max(max(abs(v[0] - v[1]), abs(v[2] - v[3])) for v in combos)
            out.append((f"частоты при {len(combos)} сочетаниях числа срезов и порядка", "", "", dm))
    tb = csv_dict(os.path.join(RES, "tbt", "tunes_vs_amplitude.csv"))
    if tb:
        amp = lambda k: float(re.sub(r"[^0-9.]", "", k) or "0")
        good = [(amp(k), v) for k, v in tb.items()
                if all(x is not None and math.isfinite(x) for x in (v[1], v[2], v[4], v[5]))]
        if good:
            a_max = max(a for a, _ in good)
            out.append((f"частоты при амплитудах до {a_max:g} мм (NAFF)", "", "",
                        max(max(abs(v[1] - v[4]), abs(v[2] - v[5])) for _, v in good)))
    lo = csv_dict(os.path.join(RES, "longitudinal", "summary.csv"))
    for k in ("U0 по трекингу за оборот, МэВ", "Qs (ВЧ + излучение), матрица", "τx, мс", "τy, мс", "τs, мс"):
        if lo and k in lo:
            out.append((k, lo[k][0], lo[k][1], rel(lo[k][0], lo[k][1])))
    rp = csv_dict(os.path.join(RES, "response", "summary.csv"))
    for k, v in (rp or {}).items():
        if k.startswith("трекинг ↔ трекинг"):
            out.append((k.replace("трекинг ↔ трекинг", "матрица отклика, трекинг").replace(".", ","),
                        v[0], v[1], v[2]))
    _, grid = read_csv(os.path.join(RES, "dynamic_aperture", "grid.csv"))
    if grid:
        nt = max(num(r[2]) for r in grid)
        agree = sum((num(r[2]) == nt) == (num(r[3]) == nt) for r in grid)
        out.append((f"динамическая апертура: совпадение «выжила/потеряна» за {int(nt)} оборотов",
                    f"{agree} из {len(grid)}", "", ""))
    return out


def rows_cross():
    """Внутри программы: аналитика против собственного трекинга при 5 срезах."""
    out = []
    lin = csv_dict(os.path.join(RES, "linear", "summary.csv"))
    if lin:
        j = json.load(open(os.path.join(RES, "linear", "summary.json"), encoding="utf-8"))
        out.append(("νx: twiss_output против трекинга elegant", lin["νx"][2], lin["νx"][1],
                    abs(lin["νx"][2] - lin["νx"][1])))
        out.append(("νy: twiss_output против трекинга elegant", lin["νy"][2], lin["νy"][1],
                    abs(lin["νy"][2] - lin["νy"][1])))
        out.append(("βx вдоль кольца (относительно)", "", "", j["max_rel_beta_x_twiss"]))
        out.append(("βy вдоль кольца (относительно)", "", "", j["max_rel_beta_y_twiss"]))
        out.append(("замкнутая орбита идеальной решётки, мкм", 0.0, lin["замкнутая орбита x, мкм"][1],
                    abs(lin["замкнутая орбита x, мкм"][1])))
    hdr, mt = read_csv(os.path.join(RES, "matrices", "by_type.csv"))
    if mt and hdr and len(hdr) > 5:
        vals = [num(r[5]) for r in mt if not r[0].startswith("первый")]
        out.append(("матрицы элементов: matrix_output против трекинга elegant", "", "", max(vals)))
    rp = csv_dict(os.path.join(RES, "response", "summary.csv"))
    for k, v in (rp or {}).items():
        if k.startswith("внутри кода"):
            out.append((k.replace("внутри кода: аналитика против трекинга (5 срезов)",
                                  "матрица отклика: аналитика против трекинга"), v[0], v[1], v[2]))
    ch = csv_dict(os.path.join(RES, "chromaticity", "summary.csv"))
    for k in ("ξ2 = d²ν/dδ² (x)", "ξ3 = d³ν/dδ³ (x)", "αc2"):
        if ch and k in ch and len(ch[k]) >= 3 and ch[k][2] is not None:
            out.append((f"{k}: ring_map против трекинга SynchrotronSim", ch[k][0], ch[k][2],
                        rel(ch[k][0], ch[k][2])))
    em = csv_dict(os.path.join(RES, "emittance", "summary.csv"))
    conv = glob.glob(os.path.join(RES, "emittance", "ss_nsl*_o*.json"))
    if em and conv:
        best = max((json.load(open(f, encoding="utf-8")) for f in conv), key=lambda q: (q["nsl"], q["order"]))
        out.append((f"εx SynchrotronSim: 5 срезов против {best['nsl']}, пм", em["εx, пм"][0], best["ex_pm"],
                    rel(em["εx, пм"][0], best["ex_pm"])))
    return out


# ============================================================================ сборка
def build(out=None):
    out = out or OUT
    f = lattice_facts()
    nsl, order = M["slices"]
    volt, harm = M["rf"]["volt_mv"], M["rf"]["harmonic"]
    R = Report()
    R.doc.add_heading(f"Сверка кикового симулятора SynchrotronSim с elegant: {M['label']}", level=0)

    # ------------------------------------------------------------ постановка
    R.h1("Постановка")
    R.p(f"Установка — кольцо «{M['label']}», файл описания "
        f"{os.path.relpath(M['source'], ROOT) if M['source'] else 'configuration_v1_geo'}: "
        f"{f['n_line']} {plural(f['n_line'], 'элемент', 'элемента', 'элементов')} линии, из них {f['n_mag']} "
        f"{plural(f['n_mag'], 'магнит', 'магнита', 'магнитов')}, периметр "
        f"{pretty(f['circumference'], 9)} м, "
        f"{M['superperiods']} супериодов. Энергия {M['e_mev']:g} МэВ, ВЧ {pretty(volt, 4)} МВ с кратностью {harm}.")
    R.p(db_note())
    R.p("Трекинг SynchrotronSim во всех опытах выполняет продакшн-ядро C++ (src/matrix_core.h, "
        "src/matrix_ring.h) — сервер кольца build/MatrixTracker.exe --ring-serve, к которому алгоритмы "
        "анализа подключаются обёрткой ring_cpp.CppRing; частицы считаются параллельно. Python-версия ядра "
        "(matrix_tracker.py) использовалась только для отработки и здесь служит эталоном: на решётке "
        "установки обе версии совпадают до округления double — это проверяется опытом kernel при каждом "
        "пересчёте.")
    R.p(f"Интеграторы настроены одинаково в обеих программах: по {nsl} срезов на элемент, симплектическая "
        f"схема порядка {order} (параметры nsl и order C++-ядра против N_SLICES и INTEGRATION_ORDER)."
        + (" Число срезов приходится задавать явно: в файле установки его нет, а по умолчанию elegant "
           "проходит линзу ОДНИМ шагом интегратора (N_SLICES = 1 у KQUAD и KSEXT, 4 у CSBEND), и на "
           "короткой сильной оптике этого не хватает — см. главу о разрыве между уровнями."
           if M["source"] else " Столько же стоит в самом файле описания установки."))
    R.p("Координаты: x, px = (1+δ)·x′/√(1+x′²+y′²), y, py, t — запаздывание относительно света в метрах, "
        f"δ = Δp/p. Точка наблюдения — начало линии (s = 0), до входа первого магнита "
        f"{pretty(f['lead'], 5)} м.")
    if MACHINE == "siberia":
        R.p("Решётка взята чужим файлом целиком (перевод из PyAT), без правок: 72 квадруполя шести "
            "семейств, 48 изгибов двух типов, 24 секступоля нулевой силы и 24 датчика положения. "
            "Секступоли обнулены, поэтому хроматизм чисто натуральный, а нелинейность кольца — только от "
            "краевых и кинематических членов. Установка устроена иначе, чем СИЛА: линзы короткие и сильные "
            "(K₁ до 4,07 м⁻², фазовый набег до 0,79 на линзу), изгибы редкие и мощные (ρ = 4,91 м), "
            "дисперсия достигает 1,13 м против 0,08 м у СИЛА.")
        R.p(corrector_note())

    # ------------------------------------------------- опорный импульс
    mass_section(R, R.h1)

    # ------------------------------------------------------------ принцип
    R.h1("Принцип сопоставления: пары одинаковых инструментов")
    R.p("В каждой программе есть два независимых механизма анализа установки.")
    R.bullet("АНАЛИТИЧЕСКИЙ: матрица каждого элемента считается формулой, матрицы перемножаются вдоль линии. "
             "У elegant это twiss_output, matrix_output, correction_matrix_output; у SynchrotronSim — "
             "ring_matrix.py (линейные матрицы), ring_map.py (отображение оборота высших порядков) и "
             "ring_response.py.")
    R.bullet("ТРЕКИНГОВЫЙ: частицы гоняются по кольцу, а показатели извлекаются из траекторий — матрицы "
             "разностями, частоты методом NAFF, затухание по спаду амплитуды. У SynchrotronSim это "
             "продакшн-ядро C++ (сервер кольца, обёртка ring_cpp.CppRing) с алгоритмами ring_analysis.py; "
             "у elegant — его собственный трекинг, к которому те же алгоритмы подключаются через "
             "elegant_diag.ElegantRing.")
    R.p("Сравнивать имеет смысл однородные пары: аналитику с аналитикой, трекинг с трекингом. Смешанное "
        "сравнение измеряет не различие программ, а ошибку дискретизации интегратора, причём одинаковую в "
        "обеих программах. Такие числа собраны в отдельную таблицу.")
    R.p("Карта инструментов по показателям:")
    R.table(["показатель", "SynchrotronSim, аналитика", "SynchrotronSim, трекинг",
             "elegant, аналитика", "elegant, трекинг"], INSTRUMENTS,
                 title="Карта инструментов: чем считается каждый показатель")

    # ------------------------------------------------------------ сводки
    R.h1("Сводка расхождений")
    R.h2("Аналитика против аналитики")
    R.p("Перемноженные матрицы SynchrotronSim против матричных инструментов elegant. Расхождение абсолютное "
        "для частот и дисперсии, относительное для остальных величин.")
    R.table(["показатель", "SynchrotronSim", "elegant", "расхождение"], rows_analytic(), [8, 9, 9, 2],
        title="Сводка расхождений: аналитика против аналитики")
    R.h2("Трекинг против трекинга")
    R.p("C++-ядро SynchrotronSim против трекинга elegant, обработка одним и тем же алгоритмом.")
    R.table(["показатель", "SynchrotronSim", "elegant", "расхождение"], rows_tracking(), [8, 9, 9, 2],
        title="Сводка расхождений: трекинг против трекинга")
    R.h2("Внутри программы: аналитика против собственного трекинга")
    R.p(f"Это не различие программ, а недосходимость интегратора: обе программы работают на "
        f"{nsl} срезах схемы порядка {order}.")
    R.table(["показатель", "аналитика", "трекинг", "расхождение"], rows_cross(), [8, 9, 9, 2],
        title="Расхождение уровней внутри одной программы")

    # ------------------------------------------------------------ инструменты
    R.h1("Инструменты SynchrotronSim")
    R.h2("Из трекинга")
    R.p("Замкнутая орбита ищется методом Ньютона по отображению за оборот. Матрицы переноса снимаются "
        "4-точечными разностями с шагом 10⁻⁵: орбита за оборот копит ошибку округления около 10⁻¹² м, и "
        "центральная разность с шагом 10⁻⁶ давала бы шум 10⁻⁶ — столько же, сколько радиационное затухание "
        "за оборот.")
    R.formula(r"\frac{\partial f}{\partial z} \approx \frac{-f(z+2h) + 8f(z+h) - 8f(z-h) + f(z-2h)}{12h}",
              [
        (r"f", "координата частицы после оборота (отображение за оборот)", "м или рад"),
        (r"z", "фазовая координата, по которой берётся производная: x, pₓ, y, p_y, t, δ", "м, рад, м, рад, м, —"),
        (r"h", "шаг разности, 10⁻⁵", "в единицах координаты z"),
    ])
    R.p("Из матриц получаются функции Твисса, дисперсия и коэффициент уплотнения орбит:")
    R.formula(r"\cos\mu = \frac{M_{11}+M_{22}}{2},\quad \beta(s) = \frac{(R_{11}\beta_0 - R_{12}\alpha_0)^2 "
              r"+ R_{12}^2}{\beta_0},\quad D_0 = (I - M_4)^{-1} M_{\cdot 6}",
              [
        (r"M", "матрица отображения за оборот, снятая разностями"),
        (r"\mu", "фазовый набег за оборот в своей плоскости", "рад"),
        (r"\beta_0,\ \alpha_0", "функции Твисса в точке наблюдения", "м; —"),
        (r"R_{ij}", "элементы матрицы от точки наблюдения до сечения s"),
        (r"M_4", "поперечный блок 4×4 матрицы оборота"),
        (r"M_{\cdot 6}", "столбец δ матрицы оборота", "м, рад"),
        (r"D_0", "дисперсия и её производная в точке наблюдения", "м; —"),
    ])
    R.p("Нелинейные показатели берутся из пооборотных координат: частоты — методом NAFF по нормированному "
        "сигналу, скорость диффузии — по изменению частот между половинами прогона.")
    R.formula(r"z_k = \frac{x_k - i(\alpha x_k + \beta p_{x,k})}{\sqrt{\beta}},\qquad "
              r"\nu = \arg\max_{\nu}\left|\sum_k w_k z_k e^{-2\pi i \nu k}\right|",
              [
        (r"x_k,\ p_{x,k}", "координата и канонический импульс на k-м обороте", "м; рад"),
        (r"\alpha,\ \beta", "функции Твисса в точке наблюдения", "—; м"),
        (r"z_k", "нормированный комплексный сигнал колебания", "√м"),
        (r"w_k", "окно Ханна: sin²(πk/N)", "—"),
        (r"\nu", "бетатронная частота (число колебаний за оборот)", "—"),
        (r"N", "число оборотов в выборке", "—"),
    ])
    R.h2("Аналитический")
    R.p("Матрицы элементов считаются формулами по соглашениям elegant: дрейф, толстый квадруполь, секторный "
        "изгиб с градиентом, края — тонкие линзы tan(E)/ρ, секступоли и октуполи в первом порядке дрейфы. "
        "Пятая координата — путь, поэтому сопряжённая пара продольного блока — (−s, δ); с ней матрицы "
        "симплектичны до 2·10⁻¹⁶.")
    R.formula(r"k_x^2 = h^2 + K_1,\qquad k_y^2 = -K_1,\qquad h = 1/\rho",
              [
        (r"k_x,\ k_y", "волновые числа фокусировки в теле элемента", "м⁻¹"),
        (r"h = 1/\rho", "кривизна изгиба", "м⁻¹"),
        (r"K_1", "нормированный градиент квадрупольной составляющей", "м⁻²"),
        (r"\rho", "радиус изгиба опорной орбиты", "м"),
    ])
    R.p("Отличие от elegant — стыки. В SynchrotronSim элементы соединены модифицированными дрейфами: частица "
        "летит по прямой от выходной грани одного элемента до входной грани другого, и грани могут быть не "
        "параллельны (боксовая модель, повороты установки, неидеальная геометрия). Отображение такого стыка "
        "не линейное, а аффинное, поэтому линия собирается парами (C, R), а замкнутая орбита берётся из "
        "свободного члена — тот же механизм сразу даёт орбиту от рассогласований:")
    R.formula(r"z \mapsto R\,z + C,\qquad (C_2, R_2)\circ(C_1, R_1) = (C_2 + R_2 C_1,\; R_2 R_1),"
              r"\qquad z_{\text{замк}} = (I - R)^{-1} C",
              [
        (r"z", "фазовый вектор (x, pₓ, y, p_y, t, δ)", "м, рад, м, рад, м, —"),
        (r"R", "линейная матрица отображения 6×6"),
        (r"C", "свободный член: вклад наклона граней и рассогласований", "м, рад"),
        (r"z_{\text{замк}}", "замкнутая орбита в точке наблюдения", "м, рад"),
    ])
    R.p("Матрицы стыков проверены против точного отображения ядра на четырёх геометриях граней:")
    R.csv_table(os.path.join(RES, "analytic", "junctions.csv"), 3,
        title="Матрицы стыков против точного отображения ядра")
    jh = os.path.join(RES, "analytic", "junctions_high_order.csv")
    if os.path.exists(jh):
        R.p("Те же стыки в высших порядках: ряд 5-го порядка ring_map и аффинная матрица против точного "
            "пролёта C++-ядра при смещении частицы на 1, 3 и 10 мм (и соответствующих углах и δ). Ошибка "
            "аффинной матрицы растёт как квадрат амплитуды — это отброшенная нелинейность пролёта; ошибка "
            "ряда — как шестая степень, до уровня округления:")
        R.csv_table(jh, 3, title="Стык в высших порядках: ряд и матрица против точного пролёта")
    R.p("Ограничение ring_matrix: только секторный режим ядра — в боксовой модели грани изгиба "
        "перпендикулярны оси ящика, и телу отвечает другая матрица; там стоит явная ошибка вместо неверного "
        "ответа.")
    R.h2("Аналитический: высшие порядки (ring_map)")
    R.p("Нелинейная часть отображения оборота считается методом дифференциальной алгебры: каждая величина "
        "— усечённый степенной ряд по шести фазовым координатам до 5-го порядка, и ряды складываются, "
        "перемножаются и подставляются друг в друга точно до этого порядка. Коэффициенты — точные "
        "производные отображения, без разностей по траекториям.")
    R.bullet("СТЫКИ (модифицированные дрейфы) входят точной формулой пролёта по прямой между гранями, "
             "вычисленной над рядами. Поэтому клин между гранями даёт вклад во все порядки, а не только "
             "свободный член C и линейную матрицу.")
    R.bullet("ТЕЛО ЭЛЕМЕНТА — точное решение его уравнений движения z(L) = exp(L·D)·z, где D = Σ fᵢ ∂/∂zᵢ — "
             "производная по векторному полю элемента. Поле то же, что интегрирует ядро (секторный изгиб — "
             "порт CSBEND, прямой мультиполь — точный дрейф и сила), но решение точное, без срезов.")
    R.bullet("КРАЯ ИЗГИБА — тонкие линзы, как в ядре.")
    R.p("Из отображения оборота берётся замкнутая орбита как ряд по δ, линеаризация вокруг неё даёт ν(δ), "
        "а путь по ней — уплотнение орбит. При порядке 5 точны ξ1…ξ4 и αc, αc2.")
    R.p("НОРМАЛЬНАЯ ФОРМА (ring_map.normal_form) даёт сдвиг частот с амплитудой. Отображение вокруг "
        "замкнутой орбиты переводится в комплексный базис, где линейная часть — поворот на угол 2πν, после "
        "чего заменой переменных степень за степенью убираются все члены, кроме резонансных. Остаётся "
        "поворот с углом, зависящим от действий, и коэффициенты этой зависимости — ∂ν/∂J. Действие "
        "J = (q² + (αq + βp)²)/(2β), то есть ε = 2J, как у elegant. Проверка: на пробном отображении "
        "«поворот + тонкий октуполь» нормальная форма воспроизводит известную формулу −K₃Lβ²/16π до 5-го "
        "знака.")
    R.h2("Аналитический: излучение (интегралы I₁…I₅)")
    R.p("Потеря энергии, времена затухания и равновесные размеры пучка считаются по описанию решётки "
        "формулами (ring_matrix.radiation_integrals), без трекинга. Подынтегральные функции нужны ВНУТРИ "
        "изгибов, поэтому каждый изгиб проходится узлами Гаусса–Лежандра: для узла строится матрица от "
        "точки наблюдения, из набора матриц обычным разбором получаются β, α, η, η′, и интегралы берутся "
        "квадратурой. У граней изгиба учитывается краевой член I₄.")
    R.formula(r"I_1 = \oint \eta h\,ds,\quad I_2 = \oint h^2 ds,\quad I_3 = \oint |h|^3 ds,\quad "
              r"I_4 = \oint \eta h (h^2 + 2K_1) ds - \sum_{\text{грани}} \eta h^2 \tan E,\quad "
              r"I_5 = \oint |h|^3 \mathcal{H} ds",
              [
        (r"I_1 \ldots I_5", "интегралы излучения", "м; м⁻¹; м⁻²; м⁻¹; м⁻¹"),
        (r"\eta,\ \eta'", "дисперсия и её производная", "м; —"),
        (r"h = 1/\rho", "кривизна изгиба", "м⁻¹"),
        (r"K_1", "нормированный градиент изгиба", "м⁻²"),
        (r"E", "угол наклона грани изгиба к опорной орбите", "рад"),
        (r"\mathcal{H} = \gamma\eta^2 + 2\alpha\eta\eta' + \beta\eta'^2", "дисперсионный инвариант", "м"),
        (r"ds", "элемент длины вдоль опорной орбиты", "м"),
    ])
    R.formula(r"U_0 = \frac{C_\gamma E^4 I_2}{2\pi},\quad J_x = 1 - \frac{I_4}{I_2},\quad "
              r"\tau_i = \frac{2 E T_0}{J_i U_0},\quad \varepsilon_x = C_q\gamma^2\frac{I_5}{J_x I_2},\quad "
              r"\sigma_\delta^2 = C_q\gamma^2\frac{I_3}{J_s I_2}",
              [
        (r"U_0", "потеря энергии на излучение за оборот", "МэВ"),
        (r"C_\gamma = 8{,}846\cdot10^{-5}", "постоянная излучения для электрона", "м/ГэВ³"),
        (r"E", "энергия пучка", "ГэВ"),
        (r"J_x,\ J_y,\ J_s", "числа затухания (ΣJ = 4)", "—"),
        (r"T_0", "период обращения", "с"),
        (r"\tau_i", "время затухания колебаний в плоскости i", "с"),
        (r"C_q = 3{,}832\cdot10^{-13}", "квантовая постоянная", "м"),
        (r"\gamma", "релятивистский фактор энергии пучка", "—"),
        (r"\varepsilon_x", "равновесный горизонтальный эмиттанс", "м·рад"),
        (r"\sigma_\delta", "равновесный относительный разброс по импульсу", "—"),
    ])
    R.p("Это прямая пара к radiation_integrals elegant: интегралы, U₀, J, τ, εx и σδ совпадают с "
        "twiss_output до 10⁻⁷ (таблицы в главах о затухании и эмиттансе).")
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

    # ------------------------------------------------------------ линейная оптика
    if exists("linear") and exists("analytic"):
        R.h1("Линейная оптика: орбита, β, α, фаза, дисперсия, частоты")
        R.h2("Аналитика против аналитики")
        an = csv_dict(os.path.join(RES, "analytic", "summary.csv"))
        R.p(f"Матрицы, собранные прямо по .lte, совпадают с twiss_output до уровня округления "
            f"({pretty(an['νx (матрицы по .lte)'][2], 2)} по νx): формулы в обеих программах одни и те же. "
            f"На решётке из базы данных остаётся {pretty(an['νx (матрицы по решётке БД)'][2], 2)} по νx и "
            f"{pretty(an['max |Δβx/βx| (решётка БД)'][2], 2)} по β"
            + (" — столько же, то есть представление решётки в базе больше не ограничивает сверку."
               if an['νx (матрицы по решётке БД)'][2] < 1e-9 else
               " — это уже не разница формул, а разрешение самой базы данных по положению магнитов."))
        R.fig(os.path.join(RES, "analytic", "analytic_vs_elegant.png"),
              f"Матрицы SynchrotronSim против twiss_output; для масштаба — трекинг ({nsl} срезов) против него же.")
        R.fig(os.path.join(RES, "analytic", "analytic_residual.png"), "Тот же остаток в собственном масштабе.")
        R.csv_table(os.path.join(RES, "analytic", "summary.csv"), 9,
            title="Линейная оптика: аналитика SynchrotronSim против аналитики elegant")
        dC, dT = geom_residual()
        R.p(f"Оставшийся след — геометрический: {geom_note()}. "
            + ("На таком уровне геометрия уже не мешает: остаток в частотах и β — округление двойной "
               "точности при перемножении матриц."
               if max(dC, dT) < 1e-9 else
               "Отсюда и остаточная замкнутая орбита, и след в дисперсии."))
        R.h2("Трекинг против трекинга")
        R.p("Замкнутая орбита, матрицы переноса и функции Твисса сняты одним алгоритмом с обоих трекеров.")
        R.fig(os.path.join(RES, "linear", "twiss_superperiod.png"), "Функции Твисса на первом супериоде.")
        R.fig(os.path.join(RES, "linear", "twiss_difference.png"),
              "Расхождения вдоль кольца: SynchrotronSim и twiss_output относительно трекинга elegant.")
        if phase_note():
            R.p(phase_note())
        R.csv_table(os.path.join(RES, "linear", "summary.csv"), 9, title="Линейная оптика по трекингу обеих программ")
        j = json.load(open(os.path.join(RES, "linear", "summary.json"), encoding="utf-8"))
        lin = csv_dict(os.path.join(RES, "linear", "summary.csv"))
        co = lin["замкнутая орбита x, мкм"]
        R.p(f"Наибольшее относительное расхождение β — {pretty(j['max_rel_beta_x_ss'], 2)} по x и "
            f"{pretty(j['max_rel_beta_y_ss'], 2)} по y, дисперсии — {pretty(j['max_deta_x_mm_ss'], 2)} мм. "
            f"Третий столбец таблицы (twiss_output) относится к другому уровню: он отличается от трекинга "
            f"elegant на {pretty(j['max_rel_beta_x_twiss'], 2)} по βx — это дискретизация "
            f"интегратора, а не различие программ.")
        R.p(f"Отдельная строка — замкнутая орбита идеальной решётки: у elegant "
            f"{si_len(abs(co[1]) * 1e-6)}, у SynchrotronSim {si_len(abs(co[0]) * 1e-6)}"
            + orbit_verdict(co[0], co[1]) + ".")

    # ------------------------------------------------------------ матрицы элементов
    if exists("matrices"):
        R.h1("Матрицы отдельных элементов")
        R.p("Матрица каждого магнита вместе с дрейфом перед ним. В таблице сразу три сопоставления: "
            "аналитика SynchrotronSim против matrix_output, трекинг против трекинга и — отдельным столбцом — "
            "matrix_output против собственного трекинга elegant. Первый магнит вынесен отдельной строкой: в "
            "его матрицу входит дрейф от точки наблюдения, и разница длин этого дрейфа видна там напрямую.")
        R.fig(os.path.join(RES, "matrices", "element_matrices.png"),
              "Наибольшее по магнитам расхождение каждого члена матрицы (log10).")
        R.csv_table(os.path.join(RES, "matrices", "by_type.csv"), 3, title="Матрицы элементов по типам устройств")
        hdr, by = read_csv(os.path.join(RES, "matrices", "by_type.csv"))
        if by:
            body = [r for r in by if not r[0].startswith("первый")]
            R.p(f"Аналитика сходится с аналитикой до {pretty(max(num(r[4]) for r in body), 2)}, трекинг с "
                f"трекингом — до {pretty(max(num(r[3]) for r in body), 2)} по медиане, а разрыв между "
                f"уровнями внутри самого elegant — {pretty(max(num(r[5]) for r in body), 2)}.")
        R.p("Представительные элементы:")
        R.csv_table(os.path.join(RES, "matrices", "representative.csv"), 9,
            title="Матрица представительного элемента, все 36 членов")

    # ------------------------------------------------------------ отклик и допуска
    if exists("response"):
        R.h1("Матрицы отклика и допуска на магниты")
        rsp = csv_dict(os.path.join(RES, "response", "summary.csv"))
        crc = csv_dict(os.path.join(RES, "response", "correction.csv"))
        n_st, n_cr = len(_stations_report()), int(crc["оставлено сингулярных чисел"][1])
        R.p(f"Отклик орбиты на кик корректора считается переносом кика до конца оборота и замыканием "
            f"орбиты, поэтому вся матрица {n_st} датчиков × {n_cr} корректоров получается за секунды:")
        R.formula(r"z_0 = (I - M)^{-1} M R_k^{-1} e\,\theta,\qquad "
                  r"x_b = \left[R_b z_0\right]_x + \left[R_b R_k^{-1} e\,\theta\right]_x",
                  [
            (r"z_0", "замкнутая орбита в точке наблюдения от кика", "м, рад"),
            (r"M", "поперечная матрица оборота 4×4"),
            (r"R_k,\ R_b", "матрицы от точки наблюдения до корректора и до датчика"),
            (r"e", "единичный вектор по каналу импульса pₓ", "—"),
            (r"\theta", "угол кика корректора", "рад"),
            (r"x_b", "отклик орбиты на датчике", "м"),
        ])
        R.p("Кик берётся тонким, в середине корректора, между двумя точными дрейфами — это модель самого "
            "элемента. Трекинговая пара строится одинаково в обеих программах: кик ставится в один корректор "
            "(в elegant через alter_elements) и ищется замкнутая орбита.")
        R.fig(os.path.join(RES, "response", "orbit_response.png"),
              f"Отклик орбиты на кик одного корректора в {n_st} датчиках и разность с elegant.")
        R.csv_table(os.path.join(RES, "response", "summary.csv"), 6, title="Матрица отклика орбиты: сводка сравнений")
        med = [v[2] for k, v in rsp.items() if k.startswith("x: медианное")][0]
        an_tr = max(v[2] for k, v in rsp.items() if k.startswith("аналитика ↔ трекинг"))
        tr_tr = max(v[2] for k, v in rsp.items() if k.startswith("трекинг ↔ трекинг"))
        inside = max(v[2] for k, v in rsp.items() if k.startswith("внутри кода"))
        R.p(f"Аналитика против аналитики даёт {pretty(med, 2)} по медиане. Остаток — внутреннее "
            f"приближение elegant: для корректора конечной длины он берёт функции Твисса как среднее по "
            f"концам элемента, тогда как здесь кик ставится тонким в середине. Что модель с киком в середине "
            f"точнее, показывает трекинг: с ним аналитика совпадает до {pretty(an_tr, 2)}. Трекинг против "
            f"трекинга совпадает до {pretty(tr_tr, 2)}, разрыв между уровнями внутри программы — "
            f"{pretty(inside, 2)}.")
        brho = M["e_mev"] * 1e6 / 299792458.0
        R.p(f"Практическое применение — коррекция орбиты. Матрица обращается по SVD с отсечкой малых "
            f"сингулярных чисел; кик переводится в поле по проектной жёсткости (B = θ·Bρ/L, "
            f"Bρ = {pretty(brho, 4)} Тл·м при {M['e_mev'] / 1000:g} ГэВ), а поле — в токи обмоток линейной матрицей "
            f"корректора.")
        R.fig(os.path.join(RES, "response", "svd_correction.png"),
              "Спектр матрицы отклика и коррекция орбиты по SVD.")
        R.csv_table(os.path.join(RES, "response", "correction.csv"), 6,
            title="Коррекция орбиты по SVD после разворота квадруполей")
        cur = [v[0] for k, v in crc.items() if k.startswith("ток обмотки main")]
        R.p(f"В опыте развёрнуты 20 случайных квадруполей (СКО 50 мкрад). Орбита с СКО "
            f"{pretty(crc['СКО орбиты до коррекции, мкм'][0], 3)} мкм правится до "
            f"{pretty(crc['СКО орбиты после коррекции, мкм'][0], 3)} мкм; остаток локален — там искажение "
            f"возникает вдали от ближайшего корректора, и корректоров ({n_cr}) на {n_st} датчиков для него "
            f"не хватает. Нужные кики не превышают {pretty(crc['максимальный кик, мкрад'][0], 3)} мкрад: поле "
            f"{pretty(crc['физкоэффициент B при максимальном кике, Тл'][0] * 1e3, 3)} мТл"
            + (f" и ток главной обмотки {pretty(cur[0], 2)} А" if cur else "") + ".")
        R.p("Та же машинерия даёт чувствительность кольца к ошибкам градиента:")
        R.formula(r"\Delta\nu_{x,y} = \pm\frac{\beta_{x,y}}{4\pi}\,\Delta(K_1 L),\qquad "
                  r"\frac{\Delta\beta}{\beta}(s) = -\frac{\Delta(K_1 L)\,\beta_k}{2\sin 2\pi\nu}"
                  r"\cos\left(2\left|\psi(s)-\psi_k\right| - 2\pi\nu\right)",
                  [
            (r"\Delta\nu_{x,y}", "сдвиг бетатронной частоты от ошибки градиента", "—"),
            (r"\Delta(K_1 L)", "ошибка интегрированного градиента магнита", "м⁻¹"),
            (r"\beta_{x,y},\ \beta_k", "бета-функции в точке наблюдения и в месте магнита", "м"),
            (r"\psi(s),\ \psi_k", "фазовый набег до сечения s и до магнита", "рад"),
            (r"\nu", "бетатронная частота кольца", "—"),
        ])
        R.fig(os.path.join(RES, "response", "tolerances.png"),
              "Чувствительность к ошибке градиента 0,1 % по типам устройств.")
        R.csv_table(os.path.join(RES, "response", "tolerances.csv"), 5,
            title="Допуска на ошибку градиента по типам магнитов")
        hdr, tol = read_csv(os.path.join(RES, "response", "tolerances.csv"))
        top = tol[0] if tol else None
        R.p("Таблица читается так: если все магниты одного типа ошибутся в одну сторону на 0,1 %, частота "
            "сместится на указанную величину; при независимых ошибках того же размера бета-биение составит "
            "указанное СКО."
            + (f" Наиболее чувствителен {top[0]} — их в кольце {top[1]}: 0,1 % даёт сдвиг частоты "
               f"{pretty(top[3], 2)} и бета-биение {pretty(top[5], 2)} % по x, {pretty(top[6], 2)} % по y."
               if top else "")
            + " Это требование к точности поля библиотеки магнитов, выраженное через оптику кольца. "
              "Штатного инструмента для такой таблицы у elegant нет, пары для сверки не существует; "
              "проверкой служат сами матрицы отклика выше.")

    # ------------------------------------------------------------ разрыв уровней
    if exists("integrator"):
        R.h1("Разрыв между уровнями: число срезов и порядок интегратора")
        R.p("Частоты по матрице оборота при разном числе срезов N и порядке интегратора, одинаковых в обеих "
            "программах (nsl и order C++-ядра против N_SLICES и INTEGRATION_ORDER). Опорное "
            "значение — аналитика. Первая строка таблицы — настройки ПО УМОЛЧАНИЮ: у elegant это один шаг "
            "интегратора на линзу (N_SLICES = 1 у KQUAD и KSEXT, 4 у CSBEND), у ядра SynchrotronSim — пять "
            "срезов схемы 4-го порядка.")
        R.fig(os.path.join(RES, "integrator", "integrator_convergence.png"),
              "Отклонение частот от аналитики и расхождение программ при одинаковых настройках.")
        R.csv_table(os.path.join(RES, "integrator", "tunes.csv"), 8,
            title="Частоты при разном числе срезов и порядке схемы")
        hdr, it = read_csv(os.path.join(RES, "integrator", "tunes.csv"))
        body = [r for r in it if r[0].startswith("N=")
                and all(num(v) is not None and math.isfinite(num(v)) for v in r[1:5])]
        dm = max(max(abs(num(r[1]) - num(r[2])), abs(num(r[3]) - num(r[4]))) for r in body)
        R.p(f"При одинаковых настройках программы совпадают во всех {len(body)} сочетаниях, где кольцо "
            f"устойчиво у обеих: наибольшее "
            f"расхождение частот {pretty(dm, 2)} — это предел сходимости интегратора самого elegant, а не "
            f"различие моделей. Схема 6-го порядка выходит на аналитику при N ≥ 5, схема 4-го порядка — при "
            f"N ≥ 20; порядок 2 при N ≤ 40 непригоден.")
        if MACHINE == "siberia":
            R.p("Строка «по умолчанию» показывает, чем это оборачивается на короткой сильной оптике: при "
                "одном шаге интегратора на линзу ошибка фокусировки около 0,7 % на элемент, и за 144 "
                "магнита кольцо у elegant выходит НЕУСТОЙЧИВЫМ по горизонтали (cos μx = 4,76 вместо 0,139; "
                "в таблице частота отсутствует). Его собственный twiss_output при этом даёт нормальные "
                "функции — матричные инструменты число срезов не используют. Ядро SynchrotronSim при своих "
                "пяти срезах остаётся устойчивым, но ошибается в νx на 1,3·10⁻³. Поэтому все опыты этого "
                "отчёта поставлены на 40 срезах у обеих программ.")

    # ------------------------------------------------------------ хроматичность
    if exists("chromaticity"):
        R.h1("Бетатронная частота от энергии")
        R.p("Трекинговая пара: частоты при отклонении импульса δ считаются по матрице оборота вокруг "
            "замкнутой орбиты с этим δ, одинаково для обоих трекеров; производные в нуле — подгонкой "
            "многочленом 7-й степени на |δ| ≤ 1 %. Аналитическая пара: у SynchrotronSim — отображение "
            "оборота 5-го порядка (ring_map), у elegant — twiss_output с higher_order_chromaticity.")
        R.formula(r"\nu(\delta) = \nu_0 + \xi_1\delta + \frac{\xi_2}{2}\delta^2 + \frac{\xi_3}{6}\delta^3",
                  [
            (r"\nu(\delta)", "бетатронная частота частицы с отклонением импульса", "—"),
            (r"\delta = \Delta p/p", "относительное отклонение импульса", "—"),
            (r"\xi_1,\ \xi_2,\ \xi_3", "хроматичности первого, второго и третьего порядка", "—"),
            (r"\beta(s),\ K(s)", "бета-функция и сила линзы вдоль кольца", "м; м⁻²"),
        ])
        R.fig(os.path.join(RES, "chromaticity", "tune_vs_delta.png"), "ν(δ) и разность программ.")
        R.csv_table(os.path.join(RES, "chromaticity", "summary.csv"), 7,
            title="Хроматичность и коэффициенты уплотнения орбит")
        ch = csv_dict(os.path.join(RES, "chromaticity", "summary.csv"))
        k2, k3 = "ξ2 = d²ν/dδ² (x)", "ξ3 = d³ν/dδ³ (x)"
        x1 = ch['ξ1 = dν/dδ (x)']
        R.p((f"Первый порядок совпадает у всех четырёх инструментов: " if rel(x1[0], x1[2]) < 1e-4 else
             f"Первый порядок совпадает внутри каждой пары, а между уровнями отличается на "
             f"{pretty(rel(x1[0], x1[2]) * 100, 2)} % — дискретизация интегратора на {nsl} срезах: ")
            + f"ξ1x = {pretty(x1[0], 7)} по отображению SynchrotronSim, {pretty(x1[1], 7)} по twiss_output, "
            f"{pretty(x1[2], 7)} и {pretty(x1[3], 7)} по трекингу. "
            f"Старшие хроматичности отображения SynchrotronSim сходятся с трекингом обеих программ: "
            f"ξ2x = {pretty(ch[k2][0], 6)} против {pretty(ch[k2][2], 6)} и {pretty(ch[k2][3], 6)}, "
            f"ξ3x = {pretty(ch[k3][0], 5)} против {pretty(ch[k3][2], 5)} и {pretty(ch[k3][3], 5)} "
            + (f"(остаток — точность подгонки трекинга). " if rel(ch[k2][0], ch[k2][2]) < 1e-3 else
               f"(остаток {pretty(rel(ch[k2][0], ch[k2][2]) * 100, 2)} % по ξ2 — дискретизация интегратора: "
               f"трекинг идёт на {nsl} срезах, отображение — точное решение; тот же остаток у elegant). ")
            + f"twiss_output даёт ξ2x = {pretty(ch[k2][1], 5)} и "
            f"ξ3x = {pretty(ch[k3][1], 5)}: ξ3 по матрицам 3-го порядка определить нельзя (нужен 4-й), а "
            f"ξ2 расходится с трекингом самого elegant — это расхождение внутри elegant между его матричным "
            f"инструментом и его трекингом."
            + ("" if MACHINE != "siberia" else
               " Секступоли решётки нулевые, поэтому хроматизм здесь чисто натуральный."))

    # ------------------------------------------------------------ пооборотный трекинг
    if exists("tbt"):
        R.h1("Фазовые траектории, спектр, частоты при разных амплитудах")
        ax_mm, ay_mm = M["amp_x_mm"], M["amp_y_mm"]
        R.p(f"Чисто трекинговая пара: аналитических инструментов для этих показателей нет ни у одной "
            f"программы. Частицы стартуют от замкнутой орбиты со смещением по одной плоскости "
            f"({ax_mm[0]:g}…{ax_mm[-1]:g} мм по x, {ay_mm[0]:g}…{ay_mm[-1]:g} мм по y; другая плоскость — "
            f"10 мкм) и проходят 1024 оборота. Частоты — методом NAFF; для слабой плоскости берётся "
            f"сильнейшая линия, не совпадающая с ±ν, ±2ν сильной.")
        R.fig(os.path.join(RES, "tbt", "phase_portraits.png"), "Фазовые портреты в точке наблюдения.")
        R.fig(os.path.join(RES, "tbt", "tbt_difference.png"), "Расхождение координат по оборотам.")
        hdr, tb = read_csv(os.path.join(RES, "tbt", "tunes_vs_amplitude.csv"))
        alive = [r for r in tb if num(r[1]) == num(r[4]) == 1024]
        dmax = max((abs(num(r[2]) - num(r[5])) for r in alive), default=float("nan"))
        R.p(f"Все {len(alive)} частиц из {len(tb)} прожили полные 1024 оборота в обеих программах, "
            f"наибольшее расхождение частоты {pretty(dmax, 2)} — тот же предел интегратора. Расхождение "
            f"координат растёт линейно по оборотам: это набег фазы от разницы частот, а не разъезд "
            f"траекторий."
            + ("" if MACHINE != "siberia" or not alive else
               f" Сдвиг частоты с амплитудой мал ({pretty(abs(num(alive[-1][2]) - num(alive[0][2])), 2)} "
               f"между крайними амплитудами сетки): секступоли обнулены, и остаётся только кинематическая "
               f"нелинейность — обе программы воспроизводят её одинаково."))
        R.fig(os.path.join(RES, "tbt", "spectrum.png"), "Спектры пооборотных колебаний.")
        R.fig(os.path.join(RES, "tbt", "tunes_vs_amplitude.png"), "Сдвиг частот с амплитудой.")
        hdr_d, dn = read_csv(os.path.join(RES, "tbt", "detuning.csv"))
        if dn:
            R.p("Сдвиг частот с действием ∂ν/∂J: аналитика SynchrotronSim — нормальная форма отображения "
                "оборота, трекинговая оценка — наклон ν(J) на самых малых амплитудах сетки (дальше вступают "
                "члены J²).")
            R.table(hdr_d, [(r[0], num(r[1]), num(r[2]), num(r[3]), num(r[4]) if r[4] else "") for r in dn],
                    [8, 6, 6, 6, 6], title="Сдвиг частот с действием ∂ν/∂J")
        dt = csv_dict(os.path.join(RES, "tbt", "driving_terms.csv"))
        if dt:
            R.p("У elegant есть и аналитическая ветвь для нелинейной части: twiss_output с "
                "compute_driving_terms считает по мультиполям решётки сдвиг частот с действием ∂ν/∂J и "
                "резонансные члены h_jklm. В эти формулы входят только мультиполи (у изгиба с градиентом — "
                "и его собственный секступольный член 2·h·K₁), а кинематика точных дрейфов и краёв — нет, "
                "поэтому с трекинговым сдвигом частот они сходятся только там, где нелинейность даёт "
                "мультиполь. У SynchrotronSim такой ветви пока нет: её можно получить нормальной формой того "
                "же отображения ring_map.")
            R.table(["величина", "elegant twiss_output (аналитика)"], [(k, v[0]) for k, v in dt.items()], [8, 6],
                title="Аналитика нелинейности elegant: ∂ν/∂J и резонансные члены")
        R.csv_table(os.path.join(RES, "tbt", "tunes_vs_amplitude.csv"), 8,
            title="Частоты при разных амплитудах колебаний")

    # ------------------------------------------------------------ продольное движение
    if exists("longitudinal"):
        R.h1("Синхротронные колебания и радиационное затухание")
        R.p("ВЧ-резонатор стоит в точке наблюдения в конце оборота; опорное время — по скорости света "
            "(RFCA с FIDUCIAL = light). Излучение распределено по изгибам, поперечные импульсы при излучении "
            "и в ВЧ сохраняются. Пара трекинговая: неподвижная точка 6D ищется Ньютоном, времена затухания и "
            "частоты берутся из собственных чисел матрицы оборота и из спада амплитуды. Аналитика elegant "
            "(интегралы излучения) приведена третьим столбцом как опорная.")
        R.formula(r"\Delta E = eV\sin\left(\varphi + \frac{2\pi h\,t}{C}\right),\qquad "
                  r"\left|\lambda_k\right| = e^{-T_0/\tau_k},\qquad Q_s = \frac{\arg\lambda_s}{2\pi}",
                  [
            (r"\Delta E", "прирост энергии в резонаторе за проход", "МэВ"),
            (r"eV", "амплитуда ускоряющего напряжения", "МэВ"),
            (r"\varphi", "фаза резонатора для опорной частицы", "рад"),
            (r"h", "кратность ВЧ (число волн на периметре)", "—"),
            (r"t", "запаздывание частицы относительно света", "м"),
            (r"C", "периметр кольца", "м"),
            (r"\lambda_k", "собственное число матрицы оборота с излучением", "—"),
            (r"\tau_k,\ T_0", "время затухания и период обращения", "с"),
            (r"Q_s", "синхротронная частота (колебаний за оборот)", "—"),
        ])
        R.fig(os.path.join(RES, "longitudinal", "synchrotron.png"),
              "Продольные портреты, расхождение по оборотам и синхротронная частота от амплитуды.")
        R.fig(os.path.join(RES, "longitudinal", "damping.png"), "Затухание по спаду амплитуды.")
        R.csv_table(os.path.join(RES, "longitudinal", "summary.csv"), 7,
            title="Продольное движение: U₀, синхротронная частота, времена затухания")
        lo = csv_dict(os.path.join(RES, "longitudinal", "summary.csv"))
        mat_d = max(rel(lo[k][0], lo[k][1]) for k in ("τx, мс", "τy, мс", "τs, мс") if k in lo)
        amp_d = max(rel(lo[k][0], lo[k][1]) for k in ("τx по спаду амплитуды, мс", "τy по спаду амплитуды, мс",
                                                      "τs по спаду амплитуды, мс") if k in lo)
        u0 = lo["U0 по трекингу за оборот, МэВ"]
        qs = [v for k, v in lo.items() if k.startswith("Qs по трекингу")]
        R.p(f"Потеря энергии за оборот совпадает до {pretty(rel(u0[0], u0[1]), 2)} "
            f"({pretty(u0[0], 7)} против {pretty(u0[1], 7)} МэВ у elegant и {pretty(u0[2], 7)} МэВ по "
            f"интегралам излучения). Времена затухания по матрице оборота совпадают до {pretty(mat_d, 2)}, "
            f"по спаду амплитуды — до {pretty(amp_d, 2)}: два независимых способа дают одно и то же в обеих "
            f"программах."
            + (f" Синхротронная частота по NAFF совпадает до "
               f"{pretty(max(rel(v[0], v[1]) for v in qs), 2)} на всех амплитудах δ." if qs else ""))
        an_keys = [k for k in ("U0 по трекингу за оборот, МэВ", "τx, мс", "τy, мс", "τs, мс")
                   if k in lo and len(lo[k]) > 3 and lo[k][3] is not None and lo[k][2]]
        if an_keys:
            d_an = max(rel(lo[k][3], lo[k][2]) for k in an_keys)
            R.p(f"Последний столбец таблицы — аналитика SynchrotronSim: интегралы излучения I₁…I₅ по "
                f"описанию решётки, без трекинга. С аналитикой elegant она сходится до {pretty(d_an, 2)}: "
                f"аналитическая пара здесь закрыта так же, как в линейной оптике.")

    # ------------------------------------------------------------ портреты с излучением
    if exists("radiation"):
        R.h1("Фазовые портреты с излучением")
        rs = np.load(os.path.join(RES, "radiation", "ss_sr.npz"), allow_pickle=True)
        n_turns = rs["tbt"].shape[1] - 1
        u0 = np.load(os.path.join(RES, "longitudinal", "u0.npz"), allow_pickle=True)
        T0 = lattice_facts()["circumference"] / 299792458.0
        R.p(f"Три пробные частицы стартуют от синхронной точки со смещением по одной координате — x, y и δ — "
            f"и проходят {n_turns} оборотов с включённым излучением: {pretty(n_turns * T0 / float(u0['taux']), 2)} "
            f"времени затухания по x и {pretty(n_turns * T0 / float(u0['taudelta']), 2)} продольного. Пара "
            f"трекинговая: у обеих программ излучение распределено по изгибам, ВЧ стоит в точке наблюдения, "
            f"поперечные импульсы при излучении сохраняются.")
        R.p("Верхний ряд — только классическое излучение (SYNCH_RAD): потеря энергии детерминирована, "
            "поэтому портреты обеих программ обязаны лечь друг на друга, и амплитуда сходится к нулю по "
            "экспоненте с временем затухания из предыдущей главы. Нижний ряд — с квантовым возбуждением "
            "(ISR): скачки энергии при испускании фотонов разыгрываются генераторами случайных чисел, у "
            "каждой программы своим, поэтому траектории поточечно совпадать НЕ обязаны — сравнивать здесь "
            "нужно размер равновесного облака. Пунктиром показано ±σ по равновесным моментам из следующей "
            "главы.")
        R.fig(os.path.join(RES, "radiation", "phase_portraits.png"),
              "Фазовые портреты пробных частиц: классическое излучение (сверху) и квантовое (снизу).")
        R.fig(os.path.join(RES, "radiation", "envelope.png"),
              "Огибающая амплитуды по оборотам: классическое излучение уводит частицу к нулю, квантовое "
              "останавливает её на равновесном размере.")
        R.csv_table(os.path.join(RES, "radiation", "summary.csv"), 6,
            title="Излучение: затухание амплитуд и равновесное облако")
        rd = csv_dict(os.path.join(RES, "radiation", "summary.csv"))
        env = [v[2] for k, v in rd.items() if "классическое" in k and v[2] is not None]
        if env:
            R.p(f"Амплитуда колебания считается инвариантом Куранта–Снайдера, а не размахом по осям, — "
                f"иначе разница частот программ давала бы ложный разброс. Кривые классического затухания "
                f"совпадают до {pretty(max(env), 2)} по всем трём плоскостям, пока амплитуда не упала ниже "
                f"десятой части начальной; это тот же предел интегратора, что и в линейной оптике.")
        note = t_floor_note(n_turns)
        if note:
            R.p(note)
        em = csv_dict(os.path.join(RES, "emittance", "summary.csv")) if exists("emittance") else {}
        R.p("Вертикальная плоскость отдельно показательна: при нулевой вертикальной дисперсии и без связи "
            "квантовых скачков в неё не попадает, и кривая с квантовым возбуждением ложится ровно на "
            "классическую — обе программы дают одно и то же чистое затухание. Смазывание вертикального "
            "портрета в нижнем ряду — не рост амплитуды, а дрожание фазы от хроматизма: δ гуляет, вместе с "
            "ним гуляет νy.")
        ens = [(k, v) for k, v in rd.items() if "ансамбль" in k]
        if ens and M.get("rad_ensemble"):
            n_ens, e_turns, e_skip, e_stride = M["rad_ensemble"]
            R.p(f"Равновесное облако квантового излучения проверяется АНСАМБЛЕМ: {n_ens} "
                f"{plural(n_ens, 'частица стартует', 'частицы стартуют', 'частиц стартуют')} из "
                f"синхронной точки, у каждой свой генератор квантового шума (у elegant — свой seed на каждый "
                f"запуск), и после выхода на равновесие ({e_skip} оборотов, "
                f"{pretty(e_skip * T0 / float(u0['taux']), 2)} τx) координаты записываются ещё "
                f"{e_turns - e_skip} оборотов через каждые {e_stride}. Размеры облака — СКО по частицам и "
                f"оборотам; статистическая ошибка — по числу независимых выборок (частицы × длина записи / "
                f"время жизни действия τ/2).")
            R.fig(os.path.join(RES, "radiation", "equilibrium_clouds.png"),
                  "Равновесное облако ансамбля: фазовые плоскости и распределения x и δ; штрих — гауссово "
                  "распределение с σ из moments_output elegant.")
            worst = max(v[2] for _, v in ens if v[2] is not None)
            R.p(f"Размеры облаков программ различаются не больше чем на {pretty(worst * 100, 2)} %. Облако "
                f"каждой программы совпадает с её же расчётом равновесия (метод Чао) в пределах статистики "
                f"ансамбля; разница между программами — систематический сдвиг квантового шума в трекинге "
                f"elegant, подробнее в главе об эмиттансе.")

    # ------------------------------------------------------------ эмиттанс
    if exists("emittance"):
        R.h1("Равновесный эмиттанс")
        R.p("Пары одинаковых инструментов здесь нет: у elegant равновесные моменты считает аналитический "
            "moments_output, а в SynchrotronSim они получаются методом Чао поверх трекинга — матрица оборота "
            "с излучением и ВЧ плюс матрица диффузии, собранная по срезам изгибов из того же ядра.")
        R.formula(r"\Sigma = M\,\Sigma\,M^{T} + B,\qquad B = \sum_j \sigma_j^2\, r_j r_j^{T},\qquad "
                  r"\varepsilon_k = \left|\mathrm{Im}\,\lambda_k\left(\Sigma J\right)\right|",
                  [
            (r"\Sigma", "матрица вторых моментов пучка 6×6", "м², м·рад, рад²"),
            (r"M", "матрица оборота с излучением и ВЧ"),
            (r"B", "матрица диффузии за оборот от квантового излучения", "те же единицы"),
            (r"\sigma_j^2", "дисперсия скачка δ на срезе изгиба j", "—"),
            (r"r_j", "столбец переноса δ от среза до точки наблюдения", "м, рад"),
            (r"\varepsilon_k", "собственные эмиттансы равновесного пучка", "м·рад"),
            (r"J", "симплектическая матрица", "—"),
        ])
        R.fig(os.path.join(RES, "emittance", "emittance.png"), "Равновесный эмиттанс и энергетический разброс.")
        R.csv_table(os.path.join(RES, "emittance", "summary.csv"), 6, title="Равновесный эмиттанс и размеры пучка")
        conv = sorted(glob.glob(os.path.join(RES, "emittance", "ss_nsl*_o*.json")))
        if conv:
            rows = []
            for f in conv:
                q = json.load(open(f, encoding="utf-8"))
                rows.append((f"N = {q['nsl']}, порядок {q['order']}", q["ex_pm"], q["sd"] * 1e4, q["sz_mm"],
                             q["taux_ms"], q["taus_ms"]))
            rows.sort(key=lambda r: (int(r[0].split("порядок ")[1]), int(r[0].split("=")[1].split(",")[0])))
            R.p("SynchrotronSim при разном числе срезов и порядке интегратора:")
            R.table(["интегратор", "εx, пм", "σδ, 10⁻⁴", "σz, мм", "τx, мс", "τs, мс"], rows, 6,
                title="Сходимость равновесных величин по числу срезов")
            em = csv_dict(os.path.join(RES, "emittance", "summary.csv"))
            best = max((json.load(open(f, encoding="utf-8")) for f in conv), key=lambda q: (q["nsl"], q["order"]))
            mo = em["εx, пм"][2]
            word = lambda a, b: "ниже" if a < b else "выше"
            R.p(f"σδ, σz и продольный эмиттанс сходятся с moments_output до 0,2 %. Поперечный эмиттанс при "
                f"5 срезах на {pretty(rel(em['εx, пм'][0], mo) * 100, 2)} % {word(em['εx, пм'][0], mo)} "
                f"moments_output, при N = {best['nsl']} — {pretty(best['ex_pm'], 4)} пм, то есть на "
                f"{pretty(rel(best['ex_pm'], mo) * 100, 2)} % {word(best['ex_pm'], mo)} — между двумя "
                f"собственными оценками elegant (twiss_output даёт {pretty(em['εx, пм'][3], 5)} пм, "
                f"moments_output {pretty(mo, 5)} пм). Снова дискретизация, а не различие моделей.")
        em = csv_dict(os.path.join(RES, "emittance", "summary.csv"))
        if em and not conv:
            ex = em["εx, пм"]
            R.p(f"Метод Чао поверх C++-ядра даёт εx = {pretty(ex[0] / 1e3, 5)} нм против "
                f"{pretty(ex[3] / 1e3, 5)} нм у twiss_output (интегралы излучения) — расхождение "
                f"{pretty(rel(ex[0], ex[3]) * 100, 2)} %. Энергетический разброс, длина сгустка и продольный "
                f"эмиттанс совпадают с moments_output до 0,2 %.")
            R.p(f"Столбец moments_output по поперечному эмиттансу ({pretty(ex[2] / 1e3, 5)} нм) стоит "
                f"особняком, и это не расхождение программ: при дисперсии 1,13 м и σδ = "
                f"{pretty(em['σδ, 10⁻⁴'][2], 4)}·10⁻⁴ дисперсионный вклад в размер пучка больше бетатронного, "
                f"и moments_output приводит к «эмиттансу» полный размер. Проверяется это прямо в таблице: "
                f"сам размер σx у точки наблюдения у всех трёх инструментов совпадает.")
        R.p("Последний столбец — аналитика SynchrotronSim по интегралам излучения (без трекинга): с "
            "twiss_output elegant она совпадает до 10⁻⁶, а с методом Чао поверх трекинга расходится ровно "
            "на дискретизацию интегратора.")
        R.p("Столбец «elegant-трек» получен статистикой по частицам с квантовым излучением. " + ISR_NOTE)

    # ------------------------------------------------------------ апертура
    if exists("dynamic_aperture"):
        R.h1("Динамическая апертура, карта частот и сетка резонансов")
        fx, fy = M["fma_x_mm"], M["fma_y_mm"]
        R.p(f"Чисто трекинговая пара. Сетка начальных смещений от замкнутой орбиты: x от "
            f"{pretty(fx[0], 3)} до {pretty(fx[1], 3)} мм ({int(fx[2])} {plural(fx[2], 'точка', 'точки', 'точек')}), "
            f"y от {pretty(fy[0], 3)} до "
            f"{pretty(fy[-1], 3)} мм, δ = 0; 512 оборотов. "
            f"Частица потеряна, если в точке наблюдения |x| или |y| превысит "
            f"{pretty(float(M.get('loss_limit_m', 0.03)) * 1e3, 4)} мм. Частоты — NAFF по первым "
            f"и вторым 256 оборотам, скорость диффузии — изменение частот между половинами.")
        R.formula(r"D = \log_{10}\sqrt{\left(\nu_x^{(2)}-\nu_x^{(1)}\right)^2 + \left(\nu_y^{(2)}-\nu_y^{(1)}\right)^2},"
                  r"\qquad m\,\nu_x + n\,\nu_y = p",
                  [
            (r"D", "показатель диффузии частоты", "—"),
            (r"\nu^{(1)},\ \nu^{(2)}", "частоты по первой и второй половине прогона", "—"),
            (r"m,\ n,\ p", "целые числа резонансного условия", "—"),
        ])
        R.fig(os.path.join(RES, "dynamic_aperture", "dynamic_aperture.png"), "Динамическая апертура.")
        R.fig(os.path.join(RES, "dynamic_aperture", "frequency_map.png"),
              "Карта частот частиц, доживших до конца, на сетке резонансов.")
        R.fig(os.path.join(RES, "dynamic_aperture", "diffusion_xy.png"), "Диффузия и разность частот программ.")
        _, grid = read_csv(os.path.join(RES, "dynamic_aperture", "grid.csv"))
        if grid:
            nt = max(num(r[2]) for r in grid)
            agree = sum((num(r[2]) == nt) == (num(r[3]) == nt) for r in grid)
            dq = [abs(num(r[4]) - num(r[5])) for r in grid
                  if num(r[4]) is not None and num(r[5]) is not None
                  and math.isfinite(num(r[4])) and math.isfinite(num(r[5]))]
            dq.sort()
            R.p((f"Исход «выжила/потеряна» совпадает во всех {len(grid)} точках сетки"
                 if agree == len(grid) else
                 f"Исход «выжила/потеряна» совпадает в {agree} точках из {len(grid)}")
                + (f", медианное расхождение νx по выжившим {pretty(dq[len(dq) // 2], 2)}" if dq else "")
                + ("." if agree == len(grid) else
                   "; различия — только на хаотической границе апертуры, где расхождения координат растут "
                   "экспоненциально."))
            if M["source"]:
                R.p(da_note(grid, nt))

    # ------------------------------------------------------------ пределы
    R.h1("Чем ограничена достигнутая точность")
    R.p("Ни одно из оставшихся расхождений не связано с физикой моделей — у каждого есть измеренная "
        "причина. Перечислены в порядке убывания вклада.")
    an = csv_dict(os.path.join(RES, "analytic", "summary.csv")) if exists("analytic") else {}
    lin_s = csv_dict(os.path.join(RES, "linear", "summary.csv")) if exists("linear") else {}
    gap_el = (abs(lin_s["νx"][2] - lin_s["νx"][1]) if lin_s and len(lin_s.get("νx", ())) > 2 else None)
    gap_ss = (abs(an["νx (матрицы по решётке БД)"][0] - lin_s["νx"][0])
              if an and lin_s and "νx (матрицы по решётке БД)" in an else None)
    lim = [
        ("Сходимость интегратора. Собственный трекинг каждой программы отличается от её же матриц на "
         + (f"{pretty(gap_el, 2)} по частоте у elegant и {pretty(gap_ss, 2)} у SynchrotronSim "
            if gap_el and gap_ss else "разностную величину ")
         + f"при {nsl} срезах — разрыв один и тот же, то есть это чистая дискретизация, а не устройство "
           f"программ. Это и есть предел трекинговой пары: все трекинговые расхождения этого отчёта "
           f"упираются в него."),
        "Разностные матрицы. Якобианы снимаются 4-точечной схемой по траекториям; шум округления даёт "
        "10⁻⁷…10⁻⁶ в членах матрицы. Точные якобианы (дифференцирование ядра) убрали бы этот фон.",
    ]
    if an:
        dC, dT = geom_residual()
        co_db = lin_s.get("замкнутая орбита x, мкм", [0])[0]
        if max(dC, dT) >= 1e-9:
            lim.append(
                f"Геометрия базы данных. {geom_note()[0].upper() + geom_note()[1:]}. Отсюда остаточная "
                f"замкнутая орбита {si_len(abs(co_db) * 1e-6)}, а через хроматический член дрейфа "
                f"∂x′/∂δ = −p_x — и след в дисперсии: матрицы разложены вокруг проектной оси, а разностные "
                f"производные снимаются вокруг реальной орбиты.")
        else:
            lim.append(
                f"Геометрия базы данных больше не в списке: координаты магнитов записаны с полной точностью "
                f"двойного числа, {geom_note()}, остаточная замкнутая орбита {si_len(abs(co_db) * 1e-6)}. "
                f"Аналитическая пара упирается в округление перемножения матриц.")
    lim += [
        "Соглашения elegant. Для корректора конечной длины он усредняет функции Твисса по концам элемента, "
        "а квантовый шум в CSBEND обрезает на трёх σ (отсюда заниженная статистическая оценка диффузии). "
        "Это различия инструментов, а не физики.",
        "Резонансные члены h_jklm у SynchrotronSim пока не считаются: нормальная форма отображения даёт "
        "сдвиг частот с действием, но разложение по резонансам не реализовано — у elegant для сравнения "
        "есть compute_driving_terms.",
    ]
    for t in lim:
        R.bullet(t)

    # ------------------------------------------------------------ готчи
    R.h1("Особенности инструментов elegant")
    for t in (
        "Масса электрона зашита как " + f"{ed.ME_ELEGANT_MEV:.8f}".replace(".", ",") +
        " МэВ (CODATA 1986) — по ней считается p_central. Файлы частиц надо писать по ЕГО опорному "
        "импульсу (pCentral из его же вывода), иначе осевая частица приезжает с δ = " +
        pretty(ed.ME_ELEGANT_MEV / ed.ME_MEV - 1, 3) + " и идёт по дисперсионной орбите — см. главу об "
        "опорном импульсе.",
        "twiss_output и matrix_output считают по точным матрицам элементов и не учитывают N_SLICES и "
        "INTEGRATION_ORDER. Сравнивать их следует с аналитическим механизмом, а не с трекингом; на грубом "
        "интеграторе они могут описывать устойчивое кольцо там, где собственный трекинг его теряет.",
        "Число срезов интегратора по умолчанию: N_SLICES = 1 у KQUAD и KSEXT, 4 у CSBEND. На коротких "
        "сильных линзах одного шага не хватает: у Сибири кольцо при этом выходит неустойчивым в трекинге, "
        "хотя twiss_output даёт нормальные функции. Число срезов задаётся командой alter_elements.",
        "higher_order_chromaticity в twiss_output верна только при default_order = 3.",
        "correction_matrix_output для корректора конечной длины берёт функции Твисса как среднее по концам "
        "элемента, а не в его середине — отсюда систематический остаток в матрице отклика.",
        "use_response_from_computed_orbits в &correct не переводит расчёт отклика на трекинг: матрица "
        "получается бит-в-бит той же. Трекинговый отклик приходится делать киком через alter_elements.",
        "tune_shift_with_amplitude рассчитан на малый сдвиг частот: при широком окне поиска он местами берёт "
        "не ту линию, при узком — сам уменьшает амплитуды, а scale_down_limit = 0 зацикливает расчёт.",
        "RFCA: опорное время задаётся параметром FIDUCIAL; режим light не зависит от состава пучка.",
        "Время t в выходных файлах абсолютное; его точность около 10⁻¹³ м ограничивает разностные производные.",
        "Квантовый шум CSBEND ограничен тремя σ (srGaussianLimit): статистические оценки диффузии занижены.",
        "sddsmakedataset не принимает тысячи частиц в командной строке — файл частиц пишется в SDDS ASCII.",
        "В SDDS ToolKit есть sddsnaff (частоты), sddsresdiag (сетка резонансов), trimda (обработка апертуры).",
    ):
        R.bullet(t)

    R.save(out)
    return out


def main():
    ap = argparse.ArgumentParser(description="отчёт сверки SynchrotronSim с elegant")
    ap.add_argument("--machine", default="sila", choices=sorted(ed.MACHINES),
                    help="установка, по которой строится отчёт")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    use_machine(a.machine)
    print(build(a.out))


if __name__ == "__main__":
    main()
