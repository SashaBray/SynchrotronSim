# -*- coding: utf-8 -*-
"""Запуск elegant и чтение SDDS для сверки SynchrotronSim с elegant (ring_analysis ↔ elegant).

Рабочая папка — elegant_run/diag. Решётка — syla_closed.lte, собирается из той же конфигурации, что
загружена в решётку 1 БД (workspace/lte/configuration_v1_closed): описание устройств + линия RING.
К линии можно добавить точку наблюдения (WATCH) и ВЧ-резонатор (RFCA) — они не меняют оптику.

Координаты elegant: x, xp (наклон), y, yp, t (с), p (βγ). У SynchrotronSim px — канонический импульс
(1+δ)·x'/√…; в линейном приближении на опорной орбите это одно и то же, для больших амплитуд
переводится через to_elegant/from_elegant.
"""
import os, re, math, shutil, subprocess
import threading, itertools
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
DIAG = os.path.join(ROOT, "elegant_run", "diag")
# configuration_v1_geo = configuration_v1_closed + замыкание по ПОЛОЖЕНИЮ (слабина 0.1125 мм убрана в
# длинный промежуток): без неё решётка БД была на 0.11 мм длиннее линии .lte, и это давало фон 1e-5.
CONFIG = os.environ.get("SS_ELEGANT_CONFIG") or os.path.join(ROOT, "workspace", "lte", "configuration_v1_geo")
E_MEV = 6000.0                                    # проектная энергия СИЛА (radiation_params.csv)
ME_MEV = 0.51099895069                            # энергия покоя электрона, CODATA 2018
ME_ELEGANT_MEV = 0.51099906                       # та же постоянная ВНУТРИ elegant (CODATA 1986)
C_LIGHT = 299792458.0

# Установки сверки. СИЛА собирается из описания+линии (можно подменять устройства и интегратор), чужие
# решётки берутся файлом целиком: у них своя линия, которую остаётся обернуть точкой наблюдения.
#   lattice — номер решётки в БД, e_mev — проектная энергия, source/beamline — готовый .lte и имя линии
#   в нём, use_beamline — имя линии, которую собирает elegant_diag (у чужого файла своё имя RING может быть
#   занято), lte_name — файл для встроенных инструментов elegant, superperiods — сколько супериодов в кольце,
#   slices — (число срезов интегратора, порядок) для ОБОИХ трекеров,
#   rf — ВЧ (кратность и напряжение), correctors/corrector_len — чем киковать в опыте отклика орбиты,
#   rad_amp — (x, мм; y, мм; δ) начальные отклонения пробных частиц в опыте с излучением,
#   amp_*_mm и fma_*_mm — амплитуды пооборотного трекинга и сетка апертуры.
MACHINES = {
    "sila": dict(label="СИЛА", lattice=1, e_mev=6000.0, config=CONFIG, source=None, beamline=None,
                 use_beamline="RING", lte_name="syla_closed.lte", superperiods=40, slices=(5, 4),
                 rf=dict(volt_mv=8.0, harmonic=1296), correctors=("SH1A", "SH2B", "SH3E"),
                 corrector_len=0.1, rad_amp=(1.0, 1.0, 2e-3),
                 amp_x_mm=(0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
                 amp_y_mm=(0.5, 1.0, 2.0, 3.0, 4.0), fma_x_mm=(-10.0, 10.0, 21),
                 # сетка симметрична по y: у решётки без скью-элементов движение строго симметрично
                 # относительно медианной плоскости, и нижняя половина карты — зеркало верхней
                 fma_y_mm=tuple(sorted({s * y for s in (-1, 1)
                                        for y in (0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)})),
                 loss_limit_m=0.03,
                 rad_turns=8192, rad_ensemble=None),
    # СИЛА, редакция 2026.06.11 — та же машина 1110 м, но 320 квадруполей стали комбинированными
    # (антибенды QF4*/QF6*/QD5*): εx 71 → 42 пм, Jx 1.52 → 2.18. Решётка 26 БД, импорт того же файла.
    "sila2026": dict(label="СИЛА 2026.06.11", lattice=26, e_mev=6000.0, config=None, beamline="ring",
                     use_beamline="SLRING",     # имя RING занято линией самого файла (elegant без регистра)
                     source=os.path.join(ROOT, "workspace", "lte", "configuration_v5_sila2026",
                                         "sila2026.lte"),
                     lte_name="sila2026.lte", superperiods=40,
                     # 5 срезов (умолчание прежней СИЛА) даёт разрыв уровней 1.3e-2 по νx — на этой
                     # оптике мало; 40 срезов дают 4.3e-6 и сходятся как N⁻⁴, поэтому обоим кодам 40
                     slices=(40, 4),
                     rf=dict(volt_mv=8.0, harmonic=1296),
                     correctors=("SH1A", "SH2B", "SH3E"), corrector_len=0.134,
                     rad_amp=(1.0, 1.0, 2e-3),
                     amp_x_mm=(0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0),
                     amp_y_mm=(0.5, 1.0, 2.0, 3.0, 4.0), fma_x_mm=(-10.0, 10.0, 21),
                     fma_y_mm=tuple(sorted({s * y for s in (-1, 1)
                                            for y in (0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0)})),
                     loss_limit_m=0.03,
                     rad_turns=8192, rad_ensemble=None),
    "siberia": dict(label="Сибирь", lattice=2, e_mev=3000.0, config=None, beamline="ring",
                    use_beamline="SSRING",      # имя RING занято линией самого файла (elegant без регистра)
                    source=os.path.join(ROOT, "workspace", "lte", "configuration_v3_syberia", "s2lattice.lte"),
                    lte_name="siberia.lte", superperiods=6,
                    # линзы короткие и сильные: одного шага интегратора (умолчание elegant) мало —
                    # его собственный трекинг при нём неустойчив, поэтому обоим кодам ставится 40
                    slices=(40, 4),
                    # ВЧ Сибири-2: при периметре 124.13 м кратность 75 даёт 181.1 МГц, U0 = 1.42 МэВ
                    rf=dict(volt_mv=2.5, harmonic=75),
                    # корректоров в файле нет: пробные корректоры (KICKER той же длины) ставятся в пустые
                    # промежутки SX и SZ — это места секступолей, у которых в файле K2 = 0 (дрейфы 0.113 м)
                    correctors=("SX", "SZ"), corrector_len=0.113, rad_amp=(1.0, 1.0, 3e-3),
                    # секступоли нулевые: кольцо почти линейно, и амплитуды берутся крупные — до них
                    # доживают все частицы, а сверяется кинематическая нелинейность (точные дрейфы, края)
                    amp_x_mm=(1.0, 2.0, 4.0, 6.0, 9.0, 12.0, 15.0), amp_y_mm=(1.0, 2.0, 4.0, 6.0, 8.0),
                    # устойчивость теряется только от кинематики — на 0.4–0.5 м (далеко за камерой); сетка
                    # апертуры берётся до 0.6 м, порог потери 1 м — теряются лишь разошедшиеся частицы
                    fma_x_mm=(-600.0, 600.0, 41),
                    fma_y_mm=tuple(sorted({s * y for s in (-1, 1)
                                           for y in (10.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0,
                                                     400.0, 450.0, 500.0, 550.0, 600.0)})),
                    loss_limit_m=1.0,
                    # излучение: пробные частицы ≈ 9 τx и ≈ 19 τs; ансамбль с ISR — (частиц, обороты,
                    # выход на равновесие ≈ 4 τx, шаг выборки) — для равновесного облака
                    rad_turns=40000, rad_ensemble=(64, 40000, 18000, 20)),
}


def machine(name="sila"):
    """Описание установки по имени (см. MACHINES)."""
    try:
        return MACHINES[name]
    except KeyError:
        raise KeyError(f"установка {name!r} не описана; есть {', '.join(MACHINES)}")


_APS = [r"C:\Program Files\APS\Elegant x64", r"C:\Program Files\APS\SDDS ToolKit x64"]


def _exe(name):
    p = shutil.which(name)
    if p:
        return p
    for d in _APS:
        q = os.path.join(d, name + ".exe")
        if os.path.exists(q):
            return q
    raise FileNotFoundError(f"не найден {name} (elegant / SDDS ToolKit)")


# ------------------------------------------------------------------ решётка
def _desc_text(config=None, source=None):
    """Текст с ОПРЕДЕЛЕНИЯМИ устройств: описание СИЛА или чужой .lte целиком."""
    if source:
        return open(source, encoding="utf-8").read()
    return open(os.path.join(config or CONFIG, "SYLA_desc.lte"), encoding="utf-8").read()


def line_names(config=None, source=None, beamline=None):
    """Последовательность имён устройств линии. У СИЛА это отдельный файл, у чужой решётки —
    определение `<beamline>: LINE=(...)` в самом .lte (переносы через «&», вложенных линий нет)."""
    if source:
        txt = re.sub(r"!.*", "", open(source, encoding="utf-8").read())
        m = re.search(rf"^[ \t]*{beamline}[ \t]*:[ \t]*LINE[ \t]*=[ \t]*\((.*?)\)", txt, re.M | re.S | re.I)
        if not m:
            raise KeyError(f"линия {beamline!r} не найдена в {os.path.basename(source)}")
        return [n for n in re.split(r"[,&\s]+", m.group(1)) if n]
    txt = open(os.path.join(config or CONFIG, "SYLA_line.lte"), encoding="utf-8").read()
    return [n for n in re.split(r"[,\s]+", re.sub(r"!.*", "", txt)) if n]


def element_params(config=None, source=None):
    """{имя (в верхнем регистре): {type, L, ANGLE, K1, …}} из описания устройств."""
    out = {}
    txt = _desc_text(config, source)
    for m in re.finditer(r"^([A-Za-z0-9_]+)[ \t]*:[ \t]*([A-Za-z]+)[ \t]*(.*)$", txt, re.M):
        d = {"type": m.group(2)}
        for kv in m.group(3).split(","):
            if "=" in kv:
                k, v = kv.split("=", 1)
                try:
                    d[k.strip()] = float(v)
                except ValueError:
                    d[k.strip()] = v.strip().strip('"')
        out[m.group(1).upper()] = d
    return out


def line_length(config=None, source=None, beamline=None):
    """Длина линии по описанию устройств, м (периметр для частоты ВЧ)."""
    par = element_params(config, source)
    return float(sum(par[n.upper()].get("L", 0.0) for n in line_names(config, source, beamline)))


def machine_lattice(m):
    """(config, source, beamline) описания установки — аргументы line_names/element_params/lead_drift."""
    return m.get("config"), m.get("source"), m.get("beamline")


def lead_drift(config=None, source=None, beamline=None):
    """Длина промежутка от начала линии (s = 0) до входа ПЕРВОГО магнита, м.

    Точка наблюдения SynchrotronSim ставится ровно сюда: тогда продольная координата s совпадает с s
    elegant поэлементно. Половина замыкающего промежутка для этого не годится — шов несимметричен."""
    par = element_params(config, source)
    s = 0.0
    for n in line_names(config, source, beamline):
        q = par[n.upper()]
        if q["type"].upper() not in ("EDRIFT", "DRIF", "MONI", "MARK", "WATCH"):
            return s
        s += q.get("L", 0.0)
    raise ValueError("в линии нет магнитов")


_INTEGRATED = "CSBEND|KQUAD|KSEXT|KOCT"


def _patch_defs(desc, overrides=None, integrator=None, csbend_extra=None):
    """Правки в ТЕКСТЕ определений устройств: подмена элементов, число срезов и порядок интегратора,
    добавка ко всем изгибам. Если N_SLICES в определении нет (чужая решётка) — дописывается."""
    if overrides:
        for nm, definition in overrides.items():
            desc, n = re.subn(rf"^{nm}[ \t]*:.*$", f"{nm}: {definition}", desc, flags=re.M | re.I)
            if not n:
                raise KeyError(nm)
    if integrator:                                   # (N срезов, порядок) для всех интегрируемых элементов
        nsl, order = integrator

        def fix(m):
            body = m.group(1)
            add = f",INTEGRATION_ORDER={int(order)}"
            if re.search(r"N_SLICES[ \t]*=", body, re.I):
                body = re.sub(r"N_SLICES[ \t]*=[ \t]*\d+", f"N_SLICES={int(nsl)}", body, flags=re.I)
            else:
                add = f",N_SLICES={int(nsl)}" + add
            return body + add

        desc = re.sub(rf"^([A-Za-z0-9_]+[ \t]*:[ \t]*(?:{_INTEGRATED})\b.*)$", fix, desc, flags=re.M | re.I)
    if csbend_extra:                                 # напр. "SYNCH_RAD=1,ISR=1" для всех изгибов
        desc = re.sub(r"^([A-Za-z0-9_]+[ \t]*:[ \t]*CSBEND\b.*)$", lambda m: m.group(1) + "," + csbend_extra,
                      desc, flags=re.M | re.I)
    return desc


def _write_foreign(name, source, beamline, watch=None, rf=None, extra_elements="", use_beamline="SSRING",
                   overrides=None, integrator=None, csbend_extra=None):
    """Чужой .lte целиком + линия use_beamline = (WATCH, <линия файла>, RFCA)."""
    txt = _patch_defs(open(source, encoding="utf-8").read().rstrip(), overrides, integrator, csbend_extra)
    defs = [f"! собрано elegant_diag.py из {os.path.relpath(source, ROOT)}", txt]
    head = []; tail = []
    if watch:
        defs.append(f'W0: WATCH, FILENAME="{watch["filename"]}", MODE="{watch.get("mode", "coordinate")}", '
                    f'INTERVAL={watch.get("interval", 1)}')
        head.append("W0")
    if rf:
        c = rf.get("circumference") or line_length(source=source, beamline=beamline)
        freq = rf["harmonic"] * C_LIGHT / c
        defs.append(f'RF0: RFCA, L=0, VOLT={rf["volt_mv"] * 1e6:.12g}, FREQ={freq:.15g}, '
                    f'PHASE={rf["phase_deg"]:.12g}, FIDUCIAL="light"')
        tail.append("RF0")
    if extra_elements:
        defs.append(extra_elements)
    defs.append(f"{use_beamline}: LINE=(" + ",".join(head + [beamline] + tail) + ")")
    path = os.path.join(DIAG, name)
    open(path, "w", encoding="utf-8", newline="\n").write("\n".join(defs) + "\n")
    return path


def write_lattice(name="syla_closed.lte", config=CONFIG, watch=None, rf=None, extra_elements="",
                  overrides=None, integrator=None, csbend_extra=None, source=None, beamline=None,
                  use_beamline="RING"):
    """Собрать .lte: описание + RING. watch — dict(filename, mode, interval) ⇒ WATCH в начале линии;
    rf — dict(volt_mv, harmonic, phase_deg) ⇒ RFCA в конце линии (FREQ = h·c/C, FIDUCIAL="light");
    integrator — (N_SLICES, INTEGRATION_ORDER) для всех интегрируемых элементов; csbend_extra — добавка к
    каждому CSBEND (напр. "SYNCH_RAD=1,ISR=1"); overrides — {имя: 'строка определения'} заменяет устройства.

    source/beamline — ЧУЖАЯ решётка: файл берётся целиком, RING оборачивает его линию beamline."""
    os.makedirs(DIAG, exist_ok=True)
    if source:
        return _write_foreign(name, source, beamline, watch, rf, extra_elements, use_beamline,
                              overrides, integrator, csbend_extra)
    desc = _patch_defs(open(os.path.join(config, "SYLA_desc.lte"), encoding="utf-8").read().strip(),
                       overrides, integrator, csbend_extra)
    names = line_names(config)
    head = []; tail = []
    defs = [f"! собрано elegant_diag.py из {os.path.relpath(config, ROOT)}", desc]
    if watch:
        defs.append(f'W0: WATCH, FILENAME="{watch["filename"]}", MODE="{watch.get("mode", "coordinate")}", '
                    f'INTERVAL={watch.get("interval", 1)}')
        head.append("W0")
    if rf:                                           # в КОНЦЕ линии: запись WATCH до ВЧ, выход — после (как SynchrotronSim)
        freq = rf["harmonic"] * C_LIGHT / line_length(config)
        defs.append(f'RF0: RFCA, L=0, VOLT={rf["volt_mv"] * 1e6:.12g}, FREQ={freq:.15g}, '
                    f'PHASE={rf["phase_deg"]:.12g}, FIDUCIAL="light"')
        tail.append("RF0")
    if extra_elements:
        defs.append(extra_elements)
    seq = head + names + tail
    rows = [",".join(seq[i:i + 12]) for i in range(0, len(seq), 12)]
    defs.append("RING: LINE=(" + ",&\n".join(rows) + ")")
    path = os.path.join(DIAG, name)
    open(path, "w", encoding="utf-8", newline="\n").write("\n".join(defs) + "\n")
    return path


# ------------------------------------------------------------------ запуск и SDDS
def run(ele_name, text, timeout=None):
    """Записать .ele и запустить elegant в DIAG. Возврат stdout (бросает при ошибке)."""
    path = os.path.join(DIAG, ele_name)
    open(path, "w", encoding="utf-8", newline="\n").write(text)
    r = subprocess.run([_exe("elegant"), ele_name], cwd=DIAG, capture_output=True, text=True,
                       timeout=timeout, encoding="utf-8", errors="replace")
    log = os.path.splitext(path)[0] + ".log"
    open(log, "w", encoding="utf-8").write(r.stdout + r.stderr)
    if r.returncode != 0 or "Thanks for using elegant" not in r.stdout:
        raise RuntimeError(f"elegant {ele_name}: код {r.returncode}\n{r.stdout[-2000:]}\n{r.stderr[-1000:]}")
    return r.stdout


def columns(sdds, cols, page=None):
    """Столбцы SDDS → dict массивов (числа — float, строки — str). page=None ⇒ все страницы подряд."""
    args = [_exe("sdds2stream"), sdds, "-col=" + ",".join(cols)]
    if page is not None:
        args.append(f"-page={page}")
    txt = subprocess.run(args, cwd=DIAG, capture_output=True, text=True, check=True).stdout
    rows = [ln.split() for ln in txt.splitlines() if ln.strip()]
    out = {}
    for k, c in enumerate(cols):
        vals = [r[k] for r in rows]
        try:
            out[c] = np.array([float(v) for v in vals])
        except ValueError:
            out[c] = vals
    return out


def parameters(sdds, pars, page=1):
    args = [_exe("sdds2stream"), sdds, "-par=" + ",".join(pars), f"-page={page}"]
    txt = subprocess.run(args, cwd=DIAG, capture_output=True, text=True, check=True).stdout.split()
    out = {}
    for k, v in zip(pars, txt):
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def n_pages(sdds):
    txt = subprocess.run([_exe("sdds2stream"), sdds, "-npages"], cwd=DIAG, capture_output=True,
                         text=True, check=True).stdout.split()
    return int(txt[0])


def write_particles(name, xs, e_mev=None, pc=None):
    """Начальные частицы для &sdds_beam: xs — (n, 6) в координатах elegant (x, xp, y, yp, t, δ).
    Файл SDDS пишется напрямую в ASCII (командная строка sddsmakedataset на тысячах частиц не помещается).

    pc — ОПОРНЫЙ ИМПУЛЬС elegant (pCentral из его же вывода). Его надо передавать: elegant считает
    p_central по своей массе электрона (0.51099906 МэВ), а ME_MEV здесь — современное значение
    (0.51099895069). Разница 2.1e-7 по импульсу, и частица с δ = 0 у нас приезжает к elegant как
    частица с δ = 2.1e-7: получается дисперсионная «орбита» η·δ и сдвиг частот ξ·δ."""
    xs = np.atleast_2d(np.asarray(xs, float))
    pc = pc if pc else (e_mev if e_mev else E_MEV) / ME_MEV
    path = os.path.join(DIAG, name)
    head = ["SDDS1",
            "&column name=x, units=m, type=double, &end",
            "&column name=xp, type=double, &end",
            "&column name=y, units=m, type=double, &end",
            "&column name=yp, type=double, &end",
            "&column name=t, units=s, type=double, &end",
            "&column name=p, units=\"m$be$nc\", type=double, &end",
            "&column name=particleID, type=long, &end",
            "&data mode=ascii, &end",
            "! page number 1",
            str(len(xs))]
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write("\n".join(head) + "\n")
        for k, r in enumerate(xs):
            f.write(f"{r[0]:.17g} {r[1]:.17g} {r[2]:.17g} {r[3]:.17g} {r[4]:.17g} {pc * (1 + r[5]):.17g} {k + 1}\n")
    return path


# ------------------------------------------------------------------ координаты
def to_elegant(st):
    """SynchrotronSim (x, px, y, py, t_путь, δ) → elegant (x, xp, y, yp, t_путь, δ) — наклоны вместо импульсов."""
    x, px, y, py, t, d = np.moveaxis(np.asarray(st, float), -1, 0)
    ps = np.sqrt((1 + d) ** 2 - px ** 2 - py ** 2)
    return np.stack([x, px / ps, y, py / ps, t, d], axis=-1)


def from_elegant(st):
    x, xp, y, yp, t, d = np.moveaxis(np.asarray(st, float), -1, 0)
    f = (1 + d) / np.sqrt(1 + xp ** 2 + yp ** 2)
    return np.stack([x, xp * f, y, yp * f, t, d], axis=-1)


_ME_ELEGANT = None                      # масса электрона ВНУТРИ elegant (МэВ) — узнаётся по его же pCentral
_ME_LOCK = threading.Lock()
_ME_SEQ = itertools.count()


def elegant_me(lte, beamline, e_mev):
    """Масса электрона, с которой elegant считает опорный импульс (p_central = E/m).

    У elegant она 0.51099906 МэВ (старое значение CODATA), у нас ME_MEV = 0.51099895069. Если писать
    частицы по нашей массе, elegant видит их с δ = 2.1e-7: появляется дисперсионная «орбита» η·δ (0.24 мкм
    у Сибири) и сдвиг частот ξ·δ (3.6e-6). Поэтому опорный импульс берётся у самого elegant — один
    служебный запуск twiss_output на процесс."""
    global _ME_ELEGANT
    with _ME_LOCK:                                   # запуски идут потоками: служебный файл должен быть свой
        if _ME_ELEGANT is not None:
            return _ME_ELEGANT
        stem = f"pcentral{os.getpid()}_{next(_ME_SEQ)}"
        run(stem + ".ele",
            f'&run_setup lattice="{lte}", use_beamline="{beamline}", p_central_mev={e_mev:g}, '
            f'default_order=1, &end\n&twiss_output filename="{stem}.twi", matched=1, &end\n')
        _ME_ELEGANT = e_mev / float(parameters(stem + ".twi", ["pCentral"])["pCentral"])
        if abs(_ME_ELEGANT - ME_ELEGANT_MEV) > 1e-9 * ME_ELEGANT_MEV:
            print(f"ВНИМАНИЕ: масса электрона у elegant {_ME_ELEGANT!r}, ожидалась {ME_ELEGANT_MEV!r}")
        for ext in (".ele", ".twi", ".log"):
            f = os.path.join(DIAG, stem + ext)
            try:
                os.remove(f)
            except OSError:
                pass
        return _ME_ELEGANT


# ============================================================================ трекер elegant
_MAGNET_TYPES_SKIP = {"EDRIFT", "DRIF", "MONI", "MARK", "WATCH", "RFCA", "MAXAMP"}


class ElegantRing:
    """Трекер elegant с интерфейсом ring_analysis.Ring: track / track_elements / track_tbt.

    Каждый вызов — отдельный запуск elegant в DIAG. Координаты на входе и выходе канонические
    (как у SynchrotronSim); t — путь в метрах (elegant хранит время, пересчёт через β0·c)."""

    label = "elegant"

    def __init__(self, tag="er", lattice_kwargs=None, run_setup_extra="", sr=False, isr=False, rf=None,
                 pre_commands="", machine_name="sila"):
        self.tag = f"{tag}{os.getpid()}"
        self.run_setup_extra = run_setup_extra
        self.pre_commands = pre_commands          # команды после &run_setup (напр. &alter_elements — кик в один элемент)
        self.machine = machine(machine_name)
        self.e_mev = float(self.machine["e_mev"])
        kw = dict(lattice_kwargs or {})
        self.use_beamline = self.machine["use_beamline"]
        if self.machine["source"]:                # чужая решётка — файл целиком и своя линия
            kw.setdefault("source", self.machine["source"])
            kw.setdefault("beamline", self.machine["beamline"])
            kw.setdefault("use_beamline", self.use_beamline)
        if sr or isr:
            kw["csbend_extra"] = "SYNCH_RAD=1" + (",ISR=1" if isr else "")
        if rf:
            kw["rf"] = rf
        self.watch_file = f"{self.tag}.w0"
        kw.setdefault("watch", dict(filename=self.watch_file, mode="coordinate", interval=1))
        self.lte = os.path.basename(write_lattice(f"{self.tag}.lte", **kw))
        self.pc = self.e_mev / elegant_me(self.lte, self.use_beamline, self.e_mev)
        self.beta0 = self.pc / math.sqrt(1 + self.pc ** 2)
        # опорные s и имена магнитов — по центроиду одной осевой частицы (там же читается pCentral)
        c = self._centroid([np.zeros(6)])[0]
        mask = self._mag_mask
        self.s_out = c["s"][mask][:-1]                            # без последней строки (конец оборота)
        self.names = [n for n, k in zip(c["ElementName"], mask) if k][:-1]
        self.circumference = float(c["s"][-1])
        self._n = 0

    # --- ввод/вывод
    def _fresh(self, stem):
        self._n = getattr(self, "_n", 0) + 1
        return f"{self.tag}_{stem}{self._n}"

    def _to_input(self, states):
        e = to_elegant(np.asarray(states, float))
        e[:, 4] = e[:, 4] / C_LIGHT                              # запаздывание, м → с
        return e

    def _from_cols(self, c, n_turns):
        """t elegant (с, абсолютное) → запаздывание относительно света за n_turns оборотов, м."""
        pc = self.pc
        st = np.stack([c["x"], c["xp"], c["y"], c["yp"], c["t"] * C_LIGHT - n_turns * self.circumference,
                       c["p"] / pc - 1], 1)
        return from_elegant(st)

    def _setup(self, extra=""):
        return (f'&run_setup lattice="{self.lte}", use_beamline="{self.use_beamline}", p_central_mev={self.e_mev:g}, '
                f'default_order=2, {self.run_setup_extra} {extra} &end\n' + (self.pre_commands or ""))

    def track(self, states, n_turns=1):
        stem = self._fresh("trk")
        write_particles(stem + ".in", self._to_input(states), pc=self.pc)
        run(stem + ".ele", self._setup(f'output="{stem}.out",') +
            f"&run_control n_passes={int(n_turns)}, &end\n"
            f'&sdds_beam input="{stem}.in", &end\n&track &end\n')
        c = columns(stem + ".out", ["x", "xp", "y", "yp", "t", "p", "particleID"])
        out = np.full((len(states), 6), np.nan)
        if len(c["x"]):
            out[c["particleID"].astype(int) - 1] = self._from_cols(c, n_turns)
        self._cleanup(stem, (".in", ".out", ".ele", ".log"))
        if os.path.exists(os.path.join(DIAG, self.watch_file)):
            os.remove(os.path.join(DIAG, self.watch_file))
        return out

    def _centroid(self, states):
        stem = self._fresh("cen")
        files = []
        for i, st in enumerate(states):
            write_particles(f"{stem}_{i}.pg", self._to_input([st]), pc=self.pc); files.append(f"{stem}_{i}.pg")
        if len(files) > 1:
            subprocess.run([_exe("sddscombine")] + files + [stem + ".in", "-overwrite"], cwd=DIAG, check=True)
        else:
            os.replace(os.path.join(DIAG, files[0]), os.path.join(DIAG, stem + ".in"))
        run(stem + ".ele", self._setup(f'centroid="{stem}.cen",') +
            f"&run_control n_passes=1, n_steps={len(states)}, &end\n"
            f'&sdds_beam input="{stem}.in", track_pages_separately=1, &end\n&track &end\n')
        pages = []
        for pg in range(1, len(states) + 1):
            pages.append(columns(stem + ".cen", ["s", "ElementName", "ElementType", "Cx", "Cxp", "Cy",
                                                  "Cyp", "Cs", "Cdelta", "Particles"], page=pg))
        self._mag_mask = np.array([t not in _MAGNET_TYPES_SKIP for t in pages[0]["ElementType"]])
        self._mag_mask[-1] = True                                 # последняя строка — конец оборота
        for f in files:
            if os.path.exists(os.path.join(DIAG, f)):
                os.remove(os.path.join(DIAG, f))
        self._cleanup(stem, (".in", ".cen", ".ele", ".log"))
        return pages

    def track_elements(self, states):
        pages = self._centroid(states)
        mask = self._mag_mask
        out = np.full((len(states), int(mask.sum()), 6), np.nan)
        for i, c in enumerate(pages):
            if c["Particles"][-1] < 1:
                continue
            st = np.stack([c["Cx"], c["Cxp"], c["Cy"], c["Cyp"], c["Cs"], c["Cdelta"]], 1)[mask]
            out[i] = from_elegant(st)
        return out

    def track_tbt(self, states, n_turns):
        """Пооборотные координаты у начала линии (WATCH): (n, n_turns+1, 6), NaN после потери."""
        stem = self._fresh("tbt")
        write_particles(stem + ".in", self._to_input(states), pc=self.pc)
        wf = self.watch_file
        run(stem + ".ele", self._setup(f'output="{stem}.out",') +
            f"&run_control n_passes={int(n_turns)}, &end\n"
            f'&sdds_beam input="{stem}.in", &end\n&track &end\n')
        out = np.full((len(states), n_turns + 1, 6), np.nan)
        # номер прохода — параметр страницы; переносим его в столбец, чтобы прочитать всё одним вызовом
        subprocess.run([_exe("sddsprocess"), wf, stem + ".wp", "-define=column,pass,Pass,type=long"],
                       cwd=DIAG, check=True, capture_output=True)
        c = columns(stem + ".wp", ["x", "xp", "y", "yp", "t", "p", "particleID", "pass"])
        if len(c["x"]):
            ps = c["pass"].astype(int)
            st = self._from_cols(c, 0)
            st[:, 4] -= ps * self.circumference
            out[c["particleID"].astype(int) - 1, ps] = st
        c = columns(stem + ".out", ["x", "xp", "y", "yp", "t", "p", "particleID"])
        if len(c["x"]):
            out[c["particleID"].astype(int) - 1, n_turns] = self._from_cols(c, n_turns)
        self._cleanup(stem, (".in", ".out", ".ele", ".log", ".wp"))
        if os.path.exists(os.path.join(DIAG, wf)):
            os.remove(os.path.join(DIAG, wf))
        return out

    def __del__(self):
        """Решётка и файл WATCH экземпляра — временные: убрать вместе с объектом."""
        if os.environ.get("SS_ELEGANT_KEEP"):
            return
        for f in (getattr(self, "lte", None), getattr(self, "watch_file", None)):
            try:
                if f and os.path.exists(os.path.join(DIAG, f)):
                    os.remove(os.path.join(DIAG, f))
            except OSError:
                pass

    def _cleanup(self, stem, exts):
        if os.environ.get("SS_ELEGANT_KEEP"):
            return
        for e in exts:
            p = os.path.join(DIAG, stem + e)
            if os.path.exists(p):
                os.remove(p)
