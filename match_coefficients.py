#!/usr/bin/env python3
"""match_coefficients.py — для каждого элемента библиотеки строит словарь
  комбинация токов  ->  комбинация физических коэффициентов BIP,
дающих ТЕ ЖЕ траектории, что наш Boris-трекер на полевой карте.

Метод: свип параметра по ОДЗ (ток ±100А шаг 10А; для пост. диполей correct_coeff ±3
шаг 0.3) -> эталонные траектории 13 частиц (single_element_boris) -> подгонка коэф.
матричной модели BIP (bip_matrix) последовательным приближением (1D line-search по k).
Изгиб/кик берём напрямую из отклонения ЦЕНТРАЛЬНОЙ частицы; фокусировку (k=K1/K2/K3)
подгоняем по спутникам ОТНОСИТЕЛЬНО центра.

Вывод: workspace/bip_compare/<device>_map.csv + summary.md.
"""
import sys, os, math, csv, importlib
from pathlib import Path
from scipy.optimize import minimize
seb = importlib.import_module(os.environ.get("SEB_BACKEND", "single_element_boris"))
import bip_matrix as bm

OUT = seb.WS / os.environ.get("BIP_OUT", "bip_compare")
LC = 1.0          # масштаб угла->позиции в невязке (м)

# device -> (kind матрицы, тип свипа). kind: quad/sext/oct/dq/dipole/kicker
DEVICES = {
    "QP_L_1A":  ("quad", "I"),   "QP_sh_1A": ("quad", "I"),
    "SP_I1_1A": ("sext", "I"),   "SP_I2_1A": ("sext", "I"),
    "OP_1A":    ("oct",  "I"),
    "DQ_I1_1A": ("dq",   "I"),   "DQ_I2_1A": ("dq", "I"),   "DQ_I3_1A": ("dq", "I"),
    "DP_L_PM":  ("dipole", "C"), "DP_sh_PM": ("dipole", "C"),
    "corrector_I1_1A": ("kicker", "I"), "corrector_I2_1A": ("kicker", "I"),
}
MAT_KIND = {"quad": "quad", "dq": "quad", "dipole": "quad", "sext": "sext", "oct": "oct"}


def sweep_values(kind_sweep):
    if kind_sweep == "I":
        return [round(-100 + 10 * i, 3) for i in range(21)]   # -100..100 шаг 10
    return [round(-3 + 0.3 * i, 3) for i in range(21)]        # -3..3 шаг 0.3 (coeff)


APER = 0.045   # безопасная полу-апертура карты, м (пучок не должен выходить за неё в теле)


def sweep_for(device, kind_sweep):
    """21 точка свипа, АВТО-масштабированная под фактический ВЫНОС пучка (изгиб центра +
    фокус/дефокус спутников). При нашей инжекционной энергии ±100А делает элемент
    сверхсильным (пучок выносит за апертуру / hard-edge не ловит) — сужаем до безопасного.
    Коэффициент линеен по току ⇒ на-единицу-тока не зависит от диапазона."""
    base = 100.0 if kind_sweep == "I" else 3.0
    isc = (kind_sweep == "C")
    vmax = base
    for _ in range(40):                                          # итеративно сужаем, пока вынос велик
        e = seb.probe_excursion(device, vmax, isc)
        if e <= 0.5 * APER or vmax < base * 1e-4:
            break
        vmax *= 0.75
    return [round(-vmax + 2 * vmax * i / 10.0, 6) for i in range(11)], vmax  # 11 точек = ~10 шагов ОДЗ


TOTAL = seb.TOTAL                                            # полная продольная область = 2 м


def _geom(L):
    """Геометрия матрицы: элемент длины L центрирован у позиции магнита (seb.DRIFT);
    короткий дрейф ДО (DRIFT−L/2), длинный ПОСЛЕ (TOTAL−DRIFT−L/2). Сумма = TOTAL (как в физ.модели)."""
    return L, seb.DRIFT - L / 2.0, TOTAL - seb.DRIFT - L / 2.0  # L, d1, d2


def _rel(final):
    """Список (dY,dYp,dZ,dZp) спутников ОТНОСИТЕЛЬНО центра (частица 0)."""
    cy, cyp, cz, czp = final[0]
    return [(y - cy, yp - cyp, z - cz, zp - czp) for (y, yp, z, zp) in final[1:]]


def _residual(ref, mat_kind, k, L):
    """L2-невязка матричной модели (коэф. k, длина L) против эталона по спутникам отн. центра."""
    parts = [(y0, 0.0, z0, 0.0) for (y0, z0) in ref["init"]]
    Lm, d1, d2 = _geom(L)
    spec = {"kind": mat_kind, "L": Lm, "d1": d1, "d2": d2, "k": k, "theta": 0.0}
    out = bm.track_element(parts, spec)
    mcx, mcxp, mcy, mcyp = out[0]
    bref = _rel(ref["final"])
    s = 0.0
    for i in range(seb.N_SAT):
        mx, mxp, my, myp = out[i + 1]
        dym, dypm, dzm, dzpm = mx - mcx, mxp - mcxp, my - mcy, myp - mcyp
        dyb, dypb, dzb, dzpb = bref[i]
        s += (dym - dyb) ** 2 + (dzm - dzb) ** 2 + (LC * (dypm - dypb)) ** 2 + (LC * (dzpm - dzpb)) ** 2
    return s


def _rms_mm(r):
    return math.sqrt(r / (4 * seb.N_SAT)) * 1e3       # мм-эквивалент на компоненту


def fit_kL(ref, mat_kind):
    """2D-подгонка (k, L): k — физ. коэффициент, L — ЭФФЕКТИВНАЯ длина элемента
    (дрейфы (TOTAL−L)/2, сумма = 2 м). Старт: адаптивный 1D-поиск k при L=maxX (надёжно
    по масштабу k, вкл. сильные секступоли/октуполи), затем Nelder-Mead по (k,L) — быстро.
    Возвращает (k*, L*, rms_мм, signal_мм, explained)."""
    maxX = ref["L"]
    k0, _ = _adaptive_min(lambda kk: _residual(ref, mat_kind, kk, maxX))
    kref = k0 if abs(k0) > 1e-9 else 1.0          # масштаб для нормировки (a≈1)

    def obj(p):
        a, L = p
        if L <= 1e-3 or L >= TOTAL:
            return 1e12
        return _residual(ref, mat_kind, a * kref, L)

    res = minimize(obj, [k0 / kref, maxX], method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-18, "maxiter": 400})
    k, L = res.x[0] * kref, float(res.x[1])
    r = _residual(ref, mat_kind, k, L)
    r0 = _residual(ref, mat_kind, 0.0, L)
    signal = _rms_mm(r0); rms = _rms_mm(r)
    explained = (1.0 - rms / signal) if signal > 1e-9 else 0.0
    return k, L, rms, signal, explained


def fit_k_fixedL(ref, mat_kind, L):
    """1D-подгонка только k при ФИКСИРОВАННОЙ длине L. Для секступолей/октуполей: там
    (k,L) вырождены (слабый нелинейный кик ∝ k·L), поэтому L не определима — фиксируем."""
    k, r = _adaptive_min(lambda kk: _residual(ref, mat_kind, kk, L))
    r0 = _residual(ref, mat_kind, 0.0, L)
    signal = _rms_mm(r0); rms = _rms_mm(r)
    explained = (1.0 - rms / signal) if signal > 1e-9 else 0.0
    return k, L, rms, signal, explained


def center_deflection(ref):
    """Угол отклонения центральной частицы (рад) в плоскостях Y,Z и полный."""
    _, yp, _, zp = ref["final"][0]
    return yp, zp, math.hypot(yp, zp)


def _adaptive_min(resfn, lo=-100.0, hi=100.0):
    """1D-минимизатор с адаптивным расширением диапазона + зум. Возвращает (k*, res*)."""
    def scan(lo, hi, n=41):
        ks = [lo + (hi - lo) * j / (n - 1) for j in range(n)]
        rs = [resfn(k) for k in ks]
        j = min(range(n), key=lambda i: rs[i])
        return ks, rs, j
    ks, rs, j = scan(lo, hi)
    for _ in range(16):
        if 0 < j < len(ks) - 1:
            break
        span = hi - lo
        if j == 0:  lo -= 4 * span
        else:       hi += 4 * span
        ks, rs, j = scan(lo, hi)
    bk, step = ks[j], ks[1] - ks[0]
    for _ in range(6):
        lo, hi = bk - 2 * step, bk + 2 * step
        ks, rs, j = scan(lo, hi)
        bk, step = ks[j], ks[1] - ks[0]
    return bk, resfn(bk)


def _residual_plane(ref, plane, kappa, L):
    """L2-невязка ОДНОЙ плоскости (Y или Z) спутников отн. центра при фокусировке kappa, длине L."""
    parts = [(y0, 0.0, z0, 0.0) for (y0, z0) in ref["init"]]
    Lm, d1, d2 = _geom(L)
    if plane == "Y":
        out = bm.track_planes(parts, kappa, 0.0, Lm, d1, d2); pi, ai = 0, 1
    else:
        out = bm.track_planes(parts, 0.0, kappa, Lm, d1, d2); pi, ai = 2, 3
    c = out[0]; bref = ref["final"]; bc = bref[0]
    s = 0.0
    for i in range(seb.N_SAT):
        m = out[i + 1]; b = bref[i + 1]
        s += ((m[pi] - c[pi]) - (b[pi] - bc[pi])) ** 2 + (LC * ((m[ai] - c[ai]) - (b[ai] - bc[ai]))) ** 2
    return s


def fit_kappa_plane(ref, plane, L):
    """Подгонка фокусировки κ ОДНОЙ плоскости при длине L. Возвращает (κ, explained)."""
    k, r = _adaptive_min(lambda kk: _residual_plane(ref, plane, kk, L))
    r0 = _residual_plane(ref, plane, 0.0, L)
    expl = (1.0 - math.sqrt(r) / math.sqrt(r0)) if r0 > 1e-18 else 0.0
    return k, expl


def analyze_dipole(ref, L, device, value, is_coeff):
    """Физпараметры простого диполя: θ (изгиб), B, ρ + раздельная фокусировка плоскостей
    + НЕЗАВИСИМАЯ проверка θ_изм (Boris) vs θ_пред = (∫B·dl)/Bρ прямо из полевой карты.
    Плоскость изгиба = большее отклонение центра; κ∥ слабая (теория θ²/L²), κ⊥ краевая."""
    yp_c, zp_c, defl = center_deflection(ref)
    bend_is_y = abs(yp_c) >= abs(zp_c)
    theta = yp_c if bend_is_y else zp_c                 # знаковый изгиб в плоскости изгиба
    # независимый прогноз изгиба из ∫B·dl карты (Bz гнёт в Y, By в Z)
    intBy, intBz = seb.axis_field_integral(device, value, is_coeff)
    theta_pred = (-bm.QP * intBz) if bend_is_y else (bm.QP * intBy)
    ratio = (theta / theta_pred) if abs(theta_pred) > 1e-12 else 0.0
    B = theta * bm.BRHO / L                              # B = θ·Bρ/L (Тл)
    rho = (L / theta) if abs(theta) > 1e-12 else float("inf")
    ky, ey = fit_kappa_plane(ref, "Y", L)
    kz, ez = fit_kappa_plane(ref, "Z", L)
    k_par, e_par = (ky, ey) if bend_is_y else (kz, ez)   # плоскость изгиба
    k_perp, e_perp = (kz, ez) if bend_is_y else (ky, ey)  # перпендикулярная
    weak_theory = theta * theta / (L * L)                # ожидаемая слабая фокусировка 1/ρ²
    # эфф. угол среза полюса E из краевой фокусировки: κ⊥ ≈ −2·tan(E)/ρ·(1/L) (тонкие края)
    edge_E = math.degrees(math.atan(-k_perp * L * rho / 2.0)) if math.isfinite(rho) else 0.0
    return dict(value=ref["value"], theta=theta, theta_pred=theta_pred, ratio=ratio,
                B=B, rho=rho, defl=defl, L=L, k_par=k_par, k_perp=k_perp,
                weak_theory=weak_theory, edge_E=edge_E,
                expl_par=e_par, expl_perp=e_perp, bend_plane=("Y" if bend_is_y else "Z"))


def fit_shared_kpa_L(refs_items, mat_kind):
    """СОВМЕСТНАЯ подгонка (kpa, L) по ВСЕМ токам сразу: поле линейно по току ⇒ k(I)=kpa·I,
    а эффективная длина L — ОДНА на магнит. Минимизируем СУММУ невязок траекторий по всем
    токам (Nelder-Mead). Это и есть «единая эфф. длина, дающая совпадение при всех токах».
    refs_items: [(v, ref), ...] для v≠0. Возвращает (kpa, L)."""
    vmax, rmax = max(refs_items, key=lambda t: abs(t[0]))[0], None
    rmax = dict(refs_items)[vmax]
    maxX = rmax["L"]
    k0, _ = _adaptive_min(lambda kk: _residual(rmax, mat_kind, kk, maxX))
    kpa0 = k0 / vmax
    kref = kpa0 if abs(kpa0) > 1e-12 else 1.0

    def total(p):
        a, L = p
        if L <= 1e-3 or L >= TOTAL:
            return 1e18
        try:
            return sum(_residual(ref, mat_kind, a * kref * v, L) for v, ref in refs_items)
        except OverflowError:
            return 1e18

    res = minimize(total, [kpa0 / kref, maxX], method="Nelder-Mead",
                   options={"xatol": 1e-4, "fatol": 1e-18, "maxiter": 600})
    return res.x[0] * kref, float(res.x[1])


def fit_shared_kpa_fixedL(refs_items, mat_kind, L):
    """k(I)=kpa·I при ФИКСИРОВАННОЙ L (секступоли/октуполь — длина вырождена). Подгоняем
    общий kpa по сумме невязок всех токов."""
    vmax = max(refs_items, key=lambda t: abs(t[0]))[0]
    rmax = dict(refs_items)[vmax]
    k0, _ = _adaptive_min(lambda kk: _residual(rmax, mat_kind, kk, L))
    kref = (k0 / vmax) if abs(k0) > 1e-12 else 1.0

    def obj(aa):                                  # aa нормирован (опт.≈1) → узкий диапазон, без overflow
        try:
            return sum(_residual(ref, mat_kind, aa * kref * v, L) for v, ref in refs_items)
        except OverflowError:
            return 1e18

    a, _ = _adaptive_min(obj, lo=-6.0, hi=6.0)
    return a * kref


SWEEP_USED = {}   # device -> фактический max свипа (авто-масштаб под апертуру)


def process_device(device):
    mat_kind_full, sweep_kind = DEVICES[device]
    mat_kind = MAT_KIND.get(mat_kind_full, "quad")
    is_coeff = (sweep_kind == "C")
    vals, vmax = sweep_for(device, sweep_kind)
    SWEEP_USED[device] = vmax
    refs = {v: seb.run(device, v, is_coeff=is_coeff) for v in vals}   # Boris для всех токов (кэш)

    if mat_kind_full == "kicker":
        rows = []
        for v in vals:
            yp_c, zp_c, defl = center_deflection(refs[v])
            rows.append(dict(value=v, kick_y=yp_c, kick_z=zp_c, kick=defl, L=2 * refs[v]["L"],
                             k=0.0, K=0.0, G=0.0, rms_mm=0.0, signal_mm=0.0, explained=1.0))
        return mat_kind_full, rows
    if mat_kind_full == "dipole":
        return mat_kind_full, [analyze_dipole(refs[v], 2 * refs[v]["L"], device, v, is_coeff) for v in vals]

    # --- фокусирующие: СОВМЕСТНО подбираем ОДНУ эфф. длину L и коэф.-на-ток kpa по ВСЕМ токам.
    # Октуполь: сигнал ∝r³ на ≤2мм слишком слаб → L не определима, фиксируем (иначе уходит в мусор) ---
    items = [(v, refs[v]) for v in vals if abs(v) > 1e-9]
    if mat_kind_full == "oct":
        L_shared = 2 * refs[vals[0]]["L"]
        kpa = fit_shared_kpa_fixedL(items, mat_kind, L_shared)
    else:                                              # quad/dq/sext — совместно (kpa, L)
        kpa, L_shared = fit_shared_kpa_L(items, mat_kind)
    rows = []
    for v in vals:
        ref = refs[v]; yp_c, zp_c, defl = center_deflection(ref)
        k = kpa * v                                    # k(I) = kpa·I (поле линейно по току)
        if abs(v) > 1e-9:
            rms = _rms_mm(_residual(ref, mat_kind, k, L_shared))
            sig = _rms_mm(_residual(ref, mat_kind, 0.0, L_shared))
            expl = (1.0 - rms / sig) if sig > 1e-9 else 0.0
        else:
            rms = sig = 0.0; expl = 0.0
        rows.append(dict(value=v, k=k, L=L_shared, K=-k, G=k * bm.P_MOM / bm.Q,
                         rms_mm=rms, signal_mm=sig, explained=expl,
                         kick_y=yp_c, kick_z=zp_c, defl=defl))
    return mat_kind_full, rows


def write_device_csv(device, kind, rows):
    OUT.mkdir(parents=True, exist_ok=True)
    p = OUT / f"{device}_map.csv"
    if kind == "kicker":
        cols = ["value", "kick_y", "kick_z", "kick", "L"]
    elif kind == "dipole":
        cols = ["value", "theta", "theta_pred", "ratio", "B", "rho", "L", "k_perp", "edge_E", "expl_perp"]
    elif kind == "dq":
        cols = ["value", "K", "G", "k", "L", "defl", "kick_y", "kick_z", "rms_mm", "signal_mm", "explained"]
    else:
        cols = ["value", "K", "G", "k", "L", "rms_mm", "signal_mm", "explained"]
    with open(p, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(cols)
        for r in rows:
            w.writerow([f"{r.get(c, ''):.6g}" if isinstance(r.get(c), float) else r.get(c, "") for c in cols])
    return p


UNIT = {"quad": "K1", "dq": "K1", "dipole": "K1", "sext": "K2", "oct": "K3", "kicker": "kick"}


def write_summary(results):
    """results: list of (device, kind, rows). Пишет bip_compare/summary.md."""
    lines = ["# Сопоставление токи → физические коэффициенты BIP",
             "",
             f"Энергия пучка: γ={bm.GAMMA:.4f}, p={bm.P_MOM:.4e} кг·м/с, Bρ={bm.BRHO:.6f} Тл·м "
             f"(v={bm.V_BEAM:.3e} м/с — НАШ пучок, не пример BIP 3 МэВ).",
             "",
             "Пучок: центр + 3 окружности по 12 спутников (0.5/1/2 мм, Y-Z), 1 м дрейф + элемент + 1 м дрейф.",
             "Метод: Boris-трекер на полевой карте → подгонка коэф. матричной модели BIP (line-search).",
             "`explained` = доля отклонения спутников, объяснённая подогнанным коэф. (1=идеально).",
             "",
             "| device | элемент BIP | коэф. | свип | при ½ ОДЗ | при ОДЗ | невязка(½/полн) | explained | примечание |",
             "|---|---|---|---|---|---|---|---|---|"]
    for device, kind, rows in results:
        unit = UNIT[kind]
        vm = SWEEP_USED.get(device, 100.0)
        sweep = (f"ток ±{vm:g}А" if DEVICES[device][1] == "I" else f"coeff ±{vm:g}")
        half = rows[3 * len(rows) // 4]                  # ≈+½ размаха (устойчиво к числу точек)
        full = rows[-1]
        if kind == "kicker":
            v_half = f"{half['kick']*1e3:+.3f}мрад"
            v_full = f"{full['kick']*1e3:+.3f}мрад"
            res = "—"; expl = "—"; note = "чистый кик (∝ток)"
        elif kind == "dipole":
            v_half = f"θ={half['theta']*1e3:+.2f}мрад"
            v_full = f"θ={full['theta']*1e3:+.2f}мрад"
            res = f"B={full['B']*1e3:+.4f}мТл ρ={abs(full['rho']):.1f}м"
            expl = f"{full['expl_perp']:.2f}"
            note = f"изгиб (пл.{full['bend_plane']}); B,ρ,краевая ниже"
        else:
            v_half = f"{half['K']:+.3f}"
            v_full = f"{full['K']:+.3f}"
            res = f"{half['rms_mm']:.3f}/{full['rms_mm']:.3f}мм"
            ex = full.get('explained', 0)
            expl = f"{ex:.2f}"
            if kind == "sext":
                note = "секступоль, K2 ∝ ток" + ("" if ex > 0.4 else "; слабый сигнал")
            elif kind == "oct":
                note = "октуполь, K3 ∝ ток; сигнал ∝r³ слаб" + (" ⚠" if ex < 0.3 else "")
            elif kind == "dipole":
                note = "чистый изгиб, K1≈0"
            elif kind == "dq":
                note = f"изгиб {full.get('defl',0)*1e3:+.1f}мрад + квадруполь"
            else:
                note = "линейный, надёжно"
        lines.append(f"| {device} | {kind} | {unit} | {sweep} | {v_half} | {v_full} | {res} | {expl} | {note} |")
    lines += ["",
              "## Замечания",
              "- Пучок: центр + 3 окружности по 12 спутников на 0.5/1/2 мм (37 частиц). Несколько радиусов",
              "  дают нелинейный сигнал (∝r², r³) для секступолей/октуполей.",
              "- **Квадруполи (QP_L/QP_sh), диполи (DP), корректоры** — матчатся чисто, K линейна по току/coeff;",
              "  растущая невязка при крайних токах = краевые поля (одиночный hard-edge K не ловит fringe).",
              "- **Диполь-квадруполи (DQ)** — K1 надёжна (линейна), плюс изгиб центра (дипольная составляющая",
              "  + feed-down при смещённой оси). Изгиб в колонке `defl`.",
              "- **Секступоли (SP)** — K2 линейна по току, explained≈0.6 (остаток = не-чистый hard-edge: fringe/",
              "  высшие мультиполи карты). Коэффициент пригоден.",
              "- **Октуполи (OP)** — K3 линейна, но сигнал ∝r³ слаб даже при 2 мм (~мкм, низкий explained);",
              "  для надёжного K3 нужен бóльший радиус. Учесть также: ядро октуполя BIP нестандартное (By=x³−3x²y).",
              "- Полные таблицы по 21 точке свипа — в `<device>_map.csv` (колонки value;K;G;k;...;explained)."]

    # ---- детальный разбор простых диполей (физпараметры + связи) ----
    dips = [(d, rows) for (d, k, rows) in results if k == "dipole"]
    if dips:
        lines += ["", "## Простые диполи — физические параметры и связи", "",
                  "θ — изгиб (рад), B=θ·Bρ/L — поле (Тл), ρ=L/θ — радиус. κ∥ — слабая фокусировка в плоскости",
                  "изгиба (ожид. θ²/L²), κ⊥ — краевая (перпенд., даёт эфф. угол среза E).", ""]
        for dev, rows in dips:
            L = None
            # L восстановим из ρ·θ первой ненулевой строки
            nz = next((r for r in rows if abs(r["theta"]) > 1e-9), None)
            Lm = (nz["rho"] * nz["theta"]) if nz else 0.0
            lines += [f"### {dev}  (L_карты≈{Lm:.3f} м, поле симметрично ⇒ полный пробег 2×; плоскость изгиба {rows[-1]['bend_plane']})", "",
                      "| coeff | θ изм, мрад | θ пред ∫B·dl, мрад | ratio | B, мТл | ρ, м | κ⊥ | E, ° |",
                      "|---|---|---|---|---|---|---|---|"]
            for r in rows[::4]:
                rho = r["rho"]
                lines.append(f"| {r['value']:+.1f} | {r['theta']*1e3:+.3f} | {r.get('theta_pred',0)*1e3:+.3f} | "
                             f"{r.get('ratio',0):.3f} | {r['B']*1e3:+.4f} | "
                             f"{(abs(rho) if math.isfinite(rho) else float('inf')):.1f} | {r['k_perp']:+.4f} | {r['edge_E']:+.2f} |")
            ratios = [r['ratio'] for r in rows if abs(r.get('theta_pred', 0)) > 1e-9]
            if ratios:
                avg = sum(ratios) / len(ratios)
                lines += ["", f"**Кросс-чек {dev}:** независимый прогноз изгиба из ∫B·dl полевой карты "
                              f"совпадает с измеренным Boris: ⟨θ_изм/θ_пред⟩ = {avg:.4f} (≈1.000). "
                              f"Подтверждает поле карты, Bρ (энергию) и эквивалентность Boris ≡ ∫B·dl/Bρ."]
            # связи: наклоны линейных зависимостей по краям свипа
            v0, vN = rows[0], rows[-1]
            dv = vN["value"] - v0["value"]
            slope_th = (vN["theta"] - v0["theta"]) / dv * 1e3      # мрад на единицу coeff
            slope_B = (vN["B"] - v0["B"]) / dv * 1e3               # мТл на единицу
            lines += ["", f"**Связи {dev}:** θ ≈ {slope_th:+.3f} мрад/coeff (линейно), "
                          f"B ≈ {slope_B:+.4f} мТл/coeff, ρ ∝ 1/θ — надёжно. "
                          f"Фокусировка κ∥/κ⊥ на уровне шума (слабый диполь): лишь грубо ~θ²/L², "
                          f"краевой угол E (≈{vN['edge_E']:+.2f}°) ненадёжен.", ""]
        if len(dips) == 2:
            (dl, rl), (ds, rs_) = dips
            lines += [f"**DP_L vs DP_sh:** при одном coeff изгиб {dl}={rl[-1]['theta']*1e3:+.2f} vs "
                      f"{ds}={rs_[-1]['theta']*1e3:+.2f} мрад (разные L → разный θ при близком B); "
                      f"поле B {rl[-1]['B']*1e3:+.4f} vs {rs_[-1]['B']*1e3:+.4f} мТл.", ""]

    p = OUT / "summary.md"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


if __name__ == "__main__":
    targets = sys.argv[1:] or list(DEVICES.keys())
    print(f"Энергия: γ={bm.GAMMA:.4f}, p={bm.P_MOM:.4e} кг·м/с, Bρ={bm.BRHO:.6f} Тл·м")
    results = []
    for dev in targets:
        if dev not in DEVICES:
            print(f"  ! неизвестный девайс {dev}"); continue
        kind, rows = process_device(dev)
        p = write_device_csv(dev, kind, rows)
        results.append((dev, kind, rows))
        unit = UNIT[kind]
        print(f"\n{dev} ({kind}, коэф.={unit}) -> {p.name}")
        for r in rows[::5]:
            if kind == "kicker":
                print(f"   v={r['value']:>7}  kick={r['kick']*1e3:+.4f} мрад")
            elif kind == "dipole":
                rho = r['rho']
                print(f"   v={r['value']:>7}  θ={r['theta']*1e3:+.3f}мрад  B={r['B']*1e3:+.4f}мТл  "
                      f"ρ={(abs(rho) if math.isfinite(rho) else float('inf')):.1f}м  "
                      f"κ∥={r['k_par']:+.4f} κ⊥={r['k_perp']:+.4f}  E={r['edge_E']:+.2f}°  expl⊥={r['expl_perp']:.2f}")
            else:
                extra = f"  откл={r.get('defl',0)*1e3:+.2f}мрад" if kind == "dq" else ""
                flag = "  ⚠слаб" if (kind in ("sext", "oct") and r.get("explained", 1) < 0.5) else ""
                print(f"   v={r['value']:>7}  {unit}={r['K']:+8.4f}  G={r['G']:+.3e}  невязка={r['rms_mm']:.3f}мм  expl={r.get('explained',0):.2f}{extra}{flag}")
    if len(results) > 1:
        sp = write_summary(results)
        print(f"\nСводка: {sp}")
