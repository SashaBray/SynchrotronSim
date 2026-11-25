# -*- coding: utf-8 -*-
"""Продольная динамика в ВЧ-резонаторе: синхронная фаза, сепаратриса, границы по импульсу и фазе, частота
синхротронных колебаний — для окна «Резонатор ВЧ» кокпита и анализа (чистые функции, без Qt).

Соглашения — как в ядре SynchrotronSim (src/matrix_ring.h, RfSpec): резонатор в конце пролёта; фаза частицы на
нём φ = φ_прог + 2πh·t/C_RF (t — продольная координата частицы, м; C_RF = h·c/f — длина ВЧ-оборота); прибавка
энергии ΔE = V·sin φ (МэВ при V в МВ). За оборот (ультрарелятивистски, ΔE/E = δ):
    φ ← φ + 2πh·η·δ           (η = αc − 1/γ², выше критической энергии η > 0)
    δ ← δ + (V·sin φ − U₀)/E  (U₀ — потеря на излучение за оборот)
Гамильтониан H(φ, δ) = πhη·δ² + (V/E)·(cos φ + φ·sin φs); устойчивая точка φs: sin φs = U₀/V, при η > 0 —
cos φs < 0, т.е. φs = π − asin(U₀/V). Сепаратриса проходит через неустойчивую точку φu = π − φs.
Затухание (излучение ∝ E⁴) в сепаратрису не входит: за оборот оно ~1/τ ≈ 4·10⁻⁴ (τ ≈ 2650 оборотов).
"""
from __future__ import annotations

import math
from typing import Any, Dict, Optional

import numpy as np

MC2_MEV = 0.51099895069
C_LIGHT = 299792458.0


def _U(phi, sphis):
    """Потенциал по фазе (без множителя V/E): cos φ + φ·sin φs."""
    return np.cos(phi) + phi * sphis


def bucket(V_mv: float, U0_mev: float, E_mev: float, h: float, alphac: float, phase_deg: float,
           freq_hz: float, f0_hz: float) -> Dict[str, Any]:
    """Параметры ВЧ-ямы (bucket). Возврат — словарь (числа; фазы в градусах, δ — доли):
    ok — яма есть (V > U₀ и η ≠ 0); phis — синхронная фаза; phiu, phi2 — границы по фазе (неустойчивая точка и
    второй край сепаратрисы); delta_max — полувысота ямы (граница по импульсу от центра); delta_s — энергия
    синхронной частицы от сдвига частоты Δf = f − f0 (δs ≈ −(Δf/f)/η); nus — частота синхротронных колебаний,
    на оборот; t_s — где по t стоит синхронная частица (м); eta, gamma."""
    gamma = E_mev / MC2_MEV
    eta = alphac - 1.0 / gamma ** 2
    out = {"V_mv": V_mv, "U0_mev": U0_mev, "E_mev": E_mev, "h": h, "alphac": alphac, "eta": eta, "gamma": gamma,
           "phase_deg": phase_deg, "freq_hz": freq_hz, "f0_hz": f0_hz, "df_hz": freq_hz - f0_hz,
           "c_rf": h * C_LIGHT / freq_hz if freq_hz > 0 else float("nan")}
    out["delta_s"] = -((freq_hz - f0_hz) / freq_hz) / eta if (eta and freq_hz > 0) else 0.0
    if not (V_mv > abs(U0_mev)) or eta == 0.0:
        out.update(ok=False, reason=("напряжения не хватает на потерю за оборот: V ≤ U₀" if not (V_mv > abs(U0_mev))
                                     else "η = 0 (критическая энергия)"))
        return out
    s = U0_mev / V_mv
    phis = (math.pi - math.asin(s)) if eta > 0 else math.asin(s)
    phiu = math.pi - phis
    Uu = float(_U(phiu, s))
    k = V_mv / E_mev / (math.pi * h * eta)                  # δ² = k·(U(φu) − U(φ)) на сепаратрисе
    dmax2 = k * (Uu - float(_U(phis, s)))
    # второй край ямы по фазе: U(φ2) = U(φu) по другую сторону от φs
    lo, hi = (phis, phis + 2.0 * math.pi) if phiu < phis else (phis - 2.0 * math.pi, phis)
    g = lambda p: float(_U(p, s)) - Uu
    a, b = (phis, phis + 0.5) if phiu < phis else (phis - 0.5, phis)
    step = 0.05 if phiu < phis else -0.05
    p = phis
    while lo - 1e-9 <= p <= hi + 1e-9 and g(p) < 0:
        p += step
    a, b = (p - step, p) if step > 0 else (p, p - step)
    for _ in range(60):                                     # бисекция
        m = 0.5 * (a + b)
        if (g(m) < 0) == (g(a) < 0):
            a = m
        else:
            b = m
    phi2 = 0.5 * (a + b)
    nus2 = -h * eta * V_mv * math.cos(phis) / (2.0 * math.pi * E_mev)
    t_s = (phis - math.radians(phase_deg)) * out["c_rf"] / (2.0 * math.pi * h)
    out.update(ok=True, phis_deg=math.degrees(phis), phiu_deg=math.degrees(phiu), phi2_deg=math.degrees(phi2),
               delta_max=math.sqrt(max(dmax2, 0.0)), nus=math.sqrt(nus2) if nus2 > 0 else float("nan"), t_s=t_s,
               k=k)
    return out


def separatrix(b: Dict[str, Any], n: int = 361):
    """Сепаратриса ямы: (φ, град; δ верхняя ветвь) от φu до φ2, вокруг δs; нижняя — δs − (δ − δs)."""
    if not b.get("ok"):
        return np.array([]), np.array([])
    s = b["U0_mev"] / b["V_mv"]
    p1, p2 = math.radians(b["phiu_deg"]), math.radians(b["phi2_deg"])
    phi = np.linspace(min(p1, p2), max(p1, p2), int(n))
    d2 = b["k"] * (float(_U(math.radians(b["phiu_deg"]), s)) - _U(phi, s))
    d = np.sqrt(np.clip(d2, 0.0, None))
    return np.degrees(phi), b["delta_s"] + d


def particle_phase_deg(t, phase_deg: float, h: float, c_rf: float):
    """Фаза ВЧ, на которую приходит частица с продольной координатой t (м): φ_прог + 360°·h·t/C_RF."""
    return phase_deg + 360.0 * h * np.asarray(t, float) / c_rf


def inside(b: Dict[str, Any], phi_deg, delta) -> np.ndarray:
    """Частицы внутри сепаратрисы (H < H_sep, по фазе — в пределах ямы)."""
    if not b.get("ok"):
        return np.zeros(np.shape(phi_deg), bool)
    s = b["U0_mev"] / b["V_mv"]
    phi = np.radians(np.asarray(phi_deg, float))
    lo, hi = sorted((math.radians(b["phiu_deg"]), math.radians(b["phi2_deg"])))
    phi = lo + np.mod(phi - lo, 2.0 * math.pi)              # к периоду, начинающемуся у края ямы
    dd = np.asarray(delta, float) - b["delta_s"]
    return (phi <= hi) & (dd ** 2 < b["k"] * (float(_U(math.radians(b["phiu_deg"]), s)) - _U(phi, s)))


# ------------------------------------------------------------------ все параметры системы «резонатор — пучок»
C_GAMMA = 8.846273435551005e-5      # м/ГэВ³ — потеря энергии (как ring_matrix)
C_Q = 3.8319402812e-13              # м — квантовое возбуждение
E_CHARGE = 1.602176634e-19          # Кл


def derived(b: Dict[str, Any], circ_m: float, rad: Optional[Dict[str, Any]] = None,
            current_a: Optional[float] = None, n_bunches: Optional[int] = None) -> Dict[str, Any]:
    """Дополнить словарь ямы (`bucket`, + пучок модели) производными параметрами для окна «Резонатор ВЧ».

    circ_m — периметр (длина проектной орбиты); rad — излучение текущей машины из модели установки
    (`TbtPlant.radiation_now`: I2, I3, I2_id, n_id, damping.s — продольная мода матрицы оборота 6D);
    current_a, n_bunches — ток и число сгустков заполнения (в данных проекта их нет — допущения вызывающего).
    Добавляет ключи (единицы СИ, фазы в градусах, δ — доли); чего не из чего посчитать — не добавляет."""
    out = b
    E, h, f = float(b["E_mev"]), float(b["h"]), float(b["freq_hz"])
    gamma = float(b["gamma"]); beta = math.sqrt(max(0.0, 1.0 - 1.0 / gamma ** 2))
    out["circ_m"] = float(circ_m)
    out["f_rev_hz"] = beta * C_LIGHT / circ_m
    out["T0_s"] = 1.0 / out["f_rev_hz"]
    out["n_buckets"] = h
    out["lambda_rf_m"] = float(b["c_rf"]) / h
    out["T_rf_s"] = 1.0 / f if f > 0 else float("nan")
    ac = float(b["alphac"])
    if ac > 0:
        out["E_t_mev"] = MC2_MEV / math.sqrt(ac)
    U0 = float(b["U0_mev"])
    if U0 > 0:
        out["q_over"] = float(b["V_mv"]) / U0
    if b.get("ok"):
        out["width_deg"] = abs(float(b["phi2_deg"]) - float(b["phiu_deg"]))
        sx, su = separatrix(b, 721)
        area = float(np.trapezoid(2.0 * (su - float(b["delta_s"])), np.radians(sx)))   # рад × δ
        out["area_evs"] = area * E * 1e6 / (2.0 * math.pi * f)                       # (t, ΔE): эВ·с
        nus = float(b["nus"])
        if nus == nus and nus > 0:
            out["fs_hz"] = nus * out["f_rev_hz"]
            out["Ts_turns"] = 1.0 / nus
    eta = float(b["eta"])
    if rad:
        E_gev = E / 1e3
        I2, I3 = float(rad.get("I2", 0.0)), float(rad.get("I3", 0.0))
        out["I2"], out["I3"] = I2, I3
        out["U0_i2_mev"] = C_GAMMA * E_gev ** 4 * I2 / (2.0 * math.pi) * 1e3
        out["U0_id_mev"] = C_GAMMA * E_gev ** 4 * float(rad.get("I2_id", 0.0)) / (2.0 * math.pi) * 1e3
        if rad.get("n_id") is not None:
            out["n_id"] = int(rad["n_id"])
        ds = (rad.get("damping") or {}).get("s")          # продольная мода матрицы оборота 6D
        if ds:
            tau = float(ds["tau_turns"])
            if U0 > 0 and math.isfinite(tau) and tau > 0:
                out["tau_s_turns"] = tau
                out["Js"] = 2.0 * E / (U0 * tau)
                # плоское кольцо (нет вертикальных изгибов): Jy = 1, теорема Робинсона Jx + Jy + Js = 4
                out["Jy"] = 1.0
                out["Jx"] = 3.0 - out["Js"]
                out["tau_y_turns"] = 2.0 * E / U0
                if out["Jx"] > 0:
                    out["tau_x_turns"] = 2.0 * E / (U0 * out["Jx"])
            nu = float(ds["tune"])
            if nu > 0:
                out["nus_m6"] = nu
        Js = out.get("Js")
        if Js and Js > 0 and I2 > 0:
            out["sdelta0"] = math.sqrt(C_Q * gamma ** 2 * I3 / (Js * I2))
            nus = out.get("nus_m6") or b.get("nus")
            if nus and nus == nus and nus > 0:
                out["sz0_m"] = circ_m * abs(eta) * out["sdelta0"] / (2.0 * math.pi * float(nus))
    if current_a is not None:
        nb = int(n_bunches or h)
        out["I_a"], out["n_bunches"] = float(current_a), nb
        out["q_bunch_c"] = float(current_a) * out["T0_s"] / nb
        out["ne_bunch"] = out["q_bunch_c"] / E_CHARGE
    if b.get("n_part") and b.get("ok") and b.get("phi_mean") is not None:
        out["dphi_deg"] = float(b["phi_mean"]) - float(b["phis_deg"])
        nus = float(b["nus"])
        if nus == nus and nus > 0:                           # малые колебания: δ-амплитуда = A·νs/(hη)
            dd = (float(b["delta_mean"]) - float(b["delta_s"])) * h * eta / nus
            out["amp_deg"] = math.hypot(out["dphi_deg"], math.degrees(dd))
    return out


def _ms(v, b):
    """Обороты → миллисекунды."""
    return f" · {v * b['T0_s'] * 1e3:.3f} мс" if b.get("T0_s") else ""


def _phase_time(v, b):
    """Градусы ВЧ-фазы → время и длина."""
    if not b.get("freq_hz"):
        return ""
    t = abs(v) / (360.0 * float(b["freq_hz"]))
    return f" · {t * 1e12:.1f} пс · {t * C_LIGHT * 1e3:.2f} мм"


def _n_id(v, b):
    return f" · {b['n_id']} ИД по кольцу" if b.get("n_id") is not None else ""


# каталог окна «Резонатор ВЧ»: группа → [(ключ, название, единицы, множитель, формат, откуда/формула, доп. текст)]
PARAMS = [
    ("исходные (берутся в расчёт)", [
        ("E_mev", "энергия пучка E", "МэВ", 1.0, "{:.1f}", "паспорт машины", None),
        ("gamma", "лоренц-фактор γ", "", 1.0, "{:.1f}", "E/mc²", None),
        ("circ_m", "периметр C", "м", 1.0, "{:.4f}", "длина проектной орбиты кольца", None),
        ("f_rev_hz", "частота обращения f₀", "кГц", 1e-3, "{:.4f}", "βc/C", None),
        ("T0_s", "период обращения T₀", "мкс", 1e6, "{:.5f}", "1/f₀", None),
        ("h", "кратность ВЧ h", "", 1.0, "{:.0f}", "паспорт машины", None),
        ("f0_hz", "номинальная частота ВЧ", "МГц", 1e-6, "{:.6f}", "h·c/C", None),
        ("loop_df_hz", "сдвиг частоты от контура ВЧ", "Гц", 1.0, "{:+.3f}", "контур rf (энергия δ̂)", None),
        ("manual_df_hz", "ручной сдвиг Δf", "Гц", 1.0, "{:+.2f}", "ручка окна", None),
        ("freq_hz", "частота ВЧ f", "МГц", 1e-6, "{:.6f}", "номинал + контур + ручной (+ΔC кольца)", None),
        ("V_mv", "амплитуда V", "МВ", 1.0, "{:.4g}", "ручка окна", None),
        ("phase_deg", "программная фаза φ_прог", "°", 1.0, "{:.3f}", "ручка окна", None),
        ("alphac", "коэффициент уплотнения орбиты αc", "", 1.0, "{:.4e}", "линейная оптика кольца", None),
    ]),
    ("излучение и затухание (текущая машина)", [
        ("U0_mev", "потеря энергии за оборот U₀", "кэВ", 1e3, "{:.2f}",
         "энергобаланс оборота у 6D-орбиты: V·sin φ − ΔE", None),
        ("U0_i2_mev", "U₀ по изгибам (сверка)", "кэВ", 1e3, "{:.2f}", "Cγ·E⁴·I₂/2π, I₂ — изгибы при текущих зазорах",
         None),
        ("U0_id_mev", "в т.ч. ондуляторы", "кэВ", 1e3, "{:.2f}", "полюса ИД при текущих зазорах", _n_id),
        ("tau_s_turns", "время затухания τs (продольное)", "об", 1.0, "{:.0f}",
         "|λ| продольной моды матрицы оборота 6D", _ms),
        ("Js", "число затухания Js", "", 1.0, "{:.4f}", "2E/(U₀·τs)", None),
        ("Jx", "число затухания Jx", "", 1.0, "{:.4f}", "Робинсон: 4 − Jy − Js (плоское кольцо, Jy = 1)", None),
        ("Jy", "число затухания Jy", "", 1.0, "{:.4f}", "плоское кольцо: вертикальных изгибов нет", None),
        ("tau_x_turns", "время затухания τx", "об", 1.0, "{:.0f}", "2E/(Jx·U₀)", _ms),
        ("tau_y_turns", "время затухания τy", "об", 1.0, "{:.0f}", "2E/(Jy·U₀)", _ms),
        ("sdelta0", "естественный энергоразброс σδ₀", "10⁻³", 1e3, "{:.4f}", "√(Cq·γ²·I₃/(Js·I₂))", None),
        ("sz0_m", "естественная длина сгустка σz₀", "мм", 1e3, "{:.3f}", "C·|η|·σδ₀/(2π·νs)",
         lambda v, b: f" · {v / C_LIGHT * 1e12:.2f} пс"),
    ]),
    ("яма (сепаратриса)", [
        ("eta", "коэффициент скольжения η", "", 1.0, "{:.4e}", "αc − 1/γ²", None),
        ("E_t_mev", "критическая энергия E_t", "МэВ", 1.0, "{:.1f}", "mc²/√αc", None),
        ("q_over", "перенапряжение q", "", 1.0, "{:.4f}", "V/U₀", None),
        ("phis_deg", "синхронная фаза φs", "°", 1.0, "{:.3f}", "sin φs = U₀/V; η > 0 ⇒ φs = 180° − asin(U₀/V)", None),
        ("delta_s", "энергия синхронной частицы δs", "10⁻³", 1e3, "{:+.4f}", "−(Δf/f)/η", None),
        ("t_s", "запаздывание синхронной частицы t_s", "мм", 1e3, "{:+.4f}", "(φs − φ_прог)·C_RF/(2πh)", None),
        ("n_buckets", "сепаратрис (мест для сгустков) по кольцу", "", 1.0, "{:.0f}", "= h", None),
        ("lambda_rf_m", "шаг между сепаратрисами λ_RF", "м", 1.0, "{:.5f}", "C_RF/h = c/f",
         lambda v, b: f" · {v / C_LIGHT * 1e9:.4f} нс"),
        ("phiu_deg", "край ямы φu (неустойчивая точка)", "°", 1.0, "{:.3f}", "180° − φs", None),
        ("phi2_deg", "второй край ямы φ2", "°", 1.0, "{:.3f}", "U(φ2) = U(φu)", None),
        ("width_deg", "ширина ямы по фазе", "°", 1.0, "{:.3f}", "|φ2 − φu|", _phase_time),
        ("delta_max", "высота ямы δmax (аксептанс ВЧ)", "%", 100.0, "{:.4f}", "√(V/(E·πhη)·(U(φu) − U(φs)))",
         lambda v, b: f" · {v * float(b['E_mev']):.2f} МэВ"),
        ("area_evs", "площадь ямы", "эВ·с", 1.0, "{:.4e}", "∮ ΔE dt по сепаратрисе", None),
    ]),
    ("синхротронное движение", [
        ("nus", "синхротронная частота νs (формула)", "", 1.0, "{:.5e}", "√(−hηV·cos φs/(2πE)), малые колебания", None),
        ("nus_m6", "νs по матрице оборота 6D", "", 1.0, "{:.5e}", "arg λ/2π продольной моды", None),
        ("fs_hz", "синхротронная частота fs", "кГц", 1e-3, "{:.4f}", "νs·f₀", None),
        ("Ts_turns", "период синхротронных колебаний", "об", 1.0, "{:.1f}", "1/νs", _ms),
    ]),
    ("заполнение и ток (допущения)", [
        ("I_a", "ток пучка I", "мА", 1e3, "{:.1f}", "допущение (undulator_sr.P.I)", None),
        ("n_bunches", "сгустков в заполнении", "", 1.0, "{:.0f}", "допущение: заняты все h мест", None),
        ("q_bunch_c", "заряд сгустка", "нКл", 1e9, "{:.4f}", "I·T₀/n_b", None),
        ("ne_bunch", "электронов в сгустке", "", 1.0, "{:.3e}", "заряд / e", None),
    ]),
    ("пучок модели (один сгусток, макрочастицы)", [
        ("n_part", "макрочастиц", "", 1.0, "{:.0f}", "опыт кокпита", None),
        ("n_inside", "внутри ямы", "", 1.0, "{:.0f}", "H < H сепаратрисы", None),
        ("n_lost", "потеряно", "", 1.0, "{:.0f}", "вылетели из апертуры (NaN)", None),
        ("dphi_deg", "⟨φ⟩ − φs", "°", 1.0, "{:+.3f}", "центр сгустка", _phase_time),
        ("phi_std", "σφ (длина сгустка)", "°", 1.0, "{:.4f}", "разброс фаз частиц", _phase_time),
        ("delta_mean", "⟨δ⟩", "10⁻³", 1e3, "{:+.4f}", "центр сгустка", None),
        ("delta_std", "σδ", "10⁻³", 1e3, "{:.4f}", "разброс энергий частиц", None),
        ("amp_deg", "амплитуда синхротронных колебаний центра", "°", 1.0, "{:.3f}", "√(Δφ² + (Δδ·hη/νs)²), малые",
         _phase_time),
    ]),
]
PARAM_INDEX = {row[0]: (g, row[1], row[2], row[3]) for g, rows in PARAMS for row in rows}


def format_param(key: str, b: Dict[str, Any]) -> str:
    """Текст значения параметра для таблицы (число в единицах каталога + доп. перевод) или «—»."""
    for _g, rows in PARAMS:
        for k, _name, un, sc, fmt, _src, extra in rows:
            if k != key:
                continue
            v = b.get(k)
            if v is None or not isinstance(v, (int, float)) or v != v:
                return "—"
            s = fmt.format(float(v) * sc) + (f" {un}" if un else "")
            return s + (extra(float(v), b) if extra else "")
    return "—"
