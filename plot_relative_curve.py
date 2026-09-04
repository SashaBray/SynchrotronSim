# -*- coding: utf-8 -*-
"""Кривая градиента против кривой намагничивания: наложение и их отношение.

Замысел (предложен пользователем): подобрать ДВА масштаба — по горизонтали и по вертикали — так,
чтобы кривая «градиент от ампер-витков» легла на кривую намагничивания материала. Что останется
после деления одной на другую, и есть собственная характеристика МАГНИТА, очищенная от материала:

    G(NI)  ≈  β · B(α · NI)          α — растяжение по току, β — по величине
    отношение(NI) = G(NI) / (β·B(α·NI))

Если отношение определяется геометрией, а не материалом, то у двух серий карт одного магнита,
посчитанных на РАЗНЫХ материалах, эти кривые обязаны совпасть. У QP_sh такие серии есть: сталь 1008
(B_max 2.59 Тл) и низконасыщающийся материал (1.70 Тл, папка `qp_sh_mu_steps`, в реестре его карт
НЕТ — читаются напрямую из .fld).

Градиент берётся У ОСИ (радиус пробы траекторий), а не на 8 мм: при глубоком насыщении профиль
поперёк апертуры уплощается, и на разных радиусах градиент ведёт себя по-разному.

    py -3.14 plot_relative_curve.py
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm

CACHE = os.path.join("workspace", "results", "relative_curve_cache.json")
MU_DIR = os.path.join("C:" + os.sep, "Sasha", "Work", "SILA", "ansis approximations", "KD",
                      "QP_sh", "qp_sh_mu_steps")


def gradient_at_axis(path, radius, order=2, nth=64):
    """Градиент у оси: гармоника порядка n, проинтегрированная вдоль пучка, делённая на r^(n-1)."""
    from field_quality import load, _circle
    h, B = load(path)
    xs = np.linspace(h["mnx"], h["mxx"], h["nx"])
    v = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, radius, nth=nth)
        sp = np.fft.fft(by + 1j * bz) / nth
        v.append(sp[0] if order == 1 else sp[-(order - 1)])
    v = np.array(v)
    c = complex(2 * np.trapezoid(v.real, xs), 2 * np.trapezoid(v.imag, xs))
    return abs(c) / radius ** (order - 1)


def series(element="QP_sh", radius=1e-3):
    """Две серии одного магнита: (ампер-витки, градиент) на стали и на другом материале."""
    cache = {}
    if os.path.exists(CACHE):
        try:
            cache = json.load(open(CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    dirty = False
    steel = []
    for r in sm.read_registry():
        if r["element"] != element or ";" in r["excitation"] or "+" in r["excitation"]:
            continue
        p = os.path.join("workspace", "fields_bin", r["bin_path"])
        try:
            ni = abs(float(r["amp_turns"]))
        except (TypeError, ValueError):
            continue
        if ni <= 0 or not os.path.exists(p):
            continue
        key = f"steel|{r['map_id']}|{radius}"
        if key not in cache:
            cache[key] = gradient_at_axis(p, radius); dirty = True
        steel.append((ni, cache[key]))
    other = []
    for f in sorted(glob.glob(os.path.join(MU_DIR, f"{element}_*A.fld"))):
        m = re.search(r"_(\d+)A\.fld$", os.path.basename(f))
        if not m:
            continue
        key = f"mu|{os.path.basename(f)}|{radius}"
        if key not in cache:
            cache[key] = gradient_at_axis(f, radius); dirty = True
        other.append((int(m.group(1)), cache[key]))
    if dirty:
        os.makedirs(os.path.dirname(CACHE) or ".", exist_ok=True)
        json.dump(cache, open(CACHE, "w", encoding="utf-8"))
    return sorted(steel), sorted(other)


def bh(path):
    """Кривая намагничивания B(H), интерполированная ГЛАДКО и монотонно в лог-лог координатах.

    Линейная интерполяция здесь не годится: в `steel_1008_bh.csv` всего 19 точек, и её изломы
    прорастают прямо в кривую отношения — у стали получался провал до 0.80 и всплеск до 1.28,
    целиком от таблицы, а не от магнита. Монотонный кубический сплайн (PCHIP) в логарифмах эти
    изломы убирает и не создаёт ложных перегибов."""
    from scipy.interpolate import PchipInterpolator
    d = np.loadtxt(path, delimiter=";", skiprows=1)
    H, B = d[:, 0], d[:, 1]
    ok = (H > 0) & (B > 0)
    f = PchipInterpolator(np.log(H[ok]), np.log(B[ok]), extrapolate=False)
    mu0 = 4e-7 * np.pi
    Hmin, Hmax, Bmax = H[ok][0], H[-1], B[-1]

    def g(x):
        # ЗАЖИМ в ноль, а не складывание по модулю: при отрицательном сдвиге x₀ аргумент уходит
        # ниже нуля, и abs() отражал бы кривую, порождая нефизичный всплеск вне области данных.
        x = np.maximum(np.asarray(x, float), 0.0)
        out = np.empty_like(x)
        lo = x < Hmin
        hi = x > Hmax
        mid = ~(lo | hi)
        out[mid] = np.exp(f(np.log(np.maximum(x[mid], Hmin))))
        out[lo] = B[ok][0] / Hmin * x[lo]              # начальный участок линеен
        out[hi] = Bmax + (x[hi] - Hmax) * mu0          # за таблицей — по воздуху
        return out
    return g


def model(q, ni, Bf, shifts):
    """β·B(α·NI + x₀) + y₀. Без сдвигов x₀ = y₀ = 0."""
    a, b = np.exp(q[0]), np.exp(q[1])
    x0, y0 = (q[2], q[3]) if shifts else (0.0, 0.0)
    return b * Bf(a * ni + x0) + y0


def fit_scales(ni, g, Bf, shifts=False):
    """Масштабы (и, по желанию, СДВИГИ) наложения кривой материала на кривую градиента.

    Подгонка по ОТНОСИТЕЛЬНОЙ невязке — иначе слабые токи, где градиент на два порядка меньше,
    в подгонке не участвуют вовсе. Сдвиги физически осмысленны: x₀ похож на коэрцитивную силу,
    y₀ — на остаточную индукцию, хотя наша B(H) однозначна и гистерезиса не содержит."""
    def res(q):
        return (model(q, ni, Bf, shifts) - g) / g

    best = None
    for s in range(60):
        rng = np.random.default_rng(s)
        q0 = [np.log(rng.lognormal(np.log(1.0), 2.0)),
              np.log(rng.lognormal(np.log(g.max() / 1.5), 1.0))]
        if shifts:
            q0 += [rng.normal(0, 30.0), rng.normal(0, 0.1 * g.max())]
        try:
            r = least_squares(res, q0, xtol=1e-15, ftol=1e-15, max_nfev=40000)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    q = best.x
    a, b = np.exp(q[0]), np.exp(q[1])
    x0, y0 = (float(q[2]), float(q[3])) if shifts else (0.0, 0.0)
    return float(a), float(b), x0, y0, float(np.sqrt(np.mean(res(q) ** 2)) * 100), q


def main(argv=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", default="QP_sh")
    ap.add_argument("--radius", type=float, default=1e-3)
    ap.add_argument("--out", default=os.path.join("workspace", "results", "relative_curve.png"))
    a = ap.parse_args(argv)

    steel, other = series(a.element, a.radius)
    sets = [("сталь 1008", steel, "workspace/materials/steel_1008_bh.csv", "#1f4e79", "o"),
            ("низконасыщающийся", other, "workspace/materials/lowsat_custom_bh.csv", "#c0392b", "s")]

    fig, ax = plt.subplots(1, 2, figsize=(14.2, 5.6))
    print(f"{'серия':20} {'α':>11} {'β':>11} {'x0':>11} {'y0':>11} "
          f"{'со сдвигами':>10} {'без сдвигов':>14}")
    for name, pts, mat, col, mk in sets:
        if not pts:
            continue
        ni = np.array([p[0] for p in pts], float)
        g = np.array([p[1] for p in pts], float)
        Bf = bh(mat)
        a2, b2, _, _, e2, _ = fit_scales(ni, g, Bf, shifts=False)
        alpha, beta, x0, y0, err, q = fit_scales(ni, g, Bf, shifts=True)
        print(f"{name:20} {alpha:11.5g} {beta:11.5g} {x0:11.5g} {y0:11.5g} "
              f"{err:9.2f}% {e2:13.2f}%")

        # ── наложение: измеренный градиент и приведённая кривая материала
        ax[0].plot(ni, g, mk, color=col, ms=6, zorder=3, label=f"{name}: градиент по картам")
        xx = np.geomspace(ni.min(), ni.max(), 400)
        ax[0].plot(xx, model(q, xx, Bf, True), "-", color=col, lw=1.5, alpha=0.65, zorder=2,
                   label=f"{name}: β·B(α·NI+x₀)+y₀,  α={alpha:.3g} β={beta:.3g} "
                         f"x₀={x0:.3g} y₀={y0:.3g}")

        # ── отношение: что осталось после деления одной кривой на другую
        # По оси — ПРИВЕДЁННАЯ напряжённость α·NI, то есть глубина захода в кривую материала.
        # Сравнивать по току нельзя: серии покрывают разные диапазоны (сталь до 19000 А·в,
        # второй материал до 5500), и кривые оказались бы на разных участках.
        ratio = g / model(q, ni, Bf, True)
        ax[1].plot(alpha * ni, ratio, "-" + mk, color=col, ms=6, lw=1.4,
                   label=f"{name}  (наложение ±{err:.1f}%)")

    ax[0].set_xscale("log"); ax[0].set_yscale("log")
    ax[0].set_xlabel("ампер-витки"); ax[0].set_ylabel("градиент (у оси)")
    ax[0].set_title("Наложение: кривая градиента и приведённая кривая намагничивания")
    ax[0].grid(alpha=0.25, which="both"); ax[0].legend(fontsize=8.5, loc="lower right")

    ax[1].axhline(1.0, color="#888", lw=1.0, ls=":")
    ax[1].set_xscale("log")
    ax[1].set_xlabel("приведённая напряжённость α·NI, А/м — глубина захода в кривую материала")
    ax[1].set_ylabel("градиент / приведённая кривая материала")
    ax[1].set_title("Отношение: если бы оно задавалось геометрией, кривые совпали бы")
    ax[1].grid(alpha=0.25, which="both"); ax[1].legend(fontsize=8.5)

    fig.suptitle(f"{a.element}: если отношение задаётся ГЕОМЕТРИЕЙ, две кривые справа обязаны "
                 f"совпасть (градиент у оси, радиус {a.radius * 1e3:.0f} мм)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or ".", exist_ok=True)
    fig.savefig(a.out, dpi=140)
    print("график ->", a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
