# -*- coding: utf-8 -*-
"""Аналитическая модель «ток → кик» для ОДНООБМОТОЧНЫХ серий: магнитная цепь на кривой материала.

Модель записана как закон Ампера по замкнутому пути потока:

    NI = A·y + Σᵢ ℓᵢ · H(y / sᵢ)

где `y` — кик (интеграл момента вдоль элемента), `A` — воздушный член (зазор, линеен по потоку),
`ℓᵢ` и `sᵢ` — длина и эффективное сечение i-го участка железа, `H(·)` — кривая материала.

Главное здесь — что материал входит ОТДЕЛЬНЫМ множителем и ничем не подогнан: `A`, `ℓᵢ`, `sᵢ` это
геометрия, она от материала не зависит. Поэтому модель можно перенести на другой материал, просто
подменив `H`, и это проверяется (`--transfer`), а не постулируется.

Несколько участков нужны потому, что железо насыщается не разом: наконечник полюса раньше ярма,
ярмо раньше обратного пути. Каждый участок входит в кривую материала на своей глубине `y/sᵢ`.

Записано ТОКОМ ЧЕРЕЗ КИК, а не наоборот: так модель явная (решать уравнение не нужно) и совпадает
с метрикой — ошибкой по току.

    py -3.14 single_winding_model.py                      # все серии, участков 1..3
    py -3.14 single_winding_model.py --element QP_sh
    py -3.14 single_winding_model.py --transfer           # перенос на другой материал
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import re
import sys

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

DATA = os.path.join("workspace", "datasets", "currents_to_kicks")
MATERIALS = os.path.join("workspace", "materials")
# Серия карт того же QP_sh, посчитанная на ДРУГОМ материале — единственная в проекте пара
# «один магнит, два материала», и потому единственная прямая проверка переносимости.
ALT_DIR = os.path.join("C:" + os.sep, "Sasha", "Work", "SILA", "ansis approximations", "KD",
                       "QP_sh", "qp_sh_mu_steps")
LEAD = {"QP_sh": "G", "QP_L": "G", "OP": "O", "SP": "S", "DQ": "G", "corrector": "B"}
PROBE_R = {"B": 1e-3, "G": 5e-3, "S": 8e-3, "O": 14e-3}
ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}


# ───────────────────────────── материал ─────────────────────────────

def material(path):
    """H(B) — напряжённость по индукции, монотонно и с продолжением за таблицу.

    За последней точкой кривой железо ведёт себя как воздух: dB/dH = μ₀. Без этого продолжения
    модель обрывалась бы ровно там, где начинается глубокое насыщение, то есть в самом
    интересном месте."""
    d = np.loadtxt(path, delimiter=";", skiprows=1)
    H, B = d[:, 0], d[:, 1]
    mu0 = 4e-7 * np.pi

    def Hof(b):
        a = np.abs(np.asarray(b, float))
        out = np.where(a <= B[-1], np.interp(a, B, H), H[-1] + (a - B[-1]) / mu0)
        return np.sign(np.asarray(b, float)) * out
    return Hof


# ───────────────────────────── данные ─────────────────────────────

def series(element, dataset=DATA):
    """{обмотка: (ампер-витки, кик)} по ОДИНОЧНЫМ картам серии — обе величины по модулю.

    Комбинационные карты сюда не входят: одна обмотка — одна цепь, и смешивать их можно только
    после того, как одиночная модель заработает."""
    rows = list(csv.DictReader(open(os.path.join(dataset, f"{element}.csv"), encoding="utf-8"),
                               delimiter=";"))
    meta = json.load(open(os.path.join(dataset, f"{element}.meta.json"), encoding="utf-8"))
    mom = LEAD[element]
    out = {}
    for w in meta["windings"]:
        pts = []
        for r in rows:
            if ";" in r["excitation"]:
                continue
            at = float(r[f"AT_{w}"])
            if abs(at) < 1e-9 or any(abs(float(r[f"AT_{v}"])) > 1e-9
                                     for v in meta["windings"] if v != w):
                continue
            pts.append((abs(at), abs(float(r[f"int_{mom}"]))))
        if len(pts) >= 5:
            pts.sort()
            out[w] = (np.array([p[0] for p in pts]), np.array([p[1] for p in pts]))
    return out, meta


def alt_series(element="QP_sh", moment="G", folder=ALT_DIR):
    """Серия того же магнита на другом материале: (ампер-витки, кик) прямо из .fld."""
    from field_quality import load, _circle
    import _qp_scale as qs
    n, r = ORDER[moment], PROBE_R[moment]
    out = []
    for f in sorted(glob.glob(os.path.join(folder, f"{element}_*A.fld"))):
        m = re.search(r"_(\d+)A\.fld$", os.path.basename(f))
        if not m:
            continue
        h, B = load(f)
        if h["mnx"] >= -1e-9:
            h, B = qs.mirror_full(h, B)
        xs = h["mnx"] + np.arange(h["nx"]) * h["sx"]
        v = []
        for ix in range(h["nx"]):
            by, bz = _circle(h, B, ix, r, nth=64)
            sp = np.fft.fft(by + 1j * bz) / 64
            v.append(sp[0] if n == 1 else sp[-(n - 1)])
        v = np.array(v)
        out.append((float(m.group(1)), abs(np.trapezoid(v.imag, xs)) / r ** (n - 1)))
    out.sort()
    return np.array([p[0] for p in out]), np.array([p[1] for p in out])


# ───────────────────────────── модель ─────────────────────────────

def predict(theta, y, Hof, nseg):
    """Ампер-витки по кику: NI = A·y + Σ ℓᵢ·H(y/sᵢ). Параметры в логарифмах — все положительны."""
    A = np.exp(theta[0])
    ni = A * np.asarray(y, float)
    for i in range(nseg):
        ell, s = np.exp(theta[1 + 2 * i]), np.exp(theta[2 + 2 * i])
        ni = ni + ell * Hof(np.asarray(y, float) / s)
    return ni


def fit(ni, y, Hof, nseg, seeds=40, rng0=0):
    """Подгонка по ОТНОСИТЕЛЬНОЙ ошибке тока — это и есть метрика, которой модель меряется.

    Многостартово: задача невыпуклая, и участки железа легко сваливаются в один и тот же масштаб.
    Масштабы `sᵢ` упорядочены по построению (каждый следующий множителем > 1), иначе решения
    отличались бы только перестановкой участков и оптимизатор блуждал бы между ними."""
    def unpack(q):
        th = [q[0]]
        s = q[2]
        for i in range(nseg):
            th += [q[1 + 2 * i], s if i == 0 else s + np.log1p(np.exp(q[2 + 2 * i]))]
            s = th[-1]
        return np.array(th)

    def res(q):
        return (predict(unpack(q), y, Hof, nseg) - ni) / ni

    best = None
    ymax = float(np.max(y))
    for k in range(seeds):
        g = np.random.default_rng(rng0 + k)
        q0 = [np.log(np.mean(ni / y) * g.lognormal(0, 0.5))]
        for i in range(nseg):
            q0 += [np.log(0.05 * float(np.mean(ni)) * g.lognormal(0, 1.5)),
                   np.log(ymax / 1.5 * g.lognormal(0, 1.0)) if i == 0 else g.normal(0, 1.0)]
        try:
            t = least_squares(res, q0, xtol=1e-14, ftol=1e-14, max_nfev=20000)
        except Exception:
            continue
        if best is None or t.cost < best.cost:
            best = t
    return unpack(best.x), float(np.sqrt(np.mean(res(best.x) ** 2)) * 100)


def branch(theta, Hof, nseg, leak, umax, npts=4000):
    """Параметрическая кривая цепи С УТЕЧКОЙ: поток в железе u → (ампер-витки NI, кик y).

    Утечку гонит ВСЯ МДС катушки: поток замыкается вокруг обмотки, минуя рабочий зазор. Тогда
    железо несёт u = y + κ·NI, и из закона Ампера

        NI = A·y + Σ ℓᵢ·H(u/sᵢ),   y = u − κ·NI    ⇒    NI = [A·u + Σ ℓᵢ·H(u/sᵢ)] / (1 + A·κ)

    обе величины выражаются через u ЯВНО. Решать неявное уравнение не нужно — и не нужно гадать,
    сойдётся ли итерация: у жёсткого насыщения простая итерация расходится, а это ровно та область,
    ради которой утечка и вводится.

    Смысл прост: когда железо насыщено, добавка МДС почти целиком уходит в утечку, и кик встаёт.
    Последовательная цепь так не умеет — её наклон снизу ограничен воздушным членом."""
    A = np.exp(theta[0])
    u = np.linspace(umax / npts, umax, npts)
    fe = np.zeros_like(u)
    for i in range(nseg):
        ell, s = np.exp(theta[1 + 2 * i]), np.exp(theta[2 + 2 * i])
        fe = fe + ell * Hof(u / s)
    ni = (A * u + fe) / (1.0 + A * leak)
    return ni, u - leak * ni


def mmf(theta, y_query, Hof, nseg, leak=0.0):
    """Ампер-витки по кику. Без утечки — прямая формула, с утечкой — по параметрической кривой.

    Берётся ВОСХОДЯЩАЯ ветвь: при сильной утечке кик, пройдя максимум, начинает падать (лишняя МДС
    только добавляет утечки), и обратная задача на спадающей ветви смысла не имеет."""
    y_query = np.asarray(y_query, float)
    if leak <= 0:
        return predict(theta, y_query, Hof, nseg)
    umax = float(np.max(y_query)) * 4.0
    for _ in range(8):                       # растим сетку, пока кривая не накроет запрос
        ni, y = branch(theta, Hof, nseg, leak, umax)
        top = int(np.argmax(y))
        if y[top] >= y_query.max() or top < len(y) - 1:
            break
        umax *= 3.0
    ni, y = ni[:top + 1], y[:top + 1]
    if len(y) < 2:
        return np.full_like(y_query, 1e12)
    out = np.interp(y_query, y, ni)
    # Запрошенный кик может лежать ВЫШЕ достижимого — тогда продолжаем по последнему наклону, а не
    # ставим обрыв: на обрыве градиентная подгонка застревает, и утечка не может быть подобрана
    # вместе с геометрией (сечения железа обязаны расти вместе с ней).
    hi = y_query > y[-1]
    if np.any(hi):
        sl = (ni[-1] - ni[-2]) / max(y[-1] - y[-2], 1e-30)
        out = np.where(hi, ni[-1] + sl * (y_query - y[-1]), out)
    return out


def _pack(theta, nseg):
    """Параметры → свободные переменные оптимизатора (масштабы участков упорядочены)."""
    q = [theta[0]]
    for i in range(nseg):
        gap = theta[2 + 2 * i] - (theta[2 * i] if i else theta[2])
        q += [theta[1 + 2 * i], theta[2] if i == 0 else np.log(np.expm1(max(gap, 1e-6)))]
    return q


def _unpack(q, nseg):
    """Свободные переменные → параметры. Масштаб каждого следующего участка БОЛЬШЕ предыдущего:
    иначе решения отличались бы лишь перестановкой участков, и оптимизатор блуждал бы между ними."""
    th = [q[0]]
    s = q[2]
    for i in range(nseg):
        th += [q[1 + 2 * i], s if i == 0 else s + np.log1p(np.exp(q[2 + 2 * i]))]
        s = th[-1]
    return np.array(th)


def fit_leak(ni, y, Hof, nseg, seeds=30, kappas=None, rng0=0):
    """Подгонка цепи с утечкой: (параметры, κ, СКО ошибки тока, %).

    Подгонка ВЛОЖЕННАЯ: решение для N участков стартует с решения для N−1 плюс слабый новый
    участок. Без этого модель с бóльшим числом участков выходила ХУЖЕ модели с меньшим, хотя
    содержит её как частный случай, — верный признак того, что оптимизатор застревает, а не что
    участков много.

    По κ идёт скан с полной переподгонкой геометрии на каждом значении: сечения железа обязаны
    расти вместе с утечкой (часть потока в кик уже не попадает), поэтому подбирать κ отдельно от
    геометрии бессмысленно."""
    kappas = kappas if kappas is not None else [0.0] + list(np.geomspace(1e-7, 3e-3, 22))

    def solve(nseg_, kappa, starts):
        def res(q):
            return (mmf(_unpack(q, nseg_), y, Hof, nseg_, kappa) - ni) / ni
        best = None
        for q0 in starts:
            try:
                t = least_squares(res, q0, xtol=1e-14, ftol=1e-14, max_nfev=20000)
            except Exception:
                continue
            if np.all(np.isfinite(t.fun)) and (best is None or t.cost < best.cost):
                best = t
        return best

    ymax, nimean = float(np.max(y)), float(np.mean(ni))
    starts = []
    for k in range(seeds):
        g = np.random.default_rng(rng0 + k)
        q0 = [np.log(np.mean(ni / y) * g.lognormal(0, 0.5)),
              np.log(0.05 * nimean * g.lognormal(0, 1.5)),
              np.log(ymax / 1.5 * g.lognormal(0, 1.0))]
        starts.append(q0)

    best_all = None
    for n_ in range(1, nseg + 1):
        if n_ > 1:                                   # вложенный старт: прежнее решение + слабый участок
            th = _unpack(best_all[1], n_ - 1)
            starts = [list(_pack(th, n_ - 1)) + [th[1] - 3.0, 0.0],
                      list(_pack(th, n_ - 1)) + [th[1] - 1.0, 1.0]]
        for kappa in kappas:
            t = solve(n_, kappa, starts)
            if t is None:
                continue
            rms = float(np.sqrt(np.mean(t.fun ** 2)) * 100)
            if best_all is None or rms < best_all[0] or n_ == best_all[3]:
                if best_all is None or rms < best_all[0]:
                    best_all = (rms, t.x, kappa, n_)
        if best_all is not None and best_all[3] == n_:
            starts = [list(best_all[1])] + starts[:4]
    rms, q, kappa, n_ = best_all
    return _unpack(q, nseg if n_ == nseg else n_), kappa, rms, n_


def loo(ni, y, Hof, nseg):
    """Ошибка тока с исключением по одной карте: подгонка на всех точках о предсказании молчит."""
    err = []
    for i in range(len(ni)):
        m = np.ones(len(ni), bool)
        m[i] = False
        if m.sum() < 2 * nseg + 2:
            continue
        th, _ = fit(ni[m], y[m], Hof, nseg, seeds=12, rng0=100 * i)
        err.append(abs(predict(th, y[i], Hof, nseg) - ni[i]) / ni[i] * 100)
    return float(np.mean(err)) if err else float("nan")


def describe(theta, nseg, ymax, Hof):
    """Параметры в читаемом виде: доля МДС на воздух и глубина захода каждого участка в железо."""
    A = np.exp(theta[0])
    parts = [("воздух", A * ymax, None)]
    for i in range(nseg):
        ell, s = np.exp(theta[1 + 2 * i]), np.exp(theta[2 + 2 * i])
        parts.append((f"железо {i + 1}", float(ell * Hof(ymax / s)), ymax / s))
    tot = sum(p[1] for p in parts)
    return [(nm, v, 100 * v / tot, b) for nm, v, b in parts]


def plot_fits(elements, Hof, nseg=2, out=os.path.join("workspace", "results",
                                                      "single_winding_model.png")):
    """Кривая «момент от тока»: измерения против двух моделей — последовательной цепи и цепи с
    утечкой. Внизу — ошибка ПО ТОКУ, та самая метрика: по вертикальному расхождению кривых её не
    видно, потому что на плато кривая почти горизонтальна и любой промах по моменту превращается
    в огромный промах по току."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = []
    for el in elements:
        ser, meta = series(el)
        for w, (ni, y) in ser.items():
            panels.append((el, w, ni, y, meta))
    n = len(panels)
    fig, ax = plt.subplots(2, n, figsize=(4.4 * n, 7.6), squeeze=False)
    for k, (el, w, ni, y, meta) in enumerate(panels):
        tpc = meta["turns_per_coil"][w]
        nom = (meta["nominal_current_a"].get(w) or 0)
        th0, r0 = fit(ni, y, Hof, nseg, seeds=40)
        th1, kap, r1, _ = fit_leak(ni, y, Hof, nseg)
        ygrid = np.geomspace(y.min() * 0.5, y.max() * 1.02, 400)
        a = ax[0][k]
        a.plot(ni / tpc, y, "o", color="#111", ms=5, zorder=5, label="измерения")
        a.plot(predict(th0, ygrid, Hof, nseg) / tpc, ygrid, "-", color="#c0392b", lw=1.7,
               label=f"последовательная цепь, {r0:.1f}%")
        a.plot(mmf(th1, ygrid, Hof, nseg, kap) / tpc, ygrid, "-", color="#1f4e79", lw=1.7,
               label=f"с утечкой, {r1:.1f}%")
        if nom:
            a.axvline(nom, color="#888", lw=1.0, ls="--")
            a.annotate("номинал", xy=(nom, 0.02), xycoords=("data", "axes fraction"),
                       fontsize=8, color="#666", rotation=90, va="bottom", ha="right")
        a.set_xlim(0, ni.max() / tpc * 1.05)
        a.set_ylim(0, y.max() * 1.15)
        a.set_title(f"{el} {w}", fontsize=11)
        a.set_ylabel(f"{LEAD[el]}, интеграл" if k == 0 else "")
        a.grid(alpha=0.25)
        a.legend(fontsize=8, loc="lower right")
        b = ax[1][k]
        b.axhline(0, color="#888", lw=1.0)
        b.plot(ni / tpc, 100 * (predict(th0, y, Hof, nseg) - ni) / ni, "-o", color="#c0392b",
               ms=4, lw=1.3)
        b.plot(ni / tpc, 100 * (mmf(th1, y, Hof, nseg, kap) - ni) / ni, "-o", color="#1f4e79",
               ms=4, lw=1.3)
        if nom:
            b.axvline(nom, color="#888", lw=1.0, ls="--")
        b.set_xlabel("ток обмотки, А")
        b.set_ylabel("ошибка по току, %" if k == 0 else "")
        b.grid(alpha=0.25)
    fig.suptitle("Ток → кик: измерения и аналитическая модель на кривой материала "
                 f"(участков железа {nseg})", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ───────────────────── граф магнита: 4 узла, две катушки ─────────────────────
#
#          ┌───────── утечка ─────────┐
#   (1) ───┴──[ железо, катушка I1 ]──┴─── (2)
#    │                                      │
#  зазор                                  ярмо
#    │                                      │
#   (4) ───┬──[ железо, катушка I1 ]──┬─── (3)
#          └───────── утечка ─────────┘
#
# Контур 1→2→3→4→1: две ОДИНАКОВЫЕ катушки на ветвях 1–2 и 3–4 (одна обмотка, две катушки),
# ярмо 2–3 со своим сечением, рабочий зазор 4–1. Параллельно каждой катушечной ветви — воздушная
# утечка. Поток через зазор и есть кик (связь вырождается в коэффициент, он поглощён единицами).

# ───────────────────── граф магнита: контур с двумя катушками ─────────────────────
#
#              ┌──── утечка ────┐
#   (1)───[уч.]┴(2)──[ I1 ]─────┴(3)
#    │                            │
#  зазор                        ярмо
#    │                            │
#   (6)───[уч.]┬(5)──[ I1 ]─────┬(4)
#              └──── утечка ────┘
#
# Замкнутый контур 1→2→3→4→5→6→1. Две ОДИНАКОВЫЕ катушки одной обмотки на ветвях 2–3 и 4–5,
# параллельно каждой — воздушная утечка. Ветви 1–2 и 5–6 — участки между полюсами и ярмом, 3–4 —
# ярмо, 6–1 — рабочий зазор. Поток зазора и есть кик.
#
# Ключ к замкнутой форме: все ветви, кроме полюсов, несут ОДИН И ТОТ ЖЕ поток — поток зазора.
# Полюса несут свой, потому что часть потока уходит в утечку рядом с катушкой. Поэтому сколько бы
# последовательных участков ни добавляли (наконечник, ярмо, торцевая плита), они лишь пополняют
# сумму Q, и структура решения не меняется.


def graph_mmf(par, y, Hof, nseries=2):
    """Ампер-витки одной катушки по кику — ЗАМКНУТОЙ ФОРМУЛОЙ, без итераций.

        Q(y) = Σₖ mₖ·ℓₖ·H(y/Sₖ) + Rз·y        падения на последовательных ветвях и зазоре
        p(y) = y + (P/2)·Q(y)                  поток полюса: полезный плюс утечка
        F(y) = ℓп·H(p/Sп) + Q(y)/2             ампер-витки катушки

    Вывод: закон Ампера по контуру даёт 2F = 2·ℓп·H(p/Sп) + Q, откуда падение на полюсе равно
    F − Q/2. Утечка стоит параллельно катушечной ветви и видит ту же разность потенциалов, значит
    её поток равен P·(D − F) = −P·Q/2, а узел на входе в полюс даёт p = y + P·Q/2.

    Множитель mₖ — сколько РАВНЫХ ветвей этого класса в контуре: участки полюс–ярмо входят дважды
    (сверху и снизу), ярмо один раз. Считать их порознь незачем — по симметрии они одинаковы, а
    лишние свободные числа только размывают подгонку.

    Отсюда же плато: полюс несёт p = y + (P/2)·Q, то есть с ростом МДС нагружается сверх полезного
    потока и, насытившись, забирает всю добавку — кик встаёт. Последовательная цепь так не умеет.

    `par` в логарифмах: [ℓп, Sп, Rз, P, (ℓ₁, S₁), (ℓ₂, S₂), …]."""
    v = np.exp(np.asarray(par, float))
    lp, Sp, Rg, P = v[0], v[1], v[2], v[3]
    y = np.asarray(y, float)
    Q = Rg * y
    for k in range(nseries):
        mult = SERIES_MULT[k] if k < len(SERIES_MULT) else 1.0
        Q = Q + mult * v[4 + 2 * k] * Hof(y / v[5 + 2 * k])
    p = y + 0.5 * P * Q
    return lp * Hof(p / Sp) + 0.5 * Q


# Сколько равных ветвей каждого последовательного класса в контуре: сперва ярмо (одно), затем
# участки полюс–ярмо (их два — сверху и снизу).
SERIES_MULT = (1.0, 2.0, 1.0, 1.0)


# Индукция в железе физически лежит в этих пределах — отсюда и коробка на сечения: S ≈ y_max/B.
B_RANGE = (0.4, 3.0)


def _lengths(y, ni, Hof, S, P, nseries, passes=6):
    """Длины ветвей и сопротивление зазора при ЗАДАННЫХ сечениях и утечке — линейным МНК.

    При фиксированных сечениях ампер-витки линейны по длинам и по сопротивлению зазора:

        F = ℓп·H(p/Sп) + Rз·(y/2) + Σ mₖ·ℓₖ·H(y/Sₖ)/2

    поэтому искать их перебором незачем — они решаются точно. Нелинейным остаётся только поток
    полюса p, который сам зависит от найденных длин через утечку; он берётся простой итерацией
    (обычно хватает трёх проходов, начиная с p = y).

    Неотрицательность обязательна: отрицательная длина или сечение — не «немного хуже», а другая
    физика, и подгонка с ними уходит в бессмысленные минимумы."""
    from scipy.optimize import nnls
    y = np.asarray(y, float)
    p = y.copy()
    Sp, Sk = S[0], S[1:]
    coef = None
    for _ in range(passes):
        cols = [Hof(p / Sp), 0.5 * y]
        for k in range(nseries):
            mult = SERIES_MULT[k] if k < len(SERIES_MULT) else 1.0
            cols.append(0.5 * mult * Hof(y / Sk[k]))
        A = np.column_stack(cols) / ni[:, None]          # веса: метрика — ОТНОСИТЕЛЬНАЯ ошибка тока
        coef, _ = nnls(A, np.ones_like(ni))
        Q = coef[1] * y + sum(coef[2 + k] * (SERIES_MULT[k] if k < len(SERIES_MULT) else 1.0)
                              * Hof(y / Sk[k]) for k in range(nseries))
        p_new = y + 0.5 * P * Q
        if np.allclose(p_new, p, rtol=1e-10, atol=0):
            p = p_new
            break
        p = p_new
    return coef


def _pack_from(S, P, coef, nseries):
    """Сечения + утечка + найденные длины → вектор параметров `graph_mmf`."""
    par = [np.log(max(coef[0], 1e-300)), np.log(S[0]), np.log(max(coef[1], 1e-300)),
           np.log(max(P, 1e-300))]
    for k in range(nseries):
        par += [np.log(max(coef[2 + k], 1e-300)), np.log(S[1 + k])]
    return np.array(par)


def fit_graph_sep(ni, y, Hof, nseries=2, seeds=8, nkappa=17, rng0=0):
    """Подгонка графа РАЗДЕЛЁННЫМ методом: (параметры, СКО %).

    Длины и сопротивление зазора решаются линейным МНК внутри, наружу выходит поиск только по
    СЕЧЕНИЯМ — с 8 свободных чисел до 2–3. Сечения ищутся в физической коробке (индукция железа
    0.4…3 Тл), а не в бесконечном логарифме, поэтому случайные старты сразу попадают в осмысленную
    область.

    Утечка не подгоняется вместе с ними, а СКАНИРУЕТСЯ по сетке с полной переподгонкой остального
    на каждом значении. Причина: она входит в задачу иначе, чем сечения — через поток полюса, — и
    при совместном поиске оптимизатор устойчиво сваливал её в ноль, вырождая граф в последовательную
    цепь (те самые 18.5%).

    Наружу отдаётся ВЕКТОР невязок, а не их среднеквадратичное: по скаляру least_squares работает
    вслепую, теряя всю структуру задачи."""
    ymax, nimax = float(np.max(y)), float(np.max(ni))
    lo, hi = np.log(ymax / B_RANGE[1]), np.log(ymax / B_RANGE[0])
    kaps = [0.0] + list(np.geomspace(ymax / nimax * 1e-3, ymax / nimax * 30, nkappa - 1))

    def solve(q, P):
        S = np.exp(q)
        coef = _lengths(y, ni, Hof, S, P, nseries)
        par = _pack_from(S, P, coef, nseries)
        return par, (graph_mmf(par, y, Hof, nseries) - ni) / ni

    best = None
    for P in kaps:
        for k in range(seeds):
            g = np.random.default_rng(rng0 + k)
            q0 = list(g.uniform(lo, hi, 1 + nseries))
            try:
                t = least_squares(lambda q: solve(q, P)[1], q0, xtol=1e-12, ftol=1e-12,
                                  diff_step=1e-5, max_nfev=2000,
                                  bounds=([lo] * (1 + nseries), [hi] * (1 + nseries)))
            except Exception:
                continue
            par, r = solve(t.x, P)
            rms = float(np.sqrt(np.mean(r ** 2)) * 100)
            if np.isfinite(rms) and (best is None or rms < best[0]):
                best = (rms, par)
    return best[1], best[0]


def fit_graph(ni, y, Hof, nseries=2, seeds=60, rng0=0):
    """Подгонка графа по ошибке тока: (параметры, СКО %).

    Подгонка ВЛОЖЕННАЯ по числу последовательных участков: решение для N стартует с решения для
    N−1 плюс слабый новый участок. Без этого модель с бóльшим числом ветвей выходила бы хуже
    модели с меньшим, хотя содержит её как частный случай."""
    def res_n(q, n):
        return (graph_mmf(q, y, Hof, n) - ni) / ni

    ymax, nimean = float(np.max(y)), float(np.mean(ni))
    best = None
    for n in range(0 if nseries == 0 else 1, nseries + 1):
        starts = []
        if best is not None:
            starts.append(list(best[1]) + [best[1][4] - 3.0, best[1][5]])
            starts.append(list(best[1]) + [best[1][4] - 1.0, best[1][5] + 1.0])
        for k in range(seeds if best is None else max(seeds // 3, 8)):
            g = np.random.default_rng(rng0 + 1000 * n + k)
            q0 = [np.log(0.05 * nimean * g.lognormal(0, 1.5)),
                  np.log(ymax / 1.5 * g.lognormal(0, 1.0)),
                  np.log(float(np.mean(ni / y)) * g.lognormal(0, 0.5)),
                  np.log(ymax / float(np.max(ni)) * g.lognormal(0, 2.0))]
            for _ in range(n):
                q0 += [np.log(0.05 * nimean * g.lognormal(0, 1.5)),
                       np.log(ymax / 1.5 * g.lognormal(0, 1.0))]
            starts.append(q0)
        cur = None
        for q0 in starts:
            try:
                t = least_squares(lambda q: res_n(q, n), q0, xtol=1e-14, ftol=1e-14,
                                  max_nfev=20000)
            except Exception:
                continue
            if np.all(np.isfinite(t.fun)) and (cur is None or t.cost < cur.cost):
                cur = t
        if cur is not None:
            rms = float(np.sqrt(np.mean(cur.fun ** 2)) * 100)
            if best is None or rms <= best[0] + 1e-9:
                best = (rms, cur.x, n)
    return best[1], best[0]


GRAPH_LEVELS = ((0, "простейший: сердечник + катушка + утечка", "#c0392b"),
                (1, "+ ярмо отдельной ветвью", "#e8a33d"),
                (2, "+ участки полюс–ярмо", "#1f4e79"))


def plot_graph_levels(elements, Hof, seeds=60,
                      out=os.path.join("workspace", "results", "graph_levels.png")):
    """Насколько граф вообще нужно усложнять: измерения против трёх уровней подробности.

    Простейший уровень (один сердечник, одна катушка, утечка на ней) приведён не для полноты: в нём
    утечка ТОЖДЕСТВЕННО вырождается в сечение железа. Разность потенциалов на катушечной ветви там
    равна падению на зазоре, то есть линейна по кику, поэтому поток полюса выходит строго
    пропорционален кику: p = y·(1 + P·Rз/2). Утечка целиком поглощается заменой S → S/(1+P·Rз/2),
    и граф вырождается в последовательную цепь. Плато он дать не может ни при каких параметрах."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = []
    for el in elements:
        ser, meta = series(el)
        for w, (ni, y) in ser.items():
            panels.append((el, w, ni, y, meta))
    n = len(panels)
    fig, ax = plt.subplots(2, n, figsize=(4.6 * n, 7.8), squeeze=False)
    for k, (el, w, ni, y, meta) in enumerate(panels):
        tpc = meta["turns_per_coil"][w]
        nom = meta["nominal_current_a"].get(w) or 0
        g = np.geomspace(y.min() * 0.5, y.max() * 1.02, 400)
        a, b = ax[0][k], ax[1][k]
        a.plot(ni / tpc, y, "o", color="#111", ms=5.5, zorder=6, label="измерения")
        b.axhline(0, color="#888", lw=1.0)
        for ns, name, col in GRAPH_LEVELS:
            par, rms = fit_graph(ni, y, Hof, nseries=ns, seeds=seeds)
            a.plot(graph_mmf(par, g, Hof, ns) / tpc, g, "-", color=col, lw=1.8,
                   label=f"{name} — {rms:.1f}%")
            b.plot(ni / tpc, 100 * (graph_mmf(par, y, Hof, ns) - ni) / ni, "-o", color=col,
                   ms=4, lw=1.3)
        for c in (a, b):
            if nom:
                c.axvline(nom, color="#888", lw=1.0, ls="--")
            c.grid(alpha=0.25)
        a.set_xlim(0, ni.max() / tpc * 1.05)
        a.set_ylim(0, y.max() * 1.15)
        a.set_title(f"{el} {w}", fontsize=11)
        a.set_ylabel("кик (интеграл момента)" if k == 0 else "")
        a.legend(fontsize=8, loc="lower right")
        b.set_xlabel("ток обмотки, А")
        b.set_ylabel("ошибка по току, %" if k == 0 else "")
    fig.suptitle("Насколько подробным должен быть граф магнита", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def best_fit(ni, y, Hof, nseries=2, tries=5, seeds=60):
    """Лучшая из нескольких подгонок с разными зёрнами.

    Ландшафт изрезан вырожденностями (ветви, несущие один и тот же поток, различимы только суммой),
    и одиночный многостарт от прогона к прогону садится в разные минимумы — у QP_sh встречались и
    6.6%, и 10.2% на одних данных. Поэтому берём лучшее из нескольких независимых прогонов."""
    best = None
    for s in range(tries):
        par, rms = fit_graph(ni, y, Hof, nseries=nseries, seeds=seeds, rng0=s * 9173)
        if best is None or rms < best[1]:
            best = (par, rms)
    return best


MOMENT_ORDER = {"B": 1, "G": 2, "S": 3, "O": 4}
MOMENT_UNIT = {"B": "Тл", "G": "Тл/м", "S": "Тл/м²", "O": "Тл/м³"}


def coef_scale(element, dataset=DATA):
    """Множитель «интеграл момента → физкоэффициент» для этого элемента.

    Датасет хранит интеграл вдоль элемента (длина в него не входит). Привычная величина — сам
    коэффициент, то есть интеграл, делённый на длину и на энергетический множитель, с (n−1)! по
    соглашению столбцов mat_*."""
    import math
    meta = json.load(open(os.path.join(dataset, f"{element}.meta.json"), encoding="utf-8"))
    n = MOMENT_ORDER[LEAD[element]]
    return (math.factorial(n - 1) / (meta["L_used_for_mat_m"] * meta["energy_factor"]),
            meta["L_used_for_mat_m"])


def plot_model(elements, Hof, nseries=2, tries=5,
               out=os.path.join("workspace", "results", "model_fit.png")):
    """Измерения и модель «участки + ярмо» по обмоткам, с подобранной геометрией в подписи.

    По вертикали — физкоэффициент в привычных единицах (у секступоля Тл/м²), а не интеграл: модель
    подгоняется на интеграле, но читать её удобнее в том же виде, что и проектные значения."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    panels = []
    for el in elements:
        ser, meta = series(el)
        for w, (ni, y) in ser.items():
            panels.append((el, w, ni, y, meta))
    n = len(panels)
    fig, ax = plt.subplots(2, n, figsize=(5.0 * n, 7.8), squeeze=False)
    for k, (el, w, ni, y, meta) in enumerate(panels):
        tpc = meta["turns_per_coil"][w]
        nom = meta["nominal_current_a"].get(w) or 0
        par, rms = best_fit(ni, y, Hof, nseries, tries)
        v = np.exp(par)
        g = np.geomspace(y.min() * 0.5, y.max() * 1.02, 500)
        a, b = ax[0][k], ax[1][k]
        sc, Lm = coef_scale(el)
        a.plot(ni / tpc, y * sc, "o", color="#111", ms=6, zorder=5, label="измерения")
        a.plot(graph_mmf(par, g, Hof, nseries) / tpc, g * sc, "-", color="#1f4e79", lw=2.0,
               label=f"модель «участки + ярмо» — {rms:.1f}%")
        b.axhline(0, color="#888", lw=1.0)
        b.plot(ni / tpc, 100 * (graph_mmf(par, y, Hof, nseries) - ni) / ni, "-o",
               color="#1f4e79", ms=4.5, lw=1.4)
        for c in (a, b):
            if nom:
                c.axvline(nom, color="#888", lw=1.0, ls="--")
                c.annotate("номинал", xy=(nom, 0.02), xycoords=("data", "axes fraction"),
                           fontsize=8, color="#666", rotation=90, va="bottom", ha="right")
            c.grid(alpha=0.25)
        ymax = float(np.max(y))
        txt = (f"полюс  B={ymax / v[1]:.2f} Тл\n"
               f"ярмо   B={ymax / v[5]:.2f} Тл\n"
               + (f"участки B={ymax / v[7]:.2f} Тл\n" if nseries > 1 else "")
               + f"утечка P={v[3]:.2e}")
        a.annotate(txt, xy=(0.02, 0.97), xycoords="axes fraction", fontsize=8.5, va="top",
                   family="monospace", bbox=dict(fc="white", ec="#ccc", alpha=0.9, pad=3))
        a.set_xlim(0, ni.max() / tpc * 1.05)
        a.set_ylim(0, ymax * sc * 1.15)
        a.set_title(f"{el}, обмотка {w}   (L_eff = {Lm:.3f} м)", fontsize=11)
        a.set_ylabel(f"{LEAD[el]}, {MOMENT_UNIT[LEAD[el]]}" if k == 0 else "")
        a.legend(fontsize=9, loc="lower right")
        b.set_xlabel("ток обмотки, А")
        b.set_ylabel("ошибка по току, %" if k == 0 else "")
    fig.suptitle("Модель магнитной цепи на кривой материала: измерения и подгонка", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


def plot_transfer(nseg=2, base="steel_1008_bh.csv", alt="lowsat_custom_bh.csv",
                  out=os.path.join("workspace", "results", "material_transfer.png")):
    """Перенос материала на графике: подгонка на одной стали, предсказание на другой.

    Смысл рисунка в том, что синяя кривая — НЕ подгонка. Её параметры целиком взяты у красной,
    заменена только кривая материала. Пунктир показывает, что было бы без замены."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    Hb, Ha = (material(os.path.join(MATERIALS, base)), material(os.path.join(MATERIALS, alt)))
    ser, meta = series("QP_sh")
    ni, y = ser["I"]
    ni_a, y_a = alt_series("QP_sh", LEAD["QP_sh"])
    tpc = meta["turns_per_coil"]["I"]
    m = ni <= ni_a.max() * 1.6
    th, kap, rms, _ = fit_leak(ni[m], y[m], Hb, nseg)
    pred = mmf(th, y_a, Ha, nseg, kap)
    ctrl = mmf(th, y_a, Hb, nseg, kap)

    fig, ax = plt.subplots(1, 3, figsize=(15.2, 4.8))
    g = np.geomspace(min(y[m].min(), y_a.min()) * 0.5, max(y[m].max(), y_a.max()) * 1.05, 400)
    a = ax[0]
    a.plot(ni[m] / tpc, y[m], "o", color="#c0392b", ms=6, zorder=5, label="сталь 1008, измерения")
    a.plot(mmf(th, g, Hb, nseg, kap) / tpc, g, "-", color="#c0392b", lw=1.8,
           label=f"модель, подогнана ({rms:.2f}%)")
    a.plot(ni_a / tpc, y_a, "s", color="#1f4e79", ms=6, zorder=5,
           label="низконасыщ., измерения")
    a.plot(mmf(th, g, Ha, nseg, kap) / tpc, g, "-", color="#1f4e79", lw=1.8,
           label=f"модель, ПРЕДСКАЗАНА ({np.mean(np.abs(pred - ni_a) / ni_a) * 100:.1f}%)")
    a.plot(mmf(th, g, Hb, nseg, kap) / tpc, g, ":", color="#888", lw=1.6,
           label="если материал не менять")
    a.set_xlim(0, max(ni[m].max(), ni_a.max()) / tpc * 1.1)
    a.set_ylim(0, max(y[m].max(), y_a.max()) * 1.15)
    a.set_xlabel("ток обмотки, А")
    a.set_ylabel("G, интеграл")
    a.set_title("QP_sh: один магнит, два материала", fontsize=11)
    a.grid(alpha=0.25)
    a.legend(fontsize=8.5, loc="lower right")

    a = ax[1]
    a.axhline(0, color="#888", lw=1.0)
    a.plot(ni_a / tpc, 100 * (pred - ni_a) / ni_a, "-s", color="#1f4e79", ms=6, lw=1.5,
           label="со сменой материала")
    a.plot(ni_a / tpc, 100 * (ctrl - ni_a) / ni_a, "-o", color="#888", ms=5, lw=1.4, ls=":",
           label="без смены (контроль)")
    a.set_xlabel("ток обмотки, А")
    a.set_ylabel("ошибка по току, %")
    a.set_title("Ошибка предсказания на новом материале", fontsize=11)
    a.grid(alpha=0.25)
    a.legend(fontsize=9)

    a = ax[2]
    for nm, path, col in ((base, base, "#c0392b"), (alt, alt, "#1f4e79")):
        d = np.loadtxt(os.path.join(MATERIALS, path), delimiter=";", skiprows=1)
        a.plot(d[:, 0], d[:, 1], "-o", color=col, ms=3.5, lw=1.6, label=nm)
    for i in range(nseg):                       # где по кривой работают участки железа
        s = np.exp(th[2 + 2 * i])
        a.axhline(float(np.max(y[m])) / s, color="#444", lw=1.0, ls="--")
        a.annotate(f"железо {i + 1}: {np.max(y[m]) / s:.2f} Тл", xy=(1.2, np.max(y[m]) / s),
                   fontsize=8, color="#444", va="bottom")
    a.set_xscale("log")
    a.set_xlabel("H, А/м")
    a.set_ylabel("B, Тл")
    a.set_title("Кривые материалов — единственное, что менялось", fontsize=11)
    a.grid(alpha=0.25, which="both")
    a.legend(fontsize=9, loc="lower right")

    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(out)) or ".", exist_ok=True)
    fig.savefig(out, dpi=145)
    plt.close(fig)
    return out


# ───────────────────────────── прогоны ─────────────────────────────

def run_fits(elements, Hof, segs, do_loo=True):
    """Таблица «участков железа → ошибка тока» по всем однообмоточным сериям."""
    print(f"{'серия':14}{'карт':>5}" + "".join(f"{f'участков {n}':>22}" for n in segs))
    print(f"{'':14}{'':>5}" + "".join(f"{'подгонка / без точки':>22}" for _ in segs))
    out = {}
    for el in elements:
        ser, _ = series(el)
        for w, (ni, y) in ser.items():
            cells = []
            for n in segs:
                th, rms = fit(ni, y, Hof, n)
                lo = loo(ni, y, Hof, n) if do_loo else float("nan")
                cells.append(f"{rms:8.2f}% /{lo:8.2f}%")
                out[(el, w, n)] = (th, rms, lo)
            print(f"{el + ' ' + w:14}{len(ni):5d}" + "".join(f"{c:>22}" for c in cells))
    return out


def run_transfer(nseg=2, base="steel_1008_bh.csv", alt="lowsat_custom_bh.csv"):
    """Перенос на ДРУГОЙ МАТЕРИАЛ: геометрия подобрана на одном, предсказание — на другом.

    Это и есть проверка требования «менять материал». Параметры A, ℓᵢ, sᵢ и κ берутся с серии на
    стали и НЕ подгоняются заново — меняется только H. Если ошибка на второй серии сравнима с
    ошибкой на первой, материал в модели действительно отделён от геометрии; если нет, значит
    что-то из «геометрии» на самом деле впитало в себя свойства стали."""
    Hb = material(os.path.join(MATERIALS, base))
    Ha = material(os.path.join(MATERIALS, alt))
    ser, _ = series("QP_sh")
    ni, y = ser["I"]
    ni_a, y_a = alt_series("QP_sh", LEAD["QP_sh"])
    m = ni <= ni_a.max() * 1.6            # сопоставимый диапазон: вторая серия кончается раньше
    print(f"перенос QP_sh: {base} → {alt};  участков железа {nseg}")
    print(f"  подгонка на {base}: {m.sum()} карт до {ni[m].max():.0f} А·в")
    th, kap, rms, _ = fit_leak(ni[m], y[m], Hb, nseg)
    print(f"  ошибка тока на своём материале: {rms:.2f}%   κ={kap:.3g}")
    for nm, v, pct, b in describe(th, nseg, float(np.max(y[m])), Hb):
        print(f"     {nm:10} МДС {v:10.0f} А·в ({pct:5.1f}%)"
              + (f", B = {b:.3f} Тл" if b is not None else ""))
    pred = mmf(th, y_a, Ha, nseg, kap)
    err = np.abs(pred - ni_a) / ni_a * 100
    same = np.abs(mmf(th, y_a, Hb, nseg, kap) - ni_a) / ni_a * 100
    print(f"\n  предсказание на {alt} БЕЗ подгонки: ошибка тока {err.mean():.2f}% "
          f"(медиана {np.median(err):.2f}%, максимум {err.max():.2f}%)")
    print(f"     {'А·в факт':>10}{'модель':>10}{'ошибка':>9}{'кик':>12}")
    for a_, b_, e, yy in zip(ni_a, pred, err, y_a):
        print(f"     {a_:10.0f}{b_:10.0f}{e:8.1f}%{yy:12.4g}")
    print(f"  контроль — та же геометрия, но материал НЕ подменён: {same.mean():.2f}%; "
          f"подмена улучшает в {same.mean() / max(err.mean(), 1e-9):.1f} раза")
    return th, kap


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", help="только одна серия")
    ap.add_argument("--material", default="steel_1008_bh.csv")
    ap.add_argument("--segments", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--transfer", action="store_true", help="проверка переноса на другой материал")
    ap.add_argument("--no-loo", action="store_true")
    a = ap.parse_args(argv)

    Hof = material(os.path.join(MATERIALS, a.material))
    els = [a.element] if a.element else ["QP_sh", "QP_L", "OP", "SP", "DQ", "corrector"]
    print(f"материал: {a.material};  метрика — ошибка ТОКА, %\n")
    run_fits(els, Hof, a.segments, do_loo=not a.no_loo)
    if a.transfer:
        print()
        run_transfer(max(a.segments))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
