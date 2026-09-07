# -*- coding: utf-8 -*-
"""Модель «токи обмоток → кики» с МАТРИЦЕЙ СВЯЗИ вместо геометрических направлений трубок.

Отличается от МЭС-σ (`mec_sigma`) ровно одним звеном. Там связь «потоки → кик» строится из
КООРДИНАТ наконечников: каждому зазору приписано направление трубки, и кик есть проекция потока на
это направление. Здесь координаты из модели убраны: кик получается ПОДБИРАЕМОЙ матрицей из потоков
наконечников. Всё остальное — граф из чертежей, решение цепи с кривой намагничивания, доля
насыщения σ, общие знаки намотки — то же самое.

Наконечник — узел графа, к которому подходит зазор кик-класса; поток наконечника — алгебраическая
сумма потоков сходящихся в нём зазоров (топология, не геометрия). Признаки на одно решение цепи:

    x = [ q, σ·q ],   q — потоки наконечников, σ — доля насыщения железа

Матрица связи отображает x в компоненты кика (нормальные и скошенные) СРАЗУ ДЛЯ ВСЕХ ГАРМОНИК —
одно решение цепи от гармоники не зависит. Матрица подбирается точным МНК на одиночных обмотках;
комбинации остаются независимой проверкой.

Маска симметрии. У матрицы связи, в отличие от геометрии, нет встроенного знания о том, какая
обмотка какой кик создать МОЖЕТ: главная обмотка корректора сидит на ярме симметрично и скошенных
составляющих не даёт, но их шум на её картах матрица подгоняла кривой. `MASK` возвращает это знание
топологически — у запрещённой компоненты цель подгонки заменяется точным нулём, и матрица вынуждена
обнулять её на характерном для обмотки распределении потоков (проверено: скош. S у main падает с
рабочих 1e-2 до 1e-7).

    py -3.14 coupling_model.py --acc          # точность по всем элементам
    py -3.14 coupling_model.py --free         # свободные числа
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import branch_network as bn
import mec_sigma as ms

ORDERS_OF = {"DQ": (1, 2, 3), "corrector": (1, 2, 3)}
PARTS_OF = {"DQ": ("норм.",), "corrector": ("норм.", "скош.")}
_PF = {"норм.": (lambda z: z.imag), "скош.": (lambda z: z.real)}

# Маска «обмотка → какие компоненты кика она может создавать» (значение — множество разрешённых
# частей; обмотки, которых нет, могут создавать всё). Несёт априорное знание симметрии.
MASK = {"corrector": {"main": {"норм."}}}
_FIT = {}


def _allowed(element, wname, part):
    m = MASK.get(element, {}).get(wname)
    return True if m is None else part in m


def allowed(element, winding, comp):
    """Может ли обмотка создать эту компоненту: comp=0 — нормальную, 1 — скошенную.

    Симметрия обмотки запрещает часть гармоник (`MASK`); модель их зануляет, а в измерениях на их
    месте стоит численный шум карты. Отчёту это нужно, чтобы не рисовать шум как величину."""
    return _allowed(element, winding, "скош." if comp else "норм.")


def tips(cir):
    """Наконечники (узлы, к которым подходят зазоры кик-класса) и сами эти зазоры."""
    gaps = [b for b in cir.gaps if cir.sym[b] in cir.kick_classes]
    return sorted({n for b in gaps for n in cir.edge[b]}), gaps


def features(cir, currents):
    """Признаки без координат: потоки наконечников и они же, умноженные на долю насыщения σ."""
    lens, perms, secs = bn.unpack(np.zeros(1), cir)
    phi, _ = cir.solve(currents, lens, perms, secs=secs)
    ts, gaps = tips(cir)
    q = {t: 0.0 for t in ts}
    for b in gaps:
        a_, c_ = cir.edge[b]
        q[a_] += phi[b]
        q[c_] -= phi[b]
    num = den = 0.0
    for x in cir.iron:
        S = secs.get(cir.sym[x], 1.0)
        a = abs(phi[x])
        num += a * min(a / S / cir.Barr[-1], 1.0)
        den += a
    qv = np.array([q[t] for t in ts])
    return np.concatenate([qv, (num / den if den else 0.0) * qv])


def columns(element):
    return [(n, p) for n in ORDERS_OF[element] for p in PARTS_OF[element]]


def train(element):
    """Матрица связи, подогнанная на одиночных картах. Возврат (цепь, матрица W, список колонок)."""
    if element in _FIT:
        return _FIT[element]
    cir = bn.Circuit(element)
    if cir.geom is None:
        raise SystemExit(f"{element}: геометрия графа не задана")
    sen, _ = ms.joint_senses(element, ORDERS_OF[element])
    cir.sense.update(sen)
    cols = columns(element)
    maps = {n: bn.read_maps(element, n) for n in ORDERS_OF[element]}
    ids = [t[3] for t in maps[ORDERS_OF[element][0]]]
    cur = {t[3]: t[0] for t in maps[ORDERS_OF[element][0]]}
    single = [i for i in ids if len(cur[i]) == 1]
    X = np.array([features(cir, cur[i]) for i in single])
    Y = np.array([[_PF[p](complex({t[3]: t[1] for t in maps[n]}[i])) for n, p in cols]
                  for i in single])
    for r, i in enumerate(single):                     # маска: запрещённой компоненте цель — ноль
        w = list(cur[i])[0]
        for c, (n, p) in enumerate(cols):
            if not _allowed(element, w, p):
                Y[r, c] = 0.0
    # Подгонка по ОТНОСИТЕЛЬНОЙ ошибке: вес карты 1/|цель|. Обычный МНК минимизирует абсолютную
    # невязку, и колонка целиком уходит под самые сильные карты — у main корректора секступоль при
    # 100 А·в (0.2) против него же при 13000 (16.7) не весит ничего, и модель занижала его втрое.
    # Пол в 2% от максимума не даёт весу уйти в бесконечность на околонулевых картах; строки,
    # занулённые маской, получают наибольший вес — запрет должен выполняться точно.
    # Выравнивать веса ПО ОБМОТКАМ нельзя: это ломает их относительный масштаб, а именно его
    # проверяют комбинации (ошибка на них взлетает с 24% до 3800%).
    W = np.zeros((X.shape[1], Y.shape[1]))
    for c in range(Y.shape[1]):
        y = Y[:, c]
        top = float(np.max(np.abs(y)))
        if top <= 0:
            continue
        wt = 1.0 / np.maximum(np.abs(y), 0.02 * top)
        W[:, c] = np.linalg.lstsq(X * wt[:, None], y * wt, rcond=None)[0]
    _FIT[element] = (cir, W, cols)
    return _FIT[element]


def predict(element, order, currents):
    """Комплексный кик гармоники: мнимая часть — нормальная, вещественная — скошенная."""
    cir, W, cols = train(element)
    v = features(cir, currents) @ W
    d = {(n, p): v[k] for k, (n, p) in enumerate(cols)}
    return complex(d.get((order, "скош."), 0.0), d.get((order, "норм."), 0.0))


def feature_breakdown(element):
    """Разбор потоковых признаков: (наконечники полюсов, шунты).

    Признак — поток узла, к которому подходит зазор кик-класса. Такие узлы двух родов: наконечники
    полюсов (имя на «p», выходят к пучку) и углы ярма, куда приходит воздушный шунт полюса. Число
    признаков есть их сумма; называть все наконечниками неверно — у DQ 4 полюса, а признаков 8
    (4 наконечника + 4 шунта)."""
    cir = bn.Circuit(element)
    ts, _ = tips(cir)
    poles = [t for t in ts if t.startswith("p")]
    return len(poles), len(ts) - len(poles)


def free_numbers(element):
    """Свободных чисел матрицы связи: (2 · число потоковых признаков) × число компонент кика."""
    cir = bn.Circuit(element)
    ts, _ = tips(cir)
    return 2 * len(ts) * len(columns(element)), len(ts), len(columns(element))


def accuracy(element):
    """Расхождение по кику на одиночных и комбинациях, по компонентам. Комбинации — проверка."""
    cir, W, cols = train(element)
    out = {}
    for n in ORDERS_OF[element]:
        maps = bn.read_maps(element, n)
        for grp, lbl in ((1, "одиночные"), (2, "комбинации")):
            rows = [(c, m) for c, m, _, _ in maps if (len(c) == 1) == (grp == 1)]
            num = sum(abs(predict(element, n, c) - complex(m)) for c, m in rows)
            den = sum(abs(complex(m)) for c, m in rows)
            out[(n, lbl)] = (100 * num / den if den else float("nan"), len(rows))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--acc", action="store_true")
    ap.add_argument("--free", action="store_true")
    a = ap.parse_args(argv)
    for el in ("DQ", "corrector"):
        if a.free:
            f, t, c = free_numbers(el)
            print(f"{el}: наконечников {t}, признаков {2 * t}, компонент {c}, свободных чисел {f}")
        if a.acc:
            acc = accuracy(el)
            for n in ORDERS_OF[el]:
                s, c = acc[(n, "одиночные")], acc[(n, "комбинации")]
                print(f"  {el:9} гарм {n}: одиночные {s[0]:5.1f}% ({s[1]})  "
                      f"комбинации {c[0]:5.1f}% ({c[1]})")
    if not (a.acc or a.free):
        ap.print_help()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
