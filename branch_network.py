# -*- coding: utf-8 -*-
"""Сеть ветвей с насыщением: подгонка на ОДИНОЧНЫХ обмотках → предсказание КОМБИНАЦИЙ.

Отличие от `saturation_model.py`: там каждая обмотка описывается своим скаляром и комбинации не
предсказываются вовсе. Здесь магнит — магнитная цепь целиком (граф из `saturation_branches.csv`,
включая зазорные рёбра), нелинейность только в железе, а «потоки в наконечниках → кик» — линейная
свёртка. Это проверка гипотезы о линейности внешней части модели.

ЧИСЛО ПАРАМЕТРОВ СНИЖЕНО СИММЕТРИЕЙ. Ветви одного класса (`sym_class`) делят длину пути; сечение
одно на весь магнит. Заметим сразу: при ОДНОМ сечении масштаб (S, проводимость зазора, вес кика)
вырожден — умножение S и проводимости на k с делением веса на k не меняет ничего. Поэтому S = 1
фиксировано, а проводимость зазора и вес кика меряются в этих единицах.

Расчёт цепи — узловой: неизвестные это магнитные потенциалы узлов U. Для ребра a→c с ампер-витками
своей катушки F: напор = (U_a − U_c) + F; железо даёт Φ = S·B(напор/ℓ) по кривой намагничивания,
зазор — Φ = P·напор (воздух линеен). В каждом узле сумма потоков нулевая.

    py -3.14 branch_network.py corrector      # подгонка на одиночных, проверка на комбинациях
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import types

import numpy as np
from scipy.optimize import least_squares

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import saturation_model as sm
from field_quality import load, _body_slice, _circle

R_PROBE, NTH = 8e-3, 128
REGISTRY = sm.REGISTRY


# ───────────────────────────── измерение по картам ─────────────────────────────

def harmonic_complex(path, n, r=R_PROBE, integrate=True):
    """КОМПЛЕКСНАЯ гармоника n карты. По умолчанию — ИНТЕГРАЛ вдоль пучка (это и есть кик).

    Гармоника в одном срезе смещена: эффективная длина поля разная у разных обмоток (у корректора
    160 мм у c1/c3/main против 176 мм у c2/c4 — они дальше от ярма), и сравнение по центральному
    сечению вносит систематику до 9%. Карта — половина магнита, поэтому интеграл удваивается."""
    h, B = load(path)
    if not integrate:
        by, bz = _circle(h, B, _body_slice(h, B), r, nth=NTH)
        sp = np.fft.fft(by + 1j * bz) / NTH
        return sp[0] if n == 1 else sp[-(n - 1)]
    xs = np.linspace(h["mnx"], h["mxx"], h["nx"])
    vals = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, r, nth=NTH)
        sp = np.fft.fft(by + 1j * bz) / NTH
        vals.append(sp[0] if n == 1 else sp[-(n - 1)])
    v = np.array(vals)
    return complex(2.0 * np.trapezoid(v.real, xs), 2.0 * np.trapezoid(v.imag, xs))


CACHE = os.path.join("workspace", "results", "saturation_model", "harmonics_cache.json")


def read_maps(element, order=1):
    """Карты элемента: [(словарь токов по обмоткам, знаковая гармоника, split, map_id)].

    Знак берётся проекцией на направление самой сильной карты — модуль гармоники его теряет."""
    import json
    cache = {}
    if os.path.exists(CACHE):
        try:
            cache = json.load(open(CACHE, encoding="utf-8"))
        except (OSError, ValueError):
            cache = {}
    rows = [r for r in sm.read_registry() if r["element"] == element]
    out, dirty = [], False
    for r in rows:
        p = os.path.join("workspace", "fields_bin", r["bin_path"])
        if not os.path.exists(p):
            continue
        cur = {}
        for part in r["excitation"].replace("+", ";").split(";"):
            part = part.strip()
            if "=" not in part:
                continue
            w, v = part.split("=", 1)
            try:
                cur[w.strip()] = float(v)
            except ValueError:
                pass
        if not cur:
            continue
        key = f"{r['map_id']}|{order}|int"
        if key in cache:
            h = complex(cache[key][0], cache[key][1])
        else:
            h = harmonic_complex(p, order)
            cache[key] = [h.real, h.imag]; dirty = True
        out.append([cur, h, r["split"], r["map_id"]])
    if dirty:
        os.makedirs(os.path.dirname(CACHE), exist_ok=True)
        json.dump(cache, open(CACHE, "w", encoding="utf-8"))
    return out                                              # гармоника остаётся КОМПЛЕКСНОЙ: кик
                                                            # двумерен (верт. и гориз. отклонение),
                                                            # проекция на одно направление теряет половину


# ───────────────────────────── магнитная цепь ─────────────────────────────

class Circuit:
    """Магнитная цепь элемента: граф из БД + материал. Параметры связаны классами симметрии."""

    def __init__(self, element):
        g = sm.BRANCH_GRAPH[element]
        self.el = element
        self.branches = list(g["branches"])
        self.edge, self.kind, self.sym = g["edge"], g["kind"], g["sym"]
        self.turns = g["turns"]
        self.coil_branches = g["coils"]                     # обмотка → список рёбер
        self.share = g["share"]                             # доля МДС обмотки на каждом её ребре
        self.nodes = []
        for b in self.branches:
            for n in self.edge[b]:
                if n not in self.nodes:
                    self.nodes.append(n)
        self.iron = [b for b in self.branches if self.kind[b] == "iron"]
        self.gaps = [b for b in self.branches if self.kind[b] == "gap"]
        self.len_classes = sorted({self.sym[b] for b in self.iron})
        self.gap_classes = sorted({self.sym[b] for b in self.gaps})
        self.mirror = g["mirror"]
        self.pos = g["pos"]                                 # поперечные положения концов ребра
        self._node_xy = {}                                  # координаты узлов из БД, если заданы
        # Координаты берутся у ВСЕХ узлов, у которых они заданы, а не только у наконечников.
        # Прежде фильтр по имени `p*` молча обнулял направление трубки у любой ветви, идущей не
        # между наконечниками: её вклад в кик выходил тождественно нулевым, и проверить гипотезу
        # «межполюсная утечка тоже создаёт поле» было нельзя — она давала ровно ноль по построению.
        for b_, ((ya, za), (yc, zc)) in (g.get("xy") or {}).items():
            a_, c_ = self.edge[b_]
            self._node_xy.setdefault(a_, (ya, za))
            self._node_xy.setdefault(c_, (yc, zc))
        self.mirror_classes = sorted({self.sym[b] for b in self.gaps if self.mirror[b] != 0})
        # поле у пучка создают только ветви МЕЖДУ НАКОНЕЧНИКАМИ; рассеяние через окно влияет на цепь,
        # но собственного вклада в кик не даёт (оно далеко от пучка)
        # Поле у пучка создают только ветви, идущие ЧЕРЕЗ АПЕРТУРУ и вдоль неё (между наконечниками).
        # Междуполюсная утечка идёт по бокам ног, в сотне миллиметров от пучка, и вклада не даёт;
        # рассеяние через окно — тем более. Иначе базис почти вырожден и МНК гасит главную обмотку.
        self.kick_classes = sorted({self.sym[b] for b in self.gaps
                                    if self.sym[b].startswith(("gap_", "tip_row"))})
        self.sense = {w: 1.0 for w in self.coil_branches}   # уточняется по данным
        # Если у ВСЕХ ветвей в БД проставлены сечение и длина — цепь задана геометрией, и подбирать
        # в ней нечего: свободной остаётся только матрица «потоки → кик» (метод МЭС в чистом виде).
        self.geom = None
        if all(b_ in g["geom"] for b_ in self.branches):
            mu0 = 4e-7 * np.pi
            lens, secs, perms = {}, {}, {}
            for b_ in self.branches:
                S, L = g["geom"][b_]
                if self.kind[b_] == "iron":
                    lens[self.sym[b_]], secs[self.sym[b_]] = L, S
                else:
                    perms[self.sym[b_]] = mu0 * S / L
            self.geom = (lens, perms, secs)
        bh = np.loadtxt(sm.MATERIAL, delimiter=";", skiprows=1)
        self.Harr, self.Barr = bh[:, 0], bh[:, 1]
        self.mu0 = 4e-7 * np.pi

    def dB_dH(self, H):
        """Наклон кривой намагничивания (для якобиана): кусочно-постоянный по таблице, μ0 за ней."""
        a = np.abs(H)
        slopes = np.diff(self.Barr) / np.maximum(np.diff(self.Harr), 1e-30)
        mid = 0.5 * (self.Harr[:-1] + self.Harr[1:])
        return np.where(a <= self.Harr[-1], np.interp(a, mid, slopes), self.mu0)

    # Физические координаты наконечников в поперечной плоскости (мм): y — шаг полюсов,
    # z — половина своей апертуры (средние полюса отодвинуты и потому дальше от пучка).
    PITCH_MM, HALF_GAP_MM = 30.0, {"o": 12.5, "c": 16.5}

    def tip_xy(self, node):
        """(y, z) наконечника в мм. None для узлов, не выходящих к апертуре.

        Если координаты заданы в БД (столбцы pos_*/z_*), берём их; иначе — правило по имени узла
        (так устроен корректор, где наконечники расположены регулярно)."""
        if node in self._node_xy:
            return self._node_xy[node]
        if len(node) < 3 or node[0] != "p":
            return None
        y = self.PITCH_MM * {"1": 1.0, "0": 0.0, "2": -1.0}[node[2]]
        z = self.HALF_GAP_MM["c" if node[2] == "0" else "o"] * (1.0 if node[1] == "T" else -1.0)
        return y, z

    def direction(self, b):
        """Единичное направление трубки потока ребра, как комплексное u_y + i·u_z.

        Поле, создаваемое потоком у пучка, смотрит ВДОЛЬ трубки: прямой зазор даёт чистую вертикаль,
        дальняя диагональ — почти горизонталь. Направление берётся из геометрии, не подбирается."""
        a_, c_ = self.edge[b]
        pa, pc = self.tip_xy(a_), self.tip_xy(c_)
        if pa is None or pc is None:
            return 0j
        dy, dz = pc[0] - pa[0], pc[1] - pa[1]
        n = math.hypot(dy, dz)
        return complex(dy / n, dz / n) if n else 0j

    def crosses(self, b):
        """Пересекает ли воздушная ветвь апертуру (верхний наконечник ↔ нижний)."""
        a_, c_ = self.edge[b]
        return (a_[:2] == "pT") != (c_[:2] == "pT")

    def row_sign(self, b):
        """Знак вклада потока ВДОЛЬ ряда: верхний ряд и нижний дают противоположный вклад."""
        a_, c_ = self.edge[b]
        return (1.0 if a_[:2] == "pT" else -1.0) * np.sign(self.pos[b][1] - self.pos[b][0])

    def B_of_H(self, H):
        """Индукция по напряжённости (кривая материала; за таблицей — по воздуху)."""
        a = np.abs(H)
        return np.sign(H) * np.where(a <= self.Harr[-1], np.interp(a, self.Harr, self.Barr),
                                     self.Barr[-1] + (a - self.Harr[-1]) * self.mu0)

    def mmf(self, currents):
        """Ампер-витки на каждом ребре, со знаком намотки.

        ГОТЧА: столбец `excitation` реестра УЖЕ в ампер-витках (не в амперах) — домножать на витки
        второй раз нельзя, у корректора это давало ошибку в 548 раз."""
        F = {b: 0.0 for b in self.branches}
        for w, brs in self.coil_branches.items():
            if w not in currents:
                continue
            f = currents[w] * self.sense.get(w, 1.0)
            for b in brs:
                F[b] += f * self.share.get((w, b), 1.0)
        return F

    def _prep(self):
        """Разовая подготовка векторов ветвей и инцидентной матрицы (иначе питоновский цикл на
        каждой итерации Ньютона — 240 мс на расчёт вместо единиц миллисекунд)."""
        if getattr(self, "_ready", False):
            return
        self.idx = {n: i for i, n in enumerate(self.nodes[1:])}   # первый узел опорный, U = 0
        nb, nn = len(self.branches), len(self.idx)
        self.ia = np.array([self.idx.get(self.edge[b][0], -1) for b in self.branches])
        self.ic = np.array([self.idx.get(self.edge[b][1], -1) for b in self.branches])
        self.is_gap = np.array([self.kind[b] == "gap" for b in self.branches])
        self.A = np.zeros((nn, nb))                              # узел × ветвь: поток из узла
        for j in range(nb):
            if self.ia[j] >= 0:
                self.A[self.ia[j], j] -= 1.0
            if self.ic[j] >= 0:
                self.A[self.ic[j], j] += 1.0
        self._ready = True

    def solve(self, currents, lens, perms, u0=None, secs=None):
        """Потоки ветвей при заданных токах — узловой метод, Ньютон с аналитическим якобианом.

        lens[класс] = длина пути по железу (м), perms[класс] = проводимость воздушной ветви."""
        self._prep()
        Fd = self.mmf(currents)                                  # МДС считаем ОДИН раз, не на ветвь
        F = np.array([Fd[b] for b in self.branches])
        Lb = np.array([lens[self.sym[b]] if self.kind[b] != "gap" else 1.0 for b in self.branches])
        Pb = np.array([perms[self.sym[b]] if self.kind[b] == "gap" else 0.0 for b in self.branches])
        Sb = np.array([(secs or {}).get(self.sym[b], 1.0) if self.kind[b] != "gap" else 1.0
                       for b in self.branches])
        nn = len(self.idx)

        def drops(u):
            ua = np.where(self.ia >= 0, u[np.maximum(self.ia, 0)], 0.0)
            uc = np.where(self.ic >= 0, u[np.maximum(self.ic, 0)], 0.0)
            return ua - uc + F

        def phi_of(d):
            return np.where(self.is_gap, Pb * d, Sb * self.B_of_H(d / Lb))

        def cond(d):                                             # ∂Φ/∂напор по ветвям
            return np.where(self.is_gap, Pb, Sb * self.dB_dH(d / Lb) / Lb)

        u = np.zeros(nn) if u0 is None else np.array(u0, float)
        # допуск ОТНОСИТЕЛЬНЫЙ: потоки порядка единиц, абсолютный 1e-14 недостижим и цикл крутится зря
        tol = 1e-12 * max(1.0, float(np.max(np.abs(phi_of(drops(u))))))
        for _ in range(40):
            d = drops(u)
            r = self.A @ phi_of(d)
            if np.max(np.abs(r)) < tol:
                break
            J = -(self.A @ (cond(d)[:, None] * self.A.T))   # знак: ∂напор/∂U = −Aᵀ
            try:
                step = np.linalg.solve(J + 1e-20 * np.eye(nn), -r)
            except np.linalg.LinAlgError:
                break
            t, r0 = 1.0, float(r @ r)
            for _ in range(40):
                rn = self.A @ phi_of(drops(u + t * step))
                if float(rn @ rn) < r0:
                    break
                t *= 0.5
            u = u + t * step
        return {b: float(v) for b, v in zip(self.branches, phi_of(drops(u)))}, u


# ───────────────────────────── подгонка ─────────────────────────────

# Сечения ЖЕЛЕЗА подбираются у всех классов, кроме одного: ярмо = 1 задаёт общий масштаб (умножение
# всех сечений и проводимостей на k с делением весов на k ничего не меняет — одну степень свободы
# обязательно закрепить, иначе подгонка бродит вдоль неё).
#
# У ВОЗДУХА длина и сечение по отдельности НЕразличимы: воздух линеен, Φ = μ0·(S/ℓ)·F, и они входят
# только частным. Поэтому проводимость зазора — это и есть «длина и сечение» зазора.
SEC_FIXED = "yoke"


def unpack(p, cir):
    """Нелинейные параметры: длины путей, проводимости зазоров и относительные сечения ног.

    Если цепь задана геометрией из БД — возвращаем её, вектор `p` не участвует."""
    if getattr(cir, "geom", None) is not None:
        lens, perms, secs = cir.geom
        return lens, perms, secs
    n1, n2 = len(cir.len_classes), len(cir.gap_classes)
    lens = {c: float(np.exp(p[i])) for i, c in enumerate(cir.len_classes)}
    perms = {c: float(np.exp(p[n1 + i])) for i, c in enumerate(cir.gap_classes)}
    free = [c for c in cir.len_classes if c != SEC_FIXED]
    secs = {c: 1.0 for c in cir.len_classes}
    for i, c in enumerate(free):
        secs[c] = float(np.exp(p[n1 + n2 + i]))
    return lens, perms, secs


def mirror_sign(cir, b):
    """Знак ветви при зеркальном отражении y → −y.

    Отражение переставляет наконечники «+y» и «−y»: ветвь (T1,B1) переходит в (T2,B2), (T1,B0) в
    (T2,B0) и так далее. Пара получает противоположные знаки, самосимметричная ветвь — ноль по сумме
    положений, и тогда знак берётся по разности (так (T1,B2) и (T2,B1) остаются парой)."""
    a, c = cir.pos[b]
    return float(np.sign(a + c)) if (a + c) != 0 else float(np.sign(a - c))


def tip_charges(cir, phi):
    """Поток, отдаваемый КАЖДЫМ наконечником полюса в воздух (магнитный «заряд» наконечника).

    Сумма по всем наконечникам тождественно нулевая — сколько потока вышло в зазор, столько и вошло."""
    q = {}
    for b in cir.gaps:
        a_, c_ = cir.edge[b]
        for n, sgn in ((a_, +1.0), (c_, -1.0)):
            if cir.tip_xy(n) is not None:
                q[n] = q.get(n, 0.0) + sgn * phi[b]
    return q


def multipole_basis(cir, currents, p, orders=(1,)):
    """Мультиполи ИЗ ГЕОМЕТРИИ: наконечник с потоком q в точке ζ даёт вклад q/ζⁿ в гармонику n.

    Это двумерное разложение поля точечного источника вокруг оси пучка. Никаких классов и свободных
    направлений: положение наконечника задаёт и величину, и фазу вклада во ВСЕ гармоники сразу.
    Свободным остаётся один масштаб на порядок — он вбирает в себя замену распределённого полюсного
    башмака точечным зарядом."""
    lens, perms, secs = unpack(p, cir)
    phi, _ = cir.solve(currents, lens, perms, secs=secs)
    q = tip_charges(cir, phi)
    out = []
    for n in orders:
        acc = 0j
        for node, qq in q.items():
            y, z = cir.tip_xy(node)
            zeta = complex(y, z) * 1e-3                     # координаты наконечников заданы в мм
            if abs(zeta) > 1e-9:
                acc += qq / zeta ** n
        out.append(acc)
    return np.array(out, dtype=complex)


def gap_basis(cir, currents, p):
    """Базис линейной части: на класс воздушных ветвей — ОДИН вещественный вес.

    Вклад ребра в кик = поток × направление его трубки (комплексное, из координат концов). Раньше
    направление задавалось двумя свободными весами на класс (симметричная и знакопеременная суммы);
    теперь оно геометрическое, а подбирается только величина. Это и строже, и вдвое дешевле."""
    lens, perms, secs = unpack(p, cir)
    phi, _ = cir.solve(currents, lens, perms, secs=secs)
    return np.array([sum(phi[b] * cir.direction(b) for b in cir.gaps if cir.sym[b] == c)
                     for c in cir.kick_classes], dtype=complex)


def solve_weights(cir, rows, p, wts=None):
    """Веса «поток → кик» — точным МНК при фиксированных нелинейных параметрах (веса линейны).

    Так оптимизатор ищет только 6 чисел вместо 9, и задача перестаёт быть плохо обусловленной."""
    A = np.array([gap_basis(cir, cur, p) for cur, _, _, _ in rows])
    y = np.array([meas for _, meas, _, _ in rows], dtype=complex)
    den = np.abs(y) + 0.05 * float(np.max(np.abs(y)))       # ОТНОСИТЕЛЬНАЯ метрика, как в resid
    if wts:                                                  # та же развесовка по обмоткам
        den = den / np.array([wts.get(tuple(sorted(c)), 1.0) for c, _, _, _ in rows])
    M = np.vstack([np.real(A) / den[:, None], np.imag(A) / den[:, None]])   # веса ВЕЩЕСТВЕННЫЕ
    b = np.concatenate([np.real(y) / den, np.imag(y) / den])
    # Регуляризация: базисные векторы почти коллинеарны, и без неё МНК находит гигантские
    # взаимно сокращающиеся веса (10^5), при которых главная обмотка предсказывается нулём.
    lam = 1e-3 * float(np.sqrt((M ** 2).sum() / M.shape[0]))
    M = np.vstack([M, lam * np.eye(M.shape[1])])
    b = np.concatenate([b, np.zeros(M.shape[1])])
    w, *_ = np.linalg.lstsq(M, b, rcond=None)
    return w


def predict(cir, currents, p, w):
    """Кик как комплексное число: вертикальная и горизонтальная составляющие сразу."""
    return complex(gap_basis(cir, currents, p) @ w)


def calibrate_senses(cir, single, p):
    """Знак намотки каждой обмотки — ИЗ ДАННЫХ: цепь с одной катушкой нечётна по её току, поэтому
    достаточно сравнить знак предсказания при sense=+1 со знаком измерения."""
    for wd in cir.coil_branches:
        rows = [t for t in single if list(t[0]) == [wd]]
        if not rows:
            continue
        t = max(rows, key=lambda t: abs(t[1]))
        cir.sense[wd] = 1.0
        base = gap_basis(cir, t[0], p)
        ref = complex(base.sum())                            # знак отклика при единичных весах
        if float(np.real(ref * np.conj(t[1]))) < 0:
            cir.sense[wd] = -1.0


def fit(element, order=1, verbose=True):
    cir = Circuit(element)
    maps = read_maps(element, order)
    if not maps:
        raise SystemExit(f"нет карт для {element}")
    single = [t for t in maps if len(t[0]) == 1]
    combo = [t for t in maps if len(t[0]) > 1]
    # Подгонка — по ВСЕМ одиночным обмоткам, проверка — по комбинациям, которые в подгонку не входят
    # вообще. Делить одиночные на train/test незачем: честность даёт именно отложенный набор комбинаций.
    train = single
    Fmax = max(abs(v) for t in maps for v in t[0].values())

    # Начальное приближение по смыслу: проводимость зазора такая, чтобы на САМОМ СИЛЬНОМ возбуждении
    # железо выходило на колено (B ≈ 1.7 Тл при S = 1), иначе модель стартует в линейном режиме и
    # оптимизатору нечем объяснить насыщение.
    n1, n2 = len(cir.len_classes), len(cir.gap_classes)
    nsec = len([c for c in cir.len_classes if c != SEC_FIXED])
    p0 = np.concatenate([np.log(np.full(n1, 0.2)), np.log(np.full(n2, 1.7 / max(Fmax, 1.0))),
                         np.zeros(nsec)])
    calibrate_senses(cir, single, p0)
    scale = max(abs(t[1]) for t in train)

    # Нормировка ОТНОСИТЕЛЬНАЯ, по каждой карте: при общем масштабе слабые карты (|C1| ~ 0.001 против
    # 0.256 у main) в оптимизации почти не участвовали, и их относительная ошибка уходила в сотни
    # процентов при формально низкой стоимости. Мягкий пол — чтобы не делить на почти ноль.
    # Вклад уравнивается ПО ОБМОТКАМ, а не по картам: у боковых по 4-6 карт, у main всего 5 из 17,
    # и при равном весе карт оптимизатор жертвовал главной обмоткой ради большинства.
    from collections import Counter
    cnt = Counter(tuple(sorted(t[0])) for t in train)
    wts = {k: 1.0 / np.sqrt(v) for k, v in cnt.items()}

    def resid(p):
        w = solve_weights(cir, train, p, wts)
        out = []
        for cur, meas, _, _ in train:
            k = wts[tuple(sorted(cur))]
            d = k * (predict(cir, cur, p, w) - meas) / (abs(meas) + 0.05 * scale)
            out += [float(np.real(d)), float(np.imag(d))]
        return np.array(out)

    # Границы физические: длина пути по железу 5 мм … 2 м, проводимость зазора — в пределах трёх
    # порядков от оценки. Без них оптимизатор уходит в переполнение (длины в тысячи км, нули).
    # ВАЖНО: сечение нормировано в единицу, поэтому параметр «длины» — это на самом деле ℓ/S.
    # Тонкое ярмо законно выглядит как очень длинный путь, и прежний потолок 2 м ЗАПРЕЩАЛ режим,
    # в котором боковые обмотки замыкаются через апертуру, а не через ярмо (именно он и нужен).
    lo = np.concatenate([np.full(n1, np.log(5e-3)), p0[n1:n1 + n2] - np.log(100.0),
                         np.full(nsec, np.log(0.02))])
    hi = np.concatenate([np.full(n1, np.log(500.0)), p0[n1:n1 + n2] + np.log(100.0),
                         np.full(nsec, np.log(50.0))])
    # Мультистарт: задача невыпуклая, с одной точки оптимизатор садится в локальный минимум
    # (длина упиралась в границу). Пробуем несколько разумных наборов длин и берём лучший.
    best = None
    for ell in ((0.25,) if cir.geom is None else ()):
        for kp in (0.3,):
            st = np.concatenate([np.full(n1, np.log(ell)), p0[n1:n1 + n2] + np.log(kp),
                                 np.zeros(nsec)])
            st = np.clip(st, lo + 1e-9, hi - 1e-9)
            try:
                cand = least_squares(resid, st, bounds=(lo, hi), x_scale="jac",
                                     xtol=1e-9, ftol=1e-9, max_nfev=120)
            except Exception:
                continue
            if verbose:
                print(f"    старт ℓ={ell:.2f} k={kp:.1f}: невязка {cand.cost:.5g}", flush=True)
            if best is None or cand.cost < best.cost:
                best = cand
    sol = best
    if cir.geom is not None:                                # геометрия задана — оптимизировать нечего
        sol = types.SimpleNamespace(x=p0, cost=float(np.sum(resid(p0) ** 2)))
    w = solve_weights(cir, train, sol.x)
    lens, perms, secs = unpack(sol.x, cir)
    if verbose:
        print(f"\n=== {element}: сеть ветвей, гармоника C{order} ===")
        print(f"  свободных чисел {len(sol.x) + len(w)} (нелинейных {len(sol.x)} + весов {len(w)}), "
              f"обучающих карт {len(train)}, комбинаций для проверки {len(combo)}")
        print("  длины путей по железу, мм: " +
              ", ".join(f"{c}={1000*lens[c]:.1f}" for c in cir.len_classes))
        print("  проводимость зазора: " + ", ".join(f"{c}={perms[c]:.4g}" for c in cir.gap_classes))
        print("  сечения (к ярму): " + ", ".join(f"{c}={secs[c]:.4g}" for c in cir.len_classes))
        print("  веса «поток → кик»: " + ", ".join(f"{v:+.4g}" for v in w))
        print("  знаки намотки: " + ", ".join(f"{k}{'+' if v > 0 else '−'}" for k, v in cir.sense.items()))

    def report(rows, title):
        if not rows:
            print(f"  {title}: карт нет"); return None
        errs = np.array([100.0 * abs(predict(cir, cur, sol.x, w) - meas) / (abs(meas) + 1e-30)
                         for cur, meas, _, _ in rows])
        print(f"  {title}: средняя {errs.mean():.1f}%, макс {errs.max():.1f}%  (n={len(errs)})")
        return errs

    report(train, "подгонка (одиночные обмотки)")
    report(combo, "ПРОВЕРКА на комбинациях")
    return cir, sol.x, w


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("element", nargs="?", default="corrector")
    ap.add_argument("--order", type=int, default=1, help="номер гармоники-наблюдаемой (1=диполь)")
    a = ap.parse_args(argv)
    fit(a.element, a.order)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
