# -*- coding: utf-8 -*-
"""Экспорт траекторий МАТРИЧНОЙ модели в формат стандартного визуализатора.

Пишет results/trajectories_exp{ID}.bin (на частицу: uint32 N + N×3×float64 глоб. X,Y,Z) — ровно то,
что читает visualize.load_trajectories. Глобальные координаты берутся из тех же рамок магнитов
(lattice_configs x,y,z,yaw), что рисует визуализатор ⇒ траектория ложится в кольцо. Заодно заводит
строку эксперимента (exp_id→lattice_id=107) в experiments.csv, если её нет.

Запуск:  py -3.14 matrix_traj_export.py [exp_id=40] [lattice=107]
Просмотр: py -3.14 visualize.py 40
"""
import os, sys, csv, struct, shutil
import numpy as np
import matrix_tracker as mt
import matrix_optimizer as mo

RESULTS = os.path.join("workspace", "results")
TB = mt.TB


def _g(m, s_long, x, y):
    """Локальная (продольная s_long вдоль оси + поперечные x,y) → глобальная точка (в рамке магнита)."""
    P = m.entry + s_long * m.axis + x * m.uhat + y * m.vhat
    return (float(P[0]), float(P[1]), float(P[2]))


def body_global(m, st):
    """Тело магнита: (st_после, [глоб. точки вдоль тела]). Может бросить ParticleLost.
    БОКСОВАЯ модель (прямая ось + изгиб как кик) — для ГЛОБАЛЬНОГО пути .bin. Локальные координаты из неё
    НЕ берём: у диполей её x = отклонение от прямой оси (растёт с изгибом), а не бетатрон от опорной орбиты —
    для локальных используется секторная m.body() (см. local_trajectories)."""
    x, px, y, py, t, d = st
    if not m.nonzero:
        pts = [_g(m, 0.0, x, y)]
        st = mt._drift(st, m.L)
        pts.append(_g(m, m.L, st[0], st[2]))
        return st, pts
    N = m.nsl
    dsg = m.L / N; dsf = m.Leff / N
    th, k1, k1v, k2, k3, kE = m.theta0, m.k1, m.k1v, m.k2, m.k3, m.kEdge
    xr = m.xref                                                 # прогиб опорной дуги (см. Magnet._ref_bow)
    s = 0.0; pts = [_g(m, s, x, y)]; st = (x, px, y, py, t, d)
    if m.kEdge1:                                                # краевая линза на ВХОДЕ (гор.дефокус/верт.фокус)
        x, px, y, py, t, d = st
        u = x - m.xref_in if xr is not None else x
        st = (x, px + m.kEdge1 * u, y, py - m.kEdge1 * y, t, d)
    for i in range(N):
        st = mt._drift(st, dsg / 2.0); s += dsg / 2.0
        x, px, y, py, t, d = st
        u = x - xr[i] if xr is not None else x                  # отклонение от опорной дуги (сектор), не от оси
        px += th / N - k1 * u * dsf - 0.5 * k2 * (u * u - y * y) * dsf - (k3 / 6.0) * (u ** 3 - 3 * u * y * y) * dsf
        py += k1v * y * dsf + k2 * u * y * dsf + (k3 / 6.0) * (3 * u * u * y - y ** 3) * dsf - (kE / N) * y
        st = mt._drift((x, px, y, py, t, d), dsg / 2.0); s += dsg / 2.0
        pts.append(_g(m, s, st[0], st[2]))
    if m.kEdge2:                                                # краевая линза на ВЫХОДЕ
        x, px, y, py, t, d = st
        u = x - m.xref_out if xr is not None else x
        st = (x, px + m.kEdge2 * u, y, py - m.kEdge2 * y, t, d)
    return st, pts


def drift_global(st, A, B):
    """Модиф. дрейф A→B: (st_в_рамке_B, [глоб. точка входа B]). Прямой отрезок от выхода A. Может бросить."""
    st2 = mt.modified_drift(st, A, B)
    Q = B.entry + st2[0] * B.uhat + st2[2] * B.vhat
    return st2, [(float(Q[0]), float(Q[1]), float(Q[2]))]


def track_global(mags, particle):
    """Одна частица через кольцо → ndarray (M,3) глоб. точек; обрывается при потере."""
    pts = []; st = tuple(particle)
    try:
        st, pp = body_global(mags[0], st); pts += pp
        for j in range(1, len(mags)):
            st, pp = drift_global(st, mags[j - 1], mags[j]); pts += pp
            st, pp = body_global(mags[j], st); pts += pp
    except mt.ParticleLost:
        pass
    return np.asarray(pts, dtype=np.float64)


def _perp_dist(p, a, b):
    """Перпендикулярное расстояние точки p от прямой a→b (3D, на чистых float — быстро)."""
    (px, py, pz), (ax, ay, az), (bx, by, bz) = p, a, b
    abx, aby, abz = bx - ax, by - ay, bz - az
    L2 = abx * abx + aby * aby + abz * abz
    apx, apy, apz = px - ax, py - ay, pz - az
    if L2 < 1e-24:
        return (apx * apx + apy * apy + apz * apz) ** 0.5
    cx = aby * apz - abz * apy; cy = abz * apx - abx * apz; cz = abx * apy - aby * apx
    return ((cx * cx + cy * cy + cz * cz) / L2) ** 0.5


class Thinner:
    """Потоковое прореживание полилинии по кривизне (Reumann–Witkam): точка сохраняется, только если
    её удаление сдвинуло бы линию > eps. Прямые участки (дрейфы) → 2 точки; изгибы/бетатрон — сохранены.
    eps в метрах; eps<=0 ⇒ прореживания нет (все точки). Держит в памяти только сохранённые точки.

    add(p, pay) — необязательный payload на точку (напр. локальная запись); result_local() вернёт payload'ы
    ВЫРОВНЕННО с result() (те же сохранённые точки). Решение о прореживании — по ГЛОБ. геометрии p."""
    def __init__(self, eps):
        self.eps = eps; self.out = []; self.pay = []; self.anchor = None; self.prev = None
        self.prev_pay = None; self._fin = False
    def add(self, p, pay=None):
        p = (float(p[0]), float(p[1]), float(p[2]))
        if self.eps <= 0:
            self.out.append(p); self.pay.append(pay); return
        if self.anchor is None:
            self.out.append(p); self.pay.append(pay); self.anchor = p; self.prev = p; self.prev_pay = pay; return
        if _perp_dist(self.prev, self.anchor, p) > self.eps:
            self.out.append(self.prev); self.pay.append(self.prev_pay); self.anchor = self.prev
        self.prev = p; self.prev_pay = pay
    def _finalize(self):
        if self._fin:
            return
        self._fin = True
        if self.eps > 0 and self.prev is not None and (not self.out or self.out[-1] != self.prev):
            self.out.append(self.prev); self.pay.append(self.prev_pay)
    def result(self):
        self._finalize()
        return np.asarray(self.out, dtype=np.float64)
    def result_local(self):
        """Payload'ы сохранённых точек (M, K) — или None, если payload не передавали."""
        self._finalize()
        return np.asarray(self.pay, dtype=np.float64) if (self.pay and self.pay[0] is not None) else None


def deviation_paths(lattice, group, radius=1e-3, source="physical", settings=None):
    """Отклонение ПРОБНОГО пучка от ОПОРНОЙ орбиты — 3D-траектория (глоб.) + 2D-графики x(s), y(s).

    Пучок [опорная(0), +x, −x, +y, −y] инжектируется В ПЛОСКОСТИ ИНЖЕКЦИИ опыта (как реальный пучок,
    bi.inject), затем ДРЕЙФ до первого элемента («пустой пробег»), кольцо, у замкнутого — возврат к
    плоскости (замыкание). Опорная (0) = центроид группы ≈ проектная орбита ⇒ лежит на 0. На точку
    ndarray (M,6): (s, x, y, X, Y, Z) — локальные x,y (отклонение, 2D) + глобальные X,Y,Z (3D + связка).
    Трекинг m.body()+modified_drift (секторный при SS_SECTOR=1). Возврат (пути, метки, meta, mags_meta):
    meta=dict(name/idx/sp) точек, mags_meta=dict(center/name/idx/sp) магнитов. СП по длинным прямым (>2м).

    Опорная = центроид пучка на плоскости (для наших опытов ≈ проектная орбита ⇒ график опорной плоский);
    смещённые ±x/±y осциллируют бетатронно от плоскости инжекции."""
    import beam_injection as bi
    plane, ordered, _st, info = bi.inject(lattice, group, source=source, settings=settings)
    finish = info.get("finish", plane); closed = bool(info.get("closed", True))
    r = float(radius)
    beam = [(0, 0, 0, 0, 0, 0), (r, 0, 0, 0, 0, 0), (-r, 0, 0, 0, 0, 0),
            (0, 0, r, 0, 0, 0), (0, 0, -r, 0, 0, 0)]
    labels = ["опорная", "+x", "−x", "+y", "−y"]

    def gp(fr, x, y):                                                 # локальные (x,y) в рамке fr → глоб. XYZ
        p = fr[0] + x * fr[2] + y * fr[3]
        return float(p[0]), float(p[1]), float(p[2])

    gap_in = float((ordered[0].fr_in[0] - plane.fr_out[0]) @ plane.fr_out[1])          # плоскость → 1-й элемент
    ord_gaps = [float((ordered[j].fr_in[0] - ordered[j - 1].fr_out[0]) @ ordered[j - 1].fr_out[1])
                for j in range(1, len(ordered))]
    gap_out = float((finish.fr_in[0] - ordered[-1].fr_out[0]) @ ordered[-1].fr_out[1]) if closed else 0.0

    # МЕТАДАННЫЕ по точкам (общие): плоскость инжекции, затем на магнит (вход fr_in, выход fr_out), у
    # замкнутого — финиш=плоскость. СП по длинным прямым между ordered-элементами.
    frames = [plane.fr_out, ordered[0].fr_in, ordered[0].fr_out]
    name = ["плоскость инжекции", ordered[0].dev, ordered[0].dev]; idx = [0, 1, 1]; sp = [1, 1, 1]; cur = 1
    for j in range(1, len(ordered)):
        if ord_gaps[j - 1] > 2.0:
            cur += 1
        frames += [ordered[j].fr_in, ordered[j].fr_out]
        name += [ordered[j].dev, ordered[j].dev]; idx += [j + 1, j + 1]; sp += [cur, cur]
    if closed:
        frames.append(finish.fr_in); name.append("плоскость инжекции (замыкание)"); idx.append(0); sp.append(cur)
    meta = dict(name=name, idx=idx, sp=sp)
    # МЕТАДАННЫЕ МАГНИТОВ (по ordered-порядку) — для ховера/маркеров магнитов в 3D.
    mc = [tuple(float(v) for v in ordered[0].center)]; mn = [ordered[0].dev]; mi = [1]; ms = [1]; c2 = 1
    for j in range(1, len(ordered)):
        if ord_gaps[j - 1] > 2.0:
            c2 += 1
        mc.append(tuple(float(v) for v in ordered[j].center)); mn.append(ordered[j].dev); mi.append(j + 1); ms.append(c2)
    mags_meta = dict(center=mc, name=mn, idx=mi, sp=ms)

    out = []
    for st0 in beam:                                                  # на точку: (s, x, y, X, Y, Z) — локал.+глоб.
        st = tuple(map(float, st0)); s = 0.0
        g = gp(frames[0], st[0], st[2]); dev = [(0.0, st[0], st[2], *g)]
        try:
            st = mt.modified_drift(st, plane, ordered[0]); s += gap_in       # ВЕДУЩИЙ дрейф (пустой пробег)
            g = gp(frames[1], st[0], st[2]); dev.append((s, st[0], st[2], *g))
            st = ordered[0].body(st); s += ordered[0].L
            g = gp(frames[2], st[0], st[2]); dev.append((s, st[0], st[2], *g))
            for j in range(1, len(ordered)):
                st = mt.modified_drift(st, ordered[j - 1], ordered[j]); s += ord_gaps[j - 1]
                g = gp(frames[1 + 2 * j], st[0], st[2]); dev.append((s, st[0], st[2], *g))
                st = ordered[j].body(st); s += ordered[j].L
                g = gp(frames[2 + 2 * j], st[0], st[2]); dev.append((s, st[0], st[2], *g))
            if closed:                                                       # финиш → плоскость инжекции (замыкание)
                st = mt.modified_drift(st, ordered[-1], finish); s += gap_out
                g = gp(frames[-1], st[0], st[2]); dev.append((s, st[0], st[2], *g))
        except mt.ParticleLost:
            pass
        out.append(np.asarray(dev, float))
    return out, labels, meta, mags_meta


def experiment_closure(lattice, group, source="physical", settings=None):
    """Сводка замыкания опыта для подписи визуализации: инжектированный пучок, 1 оборот, СЕКТОРНЫМ телом
    (m.body()+modified_drift — как бит-в-бит). Возврат dict(closed, turns, n, survived, residual_mm).
    closed из инжекции (замкнутое кольцо vs открытый сегмент); residual_mm — макс. поперечн. невязка
    старт↔финиш за оборот (только для замкнутого)."""
    import beam_injection as bi
    plane, ordered, states, info = bi.inject(lattice, group, source=source, settings=settings)
    finish = info.get("finish", plane); closed = info.get("closed", True)
    res = []; surv = 0
    for st0 in states:
        st = tuple(st0); x0, y0 = st[0], st[2]
        try:
            st = mt.modified_drift(st, plane, ordered[0]); st = ordered[0].body(st)
            for j in range(1, len(ordered)):
                st = mt.modified_drift(st, ordered[j - 1], ordered[j]); st = ordered[j].body(st)
            st = mt.modified_drift(st, ordered[-1], finish)
            surv += 1
            if closed:
                res.append(((st[0] - x0) ** 2 + (st[2] - y0) ** 2) ** 0.5)
        except mt.ParticleLost:
            pass
    return dict(closed=closed, turns=1, n=len(states), survived=surv,
                residual_mm=(max(res) * 1e3 if res else None))


def local_trajectories(lattice, group, source="physical", settings=None, n_turns=1, target=30000,
                       states_override=None):
    """ЛОКАЛЬНЫЕ координаты СЕКТОРНОЙ модели (m.body при SS_SECTOR=1; x = отклонение от ОПОРНОЙ ОРБИТЫ,
    бит-в-бит elegant) вдоль трека инжектированного пучка. На точку (границы элементов, вход/выход):
    S (опорная арк-длина, м), x, x'≡px, y, y'≡py, ct, δ + индекс элемента; плюс глоб. XYZ (сверка/3D).
    states_override — если задан список 6-состояний (x,x',y,y',ct,δ) у плоскости инжекции, трекаются ОНИ
    вместо пучка группы (для произвольной частицы / сканов апертуры).

    n_turns>1: ПОЛНАЯ траектория за все n_turns оборотов (S копится: оборот k на S=k·кольцо..(k+1)·кольцо);
    равномерное прореживание шагом до ≤target точек/частицу (сохраняет бетатрон и FFT, ограничивает вес).
    surv_turns — сколько оборотов частица прошла до потери (=n_turns у выжившей).

    NB: боксовый body_global тут НЕ годится — у диполей его x = отклонение от ПРЯМОЙ оси (растёт с изгибом,
    ~сотни мм), а не бетатрон от опорной орбиты. Поэтому трекаем через m.body()+modified_drift (сектор).
    Требует SS_SECTOR=1 (иначе m.body() = боксовая модель — предупреждение в вызывающем коде)."""
    import beam_injection as bi
    plane, ordered, states, info = bi.inject(lattice, group, source=source, settings=settings)
    if states_override is not None:
        states = [tuple(map(float, s)) for s in states_override]   # произвольные частицы вместо пучка группы
    finish = info.get("finish", plane); closed = bool(info.get("closed", True))
    gap_in = float((ordered[0].fr_in[0] - plane.fr_out[0]) @ plane.fr_out[1])
    ord_gaps = [float((ordered[j].fr_in[0] - ordered[j - 1].fr_out[0]) @ ordered[j - 1].fr_out[1])
                for j in range(1, len(ordered))]
    gap_out = float((finish.fr_in[0] - ordered[-1].fr_out[0]) @ ordered[-1].fr_out[1]) if closed else 0.0
    entry_arc = [gap_in]                                          # опорная арк-длина входа каждого элемента
    for j in range(1, len(ordered)):
        entry_arc.append(entry_arc[-1] + ordered[j - 1].L + ord_gaps[j - 1])
    ring_len = entry_arc[-1] + ordered[-1].L + gap_out

    def gp(fr, x, y):                                            # локальные (x,y) в рамке fr → глоб. XYZ
        p = fr[0] + x * fr[2] + y * fr[3]
        return float(p[0]), float(p[1]), float(p[2])

    elem_names = ["плоскость инжекции"] + [m.dev for m in ordered]
    n_turns = max(1, int(n_turns))
    if n_turns > 1 and not closed:
        n_turns = 1                                              # открытый сегмент не зациклить
    locals_all = []; globals_all = []; surv_turns = []
    for st0 in states:
        st = tuple(map(float, st0)); done = 0
        loc = [(0.0, st[0], st[1], st[2], st[3], st[4], st[5], 0)]; glo = [gp(plane.fr_out, st[0], st[2])]
        try:
            for turn in range(n_turns):
                base = turn * ring_len
                st = mt.modified_drift(st, plane, ordered[0])
                loc.append((base + entry_arc[0], st[0], st[1], st[2], st[3], st[4], st[5], 1)); glo.append(gp(ordered[0].fr_in, st[0], st[2]))
                st = ordered[0].body(st)
                loc.append((base + entry_arc[0] + ordered[0].L, st[0], st[1], st[2], st[3], st[4], st[5], 1)); glo.append(gp(ordered[0].fr_out, st[0], st[2]))
                for j in range(1, len(ordered)):
                    st = mt.modified_drift(st, ordered[j - 1], ordered[j])
                    loc.append((base + entry_arc[j], st[0], st[1], st[2], st[3], st[4], st[5], j + 1)); glo.append(gp(ordered[j].fr_in, st[0], st[2]))
                    st = ordered[j].body(st)
                    loc.append((base + entry_arc[j] + ordered[j].L, st[0], st[1], st[2], st[3], st[4], st[5], j + 1)); glo.append(gp(ordered[j].fr_out, st[0], st[2]))
                if closed:                                       # финиш = плоскость (замыкание) → следующий оборот
                    st = mt.modified_drift(st, ordered[-1], finish)
                    loc.append((base + ring_len, st[0], st[1], st[2], st[3], st[4], st[5], 0)); glo.append(gp(finish.fr_in, st[0], st[2]))
                done = turn + 1
        except mt.ParticleLost:
            pass
        L = np.asarray(loc, float); G = np.asarray(glo, float)
        if target and len(L) > target:                          # равномерное прореживание (сохраняет форму и FFT)
            step = int(np.ceil(len(L) / target)); L = L[::step]; G = G[::step]
        locals_all.append(L); globals_all.append(G); surv_turns.append(done)
    return dict(mode="full", locals=locals_all, globals=globals_all, surv_turns=surv_turns, n_turns=n_turns,
                elem_names=elem_names, closed=closed, ring_len=ring_len, lattice=lattice, group=group)


def long_trajectories(lattice, group, source="physical", settings=None, n_turns=1000, states_override=None, rf=None):
    """ДЛИННЫЙ трек (много оборотов, 10к+): вместо полной траектории пишет состояние (x,x',y,y') у ПЛОСКОСТИ
    ИНЖЕКЦИИ РАЗ ЗА ОБОРОТ (turn-by-turn, ≤n_turns+1 точек/частицу) — для анализа огибающей/фазового портрета/
    затухания. surv_turns — сколько оборотов частица прошла (=n_turns у выжившей). Сектор m.body (SS_SECTOR=1).

    rf: ВЧ-резонатор (компенсация потерь SS_SR). None — выключен (δ дрейфует вниз от излучения). 'ideal' —
    ИДЕАЛЬНЫЙ ВЧ: восстанавливает δ до значения в начале оборота (точная компенсация потери каждой частицы,
    без синхротронных колебаний/энергоразброса). Именно ВЧ переводит падение канонического импульса в диполях
    (радиация) в затухание поперечного УГЛА — см. matrix_tracker._body_sector_symplectic."""
    import beam_injection as bi
    plane, ordered, states, info = bi.inject(lattice, group, source=source, settings=settings)
    if states_override is not None:
        states = [tuple(map(float, s)) for s in states_override]
    finish = info.get("finish", plane); closed = bool(info.get("closed", True))
    gap_in = float((ordered[0].fr_in[0] - plane.fr_out[0]) @ plane.fr_out[1])
    ord_gaps = [float((ordered[j].fr_in[0] - ordered[j - 1].fr_out[0]) @ ordered[j - 1].fr_out[1])
                for j in range(1, len(ordered))]
    gap_out = float((finish.fr_in[0] - ordered[-1].fr_out[0]) @ ordered[-1].fr_out[1]) if closed else 0.0
    ring_len = gap_in + sum(ord_gaps) + sum(m.L for m in ordered) + gap_out
    n_turns = max(1, int(n_turns))
    if not closed:
        n_turns = 1
    tt_all = []; surv_turns = []
    for st0 in states:
        st = tuple(map(float, st0)); done = 0
        tt = [(st[0], st[1], st[2], st[3])]                     # состояние у плоскости на старте
        try:
            for _turn in range(n_turns):                       # оборот: плоскость→кольцо→финиш(=плоскость)
                d0 = st[5]                                      # δ в начале оборота (для идеального ВЧ)
                st = mt.modified_drift(st, plane, ordered[0]); st = ordered[0].body(st)
                for j in range(1, len(ordered)):
                    st = mt.modified_drift(st, ordered[j - 1], ordered[j]); st = ordered[j].body(st)
                st = mt.modified_drift(st, ordered[-1], finish)
                if rf == "ideal":                              # ВЧ: вернуть энергию, потерянную за оборот на SR
                    st = (st[0], st[1], st[2], st[3], st[4], d0)
                tt.append((st[0], st[1], st[2], st[3])); done += 1
        except mt.ParticleLost:
            pass
        tt_all.append(np.asarray(tt, float)); surv_turns.append(done)
    return dict(mode="long", turnbyturn=tt_all, surv_turns=surv_turns, n_turns=n_turns,
                elem_names=["плоскость инжекции"] + [m.dev for m in ordered], closed=closed,
                ring_len=ring_len, lattice=lattice, group=group)


def write_bin(particles_pts, path):
    with open(path, "wb") as f:
        for pts in particles_pts:
            f.write(struct.pack("I", len(pts)))
            np.asarray(pts, np.float64).tofile(f)


# --- ЛОКАЛЬНЫЕ координаты матричного ядра (.npz) — для анализатора траекторий ---------------------
# На точку: S (опорная арк-длина, м) · x · x'≡px · y · y'≡py · ct · δ · elem (индекс элемента).
# x',y' — КАНОНИЧЕСКИЕ импульсы матмодели (≈ угол наклона; отличие ~x'^3/2). elem: 0=плоскость инжекции,
# j+1=ordered[j]. Хранится concat по частицам + counts (длины) ⇒ без ragged-массивов и без pickle.
LOCAL_COLS = ("S", "x", "xp", "y", "yp", "ct", "delta", "elem")


def write_local_npz(path, r):
    """r — dict из local_trajectories (mode='full': полная траектория locals/globals) или long_trajectories
    (mode='long': turn-by-turn (x,x',y,y') у плоскости на оборот). Concat по частицам + counts, без ragged/pickle."""
    mode = r.get("mode", "full")
    meta = dict(mode=str(mode), elem_names=np.asarray(r["elem_names"]),
                surv_turns=np.asarray(r["surv_turns"], np.int64), n_turns=int(r["n_turns"]),
                closed=bool(r["closed"]), ring_len=float(r["ring_len"]),
                lattice=int(r["lattice"]), group=int(r["group"]))
    if mode == "long":
        tt = r["turnbyturn"]; tt_counts = np.array([len(a) for a in tt], np.int64)
        T = np.concatenate(tt, axis=0) if tt else np.zeros((0, 4))
        np.savez(path, tt_counts=tt_counts, tt_x=T[:, 0], tt_xp=T[:, 1], tt_y=T[:, 2], tt_yp=T[:, 3], **meta)
    else:
        loc, glo = r["locals"], r["globals"]
        counts = np.array([len(a) for a in loc], np.int64)
        L = np.concatenate(loc, axis=0) if loc else np.zeros((0, 8))
        G = np.concatenate(glo, axis=0) if glo else np.zeros((0, 3))
        np.savez(path, counts=counts,
                 S=L[:, 0], x=L[:, 1], xp=L[:, 2], y=L[:, 3], yp=L[:, 4], ct=L[:, 5], delta=L[:, 6],
                 elem=L[:, 7].astype(np.int64), X=G[:, 0], Y=G[:, 1], Z=G[:, 2], **meta)


def read_local_npz(path):
    """Загрузка .npz → dict(mode, parts=[...], elem_names, n_turns, ring_len, closed, lattice, group).
    mode='full': parts с S,x,xp,y,yp,ct,delta,elem,X,Y,Z; mode='long': parts с tt_x,tt_xp,tt_y,tt_yp."""
    d = np.load(path, allow_pickle=False)
    mode = str(d["mode"]) if "mode" in d else "full"
    n_turns = int(d["n_turns"]); surv = d["surv_turns"]
    out = dict(mode=mode, elem_names=[str(s) for s in d["elem_names"]], closed=bool(d["closed"]),
               ring_len=float(d["ring_len"]), n_turns=n_turns, lattice=int(d["lattice"]), group=int(d["group"]))
    if mode == "long":
        tt_counts = d["tt_counts"]; toff = np.concatenate([[0], np.cumsum(tt_counts)]).astype(int)
        parts = []
        for i in range(len(tt_counts)):
            a, b = toff[i], toff[i + 1]
            parts.append(dict(tt_x=d["tt_x"][a:b], tt_xp=d["tt_xp"][a:b], tt_y=d["tt_y"][a:b], tt_yp=d["tt_yp"][a:b],
                              surv_turns=int(surv[i]), survived=bool(surv[i] >= n_turns)))
    else:
        counts = d["counts"]; off = np.concatenate([[0], np.cumsum(counts)]).astype(int)
        cols = ("S", "x", "xp", "y", "yp", "ct", "delta", "elem", "X", "Y", "Z")
        parts = []
        for i in range(len(counts)):
            a, b = off[i], off[i + 1]; p = {c: d[c][a:b] for c in cols}
            p["surv_turns"] = int(surv[i]); p["survived"] = bool(surv[i] >= n_turns)
            parts.append(p)
    out["parts"] = parts
    return out


def ensure_experiment(exp_id, lattice):
    """Гарантировать строку exp_id→lattice в experiments.csv (шаблон — существующий exp на этой решётке)."""
    path = os.path.join(TB, "experiments.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(rows[0].keys())
    if any(r["exp_id"] == str(exp_id) for r in rows):
        # обновим lattice на всякий
        for r in rows:
            if r["exp_id"] == str(exp_id):
                r["lattice_id"] = str(lattice)
    else:
        tmpl = next((r for r in rows if r["lattice_id"] == str(lattice)), rows[0]).copy()
        tmpl["exp_id"] = str(exp_id); tmpl["lattice_id"] = str(lattice)
        rows.append(tmpl)
        shutil.copy(path, os.path.join(TB, "backups", "experiments_pre_matrixtraj.csv")) \
            if os.path.isdir(os.path.join(TB, "backups")) else None
    rows.sort(key=lambda r: int(r["exp_id"]))
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows)


def main(exp_id=40, lattice=107, radius=2e-3):
    os.makedirs(RESULTS, exist_ok=True)
    mags = mo.build_ring(lattice)                       # из НАСТРОЕННЫХ mat_* (source=physical)
    # пучок под метки визуализатора: [опорная, +x, -x, +y(верт), -y] (в пределах динам. апертуры)
    beam = [(0, 0, 0, 0, 0, 0), (radius, 0, 0, 0, 0, 0), (-radius, 0, 0, 0, 0, 0),
            (0, 0, radius, 0, 0, 0), (0, 0, -radius, 0, 0, 0)]
    trajs = [track_global(mags, p) for p in beam]
    path = os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin")
    write_bin(trajs, path)
    ensure_experiment(exp_id, lattice)
    print(f"Матричная модель, кольцо {lattice} ({len(mags)} магнитов), пучок {len(beam)} частиц "
          f"(смещения ±{radius*1e3:.0f}мм):")
    for i, tr in enumerate(trajs):
        tag = ["опорная", "+x", "-x", "+y(верт)", "-y"][i]
        print(f"  P{i+1} {tag:9}: {len(tr):5} точек" + ("" if len(tr) > 5000 else "  (оборвана — потеря)"))
    print(f"Записано: {os.path.relpath(path)}")
    print(f"Строка эксперимента exp_id={exp_id}→lattice {lattice} в experiments.csv.")
    print(f"Просмотр:  py -3.14 visualize.py {exp_id}")


def main_injected(exp_id, lattice, group_id):
    """Экспорт траектории ИНЖЕКТИРОВАННОГО пучка (particle_groups) из его плоскости — полный оборот."""
    import beam_injection as bi
    os.makedirs(RESULTS, exist_ok=True)
    plane, ordered, states, info = bi.inject(lattice, group_id)
    finish = info.get("finish", plane)
    trajs = []
    for st in states:
        s = tuple(st); pts = [_g(plane, 0.0, st[0], st[2])]      # старт: точка на плоскости
        try:
            s, pp = drift_global(s, plane, ordered[0]); pts += pp
            s, pp = body_global(ordered[0], s); pts += pp
            for j in range(1, len(ordered)):
                s, pp = drift_global(s, ordered[j - 1], ordered[j]); pts += pp
                s, pp = body_global(ordered[j], s); pts += pp
            s, pp = drift_global(s, ordered[-1], finish); pts += pp   # → ФИНИШ-плоскость (замкнут: = старт)
        except mt.ParticleLost:
            pass
        trajs.append(np.asarray(pts, float))
    path = os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin")
    write_bin(trajs, path); ensure_experiment(exp_id, lattice)
    print(f"Инжекция группы {group_id} в решётку {lattice}: плоскость между "
          f"'{info['after_elem']}'→'{info['before_elem']}', центр {info['center'].round(3)}")
    for i, tr in enumerate(trajs):
        full = len(tr) > 5000
        print(f"  частица {i}: {len(tr):6} точек" + ("" if full else "  (оборвана — потеря)"))
    print(f"Записано: {os.path.relpath(path)}.  Просмотр: py -3.14 visualize.py {exp_id}")


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "inject":
        # py -3.14 matrix_traj_export.py inject <exp_id> <lattice> <group_id>
        main_injected(int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
    else:
        eid = int(sys.argv[1]) if len(sys.argv) > 1 else 40
        lat = int(sys.argv[2]) if len(sys.argv) > 2 else 107
        main(eid, lat)
