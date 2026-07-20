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
    """Тело магнита: (st_после, [глоб. точки вдоль тела]). Может бросить ParticleLost."""
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
    eps в метрах; eps<=0 ⇒ прореживания нет (все точки). Держит в памяти только сохранённые точки."""
    def __init__(self, eps):
        self.eps = eps; self.out = []; self.anchor = None; self.prev = None
    def add(self, p):
        p = (float(p[0]), float(p[1]), float(p[2]))
        if self.eps <= 0:
            self.out.append(p); return
        if self.anchor is None:
            self.out.append(p); self.anchor = p; self.prev = p; return
        if _perp_dist(self.prev, self.anchor, p) > self.eps:
            self.out.append(self.prev); self.anchor = self.prev
        self.prev = p
    def result(self):
        if self.eps > 0 and self.prev is not None and (not self.out or self.out[-1] != self.prev):
            self.out.append(self.prev)
        return np.asarray(self.out, dtype=np.float64)


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


def write_bin(particles_pts, path):
    with open(path, "wb") as f:
        for pts in particles_pts:
            f.write(struct.pack("I", len(pts)))
            np.asarray(pts, np.float64).tofile(f)


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
