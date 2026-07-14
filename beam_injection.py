# -*- coding: utf-8 -*-
"""Инжекция пучка в матричную структуру из ПРОИЗВОЛЬНОЙ плоскости (как в диффурной модели).

Пучок задаётся в particle_groups.csv (та же таблица, что для дифф. ядра SynchrotronSim). При выборе группы:
  1. ПРОВЕРКА: все частицы должны лежать в одной плоскости (⊥ средней скорости). Иначе — ошибка.
  2. Плоскость (точка+нормаль=направление скорости) встраивается в матричную структуру.
  3. Центр координат плоскости = пересечение прямой [выход ближайшего элемента ДО → вход ближайшего
     ПОСЛЕ] с плоскостью (опорная орбита в дрейфе). Ось локальной рамки = направление этой прямой.
  4. Промежутки до/после плоскости — модифицированные пустые дрейфы ⇒ ЗАМКНУТАЯ структура: плоскость
     старт И финиш; трек идёт полный круг и останавливается при возврате в плоскость.
  5. Направление обхода — по вектору скорости частиц (совпадает с опорной орбитой в дрейфе).

Глоб. состояние частицы → каноническое (x,pₓ,y,p_y,t,δ) в рамке плоскости: pₓ,p_y=(1+δ)·(dir·орт),
δ=(p−p₀)/p₀ по её скорости, x,y — поперечные смещения от центра.
"""
import os, csv, math
import numpy as np
import matrix_tracker as mt
import bip_matrix as bm

TB = mt.TB
ME, C = bm.ME, bm.C_LIGHT
P0 = bm.P_MOM                                    # опорный импульс γ₀·m·v₀


def _gamma(v):
    return 1.0 / math.sqrt(1.0 - (v / C) ** 2)


# ---------------- загрузка пучка ----------------
def load_group(group_id):
    """particle_groups.csv → список dict(pos[3], vel[3], m, q) для группы."""
    out = []
    for r in csv.DictReader(open(os.path.join(TB, "particle_groups.csv"), encoding="utf-8"), delimiter=";"):
        if r["group_id"] == str(group_id):
            out.append(dict(pos=np.array([float(r["x"]), float(r["y"]), float(r["z"])]),
                            vel=np.array([float(r["vx"]), float(r["vy"]), float(r["vz"])]),
                            m=float(r["m"]), q=float(r["q"])))
    if not out:
        raise ValueError(f"группа {group_id} не найдена в particle_groups.csv")
    return out


# ---------------- плоскость пучка (проверка компланарности) ----------------
def fit_plane(parts, tol=1e-4):
    """Плоскость пучка: центр=центроид позиций, нормаль=направление средней скорости. Проверка: все
    частицы в одной плоскости (продольный разброс < tol). Возврат (centroid, normal_dir)."""
    pos = np.array([p["pos"] for p in parts])
    vel = np.array([p["vel"] for p in parts])
    vmean = vel.mean(axis=0)
    n = vmean / np.linalg.norm(vmean)            # нормаль плоскости = направление пучка
    centroid = pos.mean(axis=0)
    long_spread = np.abs((pos - centroid) @ n)   # продольные смещения (вдоль нормали)
    if long_spread.max() > tol:
        raise ValueError(f"частицы НЕ компланарны: продольный разброс {long_spread.max()*1e3:.3f}мм > "
                         f"{tol*1e3:.3f}мм. Все частицы пучка должны быть в одной плоскости ⊥ скорости.")
    return centroid, n


# ---------------- псевдо-магнит «плоскость» ----------------
class PlaneMagnet:
    """Нулевой элемент-плоскость: рамка (axis=опорная орбита, uhat=гориз., vhat=верт.), entry=exit=center."""
    def __init__(self, center, axis):
        a = np.array(axis, float); a = a / np.linalg.norm(a)
        u = np.array([-a[1], a[0], 0.0])
        nu = np.linalg.norm(u); u = u / nu if nu > 1e-9 else np.array([0.0, 1.0, 0.0])
        self.center = np.array(center, float); self.axis = a; self.uhat = u
        self.vhat = np.cross(a, u)               # вертикаль (для планарного кольца ≈ +z)
        self.entry = self.center; self.exit = self.center
        self.L = self.Leff = 0.0; self.dev = "INJ_PLANE"; self.devs = ["INJ_PLANE"]
        self.yaw = math.atan2(a[1], a[0]); self.nsl = 1
        self.theta0 = self.k1 = self.k2 = self.k3 = 0.0; self.nonzero = False

    def body(self, st):                          # плоскость длины 0 — состояние не меняется
        return st


# ---------------- поиск места инжекции ----------------
def find_injection_site(mags, plane_pt, normal, vel_dir):
    """Найти дрейф (между эл-тами, включая шов), где опорная орбита пересекает плоскость. Возврат
    (i, center, axis): плоскость между mags[i] и mags[(i+1)%N]; center — пересечение; axis — направление
    орбиты в дрейфе (согласовано с направлением скорости)."""
    N = len(mags); best = None
    for i in range(N):
        A = mags[i].exit; B = mags[(i + 1) % N].entry   # дрейф (шов при i=N-1)
        AB = B - A; L = np.linalg.norm(AB)
        if L < 1e-9:
            continue
        axis = AB / L
        denom = float(normal @ AB)
        if abs(denom) < 1e-12:
            continue
        t = float(normal @ (plane_pt - A)) / denom       # параметр пересечения вдоль A→B
        if -0.02 <= t <= 1.02:
            center = A + t * AB
            if float(vel_dir @ axis) <= 0:               # скорость должна идти по орбите (вперёд)
                continue
            dist = np.linalg.norm(center - plane_pt)
            if best is None or dist < best[-1]:
                best = (i, center, axis, dist)
    if best is None:
        raise ValueError("не найден дрейф, пересекающий плоскость пучка по направлению скорости "
                         "(проверь, что пучок в пустом промежутке кольца и летит по орбите)")
    return best[0], best[1], best[2]


# ---------------- глоб. состояние → каноническое в рамке плоскости ----------------
def to_canonical(part, plane):
    rel = part["pos"] - plane.center
    x = float(rel @ plane.uhat); y = float(rel @ plane.vhat)
    v = np.linalg.norm(part["vel"]); d = part["vel"] / v
    delta = (_gamma(v) * part["m"] * v - P0) / P0
    px = (1.0 + delta) * float(d @ plane.uhat); py = (1.0 + delta) * float(d @ plane.vhat)
    return (x, px, y, py, 0.0, delta)


# ---------------- сборка инжектированной замкнутой структуры ----------------
def inject(lattice, group_id, source="physical"):
    """Полная инжекция: пучок группы → (плоскость-магнит, упорядоченное кольцо, канон. состояния).
    Возврат (plane, ordered_mags, states, info)."""
    eff, phys, coup, trans = mt.load_db()
    mags = mt.build_line(lattice, eff, phys, coup, trans, source=source)
    parts = load_group(group_id)
    centroid, normal = fit_plane(parts)
    vel_dir = np.mean([p["vel"] / np.linalg.norm(p["vel"]) for p in parts], axis=0)
    vel_dir /= np.linalg.norm(vel_dir)
    i, center, axis = find_injection_site(mags, centroid, normal, vel_dir)
    plane = PlaneMagnet(center, axis)
    ordered = mags[i + 1:] + mags[:i + 1]        # кольцо стартует ПОСЛЕ mags[i], возвращается к нему
    states = [to_canonical(p, plane) for p in parts]
    info = dict(nparts=len(parts), after_elem=mags[i].dev, before_elem=mags[(i + 1) % len(mags)].dev,
                center=center, axis=axis, closed=bool(np.linalg.norm(center - centroid) < 1e-3))
    return plane, ordered, states, info


# ---------------- трек одного оборота (плоскость → кольцо → плоскость) ----------------
def track_turn(plane, ordered, states):
    """Каждая частица: инжекция в плоскости → полный круг → возврат в плоскость (финиш). Возврат
    [(final_state, lost)]. lost=None если долетела."""
    res = []
    for st in states:
        s = tuple(st); lost = None
        try:
            s = mt.modified_drift(s, plane, ordered[0])          # плоскость → 1-й элемент
            s = ordered[0].body(s)
            for j in range(1, len(ordered)):
                s = mt.modified_drift(s, ordered[j - 1], ordered[j])
                s = ordered[j].body(s)
            s = mt.modified_drift(s, ordered[-1], plane)         # последний → плоскость = полный круг
        except mt.ParticleLost as e:
            lost = str(e)
        res.append((s, lost))
    return res


if __name__ == "__main__":
    import sys
    lat = int(sys.argv[1]) if len(sys.argv) > 1 else 101
    gid = int(sys.argv[2]) if len(sys.argv) > 2 else 1
    plane, ordered, states, info = inject(lat, gid)
    print(f"Группа {gid}, решётка {lat}: {info['nparts']} частиц, плоскость КОМПЛАНАРНА ✓")
    print(f"  инжекция в дрейфе между '{info['after_elem']}' и '{info['before_elem']}'")
    print(f"  центр плоскости = {info['center'].round(4)}  ось = {info['axis'].round(4)}")
    print(f"  замкнута (центр≈центроид пучка): {info['closed']}, магнитов в круге: {len(ordered)}")
    print("Полный оборот (плоскость → кольцо → плоскость):")
    surv = 0
    for k, (st, lost) in enumerate(track_turn(plane, ordered, states)):
        if lost:
            print(f"  частица {k}: ПОТЕРЯНА — {lost}")
        else:
            print(f"  частица {k}: круг пройден, финиш x={st[0]*1e3:+.4f}мм y={st[2]*1e3:+.4f}мм t={st[4]:.4f}м"); surv += 1
    print(f"ВЫЖИЛО {surv}/{info['nparts']} полный оборот")
