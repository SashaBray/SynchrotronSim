# -*- coding: utf-8 -*-
"""Модель возмущений во времени для динамики стабилизации (SOFB/FOFB) — задаётся таблицей БД.

Каждая компонента = (мишень, плоскость, форма, корреляция, амплитуда, частота, фаза):
  ФОРМА:   tone  — A·sin(2πf·t+φ);   walk — случайная ломаная (кусочно-пост. скорость, рывок с частотой f).
  КОРРЕЛЯЦИЯ: independent — независимо по магнитам; girder — общая на опору (магниты по s-близости);
              ground — гладкая пространственная волна (низкие прострaн. частоты, корреляция по кольцу).
Возмущает поперечные позиции магнитов (dx,dy); вызывающий код вносит их в орбиту feed-down киком.
Спек: workspace/tables/disturbance_spec.csv (разделитель «;», точка — десятичная). См. load_spec.
"""
import csv, os
import numpy as np

TB = os.path.join("workspace", "tables")
_NUM = ("amplitude", "freq_hz", "phase_deg")


def load_spec(path=None):
    path = path or os.path.join(TB, "disturbance_spec.csv")
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
        if not (r.get("name") or "").strip():
            continue
        d = {k: (r.get(k) or "").strip() for k in ("name", "target", "plane", "form", "correlation")}
        d.update({k: float(r[k]) for k in _NUM})
        out.append(d)
    return out


def girder_groups(ring, max_gap=0.6):
    """Группировка магнитов в опоры/гирдеры по s-близости: зазор > max_gap начинает новую опору."""
    s_out = np.array(ring.s_out); Lel = np.array([float(getattr(m, "Leff", 0) or 0) for m in ring.mags])
    start = s_out - Lel; gid = np.zeros(len(start), int); g = 0
    for i in range(1, len(start)):
        if start[i] - s_out[i - 1] > max_gap:
            g += 1
        gid[i] = g
    return gid, g + 1


def _ground_wave(scenter, n_modes, rng):
    """Гладкая пространственная волна (сумма низких гармоник по s) — моделирует движение грунта."""
    L = scenter.max() - scenter.min() + 1e-9; p = np.zeros(len(scenter))
    for m in range(1, n_modes + 1):
        p += rng.standard_normal() * np.sin(2 * np.pi * m * (scenter - scenter.min()) / L + rng.uniform(0, 2 * np.pi))
    return p


class Disturbance:
    """Модель возмущений. offsets(t, dt) → (dx, dy) [м] по магнитам в момент t (walk продвигается на dt)."""

    def __init__(self, spec, ring, seed=0):
        self.n = ring.n_mag
        kinds = np.array(ring.kinds)
        s_out = np.array(ring.s_out); Lel = np.array([float(getattr(m, "Leff", 0) or 0) for m in ring.mags])
        scenter = s_out - Lel / 2.0
        gid, ng = girder_groups(ring)
        self.C = []
        for k, c in enumerate(spec):
            rng = np.random.default_rng(seed + 1000 * (k + 1))     # детерминированно по индексу компоненты
            tgt = np.ones(self.n, bool) if c["target"] == "all" else (kinds == c["target"])
            if c["correlation"] == "girder":
                p = rng.standard_normal(ng)[gid]
            elif c["correlation"] == "ground":
                p = _ground_wave(scenter, 4, rng)
            else:
                p = rng.standard_normal(self.n)
            p = p * tgt
            rms = np.sqrt(np.mean(p[tgt] ** 2)) if tgt.any() else 0.0
            if rms > 0:
                p = p / rms                                        # нормировка RMS=1 на мишени
            # единицы временной эволюции: walk independent → по магниту; walk girder → по опоре; иначе 1 (ког.)
            if c["form"] == "walk" and c["correlation"] == "independent":
                idx = np.where(tgt)[0]; B = np.zeros((self.n, len(idx)))
                for u, i in enumerate(idx):
                    B[i, u] = p[i]
                nu = len(idx)
            elif c["form"] == "walk" and c["correlation"] == "girder":
                B = np.zeros((self.n, ng))
                for i in range(self.n):
                    if tgt[i]:
                        B[i, gid[i]] = p[i]
                nu = ng
            else:
                B = p.reshape(-1, 1); nu = 1
            self.C.append(dict(name=c["name"], plane=c["plane"], form=c["form"], amp=c["amplitude"],
                               f=c["freq_hz"], ph=np.radians(c["phase_deg"]), B=B, nu=nu,
                               w=np.zeros(nu), v=(rng.normal(0, c["amplitude"], nu) if c["form"] == "walk" else None),
                               rng=rng))

    def offsets(self, t, dt):
        dx = np.zeros(self.n); dy = np.zeros(self.n)
        for c in self.C:
            if c["form"] == "tone":
                u = c["amp"] * np.sin(2 * np.pi * c["f"] * t + c["ph"]) * np.ones(c["nu"])
            else:
                c["w"] += c["v"] * dt                              # случайная ломаная
                reroll = c["rng"].random(c["nu"]) < (dt * c["f"])  # рывок: смена скорости
                if reroll.any():
                    c["v"][reroll] = c["rng"].normal(0, c["amp"], int(reroll.sum()))
                u = c["w"]
            off = c["B"] @ u
            if c["plane"] != "y":
                dx += off
            if c["plane"] != "x":
                dy += off
        return dx, dy
