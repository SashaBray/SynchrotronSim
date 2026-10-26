# -*- coding: utf-8 -*-
"""Модель возмущений во времени для динамики стабилизации (SOFB/FOFB) — задаётся таблицей БД.

Каждая компонента = (имя, мишень, КАНАЛ, форма, корреляция, амплитуда, частота, фаза):
  КАНАЛ (quantity):
    position           — поперечные позиции магнитов (dx,dy) [м]; вызывающий код вносит feed-down киком.
    field:b1/k1/k2/k3  — РЯБЬ ПОЛЯ (имитация тока): ΔK = amplitude(отн.)·K_номин·осц. b1→kicks (изгиб),
                         k1→dks (градиент), k2→dk2s, k3→dk3s. field:skew→dsks (амплитуда АБСОЛЮТНАЯ, K_ном=1).
    У КАЖДОЙ пары (тип магнита × компонента) свой закон — разные обмотки/источники шумят по-разному.
  ФОРМА:   tone — A·sin(2πf·t+φ);   walk — случайная ломаная (кусочно-пост. скорость, рывок с частотой f).
  КОРРЕЛЯЦИЯ: independent — независимо по магнитам; girder — общая на опору (s-близость); ground — гладкая
    пространственная волна (сейсмика); common — синфазно по всем магнитам мишени (общий источник питания).
Амплитуда: position — метры; field — ОТНОСИТЕЛЬНАЯ (доля номинала, = доля ряби тока). Спек:
workspace/tables/disturbance_spec.csv (разделитель «;», точка — десятичная). См. load_spec.
"""
import csv, os
import numpy as np

TB = os.path.join("workspace", "tables")
_NUM = ("amplitude", "freq_hz", "phase_deg")
_FKEY = {"b1": "kicks", "k1": "dks", "skew": "dsks", "k2": "dk2s", "k3": "dk3s"}   # компонента → арг apply_perturb
_FNOM = {"b1": "theta0", "k1": "k1", "k2": "k2", "k3": "k3"}                        # компонента → номинал магнита


def load_spec(path=None):
    path = path or os.path.join(TB, "disturbance_spec.csv")
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
        if not (r.get("name") or "").strip():
            continue
        d = {k: (r.get(k) or "").strip() for k in ("name", "target", "plane", "form", "correlation", "quantity")}
        if not d["quantity"]:
            d["quantity"] = "position"                       # обратная совместимость
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
    """Модель возмущений. offsets(t,dt) → (dx,dy) [м] позиций; field(t,dt) → {арг apply_perturb: вектор ΔK}."""

    def __init__(self, spec, ring, seed=0):
        self.n = ring.n_mag
        kinds = np.array(ring.kinds)
        s_out = np.array(ring.s_out); Lel = np.array([float(getattr(m, "Leff", 0) or 0) for m in ring.mags])
        scenter = s_out - Lel / 2.0
        gid, ng = girder_groups(ring)
        nom = {c: np.array([float(getattr(ring.mags[i], a, 0.0) or 0.0) for i in range(self.n)])
               for c, a in _FNOM.items()}                    # номиналы b1/k1/k2/k3 по магнитам
        self.pos = []; self.fld = []
        for k, c in enumerate(spec):
            rng = np.random.default_rng(seed + 1000 * (k + 1))
            tgt = np.ones(self.n, bool) if c["target"] == "all" else (kinds == c["target"])
            # пространственный паттерн (RMS=1 на мишени)
            if c["correlation"] == "girder":
                p = rng.standard_normal(ng)[gid]
            elif c["correlation"] == "ground":
                p = _ground_wave(scenter, 4, rng)
            elif c["correlation"] == "common":
                p = np.ones(self.n)                          # синфазно (общий источник)
            else:
                p = rng.standard_normal(self.n)
            p = p * tgt
            rms = np.sqrt(np.mean(p[tgt] ** 2)) if tgt.any() else 0.0
            if rms > 0:
                p = p / rms
            q = c["quantity"]
            if q.startswith("field:"):                       # масштаб = номинал компоненты (skew → абс., ном=1)
                comp = q.split(":", 1)[1]
                scale = p * (nom[comp] if comp in nom else np.ones(self.n))
                apply_key = _FKEY.get(comp, "dks")
            else:
                scale = p; apply_key = None                  # позиция (метры)
            # единицы временной эволюции: walk independent → по магниту; walk girder → по опоре; иначе 1 (ког.)
            if c["form"] == "walk" and c["correlation"] == "independent":
                idx = np.where(tgt)[0]; B = np.zeros((self.n, len(idx)))
                for u, i in enumerate(idx):
                    B[i, u] = scale[i]
                nu = len(idx)
            elif c["form"] == "walk" and c["correlation"] == "girder":
                B = np.zeros((self.n, ng))
                for i in range(self.n):
                    if tgt[i]:
                        B[i, gid[i]] = scale[i]
                nu = ng
            else:
                B = scale.reshape(-1, 1); nu = 1
            comp_d = dict(name=c["name"], plane=c["plane"], form=c["form"], amp=c["amplitude"],
                          f=c["freq_hz"], ph=np.radians(c["phase_deg"]), B=B, nu=nu, apply_key=apply_key,
                          w=np.zeros(nu), v=(rng.normal(0, c["amplitude"], nu) if c["form"] == "walk" else None),
                          rng=rng)
            (self.fld if apply_key else self.pos).append(comp_d)

    def _u(self, c, t, dt):
        if c["form"] == "tone":
            return c["amp"] * np.sin(2 * np.pi * c["f"] * t + c["ph"]) * np.ones(c["nu"])
        c["w"] += c["v"] * dt                                # случайная ломаная
        reroll = c["rng"].random(c["nu"]) < (dt * c["f"])    # рывок: смена скорости
        if reroll.any():
            c["v"][reroll] = c["rng"].normal(0, c["amp"], int(reroll.sum()))
        return c["w"]

    def offsets(self, t, dt):
        dx = np.zeros(self.n); dy = np.zeros(self.n)
        for c in self.pos:
            off = c["B"] @ self._u(c, t, dt)
            if c["plane"] != "y":
                dx += off
            if c["plane"] != "x":
                dy += off
        return dx, dy

    def field(self, t, dt):
        """{арг apply_perturb: вектор ΔK по магнитам} для компонент поля (рябь тока)."""
        out = {}
        for c in self.fld:
            out[c["apply_key"]] = out.get(c["apply_key"], np.zeros(self.n)) + c["B"] @ self._u(c, t, dt)
        return out
