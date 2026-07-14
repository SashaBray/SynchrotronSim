# -*- coding: utf-8 -*-
"""Построение СЕКТОРНОГО ПРОФИЛЯ combined-function элемента (DQ) из полевых карт при заданных токах —
для секторного кика в производственном ядре (matrix_tracker). Рецепт из dq_sector.py (0.62мм ≡ диффур):
  1) суперпозиция полей обмоток: B(r) = Σ_w (I_w/nom_w)·map_w(r);
  2) опорная дуга y_ref(X) — трек осевой частицы этим полем (Boris);
  3) мультиполи ВДОЛЬ орбиты b0(X),Gu(X),Gv(X),S(X) в точке (X, y_ref(X), 0).
Возврат: dict(Xg, y_ref, b0, Gu, Gv, S, half). Кэш по (устройства, токи). Профиль линеен по полю, но
y_ref нелинеен ⇒ считаем на конкретный набор токов (для кольца — по одному разу на инстанс)."""
import os
import numpy as np
from field_quality import load
import _qp_scale as qs
import dq_sector as ds
from bip_matrix import QP

FB = os.path.join("workspace", "fields_bin")
FR = os.path.join("workspace", "fields_raw")
_CACHE = {}
_MAPCACHE = {}                                                 # загруженные карты по устройству (для Newton)


def _load_full(dev, hw_row):
    """Полная (зеркалированная) карта устройства + nominal_arg (кэш по устройству)."""
    if dev in _MAPCACHE:
        return _MAPCACHE[dev]
    mb = (hw_row.get("map_bin_path") or "-").strip()
    path = os.path.join(FB, mb) if mb not in ("-", "") and os.path.exists(os.path.join(FB, mb)) else None
    if path is None:
        for c in (os.path.join(FR, dev + ".fld"), os.path.join(FB, dev + ".bin")):
            if os.path.exists(c):
                path = c; break
    h, B = load(path)
    if h["mnx"] >= -1e-9:                                      # полукарта X≥0 → полная (зеркало)
        h, B = qs.mirror_full(h, B)
    nom = float(hw_row.get("nominal_arg") or 1.0) if (hw_row.get("nominal_arg") not in (None, "-", "")) else 1.0
    res = (h, np.nan_to_num(np.asarray(B, float)), nom); _MAPCACHE[dev] = res
    return res


def _field_fn(h, Bsup):
    """Функция поля pos→(Bx,By,Bz) трилинейной интерполяцией по суперпонированной карте (нуль вне)."""
    xg = h["mnx"] + np.arange(h["nx"]) * h["sx"]; yg = h["mny"] + np.arange(h["ny"]) * h["sy"]
    zg = h["mnz"] + np.arange(h["nz"]) * h["sz"]

    def f(pos):
        X, Y, Z = pos
        if X < xg[0] or X > xg[-1] or Y < yg[0] or Y > yg[-1] or Z < zg[0] or Z > zg[-1]:
            return np.zeros(3)
        i = min(int((X - xg[0]) / h["sx"]), h["nx"] - 2); j = min(int((Y - yg[0]) / h["sy"]), h["ny"] - 2)
        k = min(int((Z - zg[0]) / h["sz"]), h["nz"] - 2)
        tx = (X - xg[i]) / h["sx"]; ty = (Y - yg[j]) / h["sy"]; tz = (Z - zg[k]) / h["sz"]
        c = 0.0
        out = np.zeros(3)
        for di in (0, 1):
            for dj in (0, 1):
                for dk in (0, 1):
                    w = (tx if di else 1 - tx) * (ty if dj else 1 - ty) * (tz if dk else 1 - tz)
                    out = out + w * Bsup[i + di, j + dj, k + dk]
        return out
    return f


def _orbit(field, half, XSTART, ds_step, nstep):
    """Опорная дуга y_ref(X): трек осевой частицы (Y=Z=0, вдоль X) полем — простой Boris (монохром)."""
    x = XSTART; u = v = up = vp = 0.0; Xs = [x]; Ys = [u]
    for _ in range(nstep):
        ps = np.sqrt(max(1e-12, 1 - up * up - vp * vp))
        x += ds_step / 2; u += up * ds_step / (2 * ps); v += vp * ds_step / (2 * ps)
        if abs(x) <= half:
            B = field(np.array([x, u, v]))
            up += -QP * B[2] * ds_step; vp += QP * B[1] * ds_step
        ps = np.sqrt(max(1e-12, 1 - up * up - vp * vp))
        x += ds_step / 2; u += up * ds_step / (2 * ps); v += vp * ds_step / (2 * ps)
        Xs.append(x); Ys.append(u)
    return np.array(Xs), np.array(Ys)


def build(devs, currents, hw_rows, ds_step=1e-3):
    """Секторный профиль по (устройства, токи). Кэш. hw_rows: {dev: hardware_library-строка}."""
    key = (tuple(devs), tuple(round(c, 6) for c in currents))
    if key in _CACHE:
        return _CACHE[key]
    Bsup = None; h0 = None
    for dev, cur in zip(devs, currents):
        h, B, nom = _load_full(dev, hw_rows[dev])
        Bsup = (cur / nom) * B if Bsup is None else Bsup + (cur / nom) * B
        h0 = h
    half = float(h0["mxx"])                                    # полу-размах зоны поля (карта симм. −mxx..mxx)
    field = _field_fn(h0, Bsup)
    XSTART = -half - 0.05; nstep = int(round((2 * half + 0.1) / ds_step))
    Xo, Yo = _orbit(field, half, XSTART, ds_step, nstep)
    Xg, b0, Gu, Gv, S = ds.orbit_multipoles(field, Xo, Yo, half)
    prof = dict(Xg=Xg, y_ref=(Xo, Yo), b0=b0, Gu=Gu, Gv=Gv, S=S, half=half)
    _CACHE[key] = prof
    return prof
