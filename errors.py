# -*- coding: utf-8 -*-
"""Ошибки установки: рассогласование положения/ориентации магнитов и ошибки силы поля.

Ядро (C++ и Python) уже трекает решётку по геометрии x/y/z/yaw/pitch/roll (Tait–Bryan, pitch/roll —
установочные наклоны) и полям mat_* из lattice_configs. Значит смоделировать ошибки = сгенерировать
ВОЗМУЩЁННУЮ lattice_configs и стравить её ядру. Здесь:
  • load_tolerances — допуски (СКО) по типу магнита из workspace/tables/element_errors.csv;
  • perturbed_workspace — на заданном seed рисует гауссовы отклонения на КАЖДЫЙ экземпляр, сдвигает
    его в ЛОКАЛЬНОЙ рамке (u,v,axis), поворачивает, масштабирует поля, и пишет отдельную рабочую папку
    (tables/lattice_configs.csv + копии element_connectivity/particle_groups), путь которой передаётся
    в ring_cpp.CppRing(ws=...). Реальная БД не трогается.

Допуски задаются СКО (амплитуда): сдвиги dx/dy/dz [м], повороты dyaw/dpitch/droll [рад], ошибка поля
field_pct [%]. Распределение — гаусс с обрезкой по cutoff·σ. Ошибка поля — ОДИН множитель (1+δ) на
магнит (ошибка возбуждения: все мультиполи экземпляра меняются согласованно).
"""
import csv
import os
import shutil

import numpy as np

# СЕКТОРНАЯ модель обязательна (дуга диполей, бит-в-бит elegant/C++). matrix_tracker читает SS_SECTOR на
# импорте, поэтому выставляем ДО импорта — иначе включится боксовая модель (хорда) и геометрия разойдётся
# с C++-ядром (страж «геометрия разошлась»). То же делает ring_analysis.
os.environ.setdefault("SS_SECTOR", "1")
import matrix_tracker as mt                                   # frame(yaw,pitch,roll) — единый источник рамки

TB = os.path.join("workspace", "tables")
TOL_COLS = ("dx", "dy", "dz", "dyaw", "dpitch", "droll", "field_pct")
# колонки геометрии/полей в lattice_configs, которые возмущаем
_GEO = ("x", "y", "z", "yaw", "pitch", "roll")
_FIELD = ("mat_B", "mat_G", "mat_S", "mat_O", "mat_Gv")


def load_tolerances(path=None):
    """Допуски по типу магнита: {kind: {dx, dy, dz, dyaw, dpitch, droll, field_pct}} (float, СИ/%)."""
    path = path or os.path.join(TB, "element_errors.csv")
    out = {}
    with open(path, encoding="utf-8-sig") as f:
        for r in csv.DictReader(f, delimiter=";"):
            k = (r.get("kind") or "").strip()
            if not k:
                continue
            out[k] = {c: float(r[c]) for c in TOL_COLS if r.get(c) not in (None, "")}
    return out


def _trunc_normal(rng, sigma, cutoff):
    """Одно гауссово число с СКО sigma, обрезанное на ±cutoff·sigma (0 при sigma=0)."""
    if not sigma:
        return 0.0
    while True:
        v = rng.normal(0.0, sigma)
        if abs(v) <= cutoff * sigma:
            return v


def _rotate(vec, axis, ang):
    """Поворот вектора vec вокруг единичной оси axis на угол ang (Родригес)."""
    import math
    vec = np.asarray(vec, float); axis = np.asarray(axis, float)
    c, s = math.cos(ang), math.sin(ang)
    return vec * c + np.cross(axis, vec) * s + axis * float(np.dot(axis, vec)) * (1.0 - c)


def _angles_from_frame(axis, uhat):
    """Глобальные (yaw, pitch, roll), согласованные с mt.frame, из ортонормальной тройки.
    Ось зависит только от yaw/pitch; roll берём по развороту uhat относительно бескреновой тройки."""
    import math
    yaw = math.atan2(axis[1], axis[0])
    pitch = -math.asin(max(-1.0, min(1.0, axis[2])))
    _, u2, v2 = mt.frame(yaw, pitch, 0.0)                       # бескреновая тройка тех же yaw/pitch
    roll = math.atan2(float(np.dot(uhat, v2)), float(np.dot(uhat, u2)))
    return yaw, pitch, roll


def perturbed_workspace(lattice, seed, tol, out_dir, cutoff=2.5):
    """Записать возмущённую рабочую папку для одной случайной машины (seed). Возврат пути (для ws=...).

    tol — из load_tolerances. Для каждого магнита решётки с типом из tol: сдвиг в локальной рамке
    (dx·u+dy·v+dz·axis), добавка к yaw/pitch/roll, множитель поля (1+δ). Прочие строки — без изменений."""
    rng = np.random.default_rng(seed)
    src_lc = os.path.join(TB, "lattice_configs.csv")
    with open(src_lc, encoding="utf-8-sig", newline="") as f:
        rows = list(csv.reader(f, delimiter=";"))
    hdr = rows[0]; ci = {c: i for i, c in enumerate(hdr)}
    li = ci["lattice_id"]; ki = ci["kind"]
    truth = []                                                 # ИСТИНА по каждому магниту (в порядке кольца) — только для валидации
    for r in rows[1:]:
        if not r or r[li] != str(lattice):
            continue
        kind = r[ki].strip()
        rec = dict(kind=kind, dx=0.0, dy=0.0, dz=0.0, dyaw=0.0, dpitch=0.0, droll=0.0, field=0.0)
        truth.append(rec)                                      # запись на КАЖДЫЙ магнит решётки (нули, если без допуска)
        t = tol.get(kind)
        if not t:
            continue
        # НОМИНАЛЬНАЯ орбитная рамка элемента (ось — вдоль пучка, uhat/vhat — поперечные)
        yaw, pitch, roll = (float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))
        axis, uhat, vhat = mt.frame(yaw, pitch, roll)
        # СДВИГ в локальной рамке: dx·uhat (попереч. гор.) + dy·vhat (попереч. верт.) + dz·axis (продольн.)
        dx = _trunc_normal(rng, t.get("dx", 0), cutoff)
        dy = _trunc_normal(rng, t.get("dy", 0), cutoff)
        dz = _trunc_normal(rng, t.get("dz", 0), cutoff)
        rec["dx"], rec["dy"], rec["dz"] = dx, dy, dz
        shift = dx * np.asarray(uhat) + dy * np.asarray(vhat) + dz * np.asarray(axis)
        for a, s in zip(("x", "y", "z"), shift):
            r[ci[a]] = repr(float(r[ci[a]] or 0) + float(s))
        # ОРИЕНТАЦИЯ: повороты в ЛОКАЛЬНОЙ рамке (рыскание вокруг vhat, тангаж вокруг uhat, крен вокруг оси),
        # скомпонованы с номинальной тройкой → извлекаем новые глобальные yaw/pitch/roll для симулятора.
        dyaw = _trunc_normal(rng, t.get("dyaw", 0), cutoff)
        dpitch = _trunc_normal(rng, t.get("dpitch", 0), cutoff)
        droll = _trunc_normal(rng, t.get("droll", 0), cutoff)
        rec["dyaw"], rec["dpitch"], rec["droll"] = dyaw, dpitch, droll
        if dyaw or dpitch or droll:
            ax, uu, vv = np.asarray(axis, float), np.asarray(uhat, float), np.asarray(vhat, float)
            if dyaw:
                ax, uu = _rotate(ax, vv, dyaw), _rotate(uu, vv, dyaw)
            if dpitch:
                ax, vv = _rotate(ax, uu, dpitch), _rotate(vv, uu, dpitch)
            if droll:
                uu, vv = _rotate(uu, ax, droll), _rotate(vv, ax, droll)
            ny, npi, nro = _angles_from_frame(ax, uu)
            r[ci["yaw"]] = repr(ny); r[ci["pitch"]] = repr(npi); r[ci["roll"]] = repr(nro)
        # поле: один множитель возбуждения на магнит
        fpct = t.get("field_pct", 0)
        if fpct:
            delta = _trunc_normal(rng, fpct / 100.0, cutoff)
            rec["field"] = delta
            gain = 1.0 + delta
            for c in _FIELD:
                if c in ci and r[ci[c]] not in (None, "", "0"):
                    r[ci[c]] = repr(float(r[ci[c]]) * gain)
    # запись возмущённой рабочей папки: возмущённый lattice_configs + копии двух таблиц
    dst_tb = os.path.join(out_dir, "tables")
    os.makedirs(dst_tb, exist_ok=True)
    with open(os.path.join(dst_tb, "lattice_configs.csv"), "w", encoding="utf-8", newline="") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    for t2 in ("element_connectivity.csv", "particle_groups.csv"):
        shutil.copyfile(os.path.join(TB, t2), os.path.join(dst_tb, t2))
    with open(os.path.join(dst_tb, "errors_truth.csv"), "w", encoding="utf-8", newline="") as f:   # эталон для валидации
        w = csv.writer(f, delimiter=";")
        w.writerow(["idx", "kind", "dx", "dy", "dz", "dyaw", "dpitch", "droll", "field"])
        for i, rec in enumerate(truth):
            w.writerow([i, rec["kind"], repr(rec["dx"]), repr(rec["dy"]), repr(rec["dz"]),
                        repr(rec["dyaw"]), repr(rec["dpitch"]), repr(rec["droll"]), repr(rec["field"])])
    return out_dir


def error_summary(tol):
    """Короткое человекочитаемое описание заданных допусков (для подписи/шапки отчёта)."""
    parts = []
    for k, t in tol.items():
        if any(t.get(c) for c in TOL_COLS):
            xy = t.get("dx", 0) * 1e6
            parts.append(f"{k}: сдвиг {xy:.0f} мкм, крен {t.get('droll', 0)*1e3:.2f} мрад, "
                         f"поле {t.get('field_pct', 0):.2f} %")
    return "; ".join(parts)
