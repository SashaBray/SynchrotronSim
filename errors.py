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
    for r in rows[1:]:
        if not r or r[li] != str(lattice):
            continue
        kind = r[ki].strip()
        t = tol.get(kind)
        if not t:
            continue
        # геометрия: локальная рамка по НОМИНАЛЬНЫМ yaw/pitch/roll, сдвиг в глобальные x/y/z
        yaw, pitch, roll = (float(r[ci[c]] or 0) for c in ("yaw", "pitch", "roll"))
        axis, uhat, vhat = mt.frame(yaw, pitch, roll)
        dx = _trunc_normal(rng, t.get("dx", 0), cutoff)
        dy = _trunc_normal(rng, t.get("dy", 0), cutoff)
        dz = _trunc_normal(rng, t.get("dz", 0), cutoff)
        shift = dx * np.asarray(uhat) + dy * np.asarray(vhat) + dz * np.asarray(axis)
        for a, s in zip(("x", "y", "z"), shift):
            r[ci[a]] = repr(float(r[ci[a]] or 0) + float(s))
        for ang in ("yaw", "pitch", "roll"):
            d = _trunc_normal(rng, t.get("d" + ang, 0), cutoff)
            if d:
                r[ci[ang]] = repr(float(r[ci[ang]] or 0) + d)
        # поле: один множитель возбуждения на магнит
        fpct = t.get("field_pct", 0)
        if fpct:
            gain = 1.0 + _trunc_normal(rng, fpct / 100.0, cutoff)
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
