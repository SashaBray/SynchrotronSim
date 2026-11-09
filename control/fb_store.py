# -*- coding: utf-8 -*-
"""Писатель результатов эксперимента ОС (режим 4) В БД.

Идемпотентность по ключу `fb_exp_id` — калька `BPMRegistry.flush` (`matrix_tracker.py:731-743`): при
записи СОХРАНЯЕМ строки ЧУЖИХ экспериментов, переписываем заголовок + чужие + строки этого эксперимента.
Повторный прогон того же `fb_exp_id` не плодит строк. Десятичная — ТОЧКА (машинные таблицы, формат .9g).

Тяжёлые массивы (интерференционная карта, фазовый портрет) — сайдкар `np.savez` в
`workspace/results/fb/<id>/<signal>.npz` + строка-манифест в `fb_artifacts.csv` (паттерн `field_maps.csv`).
"""
import csv
import os

import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TB = os.path.join(ROOT, "workspace", "tables")
RES = os.path.join(ROOT, "workspace", "results", "fb")

RESULT_COLS = ["fb_exp_id", "signal", "step", "value"]
ARTIFACT_COLS = ["fb_exp_id", "signal", "npz_path", "shape", "meta"]


def _read_raw(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.reader(f, delimiter=";"))


def _fmt(v):
    return f"{float(v):.9g}"


def _clean(s):
    return str(s).replace(";", ",").replace("\n", " ").strip()


def write_results(fb_exp_id, rows):
    """rows = список (signal, step, value). Перезапись по ключу fb_exp_id в fb_results.csv."""
    path = os.path.join(TB, "fb_results.csv")
    keep = [r for r in _read_raw(path)[1:] if r and r[0] != str(fb_exp_id)]   # чужие эксперименты
    out = [[str(fb_exp_id), str(sig), str(int(step)), _fmt(val)] for sig, step, val in rows]
    os.makedirs(TB, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(RESULT_COLS); w.writerows(keep); w.writerows(out)
    return len(out)


def write_artifact(fb_exp_id, signal, arrays, meta=""):
    """arrays = dict имя→ndarray; сохраняем .npz сайдкар + строку в fb_artifacts.csv
    (перезапись строки по ключу fb_exp_id+signal)."""
    d = os.path.join(RES, str(fb_exp_id))
    os.makedirs(d, exist_ok=True)
    npz = os.path.join(d, f"{signal}.npz")
    np.savez(npz, **{k: np.asarray(v) for k, v in arrays.items()})
    rel = os.path.relpath(npz, ROOT).replace(os.sep, "/")
    shape = _clean(", ".join(f"{k}{tuple(np.asarray(v).shape)}" for k, v in arrays.items()))
    path = os.path.join(TB, "fb_artifacts.csv")
    keep = [r for r in _read_raw(path)[1:]
            if r and not (r[0] == str(fb_exp_id) and len(r) > 1 and r[1] == signal)]
    os.makedirs(TB, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(ARTIFACT_COLS); w.writerows(keep)
        w.writerow([str(fb_exp_id), str(signal), rel, shape, _clean(meta)])
    return rel


def set_status(fb_exp_id, status):
    """Обновить колонку status у строки fb_exp_id в fb_experiments.csv (прерывание-безопасность)."""
    path = os.path.join(TB, "fb_experiments.csv")
    if not os.path.exists(path):
        return
    with open(path, encoding="utf-8-sig") as f:
        rd = csv.DictReader(f, delimiter=";")
        fields = list(rd.fieldnames or [])
        rows = list(rd)
    if "status" not in fields:
        fields.append("status")
    for r in rows:
        if str(r.get("fb_exp_id")) == str(fb_exp_id):
            r["status"] = status
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";")
        w.writeheader(); w.writerows(rows)


def read_results(fb_exp_id):
    """Чтение результатов по ключу → dict signal → (steps[], values[]). Для отчёта."""
    path = os.path.join(TB, "fb_results.csv")
    out = {}
    for r in _read_raw(path)[1:]:
        if not r or r[0] != str(fb_exp_id):
            continue
        sig = r[1]; out.setdefault(sig, ([], []))
        out[sig][0].append(int(r[2])); out[sig][1].append(float(r[3]))
    return out


def read_artifacts(fb_exp_id):
    """Чтение манифеста артефактов по ключу → dict signal → npz-данные (np.load). Для отчёта."""
    path = os.path.join(TB, "fb_artifacts.csv")
    out = {}
    for r in _read_raw(path)[1:]:
        if not r or r[0] != str(fb_exp_id) or len(r) < 3:
            continue
        full = os.path.join(ROOT, r[2])
        if os.path.exists(full):
            out[r[1]] = dict(np.load(full))
    return out
