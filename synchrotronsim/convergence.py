# -*- coding: utf-8 -*-
"""Сходимость срезов по типам устройств (шаг 6 docs/LIBRARY_ARCHITECTURE.md): сколько срезов сектор/кик-тела
нужно каждому типу, чтобы модель укладывалась в заданную точность, — и сколько это даёт скорости.

Эталон — всё кольцо на ref срезах (40 — точность паспорта). Для каждого типа устройства по отдельности его
срезы уменьшаются до кандидата n (на лету, `Ring.set_nsl`, остальные — эталон) и измеряются отклонения от
эталона: дробных тюнов (по матрице оборота вокруг нулевой орбиты идеального кольца) и нелинейной карты оборота
пробных частиц (|Δx|, |Δy| за оборот при амплитудах ~1 мм). Вклады типов складываются по модулю (консервативно),
и жадный подбор повышает срезы у того типа, где прирост точности на единицу цены (элементы × срезы) больше всего,
пока суммарная ошибка не уложится в допуск. Итоговый набор `nsl_dev` проверяется целиком.

    from synchrotronsim.convergence import scan, recommend, validate
    S = scan(26)                                   # ~1–2 мин на 16 потоках
    rec = recommend(S, tol_nu=1e-3, tol_map=1e-6)  # {device_type_id: срезов}
    validate(26, rec, S)                           # фактические ошибки и скорость

Это численная точность интегратора (схема Йошиды 4-го порядка, ошибка ∝ 1/N⁴), а не физика: менять дефолты
кокпита — решение пользователя (раздел «настройки точности»).
"""
from __future__ import annotations

import math
import queue
import threading
import time
from typing import Dict, Iterable, Optional

import numpy as np

REF_NSL = 40
CANDIDATES = (1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24)
# пробные частицы (x, px, y, py, t, δ): амплитуды ~1 мм и δ = 1e-3 — нелинейная карта оборота
PROBES = np.array([[1e-3, 0, 0, 0, 0, 0], [0, 0, 5e-4, 0, 0, 0], [1e-3, 0, 5e-4, 0, 0, 1e-3]], float)


def _ring(lattice, nsl, nsl_dev=None, **kw):
    from .engine import ExtRing
    from .db import open_db
    m = open_db().machine(lattice)
    mev = m.meta.design_energy_mev if m.meta else 6000.0
    return ExtRing(lattice=lattice, design_mev=mev, nsl=nsl, nsl_dev=nsl_dev or {}, machine=m, **kw)


def measure(ring):
    """(νx, νy дробные, карта оборота пробных частиц [3×4]) — вокруг нулевой орбиты идеального кольца."""
    import ring_analysis as ra
    _, M = ra.one_turn_matrix(ring, orbit=np.zeros(6))
    nux, nuy = ra.tunes_from_matrix(M, strict=False)
    pr = ring.track(PROBES, 1)[:, :4]
    return float(nux), float(nuy), pr


def slice_types(ring) -> Dict[str, list]:
    """Типы устройств с сектор/кик-телом (у остальных срезов нет) → индексы их элементов."""
    import ring_cpp
    d = ring_cpp.parse_describe(ring._r.describe_flat(), ring.n_mag)
    out: Dict[str, list] = {}
    for j, (nm, s, k) in enumerate(zip(ring.names, d["sector"], d["kick"])):
        if s or k:
            out.setdefault(nm, []).append(j)
    return out


def scan(lattice: int = 26, ref: int = REF_NSL, candidates: Iterable[int] = CANDIDATES, workers: Optional[int] = None,
         progress=print):
    """Ошибки по типам: {dev: {n: (|Δνx|, |Δνy|, max|Δкарты| м)}} + эталон, число элементов, вид."""
    base = _ring(lattice, ref)
    types = slice_types(base)
    kinds = {}
    for dev, idx in types.items():
        kinds[dev] = base.kinds[idx[0]]
    ref_nux, ref_nuy, ref_pr = measure(base)
    cand = [int(n) for n in candidates if int(n) < ref]
    tasks = [(dev, n) for dev in types for n in cand]
    nw = max(1, workers or 8)
    clones = [base.clone() for _ in range(nw)]
    todo = queue.Queue()
    for t in tasks:
        todo.put(t)
    res: Dict[str, Dict[int, tuple]] = {dev: {} for dev in types}
    lock = threading.Lock()
    done = [0]
    t0 = time.perf_counter()

    def work(w):
        r = clones[w]
        while True:
            try:
                dev, n = todo.get_nowait()
            except queue.Empty:
                return
            idx = types[dev]
            r._r.set_nsl(idx, n)
            try:
                nux, nuy, pr = measure(r)
                dmap = float(np.max(np.abs(pr - ref_pr)))
                e = (abs(nux - ref_nux), abs(nuy - ref_nuy), dmap)
                if not all(math.isfinite(v) for v in e):
                    e = (math.inf, math.inf, math.inf)
            except Exception:
                e = (math.inf, math.inf, math.inf)                 # частица потеряна / неустойчиво
            r._r.set_nsl(idx, ref)
            with lock:
                res[dev][n] = e
                done[0] += 1
                if progress and done[0] % 100 == 0:
                    progress(f"  {done[0]}/{len(tasks)} ({time.perf_counter() - t0:.0f} с)")

    th = [threading.Thread(target=work, args=(w,)) for w in range(nw)]
    [x.start() for x in th]
    [x.join() for x in th]
    for c in clones:
        c.close()
    return {"lattice": lattice, "ref": ref, "candidates": cand, "errors": res,
            "count": {d: len(i) for d, i in types.items()}, "kind": kinds, "indices": types,
            "ref_tunes": (ref_nux, ref_nuy), "seconds": time.perf_counter() - t0}


def _err(S, dev, n):
    if n >= S["ref"]:
        return 0.0, 0.0
    e = S["errors"][dev].get(n)
    if e is None:
        return math.inf, math.inf
    return max(e[0], e[1]), e[2]


def _envelope(S, dev):
    """Кандидаты и НЕВОЗРАСТАЮЩАЯ огибающая ошибок типа: e_env(n) = max_{m ≥ n} e(m). На очень грубых срезах
    кривые немонотонны (у QF8D карта на 2 срезах хуже, чем на 1), и решения по сырой кривой застревают."""
    cand = list(S["candidates"]) + [S["ref"]]
    vals = [_err(S, dev, n) for n in cand]
    en = [max(v[0] for v in vals[i:]) for i in range(len(cand))]
    em = [max(v[1] for v in vals[i:]) for i in range(len(cand))]
    return cand, en, em


def recommend(S, tol_nu: float = 1e-3, tol_map: float = 1e-6) -> Dict[str, int]:
    """Жадный подбор срезов: суммы |Δν| и |Δкарты| по типам (по огибающим) ≤ tol_nu, tol_map при наименьшей
    цене Σ (элементы × срезы). На каждом шаге — переход типа на ЛЮБОГО следующего кандидата с лучшим отношением
    «выигрыш / цена» (плато кривой не останавливает подбор). Возврат {device_type_id: срезов}."""
    env = {dev: _envelope(S, dev) for dev in S["errors"]}
    pos = {dev: 0 for dev in S["errors"]}
    cnt = S["count"]

    def sums():
        a = b = 0.0
        for dev, p in pos.items():
            _, en, em = env[dev]
            a += en[p]; b += em[p]
        return a, b

    while True:
        snu, smap = sums()
        if snu <= tol_nu and smap <= tol_map:
            break
        best, best_score = None, 0.0
        for dev, p in pos.items():
            cand, en, em = env[dev]
            for q in range(p + 1, len(cand)):
                gain = 0.0
                for e, s, tol in ((en, snu, tol_nu), (em, smap, tol_map)):
                    if s > tol:
                        if not math.isfinite(e[p]):
                            gain += 1e30 if math.isfinite(e[q]) else 0.0
                        else:
                            gain += (e[p] - e[q]) / tol
                score = gain / (cnt[dev] * (cand[q] - cand[p]))
                if score > best_score:
                    best, best_score = (dev, q), score
        if best is None:
            break                                                   # улучшать нечем (всё на эталоне/шум)
        pos[best[0]] = best[1]
    return {dev: int(env[dev][0][p]) for dev, p in pos.items()}


def turn_ms(ring, n=1, turns=4):
    """мс на оборот n частиц (минимум из 3 замеров)."""
    st = np.zeros((n, 6)); st[:, 0] = np.linspace(-1e-5, 1e-5, n) if n > 1 else 1e-5
    ring.track(st, 1)
    ts = []
    for _ in range(3):
        t = time.perf_counter(); ring.track(st, turns); ts.append((time.perf_counter() - t) / turns)
    return min(ts) * 1e3


def validate(lattice: int, nsl_dev: Dict[str, int], S=None, nsl: int = 8):
    """Набор срезов целиком: фактические |Δν|, |Δкарты| против эталона и время оборота одной частицы
    (без излучения и с СИ+ВЧ, как в кокпите). nsl — для типов вне nsl_dev."""
    ref = (S or {}).get("ref", REF_NSL)
    base = _ring(lattice, ref)
    rnux, rnuy, rpr = measure(base)
    r = _ring(lattice, nsl, nsl_dev)
    nux, nuy, pr = measure(r)
    out = {"dnux": abs(nux - rnux), "dnuy": abs(nuy - rnuy), "dmap": float(np.max(np.abs(pr - rpr))),
           "ms_turn": turn_ms(r), "ms_turn_ref": turn_ms(base)}
    from .db import open_db
    m = open_db().machine(lattice)
    if m.meta:
        rf = {"volt_mv": m.meta.rf_volt_mv, "harmonic": float(m.meta.rf_harmonic), "phase_deg": 165.0}
        rs = _ring(lattice, nsl, nsl_dev, sr=True, rf=rf)
        out["ms_turn_sr"] = turn_ms(rs)
        out["ms_turn_sr_12"] = turn_ms(rs, n=12)
    return out


def load_recommendation(lattice: int, which: str = "tight", root: Optional[str] = None) -> Dict[str, int]:
    """Набор срезов из workspace/results/nsl_convergence_<lat>.json (nsl_convergence.py): tight — самый строгий
    допуск, loose — самый мягкий; или путь к json (рекомендация или словарь {тип: срезов})."""
    import json
    import os
    if which not in ("tight", "loose"):
        with open(which, encoding="utf-8") as f:
            d = json.load(f)
        recs = d.get("recommendations") if isinstance(d, dict) else None
        if not recs:
            return {str(k): int(v) for k, v in d.items()}
    else:
        from .db import PROJECT_ROOT
        path = os.path.join(root or PROJECT_ROOT, "workspace", "results", f"nsl_convergence_{int(lattice)}.json")
        if not os.path.exists(path):
            raise FileNotFoundError(f"нет {path}: сначала py -3.14 nsl_convergence.py --lattice {lattice}")
        with open(path, encoding="utf-8") as f:
            recs = json.load(f)["recommendations"]
    recs = sorted(recs, key=lambda r: r["tol_nu"])
    r = recs[-1] if which == "loose" else recs[0]                 # путь к json с несколькими — самый строгий
    return {str(k): int(v) for k, v in r["nsl_dev"].items()}


# ---- набор срезов кокпита по умолчанию: workspace/configs/slices_<lat>.json (отслеживается git) ----------
def config_path(lattice: int, root: Optional[str] = None) -> str:
    import os
    from .db import PROJECT_ROOT
    return os.path.join(root or PROJECT_ROOT, "workspace", "configs", f"slices_{int(lattice)}.json")


def save_config(lattice: int, which: str = "tight", root: Optional[str] = None) -> str:
    """Записать выбранную рекомендацию (tight | loose) из nsl_convergence_<lat>.json в конфиг кокпита
    вместе с допусками, проверкой набора целиком и отпечатком машины (чтобы заметить смену решётки)."""
    import datetime
    import json
    import os
    from .db import PROJECT_ROOT, open_db
    src = os.path.join(root or PROJECT_ROOT, "workspace", "results", f"nsl_convergence_{int(lattice)}.json")
    with open(src, encoding="utf-8") as f:
        d = json.load(f)
    recs = sorted(d["recommendations"], key=lambda r: r["tol_nu"])
    r = recs[-1] if which == "loose" else recs[0]
    out = {"lattice": int(lattice), "preset": which, "ref": d.get("ref"), "tol_nu": r["tol_nu"],
           "tol_map_um": r["tol_map_um"], "validate": r["validate"],
           "fingerprint": open_db(root).machine(int(lattice)).fingerprint,
           "created": datetime.date.today().isoformat(), "source": "nsl_convergence.py", "nsl_dev": r["nsl_dev"]}
    path = config_path(lattice, root)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    return path


def cockpit_slices(lattice: int, root: Optional[str] = None):
    """(nsl_dev | None, пояснение) — набор срезов кокпита по умолчанию из workspace/configs/slices_<lat>.json.
    Нет конфига — None (как раньше: nsl + полюса ИД на 1 срезе). Отпечаток машины изменился — набор всё равно
    применяется (типы, которых в нём нет, идут на nsl), но пояснение предупреждает: перегнать nsl_convergence.py."""
    import json
    import os
    path = config_path(lattice, root)
    if not os.path.exists(path):
        return None, f"срезы: как раньше (нет {os.path.basename(path)})"
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    note = f"срезы: {cfg.get('preset', '?')} (допуск |Δν| ≤ {cfg.get('tol_nu')}, {len(cfg['nsl_dev'])} типов)"
    try:
        from .db import open_db
        if cfg.get("fingerprint") and open_db(root).machine(int(lattice)).fingerprint != cfg["fingerprint"]:
            note += " — решётка изменилась после подбора: перегнать nsl_convergence.py"
    except Exception:
        pass
    return {str(k): int(v) for k, v in cfg["nsl_dev"].items()}, note
