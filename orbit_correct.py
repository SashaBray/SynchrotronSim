# -*- coding: utf-8 -*-
"""Коррекция замкнутой орбиты корректорами (глобальная, по матрице отклика + SVD).

Погрешности установки (рассогласование магнитов, ошибки поля — см. errors.py) искажают замкнутую орбиту
всего кольца. Корректоры — это магниты-кикеры (в решётках СИЛА роль играют слоты SH1A/SH2B/SH3E, заведённые
как kind=kicker), которые задают локальные углы (гор. и верт.), чтобы выправить орбиту глобально.

Схема (стандартная для накопителей):
  1. Матрица отклика R (BPM × корректоры) — на ИДЕАЛЬНОЙ решётке, ОДИН раз: единичный кик каждого
     корректора → смещение орбиты во всех точках наблюдения. Отдельно для x (KICK) и y (VKICK).
  2. На каждую машину с ошибками: измерить искажённую орбиту o в точках наблюдения, решить θ = −R⁺·o
     (псевдообратная через SVD с обрезкой малых сингулярных чисел — гасит неэффективные комбинации
     корректоров, иначе они дают гигантские углы), применить кики, повторить 2–4 раза (отклик линейный
     лишь приближённо на искажённой орбите). Настройка СТАЦИОНАРНА для данной машины.

Точки наблюдения (BPM) — выходы всех магнитов (track_elements): least-squares минимизирует орбиту по
всему кольцу, а не в отдельных мониторах. Один сервер ядра переиспользуется через CppRing.apply_perturb
(RESET сохраняет ошибки в ws, сбрасывает лишь кики), так что ни матрица отклика, ни итерации не плодят
процессы.

БУТСТРАП ПРОТЯЖКОЙ (threading): при больших ошибках замкнутой орбиты нет вовсе (Ньютон с нуля расходится,
пучок теряется на первом же обороте). Тогда сначала выправляем ОДНООБОРОТНУЮ траекторию причинной
матрицей отклика T (кик влияет лишь ниже по ходу) — сектор за сектором, пока пучок не пройдёт оборот и не
появится замкнутая орбита, — и лишь затем переходим к штатному SVD по замкнутой орбите (см. _thread).
"""
import os
import pickle

import numpy as np

import ring_analysis as ra

RESP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workspace", "results")


def corrector_indices(ring):
    """Индексы корректоров (kind=kicker) среди магнитов кольца."""
    return [i for i, k in enumerate(ring.kinds) if k == "kicker"]


def _orbit(ring, bpm):
    """Замкнутая орбита в точках наблюдения: (ox, oy) или None, если машина неустойчива."""
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return None, None
    orb = ring.track_elements([co])[0]
    return co, orb


def response_matrix(ideal, corr=None, bpm=None, theta=1e-5):
    """Матрицы отклика на ИДЕАЛЬНОЙ решётке: Rx, Ry [м/рад] размера (BPM × корректоры).
    ideal мутируется (apply_perturb: RESET + единичный кик), в конце сбрасывается. bpm — индексы строк
    в выводе track_elements (по умолчанию выходы всех магнитов, строки 1..n_mag)."""
    corr = corrector_indices(ideal) if corr is None else corr
    bpm = np.arange(1, ideal.n_mag + 1) if bpm is None else np.asarray(bpm)
    ideal.apply_perturb()                                       # чистая решётка, без киков
    co0, base = _orbit(ideal, bpm)
    if co0 is None:
        raise RuntimeError("идеальная решётка неустойчива — матрица отклика не считается")
    bx0, by0 = base[bpm, 0], base[bpm, 2]
    Rx = np.zeros((len(bpm), len(corr)))
    Ry = np.zeros((len(bpm), len(corr)))
    for j, c in enumerate(corr):
        ideal.apply_perturb(kicks=[(c, theta)])
        co, orb = _orbit(ideal, bpm)
        if co is not None:
            Rx[:, j] = (orb[bpm, 0] - bx0) / theta
        ideal.apply_perturb(vkicks=[(c, theta)])
        co, orb = _orbit(ideal, bpm)
        if co is not None:
            Ry[:, j] = (orb[bpm, 2] - by0) / theta
    ideal.apply_perturb()                                       # вернуть в исходное состояние
    return dict(corr=list(corr), bpm=bpm, Rx=Rx, Ry=Ry)


def threading_response(ideal, corr, bpm, theta=1e-5):
    """Траекторные (причинные) матрицы отклика Tx, Ty [м/рад]: смещение ОДНООБОРОТНОЙ траектории
    частицы из 0 в точках наблюдения от единичного кика корректора. В отличие от замкнутой орбиты, кик
    влияет лишь на BPM НИЖЕ по ходу (без обхода кольца) — этим и пользуется протяжка (threading), когда
    замкнутой орбиты ещё нет. ideal мутируется и в конце сбрасывается."""
    ideal.apply_perturb()
    base = ideal.track_elements([np.zeros(6)])[0]
    bx0, by0 = base[bpm, 0], base[bpm, 2]
    Tx = np.zeros((len(bpm), len(corr)))
    Ty = np.zeros((len(bpm), len(corr)))
    for j, c in enumerate(corr):
        ideal.apply_perturb(kicks=[(c, theta)])
        tr = ideal.track_elements([np.zeros(6)])[0]
        Tx[:, j] = (tr[bpm, 0] - bx0) / theta
        ideal.apply_perturb(vkicks=[(c, theta)])
        tr = ideal.track_elements([np.zeros(6)])[0]
        Ty[:, j] = (tr[bpm, 2] - by0) / theta
    ideal.apply_perturb()
    return Tx, Ty


def cached_response_matrix(ideal, lattice, theta=1e-5, rebuild=False):
    """Матрица отклика с кэшем на диске (workspace/results/orbit_resp_exp{lattice}_nsl{nsl}.pkl).
    Кэш свой на каждое число срезов nsl (отклик слегка зависит от точности интегратора). Инвалидируется,
    если сменились геометрия (n_mag/периметр) или число корректоров — считается заново (≈2-3 мин)."""
    nsl = int(getattr(ideal, "nsl", 40))
    path = os.path.join(RESP_DIR, f"orbit_resp_exp{int(lattice)}_nsl{nsl}.pkl")
    corr = corrector_indices(ideal)
    if not rebuild and os.path.exists(path):
        try:
            with open(path, "rb") as f:
                r = pickle.load(f)
            if (r.get("n_mag") == ideal.n_mag and abs(r.get("circ", -1) - ideal.circumference) < 1e-6
                    and r.get("corr") == list(corr) and r.get("nsl") == nsl and "Tx" in r):
                return r
        except Exception:
            pass
    r = response_matrix(ideal, corr=corr, theta=theta)
    r["Tx"], r["Ty"] = threading_response(ideal, corr, r["bpm"], theta=theta)   # для протяжки
    r["n_mag"] = ideal.n_mag; r["circ"] = ideal.circumference; r["nsl"] = nsl
    os.makedirs(RESP_DIR, exist_ok=True)
    with open(path, "wb") as f:
        pickle.dump(r, f)
    return r


def _closed_orbit_seeded(ring, seed=None, eps=1e-8, tol=1e-12, max_iter=20):
    """Замкнутая орбита с ПОСЕВОМ (Ньютон из seed, а не из 0) — нужна после протяжки, когда орбита есть,
    но далека от нуля и Ньютон с нуля расходится. seed — 4-вектор (x,px,y,py). Возврат 6-вектор или None."""
    z = np.zeros(4) if seed is None else np.asarray(seed, float)[:4].copy()
    for _ in range(max_iter):
        pts = [np.r_[z, 0.0, 0.0]]
        for k in range(4):
            q = z.copy(); q[k] += eps; pts.append(np.r_[q, 0.0, 0.0])
        fin = ring.track(pts)[:, :4]
        if not np.all(np.isfinite(fin)):
            return None
        F = fin[0] - z
        if np.max(np.abs(F)) < tol:
            break
        J = ((fin[1:] - fin[0]) / eps).T
        try:
            z = z - np.linalg.solve(J - np.eye(4), F)
        except np.linalg.LinAlgError:
            return None
    co = np.r_[z, 0.0, 0.0]
    return co if np.all(np.isfinite(co)) else None


def _thread(ring, resp, aperture, passes=15, rcond=5e-2, relax=0.85):
    """Протяжка (threading): гасить корректорами ОДНООБОРОТНУЮ траекторию частицы из 0, пока не появится
    замкнутая орбита. Причинная матрица отклика T: решаем θ = −relax·T⁺·(траектория) по ДОСТИГНУТОЙ части
    (BPM в пределах апертуры) — там работают лишь корректоры выше по ходу (ниже нулевые столбцы, SVD их не
    трогает). Важно: гасим до МАЛОЙ траектории, а не «до апертуры» — большая орбита через секступоли даёт
    feed-down, тюн уезжает на резонанс и устойчивой орбиты нет. После каждого прохода пробуем замкнуть.
    Возврат (thx, thy, co) при успехе или (None, None, None). relax — недорелаксация для устойчивости
    итераций при больших ошибках; aperture — апертура потери [м]."""
    corr, bpm = resp["corr"], resp["bpm"]
    Tx, Ty = resp["Tx"], resp["Ty"]
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))
    for _ in range(passes):
        ring.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
        tr = ring.track_elements([np.zeros(6)])[0]
        x, y = tr[bpm, 0], tr[bpm, 2]
        within = np.isfinite(x) & np.isfinite(y) & (np.abs(x) <= aperture) & (np.abs(y) <= aperture)
        reach = int(np.argmin(within)) if not within.all() else len(bpm)
        if reach == 0:                                          # теряется сразу — корректоры не помогут
            return None, None, None
        rr = np.arange(reach)
        co = _closed_orbit_seeded(ring, np.array([np.mean(x[rr]), 0.0, np.mean(y[rr]), 0.0]))
        if co is not None:                                      # орбита замкнулась — протяжка окончена
            return thx, thy, co
        thx = thx - relax * (np.linalg.pinv(Tx[rr], rcond=rcond) @ x[rr])
        thy = thy - relax * (np.linalg.pinv(Ty[rr], rcond=rcond) @ y[rr])
    return None, None, None


def correct(ring, resp, iters=4, rcond=1e-2, aperture=0.012):
    """Скорректировать орбиту кольца ring (с ошибками в ws) корректорами по матрице отклика resp.
    ring мутируется: в конце на нём применены итоговые кики (готов для трекинга ДА/диффузии). Возврат
    dict(co, thx, thy, orb_before, orb_after, threaded) или None, если орбиту вернуть не удалось.

    БУТСТРАП: если замкнутой орбиты из 0 нет (сильные ошибки), сначала ПРОТЯЖКА (_thread) выправляет
    однооборотную траекторию, пока не появится замкнутая орбита, затем — штатный SVD по замкнутой орбите.
      thx, thy — углы корректоров [рад]; orb_before/after — (ox, oy) в точках наблюдения [м];
      threaded — понадобилась ли протяжка."""
    corr, bpm = resp["corr"], resp["bpm"]
    Px = np.linalg.pinv(resp["Rx"], rcond=rcond)
    Py = np.linalg.pinv(resp["Ry"], rcond=rcond)
    thx = np.zeros(len(corr))
    thy = np.zeros(len(corr))
    threaded = False
    ring.apply_perturb()                                        # чистое кольцо (только ошибки ws)
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):                             # замкнутой орбиты нет → протяжка
        if "Tx" not in resp:
            return None
        thx, thy, co = _thread(ring, resp, aperture)
        if co is None:
            return None                                         # орбиту вернуть не удалось
        threaded = True
    seed = co[:4]                                               # посев Ньютона в цикле — из найденной орбиты
    orb_before = None
    for _ in range(iters):                                      # штатный SVD по замкнутой орбите
        ring.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
        co = _closed_orbit_seeded(ring, seed)
        if co is None:
            return None
        seed = co[:4]                                           # следующий Ньютон — из текущей орбиты
        orb = ring.track_elements([co])[0]
        ox, oy = orb[bpm, 0], orb[bpm, 2]
        if orb_before is None:
            orb_before = (ox.copy(), oy.copy())                 # орбита до SVD (после протяжки, если была)
        thx = thx - Px @ ox
        thy = thy - Py @ oy
    ring.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
    co = _closed_orbit_seeded(ring, seed)
    if co is None:
        return None
    orb = ring.track_elements([co])[0]
    return dict(co=co, thx=thx, thy=thy, orb_before=orb_before,
                orb_after=(orb[bpm, 0], orb[bpm, 2]), threaded=threaded)


def rms(a):
    a = np.asarray(a, float)
    return float(np.sqrt(np.mean(a * a))) if a.size else float("nan")
