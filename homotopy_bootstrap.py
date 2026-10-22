# -*- coding: utf-8 -*-
"""ПУНКТ A статики: пробить СТЕНУ БУТСТРАПА гомотопией по амплитуде ошибки.

Штатная коррекция (oc.correct) при больших ошибках часто не находит замкнутую орбиту даже протяжкой
(стена ~20–30 мкм, ~19% машин). Гипотеза: стена алгоритмическая. Гомотопия: включаем ошибку по долям
f=0.5→0.75→1.0; на каждом шаге орбита ещё существует, корректируем и ПЕРЕНОСИМ уставки корректоров
тёплым стартом → доводим до полной ошибки. perturbed_workspace(seed, tol×f) = ТА ЖЕ машина ×f
(один RNG-поток, обрезка относительна σ).

Сравнение: доля восстановленных (замкнутая орбита при f=1) — штатно vs гомотопия, на tol×{2,3,4}.
Запуск: py -3.14 homotopy_bootstrap.py
"""
import os, sys, tempfile, shutil
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, errors as errmod, orbit_correct as oc, ring_analysis as ra

NSEED, CUTOFF = 15, 2.5
inf = dict(fr.info_of(26)); inf["_id"] = 26
ideal = fr._ring(inf, nsl=8); resp = oc.cached_response_matrix(ideal, 26); ideal.close()
corr, bpm = resp["corr"], resp["bpm"]
Px = np.linalg.pinv(resp["Rx"], rcond=1e-2); Py = np.linalg.pinv(resp["Ry"], rcond=1e-2)
TOL0 = errmod.load_tolerances()


def scale(tol, m):
    return {k: {c: v * m for c, v in d.items()} for k, d in tol.items()}


def refine(r, thx, thy, seed_co, iters=4):
    """SVD-доводка орбиты (абсолютные уставки thx,thy), посев Ньютона seed_co."""
    for _ in range(iters):
        r.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
        co = oc._closed_orbit_seeded(r, seed_co)
        if co is None:
            return None, thx, thy
        seed_co = co[:4]
        orb = r.track_elements([co])[0]
        thx = thx - Px @ orb[bpm, 0]; thy = thy - Py @ orb[bpm, 2]
    r.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
    return oc._closed_orbit_seeded(r, seed_co), thx, thy


def recover_baseline(seed, tol):
    root = tempfile.mkdtemp(prefix="bl_")
    try:
        errmod.perturbed_workspace(26, seed, tol, root, CUTOFF)
        r = fr._ring(inf, ws=root, nsl=8)
        try:
            res = oc.correct(r, resp, aperture=inf["loss_m"])
        except Exception:
            res = None
        r.close()
        return res is not None
    finally:
        shutil.rmtree(root, ignore_errors=True)


def recover_homotopy(seed, tol, fracs=(0.5, 0.75, 1.0)):
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr)); seed_co = np.zeros(4)
    root = tempfile.mkdtemp(prefix="ho_")
    try:
        for i, f in enumerate(fracs):
            ws = os.path.join(root, f"f{i}")
            errmod.perturbed_workspace(26, seed, scale(tol, f), ws, CUTOFF)
            r = fr._ring(inf, ws=ws, nsl=8)
            r.apply_perturb(kicks=list(zip(corr, thx.tolist())), vkicks=list(zip(corr, thy.tolist())))
            co = oc._closed_orbit_seeded(r, seed_co)            # тёплый старт
            if co is None:
                co = ra.closed_orbit(r, dims=4)                 # холодный из нуля
            if co is None:                                      # шаг велик → протяжка на этой доле
                tx, ty, co = oc._thread(r, resp, inf["loss_m"])
                if co is None:
                    r.close(); return False
                thx, thy = tx, ty
            co, thx, thy = refine(r, thx, thy, co[:4])
            if co is None:
                r.close(); return False
            seed_co = co[:4]; r.close()
        return True
    finally:
        shutil.rmtree(root, ignore_errors=True)


print(f"СТЕНА БУТСТРАПА: восстановление ({NSEED} машин), по-машинно. "
      f"warm-только = тёплый старт взял машину, которую штатно НЕ берётся.\n", flush=True)
print(f"{'tol×':>5} | {'штатно':>7} | {'warm':>6} | {'ОБЪЕДИН':>8} | {'warm-только':>11}", flush=True)
for m in (2, 3, 4):
    tol = scale(TOL0, m)
    bl = [recover_baseline(7000 + s, tol) for s in range(NSEED)]
    wm = [recover_homotopy(7000 + s, tol) for s in range(NSEED)]
    union = sum(b or w for b, w in zip(bl, wm))
    warmonly = sum(w and not b for b, w in zip(bl, wm))
    baseonly = sum(b and not w for b, w in zip(bl, wm))
    print(f"  ×{m:<2} | {sum(bl):3}/{NSEED} | {sum(wm):2}/{NSEED} | {union:3}/{NSEED} | "
          f"warm-только {warmonly}, штат-только {baseonly}", flush=True)
print("\n[готово]", flush=True)
