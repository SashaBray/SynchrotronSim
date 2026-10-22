# -*- coding: utf-8 -*-
"""Пересборка LOCO-якобиана J = dM/dk1 для решётки 26 (воссоздан из diag_rm; билдера в репо не было).

M — вектор матрицы отклика «корректор→орбита в BPM» (concat Rxx,Ryy), на который опирается оптическая
(LOCO) ступень коррекции в свипах чувствительности. Выборка — ТА ЖЕ, что в прежнем якобиане (обратно
восстановлена из базлайна): корректоры = каждый 5-й кикер (24 из 120), BPM = каждый 3-й монитор (120 из
400), квадруполи = kind=quad (322). Нужна после правки геометрии диполей (слияние SBM сдвинуло индексы).

J[:,k] = (M(k1 квадрупола k сдвинут на δ) − M0)/δ. Пост-обработка в свипе: dk1 = −pinv(J)·(M_машины−M0)
(δ не влияет на результат: J∝1/δ, pinv∝δ). Резюмируемо: промежуточный J пишется на диск по чанкам.
Запуск: py -3.14 build_loco_jac.py   (NSL=8)
"""
import os, sys, pickle, time
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

LAT, NSL = 26, 8
TH = 1e-5            # пробный кик корректора [рад]
DELTA = 1e-3         # сдвиг k1 квадрупола для производной (величина не влияет на итог коррекции)
OUT = os.path.join("workspace", "results", f"loco_jac_exp{LAT}_nsl{NSL}.pkl")
SAVE_EVERY = 20

inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
ring = fr._ring(inf, nsl=NSL)
corr = oc.corrector_indices(ring)[::5]                 # 24 из 120 (каждый 5-й)
bpm = oi.monitor_rows(ring)[::3][:120]                  # 120 (каждый 3-й, первые 120 — как в базлайне)
quads = oi.quad_indices(ring, kinds=("quad",))         # 322
nrow = 2 * len(bpm) * len(corr)
print(f"LOCO-якобиан реш.{LAT} nsl{NSL}: corr={len(corr)} bpm={len(bpm)} quads={len(quads)} → J {nrow}×{len(quads)}", flush=True)


def rm(dks=()):
    """Вектор матрицы отклика (Rxx,Ryy) при опциональном сдвиге k1 (dks=[(quad,δ)])."""
    def ap(tk=None, tv=None):
        kk = {}; vv = {}
        if tk is not None: kk[tk] = TH
        if tv is not None: vv[tv] = TH
        ring.apply_perturb(dks=list(dks), kicks=list(kk.items()), vkicks=list(vv.items()))
    ap(); co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)): return None
    o = ring.track_elements([co])[0]; bx0, by0 = o[bpm, 0], o[bpm, 2]
    Rxx = np.zeros((len(bpm), len(corr))); Ryy = np.zeros((len(bpm), len(corr)))
    for j, c in enumerate(corr):
        ap(tk=c); co = ra.closed_orbit(ring, dims=4)
        if np.all(np.isfinite(co)): Rxx[:, j] = (ring.track_elements([co])[0][bpm, 0] - bx0) / TH
        ap(tv=c); co = ra.closed_orbit(ring, dims=4)
        if np.all(np.isfinite(co)): Ryy[:, j] = (ring.track_elements([co])[0][bpm, 2] - by0) / TH
    return np.concatenate([Rxx.ravel(), Ryy.ravel()])


# резюм
J = np.zeros((nrow, len(quads))); M0 = None; done = 0
if os.path.exists(OUT):
    try:
        prev = pickle.load(open(OUT, "rb"))
        if prev.get("J", np.empty(0)).shape == J.shape and list(prev.get("quads", [])) == list(quads) \
           and list(prev.get("corr", [])) == list(corr) and len(prev.get("bpm", [])) == len(bpm):
            J = prev["J"]; M0 = prev["M0"]; done = int(prev.get("done", 0))
            print(f"  резюм: уже готово {done}/{len(quads)} столбцов", flush=True)
    except Exception:
        pass

if M0 is None:
    M0 = rm()
    if M0 is None:
        print("[FAIL] M0 не посчитался (идеал неустойчив)"); sys.exit(1)


def save(complete):
    with open(OUT, "wb") as f:
        pickle.dump(dict(J=J, M0=M0, quads=list(quads), corr=list(corr),
                         bpm=np.asarray(bpm), done=done, complete=complete), f)


t0 = time.time()
for k in range(done, len(quads)):
    Mk = rm(dks=[(quads[k], DELTA)])
    J[:, k] = (Mk - M0) / DELTA if Mk is not None else 0.0
    done = k + 1
    if done % SAVE_EVERY == 0 or done == len(quads):
        save(done == len(quads))
        dt = time.time() - t0; rate = (done - 0) / dt if dt else 0
        print(f"  {done}/{len(quads)}  ({rate*60:.0f} столб/мин, ост.~{(len(quads)-done)/max(rate,1e-9)/60:.1f} мин)", flush=True)
ring.close()
print(f"[OK] {OUT}  ({done} столбцов за {(time.time()-t0)/60:.1f} мин)", flush=True)
