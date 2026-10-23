# -*- coding: utf-8 -*-
"""ЕДИНЫЙ МОДУЛЬ НАСТРОЙКИ (слепая коррекция) СИЛА — один источник истины вместо инлайн-дублей в свипах.

Процесс (всё слепое, по наблюдаемому):
  1. ОРБИТА  — oc.correct (SVD по матрице отклика + протяжка-бутстрап). Робастный линейный пол.
  2. LOCO    — градиентная коррекция квадруполей по матрице отклика (dks). Главный вклад в ДА.
  3. ОКТУПОЛЬ— транзитный скан силы k3 (dk3s) вдоль Парето ДА⊥↔MA. По умолчанию g=0 (проект);
               подтверждённый оптимум g=+0.25 (октуполи +25%) доминирует проект по ОБЕИМ осям
               (см. [[static-correction-closeout]], pareto_octupole.py).

API:
  ctx = load_loco(lattice, nsl)            # LOCO-контекст (якобиан, Pl, наборы corr/bpm/quad)
  oct = octupole_context(ideal_ring)       # индексы октуполей и k3₀ (ТОЛЬКО с идеального кольца!)
  res = tune(ring, resp, aperture, loco=ctx, oct_ctx=oct, octupole_g=0.0)
        → dict(ok, co, bk, bv, bd, dk3, stages). ring остаётся со всеми уставками (готов к трекингу).

Гомотопия-бутстрап (union протяжка∪гомотопия для «стенных» машин) — отдельный решатель замыкания,
см. homotopy_bootstrap.py; применяется ПЕРЕД tune(), когда oc.correct не замкнул орбиту.
"""
import os, pickle
import numpy as np
import orbit_correct as oc
import ring_analysis as ra

RESP_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "workspace", "results")
TH = 1e-5                    # пробный кик корректора для LOCO-отклика [рад]
LOCO_RCOND = 1e-2
OCTUPOLE_PARETO_G = 0.25     # подтверждённый Парето-оптимум (октуполи +25%), не включён по умолчанию


def load_loco(lattice=26, nsl=8, rcond=LOCO_RCOND):
    """LOCO-контекст из кэша loco_jac_exp{lat}_nsl{nsl}.pkl (J, M0, quads, corr, bpm) + Pl=pinv(J)."""
    path = os.path.join(RESP_DIR, f"loco_jac_exp{int(lattice)}_nsl{int(nsl)}.pkl")
    JD = pickle.load(open(path, "rb"))
    return dict(Pl=np.linalg.pinv(JD["J"], rcond=rcond), M0=JD["M0"],
                quads=JD["quads"], corr=JD["corr"], bpm=np.asarray(JD["bpm"]))


def octupole_context(ideal_ring):
    """Индексы октуполей и их проектная k3₀. ТОЛЬКО на идеальном кольце: на возмущённом r.mags падает
    стражем «геометрия разошлась»."""
    octs = [i for i, k in enumerate(ideal_ring.kinds) if k == "oct"]
    k30 = [float(ideal_ring.mags[i].k3) for i in octs]
    return dict(octs=octs, k30=k30)


def _loco_response(ring, loco, bk, bv):
    """Вектор матрицы отклика «корректор→орбита в BPM» (concat Rxx,Ryy) при базовых киках bk,bv.
    ВАЖНО: тестовый кик СУММИРУЕТСЯ к базовому через словарь (серверный KICK перезаписывает, не +=)."""
    corr, bpm = loco["corr"], loco["bpm"]
    bkd, bvd = dict(bk), dict(bv)

    def ap(tk=None, tv=None):
        kk, vv = dict(bkd), dict(bvd)
        if tk is not None: kk[tk] = kk.get(tk, 0.0) + TH
        if tv is not None: vv[tv] = vv.get(tv, 0.0) + TH
        ring.apply_perturb(kicks=list(kk.items()), vkicks=list(vv.items()))

    ap(); co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return None
    o = ring.track_elements([co])[0]; bx0, by0 = o[bpm, 0], o[bpm, 2]
    Rxx = np.zeros((len(bpm), len(corr))); Ryy = np.zeros((len(bpm), len(corr)))
    for j, c in enumerate(corr):
        ap(tk=c); co = ra.closed_orbit(ring, dims=4)
        if np.all(np.isfinite(co)): Rxx[:, j] = (ring.track_elements([co])[0][bpm, 0] - bx0) / TH
        ap(tv=c); co = ra.closed_orbit(ring, dims=4)
        if np.all(np.isfinite(co)): Ryy[:, j] = (ring.track_elements([co])[0][bpm, 2] - by0) / TH
    return np.concatenate([Rxx.ravel(), Ryy.ravel()])


def tune(ring, resp, aperture, loco=None, oct_ctx=None, octupole_g=0.0):
    """Настроить возмущённое кольцо: орбита → LOCO → октуполь. ring мутируется (в конце все уставки
    применены). Возврат dict(ok, co, bk, bv, bd, dk3, stages) или dict(ok=False) если орбиту не вернуть.
    stages — какие шаги реально сработали."""
    stages = []
    # 1. ОРБИТА
    try:
        res = oc.correct(ring, resp, aperture=aperture)
    except Exception:
        res = None
    if res is None:
        return dict(ok=False, co=None, stages=stages)
    corr = resp["corr"]
    bk = list(zip(corr, res["thx"].tolist())); bv = list(zip(corr, res["thy"].tolist()))
    stages.append("orbit" + ("+thread" if res.get("threaded") else ""))
    # 2. LOCO
    bd = []
    if loco is not None:
        M = _loco_response(ring, loco, bk, bv)
        if M is not None:
            bd = list(zip(loco["quads"], (-(loco["Pl"] @ (M - loco["M0"]))).tolist()))
            stages.append("loco")
    # 3. ОКТУПОЛЬ
    dk3 = []
    if oct_ctx is not None and octupole_g != 0.0:
        octs, k30 = oct_ctx["octs"], oct_ctx["k30"]
        dk3 = [(octs[j], octupole_g * k30[j]) for j in range(len(octs))]
        stages.append(f"oct(g={octupole_g:+.2f})")
    # применить всё вместе
    ring.apply_perturb(dks=bd, kicks=bk, vkicks=bv, dk3s=dk3)
    co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)):
        return dict(ok=False, co=None, bk=bk, bv=bv, bd=bd, dk3=dk3, stages=stages)
    return dict(ok=True, co=co, bk=bk, bv=bv, bd=bd, dk3=dk3, stages=stages)
