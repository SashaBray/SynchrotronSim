# -*- coding: utf-8 -*-
"""ПУНКТ D: переподтвердить ПАРЕТО-ФРОНТ ДА⊥↔MA октупольной силой на ИСПРАВЛЕННОЙ решётке (железо+зазоры).

На старой («сплошное железо»): проектные октуполи (g=0) ДОМИНИРУЮТСЯ — g=+0.25 (октуполи +25%) лучше и по
ДА+x, и по MA; фронт +0.25…−1.0 (октуполи выкл = макс MA). Октуполи ξ1-нейтральны, двигают малоамплитудный
MA только через ξ2. Проверяем, переносится ли вывод после правки геометрии диполей.

На каждой машине: ошибки (проектные) → орбита+LOCO → скан g (dk3=g·k3₀ транзитно) → ДА+x (звезда) и
MA=min(MA+,MA−). Медиана по ансамблю. Запуск: py -3.14 pareto_octupole.py
"""
import os, sys, pickle, tempfile, shutil
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, errors as errmod, orbit_correct as oc, ring_analysis as ra

NSEED, NSL, TURNS, NDIR, CUTOFF = 8, 8, 384, 9, 2.5
GVALS = [0.25, 0.0, -0.5, -0.75, -1.0]                 # dk3=g·k3₀ ⇒ k3/k3₀ = 1+g
inf = dict(fr.info_of(26)); inf["_id"] = 26
JD = pickle.load(open("workspace/results/loco_jac_exp26_nsl8.pkl", "rb"))
Jl, M0, quads, corr_rm, bpm_rm = JD["J"], JD["M0"], JD["quads"], JD["corr"], JD["bpm"]
Pl = np.linalg.pinv(Jl, rcond=1e-2)
ideal = fr._ring(inf, nsl=NSL); resp = oc.cached_response_matrix(ideal, 26)
co0, xe, ye, _, _ = fr.auto_ranges(ideal, inf, 512)
OCTS = [i for i, k in enumerate(ideal.kinds) if k == "oct"]     # октуполи и k3₀ — проектные, с идеала
K30 = [float(ideal.mags[i].k3) for i in OCTS]
ideal.close()
corr = resp["corr"]; LOSS = inf["loss_m"]; TH = 1e-5


def diag_rm(ring, bk, bv):
    bkd, bvd = dict(bk), dict(bv)
    def ap(tk=None, tv=None):
        kk, vv = dict(bkd), dict(bvd)
        if tk is not None: kk[tk] = kk.get(tk, 0) + TH
        if tv is not None: vv[tv] = vv.get(tv, 0) + TH
        ring.apply_perturb(kicks=list(kk.items()), vkicks=list(vv.items()))
    ap(); co = ra.closed_orbit(ring, dims=4)
    if not np.all(np.isfinite(co)): return None
    o = ring.track_elements([co])[0]; bx0, by0 = o[bpm_rm, 0], o[bpm_rm, 2]
    Rxx = np.zeros((len(bpm_rm), len(corr_rm))); Ryy = np.zeros((len(bpm_rm), len(corr_rm)))
    for j, c in enumerate(corr_rm):
        ap(tk=c); co = ra.closed_orbit(ring, dims=4); Rxx[:, j] = (ring.track_elements([co])[0][bpm_rm, 0] - bx0) / TH
        ap(tv=c); co = ra.closed_orbit(ring, dims=4); Ryy[:, j] = (ring.track_elements([co])[0][bpm_rm, 2] - by0) / TH
    return np.concatenate([Rxx.ravel(), Ryy.ravel()])


def da_x(ring, co, coarse=8, fine=6):
    """ДА+x (луч θ=0) бинпоиском-звездой."""
    if co is None or not np.all(np.isfinite(co)): return 0.0
    Rx = 1.5 * xe
    def alive(rs):
        st = [co.copy() for _ in rs]
        for s, rr in zip(st, rs): s[0] += rr * Rx
        return ra.survived_turns(ring.track_tbt(np.asarray(st), TURNS, limit=LOSS), limit=LOSS) >= TURNS
    tc = np.linspace(1.0 / coarse, 1.0, coarse); ok = alive(tc); lost = np.where(~ok)[0]
    if not lost.size: return float(tc[-1] * Rx * 1e3)
    a = tc[lost[0] - 1] if lost[0] > 0 else 0.0; b = tc[lost[0]]
    tf = np.linspace(a, b, fine); okf = alive(tf); lf = np.where(~okf)[0]
    rb = (tf[lf[0] - 1] if lf.size and lf[0] > 0 else (a if lf.size else b))
    return float(rb * Rx * 1e3)


def ma_min(ring, co, dmax=0.08, coarse=12, fine=6):
    if co is None or not np.all(np.isfinite(co)): return 0.0
    def mk(ds):
        s = np.tile(co, (len(ds), 1)).astype(float); s[:, 5] += ds; s[:, 0] += 1e-5
        return ra.survived_turns(ring.track_tbt(s, TURNS, limit=LOSS), limit=LOSS) >= TURNS
    def bs(sgn):
        tc = np.linspace(sgn * dmax / coarse, sgn * dmax, coarse); ok = mk(tc); lost = np.where(~ok)[0]
        if not lost.size: return dmax
        a = tc[lost[0] - 1] if lost[0] > 0 else 0.0; b = tc[lost[0]]
        tf = np.linspace(a, b, fine); okf = mk(tf); lf = np.where(~okf)[0]
        return abs(tf[lf[0] - 1]) if lf.size and lf[0] > 0 else abs(a if lf.size else b)
    return float(min(bs(+1), bs(-1)) * 100)


res = {g: {"da": [], "ma": []} for g in GVALS}
root = tempfile.mkdtemp(prefix="par_")
try:
    for s in range(NSEED):
        ws = os.path.join(root, str(s)); errmod.perturbed_workspace(26, 8000 + s, errmod.load_tolerances(), ws, CUTOFF)
        r = fr._ring(inf, ws=ws, nsl=NSL)
        rr = oc.correct(r, resp, aperture=LOSS)
        if rr is None:
            r.close(); print(f"seed {s}: не скорректирована (стена)", flush=True); continue
        bk = list(zip(corr, rr["thx"].tolist())); bv = list(zip(corr, rr["thy"].tolist()))
        M = diag_rm(r, bk, bv)
        bd = list(zip(quads, (-(Pl @ (M - M0))).tolist())) if M is not None else []
        for g in GVALS:
            dk3 = [(OCTS[j], g * K30[j]) for j in range(len(OCTS))]
            r.apply_perturb(dks=bd, kicks=bk, vkicks=bv, dk3s=dk3)
            co = ra.closed_orbit(r, dims=4)
            res[g]["da"].append(da_x(r, co)); res[g]["ma"].append(ma_min(r, co))
        r.close(); print(f"seed {s}: готово", flush=True)
finally:
    shutil.rmtree(root, ignore_errors=True)

med = lambda a: float(np.median(a)) if a else 0.0
print(f"\nПАРЕТО ДА⊥↔MA (исправленная решётка, медиана по {NSEED} машинам):")
print(f"{'g':>6} {'k3/k3₀':>7} | {'ДА+x, мм':>9} | {'MA, %':>7}")
for g in GVALS:
    tag = "  ← проект" if g == 0.0 else ("  (окт выкл)" if g == -1.0 else "")
    print(f"{g:+6.2f} {1+g:7.2f} | {med(res[g]['da']):9.2f} | {med(res[g]['ma']):7.2f}{tag}")
print("\n[готово]", flush=True)
