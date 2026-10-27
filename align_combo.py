# -*- coding: utf-8 -*-
"""ПОВТОР combo_explorer, НО С ЮСТИРОВКОЙ. Вопрос пользователя: насколько БОЛЬШЕ допуски (dx/dy/dz,
крен/тангаж/рыскание, ошибка поля), чем в чистой протяжке+корректорах, если РАЗРЕШИТЬ ещё и ДВИГАТЬ
магниты (итеративная слепая юстировка) — алгоритмы работают СОВМЕСТНО.

Две руки на КАЖДОЙ машине:
  BASELINE = combo_explorer  : протяжка (_thread) + SVD-коррекция орбиты корректорами. Магниты ЗАМОРОЖЕНЫ.
  COMBO    = +ЮСТИРОВКА      : до коррекции — слепо оценить поперечные сдвиги квадруполей/DQ по СЫРОЙ орбите
                               (корректоры ВЫКЛ → сигнал сдвига не вырожден с корректорами), ПОДВИНУТЬ магниты
                               (feed-down кик K1L·δ в месте магнита = сдвиг на δ, к 1-му порядку), повторить;
                               затем протяжка+корректоры добирают остаток. Юстировка снимает ИСТОЧНИК у сильных
                               градиентов, корректоры — симптом.

Метрика envelope: доля ПОДНЯТЫХ (замкнутая орбита И ДА+x>0.3мм) и медианная ДА+x, по уровню допуска —
BASELINE vs COMBO. Где у COMBO граница (где доля падает) дальше, чем у BASELINE, — там юстировка расширила
допуск. Юстировка слепа и видит только ПОПЕРЕЧНЫЙ сдвиг (крен/рыскание/поле по орбите не наблюдаемы её
базисом) → ждём рост поперечного envelope, прочие ~без изменений (честно показываем DOF-специфичность).

Запуск: py -3.14 align_combo.py           (полный свип, фон)
        py -3.14 align_combo.py smoke      (1 машина: проверка знака кика + тайминг)
"""
import os, sys, tempfile, shutil
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, errors as errmod, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

LAT, NSL, CUTOFF, TURNS, NSEED = 26, 8, 2.5, 384, 6
APER = 0.012
NALIGN = 5                                   # итераций юстировки (двигаем магниты)
RELAX = 0.5                                  # недорелаксация шага юстировки (устойчивость на недоопредел. базисе)
MOVE_CLIP = 300e-6                           # предохранитель: макс. сдвиг магнита за итерацию [м]
BPM_NOISE = 0.1e-6                           # шум датчика (теор. граница)
SGN_X, SGN_Y = -1.0, +1.0                    # знак компенсации feed-down (валидируется smoke; остаточная орбита должна падать)
DA_MIN = 0.3e-3                              # пригодная машина: ДА+x > 0.3 мм

# свип: поперечный (главный, несёт стену бутстрапа) гуще и дальше старой стены ~25-80 мкм; поле/рыскание/крен
# по нескольку уровней — показать, что слепая юстировка их НЕ трогает (DOF-специфичность)
SMOKE = len(sys.argv) > 1 and sys.argv[1] == "smoke"
if "other" in sys.argv:                                      # #2: не-поперечные DOF (ждём combo≈baseline)
    LEVELS = {
        "field": [0.15, 0.5, 1.0],                           # % (ошибка поля)
        "yaw":   [3e-4, 7.5e-4, 1.5e-3],                     # рад (рыскание)
        "roll":  [5e-4, 1.5e-3, 3e-3],                       # рад (крен)
    }
else:
    LEVELS = {
        "trans": [20e-6, 40e-6, 60e-6, 80e-6],
    }
COL = {"trans": ("dx", "dy"), "field": ("field_pct",), "yaw": ("dyaw",), "roll": ("droll",)}
KINDS = ("dipole", "dq", "quad", "sext", "oct")   # допуск вешаем на все магниты (как в tol_sweep3)

inf = dict(fr.info_of(LAT)); inf["_id"] = LAT
ideal = fr._ring(inf, nsl=NSL)
resp = oc.cached_response_matrix(ideal, LAT)
corr = resp["corr"]; bpm_all = np.asarray(resp["bpm"])
Px = np.linalg.pinv(resp["Rx"], rcond=1e-2); Py = np.linalg.pinv(resp["Ry"], rcond=1e-2)
IDEAL_DIAG = "ideal" in sys.argv                                  # ideal → орбита известна на КАЖДОМ элементе
HOMO = "homo" in sys.argv                                         # homo → гомотопия-рамп на стене (60-100 мкм)
bpm_mon = oi.monitor_rows(ideal)                                  # реальные датчики (kind=monitor)
abpm = bpm_all if IDEAL_DIAG else bpm_mon                         # где юстировка «видит» орбиту
aelems_idx = oi.quad_indices(ideal, kinds=("quad", "dq"))         # сильные градиенты: их двигаем
abasis = oi.shift_basis(ideal, elems=aelems_idx, bpm=abpm, lattice=LAT, nsl=NSL)
aelems = abasis["elems"]; k1l = abasis["k1l"]
print(f"юстируем {len(aelems)} магнитов (quad+dq); диагностика {'ИДЕАЛ (все элементы)' if IDEAL_DIAG else f'{len(bpm_mon)} датчиков'}; "
      f"корректоров {len(corr)}", flush=True)
SCR = os.path.join("workspace", "results", "_ws_combo")


def tol_for(etype, v):
    base = dict(dx=0.0, dy=0.0, dz=0.0, dyaw=0.0, dpitch=0.0, droll=0.0, field_pct=0.0)
    for c in COL[etype]:
        base[c] = v
    return {k: dict(base) for k in KINDS}


def apply_all(r, KK, VV, DKS):
    r.apply_perturb(kicks=list(KK.items()), vkicks=list(VV.items()), dks=DKS)


def correct_orbit(r, KK, VV, DKS, seed, iters=4):
    """SVD-коррекция орбиты корректорами, СОХРАНЯЯ кики юстировки в KK/VV (в отличие от oc.correct, который
    сбрасывает). Накапливает корректорные углы в KK/VV. Возврат (co, seed) или (None, seed)."""
    co = None
    for _ in range(iters):
        apply_all(r, KK, VV, DKS)
        co = oc._closed_orbit_seeded(r, seed)
        if co is None:
            return None, seed
        seed = co[:4]
        orb = r.track_elements([co])[0]
        dthx = -(Px @ orb[bpm_all, 0]); dthy = -(Py @ orb[bpm_all, 2])
        for j, c in enumerate(corr):
            KK[c] = KK.get(c, 0.0) + float(dthx[j]); VV[c] = VV.get(c, 0.0) + float(dthy[j])
    apply_all(r, KK, VV, DKS)
    co = oc._closed_orbit_seeded(r, seed)
    return co, (co[:4] if co is not None else seed)


def _res_mon(r, co):
    o = r.track_elements([co])[0]
    return float(1e6 * np.sqrt(np.mean(o[bpm_mon, 0] ** 2 + o[bpm_mon, 2] ** 2)))


def bring_up(r, align, rng, trace=False):
    """Поднять машину. align=False → combo_explorer (протяжка+SVD-коррекция, магниты заморожены).
    align=True → +ЮСТИРОВКА: после коррекции слепо оценить сдвиги по ОСТАТОЧНОЙ орбите (корректоры сняли свой
    достижимый подпространство → остаток лежит в недостижимом 120-корректорами и несёт сигнал сдвига, не
    вырожден), ПОДВИНУТЬ магниты с недорелаксацией, перекорректировать; защита — откат, если потеряли орбиту.
    Возврат (co, KK, VV) или (None,·,·)."""
    KK, VV, DKS = {}, {}, []
    r.apply_perturb()
    co = ra.closed_orbit(r, dims=4)
    if not np.all(np.isfinite(co)):                               # бутстрап протяжкой
        thx, thy, co = oc._thread(r, resp, APER)
        if co is None:
            return None, KK, VV
        for c, t in zip(corr, thx): KK[c] = float(t)
        for c, t in zip(corr, thy): VV[c] = float(t)
    seed = co[:4]
    co, seed = correct_orbit(r, KK, VV, DKS, seed)                # штатная коррекция (= BASELINE)
    if co is None:
        return None, KK, VV
    if align:
        if trace:
            print(f"  [после коррекции] остаток орбиты {_res_mon(r, co):.2f} мкм", flush=True)
        for it in range(NALIGN):
            meas = oi.measure_orbit(r, abpm, BPM_NOISE, None, rng)      # ОСТАТОЧНАЯ орбита (корректоры включены)
            if meas is None:
                break
            ex, ey = oi.infer_shifts(abasis, meas[0], meas[1])
            dx = np.clip(SGN_X * RELAX * k1l * ex, -MOVE_CLIP, MOVE_CLIP)   # relaxed+клип (устойчивость)
            dy = np.clip(SGN_Y * RELAX * k1l * ey, -MOVE_CLIP, MOVE_CLIP)
            KKs, VVs = dict(KK), dict(VV)                          # снимок для отката
            for i, j in enumerate(aelems):                        # ДВИГАЕМ магниты (feed-down кик)
                KK[j] = KK.get(j, 0.0) + float(dx[i]); VV[j] = VV.get(j, 0.0) + float(dy[i])
            co2, seed2 = correct_orbit(r, KK, VV, DKS, seed)      # корректоры добирают остаток после сдвига
            if co2 is None:                                       # сдвиг разрушил орбиту → откат и стоп
                KK, VV = KKs, VVs; apply_all(r, KK, VV, DKS)
                co = oc._closed_orbit_seeded(r, seed); break
            co, seed = co2, seed2
            if trace:
                print(f"  [юст {it}] остаток орбиты {_res_mon(r, co):.2f} мкм", flush=True)
    return co, KK, VV


def da_x(r, co, nturns=TURNS, amax=APER):
    base = np.array([co[0], co[1], co[2], co[3], 0.0, 0.0])
    def alive(x):
        st = base.copy(); st[0] += x
        tb = r.track_tbt([st.tolist()], nturns, amax)[0]
        return np.all(np.isfinite(tb[-1]))
    if not alive(0.2e-3):
        return 0.0
    lo, hi = 0.0, amax
    for _ in range(14):
        mid = 0.5 * (lo + hi)
        if alive(mid): lo = mid
        else: hi = mid
    return lo


def orbit_rms(r, co):
    o = r.track_elements([co])[0]
    return float(1e6 * np.sqrt(np.mean(o[bpm_mon, 0] ** 2 + o[bpm_mon, 2] ** 2)))


def run_machine(etype, v, seed, align):
    ws = os.path.join(SCR, f"{etype}_{int(round(v*1e9))}_{seed}")
    errmod.perturbed_workspace(LAT, 8000 + seed, tol_for(etype, v), ws, CUTOFF)
    r = fr._ring(inf, ws=ws, nsl=NSL)
    rng = np.random.default_rng(100 + seed)
    try:
        co, KK, VV = bring_up(r, align, rng)
        if co is None:
            return None
        return dict(da=da_x(r, co), orb=orbit_rms(r, co),
                    load=float(1e6 * np.sqrt(np.mean([x**2 for x in list(KK.values()) + list(VV.values())] or [0.0]))))
    finally:
        r.close(); shutil.rmtree(ws, ignore_errors=True)


FR = [0.3, 0.5, 0.7, 0.85, 1.0]                              # доли ошибки для гомотопии (мягкий рамп)
_CORR_SET = set(corr)


def thread_base(r, KK, VV, passes=15, rcond=5e-2, relax=0.85):
    """Протяжка корректорами, СОХРАНЯЯ кики юстировки (индексы quad/dq) как фикс. базу (корректоры — индексы
    kicker, непересекающиеся). Копия oc._thread с базой. Возврат co (и корректорные углы записаны в KK/VV) или None."""
    base_k = [(j, x) for j, x in KK.items() if j not in _CORR_SET]
    base_v = [(j, x) for j, x in VV.items() if j not in _CORR_SET]
    Tx, Ty, bpm = resp["Tx"], resp["Ty"], resp["bpm"]
    thx = np.zeros(len(corr)); thy = np.zeros(len(corr))
    for _ in range(passes):
        r.apply_perturb(kicks=base_k + list(zip(corr, thx.tolist())), vkicks=base_v + list(zip(corr, thy.tolist())))
        tr = r.track_elements([np.zeros(6)])[0]
        x, y = tr[bpm, 0], tr[bpm, 2]
        within = np.isfinite(x) & np.isfinite(y) & (np.abs(x) <= APER) & (np.abs(y) <= APER)
        reach = int(np.argmin(within)) if not within.all() else len(bpm)
        if reach == 0:
            return None
        rr = np.arange(reach)
        co = oc._closed_orbit_seeded(r, np.array([np.mean(x[rr]), 0.0, np.mean(y[rr]), 0.0]))
        if co is not None:
            for i, c in enumerate(corr):
                KK[c] = float(thx[i]); VV[c] = float(thy[i])
            return co
        thx = thx - relax * (np.linalg.pinv(Tx[rr], rcond=rcond) @ x[rr])
        thy = thy - relax * (np.linalg.pinv(Ty[rr], rcond=rcond) @ y[rr])
    return None


def bring_up_homotopy(etype, v, seed, align, rng):
    """ГОМОТОПИЯ-с-юстировкой: включаем ошибку долями f (та же машина ×f), на каждой доле — тёплый старт
    корректоров+юстировки, коррекция (протяжка при нужде, с сохранением киков юстировки), [юстировка];
    к f=1.0 магниты уже частично поставлены и орбита жива там, где холодная протяжка падала. Возврат
    (r, co, KK, VV) на полной машине (co=None если не поднялась). r НЕ закрыт — закрыть у вызывающего."""
    KK, VV, DKS = {}, {}, []
    r = None; seed_co = np.zeros(4)
    for fi, f in enumerate(FR):
        ws = os.path.join(SCR, f"homo_{etype}_{int(round(v*1e9))}_{seed}_{fi}")
        errmod.perturbed_workspace(LAT, 8000 + seed, tol_for(etype, v * f), ws, CUTOFF)
        if r is not None:
            r.close()
        r = fr._ring(inf, ws=ws, nsl=NSL)
        apply_all(r, KK, VV, DKS)
        co = oc._closed_orbit_seeded(r, seed_co)
        if co is None:                                       # тёплый старт не дал орбиты → протяжка (с базой юстировки)
            co = thread_base(r, KK, VV)
            if co is None:
                shutil.rmtree(ws, ignore_errors=True); return r, None, KK, VV
        co, seed_co = correct_orbit(r, KK, VV, DKS, co[:4])
        if co is None:
            shutil.rmtree(ws, ignore_errors=True); return r, None, KK, VV
        if align:
            for _ in range(NALIGN):
                meas = oi.measure_orbit(r, abpm, BPM_NOISE, None, rng)
                if meas is None:
                    break
                ex, ey = oi.infer_shifts(abasis, meas[0], meas[1])
                dx = np.clip(SGN_X * RELAX * k1l * ex, -MOVE_CLIP, MOVE_CLIP)
                dy = np.clip(SGN_Y * RELAX * k1l * ey, -MOVE_CLIP, MOVE_CLIP)
                KKs, VVs = dict(KK), dict(VV)
                for i, j in enumerate(aelems):
                    KK[j] = KK.get(j, 0.0) + float(dx[i]); VV[j] = VV.get(j, 0.0) + float(dy[i])
                co2, seed2 = correct_orbit(r, KK, VV, DKS, seed_co)
                if co2 is None:
                    KK, VV = KKs, VVs; apply_all(r, KK, VV, DKS)
                    co = oc._closed_orbit_seeded(r, seed_co); break
                co, seed_co = co2, seed2
        shutil.rmtree(ws, ignore_errors=True)
    return r, co, KK, VV


def run_machine_homo(etype, v, seed, align):
    rng = np.random.default_rng(100 + seed)
    r, co, KK, VV = bring_up_homotopy(etype, v, seed, align, rng)
    try:
        if r is None or co is None:
            return None
        return dict(da=da_x(r, co), orb=orbit_rms(r, co))
    finally:
        if r is not None:
            r.close()


if SMOKE:
    print("\n[smoke] проверка знака юстировки (остаток орбиты должен ПАДАТЬ) + одна машина baseline vs combo", flush=True)
    ws = os.path.join(SCR, "smoke")
    errmod.perturbed_workspace(LAT, 8000, tol_for("trans", 20e-6), ws, CUTOFF)
    r = fr._ring(inf, ws=ws, nsl=NSL); rng = np.random.default_rng(100)
    co, KK, VV = bring_up(r, True, rng, trace=True)
    db = da_x(r, co) if co is not None else 0.0; ob = orbit_rms(r, co) if co is not None else float("nan")
    r.close()
    r = fr._ring(inf, ws=ws, nsl=NSL)
    co0, _, _ = bring_up(r, False, np.random.default_rng(100))
    d0 = da_x(r, co0) if co0 is not None else 0.0; o0 = orbit_rms(r, co0) if co0 is not None else float("nan")
    r.close(); shutil.rmtree(ws, ignore_errors=True); ideal.close()
    print(f"\n[smoke trans 20мкм] BASELINE ДА {d0*1e3:.2f} мм, орбита {o0:.1f} мкм | "
          f"COMBO ДА {db*1e3:.2f} мм, орбита {ob:.1f} мкм", flush=True)
    sys.exit(0)

os.makedirs(SCR, exist_ok=True)

if HOMO and "quick" in sys.argv:                                 # быстрая проверка пути: 1 машина, 60 мкм
    print("\n[homo quick] trans 60 мкм, seed0", flush=True)
    import time as _t
    for al in (False, True):
        t0 = _t.time(); r = run_machine_homo("trans", 60e-6, 0, al)
        tag = "ГОМО+ЮСТ" if al else "ГОМО-только"
        if r is None:
            print(f"  {tag}: НЕ поднято ({_t.time()-t0:.0f} c)", flush=True)
        else:
            print(f"  {tag}: ДА {r['da']*1e3:.2f} мм, орбита {r['orb']:.0f} мкм ({_t.time()-t0:.0f} c)", flush=True)
    ideal.close(); sys.exit(0)

if HOMO:                                                         # СТЕНА: гомотопия-только vs гомотопия+юстировка
    HLEV = [60e-6, 70e-6, 80e-6]; NH = 5
    print(f"\n[ГОМОТОПИЯ-РАМП {FR}] диагностика {'ИДЕАЛ' if IDEAL_DIAG else f'{len(bpm_mon)} датчиков'}", flush=True)
    print(f"\n{'σ, мкм':>6} | {'ГОМО-только подн·ДА':>20} | {'ГОМО+ЮСТ подн·ДА':>20} | {'орбита Г→ГЮ мкм':>16}", flush=True)
    for v in HLEV:
        ho = [run_machine_homo("trans", v, s, False) for s in range(NH)]
        hj = [run_machine_homo("trans", v, s, True) for s in range(NH)]
        def aggh(rr):
            up = [x for x in rr if x is not None and x["da"] > DA_MIN]
            return len(up), (np.median([x["da"] for x in up]) * 1e3 if up else 0.0), \
                   (np.median([x["orb"] for x in up]) if up else float("nan"))
        uo, do, oo = aggh(ho); uj, dj, oj = aggh(hj)
        print(f"{v*1e6:6.0f} | {f'{uo}/{NH}':>8} {do:11.2f} | {f'{uj}/{NH}':>8} {dj:11.2f} | {oo:6.0f}→{oj:<6.0f}", flush=True)
    ideal.close()
    print("\n[готово] где ГОМО+ЮСТ поднято>ГОМО-только — там юстировка пробивает стену дальше рампа", flush=True)
    sys.exit(0)

print(f"\n{'тип':>6} {'уровень':>9} | {'BASELINE поднято ДА+x мед':>26} | {'COMBO поднято ДА+x мед':>24} | {'орбита B→C мкм':>14}", flush=True)
for etype, levs in LEVELS.items():
    for v in levs:
        bl = [run_machine(etype, v, s, False) for s in range(NSEED)]
        cm = [run_machine(etype, v, s, True) for s in range(NSEED)]
        def agg(rr):
            up = [x for x in rr if x is not None and x["da"] > DA_MIN]
            da = np.median([x["da"] for x in up]) * 1e3 if up else 0.0
            orb = np.median([x["orb"] for x in up]) if up else float("nan")
            return len(up), da, orb
        ub, db, obm = agg(bl); uc, dc, ocm = agg(cm)
        unit = "мкм" if etype in ("trans",) else ("%" if etype == "field" else "мрад")
        val = v * 1e6 if etype == "trans" else (v if etype == "field" else v * 1e3)
        print(f"{etype:>6} {val:7.1f}{unit:>2} | {f'{ub}/{NSEED}':>8} {db:17.2f} | {f'{uc}/{NSEED}':>8} {dc:15.2f} | "
              f"{obm:6.0f}→{ocm:<6.0f}", flush=True)
ideal.close()
print("\n[готово] где COMBO поднято>BASELINE поднято или ДА выше — там юстировка расширила допуск", flush=True)
