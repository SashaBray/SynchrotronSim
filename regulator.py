# -*- coding: utf-8 -*-
"""ДИНАМИКА: оптимальный регулятор для НЕПРЕДСКАЗУЕМЫХ помех (стохастика + шум BPM).

Непредсказуемую помеху (дрейф/широкополосную) feed-forward не возьмёшь — только обратная связь, и
«совершенствовать» = настроить регулятор под баланс «помеха vs шум BPM». Исследование ЛИНЕЙНОЕ (орбита
µm-уровня линейна по кикам), поэтому строим ОДИН РАЗ матрицы отклика кик→орбита-в-мониторах и кик→источник
(для корректоров и для квадруполей-feed-down), и вся петля — матрично-векторная (микросекунды/такт).
Кэш откликов: workspace/results/regulator_resp_exp26.npz.

Горизонтальная плоскость (математика регулятора по плоскостям одинакова). Помеха: независимый случайный
дрейф позиций квадруполей (random walk) + широкополосная (white). Регуляторы: наивный pinv(rcond) против
помодового Тихонова (мягкий спад слабых мод, шумо-зависимая регуляризация λ). Метрика — СКО источника.
Запуск: py -3.14 regulator.py   |   пересобрать отклик: py -3.14 regulator.py rebuild
"""
import os, sys
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, orbit_correct as oc, orbit_infer as oi, ring_analysis as ra

CACHE = os.path.join("workspace", "results", "regulator_resp_exp26.npz")
EPS = 1e-6                                              # пробный кик, рад


def build_response():
    inf = dict(fr.info_of(26)); inf["_id"] = 26
    r = fr._ring(inf, nsl=8)
    corr = oc.cached_response_matrix(r, 26)["corr"]; mon = oi.monitor_rows(r)
    quads = oi.quad_indices(r, kinds=("quad",))
    K1L = np.array([float(r.mags[q].k1) * float(r.mags[q].Leff) for q in quads])
    s_out = np.array(r.s_out); Lel = np.array([float(getattr(m, "Leff", 0.0) or 0.0) for m in r.mags])
    gaps = np.append((s_out - Lel)[1:] - s_out[:-1], 0.0)
    SIDX = np.where(gaps > 5.0)[0]; DVEC = gaps[SIDX] / 2.0

    def resp_of(indices):
        H = np.zeros((len(mon), len(indices))); hs = np.zeros(len(indices))
        for j, idx in enumerate(indices):
            r.apply_perturb(kicks=[(int(idx), EPS)])
            o = r.track_elements([ra.closed_orbit(r, dims=4)])[0]
            H[:, j] = o[mon, 0] / EPS
            hs[j] = float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC)) / EPS
            if j % 50 == 0:
                print(f"   отклик {j}/{len(indices)}", flush=True)
        return H, hs

    print("строю отклик корректоров ...", flush=True); Hc, hc = resp_of(corr)
    print("строю отклик квадруполей (feed-down) ...", flush=True); Hq, hq = resp_of(list(quads))
    r.apply_perturb(); r.close()
    np.savez(CACHE, Hc=Hc, hc=hc, Hq=Hq, hq=hq, K1L=K1L)
    print(f"сохранено: {CACHE}", flush=True)
    return dict(Hc=Hc, hc=hc, Hq=Hq, hq=hq, K1L=K1L)


if __name__ == "__main__":
    if (len(sys.argv) > 1 and sys.argv[1] == "rebuild") or not os.path.exists(CACHE):
        R = build_response()
    else:
        R = dict(np.load(CACHE))
    Hc, hc, Hq, hq, K1L = R["Hc"], R["hc"], R["Hq"], R["hq"], R["K1L"]
    nmon, ncorr = Hc.shape; nquad = Hq.shape[1]
    print(f"отклик: {nmon} мониторов, {ncorr} корректоров, {nquad} квадруполей; "
          f"усиление орбиты в BPM на кик квадруполя ~×{np.linalg.norm(Hq)/np.linalg.norm(Hc)*np.median(np.abs(K1L)):.0f}", flush=True)
    # SVD отклика корректор→монитор (пространственные моды)
    U, S, Vt = np.linalg.svd(Hc, full_matrices=False)
    print(f"сингулярные числа Hc: макс {S[0]:.2e}, мин {S[-1]:.2e}, число обусл {S[0]/S[-1]:.1e}", flush=True)

    # помодовый регулятор в пространстве мод SVD (каналы развязаны): eps=Uᵀm/S, θ=V·PID(eps)
    amask = S >= 1e-2 * S[0]
    Uc = U[:, amask]; Sc = S[amask]; Vc = Vt[amask, :]; na = int(amask.sum())

    def run_pid(Kp, Ki, Kd, T=3000, rho=0.99, sigd=2e-8, sign=0.1e-6, seed=1,
                slew=np.inf, tau=0.0, Ld=0):
        """Помодовый PID на eps=Uᵀm/S, θ=V·PID. ПРИВОД (поле не мгновенно): команда → задержка Ld тактов
        → slew-лимит (|Δθ|≤slew за такт) → лаг 1-го порядка (пост. времени tau тактов). Пучок видит
        ФАКТИЧЕСКОЕ поле th_act. slew=inf, tau=0, Ld=0 — идеальный мгновенный привод."""
        rng = np.random.default_rng(seed)
        drift = np.zeros(nquad); a = np.sqrt(1 - rho ** 2)
        integ = np.zeros(na); eprev = np.zeros(na); src = np.empty(T)
        th_slew = np.zeros(ncorr); th_act = np.zeros(ncorr)
        buf = [np.zeros(ncorr) for _ in range(Ld)]
        alpha = np.exp(-1.0 / tau) if tau > 0 else 0.0
        for k in range(T):
            drift = rho * drift + a * sigd * rng.standard_normal(nquad)
            kq = K1L * drift
            src[k] = hc @ th_act + hq @ kq                    # пучок видит ФАКТИЧЕСКОЕ поле привода
            m = Hc @ th_act + Hq @ kq + sign * rng.standard_normal(nmon)
            eps = (Uc.T @ m) / Sc
            integ += eps
            u = Vc.T @ (-(Kp * eps + Ki * integ + Kd * (eps - eprev)))   # команда PID (корректорное простр.)
            eprev = eps
            if Ld:
                buf.append(u); u = buf.pop(0)                 # задержка обработки
            th_slew = th_slew + np.clip(u - th_slew, -slew, slew)        # slew-лимит
            th_act = alpha * th_act + (1 - alpha) * th_slew   # лаг 1-го порядка (полоса БП/вихр.токи)
            if not np.all(np.isfinite(th_act)) or np.max(np.abs(th_act)) > 1e-2:
                return float("inf")
        return float(np.std(src[T // 4:]))

    KI = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]; KP = [0.0, 0.2, 0.4, 0.6, 0.8]; KD = [0.0, 0.1, 0.2, 0.4]

    def best(gen, **kw):
        b = (float("inf"), None)
        for g in gen:
            v = run_pid(*g, **kw)
            if v < b[0]:
                b = (v, g)
        return b

    if "pid" in sys.argv or "act" in sys.argv:
        print(f"активных мод SVD: {na}/{len(S)}", flush=True)
        if "act" in sys.argv:                                 # профили привода (дрейф ρ=0.99)
            cases = [("идеальный привод (мгновенный)", dict()),
                     ("лаг τ=1 такт + задержка 1", dict(tau=1.0, Ld=1)),
                     ("лаг τ=2 + задержка 2 (медленный БП/вихр.токи)", dict(tau=2.0, Ld=2)),
                     ("slew-лимит 2e-4/такт + лаг τ=1", dict(slew=2e-4, tau=1.0, Ld=1))]
        else:                                                 # идеальный привод, разные помехи/шум
            cases = [("дрейф ρ=0.99, шум BPM 0.1 мкм", dict()),
                     ("быстрее ρ=0.90", dict(rho=0.90)),
                     ("дрейф, ШУМ ×5 = 0.5 мкм", dict(sign=0.5e-6))]
        for name, kw in cases:
            ol = run_pid(0, 0, 0, **kw)
            bI = best(((0, ki, 0) for ki in KI), **kw)
            bPI = best(((kp, ki, 0) for kp in KP for ki in KI), **kw)
            bPID = best(((kp, ki, kd) for kp in KP for ki in KI for kd in KD), **kw)
            print(f"\n=== {name} ===  (разомкнуто {ol*1e6:.3f} мкм)", flush=True)
            print(f"  I   : {bI[0]*1e6:.3f} мкм   Ki={bI[1][1]}", flush=True)
            print(f"  PI  : {bPI[0]*1e6:.3f} мкм   Kp={bPI[1][0]}, Ki={bPI[1][1]}", flush=True)
            print(f"  PID : {bPID[0]*1e6:.3f} мкм   Kp={bPID[1][0]}, Ki={bPID[1][1]}, Kd={bPID[1][2]}", flush=True)
        print("\n[готово]", flush=True); sys.exit(0)

    Pinv2 = np.linalg.pinv(Hc, rcond=1e-2)               # наивный, консервативное отсечение мод
    Pinv3 = np.linalg.pinv(Hc, rcond=1e-3)               # наивный, жадное (слабые моды внутри → шум/расходимость)

    def run(reg, T=4000, rho=0.99, sigd=2e-8, sign=0.1e-6, seed=1):
        """Петля в модельном такте с ЗАДЕРЖКОЙ 1 такт (измеряем такт k по θ прошлого такта, коррекция
        вступает на k+1 — физично: белый шум задавить нельзя). Помеха — AR(1) на позициях квадруполей
        (ρ: 1=дрейф, 0=белый). sigd — стац. СКО позиции квадруполя (м); sign — шум BPM (м). inf при расходимости."""
        rng = np.random.default_rng(seed)
        theta = np.zeros(ncorr); drift = np.zeros(nquad); a = np.sqrt(1 - rho ** 2)
        src = np.empty(T)
        for k in range(T):
            drift = rho * drift + a * sigd * rng.standard_normal(nquad)
            kq = K1L * drift                             # feed-down кики квадруполей (помеха этого такта)
            src[k] = hc @ theta + hq @ kq                # источник ПО θ прошлого такта (задержка)
            m = Hc @ theta + Hq @ kq + sign * rng.standard_normal(nmon)   # измерение (с шумом BPM)
            theta = theta + reg(m)                       # коррекция — со следующего такта
            if not np.all(np.isfinite(theta)) or np.max(np.abs(theta)) > 1e-2:
                return float("inf")                      # разошлось
        return float(np.std(src[T // 4:]))

    reg_p2 = lambda g: (lambda m: -g * (Pinv2 @ m))
    reg_p3 = lambda g: (lambda m: -g * (Pinv3 @ m))
    reg_tk = lambda g, lam: (lambda m: -g * (Vt.T @ ((S / (S ** 2 + lam)) * (U.T @ m))))
    lam_grid = (S[0] ** 2) * np.array([1e-4, 1e-3, 1e-2, 1e-1, 3e-1])

    def best_tikh(g, **kw):
        return min(run(reg_tk(g, lam), **kw) for lam in lam_grid)

    fmt = lambda v: "расх." if not np.isfinite(v) else f"{v*1e6:.3f}"
    import csv
    print("\n=== A. медленный дрейф ρ=0.99, шум BPM 0.1 мкм — развёртка усиления ===", flush=True)
    ol = run(lambda m: np.zeros(ncorr)); print(f"без петли: СКО источника {ol*1e6:.3f} мкм", flush=True)
    print(f"{'усил':>5} | pinv(1e-2) | pinv(1e-3) | Тихонов(λ*), мкм")
    gains = [0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6, 1.9]
    rowsA = []
    for g in gains:
        p2 = run(reg_p2(g)); p3 = run(reg_p3(g)); tk = best_tikh(g)
        rowsA.append((g, p2, p3, tk))
        print(f"{g:5.2f} | {fmt(p2):>10} | {fmt(p3):>10} | {fmt(tk)}", flush=True)
    with open(os.path.join("workspace", "results", "regulator_sweep_exp26.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter=";"); w.writerow(["gain", "pinv1e2_um", "pinv1e3_um", "tikhonov_um"])
        for g, p2, p3, tk in rowsA:
            w.writerow([f"{g:.2f}", fmt(p2), fmt(p3), fmt(tk)])

    print("\n=== B. развёртка ρ (дрейф→широкополосная), каждый регулятор в своём лучшем усилении ===", flush=True)
    print(f"{'ρ':>5} | {'без петли':>9} | pinv(1e-2) | Тихонов | трактовка")
    for rho in [0.99, 0.95, 0.9, 0.7, 0.4, 0.0]:
        ol = run(lambda m: np.zeros(ncorr), rho=rho)
        bp2 = min(run(reg_p2(g), rho=rho) for g in gains)
        btk = min(best_tikh(g, rho=rho) for g in gains)
        tag = "дрейф: берётся" if rho >= 0.95 else ("широкополосн.: петля слабеет" if rho >= 0.4 else "белый: петля не помогает")
        print(f"{rho:5.2f} | {ol*1e6:9.3f} | {fmt(bp2):>10} | {fmt(btk):>7} | {tag}", flush=True)
    print("\n[готово]", flush=True)
