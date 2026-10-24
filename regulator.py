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

    def src_of(o):
        return float(np.mean(o[SIDX, 0] + o[SIDX, 1] * DVEC))

    r.apply_perturb()                                               # базовая (невозмущённая) орбита
    o0 = r.track_elements([ra.closed_orbit(r, dims=4)])[0]
    base_mon = o0[mon, 0].copy(); base_src = src_of(o0)
    print(f"базовая орбита в мониторах: норма {np.linalg.norm(base_mon):.2e} м — вычитаем", flush=True)

    def resp_of(indices):
        H = np.zeros((len(mon), len(indices))); hs = np.zeros(len(indices))
        for j, idx in enumerate(indices):
            r.apply_perturb(kicks=[(int(idx), EPS)])
            o = r.track_elements([ra.closed_orbit(r, dims=4)])[0]
            H[:, j] = (o[mon, 0] - base_mon) / EPS                  # ВЫЧЕТ базовой орбиты (иначе o(0)/ε забивает отклик)
            hs[j] = (src_of(o) - base_src) / EPS
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

    def rpinv(m):
        return Vc.T @ ((Uc.T @ m) / Sc)                       # помодовая псевдообратная R⁺·m

    if "realistic" in sys.argv:
        # РЕАЛИЗМ: типовые значения 4-го поколения (6 ГэВ, класс СИЛА-26) — ПРЕДПОЛОЖЕНИЯ.
        # Дрейф (тепло) + вибрация (грунт/гирдеры ~50 Гц) + шум BPM; FOFB с лагом привода и ходом, SOFB offload.
        FREV = 2.99792458e8 / 1110.0; TACT_F = 135 / FREV           # FOFB такт ≈ 0.5 мс (2 кГц)
        SIG_S, RHO_S = 0.5e-6, 0.99999                              # дрейф квадруполей ~0.5 мкм, минуты
        FVIB = 50.0; SIG_F = 30e-9; RHO_F = float(np.cos(2 * np.pi * FVIB * TACT_F))  # вибрация ~30 нм, 50 Гц
        NOISE = 0.1e-6; MSOFB = 200; GF = 0.5; ALPHA = 0.1; LD = 1  # BPM шум; SOFB 10 Гц; лаг БП+вихр.токи
        TGT_Y = 0.1e-6; T = 40000                                   # цель верт. ~10% σ_y(1мкм)
        print(f"активных мод {na}/{len(S)}; вибрация 50 Гц → ρ_f={RHO_F:.4f}; цель {TGT_Y*1e6:.2f} мкм", flush=True)

        GS = 0.5                                                    # усиление SOFB (медленный интеграл)

        def run_full(sofb, fofb, stroke=np.inf, seed=1):
            rng = np.random.default_rng(seed)
            ds = np.zeros(nquad); dv = np.zeros(nquad)
            th_s = np.zeros(ncorr); thf = np.zeros(ncorr); thf_act = np.zeros(ncorr)
            qbuf = [np.zeros(ncorr) for _ in range(LD)]; acc = np.zeros(nmon); nacc = 0; sat = 0.0
            src = np.empty(T)
            for k in range(T):
                ds = RHO_S * ds + np.sqrt(1 - RHO_S ** 2) * SIG_S * rng.standard_normal(nquad)
                dv = RHO_F * dv + np.sqrt(1 - RHO_F ** 2) * SIG_F * rng.standard_normal(nquad)
                kq = K1L * (ds + dv); th = thf_act + th_s
                src[k] = hc @ th + hq @ kq
                m = Hc @ th + Hq @ kq + NOISE * rng.standard_normal(nmon)
                if fofb:                                            # быстрый интеграл + лаг привода + ход
                    thf = thf - GF * rpinv(m); ud = thf
                    if LD:
                        qbuf.append(thf.copy()); ud = qbuf.pop(0)
                    thf_act = ALPHA * thf_act + (1 - ALPHA) * ud
                    if np.isfinite(stroke):
                        cl = np.clip(thf_act, -stroke, stroke); sat += float(np.mean(cl != thf_act)); thf_act = cl
                if sofb:                                            # медленный интеграл по усреднённой орбите
                    acc += m; nacc += 1
                    if (k + 1) % MSOFB == 0:
                        th_s = th_s - GS * rpinv(acc / nacc); acc[:] = 0; nacc = 0
            return src, sat / max(T, 1)

        def rms(s, tau=30.0):
            e = s[0]; b = np.exp(-1.0 / tau); ema = np.empty_like(s)
            for k, v in enumerate(s):
                e = b * e + (1 - b) * v; ema[k] = e
            sl = slice(T // 4, None)
            return np.std(s[sl]) * 1e6, np.std(ema[sl]) * 1e6, np.std((s - ema)[sl]) * 1e6

        # ход быстрого корректора: 35% от спроса (замер по FOFB-only ход∞)
        _, _ = run_full(False, True)
        demand = []
        rng = np.random.default_rng(1); thf = np.zeros(ncorr); ds = np.zeros(nquad); dv = np.zeros(nquad)
        for k in range(T):
            ds = RHO_S * ds + np.sqrt(1 - RHO_S ** 2) * SIG_S * rng.standard_normal(nquad)
            dv = RHO_F * dv + np.sqrt(1 - RHO_F ** 2) * SIG_F * rng.standard_normal(nquad)
            thf = thf - GF * rpinv(Hc @ thf + Hq @ (K1L * (ds + dv)) + NOISE * rng.standard_normal(nmon))
            if k > T // 2:
                demand.append(np.max(np.abs(thf)))
        stroke = 0.35 * float(np.median(demand))
        print(f"ход быстрого корректора: лимит {stroke:.2e} рад (35% спроса)", flush=True)
        print(f"\n{'режим':>22} | {'источник σ':>10} | {'дрейф':>6} | {'вибр':>6} | насыщ | цель 0.1")
        for name, so, fo in [("выкл", False, False), ("только SOFB", True, False),
                             ("только FOFB", False, True), ("SOFB+FOFB", True, True)]:
            s, sat = run_full(so, fo, stroke=stroke if fo else np.inf)
            tot, dr, vb = rms(s)
            ok = "✓" if tot < TGT_Y * 1e6 else "✗"
            print(f"{name:>22} | {tot:10.3f} | {dr:6.3f} | {vb:6.3f} | {sat*100:4.1f}% | {ok}", flush=True)
        print("\n[готово realistic]", flush=True); sys.exit(0)

    if "twoloop" in sys.argv:
        # МЕХАНИКА #1: двухтемповые контуры на спектре дрейф+вибрация. Контур на темпе M тактов,
        # интегральный по средней за такт орбите. M=1 — FOFB (быстрый), M=100 — SOFB (медленный).
        RHO_S, SIG_S = 0.9999, 2e-8      # дрейф: корр. время ~10⁴ тактов
        RHO_F, SIG_F = 0.90, 1e-8        # вибрация: декоррелирует ~10 тактов
        NOISE = 0.1e-6; T = 20000

        def run_loop(M, g, seed=1):
            rng = np.random.default_rng(seed)
            th = np.zeros(ncorr); ds = np.zeros(nquad); df = np.zeros(nquad)
            acc = np.zeros(nmon); cnt = 0; src = np.empty(T)
            for k in range(T):
                ds = RHO_S * ds + np.sqrt(1 - RHO_S ** 2) * SIG_S * rng.standard_normal(nquad)
                df = RHO_F * df + np.sqrt(1 - RHO_F ** 2) * SIG_F * rng.standard_normal(nquad)
                kq = K1L * (ds + df)
                src[k] = hc @ th + hq @ kq
                m = Hc @ th + Hq @ kq + NOISE * rng.standard_normal(nmon)
                acc += m; cnt += 1
                if (k + 1) % M == 0:
                    th = th - g * rpinv(acc / cnt); acc[:] = 0; cnt = 0
            return src

        def bands(src, tau=30.0):
            e = src[0]; b = np.exp(-1.0 / tau); ema = np.empty_like(src)
            for k, v in enumerate(src):
                e = b * e + (1 - b) * v; ema[k] = e
            sl = slice(len(src) // 4, None)
            return np.std(src[sl]) * 1e6, np.std(ema[sl]) * 1e6, np.std((src - ema)[sl]) * 1e6

        print(f"активных мод SVD: {na}/{len(S)};  помеха: дрейф(ρ={RHO_S}) + вибрация(ρ={RHO_F})", flush=True)
        print(f"\n{'режим':>18} | {'итог σ':>7} | {'дрейф-полоса':>12} | {'вибр-полоса':>11}  (мкм)")
        for name, M, g in [("разомкнуто", 1, 0.0), ("только SOFB (M=100)", 100, 1.0),
                           ("только FOFB (M=1)", 1, 0.5)]:
            tot, slow, fast = bands(run_loop(M, g))
            print(f"{name:>18} | {tot:7.3f} | {slow:12.3f} | {fast:11.3f}", flush=True)
        print("\n[готово twoloop]", flush=True); sys.exit(0)

    if "offload" in sys.argv:
        # МЕХАНИКА #2: быстрые корректоры ограничены ПО ХОДУ. FOFB копит DC на дрейф → насыщается →
        # дрейф прорывается. SOFB (медл.) разгружает DC FOFB в сильные корректоры, держа его в диапазоне.
        RHO_S, SIG_S = 0.9999, 2e-8; RHO_F, SIG_F = 0.90, 1e-8
        NOISE = 0.1e-6; T = 20000; MS = 100; GF = 0.5

        def run_two(stroke=np.inf, offload=False, seed=1):
            rng = np.random.default_rng(seed)
            th_f = np.zeros(ncorr); th_s = np.zeros(ncorr); ds = np.zeros(nquad); df = np.zeros(nquad)
            src = np.empty(T); thf_sum = np.zeros(ncorr); sat = 0.0
            for k in range(T):
                ds = RHO_S * ds + np.sqrt(1 - RHO_S ** 2) * SIG_S * rng.standard_normal(nquad)
                df = RHO_F * df + np.sqrt(1 - RHO_F ** 2) * SIG_F * rng.standard_normal(nquad)
                kq = K1L * (ds + df); th = th_f + th_s
                src[k] = hc @ th + hq @ kq
                m = Hc @ th + Hq @ kq + NOISE * rng.standard_normal(nmon)
                th_f = th_f - GF * rpinv(m)                        # FOFB быстрый интеграл
                if np.isfinite(stroke):
                    cl = np.clip(th_f, -stroke, stroke); sat += float(np.mean(cl != th_f)); th_f = cl
                thf_sum += th_f
                if offload and (k + 1) % MS == 0:
                    mean_f = thf_sum / MS                         # DC-часть FOFB → в сильные корректоры (орбита-нейтрально)
                    th_s = th_s + mean_f; th_f = th_f - mean_f
                    thf_sum[:] = 0.0
            return src, sat / T

        def bands(src, tau=30.0):
            e = src[0]; b = np.exp(-1.0 / tau); ema = np.empty_like(src)
            for k, v in enumerate(src):
                e = b * e + (1 - b) * v; ema[k] = e
            sl = slice(len(src) // 4, None)
            return np.std(src[sl]) * 1e6, np.std(ema[sl]) * 1e6, np.std((src - ema)[sl]) * 1e6

        su, _ = run_two(stroke=np.inf)                            # ход неограничен — замер спроса
        rms_f = 0.0
        # оценим потребный ход по RMS уставок FOFB в неограниченном прогоне
        rng = np.random.default_rng(1); th_f = np.zeros(ncorr); ds = np.zeros(nquad); df = np.zeros(nquad); acc = []
        for k in range(T):
            ds = RHO_S * ds + np.sqrt(1 - RHO_S ** 2) * SIG_S * rng.standard_normal(nquad)
            df = RHO_F * df + np.sqrt(1 - RHO_F ** 2) * SIG_F * rng.standard_normal(nquad)
            m = Hc @ th_f + Hq @ (K1L * (ds + df)) + NOISE * rng.standard_normal(nmon)
            th_f = th_f - GF * rpinv(m)
            if k > T // 2:
                acc.append(np.max(np.abs(th_f)))
        demand = float(np.median(acc)); stroke = 0.35 * demand
        print(f"активных мод SVD: {na}/{len(S)}; спрос хода FOFB (медиана max|θ|) {demand:.2e} рад; "
              f"лимит хода = {stroke:.2e} рад (35%)", flush=True)
        print(f"\n{'режим':>26} | {'итог':>6} | {'дрейф':>6} | {'вибр':>6} | насыщ.FOFB")
        for name, stk, off in [("FOFB, ход ∞", np.inf, False), ("FOFB, ход огранич.", stroke, False),
                               ("SOFB+FOFB (offload)", stroke, True)]:
            src, sat = run_two(stroke=stk, offload=off)
            tot, slow, fast = bands(src)
            print(f"{name:>26} | {tot:6.3f} | {slow:6.3f} | {fast:6.3f} | {sat*100:5.1f}%", flush=True)
        print("\n[готово offload]", flush=True); sys.exit(0)

    if "lqg" in sys.argv:
        # МЕХАНИКА #3 (капстоун): помодовый LQG (steady-state Калман + LQR) против I/PI. Модель на моду
        # (орбитные ед.): помеха d — AR(1,ρ); привод a_{k+1}=α·a+(1−α)·u_{k−L}; y=d+a; z=d+a+ν.
        from scipy.linalg import solve_discrete_are as dare
        RHO, SIG, NOISE, T = 0.99, 2e-8, 0.1e-6, 10000
        hcm = hc @ Vc.T
        HqU = Hq.T @ Uc
        s2d_unit = ((K1L[:, None] * HqU) ** 2).sum(axis=0)        # дисперсия δ_i на единичную SIG

        def lqg_gains(alpha, L, Rctrl):
            n = 2 + L
            A = np.zeros((n, n)); A[0, 0] = RHO; A[1, 1] = alpha
            B = np.zeros((n, 1)); G = np.zeros((n, 1)); G[0, 0] = 1.0
            if L >= 1:
                A[1, n - 1] = 1 - alpha
                for j in range(1, L):
                    A[2 + j, 2 + j - 1] = 1.0
                B[2, 0] = 1.0
            else:
                B[1, 0] = 1 - alpha
            C = np.zeros((1, n)); C[0, 0] = 1.0; C[0, 1] = 1.0
            Sx = dare(A, B, C.T @ C, np.array([[Rctrl]]))
            Klqr = np.linalg.solve(np.array([[Rctrl]]) + B.T @ Sx @ B, B.T @ Sx @ A).ravel()
            Lf = np.zeros((na, n))
            for i in range(na):
                W = (1 - RHO ** 2) * (s2d_unit[i] * SIG ** 2)
                M = dare(A.T, C.T, G @ np.array([[W]]) @ G.T, np.array([[NOISE ** 2]]))
                Lf[i] = (M @ C.T / (C @ M @ C.T + NOISE ** 2)).ravel()           # фильтрационная форма
            return A, B.ravel(), Klqr, Lf

        def sim(mode, alpha=0.0, L=0, gi=1.0, kp=0.0, Rctrl=1e-6, seed=1):
            rng = np.random.default_rng(seed)
            drift = np.zeros(nquad); a = np.zeros(na); qbuf = [np.zeros(na) for _ in range(L)]
            integ = np.zeros(na); src = np.empty(T)
            if mode == "LQG":
                A, Bf, Kl, Lf = lqg_gains(alpha, L, Rctrl); xhat = np.zeros((na, 2 + L))
            for k in range(T):
                drift = RHO * drift + np.sqrt(1 - RHO ** 2) * SIG * rng.standard_normal(nquad)
                delta = Uc.T @ (Hq @ (K1L * drift))
                src[k] = hcm @ (a / Sc) + hq @ (K1L * drift)
                z = a + delta + NOISE * (Uc.T @ rng.standard_normal(nmon))
                if mode == "LQG":
                    innov = z - (xhat[:, 0] + xhat[:, 1]); xf = xhat + Lf * innov[:, None]
                    u = -(xf @ Kl); xhat = xf @ A.T + np.outer(u, Bf)
                else:
                    integ += z; u = -(gi * integ + kp * z)
                ud = u
                if L:
                    qbuf.append(u.copy()); ud = qbuf.pop(0)
                a = alpha * a + (1 - alpha) * ud
                if not np.all(np.isfinite(a)) or np.max(np.abs(a)) > 1.0:
                    return float("inf")
            return float(np.std(src[T // 4:])) * 1e6

        def bestc(mode, **kw):
            if mode == "I":
                return min(sim("I", gi=g, **kw) for g in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2])
            if mode == "PI":
                return min(sim("PI", gi=g, kp=p, **kw)
                           for g in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2] for p in [0.0, 0.2, 0.4, 0.6, 0.8])
            return min(sim("LQG", Rctrl=rc, **kw) for rc in [1e-6, 1e-5, 1e-4, 1e-3, 1e-2])

        print(f"активных мод SVD: {na}/{len(S)}; помеха AR(1) ρ={RHO}", flush=True)
        print(f"\n{'профиль привода':>30} | {'разомкн':>7} | {'I':>6} | {'PI':>6} | {'LQG':>6}  (мкм)")
        for name, alpha, L in [("идеальный", 0.0, 0), ("лаг τ=2 (α=0.607)+задержка 2", float(np.exp(-0.5)), 2)]:
            ol = sim("I", gi=0.0, alpha=alpha, L=L)
            print(f"{name:>30} | {ol:7.3f} | {bestc('I',alpha=alpha,L=L):6.3f} | "
                  f"{bestc('PI',alpha=alpha,L=L):6.3f} | {bestc('LQG',alpha=alpha,L=L):6.3f}", flush=True)
        print("\n[готово lqg]", flush=True); sys.exit(0)

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
