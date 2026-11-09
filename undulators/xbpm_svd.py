# -*- coding: utf-8 -*-
"""ВАРИАНТ А — XBPM как дополнительные строки в ГЛОБАЛЬНОМ SVD коррекции орбиты.

Идея: один регулятор стабилизирует СОВМЕСТНЫЙ вектор измерений [орбита e⁻ в мониторах ; положение ФОТОННОГО
пучка в XBPM]. К электронной матрице отклика R_e (BPM×корректоры) снизу подшиваются фотонные строки R_ph
(XBPM×корректоры). Псевдообратная augment-матрицы даёт закон θ = −P·[o_e ; w·o_ph].

Фотонная строка берётся ПРЯМО из модели экранов: центроид светового пятна на XBPM, стоящем на расстоянии d
ниже источника, в линейном приближении = положение + d·угол ЭЛЕКТРОННОЙ орбиты в центре ондулятора:
    o_ph(b) = x_ист(b) + d_b·x'_ист(b).
Значит R_ph(b,·) = (∂x_ист/∂θ) + d_b·(∂x'_ист/∂θ) — линейная комбинация откликов ПОЛОЖЕНИЯ и УГЛА орбиты
в точке источника. Положение берём из кэша R_e (в нём есть все строки), угол ∂x'/∂θ досчитываем прогоном.

w — ВЕС фотонных строк (главная ручка варианта А): балансирует СКО орбиты e⁻ против СКО фотона в общем МНК.

Печатает: размерности augment-матрицы, спектр сингулярных чисел и число обусловленности vs вес w, норму
оператора P (усиление «хода» корректоров), и ВЫЧИСЛИТЕЛЬНУЮ НАГРУЗКУ (время MVM/SVD, FLOP/с, память) против
базовой (только e-BPM) на тактах 2/10/100 кГц. График: workspace/results/screens/xbpm_svd_L{n}.png

ЗАПУСК: py -3.14 undulators/xbpm_svd.py --lattice 26 --dist 15 --rcond 1e-2
"""
import os, sys, math, time, argparse, pickle
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")
import numpy as np
import ring_cpp, facility_report as fr, orbit_correct as oc, orbit_infer as oi
import ring_analysis as ra
import errors as er
from screen_analysis import RESDIR

CLIGHT = 2.99792458e8


def source_rows(ring, lattice):
    """Строка track_elements в центре каждого ондулятора (по ближайшему mag.center к pivot группы).
    Возврат: (gids, obs_rows) — список групп и индексы строк-наблюдателей (row = индекс элемента +1)."""
    cen = np.array([m.center for m in ring.mags])
    gids, rows = [], []
    for g in sorted(er.load_groups(lattice), key=lambda x: x["group_id"]):
        piv = np.asarray(g["pivot"], float)
        i = int(np.argmin(np.linalg.norm(cen - piv, axis=1)))
        gids.append(g["group_id"]); rows.append(i + 1)
    return gids, np.array(rows)


def lever_arms(gids, lattice, default_d):
    """Плечо источник→XBPM d_b [м] из screens.csv (dist_m по source_ref), иначе проектное default_d."""
    import csv
    d = {gid: default_d for gid in gids}
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "workspace", "tables", "screens.csv")
    if os.path.exists(path):
        for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
            if r.get("source_ref") in d and r.get("dist_m"):
                d[r["source_ref"]] = float(r["dist_m"])
    return np.array([d[g] for g in gids])


def angle_response(ring, corr, obs, theta=1e-5, lattice=None, nsl=None, rebuild=False):
    """Отклик УГЛА орбиты в точках-источниках на единичный кик корректора: Ax,Ay [рад/рад] = ∂x'/∂θ.
    Строится прогоном (closed_orbit + track_elements для каждого корректора). Кэш на диске."""
    nsl = int(getattr(ring, "nsl", 8)) if nsl is None else nsl
    path = (os.path.join(oc.RESP_DIR, f"xbpm_angresp_exp{int(lattice)}_nsl{nsl}.pkl")
            if lattice is not None else None)
    if path and not rebuild and os.path.exists(path):
        try:
            r = pickle.load(open(path, "rb"))
            if r.get("corr") == list(corr) and np.array_equal(r.get("obs"), obs):
                return r["Ax"], r["Ay"]
        except Exception:
            pass
    ring.apply_perturb()
    co0 = ra.closed_orbit(ring, dims=4)
    base = ring.track_elements([co0])[0]
    pxp0, pyp0 = base[obs, 1], base[obs, 3]
    Ax = np.zeros((len(obs), len(corr))); Ay = np.zeros((len(obs), len(corr)))
    for j, c in enumerate(corr):
        ring.apply_perturb(kicks=[(c, theta)])
        co = ra.closed_orbit(ring, dims=4); orb = ring.track_elements([co])[0]
        Ax[:, j] = (orb[obs, 1] - pxp0) / theta
        ring.apply_perturb(vkicks=[(c, theta)])
        co = ra.closed_orbit(ring, dims=4); orb = ring.track_elements([co])[0]
        Ay[:, j] = (orb[obs, 3] - pyp0) / theta
    ring.apply_perturb()
    if path:
        os.makedirs(oc.RESP_DIR, exist_ok=True)
        pickle.dump(dict(corr=list(corr), obs=obs, Ax=Ax, Ay=Ay), open(path, "wb"))
    return Ax, Ay


def cond_vs_weight(Re, Rph, weights, rcond):
    """Число обусл. и норма оператора P=pinv([Re; w·Rph]) как функция веса w фотонных строк."""
    out = []
    for w in weights:
        Raug = np.vstack([Re, w * Rph])
        sv = np.linalg.svd(Raug, compute_uv=False)
        svk = sv[sv > rcond * sv[0]]
        P = np.linalg.pinv(Raug, rcond=rcond)
        out.append((w, sv[0], svk[-1], sv[0] / svk[-1], float(np.linalg.norm(P, 2))))
    return out


def bench_mvm(P, nrow, reps=2000):
    """Среднее время матрично-векторного умножения θ=P·o (один такт FOFB, одна плоскость), с."""
    o = np.random.default_rng(0).normal(size=nrow)
    P @ o                                                        # прогрев
    t0 = time.perf_counter()
    for _ in range(reps):
        _ = P @ o
    return (time.perf_counter() - t0) / reps


def main():
    ap = argparse.ArgumentParser(description="вариант А: XBPM в глобальном SVD + вычислительная нагрузка")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--dist", type=float, default=15.0, help="проектное плечо источник→XBPM, м (если нет в screens.csv)")
    ap.add_argument("--rcond", type=float, default=1e-2)
    ap.add_argument("--nsl", type=int, default=8)
    ap.add_argument("--rebuild", action="store_true")
    a = ap.parse_args()

    inf = dict(fr.info_of(a.lattice)); inf["_id"] = a.lattice
    ring = fr._ring(inf, nsl=a.nsl)
    corr = oc.corrector_indices(ring)
    mon = oi.monitor_rows(ring)
    gids, obs = source_rows(ring, a.lattice)
    dists = lever_arms(gids, a.lattice, a.dist)
    C = float(ring.circumference); frev = CLIGHT / C

    print(f"=== ВАРИАНТ А: XBPM в глобальном SVD — решётка {a.lattice} ===")
    print(f"периметр {C:.1f} м, f_об {frev/1e3:.1f} кГц; "
          f"корректоров {len(corr)}, e-BPM мониторов {len(mon)}, ондуляторных пучков {len(gids)}")
    print(f"плечо источник→XBPM: {dists.min():.1f}–{dists.max():.1f} м (из screens.csv, иначе {a.dist:.0f})")

    # электронная матрица отклика (кэш): положение во всех строках → берём мониторы и источники
    resp = oc.cached_response_matrix(ring, a.lattice, rebuild=a.rebuild)
    Rx, Ry = resp["Rx"], resp["Ry"]
    Re_x, Re_y = Rx[mon - 1], Ry[mon - 1]                        # [n_mon × n_corr], положение e⁻ в мониторах
    Pe_x = Rx[obs - 1], Ry[obs - 1]                             # положение в источниках
    # угол в источниках (прогон + кэш)
    print("строю отклик УГЛА орбиты в источниках (прогон по корректорам)…", flush=True)
    Ax, Ay = angle_response(ring, corr, obs, lattice=a.lattice, nsl=a.nsl, rebuild=a.rebuild)
    # фотонные строки: центроид = положение + d·угол
    Rph_x = Pe_x[0] + dists[:, None] * Ax
    Rph_y = Pe_x[1] + dists[:, None] * Ay
    ring.apply_perturb(); ring.close()

    # нормировка веса: w=1 соответствует «фотонная строка весит столько же, сколько строка e-BPM»
    # (масштабируем так, чтобы типичная норма строки совпадала — иначе w безразмерно неинтерпретируем)
    se = np.sqrt((Re_x ** 2).sum(axis=1).mean()); sp = np.sqrt((Rph_x ** 2).sum(axis=1).mean())
    norm = se / sp
    Rph_xn, Rph_yn = Rph_x * norm, Rph_y * norm
    print(f"масштаб фотонной строки к e-строке: ×{norm:.2f} (плечо d увеличивает чувствительность к углу)")

    # ---- КОНДИЦИОНИРОВАНИЕ vs вес w (гор. плоскость) ----
    weights = [0.0, 0.5, 1.0, 2.0, 5.0, 10.0]
    print(f"\n--- обусловленность augment-SVD vs вес фотонных строк w (gx, rcond={a.rcond:g}) ---")
    print(f"{'w':>5} | {'σ_max':>9} | {'σ_min(обр)':>10} | {'cond':>7} | {'‖P‖₂ (ход)':>10} | фот.строк")
    cw = cond_vs_weight(Re_x, Rph_xn, weights, a.rcond)
    for (w, smax, smin, cond, pn) in cw:
        print(f"{w:5.1f} | {smax:9.3e} | {smin:10.3e} | {cond:7.1f} | {pn:10.3e} | {len(gids)}")

    # ---- ВЫЧИСЛИТЕЛЬНАЯ НАГРУЗКА ----
    nc = len(corr)
    cases = [("база: только e-BPM", len(mon), 1),
             ("A: e-BPM + XBPM×1/пучок", len(mon) + len(gids), 1),
             ("A: e-BPM + XBPM×2/пучок", len(mon) + 2 * len(gids), 1)]
    print(f"\n--- ВЫЧИСЛИТЕЛЬНАЯ НАГРУЗКА (одна плоскость; закон θ=P·o — матрично-векторное умножение) ---")
    print(f"корректоров (столбцов P) = {nc};  2 плоскости (x,y) → всё ×2\n")
    print(f"{'конфигурация':<26} | {'строк':>5} | {'P [cols×rows]':>13} | {'t_MVM':>8} | {'память P':>9} | "
          f"{'GFLOP/с @10кГц(2пл)':>20}")
    rates = {2e3: "2 кГц", 1e4: "10 кГц", 1e5: "100 кГц"}
    load_rows = []
    for name, nrow, _ in cases:
        P = np.random.default_rng(1).normal(size=(nc, nrow))     # оператор того же размера (время зависит от размера)
        t = bench_mvm(P, nrow)
        flop_cycle = 2 * nc * nrow                               # mul+add
        mem = P.nbytes / 1024                                    # КиБ (float64)
        g10 = 2 * flop_cycle * 1e4 / 1e9                         # 2 плоскости @10кГц
        print(f"{name:<26} | {nrow:5d} | {nc:4d}×{nrow:<6d} | {t*1e6:6.2f}мкс | {mem:7.1f}КиБ | {g10:18.3f}")
        load_rows.append((name, nrow, t, flop_cycle, mem))

    print(f"\nFLOP/с всей петли (2 плоскости) по тактам — augment ×1 XBPM/пучок ({len(mon)+len(gids)} строк):")
    nrowA = len(mon) + len(gids); flopA = 2 * nc * nrowA
    for rate, lbl in rates.items():
        print(f"   {lbl:>7}: {2*flopA*rate/1e9:7.3f} GFLOP/с   (MVM-бюджет такта {1e6/rate:7.2f} мкс, "
              f"факт {load_rows[1][2]*2*1e6:5.2f} мкс)")

    # SVD/pinv — офлайн, один раз на состояние машины
    Raug = np.vstack([Re_x, Rph_xn])
    t0 = time.perf_counter(); _ = np.linalg.pinv(Raug, rcond=a.rcond); tsvd = time.perf_counter() - t0
    print(f"\nофлайн-пересборка псевдообратной (SVD [{Raug.shape[0]}×{Raug.shape[1]}]): {tsvd*1e3:.1f} мс "
          f"— делается ОДИН раз на состояние машины (смена зазора ИД → плечо/чувствительность), не в такте")

    # ---- график: спектр SV и cond vs w ----
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 4.0))
    sv_e = np.linalg.svd(Re_x, compute_uv=False)
    sv_a = np.linalg.svd(np.vstack([Re_x, Rph_xn]), compute_uv=False)
    ax1.semilogy(range(1, len(sv_e) + 1), sv_e, "o-", ms=3, color="#185fa5", label="только e-BPM")
    ax1.semilogy(range(1, len(sv_a) + 1), sv_a, "s-", ms=3, color="#c0392b", label=f"+XBPM (w=1)")
    ax1.axhline(a.rcond * sv_e[0], ls="--", color="#888", lw=0.8, label=f"порог rcond={a.rcond:g}")
    ax1.set_xlabel("номер моды"); ax1.set_ylabel("сингулярное число"); ax1.set_title("Спектр augment-матрицы")
    ax1.legend(fontsize=8); ax1.grid(alpha=0.25, which="both")
    ww = np.array([c[0] for c in cw]); cc = np.array([c[3] for c in cw]); pp = np.array([c[4] for c in cw])
    ax2.plot(ww, cc, "o-", color="#185fa5", label="число обусловленности")
    ax2b = ax2.twinx(); ax2b.plot(ww, pp, "s--", color="#c0392b", label="‖P‖₂ (ход корректоров)")
    ax2.set_xlabel("вес фотонных строк w"); ax2.set_ylabel("cond", color="#185fa5")
    ax2b.set_ylabel("‖P‖₂", color="#c0392b"); ax2.set_title("Обусловленность и ход vs вес w")
    ax2.grid(alpha=0.25)
    os.makedirs(RESDIR, exist_ok=True)
    png = os.path.join(RESDIR, f"xbpm_svd_L{a.lattice}.png")
    fig.tight_layout(); fig.savefig(png, dpi=120); plt.close(fig)
    print(f"\n[OK] {png}")


if __name__ == "__main__":
    main()
