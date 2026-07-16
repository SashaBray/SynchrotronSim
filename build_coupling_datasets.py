# -*- coding: utf-8 -*-
"""Датасеты «ток → подобранные физкоэффициенты» для ВСЕХ элементов (по аналогии elegant_dq_dataset)
и перерасчёт матрицы связи по ним. Через КОРРЕКТНУЮ модель element_bench (осевой кик matrix_tracker,
дрейф по L_phys — БЕЗ бага dsg=2·half, искажавшего старый elegant_dq).

Для каждого значения тока: диффур (Boris по картам, эталон) + 5-м хвост; least_squares подбирает
активные физкоэф. (B/G/S/O), чтобы КИК повторил траекторию диффура на probe-пучке ±4мм. Датасет —
ток → коэф + невязка. Матрица связи per_* = наклон коэф(ток) линейной регрессией через 0.

Одиночные обмотки (quad/sext/oct/корректор) — свип 1 тока; DQ (3 обмотки) — свип по каждой отдельно
(для матрицы связи per-обмотка) + полный 3D-датасет. Выход: workspace/research/coupling_datasets/.
Запуск: py -3.14 build_coupling_datasets.py [--no-dq3d]"""
import os, csv, sys
import numpy as np
from scipy.optimize import least_squares
import element_bench as eb
import matrix_tracker as mt
from bip_matrix import QP

OUT = os.path.join("workspace", "research", "coupling_datasets"); os.makedirs(OUT, exist_ok=True)
TB = os.path.join("workspace", "tables")
ENTRY, EXIT = 0.5, 10.0
# нормальные элементы — probe только по u (сигнал в гор. плоскости); скью-элементы (корректоры) —
# ещё по v и по диагонали, иначе кросс-плоскостные скью-мультиполи (Gs/Ss/Os) не возбуждаются
PROBES = [[0, 0, 0, 0], [2e-3, 0, 0, 0], [-2e-3, 0, 0, 0], [4e-3, 0, 0, 0], [-4e-3, 0, 0, 0]]
PROBES_SKEW = PROBES + [[0, 2e-3, 0, 0], [0, -2e-3, 0, 0], [0, 4e-3, 0, 0], [0, -4e-3, 0, 0],
                        [3e-3, 3e-3, 0, 0], [-3e-3, 3e-3, 0, 0], [3e-3, -3e-3, 0, 0], [-3e-3, -3e-3, 0, 0]]
GRID = list(range(-20, 21, 5))                          # УЗКИЙ ОДЗ ток ±20 А, шаг 5 (ближе к рабочим 1–4 А)
# КВАДРУПОЛИ (новые карты сильные, рабочий ток ~1–2 А): при |I|≫ рабочего дефокус-направление
# разносит пучок по 10-м хвосту, а на фокус-стороне толстая линза занижает per_G. ±20 А — тот же ОДЗ.
GRID_QUAD = list(range(-20, 21, 5))                     # ±20 А
SCALE = {"B": 1e-3, "Bs": 1e-3, "G": 1e-1, "Gs": 1e-1, "S": 1e0, "Ss": 1e0, "O": 1e1, "Os": 1e1}
PERCOL = {"B": "per_B_Tl", "Bs": "per_Bs_Tl", "G": "per_G_Tlm", "Gs": "per_Gs_Tlm",
          "S": "per_S_Tlm2", "Ss": "per_Ss_Tlm2", "O": "per_O_Tlm3", "Os": "per_Os_Tlm3"}
CM = {r["device_type_id"]: r for r in mt._read("coupling_matrix.csv")}
LEFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}

# (имя, обмотки, активные коэф., режим). Нормальные линзы (quad/sext/oct) и DQ — хвост-фит (tail).
# Корректор — 5-обм. (main+c1..c4): режим "mp" = ОДНОПРОХОДНОЕ мультипольное разложение (все 8 коэф.
# B/Bs/G/Gs/S/Ss/O/Os сразу, линейным МНК по сетке probe u×v). Хвост-фит для 8 коэф. вырожден
# (октуполь ∝r³ тонет, много комбинаций дают ту же траекторию); однопроходное разложение —
# корректно обусловлено и совпадает с полевым фитом v3 в coupling_matrix. Синтетика L100 (есть .bin).
CORR_ACT = ["B", "Bs", "G", "Gs", "S", "Ss", "O", "Os"]
SPEC = [
    ("QP_L_1A", ["QP_L_1A"], ["G"], "1d"), ("QP_sh_1A", ["QP_sh_1A"], ["G"], "1d"),
    ("SP_I1_1A", ["SP_I1_1A"], ["S"], "1d"), ("SP_I2_1A", ["SP_I2_1A"], ["S"], "1d"),
    ("OP_1A", ["OP_1A"], ["O"], "1d"),
    ("corrector_main", ["corrector_main_100A_L100"], CORR_ACT, "mp"),
    ("corrector_c1", ["corrector_c1_100A_L100"], CORR_ACT, "mp"),
    ("corrector_c2", ["corrector_c2_100A_L100"], CORR_ACT, "mp"),
    ("corrector_c3", ["corrector_c3_100A_L100"], CORR_ACT, "mp"),
    ("corrector_c4", ["corrector_c4_100A_L100"], CORR_ACT, "mp"),
    ("DQ", ["DQ_I1_1A", "DQ_I2_1A", "DQ_I3_1A"], ["B", "G", "S"], "dq"),
]
SKEW = {"Bs", "Gs", "Ss", "Os"}                         # наличие → probe по u+v (иначе не возбудить)
# ---- однопроходное мультипольное разложение (режим "mp"): сетка probe и базис кик-модели ----
MPROBES = [(u, v) for u in (-4e-3, -2e-3, 0.0, 2e-3, 4e-3) for v in (-4e-3, -2e-3, 0.0, 2e-3, 4e-3)]


def mp_basis(u, v):
    """Строки базиса кик-модели (интеграл.) для up, vp в порядке CORR_ACT [B,Bs,G,Gs,S,Ss,O,Os].
    up/(QP·Leff) = B − G·u − ½S(u²−v²) − (O/6)(u³−3uv²) + Gs·v + Ss·uv + (Os/6)(3u²v−v³)
    vp/(QP·Leff) = Bs + G·v + S·uv + (O/6)(3u²v−v³) + Gs·u + ½Ss(u²−v²) + (Os/6)(u³−3uv²)"""
    up = [1.0, 0.0, -u, v, -0.5 * (u * u - v * v), u * v, -(u**3 - 3 * u * v * v) / 6, (3 * u * u * v - v**3) / 6]
    vp = [0.0, 1.0, v, u, u * v, 0.5 * (u * u - v * v), (3 * u * u * v - v**3) / 6, (u**3 - 3 * u * v * v) / 6]
    return up, vp


def defl_slope(windings, u0, v0):
    """Однопроходный наклон (up,vp) после элемента — из выходного дрейфа track_diff (наклон постоянен)."""
    pts = eb.track_diff(windings, ENTRY, EXIT, [u0, v0, 0.0, 0.0])
    (s2, u2, v2), (s1, u1, v1) = pts[-1], pts[-2]
    return (u2 - u1) / (s2 - s1), (v2 - v1) / (s2 - s1)


def extract_mp(dev, cur, Leff):
    """Все 8 мультиполей одной обмотки при токе cur — линейный МНК Boris-отклонений по сетке u×v."""
    if cur == 0:
        return {c: 0.0 for c in CORR_ACT}, 0.0
    A, b = [], []
    for (u, v) in MPROBES:
        up, vp = defl_slope({dev: float(cur)}, u, v)
        ru, rv = mp_basis(u, v)
        A += [ru, rv]; b += [up / (QP * Leff), vp / (QP * Leff)]
    A = np.array(A); b = np.array(b)
    x, *_ = np.linalg.lstsq(A, b, rcond=None)
    rms = float(np.sqrt(np.mean(((A @ x - b) * QP * Leff) ** 2))) * 1e6   # невязка отклонения, мкрад
    return dict(zip(CORR_ACT, x)), rms


def sweep_mp(name, dev):
    """Датасет ток→8 мультиполей (однопроходное разложение). per_* = наклон через 0."""
    Leff = LEFF[dev]; rows = []
    for cur in GRID:
        coeffs, rms = extract_mp(dev, cur, Leff)
        rows.append([cur] + [coeffs[c] for c in CORR_ACT] + [rms])
    with open(os.path.join(OUT, f"{name}.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(["I_A"] + [f"{c}_phys" for c in CORR_ACT] + ["fit_urad"])
        for r in rows:
            w.writerow([r[0]] + [f"{v:.6e}" for v in r[1:-1]] + [f"{r[-1]:.2f}"])
    per = {c: slope0([r[0] for r in rows], [r[1 + i] for r in rows]) for i, c in enumerate(CORR_ACT)}
    return per, rows


def x0_coupling(windings, active):
    return [sum(float(CM[d].get(PERCOL[c]) or 0) * cur for d, cur in windings.items()) for c in active]


def fit(windings, active, probes=PROBES):
    """Подбор active-коэф.: КИК (осевой, element_bench) повторяет ДИФФУР на хвосте. → (coeffs, rms_мм)."""
    diffs = [eb.track_diff(windings, ENTRY, EXIT, p)[-1] for p in probes]

    def resid(vals):
        c = dict(zip(active, vals)); r = []
        for p, d in zip(probes, diffs):
            try:
                k = eb.track_kick(windings, c, ENTRY, EXIT, p)[-1]
                r += [(k[1] - d[1]) * 1e3, (k[2] - d[2]) * 1e3]
            except Exception:
                r += [1e6, 1e6]
        return r

    x0 = x0_coupling(windings, active)
    sol = least_squares(resid, x0, x_scale=[SCALE[c] for c in active], method="trf")
    return dict(zip(active, sol.x)), float(np.sqrt(np.mean(np.array(sol.fun) ** 2)))


def slope0(xs, ys):
    """Наклон через 0 (per_* = коэф/ток): Σxy/Σx²."""
    xs = np.array(xs, float); ys = np.array(ys, float)
    return float(np.sum(xs * ys) / (np.sum(xs * xs) + 1e-30))


def sweep_1d(name, dev, active, grid=GRID):
    probes = PROBES_SKEW if (set(active) & SKEW) else PROBES     # скью-набор возбуждают только u+v probe
    rows = []
    for cur in grid:
        (coeffs, rms) = fit({dev: float(cur)}, active, probes) if cur != 0 else ({a: 0.0 for a in active}, 0.0)
        rows.append([cur] + [coeffs[a] for a in active] + [rms])
    hdr = ["I_A"] + [f"{a}_phys" for a in active] + ["fit_rms_mm"]
    with open(os.path.join(OUT, f"{name}.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr)
        for r in rows:
            w.writerow([r[0]] + [f"{v:.6e}" for v in r[1:-1]] + [f"{r[-1]:.4f}"])
    per = {a: slope0([r[0] for r in rows], [r[1 + i] for r in rows]) for i, a in enumerate(active)}
    return per, rows


def sweep_dq(windings, active, do_3d=True):
    """DQ: по каждой обмотке отдельно (per-обмотка матрица связи) + опц. полный 3D-датасет."""
    per = {}                                                       # per[dev][coef]
    single = {}
    for d in windings:                                            # свип одной обмотки, остальные 0
        rows = []
        for cur in GRID:
            w = {x: (float(cur) if x == d else 0.0) for x in windings}
            (coeffs, rms) = fit(w, active) if cur != 0 else ({a: 0.0 for a in active}, 0.0)
            rows.append([cur] + [coeffs[a] for a in active] + [rms])
        single[d] = rows
        per[d] = {a: slope0([r[0] for r in rows], [r[1 + i] for r in rows]) for i, a in enumerate(active)}
    with open(os.path.join(OUT, "DQ_single_windings.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(["winding", "I_A"] + [f"{a}_phys" for a in active] + ["fit_rms_mm"])
        for d, rows in single.items():
            for r in rows:
                w.writerow([d, r[0]] + [f"{v:.6e}" for v in r[1:-1]] + [f"{r[-1]:.4f}"])
    if do_3d:
        with open(os.path.join(OUT, "DQ_full3d.csv"), "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f, delimiter=";"); w.writerow(["i1_A", "i2_A", "i3_A"] + [f"{a}_phys" for a in active] + ["fit_rms_mm"])
            for i1 in GRID:
                for i2 in GRID:
                    for i3 in GRID:
                        wd = {"DQ_I1_1A": float(i1), "DQ_I2_1A": float(i2), "DQ_I3_1A": float(i3)}
                        (coeffs, rms) = fit(wd, active) if any((i1, i2, i3)) else ({a: 0.0 for a in active}, 0.0)
                        w.writerow([i1, i2, i3] + [f"{coeffs[a]:.6e}" for a in active] + [f"{rms:.4f}"])
    return per


def main():
    do_3d = "--no-dq3d" not in sys.argv
    newper = {}                                                   # dev → {coef: per}
    print(f"{'элемент':16} {'коэф':>4} {'per_* (наклон)':>15}  среднее невязки, мм")
    for name, windings, active, mode in SPEC:
        if mode == "1d":
            dev = windings[0]; grid = GRID_QUAD if active == ["G"] else GRID   # квадруполи — узкий диапазон
            per, rows = sweep_1d(name, dev, active, grid)
            newper[dev] = per
            rms = np.mean([r[-1] for r in rows])
            for a in active:
                print(f"  {name:14} {a:>4} {per[a]:>+15.5e}  {rms:.3f}")
        elif mode == "mp":                                       # корректор: однопроходное 8-мультип. разложение
            dev = windings[0]; per, rows = sweep_mp(name, dev)
            newper[dev] = per
            urad = np.mean([r[-1] for r in rows])
            print(f"  {name:14} {'/'.join(CORR_ACT)}  невязка {urad:.2f} мкрад")
            print("    " + "  ".join(f"{c}={per[c]:+.2e}" for c in CORR_ACT))
        else:                                                    # DQ
            per = sweep_dq(windings, active, do_3d)
            for d in windings:
                newper[d] = per[d]
                print(f"  {d:14} B/G/S {per[d]['B']:+.3e}/{per[d]['G']:+.3e}/{per[d]['S']:+.3e}")
    # СРАВНЕНИЕ per_* из датасета (наклон через 0) с ТЕКУЩЕЙ матрицей связи (near-axis recalibrate).
    # НЕ перезаписываем coupling_matrix: наклон хвост-фита СМЕЩЁН у сильных элементов (толстая линза),
    # near-axis recalibrate — корректная линейная связь. Файл — для наглядности.
    comp = os.path.join(OUT, "coupling_from_datasets.csv")
    with open(comp, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["device", "coef", "per_датасет_наклон", "per_матрица_near_axis", "отношение"])
        for dev, per in newper.items():
            for coef, val in per.items():
                cur = float(CM.get(dev, {}).get(PERCOL[coef]) or 0)
                w.writerow([dev, coef, f"{val:.5e}", f"{cur:.5e}", f"{val/cur:.2f}" if cur else "-"])
    print(f"\n[OK] датасеты → {OUT}/ ; сравнение per_* → {comp} (coupling_matrix НЕ тронут)")


if __name__ == "__main__":
    main()
