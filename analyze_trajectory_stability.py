"""
Анализ устойчивости траекторий относительно ОПОРНОЙ ОРБИТЫ.

Опорная (дизайн-) орбита определяется ГЕОМЕТРИЕЙ решётки -- линией, проходящей
через центры всех элементов lattice (из lattice_configs.csv), отсортированных
по instance_id. Это аппроксимация дизайн-траектории, не зависящая от настройки
токов магнитов.

Для каждой частицы (включая референсную P1) считаем:
  A(s) = расстояние от частицы до полилинии-орбиты в момент пути s.

Это даёт ЧЕСТНУЮ оценку отклонения: если P1 сама уплывает, мы это увидим.

Анализ:
  1. A(s) -- мгновенное отклонение от орбиты
  2. Линейный фит:        A(s) = A0 + alpha*s    -- средний рост [мм/м]
  3. Экспоненциальный:    A(s) = A0 * exp(gamma*s)
  4. R^2 фита -- доверие к gamma

Verdict-логика учитывает R^2:
  R^2 < 0.3:  любой gamma считается ненадёжным
  gamma > 0.01 при R^2 >= 0.3:  GROWING
  gamma < -0.01 при R^2 >= 0.3: DAMPING
  A_max > 3*A_init:             OSCILLATING (биения)
  иначе:                         STABLE

Использование:
  python analyze_trajectory_stability.py [exp_id]
  default exp_id=4
"""
import struct
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR / "workspace"


# ---------------- TRAJECTORY LOADING ----------------

def load_trajectories(path):
    parts = []
    with open(path, "rb") as f:
        while True:
            h = f.read(4)
            if not h:
                break
            n = struct.unpack("I", h)[0]
            pts = np.fromfile(f, dtype=np.float64, count=n*3).reshape(-1, 3)
            parts.append(pts)
    return parts


def cumulative_arc_length(traj):
    if len(traj) < 2:
        return np.array([0.0])
    diffs = np.diff(traj, axis=0)
    seg = np.linalg.norm(diffs, axis=1)
    return np.concatenate([[0.0], np.cumsum(seg)])


# ---------------- DESIGN ORBIT ----------------

def load_design_orbit(exp_id):
    """Возвращает (positions, yaws, first_tangent, last_tangent, lat_id), где
       positions      -- (N, 3) центры элементов lattice по instance_id;
       yaws           -- (N,) yaw каждого элемента (для арочной интерполяции);
       first_tangent  -- касательная у первого элемента (для экстраполяции назад,
                         в участок старта до первого магнита);
       last_tangent   -- касательная у последнего элемента (для экстраполяции
                         вперёд, в drift-секцию за концом lattice)."""
    # 1. exp_id -> lattice_id
    exp_df = pd.read_csv(WS / "tables" / "experiments.csv", sep=";", skipinitialspace=True)
    row = exp_df[exp_df["exp_id"] == exp_id]
    if row.empty:
        raise ValueError(f"experiment {exp_id} not found in experiments.csv")
    lat_id = int(row["lattice_id"].iloc[0])

    # 2. lattice_id -> упорядоченные позиции и yaw элементов
    lat_df = pd.read_csv(WS / "tables" / "lattice_configs.csv", sep=";", skipinitialspace=True)
    sub = lat_df[lat_df["lattice_id"] == lat_id].sort_values("instance_id").reset_index(drop=True)
    if sub.empty:
        raise ValueError(f"lattice {lat_id} has no rows in lattice_configs.csv")
    positions = sub[["x", "y", "z"]].to_numpy()
    yaws = sub["yaw"].to_numpy()
    first_tan = np.array([np.cos(yaws[0]),  np.sin(yaws[0]),  0.0])
    last_tan  = np.array([np.cos(yaws[-1]), np.sin(yaws[-1]), 0.0])
    return positions, yaws, first_tan, last_tan, lat_id


def refine_polyline_with_arcs(positions, yaws, n_subdivide=10):
    """Сглаживает полилинию: между каждыми двумя соседними элементами вставляет
    n_subdivide-1 точек.
    - Если yaw одинаков (drift) -- линейная интерполяция.
    - Если yaw меняется -- дуга окружности, касательная к yaw_i в P_i и
      к yaw_{i+1} в P_{i+1}.

    Это сокращает «sagitta-артефакт» полилинии в области изгибов (диполей)."""
    if len(positions) < 2 or n_subdivide < 1:
        return positions
    refined = [positions[0]]
    for i in range(len(positions) - 1):
        p_a = positions[i]
        p_b = positions[i+1]
        yaw_a = yaws[i]
        yaw_b = yaws[i+1]
        dyaw = yaw_b - yaw_a
        # Нормализация к [-pi, pi] на случай wrap-around
        while dyaw >  np.pi: dyaw -= 2*np.pi
        while dyaw < -np.pi: dyaw += 2*np.pi

        chord = p_b - p_a
        chord_len = np.linalg.norm(chord)

        if abs(dyaw) < 1e-6 or chord_len < 1e-12:
            # Прямой сегмент: линейная интерполяция
            for k in range(1, n_subdivide + 1):
                t = k / n_subdivide
                refined.append(p_a + t * chord)
        else:
            # Дуга: вычисляем kappa = curvature и параметризуем по s = arc length
            arc_len = chord_len * abs(dyaw) / (2 * np.sin(abs(dyaw)/2))
            kappa = dyaw / arc_len
            sin_a, cos_a = np.sin(yaw_a), np.cos(yaw_a)
            for k in range(1, n_subdivide + 1):
                s = k * arc_len / n_subdivide
                yaw_s = yaw_a + kappa * s
                x_off = (np.sin(yaw_s) - sin_a) / kappa
                y_off = (cos_a - np.cos(yaw_s)) / kappa
                z_off = (p_b[2] - p_a[2]) * (k / n_subdivide)
                refined.append(np.array([p_a[0] + x_off,
                                          p_a[1] + y_off,
                                          p_a[2] + z_off]))
    return np.array(refined)


def extend_polyline(polyline, first_tan, last_tan,
                    pre_length=10.0, post_length=20.0):
    """Дописывает по 1 точке с каждого конца:
       перед   = first_point - first_tan * pre_length
       после   = last_point  + last_tan  * post_length
    Этим закрываем prelude (от точки старта до первого магнита) и
    дрейф-секцию за последним элементом."""
    first = polyline[0]
    last  = polyline[-1]
    pre  = first - first_tan * pre_length
    post = last  + last_tan  * post_length
    return np.vstack([pre[np.newaxis, :], polyline, post[np.newaxis, :]])


# ---------------- DISTANCE TO POLYLINE ----------------

def distance_to_polyline(points, polyline, batch=512):
    """Для каждой точки в `points` (M, 3) считает мин расстояние до полилинии
    `polyline` (N, 3). Полилиния = ломаная из (N-1) отрезков.

    Векторизованный батчами вариант, чтобы вместиться в память для M*N ~ 10M."""
    M = len(points)
    N = len(polyline)
    if N < 2:
        return np.full(M, np.inf)
    starts = polyline[:-1]              # (N-1, 3)
    dirs   = polyline[1:] - polyline[:-1]  # (N-1, 3)
    seg_len2 = np.maximum(np.sum(dirs * dirs, axis=1), 1e-30)  # (N-1,)
    out = np.zeros(M)
    for i0 in range(0, M, batch):
        i1 = min(i0 + batch, M)
        P = points[i0:i1, np.newaxis, :]              # (b, 1, 3)
        rel = P - starts[np.newaxis, :, :]            # (b, N-1, 3)
        t = np.sum(rel * dirs[np.newaxis, :, :], axis=2) / seg_len2[np.newaxis, :]
        t = np.clip(t, 0.0, 1.0)                      # (b, N-1)
        closest = starts[np.newaxis, :, :] + t[..., np.newaxis] * dirs[np.newaxis, :, :]
        d = np.linalg.norm(P - closest, axis=2)        # (b, N-1)
        out[i0:i1] = d.min(axis=1)
    return out


# ---------------- ENVELOPE FIT & VERDICT ----------------

def fit_envelope(s, a):
    """Возвращает (alpha_lin, gamma_exp, r2_lin, r2_exp)."""
    if len(s) < 5:
        return 0.0, 0.0, 0.0, 0.0
    p = np.polyfit(s, a, 1)
    alpha = p[0]
    pred = np.polyval(p, s)
    ss_tot = np.sum((a - a.mean()) ** 2)
    r2_lin = 1.0 - np.sum((a - pred) ** 2) / ss_tot if ss_tot > 1e-30 else 0.0

    eps = max(1e-9, a.max() * 1e-6)
    mask = a > eps
    if mask.sum() < 5:
        return alpha, 0.0, r2_lin, 0.0
    pe = np.polyfit(s[mask], np.log(a[mask]), 1)
    gamma = pe[0]
    pred_e = np.exp(np.polyval(pe, s[mask]))
    ss_tot_e = np.sum((a[mask] - a[mask].mean()) ** 2)
    r2_exp = 1.0 - np.sum((a[mask] - pred_e) ** 2) / ss_tot_e if ss_tot_e > 1e-30 else 0.0
    return alpha, gamma, r2_lin, r2_exp


def verdict(a_init, a_final, a_max, gamma, r2_exp):
    """Учитывает R^2: при низком фите gamma считается ненадёжным.
    Особый случай: a_init ~ 0 (референсная частица на орбите) -- любой
    значимый рост = DRIFT, а не осцилляция."""
    # Спецкейс для частицы, стартующей строго на орбите (P1).
    # Тогда «осцилляция вокруг 0» бессмысленна -- любое отклонение это drift.
    if a_init < 1e-4:  # <0.1 mm -- de facto на орбите
        if a_final > 1e-3:    # >1 мм в конце -- ушла
            return "DRIFTED"
        if a_max > 5e-3:      # >5 мм в середине -- биения, но вернулась
            return "OSCILLATING"
        return "STABLE"
    # Обычный случай для спутника (a_init ~ 1 мм)
    a_eff = a_init
    if r2_exp >= 0.3 and gamma > 0.01:
        return "GROWING"
    if r2_exp >= 0.3 and gamma < -0.01:
        return "DAMPING"
    # DRIFTED: конец большой И БЛИЗОК к максимуму (= не вернулась к норме).
    # Если A_final < 0.5*A_max, то частица БЫЛА далеко в середине, но
    # вернулась -- это OSCILLATING, не DRIFTED.
    if a_final > 2.0 * a_eff and a_final > 0.5 * a_max:
        return "DRIFTED"
    # OSCILLATING: пик в середине, но возврат к норме
    if a_max > 3.0 * a_eff:
        return "OSCILLATING"
    return "STABLE"


# ---------------- MAIN ----------------

def analyze(exp_id, show=True):
    traj_path = WS / "results" / f"trajectories_exp{exp_id}.bin"
    if not traj_path.exists():
        print(f"ERROR: {traj_path} not found.")
        print(f"  Run first: build\\SynchrotronTracker.exe {exp_id}")
        sys.exit(1)

    # Опорная орбита: загрузить, сгладить арками, добавить prelude/postlude
    positions, yaws, first_tan, last_tan, lat_id = load_design_orbit(exp_id)
    N_SUBDIVIDE = 10  # точек между соседними элементами
    polyline = refine_polyline_with_arcs(positions, yaws, n_subdivide=N_SUBDIVIDE)
    polyline_ext = extend_polyline(polyline, first_tan, last_tan,
                                   pre_length=10.0, post_length=20.0)
    print(f"Experiment {exp_id} (lattice {lat_id}):")
    print(f"  raw elements: {len(positions)}, refined polyline: {len(polyline)} "
          f"(x{N_SUBDIVIDE} subdivision), with extensions: {len(polyline_ext)}")
    print(f"  first tangent: ({first_tan[0]:+.4f}, {first_tan[1]:+.4f}, {first_tan[2]:+.4f})")
    print(f"  last  tangent: ({last_tan[0]:+.4f}, {last_tan[1]:+.4f}, {last_tan[2]:+.4f})")

    # Траектории
    parts = load_trajectories(traj_path)
    print(f"  particles loaded:      {len(parts)} (lengths: {[len(p) for p in parts]})")
    if not parts:
        return

    labels = ["P1 (ref)", "P2 (+Y)", "P3 (-Y)", "P4 (+Z)", "P5 (-Z)"]
    colors = ["#9467bd", "#1f77b4", "#2ca02c", "#ff7f0e", "#d62728"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)

    print()
    print(f"{'Particle':<10}  {'A_init [mm]':>11}  {'A_final [mm]':>12}"
          f"  {'A_max [mm]':>11}  {'alpha[mm/m]':>11}  {'gamma_exp[1/m]':>14}"
          f"  {'R2_exp':>7}  {'verdict':>13}")
    print("-" * 113)

    rows = []
    for i, traj in enumerate(parts):
        if len(traj) < 5:
            print(f"{labels[i] if i < len(labels) else f'P{i+1}':<10}  too few points")
            continue
        s_loc = cumulative_arc_length(traj)
        dist  = distance_to_polyline(traj, polyline_ext)
        a_init = dist[0]
        a_final = dist[-1]
        a_max = dist.max()
        alpha, gamma, r2_lin, r2_exp = fit_envelope(s_loc, dist)
        v = verdict(a_init, a_final, a_max, gamma, r2_exp)
        lab = labels[i] if i < len(labels) else f"P{i+1}"
        rows.append((lab, a_init, a_final, a_max, alpha, gamma, r2_exp, v))

        print(f"{lab:<10}  {a_init*1000:>11.4f}  {a_final*1000:>12.4f}"
              f"  {a_max*1000:>11.4f}  {alpha*1000:>11.3f}  {gamma:>14.5f}"
              f"  {r2_exp:>7.3f}  {v:>13s}")

        color = colors[i % len(colors)]
        ax1.plot(s_loc, dist*1000, color=color, lw=1.0, alpha=0.85,
                 label=f"{lab}: gamma={gamma:.4f}/m, R2={r2_exp:.2f} ({v})")

        ax2.semilogy(s_loc, np.maximum(dist, 1e-9)*1000, color=color, lw=1.0, alpha=0.85)

    ax1.set_ylabel("Distance to design orbit (mm)")
    ax1.set_title(f"Trajectory stability vs DESIGN ORBIT -- experiment {exp_id} (lattice {lat_id})")
    ax1.legend(loc="upper left", fontsize=8)
    ax1.grid(True, alpha=0.3)

    ax2.set_xlabel("Path length s (m)")
    ax2.set_ylabel("|particle - orbit| (mm, log)")
    ax2.grid(True, which="both", alpha=0.3)

    out_path = WS / "results" / f"trajectory_stability_exp{exp_id}.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    print(f"\nPlot saved: {out_path}")

    # ============== FFT-АНАЛИЗ БЕТАТРОННЫХ ЧАСТОТ ==============
    # Для каждой частицы делаем FFT от A(s) - mean(A(s)) (центрируем).
    # Получаем спектр мощности |F(k)|^2 по волновому числу k [1/m].
    # Переводим в "бетатронный тюн" Q = k * SP_length (циклов на SP).
    print()
    print("=== FFT analysis (betatron oscillation frequencies) ===")
    # SP_length вычисляем из layout: расстояние SP1.first -> SP2.first
    # (тождественно длине одного супериода для стандартной части кольца)
    try:
        layout = pd.read_csv(WS / "tables" / "lattice_layout.csv", sep=";", skipinitialspace=True)
        sp1f = layout[layout["global_id"] == 1].iloc[0]
        sp2f = layout[layout["global_id"] == 70].iloc[0]
        sp_length = float(np.linalg.norm([sp2f["x"]-sp1f["x"], sp2f["y"]-sp1f["y"], sp2f["z"]-sp1f["z"]]))
    except Exception:
        sp_length = 27.71  # fallback
    fig_fft, axf = plt.subplots(figsize=(13, 5))
    for i, traj in enumerate(parts):
        if len(traj) < 64:
            continue
        s_loc = cumulative_arc_length(traj)
        dist  = distance_to_polyline(traj, polyline_ext)
        # Унифицируем шаг s (FFT требует равномерной сетки)
        # На самом деле s уже почти uniform т.к. |v| const (Boris сохраняет)
        ds = (s_loc[-1] - s_loc[0]) / (len(s_loc) - 1)
        signal = dist - dist.mean()  # центрируем (убираем DC)
        # Hann window для уменьшения spectral leakage
        window = np.hanning(len(signal))
        signal_w = signal * window
        spec = np.fft.rfft(signal_w)
        freqs = np.fft.rfftfreq(len(signal_w), d=ds)  # cycles per meter
        power = np.abs(spec)
        # Топ-3 пика по мощности
        valid = (freqs > 0.001)  # отбрасываем DC
        if not valid.any():
            continue
        sorted_idx = np.argsort(-power * valid)[:5]
        top = [(freqs[k], power[k]) for k in sorted_idx if power[k] > 0]
        lab = labels[i] if i < len(labels) else f"P{i+1}"
        peak_str = ", ".join(f"f={f:.4f}/m (Q={f*sp_length:.3f}/SP)" for f, _ in top[:3])
        print(f"  {lab}: top freqs: {peak_str}")
        color = colors[i % len(colors)]
        axf.plot(freqs * sp_length, power, color=color, lw=0.8, alpha=0.85,
                 label=f"{lab}")
    # Метки рациональных тюнов (хорошие/плохие резонансы)
    for q in [0.25, 0.333, 0.5, 0.667, 0.75]:
        axf.axvline(q, color="gray", lw=0.4, alpha=0.4, ls="--")
        axf.text(q, axf.get_ylim()[1]*0.95, f"{q:.2f}", fontsize=7,
                 color="gray", rotation=90, va="top")
    axf.set_xlabel("Tune Q (cycles per superperiod)")
    axf.set_ylabel("|FFT(A - <A>)|  (arb units)")
    axf.set_title(f"Betatron frequency spectrum  -- experiment {exp_id}\n"
                  f"(Hann window; dashed lines = simple rationals)")
    axf.set_xlim(0, 1.5)
    axf.legend(fontsize=8)
    axf.grid(True, alpha=0.3)
    fft_path = WS / "results" / f"betatron_spectrum_exp{exp_id}.png"
    plt.tight_layout()
    plt.savefig(fft_path, dpi=120)
    print(f"FFT plot saved: {fft_path}")

    # Резюме
    n_growing = sum(1 for r in rows if r[-1] == "GROWING")
    n_damping = sum(1 for r in rows if r[-1] == "DAMPING")
    n_drifted = sum(1 for r in rows if r[-1] == "DRIFTED")
    n_osc     = sum(1 for r in rows if r[-1] == "OSCILLATING")
    n_stable  = sum(1 for r in rows if r[-1] == "STABLE")
    print()
    print(f"SUMMARY: growing={n_growing}, damping={n_damping}, drifted={n_drifted},"
          f" oscillating={n_osc}, stable={n_stable}")
    if n_growing or n_drifted:
        print("  -> Решётка показывает неустойчивость / уплытие от опорной орбиты.")
        print("     Рекомендуется усилить envelope_loss или добавить дополнительные")
        print("     контрольные точки для оптимизатора.")
    elif n_osc and not n_stable:
        print("  -> Сильные бетатронные биения. Огибающая ограничена, но широкая.")
    elif n_stable == len([r for r in rows if r[-1] != "DAMPING"]):
        print("  -> Все частицы стабильны относительно опорной орбиты.")

    if show:
        plt.show()


if __name__ == "__main__":
    exp_id = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    analyze(exp_id)
