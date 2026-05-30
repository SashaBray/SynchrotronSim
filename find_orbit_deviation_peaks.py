"""
Диагностика: где конкретно траектория P1 отклоняется от дизайн-орбиты?

1. Считаем A(s) = расстояние от P1 до полилинии вдоль трассы.
2. Находим локальные максимумы.
3. Для каждого пика определяем ближайший lattice-элемент (по расстоянию
   в 3D пространстве) и печатаем (gid, тип, позиция, отклонение).
4. Показывает зависимость A(s) с метками магнитов на пиках.

Использование:
    python find_orbit_deviation_peaks.py [exp_id]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_trajectory_stability as ats

WS = ats.WS


def find_local_maxima(a, min_separation=10, top_n=15):
    """Возвращает индексы локальных максимумов A(s).
       min_separation -- мин дистанция между пиками (в индексах массива)
       top_n          -- вернуть топ-N по величине."""
    n = len(a)
    peaks = []
    for i in range(1, n - 1):
        if a[i] > a[i-1] and a[i] > a[i+1]:
            peaks.append((i, a[i]))
    # Сортируем по убыванию амплитуды
    peaks.sort(key=lambda x: -x[1])
    # Жадно выбираем top_n с учётом min_separation
    selected = []
    for idx, val in peaks:
        if all(abs(idx - j) >= min_separation for j, _ in selected):
            selected.append((idx, val))
        if len(selected) >= top_n:
            break
    return sorted(selected, key=lambda x: x[0])  # вернуть в порядке s


def find_nearest_element(point, lat_df):
    """Возвращает строку lat_df (как dict) -- ближайший элемент по 3D-расстоянию."""
    dx = lat_df["x"].to_numpy() - point[0]
    dy = lat_df["y"].to_numpy() - point[1]
    dz = lat_df["z"].to_numpy() - point[2]
    d2 = dx*dx + dy*dy + dz*dz
    i = int(np.argmin(d2))
    return lat_df.iloc[i], float(np.sqrt(d2[i]))


def analyze(exp_id):
    # 1. Опорная орбита (refined arc-interpolated)
    positions, yaws, first_tan, last_tan, lat_id = ats.load_design_orbit(exp_id)
    polyline = ats.refine_polyline_with_arcs(positions, yaws, n_subdivide=10)
    polyline_ext = ats.extend_polyline(polyline, first_tan, last_tan, 10.0, 20.0)

    # 2. Траектория P1
    parts = ats.load_trajectories(WS / "results" / f"trajectories_exp{exp_id}.bin")
    if not parts:
        print("No trajectories")
        return
    ref = parts[0]

    # 3. A(s) для P1
    s = ats.cumulative_arc_length(ref)
    A = ats.distance_to_polyline(ref, polyline_ext)
    print(f"P1 trajectory: {len(ref)} points, total path {s[-1]:.2f} m")
    print(f"A_max = {A.max()*1000:.2f} mm at s = {s[int(np.argmax(A))]:.2f} m, "
          f"position {ref[int(np.argmax(A))]}")

    # 4. Локальные максимумы
    print()
    print(f"Top local maxima of A(s) (P1 deviation from design orbit):")
    print(f"{'#':<3}  {'s [m]':>8}  {'A [mm]':>8}  {'P1 position (x,y,z) [m]':<35}  "
          f"{'nearest elem (gid, type, sp)':<40}  {'dist [m]':>9}")
    print("-" * 130)

    # Загрузить таблицы для поиска ближайшего магнита
    lat_df = pd.read_csv(WS / "tables" / "lattice_configs.csv", sep=";", skipinitialspace=True)
    lat_sub = lat_df[lat_df["lattice_id"] == lat_id].sort_values("instance_id").reset_index(drop=True)
    layout_df = pd.read_csv(WS / "tables" / "lattice_layout.csv", sep=";", skipinitialspace=True)
    # layout даёт sp_number и sp_element_id

    peaks = find_local_maxima(A, min_separation=int(len(A) * 0.01), top_n=20)
    rows = []
    for k, (i, val) in enumerate(peaks):
        pos = ref[i]
        elem, dist = find_nearest_element(pos, lat_sub)
        # Найти соответствующий элемент в layout
        iid = int(elem["instance_id"])
        layout_row = layout_df[layout_df["global_id"] == iid]
        sp_n = int(layout_row["sp_number"].iloc[0]) if not layout_row.empty else -1
        sp_eid = int(layout_row["sp_element_id"].iloc[0]) if not layout_row.empty else -1
        dev = elem['device_type_id']
        rows.append((k+1, s[i], val*1000, pos, iid, dev, sp_n, sp_eid, dist))
        print(f"{k+1:<3}  {s[i]:>8.2f}  {val*1000:>8.2f}  "
              f"({pos[0]:+7.3f}, {pos[1]:+7.3f}, {pos[2]:+6.3f})  "
              f"gid={iid:>4} {dev:<18} SP{sp_n:>2}.{sp_eid:>2}  {dist:>9.4f}")

    # 5. Plot
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(s, A*1000, color="purple", lw=1.0, alpha=0.85, label="A(s) = |P1 - design orbit|")
    for k, (i, val) in enumerate(peaks):
        ax.scatter([s[i]], [val*1000], color="red", s=30, zorder=5)
        ax.annotate(f"#{k+1}", (s[i], val*1000), xytext=(3, 3),
                    textcoords="offset points", fontsize=8)
    # Метки границ SP (по yaw элементов)
    yaws_arr = yaws
    # Поместим метку при переходе SP (yaw делает прыжок назад или большой шаг)
    # Проще: разделить по sp_number из layout
    sp_layout = layout_df[layout_df["global_id"] <= len(positions)]
    sp_changes = sp_layout.groupby("sp_number").first().reset_index()
    s_at_elements = ats.cumulative_arc_length(positions)
    for _, row in sp_changes.iterrows():
        gid = int(row["global_id"])
        if gid > 0 and gid <= len(s_at_elements):
            s_sp = s_at_elements[gid - 1]
            ax.axvline(s_sp, color="gray", lw=0.3, alpha=0.5)
            ax.text(s_sp, ax.get_ylim()[1]*0.9, f"SP{int(row['sp_number'])}",
                    fontsize=7, color="gray", rotation=90, va="top")

    ax.set_xlabel("Path length s (m)")
    ax.set_ylabel("|P1 - design orbit| (mm)")
    ax.set_title(f"Experiment {exp_id}: Where does P1 deviate from design orbit?\n"
                 f"Top {len(peaks)} local maxima labeled (gray lines = SP boundaries)")
    ax.legend()
    ax.grid(True, alpha=0.3)

    out_path = WS / "results" / f"orbit_deviation_peaks_exp{exp_id}.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=120)
    print(f"\nPlot saved: {out_path}")
    plt.show()


if __name__ == "__main__":
    exp_id = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    analyze(exp_id)
