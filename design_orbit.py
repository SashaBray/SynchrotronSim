"""Сборка опорной (design) орбиты для конкретного lattice_id.

Состав:
    Segment(kind, length, angle) -- один геометрический примитив в SP.
    load_segments_db(ws)         -> dict[sp_type] -> list[Segment]
    load_sublattices(ws)         -> dict[lattice_id] -> max_gid
    sp_sequence_for_lattice(...) -> list[(sp_number, sp_type)] для lattice
    build_orbit(seq, segs_db, points_per_segment=10)
                                 -> np.ndarray (N, 3) точек орбиты
    build_orbit_for_lattice(lattice_id, workspace=None, points_per_segment=10)
                                 -> np.ndarray (N, 3) -- удобный шорткат.

Стартовая позиция первого SP в lattice = (0, 0, 0), yaw = 0
(те же конвенции, что использует lte_to_lattice.py при построении
lattice_layout.csv).

Используется:
    - visualize.py -- отрисовка design orbit вместо зигзага через центры железа;
    - prepare_optim_config.py -- генерация design_orbit_lat{id}.csv для C++;
    - C++ optimizer.cpp -- через предсобранный CSV.
"""
from __future__ import annotations
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np

from _common import project_root, workspace_dir


@dataclass(frozen=True)
class Segment:
    kind: str           # "line" | "arc"
    length: float       # м
    angle: float        # рад (0 для line)
    lte_name: str = ""  # имя LTE-элемента (для связки с lattice_layout через map_to_hw)


@dataclass(frozen=True)
class MagnetEndpoints:
    """Вход/выход одного магнита (global_id) на эталонной орбите."""
    global_id: int
    device_type_id: str
    lte_name: str
    entry: tuple[float, float, float]
    exit:  tuple[float, float, float]
    entry_yaw: float
    exit_yaw:  float


# ---------------------------------------------------------------------------
# Загрузка БД
# ---------------------------------------------------------------------------
def load_segments_db(ws: Path | None = None) -> dict[str, list[Segment]]:
    """Читает workspace/tables/design_orbit_segments.csv."""
    ws = ws or workspace_dir()
    path = ws / "tables" / "design_orbit_segments.csv"
    db: dict[str, list[tuple[int, Segment]]] = {}
    with open(path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";")
        for row in reader:
            sp_type = row["sp_type"].strip()
            idx     = int(row["seg_idx"])
            seg = Segment(kind=row["kind"].strip(),
                          length=float(row["length"]),
                          angle=float(row["angle"]),
                          lte_name=row.get("lte_name", "").strip())
            db.setdefault(sp_type, []).append((idx, seg))
    # Сортируем по seg_idx, оставляем только сегменты
    return {t: [s for _, s in sorted(lst, key=lambda p: p[0])]
            for t, lst in db.items()}


def load_sublattices(ws: Path | None = None) -> dict[int, int]:
    """Возвращает {lattice_id: max_gid}. Для главного кольца (101) max_gid=∞."""
    ws = ws or workspace_dir()
    path = ws / "configs" / "sub_lattices.json"
    out: dict[int, int] = {}
    if path.exists():
        cfg = json.loads(path.read_text(encoding="utf-8"))
        for entry in cfg.get("sub_lattices", []):
            out[int(entry["lattice_id"])] = int(entry["max_gid"])
    return out


# ---------------------------------------------------------------------------
# Какие SP входят в lattice
# ---------------------------------------------------------------------------
def sp_sequence_for_lattice(lattice_id: int,
                            ws: Path | None = None) -> list[tuple[int, str]]:
    """Возвращает [(sp_number, sp_type), ...] по lattice_layout.csv,
    отфильтровав по max_gid из sub_lattices.json (для sub-lattice).
    Главный lattice = 101 не имеет ограничения -- берём все SP.
    """
    ws = ws or workspace_dir()
    layout_path = ws / "tables" / "lattice_layout.csv"

    sub = load_sublattices(ws)
    max_gid = sub.get(lattice_id, 10**12)  # ∞ если нет в sub_lattices

    seen: dict[int, str] = {}    # sp_number -> sp_type (первый встреченный)
    with open(layout_path, encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=";", skipinitialspace=True)
        for row in reader:
            gid = int(row["global_id"])
            if gid > max_gid:
                break
            sp_num = int(row["sp_number"])
            if sp_num not in seen:
                seen[sp_num] = row["sp_type"].strip()
    return [(n, seen[n]) for n in sorted(seen)]


# ---------------------------------------------------------------------------
# Сборка точек орбиты по сегментам
# ---------------------------------------------------------------------------
_DYAW_EPS  = 1e-9
_CHORD_EPS = 1e-12


def _apply_line(pts: list[tuple[float, float, float]],
                x: float, y: float, yaw: float,
                L: float, n_sub: int) -> tuple[float, float, float]:
    """Прямой отрезок длиной L. Добавляет n_sub точек (без стартовой)."""
    if L <= 0.0:
        return x, y, yaw
    cx, cy = math.cos(yaw), math.sin(yaw)
    for k in range(1, n_sub + 1):
        t = k / n_sub
        pts.append((x + t * L * cx, y + t * L * cy, 0.0))
    return x + L * cx, y + L * cy, yaw


def _apply_arc(pts: list[tuple[float, float, float]],
               x: float, y: float, yaw: float,
               L: float, angle: float, n_sub: int) -> tuple[float, float, float]:
    """Дуга длиной L с поворотом angle. Добавляет n_sub точек (без стартовой).
    Аналитика: kappa=angle/L; на пути s касательная = yaw + kappa*s;
    смещения вычислены интегралом cos/sin от касательной."""
    if L <= 0.0:
        return x, y, yaw
    if abs(angle) < _DYAW_EPS:
        return _apply_line(pts, x, y, yaw, L, n_sub)
    kappa = angle / L
    sin_y0 = math.sin(yaw)
    cos_y0 = math.cos(yaw)
    for k in range(1, n_sub + 1):
        s = k * L / n_sub
        yaw_s = yaw + kappa * s
        x_off = (math.sin(yaw_s) - sin_y0) / kappa
        y_off = (cos_y0 - math.cos(yaw_s)) / kappa
        pts.append((x + x_off, y + y_off, 0.0))
    # Конечные значения
    yaw_end = yaw + angle
    x_end = x + (math.sin(yaw_end) - sin_y0) / kappa
    y_end = y + (cos_y0 - math.cos(yaw_end)) / kappa
    return x_end, y_end, yaw_end


def _advance_pose(x: float, y: float, yaw: float,
                  L: float, angle: float) -> tuple[float, float, float]:
    """Сдвиг (x, y, yaw) после прохождения одного сегмента. Без накопления точек.
    Возвращает (x_end, y_end, yaw_end). Используется в compute_magnet_endpoints,
    где промежуточные точки не нужны."""
    if L <= 0.0:
        return x, y, yaw
    if abs(angle) < _DYAW_EPS:
        return x + L * math.cos(yaw), y + L * math.sin(yaw), yaw
    kappa = angle / L
    sin_y0 = math.sin(yaw)
    cos_y0 = math.cos(yaw)
    yaw_end = yaw + angle
    x_end = x + (math.sin(yaw_end) - sin_y0) / kappa
    y_end = y + (cos_y0 - math.cos(yaw_end)) / kappa
    return x_end, y_end, yaw_end


def compute_magnet_endpoints(lattice_id: int,
                             workspace: Path | None = None
                             ) -> tuple[list[MagnetEndpoints], float]:
    """Возвращает (endpoints, final_yaw) для всех global_id в lattice.

    Для каждого LTE-сегмента в SP вызывается lte.map_to_hw(): если сегмент
    эмитит N строк lattice_layout (N=1..3, для DQ/KSEXT), все N получают
    одинаковые entry/exit на эталонной орбите — потому что физически это
    несколько обмоток ОДНОГО магнита, и их позиция/yaw в lattice_configs
    совпадают.

    Список упорядочен по global_id (соответствует lattice_layout.csv).
    Драйфты (EDRIFT) и BPM (MONI) не попадают в результат -- они не имеют
    строк в lattice_layout (map_to_hw возвращает []).
    """
    ws = workspace or workspace_dir()

    # Загружаем определения LTE-устройств и применяем close_scale (как в
    # lte_to_design_orbit.py), чтобы углы и L согласовались с segments_db.
    import lte_to_lattice as lte
    devices = lte.parse_desc(str(ws / "lte" / "SYLA_desc.lte"))
    # close_scale тут не нужен для самого endpoints (мы берём углы из
    # segments_db, где они УЖЕ масштабированы), но devices используем
    # для map_to_hw — которому close_scale неважен (он смотрит на L, K1, K2).

    seq = sp_sequence_for_lattice(lattice_id, ws)
    db  = load_segments_db(ws)

    out: list[MagnetEndpoints] = []
    x, y, yaw = 0.0, 0.0, 0.0
    gid = 0
    for sp_num, sp_type in seq:
        if sp_type not in db:
            raise KeyError(f"SP {sp_num}: sp_type='{sp_type}' отсутствует в БД сегментов")
        for seg in db[sp_type]:
            entry = (x, y, 0.0)
            entry_yaw = yaw
            x_end, y_end, yaw_end = _advance_pose(x, y, yaw, seg.length, seg.angle)
            exit_pt   = (x_end, y_end, 0.0)
            exit_yaw  = yaw_end

            # Эмитит ли этот сегмент строки в lattice_layout?
            dev = devices.get(seg.lte_name)
            if dev is not None:
                hw_list = lte.map_to_hw(dev)
                for hw_id, _default_arg in hw_list:
                    gid += 1
                    out.append(MagnetEndpoints(
                        global_id      = gid,
                        device_type_id = hw_id,
                        lte_name       = seg.lte_name,
                        entry          = entry,
                        exit           = exit_pt,
                        entry_yaw      = entry_yaw,
                        exit_yaw       = exit_yaw,
                    ))
            # Иначе (драйфт / BPM): просто продвигаем pose без эмита gid
            x, y, yaw = x_end, y_end, yaw_end

    return out, yaw


def build_orbit(sp_sequence: Iterable[tuple[int, str]],
                segments_db: dict[str, list[Segment]],
                points_per_segment: int = 10,
                start_xyz: tuple[float, float, float] = (0.0, 0.0, 0.0),
                start_yaw: float = 0.0) -> tuple[np.ndarray, float]:
    """Собирает плотную точечную полилинию опорной орбиты.

    Возвращает (points, final_yaw):
        points    -- np.ndarray (N, 3), N = 1 + Σ_segments points_per_segment.
        final_yaw -- касательная орбиты на конце (полезно для postlude).
    """
    n_sub = max(1, int(points_per_segment))
    x, y, _ = start_xyz
    yaw = start_yaw
    pts: list[tuple[float, float, float]] = [(x, y, 0.0)]
    for sp_num, sp_type in sp_sequence:
        if sp_type not in segments_db:
            raise KeyError(
                f"SP {sp_num}: sp_type='{sp_type}' отсутствует в БД сегментов "
                f"(design_orbit_segments.csv). Перегенерируйте её "
                f"командой `python lte_to_design_orbit.py`.")
        for seg in segments_db[sp_type]:
            if seg.kind == "line":
                x, y, yaw = _apply_line(pts, x, y, yaw, seg.length, n_sub)
            elif seg.kind == "arc":
                x, y, yaw = _apply_arc (pts, x, y, yaw, seg.length, seg.angle, n_sub)
            else:
                raise ValueError(f"Unknown segment kind: {seg.kind!r}")
    return np.asarray(pts, dtype=np.float64), yaw


def build_orbit_for_lattice(lattice_id: int,
                            workspace: Path | None = None,
                            points_per_segment: int = 10
                            ) -> tuple[np.ndarray, float]:
    """Шорткат: загрузить БД + последовательность SP для lattice_id, собрать орбиту."""
    ws = workspace or workspace_dir()
    seq = sp_sequence_for_lattice(lattice_id, ws)
    db  = load_segments_db(ws)
    return build_orbit(seq, db, points_per_segment=points_per_segment)


# ---------------------------------------------------------------------------
# Пред- и постлюдия (для distanceToPolyline / устойчивости envelope_loss)
# ---------------------------------------------------------------------------
def extend_with_pre_post(points: np.ndarray,
                         final_yaw: float,
                         start_yaw: float = 0.0,
                         pre_length: float = 10.0,
                         post_length: float = 20.0) -> np.ndarray:
    """Добавляет по 1 точке с каждого конца -- продолжая касательную."""
    if len(points) < 1:
        return points
    cs, ss = math.cos(start_yaw), math.sin(start_yaw)
    cf, sf = math.cos(final_yaw), math.sin(final_yaw)
    pre  = np.array([points[0, 0] - pre_length  * cs,
                     points[0, 1] - pre_length  * ss,
                     points[0, 2]], dtype=np.float64)
    post = np.array([points[-1, 0] + post_length * cf,
                     points[-1, 1] + post_length * sf,
                     points[-1, 2]], dtype=np.float64)
    return np.vstack([pre, points, post])


# ---------------------------------------------------------------------------
# CLI / отладка
# ---------------------------------------------------------------------------
def _cli() -> None:
    import sys
    if len(sys.argv) < 2:
        print("Usage: python design_orbit.py <lattice_id> [points_per_seg]   -- orbit")
        print("       python design_orbit.py <lattice_id> endpoints           -- magnet entries/exits")
        sys.exit(1)
    lid = int(sys.argv[1])

    if len(sys.argv) > 2 and sys.argv[2] == "endpoints":
        eps, fy = compute_magnet_endpoints(lid)
        print(f"lattice {lid}: {len(eps)} magnets (lattice_layout rows)")
        print(f"  final yaw  = {math.degrees(fy):.4f} deg")
        print(f"  first 6:")
        for e in eps[:6]:
            print(f"    gid={e.global_id:>3}  {e.device_type_id:<16} {e.lte_name:<12} "
                  f"entry=({e.entry[0]:>9.4f},{e.entry[1]:>9.4f})  "
                  f"exit=({e.exit[0]:>9.4f},{e.exit[1]:>9.4f})  "
                  f"d_yaw={math.degrees(e.exit_yaw-e.entry_yaw):>+7.4f} deg")
        print(f"  last 6:")
        for e in eps[-6:]:
            print(f"    gid={e.global_id:>3}  {e.device_type_id:<16} {e.lte_name:<12} "
                  f"entry=({e.entry[0]:>9.4f},{e.entry[1]:>9.4f})  "
                  f"exit=({e.exit[0]:>9.4f},{e.exit[1]:>9.4f})  "
                  f"d_yaw={math.degrees(e.exit_yaw-e.entry_yaw):>+7.4f} deg")
        return

    nsub = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    pts, final_yaw = build_orbit_for_lattice(lid, points_per_segment=nsub)
    print(f"lattice {lid}: {len(pts)} points")
    print(f"  start = {tuple(pts[0])}")
    print(f"  end   = {tuple(pts[-1])}  final_yaw = {math.degrees(final_yaw):.4f} deg")
    # Длина полилинии
    diffs = np.diff(pts, axis=0)
    total_len = np.linalg.norm(diffs, axis=1).sum()
    print(f"  polyline length = {total_len:.4f} m")


if __name__ == "__main__":
    _cli()
