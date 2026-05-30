"""
Сканирование бинарных карт полей (workspace/fields_bin/*.bin) на предмет провалов.

Провал (dropout) -- узел сетки, где |B| ≈ 0, хотя СОСЕДНИЕ узлы имеют значимое
поле. Бывает из-за дефектов экспорта Ansys или повреждений файла.

Алгоритм детектора:
  1. Для каждого узла (i,j,k) считаем |B| и |B|_neighbors (6 face-neighbors).
  2. Если у соседей среднее |B| значимое (> neighbor_eps), а |B|_node
     значительно меньше (ratio < dropout_ratio), -- помечаем как провал.

Использование:
    python scan_field_maps.py                    # сканировать все карты, только отчёт
    python scan_field_maps.py DQ_I1_1A           # сканировать конкретную карту
    python scan_field_maps.py --fix              # автоматически починить провалы
                                                 # (бэкап -> *.bak, перезапись)
    python scan_field_maps.py --threshold 0.05   # ratio для отнесения к провалу
                                                 # (по умолчанию 0.1, т.е. 10%)

Замена при `--fix`:
  Дроп-узел заменяется средним значением валидных соседей (6-connectivity,
  отбрасываются другие провалы и пустые края).
"""
import struct
import sys
import shutil
import os
from pathlib import Path

import numpy as np


SCRIPT_DIR     = Path(__file__).resolve().parent
FIELDS_BIN_DIR = SCRIPT_DIR / "workspace" / "fields_bin"

# Параметры детектора (можно переопределить через CLI)
DEFAULT_DROPOUT_RATIO = 0.1     # узел считается провалом если |B| < ratio × mean(|B|_neighbors)
DEFAULT_NEIGHBOR_EPS  = 1e-9    # минимальное среднее |B|_neighbors, чтобы вообще проверять


# C++ GridHeader struct layout (с учётом 4-байтного padding перед double nominalArg):
#   9 doubles (min/max/step xyz) + 3 int32 (nx, ny, nz) + 4 bytes padding + 1 double = 96 bytes
HEADER_FMT  = "<9d3i4xd"
HEADER_SIZE = struct.calcsize(HEADER_FMT)
assert HEADER_SIZE == 96, f"Unexpected header size: {HEADER_SIZE}"


def read_header(path):
    """Прочитать заголовок: возвращает dict с границами, шагами, размерами, nominal_arg."""
    with open(path, "rb") as f:
        raw = f.read(HEADER_SIZE)
    vals = struct.unpack(HEADER_FMT, raw)
    return {
        "minX": vals[0], "minY": vals[1], "minZ": vals[2],
        "maxX": vals[3], "maxY": vals[4], "maxZ": vals[5],
        "stepX": vals[6], "stepY": vals[7], "stepZ": vals[8],
        "nx": vals[9], "ny": vals[10], "nz": vals[11],
        "nominalArg": vals[12],
    }


def read_field(path):
    """Загрузить заголовок + всё поле как np.array (nx, ny, nz, 3) в SI (Тесла)."""
    h = read_header(path)
    n_total = h["nx"] * h["ny"] * h["nz"]
    with open(path, "rb") as f:
        f.seek(HEADER_SIZE)
        data = np.fromfile(f, dtype=np.float64, count=n_total * 3)
    if data.size != n_total * 3:
        raise IOError(f"Truncated file {path}: expected {n_total * 3} doubles, got {data.size}")
    # ВНИМАНИЕ: важно понять порядок индексирования.
    # В коде AnsysConverter.h данные записываются в порядке: внешний цикл по x,
    # внутри по y, внутри по z. То есть индекс = (i_x * ny + i_y) * nz + i_z.
    # Соответствующая форма numpy: (nx, ny, nz, 3) с reshape по этому порядку.
    field = data.reshape(h["nx"], h["ny"], h["nz"], 3)
    return h, field


def write_field(path, header, field):
    """Перезапись бинарника (заголовок неизменён, данные заменены)."""
    n_total = header["nx"] * header["ny"] * header["nz"]
    assert field.shape == (header["nx"], header["ny"], header["nz"], 3)
    # Восстановим бинарный заголовок (читаем оригинал и берём первые 96 байт)
    with open(path, "rb") as f:
        raw_header = f.read(HEADER_SIZE)
    with open(path, "wb") as f:
        f.write(raw_header)
        f.write(field.astype(np.float64, copy=False).tobytes())


def find_dropouts(field, dropout_ratio=DEFAULT_DROPOUT_RATIO,
                  neighbor_eps=DEFAULT_NEIGHBOR_EPS):
    """Вернуть mask (nx, ny, nz) с True на узлах-провалах.

    Провал: |B_node| < ratio * mean(|B|_neighbors), при том что
    mean(|B|_neighbors) > eps (исключаем true-zero regions).
    """
    Bmag = np.linalg.norm(field, axis=-1)  # (nx, ny, nz)
    nx, ny, nz = Bmag.shape

    # Среднее по 6 face-neighbors. Используем np.roll и маску валидности.
    # Для краёв (i=0/nx-1 и т.д.) соответствующий сосед "выпадает" -- считаем
    # только реальных соседей.
    neighbor_sum   = np.zeros_like(Bmag)
    neighbor_count = np.zeros_like(Bmag, dtype=np.int32)

    for axis in range(3):
        for shift in (-1, +1):
            shifted = np.roll(Bmag, shift, axis=axis)
            # Маска валидности: узел не должен оборачиваться за край
            valid = np.ones_like(Bmag, dtype=bool)
            if shift == -1:
                idx = [slice(None)] * 3; idx[axis] = -1
                valid[tuple(idx)] = False
            else:
                idx = [slice(None)] * 3; idx[axis] = 0
                valid[tuple(idx)] = False
            neighbor_sum   += np.where(valid, shifted, 0.0)
            neighbor_count += valid.astype(np.int32)

    neighbor_mean = np.where(neighbor_count > 0, neighbor_sum / np.maximum(neighbor_count, 1), 0.0)

    # Маска провалов
    significant_around = neighbor_mean > neighbor_eps
    is_dropout = significant_around & (Bmag < dropout_ratio * neighbor_mean)
    return is_dropout, Bmag, neighbor_mean


def fix_dropouts(field, dropout_mask):
    """Заменить значения в узлах-провалах средним по 6 face-neighbors.
    Соседи, которые САМИ помечены как провал, при усреднении игнорируются.
    Возвращает (field_fixed, count_fixed)."""
    fixed = field.copy()
    nx, ny, nz = field.shape[:3]
    drop_idx = np.argwhere(dropout_mask)
    fixed_count = 0
    valid_mask = ~dropout_mask
    for (i, j, k) in drop_idx:
        neigh_vals = []
        for axis, di, dj, dk in [(0, -1, 0, 0), (0, +1, 0, 0),
                                 (1, 0, -1, 0), (1, 0, +1, 0),
                                 (2, 0, 0, -1), (2, 0, 0, +1)]:
            ni, nj, nk = i + di, j + dj, k + dk
            if 0 <= ni < nx and 0 <= nj < ny and 0 <= nk < nz:
                if valid_mask[ni, nj, nk]:
                    neigh_vals.append(field[ni, nj, nk])
        if neigh_vals:
            fixed[i, j, k] = np.mean(neigh_vals, axis=0)
            fixed_count += 1
        # Если все соседи -- тоже провалы, оставляем как есть (редкий случай)
    return fixed, fixed_count


def scan_one(path, dropout_ratio, neighbor_eps, fix=False):
    """Возвращает словарь со статистикой по одной карте."""
    print(f"\n>>> Scanning {path.name} ...")
    h, field = read_field(path)
    n_total = h["nx"] * h["ny"] * h["nz"]
    Bmag = np.linalg.norm(field, axis=-1)
    nonzero_nodes = int(np.sum(Bmag > neighbor_eps))
    nonzero_pct   = 100.0 * nonzero_nodes / n_total
    print(f"  grid: {h['nx']}x{h['ny']}x{h['nz']} = {n_total:,} nodes  "
          f"(bbox: [{h['minX']:+.3f}..{h['maxX']:+.3f}, "
          f"{h['minY']:+.3f}..{h['maxY']:+.3f}, "
          f"{h['minZ']:+.3f}..{h['maxZ']:+.3f}] m)")
    print(f"  nominal_arg: {h['nominalArg']}  |B| range: [{Bmag.min():.3e} .. {Bmag.max():.3e}] T")
    print(f"  nonzero nodes: {nonzero_nodes:,} ({nonzero_pct:.2f}%)")

    stats = {
        "name": path.name,
        "n_total": n_total,
        "n_nonzero": nonzero_nodes,
        "n_drops": 0,
        "drops_pct": 0.0,
        "fixed": 0,
        "min_ratio": None,
        "n_edge": 0,
        "n_interior": 0,
        "mean_ratio": None,
        "Bmax": float(Bmag.max()),
    }

    if nonzero_nodes == 0:
        print(f"  (карта полностью нулевая - skip)")
        return stats

    drop_mask, _, neigh_mean = find_dropouts(field, dropout_ratio, neighbor_eps)
    drop_count = int(np.sum(drop_mask))
    drops_pct  = 100.0 * drop_count / nonzero_nodes
    print(f"  >>> DROPOUTS found: {drop_count:,} ({drops_pct:.3f}% of non-zero nodes)")

    stats["n_drops"]  = drop_count
    stats["drops_pct"] = drops_pct

    if drop_count > 0:
        nx, ny, nz = h["nx"], h["ny"], h["nz"]
        idx = np.argwhere(drop_mask)
        # Распределение: edge (любой индекс на границе) vs interior
        is_edge = (idx[:, 0] == 0) | (idx[:, 0] == nx - 1) | \
                  (idx[:, 1] == 0) | (idx[:, 1] == ny - 1) | \
                  (idx[:, 2] == 0) | (idx[:, 2] == nz - 1)
        n_edge = int(np.sum(is_edge))
        n_interior = drop_count - n_edge
        stats["n_edge"] = n_edge
        stats["n_interior"] = n_interior

        # ratios = |B| / neighbor_mean у узлов-провалов
        ratios = []
        for (i, j, k) in idx:
            nm = neigh_mean[i, j, k]
            ratios.append(Bmag[i, j, k] / nm if nm > 0 else 0.0)
        ratios = np.asarray(ratios)
        stats["min_ratio"] = float(ratios.min())
        stats["mean_ratio"] = float(ratios.mean())

        print(f"     edge:     {n_edge:,} ({100*n_edge/drop_count:.1f}%)")
        print(f"     interior: {n_interior:,} ({100*n_interior/drop_count:.1f}%)")
        print(f"     ratio (|B|/neigh): min={ratios.min():.4f}  mean={ratios.mean():.4f}  max={ratios.max():.4f}")

        # Первые 5 примеров для контекста
        print(f"     examples (up to 5):")
        for (i, j, k) in idx[:5]:
            x = h["minX"] + i * h["stepX"]
            y = h["minY"] + j * h["stepY"]
            z = h["minZ"] + k * h["stepZ"]
            b = Bmag[i, j, k]; nm = neigh_mean[i, j, k]
            print(f"       node ({i:3d},{j:3d},{k:3d}) at ({x:+.4f},{y:+.4f},{z:+.4f}) m  "
                  f"|B|={b:.3e}  neigh_mean={nm:.3e}  ratio={b/(nm+1e-30):.3f}")

        if fix:
            print(f"  --- FIXING {drop_count:,} dropouts ---")
            field_fixed, fixed_count = fix_dropouts(field, drop_mask)
            backup = path.with_suffix(path.suffix + ".bak")
            if not backup.exists():
                shutil.copy(path, backup)
                print(f"     backup -> {backup.name}")
            write_field(path, h, field_fixed)
            stats["fixed"] = fixed_count
            print(f"     patched: {fixed_count:,} nodes -> {path.name}")
            if fixed_count < drop_count:
                print(f"     warning: {drop_count - fixed_count} dropouts had no valid neighbors (left as-is)")

    return stats


def print_summary(stats_list, fix):
    """Финальная таблица со статистикой по всем картам."""
    print(f"\n{'='*80}")
    print(f"FINAL SUMMARY")
    print(f"{'='*80}")

    name_w = max(len(s["name"]) for s in stats_list) if stats_list else 20
    name_w = max(name_w, 10)

    header = f"{'name':<{name_w}}  {'nodes':>10}  {'drops':>8}  {'%nz':>7}  {'edge':>7}  {'interior':>8}  {'min ratio':>10}"
    if fix:
        header += f"  {'fixed':>7}"
    print(header)
    print("-" * len(header))

    total_drops = 0
    total_edge = 0
    total_interior = 0
    total_fixed = 0
    total_nodes = 0
    files_clean = 0
    files_dirty = 0

    for s in stats_list:
        total_drops += s["n_drops"]
        total_edge  += s["n_edge"]
        total_interior += s["n_interior"]
        total_fixed += s["fixed"]
        total_nodes += s["n_total"]
        if s["n_drops"] == 0:
            files_clean += 1
        else:
            files_dirty += 1
        min_r = f"{s['min_ratio']:.4f}" if s["min_ratio"] is not None else "  -  "
        row = (f"{s['name']:<{name_w}}  {s['n_total']:>10,}  {s['n_drops']:>8,}  "
               f"{s['drops_pct']:>6.3f}%  {s['n_edge']:>7,}  {s['n_interior']:>8,}  {min_r:>10}")
        if fix:
            row += f"  {s['fixed']:>7,}"
        print(row)

    print("-" * len(header))
    print(f"  files scanned:  {len(stats_list)}  "
          f"({files_clean} clean, {files_dirty} with dropouts)")
    print(f"  total nodes:    {total_nodes:,}")
    print(f"  total dropouts: {total_drops:,}  "
          f"(edge={total_edge:,}, interior={total_interior:,})")
    if total_nodes > 0:
        print(f"  overall rate:   {100*total_drops/total_nodes:.5f}% of all nodes")
    if fix:
        print(f"  total fixed:    {total_fixed:,}  ({total_drops - total_fixed} unfixed)")


def main():
    args = sys.argv[1:]
    fix = "--fix" in args
    args = [a for a in args if a != "--fix"]

    threshold = DEFAULT_DROPOUT_RATIO
    if "--threshold" in args:
        i = args.index("--threshold")
        threshold = float(args[i + 1])
        del args[i:i + 2]

    if not FIELDS_BIN_DIR.exists():
        print(f"ERROR: {FIELDS_BIN_DIR} not found")
        sys.exit(1)

    if args:
        targets = []
        for a in args:
            p = FIELDS_BIN_DIR / (a if a.endswith(".bin") else f"{a}.bin")
            if not p.exists():
                print(f"WARN: {p} not found, skipped")
                continue
            targets.append(p)
    else:
        targets = sorted(FIELDS_BIN_DIR.glob("*.bin"))

    if not targets:
        print("No .bin field maps found")
        return

    print(f"Scanning {len(targets)} field map(s)")
    print(f"  dropout_ratio: {threshold}  (node flagged if |B|/neigh_mean < this)")
    print(f"  neighbor_eps:  {DEFAULT_NEIGHBOR_EPS:.0e}  (skip true-zero regions)")
    print(f"  mode: {'FIX (overwrite with patched values, backups -> *.bak)' if fix else 'scan-only (read-only)'}")

    all_stats = []
    for p in targets:
        all_stats.append(scan_one(p, threshold, DEFAULT_NEIGHBOR_EPS, fix=fix))

    print_summary(all_stats, fix)

    total_drops = sum(s["n_drops"] for s in all_stats)
    if total_drops > 0 and not fix:
        print(f"\n>>> Re-run with --fix to patch them (backups created automatically). <<<")


if __name__ == "__main__":
    main()
