"""Построение ПЛОТНОГО опорного пути для path-following loss.

Идея (предложение пользователя): loss считает отклонение не от дискретных центров
магнитов, а от ОПОРНОЙ ТРАЕКТОРИИ, заданной точками:
  * центры железа каждого элемента (везде);
  * в ПРЯМЫХ участках (дрейфы между соосными элементами) — дополнительные точки
    ВДОЛЬ прямой (это задаёт НАПРАВЛЕНИЕ, а не только позицию);
  * в КРИВЫХ участках (через диполи) — только центры (изгиб определяется полем).
  * 50-метровый дрейф после SP1 — тоже прямой (yaw=9°), уплотняем вдоль линии:
    плотная линия не даёт "выстрелить" одним магнитом (в отличие от 5 редких целей).

Точки тэгируются gid'ом (для оконного режима composite_next: точка входит в loss,
когда её gid <= gid_last окна). Дрейф-точки между элементами i и i+1 -> gid = i.
Точки 50м-дрейфа после последнего gid -> gid = последний.

Выход: список (gid, x, y, z) — печатается и/или возвращается для эмита PATHPT-строк.

Запуск:  python build_reference_path.py [lattice_id=105] [spacing_m=0.25] [drift_end_x drift_end_y]
"""
from __future__ import annotations
import sys, csv, math
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except Exception:
    pass


def load_endpoints(path: Path):
    rows = list(csv.reader(open(path, encoding='utf-8'), delimiter=';'))
    h = rows[0]
    idx = {name: h.index(name) for name in
           ('global_id', 'device_type_id', 'x_in', 'y_in', 'yaw_in',
            'x_out', 'y_out', 'yaw_out')}
    els = []
    for r in rows[1:]:
        els.append(dict(
            gid=int(r[idx['global_id']]),
            dt=r[idx['device_type_id']].strip(),
            xin=float(r[idx['x_in']]), yin=float(r[idx['y_in']]),
            yawin=float(r[idx['yaw_in']]),
            xout=float(r[idx['x_out']]), yout=float(r[idx['y_out']]),
            yawout=float(r[idx['yaw_out']]),
        ))
    els.sort(key=lambda e: e['gid'])
    return els


def ang_diff_deg(a, b):
    """Минимальная разница углов в градусах (-180..180]."""
    d = (a - b) * 180.0 / math.pi
    while d > 180.0:
        d -= 360.0
    while d <= -180.0:
        d += 360.0
    return abs(d)


def build_path(els, spacing=0.25, straight_tol_deg=1.0,
               drift_end=None, drift_gid=None,
               drift_chunk_m=0.0, drift_point_spacing=1.0):
    """Возвращает список точек пути: (gid, x, y, z).

    spacing       -- шаг уплотнения прямых участков SP1 (м).
    straight_tol_deg -- допуск 'прямизны' (град).
    drift_end     -- (x,y) конца длинного 50м-дрейфа (если задан, уплотняем его).
    drift_gid     -- базовый gid для точек дрейфа (обычно последний gid SP1).
    drift_chunk_m -- если >0: прогрессивные виртуальные gid'ы по дрейфу, по одному
                     на каждые drift_chunk_m метров (для "вползания" окна в дрейф).
                     0 = все точки дрейфа в один drift_gid (старое поведение).
    drift_point_spacing -- шаг точек вдоль дрейфа (м).
    """
    pts = []  # (gid, x, y, z)

    # 1. Центры элементов (midpoint входа/выхода). Со-локализованные схлопываем.
    centers = []  # (gid, cx, cy, yaw)
    for e in els:
        cx = 0.5 * (e['xin'] + e['xout'])
        cy = 0.5 * (e['yin'] + e['yout'])
        # средний угол (для длинных гнутых элементов = (in+out)/2)
        yaw = 0.5 * (e['yawin'] + e['yawout'])
        centers.append((e['gid'], cx, cy, yaw, e))

    for gid, cx, cy, yaw, e in centers:
        pts.append((gid, cx, cy, 0.0))

    # 2. Уплотнение прямых дрейфов между ПОСЛЕДОВАТЕЛЬНЫМИ РАЗНЕСЁННЫМИ элементами.
    for i in range(len(els) - 1):
        a, b = els[i], els[i + 1]
        ax, ay = a['xout'], a['yout']
        bx, by = b['xin'], b['yin']
        dx, dy = bx - ax, by - ay
        seglen = math.hypot(dx, dy)
        if seglen < 1e-4:
            continue  # со-локализованные/нулевой зазор -> пропускаем
        seg_yaw = math.atan2(dy, dx)
        # 'Прямой' = выход a, направление сегмента и вход b сонаправлены.
        d1 = ang_diff_deg(a['yawout'], seg_yaw)
        d2 = ang_diff_deg(b['yawin'], seg_yaw)
        if d1 > straight_tol_deg or d2 > straight_tol_deg:
            continue  # кривой участок -> только центры (уже добавлены)
        if seglen < 0.03:
            continue  # совсем короткий зазор -- направление там неважно
        # уплотняем: >=3 деления (>=2 промежуточных точки) на каждый прямой сегмент,
        # либо чаще при длинном дрейфе (шаг spacing).
        n = max(3, int(seglen / spacing))
        for k in range(1, n):
            t = k / n
            pts.append((a['gid'], ax + t * dx, ay + t * dy, 0.0))

    # 3. 50-метровый дрейф после SP1: прямая линия по касательной выхода,
    #    уплотнённая точками с шагом drift_point_spacing. gid'ы:
    #    - drift_chunk_m <= 0: ВСЕ точки = drift_gid (старое поведение, всё в финале);
    #    - drift_chunk_m  > 0: ПРОГРЕССИВНЫЕ виртуальные gid'ы (drift_gid+1, +2, ...),
    #      по одному на каждые drift_chunk_m метров дрейфа. Это даёт оконному режиму
    #      "вползать" в дрейф по drift_chunk_m за шаг (окно за окном), а не получать
    #      весь дрейф разом в финальном окне (где раньше был "выстрел").
    if drift_end is not None and drift_gid is not None:
        last = els[-1]
        ax, ay = last['xout'], last['yout']
        ex, ey = drift_end
        dx, dy = ex - ax, ey - ay
        seglen = math.hypot(dx, dy)
        if seglen > 1e-6:
            ux, uy = dx / seglen, dy / seglen
            sp = max(0.25, drift_point_spacing)
            n = int(seglen / sp)
            for k in range(1, n + 1):
                d = k * sp
                if d > seglen + 1e-9:
                    break
                if drift_chunk_m > 0:
                    g = drift_gid + 1 + int(d // drift_chunk_m)  # прогрессивный gid
                else:
                    g = drift_gid                                # всё в один gid (старое)
                pts.append((g, ax + ux * d, ay + uy * d, 0.0))

    return pts


def main():
    root = Path(__file__).resolve().parent
    lat = int(sys.argv[1]) if len(sys.argv) > 1 else 105
    spacing = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
    ep_path = root / "workspace" / "results" / f"magnet_endpoints_lat{lat}.csv"
    els = load_endpoints(ep_path)
    # 50м-дрейф: конец = последняя extra-цель (центр финиш-плоскости) 77.006, 9.9956
    drift_end = (77.00608, 9.995595)
    drift_gid = els[-1]['gid']
    pts = build_path(els, spacing=spacing, drift_end=drift_end, drift_gid=drift_gid)
    # сводка
    from collections import Counter
    by_gid = Counter(p[0] for p in pts)
    n_centers = len(els)
    n_drift_intra = sum(1 for p in pts) - n_centers
    print(f"Элементов (центров): {n_centers}")
    print(f"Всего точек пути: {len(pts)}  (центры {n_centers} + уплотнение {len(pts)-n_centers})")
    print(f"Точек на 50м-дрейфе (gid={drift_gid}): {by_gid[drift_gid] - 1}")  # минус сам центр gid69
    print("Примеры точек в первых дрейфах:")
    shown = 0
    for gid, x, y, z in pts[n_centers:n_centers+12]:
        print(f"   gid={gid:2d}  ({x:.3f}, {y:.3f})")
        shown += 1


if __name__ == "__main__":
    main()
