#!/usr/bin/env python3
"""single_element_boris.py — эталонный прогон ОДНОГО элемента через наш Boris-трекер
на полевой карте: 1 м дрейф + элемент + 1 м дрейф, прямая геометрия (ось пучка = +X).

Пучок: центр + 3 окружности по 12 спутников на радиусах 0.5/1/2 мм в поперечной (Y-Z)
плоскости (37 частиц). Несколько радиусов дают сильный нелинейный сигнал (∝r², r³) —
нужно для характеризации секступолей/октуполей.
Возвращает по каждой частице начальные (Y0,Z0) и финальные (Y,Y',Z,Z') у финиш-плоскости
(1 м за выходом). Эти данные match_coefficients.py подгоняет матричной моделью BIP.

Изолированный прогон: бэкап реальных lattice_configs/particle_groups/experiments,
запись минимальных версий, трек, восстановление (рабочее состояние НЕ трогается).
"""
import os, sys, struct, shutil, subprocess, math
from pathlib import Path

WS = Path(__file__).resolve().parent / "workspace"
TABLES = WS / "tables"
FIELDS = WS / "fields_bin"
RESULTS = WS / "results"
EXE = Path(__file__).resolve().parent / "build" / "SynchrotronTracker.exe"

EXP_ID, LAT_ID, GRP_ID = 200, 200, 99
TOTAL = 5.0                      # м — ПОЛНАЯ продольная область (старт пучка → финиш-плоскость)
DRIFT = 1.0                     # м — позиция ЦЕНТРА магнита: короткий дрейф ДО (~0.6м),
#                                ДЛИННЫЙ дрейф ПОСЛЕ (~3.6м) — усиливает погрешности → точнее фит
RADII = (0.0005, 0.001, 0.002)   # м — три окружности спутников (0.5/1/2 мм)
N_PER = 12                       # частиц на окружность
N_SAT = len(RADII) * N_PER       # всего спутников (36); +1 центр = 37 частиц
V_BEAM = 2.967e8
ME, QE = 9.109e-31, -1.602e-19
HDR = 96                         # размер GridHeader, байт

_HW = {}                         # device_type -> (active, field_dir, hw_corr, bin, nominal)


def _load_hw():
    if _HW: return
    for ln in (TABLES / "hardware_library.csv").read_text(encoding="utf-8").splitlines()[1:]:
        c = ln.split(";")
        if len(c) < 10 or not c[0]: continue
        _HW[c[0]] = dict(active=(c[1].lower() in ("true", "1")), binf=c[4],
                         nominal=(float(c[7]) if c[7] not in ("-", "") else 1.0),
                         field_dir=float(c[8]), hw_corr=float(c[9]))


def map_length(device):
    """Продольная (X) протяжённость полевой карты, м (maxX)."""
    _load_hw()
    p = FIELDS / _HW[device]["binf"]
    with open(p, "rb") as f:
        blob = f.read(HDR)
    minX, minY, minZ, maxX = struct.unpack_from("<4d", blob, 0)
    return maxX


def axis_field_integral(device, value, is_coeff):
    """∫By·dl, ∫Bz·dl (Тл·м) ПРИЛОЖЕННОГО поля вдоль оси (y=z=0) по ПОЛНОМУ симметричному
    пробегу [−maxX,maxX] = 2×[0,maxX]: поле чётно по x (зеркалится только Bx, см. FieldMap.h).
    Масштаб как у трекера: hw_corr·instance_coeff·field_dir·(ток/nominal если active)."""
    _load_hw(); hw = _HW[device]
    with open(FIELDS / hw["binf"], "rb") as f:
        blob = f.read(HDR)
        mnx, mny, mnz, mxx, mxy, mxz, sx, sy, sz = struct.unpack_from("<9d", blob, 0)
        nx, ny, nz = struct.unpack_from("<3i", blob, 72)
        j0, k0 = round((0 - mny) / sy), round((0 - mnz) / sz)   # узлы y=0, z=0
        By, Bz = [], []
        for i in range(nx):
            f.seek(HDR + (i * (ny * nz) + j0 * nz + k0) * 24)
            _, by, bz = struct.unpack("<3d", f.read(24))
            By.append(by); Bz.append(bz)
    Ihy = sum((By[i] + By[i + 1]) / 2 * sx for i in range(nx - 1))   # ∫[0,maxX]
    Ihz = sum((Bz[i] + Bz[i + 1]) / 2 * sx for i in range(nx - 1))
    if hw["active"] and not is_coeff:
        scale = hw["hw_corr"] * hw["field_dir"] * (value / hw["nominal"])
    else:
        scale = hw["hw_corr"] * value * hw["field_dir"]
    return 2 * Ihy * scale, 2 * Ihz * scale                        # ×2 — симметрия


def _satellites():
    """(Y,Z) всех 36 спутников: по 12 на окружность для каждого радиуса из RADII."""
    sats = []
    for r in RADII:
        for k in range(N_PER):
            th = 2 * math.pi * k / N_PER
            sats.append((r * math.cos(th), r * math.sin(th)))
    return sats


def _backup_restore(paths, restore=False):
    for p in paths:
        bak = p.with_suffix(p.suffix + ".sebak")
        if restore:
            if bak.exists(): shutil.move(str(bak), str(p))
        else:
            shutil.copy2(str(p), str(bak))


def _write_inputs(device, value, is_coeff, dt=6.74e-13):
    """Минимальные lattice_configs / particle_groups / experiments для одного элемента."""
    L = map_length(device)
    x_en = DRIFT                          # ЦЕНТР элемента (поле карты симметрично вокруг локального X=0)
    x_fin = TOTAL                         # финиш-плоскость на 2 м (фикс. полная продольная область)
    # --- lattice_configs: одна строка магнита, yaw=0 (продольная ось карты = +X) ---
    if is_coeff:                          # постоянный диполь: arg_val="-", свип correct_coeff
        arg, corr = "-", f"{value:.6f}"
    else:
        arg, corr = f"{value:.6f}", "1.0"
    lat = ["lattice_id;instance_id;device_type_id;x;y;z;yaw;pitch;roll;arg_val;correct_coeff",
           f"{LAT_ID};1;{device};{x_en:.6f};0;0;0;0;0;{arg};{corr}"]
    (TABLES / "lattice_configs.csv").write_text("\n".join(lat) + "\n", encoding="utf-8")
    # --- particle_groups: центр + 36 спутников (3 окружности 0.5/1/2 мм) в Y-Z ---
    pg = ["group_id;part_id;m;q;x;y;z;vx;vy;vz",
          f"{GRP_ID};1;{ME};{QE};0;0;0;{V_BEAM};0;0"]
    for idx, (y, z) in enumerate(_satellites()):
        pg.append(f"{GRP_ID};{idx + 2};{ME};{QE};0;{y:.8e};{z:.8e};{V_BEAM};0;0")
    (TABLES / "particle_groups.csv").write_text("\n".join(pg) + "\n", encoding="utf-8")
    # --- experiments: финиш-плоскость x=x_fin (нормаль +X), forward-пересечение ---
    stride = max(1, round(TOTAL / (V_BEAM * dt) / 1000))   # ~1000 точек траектории (длинный путь)
    exp = ["exp_id;lattice_id;group_id;dt;max_turns;traj_stride;p1x;p1y;p1z;p2x;p2y;p2z;p3x;p3y;p3z;dirx;diry;dirz",
           f"{EXP_ID};{LAT_ID};{GRP_ID};{dt:.8e};0;{stride};"
           f"{x_fin:.6f};-1;-1;{x_fin:.6f};1;-1;{x_fin:.6f};0;1;1;0;0"]
    (TABLES / "experiments.csv").write_text("\n".join(exp) + "\n", encoding="utf-8")
    return L, x_fin


def _read_trajectories():
    """trajectories_exp200.bin → список массивов точек [(x,y,z),...] на частицу."""
    p = RESULTS / f"trajectories_exp{EXP_ID}.bin"
    parts = []
    with open(p, "rb") as f:
        while True:
            h = f.read(4)
            if not h: break
            n = struct.unpack("I", h)[0]
            raw = f.read(n * 3 * 8)
            pts = struct.unpack(f"<{n * 3}d", raw)
            parts.append([(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]) for i in range(n)])
    return parts


def _final_phase(traj):
    """Финальные (Y, Y', Z, Z') у последней точки (наклоны по X из последних 2 точек)."""
    (x1, y1, z1), (x2, y2, z2) = traj[-2], traj[-1]
    dx = x2 - x1
    yp = (y2 - y1) / dx if abs(dx) > 1e-15 else 0.0
    zp = (z2 - z1) / dx if abs(dx) > 1e-15 else 0.0
    return (y2, yp, z2, zp)


def run(device, value, is_coeff=False, dt=6.74e-13):
    """Прогнать элемент при токе value (или correct_coeff, если is_coeff). Возвращает dict:
       L, x_fin, init=[(Y0,Z0)*13], final=[(Y,Yp,Z,Zp)*13] (частица 0 = центр).
       dt — шаг интегрирования Boris, с (по умолчанию рабочий)."""
    paths = [TABLES / "lattice_configs.csv", TABLES / "particle_groups.csv", TABLES / "experiments.csv"]
    _backup_restore(paths, restore=False)
    try:
        L, x_fin = _write_inputs(device, value, is_coeff, dt=dt)
        r = subprocess.run([str(EXE), str(EXP_ID)], cwd=str(EXE.parent.parent),
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            raise RuntimeError(f"tracker rc={r.returncode}\n{r.stdout[-500:]}\n{r.stderr[-500:]}")
        trajs = _read_trajectories()
    finally:
        _backup_restore(paths, restore=True)
    if len(trajs) != 1 + N_SAT:
        raise RuntimeError(f"ожидалось {1 + N_SAT} частиц, получено {len(trajs)}")
    init = [(0.0, 0.0)] + _satellites()
    final = [_final_phase(t) for t in trajs]
    return dict(device=device, value=value, is_coeff=is_coeff, L=L, x_fin=x_fin,
                init=init, final=final, trajs=trajs, drift=DRIFT)


def run_multi(devices, currents, dt=6.74e-13):
    """Прогон НЕСКОЛЬКИХ устройств в одной позиции (поля складываются) через Boris.
       devices — список device_type; currents — токи (А, arg_val). Возврат как run().
       Нужно для диполь-квадруполя: 3 обмотки (I1/I2/I3) сразу → корректная матрица связи."""
    paths = [TABLES / "lattice_configs.csv", TABLES / "particle_groups.csv", TABLES / "experiments.csv"]
    _backup_restore(paths, restore=False)
    try:
        L = map_length(devices[0]); x_en = DRIFT; x_fin = TOTAL
        lat = ["lattice_id;instance_id;device_type_id;x;y;z;yaw;pitch;roll;arg_val;correct_coeff"]
        for n, (dev, cur) in enumerate(zip(devices, currents), start=1):
            lat.append(f"{LAT_ID};{n};{dev};{x_en:.6f};0;0;0;0;0;{cur:.6f};1.0")
        (TABLES / "lattice_configs.csv").write_text("\n".join(lat) + "\n", encoding="utf-8")
        pg = ["group_id;part_id;m;q;x;y;z;vx;vy;vz",
              f"{GRP_ID};1;{ME};{QE};0;0;0;{V_BEAM};0;0"]
        for idx, (y, z) in enumerate(_satellites()):
            pg.append(f"{GRP_ID};{idx + 2};{ME};{QE};0;{y:.8e};{z:.8e};{V_BEAM};0;0")
        (TABLES / "particle_groups.csv").write_text("\n".join(pg) + "\n", encoding="utf-8")
        stride = max(1, round(TOTAL / (V_BEAM * dt) / 1000))
        exp = ["exp_id;lattice_id;group_id;dt;max_turns;traj_stride;p1x;p1y;p1z;p2x;p2y;p2z;p3x;p3y;p3z;dirx;diry;dirz",
               f"{EXP_ID};{LAT_ID};{GRP_ID};{dt:.8e};0;{stride};"
               f"{x_fin:.6f};-1;-1;{x_fin:.6f};1;-1;{x_fin:.6f};0;1;1;0;0"]
        (TABLES / "experiments.csv").write_text("\n".join(exp) + "\n", encoding="utf-8")
        r = subprocess.run([str(EXE), str(EXP_ID)], cwd=str(EXE.parent.parent),
                           capture_output=True, text=True, encoding="utf-8", errors="replace")
        if r.returncode != 0:
            raise RuntimeError(f"tracker rc={r.returncode}\n{r.stdout[-500:]}\n{r.stderr[-500:]}")
        trajs = _read_trajectories()
    finally:
        _backup_restore(paths, restore=True)
    if len(trajs) != 1 + N_SAT:
        raise RuntimeError(f"ожидалось {1 + N_SAT} частиц, получено {len(trajs)}")
    init = [(0.0, 0.0)] + _satellites()
    final = [_final_phase(t) for t in trajs]
    return dict(devices=devices, currents=currents, L=L, x_fin=x_fin,
                init=init, final=final, trajs=trajs, drift=DRIFT)


if __name__ == "__main__":
    dev = sys.argv[1] if len(sys.argv) > 1 else "QP_L_1A"
    val = float(sys.argv[2]) if len(sys.argv) > 2 else 50.0
    res = run(dev, val)
    print(f"{dev} @ {val}:  L_map={res['L']:.3f} м  финиш x={res['x_fin']:.3f} м")
    cy, cyp, cz, czp = res["final"][0]
    print(f"  центр финал: Y={cy*1e3:+.4f}мм Y'={cyp*1e3:+.4f}мрад  Z={cz*1e3:+.4f}мм Z'={czp*1e3:+.4f}мрад")
    print(f"  {'#':>2} {'Y0мм':>7} {'Z0мм':>7} | {'Yмм':>8} {'Yмрад':>8} {'Zмм':>8} {'Zмрад':>8}  (отн. центра)")
    for i in range(1, 1 + N_SAT):
        y0, z0 = res["init"][i]
        y, yp, z, zp = res["final"][i]
        print(f"  {i:>2} {y0*1e3:>7.3f} {z0*1e3:>7.3f} | {(y-cy)*1e3:>8.4f} {(yp-cyp)*1e3:>8.4f} "
              f"{(z-cz)*1e3:>8.4f} {(zp-czp)*1e3:>8.4f}")
