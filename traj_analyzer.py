# -*- coding: utf-8 -*-
"""АНАЛИЗАТОР ТРАЕКТОРИЙ SynchrotronSim в ЛОКАЛЬНЫХ координатах матричного ядра.

Читает results/trajectories_local_exp{N}.npz (пишется `track3.py --mode kick-phys --exp N --local`):
на точку — S (опорная арк-длина, м), x, x'≡px, y, y'≡py, ct, δ и индекс элемента. Работает В ЛОКАЛЬНЫХ
координатах (отклонение от опорной орбиты + канонические углы), а НЕ в глобальных XYZ.

По каждой частице: дошла/потеряна (и где), огибающая (max|x|,max|y| и в каком элементе), СКО, макс.
углы, невязка замыкания (для замкнутого кольца), грубая оценка числа бетатронных колебаний. Плюс
сводка «горячих» элементов (где поперечное отклонение максимально по всем частицам).

Запуск:  py -3.14 traj_analyzer.py <exp_id> [--csv] [--elem-top N] [--part K]
Сначала выгрузить локальные координаты:  py -3.14 track3.py --mode kick-phys --exp <N> --local
"""
import os
import sys
import argparse
import numpy as np
import matrix_traj_export as mte

RESULTS = os.path.join("workspace", "results")


def _rms(a):
    return float(np.sqrt(np.mean(a * a))) if len(a) else 0.0


def _osc_count(a):
    """Грубое число полупериодов колебаний: смены знака относительно среднего (для оценки бетатрона)."""
    if len(a) < 3:
        return 0
    c = a - np.mean(a)
    s = np.sign(c); s = s[s != 0]
    return int(np.sum(s[1:] != s[:-1])) if len(s) > 1 else 0


def analyze_particle(p, elem_names, closed, ring_len):
    """Метрики одной частицы в локальных координатах. Возвращает dict (величины в мм / мрад / м)."""
    S, x, xp, y, yp, elem = p["S"], p["x"], p["xp"], p["y"], p["yp"], p["elem"]
    n = len(S)
    ix, iy = int(np.argmax(np.abs(x))), int(np.argmax(np.abs(y)))
    r = dict(
        n=n, survived=p["survived"], S_end=float(S[-1]) if n else 0.0,
        max_x_mm=float(x[ix] * 1e3), max_x_S=float(S[ix]), max_x_elem=elem_names[int(elem[ix])],
        max_y_mm=float(y[iy] * 1e3), max_y_S=float(S[iy]), max_y_elem=elem_names[int(elem[iy])],
        rms_x_mm=_rms(x) * 1e3, rms_y_mm=_rms(y) * 1e3,
        max_xp_mrad=float(np.max(np.abs(xp)) * 1e3) if n else 0.0,
        max_yp_mrad=float(np.max(np.abs(yp)) * 1e3) if n else 0.0,
        osc_x=_osc_count(x), osc_y=_osc_count(y),
        end_elem=elem_names[int(elem[-1])] if n else "-",
    )
    # невязка замыкания (только для замкнутого кольца и выжившей частицы): старт vs финиш в локальных коорд.
    if closed and p["survived"] and n >= 2:
        r["clo_x_mm"] = float((x[-1] - x[0]) * 1e3); r["clo_xp_mrad"] = float((xp[-1] - xp[0]) * 1e3)
        r["clo_y_mm"] = float((y[-1] - y[0]) * 1e3); r["clo_yp_mrad"] = float((yp[-1] - yp[0]) * 1e3)
    else:
        r["clo_x_mm"] = r["clo_xp_mrad"] = r["clo_y_mm"] = r["clo_yp_mrad"] = None
    return r


def hot_elements(parts, elem_names, top):
    """По каждому элементу — макс |x| (мм) по всем частицам; вернуть top «горячих»."""
    peak = {}
    for p in parts:
        x, elem = np.abs(p["x"]), p["elem"]
        for e in np.unique(elem):
            m = float(np.max(x[elem == e]) * 1e3)
            if m > peak.get(int(e), -1):
                peak[int(e)] = m
    rows = sorted(peak.items(), key=lambda kv: -kv[1])[:top]
    return [(elem_names[e], v) for e, v in rows]


def _fmt(v, prec=3):
    return "-" if v is None else f"{v:.{prec}f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("exp", type=int, help="id эксперимента (читает results/trajectories_local_exp{N}.npz)")
    ap.add_argument("--elem-top", type=int, default=8, help="сколько «горячих» элементов показать")
    ap.add_argument("--part", type=int, default=None, help="детальный дамп одной частицы (индекс)")
    ap.add_argument("--csv", action="store_true", help="сохранить сводку в results/traj_analysis_exp{N}.csv")
    a = ap.parse_args()

    path = os.path.join(RESULTS, f"trajectories_local_exp{a.exp}.npz")
    if not os.path.exists(path):
        sys.exit(f"[!] нет {path} — сначала: py -3.14 track3.py --mode kick-phys --exp {a.exp} --local")
    d = mte.read_local_npz(path)
    parts, names, closed, ring = d["parts"], d["elem_names"], d["closed"], d["ring_len"]
    kind = "замкнутое кольцо" if closed else "открытый сегмент"
    print(f"=== Анализ траекторий (локальные координаты матядра) — эксперимент {a.exp} ===")
    print(f"установка {d['lattice']}, пучок {d['group']}, {kind}, элементов {len(names)-1}, "
          f"оборот {ring:.2f}м, частиц {len(parts)}, прореживание ε={d['thin_eps']*1e3:.2f}мм")

    res = [analyze_particle(p, names, closed, ring) for p in parts]
    surv = sum(r["survived"] for r in res)
    print(f"\nВыжило {surv}/{len(res)} (дошли до конца пути).")

    # таблица по частицам
    print("\n№  дошла  точек  S_кон,м   max|x|,мм (элемент)            СКОx  max|y|,мм            СКОy  max x',мрад  ~колеб")
    for i, r in enumerate(res):
        flag = "да " if r["survived"] else "НЕТ"
        loss = "" if r["survived"] else f"  ⟂потеря @ {r['end_elem']} (S={r['S_end']:.2f}м)"
        print(f"{i:<2} {flag}   {r['n']:>5}  {r['S_end']:>7.2f}  {r['max_x_mm']:>9.3f} "
              f"({r['max_x_elem'][:22]:<22}) {r['rms_x_mm']:>5.2f}  {r['max_y_mm']:>8.3f}  "
              f"{r['rms_y_mm']:>5.2f}  {r['max_xp_mrad']:>9.3f}  {r['osc_x']//2:>4}{loss}")

    # невязка замыкания
    clo = [r for r in res if r["clo_x_mm"] is not None]
    if clo:
        print("\nНевязка замыкания за оборот (старт vs финиш, локальные координаты):")
        for i, r in enumerate(res):
            if r["clo_x_mm"] is not None:
                print(f"  ч{i}: Δx={_fmt(r['clo_x_mm'])}мм Δx'={_fmt(r['clo_xp_mrad'])}мрад "
                      f"Δy={_fmt(r['clo_y_mm'])}мм Δy'={_fmt(r['clo_yp_mrad'])}мрад")

    # «горячие» элементы
    print(f"\n«Горячие» элементы (макс |x| по всем частицам), топ-{a.elem_top}:")
    for nm, v in hot_elements(parts, names, a.elem_top):
        print(f"  {v:>9.3f} мм  {nm}")

    # детальный дамп одной частицы
    if a.part is not None and 0 <= a.part < len(parts):
        p = parts[a.part]
        print(f"\nЧастица {a.part}: первые/последние 3 точки (S,м | x,мм | x',мрад | y,мм | y',мрад | элемент):")
        for k in list(range(min(3, len(p["S"])))) + list(range(max(0, len(p["S"]) - 3), len(p["S"]))):
            print(f"  S={p['S'][k]:8.3f}  x={p['x'][k]*1e3:9.3f}  x'={p['xp'][k]*1e3:8.3f}  "
                  f"y={p['y'][k]*1e3:9.3f}  y'={p['yp'][k]*1e3:8.3f}  {names[int(p['elem'][k])]}")

    # CSV-сводка (человеко-читаемая: ЗАПЯТАЯ-десятичная, ; разделитель — русский Excel)
    if a.csv:
        cpath = os.path.join(RESULTS, f"traj_analysis_exp{a.exp}.csv")
        with open(cpath, "w", encoding="utf-8") as f:
            f.write("частица;дошла;точек;S_конца_м;max_x_мм;элемент_max_x;СКО_x_мм;max_y_мм;СКО_y_мм;"
                    "max_xprime_мрад;max_yprime_мрад;колебаний_x;невязка_x_мм;невязка_y_мм\n")
            for i, r in enumerate(res):
                row = [i, "да" if r["survived"] else "нет", r["n"], r["S_end"], r["max_x_mm"],
                       r["max_x_elem"], r["rms_x_mm"], r["max_y_mm"], r["rms_y_mm"], r["max_xp_mrad"],
                       r["max_yp_mrad"], r["osc_x"] // 2,
                       r["clo_x_mm"] if r["clo_x_mm"] is not None else "",
                       r["clo_y_mm"] if r["clo_y_mm"] is not None else ""]
                f.write(";".join(str(c).replace(".", ",") if isinstance(c, float) else str(c) for c in row) + "\n")
        print(f"\n[CSV] {cpath}")


if __name__ == "__main__":
    main()
