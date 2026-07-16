# -*- coding: utf-8 -*-
"""ЕДИНЫЙ ДРАЙВЕР ТРЕКИНГА В 3 РЕЖИМАХ SynchrotronSim (полный оборот инжектированного пучка).

  --mode diff          режим 1: ДИФФУРНОЕ ядро (Boris по 3D-картам поля). Аргументы — токи (arg_val)
                       и коэф. диполей (correct_coeff) из lattice_configs. C++ SynchrotronTracker.exe.
  --mode kick-phys     режим 2: КИКОВОЕ (матричное) ядро, аргументы — проектные ФИЗКОЭФ. (mat_* из
                       lattice_configs — заполнены из LTE/тюна). ЕДИНСТВЕННЫЙ режим кика (токовый
                       архивирован: калибровки ток→физкоэф. предвычисляются в mat_* при сборке БД).

Режим 2 гоняет пучок глобальным треком от плоскости инжекции и обратно в неё (полный круг), режим 1
— C++ бинарником с финиш-плоскостью, поставленной в плоскость инжекции.

Запуск:  python track3.py --mode kick-phys --lattice 101 --group 1
         python track3.py --mode diff --lattice 101 --group 1 --exp 44
"""
import os
import sys
import argparse
import subprocess
import numpy as np
import beam_injection as bi
import matrix_traj_export as mte
import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
RESULTS = os.path.join("workspace", "results")


def run_kick(lattice, group, source, exp_id=None):
    """Режимы 2/3: глобальный трек инжектированного пучка полный оборот. Возврат сводки."""
    plane, ordered, states, info = bi.inject(lattice, group, source=source)
    survived = 0; residuals = []; trajs = []
    for st in states:
        s = tuple(st); x0, y0 = st[0], st[2]; lost_at = None
        pts = [mte._g(plane, 0.0, x0, y0)]
        try:
            s, pp = mte.drift_global(s, plane, ordered[0]); pts += pp
            s, pp = mte.body_global(ordered[0], s); pts += pp
            for j in range(1, len(ordered)):
                s, pp = mte.drift_global(s, ordered[j - 1], ordered[j]); pts += pp
                s, pp = mte.body_global(ordered[j], s); pts += pp
            s, pp = mte.drift_global(s, ordered[-1], plane); pts += pp   # возврат в плоскость
            survived += 1
            residuals.append(np.hypot(s[0] - x0, s[2] - y0))            # невязка замыкания (поперечн.)
        except mt.ParticleLost as e:
            lost_at = str(e)
        trajs.append(np.asarray(pts, float))
    if exp_id is not None:
        os.makedirs(RESULTS, exist_ok=True)
        mte.write_bin(trajs, os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin"))
        mte.ensure_experiment(exp_id, lattice)
    return dict(n=len(states), survived=survived, residuals=residuals, info=info,
                after=info["after_elem"], before=info["before_elem"])


def setup_diff_experiment(lattice, group, exp_id, dt=1e-11, skip_steps=300000):
    """Записать эксперимент для C++ диффур-трекера: финиш-плоскость = плоскость инжекции пучка.
    skip_steps (колонка max_turns, перепрофилирована в C++ как «пропустить N шагов до проверки
    плоскости») — пропускает старт (пучок на плоскости) и противоположное пересечение кольца,
    ловит только ВОЗВРАТ (полный оборот ≈ 3.7e5 шагов при dt=1e-11)."""
    plane, ordered, states, info = bi.inject(lattice, group)     # source не важен для геометрии плоскости
    c = plane.center; u = plane.uhat; v = plane.vhat; ax = plane.axis
    p1 = c + 1.0 * u + 1.0 * v; p2 = c - 1.0 * u + 1.0 * v; p3 = c + 0.0 * u - 1.0 * v  # 3 точки плоскости
    row = [str(exp_id), str(lattice), str(group), f"{dt:.2e}".replace("e", "E"), str(skip_steps), "50"]
    row += [f"{p1[0]:.6f}", f"{p1[1]:.6f}", f"{p1[2]:.6f}",
            f"{p2[0]:.6f}", f"{p2[1]:.6f}", f"{p2[2]:.6f}",
            f"{p3[0]:.6f}", f"{p3[1]:.6f}", f"{p3[2]:.6f}"]
    row += [f"{ax[0]:.6f}", f"{ax[1]:.6f}", f"{ax[2]:.6f}"]        # направление пересечения = ось пучка
    import csv
    path = os.path.join(TB, "experiments.csv")
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";")); head = rows[0]
    rows = [r for r in rows if not (r and r[0] == str(exp_id))]    # заменить, если был
    rows.append(row)
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    return info


def run_diff(lattice, group, exp_id):
    """Режим 1: C++ диффур-трекер по картам поля, полный оборот до плоскости инжекции."""
    info = setup_diff_experiment(lattice, group, exp_id)
    exe = os.path.join("build", "SynchrotronTracker.exe")
    if not os.path.exists(exe):
        exe = os.path.join("build", "Release", "SynchrotronTracker.exe")
    print(f"  запуск {exe} {exp_id} … (загрузка карт + интегрирование, может занять минуты)")
    r = subprocess.run([exe, str(exp_id)], capture_output=True, text=True, encoding="utf-8", errors="replace")
    tail = "\n".join((r.stdout or "").splitlines()[-12:])
    print(tail)
    # разбор траекторий: survivor — дошедший до финиш-плоскости (длинная траектория)
    binp = os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin")
    surv = None
    if os.path.exists(binp):
        trajs = _load_trajs(binp)
        lens = [len(t) for t in trajs]
        surv = sum(1 for L in lens if L > 100)
        print(f"  траекторий: {len(trajs)}, длины: {lens}")
    return dict(info=info, survived=surv, ret=r.returncode)


def _load_trajs(path):
    """results/trajectories_expN.bin → список массивов точек (uint32 N + N×3×float64 на частицу)."""
    out = []
    with open(path, "rb") as f:
        data = f.read()
    off = 0
    while off + 4 <= len(data):
        n = int(np.frombuffer(data, np.uint32, 1, off)[0]); off += 4
        pts = np.frombuffer(data, np.float64, n * 3, off).reshape(n, 3); off += n * 3 * 8
        out.append(pts)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["diff", "kick-phys"])
    ap.add_argument("--lattice", type=int, default=101)
    ap.add_argument("--group", type=int, default=1)
    ap.add_argument("--exp", type=int, default=44)
    a = ap.parse_args()

    print(f"=== SynchrotronSim: трекинг, режим '{a.mode}', решётка {a.lattice}, группа {a.group} ===")
    if a.mode == "diff":
        res = run_diff(a.lattice, a.group, a.exp)
        s = res["survived"]
        print(f"\nРежим 1 (диффур): выжило {s if s is not None else '?'} частиц "
              f"(плоскость '{res['info']['after_elem']}'→'{res['info']['before_elem']}').")
    else:
        res = run_kick(a.lattice, a.group, "physical", exp_id=a.exp)
        rr = res["residuals"]
        print(f"\nРежим 2 ({a.mode}): выжило {res['survived']}/{res['n']} "
              f"частиц полный оборот.")
        print(f"  плоскость инжекции: '{res['after']}' → '{res['before']}'")
        if rr:
            print(f"  невязка замыкания (поперечн.): макс {max(rr)*1e3:.3f} мм, средн {np.mean(rr)*1e3:.3f} мм")


if __name__ == "__main__":
    main()
