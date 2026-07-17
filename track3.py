# -*- coding: utf-8 -*-
"""ЕДИНЫЙ ДРАЙВЕР ТРЕКИНГА В 3 РЕЖИМАХ SynchrotronSim (полный оборот инжектированного пучка).

  --mode diff          режим 1: ДИФФУРНОЕ ядро (Boris по 3D-картам поля). Аргументы — токи (arg_val)
                       и коэф. диполей (correct_coeff) из lattice_configs. C++ SynchrotronTracker.exe.
  --mode kick-phys     режим 2: КИКОВОЕ (матричное) ядро, аргументы — проектные ФИЗКОЭФ. (mat_* из
                       lattice_configs — заполнены из LTE/тюна). ЕДИНСТВЕННЫЙ режим кика (токовый
                       архивирован: калибровки ток→физкоэф. предвычисляются в mat_* при сборке БД).

Режим 2 гоняет пучок глобальным треком от плоскости инжекции и обратно в неё (полный круг), режим 1
— C++ бинарником с финиш-плоскостью, поставленной в плоскость инжекции.

Установка+пучок+финиш-плоскость берутся из experiments.csv[exp] (ЕДИНЫЙ конфиг эксперимента для обоих ядер).

Запуск:  python track3.py --mode kick-phys --exp 44        (установка/пучок — из experiments.csv[44])
         python track3.py --mode diff --exp 44
         python track3.py --mode kick-phys --exp 99 --lattice 1 --group 1   (создать/переопределить эксперимент)
"""
import os
import sys
import csv
import argparse
import subprocess
import numpy as np
import beam_injection as bi
import matrix_traj_export as mte
import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
RESULTS = os.path.join("workspace", "results")


def run_kick(lattice, group, source, exp_id=None, n_turns=1, thin_eps=0.0, loss_radius=None):
    """Режим 2: глобальный трек инжектированного пучка старт→финиш, n_turns оборотов. Возврат сводки.
    Замкнут (финиш=старт) ⇒ оборот(ы) + невязка; открыт (финиш≠старт) ⇒ путь до конца сегмента (n_turns=1).
    thin_eps>0 (м) ⇒ потоковое прореживание траектории по кривизне (лёгкие длинные многооборотные треки)."""
    plane, ordered, states, info = bi.inject(lattice, group, source=source, loss_radius=loss_radius)
    finish = info.get("finish", plane); closed = info.get("closed", True)
    if n_turns > 1 and not closed:
        n_turns = 1                                        # открытый сегмент не зациклить
    survived = 0; residuals = []; trajs = []; lost_turns = []
    for st in states:
        s = tuple(st); x0, y0 = st[0], st[2]; lost_turn = None; turn = 0
        th = mte.Thinner(thin_eps); th.add(mte._g(plane, 0.0, x0, y0))
        try:
            for turn in range(n_turns):
                s, pp = mte.drift_global(s, plane, ordered[0]);            [th.add(p) for p in pp]
                s, pp = mte.body_global(ordered[0], s);                    [th.add(p) for p in pp]
                for j in range(1, len(ordered)):
                    s, pp = mte.drift_global(s, ordered[j - 1], ordered[j]); [th.add(p) for p in pp]
                    s, pp = mte.body_global(ordered[j], s);                  [th.add(p) for p in pp]
                s, pp = mte.drift_global(s, ordered[-1], finish);          [th.add(p) for p in pp]   # → финиш (замкнут: =старт)
            survived += 1
            residuals.append(np.hypot(s[0] - x0, s[2] - y0) if closed else np.nan)  # невязка только для замкнутого
        except mt.ParticleLost:
            lost_turn = turn + 1
        trajs.append(th.result()); lost_turns.append(lost_turn)
    if exp_id is not None:
        os.makedirs(RESULTS, exist_ok=True)
        mte.write_bin(trajs, os.path.join(RESULTS, f"trajectories_exp{exp_id}.bin"))
        mte.ensure_experiment(exp_id, lattice)
    return dict(n=len(states), survived=survived, residuals=residuals, info=info, closed=closed,
                turns=n_turns, lost_turns=lost_turns,
                after=info["after_elem"], before=info["before_elem"])


def author_experiment(lattice, group, exp_id, dt=1e-11, skip_steps=300000, loss_mm=12.0):
    """АВТОРИНГ конфига эксперимента в experiments.csv (установка+пучок+финиш-плоскость) — ЕДИНЫЙ
    источник, который читают ОБА ядра. Финиш-плоскость = плоскость инжекции пучка (полный оборот).
    skip_steps (max_turns, перепрофилир. в C++ как «пропустить N шагов до проверки плоскости») —
    пропускает старт и противоположное пересечение, ловит ВОЗВРАТ (≈3.7e5 шагов при dt=1e-11)."""
    plane, ordered, states, info = bi.inject(lattice, group)     # source не важен для геометрии плоскости
    fin = info.get("finish", plane)                              # в БД пишем ФИНИШ-плоскость (её ловит диффур)
    closed = info.get("closed", True)
    c = fin.center; u = fin.uhat; v = fin.vhat; ax = fin.axis
    p1 = c + 1.0 * u + 1.0 * v; p2 = c - 1.0 * u + 1.0 * v; p3 = c + 0.0 * u - 1.0 * v  # 3 точки плоскости
    skip = skip_steps if closed else 0                           # открыт: финиш≠старт, пересекается 1 раз — без пропуска
    row = [str(exp_id), str(lattice), str(group), f"{dt:.2e}".replace("e", "E"), str(skip), "50"]
    row += [f"{p1[0]:.6f}", f"{p1[1]:.6f}", f"{p1[2]:.6f}",
            f"{p2[0]:.6f}", f"{p2[1]:.6f}", f"{p2[2]:.6f}",
            f"{p3[0]:.6f}", f"{p3[1]:.6f}", f"{p3[2]:.6f}"]
    row += [f"{ax[0]:.6f}", f"{ax[1]:.6f}", f"{ax[2]:.6f}"]        # направление пересечения = ось пучка
    row += [f"{loss_mm:.1f}"]                                      # ПОРОГ ПОТЕРИ ПУЧКА (мм) — параметр опыта
    import csv
    path = os.path.join(TB, "experiments.csv")
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";")); head = rows[0]
    rows = [r for r in rows if not (r and r[0] == str(exp_id))]    # заменить, если был
    rows.append(row)
    with open(path, "w", newline="", encoding="utf-8") as f:
        csv.writer(f, delimiter=";").writerows(rows)
    return info


def load_experiment(exp_id):
    """experiments.csv[exp_id] → (lattice_id, group_id, loss_radius_m) или (None, None, None).
    ЕДИНЫЙ конфиг эксперимента (установка+пучок+финиш-плоскость+порог потери), который читают ОБА ядра.
    loss_radius_m — из столбца loss_radius_mm (мм→м); None ⇒ ядро возьмёт свой дефолт."""
    path = os.path.join(TB, "experiments.csv")
    if not os.path.exists(path):
        return None, None, None
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";"):
        if r.get("exp_id") == str(exp_id):
            lr = (r.get("loss_radius_mm") or "").strip()
            return int(r["lattice_id"]), int(r["group_id"]), (float(lr) / 1e3 if lr else None)
    return None, None, None


def run_diff(lattice, group, exp_id):
    """Режим 1: C++ диффур-трекер по картам поля, полный оборот до плоскости инжекции."""
    info = author_experiment(lattice, group, exp_id)
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
    ap.add_argument("--exp", type=int, required=True, help="id эксперимента в experiments.csv (установка+пучок+финиш-плоскость)")
    ap.add_argument("--lattice", type=int, default=None, help="переопределить установку (иначе — из experiments.csv[exp])")
    ap.add_argument("--group", type=int, default=None, help="переопределить пучок (иначе — из experiments.csv[exp])")
    ap.add_argument("--turns", type=int, default=1, help="число оборотов (kick-phys, только замкнутый; >1 ⇒ авто-прореживание)")
    ap.add_argument("--thin", type=float, default=0.0, help="прореживание траектории по кривизне, ε в мм (0=авто: выкл для 1 об., 0.5мм для многих)")
    ap.add_argument("--loss-mm", type=float, default=None, help="порог потери пучка (мм); иначе из experiments.csv[loss_radius_mm]")
    a = ap.parse_args()

    # ЕДИНЫЙ ИСТОЧНИК: установка+пучок берутся из experiments.csv[exp]; --lattice/--group переопределяют
    # (и создают новую строку, если эксперимента ещё нет). Оба ядра стартуют из ОДНОГО конфига.
    lat_db, grp_db, lr_db = load_experiment(a.exp)
    lattice = a.lattice if a.lattice is not None else lat_db
    group = a.group if a.group is not None else grp_db
    loss_r = (a.loss_mm / 1e3) if a.loss_mm is not None else lr_db   # порог потери: CLI > БД > дефолт ядра
    if lattice is None or group is None:
        ap.error(f"эксперимент {a.exp} не найден в experiments.csv — задайте --lattice и --group, чтобы создать его")

    print(f"=== SynchrotronSim: трекинг, режим '{a.mode}', эксперимент {a.exp} (установка {lattice}, пучок {group}) ===")
    if a.mode == "diff":
        res = run_diff(lattice, group, a.exp)
        s = res["survived"]
        print(f"\nРежим 1 (диффур): выжило {s if s is not None else '?'} частиц "
              f"(плоскость '{res['info']['after_elem']}'→'{res['info']['before_elem']}').")
    else:
        turns = max(1, a.turns)
        thin_eps = (a.thin * 1e-3) if a.thin > 0 else (5e-4 if turns > 1 else 0.0)   # мм→м; авто 0.5мм для многооборота
        # записать/обновить конфиг (старт+финиш-плоскость+порог потери) — чтобы оба ядра читали одно
        author_experiment(lattice, group, a.exp, loss_mm=(loss_r * 1e3 if loss_r else 12.0))
        res = run_kick(lattice, group, "physical", exp_id=a.exp, n_turns=turns, thin_eps=thin_eps,
                       loss_radius=loss_r)
        kind = "замкнут" if res["closed"] else "ОТКРЫТ (старт→финиш сегмента)"
        print(f"\nРежим 2 ({a.mode}): {res['turns']} об., выжило {res['survived']}/{res['n']} частиц — путь {kind}.")
        print(f"  плоскость старта: '{res['after']}' → '{res['before']}'")
        if res["turns"] > 1:
            lost = [f"частица {i} на об.{t}" for i, t in enumerate(res["lost_turns"]) if t]
            print(f"  потери по оборотам: {'; '.join(lost) if lost else 'нет, все дошли'}")
        rr = [r for r in res["residuals"] if not np.isnan(r)]
        if res["closed"] and rr:
            print(f"  невязка замыкания за оборот (поперечн.): макс {max(rr)*1e3:.3f} мм, средн {np.mean(rr)*1e3:.3f} мм")
        if thin_eps > 0:
            print(f"  траектория прорежена по кривизне (ε={thin_eps*1e3:.2f}мм) — лёгкий bin для визуализации")


if __name__ == "__main__":
    main()
