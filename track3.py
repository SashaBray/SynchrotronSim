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
os.environ.setdefault("SS_SECTOR", "1")     # СЕКТОРНАЯ модель (m.body, бит-в-бит elegant) — нужна для локальных
#                                             координат (x=отклонение от опорной орбиты); SS_SECTOR=0 форсит боксовую
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


def run_kick(lattice, group, source, exp_id=None, n_turns=1, thin_eps=0.0, settings=None, save_local=False):
    """Режим 2: глобальный трек инжектированного пучка старт→финиш, n_turns оборотов. Возврат сводки.
    Замкнут (финиш=старт) ⇒ оборот(ы) + невязка; открыт (финиш≠старт) ⇒ путь до конца сегмента (n_turns=1).
    thin_eps>0 (м) ⇒ потоковое прореживание траектории по кривизне (лёгкие длинные многооборотные треки).
    save_local ⇒ доп. писать trajectories_local_exp{id}.npz — ЛОКАЛЬНЫЕ координаты (СЕКТОРНАЯ модель m.body,
    x=отклонение от опорной орбиты) для анализатора; глоб. .bin остаётся боксовым (body_global)."""
    plane, ordered, states, info = bi.inject(lattice, group, source=source, settings=settings)
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
        if save_local:                                     # ЛОКАЛЬНЫЕ координаты — отдельным секторным треком (см. mte)
            r = mte.local_trajectories(lattice, group, source=source, settings=settings, n_turns=n_turns)
            mte.write_local_npz(os.path.join(RESULTS, f"trajectories_local_exp{exp_id}.npz"), r)
        mte.ensure_experiment(exp_id, lattice)
    return dict(n=len(states), survived=survived, residuals=residuals, info=info, closed=closed,
                turns=n_turns, lost_turns=lost_turns,
                after=info["after_elem"], before=info["before_elem"])


def run_long(lattice, group, exp_id, n_turns, rf=None, x0=0.0, y0=0.0, settings=None):
    """ДЛИННЫЙ турн-by-turn трек СЕКТОРНОЙ моделью (m.body) для анализа затухания/огибающей — минуя боксовый
    глобальный путь run_kick (диполи в нём раздуваются, а с SR тем более бессмысленны). Пишет
    trajectories_local_exp{id}.npz (mode='long') для traj_analyzer.py. Радиация — env SS_SR/SS_ISR (ставит main);
    rf='ideal' — компенсация потерь. x0/y0 (м) ⇒ ОДНА частица со смещением (сканы апертуры) вместо пучка группы."""
    states_override = [(x0, 0.0, y0, 0.0, 0.0, 0.0)] if (x0 or y0) else None
    r = mte.long_trajectories(lattice, group, source="physical", settings=settings,
                              n_turns=n_turns, states_override=states_override, rf=rf)
    os.makedirs(RESULTS, exist_ok=True)
    mte.write_local_npz(os.path.join(RESULTS, f"trajectories_local_exp{exp_id}.npz"), r)
    mte.ensure_experiment(exp_id, lattice)
    return r


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
    """experiments.csv[exp_id] → (lattice_id, group_id) или (None, None).
    ЕДИНЫЙ конфиг эксперимента (установка+пучок+финиш-плоскость), который читают ОБА ядра.
    ПАРАМЕТРЫ ЯДРА того же опыта читает matrix_tracker.Settings.from_experiment(exp_id)."""
    path = os.path.join(TB, "experiments.csv")
    if not os.path.exists(path):
        return None, None
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";"):
        if r.get("exp_id") == str(exp_id):
            return int(r["lattice_id"]), int(r["group_id"])
    return None, None


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
    ap.add_argument("--local", action="store_true", help="доп. писать локальные координаты матядра (S,x,x',y,y',ct,δ) в trajectories_local_exp{id}.npz — для traj_analyzer.py")
    # ДЛИННЫЙ турн-by-turn режим + РАДИАЦИЯ (синхротронное излучение и затухание)
    ap.add_argument("--long", action="store_true", help="длинный турн-by-turn трек (сектор m.body, много оборотов) → trajectories_local_exp{id}.npz (mode=long) для traj_analyzer; минует боксовый глоб. путь")
    ap.add_argument("--sr", action="store_true", help="распределённое синхротронное излучение (потери энергии в изгибах, ∝γ⁴ по ПРОЕКТНОЙ энергии --design-mev)")
    ap.add_argument("--isr", action="store_true", help="+ квантовое возбуждение ISR (гауссов шум фотонов; равновесный эмиттанс/энергоразброс). Требует --sr")
    ap.add_argument("--rf", choices=["none", "ideal"], default="none", help="ВЧ-резонатор (компенсация потерь SR): ideal — восстановление δ за оборот ⇒ поперечное затухание. Требует --long")
    ap.add_argument("--design-mev", type=float, default=None, help="ПРОЕКТНАЯ энергия пучка для SR, МэВ (СИЛА 6000; трекер гоняет сурогатный γ≈7, а SR ∝ γ⁴ берётся по проектной)")
    ap.add_argument("--x0-mm", type=float, default=0.0, help="(--long) ОДНА частица со смещением по x, мм (скан апертуры/затухания) вместо пучка группы")
    ap.add_argument("--y0-mm", type=float, default=0.0, help="(--long) смещение одной частицы по y, мм")
    # ПАРАМЕТРЫ ЯДРА: дефолт ← experiments.csv[exp] ← этот аргумент (приоритет справа налево)
    ap.add_argument("--loss-mm", type=float, default=None,
                    help="порог потери пучка, мм (иначе experiments.csv[loss_radius_mm], иначе дефолт ядра)")
    ap.add_argument("--plane-same-tol-mm", type=float, default=None,
                    help="плоскости ближе этого — одна, без пролёта, мм (иначе experiments.csv[plane_same_tol_mm])")
    ap.add_argument("--plane-parallel-tol", type=float, default=None,
                    help="|axis_A × axis_B| ниже — нормали параллельны (иначе experiments.csv[plane_parallel_tol])")
    ap.add_argument("--n-slices", type=int, default=None,
                    help="срезов тела, ФОЛБЭК (приоритет у lattice_configs[n_slices]; иначе experiments.csv[n_slices])")
    a = ap.parse_args()

    # РАДИАЦИЯ: флаги → env (ядро читает SS_SR/SS_ISR/SS_DESIGN_MEV при сборке магнитов). Ставить ДО трекинга.
    if a.sr:
        os.environ["SS_SR"] = "1"
    if a.isr:
        os.environ["SS_ISR"] = "1"
    if a.design_mev is not None:
        os.environ["SS_DESIGN_MEV"] = str(a.design_mev)
    if a.isr and not a.sr:
        ap.error("--isr (квантовое возбуждение) требует --sr")
    if a.rf != "none" and not a.long:
        ap.error("--rf требует --long (ВЧ применяется в турн-by-turn цикле)")

    # ЕДИНЫЙ ИСТОЧНИК: установка+пучок берутся из experiments.csv[exp]; --lattice/--group переопределяют
    # (и создают новую строку, если эксперимента ещё нет). Оба ядра стартуют из ОДНОГО конфига.
    lat_db, grp_db = load_experiment(a.exp)
    lattice = a.lattice if a.lattice is not None else lat_db
    group = a.group if a.group is not None else grp_db
    settings = mt.Settings.from_experiment(
        a.exp,
        loss_radius_m=(a.loss_mm / 1e3) if a.loss_mm is not None else None,
        plane_same_tol_m=(a.plane_same_tol_mm / 1e3) if a.plane_same_tol_mm is not None else None,
        plane_parallel_tol=a.plane_parallel_tol, n_slices=a.n_slices)
    if lattice is None or group is None:
        ap.error(f"эксперимент {a.exp} не найден в experiments.csv — задайте --lattice и --group, чтобы создать его")

    print(f"=== SynchrotronSim: трекинг, режим '{a.mode}', эксперимент {a.exp} (установка {lattice}, пучок {group}) ===")
    if a.mode == "diff":
        res = run_diff(lattice, group, a.exp)
        s = res["survived"]
        print(f"\nРежим 1 (диффур): выжило {s if s is not None else '?'} частиц "
              f"(плоскость '{res['info']['after_elem']}'→'{res['info']['before_elem']}').")
    elif a.long:                                           # ДЛИННЫЙ турн-by-turn (сектор m.body) + радиация/ВЧ
        turns = max(1, a.turns)
        author_experiment(lattice, group, a.exp, loss_mm=settings.loss_radius_m * 1e3)
        print(f"  параметры ядра: {settings}")
        rf = None if a.rf == "none" else a.rf
        sr_txt = ("SR " + ("+ISR " if a.isr else "") + f"[проект {os.environ.get('SS_DESIGN_MEV', '6000')} МэВ] " if a.sr else "")
        who = f"1 частица x0={a.x0_mm}мм y0={a.y0_mm}мм" if (a.x0_mm or a.y0_mm) else "пучок группы"
        print(f"  ДЛИННЫЙ режим: {turns} об., {sr_txt}ВЧ={a.rf}, {who}")
        r = run_long(lattice, group, a.exp, turns, rf=rf, x0=a.x0_mm / 1e3, y0=a.y0_mm / 1e3, settings=settings)
        for i, (tt, sv) in enumerate(zip(r["turnbyturn"], r["surv_turns"])):
            if len(tt) < 4:
                print(f"  частица {i}: погибла на об.{sv}"); continue
            w = min(len(tt), 80)                           # оконная амплитуда √2·σ на старте и в конце
            ax0 = np.sqrt(2) * np.std(tt[:w, 0]); ax1 = np.sqrt(2) * np.std(tt[-w:, 0])
            ay0 = np.sqrt(2) * np.std(tt[:w, 2]); ay1 = np.sqrt(2) * np.std(tt[-w:, 2])
            print(f"  частица {i}: дошла {sv}/{turns} об.; ампл x {ax0*1e3:.3f}→{ax1*1e3:.3f} мм,  y {ay0*1e3:.3f}→{ay1*1e3:.3f} мм")
        print(f"  → results/trajectories_local_exp{a.exp}.npz  (traj_analyzer.py {a.exp} — огибающая/фазовый портрет)")
    else:
        turns = max(1, a.turns)
        thin_eps = (a.thin * 1e-3) if a.thin > 0 else (5e-4 if turns > 1 else 0.0)   # мм→м; авто 0.5мм для многооборота
        # записать/обновить конфиг (старт+финиш-плоскость+порог потери) — чтобы оба ядра читали одно
        author_experiment(lattice, group, a.exp, loss_mm=settings.loss_radius_m * 1e3)
        print(f"  параметры ядра: {settings}")
        res = run_kick(lattice, group, "physical", exp_id=a.exp, n_turns=turns, thin_eps=thin_eps,
                       settings=settings, save_local=a.local)
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
        if a.local:
            print(f"  локальные координаты матядра → results/trajectories_local_exp{a.exp}.npz (для traj_analyzer.py)")


if __name__ == "__main__":
    main()
