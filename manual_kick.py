#!/usr/bin/env python3
"""manual_kick.py — ОДНОКРАТНЫЙ ручной кик токов/коэффициентов в sp_elements.csv.

Перетряхивает ВСЕ настраиваемые значения так же, как кик оптимизатора:
    value += value * amp * U(-1, 1)
и зажимает в ОДЗ (|ток| <= 100 A, |coeff| <= 3). Замороженные элементы
(корректоры, обмотки SP_I2/DQ_I3, ток PM-диполей) НЕ трогаются — ровно как в
parameter_groups_per_position тюнера. Делает бэкап, пишет sp_elements.csv,
пересобирает lattice_configs.csv. После этого можно запускать оптимизацию из
встряхнутого состояния.

Запуск:
    python manual_kick.py <amp>            # amp=0.01 -> ±1% от значения
    python manual_kick.py 0.02 --seed 7    # воспроизводимо
    python manual_kick.py 0.5 --abs-min 1.0  # +абсолютный пол ±1A (двигает ~нулевые токи)
"""
from __future__ import annotations
import argparse
import random
import shutil
import subprocess
import sys
from pathlib import Path

try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
from _common import workspace_dir  # noqa: E402

WS = workspace_dir()
SPE = WS / "tables" / "sp_elements.csv"
# Замороженные device_type для ТОКА (arg_val не кикаем) — как exclude_device_types тюнера.
FROZEN_CURRENT = {"DP_L_PM", "DP_sh_PM", "SP_I2_1A", "DQ_I3_1A",
                  "corrector_I1_1A", "corrector_I2_1A"}
# Для этих кикаем correct_coeff (как 2-я группа тюнера).
COEFF_DEVICES = {"DP_L_PM", "DP_sh_PM"}
CUR_LIM, COEFF_LIM = 100.0, 3.0


def _kick(v, amp, abs_min, rng, lim):
    u = rng.uniform(-1.0, 1.0)
    step = v * amp                                   # относительный кик
    if abs_min > 0.0 and abs(step) < abs_min:        # абсолютный пол (двигает ~нулевые)
        step = abs_min
    nv = v + step * u
    return max(-lim, min(lim, nv))


def main():
    ap = argparse.ArgumentParser(description="Однократный ручной кик sp_elements.")
    ap.add_argument("amp", type=float, help="амплитуда кика (доля; 0.01 = +-1%% от значения)")
    ap.add_argument("--seed", type=int, default=None, help="seed для воспроизводимости")
    ap.add_argument("--abs-min", type=float, default=0.0,
                    help="абсолютный пол кика (A) — двигает околонулевые токи (0=выкл)")
    ap.add_argument("--no-check", action="store_true",
                    help="не прогонять кольцо после кика (по умолчанию печатает новый рубеж)")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    df = pd.read_csv(SPE, sep=";", skipinitialspace=True)
    for col in ("arg_val", "correct_coeff"):          # pandas 3.x читает как str -> в float
        df[col] = pd.to_numeric(df[col], errors="coerce")
    bk = WS / "tables" / "backups" / "sp_elements_pre_manual_kick.csv"
    bk.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SPE, bk)

    n_cur, n_coeff, sum_dc = 0, 0, 0.0
    for i, r in df.iterrows():
        if str(r["sp_type"]).strip() != "standard":
            continue
        dev = str(r["device_type_id"]).strip()
        if dev not in FROZEN_CURRENT:                # ток (arg_val)
            v = float(r["arg_val"])
            nv = _kick(v, args.amp, args.abs_min, rng, CUR_LIM)
            sum_dc += abs(nv - v); n_cur += 1
            df.at[i, "arg_val"] = nv
        if dev in COEFF_DEVICES:                     # коэффициент (correct_coeff)
            c = float(r["correct_coeff"])
            nc = _kick(c, args.amp, 0.0, rng, COEFF_LIM)
            n_coeff += 1
            df.at[i, "correct_coeff"] = nc

    df.to_csv(SPE, sep=";", index=False)
    avg = sum_dc / n_cur if n_cur else 0.0
    print(f"Кик ±{args.amp*100:.1f}%"
          + (f" (абс.пол {args.abs_min}A)" if args.abs_min > 0 else "")
          + f": {n_cur} токов + {n_coeff} коэфф., средн.|Δток|={avg:.3f}A. Бэкап: {bk.name}")
    subprocess.run([sys.executable, "assemble_lattice_configs.py"], cwd=str(SCRIPT_DIR), check=True)
    if not args.no_check:                            # прогон кольца -> новый рубеж
        import continuation_tuner as ct
        ep = ct.load_endpoints()
        ct.track_full_ring(ep)
        loss, dev = ct.find_loss_gid(0.012)
        rs = f"gid {loss} ({ct.gid_to_sp(loss)})" if loss is not None else "ВСЁ кольцо"
        print(f"После кика: рубеж {rs}, max откл. {dev:.1f}мм.")
        if loss is not None and loss < 200:
            print("  ВНИМАНИЕ: рубеж рухнул -> кик слишком сильный. Откати: "
                  "copy /Y workspace\\tables\\backups\\sp_elements_pre_manual_kick.csv "
                  "workspace\\tables\\sp_elements.csv")
    print("Готово — запускай оптимизацию (или откати backups/sp_elements_pre_manual_kick.csv).")


if __name__ == "__main__":
    main()
