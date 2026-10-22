# -*- coding: utf-8 -*-
"""ПУНКТ 1 статики медленной системы: ЗАПАС КОРРЕКТОРОВ (authority / насыщение).

Сколько углов корректоров реально тратится на коррекцию орбиты и как спрос растёт с амплитудой
рассогласования — это бюджет дрейфа, который медленная ОС способна снять штатными корректорами.
Сравниваем с аппаратным максимумом кика (ток ≤2.5 А рабочая область / 125 А imax в калибровке).

Сценарии: проектные допуски (element_errors.csv) + только-поперечный свип амплитуды (прокси дрейфа).
Корректор kind=kicker: горизонт — главная обмотка (554 мкрад/А), вертикаль — тримы (~88 мкрад/А).
Запуск: py -3.14 corrector_authority.py
"""
import os, sys, tempfile, shutil
os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import facility_report as fr, errors as errmod, orbit_correct as oc, ring_analysis as ra
from bip_matrix import QP

NSEED, CUTOFF = 15, 2.5
L_C = 0.1
# Прямая матрица 5-обмоточного корректора [B;Bv] (горизонт/вертикаль поле) × [main,c1,c2,c3,c4], coef/А
F2 = np.array([
    [6.525763430073686e-05, 4.379117802766501e-06, 4.955703826510264e-06, 4.38375542362199e-06, 4.967217474297807e-06],  # B (гор)
    [1.6559167696844574e-09, -8.504657643801892e-06, 1.0385492912758759e-05, 8.485113787961437e-06, -1.0384265950213528e-05],  # Bv (вер)
])
F2PINV = np.linalg.pinv(F2)                               # min-norm токи обмоток под целевое (B,Bv)
WIND = ["main", "c1", "c2", "c3", "c4"]
KICK_H_PERA = abs(6.525763430073686e-05 * QP * L_C)       # рад/А, горизонт (главная обмотка)
KINDS = ("dipole", "dq", "quad", "sext", "oct")


def winding_currents(thx, thy):
    """Кик (thx,thy) [рад] одного корректора → токи 5 обмоток [А] (min-norm). thx,thy — массивы по инстансам."""
    Bt = thx / (QP * L_C); Bvt = thy / (QP * L_C)         # целевые поля B,Bv (θ = поле·QP·L)
    return F2PINV @ np.vstack([Bt, Bvt])                  # (5, N)
inf = dict(fr.info_of(26)); inf["_id"] = 26
ideal = fr._ring(inf, nsl=8); resp = oc.cached_response_matrix(ideal, 26); ideal.close()


def tol_trans(dxy):
    return {k: dict(dx=dxy, dy=dxy, dz=0, dyaw=0, dpitch=0, droll=0, field_pct=0) for k in KINDS}


def run(tol, label):
    thx_all, thy_all = [], []; nok = 0
    root = tempfile.mkdtemp(prefix="auth_")
    try:
        for s in range(NSEED):
            ws = os.path.join(root, str(s)); errmod.perturbed_workspace(26, 9000 + s, tol, ws, CUTOFF)
            r = fr._ring(inf, ws=ws, nsl=8)
            try:
                res = oc.correct(r, resp, aperture=inf["loss_m"])
            except Exception:
                res = None
            if res is not None:
                nok += 1; thx_all.append(res["thx"]); thy_all.append(res["thy"])   # ЗНАКОВЫЕ
            r.close()
    finally:
        shutil.rmtree(root, ignore_errors=True)
    if not thx_all:
        print(f"  {label:16}: исправлено 0/{NSEED} (стена бутстрапа)"); return
    thx = np.concatenate(thx_all); thy = np.concatenate(thy_all)
    I = winding_currents(thx, thy)                        # (5, N) токи обмоток [А]
    Iabs = np.abs(I)
    perw_max = Iabs.max(axis=1)                            # макс ток по каждой обмотке
    imax_coil = int(perw_max.argmax()); imax_val = perw_max[imax_coil]
    # для контекста — сами кики
    print(f"  {label:16}: испр {nok:2}/{NSEED} | кик ГОР max={np.abs(thx).max()*1e6:6.0f} ВЕР max={np.abs(thy).max()*1e6:6.0f} мкрад "
          f"| ТОК обмоток макс: " + " ".join(f"{WIND[k]}={perw_max[k]:.2f}" for k in range(5))
          + f" А → предел {WIND[imax_coil]} {imax_val:.2f} А ({imax_val/2.5*100:.0f}% от 2.5А)")


print(f"ЗАПАС КОРРЕКТОРОВ (реш.26, {NSEED} машин/сценарий). Ток обмоток через min-norm инверс 5-обмоточного корректора; предел 2.5 А/обмотку.\n", flush=True)
print("Проектные допуски и свип поперечного сдвига (прокси дрейфа):", flush=True)
run(errmod.load_tolerances(), "проектные")
for dxy in (5e-6, 10e-6, 20e-6, 30e-6, 50e-6):
    run(tol_trans(dxy), f"trans {int(dxy*1e6)}мкм")
print("\n[готово]", flush=True)
