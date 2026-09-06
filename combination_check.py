# -*- coding: utf-8 -*-
"""Проверка КОМБИНАЦИЙ: предсказание модели против измерения на картах с несколькими обмотками.

Комбинации в подгонку не входят — веса ищутся только по одиночным обмоткам, — поэтому это
независимая проверка. По оси абсцисс измеренное, по ординат предсказанное; идеал — прямая y = x.
Одиночные показаны серым для сравнения: на них модель подогнана и обязана лежать хорошо.

Метрика — сумма расхождений к сумме модулей измеренного, и считается она по ВСЕМ картам. Этим
она отличается от ошибки тока |λ−1|, которая на пологом отклике не имеет решения и молча
выбрасывает карту (см. `mec_sigma.current_errors`): здесь ни одна карта не теряется, и число
нельзя улучшить тем, что модель перестала отвечать.

    py -3.14 combination_check.py corrector
    py -3.14 combination_check.py DQ
"""
import json, math, os, sys, warnings
sys.path.insert(0, os.path.abspath("."))
warnings.filterwarnings("ignore")
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import branch_network as bn
import mec_sigma as ms

EL = sys.argv[1] if len(sys.argv) > 1 else "corrector"
R = 8e-3
MAT_SIGN = {1: -1.0, 2: 1.0, 3: 1.0}
ORDERS = ((1, "B", "Тл"), (2, "G", "Тл/м"), (3, "S", "Тл/м²"))
PARTS = (("норм.", lambda z: z.imag),) if EL == "DQ" else \
        (("норм.", lambda z: z.imag), ("скош.", lambda z: z.real))
meta = json.load(open(f"workspace/datasets/currents_to_kicks/{EL}.meta.json", encoding="utf-8"))
L, EF = meta["L_used_for_mat_m"], meta["energy_factor"]
saved = ms.load()

cols = [(n, name, unit, pn, pf) for n, name, unit in ORDERS for pn, pf in PARTS]
fig, ax = plt.subplots(1, len(cols), figsize=(3.5 * len(cols), 4.3), squeeze=False)
summary = []
for k, (n, name, unit, pn, pf) in enumerate(cols):
    cir, wv = saved[(EL, n)]
    sc = MAT_SIGN[n] * math.factorial(n - 1) / (R ** (n - 1) * L * EF)
    a = ax[0][k]
    for grp, col, lab, ms_ in ((1, "#bbbbbb", "одиночные (подгонка)", 5),
                               (2, "#c0392b", "комбинации (проверка)", 6)):
        xs, ys = [], []
        for c, m, _, _ in bn.read_maps(EL, n):
            if (len(c) == 1) != (grp == 1):
                continue
            xs.append(pf(complex(m)) * sc)
            ys.append(pf(complex(ms.basis(cir, c) @ wv)) * sc)
        a.plot(xs, ys, "o", color=col, ms=ms_, mec="white", mew=0.6,
               label=lab if k == 0 else None, zorder=4 if grp == 2 else 3)
        if grp == 2 and xs:
            d = np.abs(np.array(ys) - np.array(xs))
            rel = 100 * d.sum() / np.abs(xs).sum()
            summary.append((f"{name} {pn}", rel, len(xs)))
            a.set_title(f"{name} {pn}\nкомбинации: расхождение {rel:.1f}%", fontsize=10)
    lim = max(abs(np.array(a.get_xlim())).max(), abs(np.array(a.get_ylim())).max())
    a.plot([-lim, lim], [-lim, lim], "-", color="#333", lw=0.9, zorder=1)
    a.set_xlim(-lim, lim); a.set_ylim(-lim, lim)
    a.set_xlabel(f"измерено, {unit}", fontsize=9)
    if k == 0:
        a.set_ylabel(f"предсказано моделью, {unit}", fontsize=9.5)
    a.grid(alpha=0.25); a.tick_params(labelsize=8); a.set_aspect("equal")
fig.legend(loc="lower center", ncol=2, fontsize=10, frameon=False,
           bbox_to_anchor=(0.5, 0.005))
fig.suptitle(f"{EL}: комбинации обмоток — независимая проверка модели", fontsize=13)
fig.tight_layout(rect=(0, 0.12, 1, 0.94))
out = f"workspace/results/{EL}_combinations_check.png"
fig.savefig(out, dpi=135)
for nm, rel, cnt in summary:
    print(f"  {nm:10} расхождение {rel:6.1f}%  по {cnt} комбинациям")
print(out)
