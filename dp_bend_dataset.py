# -*- coding: utf-8 -*-
"""Датасет «correct_coeff → угол изгиба θ» для постоянных диполей DP (DP_L_PM, DP_sh_PM).

DP — ПОСТОЯННЫЕ магниты: поле = карта × correct_coeff × field_direction (тока нет, nominal_arg='-'),
поэтому свипаем НЕ ток, а correct_coeff. Для каждого значения гоним осевой Boris-probe через
масштабированное поле и меряем изгиб θ (горизонтальный) и паразитный вертикальный θ_v.
Изгиб линеен по correct_coeff → чувствительность per_θ = θ/correct_coeff (мрад на единицу коэф.).

Выход: workspace/research/coupling_datasets/DP_<*>_bend.csv (correct_coeff; θ_гор_мрад; θ_верт_мрад)
+ dp_bend.png (θ vs correct_coeff обе кривые). Запуск: py -3.14 dp_bend_dataset.py"""
import os, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import field_quality as fq
import dq_sector_profile as dsp
from bip_matrix import QP

TB = os.path.join("workspace", "tables")
OUT = os.path.join("workspace", "research", "coupling_datasets"); os.makedirs(OUT, exist_ok=True)
HW = {r["device_type_id"]: r for r in csv.DictReader(open(os.path.join(TB, "hardware_library.csv"), encoding="utf-8"), delimiter=";")}
DEVS = ["DP_L_PM", "DP_sh_PM"]
GRID = list(np.linspace(-1e-3, 1e-3, 9))               # correct_coeff вокруг рабочего 5e-4 (симметрично через 0)


def bend(dev, coeff):
    """Осевой Boris через поле DP × coeff × field_direction → (θ_гор, θ_верт) в мрад."""
    r = HW[dev]; h, B = fq.load(os.path.join("workspace", "fields_bin", r["map_bin_path"]))
    fdir = float(r["field_direction"]); half = float(h["mxx"]); fn = dsp._field_fn(h, B)
    u = v = up = vp = 0.0; x = -half; ds = 2e-4; n = int(round(2 * half / ds))
    for _ in range(n):
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
        Bv = fn(np.array([x, u, v])) * coeff * fdir; up += -QP * Bv[2] * ds; vp += QP * Bv[1] * ds
        ps = math.sqrt(max(1e-12, 1 - up * up - vp * vp)); x += ds / 2; u += up * ds / (2 * ps); v += vp * ds / (2 * ps)
    return up * 1e3, vp * 1e3


def slope0(xs, ys):
    xs = np.array(xs); ys = np.array(ys)
    return float(np.sum(xs * ys) / (np.sum(xs * xs) + 1e-30))


def main():
    fig, ax = plt.subplots(figsize=(9, 5.5))
    print(f"{'диполь':10} {'per_θ (мрад/коэф)':>18}  рабочий θ (при 5e-4)")
    for dev in DEVS:
        rows = [(c, *bend(dev, c)) for c in GRID]
        p = os.path.join(OUT, f"{dev}_bend.csv")
        with open(p, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.writer(f, delimiter=";"); w.writerow(["correct_coeff", "theta_hor_mrad", "theta_vert_mrad"])
            for c, th, tv in rows:
                w.writerow([f"{c:.6e}", f"{th:.5f}", f"{tv:.5f}"])
        per = slope0([r[0] for r in rows], [r[1] for r in rows])
        op = float(HW[dev]["correct_coeff"]); th_op = per * op
        print(f"  {dev:10} {per:>+18.1f}  {th_op:+.4f} мрад")
        ax.plot([r[0] for r in rows], [r[1] for r in rows], "o-", label=f"{dev} (per_θ={per:+.0f})")
    ax.axhline(0, color="gray", lw=0.6); ax.axvline(0, color="gray", lw=0.6)
    ax.axvline(5e-4, color="red", ls=":", lw=1, label="рабочий correct_coeff=5e-4")
    ax.set_xlabel("correct_coeff"); ax.set_ylabel("угол изгиба θ (гор.), мрад")
    ax.set_title("Постоянные диполи DP: изгиб θ vs correct_coeff (линейно; новые карты)")
    ax.legend(); ax.grid(alpha=0.3); fig.tight_layout()
    png = os.path.join(OUT, "dp_bend.png"); fig.savefig(png, dpi=120); plt.close(fig)
    print(f"\n[OK] датасеты DP → {OUT}/DP_*_bend.csv ; график → {png}")


if __name__ == "__main__":
    main()
