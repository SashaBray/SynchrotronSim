# -*- coding: utf-8 -*-
"""correct_coeff диполей DP в кольце — из АРОЧНОГО теста dp_arc_<dev>.csv (НЕ через устаревший per_B
coupling_matrix). Диполь DP — постоянный магнит (input_type=coeff): θ линейна по correct_coeff с наклоном
per_θ (мрад/coeff), измеренным симметричным арочным проездом на ТЕКУЩИХ картах кольца. Инверсия:
  correct_coeff = θ_проект / per_θ,   θ_проект = QP·mat_B·Leff·1e3  (мрад, из lattice_configs mat_B).

Почему не per_B: recalibrate_coupling ЯВНО пропускает диполи ⇒ их per_B устарел (старые карты, per_θ≈5980).
Новые синт.карты DP_*_L### дают per_θ≈5.8/0.9 ⇒ старый correct_coeff≈0.0016 гнёт ~0 (диффур шёл прямо).
Арочный per_θ согласован с C++-диффуром (Python-Boris арки ≡ наблюдаемый недогиб при 0.0016).

Запуск ПОСЛЕ phys_to_current (перекрывает его устаревший DP-correct_coeff):
  py -3.14 dp_correct_from_arc.py [lattice=101]"""
import os, csv, sys, shutil
import numpy as np
from collections import defaultdict
from bip_matrix import QP
import matrix_tracker as mt

TB = os.path.join("workspace", "tables")
KD = os.path.join("workspace", "research", "kickdoc")
EFF = {r["device_type_id"]: float(r["eff_length"]) for r in mt._read("hardware_library.csv")}


def per_theta(dev):
    """Наклон θ_дуга/correct_coeff (мрад/coeff) из арочного CSV устройства."""
    path = os.path.join(KD, f"dp_arc_{dev}.csv")
    if not os.path.exists(path):
        return None
    D = list(csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"))
    cc = np.array([float(r["correct_coeff"]) for r in D])
    th = np.array([float(r["theta_arc_mrad"]) for r in D])
    g = np.isfinite(cc) & np.isfinite(th)
    return float(np.polyfit(cc[g], th[g], 1)[0])


def main(lattice=101):
    path = os.path.join(TB, "lattice_configs.csv")
    rows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(rows[0].keys())
    lc = [r for r in rows if r["lattice_id"] == str(lattice)]

    # per_θ по каждому DP-типу, встречающемуся в кольце
    devs = sorted({r["device_type_id"] for r in lc if r["device_type_id"].startswith("DP_")})
    per = {}
    print(f"QP={QP:.3f}. Диполи DP в решётке {lattice}: {', '.join(devs)}")
    for d in devs:
        p = per_theta(d)
        per[d] = p
        print(f"  {d:16} Leff={EFF[d]:.3f}  per_θ={'н/д (нет арочного CSV)' if p is None else f'{p:.3f} мрад/coeff'}")

    n = 0; by = defaultdict(list)
    print(f"\n{'DP-тип':16}{'θ_проект,мрад':>14}{'correct_coeff':>15}{'ОДЗ':>6}")
    for r in lc:
        d = r["device_type_id"]
        if not d.startswith("DP_") or per.get(d) is None or abs(per[d]) < 1e-9:
            continue
        Leff = EFF[d]
        mB = float(r.get("mat_B") or 0)
        th_des = QP * mB * Leff * 1e3                        # проектный изгиб, мрад
        cc = th_des / per[d]
        r["correct_coeff"] = f"{cc:.6f}"
        r["arg_val"] = ""                                    # DP — PM, тока нет
        n += 1; by[d].append((th_des, cc))
        if len(by[d]) == 1:
            ok = "ДА" if abs(cc) <= 3 else "НЕТ!"
            print(f"{d:16}{th_des:>14.3f}{cc:>15.6f}{ok:>6}")

    shutil.copy(path, path + ".bak_dparc")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(rows)
    print(f"\n[OK] correct_coeff {n} диполей DP из арочного per_θ → lattice_configs {lattice} (бэкап .bak_dparc)")
    for d, v in by.items():
        v = np.array(v)
        print(f"    {d:16}: {len(v)} шт, θ∈[{v[:,0].min():+.2f},{v[:,0].max():+.2f}]мрад, "
              f"coeff∈[{v[:,1].min():+.3f},{v[:,1].max():+.3f}]")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 101)
