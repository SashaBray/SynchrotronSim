# -*- coding: utf-8 -*-
"""Арочный тест ДИПОЛЯ DP: correct_coeff → (θ, G, S) СИММЕТРИЧНЫМ проездом (метод как для DQ) на ТЕКУЩИХ
картах. DP — постоянный магнит (input_type=coeff): поле = карта×correct_coeff×field_direction. Свип coeff,
для каждого — симметричная дуга (влёт −θ/2=вылет +θ/2) + матчинг 13-частичного пучка Ø2мм. Даёт per_θ
(наклон θ/coeff) на ТЕКУЩИХ картах — для инверсии (проектный θ → correct_coeff) токов диполей кольца, и
сравнения с прежним dp_bend_dataset (была ли рассинхронизация после переоснастки карт).
Переиспользует ядро dq_arc_dataset (переопределяет геометрию под карту DP).
Запуск: py -3.14 dp_arc_dataset.py [--dev DP_L_PM|DP_sh_PM] [--cmax 0.005] [--n 21]"""
import os, csv, argparse
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import dq_arc_dataset as da
import dq_sector_profile as dsp

HW, EFF, QP = da.HW, da.EFF, da.QP
KD = os.path.join("workspace", "research", "kickdoc")


def setup_dp(dev):
    """Загрузить карту DP и ПЕРЕОПРЕДЕЛИТЬ геометрию dq_arc_dataset под неё (half, Leff, Lphys, lead, RMAX)."""
    h, B, nom = dsp._load_full(dev, HW[dev])
    fd = float(HW[dev].get("field_direction") or 1)
    hwcc = float(HW[dev].get("correct_coeff") or 1)
    f = dsp._field_fn(h, B); half = float(h["mxx"])
    da.HALF = half; da.LEFF = EFF[dev]; da.LPHYS = da.LEFF          # phys убран → длина по L_eff
    da.LEAD = max(0.0, (2 * half - da.LPHYS) / 2.0); da.DS_G = da.LPHYS / da.N; da.DS_F = da.LEFF / da.N
    da.RMAX = min(min(abs(h["mny"]), abs(h["mny"] + (h["ny"] - 1) * h["sy"])),
                  min(abs(h["mnz"]), abs(h["mnz"] + (h["nz"] - 1) * h["sz"])))

    def field(cc):
        s = cc * hwcc * fd                                   # поле = карта × correct_coeff × field_direction
        return lambda p: s * f(p)
    return field, fd


def main(dev="DP_L_PM", cmax=0.005, npts=21, out=None):
    field, fd = setup_dp(dev)
    base = out or f"dp_arc_{dev}"
    print(f"Арочный тест {dev}: L_eff={da.LEFF:.4f} L_phys={da.LPHYS:.3f} полузона={da.HALF:.3f} "
          f"field_dir={fd:g}, RMAX={da.RMAX*1e3:.0f}мм. Свип correct_coeff ±{cmax} ({npts} т.), пучок Ø2мм.")
    grid = np.linspace(-cmax, cmax, npts)
    rows = []
    print(f"\n{'coeff':>11}{'θ_ось,мрад':>12}{'θ_дуга,мрад':>13}{'G_дуга,Тл/м':>13}{'нев,мкм':>9}")
    for cc in grid:
        cc = round(float(cc), 8)
        if abs(cc) < 1e-9:
            rows.append([f"{cc:.8f}", "0", "0", "0", "0", "0", "0"]); continue
        fld = field(cc)
        try:
            th_ax, G_ax, S_ax = da.axis_coeffs(fld)
            phi, th_sym, _ = da.symmetrize(fld)
            g0 = G_ax if np.isfinite(G_ax) else 0.0
            (th_m, G_m, S_m, O_m), rms, sag = da.fit_matrix(fld, phi, [th_sym, g0, 0.0, 0.0])
        except Exception:
            continue
        rows.append([f"{cc:.8f}", f"{th_ax*1e3:.5f}", f"{th_m*1e3:.5f}", f"{G_m:.6f}", f"{S_m:.5f}",
                     f"{sag*1e3:.4f}", f"{rms:.5f}"])
        if abs(cc) > cmax * 0.4:                             # печать по краям
            print(f"{cc:>11.5f}{th_ax*1e3:>12.4f}{th_m*1e3:>13.4f}{G_m:>13.5f}{rms*1e3:>9.1f}")

    cc = np.array([float(r[0]) for r in rows]); tha = np.array([float(r[2]) for r in rows])
    thax = np.array([float(r[1]) for r in rows]); g = np.isfinite(tha) & np.isfinite(thax)
    per_th = float(np.polyfit(cc[g], tha[g], 1)[0])          # мрад/coeff (дуга)
    per_th_ax = float(np.polyfit(cc[g], thax[g], 1)[0])      # мрад/coeff (ось)
    print(f"\nper_θ (наклон θ/correct_coeff): ДУГА {per_th:.1f} мрад/coeff, ОСЬ {per_th_ax:.1f} мрад/coeff")
    if dev == "DP_L_PM":
        print(f"  прежний dp_bend_dataset (осевой, старые карты): DP_L per_θ=5980 мрад/coeff — "
              f"{'СОВПАДАЕТ' if abs(per_th_ax-5980)/5980<0.1 else 'РАССОГЛАСОВАНО ⇒ калибровка устарела'}")
    print(f"  инверсия: для проектного θ correct_coeff = θ[мрад]/{per_th:.1f}")

    os.makedirs(KD, exist_ok=True)
    with open(os.path.join(KD, base + ".csv"), "w", newline="", encoding="utf-8-sig") as fcsv:
        w = csv.writer(fcsv, delimiter=";")
        w.writerow(["correct_coeff", "theta_axis_mrad", "theta_arc_mrad", "G_arc_Tlm", "S_arc_Tlm2",
                    "sagitta_mm", "fit_rms_mm"])
        w.writerows(rows)

    os.makedirs(os.path.join(KD, "figs"), exist_ok=True)
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    ax[0].plot(cc, thax, "o-", color="#e15759", label="ОСЬ (прямой)")
    ax[0].plot(cc, tha, "s-", color="#2f9e44", label="ДУГА (симметр.)")
    ax[0].set_xlabel("correct_coeff"); ax[0].set_ylabel("изгиб θ, мрад")
    ax[0].set_title(f"{dev}: изгиб (per_θ_дуга={per_th:.0f})"); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    Ga = np.array([float(r[3]) for r in rows])
    ax[1].plot(cc, Ga, "s-", color="#2f9e44"); ax[1].set_xlabel("correct_coeff"); ax[1].set_ylabel("G, Тл/м")
    ax[1].set_title("паразитный градиент (дуга)"); ax[1].grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(KD, "figs", base + ".png"), dpi=130); plt.close(fig)
    print(f"\n[OK] {len(rows)} строк → {os.path.join(KD, base + '.csv')}; фигура → figs/{base}.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="арочный тест диполя DP: correct_coeff → θ,G,S")
    ap.add_argument("--dev", default="DP_L_PM", help="устройство DP (DP_L_PM / DP_sh_PM)")
    ap.add_argument("--cmax", type=float, default=0.005, help="макс |correct_coeff| свипа")
    ap.add_argument("--n", type=int, default=21, help="точек свипа (нечётное — включает 0)")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    main(a.dev, a.cmax, a.n, a.out)
