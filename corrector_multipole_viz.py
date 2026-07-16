# -*- coding: utf-8 -*-
"""Визуализация мультипольного состава корректора (5 обмоток × 8 коэф.) из датасетов
build_coupling_datasets (workspace/research/coupling_datasets/corrector_*.csv).

Два вида:
  1) corrector_<w>_vs_current.png — 8 панелей (B,Bs,G,Gs,S,Ss,O,Os) vs ток обмотки: линейность+знак.
  2) corrector_kick_spectrum.png  — «кик-спектр»: вклад каждого мультиполя в отклонение (мкрад)
     при r=4 мм, 100 А, лог-шкала, сгруппировано по обмоткам → видно, какие мультиполи
     динамически значимы у каждой обмотки (нормальные vs скью).
Запуск: py -3.14 corrector_multipole_viz.py"""
import os, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bip_matrix import QP

D = os.path.join("workspace", "research", "coupling_datasets")
WIND = ["main", "c1", "c2", "c3", "c4"]
COEF = ["B", "Bs", "G", "Gs", "S", "Ss", "O", "Os"]
UNIT = {"B": "Тл", "Bs": "Тл", "G": "Тл/м", "Gs": "Тл/м", "S": "Тл/м²", "Ss": "Тл/м²", "O": "Тл/м³", "Os": "Тл/м³"}
LEFF = {r["device_type_id"]: float(r["eff_length"])
        for r in csv.DictReader(open(os.path.join("workspace", "tables", "hardware_library.csv"), encoding="utf-8"), delimiter=";")}
RREF = 4e-3                                              # опорный радиус для кик-спектра


def load(w):
    p = os.path.join(D, f"corrector_{w}.csv")
    rows = list(csv.DictReader(open(p, encoding="utf-8-sig"), delimiter=";"))
    I = np.array([float(r["I_A"]) for r in rows])
    return I, {c: np.array([float(r[f"{c}_phys"]) for r in rows]) for c in COEF}


def kick_contrib(coef, val, Leff):
    """Интегральный вклад мультиполя в отклонение (рад) при r=RREF: dipole ∝1, quad ∝r, sext ∝r²/2, oct ∝r³/6."""
    base = abs(QP) * Leff * abs(val)                    # QP<0 → берём модуль (нужна величина вклада)
    return {"B": base, "Bs": base, "G": base * RREF, "Gs": base * RREF,
            "S": base * RREF**2 / 2, "Ss": base * RREF**2 / 2,
            "O": base * RREF**3 / 6, "Os": base * RREF**3 / 6}[coef]


def per_winding_panels():
    for w in WIND:
        I, dat = load(w)
        fig, axes = plt.subplots(2, 4, figsize=(15, 7))
        for ax, c in zip(axes.ravel(), COEF):
            skew = c.endswith("s") and c != "Bs" or c == "Bs"
            ax.plot(I, dat[c], "o-", color="#c0392b" if c in ("Bs", "Gs", "Ss", "Os") else "#2c3e50", lw=1.8, ms=5)
            ax.axhline(0, color="gray", lw=0.6); ax.axvline(0, color="gray", lw=0.6)
            ax.set_title(f"{c}  ({UNIT[c]})", fontsize=11,
                         color="#c0392b" if c in ("Bs", "Gs", "Ss", "Os") else "#2c3e50")
            ax.set_xlabel("ток обмотки, А"); ax.grid(alpha=0.3)
        fig.suptitle(f"Корректор, обмотка {w}: мультипольный состав vs ток "
                     f"(красное — скью). Чёрное — нормальные (B,G,S,O)", fontsize=13)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        out = os.path.join(D, f"corrector_{w}_vs_current.png")
        fig.savefig(out, dpi=110); plt.close(fig)
        print(f"[png] {out}")


def kick_spectrum():
    """Сгруппированный бар: вклад каждого мультиполя (мкрад) при 100 А, r=4 мм, по обмоткам, лог-шкала."""
    fig, ax = plt.subplots(figsize=(15, 6))
    x = np.arange(len(COEF)); bw = 0.16
    colors = {"main": "#2c3e50", "c1": "#2980b9", "c2": "#27ae60", "c3": "#e67e22", "c4": "#8e44ab"}
    for k, w in enumerate(WIND):
        dev = f"corrector_{w}_100A_L100"; Leff = LEFF[dev]
        I, dat = load(w)
        i100 = int(np.argmin(np.abs(I - 100)))
        vals = [kick_contrib(c, dat[c][i100], Leff) * 1e6 for c in COEF]     # мкрад
        vals = [max(v, 1e-4) for v in vals]                                  # лог-пол
        ax.bar(x + (k - 2) * bw, vals, bw, label=w, color=colors[w])
    ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(COEF)
    ax.set_ylabel("вклад в отклонение при 100 А, r=4 мм, мкрад")
    ax.set_title("Кик-спектр корректора: динамическая значимость мультиполей по обмоткам\n"
                 "(main — почти чистый гор.диполь B; c1–c4 несут полный набор нормальных и скью мультиполей)",
                 fontsize=12)
    ax.legend(title="обмотка", ncol=5); ax.grid(alpha=0.3, axis="y")
    for i, c in enumerate(COEF):                          # подсветить скью
        if c in ("Bs", "Gs", "Ss", "Os"):
            ax.get_xticklabels()[i].set_color("#c0392b")
    fig.tight_layout()
    out = os.path.join(D, "corrector_kick_spectrum.png")
    fig.savefig(out, dpi=120); plt.close(fig)
    print(f"[png] {out}")


if __name__ == "__main__":
    per_winding_panels()
    kick_spectrum()
    print("[OK] визуализации мультипольного состава корректора → " + D)
