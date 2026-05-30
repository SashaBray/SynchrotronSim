"""
Визуализатор истории оптимизации, записанной SynchrotronOptimizer.exe.

Запуск:
    python viz_optim_history.py [path/to/optim_history.json]

По умолчанию ищет workspace/results/sp1_endpoint.history.json.

Рисует:
  - loss vs iteration (текущий + best, log-y)
  - mean / max distance к магнитам vs iteration
  - распределение принятых/отклонённых
  - финальные значения параметров (bar chart)
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    default_path = Path(__file__).resolve().parent / "workspace" / "results" / "sp1_endpoint.history.json"
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path
    if not path.exists():
        print(f"ERROR: history file not found: {path}")
        sys.exit(1)
    with open(path, encoding="utf-8") as f:
        h = json.load(f)

    iters_data  = h["iterations"]
    it          = np.array([e["it"]       for e in iters_data])
    loss        = np.array([e["loss"]     for e in iters_data])
    best        = np.array([e["best"]     for e in iters_data])
    mean_mm     = np.array([e["mean_mm"]  for e in iters_data])
    max_mm      = np.array([e["max_mm"]   for e in iters_data])
    marks       = [e["mark"]              for e in iters_data]

    accept_mask = np.array(["ACCEPT"  in m for m in marks])
    reject_mask = np.array(["reject"  in m for m in marks])
    kick_mask   = np.array(["KICK"    in m for m in marks])
    best_mask   = np.array(["BEST"    in m for m in marks])

    # --- Figure ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"Optimization history -- exp {h['experiment_id']}  "
                 f"({h['n_iterations']} iter, {h['loop_sec']:.1f}s, "
                 f"{h['accepted']} acc / {h['rejected']} rej / "
                 f"{h['kicks']} kicks / {h['failed']} fail)",
                 fontsize=11)

    # Loss
    ax = axes[0][0]
    ax.semilogy(it, loss, ".", color="lightgray", markersize=2, label="iter loss")
    ax.semilogy(it, best, "-", color="C0", linewidth=1.5, label="best")
    if kick_mask.any():
        ax.semilogy(it[kick_mask], loss[kick_mask], "x", color="orange",
                    markersize=6, label="kicks")
    if best_mask.any():
        ax.semilogy(it[best_mask], loss[best_mask], "o", color="C2",
                    markersize=4, label="new best")
    ax.set_xlabel("iteration"); ax.set_ylabel("loss")
    ax.set_title(f"Loss  (final best = {h['best_loss']:.6g})")
    ax.legend(loc="best"); ax.grid(True, alpha=0.3)

    # Mean/max distance
    ax = axes[0][1]
    ax.plot(it, mean_mm, ".", color="C0", markersize=2, label="mean (mm)")
    ax.plot(it, max_mm,  ".", color="C3", markersize=2, label="max (mm)")
    ax.set_xlabel("iteration"); ax.set_ylabel("distance, mm")
    ax.set_title("Particle->magnet distance")
    ax.set_yscale("log")
    ax.legend(loc="best"); ax.grid(True, alpha=0.3)

    # Accept rate over time (rolling 50-iter)
    ax = axes[1][0]
    window = max(20, len(it) // 50)
    if len(it) > window:
        acc = accept_mask.astype(float)
        roll = np.convolve(acc, np.ones(window) / window, mode="valid")
        ax.plot(it[window-1:], roll * 100, color="C1")
        ax.set_xlabel("iteration"); ax.set_ylabel(f"% accepted ({window}-iter rolling)")
        ax.set_title("Acceptance rate")
        ax.set_ylim(0, 100)
    else:
        ax.text(0.5, 0.5, "not enough iterations for rolling stats",
                ha="center", va="center")
    ax.grid(True, alpha=0.3)

    # Final values bar chart
    ax = axes[1][1]
    names  = h["group_names"]
    values = h["best_values"]
    # Сортируем по |value| убыванию для читаемости
    order = sorted(range(len(values)), key=lambda i: -abs(values[i]))
    sorted_names  = [names[i]  for i in order]
    sorted_values = [values[i] for i in order]
    y = np.arange(len(values))
    colors = ["C2" if v > 0 else "C3" for v in sorted_values]
    ax.barh(y, sorted_values, color=colors, height=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(sorted_names, fontsize=6)
    ax.set_xlabel("value")
    ax.set_title(f"Final best parameter values ({len(values)} groups)")
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()
    out_png = path.with_suffix(".png")
    plt.savefig(out_png, dpi=120)
    print(f"Saved {out_png}")
    plt.show()


if __name__ == "__main__":
    main()
