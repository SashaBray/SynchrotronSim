"""
Визуализатор хода MC-оптимизации:
  1. Loss vs iteration (current/best, log-y), кики отмечены вертикальными линиями.
  2. Эволюция самых активных параметров — как они менялись по итерациям.
  3. Гистограмма размера принятых шагов |new - old|.
  4. ТОП-N групп по активности (число ACCEPT'ов).
  5. Сводка статистики прогона.

В отличие от viz_optim_history.py (тот показывает финальные значения параметров
как bar chart), этот скрипт показывает ДИНАМИКУ параметров в течение оптимизации.

Запуск:
    python viz_mc_evolution.py                          # последний history.json
    python viz_mc_evolution.py path/to/history.json     # явный путь
    python viz_mc_evolution.py --top 20                 # сколько параметров рисовать
    python viz_mc_evolution.py --no-show                # сохранить PNG, окно не показывать
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from _common import set_matplotlib_backend, workspace_dir
set_matplotlib_backend("TkAgg")          # fallback to Qt5Agg / Agg
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                        # noqa: E402


# --------------------------------------------------------------------------
# Поиск файла истории
# --------------------------------------------------------------------------
def find_latest_history() -> Path | None:
    """Самый свежий *.history.json в workspace/results/."""
    rdir = workspace_dir() / "results"
    if not rdir.exists():
        return None
    candidates = sorted(rdir.glob("*.history.json"),
                        key=lambda p: p.stat().st_mtime,
                        reverse=True)
    return candidates[0] if candidates else None


# --------------------------------------------------------------------------
# Реконструкция эволюции параметров
# --------------------------------------------------------------------------
def reconstruct_param_traces(iters, group_names, best_values):
    """Для каждой группы name собираем (iters_when_changed, values_after).

    Returns:
        traces: dict {name: (np.array iters, np.array values)}
        activity: dict {name: число ACCEPT'ов (включая ACCEPT+BEST)}
    """
    # Для каждого параметра отслеживаем его текущее значение.
    # Стартового значения у нас нет — берём как best_values (это значения после
    # последней BEST; для тех групп, которые ни разу не меняли best, это и есть
    # стартовое значение). Не идеально, но даёт точку отсчёта для первой
    # итерации.
    current = {name: best_values[i] for i, name in enumerate(group_names)}
    # Уточняем: если есть ACCEPT-итерации по группе, первое 'old' — стартовое.
    first_old = {}
    for e in iters:
        name = e.get("name", "")
        if name in current and name not in first_old and "ACCEPT" in e.get("mark", ""):
            first_old[name] = e["old"]
    for name, v in first_old.items():
        current[name] = v

    traces = {name: ([0], [current[name]]) for name in group_names}
    activity = {name: 0 for name in group_names}

    for e in iters:
        mark = e.get("mark", "")
        name = e.get("name", "")
        if "ACCEPT" in mark and name in current:
            current[name] = e["new"]
            traces[name][0].append(e["it"])
            traces[name][1].append(e["new"])
            activity[name] += 1
        # KICK — name == "KICK"; индивидуальные new по группам не записаны,
        # поэтому конкретные значения после кика нам неизвестны. Эволюционные
        # линии параметров между ACCEPT'ами после кика будут иметь разрыв —
        # это сознательно (мы видим, как параметр продолжает меняться, но без
        # «прыжка» от кика).

    # Преобразуем в np.array
    traces = {name: (np.array(it), np.array(vs))
              for name, (it, vs) in traces.items()}
    return traces, activity


# --------------------------------------------------------------------------
# Главная функция
# --------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("path", nargs="?", help="path to *.history.json (default: latest in workspace/results/)")
    p.add_argument("--top", type=int, default=12,
                   help="number of most-active parameters to plot (default 12)")
    p.add_argument("--no-show", action="store_true",
                   help="save PNG only, do not pop up a window")
    args = p.parse_args()

    if args.path:
        path = Path(args.path)
    else:
        path = find_latest_history()
        if path is None:
            print("ERROR: no *.history.json found in workspace/results/")
            sys.exit(1)
        print(f"Latest history file: {path}")

    if not path.exists():
        print(f"ERROR: {path} not found")
        sys.exit(1)

    with open(path, encoding="utf-8") as f:
        h = json.load(f)

    iters_data  = h["iterations"]
    if not iters_data:
        print("ERROR: empty iterations[] in history")
        sys.exit(1)

    it          = np.array([e["it"]       for e in iters_data])
    loss        = np.array([e["loss"]     for e in iters_data])
    best        = np.array([e["best"]     for e in iters_data])
    mean_mm     = np.array([e.get("mean_mm", 0)  for e in iters_data])
    marks       = [e["mark"]              for e in iters_data]

    accept_mask = np.array(["ACCEPT" in m and "BEST" not in m for m in marks])
    best_mask   = np.array(["BEST"   in m for m in marks])
    kick_mask   = np.array(["KICK"   in m for m in marks])
    fail_mask   = np.array([m in ("KICK_FAIL", "FAIL") for m in marks])

    group_names = h["group_names"]
    best_vals   = h["best_values"]

    traces, activity = reconstruct_param_traces(iters_data, group_names, best_vals)

    # Размер принятых шагов (для гистограммы)
    accept_step_sizes = []
    for e in iters_data:
        mark = e.get("mark", "")
        if "ACCEPT" in mark and e.get("name") != "KICK":
            try:
                accept_step_sizes.append(abs(e["new"] - e["old"]))
            except (KeyError, TypeError):
                pass
    accept_step_sizes = np.array(accept_step_sizes)

    # Топ-N самых активных параметров
    top_names = sorted(group_names, key=lambda n: -activity[n])[:args.top]

    # =====================================================================
    # FIGURE
    # =====================================================================
    fig = plt.figure(figsize=(16, 10))
    gs  = fig.add_gridspec(3, 3, height_ratios=[1.4, 1.0, 0.9],
                           hspace=0.45, wspace=0.30)

    expid   = h.get("experiment_id", "?")
    n_iter  = h.get("n_iterations", len(iters_data))
    loop_s  = h.get("loop_sec", 0)
    acc     = h.get("accepted", int(accept_mask.sum() + best_mask.sum()))
    rej     = h.get("rejected", 0)
    kicks   = h.get("kicks", int(kick_mask.sum()))
    failed  = h.get("failed", 0)
    best_L  = h.get("best_loss", float(best.min()))

    fig.suptitle(
        f"MC-оптимизация — exp {expid}  •  {n_iter} итер  •  {loop_s:.0f} с  •  "
        f"accept={acc}  reject={rej}  kicks={kicks}  failed={failed}  •  "
        f"best loss = {best_L:.4g}",
        fontsize=11
    )

    # ---------------- 1. LOSS vs ITER (главный) ----------------
    ax = fig.add_subplot(gs[0, :])
    # Кики — вертикальные оранжевые полосы под графиком
    if kick_mask.any():
        for k_it in it[kick_mask]:
            ax.axvline(k_it, color="orange", alpha=0.20, linewidth=0.8, zorder=0)
    ax.semilogy(it, loss, ".", color="lightgray", markersize=2, label="loss за итерацию")
    ax.semilogy(it, best, "-", color="C0", linewidth=1.6, label="best loss")
    if best_mask.any():
        ax.semilogy(it[best_mask], loss[best_mask], "o", color="C2",
                    markersize=4, markeredgewidth=0, label=f"новые BEST ({best_mask.sum()})")
    if fail_mask.any():
        ax.semilogy(it[fail_mask], loss[fail_mask], "x", color="C3",
                    markersize=6, label=f"FAIL ({fail_mask.sum()})")
    ax.set_xlabel("итерация")
    ax.set_ylabel("loss (log)")
    ax.set_title(f"Loss vs iteration  (оранжевые полосы — кики)")
    ax.legend(loc="best", fontsize=9)
    ax.grid(True, alpha=0.3)

    # ---------------- 2. ЭВОЛЮЦИЯ ПАРАМЕТРОВ ----------------
    ax = fig.add_subplot(gs[1, :2])
    cmap = plt.get_cmap("tab20")
    for i, name in enumerate(top_names):
        x, y = traces[name]
        if len(x) < 2:
            continue
        ax.plot(x, y, "-", linewidth=1.0, color=cmap(i % 20),
                label=f"{name} ({activity[name]} accept)")
    ax.set_xlabel("итерация")
    ax.set_ylabel("значение параметра")
    ax.set_title(f"ТОП-{args.top} самых активных параметров (по числу ACCEPT)")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, ncol=1, frameon=False)
    ax.grid(True, alpha=0.3)

    # ---------------- 3. ГИСТОГРАММА РАЗМЕРОВ ACCEPT-ШАГОВ ----------------
    ax = fig.add_subplot(gs[1, 2])
    if len(accept_step_sizes) > 0:
        ax.hist(accept_step_sizes, bins=40, color="C2", alpha=0.75)
        ax.set_xlabel("|new − old|")
        ax.set_ylabel("count")
        ax.set_title(f"Распределение размера принятых шагов\n"
                     f"(n={len(accept_step_sizes)}, "
                     f"медиана={np.median(accept_step_sizes):.3g})")
        ax.grid(True, alpha=0.3)
    else:
        ax.text(0.5, 0.5, "нет ACCEPT-шагов", ha="center", va="center")

    # ---------------- 4. АКТИВНОСТЬ ПО ГРУППАМ ----------------
    ax = fig.add_subplot(gs[2, :2])
    sorted_act = sorted(activity.items(), key=lambda kv: -kv[1])
    n_show = min(args.top * 2, len(sorted_act))
    show = sorted_act[:n_show]
    if show:
        names_show = [n for n, _ in show]
        counts     = [c for _, c in show]
        y_pos = np.arange(len(show))
        ax.barh(y_pos, counts, color="C0", height=0.7)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(names_show, fontsize=6)
        ax.invert_yaxis()
        ax.set_xlabel("число ACCEPT'ов")
        ax.set_title(f"Активность групп — {n_show} наиболее активных")
        ax.grid(True, alpha=0.3, axis="x")

    # ---------------- 5. СВОДКА ----------------
    ax = fig.add_subplot(gs[2, 2])
    ax.axis("off")
    total_iters = len(iters_data)
    acc_rate = 100.0 * (acc / total_iters) if total_iters else 0
    median_step = (np.median(accept_step_sizes) if len(accept_step_sizes) > 0 else 0)
    loss_drop  = float(loss[0] - best[-1])
    loss_ratio = float(loss[0] / best[-1]) if best[-1] > 0 else 0
    info_lines = [
        f"Файл:  {path.name}",
        f"",
        f"Эксперимент:     {expid}",
        f"Итераций:        {total_iters}",
        f"Время:           {loop_s:.0f} с ({loop_s/60:.1f} мин)",
        f"",
        f"ACCEPT:          {acc}  ({acc_rate:.1f} %)",
        f"REJECT:          {rej}",
        f"KICK:            {kicks}",
        f"FAILED:          {failed}",
        f"",
        f"Loss[0]:         {loss[0]:.4g}",
        f"Best loss:       {best_L:.4g}",
        f"Падение:         {loss_drop:.3g}",
        f"Отношение:       {loss_ratio:.2f}×",
        f"",
        f"Параметров:      {len(group_names)}",
        f"Из них меняли:   {sum(1 for v in activity.values() if v > 0)}",
        f"Медиан. шаг:     {median_step:.4g}",
    ]
    ax.text(0.0, 1.0, "\n".join(info_lines),
            family="monospace", fontsize=9,
            verticalalignment="top")

    # =====================================================================
    out_png = path.with_suffix(".evolution.png")
    plt.savefig(out_png, dpi=120, bbox_inches="tight")
    print(f"Saved {out_png}")
    if not args.no_show:
        plt.show()


if __name__ == "__main__":
    main()
