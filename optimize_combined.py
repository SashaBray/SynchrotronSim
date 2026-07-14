"""Combined optimization: чередование оконного (windowed) и глобального (global) режимов.

Идея:
    Оконный режим притёр каждый магнит к своей эталонной траектории (локальная
    геометрия), но фокусировку обеспечить не может -- она зависит от связки
    магнитов. Глобальный режим настраивает фокусировку (envelope) и компенсации,
    но из плохой начальной точки сходится медленно.

    Решение: чередовать фазы. Оконный даёт хороший геометрический старт, глобальный
    подкручивает фокусировку, оконный снова правит локальную геометрию, и т.д.

Конфиг в JSON:
    "combined": {
        "cycles": 5,
        "pe_passes_per_cycle": 2,
        "mc_iters_per_cycle":  1000,
        "mc_loss_mode": "composite"      # или "envelope" чтобы тюнить только огибающую
    }

Запуск:
    python optimize_combined.py <exp_id|path-to-config.json>
    opt_combo <exp_id>                   # через bat-обёртку

Между фазами sp_elements.csv обновляется на диске -- следующая фаза подхватывает
найденные параметры. Optimizer.exe загружается с нуля на каждой фазе
(~1-2 сек overhead -- незначительно по сравнению с временем фазы).
"""
from __future__ import annotations
import sys
import json
import shutil
import subprocess
from pathlib import Path

from _common import project_root, workspace_dir, optimizer_exe, set_matplotlib_backend


def find_json_by_exp_id(inputs_dir: Path, exp_id: int) -> list[Path]:
    matches: list[Path] = []
    for path in sorted(inputs_dir.glob("*.json")):
        try:
            cfg = json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            continue
        if isinstance(cfg, dict) and cfg.get("experiment_id") == exp_id:
            matches.append(path)
    return matches


def history_path_for(cfg_path: Path) -> Path:
    """Путь к optim_history.json, который пишет optimizer.exe рядом с cfg."""
    return cfg_path.with_suffix(".history.json")


def archive_history(cfg_path: Path, archive_dir: Path,
                    cycle: int, phase_letter: str, phase_method: str) -> Path | None:
    """Копирует history JSON из последней фазы в archive_dir с понятным именем.
    Возвращает путь архивной копии или None если history не существует."""
    src = history_path_for(cfg_path)
    if not src.exists():
        print(f"  WARN: история не найдена ({src.name}), фаза не сохранена в архив")
        return None
    archive_dir.mkdir(parents=True, exist_ok=True)
    dst = archive_dir / f"cycle{cycle:02d}_phase{phase_letter}_{phase_method}.json"
    shutil.copy(src, dst)
    return dst


def merge_histories(archive_files: list[Path]) -> dict:
    """Сливает истории фаз в одну, перенумеровав it_index непрерывно.
    Каждая запись получает поля cycle, phase, original_it -- для прозрачности."""
    merged_iters: list[dict] = []
    phase_meta: list[dict] = []
    groups_first = None
    it_offset = 0
    for ai, ap in enumerate(archive_files):
        try:
            with open(ap, encoding="utf-8") as f:
                h = json.load(f)
        except Exception as e:
            print(f"  WARN: не удалось прочитать {ap.name}: {e}")
            continue
        # Запоминаем groups из первой непустой истории
        if groups_first is None and "groups" in h:
            groups_first = h.get("groups")
        iters = h.get("iterations") or []
        # Извлекаем фазу из имени файла
        name_parts = ap.stem.split("_")
        cycle_str = name_parts[0] if name_parts else ""
        phase_str = name_parts[1] if len(name_parts) > 1 else ""
        method_str = "_".join(name_parts[2:]) if len(name_parts) > 2 else ""
        phase_start = it_offset
        for it_rec in iters:
            new_rec = dict(it_rec)
            new_rec["original_it"] = it_rec.get("it", 0)
            new_rec["it"] = it_offset
            new_rec["cycle"] = cycle_str
            new_rec["phase"] = phase_str
            new_rec["method"] = method_str
            merged_iters.append(new_rec)
            it_offset += 1
        phase_meta.append({
            "cycle":  cycle_str,
            "phase":  phase_str,
            "method": method_str,
            "it_start": phase_start,
            "it_end":   it_offset,
            "n_records": it_offset - phase_start,
        })
    return {
        "groups": groups_first or [],
        "phases": phase_meta,
        "iterations": merged_iters,
    }


def plot_combined_history(merged: dict, png_path: Path) -> None:
    """Рисует кривую loss vs iteration с цветными зонами PE и MC."""
    try:
        set_matplotlib_backend("Agg")
        import matplotlib.pyplot as plt
    except Exception as e:
        print(f"  WARN: matplotlib недоступен, график не построен: {e}")
        return
    iters = merged.get("iterations") or []
    phases = merged.get("phases") or []
    if not iters:
        print("  WARN: пусто merged.iterations, график не построен")
        return
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
    # Раскраска зон по методу фазы
    color_by_method = {"windowed": "#cfe2ff", "global": "#ffe1c2",
                       # legacy-метки старых архивов:
                       "per_element": "#cfe2ff", "mc": "#ffe1c2", "montecarlo": "#ffe1c2"}
    for ph in phases:
        c = color_by_method.get(ph["method"], "#eeeeee")
        for ax in (ax1, ax2):
            ax.axvspan(ph["it_start"], ph["it_end"], color=c, alpha=0.5, lw=0)
    its     = [r["it"]                   for r in iters]
    losses  = [r.get("loss", 0)          for r in iters]
    bests   = [r.get("best", r.get("bestLoss", r.get("loss", 0))) for r in iters]
    means   = [r.get("mean_mm", r.get("meanMM", 0)) for r in iters]
    maxs    = [r.get("max_mm",  r.get("maxMM",  0)) for r in iters]
    # 1) loss vs iter
    ax1.semilogy(its, losses, ".", ms=2, color="#333", label="iteration loss", alpha=0.4)
    ax1.semilogy(its, bests,  "-",  lw=1.3, color="#1f77b4", label="best so far")
    ax1.set_ylabel("Loss (log)")
    ax1.set_title("Combined history: windowed ⇄ global")
    ax1.legend(loc="best")
    ax1.grid(True, which="both", ls="--", alpha=0.3)
    # 2) mean/max mm
    ax2.plot(its, means, "-", color="#2ca02c", label="mean (mm)")
    ax2.plot(its, maxs,  "-", color="#d62728", label="max (mm)", alpha=0.6)
    ax2.set_xlabel("Iteration (across all phases)")
    ax2.set_ylabel("Deviation (mm)")
    ax2.legend(loc="best")
    ax2.grid(True, which="both", ls="--", alpha=0.3)
    # Аннотация цветовых зон
    legend_handles = [
        plt.Rectangle((0, 0), 1, 1, color="#cfe2ff", alpha=0.5, label="windowed"),
        plt.Rectangle((0, 0), 1, 1, color="#ffe1c2", alpha=0.5, label="global"),
    ]
    ax1.legend(handles=ax1.get_legend().legend_handles + legend_handles,
               loc="best", fontsize=9)
    fig.tight_layout()
    fig.savefig(png_path, dpi=130)
    plt.close(fig)
    print(f"  saved plot -> {png_path}")


def run_phase(cfg_path: Path, method: str, extra_flags: list[str] | None = None) -> None:
    """Регенерирует flat-конфиг и запускает opt.exe с указанным методом."""
    root = project_root()
    rc = subprocess.call([sys.executable, str(root / "prepare_optim_config.py"),
                          str(cfg_path)])
    if rc != 0:
        raise RuntimeError(f"prepare_optim_config.py failed (rc={rc})")
    exe = optimizer_exe(root)
    if not exe.is_file():
        raise RuntimeError(f"Optimizer exe not found: {exe}")
    if exe.stat().st_size == 0:
        raise RuntimeError(f"Optimizer exe is 0 bytes (Defender?): {exe}")
    optimconfig = cfg_path.with_suffix(".optimconfig")
    cmd = [str(exe), str(optimconfig), f"--method={method}"]
    if extra_flags:
        cmd.extend(extra_flags)
    print(f"\n[combined] >>> {' '.join(cmd)}\n", flush=True)
    rc = subprocess.call(cmd)
    if rc != 0:
        raise RuntimeError(f"opt.exe ({method}) failed (rc={rc})")


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python optimize_combined.py <exp_id|path-to-config.json>")
        print()
        print(__doc__.strip())
        sys.exit(1)

    root = project_root()
    inputs_dir = workspace_dir(root) / "optim_inputs"
    arg1 = sys.argv[1]
    arg_path = Path(arg1)
    if arg_path.is_file():
        target = arg_path
    else:
        try:
            exp_id = int(arg1)
        except ValueError:
            print(f"[combined] '{arg1}' не файл и не целое число.")
            sys.exit(1)
        matches = find_json_by_exp_id(inputs_dir, exp_id)
        if not matches:
            print(f"[combined] нет JSON с experiment_id={exp_id} в {inputs_dir}")
            sys.exit(1)
        target = matches[0]
        if len(matches) > 1:
            print(f"[combined] предупреждение: {len(matches)} JSON с этим exp_id, "
                  f"использую первый: {target.name}")

    cfg = json.loads(target.read_text(encoding="utf-8"))
    combined = cfg.get("combined", {})
    cycles    = int(combined.get("cycles", 3))
    pe_passes = int(combined.get("pe_passes_per_cycle", 2))
    mc_iters  = int(combined.get("mc_iters_per_cycle", 1000))
    mc_loss   = str(combined.get("mc_loss_mode", "composite"))

    print(f"[combined] target: {target.name}")
    print(f"[combined] cycles={cycles}, pe_passes_per_cycle={pe_passes}, "
          f"mc_iters_per_cycle={mc_iters}, mc_loss_mode={mc_loss}")

    # Временный JSON, где мы мутируем passes / n_iterations для каждой фазы.
    # Имя:  <stem>.combined_tmp.json -- чтобы не пересекаться с базовым файлом.
    tmp_path = target.with_name(target.stem + ".combined_tmp.json")

    # Архив историй каждой фазы (для сборки единой кривой в конце).
    exp_id = cfg.get("experiment_id", "X")
    archive_dir = workspace_dir(root) / "results" / f"combined_history_exp{exp_id}"
    if archive_dir.exists():
        # Чистим предыдущие истории, чтобы случайные старые JSON-файлы не попали в график
        for old in archive_dir.glob("cycle*_phase*.json"):
            try: old.unlink()
            except OSError: pass
    archived_files: list[Path] = []

    try:
        for cycle in range(1, cycles + 1):
            print(f"\n{'='*60}\n  CYCLE {cycle}/{cycles}\n{'='*60}")

            # ---- Phase A: оконный (windowed) ----
            print(f"\n[cycle {cycle}] PHASE A: windowed ({pe_passes} passes)")
            cfg_pe = json.loads(target.read_text(encoding="utf-8"))
            cfg_pe.setdefault("per_element", {})["passes"] = pe_passes
            tmp_path.write_text(json.dumps(cfg_pe, ensure_ascii=False, indent=2),
                                encoding="utf-8")
            run_phase(tmp_path, "windowed")
            arc = archive_history(tmp_path, archive_dir, cycle, "A", "windowed")
            if arc: archived_files.append(arc)

            # ---- Phase B: глобальный (global) ----
            # После оконной фазы sp_elements.csv обновлён -- следующая фаза подхватывает.
            print(f"\n[cycle {cycle}] PHASE B: global ({mc_iters} iters, "
                  f"loss={mc_loss})")
            cfg_mc = json.loads(target.read_text(encoding="utf-8"))
            cfg_mc["n_iterations"] = mc_iters
            tmp_path.write_text(json.dumps(cfg_mc, ensure_ascii=False, indent=2),
                                encoding="utf-8")
            run_phase(tmp_path, "global", extra_flags=[f"--loss={mc_loss}"])
            arc = archive_history(tmp_path, archive_dir, cycle, "B", "global")
            if arc: archived_files.append(arc)

        print(f"\n{'='*60}\n[combined] DONE: {cycles} cycles complete\n{'='*60}")

        # ---- Сборка единой истории и графика ----
        if archived_files:
            print(f"\n[combined] merging {len(archived_files)} histories ...")
            merged = merge_histories(archived_files)
            merged_path = archive_dir / "combined_history.json"
            merged_path.write_text(
                json.dumps(merged, ensure_ascii=False, indent=1), encoding="utf-8")
            print(f"  saved merged history -> {merged_path}")
            png_path = archive_dir / "combined_history.png"
            plot_combined_history(merged, png_path)
        else:
            print(f"[combined] no history archived, skipping merge/plot")
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass
        # Удалим также его .optimconfig если остался
        tmp_oc = tmp_path.with_suffix(".optimconfig")
        if tmp_oc.exists():
            try:
                tmp_oc.unlink()
            except OSError:
                pass


if __name__ == "__main__":
    main()
