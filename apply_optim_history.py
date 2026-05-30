"""
Применить best-значения из истории оптимизации (sp1_endpoint.history.json) к БД:
  - workspace/tables/sp_elements.csv      (для групп с target=CURRENT, value_col=arg_val)
  - workspace/tables/hardware_library.csv (для групп с target=COEFF,    value_col=correct_coeff)

После записи -- автоматически пересобирает lattice_configs.csv через
assemble_lattice_configs.py, чтобы трекер сразу видел итог.

Запуск:
    python apply_optim_history.py [path/to/optim_history.json]

По умолчанию: workspace/results/sp1_endpoint.history.json.

Историю генерирует C++ оптимизатор:
  - после каждых 100 итераций (чекпойнт) -- "final": false
  - в конце прогона                       -- "final": true
Скрипт работает с любым из этих состояний.
"""
import csv
import json
import subprocess
import sys
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR / "workspace"
SP_PATH    = WS / "tables" / "sp_elements.csv"
HW_PATH    = WS / "tables" / "hardware_library.csv"
ASSEMBLE   = SCRIPT_DIR / "assemble_lattice_configs.py"


def read_csv(path):
    with open(path, encoding="utf-8") as f:
        reader = csv.reader(f, delimiter=";")
        rows = list(reader)
    return rows[0], rows[1:]


def write_csv(path, header, rows, attempts=12, delay=0.5):
    """Запись с retry: Windows иногда держит CSV открытым из Excel/антивируса.
    Несколько повторов с паузой решают проблему."""
    import time
    last_err = None
    for k in range(attempts):
        try:
            with open(path, "w", encoding="utf-8", newline="") as f:
                w = csv.writer(f, delimiter=";", lineterminator="\n")
                w.writerow(header)
                w.writerows(rows)
            return
        except PermissionError as e:
            last_err = e
            if k == 0:
                print(f"  [retry] {path.name if hasattr(path, 'name') else path} locked, "
                      f"waiting (close Excel/viewer if open)...")
            time.sleep(delay)
    raise last_err


def row_matches(row, header, matchers):
    """matchers: список {"column": <name>, "values": [<val>, ...]} (OR в values)."""
    for m in matchers:
        col = m["column"]
        if col not in header:
            return False
        actual = row[header.index(col)].strip()
        if actual not in m["values"]:
            return False
    return True


def apply_group(group, sp_rows, sp_header, hw_rows, hw_header):
    """Возвращает (sp_modified_count, hw_modified_count)."""
    target  = group["target"]
    fname   = group["file"]
    val_col = group["value_col"]
    matchers = group["matchers"]
    best    = group["best"]

    sp_mod = hw_mod = 0
    if target == "CURRENT" or (fname == "sp_elements.csv" and val_col == "arg_val"):
        if val_col not in sp_header:
            return 0, 0
        val_idx = sp_header.index(val_col)
        for row in sp_rows:
            if row_matches(row, sp_header, matchers):
                new_str = f"{best:.10f}"
                if row[val_idx] != new_str:
                    row[val_idx] = new_str
                    sp_mod += 1
    elif target == "COEFF" or (fname == "hardware_library.csv" and val_col == "correct_coeff"):
        if val_col not in hw_header:
            return 0, 0
        val_idx = hw_header.index(val_col)
        for row in hw_rows:
            if row_matches(row, hw_header, matchers):
                new_str = f"{best:.10f}"
                if row[val_idx] != new_str:
                    row[val_idx] = new_str
                    hw_mod += 1
    return sp_mod, hw_mod


def main():
    default_path = WS / "results" / "sp1_endpoint.history.json"
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else default_path
    if not path.exists():
        print(f"ERROR: history not found: {path}")
        sys.exit(1)
    with open(path, encoding="utf-8") as f:
        h = json.load(f)

    if "groups" not in h:
        print("ERROR: history file does not contain 'groups' section --")
        print("       это значит, что оптимизатор был собран до изменений.")
        print("       Пересоберите: cmake --build build")
        sys.exit(1)

    final = bool(h.get("final", False))
    print(f"History: {path}")
    print(f"  final={'YES' if final else 'NO (checkpoint)'}  "
          f"best_loss={h.get('best_loss', 0):.6g}  "
          f"completed={h.get('completed_iter', 0)}/{h.get('n_iterations', 0)}")
    print(f"  groups: {len(h['groups'])}")

    sp_header, sp_rows = read_csv(SP_PATH)
    hw_header, hw_rows = read_csv(HW_PATH)

    total_sp = total_hw = 0
    for g in h["groups"]:
        s, hw = apply_group(g, sp_rows, sp_header, hw_rows, hw_header)
        total_sp += s
        total_hw += hw

    if total_sp > 0:
        write_csv(SP_PATH, sp_header, sp_rows)
        print(f"  -> sp_elements.csv:      {total_sp:4d} rows updated")
    else:
        print(f"  -> sp_elements.csv:      no changes")
    if total_hw > 0:
        write_csv(HW_PATH, hw_header, hw_rows)
        print(f"  -> hardware_library.csv: {total_hw:4d} rows updated")
    else:
        print(f"  -> hardware_library.csv: no changes")

    # Пересборка lattice_configs.csv
    print(f"\nRebuilding lattice_configs.csv ...")
    r = subprocess.run([sys.executable, str(ASSEMBLE)], capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode != 0:
        print(f"  ERROR: assemble failed: {r.stderr.strip()[:300]}")
        sys.exit(1)

    # Перегенерация trajectories_exp{id}.bin -- чтобы visualize.py видел актуальную траекторию.
    exp_id = int(h.get("experiment_id", 2))
    exe_path = SCRIPT_DIR / "build" / "SynchrotronTracker.exe"
    if exe_path.exists():
        print(f"\nRe-running tracker for exp {exp_id} ...")
        r = subprocess.run([str(exe_path), str(exp_id)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  WARN: tracker failed (visualize.py покажет старую траекторию): "
                  f"{r.stderr.strip()[:300]}")
        else:
            # Покажем только финальные строки (Part 1.. + total time)
            for line in r.stdout.strip().splitlines()[-8:]:
                print(f"  {line}")
    else:
        print(f"\nWARN: {exe_path} not found, skip re-tracking. "
              f"Запустите вручную: build\\SynchrotronTracker.exe {exp_id}")

    print(f"\nDone. Запустите python visualize.py {exp_id} чтобы увидеть финальную решётку.")


if __name__ == "__main__":
    main()
