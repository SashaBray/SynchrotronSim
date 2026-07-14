"""Запуск оптимизатора по номеру эксперимента.

Использование:
    python optimize.py <exp_id> [--method=global|windowed|gradient] [--loss=composite|orbit|envelope] ...
    python optimize.py <path-to-config.json|optimconfig> [flags...]   # legacy

Режимы (--method= или синоним --mode=):
    global   -- глобальный (вся конфигурация сразу); метод поиска задаётся
                global_search_mode = line_search | monte_carlo (случайные кики).
    windowed -- оконный (скользящее окно по магнитам).
    gradient -- градиентный (FD + Adam).
    Алиасы: mc/montecarlo -> global, per_element/pe -> windowed.

Что делает обёртка:
    1. Ищет в workspace/optim_inputs/*.json файл, у которого
       experiment_id == <exp_id>.
    2. Перегенерирует .optimconfig через prepare_optim_config.py
       (чтобы свежие правки JSON попали в flat-конфиг автоматически).
    3. Запускает build/SynchrotronOptimizer{.exe} с этим .optimconfig
       и любыми флагами, переданными дальше.

Все аргументы после первого передаются оптимизатору без изменений.

Запуск через CMD-обёртку opt.bat:
    opt 9 --method=global --loss=orbit
"""
from __future__ import annotations
import sys
import os
import json
import glob
import subprocess
from pathlib import Path

from _common import project_root, workspace_dir, optimizer_exe


def find_json_by_exp_id(inputs_dir: Path, exp_id: int) -> list[Path]:
    """Возвращает все .json в inputs_dir с experiment_id == exp_id."""
    matches: list[Path] = []
    for path in sorted(inputs_dir.glob("*.json")):
        try:
            with open(path, encoding="utf-8") as f:
                cfg = json.load(f)
        except Exception:
            continue
        if isinstance(cfg, dict) and cfg.get("experiment_id") == exp_id:
            matches.append(path)
    return matches


def usage_and_exit(rc: int = 1) -> None:
    print(__doc__.strip())
    sys.exit(rc)


def main() -> None:
    if len(sys.argv) < 2 or sys.argv[1] in ("-h", "--help"):
        usage_and_exit(0 if len(sys.argv) >= 2 else 1)

    root = project_root()
    inputs_dir = workspace_dir(root) / "optim_inputs"
    arg1 = sys.argv[1]
    forward_flags = sys.argv[2:]

    # --- 1. Определить путь к JSON/optimconfig ------------------------------
    arg_path = Path(arg1)
    if arg_path.is_file():
        # Передан явный путь к файлу — legacy-режим.
        target = arg_path
    else:
        # Должно быть число — exp_id.
        try:
            exp_id = int(arg1)
        except ValueError:
            print(f"[optimize] '{arg1}' не файл и не целое число.\n")
            usage_and_exit(1)

        matches = find_json_by_exp_id(inputs_dir, exp_id)
        if not matches:
            print(f"[optimize] В {inputs_dir} нет JSON с experiment_id={exp_id}.")
            # Подскажем доступные id.
            available: list[tuple[int, str]] = []
            for p in sorted(inputs_dir.glob("*.json")):
                try:
                    cfg = json.loads(p.read_text(encoding="utf-8"))
                    eid = cfg.get("experiment_id")
                    if isinstance(eid, int):
                        available.append((eid, p.name))
                except Exception:
                    continue
            if available:
                print("[optimize] Доступные experiment_id в optim_inputs:")
                for eid, name in sorted(available):
                    print(f"   {eid}: {name}")
            sys.exit(1)
        if len(matches) > 1:
            print(f"[optimize] Предупреждение: exp_id={exp_id} есть в нескольких JSON:")
            for m in matches:
                print(f"   {m.name}")
            print(f"[optimize] Использую первый: {matches[0].name}")
        target = matches[0]

    # --- 2. Перегенерация flat-конфига (если передан .json) -----------------
    if target.suffix.lower() == ".json":
        print(f"[optimize] Регенерация flat-конфига из {target.name} ...")
        rc = subprocess.call([sys.executable,
                              str(root / "prepare_optim_config.py"),
                              str(target)])
        if rc != 0:
            print(f"[optimize] prepare_optim_config.py упал (rc={rc})")
            sys.exit(rc)
        optimconfig = target.with_suffix(".optimconfig")
    else:
        optimconfig = target

    if not optimconfig.is_file():
        print(f"[optimize] Не найден flat-конфиг: {optimconfig}")
        sys.exit(1)

    # --- 3. Запуск оптимизатора ---------------------------------------------
    exe = optimizer_exe(root)
    if not exe.is_file():
        print(f"[optimize] Бинарник не найден: {exe}\n"
              f"           Собери проект: cmake --build build --config Release")
        sys.exit(1)
    if exe.stat().st_size == 0:
        print(f"[optimize] Бинарник имеет 0 байт (вероятно Defender): {exe}\n"
              f"           Пересобери и добавь build/ в исключения Defender.")
        sys.exit(1)

    cmd = [str(exe), str(optimconfig)] + forward_flags
    print(f"[optimize] Запуск: {' '.join(cmd)}")
    sys.exit(subprocess.call(cmd))


if __name__ == "__main__":
    main()
