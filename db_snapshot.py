# -*- coding: utf-8 -*-
"""Снимок ГЕНЕРИРУЕМЫХ таблиц БД: сохранить / восстановить / сверить.

ЗАЧЕМ. lattice_configs.csv, lattice_layout.csv, sp_elements.csv — генерируемые и гитигнорятся,
поэтому в истории проекта НЕТ записи о том, из чего считалась модель. Снимок закрывает эту дыру:
это единственный коммитабельный след состояния БД.

ПОЧЕМУ СКРИПТОМ, А НЕ РУКАМИ. Прежний снимок (lattice_configs_design.csv) сняли один раз вручную,
и он протух молча: пока чинили зеркальную укладку и грузили замкнутое кольцо, снимок остался со
СТАРОЙ геометрией (Σ поворота +2π вместо −2π) и противоположным знаком mat_B. Восстановление из
него вернуло бы уже исправленные баги. Отсюда правило: снимок только скриптом, и `check` обязан
кричать, когда он разошёлся с БД.

ЧТО ВНУТРИ. Снимок — папка с копиями таблиц + manifest.json (когда снят, из какой конфигурации,
git-коммит, суммы строк по решёткам). manifest позволяет понять, ЧТО именно лежит в снимке, не
открывая CSV.

ЗАПУСК
  py -3.14 db_snapshot.py save [--out <папка>] [--note "текст"]
  py -3.14 db_snapshot.py check [--snap <папка>]      # разошёлся ли снимок с БД (код 1 — да)
  py -3.14 db_snapshot.py restore --snap <папка> [--force]
"""
import os, csv, json, shutil, argparse, subprocess, collections
from _common import no_window

ROOT = os.path.dirname(os.path.abspath(__file__))
TB = os.path.join(ROOT, "workspace", "tables")
DEFAULT_SNAP = os.path.join(TB, "snapshot_design")

# Генерируемые таблицы, которых нет в git. ИСТОЧНИКИ ПРАВДЫ (element_connectivity, experiments,
# hardware_library) сюда не входят — они и так коммитятся.
TABLES = ["lattice_configs.csv", "lattice_layout.csv", "sp_elements.csv"]


def _git(*args):
    try:
        return subprocess.run(["git"] + list(args), cwd=ROOT, capture_output=True,
                              text=True, timeout=10, **no_window()).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def _stats(path):
    """{решётка: строк} — быстрый профиль таблицы (если в ней есть lattice_id)."""
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        rd = csv.DictReader(f, delimiter=";")
        if not rd.fieldnames or "lattice_id" not in rd.fieldnames:
            return {"rows": sum(1 for _ in rd)}
        return dict(sorted(collections.Counter(r["lattice_id"] for r in rd).items()))


def save(out, note, lte_dir=None):
    os.makedirs(out, exist_ok=True)
    man = {"note": note or "", "git_commit": _git("rev-parse", "--short", "HEAD"),
           "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
           "lte_dir": lte_dir or "", "tables": {}}
    for t in TABLES:
        src = os.path.join(TB, t)
        if not os.path.exists(src):
            print(f"  пропуск (нет файла): {t}")
            continue
        shutil.copy(src, os.path.join(out, t))
        man["tables"][t] = _stats(src)
        print(f"  сохранено: {t}  {man['tables'][t]}")
    json.dump(man, open(os.path.join(out, "manifest.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=2)
    print(f"Снимок: {os.path.relpath(out, ROOT)}  (коммит {man['git_commit'] or '?'})")
    return man


def check(snap):
    """Разошёлся ли снимок с текущей БД. Возврат: список расхождений."""
    if not os.path.isdir(snap):
        raise SystemExit(f"снимок не найден: {snap}")
    diffs = []
    for t in TABLES:
        a, b = os.path.join(snap, t), os.path.join(TB, t)
        if not os.path.exists(a):
            continue
        if not os.path.exists(b):
            diffs.append(f"{t}: в БД отсутствует"); continue
        sa, sb = _stats(a), _stats(b)
        if sa != sb:
            diffs.append(f"{t}: снимок {sa} ≠ БД {sb}")
        elif open(a, encoding="utf-8").read() != open(b, encoding="utf-8").read():
            diffs.append(f"{t}: те же размеры, но СОДЕРЖИМОЕ отличается")
    man = os.path.join(snap, "manifest.json")
    if os.path.exists(man):
        m = json.load(open(man, encoding="utf-8"))
        print(f"Снимок снят на коммите {m.get('git_commit') or '?'}"
              + (f" из {m['lte_dir']}" if m.get("lte_dir") else "")
              + (f" — {m['note']}" if m.get("note") else ""))
    if diffs:
        print("РАЗОШЁЛСЯ с текущей БД:")
        for d in diffs:
            print("  ✗ " + d)
        print("  Снимок устарел. Пересними: py -3.14 db_snapshot.py save")
    else:
        print("Снимок СОВПАДАЕТ с текущей БД ✓")
    return diffs


def restore(snap, force):
    if not os.path.isdir(snap):
        raise SystemExit(f"снимок не найден: {snap}")
    if not force:
        raise SystemExit("restore ПЕРЕЗАПИШЕТ текущие таблицы. Повтори с --force "
                         "(бэкап снимается автоматически).")
    bdir = os.path.join(TB, "backups", "pre_restore")
    os.makedirs(bdir, exist_ok=True)
    for t in TABLES:
        src = os.path.join(snap, t)
        if not os.path.exists(src):
            continue
        dst = os.path.join(TB, t)
        if os.path.exists(dst):
            shutil.copy(dst, os.path.join(bdir, t))
        shutil.copy(src, dst)
        print(f"  восстановлено: {t}")
    print(f"Бэкап прежних: {os.path.relpath(bdir, ROOT)}")


def main():
    ap = argparse.ArgumentParser(description="снимок генерируемых таблиц БД (их нет в git)")
    ap.add_argument("action", choices=["save", "check", "restore"])
    ap.add_argument("--snap", default=DEFAULT_SNAP, help="папка снимка")
    ap.add_argument("--out", default=None, help="куда сохранять (save; по умолч. --snap)")
    ap.add_argument("--note", default=None, help="пометка в manifest (save)")
    ap.add_argument("--lte-dir", default=None, help="из какой конфигурации собрана БД (save, в manifest)")
    ap.add_argument("--force", action="store_true", help="restore: подтвердить перезапись")
    a = ap.parse_args()
    if a.action == "save":
        save(a.out or a.snap, a.note, a.lte_dir)
    elif a.action == "check":
        raise SystemExit(1 if check(a.snap) else 0)
    else:
        restore(a.snap, a.force)


if __name__ == "__main__":
    main()
