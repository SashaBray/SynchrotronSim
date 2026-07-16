# -*- coding: utf-8 -*-
"""ЕДИНЫЙ загрузчик конфигураций решётки из подпапок workspace/lte/ — БЕЗ копирования .lte в корень.
Одна команда на любой формат: авто-определяет тип конфигурации и вызывает нужный пайплайн, читая .lte
ПРЯМО из папки конфигурации.

Форматы (определяются по содержимому папки):
  • sila  — СИЛА: *_desc.lte + *_line.lte (+ workspace/configs/lte_mapping.json).
            Пайплайн: lte_to_lattice → assemble → build_matrix_db → warmstart_ring  (реш. из lte_mapping).
  • flat  — сторонний PyAT/elegant: один .lte с LINE=(...).  Пайплайн: s2lattice_to_db (реш. 200 по умолч.).
  • none  — нет загружаемого .lte (напр. только xlsx/pdf) → пропуск.

По умолчанию доводит до ГОТОВОЙ К ТРЕКИНГУ решётки (mat_*). Полная загрузка СИЛА СБРАСЫВАет тюн
(build_matrix_db перезаписывает откалиброванную coupling_matrix, warmstart — mat_* проектными K); поэтому
если целевая решётка УЖЕ есть в БД — требуется --force (перед перезаписью делаются бэкапы). --parse-only —
только разбор .lte в таблицы (для СИЛА: sp_elements+lattice_layout, недеструктивно к lattice_configs).

Запуск:
  py -3.14 load_config.py --list
  py -3.14 load_config.py --config v1 [--gen-assets] [--force] [--parse-only]
  py -3.14 load_config.py --config syberia [--lattice 200]
"""
import os, sys, csv, glob, json, re, shutil, subprocess, argparse

ROOT = os.path.dirname(os.path.abspath(__file__))
LTE_ROOT = os.path.join(ROOT, "workspace", "lte")
LATTICE_CONFIGS = os.path.join(ROOT, "workspace", "tables", "lattice_configs.csv")
MAPPING = os.path.join(ROOT, "workspace", "configs", "lte_mapping.json")
FMT_RU = {"sila": "СИЛА desc+line", "flat": "плоский .lte (PyAT)", "none": "— нет решётки"}


def _has_line_def(path):
    try:
        return bool(re.search(r"\bLINE\s*=", open(path, encoding="utf-8", errors="ignore").read(), re.I))
    except OSError:
        return False


def discover():
    """Скан подпапок workspace/lte/ → [{name, dir, format, descs, lines, ltes}]."""
    out = []
    for d in sorted(glob.glob(os.path.join(LTE_ROOT, "*"))):
        if not os.path.isdir(d):
            continue
        descs = sorted(set(glob.glob(os.path.join(d, "*_desc.lte")) + glob.glob(os.path.join(d, "*desc*.lte"))))
        lines = sorted(set(glob.glob(os.path.join(d, "*_line.lte")) + glob.glob(os.path.join(d, "*line*.lte"))))
        ltes = sorted(glob.glob(os.path.join(d, "*.lte")))
        flat = [p for p in ltes if p not in descs and p not in lines and _has_line_def(p)]
        fmt = "sila" if (descs and lines) else ("flat" if flat else "none")
        out.append(dict(name=os.path.basename(d), dir=d, format=fmt,
                        descs=descs, lines=lines, ltes=ltes, flat=flat))
    return out


def resolve(name, configs):
    """Гибкое сопоставление: 'v1'→configuration_v1, 'syberia'/'v3'→configuration_v3_syberia."""
    nl = name.lower()
    for c in configs:
        if c["name"].lower() == nl:
            return c
    hits = [c for c in configs if nl in c["name"].lower()]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise SystemExit(f"неоднозначно '{name}': {[c['name'] for c in hits]}")
    raise SystemExit(f"конфигурация '{name}' не найдена. Есть: {[c['name'] for c in configs]}")


def sila_lattice_id():
    try:
        return int(json.load(open(MAPPING, encoding="utf-8")).get("lattice_id", 101))
    except (OSError, ValueError):
        return 101


def lattice_exists(lat_id):
    if not os.path.exists(LATTICE_CONFIGS):
        return False
    with open(LATTICE_CONFIGS, encoding="utf-8") as f:
        return any(r.get("lattice_id") == str(lat_id) for r in csv.DictReader(f, delimiter=";"))


def _backup(rel):
    src = os.path.join(ROOT, rel)
    if os.path.exists(src):
        bdir = os.path.join(ROOT, "workspace", "tables", "backups")
        os.makedirs(bdir, exist_ok=True)
        dst = os.path.join(bdir, os.path.basename(rel) + ".pre_loadconfig")
        shutil.copy(src, dst)
        print(f"  бэкап: {rel} → {os.path.relpath(dst, ROOT)}")


def sh(args):
    """Запуск шага пайплайна тем же интерпретатором. Останов при ошибке."""
    print(f"\n$ {os.path.basename(sys.executable)} {' '.join(args)}")
    r = subprocess.run([sys.executable] + args, cwd=ROOT)
    if r.returncode != 0:
        raise SystemExit(f"[стоп] шаг завершился с кодом {r.returncode}: {args[0]}")


def load(cfg, lattice, exp, gen_assets, parse_only, force):
    fmt = cfg["format"]
    if fmt == "none":
        raise SystemExit(f"'{cfg['name']}': нет загружаемого .lte (формат none) — грузить нечего")

    if fmt == "sila":
        sh(["lte_to_lattice.py", "--lte-dir", cfg["dir"]] + (["--gen-assets"] if gen_assets else []))
        if parse_only:
            print("\n[parse-only] sp_elements.csv + lattice_layout.csv обновлены (arg_val сохранены). "
                  "assemble/build_matrix_db/warmstart НЕ запускались.")
            return None
        lat = lattice or sila_lattice_id()
        if lattice_exists(lat) and not force:
            raise SystemExit(f"[защита] решётка {lat} уже есть в БД. Полная загрузка СБРОСИТ тюн "
                             f"(mat_*). Повтори с --force (будут бэкапы) или используй --parse-only.")
        sh(["assemble_lattice_configs.py"])
        sh(["build_matrix_db.py"])
        sh(["warmstart_ring.py", str(lat), "--lte-dir", cfg["dir"]])
        return lat

    # flat (PyAT/elegant) — self-contained one-shot
    lat = str(lattice or 200)
    ex = str(exp or lat)
    if lattice_exists(lat) and not force:
        raise SystemExit(f"[защита] решётка {lat} уже есть в БД. Повтори с --force для перезаписи.")
    sh(["s2lattice_to_db.py", "--lte", cfg["flat"][0], "--lattice", lat, "--exp", ex])
    return lat


def main():
    ap = argparse.ArgumentParser(description="единый загрузчик конфигураций решётки из подпапок workspace/lte/")
    ap.add_argument("--config", help="имя конфигурации (v1 / syberia / configuration_v3_syberia; подстрока)")
    ap.add_argument("--list", action="store_true", help="показать доступные конфигурации и их формат")
    ap.add_argument("--lattice", type=int, default=None, help="переопределить id решётки")
    ap.add_argument("--exp", type=int, default=None, help="id эксперимента (flat; по умолч. = id решётки)")
    ap.add_argument("--gen-assets", action="store_true", help="СИЛА: генерировать синтетические карты/STL")
    ap.add_argument("--parse-only", action="store_true", help="только разбор .lte в таблицы (без assemble/warmstart)")
    ap.add_argument("--force", action="store_true", help="перезаписать уже существующую решётку (с бэкапами)")
    a = ap.parse_args()

    configs = discover()
    if a.list or not a.config:
        print("Конфигурации в workspace/lte/:\n")
        print(f"  {'папка':28} {'формат':22} файлы")
        for c in configs:
            key = c["descs"] + c["lines"] if c["format"] == "sila" else (c["flat"] or c["ltes"])
            files = ", ".join(os.path.basename(p) for p in key) or "—"
            print(f"  {c['name']:28} {FMT_RU[c['format']]:22} {files}")
        print("\nЗагрузка:  py -3.14 load_config.py --config <имя> [--force] [--parse-only]")
        if not a.config:
            return

    cfg = resolve(a.config, configs)
    print(f"\n=== Загрузка «{cfg['name']}» (формат: {FMT_RU[cfg['format']]}) из {os.path.relpath(cfg['dir'], ROOT)} ===")
    lat = load(cfg, a.lattice, a.exp, a.gen_assets, a.parse_only, a.force)
    if lat is not None:
        print(f"\n[OK] Решётка {lat} загружена и готова к матричному трекингу:")
        print(f"     py -3.14 track3.py --mode kick-phys --lattice {lat} --exp <N>")
        print(f"     py -3.14 viz_html.py <N> --both        # визуализация")
        if cfg["format"] == "sila":
            print(f"     Для токов/диффура: recalibrate_coupling.py → phys_to_current.py {lat} → track3 --mode diff")
            print(f"     (полная загрузка СБРОСИЛА coupling_matrix к расчётной — бэкап в tables/backups/)")


if __name__ == "__main__":
    main()
