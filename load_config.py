# -*- coding: utf-8 -*-
"""Загрузчик конфигурации решётки в БД — по ПУТИ к .lte, без зашитых имён.

Форматы (определяются по содержимому, а не по названию):
  • sila  — описание + отдельная последовательность: *_desc.lte + *_line.lte
            Пайплайн: lte_to_lattice → assemble_lattice_configs → build_matrix_db → warmstart_ring
  • flat  — сторонний PyAT/elegant: один .lte, внутри LINE=(...)
            Пайплайн: s2lattice_to_db

СОСУЩЕСТВОВАНИЕ КОНФИГУРАЦИЙ. В таблицах живут РАЗНЫЕ установки одновременно, и загрузка одной
НЕ должна ломать другие. Раньше ломала: assemble переписывал lattice_configs целиком, а
build_matrix_db — element_connectivity, поэтому каждая загрузка СИЛЫ стирала Сибирь (build_line(2)
падал с KeyError). Теперь загрузка объявляет СВОИ решётки (сама + суб-решётки из sub_lattices.json)
и передаёт их вниз через --lattices; строки чужих установок переносятся как есть, со всеми
столбцами (mat_*, структурные).

Полная загрузка сбрасывает mat_* СВОЕЙ решётки к проектным K из .lte (warmstart). Если решётка уже
есть — нужен --force (снимаются бэкапы). --parse-only — только разбор .lte в таблицы.

Запуск:
  py -3.14 load_config.py --list
  py -3.14 load_config.py path/to/config_dir --force
  py -3.14 load_config.py path/to/desc.lte --line path/to/line.lte --force
  py -3.14 load_config.py path/to/ring.lte --lattice 2 --force      # сторонний flat .lte
  py -3.14 load_config.py --config v1 --force                        # по короткому имени (удобство)
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


def from_path(path, line=None):
    """ПУТЬ (папка или .lte-файл) → cfg. Главный способ адресации: никаких зашитых имён.

    • папка              — ищем в ней *_desc.lte + *_line.lte (формат sila) либо один .lte с LINE= (flat)
    • файл + --line      — sila: описание и последовательность разными файлами
    • файл с LINE=       — flat (PyAT/elegant, всё в одном файле)
    """
    path = os.path.abspath(path)
    if os.path.isdir(path):
        descs = sorted(set(glob.glob(os.path.join(path, "*_desc.lte")) + glob.glob(os.path.join(path, "*desc*.lte"))))
        lines = sorted(set(glob.glob(os.path.join(path, "*_line.lte")) + glob.glob(os.path.join(path, "*line*.lte"))))
        ltes = sorted(glob.glob(os.path.join(path, "*.lte")))
        flat = [p for p in ltes if p not in descs and p not in lines and _has_line_def(p)]
        fmt = "sila" if (descs and lines) else ("flat" if flat else "none")
        return dict(name=os.path.basename(path), dir=path, format=fmt,
                    descs=descs, lines=lines, ltes=ltes, flat=flat)
    if not os.path.isfile(path):
        raise SystemExit(f"путь не найден: {path}")
    d = os.path.dirname(path)
    if line:
        line = os.path.abspath(line)
        if not os.path.isfile(line):
            raise SystemExit(f"файл последовательности не найден: {line}")
        return dict(name=os.path.basename(d) or os.path.basename(path), dir=d, format="sila",
                    descs=[path], lines=[line], ltes=[path, line], flat=[])
    if _has_line_def(path):
        return dict(name=os.path.splitext(os.path.basename(path))[0], dir=d, format="flat",
                    descs=[], lines=[], ltes=[path], flat=[path])
    raise SystemExit(f"в {path} нет LINE=(...) — укажи файл последовательности через --line")


def resolve(name, configs):
    """УДОБСТВО (не основной путь): короткое имя подпапки workspace/lte/ → cfg.
    'v1'→configuration_v1, 'syberia'→configuration_v3_syberia. Основной способ — передать ПУТЬ."""
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


def managed_lattices(lat):
    """{решётки, которые ЭТА загрузка вправе перезаписывать} = сама + её суб-решётки.
    Всё остальное в таблицах — чужие установки, их трогать нельзя."""
    out = {str(lat)}
    try:
        cfg = json.load(open(os.path.join(ROOT, "workspace", "configs", "sub_lattices.json"), encoding="utf-8"))
        for e in cfg.get("sub_lattices", []):
            out.add(str(e.get("lattice_id") if isinstance(e, dict) else e[0]))
    except (OSError, ValueError, KeyError, IndexError, TypeError):
        pass
    return out


def sila_lattice_id():
    try:
        return int(json.load(open(MAPPING, encoding="utf-8")).get("lattice_id", 1))
    except (OSError, ValueError):
        return 1


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
        mine = ",".join(sorted(managed_lattices(lat)))
        print(f"  перезаписываем решётки: {mine}   (остальные в БД — чужие, не трогаем)")
        sh(["assemble_lattice_configs.py"])
        # --lattices: НЕ пересчитывать mat_*/края чужих установок — их писали свои импортёры
        sh(["build_matrix_db.py", "--lte-dir", cfg["dir"], "--lattices", mine])
        sh(["warmstart_ring.py", str(lat), "--lte-dir", cfg["dir"]])
        # ФИНАЛ: суб-решётки — копия родителя УЖЕ с физикой. assemble создаёт их до warmstart,
        # когда mat_* ещё пусты, а warmstart пишет только родителя ⇒ без этого шага суб-решётка
        # остаётся с нулевыми киками при гнущейся геометрии (так решётка 3 потеряла mat_B).
        sh(["assemble_lattice_configs.py", "--sync-subs"])
        return lat

    # flat (PyAT/elegant) — self-contained one-shot
    lat = str(lattice or 2)
    ex = str(exp or lat)
    if lattice_exists(lat) and not force:
        raise SystemExit(f"[защита] решётка {lat} уже есть в БД. Повтори с --force для перезаписи.")
    sh(["s2lattice_to_db.py", "--lte", cfg["flat"][0], "--lattice", lat, "--exp", ex])
    return lat


def main():
    ap = argparse.ArgumentParser(
        description="загрузчик конфигурации решётки в БД по ПУТИ к .lte (папка или файл)")
    ap.add_argument("path", nargs="?", default=None,
                    help="ПУТЬ: папка с .lte, либо сам .lte-файл (с --line, если последовательность отдельно)")
    ap.add_argument("--line", default=None, help="путь к файлу последовательности (формат desc+line)")
    ap.add_argument("--config", default=None,
                    help="УДОБСТВО вместо пути: короткое имя подпапки workspace/lte/ (v1 / syberia)")
    ap.add_argument("--list", action="store_true", help="показать конфигурации в workspace/lte/ и их формат")
    ap.add_argument("--lattice", type=int, default=None, help="переопределить id решётки")
    ap.add_argument("--exp", type=int, default=None, help="id эксперимента (flat; по умолч. = id решётки)")
    ap.add_argument("--gen-assets", action="store_true", help="СИЛА: генерировать синтетические карты/STL")
    ap.add_argument("--parse-only", action="store_true", help="только разбор .lte в таблицы (без assemble/warmstart)")
    ap.add_argument("--force", action="store_true", help="перезаписать уже существующую решётку (с бэкапами)")
    a = ap.parse_args()

    if a.list or (not a.path and not a.config):
        print("Конфигурации в workspace/lte/:\n")
        print(f"  {'папка':28} {'формат':22} файлы")
        for c in discover():
            key = c["descs"] + c["lines"] if c["format"] == "sila" else (c["flat"] or c["ltes"])
            files = ", ".join(os.path.basename(p) for p in key) or "—"
            print(f"  {c['name']:28} {FMT_RU[c['format']]:22} {files}")
        print("\nЗагрузка:  py -3.14 load_config.py <путь> [--line <файл>] [--force] [--parse-only]")
        print("           py -3.14 load_config.py --config v1 --force      # по короткому имени")
        if not (a.path or a.config):
            return

    cfg = from_path(a.path, a.line) if a.path else resolve(a.config, discover())
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
