# -*- coding: utf-8 -*-
"""Реестр проектов оболочки.

Проект = связка «файл(ы) решётки `.lte` ↔ решётка и опыт в базе данных». Сама база одна на всё
приложение (`workspace/tables/`), в ней конфигурации сосуществуют по `lattice_id`; реестр лишь
помнит, какие решётки откуда загружены, чтобы не искать это каждый раз руками.

Формат решётки определяет сам `load_config` (импортируем его, чтобы правила детекта жили в одном
месте): `sila` — описание + последовательность (`*_desc.lte` + `*_line.lte`), `flat` — один файл с
`LINE=(...)`.
"""
from __future__ import annotations

import csv
import datetime as dt
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path  # noqa: F401 — используется в discover_projects

from . import paths

sys.path.insert(0, str(paths.ROOT))
import load_config as lc                                    # noqa: E402 — детект форматов .lte (единый источник)

REGISTRY = paths.CONFIGS / "gui_projects.json"


@dataclass
class Project:
    """Запись реестра. `path` — папка или файл `.lte`; `line` — файл последовательности (формат sila)."""
    name: str
    path: str
    line: str = ""
    fmt: str = "sila"
    lattice: int = 1
    exp: int = 1
    group: int = 1
    note: str = ""
    created: str = field(default_factory=lambda: dt.datetime.now().strftime("%Y-%m-%d %H:%M"))

    def label(self) -> str:
        return f"{self.name}  —  решётка {self.lattice}, опыт {self.exp}"

    def describe(self) -> str:
        fmt = lc.FMT_RU.get(self.fmt, self.fmt)
        s = [f"Проект:        {self.name}",
             f"Решётка (.lte): {self.path}"]
        if self.line:
            s.append(f"Последовательность: {self.line}")
        s += [f"Формат:        {fmt}",
              f"В базе:        решётка {self.lattice}, опыт {self.exp}, пучок {self.group}",
              f"Создан:        {self.created}"]
        if self.note:
            s.append(f"Примечание:    {self.note}")
        return "\n".join(s)


class Registry:
    """Список проектов в `workspace/configs/gui_projects.json` (создаётся при первом сохранении)."""

    def __init__(self, path: Path = REGISTRY):
        self.path = path
        self.items: list[Project] = []
        self.load()

    def load(self) -> None:
        self.items = []
        if not self.path.exists():
            return
        try:
            data = json.loads(self.path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            return
        for d in data.get("projects", []):
            known = {k: d[k] for k in Project.__dataclass_fields__ if k in d}
            try:
                self.items.append(Project(**known))
            except TypeError:
                continue                                    # запись из будущей версии — просто пропускаем

    def save(self) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps({"projects": [asdict(p) for p in self.items]},
                                        ensure_ascii=False, indent=2), encoding="utf-8")

    def add(self, project: Project) -> None:
        self.items = [p for p in self.items if p.name != project.name]
        self.items.append(project)
        self.save()

    def remove(self, name: str) -> None:
        self.items = [p for p in self.items if p.name != name]
        self.save()

    def get(self, name: str) -> Project | None:
        return next((p for p in self.items if p.name == name), None)


# --- сведения из базы данных -------------------------------------------------
def lattices_in_db() -> dict[int, int]:
    """Какие решётки лежат в `lattice_configs.csv` → {lattice_id: число элементов}."""
    path = paths.table("lattice_configs.csv")
    out: dict[int, int] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            try:
                lid = int(r["lattice_id"])
            except (KeyError, TypeError, ValueError):
                continue
            out[lid] = out.get(lid, 0) + 1
    return out


def experiments() -> list[dict]:
    """Строки `experiments.csv` как список словарей (для выпадающих списков опытов)."""
    path = paths.table("experiments.csv")
    if not path.exists():
        return []
    with open(path, encoding="utf-8-sig", newline="") as f:
        return [r for r in csv.DictReader(f, delimiter=";") if (r.get("exp_id") or "").strip()]


def experiment_config(exp: int) -> tuple[int | None, int | None]:
    """(установка, пучок) из строки опыта в `experiments.csv`."""
    for e in experiments():
        if (e.get("exp_id") or "").strip() == str(exp):
            try:
                return int(e["lattice_id"]), int(e["group_id"])
            except (KeyError, TypeError, ValueError):
                return None, None
    return None, None


def particle_groups() -> dict[int, list[dict]]:
    """Пучки из `particle_groups.csv` → {group_id: [частицы]}."""
    path = paths.table("particle_groups.csv")
    out: dict[int, list[dict]] = {}
    if not path.exists():
        return out
    with open(path, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f, delimiter=";"):
            try:
                gid = int(r["group_id"])
            except (KeyError, TypeError, ValueError):
                continue
            out.setdefault(gid, []).append(r)
    return out


def next_free(ids, start: int = 1) -> int:
    """Наименьший свободный целый идентификатор, начиная со `start`."""
    used = {int(i) for i in ids}
    n = start
    while n in used:
        n += 1
    return n


def discover_configs() -> list[dict]:
    """Конфигурации в `workspace/lte/` (делегируем детект самому load_config)."""
    return lc.discover()


# --- решётки, уже лежащие в базе, как готовые проекты ------------------------
MIN_ELEMENTS = 3          # меньше — это демо-конфигурации element_compare (один магнит на решётку)


def _json(path: Path) -> dict:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _sub_lattices() -> dict[int, str]:
    """{lattice_id: описание} суб-решёток (укороченных копий основной) из sub_lattices.json."""
    out: dict[int, str] = {}
    for s in _json(paths.CONFIGS / "sub_lattices.json").get("sub_lattices", []):
        try:
            out[int(s["lattice_id"])] = str(s.get("description", ""))
        except (KeyError, TypeError, ValueError):
            continue
    return out


def _short(text: str, limit: int = 58) -> str:
    """Описание суб-решётки → короткая подпись проекта.

    Берём первое предложение; если в нём есть двоеточие («Конфигурация 3 = SP1: 1 стандартный SP…»),
    оставляем часть до него — она и есть название. Длинное режем по границе слова.
    """
    t = text.strip().split(". ")[0].strip()
    head = t.split(":")[0].strip()
    if len(head) >= 8:
        t = head
    if len(t) <= limit:
        return t
    cut = t[:limit].rsplit(" ", 1)[0]
    return (cut or t[:limit]) + "…"


def discover_projects(known: list[Project] | None = None) -> list[Project]:
    """Решётки из базы, оформленные как проекты (для тех, кто ещё не заведён в реестре).

    База и реестр — разные вещи: решётки могли быть загружены до появления оболочки или из консоли.
    Поэтому проекты не выдумываются, а восстанавливаются из фактов базы и конфигураций:
      • состав решёток — `lattice_configs.csv` (демо-конфигурации в один элемент отбрасываются);
      • суб-решётки и их описания — `sub_lattices.json` (наследуют `.lte` родителя);
      • основная решётка — `lte_mapping.json[lattice_id]`, её `.lte` — из манифеста снимка базы,
        иначе первая подходящая конфигурация в `workspace/lte/`;
      • сторонние (плоские) решётки сопоставляются плоским конфигурациям, и только когда
        соответствие однозначно (ровно одна такая решётка и ровно одна такая конфигурация);
      • опыт и пучок — первая строка `experiments.csv`, ссылающаяся на решётку.
    Если источник `.lte` определить нельзя, поле остаётся пустым — проект всё равно рабочий,
    путь можно указать позже (он нужен только для перезагрузки решётки из `.lte`).
    """
    in_db = {lid: n for lid, n in lattices_in_db().items() if n >= MIN_ELEMENTS}
    if not in_db:
        return []
    taken = {(p.lattice) for p in (known or [])}
    subs = _sub_lattices()
    main = int(_json(paths.CONFIGS / "lte_mapping.json").get("lattice_id") or 1)

    configs = discover_configs()
    sila = [c for c in configs if c["format"] == "sila"]
    flat = [c for c in configs if c["format"] == "flat"]

    snap_dir = (_json(paths.TABLES / "snapshot_design" / "manifest.json").get("lte_dir") or "").strip()
    snap_name = Path(snap_dir).name if snap_dir else ""
    main_cfg = next((c for c in sila if c["name"] == snap_name), sila[0] if sila else None)

    # плоская решётка ↔ плоская конфигурация: сопоставляем только при однозначности
    others = [lid for lid in in_db if lid != main and lid not in subs]
    flat_cfg = flat[0] if (len(flat) == 1 and len(others) == 1) else None

    exp_by_lat: dict[int, dict] = {}
    for e in experiments():
        try:
            lid = int(e["lattice_id"])
        except (KeyError, TypeError, ValueError):
            continue
        exp_by_lat.setdefault(lid, e)

    def cfg_path(cfg) -> tuple[str, str, str]:
        """(путь, файл последовательности, формат). Для desc+line отдаём ПАПКУ — load_config найдёт
        в ней оба файла сам, поэтому отдельный `--line` не нужен."""
        if not cfg:
            return "", "", "sila"
        path = cfg["dir"] if cfg["format"] == "sila" else (cfg["flat"][0] if cfg["flat"] else cfg["dir"])
        return path, "", cfg["format"]

    out: list[Project] = []
    for lid, count in sorted(in_db.items()):
        if lid in taken:
            continue                                       # эту решётку пользователь уже завёл сам
        if lid in subs:
            path, line, fmt = cfg_path(main_cfg)           # суб-решётка живёт в .lte родителя
            name = _short(subs[lid]) or f"Решётка {lid}"
        elif lid == main:
            path, line, fmt = cfg_path(main_cfg)
            name = f"Решётка {lid} — {Path(path).name}" if path else f"Решётка {lid}"
        else:
            path, line, fmt = cfg_path(flat_cfg)
            # плоская конфигурация — это файл; узнаваемее звать её по папке конфигурации, а не по файлу
            name = f"Решётка {lid} — {flat_cfg['name'] if flat_cfg else Path(path).stem}" if path \
                else f"Решётка {lid}"
        e = exp_by_lat.get(lid, {})
        out.append(Project(
            name=name, path=path, line=line, fmt=fmt, lattice=lid,
            exp=int(e.get("exp_id") or lid), group=int(e.get("group_id") or 1),
            note=f"найдено в базе: {count} элементов" + ("" if path else "; источник .lte не определён"),
        ))
    return out


def config_from_path(path: str, line: str = "") -> dict:
    """Путь (папка или файл) → описание конфигурации; SystemExit от load_config → ValueError."""
    try:
        return lc.from_path(path, line or None)
    except SystemExit as e:
        raise ValueError(str(e)) from None
