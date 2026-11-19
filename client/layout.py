# -*- coding: utf-8 -*-
"""Персистентность раскладки кокпита (вместо логина).

Раскладка = сериализованный список окон-подписчиков (каждое — {kind, params, geom, state}) + камера 3D +
видимость доков. Чистый JSON (без Qt): восстановление = пересоздать окна из описаний и ПЕРЕПОДПИСАТЬ на
ТЕКУЩИЙ прогон по адресам — работает только потому, что адреса стабильны между запусками (см. addresses.py).

Файлы: `workspace/layouts/<проект>/<имя>.json` (в .gitignore). Авто-`last` поднимается при входе и
сохраняется при выходе. Имя проекта — ключ (напр. 'sila26'); мульти-пользователь позже — ключуем ещё по имени.
"""
from __future__ import annotations

import json
import os
import re
from typing import List, Optional

_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LAYOUT_ROOT = os.path.join(_ROOT, "workspace", "layouts")
LAST = "last"
VERSION = 1


def layout_dir(project: str) -> str:
    return os.path.join(LAYOUT_ROOT, _safe(project))


_BAD = re.compile(r'[<>:"/\\|?*\x00-\x1f]+')       # запрещённые в имени файла Windows + управляющие


def _safe(s: str) -> str:
    """Безопасное имя файла/папки: чистим только недопустимые символы, Unicode-буквы (кириллица) сохраняем."""
    s = _BAD.sub("_", str(s or "default")).strip(" .")
    return s or "default"


def make_layout(windows: List[dict], camera=None, docks: Optional[dict] = None) -> dict:
    """Собрать описание раскладки. `windows` — список {kind, params, geom:[x,y,w,h], state}."""
    return {"version": VERSION, "windows": list(windows),
            "camera": list(camera) if camera is not None else None,
            "docks": dict(docks or {})}


def path_of(project: str, name: str) -> str:
    return os.path.join(layout_dir(project), f"{_safe(name)}.json")


def save_layout(project: str, name: str, data: dict) -> str:
    d = layout_dir(project)
    os.makedirs(d, exist_ok=True)
    p = path_of(project, name)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2)
    return p


def load_layout(project: str, name: str) -> Optional[dict]:
    p = path_of(project, name)
    if not os.path.exists(p):
        return None
    try:
        with open(p, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def list_layouts(project: str) -> List[str]:
    d = layout_dir(project)
    if not os.path.isdir(d):
        return []
    return sorted(os.path.splitext(f)[0] for f in os.listdir(d) if f.endswith(".json"))


def save_last(project: str, data: dict) -> str:
    return save_layout(project, LAST, data)


def load_last(project: str) -> Optional[dict]:
    return load_layout(project, LAST)
