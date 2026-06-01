"""Общие cross-platform утилиты для Python-скриптов проекта.

Назначение:
    - exe_name / exe_path -- собрать имя/путь исполняемого файла с учётом
      платформы (.exe на Windows, без суффикса на Linux/macOS).
    - tracker_exe / optimizer_exe -- удобные шорткаты на бинарники проекта.
    - project_root / workspace_dir -- ищут корень проекта, поднимаясь вверх
      от каталога вызывающего скрипта.

Скрипты на верхнем уровне репозитория должны импортировать только нужное:

    from _common import tracker_exe, project_root, workspace_dir
"""
from __future__ import annotations
import os
import sys
from pathlib import Path

__all__ = [
    "EXE_SUFFIX",
    "exe_name",
    "exe_path",
    "tracker_exe",
    "optimizer_exe",
    "project_root",
    "workspace_dir",
    "set_matplotlib_backend",
]

# ---------------------------------------------------------------------------
# Платформенные утилиты
# ---------------------------------------------------------------------------
EXE_SUFFIX = ".exe" if sys.platform == "win32" else ""


def exe_name(stem: str) -> str:
    """`SynchrotronTracker` -> `SynchrotronTracker.exe` (Windows) или без суффикса."""
    return f"{stem}{EXE_SUFFIX}"


def exe_path(stem: str, root: Path | None = None) -> Path:
    """Полный путь к build/{stem}{EXE_SUFFIX} относительно корня проекта."""
    return (root or project_root()) / "build" / exe_name(stem)


def tracker_exe(root: Path | None = None) -> Path:
    return exe_path("SynchrotronTracker", root)


def optimizer_exe(root: Path | None = None) -> Path:
    return exe_path("SynchrotronOptimizer", root)


# ---------------------------------------------------------------------------
# Поиск корня проекта / workspace
# ---------------------------------------------------------------------------
_MARKER_FILES = ("CMakeLists.txt", "README.md")


def project_root(start: Path | None = None) -> Path:
    """Подняться вверх от `start` (или от текущего скрипта) до каталога с
    маркером проекта (CMakeLists.txt или README.md). Возвращает абсолютный путь.

    Бросает FileNotFoundError, если корень не найден.
    """
    here = (start or Path(__file__)).resolve()
    cur = here.parent if here.is_file() else here
    for _ in range(15):
        if any((cur / m).exists() for m in _MARKER_FILES):
            return cur
        if cur.parent == cur:
            break
        cur = cur.parent
    raise FileNotFoundError(
        f"project_root: не найден каталог с {_MARKER_FILES} от {here}"
    )


def workspace_dir(root: Path | None = None) -> Path:
    """Абсолютный путь к каталогу workspace/."""
    return (root or project_root()) / "workspace"


# ---------------------------------------------------------------------------
# Matplotlib: безопасный выбор бэкенда
# ---------------------------------------------------------------------------
def set_matplotlib_backend(prefer: str = "TkAgg") -> str:
    """Установить интерактивный бэкенд matplotlib, если он доступен.

    Если предпочтительный бэкенд недоступен (например, на headless-Linux без
    Tk), безмолвно падаем на Agg (off-screen). Возвращает имя итого бэкенда.
    Вызывать ДО `import matplotlib.pyplot as plt`.
    """
    import matplotlib
    candidates = [prefer, "Qt5Agg", "Agg"]
    for name in candidates:
        try:
            matplotlib.use(name, force=True)
            return name
        except Exception:
            continue
    return matplotlib.get_backend()
