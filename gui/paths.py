# -*- coding: utf-8 -*-
"""Пути проекта и интерпретатор для дочерних процессов GUI."""
from __future__ import annotations

import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
WORKSPACE = ROOT / "workspace"
TABLES = WORKSPACE / "tables"
RESULTS = WORKSPACE / "results"
CONFIGS = WORKSPACE / "configs"
LTE_DIR = WORKSPACE / "lte"
HTML_DIR = RESULTS / "html"
BUILD = ROOT / "build"

EXE_SUFFIX = ".exe" if sys.platform == "win32" else ""


def python_exe() -> str:
    """Интерпретатор для дочерних процессов — тот же, что запустил GUI (в нём все зависимости)."""
    return sys.executable


def exe_path(stem: str) -> Path:
    """build/<stem>[.exe]; на MSVC-сборке бинарь может лежать в build/Release/."""
    p = BUILD / f"{stem}{EXE_SUFFIX}"
    if not p.exists():
        alt = BUILD / "Release" / f"{stem}{EXE_SUFFIX}"
        if alt.exists():
            return alt
    return p


def table(name: str) -> Path:
    return TABLES / name


def child_env() -> dict:
    """Окружение дочернего процесса: UTF-8 вывод (иначе логи падают в cp1251-консоли Windows)."""
    env = dict(os.environ)
    env["PYTHONIOENCODING"] = "utf-8"
    env["PYTHONUNBUFFERED"] = "1"
    return env
