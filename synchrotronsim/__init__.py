# -*- coding: utf-8 -*-
"""SynchrotronSim как библиотека (зачаток; план — docs/LIBRARY_ARCHITECTURE.md).

Шаг 0: C++-ядро кольца как модуль расширения `synchrotronsim._engine` (то же ядро, что сервер
`MatrixTracker --ring-serve`) и обёртка `synchrotronsim.engine.ExtRing` с интерфейсом `ring_cpp.CppRing`.
Шаг 1: модель данных `synchrotronsim.model` (Machine, InstanceTable, ElementTable, …) и единый слой доступа
к таблицам `synchrotronsim.db.Repository` (`open_db`).

    import synchrotronsim as ss
    m = ss.open_db().machine(26)

Сборка расширения (MinGW, из корня проекта):
    cmake -S . -B build-py -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release -DSS_BUILD_PYTHON=ON ^
          -DPython_EXECUTABLE=<путь к python.exe>
    cmake --build build-py --target _engine
"""
try:
    from . import _engine                     # noqa: F401  (C++-модуль; нет сборки — библиотека без движка)
    ENGINE_AVAILABLE = True
    ENGINE_ERROR = None
except ImportError as e:                      # не собрано или собрано под другой Python
    _engine = None
    ENGINE_AVAILABLE = False
    ENGINE_ERROR = str(e)

from .db import Repository, open_db           # noqa: E402
from .model import Machine                    # noqa: E402
from .site import Site                        # noqa: E402  (площадка: где хранилище и вычислитель — этап 3)

__all__ = ["ENGINE_AVAILABLE", "ENGINE_ERROR", "Repository", "open_db", "Machine", "Site"]
