# -*- coding: utf-8 -*-
"""Графическая оболочка SynchrotronSim.

    py -3.14 gui_app.py          (Windows: gui.bat или ярлык SynchrotronSim.lnk)

Окно ведёт весь конвейер: проект (файл решётки .lte ↔ решётка в базе), правка таблиц базы, матрицы
связи и проектные токи, запуск трекинга тремя ядрами, визуализация и HTML-отчёты. Считают при этом
штатные инструменты проекта — оболочка запускает их процессами и показывает вывод в общем журнале.
"""
import sys


def _fail(text: str) -> None:
    """Сообщить об ошибке запуска и в консоль, и окном.

    Ярлык запускает оболочку через `pythonw.exe` — консоли нет, и обычный `print` в никуда: без окна
    двойной щелчок по ярлыку выглядел бы как «ничего не произошло». MessageBox берём из системного
    api напрямую, потому что упасть мы могли именно на отсутствии Qt.
    """
    print(text, file=sys.stderr)
    if sys.platform == "win32":
        try:
            import ctypes
            ctypes.windll.user32.MessageBoxW(None, text, "SynchrotronSim — не удалось запустить", 0x10)
        except Exception:
            pass


if __name__ == "__main__":
    try:
        from gui.main import run
    except ImportError as e:
        _fail(f"Не найден PySide6 или модуль оболочки:\n{e}\n\n"
              f"Установить: py -3.14 -m pip install PySide6")
        sys.exit(2)
    try:
        sys.exit(run())
    except Exception as e:                       # иначе при запуске без консоли ошибка пропадёт молча
        import traceback
        _fail(f"Оболочка завершилась с ошибкой:\n\n{type(e).__name__}: {e}\n\n"
              f"{traceback.format_exc(limit=6)}")
        sys.exit(1)
