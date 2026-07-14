@echo off
REM Короткая обёртка для запуска оптимизатора по номеру эксперимента.
REM Пример:
REM   opt 9 --method=global --loss=orbit       (режимы: global | windowed | gradient)
REM
REM Все аргументы передаются в optimize.py, который сам находит JSON,
REM регенерирует flat-конфиг и запускает SynchrotronOptimizer.exe.
python "%~dp0optimize.py" %*
