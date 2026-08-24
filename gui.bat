@echo off
rem Запуск графической оболочки SynchrotronSim
setlocal
set PYTHONIOENCODING=utf-8
py -3.14 "%~dp0gui_app.py" %*
