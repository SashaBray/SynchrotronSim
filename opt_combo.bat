@echo off
REM Combined orchestrator: чередование оконной (windowed) и глобальной (global) фаз.
REM Пример:
REM   opt_combo 7              -- запустить exp 7 с параметрами из его JSON.combined
REM
REM Параметры читаются из блока "combined" в JSON-конфиге эксперимента.
python "%~dp0optimize_combined.py" %*
