#!/bin/sh
# Запуск графической оболочки SynchrotronSim (Linux / macOS).
# Интерпретатор: $SS_PYTHON, иначе python3.
export PYTHONIOENCODING=utf-8
exec "${SS_PYTHON:-python3}" "$(dirname "$0")/gui_app.py" "$@"
