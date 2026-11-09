# -*- coding: utf-8 -*-
"""Обнаружение плагинов-регуляторов: пользователь кладёт .py в каталог, мы находим классы Controller.

Так «добавить алгоритм» = положить файл. Интерфейс/планировщик затем перечисляют найденное и дают
собрать из них контуры с нужными частотами.
"""
import importlib.util
import inspect
import os

from .controller import Controller


def discover(directory):
    """{имя_класса: класс} для всех Controller-наследников в directory/*.py (кроме _*.py)."""
    found = {}
    for fn in sorted(os.listdir(directory)):
        if not fn.endswith(".py") or fn.startswith("_"):
            continue
        path = os.path.join(directory, fn)
        spec = importlib.util.spec_from_file_location(f"_plugin_{fn[:-3]}", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        for name, obj in inspect.getmembers(mod, inspect.isclass):
            if issubclass(obj, Controller) and obj is not Controller:
                found[name] = obj
    return found
