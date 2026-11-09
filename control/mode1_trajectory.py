# -*- coding: utf-8 -*-
"""Режим 1 — трек пучка несколько оборотов + 3D-траектория в ГЛОБАЛЬНЫХ (лабораторных) координатах.

Платформенная точка входа поверх зрелого конвейера визуализации: движок ведёт инжектированный пучок
старт→финиш по кольцу (mte.deviation_paths, глобальные XYZ через лаб-рамки магнитов), сцена (магниты +
мониторы + траектория) собирается в самодостаточный интерактивный HTML (viz_html.build) с привязанным
2D-графиком отклонения. По умолчанию — блочные магниты (glyph, без зависимости от STL), устойчиво везде.

Машина задаётся экспериментом (experiments.csv → решётка+инжекция); тот же вход примет и конфигурацию,
которую оставил опыт с ОС (через свой exp). ЗАПУСК: py -3.14 control/mode1_trajectory.py [exp] [--detailed]
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")                      # viz_html печатает '→' — на cp1251 иначе падает print
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import viz_html


def main():
    ap = argparse.ArgumentParser(description="режим 1 — 3D-траектория пучка в лаб-координатах")
    ap.add_argument("exp", nargs="?", type=int, default=26, help="id эксперимента (experiments.csv)")
    ap.add_argument("--detailed", action="store_true", help="детальные STL-магниты вместо блоков (нужны STL)")
    ap.add_argument("--out", default=os.path.join("workspace", "results", "mode1_trajectory.html"))
    a = ap.parse_args()

    print(f"=== РЕЖИМ 1: 3D-траектория пучка (эксперимент {a.exp}) ===")
    print(f"движок ведёт пучок по кольцу в глобальных координатах; сцена → интерактивный 3D-HTML\n", flush=True)
    viz_html.main(a.exp, decimate=-1.0, out=a.out, glyph=not a.detailed)
    print(f"\n[OK] {a.out}")


if __name__ == "__main__":
    main()
