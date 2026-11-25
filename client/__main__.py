# -*- coding: utf-8 -*-
"""Standalone-запуск живого кокпита режима 4:

    SS_SECTOR=1 PYTHONPATH=C:/Sasha/Work/SILA/SynchrotronSim py -3.14 -m client --exp 1
    py -3.14 -m client --lattice 26
    py -3.14 -m client --experiment workspace/experiments/my.json     # опыт из файла (окно «Опыт…» → сохранить)
    py -3.14 -m client --experiment fb:3                              # опыт ОС из таблиц fb_*
    py -3.14 -m client --experiment bundle.zip                        # пакет: опыт + документ машины, без БД
    py -3.14 -m client --site lab                                     # модель на вычислителе площадки lab

Тяжёлая модель кольца строится в фоне — окно появляется сразу. Из оболочки кокпит открывается кнопкой
«Режим 4»; этот модуль — тот же кокпит отдельным процессом.
"""
import argparse
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def main() -> int:
    ap = argparse.ArgumentParser(description="SynchrotronSim — живой кокпит режима 4")
    ap.add_argument("--exp", type=int, default=None, help="id опыта проекта (прежний хук оболочки; не используется)")
    ap.add_argument("--experiment", default=None,
                    help="опыт: путь к JSON (Experiment.to_json, кнопка «Опыт…» → «Сохранить опыт…»), пакет .zip "
                         "(Experiment.save_bundle: опыт + документ машины — модель считает без БД) или fb:N — опыт "
                         "ОС из таблиц fb_*; задан — параметры модели (--lattice, --nsl, --slices, …) берутся из опыта")
    ap.add_argument("--lattice", type=int, default=26, help="решётка (по умолчанию 26 — активная СИЛА)")
    ap.add_argument("--gap-scale", type=float, default=1.0,
                    help="множитель визуального зазора челюстей (1 = точно по модели)")
    ap.add_argument("--publish-every", type=int, default=3,
                    help="публиковать снимок на доски/графики раз в N тактов (меньше → свежее чаще; 1 = каждый такт)")
    ap.add_argument("--nsl", type=int, default=8,
                    help="срезов на элемент — ГЛАВНЫЙ рычаг модельной скорости (8 деф.; 4 ≈ ×2.5, 2 ≈ ×4.6; "
                         "ниже 2 ломается матрица отклика FOFB; меньше nsl — грубее оптика)")
    ap.add_argument("--n-part", type=int, default=12,
                    help="частиц в пучке (12 деф.); на скорость влияет слабо (трекинг векторизован), "
                         "но 1 убирает декогеренцию бетатрона — переходный процесс идеализируется")
    ap.add_argument("--stride", type=int, default=15,
                    help="оборотов на такт (15 деф.); на скорость ~не влияет, задаёт зернистость времени "
                         "(снимок раз в stride×publish-every оборотов)")
    ap.add_argument("--settle-turns", type=int, default=250, help="оборотов декогеренции при reset")
    ap.add_argument("--gap-slew", type=float, default=500.0,
                    help="скорость привода зазора, мм/с модельного времени (дефолт 500 ≈ компромисс: видно за "
                         "~13 с, ~540 об/мм; реальный привод ~0.5; большое значение — быстрее, но менее адиабатично)")
    ap.add_argument("--no-fofb", action="store_true",
                    help="не включать контуры обратной связи при старте (по умолчанию FOFB X/Y включены)")
    ap.add_argument("--no-optics", action="store_true",
                    help="без оптических контуров (tune/chroma/adts) — быстрый старт на одном TbtPlant "
                         "(с оптикой первый билд считает матрицы отклика нормальной формы — минуты, потом кэш)")
    ap.add_argument("--autorun", action="store_true", help="сразу запустить прогон (без нажатия «Пуск»)")
    ap.add_argument("--slices", default="auto",
                    help="срезы по типам устройств: auto (по умолчанию — workspace/configs/slices_<lat>.json, набор "
                         "tight из nsl_convergence.py; нет файла — как раньше) | none (как раньше: --nsl, полюса ИД 1) | "
                         "tight | loose (из workspace/results/nsl_convergence_<lat>.json) | путь к json")
    ap.add_argument("--site", default=None,
                    help="площадка (ss_server.py site add …): модель строится и считается на её вычислителе, "
                         "здесь — сцена, окна и ручки (шов Link по сети); машина — из хранилища площадки по ссылке")
    a = ap.parse_args()

    from PySide6.QtWidgets import QApplication
    from client.cockpit import launch_cockpit

    nsl_dev = "auto"
    if a.slices == "none":
        nsl_dev = None
    elif a.slices != "auto":
        from synchrotronsim.convergence import load_recommendation
        nsl_dev = load_recommendation(a.lattice, a.slices)
        print(f"срезы по типам ({a.slices}): {len(nsl_dev)} типов")
    app = QApplication.instance() or QApplication(sys.argv)
    app.setApplicationName("SynchrotronSim — Режим 4")
    win = launch_cockpit(lattice=a.lattice, exp=a.exp, gap_scale=a.gap_scale,
                         publish_every=a.publish_every, settle_turns=a.settle_turns,
                         nsl=a.nsl, n_part=a.n_part, stride=a.stride, fofb_on=not a.no_fofb,
                         optics=not a.no_optics, autorun=a.autorun, gap_slew_mm_s=a.gap_slew, nsl_dev=nsl_dev,
                         experiment=a.experiment, site=a.site)
    assert win is not None       # держим ссылку живой на время цикла событий
    return app.exec()


if __name__ == "__main__":
    raise SystemExit(main())
