# -*- coding: utf-8 -*-
"""ПРИМЕР ИМПЕРАТИВНОГО контура: читает из доски и пишет в свой слот ПРЯМО В КОДЕ (ctl), без возврата словаря.

Это «модуль, где пишешь что хочешь»: в update() ты сам решаешь, ЧТО прочитать из структуры текущего
состояния (датчики `view.read`, пучок у магнита `view.beam.at_magnet`, поле/положение `view.machine`…)
и ЧТО записать в структуру требуемого состояния (свой слот через `ctl.set`, уставки через `ctl.ref`).

Контракт: класс ставит `wants_ctl = True`, и планировщик вызывает `update(t, view, ctl)`, где `ctl`
(SlotWriter) пишет в ЭТОТ слот. Арбитраж сохранён: планировщик обнуляет слот ПЕРЕД update и суммирует
слоты ПОСЛЕ (с Limits/slew). Единственное правило — писать свой слот (ctl.set), а не чужие/actual.

Демо-закон НАМЕРЕННО простой (не тюнинг, а образец контракта): гасим горизонтальное положение пучка на
входе выбранного магнита одним выбранным корректором. Чтение — по ИМЕНИ магнита; запись — по id актуатора
(это и есть имя слота; имя корректора SP<nn>.COR<mm> → pid через plant.correctors, если нужно).
"""
from control.controller import Controller


class ImperativeSteer(Controller):
    wants_ctl = True                                      # → update(t, view, ctl): пишем императивно
    reads = ["ebpm.x"]

    def __init__(self, name="steer", source="SP01.und_QF1A.1", corrector="corr.000.hkick", gain=0.3):
        super().__init__(name, gain=gain)
        self.source = source                              # магнит, у входа которого смотрим пучок
        self.corrector = corrector                        # актуатор (pid слота), которым рулим

    def update(self, t, view, ctl):
        # 1) ПОЛУЧИТЬ из структуры текущего состояния — по имени, что захотели:
        if view.beam is not None:
            x = float(view.beam.at_magnet(self.source, "entry")[0])   # x пучка на входе магнита
        else:
            m = view.read("ebpm.x")                                   # запасной путь: орбита в мониторах
            x = float(m[0]) if m is not None else 0.0
        # 2) ОБРАБОТАТЬ (твой закон — тут простейший пропорциональный):
        u = -self.params["gain"] * x
        # 3) ЗАПИСАТЬ в структуру требуемого состояния — в СВОЙ слот, по id актуатора:
        ctl.set(self.corrector, u)
