# -*- coding: utf-8 -*-
"""Паспортные данные установок по номеру решётки: энергия, ВЧ (напряжение МВ, кратность), суперпериоды и
пределы проб амплитуды для паспорта. Единственный источник (перенесено из facility_report.LATTICE_INFO
2026-11-24 — facility_report импортирует отсюда). По плану (docs/LIBRARY_ARCHITECTURE.md, раздел 7)
место этих данных — таблица БД; пока это код."""

LATTICE_INFO = {
    25: dict(label="СИЛА 2025.08.29", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
    26: dict(label="СИЛА 2026.06.11", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
    2:  dict(label="Сибирь", energy_mev=3000.0, rf=(2.5, 75), superperiods=6,
             # секступоли нулевые: устойчивость теряется только от кинематики, на ~0.4–0.5 м —
             # поэтому предел пробы амплитуды сотни мм, иначе карта выживания однотонная (всё выжило)
             x_probe_mm=600.0, y_probe_mm=800.0, d_probe_pct=15.0, loss_m=1.0),
    1:  dict(label="СИЛА (реш. 1)", energy_mev=6000.0, rf=(8.0, 1296), superperiods=40,
             x_probe_mm=30.0, y_probe_mm=20.0, d_probe_pct=12.0, loss_m=0.05),
}


# ---------------------------------------------------------------------------------------------- по объекту Machine
# Модель во время счёта берёт паспорт, группы и семейства ИЗ ДОКУМЕНТА МАШИНЫ (docs/DATA_PROTOCOL.md, этап 1),
# а не из таблиц и не из facility_report (его импорт тянул цепочку модулей, читающих таблицы при импорте).

def info_of(machine) -> dict:
    """Паспорт машины в прежнем виде `facility_report.info_of` (+ `_id`): label, energy_mev, rf=(МВ, кратность),
    superperiods; прочие ключи LATTICE_INFO (пределы проб паспорта) — если решётка там есть."""
    meta = getattr(machine, "meta", None)
    if meta is None:
        raise ValueError(f"у машины {getattr(machine, 'lattice_id', '?')} нет паспорта (энергия/ВЧ/суперпериоды) — "
                         f"добавьте его в synchrotronsim/machines.py LATTICE_INFO")
    inf = dict(LATTICE_INFO.get(int(machine.lattice_id), {}))
    inf.update(label=meta.label, energy_mev=float(meta.design_energy_mev),
               rf=(float(meta.rf_volt_mv), int(meta.rf_harmonic)), superperiods=int(meta.superperiods))
    inf["_id"] = int(machine.lattice_id)
    return inf


def groups_of(machine) -> list:
    """Жёсткие группы машины в прежнем виде `errors.load_groups`: [{group_id, lat, pivot, yaw, pitch, roll, insts}]."""
    import numpy as np
    return [dict(group_id=g.group_id, lat=str(machine.lattice_id), pivot=np.array(g.pivot, float),
                 yaw=float(g.yaw), pitch=float(g.pitch), roll=float(g.roll), insts=set(g.instance_ids))
            for g in machine.groups]


def family(dev: str) -> str:
    """Имя семейства магнита по типу устройства: срезаем префикс sila26_/s2_ и хвостовой индекс среза
    (_1.._5, _1A…). (Перенесено из facility_report._family — тот ссылается сюда.)"""
    import re
    base = dev.split("_", 1)[1] if dev.split("_", 1)[0].endswith(("26", "25", "s2", "2")) and "_" in dev else dev
    return re.sub(r"_(\d+[A-Z]?|[A-Z])$", "", base)
