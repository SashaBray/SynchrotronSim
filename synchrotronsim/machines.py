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
