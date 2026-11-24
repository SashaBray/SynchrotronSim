# -*- coding: utf-8 -*-
"""Контракт адресов сигналов клиента — теперь ЕДИНЫЙ реестр библиотеки `synchrotronsim.addresses` (шаг 4
docs/LIBRARY_ARCHITECTURE.md). Этот модуль оставлен для совместимости импортов клиента (`from client import
addresses as A`): всё — реэкспорт.

Адрес стабилен между запусками: раскладка окон (`client/layout.py`) сериализует адреса и переподписывается на
новый прогон по ним. Прежние имена — псевдонимы: `drive.<aid>.<поле>` → `act.<aid>.<поле>`,
`actual.<aid>` → `act.<aid>`; старые раскладки открываются (адреса приводятся `canonical` при восстановлении
окна), а `resolve`/`meta` понимают обе формы. Грамматика и области — в шапке `synchrotronsim/addresses.py`.
"""
from synchrotronsim.addresses import (  # noqa: F401
    AREAS, ACTUATOR_TYPES, AddrMeta, DIAG_TITLES, DRIVE_FIELDS, ELEMENT_ATTRS, LOOP_FIELDS, META_NAMES,
    OPTICS_AT, OPTICS_SCALARS, Parsed, UND_ATTRS, W_GRAPH, W_INFO, W_POINT, W_SCALAR, W_SCREEN, W_VECTOR,
    WRITABLE_DRIVE, WRITABLE_LOOP, WRITABLE_UNDULATOR, _ser, canonical, enumerate_addresses, is_address, kind_of,
    meta, parse, resolve, split_reducer,
)
