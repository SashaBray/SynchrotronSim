# -*- coding: utf-8 -*-
"""Стандартные диагностические графики кокпита — расчёт на стороне ДВИЖКА.

Адрес `diag.<name>`. Значение считает EngineThread (как `screen.<id>`) и отдаёт окну словарём ТОЛЬКО из чисел
и строк (шов Link). Данные статичны для сессии (матрица отклика FOFB, модельная дисперсия), поэтому движок
считает их один раз и кэширует; ключ кэша `cache_key` — то, при смене чего надо пересчитать (rcond контура).

`orbit_modes` — «η-узор и моды FOFB». Три панели:
  1) η-узор — орбита в BPM при сдвиге энергии пучка δ: x = η·δ. Модельная η_x в BPM — та же, что у датчика
     контура ВЧ (`rf_frequency._b_eta_bpm`: оптика доски, иначе twiss кольца plant'а);
  2) бетатронная волна — столбец матрицы отклика FOFB X (орбита от кика одного корректора);
  3) спектр мод SVD матрицы отклика σ/σ_max с порогом rcond контура: сколько мод берёт FOFB и какую долю
     η-узора он может нарисовать корректорами (то, что FOFB «отнимает» у контура ВЧ).
"""
from types import SimpleNamespace

import numpy as np

DIAG_NAMES = ("orbit_modes",)
ORBIT_LOOP = "fofb_x"


def modes_summary(R, eta, rcond):
    """Моды матрицы отклика R (BPM × корректоры) и их отношение к η-узору. Чистая функция.

    sv — σ/σ_max (по убыванию); n_kept — сколько мод берёт псевдообратная с этим rcond (σ > rcond·σ_max, как
    numpy.linalg.pinv); eta_reach — доля нормы η-узора в подпространстве взятых мод (0..1); eta_reach_full —
    то же для всех мод (сколько η вообще могут нарисовать корректоры)."""
    R = np.asarray(R, float)
    eta = np.asarray(eta, float)
    U, s, _ = np.linalg.svd(R, full_matrices=False)
    if s.size == 0 or s[0] <= 0:
        return {"sv": [], "n_kept": 0, "n_modes": int(s.size), "eta_reach": 0.0, "eta_reach_full": 0.0}
    keep = s > float(rcond) * s[0]
    nrm = float(np.linalg.norm(eta))
    e = eta / nrm if nrm > 0 else eta
    return {"sv": [float(v) for v in s / s[0]],
            "n_kept": int(keep.sum()),
            "n_modes": int(s.size),
            "eta_reach": float(np.linalg.norm(U[:, keep].T @ e)),
            "eta_reach_full": float(np.linalg.norm(U.T @ e))}


def eta_at_bpm(plant, state):
    """Модельная η_x в BPM — ровно та, что у датчика контура ВЧ (одна правда для окна и контура)."""
    from control.controllers.rf_frequency import _b_eta_bpm
    return np.asarray(_b_eta_bpm(SimpleNamespace(plant=plant, state=state, ring=plant.ring)), float)


def cache_key(name, loops):
    """Что инвалидирует кэш диагностики `name` (смена rcond контура FOFB X)."""
    if name == "orbit_modes":
        spec = (loops or {}).get(ORBIT_LOOP)
        ctrl = getattr(spec, "controller", None)
        return (name, float(ctrl.params.get("rcond", 0.0)) if ctrl is not None else None)
    return (name,)


def orbit_modes(plant, state, loops, resources, loop_name=ORBIT_LOOP):
    """Данные графика «η-узор и моды FOFB» (сериализуемый dict) или {"error": …}."""
    spec = (loops or {}).get(loop_name)
    ctrl = getattr(spec, "controller", None)
    if ctrl is None:
        return {"error": f"контур {loop_name} не найден"}
    act = list(getattr(ctrl, "act", None) or [])
    if not act:
        return {"error": f"у контура {loop_name} нет актуаторов (не настроен)"}
    rc = float(ctrl.params.get("rcond", 1e-2))
    R = np.asarray(resources.response(act, [ctrl.sensor]).R, float)      # кэш Resources: та же матрица, что у FOFB
    eta = eta_at_bpm(plant, state)
    if eta.size != R.shape[0]:
        return {"error": f"η в BPM ({eta.size}) не совпадает с числом строк матрицы отклика ({R.shape[0]})"}
    out = modes_summary(R, eta, rc)
    out.update(loop=loop_name, rcond=rc, eta_bpm=[float(v) for v in eta],
               cols=[[float(v) for v in R[:, j]] for j in range(R.shape[1])], corr_ids=[str(a) for a in act])
    return out


def compute(name, plant, state, loops, resources):
    """Диспетчер: `diag.<name>` → данные. Неизвестное имя — {"error"}."""
    if name == "orbit_modes":
        return orbit_modes(plant, state, loops, resources)
    return {"error": f"нет стандартного графика «{name}»"}
