# -*- coding: utf-8 -*-
"""Кольцо на C++-ядре В ПРОЦЕССЕ Python (модуль расширения `_engine`) с интерфейсом `ring_cpp.CppRing`.

Шаг 0 плана библиотеки (docs/LIBRARY_ARCHITECTURE.md): то же ядро и те же числа, что у сервера
`MatrixTracker --ring-serve`, но без процесса и канала stdin/stdout. На время трекинга GIL отпущен,
поэтому разные кольца (`ExtRing`) можно гонять из разных потоков Python параллельно — это и есть
«параллелизм по копиям машины» для задач с одной частицей (замкнутая орбита, столбцы отклика).

Совпадение с сервером бит-в-бит — при одинаковой последовательности запросов трекинга: номер запроса
входит в зерно квантового излучения (ISR), как у сервера.

Шаг 2: решётка по умолчанию приходит из модели `Machine` (`synchrotronsim.db`) МАССИВАМИ (`line_spec`) — C++ не
читает диск ни при сборке, ни при RESET; RESET — откат к снимку свежего кольца (~6 мс против ~0,3 с пересборки).
`source="csv"` — C++ читает таблицы сам (как сервер; для сверок). Обычный вход в проект — `ring_cpp.CppRing`
(бэкенд "ext" поверх этого класса, с описаниями `mags`/`plane` и копиями `with_*`/`perturbed`).

  ExtRing(lattice, start_before_first_m, sr, isr, rf, design_mev, nsl, order, seed, ws, nsl_dev, machine, source)
    track(states, n_turns)            → (n, 6), NaN у потерянных
    track_tbt(states, n_turns, limit) → (n, n_turns+1, 6)
    track_elements(states)            → (n, n_mag+1, 6)
    track_pieces(states)              → (n, n_pieces+1, 6)
    piece_meta(), junction(states, A, B)
    apply_perturb(...), add_perturb(...), set_kicks(...), set_bend_scale(...), set_rf(...), reset()
    snapshot() / restore(snap)        — снимок и откат живого состояния (магниты с возмущениями, срезы, ВЧ)
    clone()                           — независимая копия (параллельная работа в потоках Python)
"""
import math
import os

import numpy as np

from . import _engine as _eng_mod

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
WORKSPACE = os.path.join(ROOT, "workspace")
MC2_MEV = 0.51099895069
C_LIGHT = 299792458.0


def _default_slices():
    """Как ring_cpp._default_slices: env SS_SECTOR_NSL/SS_SECTOR_ORDER или 5 и 4."""
    return (int(os.environ.get("SS_SECTOR_NSL", "") or 5), int(os.environ.get("SS_SECTOR_ORDER", "") or 4))


def _nsl_dev_from(nsl_dev):
    """Срезы по типу устройства {device_type_id: n}; None — из env SS_NSL_DEV (как CppRing)."""
    if nsl_dev is None:
        nsl_dev = dict(t.split("=", 1) for t in (os.environ.get("SS_NSL_DEV") or "").split(",") if "=" in t)
    out = {}
    for k, v in nsl_dev.items():
        try:
            n = int(v)
        except (TypeError, ValueError):
            continue
        if k and n > 0:
            out[str(k)] = n
    return out


def omp_threads():
    return _eng_mod.omp_threads()


def line_spec(machine):
    """Описание линии для C++ (`mc::LineSpec`) из `Machine.elements`: магниты в порядке пучка, обмотки сложены.
    n_slices = 0 («из настроек опыта») уходит как есть — ядро кольца берёт 40, как у сервера."""
    E = machine.elements
    f = lambda a: np.ascontiguousarray(a, dtype=np.float64)
    return {"dev": [str(s) for s in E.device_type], "kind": [str(s) for s in E.kind],
            "x": f(E.x), "y": f(E.y), "z": f(E.z), "yaw": f(E.yaw), "pitch": f(E.pitch), "roll": f(E.roll),
            "leff": f(E.eff_length), "B": f(E.mat_B), "G": f(E.mat_G), "Gv": f(E.mat_Gv), "S": f(E.mat_S),
            "O": f(E.mat_O), "e1": f(E.edge1), "e2": f(E.edge2), "lf1": f(E.lf1), "lf2": f(E.lf2),
            "nsl_box": [int(v) for v in E.n_slices]}


def cavity_specs(machine):
    """Резонаторы установки для C++: [(id, x, y, z, nx, ny, nz, доля напряжения)]."""
    return [(c.cavity_id, *map(float, c.point), *map(float, c.normal), float(c.volt_share)) for c in machine.cavities]


class ExtRing:
    label = "SynchrotronSim (расширение)"

    def __init__(self, lattice=1, start_before_first_m=None, sr=False, isr=False, rf=None,
                 design_mev=6000.0, nsl=None, order=None, seed=12345, ws=None, nsl_dev=None,
                 machine=None, source="db", _ring=None):
        """source="db" (по умолчанию) — решётка из модели `Machine` (`synchrotronsim.db`, ws ⇒ overlay):
        C++ получает массивы и не читает диск ни при сборке, ни при RESET; machine= — готовая Machine.
        source="csv" — C++ сам читает таблицы ws (как сервер; для сверок)."""
        dn, do = _default_slices()
        self.lattice = int(machine.lattice_id if machine is not None else lattice)
        self.lead = start_before_first_m
        self.sr, self.isr, self.rf, self.design_mev = bool(sr or isr), bool(isr), rf, float(design_mev)
        self.nsl = int(nsl or dn); self.order = int(order or do); self.seed = int(seed)
        self.nsl_dev = _nsl_dev_from(nsl_dev)
        self.ws = ws
        self.source = "machine" if machine is not None else source
        if _ring is not None:                                        # clone(): готовое C++-кольцо
            self._r = _ring
            self.machine = None
            self._read_info()
            return
        rf_ = rf or {}
        common = dict(lead=float("nan") if self.lead is None else float(self.lead), design_mev=self.design_mev,
                      sr=self.sr, isr=self.isr, nsl=self.nsl, order=self.order,
                      rf=bool(rf), volt_mv=float(rf_.get("volt_mv", 0.0)), harmonic=float(rf_.get("harmonic", 0.0)),
                      phase_deg=float(rf_.get("phase_deg", 0.0)), seed=self.seed,
                      nsl_dev=self.nsl_dev, nsl_dev_explicit=True)
        if self.source == "csv":
            self.machine = None
            self._r = _eng_mod.Ring(lattice=self.lattice, ws=os.path.abspath(ws or WORKSPACE), **common)
        elif self.source in ("db", "machine"):
            if machine is None:
                from .db import Repository
                machine = Repository(overlay=ws).machine(self.lattice)
            self.machine = machine
            self._r = _eng_mod.Ring(spec=line_spec(machine), cavities=cavity_specs(machine),
                                    lattice=self.lattice, **common)
        else:
            raise ValueError(f"ExtRing: source = db | csv, не {source!r}")
        self._fresh = self._r.snapshot()                             # свежее кольцо: RESET = откат к нему (мс)
        if self.rf and float(self.rf.get("freq_hz") or 0.0) > 0:     # как CppRing: частота — командой после сборки
            self._send_rf()
        self._read_info()

    # ---- копии и откат ----------------------------------------------------------------------------
    def clone(self):
        """Независимая копия кольца со всеми текущими возмущениями и ВЧ (для параллельной работы в потоках)."""
        import copy
        c = copy.copy(self)
        c.rf = dict(self.rf) if self.rf else self.rf
        c._r = self._r.clone()
        return c

    def snapshot(self):
        """Снимок живого состояния (магниты с возмущениями, срезы, ВЧ); вернуть — restore(snap)."""
        return (self._r.snapshot(), dict(self.rf) if self.rf else self.rf)

    def restore(self, snap):
        """Откат к снимку без пересборки. Номер запроса (зерно ISR) не откатывается."""
        s, rf = snap
        self._r.restore(s)
        self.rf = dict(rf) if rf else rf
        return self

    def _read_info(self):
        r = self._r
        self.n_mag, self.circumference, self.n_pieces = int(r.n_mag), float(r.circumference), int(r.n_pieces)
        inf = r.info()
        self.s_in, self.s_out = inf["s_in"], inf["s_out"]
        self.kinds = list(inf["kinds"])
        self.names = [n if n else "-" for n in inf["names"]]

    # ---- ВЧ ---------------------------------------------------------------------------------------
    def _send_rf(self):
        rf = self.rf or {}
        self._r.set_rf(True, float(rf.get("volt_mv", 0.0)), float(rf.get("harmonic", 0.0)),
                       float(rf.get("phase_deg", 0.0)), float(rf.get("freq_hz") or 0.0))

    def set_rf(self, **kw):
        if not self.rf:
            raise ValueError("set_rf: кольцо построено без ВЧ (rf=None)")
        self.rf = {**self.rf, **kw}
        self._send_rf()
        return self

    def rf_circumference(self):
        f = float((self.rf or {}).get("freq_hz") or 0.0)
        return self.rf["harmonic"] * C_LIGHT / f if f > 0 else self.circumference

    def rf_frequency(self):
        f = float((self.rf or {}).get("freq_hz") or 0.0)
        return f if f > 0 else self.rf["harmonic"] * C_LIGHT / self.circumference

    # ---- трекинг (пакеты; GIL отпущен внутри) -------------------------------------------------------
    @staticmethod
    def _states(states):
        return np.ascontiguousarray(np.asarray(states, float).reshape(-1, 6))

    def track(self, states, n_turns=1):
        return self._r.track(self._states(states), int(n_turns))

    def track_tbt(self, states, n_turns, limit=None):
        return self._r.track_tbt(self._states(states), int(n_turns), float(limit) if limit else 0.0)

    def track_elements(self, states):
        return self._r.track_elements(self._states(states))

    def track_pieces(self, states):
        return self._r.track_pieces(self._states(states))

    def piece_meta(self):
        return [tuple(t) for t in self._r.piece_meta()]

    def junction(self, states, A, B):
        fa = np.ascontiguousarray(np.asarray(A, float).reshape(12))
        fb = np.ascontiguousarray(np.asarray(B, float).reshape(12))
        return self._r.junction(self._states(states), fa, fb)

    def beta(self, delta):
        pc = self.design_mev * (1.0 + delta)
        return pc / math.sqrt(pc * pc + MC2_MEV ** 2)

    # ---- возмущения (семантика CppRing) ----------------------------------------------------------
    _CMD = {"KICK": "kick", "VKICK": "vkick", "YAW": "yaw", "DK": "dk", "DKV": "dkv", "DSK": "dsk",
            "DK2": "dk2", "DK3": "dk3", "BSCALE": "bscale"}
    _ADDITIVE = {"dks": "DK", "dkvs": "DKV", "dsks": "DSK", "dk2s": "DK2", "dk3s": "DK3"}

    def _apply(self, cmds):
        for cmd, j, v in cmds:
            getattr(self._r, self._CMD[cmd])(int(j), float(v))

    def reset(self):
        """RESET: магниты свежего кольца (откат к снимку после сборки, без пересборки), ВЧ — текущий."""
        self._r.reset_to(self._fresh)
        return self

    def apply_perturb(self, kicks=(), vkicks=(), dks=(), dsks=(), dk2s=(), dk3s=(), dkvs=(), bscales=()):
        """RESET + новый набор транзитных возмущений — как CppRing.apply_perturb (ВЧ переживает RESET)."""
        new = [(cmd, int(j), float(v)) for cmd, pairs in (("KICK", kicks), ("VKICK", vkicks), ("DK", dks),
                                                          ("DSK", dsks), ("DK2", dk2s), ("DK3", dk3s),
                                                          ("DKV", dkvs), ("BSCALE", bscales))
               for j, v in pairs]
        self.reset()
        self._apply(new)
        return self

    def add_perturb(self, **channels):
        new = []
        for name, pairs in channels.items():
            cmd = self._ADDITIVE.get(name)
            if cmd is None:
                raise ValueError(f"add_perturb: канал {name} не аддитивный (допустимы {sorted(self._ADDITIVE)})")
            new += [(cmd, int(j), float(v)) for j, v in pairs]
        self._apply(new)
        return self

    def set_kicks(self, kicks=(), vkicks=()):
        self._apply([("KICK", j, v) for j, v in kicks] + [("VKICK", j, v) for j, v in vkicks])
        return self

    def set_bend_scale(self, pairs=()):
        self._apply([("BSCALE", j, f) for j, f in pairs])
        return self

    def close(self):                           # совместимость с CppRing (процесса нет)
        self._r = None
