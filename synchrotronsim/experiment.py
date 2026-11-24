# -*- coding: utf-8 -*-
"""Опыт (`Experiment`) — что делаем с машиной: слой 2 модели данных (docs/LIBRARY_ARCHITECTURE.md, раздел 3.4;
шаг 5). Один и тот же объект строят блокнот, `mode4_run` (из таблиц fb_*) и кокпит; его исполняет стандартный
цикл `synchrotronsim.simulation` (`build_session`, `Simulation.run`), а его полная копия ложится в паспорт прогона.

    exp = Experiment.from_fb(3)                       # опыт ОС из БД (fb_experiments + fb_loops + …)
    exp = Experiment.cockpit(lattice=26, optics=True)  # живой кокпит (те же параметры, что раньше в build_session)
    exp.record.append(Watch("obs.ebpm.x:rms"))         # что записывать
    exp.to_json("my_exp.json"); Experiment.from_json("my_exp.json")
    exp.save_bundle("sila26_exp3.zip")                 # опыт + таблицы его машины одним файлом (для Colab/статьи)

Параметры конструктора plant'а хранятся ДОСЛОВНО (`plant.kwargs`): опыт воспроизводит прежнюю сборку бит-в-бит,
без «нормализации», которая могла бы незаметно поменять числа.
"""
from __future__ import annotations

import io
import json
import os
import zipfile
from dataclasses import asdict, dataclass, field
from typing import Any, Dict, List, Optional, Union

from .run import Watch, _jsonable

PLANT_KINDS = ("tbt", "optics", "composite")
LOOP_KINDS = ("integral_svd", "rf_frequency", "optics_tune", "optics_chroma", "optics_adts")
# темп контура: "auto" — раз в окно усреднения plant'а; "cycle" — раз в управляющий цикл (run.steps_per_cycle);
# "optics" — медленный оптический темп (run.optics_period_ticks); число — Гц модельного времени
RATES = ("auto", "cycle", "optics")


@dataclass
class LoopSpec:
    """Контур управления: класс (kind), датчик, актуатор, усиление, порог SVD, уставка, темп, включён ли.
    None — значение по умолчанию класса контура (у оптических контуров датчик/актуатор/уставка зашиты в класс);
    ref "" — явно без уставки (как пустая ячейка fb_loops)."""
    name: str
    kind: str = "integral_svd"
    sensor: Optional[str] = None
    actuator: Optional[str] = None
    gain: Optional[float] = None
    rcond: Optional[float] = None
    ref: Optional[str] = None
    rate: Union[str, float] = "auto"
    enabled: bool = True

    def __post_init__(self):
        if self.kind not in LOOP_KINDS:
            raise ValueError(f"контур {self.name}: вид {self.kind!r}, допустимы {LOOP_KINDS}")
        if isinstance(self.rate, str) and self.rate not in RATES:
            self.rate = float(self.rate)


@dataclass
class PlantSpec:
    """Объект управления: вид (tbt — нелинейное 6D-кольцо turn-by-turn; optics — квазистатика тюна/хроматизма;
    composite — оба: орбита + линеаризованная оптика) и ДОСЛОВНЫЕ аргументы конструктора (без lattice и
    observer'ов). observer: "loops" — наблюдатель из provides() контуров; "default" — встроенный plant'а.
    extra_sensors — датчики сверх контуров (напр. "rf.delta": записать энергию без контура ВЧ)."""
    kind: str = "tbt"
    kwargs: Dict[str, Any] = field(default_factory=dict)
    observer: str = "loops"
    extra_sensors: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.kind not in PLANT_KINDS:
            raise ValueError(f"plant: вид {self.kind!r}, допустимы {PLANT_KINDS}")


@dataclass
class RunLimits:
    """Сколько и как гонять. steps — шагов планировщика (шаг = stride оборотов у tbt), None — без конца (кокпит);
    max_turns — лимит оборотов (кокпит: стоп или новый проход при repeat); steps_per_cycle — управляющий цикл,
    тактов (темп "cycle"); optics_period_ticks — оптический темп, тактов; record_mesh — сетка экрана рекордера."""
    steps: Optional[int] = None
    max_turns: Optional[int] = None
    repeat: bool = False
    steps_per_cycle: int = 9
    optics_period_ticks: int = 90
    record_mesh: int = 121


@dataclass
class Experiment:
    """Опыт: машина (номер решётки; отпечаток — в паспорте прогона), объект управления, контуры, возмущения,
    приводы, подписка на запись, лимиты. noise — закон шума на ДОСКЕ (state.noise, как в режиме 4); шум,
    заданный самому plant'у, — в plant.kwargs["disturbance_spec"] (как в кокпите). drive_profile — профиль
    приводов (slew/лаг/ОДЗ) по типам; None — мгновенные приводы; "default" — DRIVE_PROFILE_DEFAULTS."""
    name: str
    lattice_id: int = 26
    plant: PlantSpec = field(default_factory=PlantSpec)
    loops: List[LoopSpec] = field(default_factory=list)
    noise: Optional[List[Dict[str, Any]]] = None
    drive_profile: Optional[Union[str, Dict[str, Any]]] = None
    record: List[Watch] = field(default_factory=list)
    run: RunLimits = field(default_factory=RunLimits)
    seed: int = 1
    note: str = ""
    source: str = ""

    # ---- из БД: опыт ОС режима 4 -------------------------------------------------------------------
    @classmethod
    def from_fb(cls, fb_exp_id: int, smoke: bool = False) -> "Experiment":
        """Опыт ОС из таблиц fb_* (control.fb_config.load) — та же сборка, что прежний mode4_run.
        smoke — короткий прогон (как `mode4_run smoke`): ≤ 36 шагов (optics ≤ 8), осадка ≤ 120, сетка экрана 81."""
        from control import fb_config
        cfg = fb_config.load(fb_exp_id)
        steps = min(cfg.steps, 8 if cfg.plant == "optics" else 36) if smoke else cfg.steps
        settle = min(cfg.settle, 120) if smoke else cfg.settle
        loops = [LoopSpec(name=l.name, kind="rf_frequency" if l.sensor == "rf.delta" else "integral_svd",
                          sensor=l.sensor, actuator=l.actuator_key, gain=float(l.gain), rcond=float(l.rcond),
                          ref=l.ref or "", rate=l.rate_hz, enabled=True) for l in cfg.enabled_loops()]
        if cfg.plant == "optics":
            plant = PlantSpec("optics", dict(nsl=cfg.nsl, observe=("tune", "chroma", "adts"),
                                             drift_k1=cfg.drift[0], drift_k2=cfg.drift[1], drift_k3=cfg.drift[2],
                                             rho=1.0, seed=cfg.seed), observer="default")
            noise = None
        else:
            rf_act = any(l.actuator.startswith("rf.") for l in loops)
            extra = ["rf.delta"] if any(r.signal == "rf" for r in cfg.enabled_records()) else []
            plant = PlantSpec("tbt", dict(nsl=cfg.nsl, sr=cfg.sr, stride=cfg.stride, settle_turns=settle,
                                          seed=cfg.seed, dist_seed=cfg.seed, n_photon=8,
                                          relative_orbit=cfg.relative_orbit, center_bunch=cfg.center_bunch,
                                          rf_actuator=rf_act), observer="loops", extra_sensors=extra)
            noise = list(cfg.disturbance_spec)
        from .addresses import is_address
        record = []
        for r in cfg.enabled_records():                       # имя рекордера или адрес; прочее (sf) — пропуск
            if _is_recorder(r.signal):
                record.append(Watch(r.signal, cadence=r.cadence, params=r.params, recorder=r.signal))
            elif is_address(r.signal):
                record.append(Watch(r.signal, cadence=_CADENCE.get(r.cadence, r.cadence), params=r.params))
        prof = (cfg.drive_profile or "default") if cfg.realistic_drives else None
        return cls(name=f"fb_exp_{cfg.fb_exp_id}", lattice_id=cfg.lattice_id, plant=plant, loops=loops, noise=noise,
                   drive_profile=prof, record=record,
                   run=RunLimits(steps=steps, record_mesh=81 if smoke else 121), seed=cfg.seed,
                   note=(cfg.note or "") + (" [smoke]" if smoke else ""), source=f"fb_exp_id={cfg.fb_exp_id}")

    # ---- живой кокпит ----------------------------------------------------------------------------
    @classmethod
    def cockpit(cls, lattice: int = 26, nsl: int = 8, settle_turns: int = 250, fofb: bool = True,
                disturbance_spec=None, gap_mm: float = 5.0, steps_per_cycle: int = 9, meas_noise: float = 1e-7,
                n_photon: int = 1, optics: bool = False, optics_drift=(1e-3, 1e-2, 3e-2), optics_order: int = 4,
                optics_period_ticks: int = 90, optics_nsl: int = 2, rf: bool = False, id_pole_nsl: int = 1,
                **plant_kw) -> "Experiment":
        """Опыт живого кокпита — дословно параметры прежнего `client.engine.build_session` (контуры добавлены
        ВЫКЛ; оператор включает панелью). optics=True — composite (орбита + линеаризованная оптика)."""
        loops = []
        if fofb:
            loops += [LoopSpec("fofb_x", sensor="ebpm.x", actuator="hkick", gain=0.5, rcond=1e-2, rate="cycle",
                               enabled=False),
                      LoopSpec("fofb_y", sensor="ebpm.y", actuator="vkick", gain=0.5, rcond=1e-2, rate="cycle",
                               enabled=False)]
        if optics:
            loops += [LoopSpec("tune", kind="optics_tune", rate="optics", enabled=False),
                      LoopSpec("chroma", kind="optics_chroma", rate="optics", enabled=False),
                      LoopSpec("adts", kind="optics_adts", rate="optics", enabled=False)]
        if rf:
            loops += [LoopSpec("rf", kind="rf_frequency", rate="cycle", enabled=False)]
        kw = dict(nsl=nsl, settle_turns=settle_turns, meas_noise=meas_noise, n_photon=n_photon, gap_mm=gap_mm,
                  disturbance_spec=disturbance_spec, rf_actuator=rf, id_pole_nsl=id_pole_nsl, **plant_kw)
        if optics:
            kw.update(drift=optics_drift, order=optics_order, optics_nsl=optics_nsl)
            plant = PlantSpec("composite", kw, observer="loops")
        else:
            plant = PlantSpec("tbt", kw, observer="loops" if rf else "default")
        return cls(name="cockpit", lattice_id=lattice, plant=plant, loops=loops,
                   run=RunLimits(steps_per_cycle=steps_per_cycle, optics_period_ticks=optics_period_ticks),
                   source="cockpit")

    # ---- сериализация ----------------------------------------------------------------------------
    def to_dict(self) -> Dict[str, Any]:
        return _jsonable(asdict(self))

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Experiment":
        d = dict(d)
        p = d.get("plant") or {}
        kw = dict(p.get("kwargs") or {})
        for k in ("observe", "drift"):                       # кортежи, ставшие списками в JSON
            if isinstance(kw.get(k), list):
                kw[k] = tuple(kw[k])
        return cls(name=d["name"], lattice_id=int(d.get("lattice_id", 26)),
                   plant=PlantSpec(p.get("kind", "tbt"), kw, p.get("observer", "loops"),
                                   list(p.get("extra_sensors") or [])),
                   loops=[LoopSpec(**l) for l in d.get("loops") or []], noise=d.get("noise"),
                   drive_profile=d.get("drive_profile"),
                   record=[Watch(**w) for w in d.get("record") or []], run=RunLimits(**(d.get("run") or {})),
                   seed=int(d.get("seed", 1)), note=d.get("note", ""), source=d.get("source", ""))

    def to_json(self, path: str) -> str:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False, indent=1)
        return path

    @classmethod
    def from_json(cls, path: str) -> "Experiment":
        with open(path, encoding="utf-8") as f:
            return cls.from_dict(json.load(f))

    # ---- пакет опыта: опыт + таблицы его машины ----------------------------------------------------
    BUNDLE_TABLES = ("lattice_configs", "element_connectivity", "element_groups", "rf_cavities", "pulsed_kickers",
                     "bpm_config", "screens")                       # строки ЭТОЙ решётки (по lattice_id)
    BUNDLE_LIBRARY = ("hardware_library", "element_transforms", "field_maps")   # строки использованных типов

    def save_bundle(self, path: str, root: Optional[str] = None) -> str:
        """Один zip: experiment.json + манифест (отпечаток машины) + таблицы БД только этой машины (её строки
        решётки и использованные типы устройств). Карты поля не кладутся (киковому ядру не нужны)."""
        import csv
        from .db import Repository
        db = Repository(root)
        m = db.machine(self.lattice_id)
        devs = set(m.devices)
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("experiment.json", json.dumps(self.to_dict(), ensure_ascii=False, indent=1))
            z.writestr("bundle.json", json.dumps({"format": "synchrotronsim.bundle/1", "lattice_id": m.lattice_id,
                                                  "fingerprint": m.fingerprint,
                                                  "label": m.meta.label if m.meta else None},
                                                 ensure_ascii=False, indent=1))
            for name in self.BUNDLE_TABLES + self.BUNDLE_LIBRARY:
                if not db.exists(name):
                    continue
                t = db.table(name)
                if name in self.BUNDLE_TABLES:
                    rows = t.where("lattice_id", str(m.lattice_id))
                else:
                    rows = [r for r in t.rows if t.cell(r, "device_type_id") in devs]
                buf = io.StringIO()
                w = csv.writer(buf, delimiter=";", lineterminator="\n")
                w.writerow(t.header); w.writerows(rows)
                z.writestr(f"tables/{name}.csv", buf.getvalue())
        return path

    @staticmethod
    def load_bundle(path: str, dest: str):
        """Распаковать пакет в папку dest (dest/tables/*.csv) → (Experiment, Repository(dest)). Отпечаток
        собранной из пакета машины сверяется с записанным — пакет либо воспроизводит машину, либо ошибка."""
        from .db import Repository
        os.makedirs(dest, exist_ok=True)
        with zipfile.ZipFile(path) as z:
            for n in z.namelist():
                if n.startswith("tables/") and n.endswith(".csv") and "/" not in n[len("tables/"):]:
                    z.extract(n, dest)
            exp = Experiment.from_dict(json.loads(z.read("experiment.json").decode("utf-8")))
            info = json.loads(z.read("bundle.json").decode("utf-8"))
        db = Repository(dest)
        fp = db.machine(info["lattice_id"]).fingerprint
        if fp != info["fingerprint"]:
            raise RuntimeError(f"пакет {path}: машина не воспроизвелась (отпечаток {fp[:12]} ≠ {info['fingerprint'][:12]})")
        return exp, db


_CADENCE = {"per_turn": "per_step", "snapshot": "final", "post": "final"}   # прежние слова fb_record


def _is_recorder(signal: str) -> bool:
    from control.recorders import REGISTRY
    return signal in REGISTRY
