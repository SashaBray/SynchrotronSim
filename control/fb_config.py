# -*- coding: utf-8 -*-
"""Конфиг эксперимента ОС (режим 4) ИЗ БД ПО КЛЮЧУ `fb_exp_id`.

Мастер-строка `fb_experiments.csv` + сателлиты по тому же ключу (`fb_loops`, `fb_disturbances`,
`fb_record`) — тот же паттерн, что `experiments.csv`+`radiation_params.csv`
(`matrix_tracker.Settings.from_experiment`+`load_radiation`). Спек возмущений берём из общего
`disturbance.load_spec()`, отфильтрованного по именам из `fb_disturbances`. Только чтение БД.

Схема таблиц (разделитель «;», точка — десятичная):
  fb_experiments: fb_exp_id;lattice_id;steps;stride;settle;seed;sr;nsl;note;status;plant;drift_k1/k2/k3;realistic_drives
  fb_loops:       fb_exp_id;loop;sensor;actuator_key;gain;rcond;rate_hz;ref;enabled
  fb_disturbances:fb_exp_id;name                        (джойн на disturbance_spec.csv по name)
  fb_record:      fb_exp_id;signal;params;cadence;enabled   (ПОДПИСКА: за чем следить/что писать)
  fb_drives:      fb_exp_id;type;slew;lag_tau;odz_lo;odz_hi;res   (богатый профиль приводов, opt-in; пусто=дефолт)
"""
import csv
import os
from dataclasses import dataclass, field

import disturbance as dist

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
TB = os.path.join(ROOT, "workspace", "tables")


def _read(name):
    path = os.path.join(TB, name)
    if not os.path.exists(path):
        return []
    with open(path, encoding="utf-8-sig") as f:
        return list(csv.DictReader(f, delimiter=";"))


def _rows_for(name, key):
    return [r for r in _read(name) if str((r.get("fb_exp_id") or "")).strip() == str(key)]


def _b(v):
    return str(v).strip().lower() in ("1", "true", "yes", "on", "да")


def _drive_profile(key):
    """Богатый профиль приводов из fb_drives.csv по ключу: {тип: {slew?, lag_tau?, odz:(lo,hi)?, res?}}."""
    prof = {}
    for r in _rows_for("fb_drives.csv", key):
        t = (r.get("type") or "").strip()
        if not t:
            continue
        spec = {}
        for k in ("slew", "lag_tau", "res"):
            v = (r.get(k) or "").strip()
            if v:
                try:
                    spec[k] = float(v)
                except ValueError:
                    pass
        lo = (r.get("odz_lo") or "").strip(); hi = (r.get("odz_hi") or "").strip()
        if lo and hi:
            try:
                spec["odz"] = (float(lo), float(hi))
            except ValueError:
                pass
        if spec:
            prof[t] = spec
    return prof


@dataclass
class Loop:
    name: str
    sensor: str
    actuator_key: str
    gain: float
    rcond: float
    rate_hz: str          # число или "auto" (= один апдейт на окно усреднения)
    ref: str              # имя уставки или "" (= None, следим за нулём)
    enabled: bool


@dataclass
class Record:
    signal: str           # orbit_rms / light / phase / bpm / screen / tune / chroma / sf
    params: str           # напр. screen_id для screen, channel=b1 для sf
    cadence: str          # per_turn / per_window / snapshot / final / post
    enabled: bool


@dataclass
class FbConfig:
    fb_exp_id: int
    lattice_id: int
    steps: int            # число шагов планировщика (один шаг = stride оборотов)
    stride: int
    settle: int
    seed: int
    sr: bool
    nsl: int
    note: str
    loops: list = field(default_factory=list)
    disturbance_spec: list = field(default_factory=list)   # row-dict'ы для disturbance.Disturbance
    records: list = field(default_factory=list)
    plant: str = "tbt"                                     # tbt (орбита/TbtPlant) | optics (OpticsPlant)
    drift: tuple = (0.0, 0.0, 0.0)                         # дрейф градиентов (drift_k1/k2/k3) для optics-plant
    realistic_drives: bool = False                         # opt-in: конечная скорость+лаг приводов ВСЕХ магнитов
    drive_profile: dict = field(default_factory=dict)      # богатый профиль приводов из fb_drives (пусто = дефолт)

    def enabled_records(self):
        return [r for r in self.records if r.enabled]

    def enabled_loops(self):
        return [l for l in self.loops if l.enabled]


def load(fb_exp_id):
    master = _rows_for("fb_experiments.csv", fb_exp_id)
    if not master:
        raise SystemExit(f"fb_config: нет эксперимента fb_exp_id={fb_exp_id} в fb_experiments.csv")
    m = master[0]
    loops = [Loop(r["loop"].strip(), r["sensor"].strip(), r["actuator_key"].strip(),
                  float(r["gain"]), float(r["rcond"]), (r.get("rate_hz") or "auto").strip(),
                  (r.get("ref") or "").strip(), _b(r.get("enabled", "1")))
             for r in _rows_for("fb_loops.csv", fb_exp_id)]
    names = {(r.get("name") or "").strip() for r in _rows_for("fb_disturbances.csv", fb_exp_id)}
    spec = [c for c in dist.load_spec() if c["name"] in names]   # общий спек, отфильтрованный по подписке
    records = [Record(r["signal"].strip(), (r.get("params") or "").strip(),
                      (r.get("cadence") or "per_window").strip(), _b(r.get("enabled", "1")))
               for r in _rows_for("fb_record.csv", fb_exp_id)]
    def _f(key):
        try:
            return float(m.get(key) or 0.0)
        except ValueError:
            return 0.0
    return FbConfig(
        fb_exp_id=int(m["fb_exp_id"]), lattice_id=int(m["lattice_id"]), steps=int(m["steps"]),
        stride=int(m["stride"]), settle=int(m["settle"]), seed=int(m["seed"]), sr=_b(m["sr"]),
        nsl=int(m["nsl"]), note=(m.get("note") or "").strip(),
        loops=loops, disturbance_spec=spec, records=records,
        plant=(m.get("plant") or "tbt").strip() or "tbt",
        drift=(_f("drift_k1"), _f("drift_k2"), _f("drift_k3")),
        realistic_drives=_b(m.get("realistic_drives", "0")),
        drive_profile=_drive_profile(fb_exp_id))
