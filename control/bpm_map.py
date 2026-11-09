# -*- coding: utf-8 -*-
"""Карта имён BPM ↔ индексов в векторах показаний (`state.observables["ebpm.x"]` и т.п.).

Человеку неудобно обращаться к сенсору по номеру — этот модуль даёт ИМЯ `SP<nn>.BPM<mm>` (суперпериод +
порядковый в нём), плюс ground-truth `sensor_id` (= instance_id элемента решётки, как в bpm_readings.csv)
и s-координату/лаб-позицию. Карта строится ПРЯМО из кольца (`orbit_infer.monitor_rows`), поэтому индексы
РОВНО те же, что в векторах орбиты. Экспорт в `workspace/configs/bpm_map_<lattice>.json` — для standalone-
потребителей (отчёт, клиент); `RingModel`/`TbtPlant` держат готовый `BpmMap` в памяти (`self.bpm`).

ЗАПУСК (пере)экспорта JSON:  py -3.14 control/bpm_map.py --lattice 26 [--superperiods 10]
"""
import argparse
import json
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CFG = os.path.join(ROOT, "workspace", "configs")

DEFAULT_SP = {26: 10, 25: 10}                                 # СИЛА — 10 суперпериодов (документировано)


def from_ring(ring, mon, n_superperiods=10):
    """Список записей BPM из готового кольца + индексов мониторов (monitor_rows). Без тяжёлых импортов."""
    mon = [int(m) for m in mon]
    per = max(1, len(mon) // max(1, n_superperiods))
    out = []
    for j, row in enumerate(mon):
        m = ring.mags[row - 1]
        c = [float(x) for x in getattr(m, "center", [0.0, 0.0, 0.0])]
        out.append({"index": j, "name": f"SP{j // per + 1:02d}.BPM{j % per + 1:02d}",
                    "sensor_id": str(getattr(m, "sensor_id", "") or ""),
                    "s": round(float(ring.s_out[row - 1]), 4), "x": round(c[0], 6), "y": round(c[1], 6)})
    return out


def corrector_entries(ring, corr, n_superperiods=10):
    """Записи КОРРЕКТОРОВ (актуаторов орбиты): имя SP<nn>.COR<mm> + dev (напр. und_SH3E) + sensor_id + s + x,y.
    index = порядковый номер корректора (тот же, что в corr.<NNN>.hkick/vkick). Для адресного изменения
    скорости привода (slew) конкретного корректора по имени/sensor_id."""
    corr = [int(c) for c in corr]
    per = max(1, len(corr) // max(1, n_superperiods))
    out = []
    for j, idx in enumerate(corr):
        m = ring.mags[idx]
        c = [float(x) for x in getattr(m, "center", [0.0, 0.0, 0.0])]
        out.append({"index": j, "name": f"SP{j // per + 1:02d}.COR{j % per + 1:02d}",
                    "dev": str(getattr(m, "dev", "") or ""),
                    "sensor_id": str(getattr(m, "sensor_id", "") or ""),
                    "s": round(float(ring.s_out[idx]), 4), "x": round(c[0], 6), "y": round(c[1], 6)})
    return out


class BpmMap:
    """Резолвер имя↔индекс поверх списка записей. Обращение: index(name|sensor_id|int), name(index), named(vec)."""

    def __init__(self, entries):
        self.entries = list(entries)
        self._by_name = {e["name"]: e["index"] for e in self.entries}
        self._by_sid = {e["sensor_id"]: e["index"] for e in self.entries if e["sensor_id"]}
        self._name = {e["index"]: e["name"] for e in self.entries}

    def __len__(self):
        return len(self.entries)

    @property
    def names(self):
        return [e["name"] for e in self.entries]

    def index(self, key):
        """Индекс в векторе: по имени `SP03.BPM02`, по `sensor_id`, или сам int (пропускается)."""
        if isinstance(key, int):
            return key
        if key in self._by_name:
            return self._by_name[key]
        if key in self._by_sid:
            return self._by_sid[key]
        raise KeyError(f"BPM не найден: {key!r} (есть {len(self)} шт., напр. {self.names[:3]})")

    def name(self, index):
        return self._name[int(index)]

    def entry(self, key):
        return self.entries[self.index(key)]

    def named(self, vec):
        """Вектор показаний → {имя: значение} — для человека/отчёта/клиента."""
        return {e["name"]: float(vec[e["index"]]) for e in self.entries}

    @classmethod
    def load(cls, lattice):
        """Прочитать экспортированный JSON (для тех, у кого нет живого кольца — отчёт/клиент)."""
        with open(os.path.join(CFG, f"bpm_map_{lattice}.json"), encoding="utf-8") as f:
            return cls(json.load(f)["bpms"])


def build_and_save(lattice, n_superperiods=None, nsl=8):
    import facility_report as fr
    import orbit_infer as oi
    inf = dict(fr.info_of(lattice)); inf["_id"] = lattice
    ring = fr._ring(inf, nsl=nsl)
    entries = from_ring(ring, oi.monitor_rows(ring), n_superperiods or DEFAULT_SP.get(lattice, 10))
    os.makedirs(CFG, exist_ok=True)
    path = os.path.join(CFG, f"bpm_map_{lattice}.json")
    with open(path, "w", encoding="utf-8") as f:
        json.dump({"lattice": lattice, "n_superperiods": n_superperiods or DEFAULT_SP.get(lattice, 10),
                   "n_bpm": len(entries), "bpms": entries}, f, ensure_ascii=False, indent=1)
    try:
        ring.close()
    except Exception:
        pass
    return path, len(entries)


def main():
    os.environ.setdefault("SS_SECTOR", "1")
    sys.stdout.reconfigure(encoding="utf-8")
    sys.path.insert(0, ROOT)
    ap = argparse.ArgumentParser(description="экспорт карты имён BPM ↔ индексов")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--superperiods", type=int, default=None)
    a = ap.parse_args()
    path, n = build_and_save(a.lattice, a.superperiods)
    print(f"[OK] {os.path.relpath(path, ROOT).replace(os.sep, '/')} — {n} BPM")


if __name__ == "__main__":
    main()
