# -*- coding: utf-8 -*-
"""Живая 3D-сцена установки — геометрия из РАБОТАЮЩЕГО кольца (plant), не из БД/STL.

Две части, как и у батч-визуализации (геометрия отдельно от рендера):
  • `build_geometry(plant, state)` → СЕРИАЛИЗУЕМЫЙ «узел геометрии» (числа): центры элементов в X-Y
    (планарное кольцо, z=0), ряды мониторов, якоря ондуляторов + их зазор, обвод кольца. Это и есть
    то, что поедет по шву Link в сетевом режиме (по проводу — числа, не пиксели).
  • `SceneModel(bundle)` — ЧИСТАЯ математика актёров (положения челюстей ондулятора от зазора,
    стартовая плоскость, направление пути). Без PyVista/Qt → тестируется headless.
  • `RingScene` — тонкий вид на `pyvistaqt.QtInteractor`: строит актёров из SceneModel и ДВИГАЕТ их
    по снимку (зазор ондулятора → расстояние между челюстями). PyVista/Qt импортируются ЛЕНИВО.

Ядро не трогаем: читаем только `plant.ring.mags[*].center` (геометрия решётки 26, совпадает с батчем) и
живой зазор из `state.undulators`; динамика приходит снимком подписанных адресов.
"""
from __future__ import annotations

import os
from collections import defaultdict
from typing import Dict, List, Optional

import numpy as np


# ------------------------------------------------------------------ геометрия (сериализуемая)
# элементы САМОГО ондулятора (полюса + концевые) — не STL-тела, рисуются гребёнкой ондулятора; не адресуем
POLE_DEVS = {"und_UPM", "und_UPP", "und_UE0", "und_UE1"}
# устройства, СКРЫТЫЕ в сцене по просьбе (чистый визуал; ядро и якоря ондуляторов не затрагиваются)
HIDDEN_DEVS = {"und_UFC"}                 # прежний фокус-элемент ИД (в текущих решётках его нет)
# типы, у которых рисуем ТЕЛО-блок (остальное — дрейфы/маркеры/мониторы — тел не имеют)
DRAW_KINDS = {"dipole", "dq", "quad", "sext", "oct", "corrector", "kicker"}
KIND_RU = {"dipole": "диполь", "dq": "диполь-квадруполь", "quad": "квадруполь", "sext": "секступоль",
           "oct": "октуполь", "corrector": "корректор", "kicker": "корректор", "undulator": "ондулятор"}
NAV_MERGE = 1e-3                          # м: строки в одной точке (обмотки, DQ I1..I3) — один шаг навигации
NAV_LOCK = 0.05                           # м: камера «на элементе», если фокус от него ближе
NAV_SNAP = 2.0                            # м: иначе берём ближайший к фокусу элемент, если он не дальше


def _mag_names(plant):
    nm = getattr(getattr(plant, "magnets", None), "names", None)
    if nm is None:
        return []
    return list(nm() if callable(nm) else nm)


def build_geometry(plant, state=None, *, decimate: int = 1) -> dict:
    """Собрать узел геометрии из работающего plant'а (+доска для живых зазоров). Всё — числа/строки.
    Без ring_cpp (фейк-plant) вернёт вырожденный узел: пустое кольцо + ондуляторы с доски.

    `named` — адресуемые элементы (имя из plant.magnets ↔ индекс ring.mags, кроме пол-модулей): для
    пикинга (hover/ПКМ) и адреса `element.<name>.*`. Центры в X-Y (планарное кольцо, z=0)."""
    g: dict = {"lattice": int(getattr(plant, "lattice", 0)),
               "circumference": float(getattr(getattr(plant, "ring", None), "circumference", 0.0) or 0.0),
               "centers": [], "mon_xy": [], "named": [], "magnets": [], "undulators": {}, "screens": [],
               "sensors": []}
    ring = getattr(plant, "ring", None)
    if ring is not None and getattr(ring, "mags", None) is not None:
        mags = ring.mags
        step = max(1, int(decimate))
        g["centers"] = [[float(m.center[0]), float(m.center[1])] for m in mags[::step]]
        g["mon_xy"] = [[float(m.center[0]), float(m.center[1])]
                       for m in mags if str(getattr(m, "kind", "")) == "monitor"]
        names = _mag_names(plant)
        corr_pos = {int(ci): k for k, ci in enumerate(getattr(plant, "corr", []) or [])}  # ring-idx → № актуатора
        # id семейных ГРАДИЕНТНЫХ актуаторов (контуры оптики): quad.<fam>.dk / sext.<fam>.dk2 / oct.<fam>.dk3.
        # Привязываемся по имени семейства `_family(dev)` — оно совпадает в обоих кольцах (C++ орбиты и
        # аналитич. оптики), хотя индексы разные. Только если такой актуатор есть в схеме plant'а (optics=True),
        # иначе линк мёртв. Так магнит показывает ЖИВУЮ Δk, которую реализует контур (двигается при вкл/выкл ОС).
        try:
            _sids = set(plant.schema().ids)
        except Exception:
            _sids = set()
        _FAM = {"quad": ("quad", "dk"), "sext": ("sext", "dk2"), "oct": ("oct", "dk3")}
        _fam_fn = None
        if _sids:
            try:
                import facility_report as _frm
                _fam_fn = _frm._family
            except Exception:
                _fam_fn = None
        _midx = 0                                            # порядковый номер сенсора = индекс в наблюдаемой ebpm
        for i, m in enumerate(mags):
            dev = str(getattr(m, "dev", "") or "")
            if dev in POLE_DEVS or dev in HIDDEN_DEVS:       # полюса рисуются гребёнкой; скрытые — не рисуем
                continue
            nm = names[i] if i < len(names) else f"mag.{i}"
            kind = str(getattr(m, "kind", "") or "")
            cx, cy = float(m.center[0]), float(m.center[1])
            entry = {"name": nm, "xy": [cx, cy], "kind": kind, "dev": dev, "s_idx": i}   # s_idx — порядок по пучку
            if kind in ("corrector", "kicker") and i in corr_pos:   # id реализованных H/V-киков корректора
                j = corr_pos[i]
                entry["aids"] = [f"corr.{j:03d}.hkick", f"corr.{j:03d}.vkick"]
            elif kind in _FAM and _fam_fn is not None:              # id реализованной Δk семейства (контур оптики)
                pfx, sfx = _FAM[kind]
                fid = f"{pfx}.{_fam_fn(dev)}.{sfx}"
                if fid in _sids:
                    entry["aids"] = [fid]
            g["named"].append(entry)
            if kind == "monitor":                               # сенсор BPM: экран + idx = порядок = индекс ebpm
                ax = getattr(m, "axis", None)
                a2s = [float(ax[0]), float(ax[1])] if ax is not None else _axis_xy(mags, i)
                g["sensors"].append({"name": nm, "xy": [cx, cy],
                                     "z": float(m.center[2]) if len(m.center) > 2 else 0.0,
                                     "sid": str(getattr(m, "sensor_id", "") or ""),
                                     "axis": a2s, "idx": _midx,
                                     "meas_noise": float(getattr(plant, "meas_noise", 0.0) or 0.0)})
                _midx += 1
            if kind in DRAW_KINDS:                               # тело магнита: STL-инстанс + запас (блок)
                ax = getattr(m, "axis", None)
                uh = getattr(m, "uhat", None)
                a2 = [float(ax[0]), float(ax[1])] if ax is not None else _axis_xy(mags, i)
                u2 = [float(uh[0]), float(uh[1])] if uh is not None else [-a2[1], a2[0]]
                cz = float(m.center[2]) if len(m.center) > 2 else 0.0
                t = [float(getattr(m, "yaw", 0.0) or 0.0), float(getattr(m, "pitch", 0.0) or 0.0),
                     float(getattr(m, "roll", 0.0) or 0.0), cx, cy, cz]   # трансформ для STL-инстанса
                g["magnets"].append({"name": nm, "dev": dev,
                                     "c": [cx, cy], "ax": a2, "ux": u2, "t": t,
                                     "L": float(getattr(m, "L", 0.0) or 0.0), "kind": kind})
        # якоря ондуляторов: group_id → якорь (выходной полюс) + осевое направление + ЦЕНТРЫ ПОЛЮСОВ (для отрисовки
        # самого ондулятора периодическим рядом). Полюса — из plant._und_poles; старый путь (якорь-UFC после полюсов)
        # — сбор назад по s от якоря, если plant их не отдаёт.
        poles_of = getattr(plant, "_und_poles", {}) or {}
        for gid, idx in getattr(plant, "_und_anchor", {}).items():
            idx = int(idx)
            m = mags[idx]
            axis = _axis_xy(mags, idx)
            cells = []                                       # ячейки гребёнки: [x, y, is_end] (концевые — длиннее)
            if gid in poles_of:
                seq = [int(j) for j in poles_of[gid]]
            else:
                seq, j = [], idx - 1
                while j >= 0 and str(getattr(mags[j], "dev", "") or "") in ("und_UPM", "und_UPP", "und_UE0", "und_UE1"):
                    seq.append(j); j -= 1
                seq.reverse()
            for j in seq:
                dv = str(getattr(mags[j], "dev", "") or "")
                cells.append([float(mags[j].center[0]), float(mags[j].center[1]),
                              1 if dv in ("und_UE0", "und_UE1") else 0])
            cen = (np.mean(np.asarray(cells, float)[:, :2], axis=0) if cells
                   else np.array([float(m.center[0]), float(m.center[1])]))
            g["undulators"][gid] = {"xy": [float(m.center[0]), float(m.center[1])], "axis": axis,
                                    "gap_mm": _gap_of(state, gid, default=plant.gap_mm), "poles": cells,
                                    "s_idx": int(seq[0]) if seq else idx,       # порядок по пучку (первый полюс)
                                    "center": [float(cen[0]), float(cen[1])]}   # центр гребёнки — для камеры
    # доска могла знать ондуляторы, которых нет среди якорей (фейк-plant): добавить с нулём координат
    for gid in (getattr(state, "undulators", {}) or {}):
        g["undulators"].setdefault(gid, {"xy": [0.0, 0.0], "axis": [1.0, 0.0],
                                         "gap_mm": _gap_of(state, gid, default=5.0)})
    # экраны у ВСЕХ ондуляторов (даже без записи в screens.csv — тогда дефолтная дистанция/размер).
    # Дистанция/размер из таблицы по source_ref==gid, координаты — от ЖИВОГО якоря (совпадает с решёткой).
    scr = _read_screens()
    for gid, u in g["undulators"].items():
        sc = scr.get(gid, {})
        dist = float(sc.get("dist") or 0.0) or SCREEN_DEFAULT_DIST
        ax = u.get("axis", [1.0, 0.0])
        g["screens"].append({"gid": gid, "sid": sc.get("sid", "") or f"SCR_{gid}",
                             "src": [u["xy"][0], u["xy"][1]],
                             "center": [u["xy"][0] + dist * ax[0], u["xy"][1] + dist * ax[1]],
                             "axis": [ax[0], ax[1]],
                             "w": float(sc.get("w", 0) or 0), "h": float(sc.get("h", 0) or 0)})
    # СТАРТ решётки/пучка s=0 — в глобальном начале координат (раскладка elegant). Первый элемент в ring.mags
    # (BPM) стоит ПОСЛЕ входного дрейфа DR_01, поэтому стартплоскость кладём на s=0, а не на mags[0].
    if len(g["centers"]) >= 2:
        c0 = np.asarray(g["centers"][0], float)
        c1 = np.asarray(g["centers"][1], float)
        d = c1 - c0
        n = float(np.hypot(d[0], d[1])) or 1.0
        g["beam_start"] = {"pos": [0.0, 0.0], "dir": [float(d[0] / n), float(d[1] / n)]}
    return g


def _read_screens(root: Optional[str] = None) -> dict:
    """Экраны из workspace/tables/screens.csv: {source_ref(gid): {dist(м), w,h(м), sid}}. {} если файла нет.
    source_ref — ондулятор-источник; привязываемся по нему, а координаты берём от живого якоря (решётко-незав.)."""
    import csv
    p = os.path.join(root or _project_root(), "workspace", "tables", "screens.csv")
    out: dict = {}
    if not os.path.exists(p):
        return out
    try:
        with open(p, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f, delimiter=";"):
                src = (r.get("source_ref") or "").strip()
                if not src:
                    continue
                out[src] = {"dist": float(r.get("dist_m") or 0) or 0.0,
                            "w": float(r.get("width_mm") or 0) / 1000.0,
                            "h": float(r.get("height_mm") or 0) / 1000.0,
                            "sid": (r.get("screen_id") or "").strip()}
    except (OSError, ValueError, KeyError):
        pass
    return out


def _gap_of(state, gid, default):
    u = (getattr(state, "undulators", {}) or {}).get(gid) if state is not None else None
    return float(u.gap_mm) if u is not None else float(default)


def _axis_xy(mags, idx):
    """Единичное направление пути в X-Y у элемента idx (из соседних центров; запасной — ось X)."""
    j = min(idx + 1, len(mags) - 1)
    a = np.array([float(mags[idx].center[0]), float(mags[idx].center[1])])
    b = np.array([float(mags[j].center[0]), float(mags[j].center[1])])
    d = b - a
    n = float(np.hypot(*d))
    return [float(d[0] / n), float(d[1] / n)] if n > 1e-9 else [1.0, 0.0]


# ------------------------------------------------------------------ чистая модель актёров
# тела магнитов: цвет и относит. поперечная ширина по типу (палитра как у остальных визуализаторов)
KIND_COLOR = {"dipole": "#e15759", "dq": "#f28e2b", "quad": "#4e79a7", "sext": "#59a14f",
              "oct": "#b07aa1", "kicker": "#9c755f", "corrector": "#2b9099"}
DEFAULT_KIND_COLOR = "#8a8a8a"
BODY_OPACITY = 0.5      # полупрозрачность тел магнитов (как в visualize.py — видно структуру/пучок)
EDGE_COLOR = "#3a3a3a"  # цвет рёбер ярм (show_edges=True — форма читается, не плоское пятно)
KIND_WIDTH = {"dipole": 1.0, "dq": 1.0, "quad": 0.85, "sext": 0.7, "oct": 0.6,
              "kicker": 0.8, "corrector": 0.7}
MAG_LEN_EXAGG = 1.0     # масштаб длины тела магнита (1 = реальная L_eff; срезы тел тайлятся встык в тело)
MAG_MIN_LEN = 0.12      # минимальная видимая длина блока (для элементов с крошечной L)
UND_HALF_L = 0.3        # полудлина челюсти ондулятора вдоль оси (бар ~0.6 м — адекватно реальному ~0.5 м)
UND_HALF_W = 0.1        # полуширина челюсти поперёк
UND_HALF_H = 0.03       # полувысота челюсти (тонкий бар)
STL_TRI_CAP = 20_000_000  # практически БЕЗ прореживания (полное качество STL); предохранитель от патологий
PICK_THRESH = 2.0       # порог пикинга (м): клик дальше этого от магнита — не открываем инфо (пустой клик)
SCREEN_MIN_SIZE = 0.6   # видимый размер экрана в сцене (реальный ~50 мм — спек на 12 м; увеличен для глаза)
SCREEN_DEFAULT_DIST = 12.0  # дистанция экрана от ондулятора, если нет записи в screens.csv (у ВСЕХ ондуляторов)
HOVER_SIZE = 1.0        # размер жёлтой рамки-подсветки магнита под курсором (м)
START_PLANE_SIZE = 1.2  # размер стартовой плоскости (м) — небольшой, не во всё кольцо
SENSOR_THICK = 0.02     # толщина экрана-сенсора (BPM) вдоль пучка, м (тоньше ×2)
SENSOR_SIZE = 0.216     # поперечный размер экрана-сенсора (апертура), м (+20%: шире и выше)
JAW_VIS_RATE = 10.0     # скорость ПЛАВНОЙ анимации зазора челюстей, мм/с стенного времени (правдоподобный привод)


class SceneModel:
    """Чистая геометрия актёров из узла `build_geometry`. Без PyVista/Qt — для рендера и для тестов.

    Сцена в РЕАЛЬНОМ масштабе: тела магнитов реальной длины (срезы тайлятся встык в сплошное тело),
    челюсти ондулятора разводим по Z на ±`jaw_offset(gap)` — небольшое абсолютное смещение у пучка (не
    улетает), растёт линейно с зазором. `gap_scale` — безразмерный множитель-утрирование разведения
    челюстей для наглядности (чистый визуал, на физику не влияет)."""

    def __init__(self, bundle: dict, gap_scale: float = 1.0):
        self.bundle = bundle or {}
        self.gap_scale = float(gap_scale)      # множитель-утрирование визуального разведения челюстей
        self.centers = np.asarray(self.bundle.get("centers", []), float).reshape(-1, 2)
        self.mon_xy = np.asarray(self.bundle.get("mon_xy", []), float).reshape(-1, 2)
        self.undulators = dict(self.bundle.get("undulators", {}))
        self.named = list(self.bundle.get("named", []))
        self._named_xy = np.asarray([n["xy"] for n in self.named], float).reshape(-1, 2) if self.named \
            else np.zeros((0, 2))
        self.magnets = list(self.bundle.get("magnets", []))
        self.screens = list(self.bundle.get("screens", []))
        self.sensors = list(self.bundle.get("sensors", []))
        self._sensor_xy = (np.asarray([s["xy"] for s in self.sensors], float).reshape(-1, 2)
                           if self.sensors else np.zeros((0, 2)))
        pk = DRAW_KINDS | {"monitor"}                   # пикаем только «важные» элементы (не дрейфы/маркеры)
        self._pick_mask = (np.array([str(n.get("kind", "")) in pk for n in self.named], bool)
                           if self.named else np.zeros(0, bool))
        self.nav = self._build_nav()
        self._nav_xy = (np.asarray([it["xy"] for it in self.nav], float).reshape(-1, 2) if self.nav
                        else np.zeros((0, 2)))

    # ---- навигация камеры по элементам (Ctrl/Shift + ←/→) ----
    def _build_nav(self) -> List[dict]:
        """Элементы для перелёта камеры в ПОРЯДКЕ ДВИЖЕНИЯ ПУЧКА: тела магнитов (DRAW_KINDS) + ондуляторы целиком
        (центр гребёнки полюсов). Строки в одной точке (обмотки корректора, DQ I1/I2/I3, пары секступолей)
        склеиваются в один шаг — имя и тип по первой. type — тип устройства (dev; у ондулятора "undulator")."""
        items = []
        for k, n in enumerate(self.named):
            kind = str(n.get("kind", ""))
            if kind not in DRAW_KINDS:
                continue
            items.append({"name": str(n.get("name", "")), "xy": [float(n["xy"][0]), float(n["xy"][1])],
                          "kind": kind, "type": str(n.get("dev", "") or kind), "s": float(n.get("s_idx", k))})
        for gid, u in self.undulators.items():
            if "s_idx" not in u:                       # нет положения по пучку (фейк-доска) — не в навигации
                continue
            c = u.get("center") or u.get("xy") or [0.0, 0.0]
            items.append({"name": str(gid), "xy": [float(c[0]), float(c[1])], "kind": "undulator",
                          "type": "undulator", "s": float(u["s_idx"]) - 0.5})
        items.sort(key=lambda it: it["s"])
        out = []
        for it in items:
            if out and np.hypot(it["xy"][0] - out[-1]["xy"][0], it["xy"][1] - out[-1]["xy"][1]) < NAV_MERGE:
                continue
            out.append(it)
        return out

    def nav_nearest(self, x: float, y: float):
        """(индекс, расстояние) ближайшего элемента навигации к (x, y); None — навигации нет."""
        if self._nav_xy.size == 0:
            return None
        d = np.hypot(self._nav_xy[:, 0] - float(x), self._nav_xy[:, 1] - float(y))
        j = int(np.argmin(d))
        return j, float(d[j])

    def nav_step(self, i: int, step: int, same_type: bool = False) -> int:
        """Соседний элемент от i: step=+1 — по движению пучка, −1 — навстречу; кольцо замкнуто (по кругу).
        same_type — следующий элемент того же типа устройства (напр. следующий ондулятор)."""
        n = len(self.nav)
        if n == 0:
            return i
        t = self.nav[i]["type"]
        j = i
        for _ in range(n):
            j = (j + (1 if step > 0 else -1)) % n
            if not same_type or self.nav[j]["type"] == t:
                return j
        return i

    def nav_label(self, i: int) -> str:
        """Подпись для статус-строки: имя, вид и номер среди элементов того же типа."""
        it = self.nav[i]
        same = [k for k, x in enumerate(self.nav) if x["type"] == it["type"]]
        return (f"{KIND_RU.get(it['kind'], it['kind'])} {it['name']} — {same.index(i) + 1} из {len(same)} "
                f"({it['type']})" if it["kind"] != "undulator" else
                f"ондулятор {it['name']} — {same.index(i) + 1} из {len(same)}")

    def ring_polyline(self) -> np.ndarray:
        """Точки обвода кольца (x,y,0), в порядке s. Пусто → (0,3)-массив."""
        if self.centers.size == 0:
            return np.zeros((0, 3))
        return np.column_stack([self.centers, np.zeros(len(self.centers))])

    def monitor_points(self) -> np.ndarray:
        if self.mon_xy.size == 0:
            return np.zeros((0, 3))
        return np.column_stack([self.mon_xy, np.zeros(len(self.mon_xy))])

    def start_plane(self):
        """(центр, нормаль) стартовой плоскости. s=0 решётки/пучка (bundle['beam_start'] = глоб. начало
        координат); иначе — у первого нарисованного элемента. None, если кольца нет."""
        bs = self.bundle.get("beam_start")
        if bs is not None:
            p, d = bs["pos"], bs["dir"]
            return ([float(p[0]), float(p[1]), 0.0], [float(d[0]), float(d[1]), 0.0])
        if len(self.centers) < 2:
            return None
        c = self.centers[0]
        d = self.centers[1] - self.centers[0]
        n = float(np.hypot(*d)) or 1.0
        return ([float(c[0]), float(c[1]), 0.0], [float(d[0] / n), float(d[1] / n), 0.0])

    def jaw_offset(self, gap_mm: float) -> float:
        """Z-смещение ЦЕНТРА ряда зубьев (м). Зубья построены с полувысотой UND_HALF_H вокруг своего центра,
        поэтому центр = полувысота зуба + РЕАЛЬНЫЙ полузазор (gap/2, мм→м) × gap_scale. Так ВНУТРЕННИЕ грани
        рядов отстоят ровно на зазор: при gap_scale=1 пустота между 3D-моделями ТОЧНО равна заданному зазору."""
        return UND_HALF_H + 0.5 * float(gap_mm) * 1e-3 * self.gap_scale

    def undulator_jaws(self, gid: str, gap_mm: float):
        """Центры верхней и нижней челюстей ондулятора (x,y,z) при зазоре gap_mm. Разведены по Z на
        ±jaw_offset вокруг оси пучка. Больше зазор → больше |z|-смещение (монотонно)."""
        u = self.undulators.get(gid, {"xy": [0.0, 0.0]})
        x, y = float(u["xy"][0]), float(u["xy"][1])
        off = self.jaw_offset(gap_mm)
        return ([x, y, off], [x, y, -off])

    def undulator_list(self) -> List[str]:
        return list(self.undulators)

    def named_points(self) -> np.ndarray:
        """Точки адресуемых элементов (x,y,0) для облака пикинга; (0,3), если их нет."""
        if self._named_xy.size == 0:
            return np.zeros((0, 3))
        return np.column_stack([self._named_xy, np.zeros(len(self._named_xy))])

    def magnet_boxes(self):
        """Ориентированные блоки-тела магнитов: список (center3, axis3, uhat3, hl, hw, hh, kind). Длина
        вдоль оси пучка (утрирована, чтобы читалось на кольце), ширина/высота — по типу. Для рендера тел."""
        out = []
        for mg in self.magnets:
            c = np.array([float(mg["c"][0]), float(mg["c"][1]), 0.0])
            a = np.array([float(mg["ax"][0]), float(mg["ax"][1]), 0.0])
            a = a / (float(np.hypot(a[0], a[1])) or 1.0)
            u = np.array([float(mg["ux"][0]), float(mg["ux"][1]), 0.0])
            u = u / (float(np.hypot(u[0], u[1])) or 1.0)
            kind = str(mg.get("kind", ""))
            lvis = max(float(mg.get("L", 0.0)) * MAG_LEN_EXAGG, MAG_MIN_LEN)
            out.append((c, a, u, 0.5 * lvis, 0.25 * lvis * KIND_WIDTH.get(kind, 0.8), 0.16 * lvis, kind))
        return out

    def nearest_named(self, x: float, y: float, max_dist: Optional[float] = None, tag: Optional[str] = None):
        """Ближайший «важный» элемент (магнит/монитор, НЕ дрейф) к (x,y): (index, name, kind, dist).
        None, если пусто или ближайший дальше max_dist (клик мимо магнитов). `tag` — dev/kind кликнутого
        ТЕЛА: сужает кандидатов до экземпляров этого типа (точный пик, а не ближайший центр соседа)."""
        if self._named_xy.size == 0:
            return None
        base = (self._pick_mask if (self._pick_mask.size and self._pick_mask.any())
                else np.ones(len(self.named), bool))
        if tag:                                            # только экземпляры кликнутого тела (dev или kind)
            sel = np.array([(self.named[k].get("dev", "") == tag or self.named[k].get("kind", "") == tag)
                            for k in range(len(self.named))], bool) & base
            if sel.any():
                base = sel
        idxs = np.where(base)[0]
        d = np.hypot(self._named_xy[idxs, 0] - float(x), self._named_xy[idxs, 1] - float(y))
        j = int(np.argmin(d))
        i = int(idxs[j])
        dist = float(d[j])
        if max_dist is not None and dist > max_dist:
            return None
        n = self.named[i]
        return (i, n["name"], n.get("kind", ""), dist)

    def nearest_sensor(self, x: float, y: float):
        """Ближайший сенсор (BPM) к (x,y): dict из bundle['sensors'] (name/xy/z/sid/idx) или None."""
        if self._sensor_xy.size == 0:
            return None
        d = np.hypot(self._sensor_xy[:, 0] - float(x), self._sensor_xy[:, 1] - float(y))
        return self.sensors[int(np.argmin(d))]

    def bbox(self):
        """Габариты кольца в X-Y ([xmin,xmax,ymin,ymax]); None, если пусто. Для авто-камеры."""
        if self.centers.size == 0:
            return None
        return [float(self.centers[:, 0].min()), float(self.centers[:, 0].max()),
                float(self.centers[:, 1].min()), float(self.centers[:, 1].max())]


# ------------------------------------------------------------------ настоящие STL-тела магнитов
def _project_root() -> str:
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _pre_rot(row):
    """R_map (3×3) из rot_x/y/z_deg строки hardware_library (как в viz_html), или None. Геометрию магнита
    вертим так же, как C++ FieldMap вертит поле (столбцы rot_*_deg)."""
    import visualize as vz

    def _d(k):
        try:
            return float(np.radians(float(row.get(k, 0) or 0)))
        except (TypeError, ValueError):
            return 0.0
    ang = (_d("rot_x_deg"), _d("rot_y_deg"), _d("rot_z_deg"))
    return vz.rotation_matrix_zyx(*ang) if any(ang) else None


def build_magnet_stl_meshes(magnets, *, tri_cap: int = STL_TRI_CAP, root: Optional[str] = None) -> dict:
    """НАСТОЯЩИЕ тела магнитов: {dev: (combined_mesh, kind)}. Живые трансформы (yaw/pitch/roll/xyz) из узла,
    STL+поворот по устройству из `hardware_library` (справочник геометрии, не физика), инстансинг через
    `visualize.build_combined_mesh` (тот же конвейер, что у viz_html). Прореживаем под бюджет треугольников,
    чтобы живой рендер не лагал. Возврат {} → STL/hardware недоступны (рендер откатывается на блоки).
    Чистая (pyvista, без Qt) — для RingScene и для оффскрин-проверки."""
    if not magnets:
        return {}
    try:
        import pyvista as pv
        import visualize as vz
    except Exception:
        return {}
    ws = os.path.join(root or _project_root(), "workspace")
    hwp = os.path.join(ws, "tables", "hardware_library.csv")
    if not os.path.exists(hwp):
        return {}
    try:
        import csv
        with open(hwp, encoding="utf-8") as f:
            hw = {r["device_type_id"]: r for r in csv.DictReader(f, delimiter=";")}
    except (OSError, KeyError, ValueError):
        return {}
    bydev: Dict[str, list] = defaultdict(list)                   # устройство → список трансформов инстансов
    kindof: Dict[str, str] = {}
    for mg in magnets:
        dev = str(mg.get("dev", "") or "")
        t = mg.get("t")
        if dev and t:
            bydev[dev].append(t)
            kindof.setdefault(dev, str(mg.get("kind", "")))
    bases, total = {}, 0                                         # загрузка STL раз на устройство + оценка треуг.
    for dev, tl in bydev.items():
        stl = (hw.get(dev, {}) or {}).get("stl_path", "") or ""
        if stl in ("", "-", "nan"):
            continue
        p = os.path.join(ws, "geometry", stl)
        if not os.path.exists(p):
            continue
        try:
            base = pv.read(p).scale(0.001, inplace=False).triangulate()   # мм→м
        except Exception:
            continue
        bases[dev] = base
        total += (len(np.asarray(base.faces)) // 4) * len(tl)
    if not bases:
        return {}
    red = max(0.0, min(0.95, 1.0 - float(tri_cap) / total)) if total > tri_cap else 0.0
    out = {}
    for dev, base in bases.items():
        b = base
        if red > 0:
            try:
                b = base.decimate(red)
            except Exception:
                b = base
        try:
            big = vz.build_combined_mesh(b, np.asarray(bydev[dev], float), pre_rot=_pre_rot(hw.get(dev, {})))
        except Exception:
            continue
        out[dev] = (big, kindof.get(dev, ""))
    return out


# ------------------------------------------------------------------ вид (PyVista в Qt) — ленивый импорт
def make_ring_scene(parent=None, gap_scale: float = 1.0):
    """Фабрика `RingScene` (встроенный PyVista). Вынесена, чтобы PyVista/Qt грузились только при GUI."""
    from pyvistaqt import QtInteractor
    from PySide6.QtCore import Signal, QTimer, Qt
    import pyvista as pv

    class RingScene(QtInteractor):
        """Встроенный PyVista-вид: кольцо обводом, мониторы точками, стартовая плоскость, ондуляторы
        парами блоков-челюстей (расстояние = зазор), облако адресуемых элементов для пикинга.
        `update_snapshot` двигает челюсти по снимку. Пикинг (левый клик) → сигнал `picked(addr, name)`."""

        picked = Signal(str, str)        # (адрес element.<name>, человеко-имя)
        hovered = Signal(str)            # имя под курсором (best-effort)
        picked_screen = Signal(str)      # (screen_id) — клик по плоскости экрана
        picked_undulator = Signal(str)   # (gid) — клик по модели ондулятора (гребёнка полюсов)
        picked_sensor = Signal(int)      # (idx = индекс ebpm) — клик по экрану-сенсору (BPM)
        focused = Signal(str)            # камера перелетела на элемент (Ctrl/Shift+←/→) — подпись для статуса

        def __init__(self, parent=None, gap_scale=1.0):
            super().__init__(parent)
            self.model: Optional[SceneModel] = None
            self._jaw_actors: Dict[str, tuple] = {}     # gid → (actor_up, actor_lo, base_xy)
            self._target_gap: Dict[str, float] = {}     # gid → целевой зазор (из снимка движка)
            self._vis_gap: Dict[str, float] = {}        # gid → ПОКАЗЫВАЕМЫЙ зазор (плавно догоняет целевой)
            self._anim_last = None                       # время последнего кадра анимации (perf_counter)
            self._gap_scale = float(gap_scale)
            self._rc_press = None                        # позиция ПКМ-нажатия (отличить клик от перетаскивания)
            self.set_background("white")                  # белый фон (как в visualize.py)
            self._enable_quality()

        def _enable_quality(self):
            """Качество VTK-рендера как в visualize.py: сглаживание FXAA (пост-процесс, безопасно) + свет.
            depth peeling/MSAA НЕ включаем — на части GPU «гасят» геометрию; прозрачность — per-actor opacity."""
            try:
                self.enable_anti_aliasing("fxaa")
            except Exception:
                pass
            try:
                self.enable_lightkit()                   # мягкий свет — формы читаются объёмнее
            except Exception:
                pass

        def set_geometry(self, bundle: dict):
            """Построить актёров из узла геометрии (однократно при готовности модели)."""
            self.clear()
            self._jaw_actors.clear()
            self.model = SceneModel(bundle, gap_scale=self._gap_scale)
            bb = self.model.bbox()
            span = (bb[1] - bb[0]) if bb else 1.0

            poly = self.model.ring_polyline()
            if len(poly) >= 2:                                   # кольцо: обвод линией (не пикаем)
                line = pv.lines_from_points(poly, close=True)
                self.add_mesh(line, color="#1a5fb4", line_width=2, name="ring", pickable=False)
            sensors = self.model.sensors                         # сенсоры (BPM) — тонкие полупрозрачные экраны
            if sensors:
                items = []
                for s in sensors:
                    c = np.array([float(s["xy"][0]), float(s["xy"][1]), 0.0])
                    ax = s.get("axis", [1.0, 0.0])
                    a3 = np.array([float(ax[0]), float(ax[1]), 0.0])
                    a3 = a3 / (float(np.hypot(a3[0], a3[1])) or 1.0)
                    u3 = np.array([-a3[1], a3[0], 0.0])          # тонко вдоль пучка, широко поперёк (апертура)
                    items.append((c, a3, u3, SENSOR_THICK * 0.5, SENSOR_SIZE * 0.5, SENSOR_SIZE * 0.5))
                smesh = self._boxes_mesh(items)
                if smesh is not None:
                    self.add_mesh(smesh, color="#1f9e89", opacity=0.45, name="sensors", pickable=True)
            self._add_magnet_bodies()                            # тела магнитов (STL) — ЕДИНСТВЕННЫЕ пикаемые
            sp = self.model.start_plane()
            if sp is not None:                                   # стартовая плоскость на s=0 (небольшая)
                c, n = sp
                self.add_mesh(pv.Plane(center=c, direction=n,
                                       i_size=START_PLANE_SIZE, j_size=START_PLANE_SIZE),
                              color="#59a14f", opacity=0.5, name="start", pickable=False)
            for gid in self.model.undulator_list():              # ондулятор: ДВА РЯДА ПОЛЮСОВ-зубьев (челюсти)
                u_info = self.model.undulators[gid]
                gap0 = float(u_info.get("gap_mm", 5.0))
                off0 = self.model.jaw_offset(gap0)
                poles = u_info.get("poles") or []
                if len(poles) >= 2:                              # ось — вдоль массива полюсов (точный ход)
                    a3 = np.array([poles[-1][0] - poles[0][0], poles[-1][1] - poles[0][1], 0.0])
                else:
                    ax = u_info.get("axis", [1.0, 0.0])
                    a3 = np.array([float(ax[0]), float(ax[1]), 0.0])
                a3 = a3 / (float(np.hypot(a3[0], a3[1])) or 1.0)
                u3 = np.array([-a3[1], a3[0], 0.0])              # поперёк оси пучка, в плоскости кольца
                if poles:                                        # зубья по ЦЕНТРАМ полюсов = сам ондулятор
                    span = float(np.hypot(poles[-1][0] - poles[0][0], poles[-1][1] - poles[0][1]))
                    sp = span / max(len(poles) - 1, 1) if len(poles) >= 2 else 0.01
                    reg = max(0.35 * sp, 0.002)
                    items = []
                    for cell in poles:                           # cell = [x, y, is_end]; концевые длиннее
                        px, py = float(cell[0]), float(cell[1])
                        hl = reg * (3.0 if (len(cell) > 2 and cell[2]) else 1.0)
                        items.append((np.array([px, py, 0.0]), a3, u3, hl, UND_HALF_W, UND_HALF_H))
                    m_up = self._boxes_mesh(items)
                    m_lo = self._boxes_mesh(items)
                else:                                            # запас (фейк-plant без полюсов): один бар у якоря
                    xy = u_info.get("xy", [0.0, 0.0])
                    cc = np.array([float(xy[0]), float(xy[1]), 0.0])
                    m_up = self._boxes_mesh([(cc, a3, u3, UND_HALF_L, UND_HALF_W, UND_HALF_H)])
                    m_lo = self._boxes_mesh([(cc, a3, u3, UND_HALF_L, UND_HALF_W, UND_HALF_H)])
                a_up = self.add_mesh(m_up, color="#ffb000", show_edges=True, edge_color=EDGE_COLOR,
                                     line_width=1, name=f"und_up_{gid}", pickable=True)
                a_lo = self.add_mesh(m_lo, color="#ffb000", show_edges=True, edge_color=EDGE_COLOR,
                                     line_width=1, name=f"und_lo_{gid}", pickable=True)
                a_up.SetPosition(0.0, 0.0, off0)                 # ряды построены при z=0 → разводим на ±зазор
                a_lo.SetPosition(0.0, 0.0, -off0)
                self._jaw_actors[gid] = (a_up, a_lo, None)
                self._vis_gap[gid] = gap0                        # старт анимации — с текущего зазора
                self._target_gap[gid] = gap0
            for sc in (self.model.screens or []):                 # экраны + линия ондулятор→экран (LOS)
                c = np.array([float(sc["center"][0]), float(sc["center"][1]), 0.0])
                nrm = np.array([float(sc["axis"][0]), float(sc["axis"][1]), 0.0])
                nrm = nrm / (float(np.hypot(nrm[0], nrm[1])) or 1.0)
                sz = max(float(sc.get("w", 0) or 0), float(sc.get("h", 0) or 0), SCREEN_MIN_SIZE)
                self.add_mesh(pv.Plane(center=c, direction=nrm, i_size=sz, j_size=sz), color="#17c3d6",
                              opacity=0.6, name=f"screen_{sc.get('sid', '')}", pickable=True)
                src = np.array([float(sc["src"][0]), float(sc["src"][1]), 0.0])
                self.add_mesh(pv.lines_from_points(np.vstack([src, c])), color="#ff8c00",
                              line_width=2, name=f"los_{sc.get('sid', '')}", pickable=False)
            self._enable_picking()
            self._add_grid_axes()                                 # сетка + оси (как в visualize.py)
            self.fit_view()                                       # вид сверху после построения актёров

        def _add_grid_axes(self):
            """Сетка и оси как в visualize.py: show_bounds(grid='back', location='outer') + add_axes."""
            try:
                self.show_bounds(grid="back", location="outer", all_edges=False,
                                 xtitle="X, м", ytitle="Y, м", ztitle="Z, м",
                                 font_size=10, color="black")
            except Exception:
                pass
            try:
                self.add_axes()
            except Exception:
                pass

        def fit_view(self):
            """Вписать камеру и смотреть сверху (планарное кольцо в X-Y). Зовётся и ПОСЛЕ показа окна —
            у интерактора уже валидный размер, иначе вид остаётся косым (кольцо «эллипсом»)."""
            try:
                self.reset_camera()
                self.view_xy()
                self.render()
            except Exception:
                pass

        def _add_magnet_bodies(self):
            """Тела магнитов: НАСТОЯЩИЕ STL-инстансы (как viz_html). Если STL/hardware недоступны —
            откат на ориентированные блоки, чтобы сцена не осталась пустой."""
            if self.model is None:
                return
            try:
                meshes = build_magnet_stl_meshes(self.model.magnets)
            except Exception:
                meshes = {}
            if meshes:
                for dev, (mesh, kind) in meshes.items():
                    self.add_mesh(mesh, color=KIND_COLOR.get(kind, DEFAULT_KIND_COLOR),
                                  opacity=BODY_OPACITY, show_edges=True, edge_color=EDGE_COLOR,
                                  line_width=1, name=f"mag_{dev}", pickable=True)
                return
            self._add_box_bodies()                               # запасной путь: блоки

        def _add_box_bodies(self):
            """Запасные тела магнитов — ориентированные блоки по типам (один меш на тип)."""
            bykind: Dict[str, list] = defaultdict(list)
            for (c, a, u, hl, hw, hh, kind) in self.model.magnet_boxes():
                bykind[kind].append((c, a, u, hl, hw, hh))
            for kind, items in bykind.items():
                mesh = self._boxes_mesh(items)
                if mesh is not None:
                    self.add_mesh(mesh, color=KIND_COLOR.get(kind, DEFAULT_KIND_COLOR),
                                  opacity=BODY_OPACITY, show_edges=True, edge_color=EDGE_COLOR,
                                  line_width=1, name=f"mag_{kind}", pickable=True)

        def _boxes_mesh(self, items):
            """Один объединённый меш из ориентированных ящиков. items: (c,a,u,hl,hw,hh)."""
            signs = [(sx, sy, sz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
            quads = [(0, 1, 3, 2), (4, 5, 7, 6), (0, 1, 5, 4), (2, 3, 7, 6), (0, 2, 6, 4), (1, 3, 7, 5)]
            pts, faces, base = [], [], 0
            for (c, a, u, hl, hw, hh) in items:
                z = np.cross(a, u)
                z = z / (float(np.linalg.norm(z)) or 1.0)
                for (sx, sy, sz) in signs:
                    pts.append(c + sx * hl * a + sy * hw * u + sz * hh * z)
                for q in quads:
                    faces.append([4, base + q[0], base + q[1], base + q[2], base + q[3]])
                base += 8
            if not pts:
                return None
            return pv.PolyData(np.asarray(pts, float), np.asarray(faces, int).ravel())

        def _enable_picking(self):
            """Интерактив по магнитам БЕЗ конфликта с камерой: ЛКМ-клик (без перетаскивания) — выравнивание
            камеры на элемент; ПКМ-клик (без перетаскивания) — справочное окно (магнит/экран/ондулятор/сенсор);
            HOVER — жёлтая рамка-подсветка. ЛКМ-drag = вращение, ПКМ-drag = зум (дефолт VTK).
            ПРОТИВ ЗАЛИПАНИЯ ПКМ-зума: ПКМ НЕ вешаем на VTK-наблюдатели (их САМО НАЛИЧИЕ на Right-Button
            Press/Release ломало завершение dolly у интерактора — эмпирически, правки 14/20). ПКМ-клик ловим на
            уровне Qt (`mousePressEvent`/`mouseReleaseEvent` виджета), ПОСЛЕ `super()` — VTK полностью отрабатывает
            зум сам. ЛКМ и hover остаются VTK-наблюдателями (не залипают). Hover не пикает во время
            `_dragging`/`_rmb_down`. Headless — молча; остаётся pick_nearest."""
            self._lc_press = None
            self._rc_press = None
            self._rmb_down = False                           # ПКМ зажата (гасим hover-пик во время зума)
            self._dragging = False
            self._hover_name = None
            self._hover_actor = None
            self._focus = None                           # индекс model.nav, на котором стоит камера
            self._fly_timer = None
            self.fly_steps = 10                          # кадров перелёта камеры (0 — мгновенно)
            self._prop_picker = None
            try:
                import vtk
                self._prop_picker = vtk.vtkPropPicker()      # быстрый hardware-picker (hover/клик)
            except Exception:
                self._prop_picker = None
            for ev, cb in (("LeftButtonPressEvent", self._lc_down),   # ПКМ — НЕ здесь (ловим на уровне Qt)
                           ("LeftButtonReleaseEvent", self._lc_up),
                           ("MouseMoveEvent", self._on_move)):
                try:
                    self.iren.add_observer(ev, cb)
                except Exception:
                    pass

        def _event_pos(self):
            try:
                return tuple(self.iren.interactor.GetEventPosition())
            except Exception:
                return None

        def _pick_prop(self, pos=None):
            """Имя актёра + точка под курсором (быстрый prop-picker). (actor_name, x, y) или None.
            `pos`=(px,py) экранные координаты пикселей — для ОТЛОЖЕННОГО пика (позиция захвачена на release,
            чтобы не перечитывать интерактор вне колбэка); None → текущая позиция события."""
            if self._prop_picker is None:
                return None
            if pos is None:
                pos = self._event_pos()
            if pos is None:
                return None
            try:
                self._prop_picker.Pick(float(pos[0]), float(pos[1]), 0, self.renderer)
                actor = self._prop_picker.GetActor()
                wp = self._prop_picker.GetPickPosition()
            except Exception:
                return None
            if actor is None or wp is None or len(wp) < 2:
                return None
            for nm, act in list(self.actors.items()):
                if act is actor:
                    return (str(nm), float(wp[0]), float(wp[1]))
            return None

        @staticmethod
        def _mag_tag(name):
            return name[4:] if name.startswith("mag_") else None

        def _mag_at(self, hit):
            """hit=(name,x,y) → (name_el, mx, my) магнита того же типа; None, если не магнит / не нашли."""
            if self.model is None or hit is None:
                return None
            tag = self._mag_tag(hit[0])
            if not tag:
                return None
            res = self.model.nearest_named(hit[1], hit[2], max_dist=PICK_THRESH, tag=tag)
            if res is None:
                return None
            xy = self.model.named[res[0]]["xy"]
            return (res[1], float(xy[0]), float(xy[1]))

        def _is_click(self, p0, p1):
            return bool(p0) and bool(p1) and abs(p1[0] - p0[0]) <= 3 and abs(p1[1] - p0[1]) <= 3

        # ---- левый клик (без перетаскивания): выровнять камеру на элемент (дистанция сохраняется) ----
        def _lc_down(self, *a):
            self._lc_press = self._event_pos()
            self._dragging = True

        def _lc_up(self, *a):
            p0, self._lc_press = self._lc_press, None
            self._dragging = False
            if not self._is_click(p0, self._event_pos()):     # ЛКМ-перетаскивание = вращение, не клик
                return
            hit = self._pick_prop()
            if not hit:
                return
            name = hit[0]
            nm = self._mag_at(hit)                            # магнит → центр элемента; иначе точка под курсором
            if nm:
                self._center_on(nm[1], nm[2])
                self._set_focus_at(nm[1], nm[2])
            elif name.startswith(("und_up_", "und_lo_")) and self._focus_undulator(name[len("und_up_"):]):
                pass                                          # ондулятор — центр его гребёнки (как при навигации)
            elif name.startswith(("screen_", "und_up_", "und_lo_", "sensors")):
                self._center_on(hit[1], hit[2])

        # ---- ПКМ на уровне Qt (БЕЗ VTK-наблюдателей — иначе залипает dolly) ----
        @staticmethod
        def _qt_xy(ev):
            try:
                p = ev.position()                             # Qt6 (QPointF, виджет-локальные)
            except Exception:
                p = ev.pos()                                  # запас (Qt5)
            return (float(p.x()), float(p.y()))

        def mousePressEvent(self, ev):
            try:
                if ev.button() == Qt.MouseButton.RightButton:
                    self._rc_press = self._qt_xy(ev)
                    self._rmb_down = True                     # ПКМ зажата — hover не пикает (не мешает зуму)
            except Exception:
                pass
            super().mousePressEvent(ev)                       # VTK начинает dolly как обычно

        def mouseReleaseEvent(self, ev):
            super().mouseReleaseEvent(ev)                     # СНАЧАЛА VTK полностью завершает dolly — не залипает
            try:
                if ev.button() != Qt.MouseButton.RightButton:
                    return                                    # ЛКМ/СКМ — обрабатывают VTK-наблюдатели
                self._rmb_down = False
                p0, self._rc_press = self._rc_press, None
                if not self._is_click(p0, self._qt_xy(ev)):   # ПКМ-перетаскивание = зум, не клик
                    return
                pos = self._event_pos()                       # VTK-координаты release (super выставил), захват СЕЙЧАС
                if pos is not None:
                    QTimer.singleShot(0, lambda p=pos: self._rc_info(p))   # пик+сигнал вне обработчика события
            except Exception:
                pass

        def _rc_info(self, pos):
            """Отложенный пик по захваченной позиции ПКМ-клика → сигнал справочного окна (вне обработчика мыши)."""
            hit = self._pick_prop(pos)
            if not hit:
                return
            name = hit[0]
            if name.startswith("mag_"):                       # ПКМ-клик — справочное окно
                self.pick_nearest(hit[1], hit[2], tag=name[4:])
            elif name.startswith("screen_"):
                try:
                    self.picked_screen.emit(name[len("screen_"):])
                except Exception:
                    pass
            elif name.startswith("und_up_") or name.startswith("und_lo_"):
                try:
                    self.picked_undulator.emit(name[len("und_up_"):])
                except Exception:
                    pass
            elif name == "sensors":
                s = self.model.nearest_sensor(hit[1], hit[2]) if self.model else None
                if s is not None:
                    try:
                        self.picked_sensor.emit(int(s["idx"]))
                    except Exception:
                        pass

        # ---- навигация камеры: Ctrl+←/→ — соседний элемент, Shift+←/→ — соседний того же типа ----
        def _set_focus_at(self, x, y):
            """Запомнить элемент навигации в точке (x, y) как текущий фокус камеры."""
            res = self.model.nav_nearest(x, y) if self.model is not None else None
            self._focus = res[0] if res is not None and res[1] <= NAV_LOCK else None

        def _focus_undulator(self, gid):
            if self.model is None:
                return False
            for j, it in enumerate(self.model.nav):
                if it["kind"] == "undulator" and it["name"] == gid:
                    self._center_on(*it["xy"])
                    self._focus = j
                    return True
            return False

        def nav_step(self, step: int, same_type: bool = False):
            """Перелёт камеры на соседний элемент (step=+1 — по пучку, −1 — навстречу; same_type — того же типа
            устройства). Текущий — элемент, на котором стоит камера (последний клик/перелёт, если фокус не
            сдвигали), иначе ближайший к фокусу камеры в пределах NAV_SNAP. Возврат — подпись или None."""
            if self.model is None or not self.model.nav:
                return None
            try:
                f = np.array(self.camera.focal_point, float)
            except Exception:
                return None
            cur = None
            if self._focus is not None and self._focus < len(self.model.nav):
                xy = self.model.nav[self._focus]["xy"]
                if np.hypot(xy[0] - f[0], xy[1] - f[1]) <= NAV_LOCK:
                    cur = self._focus
            if cur is None:
                res = self.model.nav_nearest(f[0], f[1])
                if res is not None and res[1] <= NAV_SNAP:
                    cur = res[0]
            if cur is None:
                msg = "камера не на элементе — щёлкните по магниту или ондулятору, затем Ctrl/Shift + ←/→"
                try:
                    self.focused.emit(msg)
                except Exception:
                    pass
                return None
            j = self.model.nav_step(cur, step, same_type)
            it = self.model.nav[j]
            self._focus = j
            self._fly_to(it["xy"][0], it["xy"][1])
            self._hover_name = None
            self._set_hover(it["name"] if it["kind"] != "undulator" else None, it["xy"][0], it["xy"][1])
            label = self.model.nav_label(j)
            try:
                self.focused.emit("камера: " + label)
            except Exception:
                pass
            return label

        def _fly_to(self, mx, my):
            """Плавно перевести фокус камеры в (mx, my, 0), сохраняя направление и дистанцию (как _center_on).
            Новый перелёт поверх незаконченного начинается с текущего положения."""
            if self._fly_timer is not None:
                self._fly_timer.stop()
                self._fly_timer = None
            n = int(self.fly_steps)
            if n <= 0:
                self._center_on(mx, my)
                return
            try:
                cam = self.camera
                f0 = np.array(cam.focal_point, float); off = np.array(cam.position, float) - f0
            except Exception:
                return
            f1 = np.array([float(mx), float(my), 0.0])
            k = [0]
            timer = QTimer(self)

            def tick():
                k[0] += 1
                t = min(1.0, k[0] / n); w = t * t * (3.0 - 2.0 * t)          # smoothstep: мягкий старт/стоп
                f = f0 + (f1 - f0) * w
                try:
                    cam.focal_point = tuple(float(v) for v in f)
                    cam.position = tuple(float(v) for v in f + off)
                    self.render()
                except Exception:
                    pass
                if k[0] >= n:
                    timer.stop()
                    if self._fly_timer is timer:
                        self._fly_timer = None
            timer.timeout.connect(tick)
            self._fly_timer = timer
            timer.start(15)

        def _center_on(self, mx, my):
            """Центрировать камеру на (mx,my,0), СОХРАНЯЯ направление и дистанцию до фокуса (чистая панорама)."""
            try:
                cam = self.camera
                p = np.array(cam.position, float)
                f = np.array(cam.focal_point, float)
                off = p - f                                   # смещение камера→фокус: дистанция + направление
                cam.focal_point = (float(mx), float(my), 0.0)
                cam.position = (float(mx) + off[0], float(my) + off[1], off[2])
                self.render()
            except Exception:
                pass

        # ---- наведение: жёлтая рамка-подсветка магнита под курсором ----
        def _on_move(self, *a):
            if self._dragging or self._rmb_down or self.model is None:   # во время вращения/зума не пикаем
                return
            nm = self._mag_at(self._pick_prop())
            name = mx = my = None
            if nm:
                name, mx, my = nm
            if name == self._hover_name:                      # тот же магнит — без лишних перерисовок
                return
            self._hover_name = name
            self._set_hover(name, mx, my)
            try:
                self.hovered.emit(name or "")
            except Exception:
                pass

        def _set_hover(self, name, mx, my):
            if self._hover_actor is not None:
                try:
                    self.remove_actor(self._hover_actor, render=False)
                except Exception:
                    pass
                self._hover_actor = None
            if name and mx is not None:
                try:
                    box = pv.Cube(center=(float(mx), float(my), 0.0),
                                  x_length=HOVER_SIZE, y_length=HOVER_SIZE, z_length=HOVER_SIZE * 0.6)
                    self._hover_actor = self.add_mesh(box, color="#ffd400", style="wireframe",
                                                      line_width=3, name="_hover", pickable=False)
                except Exception:
                    self._hover_actor = None
            try:
                self.render()
            except Exception:
                pass

        def pick_nearest(self, x: float, y: float, tag: Optional[str] = None):
            """Найти элемент у (x,y) и испустить picked(addr, name). `tag` — dev/kind кликнутого тела
            (точный пик). Возврат — имя."""
            if self.model is None:
                return None
            hit = self.model.nearest_named(x, y, max_dist=PICK_THRESH, tag=tag)
            if hit is None and tag is not None:               # фильтр по телу не дал — попробуем без него
                hit = self.model.nearest_named(x, y, max_dist=PICK_THRESH)
            if hit is None:
                return None
            _i, name, _kind, _d = hit
            self.picked.emit(f"element.{name}", name)
            self.hovered.emit(name)
            return name

        def update_snapshot(self, snap):
            """Запомнить ЦЕЛЕВОЙ зазор из снимка движка. Сами челюсти двигает animate_jaws() ПЛАВНО —
            иначе рывки: движок меняет зазор дискретно за цикл и успевает много циклов между кадрами UI."""
            if self.model is None or snap is None:
                return
            for gid in self._jaw_actors:
                gap = snap.values.get(f"undulator.{gid}.gap")
                if gap is not None:
                    self._target_gap[gid] = float(gap)

        def animate_jaws(self):
            """Плавно догонять ПОКАЗЫВАЕМЫМ зазором целевой (из движка) с правдоподобной скоростью
            `JAW_VIS_RATE` мм/с (стенного времени). Рендерим ТОЛЬКО пока идёт движение (иначе сцена покоится)."""
            if self.model is None or not self._jaw_actors:
                return
            from time import perf_counter
            now = perf_counter()
            dt = 0.0 if self._anim_last is None else min(0.1, now - self._anim_last)
            self._anim_last = now
            step = JAW_VIS_RATE * dt                               # максимум мм за кадр
            changed = False
            for gid, (a_up, a_lo, _xy) in self._jaw_actors.items():
                tgt = self._target_gap.get(gid)
                if tgt is None:
                    continue
                cur = self._vis_gap.get(gid, tgt)
                if abs(tgt - cur) < 1e-4:
                    if cur != tgt:
                        self._vis_gap[gid] = tgt                   # дощёлкнуть вплотную
                    continue
                cur += max(-step, min(step, tgt - cur))            # шаг к цели, не перелетая
                self._vis_gap[gid] = cur
                off = self.model.jaw_offset(cur)
                try:
                    a_up.SetPosition(0.0, 0.0, off)
                    a_lo.SetPosition(0.0, 0.0, -off)
                    changed = True
                except Exception:
                    pass
            if changed:
                try:
                    self.render()
                except Exception:
                    pass

    return RingScene(parent, gap_scale=gap_scale)
