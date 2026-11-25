# -*- coding: utf-8 -*-
"""Кольцо на продакшн-ядре C++ (MatrixTracker --ring-serve) с интерфейсом ring_analysis.Ring.

Трекинг — только в C++: оборот, пооборотные координаты, выходы магнитов, срезы изгибов считает сервер
кольца (src/matrix_ring_serve.h), частицы параллельно (OpenMP). Python здесь — лишь транспорт и
анализ. Для аналитического механизма (ring_matrix, ring_response, ring_map) нужны описания магнитов
(атрибуты mags/plane) — с шага 3 (docs/LIBRARY_ARCHITECTURE.md) они берутся ИЗ ЯДРА (команда DESCRIBE /
Ring.describe_flat: рамки, силы, флаги тела — текущее состояние со всеми возмущениями), без второй сборки
решётки в Python и без сверки геометрии двух моделей.

  CppRing(lattice, start_before_first_m, sr, isr, rf, design_mev, nsl, order, seed)
    track(states, n_turns)            → (n, 6), NaN у потерянных
    track_elements(states)            → (n, n_mag+1, 6)
    track_tbt(states, n_turns, limit) → (n, n_turns+1, 6)
    track_pieces(states)              → (n, n_pieces+1, 6) — после каждого среза изгиба
    with_kick(j, θ) / with_yaw(j, α)  → новое кольцо с возмущением (C++ и описание магнитов)
    perturbed(kicks, yaws)            → то же для нескольких возмущений
    apply_perturb(kicks, …, dk3s)     → транзитные возмущения на том же сервере (RESET + набор)
    set_rf(volt_mv=…, freq_hz=…)      → параметры ВЧ на лету (частота ВЧ — актуатор энергии), переживает RESET
    junction(states, A, B)            → стык двух произвольных плоскостей (проверка матриц стыков)

rf = dict(volt_mv, harmonic, phase_deg[, freq_hz]): без freq_hz f_RF = h·c/C (как всегда); с freq_hz —
расстроенный резонатор, C_RF = h·c/f_RF (см. шапку src/matrix_ring_serve.h).

Возмущения двух сортов: «базовые» (with_*/perturbed — часть кольца, переживают RESET, копии их наследуют) и
транзитные (apply_perturb — заменяются его следующим вызовом). Описание магнитов mags/plane всегда по обоим.

Протокол сервера — в шапке src/matrix_ring_serve.h.

БЭКЕНДЫ (одно C++-ядро, одни числа бит-в-бит; шаг 2 docs/LIBRARY_ARCHITECTURE.md):
  * "ext"  — модуль расширения synchrotronsim._engine В ПРОЦЕССЕ Python: решётка из модели Machine
             (synchrotronsim.db, ws ⇒ overlay), RESET — откат к снимку свежего кольца (~6 мс против ~0,3 с),
             копии with_*/perturbed — клон (~30 мс против ~1 с запуска процесса), GIL отпущен на трекинге;
  * "pipe" — сервер MatrixTracker --ring-serve в своём процессе (изоляция падений, сторож timeout).
Выбор: аргумент backend= → env SS_ENGINE=ext|pipe → по умолчанию "ext", если модуль собран и timeout не задан
(сторож убивает зависший ПРОЦЕСС — у расширения его нет), иначе "pipe".
"""
import os, copy, math, subprocess, atexit, weakref
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(ROOT, "build", "MatrixTracker.exe")
MC2_MEV = 0.51099895069
_LIVE = weakref.WeakSet()
BACKENDS = ("ext", "pipe")


# Срезы для ЛИНЕЙНОЙ ОПТИКИ (тюны, β, η) — как паспорт facility_report. Дефолт 5 (ниже) — конвенция паритета
# с elegant (N_SLICES=5 в .lte); он смещает тюн СИЛА-26 на ~0,03 (νy 33,367 вместо 33,398, ошибка ∝ 1/N⁴,
# в основном DQ1B/D с √|K1|·L = 1,43) и ηx у ондуляторов на ~4%. Кто считает оптику — передаёт nsl явно.
OPTICS_NSL = 40


def _default_slices():
    """Срезы и порядок по умолчанию — как у Python-ядра: env SS_SECTOR_NSL/SS_SECTOR_ORDER или 5 и 4.
    5 — ради побитной сверки с elegant; для тюнов/оптики передавайте nsl=OPTICS_NSL."""
    return (int(os.environ.get("SS_SECTOR_NSL", "") or 5), int(os.environ.get("SS_SECTOR_ORDER", "") or 4))


def _ext_available():
    try:
        import synchrotronsim
        return bool(synchrotronsim.ENGINE_AVAILABLE)
    except Exception:
        return False


def choose_backend(backend=None, timeout=None):
    """Бэкенд кольца: явный → env SS_ENGINE → "ext" (если собран и нет timeout) → "pipe"."""
    b = (backend or os.environ.get("SS_ENGINE") or "").strip().lower()
    if b and b not in BACKENDS:
        raise ValueError(f"backend: {BACKENDS}, не {b!r}")
    if b == "ext" and not _ext_available():
        raise RuntimeError("backend='ext': модуль synchrotronsim._engine не собран (см. docs/MANUAL.md)")
    if b:
        return b
    return "ext" if (not timeout and _ext_available()) else "pipe"


# Раскладка описания узлов — DESCRIBE_COLS в src/matrix_ring_serve.h (команда DESCRIBE, Ring.describe_flat);
# совпадение строки с C++ проверяет тест. «имя:ширина», ширина 1 опускается; затем 24 числа плоскости.
DESCRIBE_COLS = ("center:3 axis:3 uhat:3 vhat:3 entry:3 exit:3 fr_in:12 fr_out:12 Leff theta0 k1 k1v k2 k3 ks "
                 "hkick vkick edge1 edge2 kEdge1 kEdge2 psi1 psi2 rinv rad_coef isr_coef lf1 lf2 theta0_0 edge1_0 "
                 "edge2_0 sector kick same_plane_as_prev nsl nsl_box")


def parse_describe(flat, n):
    """Плоское описание узлов ядра → {столбец: массив} (векторы (n,3), рамки (n,4,3), скаляры (n,)) + рамки
    плоскости наблюдения plane_fr_in/plane_fr_out (4,3)."""
    cols = []
    for tok in DESCRIBE_COLS.split():
        name, _, w = tok.partition(":")
        cols.append((name, int(w or 1)))
    width = sum(w for _, w in cols)
    flat = np.asarray(flat, float)
    if flat.size != n * width + 24:
        raise RuntimeError(f"DESCRIBE: {flat.size} чисел, ожидалось {n}×{width}+24 — раскладка C++ и Python разошлась")
    rows = flat[: n * width].reshape(n, width)
    out, c = {}, 0
    for name, w in cols:
        a = rows[:, c:c + w]
        out[name] = a.reshape(n, 4, 3) if w == 12 else (a if w > 1 else a[:, 0].copy())
        c += w
    out["plane_fr_in"] = flat[n * width:n * width + 12].reshape(4, 3)
    out["plane_fr_out"] = flat[n * width + 12:].reshape(4, 3)
    return out


class _PipeBackend:
    """Сервер MatrixTracker --ring-serve: свой процесс, двоичный протокол stdin/stdout."""
    kind = "pipe"

    def __init__(self, ring):
        if not os.path.exists(EXE):
            raise FileNotFoundError(f"нет {EXE}: соберите cmake --build build --target MatrixTracker")
        self.timeout = ring.timeout
        rf = ring.rf or {}
        args = [EXE, "--ring-serve", str(ring.lattice),
                "nan" if ring.lead is None else repr(float(ring.lead)), repr(ring.design_mev),
                "1" if ring.sr else "0", "1" if ring.isr else "0", str(ring.nsl), str(ring.order),
                "1" if ring.rf else "0", repr(float(rf.get("volt_mv", 0.0))), repr(float(rf.get("harmonic", 0.0))),
                repr(float(rf.get("phase_deg", 0.0))), str(ring.seed)]
        if ring.ws:                             # возмущённая рабочая папка (ансамбли ошибок) — 13-й аргумент
            args.append(ring.ws)
        env = dict(os.environ)
        env.pop("SS_NSL_DEV", None)
        if ring.nsl_dev:
            env["SS_NSL_DEV"] = ring._nsl_dev_env()
        self._proc = subprocess.Popen(args, cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=subprocess.DEVNULL, bufsize=0, env=env)
        head = self._proc.stdout.readline().decode("utf-8", "replace").split()
        if not head or head[0] != "READY":
            raise RuntimeError(f"сервер кольца не поднялся: {' '.join(head)}")
        self.n_mag, self.circumference, self.n_pieces = int(head[1]), float(head[2]), int(head[3])
        info = self._text("INFO")
        self.s_in = np.array([float(l.split()[1]) for l in info])
        self.s_out = np.array([float(l.split()[2]) for l in info])
        self.kinds = [l.split()[4] if l.split()[4] != "-" else "" for l in info]
        self.names = [l.split()[5] for l in info]

    # ---- протокол ---------------------------------------------------------------------------------
    def _watch(self):
        """Сторож: таймер, убивающий процесс через timeout с, если команда зависла (защита ансамбля
        от патологичной машины — убитый процесс даёт EOF, чтение падает, ансамбль метит машину неустойчивой).
        Возврат таймера (его .cancel() по завершении команды) или None, если таймаут не задан."""
        if not self.timeout:
            return None
        import threading
        t = threading.Timer(self.timeout, lambda p=self._proc: p.kill())
        t.daemon = True; t.start()
        return t

    def _send(self, line, payload=None):
        payload = b"" if payload is None else np.ascontiguousarray(payload, dtype="<f8").tobytes()
        self._proc.stdin.write(line.encode() + b"\n" + payload); self._proc.stdin.flush()
        head = self._proc.stdout.readline().decode("utf-8", "replace").split()
        if not head or head[0] != "OK":
            raise RuntimeError(f"сервер кольца: {' '.join(head)}")
        return int(head[1])

    def _doubles(self, line, payload):
        w = self._watch()
        try:
            m = self._send(line, payload)
            need = m * 8; buf = bytearray()
            while len(buf) < need:
                chunk = self._proc.stdout.read(need - len(buf))
                if not chunk:
                    raise RuntimeError("сервер кольца оборвал ответ")
                buf += chunk
        finally:
            if w:
                w.cancel()
        return np.frombuffer(bytes(buf), dtype="<f8")

    def _text(self, line):
        w = self._watch()
        try:
            k = self._send(line)
            out = [self._proc.stdout.readline().decode("utf-8", "replace").rstrip("\n") for _ in range(k)]
        finally:
            if w:
                w.cancel()
        return out

    # ---- операции ---------------------------------------------------------------------------------
    def cmd(self, cmd, j, v):
        self._text(f"{cmd} {j} {v!r}")

    def reset(self):
        self._text("RESET")

    def rf(self, on, volt, harm, phase, freq):
        self._text(f"RF {1 if on else 0} {float(volt)!r} {float(harm)!r} {float(phase)!r} {float(freq)!r}")

    def track(self, a, n_turns):
        return self._doubles(f"TRACK {len(a)} {int(n_turns)}", a).reshape(len(a), 6)

    def track_tbt(self, a, n_turns, lim):
        return self._doubles(f"TBT {len(a)} {int(n_turns)} {lim!r}", a).reshape(len(a), int(n_turns) + 1, 6)

    def track_elements(self, a):
        return self._doubles(f"ELEM {len(a)}", a).reshape(len(a), self.n_mag + 1, 6)

    def track_pieces(self, a):
        return self._doubles(f"PIECES {len(a)}", a).reshape(len(a), self.n_pieces + 1, 6)

    def piece_meta(self):
        out = []
        for l in self._text("PIECEMETA"):
            v = l.split()
            out.append((int(v[0]), v[1] == "1", v[2] == "1", *map(float, v[3:8])))
        return out

    def junction(self, a, fr):
        m = self._send(f"JUNCTION {len(a)}", np.concatenate([fr, a.ravel()]))
        buf = bytearray()
        while len(buf) < m * 8:
            chunk = self._proc.stdout.read(m * 8 - len(buf))
            if not chunk:
                raise RuntimeError("сервер кольца оборвал ответ")
            buf += chunk
        return np.frombuffer(bytes(buf), dtype="<f8").reshape(len(a), 6)

    def describe(self):
        return parse_describe(self._doubles("DESCRIBE", None), self.n_mag)

    def fresh_copy(self, ring):
        return _PipeBackend(ring)

    def close(self):
        p = getattr(self, "_proc", None)
        if p and p.poll() is None:
            try:
                p.stdin.write(b"QUIT\n"); p.stdin.flush(); p.wait(timeout=5)
            except Exception:
                p.kill()
        self._proc = None


class _ExtBackend:
    """Модуль расширения synchrotronsim._engine в процессе Python (synchrotronsim.engine.ExtRing)."""
    kind = "ext"
    _CMD = {"KICK": "kick", "VKICK": "vkick", "YAW": "yaw", "DK": "dk", "DKV": "dkv", "DSK": "dsk",
            "DK2": "dk2", "DK3": "dk3", "BSCALE": "bscale"}

    def __init__(self, ring=None, _ext=None):
        from synchrotronsim.engine import ExtRing
        self.x = _ext if _ext is not None else ExtRing(
            lattice=ring.lattice, start_before_first_m=ring.lead, sr=ring.sr, isr=ring.isr,
            rf=dict(ring.rf) if ring.rf else None, design_mev=ring.design_mev, nsl=ring.nsl, order=ring.order,
            seed=ring.seed, ws=ring.ws, nsl_dev=dict(ring.nsl_dev), machine=ring.machine)
        x = self.x
        self.n_mag, self.circumference, self.n_pieces = x.n_mag, x.circumference, x.n_pieces
        self.s_in, self.s_out, self.kinds, self.names = x.s_in, x.s_out, list(x.kinds), list(x.names)

    def cmd(self, cmd, j, v):
        getattr(self.x._r, self._CMD[cmd])(int(j), float(v))

    def reset(self):
        self.x.reset()

    def rf(self, on, volt, harm, phase, freq):
        self.x._r.set_rf(bool(on), float(volt), float(harm), float(phase), float(freq))

    def track(self, a, n_turns):
        return self.x.track(a, n_turns)

    def track_tbt(self, a, n_turns, lim):
        return self.x.track_tbt(a, n_turns, lim)

    def track_elements(self, a):
        return self.x.track_elements(a)

    def track_pieces(self, a):
        return self.x.track_pieces(a)

    def piece_meta(self):
        return self.x.piece_meta()

    def junction(self, a, fr):
        return self.x.junction(a, fr[:12], fr[12:])

    def describe(self):
        return parse_describe(self.x._r.describe_flat(), self.n_mag)

    def fresh_copy(self, ring):
        """Новое кольцо той же конфигурации без возмущений: клон + откат к свежему (без сборки, ~30 мс)."""
        c = self.x.clone(); c.reset()
        return _ExtBackend(_ext=c)

    def close(self):
        self.x = None


class CppRing:
    label = "SynchrotronSim"

    def __init__(self, lattice=1, start_before_first_m=None, settings=None, sr=False, isr=False, rf=None,
                 design_mev=6000.0, nsl=None, order=None, seed=12345, ws=None, timeout=None, nsl_dev=None,
                 backend=None, machine=None):
        self.backend = choose_backend(backend, timeout)   # "ext" (в процессе) | "pipe" (сервер) — см. шапку
        dn, do = _default_slices()
        # машина (synchrotronsim.model.Machine): готовый документ — ядро ext строится из него, таблиц не читает;
        # None — по номеру решётки из таблиц (ws ⇒ overlay), как раньше
        if machine is not None and self.backend == "pipe" and str(machine.source).startswith("document:"):
            raise ValueError("CppRing: бэкенд pipe (сервер) читает таблицы сам; машина из документа — только ext")
        self.lattice = int(machine.lattice_id) if machine is not None else int(lattice)
        self._machine = machine
        self.lead = start_before_first_m
        self.sr, self.isr, self.rf, self.design_mev = bool(sr or isr), bool(isr), rf, float(design_mev)
        self.nsl = int(nsl or dn); self.order = int(order or do); self.seed = int(seed)
        # срезы сектор/кик-тела по типу устройства {device_type_id: n} (env SS_NSL_DEV у C++ и Python-описания);
        # не задано — как у Python-эталона: из окружения SS_NSL_DEV (по умолчанию пусто — всё как было)
        if nsl_dev is None:
            nsl_dev = dict(t.split("=", 1) for t in (os.environ.get("SS_NSL_DEV") or "").split(",") if "=" in t)
        self.nsl_dev = {}
        for k, v in nsl_dev.items():
            try:
                n = int(v)
            except (TypeError, ValueError):
                continue
            if k and n > 0:
                self.nsl_dev[str(k)] = n
        self.ws = ws                            # путь к workspace с ВОЗМУЩЁННОЙ решёткой (ансамбли ошибок); None — обычная БД
        self.timeout = timeout                  # сторож: команда дольше timeout с → убить процесс (защита от зависаний)
        self.settings = settings
        self._perturb = []                      # базовые [(cmd, j, value)] — with_*/perturbed: переживают RESET
        self._transient = []                    # транзитные — от последнего apply_perturb
        self._py = None
        self.rf_cavities = self._cavity_key() if rf else ()      # резонаторы установки — для отпечатка кэша откликов
        self._start()

    @property
    def machine(self):
        """Машина этого кольца (Machine): заданная при создании или из таблиц по номеру решётки (ws ⇒ overlay)."""
        if self._machine is None:
            from synchrotronsim.db import Repository
            self._machine = Repository(overlay=self.ws).machine(self.lattice)
        return self._machine

    def _cavity_key(self):
        """Резонаторы машины (id, точка, нормаль, доля) — для отпечатка кэша откликов; () — нет."""
        return tuple((c.cavity_id, *(repr(float(v)) for v in c.point), *(repr(float(v)) for v in c.normal),
                      repr(float(c.volt_share))) for c in self.machine.cavities)

    # ---- бэкенд ядра (ext — в процессе, pipe — сервер) ---------------------------------------------
    def _start(self):
        parent = self.__dict__.pop("_be_parent", None)          # копия (_spawn): бэкенд родителя — шаблон
        if parent is not None and parent.kind == self.backend:
            self._be = parent.fresh_copy(self)
        elif self.backend == "ext":
            self._be = _ExtBackend(self)
        else:
            self._be = _PipeBackend(self)
        be = self._be
        self.n_mag, self.circumference, self.n_pieces = be.n_mag, be.circumference, be.n_pieces
        self.s_in, self.s_out, self.kinds, self.names = be.s_in, be.s_out, list(be.kinds), list(be.names)
        if self.rf and float(self.rf.get("freq_hz") or 0.0) > 0:     # расстроенный ВЧ: частота — командой RF
            self._send_rf()
        self._replay(self._perturb + self._transient)
        _LIVE.add(self)

    def _send_rf(self):
        rf = self.rf or {}
        self._be.rf(bool(self.rf), float(rf.get("volt_mv", 0.0)), float(rf.get("harmonic", 0.0)),
                    float(rf.get("phase_deg", 0.0)), float(rf.get("freq_hz") or 0.0))

    def set_rf(self, **kw):
        """Сменить параметры ВЧ НА ЛЕТУ (volt_mv/harmonic/phase_deg/freq_hz) на том же кольце, без RESET:
        состояние магнитов и возмущения не трогаются, переживает apply_perturb. freq_hz ≤ 0/None ⇒ f_RF = h·c/C.
        Новый dict (копии кольца со старым rf не задеваются). Возврат self."""
        if not self.rf:
            raise ValueError("set_rf: кольцо построено без ВЧ (rf=None)")
        self.rf = {**self.rf, **kw}
        self._send_rf()
        return self

    def rf_circumference(self):
        """Длина ВЧ-оборота C_RF (м): h·c/f_RF при заданной частоте, иначе периметр."""
        f = float((self.rf or {}).get("freq_hz") or 0.0)
        return self.rf["harmonic"] * 299792458.0 / f if f > 0 else self.circumference

    def rf_frequency(self):
        """Текущая частота ВЧ, Гц: заданная freq_hz или h·c/C (резонатор на геометрическом периметре)."""
        f = float((self.rf or {}).get("freq_hz") or 0.0)
        return f if f > 0 else self.rf["harmonic"] * 299792458.0 / self.circumference

    def _replay(self, perturb):
        for cmd, j, v in perturb:
            self._be.cmd(cmd, j, v)

    def close(self):
        be = self.__dict__.get("_be")
        if be is not None:
            be.close()
        self._be = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    # ---- интерфейс трекера (как ring_analysis.Ring) -------------------------------------------------
    @staticmethod
    def _states(states):
        a = np.asarray(states, float)
        return np.ascontiguousarray(a.reshape(-1, 6))

    def track(self, states, n_turns=1):
        return self._be.track(self._states(states), int(n_turns))

    def track_tbt(self, states, n_turns, limit=None):
        return self._be.track_tbt(self._states(states), int(n_turns), float(limit) if limit else 0.0)

    def track_elements(self, states):
        return self._be.track_elements(self._states(states))

    def track_pieces(self, states):
        return self._be.track_pieces(self._states(states))

    def piece_meta(self):
        """[(магнит, стык перед куском, изгиб, Leff, θ0, k1, k2, k3)] — срезы изгибов и прочие магниты."""
        return self._be.piece_meta()

    def beta(self, delta):
        pc = self.design_mev * (1.0 + delta)
        return pc / math.sqrt(pc * pc + MC2_MEV ** 2)

    def _nsl_dev_env(self):
        return ",".join(f"{k}={v}" for k, v in sorted(self.nsl_dev.items()))

    # ---- описание магнитов для аналитического механизма (без трекинга) ----------------------------
    def _python_ring(self):
        """Python-описание кольца (ring_analysis.Ring с магнитами matrix_tracker.Magnet) для аналитики
        (ring_matrix, ring_response, ring_map, оптика). Шаг 3 docs/LIBRARY_ARCHITECTURE.md: состояние магнитов —
        ИЗ ЯДРА (команда DESCRIBE): рамки, силы, грани, фриндж, флаги тела и стыков — ровно то, что трекает C++,
        со всеми текущими возмущениями. Вторая сборка решётки и повторное проигрывание возмущений в Python
        больше не нужны, сверять геометрию не с чем. Из модели Machine берутся только то, чего у ядра нет:
        углы установки yaw/pitch/roll, обмотки devs, sensor_id, n_slices box-тела, edge_ff (kEdge).
        Как и прежде, описание — БЕЗ излучения и ВЧ (rad_coef = isr_coef = 0, rf = None): аналитика линейной
        оптики считается вокруг проектной орбиты."""
        if self._py is None:
            self._py = self._ring_from_engine()
        return self._py

    def _ring_from_engine(self):
        import ring_analysis as ra
        import matrix_tracker as mt
        from beam_injection import PlaneMagnet
        D = self._be.describe()
        n = self.n_mag
        M = self.machine
        E = M.elements
        if len(E) != n or [str(d) or "-" for d in E.device_type] != [str(d) or "-" for d in self.names]:
            raise RuntimeError("порядок магнитов модели Machine и C++-ядра разошёлся")   # не геометрия — индексы
        nsl_default = (self.settings.n_slices if self.settings is not None else mt.Settings().n_slices)
        weak = os.environ.get("SS_WEAK_FOCUS", "") == "1"
        ydk = mt._yoshida(self.order)
        mags = []
        for j in range(n):
            m = object.__new__(mt.Magnet)
            dev = str(E.device_type[j]); kind = str(E.kind[j])
            m.dev, m.kind = dev, kind
            sl = E.members(j)
            m.devs = [str(d) for d in M.instances.device_type[sl]]
            m.sensor_id = str(E.instance_id[j]); m.monitor = (kind == "monitor")
            m.yaw, m.pitch, m.roll = float(E.yaw[j]), float(E.pitch[j]), float(E.roll[j])
            m.center, m.axis, m.uhat, m.vhat = (D[k][j].copy() for k in ("center", "axis", "uhat", "vhat"))
            m.entry, m.exit = D["entry"][j].copy(), D["exit"][j].copy()
            m.fr_in = tuple(D["fr_in"][j][q].copy() for q in range(4))
            m.fr_out = tuple(D["fr_out"][j][q].copy() for q in range(4))
            m.Leff = m.L = float(D["Leff"][j])
            m.nsl = int(E.n_slices[j]) or nsl_default
            m.edge1, m.edge2 = float(D["edge1"][j]), float(D["edge2"][j])
            m.lf1, m.lf2 = float(D["lf1"][j]), float(D["lf2"][j])
            th = float(D["theta0"][j])
            m.k1, m.k1v, m.k2, m.k3 = (float(D[k][j]) for k in ("k1", "k1v", "k2", "k3"))
            m.ks, m.hkick, m.vkick = float(D["ks"][j]), float(D["hkick"][j]), float(D["vkick"][j])
            m._bend0 = (float(D["theta0_0"][j]), float(D["edge1_0"][j]), float(D["edge2_0"][j]))
            m.kEdge1, m.kEdge2, m.rinv = float(D["kEdge1"][j]), float(D["kEdge2"][j]), float(D["rinv"][j])
            m.psi1, m.psi2 = float(D["psi1"][j]), float(D["psi2"][j])
            dt = M.devices.get(dev)
            m.kEdge = (dt.edge_ff if dt else 0.0) * th ** 2 / m.Leff if m.Leff > 0 else 0.0
            m.wkfoc = m.rinv * m.rinv if weak else 0.0
            m.sector = bool(D["sector"][j]); m.kick_sector = bool(D["kick"][j]) and not m.sector
            m.rad_coef = m.isr_coef = 0.0                   # описание без излучения (как прежде)
            if m.sector or m.kick_sector:
                m._nsl_sector = int(D["nsl"][j]); m._yoshida_dk = ydk
            if m.sector:
                m.entry_pt, m.exit_pt = m.fr_in[0], m.fr_out[0]
            m.same_plane_as_prev = bool(D["same_plane_as_prev"][j])
            # гор. кик «ящика» (не изгиб, без поля — корректор) = его θ0 box-тела, как ra.kicked; у изгиба и
            # полевого элемента θ0 не трогается (у изгиба это кривизна h для ring_matrix), кик — в hkick/vkick
            boxed = bool(m.hkick) and not m.sector and not m.kick_sector
            m.theta0 = th + m.hkick if boxed else th
            m.xref = m.xref_in = m.xref_out = None
            if not boxed and th and any((m.k1, m.k1v, m.k2, m.k3, m.kEdge1, m.kEdge2)):
                m.xref, m.xref_in, m.xref_out = m._ref_bow()
            m.nonzero = any(abs(v) > 0 for v in (m.theta0, m.k1, m.k1v, m.k2, m.k3, m.kEdge, m.kEdge1,
                                                 m.kEdge2, m.wkfoc))
            mags.append(m)
        r = object.__new__(ra.Ring)
        r.mags = mags
        pin, pout = D["plane_fr_in"], D["plane_fr_out"]
        pl = PlaneMagnet(pin[0].copy(), pin[1].copy())
        pl.uhat, pl.vhat = pin[2].copy(), pin[3].copy()           # рамка плоскости — ровно ядра
        pl.fr_in = tuple(pin[q].copy() for q in range(4)); pl.fr_out = tuple(pout[q].copy() for q in range(4))
        r.plane = pl
        r.sr, r.isr, r.rf, r.design_mev = False, False, None, float(self.design_mev)
        r.lattice = self.lattice
        r.s_in, r.s_out = np.array(self.s_in, float), np.array(self.s_out, float)
        r.circumference = float(self.circumference)
        r.names = [m.dev for m in mags]; r.kinds = [m.kind for m in mags]
        r._cav_before = {}; r._cav_tail = []; r._end_share = None; r.rf_cavities = ()
        return r

    @property
    def mags(self):
        return self._python_ring().mags

    @property
    def plane(self):
        return self._python_ring().plane

    # ---- возмущения --------------------------------------------------------------------------------
    def _copy_with(self, cmd, j, value):
        return self._spawn([(cmd, int(j), float(value))])

    def _spawn(self, extra):
        """Новое кольцо (свой процесс) = это кольцо КАК ЕСТЬ (базовые + транзитные возмущения) + extra; всё
        это — базовые возмущения копии (порядок сохранён, переживают её RESET)."""
        r = copy.copy(self)
        r._be_parent, r._be = self._be, None    # копия не делит бэкенд (её close не тронет родителя)
        r._perturb = list(self._perturb) + list(self._transient) + list(extra)
        r._transient = []
        r._py = None
        r._start()
        return r

    def with_kick(self, j, theta):
        """Кольцо с однородным ГОР. киком θ в магните j (корректор; дрейф–кик–дрейф по n_slices из БД)."""
        return self._copy_with("KICK", j, theta)

    def with_vkick(self, j, theta):
        """Кольцо с однородным ВЕРТ. киком θ в магните j (корректор орбиты по y)."""
        return self._copy_with("VKICK", j, theta)

    def with_yaw(self, j, alpha):
        """Кольцо с магнитом j, повёрнутым на α вокруг вертикали входа (ring_analysis.yaw_element)."""
        return self._copy_with("YAW", j, alpha)

    def apply_perturb(self, kicks=(), vkicks=(), dks=(), dsks=(), dk2s=(), dk3s=(), dkvs=(), bscales=()):
        """Применить возмущения НА ТОМ ЖЕ сервере (RESET + KICK/VKICK/DK/DSK/DK2/DK3), без нового процесса. Для
        матрицы отклика, коррекции орбиты, BBA, коррекции связи и нелинейного тюнера. RESET пересобирает кольцо
        из ws (ошибки рассогласования СОХРАНЯЮТСЯ), затем возвращаются базовые возмущения кольца (with_*/
        perturbed) — сбрасываются только транзитные от прежнего apply_perturb. dks=[(j,dk)] — k1 квадруполя
        (BBA); dsks=[(j,ks)] — скью (связь); dk2s=[(j,dk2)]/dk3s=[(j,dk3)] — секступ./октуп. силы (нелинейный
        тюнер ν(δ)); dkvs=[(j,dk1v)] — ТОЛЬКО вертикальный градиент; bscales=[(j,f)] — масштаб изгиба полюса
        (угол и грани = f × исходные, рамки заново: живой зазор ондулятора).
        Описание mags/plane следует за ними. apply_perturb() — вернуть кольцо в исходное. Возврат self."""
        new = [(cmd, int(j), float(v)) for cmd, pairs in (("KICK", kicks), ("VKICK", vkicks), ("DK", dks),
                                                          ("DSK", dsks), ("DK2", dk2s), ("DK3", dk3s),
                                                          ("DKV", dkvs), ("BSCALE", bscales))
               for j, v in pairs]
        self._py = None
        self._transient = []                    # до RESET: при сбое описание не припишет кольцу чужих возмущений
        self._be.reset()
        self._replay(self._perturb)
        self._replay(new)
        self._transient = new
        return self

    _ADDITIVE = {"dks": "DK", "dkvs": "DKV", "dsks": "DSK", "dk2s": "DK2", "dk3s": "DK3"}

    def add_perturb(self, **channels):
        """Доложить АДДИТИВНЫЕ транзитные возмущения поверх текущих БЕЗ RESET: dks/dkvs/dsks/dk2s/dk3s=[(j, Δ)].
        Дёшево (одна команда на элемент, без пересборки кольца, против ~0,1 с у apply_perturb) — для живых
        медленных каналов, меняющихся каждый такт (фокус ИД от движущегося зазора). Кики (KICK/VKICK — «задать»,
        а не «прибавить») сюда не входят. Следующий apply_perturb, как обычно, сбрасывает все транзитные. Возврат self."""
        new = []
        for name, pairs in channels.items():
            cmd = self._ADDITIVE.get(name)
            if cmd is None:
                raise ValueError(f"add_perturb: канал {name} не аддитивный (допустимы {sorted(self._ADDITIVE)})")
            new += [(cmd, int(j), float(v)) for j, v in pairs]
        self._push_transient(new)
        return self

    def set_kicks(self, kicks=(), vkicks=()):
        """Задать доп. кики KICK/VKICK (гор./верт.) поверх текущих БЕЗ RESET: команда «задаёт», а не прибавляет,
        поэтому для смены корректоров (FOFB) пересборка кольца не нужна — только новые значения изменившихся
        (миллисекунды против ~0,3 с RESET у apply_perturb). Снять кик — задать 0. Возврат self."""
        self._push_transient([("KICK", int(j), float(v)) for j, v in kicks]
                             + [("VKICK", int(j), float(v)) for j, v in vkicks])
        return self

    def set_bend_scale(self, pairs=()):
        """Масштаб изгиба полюсов BSCALE j f (угол и грани = f × исходные, рамки заново) поверх текущих БЕЗ
        RESET: «задаёт», а не прибавляет (f=1 — исходный полюс). Живой зазор ондулятора: поле полюсов ∝ K,
        фокусировка граней ∝ K² и нулевая горизонталь получаются сами. Возврат self."""
        self._push_transient([("BSCALE", int(j), float(f)) for j, f in pairs])
        return self

    def _push_transient(self, new):
        """Отправить команды серверу и дописать их в транзитные (описание mags следует за ними). Длинный хвост
        сворачивается: аддитивные (DK/DKV/…) — суммой на (команда, j), «задающие» (KICK/VKICK) — последним."""
        if not new:
            return
        self._py = None
        self._replay(new)
        tr = list(self._transient) + list(new)
        if len(tr) > 4096:
            add, last = {}, {}
            for cmd, j, v in tr:
                if cmd in self._ADDITIVE.values():
                    add[(cmd, j)] = add.get((cmd, j), 0.0) + v
                else:
                    last[(cmd, j)] = v
            tr = [(c, j, v) for (c, j), v in last.items()] + [(c, j, v) for (c, j), v in add.items()]
        self._transient = tr

    def perturbed(self, kicks=(), yaws=(), vkicks=()):
        """Кольцо с несколькими возмущениями сразу: kicks=[(j, θ)] (гор.), vkicks=[(j, θ)] (верт.),
        yaws=[(j, α)] — один процесс сервера."""
        return self._spawn([("KICK", int(j), float(v)) for j, v in kicks]
                           + [("VKICK", int(j), float(v)) for j, v in vkicks]
                           + [("YAW", int(j), float(v)) for j, v in yaws])

    # ---- копии для параллельной работы (шаг 6) ---------------------------------------------------
    def clone(self):
        """Независимая копия кольца в ТЕКУЩЕМ состоянии: базовые (with_*) и транзитные (apply_perturb) возмущения
        и ВЧ те же, семантика RESET та же. ext — копия модели в памяти (~30 мс); pipe — новый процесс сервера
        с проигрыванием возмущений (~1 с). Для параллельных столбцов отклика, сканов, ансамблей."""
        r = copy.copy(self)
        r._perturb, r._transient = list(self._perturb), list(self._transient)
        r.rf = dict(self.rf) if self.rf else self.rf
        r._py = None
        if self.backend == "ext":
            r._be = _ExtBackend(_ext=self._be.x.clone())
            _LIVE.add(r)
        else:
            r._be = None
            r._start()
        return r

    def set_threads(self, n):
        """Потоков OpenMP на пакетный запрос этого кольца (0 — все). Только ext (у сервера — все потоки)."""
        if self.backend == "ext":
            self._be.x._r.threads = int(n)
        return self

    def junction(self, states, A, B):
        """Стык плоскостей C++-ядра: выход рамки A → вход рамки B (рамки — (p, axis, u, v))."""
        a = self._states(states)
        fr = np.concatenate([np.asarray(v, float).ravel() for v in (*A, *B)])
        if fr.size != 24:
            raise ValueError("рамка — четыре 3-вектора")
        return self._be.junction(a, fr)


@atexit.register
def _shutdown():
    for r in list(_LIVE):
        try:
            r.close()
        except Exception:
            pass
