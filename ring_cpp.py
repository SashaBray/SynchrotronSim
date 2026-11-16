# -*- coding: utf-8 -*-
"""Кольцо на продакшн-ядре C++ (MatrixTracker --ring-serve) с интерфейсом ring_analysis.Ring.

Трекинг — только в C++: оборот, пооборотные координаты, выходы магнитов, срезы изгибов считает сервер
кольца (src/matrix_ring_serve.h), частицы параллельно (OpenMP). Python здесь — лишь транспорт и
анализ. Для аналитического механизма (ring_matrix, ring_response) нужны описания магнитов — они берутся
из той же решётки БД питоновским загрузчиком по требованию (атрибуты mags/plane), без трекинга.

  CppRing(lattice, start_before_first_m, sr, isr, rf, design_mev, nsl, order, seed)
    track(states, n_turns)            → (n, 6), NaN у потерянных
    track_elements(states)            → (n, n_mag+1, 6)
    track_tbt(states, n_turns, limit) → (n, n_turns+1, 6)
    track_pieces(states)              → (n, n_pieces+1, 6) — после каждого среза изгиба
    with_kick(j, θ) / with_yaw(j, α)  → новое кольцо с возмущением (C++ и описание магнитов)
    perturbed(kicks, yaws)            → то же для нескольких возмущений
    apply_perturb(kicks, …, dk3s)     → транзитные возмущения на том же сервере (RESET + набор)
    junction(states, A, B)            → стык двух произвольных плоскостей (проверка матриц стыков)

Возмущения двух сортов: «базовые» (with_*/perturbed — часть кольца, переживают RESET, копии их наследуют) и
транзитные (apply_perturb — заменяются его следующим вызовом). Описание магнитов mags/plane всегда по обоим.

Протокол сервера — в шапке src/matrix_ring_serve.h.
"""
import os, copy, math, subprocess, atexit, weakref
import numpy as np

ROOT = os.path.dirname(os.path.abspath(__file__))
EXE = os.path.join(ROOT, "build", "MatrixTracker.exe")
MC2_MEV = 0.51099895069
_LIVE = weakref.WeakSet()


def _default_slices():
    """Срезы и порядок по умолчанию — как у Python-ядра: env SS_SECTOR_NSL/SS_SECTOR_ORDER или 5 и 4."""
    return (int(os.environ.get("SS_SECTOR_NSL", "") or 5), int(os.environ.get("SS_SECTOR_ORDER", "") or 4))


class CppRing:
    label = "SynchrotronSim"

    def __init__(self, lattice=1, start_before_first_m=None, settings=None, sr=False, isr=False, rf=None,
                 design_mev=6000.0, nsl=None, order=None, seed=12345, ws=None, timeout=None):
        dn, do = _default_slices()
        self.lattice = int(lattice)
        self.lead = start_before_first_m
        self.sr, self.isr, self.rf, self.design_mev = bool(sr or isr), bool(isr), rf, float(design_mev)
        self.nsl = int(nsl or dn); self.order = int(order or do); self.seed = int(seed)
        self.ws = ws                            # путь к workspace с ВОЗМУЩЁННОЙ решёткой (ансамбли ошибок); None — обычная БД
        self.timeout = timeout                  # сторож: команда дольше timeout с → убить процесс (защита от зависаний)
        self.settings = settings
        self._perturb = []                      # базовые [(cmd, j, value)] — with_*/perturbed: переживают RESET
        self._transient = []                    # транзитные — от последнего apply_perturb
        self._py = None
        self._start()

    # ---- процесс сервера --------------------------------------------------------------------------
    def _start(self):
        if not os.path.exists(EXE):
            raise FileNotFoundError(f"нет {EXE}: соберите cmake --build build --target MatrixTracker")
        rf = self.rf or {}
        args = [EXE, "--ring-serve", str(self.lattice),
                "nan" if self.lead is None else repr(float(self.lead)), repr(self.design_mev),
                "1" if self.sr else "0", "1" if self.isr else "0", str(self.nsl), str(self.order),
                "1" if self.rf else "0", repr(float(rf.get("volt_mv", 0.0))), repr(float(rf.get("harmonic", 0.0))),
                repr(float(rf.get("phase_deg", 0.0))), str(self.seed)]
        if self.ws:                             # возмущённая рабочая папка (ансамбли ошибок) — 13-й аргумент
            args.append(self.ws)
        self._proc = subprocess.Popen(args, cwd=ROOT, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                      stderr=subprocess.DEVNULL, bufsize=0)
        head = self._proc.stdout.readline().decode("utf-8", "replace").split()
        if not head or head[0] != "READY":
            raise RuntimeError(f"сервер кольца не поднялся: {' '.join(head)}")
        self.n_mag, self.circumference, self.n_pieces = int(head[1]), float(head[2]), int(head[3])
        info = self._text("INFO")
        self.s_in = np.array([float(l.split()[1]) for l in info])
        self.s_out = np.array([float(l.split()[2]) for l in info])
        self.kinds = [l.split()[4] if l.split()[4] != "-" else "" for l in info]
        self.names = [l.split()[5] for l in info]
        self._replay(self._perturb + self._transient)
        _LIVE.add(self)

    def _replay(self, perturb):
        for cmd, j, v in perturb:
            self._text(f"{cmd} {j} {v!r}")

    def close(self):
        p = getattr(self, "_proc", None)
        if p and p.poll() is None:
            try:
                p.stdin.write(b"QUIT\n"); p.stdin.flush(); p.wait(timeout=5)
            except Exception:
                p.kill()
        self._proc = None

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def _watch(self):
        """Сторож: таймер, убивающий процесс через self.timeout с, если команда зависла (защита ансамбля
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

    # ---- интерфейс трекера (как ring_analysis.Ring) -------------------------------------------------
    @staticmethod
    def _states(states):
        a = np.asarray(states, float)
        return a.reshape(-1, 6)

    def track(self, states, n_turns=1):
        a = self._states(states)
        return self._doubles(f"TRACK {len(a)} {int(n_turns)}", a).reshape(len(a), 6)

    def track_tbt(self, states, n_turns, limit=None):
        a = self._states(states)
        lim = float(limit) if limit else 0.0
        return self._doubles(f"TBT {len(a)} {int(n_turns)} {lim!r}", a).reshape(len(a), int(n_turns) + 1, 6)

    def track_elements(self, states):
        a = self._states(states)
        return self._doubles(f"ELEM {len(a)}", a).reshape(len(a), self.n_mag + 1, 6)

    def track_pieces(self, states):
        a = self._states(states)
        return self._doubles(f"PIECES {len(a)}", a).reshape(len(a), self.n_pieces + 1, 6)

    def piece_meta(self):
        """[(магнит, стык перед куском, изгиб, Leff, θ0, k1, k2, k3)] — срезы изгибов и прочие магниты."""
        out = []
        for l in self._text("PIECEMETA"):
            v = l.split()
            out.append((int(v[0]), v[1] == "1", v[2] == "1", *map(float, v[3:8])))
        return out

    def beta(self, delta):
        pc = self.design_mev * (1.0 + delta)
        return pc / math.sqrt(pc * pc + MC2_MEV ** 2)

    # ---- описание магнитов для аналитического механизма (без трекинга) ----------------------------
    def _python_ring(self):
        if self._py is None:
            import ring_analysis as ra
            keys = ("SS_SECTOR_NSL", "SS_SECTOR_ORDER")
            saved = {k: os.environ.get(k) for k in keys}
            os.environ.update({"SS_SECTOR_NSL": str(self.nsl), "SS_SECTOR_ORDER": str(self.order)})
            try:
                r = ra.Ring(self.lattice, start_before_first_m=self.lead, settings=self.settings,
                            design_mev=self.design_mev)
            finally:
                for k, v in saved.items():
                    os.environ.pop(k, None) if v is None else os.environ.update({k: v})
            r.mags = self._describe(r.mags, self._perturb + self._transient)
            if abs(r.circumference - self.circumference) > 1e-9 or len(r.mags) != self.n_mag:
                raise RuntimeError("геометрия Python-описания и C++-ядра разошлась")
            self._py = r
        return self._py

    _FIELD = {"DK": ("k1", "k1v"), "DSK": ("ks",), "DK2": ("k2",), "DK3": ("k3",)}

    def _describe(self, base, perturb):
        """Магниты Python-описания с возмущениями сервера (копии; base не трогается) — та же физика, что у
        Element в C++ (matrix_ring_serve.h): KICK/VKICK задают доп. кики hkick/vkick (последний выигрывает),
        DK/DSK/DK2/DK3 прибавляют к k1(+k1v)/ks/k2/k3, YAW поворачивает рамки (ra.yaw_element).
        Тело — ПО ПОЛЮ, как Element::precache: элемент с появившимся градиентом/мультиполем идёт кик-телом.
        Гор. кик «ящика» (не изгиб, без поля — корректор) = его θ0 box-тела, как ra.kicked; у изгиба и полевого
        элемента θ0 не трогается (у изгиба это кривизна h для ring_matrix), кик лишь записан в hkick/vkick.
        Верт. кика в Python-эталоне нет: vkick — только атрибут (ring_matrix кики корректоров не учитывает)."""
        import ring_analysis as ra
        import matrix_tracker as mt
        mags = list(base); touched = set()
        for cmd, j, v in perturb:
            m = copy.copy(mags[j])
            if cmd == "YAW":
                m = ra.yaw_element(m, v)
            elif cmd in ("KICK", "VKICK"):
                setattr(m, "hkick" if cmd == "KICK" else "vkick", v)
            else:
                for a in self._FIELD[cmd]:
                    setattr(m, a, float(getattr(m, a, 0.0) or 0.0) + v)
            mags[j] = m; touched.add(j)
        for j in touched:
            m = mags[j]
            if not m.sector:
                m.kick_sector = (mt.SECTOR and m.Leff > 0
                                 and any(abs(float(getattr(m, a, 0.0) or 0.0)) > 0 for a in ("k1", "k1v", "k2", "k3")))
                if m.kick_sector and not hasattr(m, "_nsl_sector"):      # нулевой элемент получил поле
                    m._nsl_sector = self.nsl; m._yoshida_dk = mt._yoshida(self.order)
                hk = float(getattr(m, "hkick", 0.0) or 0.0)
                if hk and not m.kick_sector:                             # «ящик» корректора — как ra.kicked
                    m.theta0 = base[j].theta0 + hk; m.xref = None
            m.nonzero = any(abs(float(getattr(m, a, 0.0) or 0.0)) > 0 for a in
                            ("theta0", "k1", "k1v", "k2", "k3", "kEdge", "kEdge1", "kEdge2", "wkfoc"))
        return mags

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

    def apply_perturb(self, kicks=(), vkicks=(), dks=(), dsks=(), dk2s=(), dk3s=()):
        """Применить возмущения НА ТОМ ЖЕ сервере (RESET + KICK/VKICK/DK/DSK/DK2/DK3), без нового процесса. Для
        матрицы отклика, коррекции орбиты, BBA, коррекции связи и нелинейного тюнера. RESET пересобирает кольцо
        из ws (ошибки рассогласования СОХРАНЯЮТСЯ), затем возвращаются базовые возмущения кольца (with_*/
        perturbed) — сбрасываются только транзитные от прежнего apply_perturb. dks=[(j,dk)] — k1 квадруполя
        (BBA); dsks=[(j,ks)] — скью (связь); dk2s=[(j,dk2)]/dk3s=[(j,dk3)] — секступ./октуп. силы (нелинейный
        тюнер ν(δ)). Описание mags/plane следует за ними. apply_perturb() — вернуть кольцо в исходное. Возврат self."""
        new = [(cmd, int(j), float(v)) for cmd, pairs in (("KICK", kicks), ("VKICK", vkicks), ("DK", dks),
                                                          ("DSK", dsks), ("DK2", dk2s), ("DK3", dk3s))
               for j, v in pairs]
        self._py = None
        self._transient = []                    # до RESET: при сбое описание не припишет кольцу чужих возмущений
        self._text("RESET")
        self._replay(self._perturb)
        self._replay(new)
        self._transient = new
        return self

    def perturbed(self, kicks=(), yaws=(), vkicks=()):
        """Кольцо с несколькими возмущениями сразу: kicks=[(j, θ)] (гор.), vkicks=[(j, θ)] (верт.),
        yaws=[(j, α)] — один процесс сервера."""
        return self._spawn([("KICK", int(j), float(v)) for j, v in kicks]
                           + [("VKICK", int(j), float(v)) for j, v in vkicks]
                           + [("YAW", int(j), float(v)) for j, v in yaws])

    def junction(self, states, A, B):
        """Стык плоскостей C++-ядра: выход рамки A → вход рамки B (рамки — (p, axis, u, v))."""
        a = self._states(states)
        fr = np.concatenate([np.asarray(v, float).ravel() for v in (*A, *B)])
        if fr.size != 24:
            raise ValueError("рамка — четыре 3-вектора")
        m = self._send(f"JUNCTION {len(a)}", np.concatenate([fr, a.ravel()]))
        buf = bytearray()
        while len(buf) < m * 8:
            chunk = self._proc.stdout.read(m * 8 - len(buf))
            if not chunk:
                raise RuntimeError("сервер кольца оборвал ответ")
            buf += chunk
        return np.frombuffer(bytes(buf), dtype="<f8").reshape(len(a), 6)


@atexit.register
def _shutdown():
    for r in list(_LIVE):
        try:
            r.close()
        except Exception:
            pass
