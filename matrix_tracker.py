# -*- coding: utf-8 -*-
"""Матричный трекер SynchrotronSim — альтернативный движок на БД проекта. Каноническое фазовое
пространство (x, pₓ, y, p_y, t, δ): pₓ,p_y — нормированные на p₀ поперечные импульсы, δ=Δp/p₀,
t — путь (задел под ВЧ), s — независимая переменная. p_полный=p₀(1+δ) сохраняется в магнитах.

Модель элемента (ТЗ): прямоугольный ящик, центр=(x,y,z) из lattice_configs, ось по yaw, торцы =
центр ± L_eff/2·ось, плоскости входа/выхода ⊥ оси. Со-локованные обмотки → один магнит (сумма
физкоэф., суперпозиция). Кики — по ОПОРНОЙ жёсткости q/p₀: θ₀=(q/p₀)·B·L, k1=(q/p₀)G, k2=(q/p₀)S,
k3=(q/p₀)O. Дисперсия/хроматизм — из δ-корректного дрейфа. Чистый диполь = только кик изгиба (по B),
без слабой фокусировки (прямой ящик, равномерное поле не фокусирует).

Модифицированный пустой промежуток: геометрический пролёт от плоскости выхода пред. магнита до
плоскости входа след. (произвольная ориентация) + проверка направления импульса. Локальные орты:
u — горизонт. (в плоскости кольца), v — вертикаль; в состоянии x≡u, y≡v."""
import os, csv, math
import numpy as np
import bip_matrix as bm

TB = os.path.join("workspace", "tables")
QP = bm.QP                                   # q/p₀ (опорная жёсткость)
# СЕКТОРНЫЙ режим DQ (мультиполи вдоль изогнутой орбиты, dq_sector). Киковая форма валидирована
# (dq_sector_kick 0.65мм ≡ диффур), но перенос в ринг-геометрию ещё требует разбора рамки/масштаба —
# по умолчанию ВЫКЛ (осевой DQ, режимы 2/3 работают). Вкл: SECTOR_DQ=True или env SS_SECTOR_DQ=1.
SECTOR_DQ = os.environ.get("SS_SECTOR_DQ", "") == "1"


def _read(name):
    return list(csv.DictReader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";"))


_HW = None


def _hw_rows():
    """hardware_library по устройству (кэш) — для секторного профиля DQ (пути карт, nominal_arg)."""
    global _HW
    if _HW is None:
        _HW = {r["device_type_id"]: r for r in _read("hardware_library.csv")}
    return _HW


def load_db():
    hw = _read("hardware_library.csv")
    eff = {r["device_type_id"]: float(r["eff_length"]) for r in hw}
    phys = {r["device_type_id"]: float(r["phys_length"]) for r in hw}
    coup = {r["device_type_id"]: r for r in _read("coupling_matrix.csv")}
    trans = {r["device_type_id"]: r for r in _read("element_transforms.csv")}
    return eff, phys, coup, trans


_EDGEFF = None                                   # кэш форм-факторов краевого фринджа (по типу)


def _edge_ff(dev):
    """Форм-фактор краевой (фриндж) верт. фокусировки FE устройства из hardware_library (edge_ff).
    По умолчанию 0 (столбец отсутствует/пуст) ⇒ край выключен, старое поведение сохраняется."""
    global _EDGEFF
    if _EDGEFF is None:
        _EDGEFF = {r["device_type_id"]: float(r.get("edge_ff") or 0) for r in _read("hardware_library.csv")}
    return _EDGEFF.get(dev, 0.0)


def frame(yaw):
    axis = np.array([math.cos(yaw), math.sin(yaw), 0.0])
    uhat = np.array([-math.sin(yaw), math.cos(yaw), 0.0])
    vhat = np.array([0.0, 0.0, 1.0])
    return axis, uhat, vhat


def _drift(st, ell):
    """Точный дрейф на длину ell вдоль локальной оси (δ-корректный, канонический)."""
    x, px, y, py, t, d = st
    arg = (1.0 + d) ** 2 - px * px - py * py
    if arg <= 0.0:                                             # p⊥ > p_полный — частица потеряна
        raise ParticleLost(f"p⊥>p: px={px:+.3f} py={py:+.3f} δ={d:+.3f} (пучок раздут/потерян)")
    ps = math.sqrt(arg)
    return (x + ell * px / ps, px, y + ell * py / ps, py, t + ell * (1.0 + d) / ps, d)


class ParticleLost(Exception):
    pass


class Magnet:
    def __init__(self, insts, eff, phys, coup, trans, source="physical"):
        r0 = insts[0]
        self.center = np.array([float(r0["x"]), float(r0["y"]), float(r0["z"])])
        self.yaw = float(r0["yaw"])
        self.axis, self.uhat, self.vhat = frame(self.yaw)
        self.dev = r0["device_type_id"]
        self.Leff = eff[self.dev]                              # полевая (сила: кики по L_eff)
        self.L = phys[self.dev]                                # по металлу (геометрия ящика)
        self.entry = self.center - 0.5 * self.L * self.axis
        self.exit = self.center + 0.5 * self.L * self.axis
        self.kind = trans[self.dev]["kind"]; self.model = trans[self.dev]["model"]
        self.nsl = int(trans[self.dev]["n_slices"])
        self.devs = [r["device_type_id"] for r in insts]
        # ДАТЧИК ПОЛОЖЕНИЯ ПУЧКА (BPM): не преобразует частицы (пропускает как дрейф L=0),
        # регистрирует центроид пучка. sensor_id = instance_id экземпляра в решётке.
        self.monitor = (self.kind == "monitor"); self.sensor_id = r0.get("instance_id", "")
        # суммарные ФИЗКОЭФ. (суперпозиция обмоток). source: 'physical' — из mat_* столбцов БД;
        # 'current' — из токов через coupling_matrix (когда матрицы связи готовы/валидны).
        B = G = Gv = S = O = 0.0                          # Gv — ВЕРТ. градиент (∂By/∂v), отдельно от G=∂Bz/∂u
        Bs = Gs = Ss = Os = 0.0                           # СКЬЮ мультиполи: верт.диполь(By на оси), скью квад/секст/окт
        for r in insts:
            if source == "physical":
                B += float(r.get("mat_B") or 0); G += float(r.get("mat_G") or 0)
                Gv += float(r.get("mat_Gv") or r.get("mat_G") or 0)   # дефолт = mat_G ⇒ гармонич. (совместимо)
                S += float(r.get("mat_S") or 0); O += float(r.get("mat_O") or 0)
                Bs += float(r.get("mat_Bs") or 0); Gs += float(r.get("mat_Gs") or 0)   # скью (дефолт 0 ⇒ старое поведение)
                Ss += float(r.get("mat_Ss") or 0); Os += float(r.get("mat_Os") or 0)
            else:
                c = coup[r["device_type_id"]]
                inp = float(r["arg_val"] or 0) if c["input_type"] == "current" else float(r["correct_coeff"] or 0)
                gcoef = float(c["per_G_Tlm"]) * inp
                B += float(c["per_B_Tl"]) * inp; G += gcoef; Gv += gcoef   # current-source: гармонич. (Gv=G)
                S += float(c["per_S_Tlm2"]) * inp; O += float(c["per_O_Tlm3"]) * inp
                Bs += float(c.get("per_Bs_Tl") or 0) * inp; Gs += float(c.get("per_Gs_Tlm") or 0) * inp   # скью из матрицы связи
                Ss += float(c.get("per_Ss_Tlm2") or 0) * inp; Os += float(c.get("per_Os_Tlm3") or 0) * inp
        # ДИПОЛЬ: mat_G/Gv/S/O заданы ОТНОСИТЕЛЬНО поля (доли mat_B, аппаратные константы) → в абсолютные.
        # Даёт 1/f = θ₀·g (фокус ∝ изгибу): при излучении θ₀ пересчит. по текущему p, отношения фиксированы.
        if self.kind == "dipole" and source == "physical":
            G *= B; Gv *= B; S *= B; O *= B; Gs *= B; Ss *= B; Os *= B   # скью-градиенты тоже относит.; Bs — абс. изгиб
        # силы по ОПОРНОЙ жёсткости (по L_eff — интеграл поля); δ-зависимость — в дрейфе.
        # k1=k1u (гор. фокус ∂Bz/∂u), k1v (верт. фокус ∂By/∂v): неголоморфн. диполь Gu≠Gv ⇒ |f_x|≠|f_y|.
        self.theta0 = QP * B * self.Leff
        self.k1 = QP * G; self.k1v = QP * Gv; self.k2 = QP * S; self.k3 = QP * O
        # СКЬЮ + ВЕРТ.ИЗГИБ (универс. мультипольный элемент): thetaY — верт. изгиб (skew-диполь, поле By на оси),
        # k1s — скью-квад (связь плоскостей), k2s/k3s — скью секступоль/октуполь. Мнимая часть Cₙ=bₙ+i·aₙ.
        self.thetaY = QP * Bs * self.Leff
        self.k1s = QP * Gs; self.k2s = QP * Ss; self.k3s = QP * Os
        # краевой (фриндж) верт. фокус 2-го порядка: 1/f_y = FE·θ₀²/L_eff (∝поле², offset/угол-независим).
        # Эквивалент углового члена Δv'=−(q/p)u'v∂ₓBz после ∫по частям = −(q/p)²∫Bz²ds·v. Диполи: θ₀≠0.
        self.kEdge = _edge_ff(self.dev) * self.theta0 ** 2 / self.Leff if self.Leff > 0 else 0.0
        self.nonzero = any(abs(v) > 0 for v in (self.theta0, self.thetaY, self.k1, self.k1v, self.k2, self.k3,
                                                self.k1s, self.k2s, self.k3s, self.kEdge))
        # СЕКТОРНЫЙ режим для combined-function DQ в ТОКОВОМ режиме: осевые мультиполи неверны на
        # сагитте (пучок на изогнутой дуге). Профиль b0(X),Gu(X),Gv(X),S(X) ВДОЛЬ орбиты из карт×токи
        # (dq_sector, 0.65мм ≡ диффур vs 4.3мм осевые). В физкоэф.-режиме карт/токов нет → осевой.
        self.sector = False; self.sec = None
        if SECTOR_DQ and self.kind == "dq" and source == "current":
            try:
                import dq_sector_profile as _dsp
                cur = [float(r.get("arg_val") or 0) for r in insts]
                if any(abs(c) > 1e-9 for c in cur):
                    self.sec = _dsp.build(self.devs, cur, _hw_rows())
                    self.sector = True
                    self.sec_xref = self._sector_ref()          # опорная дуга в сжатом железе L (для du)
            except Exception:
                self.sector = False

    def body(self, st):
        """Тело: N срезов drift-kick-drift. ДРЕЙФ по геометрии (металл, ds_g=phys/N), КИКИ по полю
        (ds_f=eff/N) ⇒ полный изгиб/фокус = k·L_eff (верно), а ящик = металл (верно). Чистый диполь =
        только кик изгиба (без фокусировки); квадруполь/секступоль/октуполь — линейн./нелин. кики."""
        if self.sector:                                        # combined-function DQ вдоль изогнутой орбиты
            return self._body_sector(st)
        if not self.nonzero:                                   # выключенный — дрейф по металлу
            return _drift(st, self.L)
        th, k1, k1v, k2, k3, kE = self.theta0, self.k1, self.k1v, self.k2, self.k3, self.kEdge
        thy, k1s, k2s, k3s = self.thetaY, self.k1s, self.k2s, self.k3s
        N = self.nsl if self.model == "composite" else 40
        dsg = self.L / N; dsf = self.Leff / N
        for _ in range(N):
            st = _drift(st, dsg / 2.0)
            x, px, y, py, t, d = st
            # НОРМАЛЬНЫЕ мультиполи (Re Cₙ): изгиб θ0, квад −k1·x, секступоль, октуполь + краевая линза
            px += th / N - k1 * x * dsf - 0.5 * k2 * (x * x - y * y) * dsf - (k3 / 6.0) * (x ** 3 - 3 * x * y * y) * dsf
            py += k1v * y * dsf + k2 * x * y * dsf + (k3 / 6.0) * (3 * x * x * y - y ** 3) * dsf - (kE / N) * y
            # СКЬЮ + ВЕРТ.ИЗГИБ (Im Cₙ): верт.диполь θy, скью-квад k1s (связь плоскостей), скью секступоль/октуполь
            if k1s or k2s or k3s or thy:
                px += k1s * y * dsf + k2s * x * y * dsf + (k3s / 6.0) * (3 * x * x * y - y ** 3) * dsf
                py += thy / N + k1s * x * dsf + 0.5 * k2s * (x * x - y * y) * dsf + (k3s / 6.0) * (x ** 3 - 3 * x * y * y) * dsf
            st = _drift((x, px, y, py, t, d), dsg / 2.0)
        return st

    def _body_sector(self, st):
        """СЕКТОРНОЕ тело DQ: всё (дрейфы+кики) укладывается в ЖЕЛЕЗО L (карта поля — габарит диффур-
        модели, к кику отношения не имеет). Поле карты (зона 2·half) СЖИМАЕТСЯ в L с сохранением
        интеграла: дрейф ds_i=L/N, но кик по ds_f=2·half/N ⇒ ∫кик = ∫поля (изгиб/фокус эквивалентны).
        Мультиполи ВДОЛЬ орбиты (du=x−x_ref, x_ref — самосогласованная дуга в сжатом L). Форма
        валидирована dq_sector_kick (0.65мм ≡ диффур). x≡u гор., y≡v верт."""
        p = self.sec; half = p["half"]; Xg, b0, Gu, Gv, S = p["Xg"], p["b0"], p["Gu"], p["Gv"], p["S"]
        xref = self.sec_xref; N = len(xref)
        ds_i = self.L / N; ds_f = 2.0 * half / N               # дрейф по железу L, кик по зоне поля 2·half
        for i in range(N):
            st = _drift(st, ds_i / 2.0)
            x, px, y, py, t, d = st
            Xf = -half + (i + 0.5) * ds_f                       # позиция в профиле (полная зона поля)
            du = x - xref[i]                                    # отклонение от изогнутой орбиты (в сжатом L)
            b = np.interp(Xf, Xg, b0); gu = np.interp(Xf, Xg, Gu)
            gv = np.interp(Xf, Xg, Gv); s = np.interp(Xf, Xg, S)
            Bz = b + gu * du + s * (du * du - y * y); By = gv * y + 2 * s * du * y
            px += -QP * Bz * ds_f; py += QP * By * ds_f        # кик по полю (интеграл сохранён), знаки matrix_tracker
            st = _drift((x, px, y, py, t, d), ds_i / 2.0)
        return st

    def _sector_ref(self):
        """Опорная дуга x_ref[i] осевой частицы (старт x=0,px=0) в СЖАТОМ железе L: видит только
        диполь b0 (du=0). x_ref[i] — позиция в точке кика среза i (для du=x−x_ref в _body_sector)."""
        p = self.sec; half = p["half"]; Xg, b0 = p["Xg"], p["b0"]
        N = max(self.nsl if self.model == "composite" else 40, 60)
        ds_i = self.L / N; ds_f = 2.0 * half / N
        xr = 0.0; pxr = 0.0; out = []
        for i in range(N):
            ps = math.sqrt(max(1e-12, 1.0 - pxr * pxr)); xr += pxr * ds_i / (2.0 * ps)
            out.append(xr)                                     # позиция опорной частицы в точке кика
            pxr += -QP * float(np.interp(-half + (i + 0.5) * ds_f, Xg, b0)) * ds_f
            ps = math.sqrt(max(1e-12, 1.0 - pxr * pxr)); xr += pxr * ds_i / (2.0 * ps)
        return np.array(out)


def build_line(lattice_id, eff, phys, coup, trans, lo=None, hi=None, source="physical"):
    conf = {(r["lattice_id"], r["instance_id"]): r for r in _read("lattice_configs.csv")}
    conn = [r for r in _read("element_connectivity.csv") if r["lattice_id"] == str(lattice_id)]
    conn.sort(key=lambda r: int(r["instance_id"]))
    mags = []; i = 0
    while i < len(conn):
        r = conf[(str(lattice_id), conn[i]["instance_id"])]; group = [r]
        while i + 1 < len(conn):
            rn = conf[(str(lattice_id), conn[i + 1]["instance_id"])]
            if (rn["x"], rn["y"], rn["z"]) == (r["x"], r["y"], r["z"]):
                group.append(rn); i += 1
            else:
                break
        mags.append(Magnet(group, eff, phys, coup, trans, source)); i += 1
    return mags[lo:hi] if lo is not None else mags


def modified_drift(st, A, B):
    """Геометрический пролёт: плоскость выхода A → плоскость входа B (произв. ориентации). δ-корректно
    (p_полный сохраняется). Проверка направления импульса. st в системе A → st в системе B, + путь."""
    x, px, y, py, t, d = st
    P = A.exit + x * A.uhat + y * A.vhat                       # глоб. позиция на выходной плоскости A
    ps = math.sqrt((1.0 + d) ** 2 - px * px - py * py)
    dir_ = (px * A.uhat + py * A.vhat + ps * A.axis) / (1.0 + d)   # единичное направление импульса (глоб.)
    denom = float(B.axis @ dir_)
    if denom <= 1e-9:
        raise ParticleLost(f"пучок не долетит до входа {B.dev}: n·d={denom:.3g} ≤ 0")
    tt = float(B.axis @ (B.entry - P)) / denom
    if tt < 0:
        raise ParticleLost(f"плоскость входа {B.dev} позади частицы: t={tt:.3g} < 0")
    Q = P + tt * dir_
    rel = Q - B.entry
    x2 = float(rel @ B.uhat); y2 = float(rel @ B.vhat)
    px2 = (1.0 + d) * float(dir_ @ B.uhat); py2 = (1.0 + d) * float(dir_ @ B.vhat)
    return (x2, px2, y2, py2, t + tt, d)


def track(particles, mags, want_path=False, monitors=None):
    """Прогон частиц. Возврат [(final_state, lost)] — lost=None если долетела, иначе (индекс магнита, причина).
    Если monitors (dict) передан — накапливает {sensor_id: [(x,y,t) по частицам]} на входе каждого
    элемента-датчика (BPM): датчик пропускает частицу (body=дрейф L=0), но фиксирует её (x,y,t)."""
    finals, paths = [], []
    for st in particles:
        s = 0.0; path = [(0.0, st[0], st[2])]; lost = None; at = 0
        try:
            if monitors is not None and mags[0].monitor:
                monitors.setdefault(mags[0].sensor_id, []).append(tuple(st[:6]))   # полный 6D-стейт частицы
            st = mags[0].body(st); s += mags[0].L; path.append((s, st[0], st[2]))
            for j in range(1, len(mags)):
                at = j
                st = modified_drift(st, mags[j - 1], mags[j]); s = st[4]; path.append((s, st[0], st[2]))
                if monitors is not None and mags[j].monitor:                    # вход датчика: зафиксировать частицу
                    monitors.setdefault(mags[j].sensor_id, []).append(tuple(st[:6]))
                st = mags[j].body(st); s += mags[j].L; path.append((s, st[0], st[2]))
        except ParticleLost as e:
            lost = (at, str(e))
        finals.append((st, lost))
        if want_path:
            paths.append(np.array(path))
    return (finals, paths) if want_path else finals


# ---- Датчики положения пучка (BPM): накопление в ОЗУ, сброс в БД ПОСЛЕ моделирования ---------------
# Фишки индустриальных кодов (Elegant MONITOR/WATCH, MAD-NG): центроид + МОМЕНТЫ пучка (σx,σy),
# TURN-BY-TURN, модель НЕСОВЕРШЕНСТВА (смещение dx/dy, наклон tilt, усиление gain, setpoint).
BPM_COLS = ["id_exp", "id_sensor", "turn", "x", "y", "xp", "yp", "t", "delta", "sigma_x", "sigma_y", "n_particles"]


def load_bpm_config(lattice_id=None, path=None):
    """Несовершенства датчиков из bpm_config.csv → {id_sensor: {dx,dy,tilt,gain_x,gain_y,setpoint_x,setpoint_y}}.
    Нет файла/строки ⇒ идеальный датчик. Фильтр по lattice_id (пусто в строке = на все решётки)."""
    if path is None:
        path = os.path.join(TB, "bpm_config.csv")
    cfg = {}
    if not os.path.exists(path):
        return cfg
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";"):
        if lattice_id is not None and (r.get("lattice_id") or "") not in ("", str(lattice_id)):
            continue
        cfg[r["id_sensor"]] = {k: float(r.get(k) or (1.0 if k.startswith("gain") else 0.0))
                               for k in ("dx", "dy", "tilt", "gain_x", "gain_y", "setpoint_x", "setpoint_y")}
    return cfg


class BPMRegistry:
    """Регистр датчиков: копит показания в ОЗУ (по оборотам), сбрасывает в БД одним flush() ПОСЛЕ прогона.
    add_turn: центроид (среднее) + σ (RMS) по частицам, к (x,y) применяется модель несовершенства.
    config = {id_sensor: {dx,dy,tilt,gain_x,gain_y,setpoint_x,setpoint_y}} (по умолч. идеальный)."""

    def __init__(self, config=None):
        self.config = config or {}
        self.records = []                                      # буфер в ОЗУ: кортежи BPM_COLS[1:]

    @staticmethod
    def _distort(mx, my, c):
        """Показание из истинного центроида (как Elegant/MAD): усиление·(центроид−смещение), поворот tilt, −setpoint."""
        cx = c.get("gain_x", 1.0) * (mx - c.get("dx", 0.0))
        cy = c.get("gain_y", 1.0) * (my - c.get("dy", 0.0))
        ct, s3 = math.cos(c.get("tilt", 0.0)), math.sin(c.get("tilt", 0.0))
        return cx * ct + cy * s3 - c.get("setpoint_x", 0.0), cy * ct - cx * s3 - c.get("setpoint_y", 0.0)

    def add_turn(self, monitors, turn=0):
        """Показания одного оборота в ОЗУ: monitors={id_sensor:[(x,px,y,py,t,d) по частицам]}."""
        for sid, states in monitors.items():
            a = np.array(states, dtype=float)                  # (n,6): x,px,y,py,t,δ
            c = a.mean(axis=0); n = len(a)                     # центроид + число частиц
            x, y = self._distort(c[0], c[2], self.config.get(str(sid), {}))
            self.records.append((str(sid), turn, x, y, c[1], c[3], c[4], c[5],
                                 float(a[:, 0].std()), float(a[:, 2].std()), n))   # + σx,σy (RMS размер)

    def flush(self, id_exp, path=None):
        """Сброс всего буфера ОЗУ → bpm_readings.csv (строки того же id_exp заменяются)."""
        if path is None:
            path = os.path.join(TB, "bpm_readings.csv")
        keep = []
        if os.path.exists(path):
            rd = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
            keep = [r for r in rd[1:] if r and r[0] != str(id_exp)]      # чужие эксперименты сохранить
        with open(path, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";"); w.writerow(BPM_COLS); w.writerows(keep)
            for rec in self.records:
                w.writerow([str(id_exp)] + [f"{v:.9g}" if isinstance(v, float) else v for v in rec])
        return len(self.records)


def record_bpm(lattice_id, particles, id_exp, source="physical", n_turns=1, config=None):
    """Полный прогон: решётка → пучок n_turns оборотов (выжившие → следующий оборот) → показания в ОЗУ →
    СБРОС в БД после последнего оборота. Возврат (registry, finals)."""
    eff, phys, coup, trans = load_db()
    mags = build_line(lattice_id, eff, phys, coup, trans, source=source)
    reg = BPMRegistry(config if config is not None else load_bpm_config(lattice_id))
    parts = [list(p) for p in particles]; finals = None
    for turn in range(n_turns):
        mon = {}
        finals = track(parts, mags, monitors=mon); reg.add_turn(mon, turn)
        parts = [list(st) for st, lost in finals if lost is None]    # выжившие → следующий оборот
        if not parts:
            break
    reg.flush(id_exp)
    return reg, finals


if __name__ == "__main__":
    import sys
    lat = int(sys.argv[1]) if len(sys.argv) > 1 else 101
    lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
    hi = int(sys.argv[3]) if len(sys.argv) > 3 else 12
    src = sys.argv[4] if len(sys.argv) > 4 else "physical"     # physical | current
    eff, phys, coup, trans = load_db()
    mags = build_line(lat, eff, phys, coup, trans, lo, hi, src)
    print(f"Решётка {lat}, магниты [{lo}:{hi}], источник='{src}', q/p₀={QP:.3f}. Вектор (x,pₓ,y,p_y,t,δ):")
    for m in mags:
        print(f"  {'+'.join(m.devs):26} L={m.L*1e3:5.0f}мм yaw={m.yaw*1e3:+6.1f}мрад  "
              f"θ₀={m.theta0*1e3:+8.2f}мрад k1={m.k1:+7.3f} k2={m.k2:+7.3g} k3={m.k3:+7.3g}")
    # опорная (на оси) + малые смещения 0.1мм + дисперсия (δ=+1% на оси)
    P = [(0, 0, 0, 0, 0, 0), (1e-4, 0, 0, 0, 0, 0), (0, 0, 1e-4, 0, 0, 0), (0, 0, 0, 0, 0, 0.01)]
    lbl = ["опорная       ", "x0=0.1мм      ", "y0=0.1мм      ", "на оси δ=+1%  "]
    print("\nФинал (x,pₓ,y,p_y мм/мрад, δ):")
    for l, (st, lost) in zip(lbl, track(P, mags)):
        if lost:
            print(f"  {l}: ПОТЕРЯНА на магните #{lost[0]} — {lost[1]}")
        else:
            print(f"  {l}: x={st[0]*1e3:+8.4f} pₓ={st[1]*1e3:+8.4f}  y={st[2]*1e3:+8.4f} p_y={st[3]*1e3:+8.4f}  δ={st[5]:+.3f}")
