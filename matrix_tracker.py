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


def _read(name):
    return list(csv.DictReader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";"))


def load_db():
    hw = _read("hardware_library.csv")
    eff = {r["device_type_id"]: float(r["eff_length"]) for r in hw}
    phys = {}                                        # phys_length УДАЛЁН из библиотеки — кик работает ТОЛЬКО по L_eff
    coup = {}                                        # coupling_matrix УДАЛЁН (токовый режим вырезан) — ядро физрежима её не читает
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
    def __init__(self, insts, eff, phys, coup, trans, source="physical", settings=None):
        r0 = insts[0]
        self.center = np.array([float(r0["x"]), float(r0["y"]), float(r0["z"])])
        self.yaw = float(r0["yaw"])
        self.axis, self.uhat, self.vhat = frame(self.yaw)
        self.dev = r0["device_type_id"]
        # ЕДИНЫЙ ИСТОЧНИК — структурные метаданные из lattice_configs (eff_length/kind/n_slices/edge1/edge2).
        # Fallback на hardware_library/element_transforms — только для неполной БД (те же значения).
        self.Leff = float(r0.get("eff_length") or eff.get(self.dev) or 0.0)   # эффективная (проектная .lte) длина
        self.L = self.Leff                                     # ОТКАЗ от металл-длины: кик работает ТОЛЬКО по L_eff
        self.entry = self.center - 0.5 * self.L * self.axis
        self.exit = self.center + 0.5 * self.L * self.axis
        _tr = trans.get(self.dev, {})
        self.kind = r0.get("kind") or _tr.get("kind") or ""
        # срезов тела (drift-kick-drift): приоритет у элемента, затем у типа, затем параметр прогона
        self.nsl = int(r0.get("n_slices") or _tr.get("n_slices") or (settings or Settings()).n_slices)
        self.edge1 = float(r0.get("edge1") or 0.0)             # краевые углы E1/E2 (.lte); физика края — Фаза 3
        self.edge2 = float(r0.get("edge2") or 0.0)
        self.devs = [r["device_type_id"] for r in insts]
        # ДАТЧИК ПОЛОЖЕНИЯ ПУЧКА (BPM): не преобразует частицы (пропускает как дрейф L=0),
        # регистрирует центроид пучка. sensor_id = instance_id экземпляра в решётке.
        self.monitor = (self.kind == "monitor"); self.sensor_id = r0.get("instance_id", "")
        # суммарные ФИЗКОЭФ. (суперпозиция обмоток). source: 'physical' — из mat_* столбцов БД;
        # 'current' — из токов через coupling_matrix (когда матрицы связи готовы/валидны).
        B = G = Gv = S = O = 0.0                          # Gv — ВЕРТ. градиент (∂By/∂v), отдельно от G=∂Bz/∂u
        if source != "physical":
            raise ValueError("токовый режим матричного ядра архивирован — используйте source='physical'; "
                             "физкоэф. mat_* берутся из lattice_configs (калибровки предвычисляются туда при сборке БД)")
        for r in insts:
            B += float(r.get("mat_B") or 0); G += float(r.get("mat_G") or 0)
            Gv += float(r.get("mat_Gv") or r.get("mat_G") or 0)   # дефолт = mat_G ⇒ гармонич. (совместимо)
            S += float(r.get("mat_S") or 0); O += float(r.get("mat_O") or 0)
        # ДИПОЛЬ: mat_G/Gv/S/O заданы ОТНОСИТЕЛЬНО поля (доли mat_B, аппаратные константы) → в абсолютные.
        # Даёт 1/f = θ₀·g (фокус ∝ изгибу): при излучении θ₀ пересчит. по текущему p, отношения фиксированы.
        if self.kind == "dipole":
            G *= B; Gv *= B; S *= B; O *= B                              # градиенты диполя заданы относительно поля
        # силы по ОПОРНОЙ жёсткости (по L_eff — интеграл поля); δ-зависимость — в дрейфе.
        # k1=k1u (гор. фокус ∂Bz/∂u), k1v (верт. фокус ∂By/∂v): неголоморфн. диполь Gu≠Gv ⇒ |f_x|≠|f_y|.
        self.theta0 = QP * B * self.Leff
        self.k1 = QP * G; self.k1v = QP * Gv; self.k2 = QP * S; self.k3 = QP * O
        # краевой (фриндж) верт. фокус 2-го порядка: 1/f_y = FE·θ₀²/L_eff (∝поле², offset/угол-независим).
        # Эквивалент углового члена Δv'=−(q/p)u'v∂ₓBz после ∫по частям = −(q/p)²∫Bz²ds·v. Диполи: θ₀≠0.
        self.kEdge = _edge_ff(self.dev) * self.theta0 ** 2 / self.Leff if self.Leff > 0 else 0.0
        # КРАЕВАЯ ЛИНЗА 1-го порядка (наклон полюса E1/E2 из .lte, box-модель): 1/f = tan(E)/ρ, ρ=L_eff/θ0.
        # Верт. ФОКУС (py−=k·y), гор. ДЕФОКУС (px+=k·x) — как rectangular bend. На ВНЕШНИХ гранях диполя.
        _rinv = self.theta0 / self.Leff if self.Leff > 0 else 0.0    # 1/ρ (знак по изгибу)
        self.kEdge1 = math.tan(self.edge1) * _rinv                  # вход
        self.kEdge2 = math.tan(self.edge2) * _rinv                  # выход
        self.nonzero = any(abs(v) > 0 for v in (self.theta0, self.k1, self.k1v, self.k2, self.k3,
                                                self.kEdge, self.kEdge1, self.kEdge2))
        # ПРОГИБ ОПОРНОЙ ОРБИТЫ (сагитта). Ящик прямой, а опорная внутри изгиба идёт по дуге и отходит
        # от оси ящика на L·θ0/8. Мультиполи и края обязаны раскладываться ВДОЛЬ ОРБИТЫ (elegant делает
        # это конструктивно: CSBEND — СЕКТОР, рама гнётся вместе с опорной, x≡0 по построению), поэтому
        # кик берёт u = x − xref, а не x. Иначе сама опорная ловит кик k1·сагитта: у DQ_L1000
        # (θ0=24.8мрад ⇒ сагитта 3.10мм, k1=−1.96) это ~6мрад ⇒ ~2мм за один элемент.
        self.xref = self.xref_in = self.xref_out = None
        if self.theta0 and any((self.k1, self.k1v, self.k2, self.k3, self.kEdge1, self.kEdge2)):
            self.xref, self.xref_in, self.xref_out = self._ref_bow()

    def _ref_bow(self):
        """Опорная орбита в раме ящика: (x в точках кика, x на входной плоскости, x на выходной).

        Ящик уложен по ХОРДЕ дуги, а его центр = СЕРЕДИНА дуги (lte_to_lattice.step ставит center в
        точку дуги на половине угла, а ось — по касательной там же) ⇒ в центре опорная имеет x=0 и px=0.
        Отсюда интегрируем в обе стороны, видя ТОЛЬКО изгиб: на самой опорной мультиполи по определению
        не действуют. К торцам орбита отходит на сагитту L·θ0/8."""
        N = self.nsl; ds = self.L / N; th = self.theta0
        out = [0.0] * N
        x = px = 0.0
        for i in range(N // 2, N):                             # вперёд от центра к выходу
            x += px * (ds / 2.0) / math.sqrt(max(1e-12, 1.0 - px * px)); out[i] = x
            px += th / N
            x += px * (ds / 2.0) / math.sqrt(max(1e-12, 1.0 - px * px))
        x_out = x
        x = px = 0.0
        for i in range(N // 2 - 1, -1, -1):                    # назад от центра ко входу
            x -= px * (ds / 2.0) / math.sqrt(max(1e-12, 1.0 - px * px)); out[i] = x
            px -= th / N
            x -= px * (ds / 2.0) / math.sqrt(max(1e-12, 1.0 - px * px))
        return out, x, x_out

    def body(self, st):
        """Тело: N срезов drift-kick-drift по ЭФФЕКТИВНОЙ длине (ds_g=ds_f=L_eff/N) ⇒ полный изгиб/
        фокус = k·L_eff (совпадает с .lte). Чистый диполь = только кик изгиба (без фокусировки);
        квадруполь/секступоль/октуполь — линейн./нелин. кики.

        У изгибающих элементов мультиполи и края берут u = x − xref (отклонение от ИЗОГНУТОЙ опорной),
        а не x от оси ящика — см. _ref_bow. Для неизгибающих xref=None ⇒ u≡x, поведение бит-в-бит."""
        if not self.nonzero:                                   # выключенный — дрейф (нулевой элемент)
            return _drift(st, self.L)
        th, k1, k1v, k2, k3, kE = self.theta0, self.k1, self.k1v, self.k2, self.k3, self.kEdge
        N = self.nsl
        dsg = self.L / N; dsf = self.Leff / N
        xr = self.xref
        if self.kEdge1:                                        # краевая линза на ВХОДЕ (гор.дефокус/верт.фокус)
            x, px, y, py, t, d = st
            u = x - self.xref_in if xr is not None else x
            st = (x, px + self.kEdge1 * u, y, py - self.kEdge1 * y, t, d)
        for i in range(N):
            st = _drift(st, dsg / 2.0)
            x, px, y, py, t, d = st
            u = x - xr[i] if xr is not None else x             # отклонение от опорной дуги (сектор), не от оси
            # НОРМАЛЬНЫЕ мультиполи (Re Cₙ): изгиб θ0, квад −k1·u, секступоль, октуполь + краевая линза
            px += th / N - k1 * u * dsf - 0.5 * k2 * (u * u - y * y) * dsf - (k3 / 6.0) * (u ** 3 - 3 * u * y * y) * dsf
            py += k1v * y * dsf + k2 * u * y * dsf + (k3 / 6.0) * (3 * u * u * y - y ** 3) * dsf - (kE / N) * y
            st = _drift((x, px, y, py, t, d), dsg / 2.0)
        if self.kEdge2:                                        # краевая линза на ВЫХОДЕ
            x, px, y, py, t, d = st
            u = x - self.xref_out if xr is not None else x
            st = (x, px + self.kEdge2 * u, y, py - self.kEdge2 * y, t, d)
        return st


def build_line(lattice_id, eff, phys, coup, trans, lo=None, hi=None, source="physical", settings=None):
    """Собрать решётку. settings — Settings прогона (None ⇒ дефолты ядра)."""
    st = settings or Settings()
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
        mags.append(Magnet(group, eff, phys, coup, trans, source, settings=st)); i += 1
    out = mags[lo:hi] if lo is not None else mags
    check_junctions(out, settings=st)               # стыки: снап совпадающих + контроль клиновых по порогу опыта
    return out


# ---- ПРОВЕРКА СТЫКОВ ПЛОСКОСТЕЙ (при сборке решётки) -------------------------------------------
# Выходная плоскость A и входная плоскость B могут: СОВПАДАТЬ (стык «в упор» — переходим без пролёта),
# быть ПАРАЛЛЕЛЬНЫМИ с зазором (обычный дрейф) или образовывать КЛИН (непараллельные). У клина плоскости
# пересекаются по линии; если эта линия ближе к опорной, чем порог потери, часть пучка попадает в
# вырожденную зону (для одних частиц вход B уже позади выхода A) — это ошибка геометрии.
class Settings:
    """Параметры ЯДРА для одного прогона: ни одного зашитого числа в расчётном коде.

    Источник — столбцы experiments.csv (у каждого опыта свои), любой можно переопределить
    аргументом. Пустая ячейка ⇒ значение из DEFAULTS. Порядок приоритета:
        явный аргумент  >  столбец experiments.csv  >  DEFAULTS

    loss_radius_m      — порог потери пучка. Им же меряется опасность клина между плоскостями.
    plane_same_tol_m   — плоскости ближе этого считаем ОДНОЙ (переход без пролёта).
    plane_parallel_tol — |axis_A × axis_B| ниже ⇒ нормали параллельны.
    n_slices           — срезов тела, ФОЛБЭК: приоритет у lattice_configs.n_slices данного элемента.
    """
    DEFAULTS = dict(loss_radius_m=0.012, plane_same_tol_m=1e-4,
                    plane_parallel_tol=1e-9, n_slices=40)
    # столбец experiments.csv → (поле, множитель к СИ)
    COLUMNS = {"loss_radius_mm": ("loss_radius_m", 1e-3),
               "plane_same_tol_mm": ("plane_same_tol_m", 1e-3),
               "plane_parallel_tol": ("plane_parallel_tol", 1.0),
               "n_slices": ("n_slices", 1)}

    def __init__(self, **kw):
        bad = set(kw) - set(self.DEFAULTS)
        if bad:
            raise ValueError(f"неизвестные параметры ядра: {sorted(bad)}; "
                             f"допустимы {sorted(self.DEFAULTS)}")
        for k, v in self.DEFAULTS.items():
            setattr(self, k, v if kw.get(k) is None else type(v)(kw[k]))

    @classmethod
    def from_experiment(cls, exp_id, **override):
        """Прочитать параметры опыта из experiments.csv; **override перебивает столбцы."""
        vals = {}
        p = os.path.join(TB, "experiments.csv")
        if exp_id is not None and os.path.exists(p):
            for r in csv.DictReader(open(p, encoding="utf-8"), delimiter=";"):
                if r.get("exp_id") == str(exp_id):
                    for col, (field, mul) in cls.COLUMNS.items():
                        raw = (r.get(col) or "").strip()
                        if raw:
                            vals[field] = float(raw) * mul
                    break
        vals.update({k: v for k, v in override.items() if v is not None})
        return cls(**vals)

    def __repr__(self):
        return ("Settings(" + ", ".join(f"{k}={getattr(self, k)!r}" for k in self.DEFAULTS) + ")")


def _plane_relation(A, B, st):
    """(kind, apex, gap): kind ∈ {'same','parallel','wedge'}; apex — расстояние линии пересечения
    плоскостей до опорной точки стыка (для 'wedge'); gap — отступ A.exit от плоскости B по её нормали."""
    nA, nB, pA, pB = A.axis, B.axis, A.exit, B.entry
    gap = abs(float(nB @ (pA - pB)))
    cr = np.cross(nA, nB); s = float(np.linalg.norm(cr))
    if s < st.plane_parallel_tol:                              # нормали параллельны
        return ("same" if gap < st.plane_same_tol_m else "parallel"), float("inf"), gap
    d = cr / s                                                 # направление линии пересечения
    X0 = np.linalg.lstsq(np.vstack([nA, nB]), np.array([float(nA @ pA), float(nB @ pB)]), rcond=None)[0]
    w = pA - X0
    return "wedge", float(np.linalg.norm(w - (w @ d) * d)), gap


def check_junctions(mags, settings=None, verbose=False):
    """Проверка стыков ПРИ СБОРКЕ. Совпадающие плоскости помечаются как одна (пролёт не нужен).
    Клин с линией пересечения ближе порога потери к опорной — ОШИБКА. Возврат: список проблем.
    settings — Settings прогона (None ⇒ дефолты)."""
    st = settings or Settings()
    lr = st.loss_radius_m
    problems = []; n = {"same": 0, "parallel": 0, "wedge": 0}
    for j in range(1, len(mags)):
        A, B = mags[j - 1], mags[j]
        kind, apex, gap = _plane_relation(A, B, st)
        n[kind] += 1
        B.same_plane_as_prev = (kind == "same")                # снап: переход без пролёта
        if kind == "wedge" and apex < lr:
            problems.append((j, A.dev, B.dev, apex, gap))
    if verbose:
        print(f"  стыки плоскостей: {n['same']} совмещённых (без пролёта), "
              f"{n['parallel']} параллельных с зазором, {n['wedge']} клиновых "
              f"(порог потери {lr * 1e3:.1f} мм)")
    for j, da, db, apex, gap in problems:
        print(f"  [ОШИБКА СТЫКА] #{j} {da} → {db}: линия пересечения плоскостей в {apex * 1e3:.2f} мм "
              f"от опорной — ближе порога потери {lr * 1e3:.1f} мм (зазор {gap * 1e3:.3f} мм). "
              f"Часть пучка попадёт в вырожденную зону.")
    return problems


def modified_drift(st, A, B):
    """Геометрический пролёт: плоскость выхода A → плоскость входа B (произв. ориентации). δ-корректно
    (p_полный сохраняется). Если при сборке плоскости признаны ОДНОЙ (B.same_plane_as_prev) — переход
    без пролёта. st в системе A → st в системе B, + путь."""
    x, px, y, py, t, d = st
    P = A.exit + x * A.uhat + y * A.vhat                       # глоб. позиция на выходной плоскости A
    ps = math.sqrt((1.0 + d) ** 2 - px * px - py * py)
    dir_ = (px * A.uhat + py * A.vhat + ps * A.axis) / (1.0 + d)   # единичное направление импульса (глоб.)
    if getattr(B, "same_plane_as_prev", False):
        Q, tt = P, 0.0                                         # стык «в упор» — пролёта нет
    else:
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
    lat = int(sys.argv[1]) if len(sys.argv) > 1 else 1
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
