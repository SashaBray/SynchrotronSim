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


def _rot(v, n, a):
    """Поворот вектора v вокруг ЕДИНИЧНОЙ оси n на угол a (формула Родрига)."""
    c, s = math.cos(a), math.sin(a)
    return v * c + np.cross(n, v) * s + n * float(n @ v) * (1.0 - c)


def frame(yaw, pitch=0.0, roll=0.0):
    """Ортонормальная тройка (ось, uhat, vhat) элемента из углов установки (для рассогласований).
    Порядок: yaw → pitch → roll, каждый вокруг УЖЕ повёрнутой оси (Tait–Bryan):
      • yaw   — вокруг вертикали vhat (горизонтальный увод); при pitch=roll=0 совпадает с прежним frame(yaw);
      • pitch — вокруг поперечной uhat (вертикальный наклон оси, «клевок»);
      • roll  — вокруг самой оси (крен: uhat/vhat вокруг луча).
    Опорная тройка: ось вдоль X, uhat поперёк по Y, vhat вверх по Z."""
    axis = np.array([1.0, 0.0, 0.0]); uhat = np.array([0.0, 1.0, 0.0]); vhat = np.array([0.0, 0.0, 1.0])
    if yaw:
        axis, uhat = _rot(axis, vhat, yaw), _rot(uhat, vhat, yaw)
    if pitch:
        axis, vhat = _rot(axis, uhat, pitch), _rot(vhat, uhat, pitch)
    if roll:
        uhat, vhat = _rot(uhat, axis, roll), _rot(vhat, axis, roll)
    return axis, uhat, vhat


SECTOR = os.environ.get("SS_SECTOR", "") == "1"    # ОПЫТ: диполь = точный сектор внутри ящика
_SECTOR_DISP = -1.0                                # знак дисперсии секторной карты (под «+x наружу»); сверить с elegant
_PARA_TOL = 1e-4                                   # |axis_A × axis_B| ниже ⇒ оси параллельны (стык дуги ~1e-6,
#                                                   дрейф =0); настоящий клин диполь→элемент ~мрад >> порога


def _foc(K, L):
    """Линейная карта фокусирующего элемента силой K на длине L: (C, S, S') для [x1;px1]=[[C,S],[S',C]]·[x;px].
    K>0 — фокус (cos/sin), K<0 — дефокус (cosh/sinh), K≈0 — дрейф. Общая для гор./верт. сектора+квадруполя."""
    if K > 1e-12:
        r = math.sqrt(K); p = r * L
        return math.cos(p), math.sin(p) / r, -r * math.sin(p)
    if K < -1e-12:
        r = math.sqrt(-K); p = r * L
        return math.cosh(p), math.sinh(p) / r, r * math.sinh(p)
    return 1.0, L, 0.0                                          # дрейф


# ---- СИМПЛЕКТИЧЕСКИЙ СЕКТОРНЫЙ ИНТЕГРАТОР (порт elegant CSBEND integrate_csbend_ordn) ------------
# Точный гамильтониан изгиба: комбинированная функция (диполь+квадруполь+…) НЕЛИНЕЙНА — кривизна×градиент
# рождает члены высших порядков. Линейная секторная карта (_body_sector) их теряет (DQ: ~8мкм остаток).
# Здесь — дословный перенос elegant: точный секторный дрейф (движение по кругу в кривой раме) + кик поля с
# фактором (1+X/ρ), расщепление Йошиды (порядок 2/4/6, N срезов). Координаты канонические (наш px ≡ elegant
# QX = x'·(1+δ)/√(1+x'²+y'²)), опорная орбита = 0 в криволинейной рамке (fr_in/fr_out), конвертация не нужна.
_YOSHIDA_BETA = 1.25992104989487316477            # 2^(1/3) — константа расщепления 4-го порядка (как в elegant)
_SECTOR_LINEAR = os.environ.get("SS_SECTOR_LINEAR", "") == "1"   # откат на ЛИНЕЙНУЮ секторную карту (для сверки)
# Срезов симплектического тела изгиба = N_SLICES из .lte (elegant единообразно ставит 5 всем элементам).
# Бит-в-бит зависит от N: наш N обязан РАВНЯТЬСЯ elegant. БД держит n_slices=40 (дефолт сборки, для box-
# модели квадов/секступолей), поэтому у сектора СВОЙ дефолт 5; env SS_SECTOR_NSL перекрывает для опытов.
_SECTOR_NSLICES = 5


def _yoshida(order):
    """Коэффициенты (дрейф, кик) симплектического расщепления Йошиды — дословно из elegant integrate_csbend_ordn.
    order 2 — leapfrog; 4 — Форест–Рут (BETA=2^(1/3)); 6 — из AOP-TN-2020-064. Последний кик каждого шага = 0."""
    if order == 2:
        return [0.5, 0.5], [1.0, 0.0]
    if order == 4:
        b = _YOSHIDA_BETA
        d = [0.5 / (2 - b), (1 - b) / (2 - b) / 2, (1 - b) / (2 - b) / 2, 0.5 / (2 - b)]
        k = [1.0 / (2 - b), -b / (2 - b), 1.0 / (2 - b), 0.0]
        return d, k
    if order == 6:
        d = [0.39225680523878, 0.5100434119184585, -0.47105338540975655, 0.0687531682525181,
             0.0687531682525181, -0.47105338540975655, 0.5100434119184585, 0.39225680523878]
        k = [0.784513610477560, 0.235573213359357, -1.17767998417887, 1.3151863206839063,
             -1.17767998417887, 0.235573213359357, 0.784513610477560, 0.0]
        return d, k
    raise ValueError(f"порядок интегрирования {order} — допустимы 2/4/6")


def _csbend_fields(x, y, b1, b2, b3, b4, h):
    """Нормированное поле (Fx, Fy) комбинированной функции до 4-го порядка — порт elegant computeCSBENDFields
    для НОРМАЛЬНЫХ мультиполей (скью = 0). b_n = k_n·ρ (b1=градиент/квад, b2=секступоль, b3=октуполь), h=1/ρ.
    Кривизна h подмешивает члены высших порядков (напр. Fy ⊃ −(h·b1)/2·y²; Fy ⊃ −(h³·b1)/2·x²y²) — вот
    источник нелинейности изгиба, который линейная карта не видит. Fy чётно по y, Fx нечётно (нормальные)."""
    h2 = h * h; h3 = h2 * h
    x2 = x * x; y2 = y * y
    Fy = (1.0
          + b1 * x + 0.5 * b2 * x2 + (b3 / 6.0) * x2 * x + (b4 / 24.0) * x2 * x2
          + (-(h * b1) / 2.0 - b2 / 2.0) * y2
          + ((h2 * b1) / 2.0 - (h * b2) / 2.0 - b3 / 2.0) * x * y2
          + (-(h3 * b1) / 2.0 + (h2 * b2) / 2.0 - (h * b3) / 4.0 - b4 / 4.0) * x2 * y2
          + ((h3 * b1) / 24.0 - (h2 * b2) / 24.0 + (h * b3) / 12.0 + b4 / 24.0) * y2 * y2)
    Fx = (b1 * y + b2 * x * y + (b3 / 2.0) * x2 * y + (b4 / 6.0) * x2 * x * y
          + ((h2 * b1) / 6.0 - (h * b2) / 6.0 - b3 / 6.0) * y2 * y
          + (-(h3 * b1) / 3.0 + (h2 * b2) / 3.0 - (h * b3) / 6.0 - b4 / 6.0) * x * y2 * y)
    return Fx, Fy


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
        self.pitch = float(r0.get("pitch") or 0.0)         # установочные наклоны (рассогласование); дефолт 0
        self.roll = float(r0.get("roll") or 0.0)
        self.axis, self.uhat, self.vhat = frame(self.yaw, self.pitch, self.roll)
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
        # ОПЫТ: СЛАБАЯ (СЕКТОРНАЯ) ФОКУСИРОВКА 1/ρ². В секторе x''+(1/ρ²)x=0 ⇒ Δpx=−(1/ρ²)·u·ds (гор.
        # фокус, верт. нет). У прямого ящика с постоянным киком изгиба R21=0 — этого члена нет. Вопрос:
        # не даёт ли его сама укладка на дугу (поворот рамки на стыках)? Тест сравнит обе геометрии.
        self.wkfoc = _rinv * _rinv if os.environ.get("SS_WEAK_FOCUS", "") == "1" else 0.0
        self.nonzero = any(abs(v) > 0 for v in (self.theta0, self.k1, self.k1v, self.k2, self.k3,
                                                self.kEdge, self.kEdge1, self.kEdge2, self.wkfoc))
        # ПРОГИБ ОПОРНОЙ ОРБИТЫ (сагитта). Ящик прямой, а опорная внутри изгиба идёт по дуге и отходит
        # от оси ящика на L·θ0/8. Мультиполи и края обязаны раскладываться ВДОЛЬ ОРБИТЫ (elegant делает
        # это конструктивно: CSBEND — СЕКТОР, рама гнётся вместе с опорной, x≡0 по построению), поэтому
        # кик берёт u = x − xref, а не x. Иначе сама опорная ловит кик k1·сагитта: у DQ_L1000
        # (θ0=24.8мрад ⇒ сагитта 3.10мм, k1=−1.96) это ~6мрад ⇒ ~2мм за один элемент.
        self.xref = self.xref_in = self.xref_out = None
        if self.theta0 and any((self.k1, self.k1v, self.k2, self.k3, self.kEdge1, self.kEdge2)):
            self.xref, self.xref_in, self.xref_out = self._ref_bow()
        # РАМКИ СТЫКОВКИ: fr_in — вход (точка, ось, uhat, vhat), fr_out — выход. Для НЕизгибающих ящиков —
        # плоскости ⊥ оси (как раньше). Для ДИПОЛЯ в секторном режиме — КРИВОЛИНЕЙНЫЕ рамки: грань ⊥
        # касательной опорной дуги, опорная орбита проходит через точку рамки (x=0). Тогда секторная карта
        # тела ложится напрямую (опорная = 0 на входе и выходе), а modified_drift стыкует к наклонной грани.
        self.sector = SECTOR and self.kind in ("dipole", "dq") and abs(self.theta0) > 0
        # ПРЯМОЙ (не изгибающий) мультиполь в секторном режиме — симплектика KQUAD/KSEXT/KOCT (порт elegant
        # integrate_kick_multipole_ordn, order 4, N=5). Даёт квады/секступоли/октуполи бит-в-бит с elegant
        # (box-модель = order 2, N=40 давала ~0.13мкм/шт). Рамки — обычные ⊥оси (прямой элемент, ветка else).
        self.kick_sector = (SECTOR and not self.sector and self.nonzero
                            and self.kind in ("quad", "sext", "oct"))
        if self.sector or self.kick_sector:
            # N срезов и порядок Йошиды = как elegant (N_SLICES=5, INTEGRATION_ORDER=4); env перекрывает.
            self._nsl_sector = int(os.environ.get("SS_SECTOR_NSL", "") or _SECTOR_NSLICES)
            self._yoshida_dk = _yoshida(int(os.environ.get("SS_SECTOR_ORDER", "") or 4))
        if self.sector:
            rho = self.Leff / self.theta0
            O = self.center + rho * self.uhat                  # центр кривизны дуги
            # изгиб идёт вокруг ЛОКАЛЬНОЙ вертикали vhat (при pitch/roll плоскость изгиба наклонена/повёрнута);
            # для yaw-only vhat=[0,0,1] ⇒ _rot(·,vhat,·) ≡ прежний _rotz, поведение бит-в-бит сохраняется.
            v = self.vhat; th2 = self.theta0 / 2.0
            self.entry_pt = O + _rot(self.center - O, v, -th2)   # точка входа НА ДУГЕ
            self.exit_pt = O + _rot(self.center - O, v, th2)
            self.fr_in = (self.entry_pt, _rot(self.axis, v, -th2), _rot(self.uhat, v, -th2), v)
            self.fr_out = (self.exit_pt, _rot(self.axis, v, th2), _rot(self.uhat, v, th2), v)
        else:
            self.fr_in = (self.entry, self.axis, self.uhat, self.vhat)
            self.fr_out = (self.exit, self.axis, self.uhat, self.vhat)

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

    def _body_sector(self, st):
        """ТОЧНЫЙ СЕКТОР (линейная карта) в криволинейной рамке: опорная=0 на входе/выходе. Тело =
        край(E1)·сектор·край(E2). КОМБИНИРОВАННАЯ функция (диполь+квадруполь): гор. фокусировка
        K_x=1/ρ²+k1, верт. K_y=−k1v, через cos/cosh. Чистый диполь = частный случай k1=0 (K_x=1/ρ²).
        Дисперсия R16=h(1−C_x)/K_x, R26=h·S_x, h=1/ρ. Края — тонкие линзы tan(E)/ρ на ⊥касат. гранях.
        Даёт R11/R12/R16/R21/R22/R26/R33/R34/R43/R44 аналитически (сверено с elegant для DL и DQ)."""
        x, px, y, py, t, d = st
        h = self.theta0 / self.Leff                            # кривизна 1/ρ
        Kx = h * h + self.k1                                   # гор. фокусировка (сектор + квадруполь)
        Ky = -self.k1v                                         # верт. фокусировка (квадруполь; сектор верт. не фокус.)
        Cx, Sx, Spx = _foc(Kx, self.Leff)
        Cy, Sy, Spy = _foc(Ky, self.Leff)
        R16 = h * (1.0 - Cx) / Kx if abs(Kx) > 1e-12 else h * self.Leff * self.Leff / 2.0
        R26 = h * Sx
        disp = _SECTOR_DISP                                     # знак дисперсии под конвенцию «+x наружу»
        px += self.kEdge1 * x; py -= self.kEdge1 * y            # вход: край E1
        x1 = Cx * x + Sx * px + disp * R16 * d
        px1 = Spx * x + Cx * px + disp * R26 * d
        y1 = Cy * y + Sy * py
        py1 = Spy * y + Cy * py
        x, px, y, py = x1, px1, y1, py1
        px += self.kEdge2 * x; py -= self.kEdge2 * y            # выход: край E2
        return (x, px, y, py, t + self.Leff, d)

    def _body_sector_symplectic(self, st):
        """ТОЧНЫЙ секторный интегратор изгиба — дословный порт elegant CSBEND integrate_csbend_ordn.

        Комбинированная функция (диполь+квадруполь) НЕЛИНЕЙНА: кривизна связывается с градиентом и рождает
        члены x²/x·y²/… которые линейная карта _body_sector теряет (DQ давал ~8мкм на выходе). Здесь на каждом
        срезе — симплектический drift-kick:
          • ДРЕЙФ точный, по кругу в криволинейной рамке (asin/sin/cos, фактор (ρ+X)) — не прямой;
          • КИК поля _csbend_fields с множителем (1+X/ρ) — вот он и делает линейный градиент нелинейным.
        Опорная орбита = 0 на входе/выходе (рамки fr_in/fr_out криволинейны), координаты уже канонические
        (px ≡ elegant QX) ⇒ конвертация slope↔canonical не нужна. N срезов и порядок Йошиды = как elegant.
        Края E1/E2 — тонкие линзы tan(E)/ρ на ⊥касат. гранях (elegant edge_order=1, edge_effects=1)."""
        x, px, y, py, t, d = st
        # ρ0 в КОНВЕНЦИИ ELEGANT (положительный): elegant rho0=L/ANGLE, а наш θ0=−ANGLE ⇒ rho0=−Leff/θ0.
        # Это КРИТИЧНО: знак ρ0 задаёт знак нелинейности изгиба (фактор (1+X/ρ), члены h·b в поле) и
        # дисперсии. С нашим θ0<0 напрямую они выходили ЗЕРКАЛЬНЫМИ (DQ давал −16мкм вместо +8). При
        # rho0=−Leff/θ0>0 всё совпадает с elegant бит-в-бит (b1=k1·rho0 ≡ elegant K1·rho0). Наши (x,px,δ)
        # уже в знаках elegant (сверено по фокус./верт. R-блоку), геометрию поворота держат рамки fr_in/fr_out.
        rho0 = -self.Leff / self.theta0                        # ρ>0 (конвенция elegant), θ0<0 у изгиба
        h = 1.0 / rho0                                         # кривизна 1/ρ = −θ0/Leff
        b1 = self.k1 * rho0; b2 = self.k2 * rho0; b3 = self.k3 * rho0; b4 = 0.0   # b_n = k_n·ρ (норм. мультиполи)
        N = self._nsl_sector; s = self.Leff / N
        dfrac, kfrac = self._yoshida_dk
        px += self.kEdge1 * x; py -= self.kEdge1 * y           # вход: краевой фокус E1
        X, QX, Y, QY, DP = x, px, y, py, d                     # elegant-обозначения (наши канонические координаты)
        dist = 0.0
        for _ in range(N):
            for j in range(len(dfrac)):
                dsh = s * dfrac[j]                             # -- точный секторный ДРЕЙФ (движение по кругу) --
                f = (1.0 + DP) ** 2 - QY * QY
                if f <= 0.0:
                    raise ParticleLost("сектор: p⊥>p (потеря в дрейфе)")
                f = math.sqrt(f)
                if abs(QX / f) > 1.0:
                    raise ParticleLost("сектор: |QX/f|>1 (потеря в дрейфе)")
                sin_phi = QX / f
                phi = math.asin(sin_phi)
                sine = math.sin(dsh / rho0 + phi); cosi = math.cos(dsh / rho0 + phi)
                if cosi == 0.0:
                    raise ParticleLost("сектор: cos=0 (потеря в дрейфе)")
                cos_phi = math.cos(phi)
                QX = f * sine
                factor = (rho0 + X) * cos_phi / f * (sine / cosi - sin_phi / cos_phi)
                Y += QY * factor
                dist += factor * (1.0 + DP)
                fr = cos_phi / cosi
                X = rho0 * (fr - 1.0) + fr * X
                if kfrac[j] == 0.0:
                    break
                ds = s * kfrac[j]                              # -- КИК поля с фактором кривизны (1+X/ρ) --
                Fx, Fy = _csbend_fields(X, Y, b1, b2, b3, b4, h)
                QX += -ds * (1.0 + X / rho0) * Fy / rho0
                QY += ds * (1.0 + X / rho0) * Fx / rho0
        x, px, y, py, d = X, QX, Y, QY, DP
        px += self.kEdge2 * x; py -= self.kEdge2 * y           # выход: краевой фокус E2
        return (x, px, y, py, t + dist, d)

    def _body_kick_symplectic(self, st):
        """Прямой (не изгибающий) мультиполь — порт elegant `integrate_kick_multipole_ordn` (KQUAD/KSEXT/KOCT).

        Дрейф по СЛОПУ (наш _drift: x+=x'·ds, x'=px/√((1+δ)²−px²−py²)) + канонический НОРМАЛЬНЫЙ мультипольный
        кик, расщепление Йошиды порядка 4, N=5 срезов — ровно как elegant (INTEGRATION_ORDER=4, N_SLICES=5).
        Наши px,py уже канонические (≡ elegant qx,qy), конвертация slope↔canonical встроена в _drift/кик. Кик
        тот же нормальный мультиполь, что в box-теле (квад −k1·x/+k1v·y, секступоль, октуполь). ОТЛИЧИЕ от box:
        порядок 4 (box = 2, leapfrog) и N=5 (box = 40) ⇒ бит-в-бит с elegant (box давал ~0.13мкм/шт на квадах).
        Прямой квад фринджа не имеет (elegant EDGE1_EFFECTS=0 по умолчанию), краевых киков нет."""
        x, px, y, py, t, d = st
        k1, k1v, k2, k3 = self.k1, self.k1v, self.k2, self.k3
        N = self._nsl_sector; ds = self.Leff / N
        dfrac, kfrac = self._yoshida_dk
        for _ in range(N):
            for j in range(len(dfrac)):
                x, px, y, py, t, d = _drift((x, px, y, py, t, d), ds * dfrac[j])   # дрейф по слопу
                if kfrac[j] == 0.0:
                    break
                kf = ds * kfrac[j]                             # доля кика на этот подшаг × длина среза
                px += (-k1 * x - 0.5 * k2 * (x * x - y * y) - (k3 / 6.0) * (x ** 3 - 3.0 * x * y * y)) * kf
                py += (k1v * y + k2 * x * y + (k3 / 6.0) * (3.0 * x * x * y - y ** 3)) * kf
        return (x, px, y, py, t, d)

    def body(self, st):
        """Тело: N срезов drift-kick-drift по ЭФФЕКТИВНОЙ длине (ds_g=ds_f=L_eff/N) ⇒ полный изгиб/
        фокус = k·L_eff (совпадает с .lte). Чистый диполь = только кик изгиба (без фокусировки);
        квадруполь/секступоль/октуполь — линейн./нелин. кики.

        У изгибающих элементов мультиполи и края берут u = x − xref (отклонение от ИЗОГНУТОЙ опорной),
        а не x от оси ящика — см. _ref_bow. Для неизгибающих xref=None ⇒ u≡x, поведение бит-в-бит."""
        if self.sector:                                        # изгиб в секторном режиме
            if _SECTOR_LINEAR:                                 # SS_SECTOR_LINEAR=1 — линейная карта (сверка/чистый диполь)
                return self._body_sector(st)
            return self._body_sector_symplectic(st)            # по умолчанию — точный симплектический интегратор
        if self.kick_sector and not _SECTOR_LINEAR:            # прямой мультиполь — симплектика KQUAD/KSEXT/KOCT
            return self._body_kick_symplectic(st)              # (при SS_SECTOR_LINEAR=1 падает на box-модель ниже)
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
            # + слабая секторная фокусировка −(1/ρ²)·u (гор., опыт SS_WEAK_FOCUS)
            px += th / N - self.wkfoc * u * dsf - k1 * u * dsf - 0.5 * k2 * (u * u - y * y) * dsf - (k3 / 6.0) * (u ** 3 - 3 * u * y * y) * dsf
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


PURE_DRIFT = os.environ.get("SS_PURE_DRIFT", "") == "1"    # ОПЫТ: elegant-подобный чистый дрейф между элементами


def _pure_drift(st, A, B):
    """ЧИСТЫЙ дрейф (как в elegant): прямой отрезок вдоль оси ВЫХОДА A, выходная плоскость ПАРАЛЛЕЛЬНА
    входной (обе ⊥ axis_A). Длина — до плоскости ⊥axis_A, проходящей через начало рамки входа B.
    Затем состояние переразлагается в рамке входа B. Наклон грани B учитывается НЕ здесь, а краевым
    киком в теле B (kEdge1) — так же, как отделяет дрейф от края elegant. Рамки — fr_out(A)/fr_in(B):
    у диполя-сектора это КРИВОЛИНЕЙНЫЕ грани (⊥ касательной опорной), у прочих — ⊥ оси ящика.
    δ-корректно (p_полный сохраняется). st в системе A → st в системе B, + путь."""
    aP, aAx, aU, aV = A.fr_out
    bP, bAx, bU, bV = B.fr_in
    Ld = 0.0 if getattr(B, "same_plane_as_prev", False) else float((bP - aP) @ aAx)
    x2, px2, y2, py2, t2, d = _drift(st, Ld)                   # прямой δ-корректный дрейф вдоль axis_A
    P = aP + Ld * aAx + x2 * aU + y2 * aV                      # глоб. позиция на выходной плоскости дрейфа
    ps = math.sqrt((1.0 + d) ** 2 - px2 * px2 - py2 * py2)
    dir_ = (px2 * aU + py2 * aV + ps * aAx) / (1.0 + d)
    rel = P - bP                                              # плоскость дрейфа проходит через bP ⇒ rel·axis_B мал
    x3 = float(rel @ bU); y3 = float(rel @ bV)
    px3 = (1.0 + d) * float(dir_ @ bU); py3 = (1.0 + d) * float(dir_ @ bV)
    return (x3, px3, y3, py3, t2, d)


def modified_drift(st, A, B):
    """Геометрический пролёт: плоскость выхода A → плоскость входа B (произв. ориентации). δ-корректно
    (p_полный сохраняется). Если при сборке плоскости признаны ОДНОЙ (B.same_plane_as_prev) — переход
    без пролёта. st в системе A → st в системе B, + путь."""
    if PURE_DRIFT:
        return _pure_drift(st, A, B)
    # рамки стыковки: выход A (fr_out) → вход B (fr_in). У диполя-сектора это НАКЛОННЫЕ грани (⊥
    # касательной опорной), у прочих — ⊥ оси (как раньше). modified_drift сам ориентирует к наклонной грани.
    aP, aAx, aU, aV = A.fr_out
    bP, bAx, bU, bV = B.fr_in
    x, px, y, py, t, d = st
    P = aP + x * aU + y * aV                                   # глоб. позиция на выходной плоскости A
    ps = math.sqrt((1.0 + d) ** 2 - px * px - py * py)
    dir_ = (px * aU + py * aV + ps * aAx) / (1.0 + d)          # единичное направление импульса (глоб.)
    if getattr(B, "same_plane_as_prev", False):
        Q = P                                                 # стык «в упор» — пролёта нет
    elif float(np.linalg.norm(np.cross(aAx, bAx))) < _PARA_TOL:
        # ОСИ ПАРАЛЛЕЛЬНЫ (стык дуговых срезов = общая касательная; прямой дрейф) ⇒ ЧИСТЫЙ ДРЕЙФ (как в
        # elegant): частица летит вдоль своего направления до плоскости ⊥aAx через B. Длина знаковая
        # (μm-зазор butted-стыка любого знака корректен). Клина нет — луч-в-плоскость тут вырождается в дрейф.
        Q = P + float((bP - P) @ aAx) / float(dir_ @ aAx) * dir_
    else:
        denom = float(bAx @ dir_)                             # НАСТОЯЩИЙ клин (диполь → элемент под углом)
        if denom <= 1e-9:
            raise ParticleLost(f"пучок не долетит до входа {B.dev}: n·d={denom:.3g} ≤ 0")
        tt = float(bAx @ (bP - P)) / denom
        if tt < -1e-4:                                        # μm-перекрытие butted-клина (⊥-ящики на дуге с
            raise ParticleLost(f"плоскость входа {B.dev} позади частицы: t={tt:.3g} < 0")  # точной длиной) ок; >0.1мм = потеря
        Q = P + tt * dir_
    tt = float((Q - P) @ dir_)                                # длина пути вдоль импульса (знаковая)
    rel = Q - bP
    x2 = float(rel @ bU); y2 = float(rel @ bV)
    px2 = (1.0 + d) * float(dir_ @ bU); py2 = (1.0 + d) * float(dir_ @ bV)
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
