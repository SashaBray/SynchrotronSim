# -*- coding: utf-8 -*-
"""Регрессионные (юнит) тесты SynchrotronSim — защита от поломок при правках.

Запуск (из корня проекта):
    py -3.14 test_synchrotronsim.py            # все тесты, подробно
    py -3.14 -m unittest test_synchrotronsim   # то же через unittest
    py -3.14 test_synchrotronsim.py TestSymplecticity   # один класс

Идея: держим здесь тесты на КАЖДУЮ продакшн-программу. Прогнав файл, убеждаемся, что
изменение в одной части не сломало другую. Только ПРОДАКШН-код (не опытные скрипты, не legacy).

Покрытие:
  • matrix_tracker.py  — кик/матричное ядро (симплектичность, бит-в-бит с elegant, конвенции);
  • конвейер БД        — build_line, парсинг .lte (gated на загруженную БД);
  • traj_compare.py    — импорт;
  • C++ SynchrotronTracker.exe — наличие бинаря и smoke-прогон (gated).

Тесты ядра — на СЕКТОРНОМ режиме (бит-в-бит с elegant): SS_SECTOR=1 включается ДО импорта.
Без БД и без elegant (эталонные значения зашиты как golden из проверенных прогонов).
"""
import os, sys, math, unittest, importlib, warnings

warnings.filterwarnings("ignore", category=ResourceWarning)   # предсущ. утечки чтения в ядре — не цель тестов
os.environ["SS_SECTOR"] = "1"                    # секторный (симплектический) режим — ДО импорта matrix_tracker
try:
    sys.stdout.reconfigure(encoding="utf-8"); sys.stderr.reconfigure(encoding="utf-8")
except Exception:
    pass
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import numpy as np
import matrix_tracker as mt
from _common import exe_name


# ------------------------------------------------------------------ фикстуры/утилиты
def make_magnet(kind, Leff, angle=0.0, k1=0.0, k2=0.0, k3=0.0, edge1=0.0, edge2=0.0, nsl=5,
                center=(0.0, 0.0, 0.0), yaw=0.0, pitch=0.0, roll=0.0, lf1=None, lf2=None):
    """Синтетический Magnet БЕЗ БД: физика прямо в mat_* (mat_B=−angle/(QP·L), mat_G=k1/QP, …).
    center/yaw/pitch/roll — положение и ориентация (для тестов рассогласования установки);
    на само тело они не влияют (тело в локальных θ0/k), но задают рамки fr_in/fr_out.
    lf1/lf2 — L_f торцов (фриндж HGAP/FINT); None ⇒ столбцов нет вовсе (как в старой БД)."""
    QP = mt.QP
    matB = -angle / (QP * Leff) if angle else 0.0
    cx, cy, cz = center
    insts = [dict(x=repr(float(cx)), y=repr(float(cy)), z=repr(float(cz)),
                  yaw=repr(float(yaw)), pitch=repr(float(pitch)), roll=repr(float(roll)),
                  device_type_id="TEST", instance_id="1",
                  eff_length=repr(float(Leff)), kind=kind, n_slices=str(nsl),
                  edge1=repr(float(edge1)), edge2=repr(float(edge2)),
                  mat_B=repr(matB), mat_G=repr(k1 / QP), mat_Gv=repr(k1 / QP),
                  mat_S=repr(k2 / QP), mat_O=repr(k3 / QP))]
    if lf1 is not None: insts[0]["lf1"] = repr(float(lf1))
    if lf2 is not None: insts[0]["lf2"] = repr(float(lf2))
    return mt.Magnet(insts, {}, {}, {}, {}, settings=mt.Settings())


# эталонные элементы кольца SILA (из .lte configuration_v1)
def dq1b():   return make_magnet("dq", 1.0, angle=0.0247680586728496, k1=-1.960219550894432,
                                 edge1=-0.0116800000099, edge2=-0.0116800000099)
def dl1a5():  return make_magnet("dipole", 0.4178302891796566, angle=0.0094691022946895,
                                 edge1=-0.01268495408, edge2=0.00307504592)
def qf1a():   return make_magnet("quad", 0.295, k1=2.462177415336321)
def qd2a():   return make_magnet("quad", 0.162, k1=-2.970346620956298)
def sd1a():   return make_magnet("sext", 0.166, k2=-136.1950802309929)


J4 = np.array([[0, 1, 0, 0], [-1, 0, 0, 0], [0, 0, 0, 1], [0, 0, -1, 0]], float)


def jacobian4(fn, base, h=1e-7):
    """4×4 матрица ∂(x,px,y,py)_вых/∂(x,px,y,py)_вх конечными разностями вокруг base (6-вектор)."""
    M = np.zeros((4, 4))
    for j in range(4):
        sp = list(base); sp[j] += h
        sm = list(base); sm[j] -= h
        op = fn(tuple(sp)); om = fn(tuple(sm))
        for i in range(4):
            M[i, j] = (op[i] - om[i]) / (2 * h)
    return M


# ================================================================== 1. импорты
class TestImportsSmoke(unittest.TestCase):
    """Все продакшн-модули импортируются (ловит синтаксис и битые импорты)."""
    PROD = ["matrix_tracker", "track3", "traj_compare", "beam_injection", "bip_matrix",
            "ring_analysis", "ring_matrix", "ring_response", "elegant_diag",
            "lte_ring", "load_config", "lte_to_lattice", "assemble_lattice_configs",
            "build_matrix_db", "warmstart_ring", "lte_warmstart", "db_snapshot",
            "matrix_optimizer", "design_orbit"]

    def test_all_import(self):
        for mod in self.PROD:
            with self.subTest(module=mod):
                importlib.import_module(mod)


# ================================================================== 2. дрейф и координаты
class TestDrift(unittest.TestCase):
    """Точный дрейф matrix_tracker._drift."""

    def test_momentum_conserved(self):
        # p⊥² + p∥² = (1+δ)²: полный импульс сохраняется
        st = (1e-3, 0.5e-3, -0.4e-3, 0.3e-3, 0.0, 2e-3)
        x, px, y, py, t, d = mt._drift(st, 0.7)
        ps = math.sqrt((1 + d) ** 2 - px * px - py * py)
        self.assertAlmostEqual(px * px + py * py + ps * ps, (1 + d) ** 2, places=12)

    def test_roundtrip_identity(self):
        # дрейф(+L)·дрейф(−L) = тождество
        st = (2e-3, 0.6e-3, 1e-3, -0.2e-3, 0.0, 1e-3)
        back = mt._drift(mt._drift(st, 1.3), -1.3)
        for a, b in zip(st, back):
            self.assertAlmostEqual(a, b, places=12)


# ================================================================== 3. симплектичность
class TestSymplecticity(unittest.TestCase):
    """Мᵀ·J·M = J для тела каждого типа — необходимый признак корректного симплектического
    интегратора. Проверяет физику БЕЗ elegant. Линеаризуем вокруг НЕнулевой точки (иначе у
    секступоля/октуполя якобиан в 0 вырождается в дрейф)."""
    BASE = (0.8e-3, 0.2e-3, 0.5e-3, -0.15e-3, 0.0, 0.0)
    TOL = 1e-6

    def _check(self, m, base=None):
        M = jacobian4(m.body, base or self.BASE)
        err = float(np.max(np.abs(M.T @ J4 @ M - J4)))
        self.assertLess(err, self.TOL, f"нарушение симплектичности: max|MᵀJM−J|={err:.2e}")

    def test_dq_sector(self):    self._check(dq1b())
    def test_dipole_sector(self): self._check(dl1a5())
    def test_quad_focus(self):   self._check(qf1a())
    def test_quad_defocus(self): self._check(qd2a())
    def test_sextupole(self):    self._check(sd1a())
    def test_octupole(self):     self._check(make_magnet("oct", 0.09, k3=-13200.0))


# ================================================================== 4. бит-в-бит с elegant
class TestElegantGolden(unittest.TestCase):
    """Эталонные элементы воспроизводят выход elegant. Golden — из проверенных прогонов
    (см. сессию 2026-07-20, порт CSBEND/KQUAD). Ловит поломку физики интеграторов."""

    def test_dq1b_matches_elegant(self):
        # вход DQ1B (elegant): x=+0.015621мм, x'=+0.368903мрад → выход elegant: +0.543111мм, +0.853215мрад
        out = dq1b().body((0.015621e-3, 0.368903e-3, 0.0, 0.0, 0.0, 0.0))
        self.assertAlmostEqual(out[0] * 1e3, 0.543111, places=2,   # мм, допуск ~0.5мкм
                               msg="DQ1B: x не совпал с elegant")
        self.assertAlmostEqual(out[1] * 1e3, 0.853215, places=2,   # мрад, допуск ~0.5мкрад
                               msg="DQ1B: px не совпал с elegant")

    def test_reference_orbit_offset(self):
        # опорная (0 на входе) через DQ уводится НЕЛИНЕЙНОСТЬЮ на несколько мкм В ПЛЮС. Линейная карта
        # держала бы ~0, а старый (ошибочный) ρ0<0 давал МИНУС — guard ловит обе регрессии разом.
        off = dq1b().body((0.0, 0.0, 0.0, 0.0, 0.0, 0.0))[0] * 1e6   # мкм
        self.assertGreater(off, 5.0, f"опорная DQ не уводится (off={off:+.2f}мкм): линейн.карта или знак ρ0")


# ================================================================== 5. конвенции (guards на ГОТЧИ)
class TestConventions(unittest.TestCase):
    """Защита от возврата ошибок, которые мы уже ловили (см. память dq-symplectic-*)."""

    def test_dispersion_sign_positive(self):
        # ГОТЧА ρ0: дисперсия R16=∂x/∂δ у DQ должна быть ПОЛОЖИТЕЛЬНОЙ (+0.0145). При ρ0<0 была −.
        m = dq1b(); h = 1e-7
        op = m.body((0, 0, 0, 0, 0, h)); om = m.body((0, 0, 0, 0, 0, -h))
        R16 = (op[0] - om[0]) / (2 * h)
        self.assertGreater(R16, 0.01, f"дисперсия DQ отрицательна ({R16:+.4f}) — регресс знака ρ0")

    def test_default_nslices_is_5(self):
        # N срезов сектора = N_SLICES elegant (=5), не self.nsl из БД (=40)
        self.assertEqual(dq1b()._nsl_sector, 5)
        self.assertEqual(mt._SECTOR_NSLICES, 5)

    def test_dq_is_nonlinear(self):
        # симплектика ≠ линейная карта у DQ (>1мкм) — иначе кто-то откатил на линейную (регресс −8мкм)
        m = dq1b(); st = (0.5e-3, 0.0, 0.0, 0.0, 0.0, 0.0)
        lin = m._body_sector(st); sym = m._body_sector_symplectic(st)
        self.assertGreater(abs(sym[0] - lin[0]) * 1e6, 1.0,
                           "DQ симплектика ≈ линейная — потеряна нелинейность")


# ================================================================== 6. поле / Йошида
class TestRadiation(unittest.TestCase):
    """Распределённое синхротронное излучение (SS_SR): секторный диполь теряет энергию (δ↓);
    выключено (дефолт) — rad_coef=0, тело бит-в-бит без излучения (симплектика сохраняется)."""

    def _dipole(self, sr):
        old = os.environ.get("SS_SR")                       # SR читается в __init__ по env — строим при нужном флаге
        os.environ.setdefault("SS_DESIGN_MEV", "6000")      # проектная энергия СИЛА
        if sr:
            os.environ["SS_SR"] = "1"
        else:
            os.environ.pop("SS_SR", None)
        try:
            return make_magnet("dipole", 0.4178302891796566, angle=0.0094691022946895)
        finally:
            os.environ.pop("SS_SR", None) if old is None else os.environ.__setitem__("SS_SR", old)

    def test_sr_off_by_default(self):
        self.assertEqual(self._dipole(False).rad_coef, 0.0)   # без SS_SR радиации нет

    def test_sr_loses_energy(self):
        st = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)                 # осевая (проектная) частица
        d_off = self._dipole(False).body(st)[5]
        d_on = self._dipole(True).body(st)[5]
        self.assertEqual(d_off, 0.0)                        # без излучения энергия сохраняется
        self.assertLess(d_on, 0.0)                          # с излучением δ падает (потеря энергии)
        self.assertTrue(1e-9 < -d_on < 1e-5)                # порядок (2/3)r_e·γ³/ρ²·L ≈ 6.5e-7 на этот диполь

    def test_sr_preserves_angle(self):
        # радиация уносит импульс вдоль скорости: физ. угол QX/(1+δ) почти неизменен (энергия падает)
        st = (1e-3, 2e-4, 0.0, 0.0, 0.0, 0.0)
        so = self._dipole(False).body(st); sr = self._dipole(True).body(st)
        ang_off = so[1] / (1 + so[5]); ang_on = sr[1] / (1 + sr[5])
        self.assertAlmostEqual(ang_off, ang_on, delta=1e-8)   # угол сохранён с точностью до высших порядков
        self.assertLess(sr[5], 0.0)                            # энергия потеряна на излучение

    def test_radiation_table_rule(self):
        # ПРАВИЛО: строка опыта в radiation_params.csv → dict(энергия/ВЧ); нет строки → None (СИ выкл)
        r = mt.load_radiation(1)                               # реш.1 (СИЛА) — строка есть
        self.assertIsNotNone(r)
        self.assertIn("design_mev", r); self.assertGreater(r["design_mev"], 0)
        self.assertIsNone(mt.load_radiation(99999))           # несуществующий опыт → None
        self.assertIsNone(mt.load_radiation(None))            # без opt_id → None

    def test_rf_spec_types(self):
        # тип резонатора из radiation_params: none → None, ideal → "ideal", sine → dict параметров
        self.assertIsNone(mt.rf_spec(None))
        self.assertIsNone(mt.rf_spec({"rf": "none"}))
        self.assertEqual(mt.rf_spec({"rf": "ideal"}), "ideal")
        self.assertEqual(mt.ideal_delta("ideal"), 0.0)                         # пусто — проектная энергия
        si = mt.rf_spec({"rf": "ideal", "rf_ideal_energy_mev": 6006.0, "design_mev": 6000.0})
        self.assertEqual(si["type"], "ideal")
        self.assertAlmostEqual(mt.ideal_delta(si), 1e-3, places=15)            # δ = E/E_design − 1
        self.assertIsNone(mt.ideal_delta(None)); self.assertIsNone(mt.ideal_delta({"type": "sine"}))
        s = mt.rf_spec({"rf": "sine", "rf_voltage_mv": 8.0, "rf_harmonic": 1296, "rf_phase_deg": None,
                        "rf_freq_hz": None, "design_mev": 6000.0})
        self.assertEqual(s["type"], "sine"); self.assertEqual(s["harmonic"], 1296.0)
        self.assertIsNone(s["phase_deg"]); self.assertIsNone(s["freq_hz"])     # auto-фаза, f = h·c/L
        with self.assertRaises(SystemExit):
            mt.rf_spec({"rf": "sine"})                                         # sine без V/h — ошибка
        with self.assertRaises(SystemExit):
            mt.rf_spec({"rf": "klystron"})                                     # неизвестный тип
        r26 = mt.load_radiation(26)                                            # позиция БД: СИЛА-26, sine
        self.assertIsNotNone(r26); self.assertEqual(r26["rf"], "sine")
        self.assertEqual(mt.rf_spec(r26)["volt_mv"], 8.0)


class TestFieldYoshida(unittest.TestCase):
    """Строительные блоки симплектического интегратора."""

    def test_yoshida_fractions_sum(self):
        for order in (2, 4, 6):
            d, k = mt._yoshida(order)
            self.assertAlmostEqual(sum(d), 1.0, places=12, msg=f"Σdrift≠1 (order {order})")
            self.assertAlmostEqual(sum(k), 1.0, places=12, msg=f"Σkick≠1 (order {order})")

    def test_field_reference(self):
        # на опорной оси нормированное поле = чистый диполь: Fy=1, Fx=0
        Fx, Fy = mt._csbend_fields(0.0, 0.0, b1=0.0, b2=0.0, b3=0.0, b4=0.0, h=0.02)
        self.assertAlmostEqual(Fx, 0.0, places=12)
        self.assertAlmostEqual(Fy, 1.0, places=12)


# ================================================================== 7. рассогласование установки
class TestMisalignment(unittest.TestCase):
    """3D-ориентация frame(yaw,pitch,roll) + корректный стык луч-в-плоскости при рассогласовании.
    Инварианты дрейфа (прямая в глобале) держатся при ЛЮБОЙ ориентации — проверка БЕЗ elegant."""

    def _global(self, st, fr):
        """Глобальные (положение, единичное направление импульса) из локального состояния и рамки."""
        x, px, y, py, t, d = st; pt, ax, u, v = fr
        ps = math.sqrt((1 + d) ** 2 - px * px - py * py)
        return pt + x * u + y * v, (px * u + py * v + ps * ax) / (1 + d)

    # --- frame() ---
    def test_frame_yaw_backward_compat(self):
        # frame(yaw) без pitch/roll = прежняя формула (бит-в-бит)
        for yaw in (0.0, 0.1, -0.37):
            ax, u, v = mt.frame(yaw)
            np.testing.assert_allclose(ax, [math.cos(yaw), math.sin(yaw), 0.0], atol=1e-14)
            np.testing.assert_allclose(u, [-math.sin(yaw), math.cos(yaw), 0.0], atol=1e-14)
            np.testing.assert_allclose(v, [0.0, 0.0, 1.0], atol=1e-14)

    def test_frame_orthonormal_righthanded(self):
        # тройка ортонормирована и правая при любых yaw/pitch/roll
        for yaw, pitch, roll in [(0.1, 0.05, -0.03), (0.3, -0.2, 0.15), (-0.5, 0.4, 0.6)]:
            ax, u, v = mt.frame(yaw, pitch, roll)
            for a in (ax, u, v):
                self.assertAlmostEqual(np.linalg.norm(a), 1.0, places=12)
            self.assertAlmostEqual(float(ax @ u), 0.0, places=12)
            self.assertAlmostEqual(float(ax @ v), 0.0, places=12)
            self.assertAlmostEqual(float(u @ v), 0.0, places=12)
            np.testing.assert_allclose(np.cross(ax, u), v, atol=1e-12)   # правая: ось×uhat=vhat

    def test_pitch_roll_tilt_out_of_plane(self):
        ax, _, _ = mt.frame(0.0, 0.1, 0.0)
        self.assertAlmostEqual(abs(ax[2]), math.sin(0.1), places=12)     # pitch наклоняет ось по вертикали
        _, u, _ = mt.frame(0.0, 0.0, 0.1)
        self.assertAlmostEqual(abs(u[2]), math.sin(0.1), places=12)      # roll выводит uhat из горизонтали

    # --- луч-в-плоскость при рассогласовании ---
    def _check_junction(self, A, B, st):
        P_A, dA = self._global(st, A.fr_out)
        stB = mt.modified_drift(st, A, B)
        P_B, dB = self._global(stB, B.fr_in)
        self.assertLess(float(np.linalg.norm(dA - dB)), 1e-12, "дрейф изменил направление импульса")
        self.assertLess(float(np.linalg.norm(np.cross(P_B - P_A, dA))), 1e-9, "выход НЕ на луче из входа")
        self.assertLess(abs(float((P_B - B.fr_in[0]) @ B.fr_in[1])), 1e-9, "выход НЕ на плоскости входа B")
        return stB

    def test_raytoplane_connects_misaligned_dipole(self):
        # стык квадруполь → диполь при рассогласовании (yaw/pitch/roll + сдвиг): инварианты держатся
        A = make_magnet("quad", 0.2, k1=1.0, center=(-0.1, 0, 0))
        st = (0.4e-3, 0.7e-3, 0.2e-3, -0.1e-3, 0.0, 1e-3)
        cases = [(0, 0, 0, 0, 0, 0),                       # на месте
                 (0.005, 0, 0, 0.003, 0.002, 0),           # yaw + гор.сдвиг
                 (0, 0.004, 0, 0, 0, 0.001),               # pitch + верт.сдвиг
                 (0, 0, 0.02, 0.001, 0, 0),                # roll (крен)
                 (0.006, -0.003, 0.015, 0.002, -0.001, 0.0015)]  # всё сразу
        for yaw, pitch, roll, dx, dy, dz in cases:
            with self.subTest(yaw=yaw, pitch=pitch, roll=roll):
                B = make_magnet("dipole", 0.5, angle=0.02,
                                center=(2 + dx, dy, dz), yaw=yaw, pitch=pitch, roll=roll)
                self._check_junction(A, B, st)

    def test_transverse_shift_enters_off_center(self):
        # сдвиг диполя на Δ⊥ ⇒ опорный луч входит в B со смещением РОВНО −Δ (по гор. и по верт.)
        A = make_magnet("quad", 0.2, k1=1.0, center=(-0.1, 0, 0))
        st0 = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        base = mt.modified_drift(st0, A, make_magnet("dipole", 0.5, angle=0.02, center=(2, 0, 0)))
        # горизонтальный сдвиг dy=+1мм (глоб Y ≡ uhat при yaw=0): Δx_local ≈ −Δ, с малым сцеплением
        # через наклонённую на θ0/2 входную грань (член ∼θ0²) — потому допуск ослаблен до ~мкм.
        sh = mt.modified_drift(st0, A, make_magnet("dipole", 0.5, angle=0.02, center=(2, 1e-3, 0)))
        self.assertAlmostEqual((sh[0] - base[0]) * 1e3, -1.0, places=3)   # ≈ −Δ
        # вертикальный сдвиг dz=+1мм (глоб Z ≡ vhat, грань по вертикали НЕ наклонена): Δy_local = −Δ точно
        sv = mt.modified_drift(st0, A, make_magnet("dipole", 0.5, angle=0.02, center=(2, 0, 1e-3)))
        self.assertAlmostEqual((sv[2] - base[2]) * 1e3, -1.0, places=6)   # ровно −Δ


# ================================================================== 8. конвейер БД (gated)
@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "БД (lattice_configs.csv) не загружена")
class TestDBPipeline(unittest.TestCase):
    """build_line на загруженной БД. Gated — БД гитигнорится, может быть не загружена."""

    def test_build_line_lattice1(self):
        eff, phys, coup, trans = mt.load_db()
        mags = mt.build_line(1, eff, phys, coup, trans)
        self.assertGreater(len(mags), 0, "решётка 1 пуста")
        kinds = {m.kind for m in mags}
        self.assertIn("dq", kinds); self.assertIn("dipole", kinds); self.assertIn("quad", kinds)


# ================================================================== 9. C++ диффурное ядро (gated)
class TestCppKernel(unittest.TestCase):
    """C++ SynchrotronTracker.exe (диффурное ядро). Наличие бинаря + smoke-прогон."""
    EXE = next((p for p in (os.path.join(ROOT, "build", exe_name("SynchrotronTracker")),
                            os.path.join(ROOT, "build", "Release", exe_name("SynchrotronTracker")))
                if os.path.exists(p)), None)   # .exe только на Windows; Release — раскладка MSVC

    def test_exe_built(self):
        self.assertIsNotNone(self.EXE, "SynchrotronTracker.exe не собран (cmake --build build)")

    @unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "experiments.csv")), "нет experiments.csv")
    def test_smoke_run(self):
        if self.EXE is None:
            self.skipTest("бинарь не собран")
        import subprocess
        r = subprocess.run([self.EXE, "1"], cwd=ROOT, capture_output=True,
                           text=True, encoding="utf-8", errors="replace", timeout=300)
        self.assertEqual(r.returncode, 0, f"exe упал: {(r.stderr or r.stdout)[-400:]}")


# ================================================================== 10. Локальный экспорт .npz (анализатор)
class TestLocalNpz(unittest.TestCase):
    """write_local_npz ↔ read_local_npz round-trip: формат локальных координат для traj_analyzer."""

    def test_npz_roundtrip(self):
        import matrix_traj_export as mte, tempfile
        loc = [np.arange(8 * 3, dtype=float).reshape(3, 8), np.arange(8 * 2, dtype=float).reshape(2, 8) + 100]
        glo = [np.arange(3 * 3, dtype=float).reshape(3, 3), np.arange(3 * 2, dtype=float).reshape(2, 3) + 50]
        r = dict(locals=loc, globals=glo, surv_turns=[4, 2], n_turns=4,
                 elem_names=["плоскость инжекции", "A", "B"], closed=True, ring_len=123.4, lattice=1, group=2)
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, "loc.npz")
            mte.write_local_npz(p, r)
            d = mte.read_local_npz(p)
        self.assertEqual(len(d["parts"]), 2)
        self.assertEqual(len(d["parts"][0]["S"]), 3)
        self.assertEqual(len(d["parts"][1]["S"]), 2)
        self.assertEqual(d["n_turns"], 4)
        self.assertTrue(d["parts"][0]["survived"])          # surv_turns 4 >= n_turns 4
        self.assertFalse(d["parts"][1]["survived"])         # surv_turns 2 < 4
        self.assertAlmostEqual(d["ring_len"], 123.4)
        self.assertEqual(d["elem_names"], r["elem_names"])
        np.testing.assert_allclose(d["parts"][1]["x"], loc[1][:, 1])    # колонка x = столбец 1 записи
        np.testing.assert_allclose(d["parts"][0]["Z"], glo[0][:, 2])    # глоб. Z сохранён


class TestCouplingTable(unittest.TestCase):
    """Таблица БД coupling_matrices: формат прямых/обратных строк, реконструкция M·M⁺≈I, PM (fixed), UPSERT."""

    def test_matrix_table_roundtrip(self):
        import element_compare as ec
        M = np.array([[1.0, 2.0], [0.5, 3.0]])
        m = dict(name="T", kind="quad", outputs=["G", "S"], labels=["I1", "I2"], L=0.3, M=M,
                 imax=10.0, i_op=2.0, r2=0.99)
        rows = ec._matrix_table_rows(m)
        fw = {(r[5], r[6]): float(r[7]) for r in rows if r[4] == "forward"}
        Mr = np.array([[fw[(o, i)] for i in ("I1", "I2")] for o in ("G", "S")])
        np.testing.assert_allclose(Mr, M)                                   # forward восстанавливает M
        inv = {(r[5], r[6]): float(r[7]) for r in rows if r[4] == "inverse"}
        Mir = np.array([[inv[(w, o)] for o in ("G", "S")] for w in ("I1", "I2")])
        np.testing.assert_allclose(M @ Mir, np.eye(2), atol=1e-9)           # M·M⁺ ≈ I

    def test_pm_fixed_rows(self):
        import element_compare as ec
        rows = ec._matrix_table_rows(dict(name="DP", kind="dipole", outputs=["B", "G"], labels=["I"],
                                          L=0.4, pm=True, fixed={"B": 1e-4, "G": 2e-5}))
        self.assertTrue(all(r[4] == "fixed" for r in rows))
        self.assertEqual({r[5] for r in rows}, {"B", "G"})

    def test_upsert_preserves_others(self):
        import element_compare as ec, tempfile, csv as _csv
        row = lambda el: [el, "quad", "1", "1", "forward", "G", "I", "1.0", "coef/A", "0.3", "2.0", "10.0", "0.99"]
        old = ec.TB
        with tempfile.TemporaryDirectory() as td:
            ec.TB = td
            try:
                ec.upsert_coupling_table({"A": [row("A")], "B": [row("B")]})
                ec.upsert_coupling_table({"A": [row("A")[:2] + ["upd"] + row("A")[3:]]})   # обновить ТОЛЬКО A
                with open(os.path.join(td, ec.COUPLING_TABLE), encoding="utf-8") as f:
                    rr = list(_csv.reader(f, delimiter=";"))[1:]
                els = [r[0] for r in rr if r]
                self.assertEqual(els.count("B"), 1)                         # чужой B не тронут
                self.assertEqual(els.count("A"), 1)                         # A заменён, не продублирован
                self.assertIn("upd", [c for r in rr for c in r])            # новое значение A записано
            finally:
                ec.TB = old


class TestDesignCurrents(unittest.TestCase):
    """design_currents.py: разбор таблицы матриц, сопоставление обмоток с устройствами решётки, решение."""

    def test_label_and_base_mapping(self):
        """Имя элемента таблицы связи и device_type_id решётки должны сходиться к ОДНОЙ базе."""
        import design_currents as dc
        self.assertEqual(dc._base_name("SP_I1_1A_L166"), "SP_1A_L166")
        self.assertEqual(dc._base_name("DQ_I2_1A_L400"), dc._base_name("DQ_I1_1A_L400"))   # имя главной обмотки
        self.assertEqual(dc._base_name("corrector_c1_100A_L100"), "corrector_100A_L100")
        self.assertEqual(dc._winding_label("SP_I2_1A_L166"), "I2")
        self.assertEqual(dc._winding_label("corrector_main_100A_L100"), "main")
        self.assertEqual(dc._winding_label("OP_1A_L90"), "I")
        self.assertEqual(dc._norm_label("I_c1"), "c1")          # в таблице токи названы I_c1/I_main
        self.assertEqual(dc._norm_label("I1"), "I1")

    def test_attach_devices(self):
        """Обмотки восстанавливаются по базе+метке, а не обратным разбором имени.

        Отдельно проверяется вариант модели (`_sext` — DQ с секступольным киком): физический магнит и
        его обмотки те же, поэтому пометка варианта не должна мешать найти устройства.
        """
        import design_currents as dc
        lcrow = {d: {} for d in ("DQ_I1_1A_L400", "DQ_I2_1A_L400", "DQ_I3_1A_L400",
                                 "corrector_c1_100A_L100", "corrector_main_100A_L100")}
        mats = [dict(name="DQ_I1_1A_L400", kind="dq", pm=False, labels=["I1", "I2", "I3"]),
                dict(name="DQ_I1_1A_L400_sext", kind="dq", pm=False, labels=["I1", "I2", "I3"]),
                dict(name="corrector_100A_L100", kind="kicker", pm=False, labels=["I_c1", "I_main"])]
        dc.attach_devices(mats, lcrow)
        dq = ["DQ_I1_1A_L400", "DQ_I2_1A_L400", "DQ_I3_1A_L400"]
        self.assertEqual(mats[0]["devs"], dq)
        self.assertEqual(mats[1]["devs"], dq, "вариант _sext — тот же магнит, те же обмотки")
        self.assertEqual(mats[2]["devs"], ["corrector_c1_100A_L100", "corrector_main_100A_L100"])

    def test_load_coupling_roundtrip(self):
        """Прямые строки таблицы собираются обратно в матрицу; PM даёт ветку fixed."""
        import design_currents as dc, tempfile, csv as _csv
        hdr = ["element", "kind", "n_windings", "n_moments", "direction", "out_label", "in_label",
               "value", "units", "eff_length", "i_op", "imax", "r2"]
        rows = [["T", "sext", "2", "1", "forward", "S", "I1", "3.0", "coef/A", "0.2", "1.0", "5.0", "0.99"],
                ["T", "sext", "2", "1", "forward", "S", "I2", "0.5", "coef/A", "0.2", "1.0", "5.0", "0.99"],
                ["P", "dipole", "0", "1", "fixed", "B", "-", "-1e-4", "coef", "0.4", "-", "-", "-"]]
        with tempfile.TemporaryDirectory() as td:
            with open(os.path.join(td, dc.COUPLING_TABLE), "w", encoding="utf-8", newline="") as f:
                w = _csv.writer(f, delimiter=";")
                w.writerow(hdr); w.writerows(rows)
            mats = {m["name"]: m for m in dc.load_coupling(td)}
        self.assertEqual(mats["T"]["labels"], ["I1", "I2"])
        np.testing.assert_allclose(mats["T"]["M"], [[3.0, 0.5]])
        self.assertTrue(mats["P"]["pm"])
        self.assertAlmostEqual(mats["P"]["fixed"]["B"], -1e-4)

    @staticmethod
    def _mini_db(td):
        """Мини-база: матрица одного квадруполя + его строки в hardware_library и lattice_configs."""
        import csv as _csv
        import design_currents as dc

        def put(name, header, rows):
            with open(os.path.join(td, name), "w", encoding="utf-8", newline="") as f:
                w = _csv.writer(f, delimiter=";")
                w.writerow(header); w.writerows(rows)

        put(dc.COUPLING_TABLE,
            ["element", "kind", "n_windings", "n_moments", "direction", "out_label", "in_label",
             "value", "units", "eff_length", "i_op", "imax", "r2"],
            [["QP_X", "quad", "1", "1", "forward", "G", "I", "2.0", "coef/A", "0.3", "1.0", "9.0", "0.99"]])
        put("hardware_library.csv",
            ["device_type_id", "map_bin_path", "nominal_arg", "use_symmetry", "field_direction",
             "correct_coeff", "rot_x_deg", "rot_y_deg", "rot_z_deg", "turns_per_coil"],
            [["QP_X", "QP_X.bin", "100", "1", "1", "1", "0", "0", "0", "85"]])
        put("lattice_configs.csv",
            ["lattice_id", "device_type_id", "eff_length", "kind"], [["1", "QP_X", "0.3", "quad"]])

    def _patch(self, td, name, col, value):
        """Заменить значение в единственной строке таблицы мини-базы."""
        import csv as _csv
        p = os.path.join(td, name)
        with open(p, encoding="utf-8") as f:
            rows = list(_csv.reader(f, delimiter=";"))
        rows[1][rows[0].index(col)] = value
        with open(p, "w", encoding="utf-8", newline="") as f:
            _csv.writer(f, delimiter=";").writerows(rows)

    def test_staleness_uses_input_fingerprint(self):
        """Устаревание определяется по САМИМ входным данным, а не по дате файла.

        Дата ненадёжна: element_compare переписывает lattice_configs по ходу работы (временные и
        демо-конфигурации), и таблица оказывалась «старее» своих входов сразу после честной пересборки.
        """
        import design_currents as dc, tempfile
        with tempfile.TemporaryDirectory() as td:
            self._mini_db(td)
            self.assertIn("нет отпечатка", dc.staleness(td))          # снято прежней версией — не проверить

            dc.write_stamps(dc.current_hashes(td), td)                # как после честного измерения
            self.assertEqual(dc.staleness(td), "", "свежие матрицы не должны помечаться устаревшими")

            self._patch(td, "hardware_library.csv", "rot_z_deg", "180")   # развернули магнит — знак G меняется
            msg = dc.staleness(td)
            self.assertIn("устарели 1 из 1", msg)
            self.assertIn("QP_X", msg)

            self._patch(td, "hardware_library.csv", "rot_z_deg", "0")
            self._patch(td, "hardware_library.csv", "turns_per_coil", "42")  # витки матрицу НЕ меняют
            self.assertEqual(dc.staleness(td), "", "смена витков не должна помечать матрицу устаревшей")

            self._patch(td, "lattice_configs.csv", "eff_length", "0.4")      # геометрия — меняет
            self.assertIn("устарели", dc.staleness(td))


class TestReportPortability(unittest.TestCase):
    """Отчёт .docx собирается без Microsoft Office; формулы деградируют, а не роняют сборку."""

    def test_formula_falls_back_to_text(self):
        """Нет XSL Office (обычное дело вне Windows) ⇒ формула идёт текстом, исключения нет."""
        import latex2omml as lo
        from docx import Document
        old_xsl, old_tr, old_warn = lo._XSL, lo._transform, lo._warned
        lo._XSL, lo._transform, lo._warned = os.path.join(ROOT, "нет-такого.xsl"), None, True
        try:
            self.assertFalse(lo.available())
            p = Document().add_paragraph("Формула: ")
            lo.add_formula(p, r"\theta_0 = QP \cdot B \cdot L")          # не должно бросить
            self.assertIn("theta_0", p.text, "формула должна остаться в тексте абзаца")
        finally:
            lo._XSL, lo._transform, lo._warned = old_xsl, old_tr, old_warn

    def test_native_formula_when_xsl_present(self):
        """Есть XSL ⇒ вставляется НАТИВНЫЙ OMML-узел, а не текст."""
        import latex2omml as lo
        if not lo.available():
            self.skipTest("XSL Office недоступен на этой машине")
        from docx import Document
        p = Document().add_paragraph("Формула: ")
        lo.add_formula(p, r"k_1 = QP \cdot G")
        self.assertEqual(len([e for e in p._p if "oMath" in e.tag]), 1)
        self.assertEqual(p.text, "Формула: ", "нативная формула не должна дублироваться текстом")

    def test_convert_reports_missing_libreoffice(self):
        """Без LibreOffice конвертация даёт ВНЯТНУЮ ошибку, а не падение по месту вызова."""
        import office_convert as oc
        self.assertIn("odt", oc.FORMATS)
        with self.assertRaises(ValueError):                              # неизвестный формат
            oc.convert(__file__, to="xyz")
        old = oc.soffice
        oc.soffice = lambda: None
        try:
            with self.assertRaises(FileNotFoundError) as e:
                oc.convert(__file__, to="odt")
            self.assertIn("LibreOffice", str(e.exception))
        finally:
            oc.soffice = old


class TestGuiShell(unittest.TestCase):
    """Оболочка: модули импортируются, команды собираются правильно (без запуска Qt-приложения)."""

    def test_command_builders(self):
        from gui.runner import exe_job, py_job
        j = py_job("трек", "track3.py", "--mode", "kick-phys", "--exp", 1)
        self.assertTrue(j.argv[1].endswith("track3.py"))
        self.assertEqual(j.argv[2:], ["--mode", "kick-phys", "--exp", "1"])   # числа приводятся к строкам
        self.assertEqual(py_job("t", "s.py", "--a", None, "--b").argv[2:], ["--a", "--b"])  # None опускается
        self.assertTrue(exe_job("кик", "MatrixTracker", "--track", 1, 10).argv[0].endswith(
            "MatrixTracker.exe" if sys.platform == "win32" else "MatrixTracker"))

    def test_interrupted_run_recovers_tables(self):
        """Убитый на полпути прогон не должен оставлять базу испорченной.

        element_compare временно правит 4 таблицы и возвращает их в конце; при убийстве процесса
        (кнопка «Остановить», закрытие окна) возврат не отрабатывает. Метка прогона, пережившая
        процесс, означает прерывание — следующий старт обязан вернуть таблицы из копии.
        """
        import element_compare as ec
        import shutil
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            old_tb, old_guard, old_mark = ec.TB, ec.GUARD_DIR, ec.GUARD_MARK
            ec.TB = td
            ec.GUARD_DIR = os.path.join(td, "_cmp_backup")
            ec.GUARD_MARK = os.path.join(ec.GUARD_DIR, "in_progress")
            try:
                for f in ec.CSV4:                                   # исходное состояние таблиц
                    with open(os.path.join(td, f), "w", encoding="utf-8") as fh:
                        fh.write("a;b;kind\n1;2;quad\n")
                ec.guard_tables()
                self.assertTrue(os.path.exists(ec.GUARD_MARK), "метка прогона должна быть выставлена")
                with open(os.path.join(td, ec.CSV4[0]), "w", encoding="utf-8") as fh:
                    fh.write("a;b\n1;2\n")                          # рабочее состояние: часть столбцов потеряна
                self.assertEqual(ec.recover_tables(), ec.CSV4)      # следующий старт возвращает всё
                self.assertEqual(open(os.path.join(td, ec.CSV4[0]), encoding="utf-8").read(),
                                 "a;b;kind\n1;2;quad\n")
                self.assertFalse(os.path.exists(ec.GUARD_MARK), "после восстановления метки быть не должно")
                ec.guard_tables()
                ec.release_tables()                                 # штатное завершение — метка снята
                self.assertEqual(ec.recover_tables(), [], "без метки восстанавливать нечего")
            finally:
                ec.TB, ec.GUARD_DIR, ec.GUARD_MARK = old_tb, old_guard, old_mark

    def test_discover_projects_from_db(self):
        """Решётки из базы должны предлагаться проектами — иначе при первом запуске список пуст."""
        from gui import projects as pr
        found = pr.discover_projects()
        lats = {p.lattice for p in found}
        in_db = {lid for lid, n in pr.lattices_in_db().items() if n >= pr.MIN_ELEMENTS}
        self.assertEqual(lats, in_db, "предложены не все содержательные решётки базы")
        self.assertTrue(all(p.exp and p.group for p in found), "у проекта должны быть опыт и пучок")
        known = [p for p in found if p.lattice == min(lats)]                   # уже заведённые — не дублируем
        self.assertNotIn(min(lats), {p.lattice for p in pr.discover_projects(known)})

    def test_sub_lattice_label(self):
        """Подпись суб-решётки берётся из её описания и режется по смыслу, а не по символу."""
        from gui.projects import _short
        self.assertEqual(_short("Конфигурация 3 = SP1: 1 стандартный SP + начало SP2 (gid 1..79)"),
                         "Конфигурация 3 = SP1")
        self.assertEqual(_short("Короткое"), "Короткое")
        self.assertTrue(_short("слово " * 40).endswith("…"))

    def test_numeric_sort_key(self):
        """Сортировка таблиц по клику: числа сравниваются как числа (иначе «16» встаёт перед «2»)."""
        from gui.csvtable import NumericSortProxy as P
        self.assertEqual(P._num("16"), 16.0)
        self.assertEqual(P._num("1.00E-11"), 1e-11)
        self.assertEqual(P._num("-3,5"), -3.5)                # запятая-десятичная тоже читается
        self.assertIsNone(P._num("QP_L_1A_L295"))
        self.assertIsNone(P._num(""))

    def test_csv_roundtrip_preserves_format(self):
        """Правка таблицы не меняет разделитель, кодировку, состав столбцов и переводы строк (CRLF)."""
        from pathlib import Path
        from gui.csvtable import read_csv, write_csv
        import gui.paths as gp
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / "t.csv"
            p.write_bytes("a;b\r\n1.5;2\r\n".encode("utf-8-sig"))
            hdr, rows, enc = read_csv(p)
            self.assertEqual((hdr, rows, enc), (["a", "b"], [["1.5", "2"]], "utf-8-sig"))
            rows[0][0] = "9.75"
            old_tables = gp.TABLES
            gp.TABLES = Path(td)                                      # бэкап пишется рядом, не в БД проекта
            try:
                bak = write_csv(p, hdr, rows, enc)
            finally:
                gp.TABLES = old_tables
            self.assertEqual(p.read_bytes(), "a;b\r\n9.75;2\r\n".encode("utf-8-sig"))   # формат сохранён
            self.assertTrue(bak.exists(), "прежняя версия таблицы должна уходить в резервную копию")

            lf = Path(td) / "lf.csv"                              # таблица с LF не должна стать CRLF
            lf.write_bytes(b"a;b\n1;2\n")
            hdr, rows, enc = read_csv(lf)
            gp.TABLES = Path(td)
            try:
                write_csv(lf, hdr, rows, enc)
            finally:
                gp.TABLES = old_tables
            self.assertEqual(lf.read_bytes(), b"a;b\n1;2\n", "переводы строк исходного файла сохраняются")


class TestMecSigma(unittest.TestCase):
    """МЭС-σ: записанные веса, round-trip таблицы и ЗАКРЕПЛЁННАЯ точность по кику.

    Точность фиксируется рубежами, а не точными числами: карты и реестр могут пополняться, но
    ухудшение модели должно ронять тест.

    Метрика — расхождение ПО КИКУ (сумма модулей разностей к сумме модулей измеренного). Раньше
    здесь стояла ошибка ТОКА |λ−1|: она строже по кику, потому что в насыщении отклик пологий. Но
    как только набор пополнился картами НА ПОЛКЕ (2026-09-07, корректор), метрика тока там
    выродилась: на плоском участке d(ток)/d(кик) → ∞, и промах по кику в 5.7% давал 246% по току
    при полном охвате. Это свойство метрики, а не модели. Ошибка по кику ограничена и осмысленна на
    всём диапазоне, поэтому рубежи закреплены по ней. Охват (сколько карт модель ДОСТАЁТ по току)
    оставлен отдельным рубежом MISS — он ловит занижение полки, которое метрика кика не видит."""

    # (одиночные, КОМБИНАЦИИ), % ошибки тока. Пороги ПРИВЯЗАНЫ К НАБОРУ КАРТ: когда набор растёт,
    # числа могут законно ухудшиться — прежние были получены на данных, которых едва хватало на
    # число свободных параметров, и модель их почти интерполировала. Так и вышло 2026-09-04, когда
    # библиотека выросла со 109 карт до 154: у DQ (23→39 карт) комбинации по диполю ушли с 22.2%
    # до 29.9%, по квадруполю с 5.7% до 6.9%. Это не регресс кода, а исчезновение самообмана,
    # поэтому пороги подняты. Ужесточать их обратно можно только вслед за РЕАЛЬНЫМ улучшением
    # модели, а не ради зелёного теста.
    #   2026-09-05: то же самое повторилось у корректора. Его набор вырос с 37 карт до 56, из них
    # комбинаций стало 33 вместо 14, и ошибка тока по комбинациям ушла с 2.7% до 5.8%. Порог 3.5%
    # держался на четырнадцати комбинациях — их было меньше, чем весов, и модель их подгоняла
    # почти точно. Новое число получено на вдвое большем наборе и потому честнее прежнего.
    #   2026-09-07: набор корректора вырос с 56 до 84 карт (сняты колено, полка и скошенный
    # квадруполь). Рубежи по кику переставлены под новый набор с запасом ~30%. Числа честнее
    # прежних: добавлены трудные точки в колене и за краем, где модель работает хуже всего.
    # (одиночные, КОМБИНАЦИИ), % расхождения ПО КИКУ.
    #   2026-09-14: корректор отсюда УБРАН. Это тонкий прибор, штатно работающий при малых токах
    # (номинал 2 А) и насыщения не достигающий, поэтому модель насыщения для него не нужна: его
    # связь «ток → кик» — линейная матрица по картам рабочей области, её точность закреплена в
    # TestLinearCorrector. Карты корректора за рабочей областью из набора выведены.
    LIMITS = {("DQ", 1): (2.5, 3.0), ("DQ", 2): (2.5, 2.0), ("DQ", 3): (2.5, 6.0)}

    # Предельное число НЕДОСТИЖИМЫХ карт (одиночные, комбинации): карты, на которых модель не
    # дотягивается до измеренного кика ни при каком токе. Зафиксировано по состоянию на 2026-09-07.
    #   У корректора недостижима ПОЛОВИНА одиночных карт (11 из 23) и две трети комбинаций: во
    # всех случаях модель не дотягивается до измеренного кика даже при токе ×50, то есть её полка
    # ниже измеренной. Это настоящий изъян модели на верхних токах, а не настройка решателя.
    # Числа ниже — фиксация текущего положения, чтобы оно не ухудшалось незаметно; уменьшать их
    # надо вслед за исправлением полки.
    #   2026-09-07: сечение собственной утечки корректора поднято 0.002 → 0.05 кв.м. Разбор
    # показал, что после упрощения графа модель упиралась в потолок железа: воздушные пути — это
    # то, чем кик продолжает расти ПОСЛЕ насыщения, и слишком малая утечка эту роль не играла.
    # Охват вырос с 12/23 до 18/23 у диполя и с 12/23 до 18-20/23 у квадруполя и секступоля.
    #   Утечки катушки ПОЛЮСА и катушки ЯРМА разведены в разные классы: длины у них разные (0.17 и
    # 0.141 м), а при одном классе проницаемость класса записывалась по последней ветви и одна и
    # та же применялась ко всем шести. Прежний охват 20/23 был этой склейкой и завышен.
    #   2026-09-07: на обновлённом наборе (84 карты) полка теперь ДОСТИГАЕТСЯ почти везде на
    # одиночных (miss 1-2 из 47), а на комбинациях недостижимость держится (15-20 из 37) — там
    # много точек на полке, ниже которых садится модель. Это тот же изъян полки, теперь на большем
    # наборе; пороги фиксируют его, не пряча.
    #   2026-09-14: корректор убран вместе с LIMITS (см. выше).
    MISS = {("DQ", 1): (4, 0), ("DQ", 2): (0, 0), ("DQ", 3): (4, 0)}

    def _models(self):
        import mec_sigma as ms
        if not os.path.exists(ms.WEIGHTS):
            self.skipTest("нет таблицы весов — запустите mec_sigma.py --fit")
        return ms, ms.load()

    def test_weights_table_roundtrip(self):
        ms, M = self._models()
        self.assertEqual(set(M), set(ms.MODELS))
        for (el, order), (cir, w) in M.items():
            # длина вектора весов = 2 x число кик-классов ПОСЛЕ разделения по столбцу
            self.assertEqual(len(w), 2 * len(ms.kick_names(cir)))
            self.assertTrue(np.all(np.isfinite(w)))

    def test_zero_current_gives_zero_kick(self):
        """Точное физическое условие: без тока нет кика. Модель обязана его выполнять тождественно."""
        ms, M = self._models()
        for el, order in M:
            cir = M[(el, order)][0]
            k = ms.kick(el, order, {w: 0.0 for w in cir.coil_branches}, M)
            self.assertLess(abs(k), 1e-30)

    def test_kick_accuracy_holds(self):
        ms, M = self._models()
        import branch_network as bn
        for key, (lim_s, lim_c) in self.LIMITS.items():
            el, order = key
            cir, w = M[key]
            maps = bn.read_maps(el, order)
            for lbl, lim in (("одиночные", lim_s), ("комбинации", lim_c)):
                grp = [t for t in maps if (len(t[0]) == 1) == (lbl == "одиночные")]
                # Расхождение по кику: сумма модулей разностей к сумме модулей измеренного.
                num = sum(abs(complex(ms.basis(cir, c) @ w) - m) for c, m, _, _ in grp)
                den = sum(abs(m) for c, m, _, _ in grp)
                val = 100 * num / den if den else float("inf")
                self.assertLess(val, lim, f"{el} n={order} {lbl}: по кику {val:.1f}% (порог {lim})")

    def test_kick_reach_coverage(self):
        """Охват: сколько карт модель ДОСТАЁТ по току. Ловит занижение полки, которого метрика
        кика не видит (там модель мажет на проценты, а по току недобирает в разы)."""
        ms, M = self._models()
        import branch_network as bn
        for key, miss_lim in self.MISS.items():
            el, order = key
            maps = bn.read_maps(el, order)
            for j, lbl in ((0, "одиночные"), (1, "комбинации")):
                grp = [t for t in maps if (len(t[0]) == 1) == (lbl == "одиночные")]
                _, miss = ms.current_errors(el, order, grp, M)
                self.assertLessEqual(miss, miss_lim[j],
                                     f"{el} n={order} {lbl}: недостижимых {miss} из {len(grp)}")



class TestSaturationSolve(unittest.TestCase):
    """Единый решатель насыщения: таблица параметров и цикл прямая↔обратная задача."""

    def _tool(self):
        if not os.path.exists(os.path.join("workspace", "tables", "saturation_models.csv")):
            self.skipTest("нет saturation_models.csv — запустите build_saturation_models.py")
        import saturation_solve as ss
        return ss

    def test_models_table_covers_electromagnets(self):
        ss = self._tool()
        M = ss._models()
        for el in ("QP_sh", "QP_L", "SP", "OP", "DQ", "corrector"):
            self.assertIn(el, M, f"нет модели насыщения для {el}")

    def test_single_roundtrip(self):
        """Однообмоточный: ток → g → ток восстанавливается."""
        ss = self._tool()
        for el in ("QP_sh", "QP_L", "OP"):
            sym = ss._models()[el]["coefficient"]
            g = ss.forward(el, {"I1": 55.0})[sym][0]
            back = ss.inverse(el, {sym: g})
            self.assertAlmostEqual(next(iter(back.values())), 55.0, delta=0.5,
                                   msg=f"{el}: обратная задача не вернула ток")

    def test_dq_length_variant_same_coefficient(self):
        """Физкоэффициент длинонезависим: у исполнений g одинаков, интеграл разный."""
        ss = self._tool()
        g_base = ss.forward("DQ", {"I1": 100.0, "I2": 0.0, "I3": 0.0})["G"][0]
        g_l400 = ss.forward("DQ", {"I1": 100.0, "I2": 0.0, "I3": 0.0},
                            device="DQ_I1_1A_L400")["G"][0]
        self.assertAlmostEqual(g_base, g_l400, delta=1e-9)

    def test_dq_inverse_forward_consistent(self):
        """Многообмоточный: обратная задача под цель, прямая её воспроизводит."""
        ss = self._tool()
        cur = ss.inverse("DQ", {"B": 0.26, "G": -31.0}, fixed={"I3": 0.0})
        res = ss.forward("DQ", cur)
        self.assertAlmostEqual(res["B"][0], 0.26, delta=0.01)
        self.assertAlmostEqual(res["G"][0], -31.0, delta=0.5)



class TestLinearCorrector(unittest.TestCase):
    """Корректор: линейная матрица «ток → кик» в рабочей области и её закреплённая точность.

    Матрица подогнана только по одиночным картам; комбинации в подгонку не входят и проверяют
    суперпозицию. Метрика — расхождение по кику Σ|Δ|/Σ|изм| по порядку, нормальная и скошенная
    часть вместе. Пороги 2026-09-14 с запасом ~30% над достигнутым (одиночные / комбинации):
    B 1.00/1.07%, G 1.21/2.69%, S 1.55/5.40% на 30 одиночных и 15 комбинациях при токах до 2.5 А.
    Комбинации секступоля хуже прочих: секступоль у корректора мал и в комбинациях частично гасится."""

    LIMITS = {(1, "одиночные"): 1.3, (1, "комбинации"): 1.5,
              (2, "одиночные"): 1.6, (2, "комбинации"): 3.5,
              (3, "одиночные"): 2.0, (3, "комбинации"): 7.0}

    def _lc(self):
        if not os.path.exists(os.path.join("workspace", "tables", "linear_coupling.csv")):
            self.skipTest("нет linear_coupling.csv — запустите linear_coupling.py corrector")
        import linear_coupling as lc
        return lc

    def test_matrix_complete(self):
        """Шесть выходов на пять обмоток, числа конечны, скошенные у main запрещены симметрией."""
        lc = self._lc()
        K = lc.load("corrector")
        self.assertEqual(set(K), {"B", "Bv", "G", "Gv", "S", "Sv"})
        for sym, row in K.items():
            self.assertEqual(set(row), {"c1", "c2", "c3", "c4", "main"}, sym)
            self.assertTrue(all(np.isfinite(v) for v in row.values()), sym)
        for sym in ("Bv", "Gv", "Sv"):
            self.assertEqual(K[sym]["main"], 0.0, f"{sym}: main не создаёт скошенных компонент")

    def test_accuracy_holds(self):
        lc = self._lc()
        acc = lc.accuracy("corrector")
        for key, lim in self.LIMITS.items():
            self.assertIn(key, acc, f"нет карт группы {key}")
            val, n = acc[key]
            self.assertLess(val, lim, f"порядок {key[0]} {key[1]}: {val:.2f}% (порог {lim}, {n} карт)")

    def test_solver_roundtrip(self):
        """Через единый решатель: токи → физкоэффициенты → токи (пять целей на пять обмоток)."""
        self._lc()
        import saturation_solve as ss
        cur = {"c1": 1.2, "c2": -0.7, "c3": 0.4, "c4": 1.5, "main": 0.9}
        fwd = ss.forward("corrector", cur)
        back = ss.inverse("corrector", {k: fwd[k][0] for k in ("B", "Bv", "G", "Gv", "S")})
        for w, v in cur.items():
            self.assertAlmostEqual(back[w], v, delta=1e-6, msg=w)


class TestRingAnalysis(unittest.TestCase):
    """Алгоритмы ring_analysis (сверка с elegant) на синтетических отображениях с известным ответом — без БД."""

    @classmethod
    def setUpClass(cls):
        import ring_analysis as ra
        cls.ra = ra

    @staticmethod
    def _rot(beta, alpha, nu):
        mu = 2 * math.pi * nu; c, s = math.cos(mu), math.sin(mu)
        return np.array([[c + alpha * s, beta * s], [-(1 + alpha * alpha) / beta * s, c - alpha * s]])

    class _MapRing:
        """Трекер-заглушка: одно и то же нелинейное отображение за оборот (интерфейс Ring.track)."""
        def __init__(self, f):
            self.f = f

        def track(self, states, n_turns=1):
            out = []
            for st in states:
                z = np.asarray(st, float)
                for _ in range(n_turns):
                    z = self.f(z)
                out.append(z)
            return np.array(out)

    def test_naff_precision(self):
        k = np.arange(512)
        z = np.exp(2j * np.pi * 0.2201813 * k) + 0.05 * np.exp(2j * np.pi * 0.4403 * k)
        f, a = self.ra.naff(z, 2)
        self.assertAlmostEqual(f[0], 0.2201813, delta=1e-9)
        self.assertAlmostEqual(f[1], 0.4403, delta=1e-7)
        self.assertAlmostEqual(abs(a[1]), 0.05, delta=1e-5)

    def test_tbt_tunes_skips_forced_line(self):
        """Слабая плоскость: собственная частота, а не вынужденная линия 1−2νy от сильной."""
        k = np.arange(1024); nx, ny = 0.16, 0.2922
        tbt = np.zeros((1024, 6))
        tbt[:, 0] = 1e-5 * np.cos(2 * np.pi * nx * k) + 5e-5 * np.cos(2 * np.pi * 2 * ny * k)
        tbt[:, 1] = -1e-5 * np.sin(2 * np.pi * nx * k) - 5e-5 * np.sin(2 * np.pi * 2 * ny * k)
        tbt[:, 2] = 1e-3 * np.cos(2 * np.pi * ny * k); tbt[:, 3] = -1e-3 * np.sin(2 * np.pi * ny * k)
        qx, qy = self.ra.tbt_tunes(tbt)
        self.assertAlmostEqual(qx, nx, delta=1e-6)
        self.assertAlmostEqual(qy, ny, delta=1e-6)

    def test_twiss_from_one_turn_matrix(self):
        Mx = self._rot(6.9, 0.3, 0.22); My = self._rot(3.1, -0.2, 0.31)
        R = np.zeros((2, 6, 6))
        for i in range(2):
            R[i] = np.eye(6)
        R[1][:2, :2] = Mx; R[1][2:4, 2:4] = My
        R[0][:2, :2] = self._rot(6.9, 0.3, 0.1); R[0][2:4, 2:4] = self._rot(3.1, -0.2, 0.1)
        tw = self.ra.twiss_from_matrices(R, [1.0, 2.0], 2.0)
        self.assertAlmostEqual(tw["betax"][-1], 6.9, delta=1e-12)
        self.assertAlmostEqual(tw["alphay"][-1], -0.2, delta=1e-12)
        self.assertAlmostEqual(tw["nux"] % 1, 0.22, delta=1e-12)
        self.assertAlmostEqual(tw["nuy"] % 1, 0.31, delta=1e-12)

    def test_resonance_lines_irreducible(self):
        keys = {l[:3] for l in self.ra.resonance_lines(3, (0.15, 0.30, 0.25, 0.40))}
        self.assertEqual(keys, {(-2, 1, 0), (-1, 1, 0), (0, 3, 1), (1, 2, 1)})

    def test_stencil_exact_on_quartic_map(self):
        """4-точечная схема точна для многочленов до 4-й степени: матрица оборота = якобиан в точке."""
        c = np.array([0.3, -1.2, 0.7, 2.0, -0.5, 1.1])
        ring = self._MapRing(lambda z: z + c * z ** 2 + 0.1 * z ** 3 + 0.01 * np.roll(z, 1) ** 4)
        z0 = np.array([1e-3, -2e-3, 5e-4, 0.0, 1e-4, 2e-4])
        _, M = self.ra.one_turn_matrix(ring, orbit=z0, eps=1e-3)
        J = np.eye(6) + np.diag(2 * c * z0 + 0.3 * z0 ** 2)
        J += np.roll(np.diag(0.04 * z0 ** 3), 1, axis=0)          # d out_i / d z_{i−1}
        np.testing.assert_allclose(M, J, atol=1e-12)

    def test_damping_and_equilibrium_1d(self):
        """Затухающий поворот в каждой плоскости: τ из |λ|, равновесный эмиттанс = b/(2(1−λ²))·β-норм."""
        lam, beta, b = math.exp(-1e-3), 5.0, 1e-12
        M = np.zeros((6, 6))
        for i, nu in ((0, 0.22), (2, 0.31), (4, 0.004)):
            M[i:i + 2, i:i + 2] = lam * self._rot(beta, 0.0, nu)
        d = self.ra.damping_from_matrix(M, 1.0)
        for pl in ("x", "y", "s"):
            self.assertAlmostEqual(d[pl]["tau_turns"], 1e3, delta=1e-6)
        self.assertAlmostEqual(d["s"]["tune"], 0.004, delta=1e-12)
        B = np.zeros((6, 6)); B[1, 1] = b
        _, em = self.ra.equilibrium_moments(M, B)
        # диффузия по углу b в β-нормированных координатах: ε = β·b/(2·(1−λ²))
        self.assertAlmostEqual(em[0] / (beta * b / (2 * (1 - lam * lam))), 1.0, delta=1e-6)
        self.assertLess(em[1], 1e-30)


class TestRingMatrix(unittest.TestCase):
    """Аналитические матрицы элементов и стыков (ring_matrix) — против ядра и против симплектичности."""

    @classmethod
    def setUpClass(cls):
        import ring_matrix as rm
        cls.rm = rm

    # пятая координата — ПУТЬ (как у elegant), поэтому сопряжённая пара продольного блока — (−s, δ)
    J6 = np.zeros((6, 6))
    for _i in (0, 2):
        J6[_i, _i + 1] = 1.0; J6[_i + 1, _i] = -1.0
    J6[4, 5] = -1.0; J6[5, 4] = 1.0

    def _sympl(self, R, tol=1e-12):
        err = float(np.max(np.abs(R.T @ self.J6 @ R - self.J6)))
        self.assertLess(err, tol, f"нарушение симплектичности: {err:.2e}")

    def test_drift_and_quad_closed_form(self):
        rm = self.rm
        R = rm.drift(1.3)
        self.assertAlmostEqual(R[0, 1], 1.3); self.assertAlmostEqual(R[4, 5], 0.0)
        K, L = 2.462177415336321, 0.295                      # QF1A из .lte СИЛА
        Q = rm.quad(L, K)
        k = math.sqrt(K)
        self.assertAlmostEqual(Q[0, 0], math.cos(k * L), places=12)
        self.assertAlmostEqual(Q[1, 0], -k * math.sin(k * L), places=12)
        self.assertAlmostEqual(Q[3, 2], k * math.sinh(k * L), places=12)   # вертикаль дефокусируется
        self._sympl(Q)

    def test_bend_matches_kernel_body(self):
        """Матрица изгиба = якобиан тела ядра при сошедшемся интеграторе (40 срезов, порядок 6)."""
        import importlib
        os.environ["SS_SECTOR_NSL"] = "40"; os.environ["SS_SECTOR_ORDER"] = "6"
        try:
            importlib.reload(mt)
            m = make_magnet("dq", 1.0, angle=0.0247680586728496, k1=-1.960219550894432,
                            edge1=-0.0116800000099, edge2=-0.0116800000099)
            _, R = self.rm.element_matrix(m)
            base = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            F = np.zeros((4, 4)); h = 1e-6
            for j in range(4):
                sp = list(base); sp[j] += h; sm = list(base); sm[j] -= h
                op, om = m.body(tuple(sp)), m.body(tuple(sm))
                for i in range(4):
                    F[i, j] = (op[i] - om[i]) / (2 * h)
            err = float(np.max(np.abs(R[:4, :4] - F)))
            self.assertLess(err, 2e-6, f"матрица изгиба vs тело ядра: {err:.2e}")
            self._sympl(R, tol=1e-10)
        finally:
            os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
            importlib.reload(mt)

    def test_junction_matches_modified_drift(self):
        """Стык с НЕпараллельными гранями: аффинная матрица = точное отображение ядра (разности)."""
        from beam_injection import PlaneMagnet
        ax = np.array([1.0, 0.0, 0.0]); v = np.array([0.0, 0.0, 1.0])
        for alpha, roll, off in ((0.0, 0.0, 0.0), (5e-3, 0.0, 0.0), (5e-2, 0.0, 1e-3), (2e-2, 3e-2, 5e-4)):
            with self.subTest(alpha=alpha, roll=roll, off=off):
                A = PlaneMagnet(np.zeros(3), ax)
                b_ax = mt._rot(ax, v, alpha)
                B = PlaneMagnet(np.array([0.8, off, 0.0]), b_ax)
                if roll:
                    u = mt._rot(B.uhat, b_ax, roll)
                    B.uhat = u; B.vhat = np.cross(b_ax, u); B.fr_in = (B.center, B.axis, u, B.vhat)
                C, R = self.rm.junction_matrix(A, B)
                f = lambda z: np.array(mt.modified_drift(tuple(z), A, B), float)
                z0 = np.zeros(6); Cn = f(z0)
                Rn = np.zeros((6, 6)); h = 1e-6
                for k in range(6):
                    p = z0.copy(); p[k] += h; m_ = z0.copy(); m_[k] -= h
                    Rn[:, k] = (f(p) - f(m_)) / (2 * h)
                self.assertLess(float(np.max(np.abs(C - Cn))), 1e-12)
                self.assertLess(float(np.max(np.abs(R - Rn))), 1e-8)
                if alpha == 0.0 and off == 0.0 and roll == 0.0:
                    self.assertLess(float(np.max(np.abs(C[:4]))), 1e-15,
                                    "у параллельных граней поперечный C обязан быть нулевым")
                    self.assertAlmostEqual(C[4], 0.8, places=12, msg="C[4] — путь опорной частицы")

    def test_affine_composition(self):
        """Композиция (C, R) двух отображений = применение их подряд."""
        rm = self.rm
        C1 = np.array([1e-3, 2e-4, -5e-4, 1e-4, 0.0, 0.0]); R1 = rm.quad(0.3, 1.7)
        C2 = np.array([-2e-3, 1e-4, 1e-3, -3e-4, 0.0, 0.0]); R2 = rm.drift(0.9)
        C, R = rm.compose((C1, R1), (C2, R2))
        z = np.array([1e-3, 2e-4, -1e-3, 5e-5, 0.0, 1e-4])
        np.testing.assert_allclose(C + R @ z, C2 + R2 @ (C1 + R1 @ z), atol=1e-15)


_EXE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "build",
                    "MatrixTracker.exe" if sys.platform == "win32" else "MatrixTracker")


@unittest.skipUnless(os.path.exists(_EXE), "нет собранного MatrixTracker")
class TestRingCpp(unittest.TestCase):
    """Сервер кольца на C++-ядре (ring_cpp.CppRing) против Python-эталона ring_analysis.Ring."""

    LAT, LEAD, MEV, NSL = 2, 1.4335, 3000.0, 10              # решётка Сибири, короткий интегратор

    @classmethod
    def setUpClass(cls):
        import ring_analysis as ra, ring_cpp as rc
        cls.ra, cls.rc = ra, rc
        os.environ["SS_SECTOR_NSL"] = str(cls.NSL); os.environ["SS_SECTOR_ORDER"] = "4"
        try:
            cls.py = ra.Ring(cls.LAT, start_before_first_m=cls.LEAD, design_mev=cls.MEV)
        finally:
            os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
        cls.cpp = rc.CppRing(cls.LAT, start_before_first_m=cls.LEAD, design_mev=cls.MEV, nsl=cls.NSL, order=4)
        rng = np.random.default_rng(3)
        cls.st = rng.normal(0, [1e-3, 1e-4, 1e-3, 1e-4, 1e-3, 1e-3], (4, 6))

    @classmethod
    def tearDownClass(cls):
        cls.cpp.close()

    def test_geometry(self):
        self.assertEqual(self.cpp.n_mag, len(self.py.mags))
        self.assertAlmostEqual(self.cpp.circumference, self.py.circumference, places=9)
        self.assertEqual(self.cpp.names, self.py.names)

    def test_turn_and_elements(self):
        np.testing.assert_allclose(self.cpp.track(self.st, 2), self.py.track(self.st, 2), atol=1e-11)
        np.testing.assert_allclose(self.cpp.track_elements(self.st), self.py.track_elements(self.st), atol=1e-11)

    def _rf_pair(self, rf):
        """(C++-кольцо с ВЧ rf, Python-эталон с тем же rf) — Python-кольцо копия без пересборки."""
        import copy
        py = copy.copy(self.py); py.rf = dict(rf)
        cpp = self.rc.CppRing(self.LAT, start_before_first_m=self.LEAD, design_mev=self.MEV, nsl=self.NSL,
                              order=4, rf=dict(rf))
        return cpp, py

    def test_rf_frequency_parity(self):
        # расстроенный резонатор (freq_hz ≠ h·c/C): C++ и Python считают одну длину ВЧ-оборота C_RF = h·c/f
        f0 = 75 * 299792458.0 / self.py.circumference
        for freq in (None, f0 + 300.0):
            rf = dict(volt_mv=2.5, harmonic=75, phase_deg=170.0, **({"freq_hz": freq} if freq else {}))
            cpp, py = self._rf_pair(rf)
            try:
                np.testing.assert_allclose(cpp.track(self.st, 3), py.track(self.st, 3), atol=1e-11)
                self.assertAlmostEqual(cpp.rf_circumference(), py.rf_circumference(), places=12)
            finally:
                cpp.close()

    def test_rf_frequency_shifts_energy(self):
        # частота ВЧ — актуатор энергии: равновесная δ из 6D-замкнутой орбиты ≈ −Δf/(η·f), η = αc − 1/γ²
        rf = dict(volt_mv=2.5, harmonic=75, phase_deg=180.0)          # без СИ синхронная фаза — 180°
        cpp, _ = self._rf_pair(rf)
        try:
            f0 = cpp.rf_frequency()
            self.assertAlmostEqual(f0, 75 * 299792458.0 / cpp.circumference, delta=1e-6)
            co0 = self.ra.closed_orbit(cpp, dims=6)
            df = 2000.0
            cpp.set_rf(freq_hz=f0 + df)
            co1 = self.ra.closed_orbit(cpp, dims=6)
            alphac = self.ra.twiss(self.cpp)["alphac"]
            gamma = self.MEV / 0.51099895069
            expect = -df / ((alphac - 1.0 / gamma ** 2) * f0)
            self.assertAlmostEqual((co1[5] - co0[5]) / expect, 1.0, delta=0.03)
            # на лету = свежее кольцо с той же частотой; переживает RESET (apply_perturb)
            fresh = self.rc.CppRing(self.LAT, start_before_first_m=self.LEAD, design_mev=self.MEV, nsl=self.NSL,
                                    order=4, rf={**rf, "freq_hz": f0 + df})
            try:
                cpp.apply_perturb(kicks=[(3, 1e-5)]); cpp.apply_perturb()
                np.testing.assert_array_equal(cpp.track(self.st, 2), fresh.track(self.st, 2))
            finally:
                fresh.close()
            cpp.set_rf(freq_hz=0.0)                                   # ≤ 0 ⇒ снова h·c/C, бит-в-бит исходное
            plain, _ = self._rf_pair(rf)
            try:
                np.testing.assert_array_equal(cpp.track(self.st, 2), plain.track(self.st, 2))
            finally:
                plain.close()
        finally:
            cpp.close()

    def test_kick_and_yaw(self):
        import copy
        j = next(i for i, k in enumerate(self.py.kinds) if k == "sext")
        pert = self.ra.kicked(self.py, j, 1e-4)
        ck = self.cpp.with_kick(j, 1e-4)
        np.testing.assert_allclose(ck.track(self.st), pert.track(self.st), atol=1e-11)
        ck.close()
        q = next(i for i, k in enumerate(self.py.kinds) if k == "quad" and i > 5)
        pert = copy.copy(self.py); pert.mags = list(self.py.mags)
        pert.mags[q] = self.ra.yaw_element(copy.copy(self.py.mags[q]), 2e-4)
        cy = self.cpp.with_yaw(q, 2e-4)
        np.testing.assert_allclose(cy.track(self.st), pert.track(self.st), atol=1e-11)
        cy.close()

    def test_field_perturbation_on_zero_element(self):
        """DK2 на элементе с НУЛЕВЫМ полем (секступоль Сибири, mat_S=0) действует и совпадает с Python-эталоном
        (раньше флаг кик-тела замерзал при сборке, и возмущение терялось)."""
        import copy
        j = next(i for i, m in enumerate(self.py.mags) if self.py.kinds[i] == "sext" and not m.k2 and m.Leff > 0)
        base = self.cpp.track(self.st)
        try:
            got = self.cpp.apply_perturb(dk2s=[(j, 50.0)]).track(self.st)
        finally:
            self.cpp.apply_perturb()                            # RESET — прочие тесты класса на чистом кольце
        self.assertGreater(float(np.max(np.abs(got - base))), 1e-9)
        m = copy.copy(self.py.mags[j]); m.k2 = 50.0; m.nonzero = True; m.kick_sector = True
        m._nsl_sector = self.NSL; m._yoshida_dk = mt._yoshida(4)
        pert = copy.copy(self.py); pert.mags = list(self.py.mags); pert.mags[j] = m
        np.testing.assert_allclose(got, pert.track(self.st), atol=1e-11)

    def test_pieces_follow_perturbation(self):
        """Срезы (track_pieces/piece_meta) видят возмущения apply_perturb и совпадают с Python-эталоном: те же
        DK/DK2/KICK на копиях магнитов, срезы изгибов по ring_analysis._bend_slices (раньше копии магнитов в
        срезах собирались один раз при сборке/RESET, и возмущения в них не попадали)."""
        import copy
        mags, kinds = self.py.mags, self.py.kinds
        jb = next(i for i, m in enumerate(mags) if m.sector)                        # изгиб: DK + DK2 (b1, b2 сектора)
        jq = next(i for i, k in enumerate(kinds) if k == "quad" and i > 5)           # квадруполь: DK (кик-тело)
        js = next(i for i, m in enumerate(mags) if kinds[i] == "sext" and not m.k2 and m.Leff > 0)   # KICK (ящик)
        dkb, dkq, dk2b, th = 0.3, 0.2, 5.0, 1e-4
        try:
            self.cpp.apply_perturb(kicks=[(js, th)], dks=[(jb, dkb), (jq, dkq)], dk2s=[(jb, dk2b)])
            meta, got = self.cpp.piece_meta(), self.cpp.track_pieces(self.st)
        finally:
            self.cpp.apply_perturb()                            # RESET — прочие тесты класса на чистом кольце
        pert = list(mags)
        for j, dk, dk2 in ((jb, dkb, dk2b), (jq, dkq, 0.0)):
            m = copy.copy(mags[j]); m.k1 += dk; m.k1v += dk; m.k2 += dk2; pert[j] = m
        pieces = []                                             # (магнит, стык перед куском, кусок)
        for j, m in enumerate(pert):
            pieces += [(j, i == 0, c) for i, c in enumerate(self.ra._bend_slices(m) if m.sector else [m])]
        self.assertEqual(len(meta), len(pieces))
        self.assertEqual([(a, b, c) for a, b, c, *_ in meta], [(j, jn, bool(c.sector)) for j, jn, c in pieces])
        np.testing.assert_allclose(np.array([r[3:] for r in meta]),
                                   np.array([(c.Leff, c.theta0, c.k1, c.k2, c.k3) for _, _, c in pieces]),
                                   rtol=1e-12, atol=1e-15)
        m = copy.copy(mags[js]); m.theta0 = th; m.nonzero = True; m.xref = None   # кик корректора, как ra.kicked
        pert[js] = m
        pieces = [(j, jn, m if j == js else c) for j, jn, c in pieces]
        ref = np.full_like(got, np.nan)
        for p, z in enumerate(self.st):
            cur = tuple(map(float, z)); ref[p, 0] = cur; prev = self.py.plane
            for q, (j, jn, c) in enumerate(pieces):
                if jn:
                    cur = mt.modified_drift(cur, prev, pert[j]); prev = pert[j]
                cur = c.body(cur); ref[p, q + 1] = cur
        np.testing.assert_allclose(got, ref, atol=1e-11)

    def test_pieces_match_elements_after_perturb(self):
        """Срезы после любых возмущений, включая KICK/VKICK на изгибе (у Python-эталона такого нет), сшиваются
        в магниты: состояние после последнего среза = выход магнита в track_elements. Доп. кик у секторного
        кладётся только на последний срез, иначе он входил бы N раз."""
        jb = next(i for i, m in enumerate(self.py.mags) if m.sector)
        jq = next(i for i, k in enumerate(self.py.kinds) if k == "quad" and i > 5)
        try:
            self.cpp.apply_perturb(kicks=[(jb, 1e-4)], vkicks=[(jb, -2e-4)], dsks=[(jq, 0.05)], dk3s=[(jq, 30.0)])
            meta, P, E = self.cpp.piece_meta(), self.cpp.track_pieces(self.st), self.cpp.track_elements(self.st)
        finally:
            self.cpp.apply_perturb()
        last = {}
        for q, r in enumerate(meta):
            last[r[0]] = q + 1
        np.testing.assert_allclose(P[:, [last[j] for j in range(self.cpp.n_mag)]], E[:, :-1], atol=1e-11)

    def test_description_follows_apply_perturb(self):
        """Описание магнитов CppRing.mags (для ring_matrix/ring_response) несёт возмущения apply_perturb: та же
        физика, что у сервера — Python-трекинг описания = C++; после RESET описание снова чистое."""
        import copy
        mags, kinds = self.py.mags, self.py.kinds
        jb = next(i for i, m in enumerate(mags) if m.sector)
        jq = next(i for i, k in enumerate(kinds) if k == "quad" and i > 5)
        js, js2 = [i for i, m in enumerate(mags) if kinds[i] == "sext" and not m.k2 and m.Leff > 0][:2]
        try:
            self.cpp.apply_perturb(kicks=[(js, 1e-4)], dks=[(jb, 0.3), (jq, 0.2)], dk2s=[(jb, 5.0), (js2, 50.0)])
            d = self.cpp.mags
            got = self.cpp.track(self.st)
        finally:
            self.cpp.apply_perturb()
        self.assertAlmostEqual(d[jb].k1, mags[jb].k1 + 0.3, places=14)
        self.assertAlmostEqual(d[jq].k1v, mags[jq].k1v + 0.2, places=14)
        self.assertEqual(d[jb].theta0, mags[jb].theta0)
        self.assertTrue(d[js2].kick_sector)                     # нулевой секступоль с k2 — кик-тело, как в C++
        self.assertEqual(d[js].theta0, 1e-4)                    # кик «ящика» = θ0 box-тела (ra.kicked)
        pr = copy.copy(self.py); pr.mags = d
        np.testing.assert_allclose(got, pr.track(self.st), atol=1e-11)
        self.assertEqual([m.k1 for m in self.cpp.mags], [m.k1 for m in mags])

    def test_base_perturbations_survive_reset(self):
        """Возмущения создания (with_*/perturbed) — часть кольца: apply_perturb() их не сбрасывает, описание их
        держит. KICK/VKICK изгиба не трогает его θ0 (кривизну h для ring_matrix), VKICK — не гор. θ0."""
        import copy
        mags, kinds = self.py.mags, self.py.kinds
        jb = next(i for i, m in enumerate(mags) if m.sector)
        q = next(i for i, k in enumerate(kinds) if k == "quad" and i > 5)
        js = next(i for i, m in enumerate(mags) if kinds[i] == "sext" and not m.k2 and m.Leff > 0)
        yawed = self.ra.yaw_element(copy.copy(mags[q]), 2e-4)
        r = self.cpp.with_yaw(q, 2e-4)
        try:
            ref = r.track(self.st)
            r.apply_perturb(kicks=[(jb, 1e-4)], vkicks=[(jb, -2e-4)])
            m = r.mags[jb]
            self.assertEqual(m.theta0, mags[jb].theta0)
            self.assertEqual((m.hkick, m.vkick), (1e-4, -2e-4))
            np.testing.assert_allclose(r.mags[q].fr_in[1], yawed.fr_in[1], atol=1e-15)
            r.apply_perturb()
            np.testing.assert_array_equal(r.track(self.st), ref)
        finally:
            r.close()
        v = self.cpp.with_vkick(js, 1e-4)
        try:
            self.assertEqual((v.mags[js].theta0, v.mags[js].vkick), (mags[js].theta0, 1e-4))
        finally:
            v.close()

    def test_junction_command(self):
        from beam_injection import PlaneMagnet
        ax = np.array([1.0, 0.0, 0.0]); v = np.array([0.0, 0.0, 1.0])
        A = PlaneMagnet(np.zeros(3), ax)
        B = PlaneMagnet(np.array([0.8, 1e-3, 0.0]), mt._rot(ax, v, 5e-2))
        got = self.cpp.junction(self.st, A.fr_out, B.fr_in)
        ref = np.array([mt.modified_drift(tuple(z), A, B) for z in self.st])
        np.testing.assert_allclose(got, ref, atol=1e-14)


class TestRingMap(unittest.TestCase):
    """Отображения высших порядков (ring_map): алгебра рядов, точные стыки, хроматичность."""

    @classmethod
    def setUpClass(cls):
        import ring_map
        cls.m = ring_map
        cls.da = ring_map.DA(6, 4)

    def _eval(self, a, z):
        return sum(c * np.prod([zi ** e for zi, e in zip(z, ex)]) for c, ex in zip(a, self.da.monos))

    def test_series_functions(self):
        da = self.da
        a = da.const(1.3) + 0.2 * da.var(0) - 0.1 * da.var(5) + 0.05 * da.mul(da.var(1), da.var(1))
        z = np.array([1e-2, -2e-2, 0, 0, 0, 3e-2])
        f = lambda z: 1.3 + 0.2 * z[0] - 0.1 * z[5] + 0.05 * z[1] ** 2
        self.assertAlmostEqual(self._eval(da.sqrt(a), z), math.sqrt(f(z)), places=10)
        self.assertAlmostEqual(self._eval(da.inv(a), z), 1 / f(z), places=10)

    def test_linear_flow_equals_closed_form(self):
        """Линейная часть точного решения квадруполя = матрица ring_matrix.quad."""
        import ring_matrix as rm
        da = self.da
        K, L = 2.462177415336321, 0.295
        M = self.m.flow(da, self.m.field_multipole(da, K, K, 0.0, 0.0), L)
        R = np.array([[M[i][1 + j] for j in range(4)] for i in range(4)])
        np.testing.assert_allclose(R, rm.quad(L, K)[:4, :4], atol=1e-13)

    def test_junction_series_matches_kernel(self):
        """Ряд клинового стыка = точный modified_drift ядра с остатком порядка N+1."""
        from beam_injection import PlaneMagnet
        da = self.da
        ax = np.array([1.0, 0.0, 0.0]); v = np.array([0.0, 0.0, 1.0])
        A = PlaneMagnet(np.zeros(3), ax)
        B = PlaneMagnet(np.array([0.8, 1e-3, 0.0]), mt._rot(ax, v, 5e-2))
        J = self.m.junction(da, da.identity(), A.fr_out, B.fr_in)
        for scale in (1e-3, 5e-4):
            z = scale * np.array([1.0, 0.5, -0.7, 0.3, 0.0, 0.8])
            ref = np.array(mt.modified_drift(tuple(z), A, B))
            got = np.array([self._eval(c, z) for c in J])
            self.assertLess(float(np.max(np.abs(got - ref))), 50 * scale ** 5)

    def test_chromaticity_matches_tracking(self):
        """ξ1 кольца Сибири по отображению = производная ν(δ) трекинга (разности на ±1e-4)."""
        import ring_analysis as ra
        os.environ["SS_SECTOR_NSL"] = "40"; os.environ["SS_SECTOR_ORDER"] = "4"
        try:
            ring = ra.Ring(2, start_before_first_m=1.4335, design_mev=3000.0)
        finally:
            os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
        ch = self.m.RingMap(ring, 3).chromaticity()
        qp, qm = ra.tunes_point(ring, 1e-4), ra.tunes_point(ring, -1e-4)
        self.assertAlmostEqual(ch["xi_x"][1], (qp[0] - qm[0]) / 2e-4, delta=2e-3)
        self.assertAlmostEqual(ch["xi_y"][1], (qp[1] - qm[1]) / 2e-4, delta=2e-3)

    def test_normal_form_octupole(self):
        """Нормальная форма на пробном отображении «поворот + тонкий октуполь» = формула −K₃Lβ²/(16π)."""
        da = self.m.DA(4, 4)
        mux, muy, bx, by, k3 = 2 * math.pi * 0.31, 2 * math.pi * 0.22, 5.0, 3.0, 100.0
        R = np.zeros((4, 4))
        R[:2, :2] = [[math.cos(mux), bx * math.sin(mux)], [-math.sin(mux) / bx, math.cos(mux)]]
        R[2:, 2:] = [[math.cos(muy), by * math.sin(muy)], [-math.sin(muy) / by, math.cos(muy)]]
        x, px, y, py = da.identity()
        kick = [x, px + (k3 / 6.0) * da.mul(da.mul(x, x), x), y, py]
        M = [sum(R[i, j] * kick[j] for j in range(4)) for i in range(4)]
        nf = self.m.normal_form_of_map(da, M)
        self.assertAlmostEqual(nf["nux"], 0.31, places=9)
        self.assertAlmostEqual(nf["nuy"], 0.22, places=9)
        self.assertAlmostEqual(nf["dnux/dJx"], -k3 * bx * bx / (16 * math.pi), delta=1e-6)
        for k in ("dnux/dJy", "dnuy/dJx", "dnuy/dJy"):
            self.assertAlmostEqual(nf[k], 0.0, places=9, msg=f"{k}: октуполь по x не трогает вертикаль")

    def test_normal_form_matches_tracking(self):
        """∂ν/∂J кольца Сибири по нормальной форме = наклон ν(J) по трекингу C++-ядра."""
        import ring_analysis as ra, ring_matrix as rm
        os.environ["SS_SECTOR_NSL"] = "40"; os.environ["SS_SECTOR_ORDER"] = "4"
        try:
            ring = ra.Ring(2, start_before_first_m=1.4335, design_mev=3000.0)
        finally:
            os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
        nf = self.m.normal_form(ring, 4)
        tw = rm.twiss(ring); bx, ax = float(tw["betax"][-1]), float(tw["alphax"][-1])
        co = ra.closed_orbit(ring)
        q = []
        for a in (2e-3, 4e-3):
            st = co.copy(); st[0] += a; st[2] += 1e-5
            tb = ring.track_tbt([st], 512)[0]
            q.append((a * a * (1 + ax * ax) / (2 * bx), ra.tbt_tunes(tb, (bx, 1.0), (ax, 0.0), orbit=co)[0]))
        slope = (q[1][1] - q[0][1]) / (q[1][0] - q[0][0])
        self.assertAlmostEqual(nf["dnux/dJx"] / slope, 1.0, delta=2e-3)
        self.assertAlmostEqual(nf["dnux/dJy"], nf["dnuy/dJx"], delta=1e-6,
                               msg="перекрёстные коэффициенты обязаны совпадать")


class TestRadiationIntegrals(unittest.TestCase):
    """Интегралы излучения (ring_matrix.radiation_integrals) — против трекинга ядра и формул."""

    @classmethod
    def setUpClass(cls):
        import ring_analysis as ra, ring_matrix as rm
        cls.ra, cls.rm = ra, rm
        os.environ["SS_SECTOR_NSL"] = "40"; os.environ["SS_SECTOR_ORDER"] = "4"
        try:
            cls.ring = ra.Ring(2, start_before_first_m=1.4335, design_mev=3000.0)
            cls.ring_sr = ra.Ring(2, start_before_first_m=1.4335, design_mev=3000.0, sr=True)
        finally:
            os.environ.pop("SS_SECTOR_NSL", None); os.environ.pop("SS_SECTOR_ORDER", None)
        cls.ri = rm.radiation_integrals(cls.ring)

    def test_u0_matches_tracking(self):
        """U₀ по интегралам = потеря энергии осевой частицы за оборот в ядре с излучением."""
        u0_track = -float(self.ring_sr.track([np.zeros(6)])[0, 5]) * 3000.0
        self.assertAlmostEqual(self.ri["U0"] / u0_track, 1.0, delta=1e-3)

    def test_alphac_matches_matrices(self):
        """αc = I₁/C совпадает с уплотнением орбит из перемножения матриц."""
        self.assertAlmostEqual(self.ri["alphac"] / float(self.rm.twiss(self.ring)["alphac"]), 1.0, delta=1e-7)

    def test_robinson_sum(self):
        """Теорема Робинсона: Jx + Jy + Js = 4."""
        self.assertAlmostEqual(self.ri["Jx"] + self.ri["Jy"] + self.ri["Jdelta"], 4.0, places=12)

    def test_edge_term_present(self):
        """Краевой член I₄: у решётки с наклонными гранями он меняет I₄ (у Сибири грани прямые)."""
        import copy
        ring = copy.copy(self.ring); ring.mags = [copy.copy(m) for m in self.ring.mags]
        i4_0 = self.ri["I4"]
        for m in ring.mags:
            if getattr(m, "sector", False):
                m.kEdge1 = m.kEdge2 = 0.05 * (-m.theta0 / m.Leff)          # tg(E) = 0.05
        i4 = self.rm.radiation_integrals(ring)["I4"]
        self.assertGreater(abs(i4 - i4_0), 1e-4, "краевой член обязан входить в I₄")


class TestDipoleFringe(unittest.TestCase):
    """Фриндж диполя (эффект второго порядка, elegant CSBEND HGAP/FINT, edge_effects=1).

    Мягкий спад поля поворачивает орбиту ВНУТРИ зоны спада ⇒ эффективный угол грани для вертикали
    E → E−ψ, ψ = L_f/ρ·(1+sin²E)/cos E, L_f = 2·HGAP·FINT; горизонталь — прежняя tg E/ρ. В треке ψ/(1+δ).
    Эталоны R43 — elegant 2026.3 (CSBEND L=1, ANGLE=0.1, N_SLICES=5, HGAP=0.01, FINT=0.5 ⇒ L_f=0.01 м)."""
    L, ANG, LF = 1.0, 0.1, 0.01
    ELEGANT_R43 = {0.0: 2.000100451268587e-04, 0.05: -9.783053758787746e-03}   # E1=E2=E → R43 elegant

    def _m(self, E, lf):
        # знак БД: edge = −E, θ0 = −ANGLE (см. s2lattice_to_db)
        return make_magnet("dipole", self.L, angle=self.ANG, edge1=-E, edge2=-E, lf1=lf, lf2=lf)

    def test_lte_fringe_lengths_elegant_defaults(self):
        """L_f = 2·HGAP·FINTi по умолчаниям elegant: HGAP=0 ⇒ 0; FINT=0.5; FINT1/2<0 ⇒ FINT; EDGE*_EFFECTS=0 ⇒ 0."""
        f = mt.lte_fringe_lengths
        self.assertEqual(f({}), (0.0, 0.0))                                   # HGAP по умолчанию 0
        self.assertEqual(f({"FINT": 0.7}), (0.0, 0.0))
        self.assertEqual(f({"HGAP": 0.02}), (0.02, 0.02))                     # FINT по умолчанию 0.5
        self.assertEqual(f({"HGAP": 0.02, "FINT": 0.25, "FINT2": 0.0}), (0.01, 0.0))
        self.assertEqual(f({"HGAP": 0.02, "FINT1": 0.4, "FINT2": -1.0}), (0.016, 0.02))
        self.assertEqual(f({"HGAP": 0.02, "EDGE1_EFFECTS": 0}), (0.0, 0.02))

    def test_importer_edge_effects_off(self):
        """Импорт .lte: EDGE1_EFFECTS=0 выключает ВСЮ входную грань (и tg E, и фриндж), как elegant."""
        import s2lattice_to_db as s2
        e = {"type": "CSBEND", "E1": 0.03, "E2": 0.04, "HGAP": 0.015, "FINT": 0.5, "EDGE1_EFFECTS": 0.0}
        self.assertEqual(s2.edges_of(e), (0.0, 0.04, 0.0, 0.015))

    def test_hard_edge_bitwise(self):
        """Без L_f (нет столбцов / нули) — тело бит-в-бит прежнее: верт. коэффициент = kEdge без пересчёта tg."""
        st = (1e-3, 2e-4, 0.5e-3, -1e-4, 0.0, 1e-3)
        for E in (0.0, 0.05):
            a, b = self._m(E, None), self._m(E, 0.0)
            self.assertEqual(a.body(st), b.body(st))
            self.assertEqual(b.psi1, 0.0)
            self.assertEqual(mt._edge_kv(b.edge1, 0.0, b.rinv, b.kEdge1, 1e-3), b.kEdge1)

    def test_r43_matches_elegant(self):
        """Вертикальная линза граней с фринджем = elegant (R43), горизонталь фринджем не меняется (R21)."""
        base = (0.0,) * 6
        for E, r43_el in self.ELEGANT_R43.items():
            M = jacobian4(self._m(E, self.LF).body, base)
            M0 = jacobian4(self._m(E, 0.0).body, base)
            self.assertAlmostEqual(M[3, 2], r43_el, delta=1e-10, msg=f"E={E}: R43")
            self.assertAlmostEqual(M[1, 0], M0[1, 0], delta=1e-12, msg=f"E={E}: R21 не зависит от фринджа")

    def test_defocus_sign_independent_of_bend_direction(self):
        """При E=0 фриндж ДЕфокусирует по вертикали (R43>0) при изгибе в любую сторону: ∝ 1/ρ²."""
        for ang in (self.ANG, -self.ANG):
            m = make_magnet("dipole", self.L, angle=ang, lf1=self.LF, lf2=self.LF)
            r43 = jacobian4(m.body, (0.0,) * 6)[3, 2]
            self.assertAlmostEqual(r43, self.ELEGANT_R43[0.0], delta=1e-10, msg=f"ANGLE={ang}")

    def test_chromatic_psi(self):
        """Хроматичность края: ψ/(1+δ). R43(δ) = −k₁(δ) − k₂(δ) + k₁k₂·L/(1+δ), kᵢ = tg(E−ψ/(1+δ))/ρ."""
        E, d = 0.05, 2e-3
        m = self._m(E, self.LF)
        k = mt._edge_kv(m.edge1, m.psi1, m.rinv, m.kEdge1, d)
        expect = -2.0 * k + k * k * self.L / (1.0 + d)
        r43 = jacobian4(m.body, (0.0, 0.0, 0.0, 0.0, 0.0, d))[3, 2]
        self.assertAlmostEqual(r43, expect, delta=1e-9)
        k0 = mt._edge_kv(m.edge1, m.psi1, m.rinv, m.kEdge1, 0.0)
        self.assertGreater(abs(k - k0), 1e-8)                                 # δ-зависимость реально есть

    def test_ring_matrix_edge(self):
        """Аналитическая матрица изгиба (ring_matrix) берёт ту же вертикальную линзу tg(E−ψ)/ρ, что тело."""
        import ring_matrix as rm
        for E in (0.0, 0.05):
            m = self._m(E, self.LF)
            _, R = rm.element_matrix(m)
            self.assertAlmostEqual(R[3, 2], self.ELEGANT_R43[E], delta=1e-9, msg=f"E={E}")

    def test_ring_map_edge_series(self):
        """Ряд края в ring_map: коэффициент y → −k(0), y·δ → −∂k/∂δ (хроматичность ψ/(1+δ))."""
        import ring_map as rmap
        da = rmap.DA(6, 3)
        m = self._m(0.05, self.LF)
        st = rmap.edge(da.identity(), m.kEdge1, da, m.edge1, m.psi1, m.rinv)
        py = st[3]
        k = lambda d: mt._edge_kv(m.edge1, m.psi1, m.rinv, m.kEdge1, d)
        h = 1e-5
        self.assertAlmostEqual(da.coef(py, (0, 0, 1, 0, 0, 0)), -k(0.0), delta=1e-15)
        self.assertAlmostEqual(da.coef(py, (0, 0, 1, 0, 0, 1)), -(k(h) - k(-h)) / (2 * h), delta=1e-9)


class TestMultipolesAnyKind(unittest.TestCase):
    """Мультиполи K1/K2/K3 — свойство ПОЛЯ, а не вида: kind задаёт только геометрию тела (изгиб = сектор,
    прочее = прямой кик). Все mat_* абсолютные (у «dipole» доли mat_B отменены). Импорт .lte переносит
    мультиполи любого типа (поле от формы полюсов) и предупреждает о неучтённых параметрах."""
    ST = (1e-3, 2e-4, 0.5e-3, -1e-4, 0.0, 1e-3)

    def test_straight_body_kind_independent(self):
        """Одинаковое поле ⇒ бит-в-бит одинаковое тело у quad/sext/oct/kicker/dq и dipole с θ0=0."""
        f = dict(k1=0.8, k2=-35.0, k3=900.0)
        ref = make_magnet("quad", 0.2, **f).body(self.ST)
        for kind in ("sext", "oct", "kicker", "dq", "dipole", "monitor"):
            m = make_magnet(kind, 0.2, **f)
            self.assertTrue(m.kick_sector, kind)
            self.assertEqual(m.body(self.ST), ref, kind)

    def test_dipole_multipoles_absolute(self):
        """У «dipole» мультиполи абсолютные, как у «dq»: k2 = K2 (не K2·mat_B); тело = dq бит-в-бит."""
        g = dict(angle=0.0094691022946895, k2=-5.0, k3=40.0, edge1=-0.0127, edge2=0.0031)
        a, b = make_magnet("dipole", 0.4178302891796566, **g), make_magnet("dq", 0.4178302891796566, **g)
        self.assertTrue(a.sector)
        self.assertAlmostEqual(a.k2, -5.0, places=12); self.assertAlmostEqual(a.k3, 40.0, places=10)
        self.assertEqual(a.body(self.ST), b.body(self.ST))

    def test_kicker_static_dipole(self):
        """mat_B прямого элемента = дипольный кик: без мультиполей — «ящик» корректора, с ними — член θ0/L
        кик-тела. У опорной частицы Δpx = θ0 в обоих случаях."""
        box = make_magnet("kicker", 0.1, angle=1e-3)
        self.assertFalse(box.kick_sector)
        self.assertAlmostEqual(box.body((0.0,) * 6)[1], box.theta0, delta=1e-15)
        kick = make_magnet("kicker", 0.1, angle=1e-3, k2=30.0)
        self.assertTrue(kick.kick_sector)
        self.assertAlmostEqual(kick.body((0.0,) * 6)[1], kick.theta0, delta=1e-8)   # секступоль на смещении θL/2

    def test_symplectic_new_paths(self):
        for m in (make_magnet("kicker", 0.1, angle=1e-3, k1=0.5, k2=30.0), make_magnet("dq", 0.4, k1=-1.5),
                  make_magnet("sext", 0.166, k1=0.8, k2=-136.0),
                  make_magnet("dq", 1.0, angle=0.0247680586728496, k1=-1.96, k2=3.0, k3=150.0,
                              edge1=-0.01168, edge2=-0.01168)):
            M = jacobian4(m.body, TestSymplecticity.BASE)
            self.assertLess(float(np.max(np.abs(M.T @ J4 @ M - J4))), 1e-6, m.kind)

    def test_ring_matrix_gradient_any_kind(self):
        """Линейная оптика: K1 у секступоля/kicker/dq с θ0=0 — толстая линза (раньше дрейф) = якобиан тела."""
        import ring_matrix as rm
        for kind in ("sext", "kicker", "dq"):
            m = make_magnet(kind, 0.3, k1=0.4)     # слабый градиент: ошибка Йошиды-4 (N=5) ≪ допуска
            _, R = rm.element_matrix(m)
            F = jacobian4(m.body, (0.0,) * 6, h=1e-6)
            self.assertLess(float(np.max(np.abs(R[:4, :4] - F))), 1e-6, kind)

    def test_ring_map_follows_field(self):
        """Ряды ring_map у kicker с полем = тело ядра (раньше kicker в рядах был дрейфом)."""
        import ring_map as rmap
        mp = rmap.RingMap(None, order=4); da = mp.da
        m = make_magnet("kicker", 0.3, k1=1.2, k2=40.0, k3=500.0)
        ser = mp.element(da.identity(), m)
        z = (1e-3, 2e-4, 0.5e-3, -1e-4, 0.0, 0.0)
        val = [s[0] for s in da.compose(ser, [da.const(v) for v in z])]
        self.assertLess(max(abs(a - b) for a, b in zip(val[:4], m.body(z)[:4])), 1e-9)

    def test_importer_multipoles_any_type(self):
        import s2lattice_to_db as s2
        QP = s2.QP
        q = {"type": "KQUSE", "L": 0.3, "K1": 2.0, "K2": 15.0}
        self.assertEqual(s2.kind_of(q), "quad")
        self.assertAlmostEqual(s2.mat_of("quad", q, 0.3)["mat_S"] * QP, 15.0, places=12)
        b = {"type": "CSBEND", "L": 1.0, "ANGLE": 0.025, "K2": 3.0, "K3": 150.0}   # K2/K3 без K1
        self.assertEqual(s2.kind_of(b), "dipole")
        mb = s2.mat_of("dipole", b, 1.0)
        self.assertAlmostEqual(mb["mat_S"] * QP, 3.0, places=12)
        self.assertAlmostEqual(mb["mat_O"] * QP, 150.0, places=10)
        self.assertEqual(s2.kind_of({"type": "KSEXT", "L": 0.1, "K1": 0.5}), "sext")   # градиент — не слот корректора
        self.assertEqual(s2.kind_of({"type": "KSEXT", "L": 0.1}), "kicker")
        self.assertEqual(s2.active_of("quad", q), "k1|k2")

    def test_importer_fractional_multipoles(self):
        """F1…F3 при XREFERENCE (как elegant): профиль поля на радиусе XREF = 1 + F1 + F2 + F3."""
        import s2lattice_to_db as s2
        e = {"type": "CSBEND", "L": 1.0, "ANGLE": 0.02, "XREFERENCE": 0.01, "F1": 0.1, "F2": 0.02, "F3": 0.003}
        K = s2.lte_multipoles(e)
        rho, xr = 1.0 / 0.02, 0.01
        prof = sum(K[n] * rho * xr ** (n + 1) / math.factorial(n + 1) for n in range(3))
        self.assertAlmostEqual(prof, 0.1 + 0.02 + 0.003, places=12)
        self.assertEqual(s2.lte_multipoles({"type": "KQUAD", "L": 0.3, "F1": 0.1, "XREFERENCE": 0.01}),
                         (0.0, 0.0, 0.0))                     # дробная форма — только у изгиба

    def test_importer_reports_unused(self):
        import s2lattice_to_db as s2
        e = {"type": "CSBEND", "L": 1.0, "ANGLE": 0.02, "K4": 10.0, "TILT": 0.0, "DX": 1e-4, "N_SLICES": 5,
             "_str": ["SYSTEMATIC_MULTIPOLES"]}
        self.assertEqual(s2.unused_params(e), ["DX", "K4", "SYSTEMATIC_MULTIPOLES"])
        self.assertEqual(s2.unused_params({"type": "KQUAD", "L": 0.3, "F1": 0.1}), ["F1"])

    def test_importer_parse_comments_and_skipped(self):
        import s2lattice_to_db as s2, tempfile
        txt = ('Q1: KQUAD, L=0.3, K1=2.0 ! K2=99\nM1: MULT, L=0.0, KNL=0.1, ORDER=2\n'
               'S1: KSEXT, L=0.1, SYSTEMATIC_MULTIPOLES="m.sdds"\nring: LINE=(Q1,M1,S1)\n')
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.lte")
            with open(p, "w", encoding="utf-8") as f:
                f.write(txt)
            old = s2.LTE; s2.LTE = p
            try:
                elems, seq, skipped = s2.parse()
            finally:
                s2.LTE = old
        self.assertNotIn("K2", elems["Q1"])                  # «!» — комментарий
        self.assertEqual(skipped["M1"][0], "MULT")
        self.assertEqual(elems["S1"]["_str"], ["SYSTEMATIC_MULTIPOLES"])
        self.assertEqual(seq, ["Q1", "M1", "S1"])

    def test_flatten_kquse_and_lost(self):
        import lte_flatten as lf, tempfile
        txt = ('Q1: QUAD, L=0.3, K1=2.0, K2=15.0\nB1: SBEN, L=1.0, ANGLE=0.02, K3=150.0, K4=7.0, '
               'XREFERENCE=0.01, F2=0.01\nR: LINE=(Q1,B1)\n')
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "t.lte")
            with open(p, "w", encoding="utf-8") as f:
                f.write(txt)
            out, rep = lf.flatten(p, "R")
        self.assertIn("Q1: KQUSE", out)                       # у KQUAD нет K2 в elegant
        self.assertIn("K3=150", out); self.assertIn("F2=0.01", out); self.assertIn("XREFERENCE=0.01", out)
        self.assertEqual(rep["не перенесено параметров (типов эл.)"], {"K4": 1})


class _LinearOpticsPlant:
    """Синтетический plant оптики: observable = t0 + M·u + d (линейная модель с известным откликом M).

    Нужен ТОЛЬКО для юнит-теста контура: даёт детерминированную проверку IntegralSVD+Scheduler+ref без
    дорогого кольца. M — «истинная» матрица отклика (её же отдаёт build_response), u — вклад контура
    (state.actual по квад-актуаторам), d — синтетический увод Δk (имитация дрейфа градиентов)."""

    def __init__(self, M, t0, d, ids):
        self.M = np.asarray(M, float); self.t0 = np.asarray(t0, float)
        self.d = np.asarray(d, float); self.ids = list(ids)

    def _u(self, state):
        return np.array([state.actual[i] for i in self.ids])

    def _measure(self, state):
        state.observables["tune"] = self.t0 + self.M @ self._u(state) + self.d

    def schema(self):
        from control.state import Schema
        return Schema(self.ids)

    def history(self):
        from control.state import History
        return History(["opt"], moments=("x",))

    def reset(self, state):
        self._measure(state)

    def step(self, state, dt):
        self._measure(state)

    def build_response(self, actuators, observables):
        from control.controller import ResponseMatrix
        cols = [self.M[:, self.ids.index(a)] for a in actuators]
        R = np.column_stack(cols)
        return ResponseMatrix(R, [f"tune.{i}" for i in range(R.shape[0])], list(actuators))

    def fingerprint(self):
        return ("linopt", tuple(self.ids))


class TestOpticsLoopControl(unittest.TestCase):
    """Контур оптики (IntegralSVD + Scheduler + ref) на синтетическом Δk — быстрая, детерминированная.

    Проверяет: интегральный SVD-регулятор с уставкой (ref) сводит observable к ЦЕЛИ (не к нулю) при
    известной матрице отклика — тот же плагин, что орбитный FOFB, с векторным сенсором `tune`. Это ядро
    трёх новых контуров (квад/секступоль/октуполь); кольцевую правильность знака/масштаба отклика
    проверяет TestOpticsPlant, а сходимость на реальной машине — control/optics_loops_demo.py."""

    def _run(self, gain):
        from control.state import ExperimentState
        from control.scheduler import Scheduler
        from control.controllers.integral_svd import IntegralSVD
        ids = ["quad.A.dk", "quad.B.dk"]
        M = np.array([[3.0, 1.2], [-0.6, -2.4]])          # обусловленная 2×2 (разные знаки — как пара QF/QD)
        t0 = np.array([0.31, 0.22]); d = np.array([0.02, -0.015])
        plant = _LinearOpticsPlant(M, t0, d, ids)
        state = ExperimentState(plant.schema(), history=plant.history())
        state.refs["tune.target"] = t0.copy()             # уставка = проектные частоты (НЕ ноль)
        sch = Scheduler(plant, state, dt=1.0)
        sch.add(IntegralSVD(name="tune", sensor="tune", actuator_key="quad.", gain=gain,
                            rcond=1e-6, ref="tune.target"), rate_hz=1.0)
        e0 = float(np.linalg.norm(state.observables["tune"] - t0))
        errs = [e0]
        for _ in range(30):
            sch.step()
            errs.append(float(np.linalg.norm(state.observables["tune"] - t0)))
        return t0, state.observables["tune"].copy(), errs

    def test_converges_to_target(self):
        """‖tune − уставка‖ падает ниже 1e-9 (синтетический увод d полностью скомпенсирован)."""
        t0, final, errs = self._run(gain=0.6)
        self.assertGreater(errs[0], 1e-3, "в начале увод должен присутствовать")
        self.assertLess(errs[-1], 1e-9, "контур обязан свести увод к машинному нулю")
        self.assertTrue(all(errs[i + 1] <= errs[i] + 1e-12 for i in range(len(errs) - 1)),
                        "ошибка должна монотонно убывать (устойчивый интегральный контур)")

    def test_tracks_reference_not_zero(self):
        """Контур ведёт observable к УСТАВКЕ (ref), а не к нулю — это и даёт каскад контуров."""
        t0, final, errs = self._run(gain=0.6)
        np.testing.assert_allclose(final, t0, atol=1e-9)
        self.assertGreater(float(np.linalg.norm(final)), 0.3,
                           "цель далека от нуля — доказывает слежение за ref, а не зануление")


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск оптического plant на кольце")
class TestOpticsPlant(unittest.TestCase):
    """OpticsPlant на кольце: аналитическая оценка оптики на ВОЗМУЩЁННОЙ копии кольца реагирует на
    градиентные актуаторы (quad.<fam>.dk → Δk1). Решётка Сибири (№2, мелкая) — быстрый фикстур.

    Сторож регрессии ключевого решения Stage A: оптику считаем на _python_ring()-КОПИИ с внесённым
    Δk, а НЕ через apply_perturb (он доходит только до C++-трекера, аналитическое кольцо его не видит).
    Если кто-то «упростит» OpticsPlant обратно на apply_perturb/живое кольцо — отклик обнулится и тест
    упадёт. (Вертикальный νy и хроматизм на №2 вырождены — аналитическое кольцо Сибири не отдаёт k2 и
    вертикальный k1 как живые ручки; эти каналы проверяет optics_loops_demo на рабочей №26.)"""

    LAT, NSL = 2, 6

    @classmethod
    def setUpClass(cls):
        from control.ring_model import RingModel
        from control.optics_plant import OpticsPlant
        cls.rm = RingModel(lattice=cls.LAT, nsl=cls.NSL)
        cls.plant = OpticsPlant(model=cls.rm, observe=("tune",), order=4)

    @classmethod
    def tearDownClass(cls):
        cls.plant.close()

    def test_nominal_tune_fractional(self):
        """Проектные частоты чистого кольца — осмысленные дробные числа в (0, 1)."""
        nux, nuy = self.plant.nominal()["tune"]
        self.assertTrue(0.0 < nux < 1.0 and 0.0 < nuy < 1.0, f"νx={nux}, νy={nuy} вне (0,1)")

    def test_quad_tune_response_positive_nux(self):
        """∂νx/∂k1 > 0 для КАЖДОГО квадрупольного семейства (∝ +βx/4π, не зависит от решётки).

        Заодно отклик НЕ вырожден (|∂νx/∂k1| ~ O(1)) — значит градиентный актуатор реально доходит до
        аналитической оценки (внесён в копию кольца), а не теряется в apply_perturb."""
        R = self.plant.build_response(self.rm.dk_ids, ["tune"]).R
        self.assertEqual(R.shape, (2, len(self.rm.dk_ids)))
        nux_row = R[0]
        self.assertTrue(np.all(nux_row > 1e-4),
                        f"∂νx/∂k1 обязана быть > 0 у всех семейств (βx>0): {nux_row}")
        self.assertGreater(float(np.max(np.abs(nux_row))), 0.1,
                           "отклик νx на Δk1 не должен быть вырожденно мал (актуатор дошёл до оценки)")


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск модели полюсов ИД на решётке 26")
class TestUndulatorHardEdgePoles(unittest.TestCase):
    """Полюса ондуляторов — эквивалент «жёсткий край» (undulators/migrate_hardedge_poles.py): длина (8/π²)·λu/2,
    грани θ/2 того же знака, что угол (вертикальный ФОКУС краёв), отдельной линзы UFC нет; фокус-элемент старых
    БД (is_id_focus) — не ручка контура тюна."""

    def test_poles_are_hard_edge_without_lens(self):
        import math
        from control.ring_model import RingModel, is_id_focus, id_devices
        self.assertTrue(is_id_focus("und_UFC")); self.assertFalse(is_id_focus("und_QF1A"))
        rm = RingModel(lattice=26, nsl=2)
        try:
            mags = rm.ring.mags
            devs = id_devices(mags)
            if not devs:
                self.skipTest("в решётке 26 нет ондуляторов")
            self.assertEqual(len(devs), 38)
            self.assertFalse(any(is_id_focus(getattr(m, "dev", "")) for m in mags), "в решётке остался UFC")
            Lp = 8.0 / math.pi ** 2 * 0.009
            for i in devs[0][1]:
                m = mags[i]
                self.assertAlmostEqual(float(m.Leff), Lp, places=9)
                self.assertAlmostEqual(float(m.edge1), 0.5 * float(m.theta0), places=12)   # CSV: 15 знаков
                self.assertGreater(float(m.kEdge1), 0.0)                  # tg(E)·θ/L > 0 — вертикальный фокус
            self.assertIn("quad.und_QF1A.dk", rm.dk_ids)
        finally:
            rm.close()


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск отклика тюна на зазор ИД на решётке 26")
class TestUndulatorGapTune(unittest.TestCase):
    """Зазор планарного ондулятора двигает только ВЕРТИКАЛЬ: νy и ξy (горизонталь — дрейф при любом зазоре: грани θ/2
    полюсов гасят 1/ρ² при любом масштабе поля), а линейный член кокпита R_gap·(fac − 1) совпадает с полным
    пересчётом аналитической модели (ring_map) на тех же зазорах. Хроматизм — через бета-биение в секступолях."""

    def test_gap_moves_only_vertical_and_cockpit_term_matches(self):
        import types
        import ring_map as rmap
        from control.optics_plant import OpticsPlant
        from control.composite_plant import CompositePlant
        op = OpticsPlant(lattice=26, nsl=2, observe=("tune", "chroma"), order=2)
        try:
            if len(op._und) < 20:
                self.skipTest("в решётке 26 нет ондуляторов")

            def optics(gaps):
                rm = rmap.RingMap(op._perturbed({}, {}, {}, gaps), 2)
                ch = rm.chromaticity(rm.one_turn())
                return np.array([ch["nux"], ch["nuy"], ch["xi_x"][1], ch["xi_y"][1]])

            sub = [op._und[0][0], op._und[19][0]]                       # два ИД на 10 мм, остальные на опорном
            gaps = {g: 10.0 for g in sub}
            d_full = optics(gaps) - optics(None)
            self.assertLess(abs(d_full[0]), 1e-6, "зазор ИД сдвинул горизонтальный тюн")
            self.assertLess(d_full[1], -1e-4, "раскрытие ИД должно снижать νy (слабее верт. фокус)")
            cp = CompositePlant.__new__(CompositePlant)               # только линейный член кокпита, без трекинга
            cp.__dict__.update(optics=op, _gap_key=None, _gap_d={})
            cp._gap_ids, cp._gap_R = op.gap_response(only=sub)
            self.assertEqual(cp._gap_ids, sub)
            state = types.SimpleNamespace(undulators={g: types.SimpleNamespace(gap_mm=v) for g, v in gaps.items()})
            d = cp._gap_shift(state)
            self.assertLess(abs(float(d["tune"][0])), 1e-9)
            self.assertLess(abs(d["tune"][1] - d_full[1]), 0.01 * abs(d_full[1]))
            self.assertLess(abs(d["chroma"][0] - d_full[2]), 1e-6)
            self.assertLess(abs(d["chroma"][1] - d_full[3]), 0.02 * abs(d_full[3]))
        finally:
            op.close()


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск канала BSCALE на решётке 26")
class TestUndulatorGapTuneTracking(unittest.TestCase):
    """Живой зазор ИД в трекинге (TbtPlant): команда BSCALE масштабирует угол и грани полюсов (рамки заново) —
    сдвиг тюна C++-кольца = аналитике Σβy·(tgE₁+tgE₂)θ/L·(f²−1)/4π, горизонталь стоит; на опорном зазоре канал
    пуст; дешёвый путь set_bend_scale = полный apply_perturb бит-в-бит. DKV (только k1v) — общий канал ядра."""

    def test_bscale_moves_tracking_tune(self):
        import types
        import ring_cpp as rc
        import ring_analysis as ra
        import ring_response as rr
        import undulator_sr as sr
        from control.tbt_plant import TbtPlant
        from control.ring_model import id_devices
        r = rc.CppRing(lattice=26, nsl=8)
        try:
            devs = dict(id_devices(r.mags))
            if not devs:
                self.skipTest("в решётке 26 нет ондуляторов")
            tp = TbtPlant.__new__(TbtPlant)                            # только логика «зазор → масштаб полюсов»
            tp._und_poles = devs
            tp._und_sr = sr; tp._und_Kref = float(sr.K_of_B0(sr.B0_of_gap(sr.P.gap_ref)))

            def st(g):
                return types.SimpleNamespace(undulators={k: types.SimpleNamespace(gap_mm=g) for k in devs})

            self.assertEqual(tp._und_bscale(st(5.0)), [])             # опорный зазор → канал пуст
            bs = tp._und_bscale(st(10.0))
            self.assertEqual(len(bs), sum(len(v) for v in devs.values()))
            f = bs[0][1]
            py0 = r._python_ring()
            _, dy, _ = rr.tune_response(py0)
            expect = sum(-dy[i] * (py0.mags[i].kEdge1 + py0.mags[i].kEdge2) * (f * f - 1.0) for i, _ in bs)
            j = bs[0][0]; th0 = float(py0.mags[j].theta0)
            tw0 = ra.twiss(r)
            r.apply_perturb(bscales=bs)
            self.assertAlmostEqual(float(r._python_ring().mags[j].theta0), f * th0, places=15)
            tw1 = ra.twiss(r)
            self.assertLess(abs(tw1["nux"] - tw0["nux"]), 1e-5, "масштаб полюсов сдвинул горизонтальный тюн")
            self.assertLess(abs((tw1["nuy"] - tw0["nuy"]) - expect), 0.02 * abs(expect))
            self.assertLess(expect, -1e-3)
            p0 = [[1e-4, 0.0, 5e-5, 0.0, 0.0, 0.0]]
            ref = np.asarray(r.track(p0, 2), float)
            r.apply_perturb(); r.set_bend_scale([(i, 0.7) for i, _ in bs]); r.set_bend_scale(bs)   # без RESET
            self.assertTrue(np.array_equal(np.asarray(r.track(p0, 2), float), ref))
        finally:
            r.apply_perturb(); r.close()

    def test_dkv_vertical_only(self):
        """DKV меняет только вертикальный градиент k1v (описание следует за сервером), DK — обе плоскости."""
        import ring_cpp as rc
        r = rc.CppRing(lattice=26, nsl=2)
        try:
            j = next(i for i, m in enumerate(r.mags) if str(getattr(m, "kind", "")) == "quad")
            k1_0, k1v_0 = float(r.mags[j].k1), float(r.mags[j].k1v)
            r.apply_perturb(dkvs=[(j, 2e-3)])
            m = r._python_ring().mags[j]
            self.assertAlmostEqual(float(m.k1), k1_0, places=14)
            self.assertAlmostEqual(float(m.k1v), k1v_0 + 2e-3, places=14)
            r.apply_perturb(); r.add_perturb(dkvs=[(j, 1e-3)]); r.add_perturb(dkvs=[(j, 1e-3)])
            self.assertAlmostEqual(float(r._python_ring().mags[j].k1v), k1v_0 + 2e-3, places=14)
            with self.assertRaises(ValueError):
                r.add_perturb(kicks=[(j, 1e-6)])                      # кик — не аддитивный канал
            r.apply_perturb(dks=[(j, 1e-3)])
            m = r._python_ring().mags[j]
            self.assertAlmostEqual(float(m.k1), k1_0 + 1e-3, places=14)
            self.assertAlmostEqual(float(m.k1v), k1v_0 + 1e-3, places=14)
        finally:
            r.apply_perturb(); r.close()

    def test_per_device_slices_parity(self):
        """SS_NSL_DEV (срезы по типу устройства — полюса ИД на 1 срезе): C++ и Python-эталон режут одинаково
        (оборот совпадает до 1e-10), кусков меньше, тюн почти не меняется."""
        import ring_cpp as rc
        import ring_analysis as ra
        from control.tbt_plant import ID_POLE_DEVS
        poles = {d: 1 for d in ID_POLE_DEVS}
        c = rc.CppRing(lattice=26, nsl=8, nsl_dev=poles)
        c8 = rc.CppRing(lattice=26, nsl=8)
        saved = {k: os.environ.get(k) for k in ("SS_SECTOR_NSL", "SS_NSL_DEV")}
        try:
            if not any(str(getattr(m, "dev", "")) in poles for m in c8.mags):
                self.skipTest("в решётке 26 нет ондуляторов")
            self.assertLess(c.n_pieces, c8.n_pieces)
            os.environ["SS_SECTOR_NSL"] = "8"; os.environ["SS_NSL_DEV"] = c._nsl_dev_env()
            p0 = [[2e-4, 1e-5, -1e-4, 2e-6, 0.0, 5e-4]]
            d = np.max(np.abs(np.asarray(c.track(p0, 1), float) - np.asarray(ra.Ring(26).track(p0, 1), float)))
            self.assertLess(d, 1e-10, f"C++ ≠ Python при SS_NSL_DEV: {d}")
            os.environ.pop("SS_NSL_DEV")
            self.assertLess(abs(ra.twiss(c)["nux"] - ra.twiss(c8)["nux"]), 1e-6)
        finally:
            for k, v in saved.items():
                os.environ.pop(k, None) if v is None else os.environ.update({k: v})
            c.close(); c8.close()

    def test_set_kicks_equals_full_apply(self):
        """set_kicks без RESET (смена корректоров FOFB) даёт ТОТ ЖЕ трекинг, что полный apply_perturb, бит в бит;
        описание mags следует за ним (последний кик выигрывает)."""
        import ring_cpp as rc
        from orbit_correct import corrector_indices
        r = rc.CppRing(lattice=26, nsl=2)
        try:
            cj = corrector_indices(r)[:3]
            if len(cj) < 3:
                self.skipTest("нет корректоров")
            p0 = [[1e-5, 0.0, -1e-5, 0.0, 0.0, 0.0]]
            A = [(cj[0], 2e-5), (cj[1], -1e-5)]; B = [(cj[0], -3e-5), (cj[2], 1e-5)]
            r.apply_perturb(kicks=B, vkicks=[(cj[1], 4e-6)])
            ref = np.asarray(r.track(p0, 3), float)
            r.apply_perturb(kicks=A)                                     # было A, затем дешёво → B
            r.set_kicks(kicks=[(cj[0], -3e-5), (cj[1], 0.0), (cj[2], 1e-5)], vkicks=[(cj[1], 4e-6)])
            got = np.asarray(r.track(p0, 3), float)
            self.assertTrue(np.array_equal(got, ref), f"set_kicks ≠ apply_perturb: {np.max(np.abs(got - ref))}")
            self.assertEqual(float(r._python_ring().mags[cj[1]].hkick), 0.0)
            self.assertEqual(float(r._python_ring().mags[cj[1]].vkick), 4e-6)
        finally:
            r.apply_perturb(); r.close()

    def test_apply_cheap_paths(self):
        """Дешёвые пути TbtPlant._apply без RESET: движение челюстей — только изменившиеся полюса (set_bend_scale,
        вернувшийся к опорному = 1), смена корректоров — только изменившиеся кики (set_kicks, снятый кик = 0);
        смена «прочих» каналов (dks/…) — полный apply_perturb с текущими киками и масштабом полюсов."""
        from control.tbt_plant import TbtPlant
        calls = []

        class _Ring:
            def apply_perturb(self, **kw): calls.append(("apply", kw))
            def set_bend_scale(self, pairs): calls.append(("bscale", list(pairs)))
            def set_kicks(self, **kw): calls.append(("set", kw))

        tp = TbtPlant.__new__(TbtPlant)
        tp.ring = _Ring(); tp._last_kick = None; tp._apply_rf = lambda s: None
        cfg = {"kk": [(1, 0.0), (2, 5e-7)], "vv": [(1, 0.0)], "bs": [(7, 0.9)], "x": {}}
        tp.kick_lists = lambda s: (list(cfg["kk"]), list(cfg["vv"]), dict(cfg["x"], bscales=list(cfg["bs"])))
        tp._apply(None)
        self.assertEqual(calls[-1], ("apply", dict(kicks=[(1, 0.0), (2, 5e-7)], vkicks=[(1, 0.0)], bscales=[(7, 0.9)])))
        n = len(calls); tp._apply(None); self.assertEqual(len(calls), n)          # ничего не изменилось
        cfg["bs"] = [(9, 0.8)]; tp._apply(None)                                  # одни едут, другие вернулись
        self.assertEqual(calls[-1], ("bscale", [(7, 1.0), (9, 0.8)]))
        cfg["kk"] = [(1, 1e-6)]; cfg["vv"] = [(1, 2e-6)]; tp._apply(None)       # FOFB сменил корректоры
        self.assertEqual(calls[-1], ("set", dict(kicks=[(1, 1e-6), (2, 0.0)], vkicks=[(1, 2e-6)])))
        cfg["x"] = {"dks": [(5, 1e-4)]}; tp._apply(None)                         # другой канал → полный путь
        self.assertEqual(calls[-1][0], "apply")
        self.assertEqual(calls[-1][1]["bscales"], [(9, 0.8)])
        self.assertEqual(calls[-1][1]["kicks"], [(1, 1e-6)])


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "fb_experiments.csv")),
                     "нет fb_experiments.csv — пропуск конфига режима 4")
class TestMode4Config(unittest.TestCase):
    """Режим 4: загрузка конфига эксперимента ОС из БД по ключу + соблюдение подписки рекордерами."""

    def test_config_roundtrip(self):
        """load(1) собирает мастер+сателлиты: решётка, контуры, возмущения, подписка; plant по умолчанию tbt."""
        from control import fb_config
        cfg = fb_config.load(1)
        self.assertEqual(cfg.lattice_id, 26)
        self.assertGreater(cfg.steps, 0)
        self.assertEqual(cfg.plant, "tbt")
        self.assertEqual({l.name for l in cfg.enabled_loops()}, {"fofb_x", "fofb_y"})
        self.assertIn("ground_vibration", {c["name"] for c in cfg.disturbance_spec})
        self.assertIn("screen", {r.signal for r in cfg.enabled_records()})

    def test_optics_experiment_config(self):
        """load(2) = оптический эксперимент: plant=optics, дрейф градиентов, контуры tune/chroma."""
        from control import fb_config
        cfg = fb_config.load(2)
        self.assertEqual(cfg.plant, "optics")
        self.assertGreater(cfg.drift[0], 0.0)                 # drift_k1 прочитан
        self.assertEqual(len(cfg.drift), 3)
        self.assertEqual({l.name for l in cfg.enabled_loops()}, {"tune", "chroma"})   # adts выключен
        self.assertEqual({l.actuator_key for l in cfg.enabled_loops()}, {"quad.", "sext."})

    def test_vector_recorder(self):
        """TuneRecorder пишет компоненты векторной наблюдаемой tune.x/tune.y из state.observables."""
        import types
        from control import recorders
        rec = recorders.TuneRecorder(params="", cadence="per_window")
        st = types.SimpleNamespace(observables={"tune": np.array([0.31, 0.22])})
        rec.sample(st, 0, None); rec.sample(st, 1, None)
        rows = rec.result_rows()
        self.assertEqual(len(rows), 4)                         # 2 компоненты × 2 окна
        d = {(s, w): v for s, w, v in rows}
        self.assertAlmostEqual(d[("tune.x", 0)], 0.31)
        self.assertAlmostEqual(d[("tune.y", 1)], 0.22)

    def test_subscription_honored(self):
        """recorders.build берёт только enabled; отключённый screen НЕ инстанцируется, неизвестный сигнал пропущен."""
        from control import fb_config, recorders
        R = fb_config.Record
        recs = recorders.build([R("orbit_rms", "", "per_window", True),
                                R("screen", "SCR_U19", "final", False),   # отключён
                                R("sf", "channel=b1", "post", True)])      # нереализован в v1
        names = [type(r).__name__ for r in recs]
        self.assertEqual(names, ["OrbitRmsRecorder"], f"подписка нарушена: {names}")
        self.assertNotIn("ScreenRecorder", names, "отключённый screen не должен считаться")

    def test_rf_experiment_config(self):
        """load(3/4) = длина кольца плывёт: возмущение circumference; 4 — + контур ВЧ-частоты rf.delta→rf.freq."""
        from control import fb_config
        c3, c4 = fb_config.load(3), fb_config.load(4)
        for c in (c3, c4):
            self.assertTrue(c.relative_orbit)
            self.assertEqual([d["quantity"] for d in c.disturbance_spec], ["circumference"])
            self.assertIn("rf", {r.signal for r in c.enabled_records()})
        self.assertEqual({l.name for l in c3.enabled_loops()}, {"fofb_x", "fofb_y"})
        rf = [l for l in c4.enabled_loops() if l.name == "rf"]
        self.assertEqual(len(rf), 1)
        self.assertEqual((rf[0].sensor, rf[0].actuator_key), ("rf.delta", "rf.freq"))
        import control.mode4_run as m4
        self.assertEqual(type(m4.make_loop(rf[0])).__name__, "RfFrequencyLoop")   # класс — по датчику
        self.assertFalse(fb_config.load(1).relative_orbit)                        # прежние опыты — как были

    def test_floor_rf_pair_config(self):
        """load(5/6) = пол горизонтали: BPM абсолютно, но пучок на старте центрирован на CO; 6 — + контур ВЧ.
        Прежние опыты center_bunch не задают (None → как relative_orbit, поведение не меняется)."""
        from control import fb_config
        c5, c6 = fb_config.load(5), fb_config.load(6)
        for c in (c5, c6):
            self.assertFalse(c.relative_orbit); self.assertIs(c.center_bunch, True)
        self.assertNotIn("rf", {l.name for l in c5.enabled_loops()})
        self.assertIn("rf", {l.name for l in c6.enabled_loops()})
        self.assertIsNone(fb_config.load(1).center_bunch)
        self.assertIsNone(fb_config.load(3).center_bunch)

    def test_orbit_dc_recorder(self):
        """OrbitRmsRecorder: orbit_dc — СКО орбиты, усреднённой по второй половине окон (колебание гасится)."""
        import types
        from control import recorders
        rec = recorders.OrbitRmsRecorder()
        base = np.array([1e-6, -1e-6, 2e-6, -2e-6])
        osc = np.array([5e-6, 5e-6, -5e-6, -5e-6])
        for w in range(8):                                     # статика + знакопеременное колебание
            o = base + (osc if w % 2 else -osc)
            rec.sample(types.SimpleNamespace(observables={"ebpm.x": o, "ebpm.y": 0.5 * o}), w, None)
        rec.finalize(None, None)
        d = {s: v for s, w, v in rec.result_rows() if s.startswith("orbit_dc")}
        self.assertAlmostEqual(d["orbit_dc.x"], float(np.sqrt(np.mean(base ** 2))) * 1e6, places=9)
        self.assertAlmostEqual(d["orbit_dc.y"], 0.5 * float(np.sqrt(np.mean(base ** 2))) * 1e6, places=9)


class TestRfFrequencyLoop(unittest.TestCase):
    """Контур ВЧ-частоты (control/controllers/rf_frequency.py) и возмущение длины кольца — без 6D-кольца."""

    def test_circumference_disturbance(self):
        import types
        import disturbance as dist
        m = types.SimpleNamespace(Leff=0.1, theta0=0.0, k1=0.0, k2=0.0, k3=0.0)
        ring = types.SimpleNamespace(n_mag=3, kinds=["quad"] * 3, s_out=[1.0, 2.0, 3.0], mags=[m, m, m])
        spec = [dict(name="tide", target="ring", plane="", form="tone", correlation="common",
                     quantity="circumference", amplitude=2e-5, freq_hz=7.5, phase_deg=0.0)]
        d = dist.Disturbance(spec, ring, seed=1)
        t = 0.01
        self.assertAlmostEqual(d.circumference(t, 1e-4), 2e-5 * np.sin(2 * np.pi * 7.5 * t), places=15)
        dx, dy = d.offsets(t, 1e-4)
        self.assertEqual(float(np.abs(dx).max() + np.abs(dy).max()), 0.0)        # не трогает магниты
        self.assertEqual(d.field(t, 1e-4), {})
        self.assertEqual(dist.Disturbance([], ring).circumference(t, 1e-4), 0.0)

    def test_delta_estimate_projection(self):
        from control.controllers.rf_frequency import delta_estimate
        eta = np.array([0.0, 0.02, 0.01, -0.005, 0.0])
        self.assertAlmostEqual(delta_estimate(eta * 3e-4, eta), 3e-4, places=15)  # чистая дисперсионная орбита
        bump = np.array([1e-6, 0.0, 0.0, 0.0, -1e-6])                            # орбита ⊥ η — энергию не меняет
        self.assertAlmostEqual(delta_estimate(eta * 1e-4 + bump, eta), 1e-4, places=15)

    def test_loop_converges_on_static_plant(self):
        """Интегратор с бокс-фильтром на статическом plant δ̂ = d + R·Δf сходится к Δf = −d/R."""
        import types
        from control import Schema, ResponseMatrix
        from control.controllers.rf_frequency import RfFrequencyLoop
        R = -5.6e-5                                                              # dδ/dΔf, 1/Гц (СИЛА-26)
        res = types.SimpleNamespace(response=lambda a, o: ResponseMatrix([[R]], ["rf.delta.0"], a))
        lp = RfFrequencyLoop(gain=0.5)
        lp.setup(res, Schema(["corr.000.hkick", "rf.freq"])); lp.reset()
        self.assertEqual(lp.act, ["rf.freq"])
        self.assertIn("rf.delta", lp.provides().observables)
        d = -3.5e-4                                                              # ΔC = 20 мкм ⇒ δ ≈ −3.5e-4
        u = 0.0
        for _ in range(80):
            view = types.SimpleNamespace(read=lambda s: np.array([d + R * u]), ref=lambda n: None)
            u = lp.update(0.0, view)["rf.freq"]
        self.assertAlmostEqual(u, -d / R, delta=1e-3 * abs(d / R))               # ≈ −6.25 Гц


class TestMode4Store(unittest.TestCase):
    """Режим 4: писатель результатов — идемпотентность по ключу (перезапись, не дублирование) во ВРЕМЕННОЙ БД."""

    def setUp(self):
        import tempfile
        from control import fb_store
        self.fb_store = fb_store
        self._tb, self._res = fb_store.TB, fb_store.RES
        self.tmp = tempfile.mkdtemp()
        fb_store.TB = self.tmp
        fb_store.RES = os.path.join(self.tmp, "res")

    def tearDown(self):
        import shutil
        self.fb_store.TB, self.fb_store.RES = self._tb, self._res
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_results_keyed_rewrite_is_idempotent(self):
        """Двойная запись одного fb_exp_id не плодит строк; чужой эксперимент сохраняется."""
        rowsA = [("orbit_rms.x", 0, 1.0), ("orbit_rms.x", 1, 2.0)]
        self.fb_store.write_results(1, rowsA)
        self.fb_store.write_results(1, rowsA)                  # повтор
        self.fb_store.write_results(2, [("light.x", 0, 9.0)])  # чужой эксперимент
        r1 = self.fb_store.read_results(1); r2 = self.fb_store.read_results(2)
        self.assertEqual(len(r1["orbit_rms.x"][1]), 2, "повторная запись не должна дублировать")
        self.assertEqual(r2["light.x"][1], [9.0], "чужой эксперимент обязан сохраниться")
        self.fb_store.write_results(1, rowsA)                  # ещё раз — exp2 не трогаем
        self.assertEqual(self.fb_store.read_results(2)["light.x"][1], [9.0])

    def test_artifact_roundtrip(self):
        """write_artifact → .npz + манифест; read_artifacts читает массив; повтор не дублирует строку."""
        self.fb_store.write_artifact(1, "screen", {"ideal": np.ones((3, 3))}, meta="t")
        self.fb_store.write_artifact(1, "screen", {"ideal": np.zeros((3, 3))}, meta="t")  # перезапись
        got = self.fb_store.read_artifacts(1)
        self.assertIn("screen", got)
        np.testing.assert_array_equal(got["screen"]["ideal"], np.zeros((3, 3)))
        import csv
        with open(os.path.join(self.tmp, "fb_artifacts.csv"), encoding="utf-8-sig") as f:
            rows = [r for r in csv.reader(f, delimiter=";")][1:]
        self.assertEqual(len([r for r in rows if r and r[0] == "1"]), 1, "манифест не должен дублировать")


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "..", "configs", "bpm_map_26.json")),
                     "нет bpm_map_26.json — пропуск карты имён BPM")
class TestBpmMap(unittest.TestCase):
    """Карта имён BPM ↔ индексов: имя SP<nn>.BPM<mm> ↔ номер в векторе орбиты, round-trip, sensor_id."""

    def test_name_index_roundtrip(self):
        from control.bpm_map import BpmMap
        m = BpmMap.load(26)
        self.assertEqual(len(m), 400)
        self.assertEqual(len(set(m.names)), len(m), "имена BPM обязаны быть уникальны")
        for name in ("SP01.BPM01", "SP01.BPM06", "SP02.BPM01", "SP10.BPM40"):
            self.assertEqual(m.name(m.index(name)), name, f"round-trip имя↔индекс для {name}")
        self.assertEqual(m.index(5), 5, "int-индекс пропускается как есть")

    def test_sensor_id_and_named(self):
        """Доступ по ground-truth sensor_id и named(vector) → {имя: значение}."""
        from control.bpm_map import BpmMap
        m = BpmMap.load(26)
        idx = m.index("SP01.BPM06")
        self.assertEqual(m.index(m.entries[idx]["sensor_id"]), idx, "sensor_id → тот же индекс, что имя")
        d = m.named(np.arange(len(m), dtype=float))
        self.assertEqual(d["SP01.BPM06"], float(idx), "named() сопоставляет имя значению по его индексу")


class TestLiveStructures(unittest.TestCase):
    """Живые структуры доски для бесконечного режима: NoiseModel (закон шума) + SensorConfig + регионы."""

    def test_noise_model_revision(self):
        """set/touch бампают rev — plant по нему подхватывает смену закона на лету."""
        from control import NoiseModel
        nm = NoiseModel([{"name": "a"}])
        self.assertEqual(nm.rev, 0)
        nm.touch()
        self.assertEqual(nm.rev, 1)
        nm.set([{"name": "b"}, {"name": "c"}])
        self.assertEqual(nm.rev, 2)
        self.assertEqual(len(nm.spec), 2)

    def test_state_has_live_regions(self):
        """ExperimentState даёт позиции sensors (настройки сенсоров) и noise (закон шума)."""
        from control import ExperimentState, Schema, SensorConfig
        st = ExperimentState(Schema(["a"]))
        self.assertEqual(st.sensors, {})
        self.assertIsNone(st.noise)
        st.sensors["ebpm"] = SensorConfig(meas_noise=1e-4, enabled=False, setpoint=0.5)
        self.assertEqual(st.sensors["ebpm"].meas_noise, 1e-4)
        self.assertFalse(st.sensors["ebpm"].enabled)
        self.assertEqual(st.sensors["ebpm"].setpoint, 0.5)

    def test_per_corrector_drive_speed_is_live(self):
        """Карта корректоров + индивидуальный slew-лимит привода: set по pid ограничивает ТОЛЬКО свой
        корректор, apply читает slew каждый такт (→ правка скорости привода подхватывается живьём)."""
        import types
        from control.bpm_map import corrector_entries, BpmMap
        from control import Limits, Schema
        n = 6
        ring = types.SimpleNamespace(
            s_out=list(range(n)),
            mags=[types.SimpleNamespace(center=[float(i), 0.0, 0.0], dev=f"COR{i}", sensor_id=str(100 + i))
                  for i in range(n)])
        cm = BpmMap(corrector_entries(ring, list(range(n)), n_superperiods=2))
        self.assertEqual(len(cm), n)
        self.assertEqual(cm.name(0), "SP01.COR01")
        self.assertEqual(cm.name(3), "SP02.COR01")            # per = n/nsp = 3
        self.assertEqual(cm.index("101"), 1)                   # адресация по sensor_id
        i = cm.index("SP02.COR01")                             # конкретный корректор по имени
        sch = Schema([f"corr.{j:03d}.hkick" for j in range(n)])
        lim = Limits(sch); dt = 1e-3
        lim.set(f"corr.{i:03d}.hkick", slew=1e-4)              # медленный привод этого корректора
        a = lim.apply(np.zeros(n), np.ones(n), dt)
        self.assertAlmostEqual(a[i], 1e-4 * dt)                # ограничен slew*dt
        self.assertEqual(a[0], 1.0)                            # сосед без лимита — мгновенно
        lim.set(f"corr.{i:03d}.hkick", slew=1e-2)              # LIVE: ускорили привод
        a2 = lim.apply(a, np.ones(n), dt)
        self.assertAlmostEqual(a2[i] - a[i], 1e-2 * dt)        # новый кап подхвачен тем же apply


class TestFacilityState(unittest.TestCase):
    """Расширенная структура установки: BeamSnapshot (пучок по линии, запрос по имени магнита) и
    MachineView (lazy положение/поле магнита). Быстрые тесты на фейковых кольце/plant'е — без трекинга."""

    def test_magnet_entries_unique_and_roundtrip(self):
        import types
        from control.bpm_map import magnet_entries, BpmMap
        n = 8
        ring = types.SimpleNamespace(
            circumference=800.0, n_mag=n, s_out=[100.0 * i for i in range(n)],
            mags=[types.SimpleNamespace(center=[float(i), 0.0, 0.0], dev=f"Q{i % 3}", kind="quad")
                  for i in range(n)])          # dev повторяются (Q0/Q1/Q2) → нужен ординал
        mm = BpmMap(magnet_entries(ring, n_superperiods=2))
        self.assertEqual(len(mm), n)
        self.assertEqual(len(set(mm.names)), n, "имена магнитов обязаны быть уникальны (ординал при повторе dev)")
        for i in range(n):
            self.assertEqual(mm.index(mm.name(i)), i)         # round-trip имя↔индекс

    def test_beam_snapshot_query_by_magnet(self):
        from collections import deque
        from control import BeamSnapshot
        from control.bpm_map import BpmMap
        mm = BpmMap([{"index": 0, "name": "m0"}, {"index": 1, "name": "m1"}, {"index": 2, "name": "m2"}])
        s1 = np.arange(4 * 6, dtype=float).reshape(4, 6)      # 4 границы (3 магнита + вход)
        s2 = s1 + 10.0
        bs = BeamSnapshot(deque([s1, s2], maxlen=9), mm)
        np.testing.assert_array_equal(bs.snap, s2)            # snap = последний
        np.testing.assert_allclose(bs.avg, (s1 + s2) / 2)     # avg = среднее по окну
        self.assertEqual(bs.plane(1, "x"), s2[1][0])
        np.testing.assert_array_equal(bs.at_magnet("m1", "entry"), s2[1])   # вход m1 = строка 1
        np.testing.assert_array_equal(bs.at_magnet("m1", "exit"), s2[2])    # выход m1 = строка 2

    def test_machine_view_position_and_field(self):
        import types
        from control.bpm_map import BpmMap
        from control.tbt_plant import MachineView
        mm = BpmMap([{"index": 0, "name": "d0"}, {"index": 1, "name": "q1"}])
        ring = types.SimpleNamespace(mags=[
            types.SimpleNamespace(center=[1.0, 0.0, 0.0], theta0=0.01, k1=0.0, k2=0.0, k3=0.0),
            types.SimpleNamespace(center=[2.0, 0.0, 0.0], theta0=0.0, k1=1.5, k2=0.0, k3=0.0)])
        fake = types.SimpleNamespace(magnets=mm, ring=ring,
                                     _applied={"dx": np.array([0.0, 3e-6]), "dy": np.array([0.0, -2e-6]),
                                               "field": {"dks": np.array([0.0, 0.1])}})
        mv = MachineView(fake)
        px, py, pz = mv.position("q1")                        # центр + сдвиг возмущения
        self.assertAlmostEqual(px, 2.0 + 3e-6); self.assertAlmostEqual(py, -2e-6)
        self.assertAlmostEqual(mv.field("q1")["k1"], 1.5 + 0.1)   # номинал k1 + Δ из dks (рябь)
        self.assertAlmostEqual(mv.field("d0")["theta0"], 0.01)    # без применённой Δ — номинал


class TestImperativeLoop(unittest.TestCase):
    """Императивный контур (wants_ctl): читает из доски, пишет в СВОЙ слот через ctl — без возврата словаря.
    Планировщик обнуляет слот перед update и суммирует слоты после → ctl.set доходит до state.actual."""

    class _FakePlant:
        def __init__(self):
            from control import Schema
            self._schema = Schema(["corr.000.hkick", "corr.001.hkick"])

        def schema(self):
            return self._schema

        def history(self):
            return None

        def reset(self, state):
            state.observables["ebpm.x"] = np.array([2.0, -1.0])

        def step(self, state, dt):
            pass

        def fingerprint(self):
            return ("fake",)

    def test_ctl_set_and_ref_reach_state(self):
        from control import ExperimentState, Scheduler, Controller

        class Imp(Controller):
            wants_ctl = True

            def update(self, t, view, ctl):
                m = view.read("ebpm.x")
                ctl.set("corr.000.hkick", -0.5 * float(m[0]))     # пишем СВОЙ слот
                ctl.ref("golden.x", np.array([7.0]))              # уставка другому контуру

        plant = self._FakePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0)                     # reset() ставит observables
        sch.add(Imp(name="imp"), rate_hz=1.0)
        sch.step()
        self.assertAlmostEqual(state.actual["corr.000.hkick"], -0.5 * 2.0)   # ctl.set → слот → actual
        self.assertEqual(state.actual["corr.001.hkick"], 0.0)               # нетронутый актуатор = 0
        self.assertEqual(float(state.refs["golden.x"][0]), 7.0)             # ctl.ref записал уставку

    def test_return_dict_path_still_works(self):
        """Обратная совместимость: контур БЕЗ wants_ctl по-прежнему возвращает словарь."""
        from control import ExperimentState, Scheduler, Controller

        class Ret(Controller):
            def update(self, t, view):
                return {"corr.001.hkick": 0.3}

        plant = self._FakePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0)
        sch.add(Ret(name="ret"), rate_hz=1.0)
        sch.step()
        self.assertAlmostEqual(state.actual["corr.001.hkick"], 0.3)


class TestBumplessEnable(unittest.TestCase):
    """Безударное включение контура: после «отпустить» интегратор стартует с нуля (без броска сверх нужной
    поправки), после «заморозить» — продолжает с удержанного вклада (без скачка)."""

    class _Plant:                                   # obs = d + R·Δk: контур гонит obs к нулю → Δk = −d/R
        d, R = 1.0, 2.0

        def __init__(self):
            from control import Schema
            self._schema = Schema(["q.dk"])

        def schema(self):
            return self._schema

        def history(self):
            return None

        def reset(self, state):
            self.step(state, 0.0)

        def step(self, state, dt):
            state.observables["tune"] = np.array([self.d + self.R * state.actual["q.dk"]])

        def fingerprint(self):
            return ("bumpless",)

    def _run(self):
        from control import ExperimentState, Scheduler
        from control.controllers.integral_svd import IntegralSVD

        class Loop(IntegralSVD):
            def setup(self, resources, schema):
                self.act = ["q.dk"]; self.writes = self.act
                self.P = np.array([[1.0 / TestBumplessEnable._Plant.R]]); self.u = np.zeros(1)

        plant = self._Plant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0)
        sch.add(Loop(name="tune", sensor="tune", actuator_key="q.", gain=0.7), rate_hz=1.0, enabled=False)
        return sch, state

    def test_release_then_enable_has_no_overshoot(self):
        sch, state = self._run()
        sch.set_enabled("tune", True)
        for _ in range(40):
            sch.step()
        target = -self._Plant.d / self._Plant.R
        self.assertAlmostEqual(state.actual["q.dk"], target, places=9)
        sch.set_enabled("tune", False)                 # «отпустить» (дефолт): вклад обнулён
        sch.step()
        self.assertEqual(state.actual["q.dk"], 0.0)
        sch.set_enabled("tune", True)                  # безударно: с нуля, как при первом включении
        seq = []
        for _ in range(20):
            sch.step(); seq.append(state.actual["q.dk"])
        self.assertAlmostEqual(seq[0], 0.7 * target, places=9)       # без старого интеграла было бы 1,7·target
        self.assertTrue(all(abs(v) <= abs(target) + 1e-12 for v in seq))   # нет перерегулирования
        self.assertTrue(all(abs(b) >= abs(a) - 1e-12 for a, b in zip(seq, seq[1:])))   # монотонно к цели

    def test_freeze_then_enable_keeps_contribution(self):
        sch, state = self._run()
        sch.set_enabled("tune", True)
        for _ in range(40):
            sch.step()
        held = state.actual["q.dk"]
        sch.set_enabled("tune", False, mode="freeze")  # «заморозить»: вклад держится
        for _ in range(3):
            sch.step()
        self.assertAlmostEqual(state.actual["q.dk"], held, places=12)
        sch.set_enabled("tune", True)
        sch.step()
        self.assertAlmostEqual(state.actual["q.dk"], held, places=9)  # продолжает с удержанного — без скачка

    def test_reenable_before_release_takes_effect_is_continuous(self):
        """Выкл и сразу вкл (без такта между ними): слот ещё не обнулён — интегратор подхватывает его, скачка нет."""
        sch, state = self._run()
        sch.set_enabled("tune", True)
        for _ in range(40):
            sch.step()
        held = state.actual["q.dk"]
        sch.set_enabled("tune", False); sch.set_enabled("tune", True)
        sch.step()
        self.assertAlmostEqual(state.actual["q.dk"], held, places=9)


class TestCompositeLinearInline(unittest.TestCase):
    """CompositePlant (линеаризованная оптика) БЕЗ кольца: объединённая схема, inline-модель
    obs=nominal+drift_off+R·u, диспетчер build_response (ebpm→tbt, оптика→кэш). Логику проверяем на
    инъекции модели (минуя дорогой __init__); сквозной прогон на реальном кольце — opt-in smoke."""

    def _plant(self):
        from control.composite_plant import CompositePlant
        from control.state import Schema
        from control.controller import ResponseMatrix
        cp = object.__new__(CompositePlant)                   # минуем __init__ (без ring_cpp)
        cp._opt_ids = ["quad.F.dk", "quad.G.dk", "sext.S.dk2", "oct.O.dk3"]
        cp._fam = {"tune": ["quad.F.dk", "quad.G.dk"], "chroma": ["sext.S.dk2"], "adts": ["oct.O.dk3"]}
        cp._nominal = {"tune": np.array([0.3, 0.4]), "chroma": np.array([5.0, 6.0]), "adts": np.array([1., 2., 3.])}
        cp._drift_off = {"tune": np.array([0.01, -0.02]), "chroma": np.array([0.5, -0.3]), "adts": np.array([0.1, 0.1, 0.1])}
        cp._R = {"tune": np.array([[2.0, 0.0], [0.0, 3.0]]), "chroma": np.array([[1.0], [1.0]]),
                 "adts": np.array([[1.0], [1.0], [1.0]])}

        class _Tbt:
            def schema(self): return Schema(["corr.000.hkick"])
            def history(self): return None
            def build_response(self, act, obs): return ResponseMatrix(np.eye(len(act)), ["ebpm"], list(act))
        cp.tbt = _Tbt()
        return cp

    def test_union_schema_and_inline_and_dispatch(self):
        from control.state import ExperimentState
        cp = self._plant()
        self.assertEqual(cp.schema().ids, ["corr.000.hkick"] + cp._opt_ids)   # ОБЪЕДИНЕНИЕ актуаторов
        state = ExperimentState(cp.schema())
        state.actual["quad.F.dk"] = 0.1; state.actual["quad.G.dk"] = -0.05
        cp._measure_inline(state)                              # obs = nominal + drift_off + R·u
        exp = cp._nominal["tune"] + cp._drift_off["tune"] + cp._R["tune"] @ np.array([0.1, -0.05])
        np.testing.assert_allclose(state.observables["tune"], exp)
        # диспетчер build_response: оптика → кэш-матрица; орбита → tbt
        self.assertEqual(cp.build_response(["quad.F.dk", "quad.G.dk"], ["tune"]).R.shape, (2, 2))
        self.assertEqual(cp.build_response(["corr.000.hkick"], ["ebpm.x"]).rows, ["ebpm"])
        with self.assertRaises(ValueError):                   # неизвестная оптическая наблюдаемая
            cp.build_response(["x"], ["unknown"])


class TestObserver(unittest.TestCase):
    """Слой Observer: подписка, кастомный экстрактор, общий промежуточный считается один раз (кэш ctx)."""

    def test_subscription_and_custom_extractor(self):
        from control import Observer, ObsSpec, ObsContext, ExperimentState, Schema, OBSERVABLES
        reg = {"a": ObsSpec(("raw_a",), lambda ctx: ctx.get("raw_a")),
               "b": ObsSpec(("raw_a",), lambda ctx: ctx.get("raw_a") * 2),   # кастомная производная
               "c": ObsSpec(("raw_a",), lambda ctx: ctx.get("raw_a") * 99)}  # не подписан
        obs = Observer(registry=reg, subscribe=["a", "b"])
        st = ExperimentState(Schema(["x"]))
        ctx = ObsContext(st, None); ctx.set("raw_a", np.array([3.0]))
        obs.update(st, ctx)
        self.assertEqual(sorted(st.observables), ["a", "b"])              # подписка соблюдена (c не считается)
        np.testing.assert_array_equal(st.observables["b"], np.array([6.0]))  # кастомный экстрактор работает
        # МОДУЛЬНОСТЬ: оптика ПЕРЕЕХАЛА из observer.py в файлы-контуры (несут вход вычисления через provides())
        self.assertNotIn("tune", OBSERVABLES)                             # observer.py больше НЕ знает про tune
        from control.obs_assembly import observer_for
        from control.controllers.optics_tune import TuneLoop
        o = observer_for([TuneLoop()], serves={"tune"})
        self.assertIn("tune", o.registry)                                # контур принёс вычисление своей наблюдаемой
        self.assertIn("chroma_raw", o.intermediates)                     # и свой промежуточный (chroma_raw)

    def test_intermediate_built_once(self):
        """Общий промежуточный строится ОДИН раз и дальше берётся из кэша ObsContext."""
        from control import ObsContext
        import control.observer as obsmod
        calls = {"n": 0}

        def _counting(ctx):
            calls["n"] += 1
            return 42

        obsmod._INTERMEDIATES["_test_shared"] = _counting
        try:
            ctx = ObsContext(None, None)
            self.assertEqual(ctx.get("_test_shared"), 42)
            self.assertEqual(ctx.get("_test_shared"), 42)
            self.assertEqual(calls["n"], 1)                               # построен единожды, остальное — кэш
        finally:
            obsmod._INTERMEDIATES.pop("_test_shared", None)


class TestArbiter(unittest.TestCase):
    """Слой Arbiter: дефолт = сумма по слотам (бит-в-бит), кастомный combine, конвейер stages, выгрузка sink,
    и что планировщик реально зовёт combine (команда) и publish (выгрузка) на каждом такте."""

    def test_default_is_sum(self):
        from control import Arbiter, ExperimentState, Schema
        st = ExperimentState(Schema(["a", "b"]))
        st.config.v[:] = [1.0, 2.0]
        st.slot("l1").v[:] = [0.1, 0.0]
        st.slot("l2").v[:] = [0.0, 0.3]
        req = Arbiter().combine(st)
        np.testing.assert_allclose(req, st.effective_target())        # дефолт == прежнее поведение
        np.testing.assert_allclose(req, [1.1, 2.3])

    def test_custom_combine_and_stage(self):
        from control import Arbiter, ExperimentState, Schema
        st = ExperimentState(Schema(["a", "b"]))
        st.config.v[:] = [5.0, 5.0]
        st.slot("l1").v[:] = [1.0, 1.0]
        arb = Arbiter(combine=lambda s: s.config.v.copy(),            # свой закон: игнорирует слоты
                      stages=[lambda s, req: req + 1.0])              # конвейер: +1 (место модели токов)
        np.testing.assert_allclose(arb.combine(st), [6.0, 6.0])

    def test_publish_sink(self):
        from control import Arbiter, ExperimentState, Schema
        st = ExperimentState(Schema(["a"]))
        seen = {"n": 0}
        Arbiter().publish(st)                                         # дефолт no-op — не падает
        Arbiter(sink=lambda s: seen.__setitem__("n", seen["n"] + 1)).publish(st)
        self.assertEqual(seen["n"], 1)

    def test_scheduler_routes_through_arbiter(self):
        from control import Arbiter, ExperimentState, Scheduler, DemoPlant
        calls = {"combine": 0, "publish": 0}

        class Spy(Arbiter):
            def combine(self, state):
                calls["combine"] += 1
                return super().combine(state)

            def publish(self, state):
                calls["publish"] += 1
                return super().publish(state)

        plant = DemoPlant(n_corr=4, n_bpm=6)
        state = ExperimentState(plant.schema(), history=plant.history())
        sch = Scheduler(plant, state, dt=1.0, arbiter=Spy())
        sch.run(3)
        self.assertEqual(calls["combine"], 3)                         # команда собрана на каждом такте
        self.assertEqual(calls["publish"], 3)                         # выгрузка вызвана на каждом такте


class TestUndulatorGap(unittest.TestCase):
    """Живой зазор ондулятора: структура данных (видят контуры), gap→K→E₁, ScreenRecorder читает живой зазор."""

    def test_gap_visible_to_control_and_K_monotonic(self):
        from control import UndulatorConfig, ExperimentState, Schema, StateView
        st = ExperimentState(Schema(["a"]))
        st.undulators["und_01"] = UndulatorConfig(5.0, name="und_01")
        u = StateView(st).undulators["und_01"]                # контур видит зазор через StateView
        self.assertEqual(u.gap_mm, 5.0)
        k5 = u.K()
        u.gap_mm = 10.0                                        # двигаем зазор на лету
        k10 = u.K()
        self.assertIsNotNone(k5); self.assertIsNotNone(k10)
        self.assertLess(k10, k5)                              # открыли зазор → K упал (поле слабее), пересчёт живой

    def test_gap_raises_harmonic_energy(self):
        """Открытие зазора (K↓) монотонно поднимает энергию 1-й гармоники E₁; E₁(5мм)≈7 кэВ."""
        import control.state as cs
        sr = cs._undulator_sr()
        if sr is None:
            self.skipTest("undulator_sr недоступен")
        e = lambda g: sr.harmonic_energy_keV(1, sr.K_of_B0(sr.B0_of_gap(g)), 0.0)
        es = [e(g) for g in (5.0, 7.0, 10.0, 15.0)]
        self.assertTrue(all(b > a for a, b in zip(es, es[1:])))   # монотонно растёт с зазором
        self.assertAlmostEqual(es[0], 7.0, delta=0.6)

    def test_screen_recorder_reads_live_gap(self):
        """ScreenRecorder берёт зазор ондулятора-источника из state.undulators — разный зазор → разная E₁."""
        from control import UndulatorConfig, ExperimentState, Schema
        from control.recorders import ScreenRecorder
        st = ExperimentState(Schema(["a"]))
        st.undulators["und_19"] = UndulatorConfig(5.0, name="und_19")
        st.observables["photon.x"] = np.array([1e-7]); st.observables["photon.y"] = np.array([1e-7])
        rec = ScreenRecorder(params="SCR_U19", mesh=31)
        rec.sample(st, 0, None)                               # зазор 5 мм
        st.undulators["und_19"].gap_mm = 15.0
        rec.sample(st, 1, None)                               # зазор 15 мм (живая правка)
        rec.finalize(st, None)
        art = rec.artifact()[1]
        self.assertAlmostEqual(float(art["gap_start_mm"][0]), 5.0)
        self.assertAlmostEqual(float(art["gap_end_mm"][0]), 15.0)
        self.assertGreater(float(art["energy_end_keV"][0]), float(art["energy_start_keV"][0]))  # E₁ выросла с зазором

    def test_gap_drive_slew_and_odz(self):
        """Привод: реализованный зазор едет к уставке со скоростью slew (не телепорт), клип в ОДЗ; уставка==старт → стоит."""
        from control import UndulatorConfig
        u = UndulatorConfig(5.0, slew_mm_s=0.5, gap_min=5.0, gap_max=20.0)
        u.advance(2.0)
        self.assertEqual(u.gap_mm, 5.0)                   # уставка==старт → привод стоит (no-op)
        u.command(15.0)
        u.advance(2.0)                                     # 0.5 мм/с × 2 с = 1 мм (не телепорт)
        self.assertAlmostEqual(u.gap_mm, 6.0)
        for _ in range(100):
            u.advance(2.0)
        self.assertAlmostEqual(u.gap_mm, 15.0)            # доехал до уставки и встал
        u.command(999.0)                                   # уставка за пределом ОДЗ
        self.assertAlmostEqual(u.gap_target, 20.0)        # клипнута в ОДЗ (ход челюстей)
        for _ in range(200):
            u.advance(2.0)
        self.assertAlmostEqual(u.gap_mm, 20.0)

    def test_orbit_kick_residual(self):
        """Остаточный орбитный кик от зазора: 0 при идеальном замыкании; ненулевой ~мкрад и убывает с открытием."""
        from control import UndulatorConfig
        self.assertEqual(UndulatorConfig(5.0).orbit_kick_rad(), 0.0)      # i1_residual=0 → идеально замкнут, кика нет
        u = UndulatorConfig(5.0, i1_residual=3e-5)
        k5 = u.orbit_kick_rad()
        self.assertGreater(abs(k5), 0.0); self.assertLess(abs(k5), 1e-5)  # есть кик, масштаб ~мкрад
        u.gap_mm = 15.0
        self.assertLess(abs(u.orbit_kick_rad()), abs(k5))                 # открыли зазор → B0↓ → кик меньше


class TestOpticsView(unittest.TestCase):
    """state.optics: OpticsView отдаёт твисс по s и по имени магнита (без кольца, фейковый твисс)."""

    def test_lookup_by_s_and_magnet(self):
        from control import OpticsView
        tw = {"s": [0.0, 1.0, 2.0, 3.0], "betax": [5, 6, 7, 8], "betay": [3, 2, 1, 0.5],
              "alphax": [0, 0.1, 0.2, 0.3], "etax": [0.01, 0.02, 0.03, 0.04],
              "nux": 95.23, "nuy": 33.40, "alphac": 5e-5}

        class _Map:
            def index(self, name):
                return {"A": 1, "B": 3}[name]

        ov = OpticsView(tw, s_out=[0.5, 1.5, 2.5, 3.5], magnets=_Map())
        self.assertEqual(ov.tunes(), (95.23, 33.40))
        self.assertAlmostEqual(ov.alphac, 5e-5)
        self.assertAlmostEqual(ov.at_s(2.0)["betax"], 7.0)        # s=2 — узел сетки
        r = ov.at_magnet("B", where="exit")                       # s_out[3]=3.5 → ближайший срез s=3
        self.assertAlmostEqual(r["betax"], 8.0)
        self.assertIn("etax", r); self.assertIn("alphax", r)


class TestDriveSpeeds(unittest.TestCase):
    """Конечная скорость поля у ВСЕХ магнитов (не только корректоров): Limits.set_drive_speeds по типу актуатора."""

    def test_set_drive_speeds_by_type_and_apply(self):
        from control import Schema, Limits, DRIVE_SPEED_DEFAULTS
        sch = Schema(["corr.000.hkick", "quad.QF1A.dk", "sext.SD.dk2", "oct.OD1B.dk3", "rf.phase"])
        lim = Limits(sch)
        self.assertTrue(np.isinf(lim.slew).all())                 # дефолт — мгновенно (inf) у всех
        n = lim.set_drive_speeds(DRIVE_SPEED_DEFAULTS)
        self.assertEqual(n, 4)                                     # corr/quad/sext/oct получили скорость, rf — нет
        self.assertEqual(lim.slew[sch.index("quad.QF1A.dk")], DRIVE_SPEED_DEFAULTS["quad"])
        self.assertEqual(lim.slew[sch.index("oct.OD1B.dk3")], DRIVE_SPEED_DEFAULTS["oct"])
        self.assertTrue(np.isinf(lim.slew[sch.index("rf.phase")]))  # тип вне словаря → остался inf
        new = lim.apply(np.zeros(len(sch)), np.ones(len(sch)), dt=1.0)   # slew реально ограничивает шаг
        self.assertAlmostEqual(new[sch.index("quad.QF1A.dk")], DRIVE_SPEED_DEFAULTS["quad"])  # шаг ≤ slew·dt
        self.assertAlmostEqual(new[sch.index("rf.phase")], 1.0)                                # inf → мгновенно

    def test_set_profile_and_driveview(self):
        """Богатый профиль set_profile (slew+lag+ОДЗ) + DriveView: инфо о приводе по имени на доске."""
        from control import Schema, Limits, ExperimentState, StateView
        sch = Schema(["corr.000.hkick", "quad.QF1A.dk", "rf.phase"])
        lim = Limits(sch)
        prof = {"corr": {"slew": 1.0, "lag_tau": 2e-4, "odz": (-1e-3, 1e-3)},
                "quad": {"slew": 0.05, "lag_tau": 0.1}}
        n = lim.set_profile(prof)
        self.assertEqual(n, 2)                                            # corr, quad; rf не в профиле
        self.assertEqual(lim.slew[sch.index("corr.000.hkick")], 1.0)
        self.assertEqual(lim.lag_tau[sch.index("quad.QF1A.dk")], 0.1)
        self.assertEqual(lim.odz_hi[sch.index("corr.000.hkick")], 1e-3)   # ОДЗ из профиля
        self.assertTrue(np.isinf(lim.odz_hi[sch.index("quad.QF1A.dk")]))  # quad без odz → остался inf
        st = ExperimentState(sch, limits=lim)
        d = StateView(st).drives.of("corr.000.hkick")                     # инфо привода на доске по имени
        self.assertAlmostEqual(d["slew"], 1.0); self.assertAlmostEqual(d["lag_tau"], 2e-4)
        self.assertEqual(d["odz"], (-1e-3, 1e-3))
        self.assertIn("corr", StateView(st).drives.types())               # сводка по типам

    def test_first_order_lag(self):
        """Лаг 1-го порядка (полоса БП): actual приближается к уставке экспоненциально; τ=0 → без лага."""
        from control import Schema, Limits
        sch = Schema(["quad.QF1A.dk"])
        lim = Limits(sch); lim.set("quad.QF1A.dk", lag_tau=1.0)   # τ=1 с, slew=inf
        a = lim.apply(np.zeros(1), np.ones(1), dt=1.0)
        self.assertAlmostEqual(a[0], 1 - np.exp(-1.0), places=6)  # α=exp(-1): actual=(1-α)·уставка
        a = lim.apply(a, np.ones(1), dt=1.0)
        self.assertAlmostEqual(a[0], 1 - np.exp(-2.0), places=6)  # ещё такт — ближе к уставке
        self.assertAlmostEqual(Limits(sch).apply(np.zeros(1), np.ones(1), dt=1.0)[0], 1.0)  # τ=0 → без лага


class _MagMap:
    """Мини-карта имён магнитов (name→index) для OpticsView.at_magnet в тестах клиента."""
    def __init__(self, m): self._m = dict(m)
    def index(self, name): return self._m[name]


def _seed_client_state():
    """Посеянное состояние для тестов клиента: 2 актуатора, наблюдаемые, ондулятор, оптика, приводы."""
    from control import ExperimentState, Schema, Limits, UndulatorConfig, OpticsView
    sch = Schema(["corr.000.hkick", "corr.000.vkick"])
    lim = Limits(sch)
    lim.set("corr.000.hkick", slew=1.5, lag_tau=2e-4, odz=(-1e-3, 1e-3), res=1e-6, act=True)
    st = ExperimentState(sch, limits=lim)
    st.observables["ebpm.x"] = np.array([1.0e-6, -2.0e-6, 3.0e-6])
    st.observables["ebpm.y"] = np.array([0.0, 4.0e-6])
    st.undulators["und_01"] = UndulatorConfig(7.0, name="und_01", slew_mm_s=10.0, gap_min=5.0, gap_max=20.0)
    tw = {"s": [0.0, 1.0, 2.0], "betax": [10.0, 20.0, 30.0], "betay": [5.0, 6.0, 7.0],
          "etax": [0.1, 0.2, 0.3], "nux": 42.125, "nuy": 15.22, "alphac": 7.0e-5}
    st.optics = OpticsView(tw, s_out=np.array([1.0]), magnets=_MagMap({"q1": 0}))
    return st


class TestClientAddresses(unittest.TestCase):
    """Этап 0: адресная схема клиента — перечисление из посеянного состояния и резолв в верные значения
    (всё сериализуемо: float/bool/list/None, никаких numpy/живых объектов через шов Link)."""

    def setUp(self):
        from control import StateView
        self.st = _seed_client_state()
        self.view = StateView(self.st)

    def test_enumerate_lists_expected_addresses(self):
        from client import addresses as A
        addrs = set(A.enumerate_addresses(self.view))
        for expect in ("obs.ebpm.x", "obs.ebpm.x:rms", "obs.ebpm.y",
                       "act.corr.000.hkick.slew", "act.corr.000.hkick.act",      # канонические (шаг 4)
                       "undulator.und_01.gap", "undulator.und_01.K", "optics.nux"):
            self.assertIn(expect, addrs)
        self.assertFalse(any(a.startswith(("drive.", "actual.")) for a in addrs))   # старые — только псевдонимы

    def test_resolve_observable_vector_and_reductions(self):
        from client import addresses as A
        self.assertEqual(A.resolve(self.view, "obs.ebpm.x"), [1.0e-6, -2.0e-6, 3.0e-6])
        self.assertAlmostEqual(A.resolve(self.view, "obs.ebpm.x:rms"),
                               float(np.sqrt(np.mean(np.array([1e-6, -2e-6, 3e-6]) ** 2))))
        self.assertAlmostEqual(A.resolve(self.view, "obs.ebpm.x:max"), 3.0e-6)
        self.assertAlmostEqual(A.resolve(self.view, "obs.ebpm.x:2"), 3.0e-6)   # индекс элемента
        self.assertIsNone(A.resolve(self.view, "obs.nope.x"))                   # нет такой наблюдаемой

    def test_resolve_drive_fields(self):
        from client import addresses as A
        self.assertAlmostEqual(A.resolve(self.view, "drive.corr.000.hkick.slew"), 1.5)
        self.assertAlmostEqual(A.resolve(self.view, "drive.corr.000.hkick.lag_tau"), 2e-4)
        self.assertAlmostEqual(A.resolve(self.view, "drive.corr.000.hkick.odz_hi"), 1e-3)
        self.assertAlmostEqual(A.resolve(self.view, "drive.corr.000.hkick.odz_lo"), -1e-3)
        self.assertIs(A.resolve(self.view, "drive.corr.000.hkick.act"), True)
        self.assertIsNone(A.resolve(self.view, "drive.corr.000.hkick.bogus"))  # неизвестное поле

    def test_resolve_actual_corrector_kick(self):
        """`actual.<id>` — реализованный кик корректора из state.actual (живой дипольный момент H/V)."""
        from client import addresses as A
        self.st.actual["corr.000.hkick"] = 1.5e-4
        self.st.actual["corr.000.vkick"] = -2.0e-4
        self.assertAlmostEqual(A.resolve(self.view, "actual.corr.000.hkick"), 1.5e-4)
        self.assertAlmostEqual(A.resolve(self.view, "actual.corr.000.vkick"), -2.0e-4)
        self.assertIsNone(A.resolve(self.view, "actual.corr.999.hkick"))        # нет такого актуатора
        m = A.meta("actual.corr.000.hkick")
        self.assertEqual(m.units, "рад")
        self.assertEqual(m.widget, A.W_GRAPH)

    def test_resolve_undulator_gap_and_target(self):
        from client import addresses as A
        self.assertAlmostEqual(A.resolve(self.view, "undulator.und_01.gap"), 7.0)
        self.st.undulators["und_01"].command(9.0)                              # ручка: уставка
        self.assertAlmostEqual(A.resolve(self.view, "undulator.und_01.gap_target"), 9.0)
        self.assertAlmostEqual(A.resolve(self.view, "undulator.und_01.gap"), 7.0)  # реализованный ещё не поехал

    def test_resolve_optics_scalar_and_at_magnet(self):
        from client import addresses as A
        self.assertAlmostEqual(A.resolve(self.view, "optics.nux"), 42.125)
        self.assertAlmostEqual(A.resolve(self.view, "optics.alphac"), 7.0e-5)
        self.assertAlmostEqual(A.resolve(self.view, "optics.beta_x@q1"), 20.0)  # s_out[0]=1 → срез β=20
        self.assertIsNone(A.resolve(self.view, "optics.beta_x@nope"))           # нет магнита → None

    def test_resolve_always_serializable(self):
        from client import addresses as A
        for a in ("obs.ebpm.x", "obs.ebpm.x:rms", "drive.corr.000.hkick.slew",
                  "drive.corr.000.hkick.act", "undulator.und_01.gap", "optics.nux"):
            v = A.resolve(self.view, a)
            self.assertIsInstance(v, (float, bool, list), f"{a}→{type(v)}")
            if isinstance(v, list):
                self.assertTrue(all(isinstance(x, float) for x in v))

    def test_meta_widget_and_units(self):
        from client import addresses as A
        self.assertEqual(A.meta("obs.ebpm.x").widget, A.W_VECTOR)
        self.assertEqual(A.meta("obs.ebpm.x:rms").widget, A.W_GRAPH)
        self.assertEqual(A.meta("undulator.und_01.gap").units, "мм")
        self.assertEqual(A.meta("drive.corr.000.hkick.slew").kind, "act")        # псевдоним → область act
        self.assertEqual(A.meta("drive.corr.000.hkick.slew").addr, "act.corr.000.hkick.slew")


class TestClientDiagnostics(unittest.TestCase):
    """Стандартный график «η-узор и моды FOFB» (diag.orbit_modes): чистая сводка мод, сборка данных на
    фейках (без ring_cpp), сериализуемость через шов Link, кэш движка с инвалидацией по rcond, метаданные адреса."""

    def _fakes(self, rcond=1e-2):
        from types import SimpleNamespace as NS
        from control.controller import ResponseMatrix
        rng = np.random.default_rng(3)
        U, _ = np.linalg.qr(rng.normal(size=(8, 3)))
        V, _ = np.linalg.qr(rng.normal(size=(3, 3)))
        R = U @ np.diag([1.0, 0.5, 0.005]) @ V.T                    # 8 BPM × 3 корректора, одна слабая мода
        eta = 0.6 * U[:, 0] + 0.8 * U[:, 2]                          # η: 60% в сильной моде, 80% — в слабой
        act = ["corr.000.hkick", "corr.001.hkick", "corr.002.hkick"]
        ctrl = NS(act=act, sensor="ebpm.x", params={"rcond": rcond})
        loops = {"fofb_x": NS(controller=ctrl)}
        res = NS(response=lambda a, o: ResponseMatrix(R, [f"b{i}" for i in range(8)], a))
        mon = np.arange(2, 10)
        etax = np.zeros(12); etax[mon] = eta
        plant = NS(mon=mon, ring=None)
        state = NS(optics=NS(etax=etax))
        return R, eta, plant, state, loops, res, ctrl

    def test_modes_summary_counts_and_eta_reach(self):
        from client.diagnostics import modes_summary
        R, eta, *_ = self._fakes()
        d = modes_summary(R, eta, 1e-2)
        self.assertEqual(d["n_modes"], 3)
        self.assertEqual(d["n_kept"], 2)                             # 0.005 < 1e-2 — слабая мода отброшена
        np.testing.assert_allclose(d["sv"], [1.0, 0.5, 0.005], rtol=1e-12)
        self.assertAlmostEqual(d["eta_reach"], 0.6, places=12)       # в подпространстве 2 мод — только 60%
        self.assertAlmostEqual(d["eta_reach_full"], 1.0, places=12)
        self.assertEqual(modes_summary(R, eta, 1e-4)["n_kept"], 3)   # «без усечения» — все моды

    def test_orbit_modes_serializable_and_engine_cache(self):
        import json
        from types import SimpleNamespace as NS
        from client import InProcLink, addresses as A
        from client.diagnostics import orbit_modes
        from client.engine import EngineThread
        R, eta, plant, state, loops, res, ctrl = self._fakes()
        d = orbit_modes(plant, state, loops, res)
        json.dumps(d)                                                 # шов Link: только числа/строки
        np.testing.assert_allclose(d["eta_bpm"], eta)                 # η в BPM — как у датчика контура ВЧ
        np.testing.assert_allclose(np.asarray(d["cols"]).T, R)        # столбцы = бетатронные волны корректоров
        self.assertEqual(d["corr_ids"], ctrl.act)
        self.assertIn("error", orbit_modes(plant, state, {}, res))    # нет контура — ошибка, не исключение
        eng = EngineThread(InProcLink(), NS(state=state, plant=plant, loops=loops, resources=res))
        v1 = eng._diag_value("diag.orbit_modes")
        self.assertIs(eng._diag_value("diag.orbit_modes"), v1)        # статично: считается один раз
        ctrl.params["rcond"] = 1e-4                                   # смена rcond → пересчёт
        v2 = eng._diag_value("diag.orbit_modes")
        self.assertIsNot(v2, v1)
        self.assertEqual(v2["n_kept"], 3)
        self.assertIn("error", eng._diag_value("diag.bogus"))
        m = A.meta("diag.orbit_modes")
        self.assertEqual((m.kind, m.widget, m.name), ("diag", A.W_INFO, "η-узор и моды FOFB"))
        self.assertIsNone(A.resolve(None, "diag.orbit_modes"))        # считает движок, не resolve


class TestClientLink(unittest.TestCase):
    """Этап 0: InProcLink проходит круг subscribe→snapshot на ФЕЙКОВОМ цикле (без реального движка).
    Проверяем транспорт: подписки ведутся, команды/контроль доходят сообщениями, снимок с числами и
    метой возвращается последним-побеждает, poll до публикации — None, unsubscribe снимает подписку."""

    def test_subscribe_command_snapshot_round_trip(self):
        from client import InProcLink, Snapshot, Message, addresses as A
        from client.link import COMMAND, CONTROL
        from control import StateView
        st = _seed_client_state()
        link = InProcLink()

        self.assertIsNone(link.poll_snapshot())                    # до публикации — пусто

        # --- клиент ---
        link.subscribe("obs.ebpm.x")
        link.subscribe("undulator.und_01.gap_target")
        link.command("undulator.und_01.gap_target", 9.0)
        link.control("run")

        # --- фейковый движок: осушить сообщения, применить, опубликовать снимок подписанного ---
        msgs = link.pull()
        kinds = [m.kind for m in msgs]
        self.assertIn(COMMAND, kinds)
        self.assertIn(CONTROL, kinds)
        self.assertEqual(next(m.value for m in msgs if m.kind == CONTROL), "run")
        for m in msgs:                                              # штатный путь ручки: команда → доска
            if m.kind == COMMAND and m.addr == "undulator.und_01.gap_target":
                st.undulators["und_01"].command(m.value)
        self.assertEqual(link.subscriptions(), {"obs.ebpm.x", "undulator.und_01.gap_target"})
        view = StateView(st)
        values = {a: A.resolve(view, a) for a in link.subscriptions()}
        link.publish(Snapshot(values=values, meta={"model_time": 0.0, "cycle": 1, "seq": 1, "speed": 0.0}))

        # --- клиент видит кадр ---
        snap = link.poll_snapshot()
        self.assertIsNotNone(snap)
        self.assertEqual(snap.values["obs.ebpm.x"], [1.0e-6, -2.0e-6, 3.0e-6])
        self.assertAlmostEqual(snap.values["undulator.und_01.gap_target"], 9.0)  # ручка отразилась
        self.assertEqual(snap.meta["cycle"], 1)

    def test_unsubscribe_and_latest_wins(self):
        from client import InProcLink, Snapshot
        link = InProcLink()
        link.subscribe("a"); link.subscribe("b")
        link.pull()
        self.assertEqual(link.subscriptions(), {"a", "b"})
        link.unsubscribe("a"); link.pull()
        self.assertEqual(link.subscriptions(), {"b"})
        link.publish(Snapshot(values={"b": 1.0}, meta={"seq": 1}))
        link.publish(Snapshot(values={"b": 2.0}, meta={"seq": 2}))   # последний побеждает
        self.assertEqual(link.poll_snapshot().values["b"], 2.0)
        self.assertEqual(link.poll_snapshot().meta["seq"], 2)


class _FakeInfinitePlant:
    """Лёгкий plant по контракту Plant для тестов движка БЕЗ ring_cpp: ведёт зазор штатно (plant.advance),
    двигает «орбиту» с каждым шагом (чтобы было что публиковать). Мгновенный привод зазора (slew велик)."""

    def __init__(self):
        from control import Schema
        self._schema = Schema(["corr.000.hkick", "corr.000.vkick"])
        self.stride = 1
        self.TREV = 1.0e-6
        self.turn = 0
        self.alive = True

    def schema(self):
        return self._schema

    def history(self):
        return None

    def reset(self, state):
        from control import UndulatorConfig
        self.turn = 0
        self.alive = True
        state.observables["ebpm.x"] = np.zeros(3)
        state.undulators.setdefault(
            "und_01", UndulatorConfig(7.0, name="und_01", slew_mm_s=1.0e6, gap_min=5.0, gap_max=20.0))

    def step(self, state, dt):
        for u in state.undulators.values():          # ШТАТНЫЙ путь: plant реализует зазор к уставке
            u.advance(dt)
        self.turn += 1
        state.observables["ebpm.x"] = np.array([1.0e-6 * self.turn, -1.0e-6, 2.0e-6])

    def fingerprint(self):
        return ("fake-inf",)


class TestClientEngine(unittest.TestCase):
    """Этап 1 (headless, фейк-plant — всегда): EngineThread гоняет поток, фейк-клиент подписан на obs.ebpm.x
    и зазор, командует смену зазора и видит отражение; оценка скорости > 0; пуск/стоп управляют циклом."""

    def _wait(self, link, cond, timeout=5.0):
        import time as _t
        end = _t.time() + timeout
        snap = None
        while _t.time() < end:
            snap = link.poll_snapshot()
            if snap is not None and cond(snap):
                return snap
            _t.sleep(0.005)
        return snap

    def test_thread_runs_commands_and_reflects_gap(self):
        import time as _t
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler

        plant = _FakeInfinitePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        link = InProcLink()
        eng = EngineThread(link, sch, steps_per_cycle=1, idle_sleep=0.002)

        link.subscribe("obs.ebpm.x")
        link.subscribe("undulator.und_01.gap")
        link.subscribe("undulator.und_01.gap_target")
        link.subscribe("meta.speed")
        link.control("run")
        link.command("undulator.und_01.gap", 9.0)      # ручка: штатный путь команда → доска → plant
        link.command("undulator.bogus.gap", 1.0)       # нераспознанное — не должно ронять поток

        eng.start()
        try:
            snap = self._wait(link, lambda s: s.meta.get("cycle", 0) >= 3
                              and abs((s.values.get("undulator.und_01.gap") or 0.0) - 9.0) < 1e-6)
            self.assertIsNotNone(snap, "движок не опубликовал ожидаемый кадр")
            self.assertGreaterEqual(snap.meta["cycle"], 3)                 # поток реально шагал
            self.assertAlmostEqual(snap.values["undulator.und_01.gap_target"], 9.0)
            self.assertAlmostEqual(snap.values["undulator.und_01.gap"], 9.0)   # реализованный зазор отразил ручку
            self.assertEqual(len(snap.values["obs.ebpm.x"]), 3)            # наблюдаемая течёт
            self.assertGreater(snap.meta["speed"], 0.0)                    # оценка скорости > 0
            self.assertEqual(snap.meta.get("unknown_command"), "undulator.bogus.gap")

            # --- стоп замораживает цикл ---
            link.control("stop")
            _t.sleep(0.05)
            c0 = link.poll_snapshot().meta["cycle"]
            _t.sleep(0.1)
            self.assertEqual(link.poll_snapshot().meta["cycle"], c0)       # на паузе цикл не растёт
            self.assertFalse(link.poll_snapshot().meta["running"])
        finally:
            eng.stop()
            eng.join(timeout=5)
        self.assertFalse(eng.is_alive())                                   # поток погашен

    def test_restart_resets_model_time(self):
        import time as _t
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler

        plant = _FakeInfinitePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        link = InProcLink()
        eng = EngineThread(link, sch, steps_per_cycle=1, idle_sleep=0.002)
        link.subscribe("meta.model_time")
        link.control("run")
        eng.start()
        try:
            self._wait(link, lambda s: s.meta.get("model_time", 0.0) > 3e-3)
            link.control("stop")                                           # пауза — чтобы рестарт не гнал мгновенно
            link.control("restart")                                        # сброс (сохраняет режим паузы)
            snap = self._wait(link, lambda s: not s.meta.get("running", True)
                              and s.meta.get("model_time", 1.0) == 0.0)
            self.assertIsNotNone(snap)
            self.assertEqual(snap.meta["model_time"], 0.0)                 # время сброшено рестартом
            self.assertFalse(snap.meta["running"])
        finally:
            eng.stop()
            eng.join(timeout=5)

    def test_publish_cadence_decoupled_from_steps(self):
        """Публикация раз в publish_every тактов, развязана от шага. Фейк-plant: turn += 1/такт, stride=1;
        в момент снимка turn == cycle·publish_every (движок трекает КАЖДЫЙ такт, а снимок — раз в K)."""
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler

        K = 4
        plant = _FakeInfinitePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        link = InProcLink()
        eng = EngineThread(link, sch, steps_per_cycle=1, publish_every=K, idle_sleep=0.002)
        self.assertEqual(eng.publish_every, K)
        link.subscribe("meta.turn")
        link.control("run")
        eng.start()
        try:
            snap = self._wait(link, lambda s: s.meta.get("cycle", 0) >= 2)
            self.assertIsNotNone(snap, "движок не опубликовал кадр")
            self.assertGreaterEqual(snap.meta["cycle"], 2)
            # инвариант развязки: трекается каждый такт, публикуется раз в K → turn кратен K·cycle
            self.assertEqual(snap.meta["turn"], snap.meta["cycle"] * K)
        finally:
            eng.stop()
            eng.join(timeout=5)
        self.assertFalse(eng.is_alive())

    def test_publish_every_clamped_to_one(self):
        """publish_every < 1 зажимается в 1 (каждый такт = снимок) — не делить на ноль / не зависнуть."""
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler

        plant = _FakeInfinitePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        eng = EngineThread(InProcLink(), sch, steps_per_cycle=1, publish_every=0)
        self.assertEqual(eng.publish_every, 1)


class _StubMag:
    def __init__(self, x, y, kind="drift"):
        self.center = np.array([float(x), float(y), 0.0])
        self.kind = kind


class _StubRing:
    def __init__(self, mags):
        self.mags = mags
        self.circumference = 10.0


class _StubRingPlant:
    """Фейк-plant с минимальным кольцом (без ring_cpp) — для headless-тестов геометрии сцены."""
    def __init__(self):
        self.ring = _StubRing([_StubMag(0, 0), _StubMag(1, 0, "monitor"),
                               _StubMag(2, 1), _StubMag(3, 1, "monitor")])
        self.lattice = 26
        self.gap_mm = 5.0
        self._und_anchor = {"und_01": 2}


class TestClientScreenBackground(unittest.TestCase):
    """Экран (тяжёлая интерференция) считается в ФОНОВОМ потоке: _screen_value не ждёт расчёта, пока он идёт —
    отдаёт прежнюю картину с флагом computing (первый раз — только флаг), при неизменном входе не пересчитывает,
    новые настройки окна запускают пересчёт. Вход фонового расчёта — снимок доски, не живое состояние."""

    def test_screen_computed_in_background(self):
        import time as _t
        import types
        from client.engine import EngineThread
        st = types.SimpleNamespace(observables={"photon.x": np.array([1e-6])}, undulators={})
        eng = EngineThread(link=None, scheduler=types.SimpleNamespace(state=st, plant=None))
        calls = []

        def slow(addr, cfg=None, inp=None):
            calls.append((addr, dict(cfg or {}), inp)); _t.sleep(0.3)
            return {"heatmap": [[1.0]], "gap": 5.0}

        eng._compute_screen = slow
        try:
            t0 = _t.perf_counter(); v = eng._screen_value("screen.SCR_X")
            self.assertLess(_t.perf_counter() - t0, 0.1, "поток модели ждал расчёта экрана")
            self.assertEqual(v, {"computing": True})
            inp = calls[0][2]
            self.assertIsNot(inp.observables["photon.x"], st.observables["photon.x"])   # снимок, а не ссылка
            _t.sleep(0.5)
            self.assertEqual(eng._screen_value("screen.SCR_X"), {"heatmap": [[1.0]], "gap": 5.0})
            n = len(calls); eng._screen_value("screen.SCR_X"); self.assertEqual(len(calls), n)  # вход тот же
            eng._screens["screen.SCR_X"] = {"res": 61}
            v = eng._screen_value("screen.SCR_X")
            self.assertTrue(v.get("computing")); self.assertEqual(v.get("gap"), 5.0)   # прежняя картина + флаг
            end = _t.time() + 2.0
            while len(calls) == n and _t.time() < end:                   # фоновый поток взял новый расчёт
                _t.sleep(0.01)
            self.assertEqual(calls[-1][1].get("res"), 61)
        finally:
            eng.stop()

    def test_engine_bend_screen_power_and_flux(self):
        """Экран у линейки диполя (gated на БД): движок отдаёт карту мощности/потока и сплошной спектр, kind=bend."""
        import types
        tb = os.path.join(ROOT, "workspace", "tables")
        if not os.path.exists(os.path.join(tb, "element_groups.csv")):
            self.skipTest("нет element_groups.csv")
        und = os.path.join(ROOT, "undulators")
        if und not in sys.path:
            sys.path.insert(0, und)
        cwd = os.getcwd(); os.chdir(ROOT)
        from client.engine import EngineThread
        st = types.SimpleNamespace(observables={}, undulators={})
        eng = EngineThread(link=None, scheduler=types.SimpleNamespace(state=st, plant=types.SimpleNamespace(lattice=26)))
        try:
            import screen_analysis as sa
            sid = sa.bend_line_screens(26)[0]["sid"]
            p = eng._compute_screen(f"screen.{sid}", {"res": 41, "quantity": "power"}, st)
            f = eng._compute_screen(f"screen.{sid}", {"res": 41, "quantity": "flux", "point": (0.0, 0.0)}, st)
            u = eng._bend_context("SCR_und_01")
        finally:
            eng.stop(); os.chdir(cwd)
        self.assertNotIn("error", p)
        self.assertEqual((p["kind"], p["quantity"], f["quantity"]), ("bend", "power", "flux"))
        self.assertGreater(p["total"], 0.0)
        self.assertEqual(np.asarray(f["heatmap"]).shape, (41, 41))
        self.assertGreater(max(f["spectrum"]["flux"]), 0.0)
        self.assertIsNone(u)                                                    # экран ондулятора — прежний путь

    def test_screen_quantity_switch_recomputes_power(self):
        """Переключатель величины окна экрана: смена flux→power — новый расчёт; power идёт в карту мощности."""
        import time as _t
        import types
        from client.engine import EngineThread, _quantity
        self.assertEqual(_quantity({}), "flux")
        self.assertEqual(_quantity({"quantity": "power"}), "power")
        st = types.SimpleNamespace(observables={}, undulators={})
        eng = EngineThread(link=None, scheduler=types.SimpleNamespace(state=st, plant=None))
        calls = []
        eng._compute_screen = lambda addr, cfg=None, inp=None: calls.append(dict(cfg or {})) or {"heatmap": [[1.0]]}
        eng._compute_power = lambda *a: {"heatmap": [[2.0]], "quantity": "power"}
        try:
            eng._screen_value("screen.SCR_X"); _t.sleep(0.2); eng._screen_value("screen.SCR_X")
            n = len(calls)
            eng._screens["screen.SCR_X"] = {"res": 121, "quantity": "power"}
            eng._screen_value("screen.SCR_X")
            end = _t.time() + 2.0
            while len(calls) == n and _t.time() < end:
                _t.sleep(0.01)
            self.assertEqual(calls[-1].get("quantity"), "power")
        finally:
            eng.stop()
        real = EngineThread(link=None, scheduler=types.SimpleNamespace(state=st, plant=None))
        try:
            real._compute_power = lambda sid, mesh, pt, s: {"quantity": "power", "sid": sid, "mesh": mesh}
            self.assertEqual(real._compute_screen("screen.SCR_Y", {"res": 61, "quantity": "power"}, st),
                             {"quantity": "power", "sid": "SCR_Y", "mesh": 61})
        finally:
            real.stop()


class TestClientScene(unittest.TestCase):
    """Этап 2 (headless): чистая геометрия сцены — узел из plant'а и модель актёров, без PyVista/Qt."""

    def test_build_geometry_from_stub_ring(self):
        from client import scene3d
        from control import ExperimentState, Schema, UndulatorConfig
        st = ExperimentState(Schema(["a"]))
        st.undulators["und_01"] = UndulatorConfig(7.0, name="und_01")
        g = scene3d.build_geometry(_StubRingPlant(), st)
        self.assertEqual(len(g["centers"]), 4)
        self.assertEqual(len(g["mon_xy"]), 2)                    # два монитора
        self.assertIn("und_01", g["undulators"])
        self.assertEqual(g["undulators"]["und_01"]["xy"], [2.0, 1.0])
        self.assertAlmostEqual(g["undulators"]["und_01"]["gap_mm"], 7.0)  # живой зазор с доски

    @staticmethod
    def _nav_bundle():
        named = [{"name": "Q1", "xy": [0.0, 0.0], "kind": "quad", "dev": "QF", "s_idx": 0},
                 {"name": "C1", "xy": [1.0, 0.0], "kind": "corrector", "dev": "COR", "s_idx": 1},
                 {"name": "C1.c1", "xy": [1.0, 0.0], "kind": "corrector", "dev": "COR_c1", "s_idx": 2},
                 {"name": "BPM1", "xy": [1.5, 0.0], "kind": "monitor", "dev": "BPM", "s_idx": 3},
                 {"name": "Q2", "xy": [5.0, 0.0], "kind": "quad", "dev": "QF", "s_idx": 10},
                 {"name": "B1", "xy": [6.0, 0.0], "kind": "dipole", "dev": "DP", "s_idx": 11}]
        und = {"u1": {"xy": [3.2, 0.0], "center": [3.0, 0.0], "s_idx": 5},
               "u2": {"xy": [7.2, 0.0], "center": [7.0, 0.0], "s_idx": 12},
               "fake": {"xy": [0.0, 0.0]}}                          # без положения по пучку — не в навигации
        return {"named": named, "undulators": und}

    def test_nearest_screen_by_kind_and_stub_has_no_bend_screens(self):
        """Экран у линейки диполя ищется по ближайшему центру среди своего вида; фейк-кольцо тестов (модулей
        линеек нет) экранов у линеек не получает."""
        from client import scene3d
        m = scene3d.SceneModel({"screens": [
            {"sid": "SCR_u", "kind": "undulator", "center": [0.0, 0.0], "src": [0, 0], "axis": [1, 0]},
            {"sid": "SCR_DP1", "kind": "bend", "center": [10.0, 0.0], "src": [0, 0], "axis": [1, 0]},
            {"sid": "SCR_DP2", "kind": "bend", "center": [20.0, 0.0], "src": [0, 0], "axis": [1, 0]}]})
        self.assertEqual(m.nearest_screen(1.0, 0.0)["sid"], "SCR_u")
        self.assertEqual(m.nearest_screen(1.0, 0.0, kind="bend")["sid"], "SCR_DP1")
        self.assertEqual(m.nearest_screen(18.0, 1.0, kind="bend")["sid"], "SCR_DP2")
        g = scene3d.build_geometry(_StubRingPlant())
        self.assertFalse([s for s in g["screens"] if s.get("kind") == "bend"])

    def test_nav_beam_order_merge_and_same_type(self):
        """Навигация камеры: порядок по пучку, ондулятор — один шаг в своём месте, строки в одной точке
        склеены, мониторы пропущены; same_type — следующий того же типа; кольцо замкнуто."""
        from client import scene3d
        m = scene3d.SceneModel(self._nav_bundle())
        self.assertEqual([it["name"] for it in m.nav], ["Q1", "C1", "u1", "Q2", "B1", "u2"])
        ix = {it["name"]: k for k, it in enumerate(m.nav)}
        self.assertEqual(m.nav[ix["u1"]]["xy"], [3.0, 0.0])              # центр гребёнки, не якорь
        self.assertEqual(m.nav_step(ix["Q1"], +1), ix["C1"])             # → по пучку
        self.assertEqual(m.nav_step(ix["Q1"], -1), ix["u2"])             # ← навстречу, через «конец» кольца
        self.assertEqual(m.nav_step(ix["u1"], +1, same_type=True), ix["u2"])
        self.assertEqual(m.nav_step(ix["u2"], +1, same_type=True), ix["u1"])   # по кругу
        self.assertEqual(m.nav_step(ix["Q1"], +1, same_type=True), ix["Q2"])   # тот же тип устройства
        self.assertEqual(m.nav_step(ix["B1"], +1, same_type=True), ix["B1"])   # единственный своего типа
        self.assertEqual(m.nav_nearest(3.01, 0.0)[0], ix["u1"])
        self.assertIn("ондулятор u1 — 1 из 2", m.nav_label(ix["u1"]))
        self.assertIn("квадруполь Q2 — 2 из 2", m.nav_label(ix["Q2"]))
        self.assertEqual(scene3d.SceneModel({}).nav, [])

    def test_scene_model_jaws_monotonic_in_gap(self):
        from client import scene3d
        model = scene3d.SceneModel({"undulators": {"u": {"xy": [1.0, 2.0]}}}, gap_scale=100.0)
        up5, lo5 = model.undulator_jaws("u", 5.0)
        up10, lo10 = model.undulator_jaws("u", 10.0)
        self.assertEqual(up5[:2], [1.0, 2.0])                    # челюсть над осью пучка (x,y)
        self.assertGreater(up10[2], up5[2])                      # шире зазор → выше верхняя челюсть
        self.assertLess(lo10[2], lo5[2])                         # и ниже нижняя
        self.assertAlmostEqual(up5[2], -lo5[2])                  # симметрично вокруг оси

    def test_scene_model_empty_ring_is_safe(self):
        from client import scene3d
        from control import ExperimentState, Schema
        g = scene3d.build_geometry(_FakeInfinitePlant(), ExperimentState(Schema(["a"])))
        model = scene3d.SceneModel(g)
        self.assertEqual(len(model.ring_polyline()), 0)          # нет кольца — пустой обвод
        self.assertIsNone(model.start_plane())
        self.assertIsNone(model.bbox())

    def test_named_elements_and_picking(self):
        from client import scene3d
        from control import ExperimentState, Schema
        g = scene3d.build_geometry(_StubRingPlant(), ExperimentState(Schema(["a"])))
        self.assertEqual(len(g["named"]), 4)                     # все 4 стаб-элемента адресуемы
        model = scene3d.SceneModel(g)
        hit = model.nearest_named(2.9, 1.1)                      # ближе всего к (3,1) — SP01.BPM.1? имя из стаба нет
        self.assertIsNotNone(hit)
        idx, name, kind, dist = hit
        self.assertEqual(model.named[idx]["xy"], [3.0, 1.0])     # верный ближайший элемент
        self.assertLess(dist, 0.3)
        self.assertEqual(len(model.named_points()), 4)

    def test_composite_plant_delegates_to_tbt(self):
        """CompositePlant делегирует ring-атрибуты в .tbt (фикс пустой 3D-сцены/меты при optics=True)."""
        import types
        from control.composite_plant import CompositePlant
        cp = object.__new__(CompositePlant)
        cp.tbt = types.SimpleNamespace(ring="RING", lattice=26, turn=7, alive=True, corr=[0, 5])
        self.assertEqual(cp.ring, "RING")
        self.assertEqual(cp.lattice, 26)
        self.assertEqual(cp.turn, 7)
        self.assertIs(cp.alive, True)
        self.assertEqual(cp.corr, [0, 5])
        with self.assertRaises(AttributeError):
            _ = cp.no_such_attr_xyz

    def test_corrector_aids_and_tbt_unwrap(self):
        """build_geometry разворачивает CompositePlant (.tbt) и даёт корректору actuator-id (actual-адрес кика)."""
        from client import scene3d
        from control import ExperimentState, Schema

        class _M:
            def __init__(self, x, kind, dev=""):
                self.center = np.array([float(x), 0.0, 0.0]); self.kind = kind; self.dev = dev

        class _Ring:
            def __init__(self, mags): self.mags = mags; self.circumference = 10.0

        class _Tbt:
            def __init__(self):
                self.ring = _Ring([_M(0, "corrector", "und_SH1A"), _M(1, "quad", "QP"), _M(2, "monitor", "BPM")])
                self.lattice = 26; self.gap_mm = 5.0; self._und_anchor = {}; self.corr = [0]   # корректор ring-idx 0 → актуатор 0

        class _Comp:                                       # имитация CompositePlant: делегирует в .tbt
            def __init__(self): self.tbt = _Tbt()
            def __getattr__(self, n): return getattr(self.__dict__["tbt"], n)

        g = scene3d.build_geometry(_Comp(), ExperimentState(Schema(["a"])))
        self.assertGreaterEqual(len(g["named"]), 3)                # сцена НЕ пустая (разворот .tbt сработал)
        corr = [n for n in g["named"] if n["kind"] == "corrector"]
        self.assertTrue(corr and corr[0].get("aids") == ["corr.000.hkick", "corr.000.vkick"])
        quad = [n for n in g["named"] if n["kind"] == "quad"]
        self.assertTrue(quad and "aids" not in quad[0])            # у не-корректора aids нет (семейного id нет в схеме)

    def test_family_gradient_aids_and_realized_dk(self):
        """Квадруполь/секступоль/октуполь получают id СЕМЕЙНОЙ Δk (контур оптики) — если он есть в схеме
        plant'а; info_fields добавляет живое поле actual.<famid>. Так магнит показывает реализованную контуром
        коррекцию (двигается при вкл/выкл ОС), а element.k* остаётся проектным."""
        from client import scene3d
        from client.windows import info_fields
        from control import ExperimentState, Schema

        class _M:
            def __init__(self, x, kind, dev):
                self.center = np.array([float(x), 0.0, 0.0]); self.kind = kind; self.dev = dev

        class _Ring:
            def __init__(self, mags): self.mags = mags; self.circumference = 10.0

        class _P:
            def __init__(self):
                self.ring = _Ring([_M(0, "quad", "und_QF1A"), _M(1, "sext", "und_SD1A"),
                                   _M(2, "oct", "und_OD1B")])
                self.lattice = 26; self.gap_mm = 5.0; self._und_anchor = {}
                # схема содержит семейные градиентные id (как у CompositePlant при optics=True)
                self._ids = ["quad.und_QF1A.dk", "sext.und_SD1A.dk2", "oct.und_OD1B.dk3"]
            def schema(self): return Schema(self._ids)

        g = scene3d.build_geometry(_P(), ExperimentState(Schema(["a"])))
        byk = {n["kind"]: n for n in g["named"]}
        self.assertEqual(byk["quad"].get("aids"), ["quad.und_QF1A.dk"])   # семейный id присвоен
        self.assertEqual(byk["sext"].get("aids"), ["sext.und_SD1A.dk2"])
        self.assertEqual(byk["oct"].get("aids"), ["oct.und_OD1B.dk3"])
        # info_fields проводит его в живое поле actual.<famid> рядом с проектным element.k*
        rows = [r[0] for r in info_fields("quad", byk["quad"]["name"], byk["quad"]["aids"])]
        self.assertIn("act.quad.und_QF1A.dk", rows)                       # ЖИВАЯ Δk контура (бывш. actual.*)
        self.assertIn(f"element.{byk['quad']['name']}.k1", rows)          # проектный градиент остаётся
        # нет семейного id в схеме → aids не присваивается (optics=False), поля actual нет
        class _P0(_P):
            def schema(self): return Schema(["corr.000.hkick"])
        g0 = scene3d.build_geometry(_P0(), ExperimentState(Schema(["a"])))
        self.assertNotIn("aids", [n for n in g0["named"] if n["kind"] == "quad"][0])

    def test_screen_live_gap_fallback_from_id(self):
        """ScreenRecorder._live_gap: экран БЕЗ строки source_ref в screens.csv берёт ондулятор из id экрана
        (SCR_<gid>→<gid>). Раньше такой экран читал зашитый GAP_MM и выглядел как включённый."""
        from control.recorders import ScreenRecorder, GAP_MM
        from control.state import UndulatorConfig

        class _St:
            undulators = {"und_01": UndulatorConfig(gid="und_01", gap_mm=5.0),
                          "und_07": UndulatorConfig(gid="und_07", gap_mm=50.0)}
        st = _St()
        self.assertAlmostEqual(ScreenRecorder(params="SCR_und_07")._live_gap(st), 50.0)   # фолбэк из id
        self.assertNotAlmostEqual(50.0, GAP_MM)                          # и это НЕ зашитый дефолт
        # неизвестный источник → зашитый GAP_MM (поведение по умолчанию сохранено)
        self.assertAlmostEqual(ScreenRecorder(params="SCR_zzz")._live_gap(st), GAP_MM)

    def test_magnet_boxes_oriented_and_typed(self):
        from client import scene3d
        g = {"centers": [[0, 0], [10, 0]],
             "magnets": [{"name": "Q1", "c": [5.0, 0.0], "ax": [1.0, 0.0], "ux": [0.0, 1.0], "L": 0.4, "kind": "quad"},
                         {"name": "B1", "c": [6.0, 0.0], "ax": [1.0, 0.0], "ux": [0.0, 1.0], "L": 1.0, "kind": "dipole"}]}
        boxes = scene3d.SceneModel(g).magnet_boxes()
        self.assertEqual(len(boxes), 2)
        c, a, u, hl, hw, hh, kind = boxes[0]
        self.assertEqual(kind, "quad")
        self.assertTrue(hl > 0 and hw > 0 and hh > 0)            # невырожденный блок
        self.assertAlmostEqual(float(a[0]), 1.0)                 # ось вдоль пучка (x)
        self.assertGreater(boxes[1][3], boxes[0][3])             # диполь длиннее квадруполя → больше блок

    def test_build_geometry_includes_magnet_bodies(self):
        from client import scene3d

        class _M:
            def __init__(self, x, kind, L):
                self.center = np.array([float(x), 0.0, 0.0]); self.kind = kind; self.L = float(L)
                self.dev = ""; self.axis = np.array([1.0, 0.0, 0.0]); self.uhat = np.array([0.0, 1.0, 0.0])

        class _R:
            def __init__(self):
                self.mags = [_M(0, "drift", 0.0), _M(1, "quad", 0.4), _M(2, "dipole", 1.0)]
                self.circumference = 5.0

        class _P:
            def __init__(self):
                self.ring = _R(); self.lattice = 26; self.gap_mm = 5.0; self._und_anchor = {}

        g = scene3d.build_geometry(_P())
        self.assertEqual(sorted(mg["kind"] for mg in g["magnets"]), ["dipole", "quad"])  # дрейф тела не имеет
        self.assertEqual(g["magnets"][0]["ax"], [1.0, 0.0])

    def test_pick_only_magnets_with_distance_gate(self):
        from client import scene3d
        g = {"centers": [[0, 0], [10, 0]],
             "named": [{"name": "Q1", "xy": [0.0, 0.0], "kind": "quad"},
                       {"name": "DR", "xy": [5.0, 0.0], "kind": "drift"}]}
        m = scene3d.SceneModel(g)
        self.assertEqual(m.nearest_named(5.0, 0.0)[1], "Q1")           # дрейф не пикаем → ближайший магнит
        self.assertIsNone(m.nearest_named(100.0, 0.0, max_dist=2.0))   # клик далеко от магнитов → None
        self.assertIsNotNone(m.nearest_named(0.3, 0.0, max_dist=2.0))  # у магнита → есть

    def test_pick_by_body_tag_returns_that_type(self):
        from client import scene3d
        g = {"centers": [[0, 0], [10, 0]],
             "named": [{"name": "QF1A", "xy": [0.0, 0.0], "kind": "quad", "dev": "und_QF1A"},
                       {"name": "DL1E", "xy": [0.3, 0.0], "kind": "dipole", "dev": "und_DL1E"}]}
        m = scene3d.SceneModel(g)
        # точка ближе к центру диполя, НО кликнули по телу квадруполя (tag) → вернётся квадруполь, не сосед
        self.assertEqual(m.nearest_named(0.25, 0.0, tag="und_QF1A")[1], "QF1A")
        self.assertEqual(m.nearest_named(0.25, 0.0)[1], "DL1E")        # без тега — ближайший центр (сосед)

    def test_build_geometry_screens_from_undulator_anchor(self):
        from client import scene3d
        if not os.path.exists(os.path.join(scene3d._project_root(), "workspace", "tables", "screens.csv")):
            self.skipTest("нет screens.csv")
        g = scene3d.build_geometry(_StubRingPlant())
        if "und_01" not in [s["gid"] for s in g["screens"]]:
            self.skipTest("в screens.csv нет источника und_01")
        sc = next(s for s in g["screens"] if s["gid"] == "und_01")
        self.assertEqual(sc["src"], [2.0, 1.0])                        # источник = живой якорь und_01 (mags[2])
        self.assertEqual(len(sc["center"]), 2)                         # центр экрана посчитан

    def test_start_plane_uses_beam_start_at_origin(self):
        from client import scene3d
        g = {"centers": [[5.0, 0.0], [6.0, 0.0]], "beam_start": {"pos": [0.0, 0.0], "dir": [1.0, 0.0]}}
        (cx, cy, _), (nx, ny, _) = scene3d.SceneModel(g).start_plane()
        self.assertEqual([cx, cy], [0.0, 0.0])                         # s=0 в начале координат, НЕ у centers[0]
        self.assertEqual([nx, ny], [1.0, 0.0])
        (cx2, cy2, _), _ = scene3d.SceneModel({"centers": [[5.0, 0.0], [6.0, 0.0]]}).start_plane()
        self.assertEqual([cx2, cy2], [5.0, 0.0])                       # без beam_start — у первого элемента

    def test_build_geometry_sensors_and_nearest(self):
        from client import scene3d
        g = scene3d.build_geometry(_StubRingPlant())
        self.assertEqual(len(g["sensors"]), 2)                        # два монитора → два сенсора
        self.assertEqual([s["idx"] for s in g["sensors"]], [0, 1])    # idx = порядок = индекс ebpm
        s = scene3d.SceneModel(g).nearest_sensor(3.0, 1.0)
        self.assertEqual(s["xy"], [3.0, 1.0])                         # ближайший сенсор к (3,1)

    @unittest.skipUnless(os.environ.get("SS_GUI_SMOKE") == "1",
                         "opt-in: НАСТОЯЩИЕ STL-тела (pyvista+visualize+STL) — SS_GUI_SMOKE=1")
    def test_stl_bodies_real_device(self):
        from client import scene3d
        try:
            import pyvista  # noqa: F401
            import visualize  # noqa: F401
        except Exception as e:
            self.skipTest(f"нет pyvista/visualize: {e}")
        hwp = os.path.join(scene3d._project_root(), "workspace", "tables", "hardware_library.csv")
        if not os.path.exists(hwp):
            self.skipTest("нет hardware_library.csv")
        mags = [{"dev": "und_UFC", "t": [0.0, 0.0, 0.0, 1.0, 2.0, 0.0], "kind": "quad"}]
        meshes = scene3d.build_magnet_stl_meshes(mags)   # STL-инстанс без построения кольца
        if not meshes:
            self.skipTest("STL для und_UFC недоступен в этой сборке")
        self.assertIn("und_UFC", meshes)
        mesh, kind = meshes["und_UFC"]
        self.assertEqual(kind, "quad")
        self.assertGreater(len(mesh.points), 0)          # меш реально собран из STL


class TestClientPanelFormat(unittest.TestCase):
    """Этап 2 (headless): чистое форматирование панели времени/скорости (без Qt-окна)."""

    def test_format_meta_strings(self):
        from client.cockpit import format_meta
        f = format_meta({"running": True, "model_time": 0.01234, "wall_time": 2.5,
                         "speed": 4.9e-3, "cycle": 7, "turn": 105, "alive": True})
        self.assertEqual(f["running"], "БЕЖИТ")
        self.assertIn("мс", f["model_time"])                     # 12.34 мс
        self.assertIn("×реал", f["speed"])
        self.assertEqual(f["cycle"], "7")
        self.assertEqual(f["alive"], "жив")
        self.assertEqual(format_meta({"running": False})["running"], "ПАУЗА")
        self.assertEqual(format_meta({"alive": False})["alive"], "ПОТЕРЯН")


class TestClientLayout(unittest.TestCase):
    """Этап 6 (headless): персистентность раскладки — make_layout + save/load/list round-trip (чистый JSON)."""

    def setUp(self):
        from client import layout as L
        self.L = L
        self.proj = "unittest_tmp_layout"
        import shutil
        shutil.rmtree(L.layout_dir(self.proj), ignore_errors=True)

    def tearDown(self):
        import shutil
        shutil.rmtree(self.L.layout_dir(self.proj), ignore_errors=True)

    def test_make_save_load_list(self):
        wins = [{"kind": "graph", "params": {"addr": "obs.ebpm.x:rms"}, "geom": [0, 0, 400, 300], "state": "normal"},
                {"kind": "info", "params": {"element": "SP01.QF1A"}, "geom": [10, 10, 300, 200], "state": "min"}]
        data = self.L.make_layout(wins, camera=[[1, 2, 3], [0, 0, 0], [0, 0, 1]], docks={"tree": True})
        self.assertEqual(data["version"], self.L.VERSION)
        p = self.L.save_layout(self.proj, "раскладка1", data)
        self.assertTrue(os.path.exists(p))
        got = self.L.load_layout(self.proj, "раскладка1")
        self.assertEqual(got["windows"][0]["params"]["addr"], "obs.ebpm.x:rms")
        self.assertEqual(got["camera"][0], [1, 2, 3])
        self.assertIn("раскладка1", self.L.list_layouts(self.proj))
        self.assertIsNone(self.L.load_layout(self.proj, "нет-такой"))
        self.L.save_last(self.proj, data)
        self.assertIsNotNone(self.L.load_last(self.proj))


@unittest.skipUnless(os.environ.get("SS_GUI_SMOKE") == "1",
                     "opt-in offscreen-smoke кокпита: SS_GUI_SMOKE=1 (PyVista+Qt offscreen)")
class TestClientCockpitOffscreen(unittest.TestCase):
    """Этапы 2–6 (opt-in): кокпит строится и обновляется offscreen без краха — панель отражает мету,
    сцена принимает снимок (зазор) и пикинг открывает инфо, окна-подписчики подписываются/обновляются,
    дерево строится, раскладка сериализуется. Гейт флагом; при ошибке Qt/VTK — skip."""

    def _env(self):
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from PySide6.QtWidgets import QApplication
        import pyvista as pv
        pv.OFF_SCREEN = True
        return QApplication.instance() or QApplication([])

    def test_cockpit_offscreen_tick_and_gap(self):
        try:
            app = self._env()
            from client import InProcLink, Snapshot
            from client.cockpit import Cockpit
            from client.scene3d import make_ring_scene
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, poll_ms=10)
        try:
            scene = make_ring_scene(gap_scale=100.0)
            scene.set_geometry({"centers": [[0, 0], [1, 0], [2, 1]], "mon_xy": [[1, 0]],
                                "named": [{"name": "SP01.QF1A", "xy": [2.0, 1.0], "kind": "quad"}],
                                "undulators": {"u": {"xy": [1.0, 0.0], "gap_mm": 5.0}}})
            win.bind_scene(scene)
            win.control.set_undulators(["u"])
            link.publish(Snapshot(values={"undulator.u.gap": 12.0},
                                  meta={"running": True, "model_time": 3e-3, "wall_time": 1.0,
                                        "speed": 1.2e-3, "cycle": 5, "seq": 1, "turn": 75, "alive": True}))
            win._tick()
            app.processEvents()
            self.assertEqual(win.control.lbl["cycle"].text(), "5")
            self.assertIn("мс", win.control.lbl["model_time"].text())
        finally:
            win.close()

    def test_bend_line_screens_drawn_as_one_actor(self):
        """Экраны у линеек диполей — одним актёром «bscreens» (+ линии визирования «blos»), не по актёру на экран."""
        try:
            self._env()
            from client.scene3d import make_ring_scene
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        scene = make_ring_scene()
        try:
            scr = [{"sid": f"SCR_DP{k}", "kind": "bend", "src": [float(k), 0.0], "center": [float(k), 20.0],
                    "axis": [0.0, 1.0], "w": 0.1, "h": 0.03} for k in range(3)]
            scene.set_geometry({"centers": [[0, 0], [1, 0], [2, 1]], "screens": scr})
            self.assertIn("bscreens", scene.actors)
            self.assertIn("blos", scene.actors)
            self.assertFalse([n for n in scene.actors if str(n).startswith("screen_SCR_DP")])
            self.assertEqual(scene.model.nearest_screen(2.1, 19.0, kind="bend")["sid"], "SCR_DP2")
        finally:
            scene.close()

    def test_windows_wait_while_camera_moves(self):
        """Пока крутят камеру (кнопка мыши, колесо, перелёт), окна графиков не перерисовываются — снимок
        откладывается и отдаётся после; сцена и панель обновляются всегда. Кривые рядов — без сглаживания."""
        try:
            self._env()
            from time import perf_counter
            from client import InProcLink, Snapshot
            from client.cockpit import Cockpit
            from client.scene3d import make_ring_scene
            from client.windows import GraphWindow
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, poll_ms=10)
        try:
            scene = make_ring_scene()
            scene.set_geometry({"centers": [[0, 0], [1, 0], [2, 1]], "undulators": {"u": {"xy": [1.0, 0.0], "gap_mm": 5.0}}})
            win.bind_scene(scene)

            class Spy:
                def __init__(self):
                    self.got = []

                def on_snapshot(self, snap):
                    self.got.append(int(snap.meta["seq"]))
            spy = Spy(); win._windows.append(spy)
            meta = lambda s: {"running": True, "model_time": 1e-3 * s, "wall_time": 1.0, "speed": 1e-3, "cycle": s,
                              "seq": s, "turn": s, "alive": True}
            self.assertFalse(scene.is_interacting())
            scene._dragging = True                                  # ЛКМ зажата — вращение
            link.publish(Snapshot(values={"undulator.u.gap": 8.0}, meta=meta(1))); win._tick()
            link.publish(Snapshot(values={"undulator.u.gap": 9.0}, meta=meta(2))); win._tick()
            self.assertEqual(spy.got, [])                           # окна ждут
            self.assertEqual(scene._target_gap["u"], 9.0)           # сцена получила снимок
            self.assertEqual(win.control.lbl["cycle"].text(), "2")  # панель тоже
            scene._dragging = False
            win._tick()
            self.assertEqual(spy.got, [2])                          # после — один, последний снимок
            scene._wheel_t = perf_counter()                         # колесо — тоже осмотр
            self.assertTrue(scene.is_interacting())
            scene._wheel_t = -1e9; scene._mmb_down = True           # СКМ (панорама)
            self.assertTrue(scene.is_interacting())
            scene._mmb_down = False
            link.publish(Snapshot(values={}, meta=meta(3))); win._tick()
            self.assertEqual(spy.got, [2, 3])                       # без осмотра — сразу
            win._windows.remove(spy)
            gw = GraphWindow("actual.quad.QF.dk")
            self.assertFalse(gw.curve.opts["antialias"])            # ряды без сглаживания (×2 дешевле перерисовка)
        finally:
            win.close()

    def test_scene_camera_navigation_keys(self):
        """Ctrl/Shift+←/→: камера перелетает на соседний элемент / соседний того же типа (без фокуса — подсказка)."""
        try:
            app = self._env()
            from client import InProcLink
            from client.cockpit import Cockpit
            from client.scene3d import make_ring_scene
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, poll_ms=10)
        try:
            scene = make_ring_scene()
            scene.set_geometry(dict(TestClientScene._nav_bundle(), centers=[[0, 0], [4, 0], [8, 0]]))
            win.bind_scene(scene)
            scene.fly_steps = 0                                     # без анимации — проверяем итог
            got = []
            scene.focused.connect(got.append)
            scene._center_on(100.0, 100.0)                          # камера не на элементе
            self.assertIsNone(scene.nav_step(+1))
            self.assertIn("не на элементе", got[-1])
            scene._center_on(0.0, 0.0); scene._set_focus_at(0.0, 0.0)   # «клик» по Q1
            win._nav(+1, False)
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (1.0, 0.0))   # → C1
            win._nav(+1, True); win._nav(+1, True)                  # C1 → C1 (единственный COR)
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (1.0, 0.0))
            win._nav(+1, False)                                     # → u1 (центр гребёнки)
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (3.0, 0.0))
            win._nav(+1, True)                                      # Shift+→ — следующий ондулятор
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (7.0, 0.0))
            self.assertIn("ондулятор u2", win.statusBar().currentMessage())
            win._nav(-1, False)                                     # Ctrl+← — навстречу пучку
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (6.0, 0.0))
            scene.fly_steps = 4                                     # с анимацией — доезжает до цели
            win._nav(+1, False)                                     # Ctrl+→ B1 → u2 плавным перелётом
            import time as _t
            end = _t.time() + 2.0
            while scene._fly_timer is not None and _t.time() < end:
                app.processEvents(); _t.sleep(0.01)
            self.assertEqual(tuple(round(v, 9) for v in scene.camera.focal_point[:2]), (7.0, 0.0))
            self.assertEqual(len(win._nav_keys), 6)
        finally:
            win.close()

    def test_subscriber_windows_and_refcount(self):
        try:
            app = self._env()
            from client import InProcLink, Snapshot
            from client.cockpit import Cockpit
        except Exception as e:
            self.skipTest(f"нет Qt offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, poll_ms=10)
        try:
            win._open_graph("obs.ebpm.x:rms")               # окно подписалось на адрес
            link.pull()
            self.assertIn("obs.ebpm.x:rms", link.subscriptions())
            link.publish(Snapshot(values={"obs.ebpm.x:rms": 1.5e-6},
                                  meta={"seq": 1, "model_time": 1e-3}))
            win._tick()
            app.processEvents()
            self.assertEqual(len(win._windows), 1)
            win._windows[0]._toolwin.close()                # крестик окна снимает подписку
            app.processEvents()
            link.pull()
            self.assertEqual(len(win._windows), 0)
            self.assertNotIn("obs.ebpm.x:rms", link.subscriptions())
        finally:
            win.close()

    def test_point_and_screen_windows_construct(self):
        try:
            app = self._env()
            from client import Snapshot
            from client.windows import PointWindow, ScreenWindow, make_window
        except Exception as e:
            self.skipTest(f"нет Qt/pyqtgraph offscreen: {e}")
        pw = PointWindow("obs.ebpm.x:0", "obs.ebpm.y:0", trail=50)
        self.assertEqual(set(pw.addresses()), {"obs.ebpm.x:0", "obs.ebpm.y:0"})
        pw.on_snapshot(Snapshot(values={"obs.ebpm.x:0": 1e-6, "obs.ebpm.y:0": -2e-6}, meta={}))
        sw = ScreenWindow("SCR1", mesh=61)
        self.assertEqual(sw.addresses(), ["screen.SCR1"])
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {"heatmap": [[1.0, 2.0], [3.0, 4.0]],
                                                         "e1_series": [[0, 5.0, 7.0]], "gap": 5.0, "e1": 7.0}},
                                meta={}))
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {"error": "нет модулей"}}, meta={}))  # ошибка — не крах
        round_trip = make_window(pw.to_dict())                 # сериализация окна → пересоздание
        self.assertIsInstance(round_trip, PointWindow)
        self.assertEqual(round_trip.trail, 50)
        from client.windows import ScreenInfoWindow            # справка экрана: кнопка → request_open("screen", sid)
        si = ScreenInfoWindow("SCR_U01")
        self.assertEqual(si.to_dict()["kind"], "screen_info")
        got = {}
        si.request_open.connect(lambda k, a: got.update(kind=k, addr=a))
        si.request_open.emit("screen", "SCR_U01")
        self.assertEqual(got, {"kind": "screen", "addr": "SCR_U01"})
        self.assertIsInstance(make_window(si.to_dict()), ScreenInfoWindow)
        from client.windows import UndulatorInfoWindow                 # окно ондулятора: живой зазор + ручка
        uw = UndulatorInfoWindow("und_01")
        self.assertEqual(uw.addresses(), ["undulator.und_01.gap"])
        uw.on_snapshot(Snapshot(values={"undulator.und_01.gap": 7.5}, meta={}))
        self.assertIn("7.5", uw.val.text())
        self.assertIsInstance(make_window(uw.to_dict()), UndulatorInfoWindow)
        from client.windows import SensorInfoWindow                   # окно сенсора: параметры + живое x/y
        sn = SensorInfoWindow(5, name="SP01.BPM.6", pos=[1.0, 2.0, 0.0], sid="42")
        self.assertEqual(sn.addresses()[:2], ["obs.ebpm.x:5", "obs.ebpm.y:5"])     # живые показания x/y
        self.assertIn("element.SP01.BPM.6.yaw", sn.addresses())                    # + углы ориентации прибора
        sn.on_snapshot(Snapshot(values={"obs.ebpm.x:5": 1e-6, "obs.ebpm.y:5": -2e-6}, meta={}))
        self.assertIn("мкм", sn.vx.text())
        self.assertIsInstance(make_window(sn.to_dict()), SensorInfoWindow)
        app.processEvents()

    def test_orbit_modes_window(self):
        """Стандартный график «η-узор и моды FOFB»: окно рисует данные diag.orbit_modes, ручки δ/корректор
        перерисовывают, ошибка — текстом, раскладка пересоздаёт окно с параметрами."""
        try:
            app = self._env()
            from client import Snapshot
            from client.windows import OrbitModesWindow, DIAG_WINDOWS, make_window
        except Exception as e:
            self.skipTest(f"нет Qt/pyqtgraph offscreen: {e}")
        self.assertIs(DIAG_WINDOWS["orbit_modes"], OrbitModesWindow)
        w = OrbitModesWindow(delta_ppm=50.0, corr=1)
        self.assertEqual(w.addresses(), ["diag.orbit_modes"])
        d = {"eta_bpm": [0.01, 0.07, 0.002], "cols": [[1.0, -2.0, 0.5], [3.0, 0.0, -1.0]],
             "corr_ids": ["corr.000.hkick", "corr.001.hkick"], "sv": [1.0, 0.004], "rcond": 1e-2,
             "n_kept": 1, "n_modes": 2, "eta_reach": 0.19, "eta_reach_full": 0.75}
        w.on_snapshot(Snapshot(values={"diag.orbit_modes": d}, meta={}))
        x, y = w.c_eta.getData()
        np.testing.assert_allclose(y, np.array(d["eta_bpm"]) * 50.0)          # η·δ, мкм
        np.testing.assert_allclose(w.c_bet.getData()[1], d["cols"][1])       # волна корректора 1 (×1 мкрад, мкм)
        self.assertIn("1</b> из 2", w.info.text())
        self.assertIn("19%", w.info.text())
        w.sp_corr.setValue(0)
        np.testing.assert_allclose(w.c_bet.getData()[1], d["cols"][0])
        w.sp_delta.setValue(200.0)
        np.testing.assert_allclose(w.c_eta.getData()[1], np.array(d["eta_bpm"]) * 200.0)
        w.on_snapshot(Snapshot(values={"diag.orbit_modes": {"error": "нет FOFB"}}, meta={}))
        self.assertIn("нет FOFB", w.info.text())
        rt = make_window(w.to_dict())
        self.assertIsInstance(rt, OrbitModesWindow)
        self.assertEqual((rt.delta_ppm, rt.corr), (200.0, 0))
        app.processEvents()

    def test_info_redesign_object_specific_links_and_undulators(self):
        """Редизайн: объекто-специфичное InfoWindow, значения-гиперссылки, цветная карта+спектр+точка,
        конфиг ондуляторов (один включён), перенос точки через Link.set_screen."""
        try:
            self._env()
            from client import Snapshot, InProcLink
            from client.link import SET_SCREEN
            from client.windows import InfoWindow, ValueLink, ScreenWindow, make_window
        except Exception as e:
            self.skipTest(f"нет Qt/pyqtgraph offscreen: {e}")

        iq = InfoWindow("Q1", kind="quad")                          # квадруполь: только k1 (не k2/k3/theta0)
        self.assertIn("element.Q1.k1", iq.addresses())
        self.assertNotIn("element.Q1.k2", iq.addresses())
        self.assertNotIn("element.Q1.theta0", iq.addresses())
        ic = InfoWindow("C1", kind="corrector", aids=["corr.000.hkick", "corr.000.vkick"])
        self.assertIn("act.corr.000.hkick", ic.addresses())         # корректор: живые H/V кики (канон. act.*, шаг 4)
        self.assertIn("act.corr.000.vkick", ic.addresses())
        rt = make_window(ic.to_dict())                              # сериализация сохраняет kind/aids
        self.assertEqual((rt.kind, rt.aids), ("corrector", ["corr.000.hkick", "corr.000.vkick"]))

        got = {}                                                    # значение-ссылка эмитит request_open
        vl = ValueLink("element.Q1.k1", emit=lambda k, a: got.update(k=k, a=a))
        vl._on_link("graph::element.Q1.k1")
        self.assertEqual(got, {"k": "graph", "a": "element.Q1.k1"})
        vl._on_link("point::obs.ebpm.x:5|obs.ebpm.y:5")             # point-адрес 'ax|ay' цел (разделитель '::')
        self.assertEqual(got, {"k": "point", "a": "obs.ebpm.x:5|obs.ebpm.y:5"})

        sw = ScreenWindow("SCR1", mesh=61)                          # цветная карта + extent + настоящий спектр
        self.assertEqual(sw._point, (0.0, 0.0))
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {
            "heatmap": [[1.0, 2.0], [3.0, 4.0]], "extent": [-1, 1, -1, 1],
            "spectrum": {"E_keV": [5.0, 7.0, 9.0], "flux": [1e3, 1e4, 1e2]}, "gap": 5.0, "e1": 7.0}}, meta={}))
        self.assertEqual(sw.params()["point"], [0.0, 0.0])
        self.assertIsInstance(make_window(sw.to_dict()), ScreenWindow)

        from client.cockpit import init_undulators                  # один ондулятор сведён, остальные разведены
        from control import UndulatorConfig, ExperimentState, Schema
        st = ExperimentState(Schema(["a"]))
        for i in (1, 2, 3):
            st.undulators[f"und_{i:02d}"] = UndulatorConfig(5.0, name=f"und_{i:02d}")
        init_undulators(st, 500.0)
        self.assertAlmostEqual(st.undulators["und_01"].gap_mm, 5.0)
        self.assertAlmostEqual(st.undulators["und_02"].gap_mm, 50.0)
        self.assertAlmostEqual(st.undulators["und_03"].gap_max, 50.0)

        lk = InProcLink()                                           # set_screen доносит точку спектра до движка
        lk.set_screen("screen.S", area=None, res=121, point=(0.3, -0.2))
        ss = [m for m in lk.pull() if m.kind == SET_SCREEN]
        self.assertTrue(ss and ss[0].value.get("point") == (0.3, -0.2))

        sw.bind_link(lk); lk.pull()                                 # переключатель «величина»: поток ↔ мощность
        self.assertEqual(sw.quantity, "flux")
        sw.qcmb.setCurrentIndex(1)
        ss = [m for m in lk.pull() if m.kind == SET_SCREEN]
        self.assertEqual(ss[-1].value.get("quantity"), "power")
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {           # кадр мощности: пик и перехват в строке
            "heatmap": [[0.0, 640.0], [10.0, 5.0]], "extent": [-25, 25, -25, 25], "quantity": "power",
            "gap": 5.0, "e1": 7.0, "peak": 640.0, "total": 3370.0, "und_W": 2978.0, "bend_W": 392.0, "n_bend": 2}},
            meta={}))
        self.assertIn("Вт/мм²", sw.info.text())
        rt = make_window(sw.to_dict())
        self.assertEqual(rt.quantity, "power")                      # раскладка помнит величину
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {"heatmap": [[9.0]], "quantity": "flux"}}, meta={}))
        self.assertIn("Вт/мм²", sw.info.text())                     # чужой кадр (до пересчёта) не рисуется

        self.assertEqual(sw._heat_value(25.0, 25.0), 5.0)           # значение карты под курсором ([u, v] = (x, y))
        self.assertIsNone(sw._heat_value(99.0, 0.0))
        sw.chk_log.setChecked(True)                                 # лог. шкала: цвет по log₁₀, данные — как были
        self.assertAlmostEqual(float(sw.item.image.max()), np.log10(640.0), places=9)
        self.assertEqual(sw._heat_value(25.0, 25.0), 5.0)
        E = np.linspace(1.0, 20.0, 191); F = 1e10 * (1 + 100 * np.exp(-((E - 7.0) / 0.2) ** 2))
        sw.on_snapshot(Snapshot(values={"screen.SCR1": {"heatmap": [[1.0, 2.0], [3.0, 4.0]], "quantity": "power",
                                                        "extent": [-1, 1, -1, 1], "spectrum": {"E_keV": E.tolist(),
                                                                                             "flux": F.tolist()}}},
                                meta={}))
        mk = sw._add_mark(7.05)                                     # «к пику»: встаёт на максимум 7.0 кэВ
        self.assertAlmostEqual(mk["E"], 7.0, places=6)
        self.assertEqual(sw.marks_tbl.rowCount(), 1)
        sw.marks_tbl.item(0, 0).setText("E₁")                       # имя правится в таблице → легенда
        sw.marks_tbl.item(0, 1).setText("12,5")                     # энергия правится (запятая тоже) → вертикаль едет
        self.assertEqual(mk["name"], "E₁")
        self.assertAlmostEqual(float(mk["line"].value()), 12.5)
        sw._add_mark(15.0, snap=False)
        rt = make_window(sw.to_dict())                              # раскладка помнит лог-шкалу и метки
        self.assertTrue(rt.log)
        self.assertEqual([(m["name"], round(m["E"], 6)) for m in rt._marks], [("E₁", 12.5), ("метка 2", 15.0)])
        sw.marks_tbl.selectRow(0); sw._delete_selected()
        self.assertEqual([m["name"] for m in sw._marks], ["метка 2"])
        sw._clear_marks()
        self.assertEqual(sw.marks_tbl.rowCount(), 0)

        import re                                                   # подписи осей без латинских SI-приставок
        from client.windows import _plot
        for units, y, log, want in (("1/м²", [-1.2e-3, 0.0], False, "(×10⁻³ 1/м²)"), ("м", [-2e-5, 1.5e-5], False, "(мкм)"),
                                    ("ф/с/мрад²/0.1%BW", [1e10, 3e17], True, "(ф/с/мрад²/0.1%BW)"),
                                    ("", [1e-6, 3e-6], False, "(×10⁻⁶)")):
            pw = _plot(); pw.resize(400, 300); pw.show()
            pw.setLabel("left", "q", units=units or None)
            if log:
                pw.setLogMode(x=False, y=True)
            pw.plot([0.0, 1.0], y)
            self._env().processEvents()
            self.assertIn(want, re.sub("<[^>]+>", "", pw.getAxis("left").labelString()))

        from client.windows import CurveHover                       # значение под курсором на кривой
        import pyqtgraph as pg
        pw = pg.PlotWidget(); cv = pw.plot([1.0, 2.0, 3.0], [10.0, 20.0, 30.0]); pw.resize(400, 300); pw.show()
        hv = CurveHover(pw, [(cv, "ряд", lambda x: f"x={x:g}", lambda y: f"y={y:g}")])
        self._env().processEvents()
        hv._move(pw.getPlotItem().vb.mapViewToScene(pg.Point(2.1, 12.0)))   # ближайшая по x — точка (2, 20)
        self.assertIn("y=20", hv.text.toHtml())

    def test_orientation_noise_and_panel_optics(self):
        """Доп. инфо: углы ориентации прибора (element.*.yaw/pitch/roll), закон+σ шума датчика, и живые
        тюн/хроматизм/adts на панели управления (клик по значению → график)."""
        try:
            self._env()
            from client import addresses as A
            from client.windows import info_fields, InfoWindow, SensorInfoWindow
            from client.cockpit import _optics_specs, ControlPanel
            from PySide6.QtWidgets import QLabel
        except Exception as e:
            self.skipTest(f"нет Qt offscreen: {e}")

        # углы ориентации резолвятся и у любого прибора есть в info_fields (в мрад, scale 1e3)
        class _MV:
            def orientation(self, dev): return {"yaw": -6.283, "pitch": 0.0, "roll": 1.5e-3}
            def position(self, dev): return (0.0, 0.0, 0.0)
            def field(self, dev): return {"theta0": 0, "k1": 0, "k2": 0, "k3": 0}
        view = type("V", (), {"machine": _MV(), "actual": None})()
        self.assertAlmostEqual(A.resolve(view, "element.Q.roll"), 1.5e-3)
        self.assertEqual(A.meta("element.Q.yaw").units, "рад")
        ori = [r for r in info_fields("quad", "Q1", None) if r[0].endswith(".yaw")]
        self.assertTrue(ori and len(ori[0]) > 3 and ori[0][3] == 1e3)          # 4-кортеж со scale рад→мрад
        self.assertIn("element.Q1.yaw", InfoWindow("Q1", kind="quad").addresses())
        self.assertTrue(any(r[0].endswith(".yaw") for r in                      # корректор — тоже прибор
                            info_fields("corrector", "C1", ["corr.000.hkick", "corr.000.vkick"])))

        # сенсор: ориентация + шум (закон + σ в мкм)
        sw = SensorInfoWindow(0, name="BPM1", pos=[0, 0, 0], sid="B01", meas_noise=1e-7)
        self.assertIn("element.BPM1.yaw", sw.addresses())
        texts = " ".join(w.text() for w in sw.findChildren(QLabel))
        self.assertIn("гауссов белый", texts)                                   # закон шума
        self.assertIn("0.100 мкм", texts)                                       # σ = 1e-7 м = 0.1 мкм
        self.assertEqual(SensorInfoWindow.from_params(sw.params()).meas_noise, 1e-7)   # раскладка хранит σ

        # панель: тюн/хроматизм/adts по компонентам, только если наблюдаемая есть
        st = type("S", (), {"observables": {"ebpm.x": 0, "tune": [0.3, 0.2],
                                            "chroma": [5., 6.], "adts": [1., 2., 3.]}})()
        specs = _optics_specs(st)
        self.assertEqual([s[0] for s in specs],
                         ["obs.tune:0", "obs.tune:1", "obs.chroma:0", "obs.chroma:1",
                          "obs.adts:0", "obs.adts:1", "obs.adts:2"])
        self.assertEqual(_optics_specs(type("S0", (), {"observables": {"ebpm.x": 0}})()), [])   # orbit-only пусто
        got = []
        cp = ControlPanel(link=type("L", (), {})())
        cp.set_optics(specs, lambda k, a: got.append((k, a)))
        self.assertEqual(len(cp._opt_links), 7)
        cp.update_optics({"obs.tune:0": 0.3141})
        self.assertIn("0.3141", cp._opt_links["obs.tune:0"].text())             # живое значение
        cp._opt_links["obs.tune:0"]._on_link("graph::obs.tune:0")
        self.assertEqual(got, [("graph", "obs.tune:0")])                        # клик → график (мониторинг)

    def test_pick_opens_info_and_tree_builds(self):
        try:
            app = self._env()
            from client import InProcLink
            from client.cockpit import Cockpit
            from client.scene3d import make_ring_scene
            from control import ExperimentState, Schema, UndulatorConfig
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, poll_ms=10)
        try:
            scene = make_ring_scene()
            scene.set_geometry({"centers": [[0, 0], [1, 0], [2, 1]],
                                "named": [{"name": "SP01.QF1A", "xy": [2.0, 1.0], "kind": "quad"}],
                                "undulators": {}})
            win.bind_scene(scene)
            scene.pick_nearest(2.0, 1.0)                    # пикинг → picked → _open_info
            app.processEvents()
            self.assertTrue(any(w.KIND == "info" for w in win._windows))
            # дерево
            st = ExperimentState(Schema(["a"]))
            st.observables["ebpm.x"] = np.zeros(3)
            st.undulators["und_01"] = UndulatorConfig(7.0, name="und_01")
            win.build_tree(st, {}, ["SCR1"])
            self.assertGreaterEqual(win.tree.topLevelItemCount(), 5)
            win.build_tree(st, {}, [("SCR_und_01", "undulator"), ("SCR_DP7", "bend")])   # экраны — по видам
            scr = next(win.tree.topLevelItem(i) for i in range(win.tree.topLevelItemCount())
                       if win.tree.topLevelItem(i).text(0) == "экраны")
            self.assertEqual(sorted(scr.child(i).text(0) for i in range(scr.childCount())),
                             ["у линеек диполей (модуль №3)", "у ондуляторов"])
        finally:
            win.close()

    def test_rmb_zoom_not_sticky_gates_hover(self):
        """Регрессия на залипание ПКМ-зума: ПКМ ловится на уровне Qt (`mousePressEvent`/`mouseReleaseEvent`),
        а НЕ через VTK-наблюдатели (их наличие на Right-Button ломало завершение dolly). Пока ПКМ зажата
        (`_rmb_down`), hover НЕ пикает. `_rmb_down` переключается синтезированными Qt-событиями ПКМ."""
        try:
            self._env()
            from client.scene3d import make_ring_scene
        except Exception as e:
            self.skipTest(f"нет Qt/PyVista offscreen: {e}")
        scene = make_ring_scene()
        try:
            scene.set_geometry({"centers": [[0, 0], [1, 0], [2, 1]],
                                "named": [{"name": "SP01.QF1A", "xy": [2.0, 1.0], "kind": "quad"}],
                                "undulators": {}})
            self.assertFalse(scene._rmb_down)                  # _enable_picking инициализировал флаг
            # ПКМ — на уровне Qt-виджета, НЕ VTK-наблюдатель (структурная гарантия от залипания):
            self.assertIn("mousePressEvent", type(scene).__dict__)
            self.assertIn("mouseReleaseEvent", type(scene).__dict__)
            # во время ПКМ-зума hover ДОЛЖЕН молчать (не пикать — иначе мешает dolly):
            scene._rmb_down = True
            scene._hover_name = "SENTINEL"
            scene._on_move()
            self.assertEqual(scene._hover_name, "SENTINEL")    # ранний выход, пик не вызван
            scene._rmb_down = False
            # синтезированные Qt-события ПКМ переключают флаг (super() отрабатывает зум VTK):
            try:
                from PySide6.QtGui import QMouseEvent
                from PySide6.QtCore import QPointF, QEvent, Qt
                def _mev(t):
                    return QMouseEvent(t, QPointF(10, 10), QPointF(10, 10), Qt.MouseButton.RightButton,
                                       Qt.MouseButton.RightButton, Qt.KeyboardModifier.NoModifier)
                scene.mousePressEvent(_mev(QEvent.Type.MouseButtonPress))
                self.assertTrue(scene._rmb_down)               # зажали ПКМ
                scene.mouseReleaseEvent(_mev(QEvent.Type.MouseButtonRelease))
                self.assertFalse(scene._rmb_down)              # отпустили — флаг снят
            except ImportError:
                pass                                           # синтез события недоступен — структурных проверок хватает
        finally:
            scene.close()

    def test_loops_panel_toggles_and_master(self):
        """Панель контуров САУ: список по именам файлов логики, переключатели вкл/выкл, мастер на все сразу,
        синхрон с живым состоянием из снимка."""
        try:
            self._env()
            from client.cockpit import ControlPanel, _loops_info
            from PySide6.QtCore import Qt
        except Exception as e:
            self.skipTest(f"нет Qt offscreen: {e}")

        class _RecLink:                                        # пишет команды, не трогая движок
            def __init__(self): self.cmds = []
            def command(self, a, v): self.cmds.append((a, v))
            def control(self, *a): pass

        # _loops_info достаёт имя файла логики контроллера из Scheduler.loops
        from control.controllers.integral_svd import IntegralSVD

        class _Spec:
            def __init__(self, name, en): self.controller = IntegralSVD(name, sensor="ebpm.x",
                                                                        actuator_key="hkick"); self.enabled = en; self.rate_hz = 2000.0
        info = _loops_info({"fofb_x": _Spec("fofb_x", True), "fofb_y": _Spec("fofb_y", False)})
        self.assertEqual({d["name"] for d in info}, {"fofb_x", "fofb_y"})
        self.assertTrue(all(d["file"].endswith(".py") and "integral_svd" in d["file"] for d in info))

        link = _RecLink()
        cp = ControlPanel(link)
        self.assertFalse(cp.box_loops.isVisible())             # пусто → панель скрыта
        cp.set_loops(info)
        self.assertEqual(set(cp._loop_chk), {"fofb_x", "fofb_y"})
        self.assertTrue(cp._loop_chk["fofb_x"].isChecked())
        self.assertFalse(cp._loop_chk["fofb_y"].isChecked())
        self.assertEqual(cp.chk_all.checkState(), Qt.PartiallyChecked)   # один вкл, один выкл

        cp._toggle_loop("fofb_y", True)                        # переключатель одного контура → команда
        self.assertIn(("loop.fofb_y.enabled", True), link.cmds)

        link.cmds.clear()                                      # мастер: сейчас все вкл → выключить все
        cp._loop_chk["fofb_y"].setChecked(True); cp._refresh_all_state()
        self.assertEqual(cp.chk_all.checkState(), Qt.Checked)
        cp._toggle_all()
        self.assertEqual(sorted(link.cmds), [("loop.fofb_x.enabled", False), ("loop.fofb_y.enabled", False)])
        self.assertEqual(cp.chk_all.checkState(), Qt.Unchecked)
        self.assertFalse(any(c.isChecked() for c in cp._loop_chk.values()))

        # синхрон с движком из снимка (без обратной команды)
        link.cmds.clear()
        cp.update_loops({"loop.fofb_x.enabled": True, "loop.fofb_y.enabled": False})
        self.assertTrue(cp._loop_chk["fofb_x"].isChecked())
        self.assertFalse(cp._loop_chk["fofb_y"].isChecked())
        self.assertEqual(link.cmds, [])                        # синхрон не шлёт команд

    def test_layout_roundtrip_reopens_windows(self):
        try:
            app = self._env()
            from client import InProcLink
            from client.cockpit import Cockpit
        except Exception as e:
            self.skipTest(f"нет Qt offscreen: {e}")
        link = InProcLink()
        win = Cockpit(link, project="unittest_tmp_cockpit", poll_ms=10)
        try:
            win._open_graph("meta.speed")
            win._open_graph("obs.ebpm.y:rms")
            app.processEvents()
            data = win.collect_layout()
            self.assertEqual(len(data["windows"]), 2)
            win.restore_layout(data)                        # закрывает и пересоздаёт окна
            app.processEvents()
            self.assertEqual(len(win._windows), 2)
            kinds = {w.KIND for w in win._windows}
            self.assertEqual(kinds, {"graph"})
        finally:
            win.close()
            import shutil
            from client import layout as L
            shutil.rmtree(L.layout_dir("unittest_tmp_cockpit"), ignore_errors=True)


@unittest.skipUnless(os.environ.get("SS_ENGINE_SMOKE") == "1",
                     "opt-in smoke на РЕАЛЬНОМ TbtPlant: SS_ENGINE_SMOKE=1 (строит 6D-кольцо, ~20 с)")
class TestClientEngineSmoke(unittest.TestCase):
    """Этап 1 (opt-in): build_session + EngineThread против НАСТОЯЩЕГО TbtPlant (ring_cpp). Поток шагает,
    модельное время/скорость идут, наблюдаемая и зазор резолвятся. Гейт флагом, чтобы не утяжелять набор."""

    def test_real_plant_runs(self):
        import time as _t
        try:
            import ring_cpp  # noqa: F401
        except Exception as e:
            self.skipTest(f"нет ring_cpp: {e}")
        from client import InProcLink
        from client.engine import build_session, EngineThread
        from client.scene3d import build_geometry

        plant, state, sch = build_session(settle_turns=30, fofb=True, gap_mm=6.0, steps_per_cycle=2)
        geom = build_geometry(plant, state)                      # геометрия сцены на реальном кольце
        self.assertGreater(len(geom["centers"]), 100)
        self.assertTrue(geom["undulators"])                      # якоря ондуляторов найдены
        link = InProcLink()
        eng = EngineThread(link, sch, steps_per_cycle=2)
        gid = next(iter(state.undulators), None)
        link.subscribe("obs.ebpm.x:rms")
        link.subscribe("meta.speed")
        link.subscribe("meta.cycle")
        if gid:
            link.subscribe(f"undulator.{gid}.gap_target")
            link.command(f"undulator.{gid}.gap", 8.0)
        link.control("run")
        eng.start()
        snap = None
        try:
            end = _t.time() + 90
            while _t.time() < end:
                snap = link.poll_snapshot()
                if snap is not None and snap.meta.get("cycle", 0) >= 2:
                    break
                _t.sleep(0.05)
        finally:
            eng.stop()
            eng.join(timeout=30)
            plant.close()
        self.assertIsNotNone(snap)
        self.assertGreaterEqual(snap.meta["cycle"], 2)
        self.assertGreater(snap.meta["speed"], 0.0)
        self.assertIsInstance(snap.values["obs.ebpm.x:rms"], float)
        if gid:
            self.assertAlmostEqual(snap.values[f"undulator.{gid}.gap_target"], 8.0)


class TestScreenHeatMap(unittest.TestCase):
    """Тепловая карта экранов: точная угловая мощность ондулятора (Льенар по периоду) и веера диполей."""

    @classmethod
    def setUpClass(cls):
        und = os.path.join(ROOT, "undulators")
        if und not in sys.path:
            sys.path.insert(0, und)
        import undulator_sr as sr
        cls.sr = sr

    def test_lienard_peak_and_total_match_closed_forms(self):
        """Таблица dP/dΩ: пик на оси = формула Кима, интеграл = 0.633·E²B0²·(Nλ_u)·I (без подгонки)."""
        sr = self.sr
        for gap in (5.0, 12.0):
            B0 = sr.B0_of_gap(gap); K = sr.K_of_B0(B0)
            tab = sr.angular_power_table(K)
            self.assertAlmostEqual(tab["peak"] / 1e9 / sr.power_density_onaxis(B0, K), 1.0, delta=5e-3)
            P_nlu = 0.633 * sr.P.E ** 2 * B0 ** 2 * sr.P.N * sr.P.lam_u * sr.P.I * 1e3
            self.assertAlmostEqual(tab["total_W"] / P_nlu, 1.0, delta=5e-3)

    def test_lienard_shape_elongated_in_oscillation_plane(self):
        """При K≈1.85 пятно шире в плоскости колебаний (±K/γ), чем поперёк (~1/γ) — не круглый гаусс."""
        sr = self.sr
        K = sr.K_of_B0(sr.B0_of_gap(5.0))
        tab = sr.angular_power_cached(K)
        F, tx, ty = tab["F"], tab["tx"], tab["ty"]

        def fwhm(a, pr):
            i = np.where(pr >= pr.max() / 2)[0]
            return a[i[-1]] - a[i[0]]
        wx = fwhm(tx, F[:, ty.size // 2]); wy = fwhm(ty, F[tx.size // 2, :])
        self.assertGreater(wx / wy, 2.0)
        self.assertAlmostEqual(float(sr.angular_power(tab, 0.0, 0.0)), tab["peak"], delta=1e-6 * tab["peak"])
        self.assertEqual(float(sr.angular_power(tab, 1.0, 0.0)), 0.0)          # вне таблицы — ноль

    def test_interference_summed_over_energy_equals_lienard(self):
        """Парсеваль: интерференционная модель (этап B), просуммированная по энергиям, на оси = Льенар.
        Поэтому интерференция в тепловую карту не входит."""
        sr = self.sr
        import undulator_field as uf
        B0 = sr.B0_of_gap(5.0); K = sr.K_of_B0(B0)
        tr = uf.trajectory(K, sr.P.lam_u, sr.P.N, sr.gamma(), ppp=128)   # ppp=64 даёт алиасинг выше ~150 кэВ
        scale, _ = uf.calibrate(tr)
        Es = np.arange(0.2, 230.0, 0.025)                           # на оси спектр жёсткий: до 120 кэВ лишь ~95%
        d = uf._dirs(np.array([0.0]), np.array([0.0]))
        F = np.concatenate([uf.amplitude(tr, Es[i:i + 1500], d)[0] for i in range(0, Es.size, 1500)]) * scale
        P_par = 1.602176634e-13 * np.trapezoid(F, Es) * 1e6          # Вт/мрад² → Вт/рад²
        P_lie = sr.angular_power_cached(K)["peak"]
        self.assertAlmostEqual(P_par / P_lie, 1.0, delta=0.03)

    def test_bend_fan_reaches_screen_edge_not_center(self):
        """reaches_screen ловит веер, задевающий лишь край экрана (illuminates по центру его пропускает)."""
        import bend_screen as bs
        g = self.sr.gamma()
        d = dict(C=np.zeros(3), axis=np.array([1.0, 0, 0]), vhat=np.array([0, 0, 1.0]), theta0=0.010)
        D, ang = 10.0, 0.0065                                       # центр экрана за краем веера (±5 мрад)
        C = D * np.array([math.cos(ang), math.sin(ang), 0.0])
        scr = dict(C=C, u=np.array([-math.sin(ang), math.cos(ang), 0.0]), v=np.array([0, 0, 1.0]),
                   n=-C / D, w=40.0, h=20.0)                       # ±2 мрад по горизонтали
        self.assertFalse(bs.illuminates(d, C, g))
        self.assertTrue(bs.reaches_screen(d, scr, g))
        scr_far = dict(scr, C=D * np.array([math.cos(0.02), math.sin(0.02), 0.0]))
        self.assertFalse(bs.reaches_screen(d, scr_far, g))

    def test_undulator_screen_gets_flanking_dipole_fans(self):
        """Реш.27, SCR_U01: ондулятор ≈ Nλ_u-мощность целиком + веера двух соседних диполей (gated на БД)."""
        tb = os.path.join(ROOT, "workspace", "tables")
        if not (os.path.exists(os.path.join(tb, "screens.csv")) and os.path.exists(os.path.join(tb, "element_groups.csv"))):
            self.skipTest("нет screens.csv/element_groups.csv")
        cwd = os.getcwd(); os.chdir(ROOT)
        try:
            import screen_analysis as sa
            scr = next((s for s in sa.load_screens(27) if s["sid"] == "SCR_U01"), None)
            if scr is None:
                self.skipTest("нет экрана SCR_U01 в реш.27")
            ci, byid = sa.load_lattice(27, None); groups = sa.load_groups(27)
            emit = sa.make_emitters(groups, ci, byid, 5.0)
            found = sa.find_illuminators(scr, emit)
            bends = sa.find_bends(scr, ci, byid, set().union(*groups.values()))
            fp = sa.footprint(scr, [emit[gid] for gid in found], 161, bends)
        finally:
            os.chdir(cwd)
        self.assertEqual(found, ["und_01"])
        self.assertEqual(fp["n_bend"], 2)
        self.assertGreater(fp["und_W"], 0.98 * 2977.9)
        self.assertTrue(200.0 < fp["bend_W"] < 600.0, fp["bend_W"])

    def test_live_power_map_for_cockpit(self):
        """Карта мощности кокпита (gated на БД): экран, которого нет в этой решётке, строится по оси ИД как в
        сцене; пятно ИД вытянуто по горизонтали; зазор живой (больше зазор — меньше мощность)."""
        tb = os.path.join(ROOT, "workspace", "tables")
        if not (os.path.exists(os.path.join(tb, "screens.csv")) and os.path.exists(os.path.join(tb, "element_groups.csv"))):
            self.skipTest("нет screens.csv/element_groups.csv")
        cwd = os.getcwd(); os.chdir(ROOT)
        try:
            import screen_analysis as sa
            if "und_01" not in sa.load_groups(26):
                self.skipTest("нет und_01 в реш.26")
            gaps = {f"und_{i:02d}": 50.0 for i in range(1, 39)}
            fp5 = sa.live_power_map(26, "SCR_U01", dict(gaps, und_01=5.0), 161, bends=False)
            fp8 = sa.live_power_map(26, "SCR_U01", dict(gaps, und_01=8.0), 161, bends=False)
            self.assertIsNone(sa.live_power_map(26, "SCR_BM08", gaps, 41))   # волноотвод — не ондуляторный экран
        finally:
            os.chdir(cwd)
        self.assertEqual([g for g, _ in fp5["illuminators"]][:1], ["und_01"])
        self.assertAlmostEqual(fp5["D"], 20.0, delta=1e-6)                   # экраны кокпита — 20 м от ИД
        S = fp5["S"]; m = S.shape[0] // 2                                   # ширина по полувысоте (СКО — хвосты)
        wu = int((S[:, m] >= S.max() / 2).sum()); wv = int((S[m, :] >= S.max() / 2).sum())
        self.assertGreater(wu, 2 * wv)                                       # ±K/γ по горизонтали, ~1/γ по вертикали
        self.assertGreater(fp5["und_W"], 0.98 * self.sr.angular_power_cached(
            self.sr.K_of_B0(self.sr.B0_of_gap(5.0)))["total_W"])
        self.assertLess(fp8["peak"], fp5["peak"])

    def test_bend_line_screens_at_central_module(self):
        """Экраны у линеек диполя (gated на БД): у каждой линейки из 5 модулей (и у разбитых вставкой) — экран
        на 20 м по касательной в центре модуля №3; картина по траектории — непрерывная полоса без стыков вееров,
        поток на ε_c в центре = локальная ε_c модуля №3, спектр сплошной с пиком у ε_c."""
        tb = os.path.join(ROOT, "workspace", "tables")
        if not os.path.exists(os.path.join(tb, "element_groups.csv")):
            self.skipTest("нет element_groups.csv")
        cwd = os.getcwd(); os.chdir(ROOT)
        try:
            import screen_analysis as sa
            import sr_trace as st
            bs = sa.bend_line_screens(26)
            if not bs:
                self.skipTest("нет линеек диполей в реш.26")
            b = next(x for x in bs if x["dev"].endswith("_3"))
            v = sa.bend_view(26, b["sid"], 81, "power")
            ctx, trs = sa._bend_traces(26, sa.live_context(26, b["sid"])["scr"])
            ok, psi, ec = st.stationary_points(trs[0], np.array([b["C"]]))
            sp = sa.bend_spectrum(26, b["sid"], 0.0, 0.0)
        finally:
            os.chdir(cwd)
        self.assertEqual(len({x["sid"] for x in bs}), len(bs))
        self.assertEqual(sum(1 for x in bs if x["dev"].startswith("und_J")), 2)   # инжекционные (4+1) — тоже
        for x in bs:
            np.testing.assert_allclose(x["C"], x["P"] + 20.0 * x["t"], atol=1e-9)
            self.assertAlmostEqual(float(np.dot(x["n"], x["t"])), -1.0, places=12)
        S = v["S"]; m = S.shape[1] // 2
        row = S[:, m]
        self.assertGreater(v["total"], 0.0)
        self.assertGreater(float(row.min()), 0.2 * float(row.max()))           # по горизонтали — без провалов
        self.assertTrue(bool(ok[0]))
        self.assertAlmostEqual(float(ec[0]) / b["ec"], 1.0, delta=2e-3)        # в центре — модуль №3
        self.assertAlmostEqual(float(psi[0]), 0.0, delta=1e-9)
        Ep = float(sp["E_keV"][int(np.argmax(sp["flux"]))])
        self.assertTrue(0.6 * b["ec"] < Ep < 1.0 * b["ec"], Ep)                # пик угловой плотности ~0,83·ε_c

    def test_fan_on_orbit_shifts_point_and_tangent(self):
        """Веер по орбите: точка отбора сдвигается на (x, y), касательная поворачивается на (x′, y′)."""
        import electron_source as es
        st = np.array([[2e-6, 3e-6, -1e-6, 4e-6, 0.0, 0.0]])          # вход = выход: орбита постоянна по дуге
        et = dict(inst2mag={"7": 0}, st_in=st, st_out=st)
        U = np.array([0.0, 1.0, 0.0]); V = np.array([0.0, 0.0, 1.0]); t = np.array([1.0, 0.0, 0.0])
        em = dict(iid="7", theta0=-0.01, alpha=0.002, U=U, vn=V, t=t, P=np.zeros(3))
        out = es.fan_on_orbit(em, et)
        np.testing.assert_allclose(out["P"], 2e-6 * U - 1e-6 * V, atol=1e-15)
        self.assertAlmostEqual(float(np.dot(out["t"], U)), 3e-6, delta=1e-12)
        self.assertAlmostEqual(float(np.dot(out["t"], V)), 4e-6, delta=1e-12)
        self.assertAlmostEqual(float(np.dot(out["vn"], out["t"])), 0.0, delta=1e-15)
        self.assertIs(es.fan_on_orbit(dict(em, iid="нет"), et)["P"], em["P"])   # магнита нет — без изменений

    def test_rigid_undulator_pitch_does_not_steer_light(self):
        """Жёсткий наклон ИД на 0.2 мрад поворачивает рамку железа, но НЕ электрон (нулевой интеграл поля):
        источник по орбите — пятно стоит; по рамке железа было бы D·θ (gated на БД и C++-ядро)."""
        import tempfile, shutil
        tb = os.path.join(ROOT, "workspace", "tables")
        if not os.path.exists(os.path.join(tb, "element_groups.csv")):
            self.skipTest("нет element_groups.csv")
        import ring_cpp as rc
        if not os.path.exists(rc.EXE):
            self.skipTest("нет C++-ядра MatrixTracker")
        cwd = os.getcwd(); os.chdir(ROOT)
        tmp = tempfile.mkdtemp(prefix="ss_pitch_")
        try:
            import misalign_demo as md, electron_source as es, screen_analysis as sa
            groups = sa.load_groups(27)
            if "und_19" not in groups:
                self.skipTest("нет und_19 в реш.27")
            md.perturbed_group_ws(27, "und_19", dict(dx=0, dy=0, dz=0, dyaw=0, dpitch=2e-4, droll=0), tmp)
            read0 = es.mt._read
            g = groups["und_19"]
            _, a0, _, v0 = sa.emitter_frame(g, *sa.load_lattice(27, None))
            _, a1, _, _ = sa.emitter_frame(g, *sa.load_lattice(27, tmp))
            _, d0 = es.device_source(es.electron_track(27), g)
            _, d1 = es.device_source(es.electron_track(27, tmp), g)
            self.assertIs(es.mt._read, read0)                         # ws_tables вернул чтение БД
        finally:
            os.chdir(cwd); shutil.rmtree(tmp, ignore_errors=True)
        self.assertAlmostEqual(float(np.dot(a1 - a0, v0)), -2e-4, delta=1e-8)   # железо наклонено
        self.assertLess(abs(float(np.dot(d1 - d0, v0))), 1e-6)                 # электрон — нет (< 1 мкрад)


class TestScreenHeatMapAccurate(unittest.TestCase):
    """Точный режим тепловой карты (sr_trace): Льенар по траектории против замкнутых форм и быстрого режима."""

    @classmethod
    def setUpClass(cls):
        und = os.path.join(ROOT, "undulators")
        if und not in sys.path:
            sys.path.insert(0, und)
        import undulator_sr as sr, sr_trace as st, bend_sr as bsr
        cls.sr, cls.st, cls.bsr = sr, st, bsr
        cls.g = sr.gamma()

    @staticmethod
    def _screen(D, w_m, h_m):
        return dict(C=np.array([D, 0.0, 0.0]), n=np.array([-1.0, 0.0, 0.0]), u=np.array([0.0, 1.0, 0.0]),
                    v=np.array([0.0, 0.0, 1.0]), w=w_m * 1e3, h=h_m * 1e3)

    def _map(self, tr, scr, mesh):
        import screen_analysis as sa
        hw, hh, UU, VV, cell, pts = sa._grid(scr, mesh)
        S = self.st.power_on_screen(tr, pts.reshape(-1, 3), scr["n"], self.g, self.sr.P.I).reshape(mesh, mesh)
        return S, UU, VV, cell, pts

    def _ideal_undulator(self, ppp=32, ends=None):
        sr = self.sr; p = sr.P
        kap0 = sr.B0_of_gap(5.0) / (p.E / 0.299792458)
        n = 2 * p.N
        sj = (np.arange(n) - (n - 1) / 2.0) * p.lam_u / 2.0
        a = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
        if ends is not None:
            a[0] *= ends; a[-1] *= ends
        return self.st.wiggle_trace(np.zeros(3), [1.0, 0, 0], [0, 1.0, 0], sj, a, kap0, p.lam_u, ppp, self.g,
                                    centered=ends is None), a

    def _module(self, theta=-9.573e-3, L=0.4178, kind="dq"):
        from build_screens import _rot
        st = self.st
        rho, B, ec = self.bsr.bend_geom(theta, L, self.sr.P.E)
        m = dict(C=np.zeros(3), axis=np.array([1.0, 0, 0]), uhat=np.array([0, 1.0, 0]), vhat=np.array([0, 0, 1.0]),
                 theta0=theta, rho=rho, B=B, L=L, dev="TEST", iid="1", kind=kind)
        O = m["C"] + rho * m["uhat"]
        m.update(O=O, P_in=O + _rot(m["C"] - O, m["vhat"], -theta / 2), P_out=O + _rot(m["C"] - O, m["vhat"], theta / 2),
                 t_in=st._unit(_rot(m["axis"], m["vhat"], -theta / 2)), t_out=st._unit(_rot(m["axis"], m["vhat"], theta / 2)))
        return m

    def test_parse_effects(self):
        st = self.st
        self.assertEqual(st.parse_effects(None), frozenset(st.EFFECTS))
        self.assertEqual(st.parse_effects("all,-beam"), frozenset(st.EFFECTS) - {"beam"})
        self.assertEqual(st.parse_effects("-fringe,-poles"), frozenset({"extended", "beam"}))
        self.assertEqual(st.parse_effects("extended,poles"), frozenset({"extended", "poles"}))
        self.assertEqual(st.parse_effects("none"), frozenset())
        with self.assertRaises(ValueError):
            st.parse_effects("all,-магия")

    def test_ideal_undulator_matches_lienard_table(self):
        """N идеальных периодов, далёкий экран: карта ≡ табличный Льенар быстрого режима, интеграл ≡ Лармор."""
        sr, st = self.sr, self.st
        K = sr.K_of_B0(sr.B0_of_gap(5.0)); tab = sr.angular_power_cached(K)
        tr, _ = self._ideal_undulator(ppp=64)
        D = 50.0
        scr = self._screen(D, 2 * (K + 8) / self.g * D, 2 * 8 / self.g * D)
        S, UU, VV, cell, pts = self._map(tr, scr, 81)
        R = np.linalg.norm(pts, axis=2)
        Ft = sr.angular_power(tab, pts[..., 1] / R, pts[..., 2] / R) * (D / R) / R ** 2
        self.assertAlmostEqual(S[40, 40] / Ft[40, 40], 1.0, delta=1e-3)
        self.assertLess(float(np.max(np.abs(S - Ft))) / float(Ft.max()), 3e-3)
        P_tr = st.trace_power(tr, self.g, sr.P.I)
        self.assertAlmostEqual(P_tr / tab["total_W"], 1.0, delta=1e-3)
        self.assertAlmostEqual(float(S.sum() * cell) / P_tr, 1.0, delta=3e-3)
        tr2, _ = self._ideal_undulator(ppp=128)                         # сходимость по шагу (дефолт CLI — 64)
        S2 = self._map(tr2, scr, 81)[0]
        self.assertLess(float(np.max(np.abs(S2 - S))) / float(S.max()), 1e-5)

    def test_end_poles_scale_power_and_close_trajectory(self):
        """Концевые полюса половинной силы: мощность ∝ Σa², электрон выходит параллельно входу, средний угол ~0."""
        st = self.st
        tr0, a0 = self._ideal_undulator()
        tr1, a1 = self._ideal_undulator(ends=0.5)
        r = st.trace_power(tr1, self.g, 0.2) / st.trace_power(tr0, self.g, 0.2)
        self.assertAlmostEqual(r, float(np.sum(a1 ** 2) / np.sum(a0 ** 2)), delta=1e-4)
        K_g = self.sr.K_of_B0(self.sr.B0_of_gap(5.0)) / self.g
        self.assertLess(abs(float(tr1["t"][-1, 1])), 1e-2 * K_g)             # выход ∥ входу (до полушага)
        self.assertLess(abs(float(np.mean(tr1["t"][:, 1]))), 2e-2 * K_g)     # конус по оси, не на K/γ

    def test_hard_edge_dipole_power_profile_and_exit(self):
        """Жёсткий диполь: выход = аналитика дуги, мощность = 1,2654·E²B²LI, верт. профиль = power_vprofile."""
        st, bsr = self.st, self.bsr
        m = self._module()
        tr = st.chain_trace([m], self.g, None, fringe=False)
        self.assertLess(float(np.linalg.norm(tr["exit_t"] - m["t_out"])), 1e-12)
        self.assertLess(float(np.linalg.norm(tr["exit_r"] - m["P_out"])), 1e-10)
        P_tr = st.trace_power(tr, self.g, 0.2)
        self.assertAlmostEqual(P_tr / (bsr.total_power_kW(m["B"], m["L"], 6.0, 0.2) * 1e3 * 1.2654 / 1.266), 1.0, delta=2e-4)
        D = 10.0
        scr = self._screen(D, abs(m["theta0"]) * D * 1.6, 30 / self.g * D)
        S, UU, VV, cell, _ = self._map(tr, scr, 121)
        self.assertAlmostEqual(float(S.sum() * cell) / P_tr, 1.0, delta=2e-3)
        vp = S[60, :]; ref = bsr.power_vprofile(VV[60, :] / D, self.g)
        self.assertLess(float(np.max(np.abs(vp / vp.max() - ref / ref.max()))), 1e-4)

    def test_soft_edge_keeps_bend_angle_and_scales_power(self):
        """Мягкий край (симметричный профиль, ∫f = L): касательная и линия выхода — как у жёсткого края,
        мощность ∝ ∫f²ds/L (меньше жёсткой)."""
        st = self.st
        m = self._module(kind="dipole")
        L, w = m["L"], 0.010
        x = np.linspace(-0.35, 0.35, 7001)
        f = 0.5 * (np.tanh((x + L / 2) / w) - np.tanh((x - L / 2) / w))
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])
        prof = dict(x=x, f=f, cum=cum, half=float(x[-1]), integral=float(cum[-1]))
        hard = st.chain_trace([m], self.g, None, fringe=False)
        soft = st.chain_trace([m], self.g, None, fringe=True, profiles=[prof])
        self.assertLess(float(np.linalg.norm(soft["exit_t"] - hard["exit_t"])), 1e-12)
        dr = soft["exit_r"] - hard["exit_r"]
        self.assertLess(float(np.linalg.norm(dr - np.dot(dr, hard["exit_t"]) * hard["exit_t"])), 1e-7)
        fs = f * (L / prof["integral"])
        ratio = float(np.trapezoid(fs ** 2, x) / L)
        self.assertAlmostEqual(st.trace_power(soft, self.g, 0.2) / st.trace_power(hard, self.g, 0.2), ratio, delta=1e-4)
        self.assertLess(ratio, 1.0)

    def test_beam_blur_adds_second_moments(self):
        """Свёртка с пучком: вторые моменты складываются, интеграл сохраняется (FFT-гауссиана, σ и меньше пикселя)."""
        st = self.st
        n = 201; d = 2e-5
        u = (np.arange(n) - n // 2) * d
        UU, VV = np.meshgrid(u, u, indexing="ij")
        s0 = 3e-4
        S = np.exp(-(UU ** 2 + VV ** 2) / (2 * s0 ** 2))
        for su, sv in ((2e-4, 1e-4), (5e-6, 3e-6)):
            B = st.gauss_blur(S, np.diag([su ** 2, sv ** 2]), d, d)
            self.assertAlmostEqual(float(B.sum() / S.sum()), 1.0, delta=1e-6)      # до вытекания за край
            mu = float((B * UU ** 2).sum() / B.sum()); mv = float((B * VV ** 2).sum() / B.sum())
            self.assertAlmostEqual(mu / (s0 ** 2 + su ** 2), 1.0, delta=1e-3)
            self.assertAlmostEqual(mv / (s0 ** 2 + sv ** 2), 1.0, delta=1e-3)

    def test_scr_u01_accurate_vs_fast(self):
        """Реш.27, SCR_U01 (gated на БД): точный режим находит тот же ИД и обе цепочки соседних диполей; ИД даёт
        Σa²/2N мощности быстрого (концевые полюса); мягкий край уменьшает мощность вееров."""
        tb = os.path.join(ROOT, "workspace", "tables")
        if not (os.path.exists(os.path.join(tb, "screens.csv")) and os.path.exists(os.path.join(tb, "element_groups.csv"))):
            self.skipTest("нет screens.csv/element_groups.csv")
        st = self.st
        cwd = os.getcwd(); os.chdir(ROOT)
        try:
            import screen_analysis as sa
            scr = next((s for s in sa.load_screens(27) if s["sid"] == "SCR_U01"), None)
            if scr is None:
                self.skipTest("нет экрана SCR_U01 в реш.27")
            ctx = st.make_context(27, None, None, 5.0, ppp=32)
            fp = st.footprint_accurate(scr, ctx, "all,-beam", 61)
            nof = st.footprint_accurate(scr, ctx, "all,-beam,-fringe", 61)
            tr = st._und_trace(ctx, "und_01", True)
        finally:
            os.chdir(cwd)
        self.assertEqual(fp["illuminators"], ["und_01"])
        self.assertEqual(len(fp["chains"]), 2)
        self.assertEqual(tr["n_poles"], 60)
        P_tr = st.trace_power(tr, self.g, self.sr.P.I)
        self.assertAlmostEqual(fp["und_W"] / P_tr, 1.0, delta=1e-2)          # экран ловит весь конус ИД
        ideal = self.sr.angular_power_cached(self.sr.K_of_B0(self.sr.B0_of_gap(5.0)))["total_W"]
        self.assertAlmostEqual(P_tr / ideal, tr["sum_a2"] / (2 * self.sr.P.N), delta=2e-3)
        self.assertLess(fp["bend_W"], nof["bend_W"])


import contextlib as _ctx


@_ctx.contextmanager
def _temp_cavities(rows):
    """Временные резонаторы установки в rf_cavities.csv: rows — [(lattice, (точка, нормаль), share|None)];
    по выходе файл возвращается как был."""
    import shutil
    import exp_planes as ep
    path = os.path.join(mt.TB, "rf_cavities.csv"); bak = path + ".testbak"
    had = os.path.exists(path)
    if had:
        shutil.copy2(path, bak)
    try:
        for k, (lat, pl, sh) in enumerate(rows, 1):
            ep.write_cavity(lat, pl, sh, cavity_id=900 + k, note="тест")
        yield
    finally:
        if had:
            shutil.move(bak, path)
        elif os.path.exists(path):
            os.remove(path)


@_ctx.contextmanager
def _temp_kick_tables():
    """Временные строки импульсных кикеров: pulsed_kickers.csv / kicker_timing.csv по выходе как были
    (не было файла — удаляется)."""
    import shutil
    saved = []
    for n in ("pulsed_kickers.csv", "kicker_timing.csv"):
        p = os.path.join(mt.TB, n); had = os.path.exists(p)
        if had:
            shutil.copy2(p, p + ".testbak")
        saved.append((p, had))
    try:
        yield
    finally:
        for p, had in saved:
            if had:
                shutil.move(p + ".testbak", p)
            elif os.path.exists(p):
                os.remove(p)


@_ctx.contextmanager
def _env(**kw):
    saved = {k: os.environ.get(k) for k in kw}
    os.environ.update({k: str(v) for k, v in kw.items()})
    try:
        yield
    finally:
        for k, v in saved.items():
            os.environ.pop(k, None) if v is None else os.environ.update({k: v})


def _mid_gap_plane(lattice, frac=0.5, min_gap=0.3):
    """Плоскость в середине промежутка ≥ min_gap после доли frac кольца → ((точка, нормаль), s, индекс магнита)."""
    import exp_planes as ep
    mags = ep._mags(lattice); ring, lead, s_in, s_out, C = ep._arcs(mags)
    for j in range(int(frac * len(mags)), len(mags) - 1):
        if s_in[j + 1] - s_out[j] > min_gap:
            s = 0.5 * (s_out[j] + s_in[j + 1])
            return ep.plane_at_s(lattice, s, mags), s, j
    raise RuntimeError("нет промежутка")


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск явных плоскостей опыта")
class TestExplicitPlanes(unittest.TestCase):
    """Явные старт/финиш-плоскости опыта (experiments.csv), пучок в рамке старта, портал разомкнутой линии,
    резонаторы установки (rf_cavities.csv) — Python-эталон beam_injection / exp_planes / elegant_diag."""

    def test_local_canonical_roundtrip(self):
        import beam_injection as bi
        loc = (1e-3, 2e-4, -5e-4, -1e-4, 0.01, 3e-3)
        back = bi.canonical_to_local(bi.local_to_canonical(loc))
        for a, b in zip(loc, back):
            self.assertAlmostEqual(a, b, places=15)
        pl = bi.PlaneMagnet(np.array([1.0, 2.0, 0.0]), np.array([0.6, 0.8, 0.0]))
        pos, vel = bi.local_to_global(pl, (0, 0, 0, 0, 0, 0))
        self.assertTrue(np.allclose(pos, pl.center))
        self.assertAlmostEqual(float(np.linalg.norm(vel)), 2.967e8, delta=1e-3)     # δ = 0 ⇒ опорная скорость
        self.assertTrue(np.allclose(vel / np.linalg.norm(vel), pl.axis))
        pos, vel = bi.local_to_global(pl, (1e-3, 2e-3, 0, 0, 0, 0))
        self.assertAlmostEqual(float((pos - pl.center) @ pl.uhat), 1e-3, places=15)
        d = vel / np.linalg.norm(vel)
        self.assertAlmostEqual(float(d @ pl.uhat) / float(d @ pl.axis), 2e-3, places=12)   # xp = dx/ds

    def test_tables_migrated(self):
        import beam_injection as bi
        rows = bi._rows("experiments.csv")
        self.assertIn("start_x", rows[0]); self.assertNotIn("p1x", rows[0])
        for r in rows:
            bi.planes_of(r)                                       # у каждого опыта задан старт
        self.assertIsNone(bi.planes_of(bi.experiment_row(1))[1])          # СИЛА — кольцо
        self.assertIsNotNone(bi.planes_of(bi.experiment_row(3))[1])       # SP1 — разомкнутая линия
        self.assertEqual([p["local"][0] for p in bi.load_group(1)][:3], [0.0, -0.00025, 0.00025])
        self.assertTrue(os.path.exists(os.path.join(mt.TB, "rf_cavities.csv")))

    def test_locate_plane_and_path(self):
        import io, contextlib
        import beam_injection as bi
        eff, phys, coup, trans = mt.load_db()
        with contextlib.redirect_stdout(io.StringIO()):
            sp1 = mt.build_line(3, eff, phys, coup, trans)
        st, fin = bi.planes_of(bi.experiment_row(3))
        self.assertEqual(bi.locate_plane(sp1, *st)[0], -1)                 # до первого элемента
        self.assertEqual(bi.locate_plane(sp1, *fin)[0], len(sp1) - 1)      # после последнего
        self.assertFalse(bi.seam_is_drift(sp1))
        m = next(x for x in sp1 if x.L > 0.1)
        with self.assertRaises(ValueError):                                 # плоскость внутри элемента
            bi.locate_plane(sp1, m.center, m.axis)
        with self.assertRaises(ValueError):                                 # не кольцо — финиш обязателен
            bi.build_path(sp1, st, None)
        with self.assertRaises(ValueError):                                 # финиш раньше старта на линии
            bi.build_path(sp1, fin, st)
        path = bi.build_path(sp1, st, fin)
        self.assertFalse(path["closed"]); self.assertEqual(len(path["ordered"]), len(sp1))

    def test_portal_open_line(self):
        """Портал: второй пролёт разомкнутой линии = первый пролёт, применённый к состоянию финиша."""
        import beam_injection as bi
        import matrix_traj_export as mte
        d = mte.long_trajectories(3, 1, n_turns=2)
        self.assertEqual(d["n_turns"], 2); self.assertFalse(d["closed"]); self.assertEqual(d["surv_turns"][0], 2)
        plane, ordered, states, info = bi.inject(3, 1)
        s1 = bi.track_turn(plane, ordered, states[:1], info["finish"])[0][0]
        s2 = bi.track_turn(plane, ordered, [s1], info["finish"])[0][0]
        self.assertTrue(np.array_equal(np.asarray(d["finals"][0]), np.asarray(s2)))

    def test_exp_planes_authoring(self):
        import exp_planes as ep
        import beam_injection as bi
        st, fin = ep.default_planes(26)
        self.assertIsNone(fin)                                              # кольцо
        self.assertLess(float(np.linalg.norm(st[0])), 1e-6)                 # середина шва СИЛА ≈ (0, 0, 0)
        st3, fin3 = ep.default_planes(3)
        self.assertIsNotNone(fin3)                                          # линия: финиш на выходе последнего
        mags = ep._mags(26); _, _, s_in, s_out, _ = ep._arcs(mags)
        j = next(k for k in range(len(mags)) if mags[k].L > 0.2)
        with self.assertRaises(ValueError):
            ep.plane_at_s(26, 0.5 * (s_in[j] + s_out[j]), mags)            # внутри элемента
        pl = ep.plane_after(26, mags[j].dev, 1e-3, mags)
        self.assertEqual(bi.locate_plane(mags, *pl)[0], j)

    def test_ideal_target_energy(self):
        """Идеальный резонатор знает энергию «после себя»: частица с любой начальной δ после него имеет заданную
        δ (пусто — проектная, δ = 0); резонатор в середине — заданная δ минус потеря второй половины пути."""
        import matrix_traj_export as mte
        st = [(0.0, 0.0, 0.0, 0.0, 0.0, 2e-3), (1e-4, 0.0, 0.0, 0.0, 0.0, -3e-3)]
        tgt = dict(type="ideal", delta=1e-3)
        with _env(SS_SR="1", SS_DESIGN_MEV="3000"):
            d0 = mte.long_trajectories(2, 1, n_turns=2, rf="ideal", states_override=st)
            d1 = mte.long_trajectories(2, 1, n_turns=2, rf=tgt, states_override=st)
            pl, s, j = _mid_gap_plane(2)
            with _temp_cavities([(2, pl, None)]):
                dm = mte.long_trajectories(2, 1, n_turns=2, rf=tgt, states_override=st)
                dz = mte.long_trajectories(2, 1, n_turns=2, rf=tgt, states_override=[(0.0,) * 5 + (1e-3,)])
        self.assertEqual([f[5] for f in d0["finals"]], [0.0, 0.0])          # прилетели с ±δ — ушли с проектной
        self.assertEqual([f[5] for f in d1["finals"]], [1e-3, 1e-3])
        loss = 1e-3 - dz["finals"][0][5]                                    # потеря от резонатора до финиша
        self.assertGreater(loss, 0.0)
        for f in dm["finals"]:                       # δ на финише от начальной не зависит (только орбита → СИ)
            self.assertAlmostEqual(f[5], dz["finals"][0][5], delta=0.05 * loss)

    def test_cavity_mid_ring_ideal(self):
        """ideal-резонатор в середине кольца: на нём δ ← заданная (проектная), к финишу — потеря второй
        половины; резонатор в плоскости старта ≡ решётка без резонаторов (бит-в-бит)."""
        import matrix_traj_export as mte
        pl, s, j = _mid_gap_plane(2)
        with _env(SS_SR="1", SS_DESIGN_MEV="3000"):
            ref = mte.long_trajectories(2, 1, n_turns=2, rf="ideal")
            self.assertEqual(ref["finals"][0][5], 0.0)                    # ВЧ в конце: δ восстановлено
            with _temp_cavities([(2, pl, None)]):
                mid = mte.long_trajectories(2, 1, n_turns=2, rf="ideal")
            with _temp_cavities([(2, (np.zeros(3), np.array([1.0, 0.0, 0.0])), None)]):
                seam = mte.long_trajectories(2, 1, n_turns=2, rf="ideal")
        self.assertEqual(len(mid["cavities"]), 1)
        self.assertAlmostEqual(mid["cavities"][0][1], s, delta=0.01)
        d_mid = mid["finals"][0][5]
        self.assertLess(d_mid, 0.0)                                         # теряется после резонатора
        self.assertTrue(np.array_equal(np.asarray(seam["finals"]), np.asarray(ref["finals"])))

    def test_ring_serve_cavity(self):
        """Сервер кольца с резонатором установки в середине: C++ ≡ Python (оборот и ELEM), скачок δ 6D-орбиты
        стоит на резонаторе, отпечаток кэша откликов меняется только при наличии резонаторов."""
        import ring_analysis as ra, ring_cpp as rc
        from control.response import ring_fingerprint
        if not os.path.exists(os.path.join(ROOT, "build", exe_name("MatrixTracker"))):
            self.skipTest("C++ ядро не собрано")
        rf = dict(volt_mv=2.5, harmonic=75, phase_deg=145.0)
        kw = dict(start_before_first_m=1.4335, sr=True, rf=rf, design_mev=3000.0)
        r0 = rc.CppRing(2, nsl=10, order=4, **kw)
        fp0 = ring_fingerprint(r0); r0.close()
        pl, s, j = _mid_gap_plane(2)
        with _temp_cavities([(2, pl, None)]):
            cr = rc.CppRing(2, nsl=10, order=4, **kw)
            with _env(SS_SECTOR_NSL="10", SS_SECTOR_ORDER="4"):
                pr = ra.Ring(2, **kw)
            try:
                st = [[1e-4, 0, 5e-5, 0, 0, 0], [0, 0, 0, 0, 1e-3, 1e-3]]
                self.assertLess(float(np.max(np.abs(np.asarray(cr.track(st, 2)) - pr.track(st, 2)))), 1e-9)
                e1 = np.asarray(cr.track_elements(st[:1]))[0]; e2 = pr.track_elements(st[:1])[0]
                self.assertLess(float(np.nanmax(np.abs(e1 - e2))), 1e-9)
                self.assertNotEqual(ring_fingerprint(cr), fp0)
                co = ra.closed_orbit(cr, dims=6)
                o = np.asarray(cr.track_elements([co]))[0]
                k = int(np.argmax(np.diff(o[:, 5])))
                self.assertEqual(k, j)                                       # скачок δ — между магнитами j и j+1
            finally:
                cr.close()

    def test_elegant_splice_rf(self):
        import elegant_diag as ed
        par = {"D1": {"type": "EDRIFT", "L": 2.0}, "Q": {"type": "KQUAD", "L": 0.5}}
        names, defs = ed._splice_rf(["D1", "Q", "D1"], par, dict(volt_mv=1.0, phase_deg=180.0, at_s=[0.5, 2.0]), 1e8)
        self.assertEqual(names, ["RFD1", "RFC0", "RFD2", "RFC1", "Q", "D1"])
        self.assertIn("RFD1: EDRIFT, L=0.5", defs); self.assertIn("RFD2: EDRIFT, L=1.5", defs)
        with self.assertRaises(ValueError):
            ed._splice_rf(["D1", "Q", "D1"], par, dict(volt_mv=1.0, phase_deg=180.0, at_s=[2.2]), 1e8)


def _kick_spec(center, normal, **kw):
    spec = dict(kicker_id="k", center=np.asarray(center, float), normal=np.asarray(normal, float), shape="half_sine",
                duration=2e-9, rise=0.0, theta_x=0.0, theta_y=0.0, delay=0.0, period=None)
    spec.update(kw)
    return spec


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет lattice_configs.csv — пропуск импульсных кикеров")
class TestPulsedKickers(unittest.TestCase):
    """Импульсные кикеры (pulsed_kickers.csv + kicker_timing.csv) — Python-эталон: форма импульса, загрузка и
    проверки ввода, место на пути (промежуток / старт / финиш), удар по времени прихода частицы, авторинг
    exp_planes, вставка BUMPER для elegant. C++ ≡ Python и аналитика по Твиссу — test_matrix_core_ring."""

    def test_pulse_shape_and_period(self):
        import beam_injection as bi
        D = 2e-9
        self.assertEqual(bi.pulse_shape("half_sine", D, 0.0, -1e-12), 0.0)
        self.assertEqual(bi.pulse_shape("half_sine", D, 0.0, 1.001 * D), 0.0)
        self.assertAlmostEqual(bi.pulse_shape("half_sine", D, 0.0, D / 2), 1.0, places=15)
        self.assertAlmostEqual(bi.pulse_shape("half_sine", D, 0.0, D / 6), 0.5, places=12)
        self.assertAlmostEqual(bi.pulse_shape("trapezoid", 4e-9, 1e-9, 0.5e-9), 0.5, places=12)
        self.assertEqual(bi.pulse_shape("trapezoid", 4e-9, 1e-9, 2e-9), 1.0)
        self.assertAlmostEqual(bi.pulse_shape("trapezoid", 4e-9, 1e-9, 3.75e-9), 0.25, places=12)
        self.assertEqual(bi.pulse_shape("trapezoid", 4e-9, 0.0, 1e-9), 1.0)          # фронт 0 — прямоугольник
        k = bi.KickerNode(np.zeros(3), np.array([1.0, 0.0, 0.0]),
                          _kick_spec(np.zeros(3), (1, 0, 0), theta_x=1e-3, theta_y=-2e-3, delay=1e-6, period=1e-6))
        T = lambda sec: sec * bi.C                                                     # время, с → c·t, м
        self.assertEqual(k.amplitude(T(0.5e-6)), 0.0)                                  # до задержки — 0
        self.assertAlmostEqual(k.amplitude(T(1e-6 + D / 2)), 1.0, places=9)
        self.assertAlmostEqual(k.amplitude(T(5e-6 + D / 2)), 1.0, places=6)            # период повторяет импульс
        self.assertEqual(k.amplitude(T(5.5e-6)), 0.0)                                  # между импульсами
        st = k.kick((1e-3, 2e-4, 0.0, 0.0, 0.3, 1e-3), T(1e-6 + D / 6))
        self.assertAlmostEqual(st[1], 2e-4 + 0.5e-3, places=12); self.assertAlmostEqual(st[3], -1e-3, places=12)
        self.assertEqual((st[0], st[2], st[4], st[5]), (1e-3, 0.0, 0.3, 1e-3))         # x, y, t, δ — не трогает

    def test_load_and_validation(self):
        import beam_injection as bi
        import exp_planes as ep
        pl = _mid_gap_plane(3, 0.5, 0.1)[0]
        with _temp_kick_tables():
            self.assertEqual(bi.load_kickers(3, 3), [])                                # нет таблиц — нет кикеров
            kid = ep.write_kicker(3, pl, 2e-9, note="тест")
            self.assertEqual(bi.load_kickers(3, 3), [])                                # без импульса не стреляет
            ep.write_timing(3, kid, 1e-4, 0.0, 1e-8)
            kk = bi.load_kickers(3, 3)
            self.assertEqual((len(kk), kk[0]["shape"], kk[0]["theta_x"], kk[0]["period"]), (1, "half_sine", 1e-4, None))
            ep.write_timing(3, kid, 2e-4, 0.0, 2e-8, add=True)
            self.assertEqual(len(bi.load_kickers(3, 3)), 2)                            # два импульса одного кикера
            ep.write_timing(3, kid, 1e-4, 0.0, 0.0, period=1e-9)
            with self.assertRaises(ValueError):
                bi.load_kickers(3, 3)                                                  # период короче импульса
            with self.assertRaises(ValueError):
                ep.write_kicker(3, pl, 2e-9, shape="gauss")
            with self.assertRaises(ValueError):
                ep.write_kicker(3, pl, 4e-9, shape="trapezoid", rise=3e-9)             # 2·rise > duration
            with self.assertRaises(ValueError):
                ep.write_timing(3, 999, 1e-4)                                          # кикера нет у решётки

    def test_place_on_path(self):
        import beam_injection as bi
        import exp_planes as ep
        st, fin = bi.planes_of(bi.experiment_row(3))
        mid = _mid_gap_plane(3, 0.5, 0.1)[0]
        path = bi.build_path(ep._mags(3), st, fin, kickers=[_kick_spec(*fin, kicker_id="f"),
                                                             _kick_spec(*st, kicker_id="s"), _kick_spec(*mid)])
        self.assertEqual([k.kicker_id for k in path["kick_end"]], ["f"])              # финиш — конец пролёта
        self.assertEqual([k.kicker_id for k in path["kick_begin"]], ["s"])            # старт линии — начало
        self.assertEqual(len(path["kickers"]), 1)
        self.assertEqual(path["kick_end"][0].s_k, path["line_len"]); self.assertEqual(path["kick_begin"][0].s_k, 0.0)
        self.assertTrue(0.0 < path["kickers"][0].s_k < path["line_len"])
        self.assertIn(path["kickers"][0], path["ordered"])
        st26, _ = bi.planes_of(bi.experiment_row(26))
        p26 = bi.build_path(ep._mags(26), st26, None, kickers=[_kick_spec(*st26)])
        self.assertEqual((len(p26["kick_end"]), p26["kick_begin"]), (1, []))           # кольцо: старт = конец пролёта

    def test_kick_by_arrival_time(self):
        """Δpx на узле = θ·f(τ), τ = T/c − delay, T = ct + Δt_raw/β0 — у каждой частицы по своему времени прихода."""
        import io
        import beam_injection as bi
        import exp_planes as ep
        import matrix_traj_export as mte
        D, TH = 2e-9, 1e-4
        states = [(0, 0, 0, 0, 0.0, 0), (0, 0, 0, 0, 0.15, 0), (0, 0, 0, 0, -0.25, 0)]
        with _temp_kick_tables():
            kid = ep.write_kicker(3, _mid_gap_plane(3, 0.5, 0.1)[0], D, note="тест")
            delay = ep.reference_time(3, kid, 0) - D / 2
            ep.write_timing(3, kid, TH, 0.0, delay)
            with _ctx.redirect_stdout(io.StringIO()):
                d = mte.local_trajectories(3, 1, n_turns=1, exp_id=3, states_override=states)
        ki = d["elem_names"].index(f"PULSE_KICKER_{kid}")
        beta0 = 6000.0 / math.hypot(6000.0, mte.MC2_MEV)
        fs = []
        for st0, L in zip(states, d["locals"]):
            ent, ext = L[L[:, 7] == ki]
            T = st0[4] + (ent[5] - st0[4]) / beta0
            f = bi.pulse_shape("half_sine", D, 0.0, T / bi.C - delay)
            self.assertAlmostEqual(ext[2] - ent[2], TH * f, delta=1e-15)               # px: удар θ·f(τ)
            self.assertEqual((ext[1], ext[3], ext[4], ext[5]), (ent[1], ent[3], ent[4], ent[5]))
            fs.append(f)
        self.assertGreater(fs[0], 0.999999)                                            # опорная — на пике
        self.assertLess(fs[1], 0.75); self.assertLess(fs[2], 0.75)                     # раньше/позже — слабее

    def test_exp_planes_kicker_authoring(self):
        import io
        import beam_injection as bi
        import exp_planes as ep
        with _temp_kick_tables():
            kid = ep.write_kicker(3, _mid_gap_plane(3, 0.5, 0.1)[0], 3e-9, "trapezoid", 1e-9, note="тест")
            self.assertEqual(ep.kicker_row(3, kid)["shape"], "trapezoid")
            t0, t1 = ep.reference_time(3, kid, 0), ep.reference_time(3, kid, 1)
            beta0 = 6000.0 / math.hypot(6000.0, ep.MC2_MEV)
            L = bi.inject(3, 1, exp_id=3)[3]["line_len"]
            self.assertAlmostEqual((t1 - t0) * bi.C, L / beta0, places=9)              # без ВЧ: пролёт = L/β0
            ep.write_timing(3, kid, 1e-4, 0.0, t0)
            out = io.StringIO()
            with _ctx.redirect_stdout(out):
                ep.show(3)
            self.assertIn("кикер : " + str(kid), out.getvalue())
            ep.remove_timing(3, kid)
            self.assertEqual(bi.load_kickers(3, 3), [])
            with self.assertRaises(ValueError):
                ep.remove_timing(3, kid)
            ep.remove_kicker(3, kid)
            with self.assertRaises(ValueError):
                ep.kicker_row(3, kid)

    def test_elegant_splice_kicker(self):
        import elegant_diag as ed
        par = {"D1": {"type": "EDRIFT", "L": 2.0}, "Q": {"type": "KQUAD", "L": 0.5}}
        kk = [dict(at_s=0.5, theta_x=1e-3, theta_y=0.0, shape="half_sine", duration=2e-9, time_offset=1e-9),
              dict(at_s=3.0, theta_x=0.0, theta_y=2e-3, shape="trapezoid", duration=4e-9, rise=1e-9, period=1e-6)]
        try:
            names, defs = ed._splice_rf(["D1", "Q", "D1"], par, dict(volt_mv=1.0, phase_deg=180.0, at_s=[0.5]), 1e8, kk)
            self.assertEqual(names, ["RFD1", "RFC0", "PK0", "RFD2", "Q", "RFD3", "PK1", "RFD4"])
            self.assertIn("RFD2: EDRIFT, L=1.5", defs); self.assertIn("RFD3: EDRIFT, L=0.5", defs)
            pk = {d.split(":")[0]: d for d in defs if ": BUMPER" in d}
            self.assertIn("ANGLE=0.001,", pk["PK0"]); self.assertIn("TILT=0,", pk["PK0"])
            self.assertIn("TIME_OFFSET=1.0000000000000001e-09", pk["PK0"]); self.assertIn("PERIODIC=0", pk["PK0"])
            self.assertIn(f"TILT={math.pi / 2:.17g}", pk["PK1"]); self.assertIn("PERIODIC=1", pk["PK1"])
            with open(os.path.join(ed.DIAG, "PK1_wf.sdds"), encoding="ascii") as f:
                rows = [ln.split() for ln in f.read().splitlines()[6:]]
            self.assertEqual(len(rows), 4002)                                          # импульс + нули до периода
            self.assertEqual((float(rows[-1][0]), float(rows[-1][1])), (1e-6, 0.0))
        finally:
            for n in ("PK0_wf.sdds", "PK1_wf.sdds"):
                p = os.path.join(ed.DIAG, n)
                if os.path.exists(p):
                    os.remove(p)


def _engine_ext_available():
    try:
        import synchrotronsim
        return synchrotronsim.ENGINE_AVAILABLE
    except Exception:
        return False


@unittest.skipUnless(_engine_ext_available() and os.path.exists(_EXE)
                     and os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет собранного расширения synchrotronsim._engine / MatrixTracker / БД")
class TestEngineExtension(unittest.TestCase):
    """Шаг 0 плана библиотеки (docs/LIBRARY_ARCHITECTURE.md): C++-ядро как модуль расширения CPython против
    сервера MatrixTracker --ring-serve — те же числа бит-в-бит (включая квантовое излучение ISR) при той же
    последовательности запросов; разные кольца в разных потоках Python считаются параллельно и одинаково."""

    LAT, MEV = 26, 6000.0
    RF = {"volt_mv": 8.0, "harmonic": 1296.0, "phase_deg": 165.0}   # ВЧ СИЛА-26 (facility_report.info_of)

    @staticmethod
    def _bits_equal(a, b):
        a = np.ascontiguousarray(a, float); b = np.ascontiguousarray(b, float)
        return a.shape == b.shape and np.array_equal(a.view(np.uint64), b.view(np.uint64))

    def _pair(self, source="db", **kw):
        """Сервер (процесс, читает CSV сам) и расширение (по умолчанию — из массивов модели Machine)."""
        import ring_cpp as rc
        from synchrotronsim.engine import ExtRing
        base = dict(lattice=self.LAT, design_mev=self.MEV, nsl=5, seed=4242)
        base.update(kw)
        s, e = rc.CppRing(backend="pipe", **base), ExtRing(source=source, **base)
        self.addCleanup(s.close)
        return s, e

    def _states(self, n):
        g = np.random.default_rng(11)
        st = np.zeros((n, 6))
        st[:, 0] = g.normal(0, 1e-4, n); st[:, 2] = g.normal(0, 3e-5, n); st[:, 5] = g.normal(0, 1e-3, n)
        return st

    def test_geometry_and_tracking_bitwise_sr_isr_rf(self):
        """СИ + ISR + ВЧ, полюса ИД 1 срез: описание кольца и все виды трекинга — бит-в-бит."""
        from control.ring_model import ID_POLE_DEVS
        s, e = self._pair(sr=True, isr=True, rf=dict(self.RF), nsl_dev={d: 1 for d in ID_POLE_DEVS})
        self.assertEqual((s.n_mag, s.n_pieces), (e.n_mag, e.n_pieces))
        self.assertTrue(self._bits_equal([s.circumference], [e.circumference]))
        self.assertTrue(self._bits_equal(s.s_in, e.s_in) and self._bits_equal(s.s_out, e.s_out))
        self.assertEqual((s.kinds, s.names), (e.kinds, e.names))
        st = self._states(12)
        self.assertTrue(self._bits_equal(s.track(st, 2), e.track(st, 2)))
        self.assertTrue(self._bits_equal(s.track_tbt(st, 6), e.track_tbt(st, 6)))
        self.assertTrue(self._bits_equal(s.track_tbt(st * 30, 4, limit=1e-3), e.track_tbt(st * 30, 4, limit=1e-3)))
        self.assertTrue(self._bits_equal(s.track_elements(st[:4]), e.track_elements(st[:4])))
        self.assertTrue(self._bits_equal(s.track_pieces(st[:2]), e.track_pieces(st[:2])))
        self.assertTrue(self._bits_equal(np.array(s.piece_meta(), float), np.array(e.piece_meta(), float)))

    def test_live_perturbations_and_rf_bitwise(self):
        """Возмущения (RESET + набор, кики без RESET, BSCALE, аддитивный DKV) и частота ВЧ на лету — бит-в-бит."""
        s, e = self._pair(sr=True, rf=dict(self.RF))
        quads = [j for j, k in enumerate(s.kinds) if k == "quad"][:3]
        bends = [j for j, k in enumerate(s.kinds) if k in ("dipole", "dq")][:1]
        pert = dict(kicks=[(quads[0], 2e-6)], vkicks=[(quads[1], -1e-6)], dks=[(quads[2], 1e-3)],
                    dk2s=[(quads[0], 0.5)], dsks=[(quads[1], 1e-4)])
        for r in (s, e):
            r.apply_perturb(**pert)
            r.set_kicks(kicks=[(quads[2], 3e-6)]); r.set_bend_scale([(bends[0], 0.97)])
            r.add_perturb(dkvs=[(quads[0], 2e-4)])
            r.set_rf(freq_hz=r.rf_frequency() + 50.0)
        st = self._states(6)
        self.assertTrue(self._bits_equal(s.track_tbt(st, 4), e.track_tbt(st, 4)))
        self.assertTrue(self._bits_equal(s.track_elements(st[:3]), e.track_elements(st[:3])))
        self.assertTrue(self._bits_equal(s.track_pieces(st[:2]), e.track_pieces(st[:2])))
        A = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
        c, sn = math.cos(0.01), math.sin(0.01)
        B = np.array([[0.5, 0.01, 0], [c, sn, 0], [-sn, c, 0], [0, 0, 1]], float)
        self.assertTrue(self._bits_equal(s.junction(st, A, B), e.junction(st, A, B)))

    def test_csv_source_bitwise(self):
        """Расширение, читающее таблицы само (source=csv), — тоже бит-в-бит с сервером (и с источником db)."""
        s, e = self._pair(source="csv", sr=True)
        st = self._states(6)
        self.assertTrue(self._bits_equal(s.track_tbt(st, 3), e.track_tbt(st, 3)))
        self.assertTrue(self._bits_equal(np.array(s.piece_meta(), float), np.array(e.piece_meta(), float)))

    def test_snapshot_restore_clone_reset(self):
        """Шаг 2: снимок/откат возвращают кольцо бит-в-бит; клон независим; RESET без пересборки сохраняет ВЧ."""
        from synchrotronsim.engine import ExtRing
        e = ExtRing(lattice=self.LAT, design_mev=self.MEV, nsl=5, sr=True, rf=dict(self.RF))
        st = self._states(4)
        ref = e.track_tbt(st, 3)
        q = [j for j, k in enumerate(e.kinds) if k == "quad"]
        snap = e.snapshot()
        e.apply_perturb(dks=[(q[0], 1e-3)]); e.set_kicks(kicks=[(q[1], 5e-6)])
        self.assertFalse(self._bits_equal(e.track_tbt(st, 3), ref))
        e.restore(snap)
        self.assertTrue(self._bits_equal(e.track_tbt(st, 3), ref))
        c = e.clone(); c.set_kicks(kicks=[(q[1], 5e-6)])
        self.assertTrue(self._bits_equal(e.track_tbt(st, 3), ref))
        self.assertFalse(self._bits_equal(c.track_tbt(st, 3), ref))
        # RESET (откат к свежему) сохраняет ВЧ на лету — как команда RESET сервера
        s, _ = self._pair(sr=True, rf=dict(self.RF))
        for r in (s, e):
            r.set_rf(freq_hz=r.rf_frequency() + 40.0)
            r.apply_perturb(kicks=[(q[2], 1e-6)])
        self.assertTrue(self._bits_equal(s.track_tbt(st, 3), e.track_tbt(st, 3)))

    def test_cppring_backends_identical(self):
        """CppRing на бэкенде ext (модуль, по умолчанию) и pipe (сервер): одинаковые числа и описание магнитов,
        копии with_*/perturbed, apply_perturb; timeout ⇒ pipe (сторож процесса)."""
        import ring_cpp as rc
        kw = dict(lattice=self.LAT, design_mev=self.MEV, nsl=5, sr=True, rf=dict(self.RF))
        p, x = rc.CppRing(backend="pipe", **kw), rc.CppRing(backend="ext", **kw)
        self.addCleanup(p.close); self.addCleanup(x.close)
        self.assertEqual((p.backend, x.backend), ("pipe", "ext"))
        t = rc.CppRing(lattice=2, design_mev=3000.0, timeout=30)
        self.addCleanup(t.close)
        self.assertEqual(t.backend, "pipe")
        st = self._states(4)
        q = [j for j, k in enumerate(p.kinds) if k == "quad"]
        for r in (p, x):
            r.apply_perturb(dks=[(q[0], 1e-3)], kicks=[(q[1], 2e-6)])
        self.assertTrue(self._bits_equal(p.track_elements(st[:2]), x.track_elements(st[:2])))
        p2, x2 = p.with_kick(q[2], 3e-6), x.with_kick(q[2], 3e-6)
        self.addCleanup(p2.close)
        self.assertEqual(x2.backend, "ext")
        self.assertTrue(self._bits_equal(p2.track_tbt(st, 3), x2.track_tbt(st, 3)))
        self.assertTrue(self._bits_equal(p.track_tbt(st, 3), x.track_tbt(st, 3)))       # родители не задеты
        y2 = x.perturbed(kicks=[(q[2], 3e-6)])
        self.assertTrue(self._bits_equal(y2.track_tbt(st, 3), x2.track_tbt(st, 3)))
        self.assertEqual(len(x.mags), len(p.mags))                                       # Python-описание
        self.assertEqual([m.k1 for m in x.mags[:50]], [m.k1 for m in p.mags[:50]])

    def test_description_from_engine(self):
        """Шаг 3: CppRing.mags/plane — из ядра (DESCRIBE): раскладка C++ = Python; ext и pipe дают одно описание;
        оно совпадает с Python-эталоном (геометрия ≤1e-15, силы и флаги точно, тюны ring_matrix одинаковы) и
        следует за возмущениями ядра без повторного проигрывания."""
        import io, contextlib
        import ring_cpp as rc, ring_analysis as ra, ring_matrix as rmx
        import synchrotronsim
        self.assertEqual(rc.DESCRIBE_COLS, synchrotronsim._engine.DESCRIBE_COLS)
        kw = dict(lattice=self.LAT, design_mev=self.MEV, nsl=5)
        p, x = rc.CppRing(backend="pipe", **kw), rc.CppRing(backend="ext", **kw)
        self.addCleanup(p.close); self.addCleanup(x.close)
        dp, dx = p._be.describe(), x._be.describe()
        for k in dp:
            self.assertTrue(self._bits_equal(dp[k], dx[k]), k)
        with contextlib.redirect_stdout(io.StringIO()):
            ref = ra.Ring(self.LAT, design_mev=self.MEV)              # Python-эталон (своя сборка из CSV)
        self.assertEqual(len(x.mags), len(ref.mags))
        geo = max(float(np.max(np.abs(np.asarray(a.fr_out[q]) - np.asarray(b.fr_out[q]))))
                  for a, b in zip(x.mags, ref.mags) for q in range(4))
        self.assertLess(geo, 1e-12)
        for a, b in zip(x.mags[:400], ref.mags[:400]):
            self.assertEqual((a.theta0, a.k1, a.k1v, a.k2, a.k3, a.kEdge1, a.kEdge2, a.sector, a.kick_sector,
                              a.nonzero, a.nsl, a.dev, a.kind),
                             (b.theta0, b.k1, b.k1v, b.k2, b.k3, b.kEdge1, b.kEdge2, b.sector, b.kick_sector,
                              b.nonzero, b.nsl, b.dev, b.kind))
        tx, tr = rmx.twiss(x._python_ring()), rmx.twiss(ref)
        self.assertAlmostEqual(tx["nux"], tr["nux"], places=11); self.assertAlmostEqual(tx["nuy"], tr["nuy"], places=11)
        # описание следует за ядром: DK, BSCALE полюса ИД, кик «ящика» корректора, поворот (YAW)
        from control.ring_model import ID_POLE_DEVS
        q = [j for j, k in enumerate(x.kinds) if k == "quad"]
        kk = [j for j, k in enumerate(x.kinds) if k == "kicker"]
        pole = [j for j, nm in enumerate(x.names) if nm in ID_POLE_DEVS]
        m0 = x.mags
        x.apply_perturb(dks=[(q[0], 1e-3)], kicks=[(kk[0], 2e-6)], bscales=[(pole[0], 0.5)] if pole else ())
        m1 = x.mags
        self.assertEqual(m1[q[0]].k1, m0[q[0]].k1 + 1e-3)
        if m1[kk[0]].kick_sector:
            self.assertEqual(m1[kk[0]].hkick, 2e-6)
        else:
            self.assertEqual(m1[kk[0]].theta0, m0[kk[0]].theta0 + 2e-6)   # «ящик»: кик = θ0 box-тела
        if pole:
            self.assertEqual(m1[pole[0]].theta0, 0.5 * m0[pole[0]].theta0)
            self.assertFalse(np.array_equal(m1[pole[0]].fr_out[1], m0[pole[0]].fr_out[1]))   # рамки заново
        y = x.with_yaw(q[5], 1e-4)
        self.assertGreater(float(np.linalg.norm(y.mags[q[5]].fr_in[1] - m1[q[5]].fr_in[1])), 9e-5)

    def test_parallel_rings_in_python_threads(self):
        """Копии машины в разных потоках Python: GIL отпущен — результаты те же, что у одиночного прогона."""
        import threading
        from synchrotronsim.engine import ExtRing
        rings = [ExtRing(lattice=self.LAT, design_mev=self.MEV, nsl=5, sr=True) for _ in range(4)]
        st = self._states(1)
        ref = rings[0].track(st, 2)
        out = [None] * len(rings)

        def run(i):
            out[i] = rings[i].track(st, 2)

        th = [threading.Thread(target=run, args=(i,)) for i in range(len(rings))]
        [t.start() for t in th]; [t.join() for t in th]
        for o in out:
            self.assertTrue(self._bits_equal(o, ref))


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")), "нет БД")
class TestDataModel(unittest.TestCase):
    """Шаг 1 плана библиотеки: модель данных `synchrotronsim.model` и слой доступа `synchrotronsim.db`.
    Магниты модели — бит-в-бит как `matrix_tracker.build_line` (склейка обмоток, фолбэки на тип устройства);
    отпечаток стабилен, локален (чужая решётка его не меняет) и чувствителен; пути не зависят от текущей папки."""

    @staticmethod
    def _bits(a):
        return np.ascontiguousarray(np.asarray(a, float)).view(np.uint64)

    def _same(self, a, b, what):
        self.assertTrue(np.array_equal(self._bits(a), self._bits(b)), what)

    def test_elements_match_build_line(self):
        import io, contextlib
        import synchrotronsim as ss
        db = ss.open_db()
        eff, phys, coup, trans = mt.load_db()
        for lat in (2, 3, 27):                     # 27 — многообмоточные магниты и ондуляторы
            if lat not in db.lattices():
                continue
            with self.subTest(lattice=lat):
                E = db.machine(lat).elements
                with contextlib.redirect_stdout(io.StringIO()):
                    mags = mt.build_line(lat, eff, phys, coup, trans)
                self.assertEqual(len(E), len(mags))
                self._same(E.centers, [q.center for q in mags], "центры")
                for a in ("yaw", "pitch", "roll", "edge1", "edge2", "lf1", "lf2"):
                    self._same(getattr(E, a), [getattr(q, a) for q in mags], a)
                self._same(E.eff_length, [q.Leff for q in mags], "L_eff")
                self._same(mt.QP * E.mat_B * E.eff_length, [q.theta0 for q in mags], "θ0")
                for col, attr in (("mat_G", "k1"), ("mat_Gv", "k1v"), ("mat_S", "k2"), ("mat_O", "k3")):
                    self._same(mt.QP * getattr(E, col), [getattr(q, attr) for q in mags], attr)
                self.assertEqual(list(np.where(E.n_slices > 0, E.n_slices, mt.Settings().n_slices)),
                                 [q.nsl for q in mags])
                self.assertEqual(list(E.kind), [q.kind for q in mags])
                m = db.machine(lat)
                self.assertEqual([list(m.instances.device_type[E.members(i)]) for i in range(len(E))],
                                 [q.devs for q in mags])
                self.assertEqual(int(E.member_count.sum()), len(m.instances))          # обмотки — разбиение строк
                self.assertTrue(np.array_equal(E.member_start[1:], np.cumsum(E.member_count)[:-1]))

    def test_fingerprint_stable_local_and_sensitive(self):
        import csv, shutil, tempfile
        import synchrotronsim as ss
        lat = 3
        fp = ss.open_db().machine(lat).fingerprint
        self.assertEqual(fp, ss.open_db().machine(lat).fingerprint)
        self.assertEqual(len(fp), 64)
        src = os.path.join(mt.TB, "lattice_configs.csv")
        with open(src, encoding="utf-8", newline="") as fh:
            rows = list(csv.reader(fh, delimiter=";"))
        h = {c: i for i, c in enumerate(rows[0])}
        tmp = tempfile.mkdtemp(prefix="ss_overlay_")
        self.addCleanup(shutil.rmtree, tmp, True)
        os.makedirs(os.path.join(tmp, "tables"))
        dst = os.path.join(tmp, "tables", "lattice_configs.csv")

        def write(mod_lat, col):
            out = [list(r) for r in rows]
            for r in out[1:]:
                if r[h["lattice_id"]] == str(mod_lat) and r[h[col]] not in ("", "0"):
                    r[h[col]] = repr(float(r[h[col]]) * (1 + 1e-9)); break
            with open(dst, "w", encoding="utf-8", newline="") as fh:
                csv.writer(fh, delimiter=";").writerows(out)

        write(2, "x")                                                      # правка ЧУЖОЙ решётки
        self.assertEqual(ss.open_db(overlay=tmp).machine(lat).fingerprint, fp)
        write(lat, "x")                                                    # правка своей — на 1e-9 отн.
        m2 = ss.open_db(overlay=tmp).machine(lat)
        self.assertNotEqual(m2.fingerprint, fp)
        self.assertEqual(m2.source, os.path.join(tmp, "tables"))
        self.assertEqual(len(m2.devices), len(ss.open_db().machine(lat).devices))   # прочие таблицы — основные

    def test_independent_of_cwd(self):
        import tempfile
        import synchrotronsim as ss
        ref = ss.open_db().machine(3)
        cwd = os.getcwd()
        tmp = tempfile.mkdtemp(prefix="ss_cwd_")
        try:
            os.chdir(tmp)
            m = ss.open_db().machine(3)
        finally:
            os.chdir(cwd); os.rmdir(tmp)
        self.assertEqual(m.fingerprint, ref.fingerprint)

    def test_units_devices_immutability(self):
        import dataclasses, csv
        import synchrotronsim as ss
        import facility_report as fr
        db = ss.open_db()
        m = db.machine(27 if 27 in db.lattices() else 26)
        if m.screens:                                                       # мм/мрад таблицы → м/рад модели
            with open(os.path.join(mt.TB, "screens.csv"), encoding="utf-8") as fh:
                raw = {r["screen_id"]: r for r in csv.DictReader(fh, delimiter=";") if r["lattice_id"] == str(m.lattice_id)}
            s = m.screens[0]
            self.assertAlmostEqual(s.width, float(raw[s.screen_id]["width_mm"]) * 1e-3, places=15)
            self.assertAlmostEqual(s.dist, float(raw[s.screen_id]["dist_m"]), places=15)
        ids = set(m.instances.instance_id.tolist())
        for g in m.groups:
            self.assertTrue(set(g.instance_ids) <= ids, g.group_id)
        for dev in m.instances.device_type:
            self.assertIn(dev, m.devices)
        err = db.error_model()
        if "quad" in err:
            with open(os.path.join(mt.TB, "element_errors.csv"), encoding="utf-8") as fh:
                pct = {r["kind"]: float(r["field_pct"]) for r in csv.DictReader(fh, delimiter=";")}
            self.assertAlmostEqual(err["quad"].field_rel, pct["quad"] / 100.0, places=15)
        import synchrotronsim.machines as smach
        self.assertIs(fr.LATTICE_INFO, smach.LATTICE_INFO)                  # один источник паспортных данных
        self.assertEqual(m.meta.design_energy_mev if m.meta else 6000.0, 6000.0)
        with self.assertRaises(ValueError):
            m.elements.mat_G[0] = 1.0                                       # столбцы только для чтения
        with self.assertRaises(dataclasses.FrozenInstanceError):
            m.lattice_id = 1


class TestAddressRegistry(unittest.TestCase):
    """Шаг 4: единый реестр адресов `synchrotronsim.addresses` — грамматика, псевдонимы прежних имён (drive.*,
    actual.*, голые id актуаторов), метаданные (писать ли, цена, владелец), восстановление старых раскладок."""

    def test_canonical_aliases_and_idempotence(self):
        from synchrotronsim import addresses as AD
        cases = {"drive.corr.000.hkick.slew": "act.corr.000.hkick.slew",
                 "actual.corr.000.hkick": "act.corr.000.hkick",
                 "corr.000.vkick": "act.corr.000.vkick", "quad.QF1A.dk": "act.quad.QF1A.dk",
                 "sext.SD1A.dk2": "act.sext.SD1A.dk2", "rf.freq": "act.rf.freq",
                 "actual.rf.freq": "act.rf.freq",
                 "obs.ebpm.x:rms": "obs.ebpm.x:rms", "undulator.und_01.gap": "undulator.und_01.gap"}
        for old, new in cases.items():
            self.assertEqual(AD.canonical(old), new, old)
            self.assertEqual(AD.canonical(new), new, new)               # каноническое не меняется
        p = AD.parse("drive.corr.000.hkick.lag_tau")
        self.assertEqual((p.area, p.obj, p.attr), ("act", "corr.000.hkick", "lag_tau"))
        self.assertEqual(AD.parse("act.quad.QF1A.dk").attr, "")          # реализованное значение — без поля
        self.assertEqual(AD.parse("optics.beta_x@QF1A_03").at, "QF1A_03")
        self.assertEqual(AD.parse("element.SP03.COR05.k1").obj, "SP03.COR05")
        self.assertEqual(AD.parse("obs.ebpm.x:rms").red, "rms")
        self.assertTrue(AD.is_address("rec.orbit_rms.x"))
        for bad in ("sf", "bogus.x", "screen", "obs"):
            self.assertFalse(AD.is_address(bad), bad)

    def test_old_and_new_names_resolve_same(self):
        from control import StateView
        from synchrotronsim import addresses as AD
        st = _seed_client_state(); view = StateView(st)
        st.actual["corr.000.hkick"] = 1.5e-4
        for old, new in (("drive.corr.000.hkick.slew", "act.corr.000.hkick.slew"),
                         ("drive.corr.000.hkick.act", "act.corr.000.hkick.act"),
                         ("actual.corr.000.hkick", "act.corr.000.hkick"),
                         ("corr.000.hkick", "act.corr.000.hkick")):
            self.assertEqual(AD.resolve(view, old), AD.resolve(view, new), old)
            self.assertIsNotNone(AD.resolve(view, new), new)

    def test_meta_writable_owner_cost(self):
        from synchrotronsim import addresses as AD
        self.assertTrue(AD.meta("undulator.und_01.gap_target").writable)
        self.assertFalse(AD.meta("undulator.und_01.K").writable)
        self.assertTrue(AD.meta("drive.corr.000.hkick.slew").writable)
        self.assertFalse(AD.meta("act.corr.000.hkick").writable)          # реализованное — только чтение
        self.assertTrue(AD.meta("loop.fofb_x.gain").writable)
        s = AD.meta("screen.SCR_U01")
        self.assertEqual((s.cost, s.owner, s.widget), ("heavy", "engine", AD.W_SCREEN))
        r = AD.meta("rec.orbit_rms.x")
        self.assertEqual((r.units, r.owner), ("мкм", "recorder"))
        self.assertEqual(AD.meta("rec.rf.df_hz").units, "Гц")

    def test_old_layout_window_is_canonicalized(self):
        """Раскладка со старыми адресами (drive.*, actual.*) открывается: окно получает каноническую форму."""
        try:
            from client import windows as W
        except Exception as e:                                           # нет Qt/pyqtgraph — пропуск
            self.skipTest(f"нет GUI-зависимостей: {e}")
        import os as _os
        _os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        from PySide6.QtWidgets import QApplication
        QApplication.instance() or QApplication([])
        w = W.make_window({"kind": "graph", "params": {"addr": "actual.corr.000.hkick"}})
        self.assertEqual(w.addresses(), ["act.corr.000.hkick"])
        pw = W.make_window({"kind": "point", "params": {"x_addr": "drive.corr.000.hkick.slew",
                                                         "y_addr": "obs.ebpm.x:rms"}})
        self.assertEqual(pw.addresses()[0], "act.corr.000.hkick.slew")


class TestRunAndWatch(unittest.TestCase):
    """Шаг 4: прогон `synchrotronsim.run.Run` (папка с паспортом, ряды, артефакты, индекс runs.csv) и подписка
    `Watch`; WatchRecorder пишет любой адрес из fb_record."""

    def setUp(self):
        import tempfile, shutil
        self.tmp = tempfile.mkdtemp(prefix="ss_runs_")
        os.makedirs(os.path.join(self.tmp, "tables"))
        self.addCleanup(shutil.rmtree, self.tmp, True)

    def test_run_roundtrip_and_provenance(self):
        import json, csv
        from synchrotronsim.run import Run, list_runs
        run = Run.create("test", config={"a": 1, "arr": np.arange(3)}, seed=7, source="unit", root=self.tmp)
        self.assertEqual(run.status, "in_progress")
        self.assertEqual(list_runs(self.tmp)[-1]["status"], "in_progress")   # прерванный — так и останется
        for i in range(5):
            run.record("obs.ebpm.x:rms", i, 1e-6 * (i + 1), t=1e-3 * i)
            run.record("obs.ebpm.x", i, [1e-6, -2e-6, 3e-6 * i])           # векторный ряд
        run.record_rows([("orbit_rms.x", 0, 3.5), ("orbit_rms.x", 1, 2.5)], t_of_step={0: 0.0, 1: 0.01})
        run.record("drive.corr.000.hkick.slew", 0, 1.5)                    # псевдоним → каноническое имя
        run.artifact("bpm", {"x": np.ones(4), "spec": {"E": np.arange(2)}, "none": None}, meta="тест")
        run.finish("done")
        r2 = Run.open(run.run_id, root=self.tmp)
        self.assertEqual(r2.status, "done")
        st, t, v = r2.series("obs.ebpm.x:rms")
        np.testing.assert_array_equal(st, np.arange(5)); np.testing.assert_allclose(v, 1e-6 * np.arange(1, 6))
        np.testing.assert_allclose(t, 1e-3 * np.arange(5))
        self.assertEqual(r2.series("obs.ebpm.x")[2].shape, (5, 3))
        self.assertEqual(list(r2.series("rec.orbit_rms.x")[2]), [3.5, 2.5])
        self.assertIn("act.corr.000.hkick.slew", r2.names())
        self.assertEqual(set(r2.load_artifact("bpm")), {"x", "spec.E"})
        man = json.load(open(os.path.join(r2.path, "manifest.json"), encoding="utf-8"))
        self.assertEqual(man["config"], {"a": 1, "arr": [0, 1, 2]})
        self.assertIn("sha", man["provenance"]["git"])
        self.assertIn("python", man["provenance"])
        self.assertEqual(man["signals"]["rec.orbit_rms.x"]["units"], "мкм")
        with open(os.path.join(r2.path, "tables", "signals.csv"), encoding="utf-8-sig") as f:
            rows = list(csv.reader(f, delimiter=";"))
        self.assertIn("2,5", [r[5] for r in rows[1:]])                     # человеку — десятичная запятая
        self.assertFalse(any(r[0] == "obs.ebpm.x" for r in rows[1:]))      # векторы — только в signals.npz
        idx = list_runs(self.tmp)
        self.assertEqual([r["run_id"] for r in idx].count(run.run_id), 1)
        self.assertEqual(idx[-1]["status"], "done")

    def test_watch_validation(self):
        from synchrotronsim.run import Watch
        w = Watch("actual.corr.000.hkick")
        self.assertEqual(w.addr, "act.corr.000.hkick")
        self.assertEqual(w.meta.units, "рад")
        for bad in (dict(addr="bogus.x"), dict(addr="obs.tune", cadence="hourly"), dict(addr="obs.tune", every=0)):
            with self.assertRaises(ValueError):
                Watch(**bad)
        self.assertEqual(Watch("screen", recorder="screen", params="SCR_U19").addr, "screen")   # рекордер

    def test_watch_recorder_from_fb_record(self):
        """fb_record: signal-адрес → WatchRecorder (ряд под каноническим адресом); имя рекордера → как раньше."""
        from control import fb_config, recorders
        R = fb_config.Record
        recs = recorders.build([R("orbit_rms", "", "per_window", True),
                                R("obs.ebpm.x:rms", "every=2", "per_window", True),
                                R("actual.corr.000.hkick", "", "final", True),
                                R("sf", "", "post", True)])
        self.assertEqual([type(r).__name__ for r in recs], ["OrbitRmsRecorder", "WatchRecorder", "WatchRecorder"])
        st = _seed_client_state(); st.actual["corr.000.hkick"] = 2e-4
        for w in range(4):
            for r in recs[1:]:
                r.sample(st, w, None)
        for r in recs[1:]:
            r.finalize(st, None)
        self.assertEqual([s for s, _, _ in recs[1].series_rows()], ["obs.ebpm.x:rms"] * 2)   # каждое второе окно
        self.assertEqual([st_ for _, st_, _ in recs[1].series_rows()], [0, 2])
        self.assertEqual(recs[2].result_rows(), [("act.corr.000.hkick", 0, 2e-4)])
        self.assertEqual(recs[0].prefix, "rec."); self.assertEqual(recs[1].prefix, "")


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "fb_experiments.csv")), "нет таблиц fb_*")
class TestExperimentAndSimulation(unittest.TestCase):
    """Шаг 5: опыт-объект `synchrotronsim.experiment.Experiment` (из fb_*, кокпит, JSON, пакет), единый сборщик
    `synchrotronsim.simulation` (контуры по LoopSpec, темпы, подписка), запись прогона и лимит в движке кокпита."""

    def test_from_fb_matches_fb_config_and_json_roundtrip(self):
        import tempfile, shutil
        from control import fb_config
        from synchrotronsim.experiment import Experiment
        cfg = fb_config.load(3)
        e = Experiment.from_fb(3)
        self.assertEqual((e.lattice_id, e.seed, e.source, e.plant.kind), (cfg.lattice_id, cfg.seed, "fb_exp_id=3", "tbt"))
        self.assertEqual(e.plant.kwargs["settle_turns"], cfg.settle)
        self.assertEqual([l.name for l in e.loops], [l.name for l in cfg.enabled_loops()])
        self.assertEqual(e.run.steps, cfg.steps)
        sm = Experiment.from_fb(3, smoke=True)
        self.assertEqual((sm.run.steps, sm.plant.kwargs["settle_turns"], sm.run.record_mesh),
                         (min(cfg.steps, 36), min(cfg.settle, 120), 81))
        self.assertIn("rf.delta", e.plant.extra_sensors)                 # подписка rf без контура ВЧ → датчик δ̂
        self.assertEqual({w.recorder for w in e.record}, {"orbit_rms", "rf"})
        tmp = tempfile.mkdtemp(prefix="ss_exp_"); self.addCleanup(shutil.rmtree, tmp, True)
        p = e.to_json(os.path.join(tmp, "e.json"))
        self.assertEqual(Experiment.from_json(p).to_dict(), e.to_dict())
        o = Experiment.from_fb(2)                                         # оптический опыт
        self.assertEqual((o.plant.kind, o.plant.kwargs["observe"]), ("optics", ("tune", "chroma", "adts")))
        self.assertEqual(Experiment.from_json(o.to_json(os.path.join(tmp, "o.json"))).plant.kwargs["observe"],
                         ("tune", "chroma", "adts"))                       # кортеж переживает JSON

    def test_cockpit_experiment_and_loops(self):
        from synchrotronsim.experiment import Experiment
        from synchrotronsim.simulation import make_loop, loop_rate_hz
        e = Experiment.cockpit(optics=True, rf=True, n_part=4)
        self.assertEqual([l.name for l in e.loops], ["fofb_x", "fofb_y", "tune", "chroma", "adts", "rf"])
        self.assertFalse(any(l.enabled for l in e.loops))                 # контуры кокпита — ВЫКЛ на старте
        self.assertEqual(e.plant.kind, "composite"); self.assertEqual(e.plant.kwargs["n_part"], 4)
        classes = [type(make_loop(l)).__name__ for l in e.loops]
        self.assertEqual(classes, ["IntegralSVD", "IntegralSVD", "TuneLoop", "ChromaLoop", "AdtsLoop", "RfFrequencyLoop"])
        t = make_loop(e.loops[2])
        self.assertEqual((t.sensor, t.actuator_key, t.ref, t.params["gain"]), ("tune", "quad.", "tune.target", 0.7))  # умолчания класса
        dt, window = 1e-4, 9
        self.assertAlmostEqual(loop_rate_hz(e.loops[0], dt, window, e.run), 1.0 / (dt * e.run.steps_per_cycle))
        self.assertAlmostEqual(loop_rate_hz(e.loops[2], dt, window, e.run), 1.0 / (dt * e.run.optics_period_ticks))
        f = Experiment.from_fb(3)
        self.assertAlmostEqual(loop_rate_hz(f.loops[0], dt, window, f.run), 1.0 / (dt * window))   # auto = окно
        self.assertEqual(make_loop(f.loops[0]).ref, None)                 # пустая ячейка ref → без уставки

    def test_bundle_reproduces_machine(self):
        import tempfile, shutil, zipfile
        from synchrotronsim.experiment import Experiment
        from synchrotronsim import open_db
        e = Experiment.cockpit(lattice=3)
        tmp = tempfile.mkdtemp(prefix="ss_bundle_"); self.addCleanup(shutil.rmtree, tmp, True)
        path = e.save_bundle(os.path.join(tmp, "b.zip"))
        names = zipfile.ZipFile(path).namelist()
        self.assertIn("tables/lattice_configs.csv", names); self.assertIn("experiment.json", names)
        e2, db2 = Experiment.load_bundle(path, os.path.join(tmp, "ws"))
        self.assertEqual(e2.to_dict(), e.to_dict())
        self.assertEqual(db2.machine(3).fingerprint, open_db().machine(3).fingerprint)
        self.assertEqual(db2.lattices(), [3])                             # в пакете — только своя машина

    def test_engine_records_run_and_turn_limit(self):
        """Кокпит: адреса «● запись» → Watch движка; запись → прогон на каждый проход; лимит оборотов
        завершает проход (done) и, при «повторять», начинает новый с нуля."""
        import tempfile, shutil, time as _t
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler
        from synchrotronsim.run import Run, list_runs
        tmp = tempfile.mkdtemp(prefix="ss_cockpit_runs_"); self.addCleanup(shutil.rmtree, tmp, True)
        os.makedirs(os.path.join(tmp, "tables"))
        plant = _FakeInfinitePlant()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        plant.reset(state)
        link = InProcLink()
        eng = EngineThread(link, sch, steps_per_cycle=1, publish_every=1, idle_sleep=0.002, run_root=tmp)
        link.watch("obs.ebpm.x:max"); link.watch("undulator.und_01.gap"); link.watch("bogus.x")
        link.record(True, max_turns=20, repeat=True)
        link.control("run")
        eng.start()
        try:
            end = _t.time() + 10
            while _t.time() < end and len([r for r in list_runs(tmp) if r["status"] == "done"]) < 2:
                _t.sleep(0.01)
            link.record(False)
            _t.sleep(0.1)
        finally:
            eng.stop(); eng.join(timeout=5)
        rows = list_runs(tmp)
        done = [r for r in rows if r["status"] == "done"]
        self.assertGreaterEqual(len(done), 2, rows)                       # два полных прохода по лимиту
        self.assertTrue(all(r["status"] in ("done", "stopped") for r in rows))
        r0 = Run.open(done[0]["run_id"], root=tmp)
        self.assertEqual(r0.names(), ["obs.ebpm.x:max", "undulator.und_01.gap"])   # bogus.x отвергнут
        st, t, v = r0.series("obs.ebpm.x:max")
        self.assertEqual(int(st[-1]), 20)                                 # проход ровно до лимита (оборот = такт)
        self.assertAlmostEqual(float(v[-1]), 2.0e-5)                      # max(1e-6·оборот, …) на 20-м обороте
        self.assertEqual(r0.manifest["source"], "cockpit")
        self.assertEqual(r0.manifest["config"]["max_turns"], 20)          # без опыта — сводка движка

    def test_cockpit_record_ui(self):
        """Кнопка «● запись» окна графика шлёт Watch; блок «запись прогона» панели шлёт record(on, лимит, повтор)
        и показывает текущий прогон из меты движка."""
        try:
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            from PySide6.QtWidgets import QApplication
            from client import windows as W
            from client.cockpit import ControlPanel
        except Exception as e:
            self.skipTest(f"нет GUI-зависимостей: {e}")
        from client import InProcLink
        from client.link import WATCH, UNWATCH, RECORD
        QApplication.instance() or QApplication([])
        got = []
        w = W.GraphWindow("obs.ebpm.x:rms")
        w.request_watch.connect(lambda a, on: got.append((a, on)))
        w.btn_watch.setChecked(True); w.btn_watch.setChecked(False)
        self.assertEqual(got, [("obs.ebpm.x:rms", True), ("obs.ebpm.x:rms", False)])
        link = InProcLink()
        p = ControlPanel(link)
        p.chk_rec.setChecked(True); p.spin_turns.setValue(5000); p.chk_repeat.setChecked(True)
        p._send_record()
        msgs = [m for m in link.pull() if m.kind == RECORD]
        self.assertEqual(msgs[-1].value, {"on": True, "max_turns": 5000, "repeat": True})
        p.update_meta({"run_id": "20261125-x", "watches": ["obs.ebpm.x:rms"], "pass": 2})
        self.assertIn("20261125-x", p.lbl_run.text())
        link.watch("obs.tune"); link.unwatch("obs.tune")
        self.assertEqual([m.kind for m in link.pull()], [WATCH, UNWATCH])

    def test_loop_gain_address_reads_and_commands(self):
        """loop.<имя>.gain — читает и задаёт Controller.params["gain"] (раньше команда молча не применялась)."""
        from client import InProcLink
        from client.engine import EngineThread
        from control import ExperimentState, Scheduler, StateView
        from control.controllers.integral_svd import IntegralSVD
        from synchrotronsim import addresses as AD
        from control.controller import ResponseMatrix

        class _P(_FakeInfinitePlant):
            def build_response(self, actuators, observables):
                return ResponseMatrix(np.ones((3, len(actuators))), ["b0", "b1", "b2"], actuators)
        plant = _P()
        state = ExperimentState(plant.schema())
        sch = Scheduler(plant, state, dt=1.0e-3)
        sch.add(IntegralSVD("fofb_x", sensor="ebpm.x", actuator_key="hkick", gain=0.5), rate_hz=100.0, enabled=False)
        self.assertEqual(AD.resolve(StateView(state), "loop.fofb_x.gain", loops=sch.loops), 0.5)
        eng = EngineThread(InProcLink(), sch)
        eng._command("loop.fofb_x.gain", 0.25)
        self.assertEqual(sch.loops["fofb_x"].controller.params["gain"], 0.25)
        self.assertIsNone(eng._unknown)


class TestCockpitOnExperiment(unittest.TestCase):
    """Клиент на библиотеке: кокпит собирает модель ПО ОПЫТУ (`cockpit_experiment` / JSON / fb:N), окно «Опыт»
    правит точность и контуры и даёт новый опыт, пересборка переносит отметки записи, лимит, зазоры и пуск."""

    def _qt(self):
        try:
            os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
            from PySide6.QtWidgets import QApplication
        except Exception as e:
            self.skipTest(f"нет Qt: {e}")
        return QApplication.instance() or QApplication([])

    def test_cockpit_experiment_matches_previous_launch_and_roundtrips(self):
        try:
            from client.cockpit import cockpit_experiment, START_ON
        except Exception as e:
            self.skipTest(f"нет GUI-зависимостей: {e}")
        from synchrotronsim.experiment import Experiment
        e, note = cockpit_experiment(lattice=26, settle_turns=200, nsl=6, n_part=5, stride=10, nsl_dev=None)
        old = Experiment.cockpit(lattice=26, settle_turns=200, steps_per_cycle=9, nsl=6, n_part=5, stride=10,
                                 relative_orbit=True, optics=True, rf=True)          # прежние аргументы сборщика
        for l in old.loops:
            l.enabled = l.name in START_ON
        self.assertEqual(e.to_dict(), old.to_dict())
        self.assertEqual({l.name for l in e.loops if l.enabled}, {"fofb_x", "fofb_y", "tune", "chroma"})
        off, _ = cockpit_experiment(fofb_on=False, optics=False, nsl_dev={"QF1": 3})
        self.assertFalse(any(l.enabled for l in off.loops))
        self.assertEqual(off.plant.kwargs["nsl_dev"], {"QF1": 3})
        for x in (e, off):                                             # обратное к Experiment.cockpit
            y = Experiment.cockpit(**x.cockpit_params())
            for a, b in zip(y.loops, x.loops):
                a.enabled = b.enabled
            self.assertEqual(y.to_dict(), x.to_dict())

    def test_load_experiment_json_and_fb(self):
        import tempfile, shutil
        try:
            from client.cockpit import load_experiment, cockpit_experiment
        except Exception as e:
            self.skipTest(f"нет GUI-зависимостей: {e}")
        e, _ = cockpit_experiment(nsl_dev=None)
        tmp = tempfile.mkdtemp(prefix="ss_cexp_"); self.addCleanup(shutil.rmtree, tmp, True)
        p = e.to_json(os.path.join(tmp, "c.json"))
        self.assertEqual(load_experiment(p).to_dict(), e.to_dict())
        self.assertIs(load_experiment(e), e)
        self.assertEqual(load_experiment("fb:3").source, "fb_exp_id=3")

    def test_experiment_dialog_form(self):
        self._qt()
        from client.experiment_dialog import ExperimentDialog, slices_for
        from client.cockpit import cockpit_experiment
        from synchrotronsim.experiment import Experiment
        e, _ = cockpit_experiment(optics=False, nsl_dev=None)
        d = ExperimentDialog(e)
        self.assertEqual(d.experiment_from_form().to_dict(), e.to_dict())   # без правок — тот же опыт
        d.sp_nsl.setValue(4); d.sp_part.setValue(3)
        d.loop_chk["fofb_y"].setChecked(False)
        d.chk_optics.setChecked(True)                                  # вид модели: + оптические контуры
        n = d.experiment_from_form()
        self.assertEqual((n.plant.kind, n.plant.kwargs["nsl"], n.plant.kwargs["n_part"]), ("composite", 4, 3))
        en = {l.name: l.enabled for l in n.loops}
        self.assertEqual(en, {"fofb_x": True, "fofb_y": False, "tune": False, "chroma": False, "adts": False,
                              "rf": False})                            # новые контуры — ВЫКЛ, как в Experiment.cockpit
        self.assertEqual((n.name, n.source, n.run.steps_per_cycle), (e.name, "cockpit", e.run.steps_per_cycle))
        # опыт из БД: вид модели не меняется, в аргументы plant'а пишутся только правки
        f = Experiment.from_fb(3)
        df = ExperimentDialog(f)
        self.assertFalse(df.chk_optics.isEnabled())
        self.assertEqual(df.experiment_from_form().to_dict(), f.to_dict())
        self.assertNotIn("n_part", f.plant.kwargs)
        df.sp_part.setValue(5)
        self.assertEqual(df.experiment_from_form().plant.kwargs["n_part"], 5)
        # срезы: «как в опыте» на другой решётке — набор новой решётки (набор прежней к ней не относится)
        self.assertEqual(slices_for("keep", 26, {"X": 2})[0], {"X": 2})
        from synchrotronsim.convergence import cockpit_slices
        self.assertEqual(slices_for("keep", 3, {"X": 2}, same_lattice=False)[0], cockpit_slices(3)[0])
        self.assertIsNone(slices_for("none", 26, {"X": 2})[0])

    def test_attach_and_rebuild_carry_over(self):
        """Сборка → движок; пересборка (detach + attach) переносит отметки «● запись», лимит, зазоры, пуск;
        подписки панели не копятся; опыт без кольца (optics) кокпит не принимает."""
        import time as _t
        from types import SimpleNamespace
        self._qt()
        try:
            from client.cockpit import Cockpit
            import client.scene3d as s3
        except Exception as e:
            self.skipTest(f"нет GUI-зависимостей: {e}")
        from client import InProcLink
        from control import ExperimentState, Scheduler
        from synchrotronsim.experiment import Experiment, PlantSpec, RunLimits

        class _Scene:                                                  # без VTK: сцена-заглушка
            model = SimpleNamespace(screens=[])
            n = 0
            def set_geometry(self, b): self.n += 1
            def fit_view(self): pass
            def update_snapshot(self, s): pass
            def animate_jaws(self): pass
            def is_interacting(self): return False

        orig = s3.build_geometry
        s3.build_geometry = lambda plant, state=None, **k: {}
        self.addCleanup(setattr, s3, "build_geometry", orig)

        def session(name, max_turns=None):
            plant = _FakeInfinitePlant()
            state = ExperimentState(plant.schema())
            plant.reset(state)
            sch = Scheduler(plant, state, dt=1.0e-3)
            sch.experiment = Experiment(name=name, plant=PlantSpec("tbt", {}), run=RunLimits(max_turns=max_turns))
            return plant, state, sch

        link = InProcLink()
        win = Cockpit(link, project="unittest_tmp_cockpit", poll_ms=100000)
        win.scene = _Scene()
        try:
            win._attach(session("t1", max_turns=777))
            e1 = win._engine
            self.assertEqual(win.experiment.name, "t1")
            self.assertEqual(e1.runner.max_turns, 777)                 # лимит — из опыта
            self.assertEqual(win.control.spin_turns.value(), 777)
            link.watch("obs.ebpm.x:max")
            link.command("undulator.und_01.gap", 9.0)
            link.control("run")
            end = _t.time() + 5
            while _t.time() < end and not (e1.running and "obs.ebpm.x:max" in e1.runner.watches
                                           and abs(win._state.undulators["und_01"].gap_mm - 9.0) < 1e-9):
                _t.sleep(0.005)
            cur = win.current_experiment()
            self.assertEqual([w.addr for w in cur.record], ["obs.ebpm.x:max"])
            self.assertEqual(cur.run.max_turns, 777)
            n_subs = len(link.subscriptions())

            win._carry = win._detach()                                 # = rebuild без фоновой сборки
            self.assertFalse(e1.is_alive())
            win._attach(session("t2"))
            e2 = win._engine
            self.assertIsNot(e2, e1)
            self.assertEqual(win.experiment.name, "t2")
            self.assertIn("obs.ebpm.x:max", e2.runner.watches)
            self.assertEqual(e2.runner.max_turns, 777)                 # лимит прежней сборки, не опыта t2
            self.assertTrue(e2._running)
            self.assertAlmostEqual(win._state.undulators["und_01"].gap_mm, 9.0)
            self.assertEqual(len(link.subscriptions()), n_subs)        # подписки панели не копятся
            self.assertEqual(win.scene.n, 2)
            with self.assertRaises(ValueError):                        # optics — без кольца
                win.build(Experiment(name="o", plant=PlantSpec("optics", {})))
        finally:
            win._detach()
            win.close()


@unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")), "нет БД")
class TestDocumentsAndStore(unittest.TestCase):
    """Протокол данных, этап 1 (docs/DATA_PROTOCOL.md): документ машины (конверт + столбцы) сохраняет отпечаток и
    ядро бит-в-бит; хранилище фиксирует ревизии снимками, опыты — ревизиями с проверкой родителя (конфликт)."""

    def test_machine_document_roundtrip(self):
        from synchrotronsim import open_db
        from synchrotronsim.documents import machine_to_bytes, machine_from_bytes, unpack, pack, DocumentError
        db = open_db()
        for lat in (2, 26, 27):
            m = db.machine(lat)
            b = machine_to_bytes(m)
            m2 = machine_from_bytes(b)
            self.assertEqual(m2.fingerprint, m.fingerprint)
            for tn in ("instances", "elements"):
                A, B = getattr(m, tn), getattr(m2, tn)
                for c in A.column_names():
                    a, bb = getattr(A, c), getattr(B, c)
                    self.assertEqual(a.dtype, bb.dtype, (lat, tn, c))
                    if a.dtype.kind == "f":
                        self.assertTrue(np.array_equal(a, bb, equal_nan=True) and
                                        np.array_equal(np.signbit(a), np.signbit(bb)), (lat, tn, c))
                    else:
                        self.assertTrue((a == bb).all(), (lat, tn, c))
            self.assertEqual((m2.groups, m2.meta, m2.screens, m2.cavities), (m.groups, m.meta, m.screens, m.cavities))
        env, arr = unpack(b)
        self.assertEqual((env["format"], env["kind"], env["id"]), ("synchrotronsim.machine/1", "machine", "27"))
        env2 = dict(env, content=dict(env["content"], lattice_id=99))          # содержимое поправлено, отпечаток — нет
        with self.assertRaises(DocumentError):
            machine_from_bytes(pack(env2, arr))
        with self.assertRaises(DocumentError):                                 # старшая версия формата новее
            machine_from_bytes(pack(dict(env, format="synchrotronsim.machine/2"), arr))
        if _engine_ext_available():                                            # ядро по документу — бит-в-бит
            from synchrotronsim.engine import ExtRing
            m = db.machine(26); m2 = machine_from_bytes(machine_to_bytes(m))
            rf = {"volt_mv": 8.0, "harmonic": 1296, "phase_deg": 165.0}
            x0 = np.array([[1e-4, 0, 5e-5, 0, 0, 1e-3], [2e-4, 1e-5, -5e-5, 0, 0, 0]])
            ta = np.asarray(ExtRing(machine=m, sr=True, rf=rf).track(x0, 5))
            tb = np.asarray(ExtRing(machine=m2, sr=True, rf=rf).track(x0, 5))
            self.assertTrue(np.array_equal(ta, tb))

    def test_store_pins_revisions_and_experiment_conflicts(self):
        import tempfile, shutil, dataclasses
        from synchrotronsim.store import LocalCsvStore, Ref, StoreError, ConflictError
        from synchrotronsim.experiment import Experiment
        tmp = tempfile.mkdtemp(prefix="ss_store_"); self.addCleanup(shutil.rmtree, tmp, True)
        st = LocalCsvStore(store_dir=tmp)
        head = st.get("machine", "2")
        self.assertIs(st.get("machine", "2", head.fingerprint[:10]), head)    # префикс ревизии = последняя
        old = dataclasses.replace(head, meta=dataclasses.replace(head.meta, label="старая ревизия"), fingerprint="")
        self.assertNotEqual(old.fingerprint, head.fingerprint)
        ref = st.pin(old)
        self.assertEqual((ref.kind, ref.id, ref.rev), ("machine", "2", old.fingerprint))
        self.assertEqual(st.pin(old), ref)                                    # повтор — та же ревизия, без записи
        got = st.resolve(str(ref))
        self.assertEqual((got.fingerprint, got.meta.label), (old.fingerprint, "старая ревизия"))
        self.assertEqual(len(st.history("machine", 2)), 1)
        with self.assertRaises(StoreError):
            st.get("machine", "2", "0123456789ab")
        with self.assertRaises(StoreError):                                   # машины пишутся ревизиями на этапе 2
            st.commit("machine", "2", head)
        self.assertEqual(Ref.parse("machine:26@abc"), Ref("machine", "26", "abc"))
        self.assertEqual(str(Ref("experiment", "x")), "experiment:x")
        e = Experiment.cockpit(lattice=2)
        r1 = st.commit("experiment", "demo", e, note="первая")
        with self.assertRaises(ConflictError):                                # поверх «ничего», а ревизия уже есть
            st.commit("experiment", "demo", dataclasses.replace(e, note="b"))
        r2 = st.commit("experiment", "demo", dataclasses.replace(e, note="b"), parent=r1.rev)
        self.assertNotEqual(r1.rev, r2.rev)
        self.assertEqual(st.get("experiment", "demo").note, "b")
        self.assertEqual(st.get("experiment", "demo", r1.rev).to_dict(), e.to_dict())
        self.assertEqual([h["rev"] for h in st.history("experiment", "demo")], [r1.rev, r2.rev])
        self.assertEqual(st.get("experiment", "fb:3").source, "fb_exp_id=3")
        self.assertEqual(e.machine_ref(), Ref("machine", "2", ""))
        h = st.blob_put(b"abc")
        self.assertEqual(st.blob_get(h), b"abc")


@unittest.skipUnless(_engine_ext_available() and os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет собранного расширения или БД")
class TestModelFromDocuments(unittest.TestCase):
    """Критерий этапа 1 протокола данных: документы опыта и машины загружены заранее — дальше пакетный прогон,
    сборка и шаги кокпита, экраны (ондулятор и линейка диполей) и сцена НЕ открывают ни одной таблицы. Паспорт
    прогона ссылается на машину с ревизией; пакет опыта считается без БД."""

    ALLOW = ("workspace/results/", "workspace/store/", "workspace/runs/", "workspace/fields_bin/")

    def _guard(self):
        import builtins, io as _io
        root = os.path.dirname(os.path.abspath(mt.__file__))
        seen = []
        o_open, o_io = builtins.open, _io.open

        def h(file, mode="r", *a, **k):
            if isinstance(file, (str, os.PathLike)):
                ap = os.path.abspath(os.fspath(file))
                if ap.lower().startswith(root.lower()) and not ap.endswith((".py", ".pyc", ".pyd")):
                    rel = os.path.relpath(ap, root).replace("\\", "/")
                    if not rel.startswith(self.ALLOW):
                        seen.append(rel)
            return o_open(file, mode, *a, **k)
        builtins.open = h; _io.open = h

        def restore():
            builtins.open = o_open; _io.open = o_io
        return seen, restore

    def test_compute_path_opens_no_tables(self):
        import tempfile, shutil
        sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(mt.__file__)), "undulators"))
        from synchrotronsim.experiment import Experiment
        from synchrotronsim.simulation import Simulation, build_session
        from synchrotronsim.store import open_store
        from client.cockpit import cockpit_experiment
        from client import InProcLink
        from client.engine import EngineThread
        from client.scene3d import build_geometry
        import screen_analysis as sa
        import control.recorders, control.tbt_plant  # noqa: F401 — импорты до проверки
        tmp = tempfile.mkdtemp(prefix="ss_noread_"); self.addCleanup(shutil.rmtree, tmp, True)
        os.makedirs(os.path.join(tmp, "tables"))
        exp = Experiment.from_fb(3, smoke=True)
        exp.run.steps = 9
        m = open_store().machine(exp.lattice_id)
        ce, _ = cockpit_experiment(optics=False, nsl_dev=None)
        bend_sid = sa.bend_line_screens(m)[0]["sid"]
        seen, restore = self._guard()
        try:
            run = Simulation(exp, root=tmp, legacy_fb=False, machine=m, progress=None).run()
            S = build_session(ce, machine=m)
            for _ in range(2):
                S.sch.step()
            eng = EngineThread(InProcLink(), S.sch)
            g = build_geometry(S.plant, S.state)
            und_sid = next(x["sid"] for x in g["screens"] if x["kind"] == "undulator")
            res = {(sid, q): eng._compute_screen(f"screen.{sid}", {"res": 21, "quantity": q})
                   for sid in (und_sid, bend_sid) for q in ("flux", "power")}
            S.close()
        finally:
            restore()
        self.assertEqual(seen, [], f"модель открыла таблицы: {sorted(set(seen))}")
        self.assertTrue(all(not r.get("error") and r.get("heatmap") for r in res.values()),
                        {k: r.get("error") for k, r in res.items()})
        ref = run.manifest["machine"]["ref"]
        self.assertEqual(ref, f"machine:26@{m.fingerprint}")                  # паспорт: машина с ревизией
        self.assertTrue(os.path.exists(os.path.join(tmp, "store", "machine", "26", f"{m.fingerprint}.ssd")))
        self.assertEqual(und_sid, "SCR_und_01")                               # экраны — из документа машины 26

    def test_bundle_runs_without_db_and_rev_mismatch_rejected(self):
        import tempfile, shutil, dataclasses
        from synchrotronsim.experiment import Experiment
        from synchrotronsim.simulation import build_session, resolve_machine
        tmp = tempfile.mkdtemp(prefix="ss_bundle2_"); self.addCleanup(shutil.rmtree, tmp, True)
        e = Experiment.cockpit(lattice=26, settle_turns=20)
        p = e.save_bundle(os.path.join(tmp, "b.zip"))
        e2, m = Experiment.open_bundle(p)
        self.assertEqual(e2.to_dict(), e.to_dict())
        self.assertTrue(m.source.startswith("document:"))
        seen, restore = self._guard()
        try:
            S = build_session(e2, machine=m)
            S.sch.step()
            S.close()
        finally:
            restore()
        self.assertEqual(seen, [], f"пакет: модель открыла таблицы: {sorted(set(seen))}")
        with self.assertRaises(ValueError):                                   # опыт на другой ревизии
            resolve_machine(dataclasses.replace(e2, machine_rev="0123456789ab"), machine=m)
        with self.assertRaises(ValueError):                                   # опыт на другой решётке
            resolve_machine(dataclasses.replace(e2, lattice_id=27), machine=m)


@unittest.skipUnless(_engine_ext_available() and os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                     "нет собранного расширения или БД")
class TestParallelColumnsAndClosedOrbit(unittest.TestCase):
    """Шаг 6: замкнутая орбита останавливается у шума округления (floor); копии кольца (`CppRing.clone`,
    потоки на кольцо) и параллельные столбцы отклика — бит-в-бит как последовательно."""

    LAT, MEV = 2, 3000.0                                    # Сибирь — маленькая решётка, быстро

    def _ring(self, **kw):
        import ring_cpp as rc
        r = rc.CppRing(lattice=self.LAT, design_mev=self.MEV, nsl=5, **kw)
        self.addCleanup(r.close)
        return r

    def test_closed_orbit_stops_at_noise_floor(self):
        """СИЛА-26 (4882 магнита): шум округления орбиты ~1e-12 м — прежний критерий tol=1e-13 недостижим, Ньютон
        шёл все 10 итераций; с floor — 4–5, то же решение в пределах шума."""
        import ring_analysis as ra
        import ring_cpp as rc
        r = rc.CppRing(lattice=26, design_mev=6000.0, nsl=5)
        self.addCleanup(r.close)
        kk = [j for j, k in enumerate(r.kinds) if k == "kicker"][:1]
        r.apply_perturb(kicks=[(kk[0], 2e-6)])
        calls = {"n": 0}
        orig = r.track

        def counted(*a, **k):
            calls["n"] += 1
            return orig(*a, **k)
        r.track = counted
        co = ra.closed_orbit(r, dims=4)
        n_new = calls["n"]; calls["n"] = 0
        co_old = ra.closed_orbit(r, dims=4, floor=0.0)           # прежнее поведение: все max_iter итераций
        self.assertEqual(calls["n"], 10)
        self.assertLess(n_new, 10)
        r.track = orig
        res = np.max(np.abs(r.track([co])[0][:4] - co[:4]))
        self.assertLess(res, ra.CO_FLOOR * 10)                   # замкнута до шума
        self.assertLess(np.max(np.abs(co - co_old)), 1e-10)      # то же решение (в пределах шума)

    def test_clone_and_threads(self):
        r = self._ring()
        q = [j for j, k in enumerate(r.kinds) if k == "quad"]
        r.apply_perturb(dks=[(q[0], 1e-3)])
        c = r.clone().set_threads(2)
        self.addCleanup(c.close)
        st = np.zeros((3, 6)); st[:, 0] = [1e-4, 0, -1e-4]
        a = r.track_tbt(st, 3)
        self.assertTrue(np.array_equal(c.track_tbt(st, 3).view(np.uint64), a.view(np.uint64)))   # то же состояние
        c.apply_perturb()                                        # RESET копии не трогает оригинал
        self.assertFalse(np.array_equal(c.track_tbt(st, 3), a))
        self.assertTrue(np.array_equal(r.track_tbt(st, 3).view(np.uint64), a.view(np.uint64)))
        self.assertEqual(c._be.x._r.threads, 2)

    def test_parallel_response_bitwise_equals_serial(self):
        from control import response as resp
        import orbit_correct as oc
        r = self._ring()
        kick = [j for j, k in enumerate(r.kinds) if k in ("kicker", "quad")][:6]
        mon = np.arange(1, r.n_mag + 1)
        applier = lambda ring, j, eps: ring.apply_perturb(kicks=[(j, eps)])
        reader = lambda o, orb: orb[mon, 0]
        R1 = resp.build_response(r, kick, ["x"], applier, reader, workers=1)["x"]
        R4 = resp.build_response(r, kick, ["x"], applier, reader, workers=4)["x"]
        self.assertTrue(np.array_equal(R1.view(np.uint64), R4.view(np.uint64)))
        os.environ["SS_RESP_WORKERS"] = "1"
        try:
            m1 = oc.response_matrix(r, corr=kick[:3])
        finally:
            os.environ.pop("SS_RESP_WORKERS", None)
        m3 = oc.response_matrix(r, corr=kick[:3])                # авто: параллельно
        for k in ("Rx", "Ry"):
            self.assertTrue(np.array_equal(m1[k].view(np.uint64), m3[k].view(np.uint64)), k)


class TestSliceConvergence(unittest.TestCase):
    """Шаг 6: инструмент сходимости срезов — смена срезов элемента на лету (= сборка с nsl_dev, бит-в-бит) и
    жадный подбор по допуску (уложиться при наименьшей цене элементы × срезы)."""

    def test_recommend_greedy_meets_tolerance_cheaply(self):
        from synchrotronsim.convergence import recommend
        cand = [1, 2, 4, 8]
        # тип A: много элементов, ошибка ∝ 1/n⁴; тип B: мало элементов, крупная ошибка
        errs = {"A": {n: (1e-3 / n ** 4, 0.0, 1e-7 / n ** 4) for n in cand},
                "B": {n: (8e-3 / n ** 4, 0.0, 8e-7 / n ** 4) for n in cand}}
        S = {"ref": 40, "candidates": cand, "errors": errs, "count": {"A": 1000, "B": 10}}
        rec = recommend(S, tol_nu=1e-3, tol_map=1e-6)
        tot = sum(max(errs[d][rec[d]][0], errs[d][rec[d]][1]) for d in rec)
        self.assertLessEqual(tot, 1e-3)
        self.assertGreater(rec["B"], rec["A"])                         # дешёвый тип поднимается первым
        self.assertEqual(recommend(S, tol_nu=1.0, tol_map=1.0), {"A": 1, "B": 1})   # допуск широк — минимум
        # немонотонная кривая (как QF8D: карта на 2 срезах хуже, чем на 1) не должна застревать на 1 срезе
        errs["C"] = {1: (0.4, 0.0, 3.7e-4), 2: (0.2, 0.0, 9.5e-4), 4: (1e-2, 0.0, 7e-5), 8: (7e-4, 0.0, 4e-6)}
        S["count"]["C"] = 40
        rec = recommend(S, tol_nu=1e-3, tol_map=1e-5)
        self.assertGreaterEqual(rec["C"], 8)

    @unittest.skipUnless(_engine_ext_available() and os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")),
                         "нет расширения или БД")
    def test_set_nsl_equals_nsl_dev_build(self):
        from synchrotronsim.engine import ExtRing
        a = ExtRing(lattice=2, design_mev=3000.0, nsl=10)
        dev = next(nm for nm, k in zip(a.names, a.kinds) if k in ("dipole", "quad"))
        idx = [j for j, nm in enumerate(a.names) if nm == dev]
        a._r.set_nsl(idx, 3)
        b = ExtRing(lattice=2, design_mev=3000.0, nsl=10, nsl_dev={dev: 3})
        st = np.zeros((2, 6)); st[:, 0] = [1e-3, -5e-4]
        self.assertTrue(np.array_equal(a.track(st, 2).view(np.uint64), b.track(st, 2).view(np.uint64)))
        a._r.set_nsl(idx, 10)                                         # вернуть — снова как исходное
        c0 = ExtRing(lattice=2, design_mev=3000.0, nsl=10)
        self.assertTrue(np.array_equal(a.track(st, 2).view(np.uint64), c0.track(st, 2).view(np.uint64)))

    @unittest.skipUnless(os.path.exists(os.path.join(mt.TB, "lattice_configs.csv")), "нет БД")
    def test_cockpit_default_slices_config(self):
        """Кокпит по умолчанию берёт набор срезов из workspace/configs/slices_<lat>.json (tight для реш. 26);
        нет конфига — как раньше (None); другой отпечаток машины — предупреждение в пояснении."""
        import json, tempfile, shutil
        from synchrotronsim import convergence as cv
        if os.path.exists(cv.config_path(26)):
            d, note = cv.cockpit_slices(26)
            self.assertGreater(len(d), 50); self.assertIn("tight", note)
            self.assertEqual(d.get("und_UPM"), 1)                        # полюса ИД — 1 срез
        self.assertEqual(cv.cockpit_slices(9999)[0], None)
        tmp = tempfile.mkdtemp(prefix="ss_slices_"); self.addCleanup(shutil.rmtree, tmp, True)
        fake = os.path.join(tmp, "slices_3.json")
        with open(fake, "w", encoding="utf-8") as f:
            json.dump({"preset": "tight", "tol_nu": 1e-4, "fingerprint": "0" * 64, "nsl_dev": {"X": 3}}, f)
        orig = cv.config_path
        cv.config_path = lambda lattice, root=None: fake
        try:
            d, note = cv.cockpit_slices(3)
        finally:
            cv.config_path = orig
        self.assertEqual(d, {"X": 3}); self.assertIn("решётка изменилась", note)


if __name__ == "__main__":
    unittest.main(verbosity=2, warnings="ignore")
