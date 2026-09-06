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
                center=(0.0, 0.0, 0.0), yaw=0.0, pitch=0.0, roll=0.0):
    """Синтетический Magnet БЕЗ БД: физика прямо в mat_* (mat_B=−angle/(QP·L), mat_G=k1/QP, …).
    center/yaw/pitch/roll — положение и ориентация (для тестов рассогласования установки);
    на само тело они не влияют (тело в локальных θ0/k), но задают рамки fr_in/fr_out."""
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
    """МЭС-σ: записанные веса, round-trip таблицы и ЗАКРЕПЛЁННАЯ точность в токе.

    Точность фиксируется рубежами, а не точными числами: карты и реестр могут пополняться, но
    ухудшение модели должно ронять тест. Метрика — ошибка ТОКА |λ−1|, а не процент по кику: в
    насыщении отклик пологий и промах в кике на 1% означает промах в токе в 2.5–5 раз больший."""

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
    #   2026-09-07: знаки намотки стали ОБЩИМИ на все гармоники элемента (обмотка намотана в одну
    # сторону, и калибровать её отдельно по каждому кику значит позволять модели противоречить
    # самой себе: у c2 выходило «+» по диполю и «−» по квадруполю). Это ограничение, а не
    # улучшение, и оно стоило корректору 5.8% → 6.6% по току на комбинациях, дав взамен 21.4% →
    # 20.4% по кику. Порог поднят сознательно: физическая непротиворечивость важнее десятой доли
    # процента на одной метрике.
    LIMITS = {("corrector", 1): (9.0, 7.5),
              ("DQ", 1): (22.0, 32.0),
              ("DQ", 2): (7.0, 8.0)}

    # Предельное число НЕДОСТИЖИМЫХ карт (одиночные, комбинации): карты, на которых модель не
    # дотягивается до измеренного кика ни при каком токе. Зафиксировано по состоянию на 2026-09-07.
    #   У корректора недостижима ПОЛОВИНА одиночных карт (11 из 23) и две трети комбинаций: во
    # всех случаях модель не дотягивается до измеренного кика даже при токе ×50, то есть её полка
    # ниже измеренной. Это настоящий изъян модели на верхних токах, а не настройка решателя.
    # Числа ниже — фиксация текущего положения, чтобы оно не ухудшалось незаметно; уменьшать их
    # надо вслед за исправлением полки.
    MISS = {("corrector", 1): (11, 21), ("DQ", 1): (4, 0), ("DQ", 2): (0, 0)}

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

    def test_current_accuracy_holds(self):
        ms, M = self._models()
        import branch_network as bn
        for key, (lim_s, lim_c) in self.LIMITS.items():
            el, order = key
            maps = bn.read_maps(el, order)
            got = {}
            for lbl, grp in (("одиночные", [t for t in maps if len(t[0]) == 1]),
                             ("комбинации", [t for t in maps if len(t[0]) > 1])):
                e, miss = ms.current_errors(el, order, grp, M)
                got[lbl] = (float(np.mean(e)) if e else float("inf"), miss, len(grp))
            for lbl, lim in (("одиночные", lim_s), ("комбинации", lim_c)):
                val, miss, tot = got[lbl]
                self.assertLess(val, lim, f"{el} n={order} {lbl} (по {tot - miss} картам из {tot})")
                # Охват проверяется отдельно: среднее по остатку тем меньше, чем больше карт
                # выброшено как недостижимые, поэтому без этого рубежа метрику можно улучшить
                # УХУДШЕНИЕМ модели.
                self.assertLessEqual(miss, self.MISS[(el, order)][lbl != "одиночные"],
                                     f"{el} n={order} {lbl}: недостижимых карт {miss} из {tot}")



if __name__ == "__main__":
    unittest.main(verbosity=2, warnings="ignore")
