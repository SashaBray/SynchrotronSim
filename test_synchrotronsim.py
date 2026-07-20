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
    EXE = next((p for p in (os.path.join(ROOT, "build", "SynchrotronTracker.exe"),
                            os.path.join(ROOT, "build", "Release", "SynchrotronTracker.exe"))
                if os.path.exists(p)), None)

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


if __name__ == "__main__":
    unittest.main(verbosity=2, warnings="ignore")
