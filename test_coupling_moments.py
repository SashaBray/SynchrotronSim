# -*- coding: utf-8 -*-
"""Имена моментов матриц связи: скью-квадруполь «Gs» ≠ нормальный вертикальный градиент ядра «mat_Gv».

`element_compare.extract_coeffs` извлекает из кика −∂x′/∂y — СКЬЮ-квадруполь (связь y→x′). Ядра
SynchrotronSim читают `mat_Gv` как НОРМАЛЬНЫЙ вертикальный градиент k1v = QP·mat_Gv (∂By/∂v, y→y′).
До 2026-11-16 скью назывался «Gv» и через `_design_ivec` мог принять mat_Gv решётки за свой проект.
Тесты держат разделение: имя, смысл извлечения, отсутствие столбца mat_* и чтение старых таблиц.

Без БД и без диффура (траектории синтетические). Запуск:
    py -3.14 -m unittest test_coupling_moments
"""
import os, sys, csv, tempfile, unittest, warnings

warnings.filterwarnings("ignore", category=ResourceWarning)
os.environ.setdefault("SS_SECTOR", "1")
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

import numpy as np

HDR = ["element", "kind", "n_windings", "n_moments", "direction", "out_label", "in_label",
       "value", "units", "eff_length", "i_op", "imax", "r2"]


def _synthetic_trajs(ec, kick):
    """Траектории (s, Y, Z) проб ec.PROBES: прямые с наклонами (x′, y′) = kick(x, y)."""
    s = np.linspace(0.0, 1.0, 21)
    out = []
    for x, y in ec.PROBES:
        xp, yp = kick(x, y)
        out.append(np.column_stack([s, x + xp * s, y + yp * s]))
    return out


class TestSkewQuadMoment(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        import element_compare as ec
        cls.ec = ec

    def _extract(self, kick, L=0.1):
        ec = self.ec
        old_run, old_clean = ec.run_diffur_windings, ec.cleanup
        ec.run_diffur_windings = lambda windings, yaw=0.0: (_synthetic_trajs(ec, kick), None)
        ec.cleanup = lambda: None
        try:
            return ec.extract_coeffs([("dev", 1.0)], L)
        finally:
            ec.run_diffur_windings, ec.cleanup = old_run, old_clean

    def test_name_is_gs(self):
        """Скью-квадруполь называется Gs; «Gv» в именах моментов матриц связи больше нет."""
        ec = self.ec
        self.assertIn("Gs", ec.MULTIPOLES)
        self.assertNotIn("Gv", ec.MULTIPOLES)
        self.assertIn("скью", ec.MOMENT_DESC["Gs"])
        self.assertNotIn("Gv", ec.MOMENT_DESC)

    def test_gs_is_skew_not_vertical_gradient(self):
        """Gs ловит связь y→x′ (скью) и НЕ реагирует на нормальный вертикальный градиент y→y′ (k1v ядра)."""
        ec, L = self.ec, 0.1
        QPL = ec.mt.QP * L
        gs = 2.5e-3
        c = self._extract(lambda x, y: (-QPL * gs * y, -QPL * gs * x), L)      # чистый скью-квадруполь
        self.assertAlmostEqual(c["Gs"] / gs, 1.0, places=6)
        self.assertLess(abs(c["G"]), 1e-9 * gs)
        gv = 4.0e-3
        c = self._extract(lambda x, y: (0.0, -QPL * gv * y), L)               # только ∂y′/∂y (k1v ядра)
        self.assertLess(abs(c["Gs"]), 1e-9 * gv, "нормальный верт. градиент не должен попадать в Gs")
        self.assertLess(abs(c["G"]), 1e-9 * gv)

    def test_gs_has_no_lattice_column(self):
        """У Gs нет столбца mat_*; в частности mat_Gv (нормальный k1v ядра) за скью не принимается."""
        ec = self.ec
        self.assertNotIn("Gs", ec.MOMENT_MATCOL)
        self.assertNotIn("mat_Gv", ec.MOMENT_MATCOL.values())

    def test_design_ivec_ignores_mat_gv(self):
        """Обратная задача: ненулевой mat_Gv решётки не становится проектом скью-момента Gs."""
        ec = self.ec
        old = ec.LCROW
        try:
            m = dict(labels=["I_c1", "I_c2"], devs=["C", "C"], outputs=["B", "Gs"],
                     M=np.array([[2.0, 0.0], [0.0, 3.0]]))
            ec.LCROW = {"C": {"mat_B": "0", "mat_G": "0", "mat_Gv": "5e-3"}}
            Ivec, ctrl, design = ec._design_ivec(m)
            self.assertEqual(design["Gs"], 0.0)
            self.assertEqual(ctrl, [])
            np.testing.assert_array_equal(Ivec, [0.0, 0.0])
            ec.LCROW = {"C": {"mat_B": "1e-3", "mat_G": "0", "mat_Gv": "5e-3"}}
            Ivec, ctrl, design = ec._design_ivec(m)
            self.assertEqual(ctrl, ["B"])                                   # управляется только B
            self.assertEqual(design["Gs"], 0.0)
            self.assertAlmostEqual(Ivec[0], 1e-3 / 2.0)
        finally:
            ec.LCROW = old


class TestLegacyGvTable(unittest.TestCase):
    """Таблица связи, снятая до переименования (скью под именем «Gv»), читается и мигрирует в «Gs»."""

    ROWS = [["K", "kicker", "2", "2", "forward", "B", "I_c1", "1.0", "coef/A", "0.1", "1", "2", "0.99"],
            ["K", "kicker", "2", "2", "forward", "B", "I_c2", "0.0", "coef/A", "0.1", "1", "2", "0.99"],
            ["K", "kicker", "2", "2", "forward", "Gv", "I_c1", "0.0", "coef/A", "0.1", "1", "2", "0.99"],
            ["K", "kicker", "2", "2", "forward", "Gv", "I_c2", "4.0", "coef/A", "0.1", "1", "2", "0.99"],
            ["K", "kicker", "2", "2", "inverse", "I_c2", "Gv", "0.25", "A/coef", "0.1", "1", "2", "0.99"]]

    def _write(self, td):
        import design_currents as dc
        with open(os.path.join(td, dc.COUPLING_TABLE), "w", encoding="utf-8", newline="") as f:
            w = csv.writer(f, delimiter=";")
            w.writerow(HDR); w.writerows(self.ROWS)

    def test_load_maps_gv_to_gs(self):
        import design_currents as dc
        self.assertEqual(dc.moment_name("Gv"), "Gs")
        self.assertEqual(dc.moment_name("G"), "G")
        with tempfile.TemporaryDirectory() as td:
            self._write(td)
            m = {x["name"]: x for x in dc.load_coupling(td)}["K"]
        self.assertEqual(m["outputs"], ["B", "Gs"])
        np.testing.assert_allclose(m["M"], [[1.0, 0.0], [0.0, 4.0]])

    def test_upsert_migrates_kept_rows(self):
        """UPSERT другого элемента переименовывает «Gv» в сохраняемых строках, значения не трогает."""
        import element_compare as ec
        old = ec.TB
        with tempfile.TemporaryDirectory() as td:
            self._write(td)
            ec.TB = td
            try:
                other = ["Q", "quad", "1", "1", "forward", "G", "I", "2.0", "coef/A", "0.3", "1", "9", "0.99"]
                ec.upsert_coupling_table({"Q": [other]})
                with open(os.path.join(td, ec.COUPLING_TABLE), encoding="utf-8") as f:
                    rr = [r for r in list(csv.reader(f, delimiter=";"))[1:] if r]
            finally:
                ec.TB = old
        kept = [r for r in rr if r[0] == "K"]
        self.assertEqual(len(kept), len(self.ROWS))
        self.assertNotIn("Gv", {c for r in kept for c in r[5:7]})
        self.assertEqual(sum(1 for r in kept if "Gs" in r[5:7]), 3)
        self.assertEqual([r[7] for r in kept], [r[7] for r in self.ROWS])  # значения те же


if __name__ == "__main__":
    unittest.main(verbosity=2)
