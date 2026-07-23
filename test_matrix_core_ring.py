# -*- coding: utf-8 -*-
"""Регрессия АВТОНОМНОГО кольца C++ против Python-ЭТАЛОНА (полный трек из БД).

C++ (build/MatrixTracker.exe --track <exp> <turns>) читает БД, инжектирует пучок, трекает N оборотов с СИ+ВЧ
(по radiation_params.csv) и печатает финальные состояния. Python-эталон — matrix_traj_export.long_trajectories
(тот же цикл: modified_drift + body + ВЧ). Сравниваются (x,px,y,py) у плоскости старта и выживание.

Экзерсизит весь геометрический слой: build_line, рамки frame/_rot, секторные дуги, modified_drift, инжекцию.
Запуск: py -3.14 test_matrix_core_ring.py    Сборка: cmake --build build --target MatrixTracker
"""
import os, sys, subprocess, unittest

os.environ["SS_SECTOR"] = "1"                    # секторный режим — ДО импорта matrix_tracker
os.environ["SS_SR"] = "1"                        # СИ (exp1/2 в radiation_params.csv включают его)
os.environ.setdefault("SS_DESIGN_MEV", "6000")   # СИЛА 6 ГэВ (как radiation_params.csv exp1)
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT); os.chdir(ROOT)
import numpy as np
import matrix_traj_export as mte

EXE = os.path.join(ROOT, "build", "MatrixTracker.exe")


def cpp_track(exp, turns):
    """C++ --track → {idx: (x,px,y,py,t,delta,surv) | ('LOST',surv)}."""
    out = subprocess.run([EXE, "--track", str(exp), str(turns)], capture_output=True, text=True).stdout
    res = {}
    for ln in out.splitlines():
        if ln.startswith("#") or not ln.strip():
            continue
        p = ln.split()
        if p[1] == "LOST":
            res[int(p[0])] = ("LOST", int(p[2]))
        else:
            res[int(p[0])] = tuple(float(v) for v in p[1:7]) + (int(p[7]),)
    return res


class TestMatrixCoreRing(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(EXE):
            raise unittest.SkipTest(f"C++ ядро не собрано: {EXE} (cmake --build build --target MatrixTracker)")

    def _compare(self, exp, lattice, group, turns, tol=1e-9):
        # Python-эталон: тот же цикл (инжекция + N оборотов + ВЧ ideal)
        d = mte.long_trajectories(lattice, group, n_turns=turns, rf="ideal")
        cpp = cpp_track(exp, turns)
        worst = 0.0
        for i, (tt, surv) in enumerate(zip(d["turnbyturn"], d["surv_turns"])):
            self.assertIn(i, cpp, f"C++ не выдал частицу {i}")
            c = cpp[i]
            lost_py = surv < turns
            lost_cpp = (c[0] == "LOST")
            self.assertEqual(lost_py, lost_cpp, f"exp{exp} ч{i}: расхождение выживания (py surv={surv}, cpp={c})")
            if lost_py:
                continue
            py = tt[-1]                          # (x, x'=px, y, y'=py) после последнего оборота
            for a, b, lab in zip(py, c[:4], ("x", "px", "y", "py")):
                dd = abs(float(a) - float(b)); worst = max(worst, dd)
                self.assertLess(dd, tol, f"exp{exp} ч{i} [{lab}]: py={a:.17g} cpp={b:.17g} Δ={dd:.2e}")
        return worst

    def test_sila_one_turn(self):
        w = self._compare(1, 1, 1, turns=1)
        print(f"\n  [СИЛА exp1, 1 оборот] макс. расхождение C++ vs Python: {w:.2e}")

    def test_sila_five_turns(self):
        w = self._compare(1, 1, 1, turns=5)
        print(f"  [СИЛА exp1, 5 оборотов] макс. расхождение C++ vs Python: {w:.2e}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
