# -*- coding: utf-8 -*-
"""Регрессия АВТОНОМНОГО кольца C++ против Python-ЭТАЛОНА (полный трек из БД).

C++ (build/MatrixTracker.exe --track <exp> <turns>) читает БД, инжектирует пучок, трекает N оборотов с СИ+ВЧ
(по radiation_params.csv) и печатает финальные состояния. Python-эталон — matrix_traj_export.long_trajectories
(тот же цикл: modified_drift + body + ВЧ). Сравниваются (x,px,y,py) у плоскости старта и выживание.

Экзерсизит весь геометрический слой: build_line, рамки frame/_rot, секторные дуги, modified_drift, инжекцию.
Запуск: py -3.14 test_matrix_core_ring.py    Сборка: cmake --build build --target MatrixTracker
"""
import os, sys, subprocess, unittest, shutil

os.environ["SS_SECTOR"] = "1"                    # секторный режим — ДО импорта matrix_tracker
os.environ["SS_SR"] = "1"                        # СИ (exp1/2 в radiation_params.csv включают его)
os.environ.setdefault("SS_DESIGN_MEV", "6000")   # СИЛА 6 ГэВ (как radiation_params.csv exp1)
ROOT = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, ROOT); os.chdir(ROOT)
import numpy as np
import matrix_traj_export as mte
from _common import exe_name

EXE = os.path.join(ROOT, "build", exe_name("MatrixTracker"))   # .exe только на Windows
NPZ1 = os.path.join(ROOT, "workspace", "results", "trajectories_local_exp1.npz")   # --track перезапишет — бэкапим


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
        cls._bak = NPZ1 + ".testbak" if os.path.exists(NPZ1) else None   # --track перезапишет exp1 — сохраним
        if cls._bak:
            shutil.copy2(NPZ1, cls._bak)

    @classmethod
    def tearDownClass(cls):
        if cls._bak and os.path.exists(cls._bak):                        # вернуть пользовательский exp1 npz
            shutil.move(cls._bak, NPZ1)

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

    def test_npz_output(self):
        # C++ пишет trajectories_local_exp1.npz — читаемый numpy/read_local_npz, tt совпадает с эталоном.
        turns = 5
        subprocess.run([EXE, "--track", "1", str(turns)], capture_output=True, text=True)
        self.assertTrue(os.path.exists(NPZ1), "C++ не создал .npz")
        np.load(NPZ1, allow_pickle=False).close()                       # формат npz валиден для numpy
        r = mte.read_local_npz(NPZ1)
        self.assertEqual(r["mode"], "long")
        self.assertEqual(r["n_turns"], turns)
        self.assertEqual(r["sr"], 1); self.assertEqual(r["rf"], "ideal")
        self.assertAlmostEqual(r["design_mev"], 6000.0)
        self.assertGreater(r["ring_len"], 1000.0)                       # СИЛА ~1110 м
        d = mte.long_trajectories(1, 1, n_turns=turns, rf="ideal")      # эталон
        worst = 0.0
        for cp, py in zip(r["parts"], d["turnbyturn"]):
            for comp in ("tt_x", "tt_xp", "tt_y", "tt_yp"):
                worst = max(worst, float(np.max(np.abs(np.asarray(cp[comp]) - py[:, ("tt_x", "tt_xp", "tt_y", "tt_yp").index(comp)]))))
        self.assertLess(worst, 1e-9, f".npz tt расходится с эталоном: {worst:.2e}")
        print(f"  [СИЛА exp1, .npz] numpy читает OK; tt vs эталон: {worst:.2e}; ring_len={r['ring_len']:.1f}м")

    def test_full_mode(self):
        # ПОЛНЫЙ режим (x(S)/спектр/тюны): C++ --full пишет mode=full npz — сверяем с local_trajectories.
        subprocess.run([EXE, "--full", "1", "1"], capture_output=True, text=True)
        r = mte.read_local_npz(NPZ1)
        self.assertEqual(r["mode"], "full")
        d = mte.local_trajectories(1, 1, n_turns=1)                    # Python эталон (SS_SR=1 в env)
        worst = 0.0
        for cp, loc, glo in zip(r["parts"], d["locals"], d["globals"]):
            self.assertEqual(len(cp["S"]), len(loc), "разное число точек (прореживание?)")
            for k, col in enumerate(("S", "x", "xp", "y", "yp", "ct", "delta")):
                worst = max(worst, float(np.max(np.abs(np.asarray(cp[col]) - loc[:, k]))))
            self.assertTrue(np.array_equal(np.asarray(cp["elem"]), loc[:, 7].astype(int)), "elem-индексы разошлись")
            for k, col in enumerate(("X", "Y", "Z")):
                worst = max(worst, float(np.max(np.abs(np.asarray(cp[col]) - glo[:, k]))))
        self.assertLess(worst, 1e-9, f"full-режим расходится с эталоном: {worst:.2e}")
        print(f"  [СИЛА exp1, full] точек/частицу {len(r['parts'][0]['S'])}; локальн.+глоб. vs эталон: {worst:.2e}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
