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

    def test_sine_rf_sila26(self):
        # РЕАЛЬНЫЙ резонатор (radiation_params: exp 26, rf=sine, 8 МВ, h=1296, фаза auto): C++ --track против
        # Python-эталона по ВСЕМ шести координатам (t — запаздывание от часов ВЧ, δ — синхротронное движение).
        import matrix_tracker as mt
        npz26 = os.path.join(ROOT, "workspace", "results", "trajectories_local_exp26.npz")
        bak = npz26 + ".testbak" if os.path.exists(npz26) else None
        if bak:
            shutil.copy2(npz26, bak)
        try:
            spec = mt.rf_spec(mt.load_radiation(26))
            self.assertIsInstance(spec, dict); self.assertEqual(spec["type"], "sine")
            turns = 3
            d = mte.long_trajectories(26, 1, n_turns=turns, rf=spec)
            cpp = cpp_track(26, turns)
            self.assertAlmostEqual(d["rf_phase_rad"], 2.8102, delta=2e-3)     # φs ≈ 161,0° (U0 ≈ 2,60 МэВ: диполи + СИ ондуляторов)
            worst = 0.0
            for i, (fin, surv) in enumerate(zip(d["finals"], d["surv_turns"])):
                self.assertEqual(surv, turns); self.assertNotEqual(cpp[i][0], "LOST")
                worst = max(worst, max(abs(a - b) for a, b in zip(fin, cpp[i][:6])))
            self.assertLess(worst, 1e-9, f"sine-ВЧ: C++ vs Python {worst:.2e}")
            self.assertEqual(mte.read_local_npz(npz26)["rf"], "sine")
            print(f"  [СИЛА-26 exp26, ВЧ sine, {turns} об.] 6D C++ vs Python: {worst:.2e}; "
                  f"U0={d['rf_u0_mev']:.4f} МэВ")
        finally:
            if bak and os.path.exists(bak):
                shutil.move(bak, npz26)

    def _npz_guard(self, exp):
        npz = os.path.join(ROOT, "workspace", "results", f"trajectories_local_exp{exp}.npz")
        bak = npz + ".testbak2"
        if os.path.exists(npz):
            shutil.copy2(npz, bak)
        return npz, bak

    def _npz_restore(self, npz, bak):
        if os.path.exists(bak):
            shutil.move(bak, npz)

    def test_portal_open_line_sp1(self):
        # РАЗОМКНУТАЯ линия (exp 3, SP1): 4 пролёта старт → финиш через портал, C++ ≡ Python (СИ у опыта 3
        # выключена — строки radiation_params нет; Python-эталон собирается без СИ).
        npz, bak = self._npz_guard(3)
        try:
            os.environ["SS_SR"] = ""
            try:
                d = mte.long_trajectories(3, 1, n_turns=4)
            finally:
                os.environ["SS_SR"] = "1"
            cpp = cpp_track(3, 4)
            self.assertEqual(d["n_turns"], 4); self.assertFalse(d["closed"])
            worst = 0.0
            for i, (fin, surv) in enumerate(zip(d["finals"], d["surv_turns"])):
                self.assertEqual(surv, 4); self.assertNotEqual(cpp[i][0], "LOST"); self.assertEqual(cpp[i][6], 4)
                worst = max(worst, max(abs(a - b) for a, b in zip(fin, cpp[i][:6])))
            self.assertLess(worst, 1e-9, f"портал SP1: C++ vs Python {worst:.2e}")
            print(f"  [SP1 exp3, портал, 4 пролёта] 6D C++ vs Python: {worst:.2e}")
        finally:
            self._npz_restore(npz, bak)

    def test_ideal_target_energy(self):
        # Идеальный резонатор с заданной энергией (radiation_params.rf_ideal_energy_mev, временно у опыта 1):
        # C++ ≡ Python по 6D, δ после резонатора (в конце пролёта) = E/E_design − 1 у каждой частицы.
        import csv
        import matrix_tracker as mt
        path = os.path.join(mt.TB, "radiation_params.csv"); bak = path + ".testbak"
        npz, nbak = self._npz_guard(1)
        shutil.copy2(path, bak)
        try:
            with open(path, encoding="utf-8") as f:
                rows = list(csv.DictReader(f, delimiter=";"))
            hdr = list(rows[0].keys())
            for r in rows:
                if r["exp_id"] == "1":
                    r["rf_ideal_energy_mev"] = "6006"
            with open(path, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=hdr, delimiter=";"); w.writeheader(); w.writerows(rows)
            spec = mt.rf_spec(mt.load_radiation(1))
            self.assertAlmostEqual(mt.ideal_delta(spec), 1e-3, places=15)
            d = mte.long_trajectories(1, 1, n_turns=3, rf=spec)
            cpp = cpp_track(1, 3)
            worst = 0.0
            for i, fin in enumerate(d["finals"]):
                self.assertEqual(fin[5], 6006.0 / 6000.0 - 1.0); self.assertEqual(cpp[i][5], fin[5])
                worst = max(worst, max(abs(a - b) for a, b in zip(fin, cpp[i][:6])))
            self.assertLess(worst, 1e-9, f"ideal с заданной энергией: C++ vs Python {worst:.2e}")
            print(f"  [exp 1, ideal → 6006 МэВ] δ после резонатора {d['finals'][0][5]:.6e}; C++ vs Python {worst:.2e}")
        finally:
            shutil.move(bak, path)
            self._npz_restore(npz, nbak)

    def test_cavity_mid_ring(self):
        # РЕЗОНАТОР УСТАНОВКИ в середине кольца (rf_cavities.csv, временная строка): sine (exp 26) и ideal (exp 1)
        # — C++ ≡ Python по 6D; резонатор в плоскости старта ≡ решётка без резонаторов (бит-в-бит в обоих ядрах).
        import numpy as np
        import matrix_tracker as mt
        import exp_planes as ep
        path = os.path.join(mt.TB, "rf_cavities.csv"); bak = path + ".testbak"
        guards = [self._npz_guard(e) for e in (1, 26)]
        shutil.copy2(path, bak)

        def mid(lat, s0):
            mags = ep._mags(lat); _ring, _lead, s_in, s_out, _C = ep._arcs(mags)
            j = next(k for k in range(len(mags) - 1) if s_out[k] >= s0 and s_in[k + 1] - s_out[k] > 0.2)
            return ep.plane_at_s(lat, 0.5 * (s_out[j] + s_in[j + 1]), mags)
        try:
            ep.write_cavity(26, mid(26, 555.0), cavity_id=901, note="тест")
            ep.write_cavity(1, mid(1, 300.0), cavity_id=902, note="тест")
            spec = mt.rf_spec(mt.load_radiation(26))
            for exp, lat, rf, turns in ((26, 26, spec, 3), (1, 1, "ideal", 3)):
                d = mte.long_trajectories(lat, 1, n_turns=turns, rf=rf)
                self.assertEqual(len(d["cavities"]), 1)
                cpp = cpp_track(exp, turns)
                worst = max(max(abs(a - b) for a, b in zip(fin, cpp[i][:6])) for i, fin in enumerate(d["finals"]))
                self.assertLess(worst, 1e-9, f"резонатор в середине, exp {exp}: C++ vs Python {worst:.2e}")
                print(f"  [exp {exp}, резонатор s={d['cavities'][0][1]:.3f} м, {d['rf']}] 6D C++ vs Python: {worst:.2e}")
            shutil.copy2(bak, path)
            ref_py = mte.long_trajectories(26, 1, n_turns=2, rf=spec); ref_c = cpp_track(26, 2)
            ep.write_cavity(26, (np.zeros(3), np.array([1.0, 0.0, 0.0])), cavity_id=903, note="тест: шов")
            seam_py = mte.long_trajectories(26, 1, n_turns=2, rf=spec); seam_c = cpp_track(26, 2)
            self.assertTrue(np.array_equal(np.asarray(seam_py["finals"]), np.asarray(ref_py["finals"])))
            self.assertEqual(seam_c, ref_c)
        finally:
            shutil.move(bak, path)
            for g in guards:
                self._npz_restore(*g)


if __name__ == "__main__":
    unittest.main(verbosity=2)
