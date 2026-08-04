# -*- coding: utf-8 -*-
"""Регрессия C++ матричного ядра против Python-ЭТАЛОНА (matrix_tracker.py).

Для набора эталонных элементов СИЛЫ (диполь/DQ/квад/секступоль/октуполь) и входных состояний прогоняет
тело элемента ДВУМЯ ядрами:
  • Python: matrix_tracker.Magnet.body (секторный интегратор + классическое СИ) — эталон;
  • C++:    build/MatrixTracker.exe --element … (movok matrix_core.h).
и требует совпадения с точностью до ~1e-11 (перенос дословный; расхождение — только последние биты libm/pow).

Запуск: py -3.14 test_matrix_core_cpp.py   (Python 3.14, pytest НЕ нужен).
Собрать ядро: cmake --build build --target MatrixTracker
"""
import os, sys, subprocess, unittest

os.environ["SS_SECTOR"] = "1"                    # секторный режим — ДО импорта matrix_tracker (как в проде)
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)
os.chdir(ROOT)
import matrix_tracker as mt

EXE = os.path.join(ROOT, "build", "MatrixTracker.exe")
QP = mt.QP
NSL = 5
# эталонные элементы: (kind, Leff, angle, k1, k2, k3, edge1, edge2)
ELEMENTS = {
    "dq1b":  ("dq", 1.0, 0.0247680586728496, -1.960219550894432, 0.0, 0.0, -0.0116800000099, -0.0116800000099),
    "dl1a5": ("dipole", 0.4178302891796566, 0.0094691022946895, 0.0, 0.0, 0.0, -0.01268495408, 0.00307504592),
    "qf1a":  ("quad", 0.295, 0.0, 2.462177415336321, 0.0, 0.0, 0.0, 0.0),
    "qd2a":  ("quad", 0.162, 0.0, -2.970346620956298, 0.0, 0.0, 0.0, 0.0),
    "sd1a":  ("sext", 0.166, 0.0, 0.0, -136.1950802309929, 0.0, 0.0, 0.0),
    "oct1":  ("oct", 0.09, 0.0, 0.0, 0.0, -13200.0, 0.0, 0.0),
}
STATES = [
    (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (1e-3, 2e-4, 0.5e-3, -1e-4, 0.0, 0.0),
    (2e-3, 0.0, -1e-3, 3e-4, 0.0, 1e-3),
    (-1.5e-3, -5e-4, 2e-3, 1e-4, 0.0, -2e-3),
]


def make_py(kind, Leff, angle, k1, k2, k3, e1, e2, sr, design):
    """Python-магнит через фикстуру (mat_* как в make_magnet теста), при sr — с включённым СИ."""
    old_sr, old_de = os.environ.get("SS_SR"), os.environ.get("SS_DESIGN_MEV")
    if sr:
        os.environ["SS_SR"] = "1"; os.environ["SS_DESIGN_MEV"] = repr(design)
    else:
        os.environ.pop("SS_SR", None)
    try:
        matB = -angle / (QP * Leff) if angle else 0.0
        insts = [dict(x="0.0", y="0.0", z="0.0", yaw="0.0", pitch="0.0", roll="0.0",
                      device_type_id="TEST", instance_id="1", eff_length=repr(float(Leff)),
                      kind=kind, n_slices=str(NSL), edge1=repr(float(e1)), edge2=repr(float(e2)),
                      mat_B=repr(matB), mat_G=repr(k1 / QP), mat_Gv=repr(k1 / QP),
                      mat_S=repr(k2 / QP), mat_O=repr(k3 / QP))]
        return mt.Magnet(insts, {}, {}, {}, {}, settings=mt.Settings()), (matB, k1 / QP, k1 / QP, k2 / QP, k3 / QP)
    finally:
        os.environ.pop("SS_SR", None) if old_sr is None else os.environ.__setitem__("SS_SR", old_sr)
        if old_de is not None: os.environ["SS_DESIGN_MEV"] = old_de


def run_cpp(kind, Leff, mats, e1, e2, design, sr, st):
    matB, matG, matGv, matS, matO = mats
    args = [EXE, "--element", kind, repr(Leff), repr(matB), repr(matG), repr(matGv), repr(matS), repr(matO),
            repr(e1), repr(e2), str(NSL), repr(design), "1" if sr else "0"] + [repr(v) for v in st]
    out = subprocess.run(args, capture_output=True, text=True).stdout.strip()
    if out == "LOST":
        return "LOST"
    return tuple(float(v) for v in out.split())


class TestMatrixCoreCpp(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        if not os.path.exists(EXE):
            raise unittest.SkipTest(f"C++ ядро не собрано: {EXE} (cmake --build build --target MatrixTracker)")

    def _check(self, name, sr, design, tol=1e-10):
        kind, Leff, angle, k1, k2, k3, e1, e2 = ELEMENTS[name]
        m, mats = make_py(kind, Leff, angle, k1, k2, k3, e1, e2, sr, design)
        worst = 0.0
        for st in STATES:
            try:
                py = m.body(tuple(map(float, st)))
            except mt.ParticleLost:
                py = "LOST"
            cpp = run_cpp(kind, Leff, mats, e1, e2, design, sr, st)
            self.assertEqual(py == "LOST", cpp == "LOST", f"{name} {st}: расхождение потери")
            if py == "LOST":
                continue
            for a, b, lab in zip(py, cpp, ("x", "px", "y", "py", "t", "delta")):
                d = abs(a - b)
                worst = max(worst, d)
                self.assertLess(d, tol, f"{name} {st} [{lab}]: py={a:.17g} cpp={b:.17g} Δ={d:.2e}")
        return worst

    def test_bodies_no_radiation(self):
        w = max(self._check(n, sr=False, design=6000.0) for n in ELEMENTS)
        print(f"\n  [без СИ] макс. расхождение C++ vs Python по 6 элементам × 4 состояния: {w:.2e}")

    def test_dipoles_with_radiation(self):
        w = max(self._check(n, sr=True, design=6000.0) for n in ("dq1b", "dl1a5"))
        print(f"  [с СИ 6ГэВ] макс. расхождение (диполи DQ/DL): {w:.2e}")


class TestKickServe(unittest.TestCase):
    """C++ `--serve` (кик-траектория дрейф→элемент(M суб-срезов L/M)→дрейф) vs Python-ЭТАЛОН — то, что
    element_compare.track_matrix теперь считает C++ ядром. Бит-в-бит (перенос дословный)."""
    M, NSUB, DI, DO, DELTA = 24, 8, 1.35, 1.35, -4.212140059489028e-05
    PROBES = [(0.0, 0.0), (1e-3, 0.0), (-1e-3, 5e-4), (0.0, -1e-3)]

    @classmethod
    def setUpClass(cls):
        import numpy as np
        cls.np = np
        if not os.path.exists(EXE):
            raise unittest.SkipTest(f"C++ ядро не собрано: {EXE}")
        cls.srv = subprocess.Popen([EXE, "--serve"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                   text=True, bufsize=1)

    @classmethod
    def tearDownClass(cls):
        try:
            cls.srv.stdin.write("QUIT\n"); cls.srv.stdin.flush(); cls.srv.wait(timeout=5)
        except Exception:
            cls.srv.kill()

    def _py(self, kind, L, mats, e1, e2, probe):                    # эталон: дрейф(8)→M×Magnet.body(L/M)→дрейф(8)
        np = self.np; mB, mG, mGv, mS, mO = mats
        insts = [dict(x="0", y="0", z="0", yaw="0", pitch="0", roll="0", device_type_id="T", instance_id="1",
                      eff_length=repr(L / self.M), kind=kind, n_slices=str(NSL), edge1=repr(e1), edge2=repr(e2),
                      mat_B=repr(mB), mat_G=repr(mG), mat_Gv=repr(mGv), mat_S=repr(mS), mat_O=repr(mO))]
        sub = mt.Magnet(insts, {}, {}, {}, {}, settings=mt.Settings())
        u0, v0 = probe; st = (u0, 0.0, v0, 0.0, 0.0, self.DELTA); pts = [(0.0, u0, v0)]; s = 0.0
        for _ in range(self.NSUB): st = mt._drift(st, self.DI / self.NSUB); s += self.DI / self.NSUB; pts.append((s, st[0], st[2]))
        for _ in range(self.M):    st = sub.body(st);                       s += L / self.M;           pts.append((s, st[0], st[2]))
        for _ in range(self.NSUB): st = mt._drift(st, self.DO / self.NSUB); s += self.DO / self.NSUB; pts.append((s, st[0], st[2]))
        return np.array(pts)

    def _cpp(self, kind, L, mats, e1, e2):
        np = self.np; mB, mG, mGv, mS, mO = mats
        req = (f"{kind} {L!r} {self.M} {self.DI!r} {self.DO!r} {NSL} 6000.0 0 {self.DELTA!r} "
               f"{mB!r} {mG!r} {mGv!r} {mS!r} {mO!r} {e1!r} {e2!r} {len(self.PROBES)}")
        for u, v in self.PROBES: req += f" {u!r} {v!r}"
        self.srv.stdin.write(req + "\n"); self.srv.stdin.flush()
        out = []
        for _ in self.PROBES:
            ln = self.srv.stdout.readline().split()
            out.append(None if ln[0] == "LOST" else np.asarray(ln[1:], float).reshape(int(ln[0]), 3))
        return out

    def test_serve_vs_python(self):
        worst = 0.0
        for name, (kind, Leff, angle, k1, k2, k3, e1, e2) in ELEMENTS.items():
            mats = (-angle / (QP * Leff) if angle else 0.0, k1 / QP, k1 / QP, k2 / QP, k3 / QP)
            cpp = self._cpp(kind, Leff, mats, e1, e2)
            for i, p in enumerate(self.PROBES):
                py = self._py(kind, Leff, mats, e1, e2, p)
                self.assertIsNotNone(cpp[i], f"{name} {p}: C++ LOST, Python нет")
                worst = max(worst, float(self.np.abs(py - cpp[i]).max()))
        self.assertLess(worst, 1e-10, f"--serve vs Python: Δ={worst:.2e}")
        print(f"\n  [--serve] макс. расхождение кик-траектории C++ vs Python (6 элем × 4 пробы): {worst:.2e}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
