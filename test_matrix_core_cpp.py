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


if __name__ == "__main__":
    unittest.main(verbosity=2)
