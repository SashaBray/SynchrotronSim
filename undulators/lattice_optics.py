# -*- coding: utf-8 -*-
"""Реальные оптические функции (β, η) в месте элемента из twiss решётки — вместо допущений.

Twiss считается по C++-ядру (ring_cpp + ring_analysis), индекс в месте устройства ищется ПО ГЛОБАЛЬНОМУ
ЦЕНТРУ (порядок трекинга ≠ instance_id, поэтому сопоставляем по mag.center). Кэш по решётке.

Исп.: optics = undulator_optics(27, "und_19") → betax, betay, etax, etaxp, s.
"""
import os, sys, functools
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")
import numpy as np
import ring_cpp as rc
import ring_analysis as ra
import errors as er


@functools.lru_cache(maxsize=8)
def _twiss(lattice):
    r = rc.CppRing(lattice=lattice)
    tw = ra.twiss(r)
    cen = np.array([m.center for m in r.mags])                 # глобальные центры в порядке трекинга
    sout = np.asarray(r.s_out)
    out = dict(cen=cen, sout=sout, s=np.asarray(tw["s"]),
               bx=np.asarray(tw["betax"]), by=np.asarray(tw["betay"]),
               ex=np.asarray(tw["etax"]), exp=np.asarray(tw["etaxp"]),
               ax=np.asarray(tw["alphax"]), ay=np.asarray(tw["alphay"]),
               nux=float(tw["nux"]), nuy=float(tw["nuy"]))
    r.close()
    return out


def optics_at_point(lattice, pivot):
    """β, η в точке глобального пространства pivot (ближайший магнит → twiss по s)."""
    tw = _twiss(lattice); piv = np.asarray(pivot, float)
    i = int(np.argmin(np.linalg.norm(tw["cen"] - piv, axis=1)))
    j = int(np.argmin(np.abs(tw["s"] - tw["sout"][i])))
    return dict(betax=float(tw["bx"][j]), betay=float(tw["by"][j]),
                etax=float(tw["ex"][j]), etaxp=float(tw["exp"][j]),
                alphax=float(tw["ax"][j]), alphay=float(tw["ay"][j]),
                s=float(tw["sout"][i]), match_mm=float(np.linalg.norm(tw["cen"][i] - piv)) * 1e3,
                nux=tw["nux"], nuy=tw["nuy"])


def undulator_optics(lattice, group_id):
    g = next((x for x in er.load_groups(lattice) if x["group_id"] == group_id), None)
    if g is None:
        raise SystemExit(f"группа {group_id} не найдена в element_groups решётки {lattice}")
    return optics_at_point(lattice, g["pivot"])


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="реальные β,η в месте ондулятора из twiss")
    ap.add_argument("--lattice", type=int, default=27)
    ap.add_argument("--group", default="und_19")
    a = ap.parse_args()
    o = undulator_optics(a.lattice, a.group)
    print(f"{a.group} (реш.{a.lattice}): β=({o['betax']:.3f},{o['betay']:.3f}) м, "
          f"αx={o['alphax']:+.3f} αy={o['alphay']:+.3f}, ηx={o['etax']*1e3:.3f} мм, "
          f"η'x={o['etaxp']*1e3:.3f} мрад, s={o['s']:.3f} м (совпад. {o['match_mm']:.2f} мм)")
    print(f"νx={o['nux']:.5f} νy={o['nuy']:.5f}")
