# -*- coding: utf-8 -*-
"""Проверка решётки с ондуляторами: тюн/оптика ВКЛ vs ВЫКЛ на ОДНОЙ базе (изоляция эффекта ИД).

ВКЛ  — решётка LAT как импортирована (полюса-«жёсткий край»: mat_B и грани θ/2 — изгиб и фокусировка).
ВЫКЛ — та же решётка в рабочей папке с ОБНУЛЁННЫМИ mat_B и гранями полюсов (ондуляторы «разведены в ноль»):
       изгиб и фокус исчезают, остаётся дрейф. Разница тюна = чистый эффект ондуляторов.
Оба прогона — продакшн C++-ядро (ring_cpp), оптика из ring_analysis.twiss (из матриц трекинга).

ЗАПУСК: py -3.14 undulators/verify_undulators.py --lattice 26
"""
import os, sys, csv, shutil, argparse, tempfile
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import ring_cpp as rc
import ring_analysis as ra

TB = os.path.join("workspace", "tables")
POLES = ("UPP", "UPM", "UE0", "UE1")


def zeroed_ws(lat, prefix, out_dir):
    """Рабочая папка с ВЫКЛЮЧЕННЫМИ ондуляторами: полюса без поля и граней (mat_B=0, edge1=edge2=0)."""
    pole_devs = {prefix + n for n in POLES}
    dst = os.path.join(out_dir, "tables")
    os.makedirs(dst, exist_ok=True)
    rows = list(csv.reader(open(os.path.join(TB, "lattice_configs.csv"), encoding="utf-8"), delimiter=";"))
    hdr = rows[0]; ci = {c: i for i, c in enumerate(hdr)}
    for r in rows[1:]:
        if not r or r[ci["lattice_id"]] != str(lat):
            continue
        if r[ci["device_type_id"]] in pole_devs:
            r[ci["mat_B"]] = "0"; r[ci["edge1"]] = "0"; r[ci["edge2"]] = "0"
    with open(os.path.join(dst, "lattice_configs.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(rows[1:])
    for t in ("element_connectivity.csv", "particle_groups.csv"):
        shutil.copyfile(os.path.join(TB, t), os.path.join(dst, t))
    return out_dir


def tunes(ring):
    tw = ra.twiss(ring)
    return tw["nux"], tw["nuy"], tw["betay"].max(), tw["alphac"]


def main():
    ap = argparse.ArgumentParser(description="тюн ВКЛ/ВЫКЛ ондуляторов на одной базе")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--prefix", default="und_")
    a = ap.parse_args()

    ron = rc.CppRing(lattice=a.lattice, nsl=rc.OPTICS_NSL)
    nxon, nyon, byon, acon = tunes(ron)
    # аналитика: Σ по полюсам β_y·(фокусировка граней tanE·θ/L на входе и выходе)/4π при РЕАЛЬНОЙ β_y (ВКЛ)
    import ring_response as rr
    py = ron._python_ring(); _, dy, _ = rr.tune_response(py)
    pole_devs = {a.prefix + n for n in POLES}
    est = sum(-dy[i] * (m.kEdge1 + m.kEdge2) for i, m in enumerate(py.mags) if str(getattr(m, "dev", "")) in pole_devs)
    ron.close()

    tmp = tempfile.mkdtemp(prefix="und_off_")
    try:
        zeroed_ws(a.lattice, a.prefix, tmp)
        roff = rc.CppRing(lattice=a.lattice, ws=tmp, nsl=rc.OPTICS_NSL)
        nxoff, nyoff, byoff, acoff = tunes(roff); roff.close()
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print(f"  решётка {a.lattice}")
    print(f"  ВЫКЛ (ондуляторы в ноль):  νx={nxoff:.5f}  νy={nyoff:.5f}  β_y^max={byoff:.3f}  αc={acoff:.3e}")
    print(f"  ВКЛ  (ондуляторы 5 мм):    νx={nxon:.5f}  νy={nyon:.5f}  β_y^max={byon:.3f}  αc={acon:.3e}")
    print(f"  Δνx={nxon-nxoff:+.5f}  Δνy={nyon-nyoff:+.5f}   (аналитика Σβ_y·(tgE₁+tgE₂)θ/L/4π = {est:+.5f})")


if __name__ == "__main__":
    main()
