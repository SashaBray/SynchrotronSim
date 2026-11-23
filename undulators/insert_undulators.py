# -*- coding: utf-8 -*-
"""ДОСТОВЕРНАЯ вставка ондуляторов в КОПИЮ выверенной решётки (сохраняет кривую оптику β).

В отличие от плоского переимпорта (import_undulators.py, искажает β вдвое), здесь берётся выверенная
решётка SRC (по умолч. 1 = СИЛА) ДОСЛОВНО (вся физика mat_*/edge/kind/геометрия — как есть, β_y=6.52 м),
и в её прямолинейные промежутки ВСТАВЛЯЮТСЯ строки полюсов-«жёсткий край» (kind=dipole, длина (8/π²)·λu/2,
грани θ/2 — фокусировка самими полюсами, без отдельной линзы). Геометрия вставки — центрами вдоль оси промежутка (прямой
дрейф: центр строки = x/y/z, ось из yaw/pitch/roll, см. matrix_tracker.Magnet). Типы und_* переиспользуются
из БД (их заводит import_undulators при первом прогоне).

ВАЖНО: build_line упорядочивает кольцо по instance_id (sorted), поэтому ВСЕ инстансы новой решётки
перенумеровываются 1..M в порядке кольца (база + вставки). Со-локаций в реш.1 нет — перенумерация безопасна.

ЗАПУСК: py -3.14 undulators/insert_undulators.py --src 1 --dst 27 --gap 5
"""
import os, sys, csv, argparse
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
os.environ.setdefault("SS_SECTOR", "1")

import numpy as np
import matrix_tracker as mt
import undulator_sr as sr
from assemble_undulators import read_device, resolve_pos

TB = os.path.join("workspace", "tables")
LC = os.path.join(TB, "lattice_configs.csv")
EC = os.path.join(TB, "element_connectivity.csv")
EX = os.path.join(TB, "experiments.csv")
GAP_MIN = 5.0


def _read(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter=";"))
    return rows[0], rows[1:]


def _rewrite_lattice(path, hdr, lat, new_rows):
    """Заменить строки решётки lat (чужие сохранить), записать новые."""
    _, rows = _read(path)
    li = hdr.index("lattice_id")
    keep = [r for r in rows if r and r[li] != str(lat)]
    import shutil; shutil.copy(path, path + ".bak_ins")
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(keep + new_rows)


def main():
    ap = argparse.ArgumentParser(description="достоверная вставка ондуляторов в копию выверенной решётки")
    ap.add_argument("--src", type=int, default=1, help="исходная выверенная решётка")
    ap.add_argument("--dst", type=int, default=27, help="id новой решётки с ондуляторами")
    ap.add_argument("--gap", type=float, default=5.0, help="зазор всех ондуляторов, мм")
    ap.add_argument("--device", default=os.path.join(os.path.dirname(__file__), "undulator_cpmu18.lte"))
    ap.add_argument("--manifest", default=os.path.join(os.path.dirname(__file__), "placements.csv"))
    a = ap.parse_args()
    QP = mt.QP; Brho = sr.P.E / 0.299792458

    hdr, all_rows = _read(LC)
    ci = {c: i for i, c in enumerate(hdr)}
    base = sorted([r for r in all_rows if r and r[ci["lattice_id"]] == str(a.src)],
                  key=lambda r: int(r[ci["instance_id"]]))
    if not base:
        raise SystemExit(f"решётка {a.src} пуста")

    delems, dline, L_dev, b0_dev = read_device(a.device)
    dseq = [t for t in dline if delems[t]["type"] != "MARK"]            # полюса и промежутки по порядку
    poles = [t for t in dseq if delems[t]["type"] == "CSBEND"]
    L_foot = L_dev
    B0 = sr.B0_of_gap(a.gap, sr.P); k = B0 / b0_dev
    Ky = 0.5 * (B0 / Brho) ** 2

    manifest = list(csv.DictReader(open(a.manifest, encoding="utf-8"), delimiter=";"))
    pos_by_order = [r["pos"] for r in sorted(manifest, key=lambda r: float(r["s_start"]))]

    def C(r): return np.array([float(r[ci["x"]] or 0), float(r[ci["y"]] or 0), float(r[ci["z"]] or 0)])
    def L(r): return float(r[ci["eff_length"]] or 0)

    # найти промежутки ≥ GAP_MIN между соседними (в порядке кольца)
    straights = []
    for i in range(len(base) - 1):
        d = float(np.linalg.norm(C(base[i + 1]) - C(base[i]))) - 0.5 * (L(base[i]) + L(base[i + 1]))
        if d >= GAP_MIN:
            straights.append((i, d))
    if len(straights) != len(pos_by_order):
        print(f"  [!] промежутков {len(straights)} ≠ слотов манифеста {len(pos_by_order)} — беру min")

    blank = {c: "" for c in hdr}

    def mkrow(dev, kind, center, yaw, pitch, roll, Leff, nsl, matB=0.0, matG=0.0, matGv=0.0, edge=0.0):
        r = dict(blank)
        r.update(lattice_id=str(a.dst), device_type_id=dev, kind=kind, n_slices=str(nsl),
                 x=f"{center[0]:.12f}", y=f"{center[1]:.12f}", z=f"{center[2]:.12f}",
                 yaw=yaw, pitch=pitch, roll=roll, arg_val="0", correct_coeff="1",
                 eff_length=f"{Leff:.12f}", mat_B=f"{matB:.15g}", mat_G=f"{matG:.15g}",
                 mat_S="0", mat_O="0", mat_Gv=f"{matGv:.15g}", edge1=f"{edge:.15g}", edge2=f"{edge:.15g}")
        return [r[c] for c in hdr]

    # собрать упорядоченный список строк новой решётки (база копией + вставки), затем перенумеровать
    out_rows = []          # строки в порядке кольца
    groups = []            # (slot_idx, [индексы вставленных строк в out_rows])
    straight_at = {i: j for j, (i, d) in enumerate(straights)}   # base-индекс → порядковый № промежутка
    L_slot_at = dict(straights)
    for i, br in enumerate(base):
        row = list(br); row[ci["lattice_id"]] = str(a.dst)
        out_rows.append(row)
        if i not in straight_at:
            continue
        sidx = straight_at[i]
        rowA, rowB = base[i], base[i + 1]
        yaw, pitch, roll = rowA[ci["yaw"]], rowA[ci["pitch"]], rowA[ci["roll"]]
        axis = C(rowB) - C(rowA); axis = axis / np.linalg.norm(axis)        # направление промежутка (к B)
        P0 = C(rowA) + 0.5 * L(rowA) * axis                                 # выход предыдущего магнита
        aoff = resolve_pos(pos_by_order[sidx] if sidx < len(pos_by_order) else "center",
                           L_slot_at[i], L_foot, 0.05)
        member = []
        s_run = 0.0                                                      # путь от входа устройства
        for t in dseq:
            Lt = delems[t].get("L", 0.0)
            if delems[t]["type"] == "CSBEND":
                th = -delems[t].get("ANGLE", 0.0) * k                    # угол ядра θ0 = QP·mat_B·L (знак БД)
                cen = P0 + (aoff + s_run + 0.5 * Lt) * axis
                out_rows.append(mkrow("und_" + t, "dipole", cen, yaw, pitch, roll, Lt, 1,
                                      matB=th / (QP * Lt), edge=0.5 * th))   # грань = θ/2 (прямоугольный полюс)
                member.append(len(out_rows) - 1)
            s_run += Lt
        groups.append((sidx, member))

    # перенумеровать instance_id = 1..M в порядке кольца
    for newid, row in enumerate(out_rows, 1):
        row[ci["instance_id"]] = str(newid)

    _rewrite_lattice(LC, hdr, a.dst, out_rows)

    # connectivity: один ряд на инстанс, порядок = instance_id
    ec_hdr = ["lattice_id", "instance_id", "device_type_id", "prev_instance_id", "next_instance_id"]
    ec_rows = []
    M = len(out_rows)
    for j, row in enumerate(out_rows):
        iid = j + 1
        ec_rows.append([str(a.dst), str(iid), row[ci["device_type_id"]],
                        str(iid - 1) if j > 0 else "", str(iid + 1) if j < M - 1 else ""])
    _rewrite_lattice(EC, ec_hdr, a.dst, ec_rows)

    # experiments: exp dst (1 оборот)
    ex_hdr, ex_rows = _read(EX)
    exi = ex_hdr.index("exp_id")
    keep = [r for r in ex_rows if r and r[exi] != str(a.dst)]
    tmpl = {c: "" for c in ex_hdr}
    tmpl.update(exp_id=str(a.dst), lattice_id=str(a.dst), group_id="1", dt="6,74E-13", max_turns="1", traj_stride="1")
    import shutil; shutil.copy(EX, EX + ".bak_ins")
    with open(EX, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(ex_hdr)
        w.writerows(keep); w.writerow([tmpl[c] for c in ex_hdr])

    # element_groups: группа на устройство (полюса), pivot = среднее положение (центр устройства)
    eg_hdr = ["group_id", "lattice_id", "kind", "pivot_x", "pivot_y", "pivot_z", "yaw", "pitch", "roll", "instance_ids"]
    eg_rows = []
    for sidx, member in sorted(groups):
        insts = [out_rows[j][ci["instance_id"]] for j in member]
        cs = np.array([[float(out_rows[j][ci["x"]]), float(out_rows[j][ci["y"]]), float(out_rows[j][ci["z"]])]
                       for j in member])
        piv = cs.mean(axis=0); mid = out_rows[member[len(member) // 2]]
        eg_rows.append([f"und_{sidx+1:02d}", str(a.dst), "undulator", f"{piv[0]:.12f}", f"{piv[1]:.12f}",
                        f"{piv[2]:.12f}", mid[ci["yaw"]], mid[ci["pitch"]], mid[ci["roll"]], "|".join(insts)])
    egp = os.path.join(TB, "element_groups.csv")
    if os.path.exists(egp):
        _rewrite_lattice(egp, eg_hdr, a.dst, eg_rows)
    else:
        with open(egp, "w", newline="", encoding="utf-8") as f:
            w = csv.writer(f, delimiter=";"); w.writerow(eg_hdr); w.writerows(eg_rows)

    import exp_planes as ep                                    # явные плоскости опыта: старт s = 0 (шов), кольцо
    st, fin = ep.default_planes(int(a.dst))
    ep.write_experiment(a.dst, start=st, finish=fin)
    print(f"  решётка {a.dst} = копия реш.{a.src} ({len(base)} строк) + {len(groups)} ондуляторов "
          f"({len(poles)} полюсов-«жёсткий край» каждый); всего {M} инстансов")
    print(f"  k={k:.4f}, справочно k_y={Ky:.4e}, ½(B0/Bρ)²·L={Ky*L_dev:.4e} (фокусируют грани полюсов)")
    print(f"[OK] lattice_configs/element_connectivity/experiments/element_groups ← решётка {a.dst}")


if __name__ == "__main__":
    main()
