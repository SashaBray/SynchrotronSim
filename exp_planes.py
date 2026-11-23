# -*- coding: utf-8 -*-
"""Плоскости опыта и резонаторы установки: авторинг и просмотр.

Опыт (experiments.csv) хранит старт- и финиш-плоскость как точку + нормаль (start_x..start_nz,
finish_x..finish_nz; пустой финиш ⇒ кольцо). Резонаторы установки — rf_cavities.csv (по lattice_id). Писать
глобальные координаты руками неудобно, поэтому плоскость задаётся положением на опорной орбите:

  s (м)               — опорная дуга. Кольцо: s = 0 — середина шва (промежуток между последним и первым
                        элементом; как ring_analysis/elegant), s берётся по модулю периметра. Разомкнутая
                        линия: s = 0 — вход первого элемента, s < 0 — до него, s > конца — после последнего.
  ЭЛЕМЕНТ[#k] + d (м) — d после выхода k-го (с 1) вхождения устройства ЭЛЕМЕНТ.

Плоскость ставится в ПРОМЕЖУТОК (дрейф): точка на опорной прямой, нормаль вдоль неё. s внутри элемента —
ошибка с подсказкой ближайших промежутков. Записываются точные числа (округление до 1e-9 м / 1e-9).

ЗАПУСК
  py -3.14 exp_planes.py --show --exp 26
  py -3.14 exp_planes.py --exp 26 --start-s 0                       # старт в шве, финиш пуст (кольцо)
  py -3.14 exp_planes.py --exp 3 --start-s -0.5 --finish-after QF1J#2 --finish-d 0.1
  py -3.14 exp_planes.py --cavity --lattice 26 --s 555.2 [--share 1] [--id 1] [--note "резонатор 1"]
  py -3.14 exp_planes.py --cavity-remove 1 --lattice 26
"""
import os, sys, csv, io, argparse, contextlib
import numpy as np
import matrix_tracker as mt
import beam_injection as bi

EXP = "experiments.csv"
CAV = "rf_cavities.csv"
ROUND = 9


def _mags(lattice):
    eff, phys, coup, trans = mt.load_db()
    with contextlib.redirect_stdout(io.StringIO()):
        return mt.build_line(lattice, eff, phys, coup, trans)


def _arcs(mags):
    """(ring, lead, s_in, s_out, C): опорные дуги входа/выхода элементов (как ring_analysis.Ring)."""
    ring = bi.seam_is_drift(mags)
    m0, m1 = mags[0], mags[-1]
    lead = float((m0.fr_in[0] - m1.fr_out[0]) @ m1.fr_out[1]) / 2.0 if ring else 0.0
    s = lead; s_in = []; s_out = []
    for j, m in enumerate(mags):
        if j:
            s += float((m.fr_in[0] - mags[j - 1].fr_out[0]) @ mags[j - 1].fr_out[1])
        s_in.append(s); s += m.L; s_out.append(s)
    C = (s + lead) if ring else s
    return ring, lead, np.array(s_in), np.array(s_out), C


def _clean(v):
    return np.round(np.asarray(v, float), ROUND) + 0.0


def plane_at_s(lattice, s, mags=None):
    """Плоскость на опорной орбите в точке s (в промежутке) → (точка, нормаль)."""
    mags = mags if mags is not None else _mags(lattice)
    ring, lead, s_in, s_out, C = _arcs(mags)
    s = float(s)
    if ring:
        s = s % C
    N = len(mags)
    for j in range(N):
        if s_in[j] + 1e-12 < s < s_out[j] - 1e-12:
            gaps = []
            if j:
                gaps.append(f"[{s_out[j - 1]:.4f}; {s_in[j]:.4f}]")
            if j + 1 < N:
                gaps.append(f"[{s_out[j]:.4f}; {s_in[j + 1]:.4f}]")
            raise ValueError(f"s = {s:.6f} м внутри элемента {mags[j].dev} [{s_in[j]:.4f}; {s_out[j]:.4f}]; "
                             f"ближайшие промежутки: {', '.join(gaps)}")
    if s <= s_in[0]:                                         # до первого элемента (у кольца — вторая половина шва)
        m = mags[0]
        return _clean(m.fr_in[0] - (s_in[0] - s) * m.fr_in[1]), _clean(m.fr_in[1])
    j = int(np.searchsorted(s_out, s, side="right")) - 1     # последний элемент, вышедший до s
    m = mags[j]
    return _clean(m.fr_out[0] + (s - s_out[j]) * m.fr_out[1]), _clean(m.fr_out[1])


def plane_after(lattice, dev, d=0.0, mags=None):
    """Плоскость на d метров после выхода устройства dev ('QF1J' или 'QF1J#2' — второе вхождение)."""
    mags = mags if mags is not None else _mags(lattice)
    name, _, k = str(dev).partition("#")
    k = int(k or 1)
    hits = [j for j, m in enumerate(mags) if m.dev == name]
    if len(hits) < k:
        raise ValueError(f"в решётке {lattice} {len(hits)} вхождений {name}, запрошено #{k}")
    _ring, _lead, _si, s_out, _C = _arcs(mags)
    return plane_at_s(lattice, float(s_out[hits[k - 1]]) + float(d), mags)


def default_start(lattice):
    """Старт по умолчанию: s = 0 (у кольца — середина шва, у линии — вход первого элемента)."""
    return plane_at_s(lattice, 0.0)


def default_planes(lattice):
    """(старт, финиш) по умолчанию: кольцо — s = 0 и пустой финиш; линия — вход первого и выход последнего."""
    mags = _mags(lattice)
    ring, _lead, _si, s_out, _C = _arcs(mags)
    return plane_at_s(lattice, 0.0, mags), (None if ring else plane_at_s(lattice, float(s_out[-1]), mags))


# ---------------- запись ----------------
def _read(name):
    p = os.path.join(mt.TB, name)
    with open(p, encoding="utf-8") as f:
        rows = list(csv.reader(f, delimiter=";"))
    return rows[0], [dict(zip(rows[0], r)) for r in rows[1:]]


def _write(name, hdr, rows):
    with open(os.path.join(mt.TB, name), "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=hdr, delimiter=";", extrasaction="ignore")
        w.writeheader(); w.writerows(rows)


def _f(v):
    return repr(float(v) + 0.0)


def plane_cells(prefix, plane):
    """(точка, нормаль) | None → ячейки столбцов prefix_x..prefix_nz."""
    keys = [f"{prefix}_{c}" for c in ("x", "y", "z", "nx", "ny", "nz")]
    if plane is None:
        return dict.fromkeys(keys, "")
    c, n = plane
    return dict(zip(keys, [_f(v) for v in list(c) + list(n)]))


def write_experiment(exp_id, lattice=None, group=None, start=None, finish="keep", **cells):
    """Создать/обновить строку опыта. start=None — оставить (новой строке — default_planes), finish: 'keep' |
    None (кольцо) | (точка, нормаль). cells — прочие столбцы (dt, max_turns, loss_radius_mm, …)."""
    hdr, rows = _read(EXP)
    row = next((r for r in rows if r["exp_id"] == str(exp_id)), None)
    if row is None:
        if lattice is None or group is None:
            raise ValueError(f"опыта {exp_id} нет — задайте решётку и группу")
        row = dict.fromkeys(hdr, ""); row.update(exp_id=str(exp_id), dt="1.00E-11", max_turns="0", traj_stride="50")
        rows.append(row)
        if start is None:
            start, fin0 = default_planes(int(lattice))
            if finish == "keep":
                finish = fin0
        if finish == "keep":
            finish = None
    if lattice is not None:
        row["lattice_id"] = str(lattice)
    if group is not None:
        row["group_id"] = str(group)
    if start is not None:
        row.update(plane_cells("start", start))
    if finish != "keep":
        row.update(plane_cells("finish", finish))
    row.update({k: str(v) for k, v in cells.items() if k in hdr})
    rows.sort(key=lambda r: int(r["exp_id"]))
    _write(EXP, hdr, rows)
    return row


def write_cavity(lattice, plane, share=None, cavity_id=None, note=""):
    p = os.path.join(mt.TB, CAV)
    if not os.path.exists(p):
        with open(p, "w", newline="", encoding="utf-8") as f:
            csv.writer(f, delimiter=";").writerow(bi.CAVITY_COLS)
    hdr, rows = _read(CAV)
    if cavity_id is None:
        cavity_id = 1 + max([int(r["cavity_id"]) for r in rows if r["cavity_id"].isdigit()] or [0])
    rows = [r for r in rows if not (r["cavity_id"] == str(cavity_id) and r["lattice_id"] == str(lattice))]
    c, n = plane
    rows.append(dict(cavity_id=str(cavity_id), lattice_id=str(lattice), x=_f(c[0]), y=_f(c[1]), z=_f(c[2]),
                     nx=_f(n[0]), ny=_f(n[1]), nz=_f(n[2]), volt_share="" if share is None else repr(float(share)),
                     note=note))
    rows.sort(key=lambda r: (int(r["lattice_id"]), r["cavity_id"]))
    _write(CAV, hdr, rows)
    return cavity_id


def remove_cavity(lattice, cavity_id):
    hdr, rows = _read(CAV)
    keep = [r for r in rows if not (r["cavity_id"] == str(cavity_id) and r["lattice_id"] == str(lattice))]
    if len(keep) == len(rows):
        raise ValueError(f"резонатора {cavity_id} у решётки {lattice} нет")
    _write(CAV, hdr, keep)


# ---------------- просмотр ----------------
def _where(mags, plane, arcs):
    ring, lead, s_in, s_out, C = arcs
    slot, along = bi.locate_plane(mags, *plane)
    N = len(mags)
    if slot < 0:
        return f"до первого элемента {mags[0].dev} (s = {s_in[0] + along:.6f} м)"
    s = s_out[slot] + along
    if ring:                                                 # кольцо: s в (−C/2; C/2] — шов около нуля
        s = s % C
        s = s - C if s > C / 2 else s
    nxt = mags[(slot + 1) % N].dev if (slot + 1 < N or ring) else "конец линии"
    return f"после {mags[slot].dev} (элемент {slot + 1}), перед {nxt}; s = {s:.6f} м"


def show(exp_id):
    row = bi.experiment_row(exp_id)
    lat, gid = int(row["lattice_id"]), int(row["group_id"])
    mags = _mags(lat); arcs = _arcs(mags)
    start, finish = bi.planes_of(row)
    print(f"опыт {exp_id}: решётка {lat} ({'кольцо' if arcs[0] else 'разомкнутая линия'}, "
          f"опорная длина {arcs[4]:.6f} м), группа {gid} ({len(bi.load_group(gid))} частиц)")
    print(f"  старт : точка {start[0].tolist()} нормаль {start[1].tolist()}\n          {_where(mags, start, arcs)}")
    if finish is None:
        print("  финиш : пусто ⇒ = старт (кольцо, полный оборот)")
    else:
        print(f"  финиш : точка {finish[0].tolist()} нормаль {finish[1].tolist()}\n          {_where(mags, finish, arcs)}")
    cav = bi.load_cavities(lat)
    path = bi.build_path(mags, start, finish, cav)
    kind = "КОЛЬЦО" if path["closed"] else "РАЗОМКНУТ (портал: финиш → старт)"
    print(f"  путь  : {kind}, элементов {sum(1 for m in path['ordered'] if not m.is_cavity)}, "
          f"длина пролёта {path['line_len']:.6f} м")
    if not cav:
        print("  ВЧ    : у решётки нет резонаторов (rf_cavities.csv) ⇒ ВЧ в конце пролёта")
    for c in path["cavities"]:
        print(f"  ВЧ    : резонатор {c.cavity_id}, s = {c.s_c:.6f} м от старта, доля {c.share:g}")
    if cav and path["end_share"]:
        print(f"  ВЧ    : резонатор(ы) в плоскости финиша, доля {path['end_share']:g}")
    for cid in path["skipped"]:
        print(f"  ВЧ    : резонатор {cid} вне пути — не работает в этом опыте")
    rad = mt.load_radiation(exp_id)
    print(f"  radiation_params: {('rf=' + rad['rf']) if rad else 'нет строки (СИ и ВЧ выкл)'}")


def _plane_arg(lat, s, after, d, mags):
    if after:
        return plane_after(lat, after, d or 0.0, mags)
    return plane_at_s(lat, s, mags) if s is not None else None


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ap = argparse.ArgumentParser(description="плоскости опыта и резонаторы установки")
    ap.add_argument("--exp", type=int)
    ap.add_argument("--show", action="store_true")
    ap.add_argument("--lattice", type=int)
    ap.add_argument("--group", type=int)
    ap.add_argument("--start-s", type=float)
    ap.add_argument("--start-after"); ap.add_argument("--start-d", type=float)
    ap.add_argument("--finish-s", type=float)
    ap.add_argument("--finish-after"); ap.add_argument("--finish-d", type=float)
    ap.add_argument("--ring", action="store_true", help="финиш пуст (= старт, кольцо)")
    ap.add_argument("--cavity", action="store_true", help="добавить/заменить резонатор установки")
    ap.add_argument("--cavity-remove", type=int)
    ap.add_argument("--s", type=float); ap.add_argument("--after"); ap.add_argument("--d", type=float)
    ap.add_argument("--share", type=float); ap.add_argument("--id", type=int); ap.add_argument("--note", default="")
    a = ap.parse_args()
    if a.cavity or a.cavity_remove is not None:
        if a.lattice is None:
            ap.error("--lattice обязателен для резонатора")
        if a.cavity_remove is not None:
            remove_cavity(a.lattice, a.cavity_remove); print(f"резонатор {a.cavity_remove} решётки {a.lattice} удалён")
            return
        mags = _mags(a.lattice)
        pl = _plane_arg(a.lattice, a.s, a.after, a.d, mags)
        if pl is None:
            ap.error("положение резонатора: --s или --after [--d]")
        cid = write_cavity(a.lattice, pl, a.share, a.id, a.note)
        print(f"резонатор {cid} решётки {a.lattice}: точка {pl[0].tolist()} нормаль {pl[1].tolist()}")
        return
    if a.exp is None:
        ap.error("--exp обязателен")
    if any(v is not None for v in (a.start_s, a.start_after, a.finish_s, a.finish_after)) or a.ring:
        row = None
        try:
            row = bi.experiment_row(a.exp)
        except ValueError:
            pass
        lat = a.lattice if a.lattice is not None else (int(row["lattice_id"]) if row else None)
        if lat is None:
            ap.error("новый опыт: задайте --lattice и --group")
        mags = _mags(lat)
        start = _plane_arg(lat, a.start_s, a.start_after, a.start_d, mags)
        fin = _plane_arg(lat, a.finish_s, a.finish_after, a.finish_d, mags)
        finish = None if a.ring else (fin if fin is not None else "keep")
        write_experiment(a.exp, a.lattice, a.group, start, finish)
    show(a.exp)


if __name__ == "__main__":
    main()
