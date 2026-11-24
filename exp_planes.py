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

ИМПУЛЬСНЫЕ КИКЕРЫ: «железо» — у установки (pulsed_kickers.csv), импульс — у опыта (kicker_timing.csv):
  py -3.14 exp_planes.py --kicker --lattice 26 --after QF1J#2 --d 0.2 --shape half_sine --duration 2e-6 [--rise 0]
  py -3.14 exp_planes.py --fire --exp 26 --kicker-id 1 --theta-x 1e-3 --on-pass 0 [--offset 0] [--period 1e-3] [--add]
  py -3.14 exp_planes.py --fire --exp 26 --kicker-id 1 --theta-x 1e-3 --delay 1.5e-6
  py -3.14 exp_planes.py --unfire --exp 26 --kicker-id 1
  py -3.14 exp_planes.py --kicker-remove 1 --lattice 26
--on-pass n ставит пик импульса на приход опорной частицы на n-м пролёте: T_ref = (n·C_pass + s_k/β0)/c,
C_pass = C_RF у sine-ВЧ опыта, иначе L/β0 (β0 — по проектной энергии опыта). --add — ещё один импульс того же
кикера (иначе строки опыта с этим кикером заменяются).
"""
import os, sys, csv, io, argparse, contextlib
import numpy as np
import matrix_tracker as mt
import beam_injection as bi

EXP = "experiments.csv"
CAV = "rf_cavities.csv"
KIK = "pulsed_kickers.csv"
TIM = "kicker_timing.csv"
ROUND = 9
C_LIGHT = 299792458.0
MC2_MEV = 0.51099895069


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


def _ensure(name, cols):
    p = os.path.join(mt.TB, name)
    if not os.path.exists(p):
        with open(p, "w", newline="", encoding="utf-8") as f:
            csv.writer(f, delimiter=";").writerow(cols)


def write_kicker(lattice, plane, duration, shape="half_sine", rise=None, kicker_id=None, note=""):
    """Импульсный кикер установки (pulsed_kickers.csv): плоскость, форма half_sine | trapezoid, длительность
    основания duration (с), фронт = спад rise (с, трапеция). kicker_id=None ⇒ следующий свободный."""
    if shape not in bi.PULSE_SHAPES:
        raise ValueError(f"форма {shape!r}: допустимо {', '.join(bi.PULSE_SHAPES)}")
    if not float(duration) > 0:
        raise ValueError("длительность импульса должна быть > 0")
    if shape == "trapezoid" and rise is not None and not (0 <= float(rise) <= float(duration) / 2):
        raise ValueError("фронт трапеции должен быть в [0; duration/2]")
    _ensure(KIK, bi.KICKER_COLS)
    hdr, rows = _read(KIK)
    if kicker_id is None:
        kicker_id = 1 + max([int(r["kicker_id"]) for r in rows if r["kicker_id"].isdigit()] or [0])
    rows = [r for r in rows if not (r["kicker_id"] == str(kicker_id) and r["lattice_id"] == str(lattice))]
    c, n = plane
    rows.append(dict(kicker_id=str(kicker_id), lattice_id=str(lattice), x=_f(c[0]), y=_f(c[1]), z=_f(c[2]),
                     nx=_f(n[0]), ny=_f(n[1]), nz=_f(n[2]), shape=shape, duration_s=repr(float(duration)),
                     rise_s="" if rise is None else repr(float(rise)), note=note))
    rows.sort(key=lambda r: (int(r["lattice_id"]), r["kicker_id"]))
    _write(KIK, hdr, rows)
    return kicker_id


def remove_kicker(lattice, kicker_id):
    hdr, rows = _read(KIK)
    keep = [r for r in rows if not (r["kicker_id"] == str(kicker_id) and r["lattice_id"] == str(lattice))]
    if len(keep) == len(rows):
        raise ValueError(f"кикера {kicker_id} у решётки {lattice} нет")
    _write(KIK, hdr, keep)


def kicker_row(lattice, kicker_id):
    for r in bi._rows(KIK):
        if r.get("kicker_id") == str(kicker_id) and r.get("lattice_id") == str(lattice):
            return r
    raise ValueError(f"кикера {kicker_id} у решётки {lattice} нет в {KIK}")


def reference_time(exp_id, kicker_id, n_pass):
    """Момент прихода опорной частицы к кикеру на пролёте n_pass (с): T_ref = (n·C_pass + s_k/β0)/c,
    C_pass = C_RF у sine-ВЧ опыта (h·c/f или L), иначе L/β0; β0 — по проектной энергии опыта."""
    row = bi.experiment_row(exp_id); lat = int(row["lattice_id"])
    h = kicker_row(lat, kicker_id)
    start, finish = bi.planes_of(row)
    kk = dict(kicker_id=str(kicker_id), center=bi._vec(h, ("x", "y", "z", "nx", "ny", "nz"))[0],
              normal=bi._vec(h, ("x", "y", "z", "nx", "ny", "nz"))[1], shape="half_sine", duration=1.0, rise=0.0,
              theta_x=0.0, theta_y=0.0, delay=0.0, period=None)
    path = bi.build_path(_mags(lat), start, finish, kickers=[kk])
    nodes = path["kick_begin"] + path["kickers"] + path["kick_end"]
    if not nodes:
        raise ValueError(f"кикер {kicker_id} вне пути старт → финиш опыта {exp_id}")
    rad = mt.load_radiation(exp_id)
    design = rad["design_mev"] if rad else 6000.0
    beta0 = design / np.hypot(design, MC2_MEV)
    L = path["line_len"]
    spec = mt.rf_spec(rad) if rad else None
    if isinstance(spec, dict) and spec.get("type") == "sine":
        c_pass = spec["harmonic"] * C_LIGHT / spec["freq_hz"] if spec.get("freq_hz") else L
    else:
        c_pass = L / beta0
    return (n_pass * c_pass + nodes[0].s_k / beta0) / C_LIGHT


def write_timing(exp_id, kicker_id, theta_x=0.0, theta_y=0.0, delay=0.0, period=None, note="", add=False):
    """Импульс кикера в опыте (kicker_timing.csv). add=False — строки (опыт, кикер) заменяются этой."""
    row = bi.experiment_row(exp_id)
    kicker_row(int(row["lattice_id"]), kicker_id)                # кикер должен быть у решётки опыта
    _ensure(TIM, bi.TIMING_COLS)
    hdr, rows = _read(TIM)
    if not add:
        rows = [r for r in rows if not (r["exp_id"] == str(exp_id) and r["kicker_id"] == str(kicker_id))]
    rows.append(dict(exp_id=str(exp_id), kicker_id=str(kicker_id), theta_x_rad=repr(float(theta_x)),
                     theta_y_rad=repr(float(theta_y)), delay_s=repr(float(delay)),
                     period_s="" if period is None else repr(float(period)), note=note))
    rows.sort(key=lambda r: (int(r["exp_id"]), r["kicker_id"]))
    _write(TIM, hdr, rows)


def remove_timing(exp_id, kicker_id):
    hdr, rows = _read(TIM)
    keep = [r for r in rows if not (r["exp_id"] == str(exp_id) and r["kicker_id"] == str(kicker_id))]
    if len(keep) == len(rows):
        raise ValueError(f"у опыта {exp_id} нет импульсов кикера {kicker_id}")
    _write(TIM, hdr, keep)


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
    kk = bi.load_kickers(lat, exp_id)
    path = bi.build_path(mags, start, finish, cav, kickers=kk)
    kind = "КОЛЬЦО" if path["closed"] else "РАЗОМКНУТ (портал: финиш → старт)"
    print(f"  путь  : {kind}, элементов {sum(1 for m in path['ordered'] if not bi.is_node(m))}, "
          f"длина пролёта {path['line_len']:.6f} м")
    if not cav:
        print("  ВЧ    : у решётки нет резонаторов (rf_cavities.csv) ⇒ ВЧ в конце пролёта")
    for c in path["cavities"]:
        print(f"  ВЧ    : резонатор {c.cavity_id}, s = {c.s_c:.6f} м от старта, доля {c.share:g}")
    if cav and path["end_share"]:
        print(f"  ВЧ    : резонатор(ы) в плоскости финиша, доля {path['end_share']:g}")
    for cid in path["skipped"]:
        print(f"  ВЧ    : резонатор {cid} вне пути — не работает в этом опыте")
    where = {id(k): "в старте (начало пролёта)" for k in path["kick_begin"]}
    where.update({id(k): "в плоскости финиша (конец пролёта)" for k in path["kick_end"]})
    for k in path["kick_begin"] + path["kickers"] + path["kick_end"]:
        per = f", период {k.period:g} с" if k.period is not None else ""
        rise = f" фронт {k.rise:g} с" if k.shape == "trapezoid" else ""
        print(f"  кикер : {k.kicker_id}, s = {k.s_k:.6f} м от старта {where.get(id(k), '')}; θ = ({k.theta_x:g}, "
              f"{k.theta_y:g}) рад, {k.shape} {k.duration:g} с{rise}, задержка {k.delay:g} с{per}")
    for kid in path["kick_skipped"]:
        print(f"  кикер : {kid} вне пути — не стреляет в этом опыте")
    idle = [r["kicker_id"] for r in bi._rows(KIK) if r.get("lattice_id") == str(lat)
            and r["kicker_id"] not in {k["kicker_id"] for k in kk}]
    if idle:
        print(f"  кикеры решётки без импульса в этом опыте: {', '.join(idle)}")
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
    ap.add_argument("--kicker", action="store_true", help="добавить/заменить импульсный кикер установки")
    ap.add_argument("--kicker-remove", type=int)
    ap.add_argument("--shape", default="half_sine", choices=bi.PULSE_SHAPES)
    ap.add_argument("--duration", type=float, help="длительность основания импульса, с")
    ap.add_argument("--rise", type=float, help="фронт = спад трапеции, с")
    ap.add_argument("--fire", action="store_true", help="импульс кикера в опыте (kicker_timing.csv)")
    ap.add_argument("--unfire", action="store_true")
    ap.add_argument("--kicker-id", type=int, help="номер кикера для --fire/--unfire")
    ap.add_argument("--theta-x", type=float, default=0.0); ap.add_argument("--theta-y", type=float, default=0.0)
    ap.add_argument("--delay", type=float, help="начало импульса, с от старта опыта")
    ap.add_argument("--on-pass", type=int, help="пик импульса — на приход опорной частицы на этом пролёте")
    ap.add_argument("--offset", type=float, default=0.0, help="сдвиг пика относительно --on-pass, с")
    ap.add_argument("--period", type=float); ap.add_argument("--add", action="store_true")
    a = ap.parse_args()
    if a.fire or a.unfire:
        if a.exp is None or a.kicker_id is None:
            ap.error("--fire/--unfire: нужны --exp и --kicker-id")
        if a.unfire:
            remove_timing(a.exp, a.kicker_id); print(f"импульсы кикера {a.kicker_id} в опыте {a.exp} удалены")
        else:
            if (a.delay is None) == (a.on_pass is None):
                ap.error("--fire: задайте ровно одно из --delay и --on-pass")
            delay = a.delay
            if a.on_pass is not None:
                lat = int(bi.experiment_row(a.exp)["lattice_id"])
                D = float(kicker_row(lat, a.kicker_id)["duration_s"])
                t_ref = reference_time(a.exp, a.kicker_id, a.on_pass)
                delay = t_ref - D / 2.0 + a.offset
                print(f"приход опорной частицы к кикеру {a.kicker_id} на пролёте {a.on_pass}: {t_ref:.12g} с ⇒ "
                      f"задержка {delay:.12g} с")
            write_timing(a.exp, a.kicker_id, a.theta_x, a.theta_y, delay, a.period, a.note, a.add)
        show(a.exp)
        return
    if a.kicker or a.kicker_remove is not None:
        if a.lattice is None:
            ap.error("--lattice обязателен для кикера")
        if a.kicker_remove is not None:
            remove_kicker(a.lattice, a.kicker_remove); print(f"кикер {a.kicker_remove} решётки {a.lattice} удалён")
            return
        if a.duration is None:
            ap.error("--kicker: нужна --duration (с)")
        mags = _mags(a.lattice)
        pl = _plane_arg(a.lattice, a.s, a.after, a.d, mags)
        if pl is None:
            ap.error("положение кикера: --s или --after [--d]")
        kid = write_kicker(a.lattice, pl, a.duration, a.shape, a.rise, a.id, a.note)
        print(f"импульсный кикер {kid} решётки {a.lattice}: точка {pl[0].tolist()} нормаль {pl[1].tolist()}, "
              f"{a.shape} {a.duration:g} с")
        return
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
