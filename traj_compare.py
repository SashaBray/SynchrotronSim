# -*- coding: utf-8 -*-
"""Поэлементная сверка траектории: SynchrotronSim vs elegant, 1 оборот. Вектор состояния на ВХОДЕ и
ВЫХОДЕ каждого элемента + разница ⇒ видно, ГДЕ именно возникает расхождение. Таблицы архивируются.

Гоняет ОБА трекера на одном луче:
  • elegant — генерирует .ele, запускает, читает centroid (Cx,Cxp,Cy,Cyp на выходе каждого элемента);
  • SynchrotronSim — build_line + трек, состояние на входе/выходе каждого магнита.
Выравнивание по продольной s (границы совпадают: наш магнит ↔ элемент elegant). Разница в мкм/мкрад.

ВАЖНО: elegant .lte и БД должны быть ОДНОЙ конфигурацией. По умолчанию — открытая (elegant_run/
syla_ring.lte, RING_Q ↔ решётка 1 open). Луч задаётся одинаковым обоим (--x0/--px0/--y0/--py0/--d0).

Секторный/иной режим ядра — через env (SS_SECTOR=1 и т.п.): скрипт трекает тем же ядром, что и трек3.

Запуск:
  py -3.14 traj_compare.py                              # x0=1мм, 1 оборот, автотег
  py -3.14 traj_compare.py --x0 5 --tag dl_5mm
  SS_SECTOR=1 py -3.14 traj_compare.py --tag sector_1mm
  py -3.14 traj_compare.py --turns 1 --lattice 1 --beamline RING_Q --stop-dq
Архив: workspace/results/traj_compare/<tag>/  (table.csv + manifest.json + summary.txt).
"""
import os, sys, csv, json, subprocess, argparse, re, datetime
import numpy as np
import matrix_tracker as mt

ROOT = os.path.dirname(os.path.abspath(__file__))
ELEG_DIR = os.path.join(ROOT, "elegant_run")
OUT_ROOT = os.path.join(ROOT, "workspace", "results", "traj_compare")


def _git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT,
                              capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return "?"


# --------------------------------------------------------------------- elegant
def run_elegant(lte, beamline, probe, turns):
    """Прогнать elegant на луче probe=(x,xp,y,yp,δ). Вернуть [(s, name, typ, x, xp, y, yp), ...]
    (centroid на выходе каждого элемента)."""
    x0, px0, y0, py0, d0 = probe
    ele = os.path.join(ELEG_DIR, "_trajcmp.ele")
    with open(ele, "w", encoding="utf-8", newline="\n") as f:
        f.write(f'''&run_setup
    lattice = "{os.path.basename(lte)}", use_beamline = "{beamline}",
    p_central_mev = 3000, centroid = "%s.cen", default_order = 1,
&end
&run_control n_passes = {turns}, &end
&bunched_beam n_particles_per_bunch = 1, emit_x=0, emit_y=0, sigma_dp=0, sigma_s=0,
    centroid[0]={x0:.12g}, centroid[1]={px0:.12g}, centroid[2]={y0:.12g},
    centroid[3]={py0:.12g}, centroid[5]={d0:.12g}, &end
&track &end
''')
    r = subprocess.run(["elegant", "_trajcmp.ele"], cwd=ELEG_DIR, capture_output=True, text=True)
    cen = os.path.join(ELEG_DIR, "_trajcmp.cen")
    if not os.path.exists(cen):
        raise SystemExit(f"elegant не создал .cen:\n{r.stdout[-800:]}\n{r.stderr[-800:]}")
    txt = os.path.join(ELEG_DIR, "_trajcmp.txt")
    subprocess.run(["sddsprintout", "_trajcmp.cen", "-col=s", "-col=ElementName", "-col=ElementType",
                    "-col=Cx", "-col=Cxp", "-col=Cy", "-col=Cyp", "-noTitle", "-noLabel", "-width=250"],
                   cwd=ELEG_DIR, stdout=open(txt, "w", encoding="utf-8"))
    out = []
    for line in open(txt, encoding="utf-8"):
        p = line.split()
        if len(p) >= 7:
            try:
                out.append((float(p[0]), p[1], p[2], float(p[3]), float(p[4]), float(p[5]), float(p[6])))
            except ValueError:
                continue
    if not out:
        raise SystemExit("elegant: пустая траектория (проверь конфиг/линию)")
    return out


# --------------------------------------------------------------- SynchrotronSim
def run_ours(lattice, probe, turns, stop_s):
    """Трек луча probe по решётке. Вернуть [(dev, s_in, st_in, s_out, st_out), ...] по магнитам."""
    x0, px0, y0, py0, d0 = probe
    eff, phys, coup, trans = mt.load_db()
    mags = mt.build_line(lattice, eff, phys, coup, trans)
    st = (x0, px0, y0, py0, 0.0, d0); s = 0.0; recs = []
    for j, m in enumerate(mags):
        if j > 0:
            st = mt.modified_drift(st, mags[j - 1], m); s = st[4]
        st_in, s_in = st, s
        st = m.body(st); s = s_in + m.L
        recs.append(("+".join(m.devs), s_in, st_in, s, st))
        if stop_s and s > stop_s:
            break
    return recs


# --------------------------------------------------------------- выравнивание
def eleg_at(eleg_s, eleg_arr, s):
    """(имя, x, xp, y, yp) elegant на продольной s — БЛИЖАЙШИЙ элемент (границы совпадают ⇒ точно)."""
    i = int(np.argmin(np.abs(eleg_s - s)))
    exact = abs(eleg_s[i] - s) < 1e-6
    nm = eleg_arr[i][1] if exact else "~" + eleg_arr[i][1]
    return nm, eleg_arr[i][3], eleg_arr[i][4], eleg_arr[i][5], eleg_arr[i][6]


def build_table(ours, eleg):
    """Таблица: строка на магнит, состояние на входе/выходе, наше и elegant, разница (мкм/мкрад)."""
    es = np.array([e[0] for e in eleg])
    hdr = ["n", "elem_SS", "elem_elegant", "s_in_m", "s_out_m", "faza",
           "x_elegant_mm", "x_SS_mm", "dx_um", "xp_elegant_mrad", "px_SS_mrad", "dpx_urad",
           "y_elegant_mm", "y_SS_mm", "dy_um", "yp_elegant_mrad", "py_SS_mrad", "dpy_urad", "dipol"]
    # ПОЖЕЛАНИЕ пользователя: запятая как десятичный разделитель (удобно в русском Excel). Поле = ';',
    # так что запятая в числах не конфликтует. Касается ТОЛЬКО table.csv (manifest/summary — точка).
    f5 = lambda v: f"{v:+.6f}".replace(".", ",")           # мм/мрад, 6 знаков
    f2 = lambda v: f"{v:+.2f}".replace(".", ",")           # Δ, мкм/мкрад
    fs = lambda v: f"{v:.6f}".replace(".", ",")            # s, м
    rows = []
    for n, (dev, s_in, sti, s_out, sto) in enumerate(ours, 1):
        for faza, s, stt in (("вход", s_in, sti), ("выход", s_out, sto)):
            nm, ex, exp_, ey, eyp = eleg_at(es, eleg, s)
            rows.append([n, dev, nm, fs(s_in), fs(s_out), faza,
                         f5(ex * 1e3), f5(stt[0] * 1e3), f2((stt[0] - ex) * 1e6),
                         f5(exp_ * 1e3), f5(stt[1] * 1e3), f2((stt[1] - exp_) * 1e6),
                         f5(ey * 1e3), f5(stt[2] * 1e3), f2((stt[2] - ey) * 1e6),
                         f5(eyp * 1e3), f5(stt[3] * 1e3), f2((stt[3] - eyp) * 1e6),
                         "ДИПОЛЬ" if ("DP_" in dev or "DQ_" in dev) else ""])
    return hdr, rows


# ------------------------------------------------------------------- архив/CLI
def main():
    ap = argparse.ArgumentParser(description="поэлементная сверка траектории SynchrotronSim vs elegant")
    ap.add_argument("--lattice", type=int, default=1, help="решётка БД (по умолч. 1)")
    ap.add_argument("--lte", default=os.path.join(ELEG_DIR, "syla_ring.lte"), help="elegant .lte")
    ap.add_argument("--beamline", default="RING_Q", help="линия elegant (по умолч. RING_Q)")
    ap.add_argument("--turns", type=int, default=1, help="оборотов")
    ap.add_argument("--x0", type=float, default=1.0, help="старт x, мм")
    ap.add_argument("--px0", type=float, default=0.0, help="старт x', мрад")
    ap.add_argument("--y0", type=float, default=0.0, help="старт y, мм")
    ap.add_argument("--py0", type=float, default=0.0, help="старт y', мрад")
    ap.add_argument("--d0", type=float, default=0.0, help="старт δ=Δp/p")
    ap.add_argument("--stop-dq", action="store_true", help="остановиться перед первым DQ (только чистые диполи)")
    ap.add_argument("--tag", default=None, help="имя архива (по умолч. авто)")
    a = ap.parse_args()

    probe = (a.x0 / 1e3, a.px0 / 1e3, a.y0 / 1e3, a.py0 / 1e3, a.d0)
    mode = "sector" if os.environ.get("SS_SECTOR") == "1" else "kick"
    tag = a.tag or f"{mode}_x{a.x0:g}mm_t{a.turns}_{_git_commit()}"
    stop_s = 8.5 if a.stop_dq else 0.0                     # ~s первого DQ (для --stop-dq)

    print(f"elegant: {os.path.basename(a.lte)} / {a.beamline}, луч x0={a.x0}мм, {a.turns} об...")
    eleg = run_elegant(a.lte, a.beamline, probe, a.turns)
    print(f"  elegant: {len(eleg)} элементов")
    print(f"SynchrotronSim: решётка {a.lattice}, режим ядра '{mode}'...")
    ours = run_ours(a.lattice, probe, a.turns, stop_s)
    print(f"  SynchrotronSim: {len(ours)} магнитов")

    hdr, rows = build_table(ours, eleg)
    outdir = os.path.join(OUT_ROOT, tag); os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "table.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f, delimiter=";"); w.writerow(hdr); w.writerows(rows)

    # сводка: макс |Δ| по не-дипольным выходам (кадры диполей разные — там Δ не показатель).
    # Δx в table.csv — с запятой (для Excel), для арифметики возвращаем к точке.
    _dx = lambda r: abs(float(r[8].replace(",", ".")))
    nd = [r for r in rows if r[5] == "выход" and not r[18]]
    dxs = [_dx(r) for r in nd]
    worst = max(nd, key=_dx) if nd else None
    summary = [f"SynchrotronSim vs elegant — {tag}",
               f"дата: {datetime.datetime.now():%Y-%m-%d %H:%M}, коммит {_git_commit()}, режим ядра {mode}",
               f"луч: x0={a.x0}мм px0={a.px0}мрад y0={a.y0}мм py0={a.py0}мрад δ0={a.d0}, оборотов {a.turns}",
               f"элементов: {len(ours)} (SynchrotronSim) / {len(eleg)} (elegant)",
               f"макс |Δx| по не-дипольным выходам: {max(dxs):.2f} мкм" if dxs else "нет данных",
               f"  на элементе: {worst[1]} (s={worst[4]}м, elegant {worst[2]})" if worst else "",
               f"СКО |Δx| не-дипольных: {np.sqrt(np.mean([v*v for v in dxs])):.2f} мкм" if dxs else ""]
    open(os.path.join(outdir, "summary.txt"), "w", encoding="utf-8").write("\n".join(summary) + "\n")
    json.dump(dict(tag=tag, commit=_git_commit(), mode=mode, lattice=a.lattice, beamline=a.beamline,
                   lte=os.path.basename(a.lte), turns=a.turns, probe_mm=[a.x0, a.px0, a.y0, a.py0, a.d0],
                   date=f"{datetime.datetime.now():%Y-%m-%d %H:%M}",
                   n_ours=len(ours), n_eleg=len(eleg),
                   max_dx_um=max(dxs) if dxs else None, rms_dx_um=float(np.sqrt(np.mean([v*v for v in dxs]))) if dxs else None),
              open(os.path.join(outdir, "manifest.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=2)

    print("\n".join(summary))
    print(f"\nАрхив: {os.path.relpath(outdir, ROOT)}  (table.csv, manifest.json, summary.txt)")


if __name__ == "__main__":
    main()
