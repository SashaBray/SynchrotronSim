# -*- coding: utf-8 -*-
"""Сборка решётки с ондуляторами: базовый .lte + устройство + манифест → составной .lte.

ЧТО ДЕЛАЕТ. Для каждого АКТИВНОГО (on=1) слота манифеста placements.csv:
  • считает положение входа устройства в промежутке из pos (start|center|end|доля|метры);
  • ПРОВЕРЯЕТ, что устройство (+ зазор безопасности) помещается в промежуток — иначе отказ;
  • режет прогон дрейфов слота [i_first..i_last] на DRIF(a) + <полюса устройства> + DRIF(b).

Полюса — ОБЩИЕ имена (UPP/UPM/UE0/UE1), определяются ОДИН раз (k=1); несут ИЗГИБ-орбиту, путь и увод
при перекосе (kind=dipole, секторный интегратор). ФОКУС устройства несёт отдельный элемент UFC (KQUAD,
один на устройство): ядро для kind=quad идёт в body_kick, где k1 (горизонт) и k1v (вертикаль) НЕЗАВИСИМЫ.
UFC несёт ВЕРТИКАЛЬНЫЙ k1v (естественный фокус ондулятора k_y=½(B0/Bρ)²) и ГОРИЗОНТАЛЬНЫЙ k1, гасящий
побочную слабую фокусировку 1/ρ² секторных полюсов ⇒ нетто-оптика планарного ондулятора: горизонталь —
дрейф, вертикаль — k_y. Сила/фокус КОНКРЕТНОГО устройства не в .lte, а в undulator_knobs.csv (по зазору),
и применяются ПО ИНСТАНСУ при импорте — каждый ондулятор управляется отдельно, типы в .lte не плодятся.

Базовый .lte НЕ меняется: определения берутся из него ДОСЛОВНО, переписывается только линия ring.

ЗАПУСК
  py -3.14 undulators/assemble_undulators.py --all-on            # все std-* установлены
  py -3.14 undulators/assemble_undulators.py --only U05          # установлен один
  py -3.14 undulators/assemble_undulators.py                     # как в манифесте (колонка on)
"""
import os, sys, argparse, re, csv
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import lte_flatten as lf
import undulator_sr as sr
import matrix_tracker as mt
_QP = mt.QP                                          # q/p₀ ядра: mat_Gv = k1v/QP

HERE = os.path.dirname(__file__)
DEFAULT_LTE = os.path.join("workspace", "lte", "configuration_v5_sila2026", "sila2026.lte")


L_FC = 0.02                                          # длина фокус-элемента UFC, м (тонкая линза в промежутке)


def read_device(path):
    """(elems, line_tokens, L_dev, B0, sum_th2) устройства. sum_th2=Σθ² полюсов (для 1/ρ²-гашения)."""
    elems, lines, order = lf.read(path)
    line = lines[order[-1]]                           # единственная LINE устройства
    L_dev = sum(elems[t].get("L", 0.0) for t in line)
    sum_th2 = sum(elems[t].get("ANGLE", 0.0) ** 2 for t in line if elems[t]["type"] != "MARK")
    m = re.search(r"B0=([\d.]+)", open(path, encoding="utf-8").read())
    return elems, line, L_dev, (float(m.group(1)) if m else sr.P.B0_ref), sum_th2


def device_defs(delems, line):
    """Определения УНИКАЛЬНЫХ элементов-ПОЛЮСОВ (по одному разу, k=1) + фокус-элемент UFC.
    Маркеры устройства НЕ включаются: при monitors=True они стали бы BPM и разорвали бы прогон группы."""
    out, seen = [], set()
    for t in line:
        if t in seen or delems[t]["type"] == "MARK":
            continue
        seen.add(t)
        e = delems[t]
        out.append(f"{t}: CSBEND, L={e.get('L',0.0):.12g}, ANGLE={e.get('ANGLE',0.0):.15e}, N_SLICES=1")
    out.append(f"UFC: KQUAD, L={L_FC:.12g}, K1=0, N_SLICES=1")   # фокус-элемент (k1/k1v по инстансу при импорте)
    return out


def resolve_pos(pos, L_slot, L_dev, margin):
    """Смещение входа устройства a от начала промежутка. pos: start|center|end|доля|метры."""
    free = L_slot - L_dev
    pos = str(pos).strip().lower().replace(",", ".")
    if pos == "start":
        return margin
    if pos == "end":
        return free - margin
    if pos == "center":
        return free / 2.0
    v = float(pos)
    return v * free if 0.0 <= v <= 1.0 else v          # доля [0..1] или явные метры


def assemble(base_lte, device_lte, manifest, out_lte, knobs_csv,
             scenario_all=False, only=None, gap_override=None, margin=0.05, line_name="ring"):
    elems, lines, order = lf.read(base_lte)
    beamline = "RING" if "RING" in lines else order[-1]
    seq = list(lf.expand(beamline, lines))

    delems, dline, L_dev, b0_dev, sum_th2 = read_device(device_lte)
    Brho = sr.P.E / 0.299792458
    L_pole = next(delems[t].get("L", 0.0) for t in dline if delems[t]["type"] != "MARK")
    L_foot = L_dev + L_FC                              # footprint устройства (полюса + фокус-элемент)

    rows = list(csv.DictReader(open(manifest, encoding="utf-8"), delimiter=";"))
    for r in rows:                                     # нормализация активности (установлен ли слот)
        on = int(r["on"])
        if only is not None:
            on = 1 if r["slot_id"] in only else 0
        elif scenario_all:
            on = 1 if r["sp_pair"] != "inj-inj" else 0
        r["_on"] = on
    active = [r for r in rows if r["_on"]]
    active.sort(key=lambda r: int(r["i_first"]), reverse=True)   # врезка с конца — индексы не съезжают

    pole_seq = [t for t in dline if delems[t]["type"] != "MARK"]   # полюса без маркеров
    npole = len(pole_seq)
    drift_defs, knobs, report = [], [], []
    for r in active:
        slot = r["slot_id"]
        i0, i1 = int(r["i_first"]), int(r["i_last"])
        L_slot = float(r["L_slot"].replace(",", "."))
        gap = gap_override if gap_override is not None else float(r["gap_mm"].replace(",", "."))
        B0slot = sr.B0_of_gap(gap, sr.P)
        k = B0slot / b0_dev                            # сила относительно устройства (gap→поле)
        Ky = 0.5 * (B0slot / Brho) ** 2                # верт. фокус k_y (∫k1v ds = Ky·L_dev)
        Kx = k ** 2 * sum_th2 / L_pole                 # побочный горизонт. фокус секторных полюсов ∫1/ρ² ds
        # фокус-элемент UFC (kind=quad → body_kick): вертикаль — естеств. фокус, горизонталь — гашение 1/ρ²
        fc_mat_Gv = (-Ky * L_dev / L_FC) / _QP         # k1v<0 (фокус); ∫k1v ds = −Ky·L_dev
        fc_mat_G = (-Kx / L_FC) / _QP                  # ∫k1 ds = −Kx (гасит горизонт. фокус полюсов)
        if L_foot + 2 * margin > L_slot + 1e-9:
            raise SystemExit(f"[ОТКАЗ] слот {slot}: устройство {L_foot:.4f} м + зазор {2*margin:.4f} "
                             f"не влезает в промежуток {L_slot:.4f} м")
        a = resolve_pos(r["pos"], L_slot, L_foot, margin)
        b = L_slot - L_foot - a
        if a < margin - 1e-9 or b < margin - 1e-9:
            raise SystemExit(f"[ОТКАЗ] слот {slot}: pos={r['pos']} даёт a={a:.4f}, b={b:.4f} м "
                             f"(промежуток {L_slot:.4f}, footprint {L_foot:.4f}, зазор {margin:.4f})")
        da, db = f"{slot}_DA", f"{slot}_DB"
        drift_defs.append(f"{da}: DRIF, L={a:.15g}")
        drift_defs.append(f"{db}: DRIF, L={b:.15g}")
        seq[i0:i1 + 1] = [da] + pole_seq + ["UFC", db]   # полюса + фокус-элемент, свои дрейфы по краям
        knobs.append(f"{slot};{r['sp_pair']};cpmu18;{gap:.3f};{k:.6f};{npole};{r['pos']};{a:.6f};{b:.6f};"
                     f"{fc_mat_G:.9e};{fc_mat_Gv:.9e};{Ky:.6e};{Kx:.6e}")
        report.append((slot, r["sp_pair"], gap, k, a, b))

    base_txt = open(base_lte, encoding="utf-8").read()
    preamble = strip_line(base_txt, line_name)
    body = [",".join(seq[i:i + 12]) for i in range(0, len(seq), 12)]
    ring_txt = f"{line_name}: LINE=(" + ",&\n".join(body) + ")"
    header = (f"! сгенерировано assemble_undulators.py из {os.path.basename(base_lte)} + "
              f"{os.path.basename(device_lte)}\n"
              f"! ондуляторов установлено: {len(active)}; устройство {L_dev:.4f} м (общие полюса, k=1); "
              f"сила/фокус — по инстансу при импорте (undulator_knobs.csv)\n")
    os.makedirs(os.path.dirname(os.path.abspath(out_lte)), exist_ok=True)
    open(out_lte, "w", encoding="utf-8", newline="\n").write(
        header + preamble.rstrip() + "\n\n! ── полюса ондулятора (общие, k=1) ──\n" + "\n".join(device_defs(delems, dline))
        + "\n\n! ── дрейфы промежутков ──\n" + "\n".join(drift_defs) + "\n\n" + ring_txt + "\n")
    open(knobs_csv, "w", encoding="utf-8", newline="\n").write(
        "slot_id;sp_pair;device;gap_mm;k;n_poles;pos;a_m;b_m;fc_mat_G;fc_mat_Gv;Ky;Kx_poles\n"
        + "\n".join(knobs) + ("\n" if knobs else ""))
    return report, L_dev, len(seq)


def strip_line(txt, line_name):
    """Текст .lte без определения линии line_name (многострочное через &)."""
    lines = txt.splitlines()
    out, i = [], 0
    pat = re.compile(rf"^\s*{re.escape(line_name)}\s*:\s*LINE\s*=", re.I)
    while i < len(lines):
        if pat.match(lines[i]):
            while i < len(lines) and ")" not in lines[i]:
                i += 1
            i += 1
            continue
        out.append(lines[i]); i += 1
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description="вставить ондуляторы в решётку по манифесту")
    ap.add_argument("--lte", default=DEFAULT_LTE, help="базовый .lte кольца")
    ap.add_argument("--device", default=os.path.join(HERE, "undulator_cpmu18.lte"), help=".lte устройства")
    ap.add_argument("--manifest", default=os.path.join(HERE, "placements.csv"))
    ap.add_argument("--out-config", default="configuration_v6_sila2026_und", help="имя нового конфига")
    ap.add_argument("--line-name", default="ring")
    ap.add_argument("--all-on", action="store_true", help="установить все слоты кроме inj-inj")
    ap.add_argument("--only", nargs="*", default=None, help="установить только эти слоты")
    ap.add_argument("--gap", type=float, default=None, help="зазор всех устройств, мм (переопределяет манифест)")
    ap.add_argument("--margin", type=float, default=0.05, help="зазор безопасности с каждой стороны, м")
    a = ap.parse_args()

    out_dir = os.path.join("workspace", "lte", a.out_config)
    out_lte = os.path.join(out_dir, "sila2026_und.lte")
    knobs_csv = os.path.join(out_dir, "undulator_knobs.csv")
    report, L_dev, n = assemble(a.lte, a.device, a.manifest, out_lte, knobs_csv,
                                scenario_all=a.all_on, only=a.only, gap_override=a.gap,
                                margin=a.margin, line_name=a.line_name)
    print(f"  устройство L={L_dev:.4f} м; элементов в линии: {n}; установлено ондуляторов: {len(report)}")
    for slot, sp, gap, k, aa, bb in sorted(report)[:6]:
        print(f"  {slot}  {sp:8s}  gap={gap:5.2f} мм  k={k:.4f}  a={aa:.3f}  b={bb:.3f} м")
    if len(report) > 6:
        print(f"  … ещё {len(report)-6}")
    print(f"[OK] {out_lte}")
    print(f"[OK] {knobs_csv}")


if __name__ == "__main__":
    main()
