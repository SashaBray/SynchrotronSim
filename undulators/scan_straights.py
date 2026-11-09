# -*- coding: utf-8 -*-
"""Сканер ПУСТЫХ ПРЯМОЛИНЕЙНЫХ ПРОМЕЖУТКОВ в кольце — предлагает манифест размещения ондуляторов.

ЗАЧЕМ. Промежуток между супериодами — это НЕПРЕРЫВНЫЙ прогон дрейфов (между двумя стандартными SP это
DR_01,DR_01 ≈ 5.48 м; между стандартным и инжекционным — DR_01,DR_45). Сканер идёт по развёрнутому
кольцу, считает продольную координату s (= накопленная сумма L), находит максимальные прогоны подряд
идущих DRIF длиннее порога и КЛАССИФИЦИРУЕТ каждый по соседним магнитам:
  std-std  — оба соседних магнита стандартные,
  std-inj  — один стандартный, один инжекционный/junction (имена семейства J или Q[FD]\\dL),
  inj-inj  — оба инжекционные (ВНУТРИ инжекционного SP — туда ондулятор НЕ ставится).

Результат — placements.csv (предложение; человек правит, какие слоты активны и как позиционированы).
Колонки: slot_id;sp_pair;s_start;s_end;L_slot;i_first;i_last;left;right;device;gap_mm;pos;on
  i_first/i_last — индексы элементов в развёрнутой линии (НАДЁЖНЫЙ ключ врезки, имя DR_01 неуникально);
  pos — положение устройства в промежутке: start|center|end | доля 0..1 | смещение входа в метрах;
  on  — 1 активен (устройство вставляется), 0 пропуск.

ЗАПУСК
  py -3.14 undulators/scan_straights.py                                  # базовый конфиг СИЛА-26
  py -3.14 undulators/scan_straights.py --lte workspace/lte/.../x.lte --min-len 3.0
"""
import os, sys, argparse, re
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import lte_flatten as lf

MAGNET = {"CSBEND", "KQUAD", "KSEXT", "KOCT"}        # «настоящие» магниты для классификации
DEFAULT_LTE = os.path.join("workspace", "lte", "configuration_v5_sila2026", "sila2026.lte")


def is_injection(name):
    """Имя магнита из инжекционного/junction-семейства СИЛА-26."""
    return ("J" in name) or bool(re.match(r"Q[FD]\dL$", name))


def scan(lte_path, beamline="RING", min_len=3.0):
    elems, lines, order = lf.read(lte_path)
    if beamline not in lines:
        beamline = order[-1]
    seq = lf.expand(beamline, lines)

    # s каждого элемента + его тип/длина
    info = []
    s = 0.0
    for idx, nm in enumerate(seq):
        e = elems.get(nm, {})
        L = e.get("L", 0.0)
        info.append(dict(idx=idx, name=nm, type=e.get("type", "?"), L=L, s0=s, s1=s + L))
        s += L
    perimeter = s

    def neighbor_magnet(i, step):
        """Ближайший настоящий магнит в направлении step (±1), пропуская дрейфы/маркеры."""
        j = i + step
        while 0 <= j < len(info):
            if info[j]["type"] in MAGNET:
                return info[j]["name"]
            j += step
        return "?"

    # максимальные прогоны подряд идущих DRIF
    runs = []
    i = 0
    while i < len(info):
        if info[i]["type"] == "DRIF":
            j = i
            while j < len(info) and info[j]["type"] == "DRIF":
                j += 1
            L_slot = sum(info[k]["L"] for k in range(i, j))
            if L_slot >= min_len:
                left = neighbor_magnet(i, -1)
                right = neighbor_magnet(j - 1, +1)
                li, ri = is_injection(left), is_injection(right)
                sp_pair = "inj-inj" if (li and ri) else ("std-inj" if (li or ri) else "std-std")
                runs.append(dict(s0=info[i]["s0"], s1=info[j - 1]["s1"], L=L_slot,
                                 i0=i, i1=j - 1, left=left, right=right, sp_pair=sp_pair))
            i = j
        else:
            i += 1
    return runs, perimeter


def write_manifest(runs, path, device="cpmu18"):
    rows = ["slot_id;sp_pair;s_start;s_end;L_slot;i_first;i_last;left;right;device;gap_mm;pos;on"]
    for k, r in enumerate(runs, 1):
        # inj-inj по правилу размещения выключен по умолчанию; остальные — пользователь включит
        on = 0
        rows.append(f"U{k:02d};{r['sp_pair']};{r['s0']:.6f};{r['s1']:.6f};{r['L']:.6f};"
                    f"{r['i0']};{r['i1']};{r['left']};{r['right']};{device};20.0;center;{on}")
    open(path, "w", encoding="utf-8", newline="\n").write("\n".join(rows) + "\n")


def main():
    ap = argparse.ArgumentParser(description="найти пустые прямолинейки и предложить манифест ондуляторов")
    ap.add_argument("--lte", default=DEFAULT_LTE, help="базовый .lte кольца")
    ap.add_argument("--beamline", default="RING", help="имя линии")
    ap.add_argument("--min-len", type=float, default=3.0, help="минимальная длина прогона дрейфов, м")
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "placements.csv"))
    a = ap.parse_args()

    runs, perimeter = scan(a.lte, a.beamline, a.min_len)
    by = {}
    for r in runs:
        by[r["sp_pair"]] = by.get(r["sp_pair"], 0) + 1
    print(f"  периметр линии, м: {perimeter:.3f}")
    print(f"  промежутков ≥{a.min_len} м: {len(runs)}  " + ", ".join(f"{k}={v}" for k, v in sorted(by.items())))
    for k, r in enumerate(runs, 1):
        print(f"  U{k:02d}  {r['sp_pair']:8s}  s={r['s0']:8.3f}…{r['s1']:8.3f}  L={r['L']:.3f} м  "
              f"[{r['left']} | {r['right']}]")
    write_manifest(runs, a.out)
    print(f"[OK] {a.out}")


if __name__ == "__main__":
    main()
