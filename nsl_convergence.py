# -*- coding: utf-8 -*-
"""Сходимость срезов по типам устройств: сколько срезов нужно каждому типу для заданной точности тюнов и
карты оборота — и сколько это даёт скорости (synchrotronsim.convergence; шаг 6 docs/LIBRARY_ARCHITECTURE.md).

    py -3.14 nsl_convergence.py --lattice 26
    py -3.14 nsl_convergence.py --lattice 26 --tol-nu 1e-3 --tol-map-um 1 --tol-nu 1e-4 --tol-map-um 0.1

Выход: workspace/results/nsl_convergence_<lat>.csv (для человека: ошибки по типам и кандидатам, десятичная
запятая) и .json (рекомендации nsl_dev для каждого допуска + проверка набора целиком и время оборота).
Ничего в конфигурацию кокпита не пишет — рекомендации применяет пользователь.
"""
import argparse
import csv
import json
import os
import sys

os.environ.setdefault("SS_SECTOR", "1")
sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, ROOT)

from synchrotronsim import convergence as cv  # noqa: E402


def _h(v):
    return "" if v is None else (f"{v:.4g}".replace(".", ",") if v == v and v != float("inf") else "потеря")


def main():
    ap = argparse.ArgumentParser(description="сходимость срезов по типам устройств")
    ap.add_argument("--lattice", type=int, default=26)
    ap.add_argument("--ref", type=int, default=cv.REF_NSL, help="эталонное число срезов (паспорт: 40)")
    ap.add_argument("--tol-nu", type=float, action="append", help="допуск суммы |Δν| (можно несколько)")
    ap.add_argument("--tol-map-um", type=float, action="append", help="допуск суммы |Δкарты| за оборот, мкм")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--reuse", action="store_true", help="не сканировать заново — взять сохранённый скан")
    ap.add_argument("--save", choices=("tight", "loose"), default=None,
                    help="записать рекомендацию в workspace/configs/slices_<lat>.json — набор кокпита по умолчанию")
    a = ap.parse_args()
    tols = list(zip(a.tol_nu or [1e-3, 1e-4], a.tol_map_um or [1.0, 0.1]))

    print(f"=== сходимость срезов: решётка {a.lattice}, эталон {a.ref} срезов ===", flush=True)
    out = os.path.join(ROOT, "workspace", "results")
    os.makedirs(out, exist_ok=True)
    base = os.path.join(out, f"nsl_convergence_{a.lattice}")
    raw = base + "_scan.json"
    if a.reuse and os.path.exists(raw):                       # подбор заново по сохранённому скану
        with open(raw, encoding="utf-8") as f:
            S = json.load(f)
        S["errors"] = {d: {int(n): tuple(v) for n, v in e.items()} for d, e in S["errors"].items()}
        print(f"скан из {os.path.relpath(raw)}", flush=True)
    else:
        S = cv.scan(a.lattice, ref=a.ref, workers=a.workers)
        with open(raw, "w", encoding="utf-8") as f:
            json.dump({k: v for k, v in S.items() if k != "indices"}, f, ensure_ascii=False,
                      default=lambda x: None if x != x or x in (float("inf"),) else x)
    print(f"скан: {len(S['errors'])} типов × {len(S['candidates'])} кандидатов за {S['seconds']:.0f} с; "
          f"эталонные тюны {S['ref_tunes'][0]:.6f} / {S['ref_tunes'][1]:.6f}\n", flush=True)

    recs = []
    for tn, tm in tols:
        rec = cv.recommend(S, tol_nu=tn, tol_map=tm * 1e-6)
        v = cv.validate(a.lattice, rec, S)
        cost = sum(S["count"][d] * n for d, n in rec.items())
        recs.append({"tol_nu": tn, "tol_map_um": tm, "nsl_dev": rec, "validate": v, "cost": cost})
        print(f"допуск |Δν| ≤ {tn:g}, |Δкарты| ≤ {tm:g} мкм: срезы × элементы = {cost}; "
              f"факт: Δνx {v['dnux']:.1e}, Δνy {v['dnuy']:.1e}, Δкарты {v['dmap'] * 1e6:.3f} мкм; "
              f"оборот 1 частицы {v['ms_turn']:.1f} мс (эталон {v['ms_turn_ref']:.1f}), "
              f"с СИ+ВЧ {v.get('ms_turn_sr', float('nan')):.1f} мс, 12 частиц {v.get('ms_turn_sr_12', float('nan')):.1f} мс",
              flush=True)

    # нынешние наборы для сравнения: кокпит (nsl 8, полюса ИД 1 срез) и nsl 8 везде
    from control.ring_model import ID_POLE_DEVS
    cmp = {}
    for name, nsl, dev in (("кокпит (nsl 8, полюса ИД 1)", 8, {d: 1 for d in ID_POLE_DEVS}), ("nsl 8 везде", 8, {}),
                           ("nsl 5 везде (паритет elegant)", 5, {})):
        v = cv.validate(a.lattice, dev, S, nsl=nsl)
        cmp[name] = v
        print(f"{name}: Δνx {v['dnux']:.1e}, Δνy {v['dnuy']:.1e}, Δкарты {v['dmap'] * 1e6:.3f} мкм; "
              f"оборот {v['ms_turn']:.1f} мс, с СИ+ВЧ {v.get('ms_turn_sr', float('nan')):.1f} мс, "
              f"12 частиц {v.get('ms_turn_sr_12', float('nan')):.1f} мс", flush=True)

    with open(base + ".csv", "w", encoding="utf-8-sig", newline="") as f:
        w = csv.writer(f, delimiter=";")
        w.writerow(["тип устройства", "вид", "элементов"]
                   + [f"|Δν| при {n}" for n in S["candidates"]] + [f"|Δкарты| мкм при {n}" for n in S["candidates"]]
                   + [f"рекомендация {r['tol_nu']:g}/{r['tol_map_um']:g} мкм" for r in recs])
        order = sorted(S["errors"], key=lambda d: -S["count"][d])
        for d in order:
            e = S["errors"][d]
            w.writerow([d, S["kind"][d], S["count"][d]]
                       + [_h(max(e[n][0], e[n][1])) for n in S["candidates"]]
                       + [_h(e[n][2] * 1e6) for n in S["candidates"]]
                       + [r["nsl_dev"][d] for r in recs])
    with open(base + ".json", "w", encoding="utf-8") as f:
        json.dump({"lattice": a.lattice, "ref": a.ref, "ref_tunes": S["ref_tunes"], "recommendations": recs,
                   "compare": cmp}, f, ensure_ascii=False, indent=1)
    print(f"\nтаблица: {os.path.relpath(base + '.csv')}; рекомендации: {os.path.relpath(base + '.json')}")
    if a.save:
        print(f"набор кокпита по умолчанию ({a.save}): {os.path.relpath(cv.save_config(a.lattice, a.save))}")


if __name__ == "__main__":
    main()
