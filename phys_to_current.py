# -*- coding: utf-8 -*-
"""Конвертер «физкоэффициенты → токи»: обратная матрица связи для запуска настроенной решётки в
режимах ТОКОВ (режим 1 — диффурное ядро, режим 3 — киковое ядро на токах).

Решётка настраивается в режиме ФИЗКОЭФ. (mat_B/mat_G/mat_S/mat_O в lattice_configs — режим 2). Чтобы
ТА ЖЕ конфигурация пошла в режимах на токах, для каждой группы со-локованных обмоток решаем обратную
задачу: найти входы (токи обмоток / коэф. диполя), которые через матрицу связи дают тот же набор
физкоэф. Для группы обмоток W и активных мультиполей M:  A·inp = target,  где
  A[m][w] = coupling[w][per_<m>]  (вклад обмотки w в мультиполь m на 1 ед. входа),
  target[m] = Σ_w mat_<m>[w]      (суммарный физкоэф. группы, режим 2),
решение — наименьшей нормы (np.linalg.lstsq): QP 1×1 точно, SP 1×2 делит поровну (per_S равны),
DQ 3×3 точно (B,G,S). Диполи (input_type=coeff): correct_coeff = mat_B/per_B; тот же коэф. верен и
для диффурного ядра (поле карты × correct_coeff даёт нужный изгиб). Остальные (current): arg_val = ток.

Идемпотентно. Бэкап lattice_configs сохраняется. Запуск: python phys_to_current.py [lattice=101]"""
import os
import sys
import csv
import shutil
import numpy as np

TB = os.path.join("workspace", "tables")
MULTS = [("per_B_Tl", "mat_B"), ("per_G_Tlm", "mat_G"),
         ("per_S_Tlm2", "mat_S"), ("per_O_Tlm3", "mat_O")]


def _read(name):
    return list(csv.DictReader(open(os.path.join(TB, name), encoding="utf-8"), delimiter=";"))


def convert(lattice=101, verbose=True, dipole_priority=False):
    """dipole_priority: ПРАВИЛО для вырожденной матрицы связи. Если обмотки не позволяют
    удовлетворить ВСЕ физкоэф. (combined-function DQ: обмотки I2≡I3, ранг 2), приоритет отдаётся
    ДИПОЛЬНОЙ составляющей — ток подбирается так, чтобы РЕАЛЬНОЕ поле (секторная модель, не линейная
    матрица связи — её per_B у DQ завышает изгиб ×4–8) дало проектный изгиб ∫b0=−mat_B·L_eff;
    градиент/секступоль становятся тем, что выходит (I2=I3=0). Цель — центроид-траектория похожа на
    проектную (изгиб — главное для орбиты). Для обычных (невырожденных) элементов — как раньше."""
    coup = {r["device_type_id"]: r for r in _read("coupling_matrix.csv")}
    dbs = None
    if dipole_priority:
        import dq_bend_solve as dbs
        _hw = {r["device_type_id"]: r for r in _read("hardware_library.csv")}
        _eff = {d: float(r["eff_length"]) for d, r in _hw.items()}
        _dqcache = {}
    path = os.path.join(TB, "lattice_configs.csv")
    allrows = list(csv.DictReader(open(path, encoding="utf-8"), delimiter=";"))
    fields = list(allrows[0].keys())
    ring = [r for r in allrows if r["lattice_id"] == str(lattice)]
    ring.sort(key=lambda r: int(r["instance_id"]))

    # группы со-локованных обмоток (как build_line трекера: совпадение x,y,z)
    groups = []
    i = 0
    while i < len(ring):
        g = [ring[i]]
        while i + 1 < len(ring) and (ring[i + 1]["x"], ring[i + 1]["y"], ring[i + 1]["z"]) == \
                (ring[i]["x"], ring[i]["y"], ring[i]["z"]):
            g.append(ring[i + 1]); i += 1
        groups.append(g); i += 1

    stats = {"current": 0, "coeff": 0, "skipped": 0, "maxI": 0.0, "maxC": 0.0, "unreached": [], "dep": {}}
    for g in groups:
        cpls = [coup.get(w["device_type_id"]) for w in g]
        if any(c is None for c in cpls):
            stats["skipped"] += len(g); continue
        itype = cpls[0].get("input_type", "current")
        # ПРАВИЛО ПРИОРИТЕТА ДИПОЛЯ: combined-function DQ (3 обмотки, вырожд. связь) — ток по РЕАЛЬНОМУ
        # изгибу (сектор), а не по завышающей матрице связи. I1 → проектный изгиб, I2=I3=0.
        if dbs is not None and len(g) == 3 and g[0]["device_type_id"].startswith("DQ_I"):
            devs = tuple(w["device_type_id"] for w in g)
            Leff = _eff.get(devs[0], 1.0); matB = float(g[0].get("mat_B") or 0)
            key = (devs, round(matB, 10))
            if key not in _dqcache:
                _dqcache[key] = dbs.solve_i1(devs, -matB * Leff, _hw, float(g[0].get("arg_val") or 0))
            I1, _res = _dqcache[key]
            g[0]["arg_val"] = f"{I1:.9g}"; g[1]["arg_val"] = "0"; g[2]["arg_val"] = "0"
            stats["current"] += 3; stats["maxI"] = max(stats["maxI"], abs(I1))
            stats["dep"]["mat_G(DQ→0)"] = stats["dep"].get("mat_G(DQ→0)", 0) + 1
            continue
        target = {mc: sum(float(w.get(mc) or 0) for w in g) for (_, mc) in MULTS}
        active = [(pc, mc) for (pc, mc) in MULTS
                  if any(abs(float(c.get(pc) or 0)) > 1e-15 for c in cpls)]
        if not active:
            inp = np.zeros(len(g))
        else:
            A = np.array([[float(c.get(pc) or 0) for c in cpls] for (pc, mc) in active])
            # РАНГ матрицы связи = сколько мультиполей реально независимо управляемы обмотками.
            # У DQ обмотки I2≡I3 (вырождение) ⇒ ранг 2: можно задать только 2 мультиполя (B,G),
            # а S (низший приоритет) — ЗАВИСИМЫЙ. Оставляем первые `rank` (порядок MULTS: B,G,S,O),
            # лишние убираем из цели, иначе солвер гасит их огромными токами (было 242А).
            sv = np.linalg.svd(A, compute_uv=False)
            rank = int((sv > 1e-2 * sv[0]).sum()) if len(sv) else 0
            dropped = []
            if 0 < rank < len(active):
                dropped = active[rank:]; active = active[:rank]
                A = np.array([[float(c.get(pc) or 0) for c in cpls] for (pc, mc) in active])
            b = np.array([target[mc] for (pc, mc) in active])
            inp, *_ = np.linalg.lstsq(A, b, rcond=None)
            # ЗАВИСИМЫЕ мультиполи (S у DQ): вычислить достигнутое значение и записать в mat_* первой
            # обмотки — тогда режимы 2 (mat_*) и 3 (ток×coupling) согласованы на реальном (не нулевом) S.
            for (pc, mc) in dropped:
                got = float(sum(float(cpls[j].get(pc) or 0) * inp[j] for j in range(len(g))))
                g[0][mc] = f"{got:.9g}"
                if mc not in stats["dep"]:
                    stats["dep"][mc] = 0
                stats["dep"][mc] += 1
        # запись входов в обмотки
        for j, w in enumerate(g):
            if itype == "coeff":
                w["correct_coeff"] = f"{inp[j]:.9g}"; w["arg_val"] = "-"
                stats["coeff"] += 1; stats["maxC"] = max(stats["maxC"], abs(inp[j]))
            else:
                w["arg_val"] = f"{inp[j]:.9g}"
                stats["current"] += 1; stats["maxI"] = max(stats["maxI"], abs(inp[j]))

    shutil.copy(path, path + ".bak")
    # обновить строки кольца в общем файле (ring — те же объекты dict, что и в allrows)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter=";"); w.writeheader(); w.writerows(allrows)

    if verbose:
        print(f"Решётка {lattice}: групп {len(groups)}. Заполнено токов (arg_val): {stats['current']}, "
              f"коэф. диполей (correct_coeff): {stats['coeff']}, пропущено: {stats['skipped']}.")
        print(f"  макс |ток|={stats['maxI']:.2f} А, макс |коэф.диполя|={stats['maxC']:.4g} "
              f"(бэкап: lattice_configs.csv.bak)")
        if stats["dep"]:
            print(f"  ЗАВИСИМЫЕ мультиполи (вырожд. связь ⇒ записаны в mat_* реальные значения): {stats['dep']}")
        if stats["unreached"]:
            print(f"  ВНИМАНИЕ: {len(stats['unreached'])} мультиполей без опоры в матрице связи "
                  f"(не воспроизводимы током): {stats['unreached'][:5]}")
    return stats


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("-")]
    lat = int(args[0]) if args else 101
    dp = "--dipole-priority" in sys.argv        # опционально: приоритет диполя у DQ (иначе min-norm)
    convert(lat, dipole_priority=dp)
