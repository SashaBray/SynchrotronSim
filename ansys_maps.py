# -*- coding: utf-8 -*-
"""Автоматическое снятие карт поля из Maxwell (AEDT): заказ из таблицы → готовые .fld.

Собирает по таблице заказа скрипт для штатного интерфейса AEDT и запускает его пакетно. PyAEDT не
нужен — используется тот же API, что стоит за «Record Script to File»; сигнатуры вызовов взяты из
`gRPCFiles/API/aedt/**/*.proto` установленной версии, а не из документации по памяти.

На каждую точку заказа: `SetVariableValue('current', ...)` → `Analyze('Setup1')` →
`FieldsReporter.ExportOnGrid(...)`. Ток вынесен в переменную дизайна во всех проектах, поэтому
геометрия не трогается вовсе.

ГЛАВНАЯ ЗАЩИТА: сетка НЕ задаётся руками, а читается из шапки уже снятой карты того же элемента
(`Grid Output Min: [...] Max: [...] Grid Size: [...]`), и после экспорта шапка новой карты
сверяется с эталонной. Карты сравнимы только если сняты одинаково; несовпадение — это остановка,
а не предупреждение.

    py -3.14 ansys_maps.py --plan                 # что будет сделано, без запуска
    py -3.14 ansys_maps.py --script               # только собрать скрипт для AEDT
    py -3.14 ansys_maps.py --run --element QP_sh  # собрать и запустить (СЧИТАЕТ, берёт лицензию)
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import subprocess
import sys

REQUEST = os.path.join("workspace", "tables", "map_request_quadrupoles.csv")
AEDT_EXE = r"C:\Program Files\ANSYS Inc\v252\AnsysEM\ansysedt.exe"
PROJECTS = r"C:\Sasha\Work\ansys_projects"
KD = r"C:\Sasha\Work\SILA\ansis approximations\KD"

# Проект, дизайн, сетап и имя токовой переменной. Проверено по .aedt: во всех проектах ток задан
# ДИЗАЙНОВОЙ переменной (без префикса $), через которую выражены все возбуждения.
DESIGN = {
    "QP_sh": dict(project="QP_sh", design="Maxwell3DDesign1", setup="Setup1", var="current"),
    "QP_L": dict(project="QP_L", design="Maxwell3DDesign1", setup="Setup1", var="current"),
}

HDR = re.compile(r"Grid Output Min:\s*\[(.*?)\]\s*Max:\s*\[(.*?)\]\s*Grid Size:\s*\[(.*?)\]")


def reference_grid(element):
    """Сетка и эталонная шапка — из уже снятой карты этого элемента. Руками не задаём."""
    pats = [os.path.join(KD, element, d, f"{element}_*.fld")
            for d in (f"{element.lower()}_steps", f"qp_{element.split('_')[-1]}_steps", "")]
    for p in pats:
        for f in sorted(glob.glob(p)):
            with open(f, encoding="utf-8", errors="ignore") as fh:
                head = fh.readline()
            m = HDR.search(head)
            if m:
                return dict(sample=f, header=head.rstrip("\n"),
                            lo=m.group(1).split(), hi=m.group(2).split(), step=m.group(3).split())
    raise SystemExit(f"{element}: не нашёл ни одной снятой карты, не с чем сверить сетку")


def jobs(element=None, path=REQUEST):
    if not os.path.exists(path):
        raise SystemExit(f"нет таблицы заказа {path}")
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8-sig"), delimiter=";"):
        if element and r["элемент"] != element:
            continue
        folder = os.path.join(KD, *r["папка"].split("/")[1:]).rstrip(os.sep)
        out.append(dict(element=r["элемент"], amp_turns=int(r["ампер_витков"]),
                        name=r["имя_файла"], folder=folder,
                        target=os.path.join(folder, r["имя_файла"]), why=r["назначение"]))
    return out


def build_script(element, todo, grid, out_path):
    """Скрипт для интерпретатора AEDT. Внутри доступен только oDesktop — библиотек не подключаем.

    Проект НЕ сохраняем и в конце возвращаем току прежнее значение: расчёт не должен менять
    состояние чужого проекта на диске."""
    d = DESIGN[element]
    lines = [
        "# -*- coding: utf-8 -*-",
        "# Сгенерировано ansys_maps.py. Запуск: ansysedt.exe -ng -runscriptandexit <этот файл>",
        "import os",
        "import ScriptEnv",
        'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")',
        # ГОТЧА: SetActiveProject не возвращает None для незагруженного проекта, а БРОСАЕТ
        # исключение. В пакетном режиме открыт вообще ничего не бывает, поэтому открываем сами.
        f'oProject = oDesktop.OpenProject(r"{os.path.join(PROJECTS, d["project"])}.aedt")',
        "if oProject is None:",
        f'    oProject = oDesktop.SetActiveProject("{d["project"]}")',
        f'oDesign = oProject.SetActiveDesign("{d["design"]}")',
        'oModule = oDesign.GetModule("FieldsReporter")',
        f'OLD = oDesign.GetVariableValue("{d["var"]}")',
        'oDesktop.AddMessage("", "", 0, "исходное значение тока: " + str(OLD))',
        "",
        f"LO = {grid['lo']!r}",
        f"HI = {grid['hi']!r}",
        f"STEP = {grid['step']!r}",
        "",
    ]
    for j in todo:
        val = f"{j['amp_turns']}A"
        lines += [
            f'oDesktop.AddMessage("", "", 0, "=== {j["name"]}: ток {val} ===")',
            f'oDesign.SetVariableValue("{d["var"]}", "{val}")',
            f'oDesign.Analyze("{d["setup"]}", True)',
            "oModule.CalcStack(\"clear\")",
            'oModule.EnterQty("B")',
            "oModule.ExportOnGrid(",
            f'    r"{j["target"]}",',
            "    LO, HI, STEP,",
            f'    "{d["setup"]} : LastAdaptive",',
            f'    ["{d["var"]}:=", "{val}"],',
            "    False, \"Cartesian\", [0, 0, 0], False)",
            "",
        ]
    lines += [
        f'oDesign.SetVariableValue("{d["var"]}", OLD)',          # вернуть как было
        'oDesktop.AddMessage("", "", 0, "=== все точки пройдены ===")',
    ]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    open(out_path, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return out_path


def probe(path, order=2, r=1e-3, nth=64):
    """Физическая проверка карты: максимум |B| и интегральная гармоника.

    Сверки шапки НЕ ДОСТАТОЧНО. Экспорт, сделанный до окончания расчёта или без пересчёта после
    смены тока, даёт файл с идеально правильной шапкой и полем от ПРЕЖНЕГО решения. Поймали такое
    на QP_sh_14500A: |B| в нём оказалось слабее, чем на карте 3294 А·в."""
    import numpy as np
    from field_quality import load, _circle
    h, B = load(path)
    mx = float(np.linalg.norm(B, axis=-1).max())
    xs = np.linspace(h["mnx"], h["mxx"], h["nx"])
    v = []
    for ix in range(h["nx"]):
        by, bz = _circle(h, B, ix, r, nth=nth)
        sp = np.fft.fft(by + 1j * bz) / nth
        v.append(sp[0] if order == 1 else sp[-(order - 1)])
    v = np.array(v)
    c = complex(2 * np.trapezoid(v.real, xs), 2 * np.trapezoid(v.imag, xs))
    return mx, abs(c) / r ** (order - 1)


def measure(element, order=2):
    """Все снятые карты элемента по возрастанию тока. Кривая насыщения обязана РАСТИ: любое
    падение при большем токе — испорченный экспорт, а не физика."""
    folder = os.path.join(KD, element, f"{element.lower()}_steps")
    if not os.path.isdir(folder):
        folder = os.path.join(KD, element, f"qp_{element.split('_')[-1]}_steps")
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, f"{element}_*A.fld"))):
        m = re.search(rf"{re.escape(element)}_(\d+)A\.fld$", os.path.basename(f))
        if not m:
            continue
        mx, hm = probe(f, order)
        rows.append((int(m.group(1)), os.path.basename(f), mx, hm))
    rows.sort()
    # Критерий — ГАРМОНИКА (то, что видит пучок), а не максимум |B|. Максимум сидит в железе, и при
    # насыщении локальный пик законно смещается и слегка проседает: у QP_sh с 9882 на 13176 А·в
    # max|B| падает 2.974 → 2.907, тогда как гармоника растёт 13.59 → 14.12. Это исправная карта.
    bad = [r for i, r in enumerate(rows) if i and r[3] < rows[i - 1][3] * 0.99]
    return rows, bad


def compare_maps(a_path, b_path):
    """Численное сравнение двух карт: шапка, форма массива и расхождение поля."""
    import numpy as np
    from field_quality import load
    ha, A = load(a_path)
    hb, B = load(b_path)
    head = [open(p, encoding="utf-8", errors="ignore").readline().rstrip("\n")
            for p in (a_path, b_path)]
    if A.shape != B.shape:
        return dict(header_same=head[0] == head[1], shape_same=False,
                    shape=(A.shape, B.shape))
    d = np.linalg.norm(A - B, axis=-1)
    scale = float(np.linalg.norm(A, axis=-1).max())
    return dict(header_same=head[0] == head[1], shape_same=True, shape=A.shape,
                max_abs=float(d.max()), rms=float(np.sqrt((d ** 2).mean())),
                max_rel=100.0 * float(d.max()) / scale if scale else float("nan"),
                identical=bool(np.array_equal(A, B)), scale=scale)


def cmd_repro(element, amp_turns, outdir):
    """Пересчитать УЖЕ СНЯТУЮ карту в отдельный файл и сравнить с ней. Оригинал не трогаем."""
    grid = reference_grid(element)
    folder = os.path.dirname(grid["sample"])
    orig = os.path.join(folder, f"{element}_{amp_turns}A.fld")
    if not os.path.exists(orig):
        raise SystemExit(f"нет исходной карты {orig}")
    os.makedirs(outdir, exist_ok=True)
    repro = os.path.join(os.path.abspath(outdir), f"{element}_{amp_turns}A__repro.fld")
    job = [dict(element=element, amp_turns=amp_turns, name=os.path.basename(repro),
                folder=outdir, target=repro, why="воспроизведение")]
    sp = build_script(element, job, grid, os.path.join(outdir, f"aedt_repro_{element}.py"))
    lock = os.path.join(PROJECTS, f"{DESIGN[element]['project']}.aedt.lock")
    if os.path.exists(lock):
        raise SystemExit(f"проект открыт в AEDT ({lock}) — закройте его")
    print(f"  оригинал: {orig}")
    print(f"  пересчёт: {repro}")
    print(f"  скрипт:   {sp}\n  запуск AEDT…")
    rc = subprocess.call([AEDT_EXE, "-ng", "-runscriptandexit", os.path.abspath(sp)])
    print(f"  AEDT вернул {rc}")
    if not os.path.exists(repro):
        raise SystemExit("файл не появился — аргументы ExportOnGrid неверны либо расчёт не прошёл")
    r = compare_maps(orig, repro)
    print(f"\n  шапка совпала: {'да' if r['header_same'] else 'НЕТ'}")
    if not r["shape_same"]:
        print(f"  РАЗНАЯ ФОРМА МАССИВА: {r['shape']}")
        return r
    print(f"  форма массива: {r['shape']}, максимум |B| {r['scale']:.4f} Тл")
    print(f"  бит-в-бит: {'ДА' if r['identical'] else 'нет'}")
    print(f"  расхождение поля: макс {r['max_abs']:.3e} Тл ({r['max_rel']:.4f}% от максимума), "
          f"СКО {r['rms']:.3e} Тл")
    return r


def verify(element, todo, grid):
    """Сверка шапок: новая карта обязана быть снята на ТОЙ ЖЕ сетке, что эталонная."""
    ok, bad, missing = [], [], []
    for j in todo:
        if not os.path.exists(j["target"]):
            missing.append(j["name"]); continue
        with open(j["target"], encoding="utf-8", errors="ignore") as fh:
            head = fh.readline().rstrip("\n")
        (ok if head.strip() == grid["header"].strip() else bad).append((j["name"], head))
    return ok, bad, missing


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", help="только этот элемент (иначе все из таблицы)")
    ap.add_argument("--plan", action="store_true", help="показать план, ничего не делать")
    ap.add_argument("--script", action="store_true", help="собрать скрипт для AEDT, не запускать")
    ap.add_argument("--run", action="store_true", help="собрать И ЗАПУСТИТЬ (считает, берёт лицензию)")
    ap.add_argument("--verify", action="store_true", help="сверить шапки уже снятых карт")
    ap.add_argument("--measure", action="store_true",
                    help="физическая проверка снятых карт: |B| обязано расти с током")
    ap.add_argument("--order", type=int, default=2, help="порядок гармоники для проверки")
    ap.add_argument("--limit", type=int, help="взять только первые N точек из оставшихся: "
                                             "первый прогон разумно делать на ОДНОЙ")
    ap.add_argument("--repro", type=int, metavar="А·ВИТКОВ",
                    help="пересчитать УЖЕ СНЯТУЮ карту и сравнить с ней (оригинал не трогается)")
    ap.add_argument("--outdir", default=os.path.join("workspace", "results", "aedt"))
    a = ap.parse_args(argv)

    els = [a.element] if a.element else sorted({j["element"] for j in jobs()})
    if a.repro:
        if len(els) != 1:
            raise SystemExit("--repro требует --element")
        print(f"=== {els[0]}: воспроизведение карты {a.repro} А·в")
        cmd_repro(els[0], a.repro, a.outdir)
        return 0
    for el in els:
        if el not in DESIGN:
            print(f"{el}: проект не описан в DESIGN — пропускаю"); continue
        grid = reference_grid(el)
        todo = jobs(el)
        done = [j for j in todo if os.path.exists(j["target"])]
        left = [j for j in todo if not os.path.exists(j["target"])]
        print(f"\n=== {el}: заказано {len(todo)}, уже снято {len(done)}, осталось {len(left)}")
        print(f"  сетка из эталона {os.path.basename(grid['sample'])}: "
              f"min {' '.join(grid['lo'])}  max {' '.join(grid['hi'])}  шаг {' '.join(grid['step'])}")
        if a.measure:
            rows, bad = measure(el, a.order)
            print(f"  {'А·витков':>9} {'файл':24} {'макс |B|, Тл':>13} {'гармоника':>12}")
            for ni, nm, mx, hm in rows:
                flag = "  ← ПАДЕНИЕ ПРИ БОЛЬШЕМ ТОКЕ" if (ni, nm, mx, hm) in bad else ""
                print(f"  {ni:9d} {nm:24} {mx:13.4f} {hm:12.5g}{flag}")
            print(f"  испорченных карт: {len(bad)}")
            continue
        if a.verify:
            ok, bad, missing = verify(el, todo, grid)
            print(f"  сверка шапок: совпало {len(ok)}, РАСХОЖДЕНИЙ {len(bad)}, нет файла {len(missing)}")
            for n, h in bad:
                print(f"    ! {n}\n      {h}")
            continue
        if a.limit:
            left = left[:a.limit]
        for j in left:
            print(f"  {j['amp_turns']:7d} А·в → {j['name']:24} {j['why']}")
        if a.plan or not left:
            continue
        sp = build_script(el, left, grid, os.path.join(a.outdir, f"aedt_job_{el}.py"))
        print(f"  скрипт: {sp}")
        if not a.run:
            continue
        if not os.path.exists(AEDT_EXE):
            raise SystemExit(f"не найден {AEDT_EXE}")
        lock = os.path.join(PROJECTS, f"{DESIGN[el]['project']}.aedt.lock")
        if os.path.exists(lock):
            raise SystemExit(f"проект открыт в AEDT ({lock}) — закройте его перед пакетным прогоном")
        print(f"  запуск AEDT ({len(left)} точек, это надолго)…")
        rc = subprocess.call([AEDT_EXE, "-ng", "-runscriptandexit", os.path.abspath(sp)])
        print(f"  AEDT вернул {rc}")
        ok, bad, missing = verify(el, todo, grid)
        print(f"  итог: снято {len(ok)}, расхождений шапки {len(bad)}, не появилось {len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
