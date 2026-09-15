# -*- coding: utf-8 -*-
"""Автоматическое снятие карт поля из Maxwell (AEDT): заказ из таблицы → готовые .fld.

Собирает по таблице заказа скрипт для штатного интерфейса AEDT и запускает его пакетно. PyAEDT не
нужен — используется тот же API, что стоит за «Record Script to File»; сигнатуры вызовов взяты из
`gRPCFiles/API/aedt/**/*.proto` установленной версии, а не из документации по памяти.

На каждую точку заказа: выставить ток КАЖДОЙ обмотки (незадействованные — в ноль), `Analyze`,
`FieldsReporter.ExportOnGrid`. Токи во всех проектах вынесены в переменные дизайна, поэтому
геометрия не трогается вовсе.

ТРИ ЗАЩИТЫ, и каждая понадобилась на практике:
  1. Сетка не задаётся руками, а читается из шапки уже снятой карты того же элемента, и шапка
     новой сверяется с эталонной.
  2. `guard_material` сверяет САМУ КРИВУЮ намагничивания в проекте, а не имя материала: имя лжёт
     (в QP_sh материал звался `steel_1008 _custom`, а содержал низконасыщающуюся кривую).
  3. `--measure`: гармоника обязана расти с током. Экспорт до окончания расчёта даёт файл с
     идеально правильной шапкой и полем от ПРЕЖНЕГО решения — формальные сверки такое не ловят.

    py -3.14 ansys_maps.py --plan                  # что будет сделано, без запуска
    py -3.14 ansys_maps.py --element OP --run      # собрать и запустить (СЧИТАЕТ, берёт лицензию)
    py -3.14 ansys_maps.py --element OP --measure  # физическая проверка снятого
    py -3.14 ansys_maps.py --element OP --register # завести в реестр карт
"""
from __future__ import annotations

import argparse
import csv
import filecmp
import glob
import os
import re
import subprocess
import sys

REQUEST = os.path.join("workspace", "tables", "map_request.csv")
AEDT_EXE = r"C:\Program Files\ANSYS Inc\v252\AnsysEM\ansysedt.exe"
PROJECTS = r"C:\Sasha\Work\ansys_projects"
# Папка, куда AEDT пишет экспорт и откуда берётся эталонная сетка. Прежняя (рядом с КД) убрана с
# диска; если её нет, экспорт идёт в папку проекта. Переопределяется переменной SS_AEDT_EXPORT.
_KD_OLD = r"C:\Sasha\Work\SILA\ansis approximations\KD"
KD = (os.environ.get("SS_AEDT_EXPORT")
      or (_KD_OLD if os.path.isdir(_KD_OLD)
          else os.path.join(os.path.dirname(os.path.abspath(__file__)), "workspace", "aedt_exports")))
# Именованное выражение калькулятора, которым сняты ВСЕ прежние карты — видно в их шапке:
# `X, Y, Z, Vector data "B_Vector"`. Сырое EnterQty("B") даёт величину с именем "<Bx,By,Bz>",
# и карта становится несопоставимой с серией. Менять нельзя.
QUANTITY = "B_Vector"
# Операция калькулятора перед экспортом. Пользователь выгружал карты СГЛАЖЕННЫМИ, поэтому делаем
# так же. ИЗМЕРЕНО: на QP_sh это ПУСТАЯ операция — результат со Smooth и без совпал до 3.5e-12 Тл,
# то есть ExportOnGrid и так отдаёт узловое поле. Оставлено ради точного совпадения с ручным
# способом. Пустая строка отключает операцию.
SMOOTH = "Smooth"

# Проект, дизайн, сетап и соответствие «обмотка → переменная дизайна». Проверено по .aedt: токи
# везде заданы ДИЗАЙНОВЫМИ переменными (без префикса $), через которые выражены все возбуждения.
# `bare` — стиль имени файла: у квадруполей `QP_sh_1647A.fld` без токена обмотки, у остальных
# `OP_I1_100A.fld` и `SP_I1_1891A_I2_800A.fld`.
DESIGN = {
    "QP_sh": dict(project="QP_sh", vars={"I": "current"}, folder="qp_sh_steps", bare=True,
                  turns={"I": 68}),
    "QP_L": dict(project="QP_L", vars={"I": "current"}, folder="qp_L_steps", bare=True,
                 turns={"I": 85}),
    "OP": dict(project="OP", vars={"I1": "I"}, folder="op_step", bare=False,
               turns={"I1": 34}),
    "SP": dict(project="SP", vars={"I1": "I1", "I2": "I2"}, folder="SP_steps", bare=False,
               turns={"I1": 41, "I2": 390}),
    "DQ": dict(project="DQ", vars={"I1": "I1", "I2": "I2", "I3": "I_cor"}, folder="DQ_steps",
               bare=False, turns={"I1": 53, "I2": 13, "I3": 90}),
    "corrector": dict(project="corrector", folder="corrector_steps", bare=False,
                      vars={"main": "I_main", "c1": "I_cor_1", "c2": "I_cor_2",
                            "c3": "I_cor_3", "c4": "I_cor_4"},
                      turns={"main": 622, "c1": 548, "c2": 548, "c3": 548, "c4": 548}),
}
for _d in DESIGN.values():
    _d.setdefault("design", "Maxwell3DDesign1")
    _d.setdefault("setup", "Setup1")

ORDER_BY_ELEMENT = {"QP_sh": 2, "QP_L": 2, "SP": 3, "OP": 4, "DQ": 2, "corrector": 1}


def probe_radius(element):
    """Радиус пробы по типу элемента — тот же, что у saturation_model.

    Брать 1 мм для всех НЕЛЬЗЯ: поле мультиполя ∝ r^(n−1), и у октуполя на 1 мм сигнал лежит ПОД
    шумовым полом карт. Сравнение двух таких карт даёт бессмыслицу — у OP получалось −12.7% и
    +71.4% там, где на правильных 14 мм расхождение всего −1.1% и +0.2%."""
    import saturation_model as sm
    return sm.PROBE_R.get(sm.KIND.get(element, "quad"), 1e-3)

HDR = re.compile(r"Grid Output Min:\s*\[(.*?)\]\s*Max:\s*\[(.*?)\]\s*Grid Size:\s*\[(.*?)\]")
BH_BLOCK = re.compile(r"Points\[\d+:\s*([0-9eE.,+\-\s]+)\]")


# ───────────────────────────── имена и папки ─────────────────────────────

def folder_of(element):
    return os.path.join(KD, element, DESIGN[element]["folder"])


def order_of(element):
    """Обмотки в порядке, в котором они идут в имени файла и в скрипте."""
    return list(DESIGN[element]["vars"])


def map_name(element, currents):
    """Имя файла по набору токов: `SP_I1_1891A_I2_800A.fld`, у квадруполей `QP_sh_1647A.fld`."""
    used = [(w, currents[w]) for w in order_of(element) if currents.get(w)]
    if DESIGN[element]["bare"]:
        return f"{element}_{int(used[0][1])}A.fld"
    return element + "".join(f"_{w}_{int(v)}A" for w, v in used) + ".fld"


def parse_name(element, fname):
    """Обратный разбор имени в набор токов; None, если имя не наше."""
    base = os.path.basename(fname)
    if not base.startswith(element + "_") or not base.endswith(".fld"):
        return None
    body = base[len(element) + 1:-4]
    if DESIGN[element]["bare"]:
        m = re.fullmatch(r"(-?\d+)A", body)
        return {order_of(element)[0]: int(m.group(1))} if m else None
    # Разбор ПО ТОКЕНАМ, а не регуляркой: в имени `main_0A_c1_c2_c3_c4_1808A` один ток относится
    # сразу к четырём обмоткам, а регулярка `(имя)_(число)A` вдобавок затягивала подчёркивание в
    # имя обмотки (`_I2` вместо `I2`) и молча отбрасывала все комбинации.
    known, out, pend = set(order_of(element)), {}, []
    for tok in body.split("_"):
        m = re.fullmatch(r"(-?\d+)A", tok)
        if m:
            if not pend:
                return None
            for w in pend:
                out[w] = int(m.group(1))
            pend = []
        elif tok in known:
            pend.append(tok)
        else:
            return None
    return out if out and not pend else None


def excitation(element, currents):
    """Строка возбуждения в конвенции реестра: `I1=1891;I2=800`."""
    return ";".join(f"{w}={int(currents[w])}" for w in order_of(element) if currents.get(w))


def map_id(element, currents):
    used = [w for w in order_of(element) if currents.get(w)]
    return f"{element}.{'+'.join(used)}." + "+".join(str(int(currents[w])) for w in used)


# ───────────────────────────── сетка и заказ ─────────────────────────────

def reference_grid(element):
    """Сетка и эталонная шапка — из уже снятой карты этого элемента. Руками не задаём."""
    for f in sorted(glob.glob(os.path.join(folder_of(element), f"{element}_*.fld"))):
        with open(f, encoding="utf-8", errors="ignore") as fh:
            head = fh.readline()
        m = HDR.search(head)
        if m:
            return dict(sample=f, header=head.rstrip("\n"),
                        lo=m.group(1).split(), hi=m.group(2).split(), step=m.group(3).split())
    raise SystemExit(f"{element}: не нашёл ни одной снятой карты, не с чем сверить сетку")


def jobs(element=None, path=REQUEST):
    """Заказ: строки `элемент;возбуждение;назначение`, возбуждение вида `I1=1477;I2=800`."""
    if not os.path.exists(path):
        raise SystemExit(f"нет таблицы заказа {path}")
    out = []
    for r in csv.DictReader(open(path, encoding="utf-8"), delimiter=";"):
        el = r["элемент"]
        if element and el != element:
            continue
        cur = {}
        for part in r["возбуждение"].split(";"):
            if "=" in part:
                w, v = part.split("=", 1)
                cur[w.strip()] = int(float(v))
        nm = map_name(el, cur)
        out.append(dict(element=el, currents=cur, name=nm, folder=folder_of(el),
                        target=os.path.join(folder_of(el), nm), why=r.get("назначение", "")))
    return out


# ───────────────────────────── скрипт для AEDT ─────────────────────────────

def build_script(element, todo, grid, out_path, reset=False):
    """Скрипт для интерпретатора AEDT. Внутри доступен только oDesktop — библиотек не подключаем.

    Проект НЕ сохраняем и в конце возвращаем токам прежние значения: расчёт не должен менять
    состояние чужого проекта на диске. Незадействованные обмотки ЯВНО зануляются — иначе точка
    унаследует ток от предыдущей."""
    d = DESIGN[element]
    lines = [
        "# -*- coding: utf-8 -*-",
        "# Сгенерировано ansys_maps.py. Запуск: ansysedt.exe -ng -runscriptandexit <этот файл>",
        "import ScriptEnv",
        'ScriptEnv.Initialize("Ansoft.ElectronicsDesktop")',
        # ГОТЧА: SetActiveProject не возвращает None для незагруженного проекта, а БРОСАЕТ
        # исключение. В пакетном режиме открыт вообще ничего не бывает, поэтому открываем сами.
        f'oProject = oDesktop.OpenProject(r"{os.path.join(PROJECTS, d["project"])}.aedt")',
        "if oProject is None:",
        f'    oProject = oDesktop.SetActiveProject("{d["project"]}")',
        f'oDesign = oProject.SetActiveDesign("{d["design"]}")',
        'oModule = oDesign.GetModule("FieldsReporter")',
        # Сброс сохранённых решений: без него Analyze молча отдаёт готовое, рапортуя «Normal
        # completion», и проверка воспроизводимости вырождается в тавтологию — ровно так впустую
        # прошёл первый пересчёт QP_L.
        *(['oDesign.DeleteFullVariation("All", False)'] if reset else []),
        f"VARS = {list(d['vars'].values())!r}",
        "OLD = dict((v, oDesign.GetVariableValue(v)) for v in VARS)",
        'oDesktop.AddMessage("", "", 0, "исходные токи: " + str(OLD))',
        "",
        f"LO = {grid['lo']!r}",
        f"HI = {grid['hi']!r}",
        f"STEP = {grid['step']!r}",
        "",
    ]
    for j in todo:
        vals = {d["vars"][w]: int(j["currents"].get(w, 0)) for w in order_of(element)}
        var_list = ", ".join(f'"{v}:=", "{a}A"' for v, a in vals.items())
        lines.append(f'oDesktop.AddMessage("", "", 0, "=== {j["name"]} ===")')
        for v, a in vals.items():
            lines.append(f'oDesign.SetVariableValue("{v}", "{a}A")')
        lines += [
            f'oDesign.Analyze("{d["setup"]}", True)',
            'oModule.CalcStack("clear")',
            # Именованное выражение, а НЕ сырое EnterQty("B"): иначе величина в файле называется
            # "<Bx,By,Bz>" вместо "B_Vector" и шапка расходится с прежними картами.
            f'oModule.CopyNamedExprToStack("{QUANTITY}")',
        ] + ([f'oModule.CalcOp("{SMOOTH}")'] if SMOOTH else []) + [
            "oModule.ExportOnGrid(",
            f'    r"{j["target"]}",',
            "    LO, HI, STEP,",
            f'    "{d["setup"]} : LastAdaptive",',
            f"    [{var_list}],",
            # exportoption=True — писать столбцы X, Y, Z. При False AEDT кладёт только компоненты
            # поля, файл выходит втрое короче и несопоставим с прежними картами.
            '    True, "Cartesian", [0, 0, 0], False)',
            "",
        ]
    lines += [
        "for v in VARS:",
        "    oDesign.SetVariableValue(v, OLD[v])",          # вернуть как было
        'oDesktop.AddMessage("", "", 0, "=== все точки пройдены ===")',
    ]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    open(out_path, "w", encoding="utf-8", newline="\n").write("\n".join(lines) + "\n")
    return out_path


# ───────────────────────────── сверка материала ─────────────────────────────

def project_bh(element):
    """Кривая намагничивания, ЛЕЖАЩАЯ В ПРОЕКТЕ прямо сейчас (а не та, что в имени материала).

    ИМЯ МАТЕРИАЛА ЛЖЁТ: в QP_sh материал звался `steel_1008 _custom`, а содержал
    низконасыщающуюся кривую (B_max 1.70 Тл вместо 2.59) — остаток от опытов со сменой материала.
    Карты, снятые на ней, несопоставимы с производственной серией, и по имени этого не видно."""
    import numpy as np
    p = os.path.join(PROJECTS, f"{DESIGN[element]['project']}.aedt")
    txt = open(p, encoding="utf-8", errors="ignore").read()
    # ГОТЧА: в проекте НЕСКОЛЬКО материалов с кривыми (у QP_sh четыре — три библиотечных плюс
    # правленый). Брать ПЕРВЫЙ блок BHCoordinates нельзя: он библиотечный и к железу отношения
    # не имеет. Ищем кривую НАЗНАЧЕННОГО материала.
    assigned = [m for m in re.findall(r"MaterialValue='\"([^\"]+)\"'", txt)
                if m.lower() not in ("copper", "vacuum", "air")]
    if not assigned:
        return None, None
    name = max(set(assigned), key=assigned.count)      # тот, которым покрыто больше всего тел
    i = txt.find(f"$begin '{name}'")
    if i < 0:
        return name, None
    j = txt.find("$begin 'BHCoordinates'", i)
    if j < 0 or j > txt.find(f"$end '{name}'", i):
        return name, None
    m = BH_BLOCK.search(txt, j)
    if not m:
        return name, None
    v = np.array([float(x) for x in m.group(1).replace("\n", " ").split(",") if x.strip()])
    return name, v.reshape(-1, 2)


def check_material(element, expected="workspace/materials/steel_1008_bh.csv", tol=1e-4):
    """Сверка кривой проекта с ожидаемой: (совпало, имя, B_max проекта, B_max файла)."""
    import numpy as np
    name, bh = project_bh(element)
    if bh is None:
        return None, name, None, None
    ref = np.loadtxt(expected, delimiter=";", skiprows=1)
    same = (len(ref) == len(bh) and float(np.abs(ref - bh).max()) < tol)
    return same, name, float(bh[:, 1].max()), float(ref[:, 1].max())


def guard_material(element, expected="workspace/materials/steel_1008_bh.csv"):
    """Не считать, пока в проекте не та кривая намагничивания.

    Молча снятая серия на чужом материале несопоставима с прежними картами, а по ИМЕНИ материала
    этого не видно. Поймали тем, что пересчёт уже снятой карты 6588 А·в дал кик на 36% слабее при
    полном совпадении в линейном режиме: геометрия и возбуждение те же, разошлось только
    насыщение."""
    same, name, bmax_p, bmax_r = check_material(element, expected)
    if same is None:
        print("  ВНИМАНИЕ: кривую из проекта прочитать не удалось — сверка материала пропущена")
        return
    if same:
        print(f"  материал: {name!r} совпал с {os.path.basename(expected)} "
              f"(B_max {bmax_p:.4f} Тл)")
        return
    print("  МАТЕРИАЛ В ПРОЕКТЕ НЕ ТОТ:")
    print(f"    в проекте: {name!r}, B_max {bmax_p:.4f} Тл")
    print(f"    ожидается: {expected}, B_max {bmax_r:.4f} Тл")
    print("    Имя материала может лгать — сверяется САМА КРИВАЯ. Восстановите её в AEDT")
    print("    либо укажите через --material файл, на котором серия задумана.")
    raise SystemExit("остановлено: несовпадение кривой намагничивания")


# ───────────────────────────── проверка снятого ─────────────────────────────

def probe(path, order=2, r=1e-3, nth=64):
    """Физическая проверка карты: максимум |B| и интегральная гармоника.

    Сверки шапки НЕ ДОСТАТОЧНО. Экспорт, сделанный до окончания расчёта или без пересчёта после
    смены тока, даёт файл с идеально правильной шапкой и полем от ПРЕЖНЕГО решения."""
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


def measure(element, order=None, winding=None):
    """Снятые карты элемента по возрастанию тока. Кривая насыщения обязана РАСТИ: падение при
    большем токе — испорченный экспорт, а не физика.

    Монотонность проверяется ТОЛЬКО по одиночным обмоткам: у комбинаций ток не одномерен, и
    сравнивать их по сумме модулей бессмысленно."""
    order = order or ORDER_BY_ELEMENT.get(element, 2)
    rows = []
    for f in sorted(glob.glob(os.path.join(folder_of(element), f"{element}_*.fld"))):
        cur = parse_name(element, f)
        if not cur:
            continue
        if winding and set(cur) != {winding}:
            continue
        mx, hm = probe(f, order, probe_radius(element))
        rows.append((sum(abs(v) for v in cur.values()), os.path.basename(f), mx, hm, cur))
    rows.sort()
    # Критерий — ГАРМОНИКА (то, что видит пучок), а не максимум |B|: максимум сидит в железе и при
    # насыщении законно проседает (QP_sh 9882→13176: max|B| 2.974→2.907, гармоника растёт).
    by_w = {}
    for r in rows:
        if len(r[4]) == 1:
            by_w.setdefault(next(iter(r[4])), []).append(r)
    bad = []
    for seq in by_w.values():
        bad += [r for i, r in enumerate(seq) if i and r[3] < seq[i - 1][3] * 0.99]
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
        return dict(header_same=head[0] == head[1], shape_same=False, shape=(A.shape, B.shape))
    d = np.linalg.norm(A - B, axis=-1)
    scale = float(np.linalg.norm(A, axis=-1).max())
    return dict(header_same=head[0] == head[1], shape_same=True, shape=A.shape,
                max_abs=float(d.max()), rms=float(np.sqrt((d ** 2).mean())),
                max_rel=100.0 * float(d.max()) / scale if scale else float("nan"),
                identical=bool(np.array_equal(A, B)), scale=scale)


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


# ───────────────────────────── регистрация ─────────────────────────────

def cmd_register(element, material="workspace/materials/steel_1008_bh.csv"):
    """Завести снятые карты в реестр: .fld → копия в fields_raw, .bin в fields_bin, строка в CSV.

    Переконвертирует и УЖЕ зарегистрированные карты, если их .fld в KD новее .bin: после
    `--recompute` файлы в KD обновились, а .bin остался прежним, и модель молча читала бы
    старые данные."""
    import shutil

    from field_quality import load
    from symmetrize_map import write_bin

    reg_path = os.path.join("workspace", "tables", "field_maps.csv")
    rows = list(csv.DictReader(open(reg_path, encoding="utf-8"), delimiter=";"))
    cols = list(rows[0].keys())
    have = {r["map_id"] for r in rows}
    proto = next(r for r in rows if r["element"] == element)
    added = reconv = 0
    for f in sorted(glob.glob(os.path.join(folder_of(element), f"{element}_*.fld"))):
        cur = parse_name(element, f)
        if not cur:
            continue
        mid = map_id(element, cur)
        raw = f"model/{element}/{mid}.fld"
        bin_ = f"model/{element}/{mid}.bin"
        raw_abs = os.path.join("workspace", "fields_raw", *raw.split("/"))
        bin_abs = os.path.join("workspace", "fields_bin", *bin_.split("/"))
        if not os.path.exists(bin_abs) or os.path.getmtime(f) > os.path.getmtime(bin_abs):
            os.makedirs(os.path.dirname(raw_abs), exist_ok=True)
            os.makedirs(os.path.dirname(bin_abs), exist_ok=True)
            shutil.copy2(f, raw_abs)
            h, B = load(f)
            write_bin(h, B, bin_abs, nominal_arg=float(sum(abs(v) for v in cur.values())))
            if mid in have:
                reconv += 1
                print(f"  переконвертирована {mid} (.fld новее .bin)")
        if mid in have:
            continue
        turns = DESIGN[element]["turns"]
        amps = ";".join(f"{w}={cur[w] / turns.get(w, 1):.2f}"
                        for w in order_of(element) if cur.get(w))
        r = {c: "" for c in cols}
        r.update({"map_id": mid, "element": element,
                  "winding": "+".join(w for w in order_of(element) if cur.get(w)),
                  "device_type_id": proto["device_type_id"],
                  "amp_turns": str(sum(abs(v) for v in cur.values())),
                  "current_a": amps, "excitation": excitation(element, cur),
                  "split": "test" if added % 3 == 2 else "train",
                  "grid": proto["grid"], "material_curve": material,
                  "raw_path": raw, "bin_path": bin_, "source_fld": f.replace(os.sep, "/"),
                  "note": "автоматическое снятие"})
        rows.append(r)
        added += 1
        print(f"  добавлена {mid}  ({excitation(element, cur)}, {r['split']})")
    if added:
        # БЕЗ BOM: read_registry открывает реестр как обычный utf-8, и с BOM первый столбец
        # превращается в "\ufeffmap_id" — вся цепочка падает на KeyError('map_id').
        with open(reg_path, "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter=";")
            w.writeheader()
            w.writerows(rows)
    print(f"  итог: добавлено {added}, переконвертировано {reconv}, всего у {element} "
          f"{sum(1 for r in rows if r['element'] == element)} карт")


# ───────────────────────────── команды счёта ─────────────────────────────

def cmd_recompute(element, outdir, backup="_original", everything=False):
    """Пересчитать УЖЕ СНЯТЫЕ карты, отложив оригиналы в сторону.

    Зачем: серия должна быть снята ОДНИМ состоянием проекта. Между прогонами разных дней сидит
    систематический сдвиг ~2.6% (замерено на шести картах QP_sh), а это порядок всей достигнутой
    точности. Равномерный масштаб модель поглощает, но неравномерность — нет."""
    import numpy as np
    grid = reference_grid(element)                       # ДО переноса: иначе не с чем сверять
    folder = folder_of(element)
    rows, _ = measure(element)
    if not rows:
        raise SystemExit(f"{element}: нечего пересчитывать")
    if not everything:
        newest = max(os.path.getmtime(os.path.join(folder, r[1])) for r in rows)
        old = [r for r in rows if os.path.getmtime(os.path.join(folder, r[1])) < newest - 86400]
        if old:
            print(f"  из {len(rows)} карт пересчитываются {len(old)} из прежних прогонов "
                  f"(остальные сняты в последние сутки и уже согласованы)")
            rows = old
    bak = os.path.join(folder, backup)
    os.makedirs(bak, exist_ok=True)
    todo = []
    for _, nm, _, _, cur in rows:
        src, dst = os.path.join(folder, nm), os.path.join(bak, nm)
        if os.path.exists(dst) and os.path.exists(src):
            os.remove(src)
        elif os.path.exists(src):
            os.replace(src, dst)
        todo.append(dict(element=element, currents=cur, name=nm, folder=folder,
                         target=src, why="пересчёт"))
    print(f"  оригиналы отложены в {bak} ({len(todo)} шт)")
    guard_material(element)
    sp = build_script(element, todo, grid, os.path.join(outdir, f"aedt_recompute_{element}.py"))
    print(f"  скрипт: {sp}\n  запуск AEDT ({len(todo)} точек)…")
    rc = subprocess.call([AEDT_EXE, "-ng", "-runscriptandexit", os.path.abspath(sp)])
    print(f"  AEDT вернул {rc}\n")
    print(f"  {'карта':34} {'прежде':>10} {'теперь':>10} {'расхождение':>12}")
    diffs, cached = [], []
    for j in todo:
        old_f, new_f = os.path.join(bak, j["name"]), j["target"]
        if not (os.path.exists(old_f) and os.path.exists(new_f)):
            print(f"  {j['name']:34}   нет файла"); continue
        _, ho = probe(old_f, ORDER_BY_ELEMENT.get(element, 2))
        _, hn = probe(new_f, ORDER_BY_ELEMENT.get(element, 2))
        d = 100.0 * (hn - ho) / ho if ho else float("nan")
        diffs.append(d)
        same = filecmp.cmp(old_f, new_f, shallow=False)
        cached.append(same)
        print(f"  {j['name']:34} {ho:10.5g} {hn:10.5g} {d:11.2f}%"
              + ("   ← БАЙТ-В-БАЙТ: решение взято из кэша" if same else ""))
    if cached and all(cached):
        # Maxwell хранит решения по вариациям и на повторный Analyze отдаёт готовое, рапортуя
        # «Normal completion». Пустой прогон выглядит как ИДЕАЛЬНЫЙ результат (разброс 0.00%),
        # хотя означает, что не пересчитано ничего. Отличаем по времени: кэш ~35 с на точку
        # против ~4 мин у настоящего счёта, и по побайтовому совпадению файлов.
        print("  ПЕРЕСЧЁТА НЕ БЫЛО: все файлы совпали байт-в-байт, Maxwell отдал решения из кэша.")
        print("  Это не «нулевой разброс», а отсутствие счёта. Зато отсюда следует, что проект")
        print("  с тех пор НЕ МЕНЯЛСЯ (иначе решения были бы сброшены), и серия самосогласована.")
        print("  Принудительный пересчёт: --fresh (удалит сохранённые решения вариаций).")
    elif diffs:
        print(f"  разброс прогон-к-прогону: средний модуль {np.mean(np.abs(diffs)):.2f}%, "
              f"максимум {np.max(np.abs(diffs)):.2f}%")
    return diffs


def cmd_reproduce(element, outdir, count=2, material="workspace/materials/steel_1008_bh.csv"):
    """Проверка ВОСПРОИЗВОДИМОСТИ: пересчитать самые старые карты и сверить с тем, что на диске.

    Оригиналы НЕ трогаются — пересчёт идёт в отдельную папку. Решения принудительно сбрасываются
    (`DeleteFullVariation`), иначе Maxwell молча отдаст готовое и проверка выродится в тавтологию:
    именно так впустую прошёл первый пересчёт QP_L.

    Смысл: если проект не менялся с тех пор, как снимались старые карты, они обязаны воспроизвестись
    точно. Так проверен DQ (0.00% на трёх картах). Если не воспроизводятся — в проекте что-то
    поправили, и карты разных дат смешивать нельзя. У QP_sh расхождение оказалось до 3.5%, и
    виновата была СЕТКА."""
    import numpy as np
    grid = reference_grid(element)
    folder = folder_of(element)
    rows, _ = measure(element)
    singles = [r for r in rows if len(r[4]) == 1]
    if not singles:
        raise SystemExit(f"{element}: нет одиночных карт")
    oldest = sorted(singles, key=lambda r: os.path.getmtime(os.path.join(folder, r[1])))[:count]
    os.makedirs(outdir, exist_ok=True)
    todo = [dict(element=element, currents=cur, name=nm, folder=outdir,
                 target=os.path.join(os.path.abspath(outdir), nm), why="воспроизводимость")
            for _, nm, _, _, cur in oldest]
    guard_material(element, material)
    sp = build_script(element, todo, grid, os.path.join(outdir, f"aedt_reproduce_{element}.py"),
                      reset=True)
    print(f"  пересчитываются {len(todo)} самых старых карт, решения сбрасываются")
    rc = subprocess.call([AEDT_EXE, "-ng", "-runscriptandexit", os.path.abspath(sp)])
    print(f"  AEDT вернул {rc}")
    order, R = ORDER_BY_ELEMENT.get(element, 2), probe_radius(element)
    print(f"  проба на {R * 1e3:.0f} мм, гармоника порядка {order}")
    print(f"  {'карта':34} {'на диске':>11} {'пересчёт':>11} {'расхождение':>13}")
    d = []
    for j in todo:
        old = os.path.join(folder, j["name"])
        if not os.path.exists(j["target"]):
            print(f"  {j['name']:34}   файл не появился"); continue
        _, ho = probe(old, order, R)
        _, hn = probe(j["target"], order, R)
        d.append(100.0 * (hn - ho) / ho if ho else float("nan"))
        print(f"  {j['name']:34} {ho:11.6g} {hn:11.6g} {d[-1]:12.2f}%")
    if d:
        m = float(np.max(np.abs(d)))
        print(f"  ИТОГ: максимум {m:.2f}% — "
              + ("проект НЕ менялся, карты разных дат совместимы" if m < 1.5
                 else "проект МЕНЯЛСЯ, карты разных дат смешивать нельзя"))
    return d


def cmd_repro(element, spec, outdir):
    """Пересчитать УЖЕ СНЯТУЮ карту в отдельный файл и сравнить с ней. Оригинал не трогаем."""
    grid = reference_grid(element)
    cur = {}
    for part in spec.split(";"):
        if "=" in part:
            w, v = part.split("=", 1)
            cur[w.strip()] = int(float(v))
        else:
            cur[order_of(element)[0]] = int(float(part))
    orig = os.path.join(folder_of(element), map_name(element, cur))
    if not os.path.exists(orig):
        raise SystemExit(f"нет исходной карты {orig}")
    os.makedirs(outdir, exist_ok=True)
    nm = map_name(element, cur).replace(".fld", "__repro.fld")
    repro = os.path.join(os.path.abspath(outdir), nm)
    job = [dict(element=element, currents=cur, name=nm, folder=outdir,
                target=repro, why="воспроизведение")]
    guard_material(element)
    sp = build_script(element, job, grid, os.path.join(outdir, f"aedt_repro_{element}.py"))
    print(f"  оригинал: {orig}\n  пересчёт: {repro}\n  скрипт:   {sp}\n  запуск AEDT…")
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


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--element", help="только этот элемент (иначе все из таблицы)")
    ap.add_argument("--plan", action="store_true", help="показать план, ничего не делать")
    ap.add_argument("--script", action="store_true", help="собрать скрипт для AEDT, не запускать")
    ap.add_argument("--run", action="store_true",
                    help="собрать И ЗАПУСТИТЬ (считает, берёт лицензию)")
    ap.add_argument("--verify", action="store_true", help="сверить шапки уже снятых карт")
    ap.add_argument("--measure", action="store_true",
                    help="физическая проверка снятых карт: гармоника обязана расти с током")
    ap.add_argument("--register", action="store_true",
                    help="завести снятые карты в реестр (.fld -> .bin -> строка в field_maps.csv)")
    ap.add_argument("--recompute", action="store_true",
                    help="пересчитать УЖЕ СНЯТЫЕ карты (оригиналы переезжают в подпапку)")
    ap.add_argument("--all", action="store_true",
                    help="с --recompute: пересчитать ВСЕ карты, а не только из прежних прогонов")
    ap.add_argument("--reproduce", type=int, nargs="?", const=2, metavar="N",
                    help="проверить воспроизводимость: пересчитать N самых старых карт со сбросом "
                         "решений и сверить с диском (оригиналы не трогаются)")
    ap.add_argument("--repro", metavar="ВОЗБУЖДЕНИЕ",
                    help="пересчитать снятую карту и сравнить, напр. 3294 или I1=1891;I2=800")
    ap.add_argument("--order", type=int, help="порядок гармоники для проверки")
    ap.add_argument("--material", default="workspace/materials/steel_1008_bh.csv",
                    help="кривая, на которой ДОЛЖНА сниматься серия")
    ap.add_argument("--limit", type=int, help="взять только первые N точек из оставшихся")
    ap.add_argument("--request", default=REQUEST, help="таблица заказа")
    ap.add_argument("--outdir", default=os.path.join("workspace", "results", "aedt"))
    a = ap.parse_args(argv)

    if a.repro or a.recompute or a.register or a.measure or a.reproduce:
        if not a.element:
            raise SystemExit("этот режим требует --element")
        el = a.element
        if a.reproduce:
            print(f"=== {el}: проверка воспроизводимости старых карт")
            cmd_reproduce(el, os.path.join(a.outdir, f"reproduce_{el}"), a.reproduce, a.material)
        elif a.repro:
            cmd_repro(el, a.repro, a.outdir)
        elif a.recompute:
            print(f"=== {el}: пересчёт снятой серии одним состоянием проекта")
            cmd_recompute(el, a.outdir, everything=a.all)
        elif a.register:
            print(f"=== {el}: регистрация карт")
            cmd_register(el, a.material)
        else:
            rows, bad = measure(el, a.order)
            print(f"=== {el}: снято {len(rows)} карт, гармоника порядка "
                  f"{a.order or ORDER_BY_ELEMENT.get(el, 2)}")
            print(f"  {'|ток|, А·в':>11} {'файл':40} {'макс |B|, Тл':>13} {'гармоника':>12}")
            for row in rows:
                flag = "  ← ПАДЕНИЕ ПРИ БОЛЬШЕМ ТОКЕ" if row in bad else ""
                print(f"  {row[0]:11d} {row[1]:40} {row[2]:13.4f} {row[3]:12.5g}{flag}")
            print(f"  подозрительных карт: {len(bad)}")
        return 0

    els = [a.element] if a.element else sorted({j["element"] for j in jobs(path=a.request)})
    for el in els:
        if el not in DESIGN:
            print(f"{el}: проект не описан в DESIGN — пропускаю"); continue
        grid = reference_grid(el)
        todo = jobs(el, a.request)
        left = [j for j in todo if not os.path.exists(j["target"])]
        print(f"\n=== {el}: заказано {len(todo)}, уже снято {len(todo) - len(left)}, "
              f"осталось {len(left)}")
        print(f"  сетка из эталона {os.path.basename(grid['sample'])}: "
              f"min {' '.join(grid['lo'])}  max {' '.join(grid['hi'])}  "
              f"шаг {' '.join(grid['step'])}")
        if a.verify:
            ok, bad, missing = verify(el, todo, grid)
            print(f"  сверка шапок: совпало {len(ok)}, РАСХОЖДЕНИЙ {len(bad)}, "
                  f"нет файла {len(missing)}")
            for n, h in bad:
                print(f"    ! {n}\n      {h}")
            continue
        if a.limit:
            left = left[:a.limit]
        for j in left:
            print(f"  {excitation(el, j['currents']):26} → {j['name']:36} {j['why']}")
        if a.plan or not left:
            continue
        sp = build_script(el, left, grid, os.path.join(a.outdir, f"aedt_job_{el}.py"))
        print(f"  скрипт: {sp}")
        if not a.run:
            continue
        lock = os.path.join(PROJECTS, f"{DESIGN[el]['project']}.aedt.lock")
        if os.path.exists(lock):
            raise SystemExit(f"проект открыт в AEDT ({lock}) — закройте его перед прогоном")
        guard_material(el, a.material)
        print(f"  запуск AEDT ({len(left)} точек, это надолго)…")
        rc = subprocess.call([AEDT_EXE, "-ng", "-runscriptandexit", os.path.abspath(sp)])
        print(f"  AEDT вернул {rc}")
        ok, bad, missing = verify(el, todo, grid)
        print(f"  итог: снято {len(ok)}, расхождений шапки {len(bad)}, "
              f"не появилось {len(missing)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
