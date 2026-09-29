# -*- coding: utf-8 -*-
"""Приведение чужого .lte к плоскому виду, который принимают наши загрузчики и elegant-сверка.

ЗАЧЕМ. Файлы установок приходят в разном синтаксисе: многострочные определения через `&`, вложенные
линии, повторители `38*CELL`, матричные типы элементов, служебные элементы (импеданс, апертура, ВЧ).
`s2lattice_to_db.parse()` ждёт ОДНО определение на строку и одну развёрнутую линию `ring: LINE=(...)`,
а сверка с elegant — КИКОВЫЕ типы (иначе elegant посчитает элемент матрицей и `N_SLICES` на него не
подействует, см. docs/MANUAL.md). Этот инструмент делает и то и другое за один проход.

ЧТО ДЕЛАЕТ
  • Склеивает продолжения строк: `&` в конце переносит ТОЛЬКО следующую строку (пустая строка обрывает
    определение — так же, как это читает elegant).
  • Разворачивает выбранную линию рекурсивно: вложенные LINE, `n*ИМЯ`, `-ИМЯ` (обратный порядок).
  • Переводит МАТРИЧНЫЕ типы в КИКОВЫЕ: SBEN/RBEN/BEND → CSBEND, QUAD → KQUAD, SEXT → KSEXT,
    OCTU/OCT → KOCT, EDRIFT/DRIFT → DRIF. Параметры (L, ANGLE, K1…K3, E1, E2, TILT) переносятся как есть.
  • Выбрасывает то, чего нет в матричной модели: ВЧ, апертуры, импеданс, заряд, WATCH. Что именно
    выброшено — печатает; ВЧ подставляет сама сверка (elegant_diag.write_lattice).
  • Ставит N_SLICES на киковые элементы (по умолчанию 5 — конвенция проекта для СИЛА).

ЧЕГО НЕ ДЕЛАЕТ. Не трогает силы и углы: замыканием по углу занимается lte_ring.py, тюном — оптимизатор.

ЗАПУСК
  py -3.14 lte_flatten.py вход.lte --beamline R --out workspace/lte/configuration_v5/sila.lte
  py -3.14 lte_flatten.py вход.lte --list                 # какие линии есть в файле
"""
import argparse
import collections
import io
import math
import os
import re

import lte_to_lattice as _geom          # ЕДИНСТВЕННЫЙ источник конвенции знака (step)

# матричный тип elegant → киковый (тот же магнит, но с интегратором по срезам)
KICKIFY = {"SBEN": "CSBEND", "SBEND": "CSBEND", "RBEN": "CSBEND", "RBEND": "CSBEND", "BEND": "CSBEND",
           "CSBEND": "CSBEND", "CCBEND": "CSBEND", "KSBEND": "CSBEND",
           "QUAD": "KQUAD", "KQUAD": "KQUAD", "QUADRUPOLE": "KQUAD",
           "SEXT": "KSEXT", "KSEXT": "KSEXT", "SEXTUPOLE": "KSEXT",
           "OCTU": "KOCT", "OCT": "KOCT", "KOCT": "KOCT", "OCTUPOLE": "KOCT",
           "DRIF": "DRIF", "DRIFT": "DRIF", "EDRIFT": "DRIF", "CSRDRIFT": "DRIF",
           "MONI": "MONI", "MONITOR": "MONI", "MARK": "MARK", "MARKER": "MARK"}
SLICED = ("CSBEND", "KQUAD", "KSEXT", "KOCT")          # у кого есть N_SLICES
KEEP = ("L", "ANGLE", "K1", "K2", "K3", "E1", "E2", "TILT", "HGAP", "FINT")
DROP_NOTE = {"RFCA": "ВЧ (подставляет сверка)", "RFCW": "ВЧ", "MAXAMP": "апертура",
             "WATCH": "точка наблюдения elegant", "CHARGE": "заряд пучка",
             "ZLONGIT": "продольный импеданс", "ZTRANSVERSE": "поперечный импеданс",
             "SCRAPER": "скрепер", "MALIGN": "сдвиг опорной частицы", "RCOL": "коллиматор",
             "ECOL": "коллиматор", "SREFFECTS": "модель излучения elegant"}


def read(path):
    """Определения и линии из .lte. Возвращает (elems, lines, порядок объявления линий)."""
    txt = io.open(path, encoding="utf-8", errors="replace").read()
    txt = re.sub(r"^\s*!.*$", "", txt, flags=re.M)             # комментарии
    txt = re.sub(r"&[ \t]*\r?\n", " ", txt)                    # продолжение — ТОЛЬКО следующая строка
    elems, lines, order = {}, {}, []
    for m in re.finditer(r"^\s*\"?([A-Za-z0-9_.]+)\"?\s*:\s*(.*)$", txt, re.M):
        name, body = m.group(1).upper(), m.group(2).strip()
        if re.match(r"(?i)line\s*=", body):
            lines[name] = [t.strip() for t in body[body.index("(") + 1:body.rindex(")")].split(",") if t.strip()]
            order.append(name)
            continue
        t = re.match(r"([A-Za-z_]+)", body)
        if not t:
            continue
        d = {"type": t.group(1).upper()}
        for kv in re.finditer(r"([A-Za-z0-9_]+)\s*=\s*(-?[\d.]+(?:[eE][-+]?\d+)?)", body):
            d[kv.group(1).upper()] = float(kv.group(2))
        elems[name] = d
    return elems, lines, order


def expand(name, lines, depth=0):
    """Развёртка линии: вложенные LINE, `n*ИМЯ`, `-ИМЯ` (обратный порядок)."""
    if depth > 50:
        raise RuntimeError(f"линия {name}: вложенность глубже 50 — похоже на цикл")
    out = []
    for tok in lines[name]:
        t, rev = tok.strip(), False
        if t.startswith("-"):
            t, rev = t[1:].strip(), True
        n = 1
        m = re.match(r"(\d+)\s*\*\s*(.+)", t)
        if m:
            n, t = int(m.group(1)), m.group(2).strip()
        t = t.strip("()").upper()
        piece = expand(t, lines, depth + 1) if t in lines else [t]
        if rev:
            piece = piece[::-1]
        out += piece * n
    return out


def num(v):
    """Число для .lte: без экспоненты там, где можно — наши парсеры читают простой формат."""
    s = f"{v:.15g}"
    return s


def flatten(src, beamline=None, n_slices=5, line_name="ring"):
    """Читает src и возвращает (текст плоского .lte, отчёт-словарь)."""
    elems, lines, order = read(src)
    if beamline is None:
        for cand in ("RING", "R", "LATTICE", "MACHINE", "ALL"):
            if cand in lines:
                beamline = cand
                break
        else:
            beamline = order[-1] if order else None
    if beamline is None or beamline.upper() not in lines:
        raise SystemExit(f"линия {beamline!r} не найдена; есть: {', '.join(sorted(lines))}")
    seq = expand(beamline.upper(), lines)

    dropped = collections.Counter()
    kept = []
    for t in seq:
        e = elems.get(t)
        if e is None:
            raise SystemExit(f"элемент {t} встречается в линии, но не объявлен")
        if e["type"] not in KICKIFY:
            dropped[e["type"]] += 1
            continue
        kept.append(t)

    used = sorted({t for t in kept}, key=kept.index)
    out = [f"! сгенерировано lte_flatten.py из {os.path.basename(src)}, линия {beamline}",
           f"! матричные типы переведены в киковые, N_SLICES={n_slices}; служебные элементы убраны"]
    for t in used:
        e = elems[t]
        typ = KICKIFY[e["type"]]
        par = [f"{k}={num(e[k])}" for k in KEEP if k in e and e[k] != 0.0 or (k == "L" and "L" in e)]
        par = list(dict.fromkeys(par))
        if typ in SLICED:
            par.append(f"N_SLICES={n_slices}")
        out.append(f"{t}: {typ}" + (", " + ", ".join(par) if par else ""))
    out.append("")
    body, row = [], []
    for t in kept:
        row.append(t)
        if len(row) == 12:
            body.append(",".join(row)); row = []
    if row:
        body.append(",".join(row))
    out.append(f"{line_name}: LINE=(" + ",&\n".join(body) + ")")

    x = y = th = 0.0
    ang = 0.0
    for t in kept:
        e = elems[t]
        _, _, _, x, y, th = _geom.step(x, y, th, e.get("L", 0.0), e.get("ANGLE", 0.0))
        ang += e.get("ANGLE", 0.0)
    rep = {"beamline": beamline, "позиций": len(seq), "оставлено": len(kept), "типов": len(used),
           "периметр, м": sum(elems[t].get("L", 0.0) for t in kept), "ΣANGLE/2π": ang / (2 * math.pi),
           "зазор кольца, мм": math.hypot(x, y) * 1e3, "выброшено": dict(dropped),
           "состав": dict(collections.Counter(KICKIFY[elems[t]["type"]] for t in kept))}
    return "\n".join(out) + "\n", rep


def main():
    ap = argparse.ArgumentParser(description="чужой .lte → плоский киковый .lte для загрузчиков и elegant")
    ap.add_argument("src", help="исходный .lte")
    ap.add_argument("--beamline", default=None, help="имя линии кольца (по умолчанию RING/R/…)")
    ap.add_argument("--out", default=None, help="куда писать; без него только отчёт")
    ap.add_argument("--line-name", default="ring", help="как назвать линию в выходном файле")
    ap.add_argument("--n-slices", type=int, default=5, help="N_SLICES у киковых элементов")
    ap.add_argument("--list", action="store_true", help="показать линии в файле и выйти")
    a = ap.parse_args()
    if a.list:
        _, lines, order = read(a.src)
        for n in order:
            print(f"  {n:20s} {len(lines[n])} токенов")
        return
    txt, rep = flatten(a.src, a.beamline, a.n_slices, a.line_name)
    for k, v in rep.items():
        print(f"  {k}: {v}")
    if a.out:
        os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
        io.open(a.out, "w", encoding="utf-8", newline="\n").write(txt)
        print(f"[OK] {a.out}")


if __name__ == "__main__":
    main()
