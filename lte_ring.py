# -*- coding: utf-8 -*-
"""Проверка и замыкание кольца в .lte — УНИВЕРСАЛЬНЫЙ инструмент, работает по ПУТИ к файлу.

Ничего не знает про конкретные установки: ни имён файлов, ни папок, ни названий линий. Любой .lte
(elegant/PyAT) годится, лишь бы в нём были определения элементов и последовательность кольца.

ЧТО ПРОВЕРЯЕТ
  • Σ ANGLE по кольцу — замыкание по УГЛУ (должно быть ровно 2π).
  • Замыкание по ПОЛОЖЕНИЮ — прогон всей линии, ВКЛЮЧАЯ ДРЕЙФЫ, и зазор конец→начало.
    Это разные условия: угол можно свести к 2π, а кольцо всё равно не сойдётся в точку.
    NB: мерить «выход последнего магнита → вход первого» БЕСПОЛЕЗНО — на стыке обычно стоит
    законный длинный дрейф, и такой замер даёт метры мнимого зазора.
  • δ = E1+E2−ANGLE у изгибающих — отклонение полюса от rectangular bend. Должна быть одинаковой
    у элементов одного типа; иначе сборка БД падает (build_matrix_db.add_edge_delta).

ЧЕМ ЗАМЫКАЕТ
  Масштабирует ANGLE (и вместе с ним E1/E2) ТОЛЬКО чистых диполей — тех, у кого нет K1, — на
  close_scale = (2π − Σ ANGLE у combined-function) / Σ ANGLE у чистых. Combined-function (K1≠0)
  не трогает: их градиент задаёт оптику, менять их изгиб опаснее.
  Края E1/E2 обязаны ехать тем же множителем: δ = E1+E2−ANGLE должна остаться КОНСТАНТОЙ типа.
  Сдвинув только ANGLE, мы меняли бы δ у каждого среза на свою величину. Общий множитель даёт
  δ → s·δ — по-прежнему константа (s глобальный), и внешние грани группы срезов остаются
  согласованы с box = (Σθ+Σδ)/2.

ЗНАК. Конвенция одна на весь проект: положительный ANGLE = поворот ВПРАВО. Съём геометрии берём
из lte_to_lattice.step, чтобы источник конвенции был единственным.

ЗАПУСК
  py -3.14 lte_ring.py path/to/ring.lte                        # проверка (один файл с LINE=)
  py -3.14 lte_ring.py path/to/desc.lte --line path/to/line.lte    # проверка (линия отдельно)
  py -3.14 lte_ring.py path/to/desc.lte --line path/to/line.lte --out path/to/fixed_desc.lte
  py -3.14 lte_ring.py ring.lte --beamline RING --out fixed.lte    # выбрать линию, если их несколько
Без --out скрипт ничего не пишет — только отчёт.
"""
import os, re, math, shutil, argparse

import lte_to_lattice as _geom          # ЕДИНСТВЕННЫЙ источник конвенции знака (step)

BEND_TYPES = ("CSBEND", "SBEND", "RBEN", "RBEND", "BEND", "KSBEND", "CCBEND")


# --------------------------------------------------------------------------- разбор
def strip_comments(text):
    """Убрать комментарии (!) и склеить продолжения строк (&)."""
    text = re.sub(r"!.*", "", text)
    return re.sub(r"&\s*\n\s*", " ", text)


def parse_elements(text):
    """{ИМЯ: dict(typ, angle, has_k1, params)} — определения элементов (не LINE).

    ВНИМАНИЕ на пробелы: между именем/типом/телом только [ \\t], НЕ \\s — иначе \\s* проглатывает
    перевод строки, и у элемента без параметров («BPM: MONI») телом становится СЛЕДУЮЩАЯ строка:
    BPM получал чужую длину, а съеденное определение исчезало из словаря."""
    out = {}
    for m in re.finditer(r'^[ \t]*"?([A-Za-z0-9_.$]+)"?[ \t]*:[ \t]*([A-Za-z_]+)[ \t]*,?([^\n]*)$',
                         strip_comments(text), re.M):
        name, typ, body = m.group(1).upper(), m.group(2).upper(), m.group(3)
        if typ == "LINE":
            continue
        def num(k):
            mm = re.search(rf'\b{k}\s*=\s*("?)(-?[\d.eE+-]+)\1', body)
            return float(mm.group(2)) if mm else None
        a, k1 = num("ANGLE"), num("K1")
        out[name] = dict(typ=typ, angle=a or 0.0, L=num("L") or 0.0,
                         e1=num("E1"), e2=num("E2"),
                         has_k1=bool(k1) and abs(k1) > 0.0)
    return out


def parse_beamlines(text):
    """{ИМЯ: [токен, ...]} — определения LINE=(...) (токены могут быть вложенными линиями)."""
    out = {}
    for m in re.finditer(r'^[ \t]*"?([A-Za-z0-9_.$]+)"?[ \t]*:[ \t]*LINE[ \t]*=[ \t]*\(([^)]*)\)',
                         strip_comments(text), re.M | re.I):
        toks = [t.strip().upper() for t in re.split(r"[,\s]+", m.group(2)) if t.strip()]
        out[m.group(1).upper()] = toks
    return out


def expand(tokens, lines, depth=0):
    """Развернуть линию в плоский список имён: вложенные линии, кратность n*X, реверс -X."""
    if depth > 50:
        raise SystemExit("слишком глубокая вложенность LINE (циклическая ссылка?)")
    out = []
    for t in tokens:
        rev = t.startswith("-")
        t = t.lstrip("-")
        mult = 1
        m = re.match(r"^(\d+)\s*\*\s*(.+)$", t)
        if m:
            mult, t = int(m.group(1)), m.group(2)
        seq = expand(lines[t], lines, depth + 1) if t in lines else [t]
        if rev:
            seq = seq[::-1]
        out.extend(seq * mult)
    return out


def read_sequence(desc_path, line_path=None, beamline=None):
    """(элементы, последовательность имён). Линия берётся, в порядке приоритета:
    из --line файла; из LINE=(...) в самом файле (--beamline или единственная/самая длинная)."""
    text = open(desc_path, encoding="utf-8", errors="ignore").read()
    elements = parse_elements(text)
    lines = parse_beamlines(text)

    if line_path:                                   # линия отдельным файлом
        ltext = strip_comments(open(line_path, encoding="utf-8", errors="ignore").read())
        elements.update(parse_elements(ltext))      # в файле линии тоже могут быть определения
        llines = parse_beamlines(ltext)
        lines.update(llines)
        if llines:
            name = beamline.upper() if beamline else max(llines, key=lambda k: len(llines[k]))
            return elements, expand(llines[name], lines)
        # голый список токенов (без LINE=) — формат СИЛА
        toks = [t.strip().upper() for t in re.split(r"[,\s]+", ltext) if t.strip()]
        return elements, expand(toks, lines)

    if not lines:
        raise SystemExit(f"в {desc_path} нет LINE=(...); укажи --line с файлом последовательности")
    if beamline:
        if beamline.upper() not in lines:
            raise SystemExit(f"линия '{beamline}' не найдена. Есть: {', '.join(sorted(lines))}")
        name = beamline.upper()
    else:
        name = max(lines, key=lambda k: len(expand(lines[k], lines)))   # самая длинная = кольцо
    return elements, expand(lines[name], lines)


# --------------------------------------------------------------------------- проверки
def bend_sums(elements, seq):
    """(Σ ANGLE чистых диполей, Σ ANGLE combined-function) по всей последовательности."""
    pure = dq = 0.0
    for nm in seq:
        d = elements.get(nm)
        if not d or d["typ"] not in BEND_TYPES or abs(d["angle"]) < 1e-12:
            continue
        if d["has_k1"]:
            dq += d["angle"]
        else:
            pure += d["angle"]
    return pure, dq


def survey(elements, seq):
    """Прогон ВСЕЙ линии (с дрейфами) → (зазор конец→начало в метрах, финальный угол в радианах).
    Конвенция знака — из lte_to_lattice.step (единственный источник)."""
    x = y = th = 0.0
    for nm in seq:
        d = elements.get(nm)
        if not d:
            continue
        _, _, _, x, y, th = _geom.step(x, y, th, d["L"], d["angle"])
    return math.hypot(x, y), th


def delta_report(elements, seq):
    """{(тип, длина_мм): {δ}} — δ=E1+E2−ANGLE по изгибающим.

    Группируем по ТИПУ И ДЛИНЕ, а не по одному типу: в .lte все изгибы обычно одного типа (CSBEND),
    но физически это разные магниты, и δ у них законно разная. Устройство в БД тоже именуется с
    длиной (DP_L_PM_L418, DQ_I1_1A_L1000) — значит инвариант «δ константна» проверяем в тех же
    группах, в каких его позже потребует build_matrix_db.add_edge_delta."""
    by_grp = {}
    for nm in set(seq):
        d = elements.get(nm)
        if not d or d["typ"] not in BEND_TYPES or abs(d["angle"]) < 1e-12:
            continue
        if d["e1"] is None and d["e2"] is None:
            continue
        dl = (d["e1"] or 0.0) + (d["e2"] or 0.0) - d["angle"]
        by_grp.setdefault((d["typ"], round(d["L"] * 1000)), set()).add(round(dl, 12))
    return by_grp


def report(elements, seq, tol):
    """Печать отчёта. Возврат (pure, dq, замкнуто_ли)."""
    pure, dq = bend_sums(elements, seq)
    total = pure + dq
    gap, th_end = survey(elements, seq)
    nb = sum(1 for nm in seq if elements.get(nm, {}).get("typ") in BEND_TYPES
             and abs(elements[nm]["angle"]) > 1e-12)
    print(f"  элементов в линии: {len(seq)}  (определений: {len(elements)}, изгибающих: {nb})")
    print(f"  Σ ANGLE = {total:.9f} рад = {total/math.pi:.9f}·π   "
          f"[чистые {pure:.6f} + combined {dq:.6f}]")
    print(f"    недобор до 2π: {(2*math.pi - total)*1e3:+.3f} мрад")
    print(f"  замыкание по ПОЛОЖЕНИЮ (вся линия с дрейфами): зазор {gap*1e3:.3f} мм")
    print(f"  финальный угол: {(th_end + 2*math.pi)*1e6:+.3f} мкрад от −2π")
    bad = {g: v for g, v in delta_report(elements, seq).items() if len(v) > 1}
    if bad:
        for (typ, lmm), vals in sorted(bad.items()):
            print(f"  ✗ δ=E1+E2−ANGLE у {typ} L={lmm}мм: {len(vals)} разных значений "
                  f"— сборка БД упадёт (add_edge_delta)")
    else:
        print("  δ=E1+E2−ANGLE: константна внутри каждой группы (тип+длина) ✓")
    closed = abs(total - 2 * math.pi) < tol
    print(f"  ИТОГ: кольцо {'ЗАМКНУТО' if closed else 'НЕ замкнуто'} по углу")
    return pure, dq, closed


# --------------------------------------------------------------------------- замыкание
def close_text(text, elements, scale):
    """Переписать текст .lte: у чистых диполей ANGLE, E1, E2 ×= scale. Возврат (текст, сколько)."""
    pure_names = {nm for nm, d in elements.items()
                  if d["typ"] in BEND_TYPES and abs(d["angle"]) > 1e-12 and not d["has_k1"]}
    n = 0
    out = []
    for raw in text.splitlines():
        m = re.match(r'\s*"?([A-Za-z0-9_.$]+)"?\s*:', raw)
        if m and m.group(1).upper() in pure_names:
            raw = re.sub(r'(\b(?:ANGLE|E1|E2)\s*=\s*)([-\d.eE+]+)',
                         lambda mm: mm.group(1) + repr(float(mm.group(2)) * scale), raw)
            n += 1
        out.append(raw)
    return "\n".join(out) + "\n", n


def main():
    ap = argparse.ArgumentParser(
        description="проверка и замыкание кольца в .lte (универсально, по пути к файлу)")
    ap.add_argument("lte", help="путь к .lte с определениями элементов")
    ap.add_argument("--line", default=None,
                    help="путь к файлу с последовательностью кольца (если она не в основном файле)")
    ap.add_argument("--beamline", default=None,
                    help="имя линии, если в файле несколько LINE=(...) (по умолч. самая длинная)")
    ap.add_argument("--out", default=None,
                    help="куда писать ИСПРАВЛЕННЫЙ файл (путь к .lte или папка). Без него — только проверка")
    ap.add_argument("--tol", type=float, default=1e-9, help="допуск замыкания по углу, рад")
    a = ap.parse_args()

    if not os.path.isfile(a.lte):
        raise SystemExit(f"файл не найден: {a.lte}")
    elements, seq = read_sequence(a.lte, a.line, a.beamline)

    print(f"ИСХОДНЫЙ: {a.lte}" + (f"  + линия {a.line}" if a.line else ""))
    pure, dq, closed = report(elements, seq, a.tol)

    if not a.out:
        return
    if closed:
        print("\nКольцо уже замкнуто — исправлять нечего, файл не записан.")
        return
    if abs(pure) < 1e-12:
        raise SystemExit("\nнет чистых диполей (все изгибы с K1) — замыкать нечем; "
                         "масштабировать combined-function небезопасно, правь .lte вручную")

    scale = (2 * math.pi - dq) / pure
    print(f"\nclose_scale (чистые диполи) = {scale:.10f}")

    out_desc = (os.path.join(a.out, os.path.basename(a.lte))
                if os.path.isdir(a.out) or a.out.endswith(("/", "\\")) else a.out)
    os.makedirs(os.path.dirname(os.path.abspath(out_desc)), exist_ok=True)
    text = open(a.lte, encoding="utf-8", errors="ignore").read()
    new_text, n = close_text(text, elements, scale)
    open(out_desc, "w", encoding="utf-8", newline="\n").write(new_text)
    print(f"Записан: {out_desc}  (масштабировано определений: {n})")

    out_line = None
    if a.line:                                        # линия не менялась — копируем рядом
        out_line = os.path.join(os.path.dirname(os.path.abspath(out_desc)), os.path.basename(a.line))
        if os.path.abspath(out_line) != os.path.abspath(a.line):
            shutil.copy(a.line, out_line)
            print(f"Скопирован: {out_line}")

    print("\nПРОВЕРКА записанного:")
    el2, seq2 = read_sequence(out_desc, out_line or a.line, a.beamline)
    report(el2, seq2, a.tol)


if __name__ == "__main__":
    main()
