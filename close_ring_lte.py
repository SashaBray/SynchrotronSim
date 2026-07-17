# -*- coding: utf-8 -*-
"""Замкнуть кольцо В ИСТОЧНИКЕ (.lte), а не хаком в парсере.

ПРОБЛЕМА. Углы диполей в исходном .lte суммируются на Σ=1.9866·π (кольцо недогнуто на
~42 мрад = 0.67%). Раньше замыкание «размазывалось» скрытно и ТОЛЬКО в геометрию
(lte_to_lattice/lte_to_design_orbit масштабируют ANGLE до 2π), а кики (warmstart→mat_B)
читали СЫРОЙ угол — отсюда рассогласование кадры(2π)↔кики(1.9866π) и уход опорной ~2.4мм.

РЕШЕНИЕ. Правим ANGLE прямо в .lte и пишем ИСПРАВЛЕННУЮ КОПИЮ конфигурации. Тогда ОБА
парсера (геометрия и warmstart) читают один и тот же замкнутый файл — один источник правды,
рассогласования нет. Closure-хак в lte_to_lattice на замкнутом файле обнуляется (scale→1.0).

МЕТОД (как в историческом хаке): масштабируем ANGLE ТОЛЬКО чистых диполей (без K1) на
close_scale = (2π − Σ ANGLE_DQ) / Σ ANGLE_pure; комбинированные DQ (K1≠0) не трогаем.
Так геометрия кадров не меняется (та же 2π), но кики теперь ей соответствуют.

ВНИМАНИЕ. Замкнутое кольцо с текущими квадруполями (mat_* из warmstart под недогнутую
решётку) НЕустойчиво — нужен рематч тюна (шаг 2). См. скан close_kicks.

Запуск:
  py -3.14 close_ring_lte.py                         # configuration_v1 → configuration_v1_closed
  py -3.14 close_ring_lte.py --src configuration_v1 --out configuration_v1_closed
Затем:
  py -3.14 load_config.py --config v1_closed --force  # собрать БД из замкнутого .lte
"""
import os, re, math, shutil, argparse, glob

ROOT = os.path.dirname(os.path.abspath(__file__))
LTE_ROOT = os.path.join(ROOT, "workspace", "lte")


def find_desc_line(folder):
    descs = glob.glob(os.path.join(folder, "*desc*.lte"))
    lines = glob.glob(os.path.join(folder, "*line*.lte"))
    if not descs or not lines:
        raise SystemExit(f"в {folder} нет *desc*.lte / *line*.lte")
    return descs[0], lines[0]


def parse_desc(path):
    """{NAME: dict(raw=строка, angle, has_k1)} — по одному определению на строку."""
    dev = {}
    for raw in open(path, encoding="utf-8", errors="ignore").read().splitlines():
        m = re.match(r'\s*"?([A-Za-z0-9_.]+)"?\s*:\s*([A-Za-z]+)\b', raw)
        if not m:
            continue
        name, typ = m.group(1).upper(), m.group(2).upper()
        angle_m = re.search(r'\bANGLE\s*=\s*([-\d.eE+]+)', raw)
        k1_m = re.search(r'\bK1\s*=\s*([-\d.eE+]+)', raw)
        dev[name] = dict(
            raw=raw, typ=typ,
            angle=float(angle_m.group(1)) if angle_m else 0.0,
            has_k1=bool(k1_m) and abs(float(k1_m.group(1))) > 0.0,
        )
    return dev


def ring_bend_sums(dev, line_names):
    """Σ ANGLE по всему кольцу, раздельно: чистые диполи (без K1) и DQ (с K1)."""
    pure = dq = 0.0
    for nm in line_names:
        d = dev.get(nm)
        if not d or d["typ"] != "CSBEND" or abs(d["angle"]) < 1e-12:
            continue
        if d["has_k1"]:
            dq += d["angle"]
        else:
            pure += d["angle"]
    return pure, dq


def main():
    ap = argparse.ArgumentParser(description="замкнуть кольцо правкой ANGLE в .lte (исправленная копия конфигурации)")
    ap.add_argument("--src", default="configuration_v1", help="исходная папка-конфиг в workspace/lte/")
    ap.add_argument("--out", default=None, help="папка-конфиг назначения (по умолч. <src>_closed)")
    a = ap.parse_args()

    src_dir = os.path.join(LTE_ROOT, a.src)
    out_name = a.out or (a.src + "_closed")
    out_dir = os.path.join(LTE_ROOT, out_name)
    desc_path, line_path = find_desc_line(src_dir)

    dev = parse_desc(desc_path)
    line_names = [t.strip().upper() for t in re.split(r"[,\s]+",
                  open(line_path, encoding="utf-8", errors="ignore").read()) if t.strip()]

    raw_pure, raw_dq = ring_bend_sums(dev, line_names)
    raw_total = raw_pure + raw_dq
    print(f"Исходный .lte: {a.src}")
    print(f"  Σ ANGLE по кольцу = {raw_total:.6f} рад = {raw_total/math.pi:.6f}·π  "
          f"(недобор {(raw_total-2*math.pi)*1e3:+.3f} мрад)")
    print(f"    чистые диполи Σ = {raw_pure:.6f} рад,  DQ (K1≠0) Σ = {raw_dq:.6f} рад")

    if abs(raw_total - 2 * math.pi) < 1e-9:
        print("  Кольцо уже замкнуто (Σ=2π) — правка не нужна.")
        return

    close_scale = (2 * math.pi - raw_dq) / raw_pure
    print(f"  close_scale (только чистые диполи) = {close_scale:.10f}")

    # какие ИМЕНА устройств — чистые диполи (масштабируем их ANGLE в тексте)
    pure_names = {nm for nm, d in dev.items()
                  if d["typ"] == "CSBEND" and abs(d["angle"]) > 1e-12 and not d["has_k1"]}

    # Переписываем desc построчно: у чистых диполей ANGLE ×= close_scale, и КРАЯ E1/E2 ×= close_scale.
    # ЗАЧЕМ края. Масштабировать ANGLE, не трогая E, нельзя: δ = E1+E2−ANGLE (отклонение полюса от
    # rectangular bend) обязана быть КОНСТАНТОЙ типа устройства — build_matrix_db.add_edge_delta это
    # проверяет и падает при неоднозначности, а add_edge_angles сверяет внешние грани группы срезов с
    # box = (Σθ+Σδ)/2. Меняя только ANGLE, мы сдвигали бы δ у каждого среза на свою величину.
    # Общий множитель обе инварианты сохраняет: δ → s·δ (по-прежнему константа типа, s глобальный),
    # box → s·box, внешние грани → s·(грани) ⇒ сверка сходится. Геометрия дуги масштабируется целиком.
    n_scaled = 0
    out_lines = []
    for raw in open(desc_path, encoding="utf-8", errors="ignore").read().splitlines():
        m = re.match(r'\s*"?([A-Za-z0-9_.]+)"?\s*:', raw)
        if m and m.group(1).upper() in pure_names:
            def repl(mm):
                return mm.group(1) + repr(float(mm.group(2)) * close_scale)
            new = re.sub(r'(\b(?:ANGLE|E1|E2)\s*=\s*)([-\d.eE+]+)', repl, raw)
            out_lines.append(new)
            n_scaled += 1
        else:
            out_lines.append(raw)

    os.makedirs(out_dir, exist_ok=True)
    out_desc = os.path.join(out_dir, os.path.basename(desc_path))
    out_line = os.path.join(out_dir, os.path.basename(line_path))
    with open(out_desc, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out_lines) + "\n")
    shutil.copy(line_path, out_line)
    # прочие файлы конфига (если есть) — копируем как есть
    for p in glob.glob(os.path.join(src_dir, "*")):
        b = os.path.basename(p)
        if os.path.isfile(p) and b not in (os.path.basename(desc_path), os.path.basename(line_path)):
            shutil.copy(p, os.path.join(out_dir, b))

    # верификация: пересчитать по замкнутому файлу
    dev2 = parse_desc(out_desc)
    p2, d2 = ring_bend_sums(dev2, line_names)
    tot2 = p2 + d2
    print(f"\nЗаписано: {os.path.relpath(out_dir, ROOT)}  (изменён ANGLE у {n_scaled} определений чистых диполей)")
    print(f"  Σ ANGLE замкнутого = {tot2:.6f} рад = {tot2/math.pi:.9f}·π  "
          f"(невязка {(tot2-2*math.pi)*1e6:+.3f} мкрад)")
    print(f"\nДалее:  py -3.14 load_config.py --config {out_name.replace('configuration_','')} --force")
    print("        (оба парсера — геометрия и warmstart — прочтут замкнутый .lte; хак самонейтрализуется)")
    print("        ЗАТЕМ рематч квадруполей — замкнутое кольцо неустойчиво со старым тюном.")


if __name__ == "__main__":
    main()
