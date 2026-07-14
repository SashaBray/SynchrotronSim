"""Сброс магнитов к дефолту (untuned): все arg_val -> 0, все correct_coeff -> 1.
Нужен для прогонов "из нуля". Запуск:  python reset_magnets.py
"""
import csv
from pathlib import Path

sp = Path(__file__).resolve().parent / "workspace" / "tables" / "sp_elements.csv"
rows = list(csv.reader(open(sp, encoding="utf-8"), delimiter=";"))
h = rows[0]
a = h.index("arg_val")
c = h.index("correct_coeff")
na = nc = 0
for r in rows[1:]:
    try:
        float(r[a]); r[a] = "0.0000000000"; na += 1
    except ValueError:
        pass  # '-' и прочие нечисловые ячейки пропускаем
    try:
        float(r[c]); r[c] = "1.0000000000"; nc += 1
    except ValueError:
        pass
with open(sp, "w", encoding="utf-8", newline="") as f:
    csv.writer(f, delimiter=";", lineterminator="\n").writerows(rows)
print(f"reset to default: arg_val=0 в {na} строках, correct_coeff=1 в {nc} строках")
