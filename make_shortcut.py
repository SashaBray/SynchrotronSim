# -*- coding: utf-8 -*-
"""Ярлык запуска графической оболочки (Windows .lnk).

Ярлык хранит АБСОЛЮТНЫЕ пути и потому машинно-зависим (в git не входит) — этот скрипт позволяет
пересоздать его в один вызов на любой машине:

    py -3.14 make_shortcut.py                     # ярлык в корне проекта
    py -3.14 make_shortcut.py --desktop           # ещё и на рабочий стол
    py -3.14 make_shortcut.py --console           # запуск с консольным окном (видно ранние ошибки)
    py -3.14 make_shortcut.py --name "СИЛА"       # своё имя ярлыка

По умолчанию ярлык запускает `pythonw.exe gui_app.py` — без консольного окна. Иконка рисуется в
`gui/synchrotronsim.ico`, если её ещё нет.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
ICON = ROOT / "gui" / "synchrotronsim.ico"
SIZES = (16, 24, 32, 48, 64, 128, 256)


def make_icon(path: Path = ICON, force: bool = False) -> Path:
    """Нарисовать иконку: кольцо синхротрона с ярким сгустком пучка. Требует PIL (Pillow)."""
    if path.exists() and not force:
        return path
    from PIL import Image, ImageDraw

    S = 1024                                                  # рисуем крупно, потом уменьшаем — сглаживание
    img = Image.new("RGBA", (S, S), (0, 0, 0, 0))
    d = ImageDraw.Draw(img)
    pad, ring_w = int(S * 0.10), int(S * 0.11)
    box = (pad, pad, S - pad, S - pad)
    d.ellipse(box, outline=(26, 95, 180, 255), width=ring_w)                    # вакуумная камера
    inner = int(ring_w * 0.34)
    d.arc(box, start=-58, end=32, fill=(245, 166, 35, 255), width=ring_w)       # сгусток пучка на орбите
    d.arc(box, start=-58, end=32, fill=(255, 235, 160, 255), width=inner)       # яркая сердцевина
    r = int(S * 0.055)                                                          # головка сгустка
    cx = S / 2 + (S / 2 - pad - ring_w / 2) * 0.9945
    cy = S / 2 - (S / 2 - pad - ring_w / 2) * 0.104
    d.ellipse((cx - r, cy - r, cx + r, cy + r), fill=(255, 248, 214, 255))

    path.parent.mkdir(parents=True, exist_ok=True)
    img.save(path, format="ICO", sizes=[(s, s) for s in SIZES])
    return path


def pythonw() -> str:
    """`pythonw.exe` рядом с текущим интерпретатором (запуск без консоли); нет — обычный python."""
    p = Path(sys.executable).with_name("pythonw.exe")
    return str(p if p.exists() else Path(sys.executable))


def create_lnk(link: Path, target: str, args: str, workdir: Path, icon: Path, note: str) -> Path:
    """Создать .lnk через WScript.Shell (штатный способ Windows, без сторонних пакетов)."""
    ps = f"""
$s = (New-Object -ComObject WScript.Shell).CreateShortcut('{link}')
$s.TargetPath = '{target}'
$s.Arguments = '{args}'
$s.WorkingDirectory = '{workdir}'
$s.Description = '{note}'
$s.IconLocation = '{icon},0'
$s.Save()
"""
    r = subprocess.run(["powershell", "-NoProfile", "-NonInteractive", "-Command", ps],
                       capture_output=True, text=True, encoding="utf-8", errors="replace")
    if r.returncode != 0 or not link.exists():
        raise SystemExit(f"не удалось создать ярлык:\n{r.stdout}\n{r.stderr}")
    return link


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Ярлык запуска графической оболочки SynchrotronSim")
    ap.add_argument("--name", default="SynchrotronSim", help="имя ярлыка (без .lnk)")
    ap.add_argument("--desktop", action="store_true", help="создать ещё и на рабочем столе")
    ap.add_argument("--console", action="store_true",
                    help="запускать с консольным окном (видно вывод самой оболочки)")
    ap.add_argument("--icon-force", action="store_true", help="перерисовать иконку, даже если она есть")
    a = ap.parse_args(argv)

    if sys.platform != "win32":
        print("ярлыки .lnk — только для Windows; на других системах запускайте py gui_app.py")
        return 1

    try:
        icon = make_icon(force=a.icon_force)
    except ImportError:
        icon = Path(sys.executable)                            # без Pillow — иконка интерпретатора
        print("Pillow не установлен — беру иконку интерпретатора")

    target = str(Path(sys.executable)) if a.console else pythonw()
    args = f'"{ROOT / "gui_app.py"}"'
    note = "Графическая оболочка SynchrotronSim: проекты, база, матрицы связи, трекинг, визуализация"

    made = [create_lnk(ROOT / f"{a.name}.lnk", target, args, ROOT, icon, note)]
    if a.desktop:
        desktop = Path(os.path.expanduser("~")) / "Desktop"
        if desktop.is_dir():
            made.append(create_lnk(desktop / f"{a.name}.lnk", target, args, ROOT, icon, note))
        else:
            print(f"папка рабочего стола не найдена: {desktop}")

    for p in made:
        print(f"ярлык создан: {p}")
    print(f"запускает: {target} {args}")
    print(f"иконка:    {icon}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
