# -*- coding: utf-8 -*-
"""Конвертация готового отчёта в другой формат силами LibreOffice.

Отчёты проекта собираются в `.docx` (`python-docx` — чистый Python, работает на любой ОС, Microsoft
Word для сборки не нужен). LibreOffice открывает `.docx` напрямую, поэтому конвертация — вопрос
удобства, а не совместимости: `.odt` для правки в LibreOffice, `.pdf` для рассылки.

Делает это сам LibreOffice в пакетном режиме (`soffice --headless --convert-to`), поэтому в проекте
нет второго генератора документов и оформление отчёта не дублируется.

    py -3.14 office_convert.py workspace/results/element_compare/element_compare_report.docx
    py -3.14 office_convert.py отчёт.docx --to pdf
    py -3.14 office_convert.py отчёт.docx --to odt --outdir куда/положить
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

from _common import no_window

FORMATS = ("odt", "pdf", "doc", "rtf", "html")

# Типовые места установки LibreOffice, если её нет в PATH
_CANDIDATES = (
    r"C:\Program Files\LibreOffice\program\soffice.exe",
    r"C:\Program Files (x86)\LibreOffice\program\soffice.exe",
    "/usr/bin/soffice", "/usr/bin/libreoffice", "/usr/local/bin/soffice",
    "/snap/bin/libreoffice", "/opt/libreoffice/program/soffice",
    "/Applications/LibreOffice.app/Contents/MacOS/soffice",
)


def soffice() -> str | None:
    """Путь к исполняемому файлу LibreOffice или None, если она не установлена."""
    env = os.environ.get("SS_SOFFICE")
    if env and os.path.exists(env):
        return env
    for name in ("soffice", "libreoffice"):
        p = shutil.which(name)
        if p:
            return p
    return next((p for p in _CANDIDATES if os.path.exists(p)), None)


def convert(src: str | Path, to: str = "odt", outdir: str | Path | None = None,
            timeout: int = 600) -> Path:
    """Конвертировать документ в формат `to`. Возврат — путь результата.

    Бросает FileNotFoundError, если нет исходного файла или самой LibreOffice, и RuntimeError,
    если конвертация не удалась (текст ошибки — из вывода LibreOffice)."""
    src = Path(src).resolve()
    if not src.exists():
        raise FileNotFoundError(f"нет файла: {src}")
    if to not in FORMATS:
        raise ValueError(f"формат '{to}' не поддержан; есть: {', '.join(FORMATS)}")
    exe = soffice()
    if exe is None:
        raise FileNotFoundError(
            "LibreOffice не найдена. Установите её либо укажите путь в переменной SS_SOFFICE "
            "(например SS_SOFFICE=/usr/bin/soffice). Сам отчёт .docx уже собран и открывается "
            "и в Word, и в LibreOffice — конвертация нужна лишь для другого формата.")
    out = Path(outdir).resolve() if outdir else src.parent
    out.mkdir(parents=True, exist_ok=True)
    r = subprocess.run([exe, "--headless", "--norestore", "--convert-to", to, "--outdir", str(out), str(src)],
                       capture_output=True, text=True, encoding="utf-8", errors="replace",
                       timeout=timeout, **no_window())
    dst = out / (src.stem + "." + to)
    if r.returncode != 0 or not dst.exists():
        raise RuntimeError(f"LibreOffice не смогла конвертировать: {(r.stdout or '') + (r.stderr or '')}".strip())
    return dst


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Конвертация отчёта силами LibreOffice")
    ap.add_argument("src", help="исходный документ (обычно .docx)")
    ap.add_argument("--to", default="odt", choices=FORMATS, help="целевой формат (по умолч. odt)")
    ap.add_argument("--outdir", default=None, help="куда положить (по умолч. рядом с исходным)")
    a = ap.parse_args(argv)
    try:
        dst = convert(a.src, a.to, a.outdir)
    except (FileNotFoundError, ValueError, RuntimeError, subprocess.TimeoutExpired) as e:
        print(f"[!] {e}", file=sys.stderr)
        return 1
    print(f"{Path(a.src).name} → {dst}  ({dst.stat().st_size / 1e6:.1f} МБ)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
