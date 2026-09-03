# -*- coding: utf-8 -*-
"""Word-отчёт: модель насыщения магнита «токи обмоток → параметры пучка».

Пошаговое описание расчёта с нативными формулами Word (`latex2omml`) и расшифровкой ВСЕХ переменных
под каждой формулой — с единицами измерения и указанием типа величины (скаляр / вектор / матрица).
Текст рассчитан на читателя без специальной подготовки: каждый физический или математический термин
раскрывается при первом появлении. В конце — три разобранных примера с измеренными числами
(простая линза QP, диполь-квадруполь DQ, корректор), по которым расчёт можно повторить руками.

    py -3.14 saturation_report.py
    py -3.14 saturation_report.py --out путь/отчёт.docx --to pdf

Источники содержания: `docs/SATURATION_MODEL.md` (формальное изложение),
`workspace/research/2026-09-01_corrector_branch_model.md` (разбор проверки),
`workspace/research/corrector_model/corr_model.json` (числа модели корректора),
`workspace/materials/saturation_fits.csv` (однованные подгонки).
"""
from __future__ import annotations

import argparse
import json
import os

OUT_DEFAULT = os.path.join("workspace", "results", "saturation_model_report.docx")
FIGDIR = os.path.join("workspace", "results", "saturation_report_figs")
MODEL_JSON = os.path.join("workspace", "research", "corrector_model", "corr_model.json")
BH_CSV = os.path.join("workspace", "materials", "steel_1008_bh.csv")


# ─────────────────────────────── оформление ────────────────────────────────

def _formula(doc, latex, lo):
    """Формула отдельным абзацем по центру."""
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    lo.add_formula(p, latex)
    return p


_KIND = (("ДИАГОНАЛЬНАЯ МАТРИЦА", "диагональная матрица"), ("СКАЛЯРНАЯ ФУНКЦИЯ", "скалярная функция"),
         ("ВЕКТОР-ФУНКЦИЯ", "вектор-функция"), ("КОМПЛЕКСНОЕ", "комплексное"),
         ("МАТРИЦА", "матрица"), ("ВЕКТОР", "вектор"), ("СКАЛЯРЫ", "скаляры"))


_SEEN = set()


def _vkey(sym):
    """Ключ переменной: буквы, цифры и штрихи, РЕГИСТР ЗНАЧИМ.
    «H(·)» и «H» — одно; «H'», «R» и «r» — разные."""
    return "".join(ch for ch in sym if ch.isalnum() or ch in "'′")


def _fresh(sym):
    """Новая ли переменная. Составное обозначение («ℓ_b, S_b») ново, если ново хоть одно из имён."""
    keys = [_vkey(s) for s in sym.split(",") if _vkey(s)]
    if not keys or all(k in _SEEN for k in keys):
        return False
    _SEEN.update(keys)
    return True


def _vars(doc, triples, title="Здесь приняты следующие обозначения"):
    """Расшифровка переменных списком: «— символ — пояснение; тип, единицы.»

    triples: (обозначение, «тип, единицы», пояснение)."""
    from docx.shared import Cm, Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    doc.add_paragraph(title + ":")
    for sym, kind, txt in triples:
        for a, b in _KIND:
            kind = kind.replace(a, b)
        p = doc.add_paragraph()
        p.paragraph_format.left_indent = Cm(1.0)
        p.paragraph_format.first_line_indent = Cm(-1.0)
        p.paragraph_format.space_after = Pt(2)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY
        p.add_run("— ")
        p.add_run(sym).italic = True
        p.add_run(f" — {txt}; {kind}.")
    doc.add_paragraph()


def _symtab(doc, rows, caption=None):
    """Компактная формальная расшифровка формулы: таблица «обозначение | смысл | тип | единицы».

    rows: (обозначение, смысл, тип, единицы). Короче и строже списка `_vars` — для отчётов, где
    формулы идут подряд и развёрнутая проза только мешает."""
    return _table(doc, ["Обозн.", "Смысл", "Тип", "Ед."], rows, caption=caption, pt=9)


def _tcap(doc, text):
    from docx.shared import Pt
    _tcap.n = getattr(_tcap, "n", 0) + 1
    p = doc.add_paragraph()
    r = p.add_run(f"Таблица {_tcap.n}. {text}")
    r.bold = True
    r.font.size = Pt(10)
    return p


def _fcap(doc, text):
    from docx.shared import Pt
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    _fcap.n = getattr(_fcap, "n", 0) + 1
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = p.add_run(f"Рисунок {_fcap.n}. {text}")
    r.bold = True
    r.font.size = Pt(10)
    return p


def _table(doc, header, rows, caption=None, pt=9):
    from docx.shared import Pt
    if caption:
        _tcap(doc, caption)
    t = doc.add_table(rows=1, cols=len(header))
    t.style = "Light Grid Accent 1"
    for i, h in enumerate(header):
        t.rows[0].cells[i].text = h
    for row in rows:
        c = t.add_row().cells
        for i, v in enumerate(row):
            c[i].text = str(v)
    for row in t.rows:
        for cell in row.cells:
            for par in cell.paragraphs:
                for run in par.runs:
                    run.font.size = Pt(pt)
    doc.add_paragraph()
    return t


def _pic(doc, path, cm=15.0):
    from docx.shared import Cm
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    doc.add_picture(path, width=Cm(cm))
    doc.paragraphs[-1].alignment = WD_ALIGN_PARAGRAPH.CENTER


def _body_format(doc):
    from docx.shared import Cm
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    for p in doc.paragraphs:
        if not p.text.strip() or p.style.name.startswith("Heading") or p.style.name == "Title":
            continue
        if p.alignment == WD_ALIGN_PARAGRAPH.CENTER:
            continue
        if p.text.startswith("— "):                       # список обозначений — свой висячий отступ
            continue
        p.paragraph_format.first_line_indent = Cm(1.25)
        p.alignment = WD_ALIGN_PARAGRAPH.JUSTIFY


# ─────────────────────────────── рисунки ────────────────────────────────

def _figures():
    """Сгенерировать иллюстрации. Возврат — словарь имя → путь."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    import numpy as np
    import pandas as pd
    os.makedirs(FIGDIR, exist_ok=True)
    out = {}

    bh = pd.read_csv(BH_CSV, sep=";")
    H, B = bh["H_A_per_m"].to_numpy(float), bh["B_T"].to_numpy(float)
    MU0 = 4e-7 * np.pi
    fig, ax = plt.subplots(1, 3, figsize=(13.2, 3.5))
    ax[0].plot(H, B, "o-", ms=4, lw=1.6)
    ax[0].set_xlim(0, 20000); ax[0].set_ylim(0, 2.1)
    ax[0].set_xticks([0, 5000, 10000, 15000, 20000])
    ax[0].set_xlabel("H — напряжённость поля, А/м"); ax[0].set_ylabel("B — индукция, Тл")
    ax[0].set_title("Колено кривой: рабочая область"); ax[0].grid(alpha=.3)
    ax[1].plot(H, B, "o-", ms=4, lw=1.6)
    ax[1].set_xlim(0, H[-1] * 1.02); ax[1].set_ylim(0, 2.8)
    ax[1].set_xticks([0, 100000, 200000, 300000, 400000])
    ax[1].set_xticklabels(["0", "100 000", "200 000", "300 000", "400 000"])
    ax[1].set_xlabel("H — напряжённость поля, А/м"); ax[1].set_ylabel("B — индукция, Тл")
    ax[1].set_title("Вся таблица целиком"); ax[1].grid(alpha=.3)
    ax[2].plot(B[1:], B[1:] / (MU0 * H[1:]), "o-", ms=4, lw=1.6, color="crimson")
    ax[2].set_xlim(0, 2.7); ax[2].set_ylim(0, 2400)
    ax[2].set_xlabel("B — индукция, Тл"); ax[2].set_ylabel("во сколько раз лучше воздуха")
    ax[2].set_title("Относительная проницаемость"); ax[2].grid(alpha=.3)
    fig.tight_layout()
    out["bh"] = os.path.join(FIGDIR, "bh.png")
    fig.savefig(out["bh"], dpi=150); plt.close(fig)

    IR = dict(fc="#c8d4e3", ec="#33475b", lw=1.2)
    CU = dict(fc="#e8a33d", ec="#8a5a12", lw=0.8)
    fig, ax = plt.subplots(figsize=(6.6, 4.6))
    ax.add_patch(Rectangle((-4.6, 2.7), 8.3, 0.85, **IR))
    ax.add_patch(Rectangle((-4.6, -3.55), 8.3, 0.85, **IR))
    ax.add_patch(Rectangle((2.8, -2.7), 0.9, 5.4, **IR))
    for x, br, coil in ((-3.9, "T+", "c1"), (-1.9, "T0", None), (0.1, "T−", "c2")):
        ax.add_patch(Rectangle((x, 1.0), 0.9, 1.7, **IR))
        ax.text(x + .45, 1.85, br, ha="center", va="center", fontsize=9.5, weight="bold")
        if coil:
            ax.add_patch(Rectangle((x - .45, 1.25), .42, 1.2, **CU))
            ax.add_patch(Rectangle((x + .93, 1.25), .42, 1.2, **CU))
            ax.text(x + .45, .72, coil, ha="center", va="center", fontsize=10.5,
                    color="#8a5a12", weight="bold")
    for x, br, coil in ((-3.9, "B+", "c3"), (-1.9, "B0", None), (0.1, "B−", "c4")):
        ax.add_patch(Rectangle((x, -2.7), 0.9, 1.7, **IR))
        ax.text(x + .45, -1.85, br, ha="center", va="center", fontsize=9.5, weight="bold")
        if coil:
            ax.add_patch(Rectangle((x - .45, -2.45), .42, 1.2, **CU))
            ax.add_patch(Rectangle((x + .93, -2.45), .42, 1.2, **CU))
            ax.text(x + .45, -.72, coil, ha="center", va="center", fontsize=10.5,
                    color="#8a5a12", weight="bold")
    ax.add_patch(Rectangle((2.34, -.8), .42, 1.6, **CU))
    ax.add_patch(Rectangle((3.74, -.8), .42, 1.6, **CU))
    ax.text(3.25, 1.25, "main", ha="center", fontsize=10.5, color="#8a5a12", weight="bold")
    ax.text(3.25, -1.35, "ярмо", ha="center", fontsize=9.5, weight="bold")
    ax.plot([-5.1, 2.1], [0, 0], color="#cc2222", lw=1.6, ls="--")
    ax.text(-5.3, 0, "пучок", ha="right", va="center", fontsize=9.5, color="#cc2222")
    ax.annotate("", xy=(-1.45, .98), xytext=(-1.45, -.98),
                arrowprops=dict(arrowstyle="<->", color="#666", lw=1))
    ax.text(-1.30, .55, "зазор", ha="left", va="center", fontsize=9, color="#555")
    ax.set_xlim(-6.4, 4.6); ax.set_ylim(-4.1, 4.1); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout()
    out["corr"] = os.path.join(FIGDIR, "graph_corrector.png")
    fig.savefig(out["corr"], dpi=170); plt.close(fig)

    fig, axs = plt.subplots(1, 3, figsize=(12.6, 3.6))
    # QP — один сердечник, одна катушка
    a = axs[0]
    a.add_patch(Rectangle((-1.5, 1.05), 3.0, .75, **IR)); a.add_patch(Rectangle((-1.5, -1.8), 3.0, .75, **IR))
    a.add_patch(Rectangle((-1.85, -1.05), .8, 2.1, **IR))
    a.add_patch(Rectangle((-.45, .25), .9, .8, **IR)); a.add_patch(Rectangle((-.45, -1.05), .9, .8, **IR))
    for y in (.4, -.95):
        a.add_patch(Rectangle((-.87, y), .4, .55, **CU)); a.add_patch(Rectangle((.47, y), .4, .55, **CU))
    a.text(1.05, .6, "I", fontsize=11, color="#8a5a12", weight="bold")
    a.text(0, 2.25, "QP — одна ветвь, одна катушка", ha="center", fontsize=10, weight="bold")
    a.plot([-2.4, 2.1], [-.1, -.1], color="#cc2222", lw=1.4, ls="--")
    a.set_xlim(-2.7, 2.3); a.set_ylim(-2.5, 2.7); a.set_aspect("equal"); a.axis("off")
    # SP — один сердечник, ДВЕ катушки
    b = axs[1]
    b.add_patch(Rectangle((-1.5, 1.05), 3.0, .75, **IR)); b.add_patch(Rectangle((-1.5, -1.8), 3.0, .75, **IR))
    b.add_patch(Rectangle((-1.85, -1.05), .8, 2.1, **IR))
    b.add_patch(Rectangle((-.45, .25), .9, .8, **IR)); b.add_patch(Rectangle((-.45, -1.05), .9, .8, **IR))
    for y, c in ((.72, "#e8a33d"), (.4, "#c47d1c")):
        b.add_patch(Rectangle((-.87, y), .4, .3, fc=c, ec="#8a5a12", lw=.8))
        b.add_patch(Rectangle((.47, y), .4, .3, fc=c, ec="#8a5a12", lw=.8))
    for y, c in ((-.72, "#e8a33d"), (-1.04, "#c47d1c")):
        b.add_patch(Rectangle((-.87, y), .4, .3, fc=c, ec="#8a5a12", lw=.8))
        b.add_patch(Rectangle((.47, y), .4, .3, fc=c, ec="#8a5a12", lw=.8))
    b.text(1.05, .85, "I1", fontsize=10.5, color="#8a5a12", weight="bold")
    b.text(1.05, .45, "I2", fontsize=10.5, color="#8a5a12", weight="bold")
    b.text(0, 2.25, "SP — одна ветвь, две катушки", ha="center", fontsize=10, weight="bold")
    b.plot([-2.4, 2.1], [-.1, -.1], color="#cc2222", lw=1.4, ls="--")
    b.set_xlim(-2.7, 2.3); b.set_ylim(-2.5, 2.7); b.set_aspect("equal"); b.axis("off")
    # DQ — общий сердечник, на каждом конце раздваивается: I1 на одной ветви, I2+I3 на другой
    c = axs[2]
    c.add_patch(Rectangle((2.05, -1.15), .85, 2.3, **IR))                     # общий сердечник
    c.add_patch(Rectangle((-1.9, 1.15), 4.8, .7, **IR))                       # верхняя развилка
    c.add_patch(Rectangle((-1.9, -1.85), 4.8, .7, **IR))                      # нижняя развилка
    for y0, y1 in ((.3, 1.15), (-1.15, -.3)):
        c.add_patch(Rectangle((-1.55, y0), .8, y1 - y0, **IR))                # ветвь I1
        c.add_patch(Rectangle((.25, y0), .8, y1 - y0, **IR))                  # ветвь I2+I3
    for y, col in ((.72, "#e8a33d"), (.4, "#c47d1c")):                    # I2 + I3 — ЛЕВАЯ ветвь
        c.add_patch(Rectangle((-1.95, y), .38, .28, fc=col, ec="#8a5a12", lw=.8))
        c.add_patch(Rectangle((-.75, y), .38, .28, fc=col, ec="#8a5a12", lw=.8))
    for y, col in ((-.72, "#e8a33d"), (-1.04, "#c47d1c")):
        c.add_patch(Rectangle((-1.95, y), .38, .28, fc=col, ec="#8a5a12", lw=.8))
        c.add_patch(Rectangle((-.75, y), .38, .28, fc=col, ec="#8a5a12", lw=.8))
    for y in (.45, -1.0):                                                 # I1 — ПРАВАЯ ветвь
        c.add_patch(Rectangle((-.15, y), .38, .5, **CU)); c.add_patch(Rectangle((1.05, y), .38, .5, **CU))
    c.text(-1.15, 2.02, "I2 + I3", ha="center", fontsize=10.5, color="#8a5a12", weight="bold")
    c.text(.65, 2.02, "I1", ha="center", fontsize=10.5, color="#8a5a12", weight="bold")
    c.text(2.48, 0, "общий сердечник", ha="center", va="center", rotation=90,
           fontsize=9, weight="bold")
    c.text(0, 2.35, "DQ — общий сердечник, раздвоенный на концах",
           ha="center", fontsize=10, weight="bold")
    c.plot([-2.4, 1.9], [-.1, -.1], color="#cc2222", lw=1.4, ls="--")
    c.set_xlim(-2.7, 3.3); c.set_ylim(-2.5, 2.9); c.set_aspect("equal"); c.axis("off")
    fig.tight_layout()
    out["qpdq"] = os.path.join(FIGDIR, "graph_qp_sp_dq.png")
    fig.savefig(out["qpdq"], dpi=170); plt.close(fig)

    # ── графы, иллюстрирующие три матрицы модели
    from matplotlib.patches import Circle, FancyArrowPatch
    Qm = np.array(json.load(open(MODEL_JSON, encoding="utf-8"))["Q_Wb_per_At"]) * 1e7
    Gm = np.array(json.load(open(MODEL_JSON, encoding="utf-8"))["G"])
    QTY = json.load(open(MODEL_JSON, encoding="utf-8"))["qty"]
    BRN = ["T+", "T0", "T−", "B+", "B0", "B−", "ярмо"]
    POS = {"T+": (-1.3, 1.0), "T0": (0, 1.0), "T−": (1.3, 1.0), "B+": (-1.3, -1.0),
           "B0": (0, -1.0), "B−": (1.3, -1.0), "ярмо": (2.6, 0)}
    def _node(ax, xy, txt, r=.30, fc="#dce6f2", fs=9):
        ax.add_patch(Circle(xy, r, fc=fc, ec="#33475b", lw=1.4, zorder=3))
        ax.text(xy[0], xy[1], txt, ha="center", va="center", fontsize=fs, zorder=4, weight="bold")
    fig, axs = plt.subplots(1, 3, figsize=(13.6, 4.3))
    a = axs[0]
    wy = {"main": 2.0, "c1": 1.0, "c2": 0.0, "c3": -1.0, "c4": -2.0}
    byy = {"ярмо": 2.0, "T+": 1.0, "T0": 0.3, "T−": -0.4, "B+": -1.1, "B0": -1.8, "B−": -2.5}
    for w, y in wy.items(): _node(a, (0, y), w, fc="#f6d9a8")
    for br, y in byy.items(): _node(a, (2.2, y), br)
    for w, br, nt in (("main", "ярмо", 622), ("c1", "T+", 548), ("c2", "T−", 548),
                      ("c3", "B+", 548), ("c4", "B−", 548)):
        a.annotate("", xy=(1.9, byy[br]), xytext=(.3, wy[w]),
                   arrowprops=dict(arrowstyle="-|>", color="#8a5a12", lw=1.6))
        a.text(1.1, (wy[w] + byy[br]) / 2 + .1, f"{nt}", fontsize=8, color="#8a5a12", ha="center")
    a.text(0, 2.75, "токи I", ha="center", fontsize=10, weight="bold")
    a.text(2.2, 2.75, "ветви", ha="center", fontsize=10, weight="bold")
    a.set_title("Матрица T — размещение обмоток", fontsize=10)
    a.text(1.1, 3.18, "рёбра подписаны числом витков", ha="center", fontsize=8.5, color="#555")
    a.set_xlim(-.8, 3.0); a.set_ylim(-3.1, 3.5); a.axis("off"); a.set_aspect("equal")
    b = axs[1]
    for nm, xy in POS.items(): _node(b, xy, nm)
    for i in range(6):
        for j in range(i + 1, 7):
            v = Qm[i, j]
            if abs(v) < 0.5: continue
            pp, qq = POS[BRN[i]], POS[BRN[j]]
            col = "#1f6fb4" if v > 0 else "#c0392b"
            rad = 0.0 if abs(pp[1] - qq[1]) < 1e-9 or abs(pp[0] - qq[0]) < 1e-9 else 0.18
            b.add_patch(FancyArrowPatch(pp, qq, arrowstyle="-", color=col,
                                        lw=0.6 + 1.1 * min(abs(v), 5.5) / 2.2, alpha=.85,
                                        connectionstyle=f"arc3,rad={rad}", zorder=1,
                                        shrinkA=17, shrinkB=17))
            t = 0.30 if BRN[j] == "ярмо" else 0.5
            mxx, myy = pp[0] + t * (qq[0] - pp[0]), pp[1] + t * (qq[1] - pp[1])
            dx = .30 if abs(pp[0] - qq[0]) < 1e-9 else 0.0
            b.text(mxx + dx, myy + .13, f"{v:+.1f}", fontsize=7.5, color=col, ha="center", zorder=5,
                   bbox=dict(fc="white", ec="none", alpha=.8, pad=.6))
    b.plot([-2.0, 2.0], [0, 0], color="#cc2222", lw=1.2, ls="--", zorder=0)
    b.text(-2.15, 0, "пучок", ha="right", va="center", fontsize=8.5, color="#cc2222")
    b.set_title("Матрица Q — проводимости ветвей", fontsize=10)
    b.text(0.2, 1.75, "синее: потоки складываются, красное: встречны", ha="center",
           fontsize=8.5, color="#555")
    b.set_xlim(-2.9, 3.4); b.set_ylim(-2.0, 2.0); b.axis("off"); b.set_aspect("equal")
    c = axs[2]
    qy = {q: 1.5 - i for i, q in enumerate(QTY)}
    gby = {"T+": 2.5, "T0": 1.5, "T−": 0.5, "B+": -0.5, "B0": -1.5, "B−": -2.5}
    for nm, y in gby.items(): _node(c, (0, y), nm)
    for q, y in qy.items(): _node(c, (2.4, y), q, fc="#cfe8cf")
    cmax = np.abs(Gm).max(axis=0)
    for i, nm in enumerate(gby):
        for j, q in enumerate(QTY):
            v = Gm[i, j]
            if abs(v) < 0.5 * cmax[j]: continue
            c.annotate("", xy=(2.1, qy[q]), xytext=(.3, gby[nm]),
                       arrowprops=dict(arrowstyle="-|>", color="#1f6fb4" if v > 0 else "#c0392b",
                                       lw=0.5 + 2.0 * abs(v) / cmax[j], alpha=.8))
    c.text(0, 3.05, "ветви", ha="center", fontsize=10, weight="bold")
    c.text(2.4, 3.05, "коэффициенты y", ha="center", fontsize=10, weight="bold")
    c.set_title("Матрица G — вклады ветвей в поле", fontsize=10)
    c.text(1.2, 3.62, "на каждый коэффициент — его сильнейшие ветви", ha="center",
           fontsize=8.5, color="#555")
    c.set_xlim(-.8, 3.2); c.set_ylim(-3.1, 3.9); c.axis("off"); c.set_aspect("equal")
    fig.tight_layout()
    out["graphs"] = os.path.join(FIGDIR, "graphs.png")
    fig.savefig(out["graphs"], dpi=165); plt.close(fig)

    # ── графы простых элементов: цепочка «токи → ветви → коэффициенты»
    XW, XB, XC = 0.0, 2.6, 5.2

    def _chain(path, title, size, wnd, brs, coef, tedges, qedges, gedges, ylim):
        fg, ax = plt.subplots(figsize=size)
        for nm, y in wnd.items(): _node(ax, (XW, y), nm, r=.40, fc="#f6d9a8", fs=9)
        for nm, y in brs.items(): _node(ax, (XB, y), nm, r=.46, fs=8.5)
        for nm, y in coef.items(): _node(ax, (XC, y), nm, r=.42, fc="#cfe8cf", fs=8.5)
        for w, br, lab in tedges:
            ax.annotate("", xy=(XB - .55, brs[br]), xytext=(XW + .45, wnd[w]),
                        arrowprops=dict(arrowstyle="-|>", color="#8a5a12", lw=1.6))
            ax.text((XW + XB) / 2, (wnd[w] + brs[br]) / 2 + .18, lab, fontsize=8.5,
                    color="#8a5a12", ha="center")
        for b1, b2 in qedges:
            ax.add_patch(FancyArrowPatch((XB, brs[b1]), (XB, brs[b2]), arrowstyle="-",
                                         color="#1f6fb4", lw=2.2, alpha=.85, shrinkA=26, shrinkB=26,
                                         connectionstyle="arc3,rad=0.40", zorder=1))
        for br, q in gedges:
            ax.annotate("", xy=(XC - .5, coef[q]), xytext=(XB + .52, brs[br]),
                        arrowprops=dict(arrowstyle="-|>", color="#1f6fb4", lw=1.3, alpha=.75))
        top = ylim[1] - .30
        for x, t in ((XW, "токи I"), (XB, "ветви"), (XC, "коэффициенты y")):
            ax.text(x, top, t, ha="center", fontsize=9.5, weight="bold")
        bot = ylim[0] + .22
        for x, t in (((XW + XB) / 2, "матрица T"), (XB, "матрица Q"), ((XB + XC) / 2, "матрица G")):
            ax.text(x, bot, t, ha="center", fontsize=8.5, color="#555", style="italic")
        ax.set_title(title, fontsize=10)
        ax.set_xlim(XW - 1.1, XC + 1.1); ax.set_ylim(*ylim); ax.axis("off"); ax.set_aspect("equal")
        fg.tight_layout(); fg.savefig(path, dpi=170); plt.close(fg)

    out["gqp"] = os.path.join(FIGDIR, "graph_qp_chain.png")
    _chain(out["gqp"], "Квадруполь: одна обмотка, одна ветвь, один коэффициент", (8.0, 2.9),
           {"I": 0.0}, {"ветвь": 0.0}, {"∫G dl": 0.0},
           [("I", "ветвь", "N витков")], [], [("ветвь", "∫G dl")], (-1.5, 1.5))

    out["gdq"] = os.path.join(FIGDIR, "graph_dq_chain.png")
    _chain(out["gdq"], "Диполь-квадруполь: I1 на своей ветви, I2 и I3 — вместе на другой", (8.6, 4.6),
           {"I1": 1.5, "I2": 0.0, "I3": -1.5},
           {"I1": 1.5, "I2 I3": -0.8, "общий": -2.9},
           {"B": 1.5, "G": 0.0, "S": -1.5},
           [("I1", "I1", "N₁"), ("I2", "I2 I3", "N₂"), ("I3", "I2 I3", "N₃")],
           [("I1", "общий"), ("I2 I3", "общий")],
           [(b, q) for b in ("I1", "I2 I3") for q in ("B", "G", "S")], (-4.3, 2.6))
    return out


# ─────────────────────────────── разделы ────────────────────────────────

def _ch_intro(doc, lo, F):
    doc.add_heading("1. Задача, которую решает модель", level=1)
    doc.add_paragraph(
        "Магнит управляет пучком заряженных частиц. Сила и характер этого управления задаются токами "
        "в его обмотках. Пока железо сердечника далеко от предела своих возможностей, связь простая: "
        "вдвое больший ток даёт вдвое более сильное поле. Но по мере роста тока железо намагничивается "
        "всё труднее, и эта пропорциональность нарушается — магнит начинает отставать от тока. Явление "
        "называют насыщением.")
    doc.add_paragraph(
        "Задача модели: по заданным токам обмоток предсказать, как магнит подействует на пучок, с учётом "
        "насыщения. При этом модель должна опираться на разумное число расчётов поля и допускать замену "
        "материала сердечника без пересчёта всего заново.")
    doc.add_paragraph(
        "Сложность появляется у магнитов с несколькими обмотками на общем сердечнике. Их поля нельзя "
        "просто сложить: обмотки соперничают за одно и то же железо, и результат зависит не только от "
        "величин токов, но и от их знаков. Именно этот случай модель и описывает.")

    doc.add_heading("2. Понятия, которыми пользуется модель", level=1)
    doc.add_paragraph(
        "Все дальнейшие рассуждения опираются на десяток понятий. Здесь они собраны вместе и объяснены "
        "обиходными словами; в тексте они больше не расшифровываются.")
    _table(doc, ["Понятие", "Что это означает"], [
        ("Магнитный поток",
         "Мера того, сколько магнитного поля проходит сквозь заданное сечение. Аналог расхода воды в трубе. "
         "Измеряется в веберах (Вб)."),
        ("Ампер-виток",
         "Произведение тока на число витков обмотки. Именно оно, а не ток сам по себе, создаёт магнитное "
         "поле: 1 ампер в обмотке из 500 витков действует так же, как 500 ампер в одном витке."),
        ("Магнитодвижущая сила",
         "То же самое, что ампер-витки, названное по роли: «напор», который гонит магнитный поток по "
         "железу. Прямой аналог напряжения в электрической цепи."),
        ("Магнитная проницаемость",
         "Насколько легко материал пропускает магнитное поле. У железа она в тысячи раз выше, чем у "
         "воздуха, поэтому поле идёт по железу, как вода по трубе, а не растекается вокруг."),
        ("Насыщение",
         "Предел, за которым железо больше не отвечает на рост тока: все его внутренние магнитики уже "
         "выстроены по полю. После этого прибавка поля идёт как в воздухе — в тысячи раз слабее."),
        ("Кривая намагничивания",
         "Таблица, отвечающая на вопрос «какой напор нужен, чтобы получить в этом материале такую "
         "плотность поля». Единственное место, где в модели упоминается материал. Обозначается H(B)."),
        ("Магнитная цепь",
         "Способ считать магнит как электрическую схему: поток играет роль тока, ампер-витки — роль "
         "напряжения, участки железа — роль сопротивлений. Приём позволяет обойтись арифметикой вместо "
         "трёхмерного расчёта."),
        ("Ветвь",
         "Участок магнитопровода, вдоль которого поток можно считать одинаковым: отдельный полюс, "
         "сердечник, кусок ярма. Разбиение магнита на ветви — единственное, что нужно задать о его "
         "устройстве."),
        ("Зазор (апертура)",
         "Промежуток между полюсами, где летит пучок. Здесь поле и работает."),
        ("Мультиполи",
         "Разложение поля в зазоре на простые составляющие: постоянная по сечению часть отклоняет пучок "
         "целиком (дипольная), линейно нарастающая фокусирует его (квадрупольная), квадратичная "
         "исправляет хроматические ошибки (секступольная) и так далее."),
        ("Кик",
         "Изменение направления полёта частицы после прохождения магнита — то, ради чего магнит и "
         "поставлен. Вычисляется из мультиполей."),
        ("Эффективная длина",
         "Длина воображаемого магнита с резкими краями и однородным полем, который действует на пучок "
         "так же, как настоящий с его размытыми краями."),
    ], pt=9)


def _ch_scheme(doc, lo, F):
    doc.add_heading("3. Общая схема расчёта", level=1)
    doc.add_paragraph(
        "Расчёт устроен как цепочка из трёх звеньев. Крайние звенья линейны и измеряются прямо из "
        "расчётов поля; всё нелинейное собрано в среднем звене, и только оно содержит подбираемые числа.")
    _table(doc, ["Звено", "Что делает", "Линейно?", "Откуда берётся"], [
        ("1", "Ампер-витки обмоток → потоки по ветвям, пока железо не мешает",
         "да", "измеряется из расчётов поля"),
        ("2", "Поправка на насыщение: потоки пересчитываются с учётом сопротивления железа",
         "нет", "подбирается: два числа на ветвь"),
        ("3", "Потоки по ветвям → поле в зазоре → мультиполи → кик",
         "да", "измеряется из тех же расчётов поля"),
    ])
    doc.add_paragraph(
        "Такое разделение принципиально. Насыщение не действует на поле в зазоре напрямую: оно "
        "перераспределяет поток между ветвями железа, а поле в зазоре складывается из потоков ветвей уже "
        "по неизменным правилам. Поэтому насыщение меняет не только силу магнита, но и состав его поля — "
        "соотношение отклоняющей, фокусирующей и прочих составляющих.")



def _ch_steps(doc, lo, F):
    doc.add_heading("4. Расчёт шаг за шагом", level=1)
    doc.add_paragraph(
        "Ниже n — число обмоток магнита, m — число ветвей магнитопровода, k — число учитываемых "
        "составляющих поля. Для корректора из примера в приложении n = 5, m = 7, k = 4.")

    doc.add_heading("Шаг 1. От токов к ампер-виткам", level=2)
    doc.add_paragraph(
        "Поле создаёт не ток сам по себе, а ток, умноженный на число витков. Поэтому первое, что делает "
        "расчёт, — переводит амперы в ампер-витки.")
    _formula(doc, r"F_w = N_w \cdot I_w", lo)
    _vars(doc, [
        ("F_w", "скаляр, А·виток", "магнитодвижущая сила обмотки w"),
        ("N_w", "скаляр, штук", "число витков обмотки w"),
        ("I_w", "скаляр, ампер", "ток в обмотке w"),
        ("w", "индекс, 1…n", "номер обмотки"),
    ])
    doc.add_paragraph(
        "Смысл шага чисто пересчётный, но пропускать его нельзя: обмотки одного магнита обычно имеют "
        "разное число витков, и одинаковые токи в них создают разный напор.")

    doc.add_heading("Шаг 2. Чем ветвь сопротивляется потоку", level=2)
    doc.add_paragraph(
        "Поток, идущий по ветви магнитопровода, встречает сопротивление и тратит на него часть напора. "
        "Расход напора складывается из двух слагаемых, устроенных по-разному.")
    _formula(doc, r"U_b = R_b \cdot \Phi_b + \ell_b \cdot H\left(\frac{\Phi_b}{S_b}\r\right)", lo)
    _vars(doc, [
        ("U_b", "скаляр, А·виток", "расход напора на ветви b"),
        ("Φ_b", "скаляр, вебер", "поток в ветви b"),
        ("R_b", "скаляр, А·виток/Вб", "сопротивление воздушной части пути"),
        ("ℓ_b", "скаляр, метр", "длина железа ветви b"),
        ("S_b", "скаляр, м²", "сечение железа ветви b"),
        ("Φ_b / S_b", "скаляр, тесла", "индукция в железе ветви b"),
        ("H(·)", "функция, Тл → А/м", "кривая намагничивания материала"),
        ("b", "индекс, 1…m", "номер ветви"),
    ])
    doc.add_paragraph(
        "Первое слагаемое пропорционально потоку: воздух отвечает на любое воздействие одинаково. Второе "
        "не пропорционально: пока железо не насыщено, оно почти ничего не забирает, а за пределом "
        "насыщения начинает забирать стремительно. Вся нелинейность модели сосредоточена во втором "
        "слагаемом, и она местная: расход на ветви зависит только от потока этой же ветви.")

    doc.add_heading("Шаг 3. Связь напора и потоков, пока железо не мешает", level=2)
    doc.add_paragraph(
        "Представим на минуту, что железо бесконечно податливо и вторым слагаемым можно пренебречь. Тогда "
        "магнит — обычная линейная схема, и потоки во всех ветвях выражаются через приложенные напоры "
        "одной матрицей.")
    _formula(doc, r"\Phi_b = \sum_{c} Q_{bc} \cdot F_c", lo)
    _vars(doc, [
        ("Φ_b", "скаляр, вебер", "поток в ветви b — искомое"),
        ("F_c", "скаляр, А·виток", "напор обмотки на ветви c; нет обмотки — ноль"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей"),
        ("b, c", "индексы, 1…m", "номера ветвей"),
    ])
    doc.add_paragraph(
        "Диагональные элементы матрицы говорят, сколько потока обмотка гонит по своей собственной ветви, "
        "внедиагональные — сколько потока попадает к соседям. Матрица целиком описывает устройство "
        "магнита и от токов не зависит.")

    doc.add_heading("Свойство, которое экономит расчёты", level=3)
    doc.add_paragraph("Матрица проводимостей симметрична: перестановка номеров ветвей ничего не меняет.")
    _formula(doc, r"Q_{bc} = Q_{cb}", lo)
    _vars(doc, [
        ("Q_bc", "скаляр, Вб/(А·виток)", "поток в ветви b от напора в ветви c"),
        ("Q_cb", "скаляр, Вб/(А·виток)", "то же в обратную сторону"),
    ])
    doc.add_paragraph(
        "Причина в том, что энергия, запасённая магнитом, выражается через ту же матрицу, а энергии "
        "безразличен порядок слагаемых:")
    _formula(doc, r"W = \frac{1}{2} \sum_{b} \sum_{c} F_b \cdot Q_{bc} \cdot F_c", lo)
    _vars(doc, [
        ("W", "скаляр, джоуль", "энергия магнитного поля"),
        ("F_b, F_c", "скаляры, А·виток", "напоры в ветвях b и c"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей"),
    ])
    doc.add_paragraph(
        "Практическая польза от симметрии решающая. У магнита есть ветви, на которых обмоток нет — их "
        "не получится возбудить и измерить напрямую. Но поток, который такая ветвь получает от соседней "
        "обмотки, измерить можно, а по симметрии это ровно то же число, что и обратное влияние. Так "
        "столбцы матрицы для ветвей без обмоток достраиваются бесплатно.")
    doc.add_paragraph(
        "Симметрия служит и проверкой правильности разбиения на ветви: если она нарушена сильнее "
        "нескольких процентов, значит ветви размечены неверно. На корректоре проверка дала расхождение "
        "4–5 процентов, то есть разбиение принято.")

    doc.add_heading("Шаг 4. Как учесть железо", level=2)
    doc.add_paragraph(
        "Теперь вернём второе слагаемое из шага 2. Ключевое наблюдение: расход напора на железе входит "
        "в уравнения совершенно так же, как напор обмотки, только с обратным знаком. Для расчёта железо "
        "неотличимо от обмотки, включённой навстречу. Значит матрицу проводимостей можно применять как "
        "есть — надо лишь вычесть из приложенного напора то, что забирает железо.")
    _formula(doc, r"\Phi_b = \sum_{c} Q_{bc} \left( F_c - D_c \r\right)", lo)
    _formula(doc, r"D_b = \ell_b \cdot H\left(\frac{\Phi_b}{S_b}\r\right)", lo)
    _vars(doc, [
        ("Φ", "вектор m, вебер", "потоки ветвей — искомое"),
        ("F", "вектор m, А·виток", "напоры обмоток по ветвям; без обмотки — ноль"),
        ("D", "вектор m, А·виток", "расход напора на железе ветвей"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей"),
        ("ℓ_b, S_b", "скаляры, м и м²", "длина и сечение железа ветви b"),
        ("H(·)", "функция, Тл → А/м", "кривая намагничивания"),
    ])
    doc.add_paragraph(
        "Эти два равенства и есть модель. Особенность в том, что они неявные: чтобы вычислить расход на "
        "железе, нужно знать поток, а чтобы вычислить поток, нужно знать расход. Такие системы решаются "
        "последовательными приближениями.")

    doc.add_heading("Шаг 5. Как решается система", level=2)
    doc.add_paragraph(
        "Начинают с ответа для линейного случая — считают потоки так, будто железо не сопротивляется. По "
        "этим потокам вычисляют расход на железе, вычитают его из напоров, получают уточнённые потоки, и "
        "так далее. Чтобы приближения сходились быстро, шаг выбирают по методу Ньютона — приёму, который "
        "на каждом шаге заменяет кривую касательной прямой и решает уже линейную задачу.")
    _formula(doc, r"\left( R + J \r\right) \cdot \Delta\Phi = -\, r", lo)
    _formula(doc, r"J_{bb} = \frac{\ell_b}{S_b} \cdot H'\left(\frac{\Phi_b}{S_b}\r\right)", lo)
    _vars(doc, [
        ("ΔΦ", "вектор m, вебер", "поправка к потокам за шаг"),
        ("ε", "вектор m, вебер", "невязка; в решении ноль"),
        ("R", "матрица m × m, А·виток/Вб", "обратная к Q"),
        ("J", "диагональная матрица m × m, А·виток/Вб", "крутизна сопротивления железа"),
        ("H'(·)", "функция, (А/м)/Тл", "производная кривой намагничивания"),
        ("ℓ_b, S_b", "скаляры, м и м²", "длина и сечение железа ветви b"),
    ])
    doc.add_paragraph(
        "У такого расчёта есть важное свойство: сумма двух матриц в левой части всегда обратима, потому "
        "что первая положительна, а вторая неотрицательна. Приближения никогда не упираются в вырожденный "
        "случай, при каком угодно уровне насыщения.")
    doc.add_paragraph(
        "Практика подтверждает: в рабочей области решение находится за два-шесть шагов. Есть, однако, "
        "оговорка вычислительного свойства. Кривая намагничивания задана таблицей и обычно достраивается "
        "отрезками прямых, поэтому её крутизна меняется скачками на узлах таблицы. Метод Ньютона от таких "
        "изломов может начать топтаться: на глубоком насыщении, при пятикратном номинальном возбуждении, "
        "сходимость теряется. Увеличение числа шагов не помогает — невязка выходит на полку и стоит "
        "на ней хоть три тысячи шагов. Лечится это заменой отрезков прямых на монотонную гладкую "
        "интерполяцию по тем же узлам: после неё тот же режим сходится за десять шагов.")
    doc.add_paragraph(
        "Важно понимать, чего эта замена НЕ даёт. Остаточная невязка на полке составляет менее "
        "процента от величины потока, поэтому ошибку модели на краю диапазона она не объясняет: "
        "переподгонка на сглаженной кривой меняет среднюю ошибку предсказания с 2.2 до 2.1 процента, "
        "то есть в пределах разброса самой подгонки. Гладкая интерполяция нужна для надёжности "
        "решателя, а не для точности; ошибка на краю диапазона — свойство модели, а не численный "
        "артефакт.")

    doc.add_heading("Почему ответ единственный", level=3)
    doc.add_paragraph(
        "Уравнения шага 4 равносильны поиску наименьшего значения одной величины — она играет роль "
        "рельефа, по которому расчёт скатывается в низшую точку.")
    _formula(doc,
             r"E(\Phi) = \frac{1}{2}\sum_{b}\sum_{c}\Phi_b R_{bc}\Phi_c"
             r" + \sum_{b} \ell_b S_b \int_{0}^{\Phi_b / S_b} H(\beta)\, d\beta"
             r" - \sum_{b} F_b \Phi_b", lo)
    _vars(doc, [
        ("E(Φ)", "скаляр, джоуль", "минимизируемая величина"),
        ("три слагаемых", "скаляры, джоуль",
         "энергия воздуха, энергия намагничивания железа, работа обмоток"),
        ("β", "переменная интегрирования, тесла", "индукция от нуля до достигнутой"),
        ("R", "матрица m × m, А·виток/Вб", "обратная к Q"),
    ])
    doc.add_paragraph(
        "Этот рельеф имеет форму чаши без побочных ямок: куда шарик ни положи, он скатится в одну и ту же "
        "точку. Формально это следует из того, что кривая намагничивания монотонна — большему напору "
        "всегда отвечает большее поле. Отсюда вывод, важный для доверия к расчёту: у задачи ровно один "
        "ответ, и он не зависит от того, с какого приближения начали.")

    doc.add_heading("Шаг 6. От потоков к полю в зазоре", level=2)
    doc.add_paragraph(
        "В зазоре нет ни токов, ни железа, поэтому поле там подчиняется простым уравнениям, решение "
        "которых полностью определяется тем, сколько потока вышло из каждого полюса. Отсюда следует, что "
        "поле в любой точке зазора — просто сумма вкладов ветвей.")
    _formula(doc, r"\mathbf{B}(\mathbf{r}) = \sum_{b} \Phi_b \cdot \mathbf{G}_b(\mathbf{r})", lo)
    _vars(doc, [
        ("B(r)", "вектор из 3 составляющих, тесла", "поле в точке r"),
        ("r", "вектор из 3 координат, метр", "точка в зазоре"),
        ("Φ_b", "скаляр, вебер", "поток ветви b"),
        ("G_b(r)", "вектор-функция, Тл/Вб", "поле от единичного потока ветви b"),
    ])
    doc.add_paragraph(
        "Насыщение в это выражение не входит вовсе. Оно уже полностью учтено тем, какими получились "
        "потоки. Отсюда главное следствие модели: когда железо садится, потоки перераспределяются между "
        "ветвями, и меняется не только величина поля, но и его форма.")

    doc.add_heading("Шаг 7. От поля к мультиполям", level=2)
    doc.add_paragraph(
        "Поле в зазоре раскладывают на простые составляющие: его снимают по окружности вокруг оси пучка "
        "и разделяют на части, отвечающие постоянной, линейно нарастающей, квадратичной и так далее "
        "зависимости от расстояния до оси.")
    _formula(doc,
             r"C_n = \frac{1}{2\pi R^{\,n-1}} \oint \left( B_v + i\,B_u \r\right)"
             r" e^{-i(n-1)\theta}\, d\theta", lo)
    _vars(doc, [
        ("C_n", "комплексное число, Тл/м^(n−1)",
         "составляющая n: 1 — отклоняющая, 2 — фокусирующая, 3 — секступольная"),
        ("R", "скаляр, метр", "радиус окружности разложения, внутри области хорошего поля"),
        ("θ", "скаляр, радиан", "угол по окружности"),
        ("B_u, B_v", "скаляры, тесла", "поперечные составляющие поля"),
        ("i", "мнимая единица", "запись двух составляющих одним числом"),
        ("n", "индекс, 1…k", "номер составляющей"),
    ])
    doc.add_paragraph(
        "Поскольку это суммирование линейно, а поле по шагу 6 есть сумма вкладов ветвей, коэффициенты "
        "тоже выражаются через потоки ветвей линейно:")
    _formula(doc, r"\int C_n \, dl = \sum_{b} G_{n b} \cdot \Phi_b", lo)
    _vars(doc, [
        ("∫ C_n dl", "комплексное число, Тл·м/м^(n−1)", "составляющая n вдоль пути пучка"),
        ("G", "матрица k × m, ед. коэф./Вб", "вклады потоков ветвей в составляющие"),
        ("Φ_b", "скаляр, вебер", "поток ветви b"),
        ("l", "скаляр, метр", "координата вдоль оси пучка"),
    ])

    doc.add_heading("Шаг 8. От мультиполей к действию на пучок", level=2)
    doc.add_paragraph(
        "Последнее звено — перевод составляющих поля в изменение направления полёта частицы. Для "
        "короткого магнита оно даётся прямой формулой.")
    _formula(doc,
             r"\Delta x' + i\,\Delta y' = \mathrm{QP} \sum_{n} \left( \int C_n \, dl \r\right) \zeta^{\,n-1}", lo)
    _vars(doc, [
        ("Δx', Δy'", "скаляры, радиан", "изменение наклона траектории"),
        ("QP", "скаляр, 1/(м·Тл)", "жёсткость пучка: заряд, делённый на импульс"),
        ("∫ C_n dl", "комплексное число, Тл·м/м^(n−1)", "составляющие поля из шага 7"),
        ("ζ", "комплексное число, метр", "положение частицы в сечении"),
        ("n", "индекс, 1…k", "номер составляющей"),
    ])
    doc.add_paragraph(
        "Полученные коэффициенты — это и есть величины mat_B, mat_G, mat_S и остальные, которыми "
        "программа описывает каждый магнит решётки. Точные знаки и нормировки приведены в документе "
        "docs/algorithms.md; здесь важна структура.")

    doc.add_heading("Шаг 9. Проверка на частном случае", level=2)
    doc.add_paragraph(
        "У магнита с одной обмоткой и одной ветвью вся система сворачивается в единственное равенство:")
    _formula(doc, r"N \cdot I = A \cdot y + L \cdot H\left(C \cdot y\r\right)", lo)
    _vars(doc, [
        ("N · I", "скаляр, А·виток", "ампер-витки обмотки"),
        ("y", "скаляр, единицы коэффициента", "коэффициент магнита"),
        ("A", "скаляр, (А·виток) на единицу y", "линейный член: зазор"),
        ("L", "скаляр, метр", "длина железного пути"),
        ("C", "скаляр, Тл на единицу y", "пересчёт коэффициента в индукцию железа"),
        ("H(·)", "функция, Тл → А/м", "кривая намагничивания"),
    ])
    doc.add_paragraph(
        "Это в точности та формула, которая ранее была проверена на квадруполях, секступоле, октуполе и "
        "диполе-квадруполе с ошибкой один-два процента. Новая модель прежнюю не отменяет, а обобщает: "
        "все накопленные подгонки переносятся без изменений как частный случай с одной ветвью.")
    doc.add_paragraph(
        "Если на одну ветвь надеты несколько катушек — как у секступоля, — в левой части этого равенства "
        "стоит сумма их ампер-витков, а всё остальное не меняется. Насыщение при этом считается от "
        "суммарного напора, что и наблюдается на опыте.")

    doc.add_heading("Модель целиком, одной строкой", level=2)
    doc.add_paragraph(
        "На входе задачи стоит вектор токов обмоток, на выходе — вектор физических коэффициентов "
        "магнита. Всё, что между ними, — это три матрицы, две диагональные матрицы и одна скалярная "
        "функция. Ничего другого в модели нет.")
    doc.add_paragraph("Первое звено переводит токи в напоры, приложенные к ветвям:")
    _formula(doc, r"\mathbf{F} = T \cdot \mathbf{I}", lo)
    doc.add_paragraph(
        "Второе звено — единственное нелинейное. Оно неявное: искомый вектор потоков стоит и слева, и "
        "справа. Функция материала применяется поэлементно, к каждой ветви отдельно.")
    _formula(doc, r"\Phi = Q \left( \mathbf{F} - \Lambda \, H\!\left( \Sigma^{-1} \Phi \right) \right)", lo)
    doc.add_paragraph("Третье звено переводит потоки ветвей в коэффициенты магнита:")
    _formula(doc, r"\mathbf{y} = G \cdot \Phi", lo)
    doc.add_paragraph("Вся модель одной строкой:")
    _formula(doc, r"\mathbf{y} = G\, \Phi, \qquad \Phi = Q\left( T\, \mathbf{I} - \Lambda\, H(\Sigma^{-1}\Phi) \right)", lo)
    doc.add_paragraph(
        "Если убрать железо, то есть положить функцию материала равной нулю, остаётся чистое "
        "произведение трёх матриц — это и есть матрица связи, которой программа пользуется сейчас:")
    _formula(doc, r"\mathbf{y} = G\,Q\,T\,\mathbf{I} = K\,\mathbf{I}", lo)
    _vars(doc, [
        ("I", "вектор n, ампер", "токи обмоток — вход"),
        ("y", "вектор k, единицы коэффициентов", "коэффициенты магнита — выход"),
        ("F", "вектор m, А·виток", "напоры по ветвям"),
        ("Φ", "вектор m, вебер", "потоки ветвей"),
        ("T", "матрица m × n, витки", "размещение обмоток по ветвям"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей, симметрична"),
        ("Λ", "диагональная матрица m × m, метр", "длины железа ветвей"),
        ("Σ", "диагональная матрица m × m, м²", "сечения железа ветвей"),
        ("H", "скалярная функция, Тл → А/м", "кривая намагничивания, поэлементно"),
        ("G", "матрица k × m, ед. коэф./Вб", "вклады ветвей в коэффициенты"),
        ("K = G Q T", "матрица k × n, ед. коэф./А", "матрица связи — линейный предел"),
        ("n, m, k", "скаляры, штук", "число обмоток, ветвей, коэффициентов"),
    ], title="Все объекты модели")
    doc.add_paragraph(
        "Три матрицы удобно читать как графы. Матрица размещения обмоток соединяет токи с ветвями; "
        "матрица проводимостей — ветви между собой; матрица вкладов — ветви с коэффициентами. "
        "Диагональные матрицы графа не имеют: они действуют на каждую ветвь по отдельности.")
    _pic(doc, F["qpdq"], cm=16)
    _fcap(doc, "Устройство трёх элементов, от которого и берётся граф ветвей. Квадруполь — одна ветвь "
               "и одна катушка. Секступоль — одна ветвь и две катушки, их ампер-витки складываются. "
               "Диполь-квадруполь — общий сердечник, раздваивающийся на концах: I1 на одной ветви, "
               "I2 и I3 на другой. Схемы показывают связность ветвей, а не истинные пропорции")
    _pic(doc, F["graphs"], cm=16.5)
    _fcap(doc, "Три матрицы модели корректора как графы. Слева — размещение обмоток: каждая обмотка "
               "сидит ровно на одной ветви. В середине — проводимости: соседние ветви одной стороны "
               "связаны встречно (красное), противолежащие через зазор — согласно (синее), ярмо связано "
               "со всеми. Справа — вклады ветвей в коэффициенты поля")
    doc.add_paragraph(
        "Из графа проводимостей сразу видно устройство магнита. Внутри верхнего и нижнего рядов связи "
        "отрицательные: поток, вышедший из одной ветви, возвращается через соседнюю. Между "
        "противолежащими ветвями связи положительные: там поток идёт через зазор в одном направлении. "
        "Ярмо связано со всеми ветвями положительно и сильно — через него замыкается основной поток.")



def _ch_params(doc, lo, F):
    doc.add_heading("5. Сколько в модели подбираемых параметров", level=1)
    doc.add_paragraph(
        "Подавляющая часть чисел модели измеряется, а не подбирается. Подбираются только те, которые "
        "нельзя получить прямым измерением.")
    _table(doc, ["Что", "Сколько", "Подбирается?"], [
        ("Столбцы матрицы проводимостей для ветвей с обмотками",
         "n столбцов", "нет — измеряются, по одной линейной карте на обмотку"),
        ("Столбцы для ветвей без обмоток", "—", "нет — достраиваются по симметрии"),
        ("Элементы матрицы, связывающие ветви без обмоток между собой",
         "u(u+1)/2, где u — число ветвей без обмоток", "да, нелинейно"),
        ("Сечение и длина железа каждой ветви", "2 на ветвь, всего 2m", "да, нелинейно"),
        ("Матрица G (вклады ветвей в мультиполи)", "k × m", "нет — обычный линейный МНК, поиска не требует"),
    ])
    _formula(doc, r"P = 2m + \frac{u\,(u+1)}{2}", lo)
    _vars(doc, [
        ("P", "скаляр, штук", "число нелинейно подбираемых параметров"),
        ("m", "скаляр, штук", "число ветвей"),
        ("u", "скаляр, штук", "число ветвей без обмоток"),
    ])
    _table(doc, ["Элемент", "Обмоток n", "Ветвей m", "Ветвей без обмоток u", "Подбираемых параметров"], [
        ("QP (квадруполь)", "1", "1", "0", "2"),
        ("DQ (диполь-квадруполь)", "3", "2 (по составляющей)", "0", "4 = 2 на составляющую"),
        ("Корректор", "5", "7", "2", "16"),
    ], caption="Число подбираемых параметров по элементам")
    doc.add_paragraph(
        "Для сравнения: у корректора эти шестнадцать чисел определяются по восемнадцати расчётам поля, "
        "каждый из которых даёт шесть измеренных потоков, то есть по 108 измерениям. Задача переопределена "
        "почти в семь раз.")
    doc.add_paragraph(
        "Отдельно стоит сказать, что не все параметры одинаково важны. Проверка чувствительности показала: "
        "ошибка в сечении ветви, несущей обмотку, на 30 процентов сдвигает результат на 3–4 процента, а "
        "такая же ошибка в сечении ветви без обмотки — всего на 0.7 процента. Труднее всего определяются "
        "именно те параметры, которые меньше всего влияют на ответ.")


def _ch_calib(doc, lo, F):
    doc.add_heading("6. Какие измерения нужны, чтобы настроить модель", level=1)
    doc.add_paragraph(
        "Под измерением здесь понимается расчёт поля магнита при заданном возбуждении — карта поля, "
        "полученная в программе конечно-элементного анализа.")

    doc.add_heading("6.1. Что задаёт человек, а что находит программа", level=2)
    _table(doc, ["Вход", "Кто даёт"], [
        ("Список ветвей магнитопровода и какая обмотка на какой ветви", "человек, по чертежу"),
        ("Число витков и номинальный ток каждой обмотки", "человек"),
        ("Кривая намагничивания материала", "человек, файл таблицы"),
        ("Три карты поля на обмотку плюс две-три проверочные", "человек"),
        ("Поверхности, по которым считается поток ветви", "программа, по картам"),
        ("Сечения и длины железа ветвей", "программа, подбором"),
        ("Матрица проводимостей и вклады ветвей в поле", "программа, измерением"),
    ])
    doc.add_paragraph(
        "Сечения магнитопровода задавать не требуется, и это не упрощение, а следствие проверки "
        "чувствительности из предыдущего раздела: параметры, которые трудно измерить, слабо влияют на "
        "результат. Две независимые настройки модели корректора дали для средней ветви сечения 6.8 и "
        "18.5 квадратных сантиметра, а предсказания совпали в пределах 0.2 процента.")

    doc.add_heading("6.2. Сколько карт поля нужно", level=2)
    _table(doc, ["Карта", "Возбуждение", "Что из неё извлекается"], [
        ("Первая", "около 10 процентов номинального тока",
         "столбец матрицы проводимостей целиком и линейный наклон"),
        ("Вторая", "около 50 процентов номинального тока", "первая точка кривой насыщения"),
        ("Третья", "около 300 процентов номинального тока, просадка 25–35 процентов",
         "вторая точка кривой насыщения"),
    ], caption="Три карты на каждую обмотку — минимальный набор для настройки")
    doc.add_paragraph(
        "Требование одно, но обязательное: каждая обмотка снимается отдельно, и у каждой должна быть "
        "глубокая карта. Возвратный поток обмотки проходит через соседние ветви, и именно глубокий "
        "прогон их нагружает и проявляет. Для корректора это 5 × 3 = 15 карт.")
    doc.add_paragraph(
        "Проверено прямым опытом, что урезать набор нельзя. Если обойтись девятью картами, возбудив "
        "только три обмотки из пяти, потоки ещё считаются верно (ошибка 5 процентов), но перевод потоков "
        "в поле рассыпается: ветви, которые ни разу не возбуждались, остаются неопределёнными, и ошибка "
        "предсказания вырастает до 24 процентов.")

    doc.add_heading("6.3. Нужны ли карты при нескольких включённых обмотках", level=2)
    doc.add_paragraph(
        "Для настройки — нет. Модель, подогнанная только на одиночных обмотках, предсказала все "
        "тринадцать комбинационных прогонов корректора с той же точностью, с какой описывает саму "
        "подгонку.")
    doc.add_paragraph(
        "Для проверки — да, и это не формальность. Согласие модели на тех данных, по которым она "
        "настраивалась, не доказывает ничего: предыдущая версия модели прекрасно согласовывалась с "
        "подгонкой и разошлась с реальностью именно на комбинациях. Достаточно двух-трёх карт на элемент "
        "при рабочих токах, выбранных там, где обмотки делят железо. Если пару можно включить в обе "
        "полярности — снять обе: знак меняет результат качественно.")

    doc.add_heading("6.4. Порядок настройки", level=2)
    for i, s in enumerate((
        "Задать список ветвей и какая обмотка на какой сидит.",
        "Программа находит полюса по картам как области наибольшего поля и строит поверхности для "
        "измерения потоков; найденную разметку она печатает, чтобы её можно было сверить с чертежом.",
        "Собрать матрицу проводимостей из линейных карт и достроить по симметрии.",
        "Проверить симметрию. Расхождение сверх десяти процентов означает, что разметка неверна, и "
        "дальше идти нельзя.",
        "Подобрать сечения и длины железа по одиночным картам.",
        "Определить вклады ветвей в составляющие поля линейным МНК.",
        "Проверить на комбинационных картах, которые в подборе не участвовали.",
    ), 1):
        doc.add_paragraph(f"{i}. {s}")


def _ch_accuracy(doc, lo, F):
    doc.add_heading("7. Точность по типам элементов", level=1)
    doc.add_paragraph(
        "Ниже сведены две разные величины. «Перекрёстное насыщение» — насколько поле при нескольких "
        "включённых обмотках отличается от суммы полей тех же обмоток по отдельности; это мера того, "
        "нужна ли сеть ветвей вообще. «Точность» — ошибка модели после настройки.")
    _table(doc, ["Элемент", "Обмоток", "Как соединены обмотки", "Правило сложения", "Точность"], [
        ("QP_sh (квадруполь короткий)", "1", "одна ветвь", "—", "1.7 %"),
        ("QP_L (квадруполь длинный)", "1", "одна ветвь", "—", "1.5 %"),
        ("OP (октуполь)", "1", "одна ветвь", "—", "0.8 %"),
        ("SP (секступоль)", "2", "обе катушки на ОДНОЙ ветви",
         "складываются ампер-витки", "0.2…0.7 %"),
        ("DQ (диполь-квадруполь)", "3", "I1 на одной ветви, I2 и I3 на другой",
         "складываются поля", "0.4…0.9 %"),
        ("Корректор", "5", "шесть ветвей, влияние через общую среднюю",
         "полная сеть ветвей", "2.1 %"),
    ], caption="Устройство элементов, правило сложения обмоток и точность модели")
    doc.add_paragraph(
        "Правило сложения обмоток не выбирается произвольно — оно прямо следует из того, как обмотки "
        "сидят на ветвях, и проверяется по картам. Если две обмотки надеты на ОДНУ ветвь, они гонят "
        "поток по одному и тому же железу, поэтому складывать надо их ампер-витки, а насыщение считать "
        "уже от суммы. Если обмотки сидят на РАЗНЫХ ветвях, каждая насыщает своё железо, и складывать "
        "надо готовые поля.")
    doc.add_paragraph(
        "Проверка на измеренных картах показывает, что подмена правила стоит дорого, особенно при "
        "встречном включении обмоток.")
    _table(doc, ["Элемент и режим", "Складываем ампер-витки", "Складываем поля", "Верное правило"], [
        ("SP: I1 = 1891, I2 = +800 ампер-витков", "0.7 %", "1.0 %", "ампер-витки"),
        ("SP: I1 = 1891, I2 = −800 ампер-витков", "0.2 %", "3.7 %", "ампер-витки"),
        ("DQ: I1 = 7747, I2 = +800 ампер-витков", "17.9 %", "0.4 %", "поля"),
        ("DQ: I1 = 4996, I2 = +882 ампер-витков", "41.8 %", "0.9 %", "поля"),
    ], caption="Ошибка каждого правила сложения; верное правило задаётся устройством элемента")
    doc.add_paragraph(
        "У корректора ни одно из двух простых правил не работает: часть его обмоток делит общую среднюю "
        "ветвь, часть стоит на своих, и результат зависит ещё и от знаков токов. Именно поэтому ему нужна "
        "полная сеть ветвей. У DQ и SP при известном правиле сложения достаточно одноветвевой формулы, и "
        "сеть не требуется.")

    doc.add_heading("7.1. Как понимать цифры точности", level=2)
    doc.add_paragraph(
        "Приведённые проценты — это ошибка, отнесённая к полному размаху величины по всем режимам, а не "
        "к значению в каждом отдельном режиме. Такая нормировка отвечает на вопрос «насколько ошибётся "
        "кик по сравнению с максимальным киком этого магнита», и именно она важна для расчёта пучка.")
    doc.add_paragraph(
        "Если же относить ошибку к значению в каждом режиме по отдельности, у корректора получается "
        "медиана 7 процентов, а в режимах, где рассматриваемая составляющая мала сама по себе, "
        "относительная ошибка доходит до 40 процентов при абсолютной ошибке 0.1 процента от размаха. Обе "
        "цифры описывают одно и то же и не противоречат друг другу; выбирать метрику следует по задаче.")

    doc.add_heading("7.2. Перекрёстное насыщение зависит от знака токов", level=2)
    doc.add_paragraph(
        "Самая жёсткая из выполненных проверок — одна и та же пара обмоток при одинаковых по величине "
        "токах, но разного знака. Величины токов одни и те же, результат противоположный:")
    _table(doc, ["Включение обмоток c1 и c2", "Осталось от c1", "Осталось от c2", "Посторонняя составляющая"], [
        ("Согласное: +1808 и +1808 ампер-витков", "0.941", "0.927", "0.065"),
        ("Встречное: +1808 и −1808 ампер-витков", "0.986", "0.992", "−0.001"),
    ], caption="Влияние знака токов на перекрёстное насыщение корректора")
    doc.add_paragraph(
        "При согласном включении обмотки перегружают общую среднюю ветвь, и лишний поток выталкивается в "
        "зазор, создавая составляющую поля, которой нет ни в одной одиночной карте. При встречном поток "
        "замыкается поверху с полюса на полюс, средняя ветвь несёт около нуля, и ничего не происходит. "
        "Предсказание встречного случая было записано до расчёта и подтвердилось по обоим числам.")
    doc.add_paragraph(
        "Измерен и порог: при половинном токе, 904 ампер-витка на обмотку, перекрёстного эффекта нет "
        "вовсе. Номинал корректора равен 1808 ампер-виткам, то есть рабочая точка находится сразу за "
        "порогом — это стоит учитывать при выборе режимов.")


def _ch_material(doc, lo, F):
    doc.add_heading("8. Смена материала сердечника", level=1)
    doc.add_paragraph(
        "В выражении для расхода напора на железе роли разделены строго. Длина и сечение — геометрия, от "
        "материала не зависят. Кривая намагничивания — материал, и он входит только сюда. Поэтому смена "
        "стали сводится к замене одной таблицы, без повторного подбора чего бы то ни было.")
    doc.add_paragraph(
        "Модель допускает и разные материалы у разных ветвей: если полюса и ярмо сделаны из разной стали, "
        "каждая ветвь получает собственную таблицу.")
    doc.add_paragraph(
        "Утверждение проверено прямым опытом. Числа для квадруполя были подобраны только на одной стали "
        "и затем без единой поправки применены к материалу с насыщением 1.198 тесла против 2.085 у "
        "исходного, то есть с контрастом в 1.74 раза. Предсказание легло в 1.1 процента на всём участке, "
        "где отдача магнита падала со ста процентов до шестидесяти девяти.")
    doc.add_heading("Оговорка, которую важно не пропустить", level=2)
    doc.add_paragraph(
        "Матрица проводимостей измеряется на реальном магните, и в неё входит проницаемость железа при "
        "слабом возбуждении. Строго не зависела бы от материала только матрица чисто воздушной цепи. "
        "Поэтому замена материала точна для насыщающейся части, но если у нового материала заметно другая "
        "начальная проницаемость, матрицу придётся измерить заново. В описанной проверке это учитывалось "
        "намеренно: у подменного материала начальная проницаемость была близка к исходной, а контраст "
        "создавался только по насыщению.")


def _ch_matrices(doc, lo, F):
    doc.add_heading("9. Как в модель входят матрица связи и матрица индуктивностей", level=1)

    doc.add_heading("9.1. Матрица связи — это линейный предел модели", level=2)
    doc.add_paragraph(
        "Матрица связи, которой программа пользуется сейчас, переводит токи обмоток в коэффициенты "
        "магнита напрямую и линейно. В новой модели она получается сама собой, если убрать расход напора "
        "на железе, то есть объединить первое и третье звенья без второго.")
    _formula(doc, r"K_{n w} = N_w \sum_{b} G_{n b} \cdot Q_{b w}", lo)
    _vars(doc, [
        ("K", "матрица k × n, ед. коэф./А", "матрица связи"),
        ("N_w", "скаляр, штук", "число витков обмотки w"),
        ("G", "матрица k × m, ед. коэф./Вб", "вклады ветвей в коэффициенты — звено 3"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей — звено 1"),
    ])
    doc.add_paragraph(
        "Отсюда следует, что матрица связи не отменяется: она остаётся верной в области малых токов и "
        "становится частным случаем модели. Это же даёт удобную проверку — модель обязана её "
        "воспроизводить.")
    doc.add_paragraph(
        "Замечание, требующее отдельного разбирательства: для корректора отношения между обмотками в "
        "нынешней таблице расходятся с картами поля на семь-восемь процентов. Причина не установлена. "
        "Разобраться следует прежде, чем сверять модель с этой таблицей.")

    doc.add_heading("9.2. Матрица индуктивностей в модель не входит", level=2)
    doc.add_paragraph(
        "Матрица индуктивностей, которую выдаёт программа расчёта поля, не входит в модель ни одним "
        "числом. Её роль оказалась разведочной, и в этом качестве она была очень полезна: по знакам и "
        "величинам её внедиагональных членов сразу читается устройство магнита — какие катушки стоят "
        "друг напротив друга через зазор, какие соседствуют на одной стороне, какие развязаны. Она же "
        "подтвердила, какая катушка расчётной модели какой обмотке соответствует.")
    doc.add_paragraph(
        "Но попытка взять из неё веса, с которыми обмотки нагружают железо, провалилась количественно: "
        "для основной обмотки она дала 0.54 против измеренных 1.08, то есть ошиблась ровно вдвое. Причина "
        "принципиальная и видна из состава самой величины.")
    _formula(doc, r"L_{jk} = N_j N_k \sum_{b} \alpha_{j b} \cdot Q_{b k} \; + \; \Lambda_{jk}", lo)
    _vars(doc, [
        ("L", "матрица n × n, генри", "индуктивности обмоток"),
        ("N_j, N_k", "скаляры, штук", "числа витков обмоток j и k"),
        ("α", "матрица n × m, безразмерная", "доля потока ветви, сцепленная с обмоткой"),
        ("Q", "матрица m × m, Вб/(А·виток)", "проводимости ветвей"),
        ("Λ", "матрица n × n, генри", "рассеяние: поток мимо насыщающегося железа"),
    ])
    doc.add_paragraph(
        "Индуктивность считает весь поток, сцепленный с витками, включая рассеяние, которое до "
        "насыщающегося железа не доходит. Насыщение же определяется потоком в конкретном сечении железа. "
        "Это разные величины, и последнее слагаемое их разводит. Поэтому матрица проводимостей берётся из "
        "карт поля, где поток по сечению ветви измеряется напрямую, а не из индуктивностей.")


def _ch_limits(doc, lo, F):
    doc.add_heading("10. Границы применимости", level=1)
    _table(doc, ["Ограничение", "Суть", "Что с этим делать"], [
        ("Потоки измерены косвенно",
         "Полюсные наконечники корректора лежат за границей выгруженной области, поэтому поток считается "
         "в горловине, и часть его уходит мимо. Сечения и длины получаются эффективными, а не чертёжными.",
         "Выгрузить те же готовые расчёты на большей области, включающей железо. Это обработка готового "
         "результата, повторный расчёт не нужен."),
        ("Край диапазона",
         "При шести-семи номинальных токах ошибка вырастает до 9–21 процента.",
         "В рабочей области до двух номиналов точность 2–5 процентов."),
        ("Гистерезис не учитывается",
         "Кривая намагничивания принята однозначной: одному напору отвечает одно поле независимо от "
         "предыстории.",
         "Для проектных расчётов это обычное допущение. При разговоре о цикле перемагничивания оно "
         "перестаёт работать."),
        ("Материал однороден внутри ветви",
         "Слоистость сердечника и зависимость свойств от направления в модель не заведены.",
         "При необходимости ветвь разбивается на несколько с разными таблицами."),
        ("Разбиение на ветви должно быть физическим",
         "Если ветви размечены неверно, симметрия матрицы проводимостей нарушается.",
         "Проверять симметрию перед подбором — встроенный контроль правильности разметки."),
    ])


# ─────────────────────────────── приложения ────────────────────────────────

def _ch_material_data(doc, lo, F):
    import pandas as pd
    doc.add_heading("Приложение А. Кривая намагничивания steel_1008", level=1)
    doc.add_paragraph(
        "Единственное место, где в модель входит материал. Таблица читается так: чтобы получить в стали "
        "индукцию B, требуется напряжённость поля H. Промежуточные значения берутся линейной "
        "интерполяцией, за последней точкой — по наклону воздуха.")
    _pic(doc, F["bh"], cm=16)
    _fcap(doc, "Кривая намагничивания steel_1008. Слева — колено, то есть рабочая область; "
               "в середине — вся таблица целиком; справа — во сколько раз железо лучше воздуха. "
               "После 2 тесла проницаемость падает почти до единицы: железо становится "
               "неотличимо от пустоты")
    bh = pd.read_csv(BH_CSV, sep=";")
    rows = [(f"{r.H_A_per_m:g}", f"{r.B_T:g}") for r in bh.itertuples()]
    half = (len(rows) + 1) // 2
    pairs = [(rows[i][0], rows[i][1],
              rows[i + half][0] if i + half < len(rows) else "",
              rows[i + half][1] if i + half < len(rows) else "") for i in range(half)]
    _table(doc, ["H, А/м", "B, Тл", "H, А/м", "B, Тл"], pairs,
           caption="Таблица steel_1008 (файл workspace/materials/steel_1008_bh.csv)")


def _ch_example_qp(doc, lo, F):
    doc.add_heading("Приложение Б. Пример 1 — квадруполь QP_sh (одна обмотка)", level=1)
    doc.add_paragraph(
        "Простейший случай: одна обмотка, одна ветвь, один коэффициент. Все числа ниже измерены, "
        "и по ним расчёт можно повторить.")
    _pic(doc, F["gqp"], cm=15)
    _fcap(doc, "Квадруполь: обмотка сидит на единственной ветви, ветвь даёт единственный коэффициент. "
               "Все три матрицы модели вырождаются в числа")

    doc.add_heading("Б.1. Модель в одну строку", level=2)
    _formula(doc, r"\mathbf{y} = G\,\Phi, \qquad \Phi = Q\left( T\,\mathbf{I} - \Lambda\, H(\Sigma^{-1}\Phi) \right)", lo)
    doc.add_paragraph(
        "При одной обмотке и одной ветви все матрицы становятся числами, и систему удобно разрешить "
        "относительно ампер-витков — тогда неизвестной остаётся сам коэффициент магнита:")
    _formula(doc, r"N I = A\,y + L\,H\!\left( C\,y \right)", lo)
    _vars(doc, [
        ("N I", "скаляр, А·виток", "ампер-витки обмотки"),
        ("y", "скаляр, тесла", "интегральный градиент ∫G dl"),
        ("A", "скаляр, (А·виток)/Тл", "линейный член: зазор"),
        ("L", "скаляр, метр", "длина железного пути"),
        ("C", "скаляр, безразмерный", "пересчёт y в индукцию железа"),
        ("H", "функция, Тл → А/м", "кривая steel_1008 из приложения А"),
    ])

    doc.add_heading("Б.2. Расчётная формула", level=2)
    _table(doc, ["Величина", "Значение", "Единицы"], [
        ("A", "458.2", "(А·виток)/Тл"),
        ("C", "0.1539", "безразмерная"),
        ("L", "0.0652", "м"),
        ("Эффективная длина", "0.162", "м"),
        ("Матрица связи K (линейный предел)", "−0.00902354", "Тл на ампер"),
        ("Рабочий ток", "3.823", "А"),
    ], caption="Параметры модели QP_sh")
    doc.add_paragraph(
        "Подставляем числа. Возбуждение 6588 ампер-витков; неизвестной остаётся только y, "
        "а единственная величина, не заменённая числом, — закон намагничивания:")
    _formula(doc, r"6588 = 458.2\,y + 0.0652 \cdot H\!\left( 0.1539\,y \right)", lo)

    doc.add_heading("Б.3. Что перебирать", level=2)
    doc.add_paragraph(
        "Уравнение неявное, поэтому решение ищется перебором одной величины. Минимизируется невязка:")
    _formula(doc, r"\mathrm{Loss}(y) = \left| 458.2\,y + 0.0652 \cdot H(0.1539\,y) - 6588 \right| \to \min_{y}", lo)
    doc.add_paragraph(
        "Перебирается один скаляр — коэффициент y. Начальное приближение даёт линейная оценка "
        "6588 / 458.2 = 14.38 тесла.")
    _table(doc, ["Величина", "Значение"], [
        ("Решение модели", "y = 12.268 Тл"),
        ("Измерено", "12.375 Тл"),
        ("Расхождение", "0.9 %"),
        ("Линейная оценка без насыщения", "14.38 Тл, ошибка +16 %"),
    ], caption="Результат расчёта QP_sh при 6588 ампер-витков")
    doc.add_paragraph(
        "На всей серии из шести карт, от 100 до 13176 ампер-витков, средняя ошибка модели 1.0 процента, "
        "максимальная 1.9 процента. Отдача магнита на этом участке падает со 100 до 51 процента.")


def _ch_example_dq(doc, lo, F):
    doc.add_heading("Приложение В. Пример 2 — диполь-квадруполь DQ (три обмотки)", level=1)
    doc.add_paragraph(
        "Есть общий сердечник; на каждом из его концов он разделяется на две ветви. На одной ветви "
        "сидит обмотка I1, на другой — обмотки I2 и I3 вместе. Снизу устроено так же. Расположение "
        "обмоток по сторонам подтверждается картами: поле обмотки I1 втрое сильнее с одной стороны "
        "апертуры, а поле I2 и I3 — втрое сильнее с противоположной.")
    _pic(doc, F["gdq"], cm=15.5)
    _fcap(doc, "Диполь-квадруполь: обмотка I1 на своей ветви, обмотки I2 и I3 — вместе на другой. "
               "Обе ветви сходятся в общий сердечник. Матрица связи переводит токи трёх обмоток "
               "в три составляющие поля")
    doc.add_paragraph(
        "Из этого устройства следуют два наблюдаемых свойства, и оба подтверждаются картами. Первое: "
        "обмотки I2 и I3 сидят на одной ветви, поэтому их поля почти неразличимы — измерение при "
        "800 ампер-витках даёт совпадение составляющих в пределах 0.2…0.9 процента. Второе: I1 и I2 "
        "сидят на разных ветвях и насыщают разное железо, поэтому складывать надо готовые поля, а не "
        "ампер-витки. Сложение полей даёт ошибку 0.4…0.9 процента, сложение ампер-витков — 18…42.")

    doc.add_heading("В.1. Матрица связи", level=2)
    doc.add_paragraph(
        "Линейный предел модели, то есть произведение трёх матриц. Строки — составляющие поля, "
        "столбцы — обмотки; значения на ампер, эффективная длина 0.400 метра:")
    _formula(doc,
             r"K = \begin{pmatrix}"
             r" 5.925 \cdot 10^{-5} & 3.892 \cdot 10^{-5} & 3.925 \cdot 10^{-5} \\"
             r" 5.322 \cdot 10^{-3} & -5.218 \cdot 10^{-3} & -5.194 \cdot 10^{-3} \\"
             r" -8.611 \cdot 10^{-2} & -5.351 \cdot 10^{-1} & -5.316 \cdot 10^{-1}"
             r" \end{pmatrix}", lo)
    _vars(doc, [
        ("K", "матрица 3 × 3, ед. коэф./А", "строки — составляющие B, G, S; столбцы — обмотки I1, I2, I3"),
        ("строки", "Тл/А, (Тл/м)/А, (Тл/м²)/А", "отклоняющая, фокусирующая, секступольная"),
    ])
    doc.add_paragraph(
        "Второй и третий столбцы почти совпадают — это и есть следствие того, что обмотки I2 и I3 "
        "сидят на одной ветви. По отклику на пучок они неразличимы.")

    doc.add_heading("В.2. Модель в одну строку", level=2)
    doc.add_paragraph(
        "Обмотки сидят на разных ветвях, поэтому вклады складываются на уровне полей. Обмотка I1 "
        "работает вблизи своего колена, и её вклад считается с насыщением; обмотки I2 и I3 в рабочей "
        "области линейны — измерение подтверждает постоянство отдачи в пределах 0.1 процента до "
        "882 ампер-витков и 2 процентов до 1764:")
    _formula(doc, r"y_j = y_j^{(1)} + K_{j2} I_2 + K_{j3} I_3, \qquad N_1 I_1 = A_j\,y_j^{(1)} + L_j\,H\!\left( C_j\,y_j^{(1)} \right)", lo)
    _vars(doc, [
        ("y_j", "скаляр, единицы составляющей j", "полная составляющая j"),
        ("y_j^(1)", "скаляр, единицы составляющей j", "вклад I1, с насыщением"),
        ("K_j2, K_j3", "скаляры, ед. составляющей/А", "элементы матрицы связи для I2 и I3"),
        ("I_2, I_3", "скаляры, ампер", "токи вспомогательных обмоток"),
        ("N_1 I_1", "скаляр, А·виток", "ампер-витки обмотки I1"),
        ("A_j, C_j, L_j", "скаляры", "параметры насыщения составляющей j"),
    ])

    doc.add_heading("В.3. Расчётная формула", level=2)
    _table(doc, ["Составляющая", "A", "C", "L, м", "Ошибка на серии"], [
        ("∫B dl (отклоняющая)", "39790 (А·виток)/(Тл·м)", "8.692", "0.2674", "0.8 % средняя, 2.5 % макс"),
        ("∫G dl (фокусирующая)", "447.1 (А·виток)/Тл", "0.1141", "0.0104", "0.9 % средняя, 1.7 % макс"),
    ], caption="Параметры насыщения DQ, подогнаны по девяти картам серии I1")
    doc.add_paragraph(
        "Задача: обмотка I1 возбуждена на 7747 ампер-витков, остальные выключены. Подставляем числа — "
        "не заменённым остаётся только закон намагничивания:")
    _formula(doc, r"7747 = 39790\,y_B + 0.2674 \cdot H\!\left( 8.692\,y_B \right)", lo)
    _formula(doc, r"7747 = 447.1\,y_G + 0.0104 \cdot H\!\left( 0.1141\,y_G \right)", lo)

    doc.add_heading("В.4. Что перебирать", level=2)
    doc.add_paragraph(
        "Каждое уравнение решается отдельно перебором своего одного скаляра:")
    _formula(doc, r"\mathrm{Loss}(y_B) = \left| 39790\,y_B + 0.2674 \cdot H(8.692\,y_B) - 7747 \right| \to \min_{y_B}", lo)
    doc.add_paragraph(
        "Перебираются две независимые величины, по одной на составляющую. Вклады обмоток I2 и I3 "
        "перебора не требуют вовсе — они прибавляются к результату умножением тока на элемент "
        "матрицы связи.")
    _table(doc, ["Составляющая", "Решение модели", "Индукция в железе", "Измерено", "Расхождение"], [
        ("∫B dl", "0.1785 Тл·м", "1.551 Тл", "0.17867 Тл·м", "0.1 %"),
        ("∫G dl", "16.889 Тл", "1.927 Тл", "16.8178 Тл", "0.4 %"),
    ], caption="Результат расчёта DQ при 7747 ампер-витков")


def _ch_example_corr(doc, lo, F, M):
    doc.add_heading("Приложение Г. Пример 3 — корректор (пять обмоток, сеть ветвей)", level=1)
    doc.add_paragraph(
        "Полный случай, ради которого модель и строилась. Основная обмотка сидит на ярме сбоку; "
        "магнитопровод расходится на три сердечника сверху и три снизу, стоящих друг напротив друга. "
        "Средние сердечники обмоток не несут, боковые несут четыре дополнительные обмотки.")
    _pic(doc, F["corr"], cm=14)
    _fcap(doc, "Устройство корректора: шесть сердечников и ярмо. Обозначения ветвей — T+, T0, T−, "
               "B+, B0, B−; обмотки c1…c4 на боковых сердечниках, main на ярме")
    doc.add_paragraph(
        "Те же три матрицы модели, изображённые графами, приведены на рисунке 1 в разделе 3.1.")

    doc.add_heading("Г.1. Матрица размещения обмоток", level=2)
    doc.add_paragraph(
        "Матрица T переводит токи в напоры, приложенные к ветвям. Строки — ветви в порядке "
        "T+, T0, T−, B+, B0, B−, ярмо; столбцы — обмотки c1, c2, c3, c4, main. В каждом столбце ровно "
        "одно ненулевое число — число витков той обмотки:")
    _formula(doc,
             r"T = \begin{pmatrix}"
             r" 548 & 0 & 0 & 0 & 0 \\"
             r" 0 & 0 & 0 & 0 & 0 \\"
             r" 0 & 548 & 0 & 0 & 0 \\"
             r" 0 & 0 & 548 & 0 & 0 \\"
             r" 0 & 0 & 0 & 0 & 0 \\"
             r" 0 & 0 & 0 & 548 & 0 \\"
             r" 0 & 0 & 0 & 0 & 622"
             r" \end{pmatrix}", lo)
    doc.add_paragraph(
        "Две нулевые строки — это средние ветви T0 и B0: обмоток на них нет, и напор к ним напрямую "
        "не прикладывается.")

    doc.add_heading("Г.2. Матрица проводимостей", level=2)
    doc.add_paragraph(
        "Измерена по картам полей. Строки — ветви, в которых считается поток; столбцы — ветви, куда "
        "приложен напор. Общий множитель вынесен, единицы — вебер на ампер-виток:")
    Q = M["Q_Wb_per_At"]
    rows = r" \\ ".join(" & ".join(f"{1e7 * Q[i][j]:+.2f}" for j in range(7)) for i in range(6))
    _formula(doc, r"Q = 10^{-7} \cdot \begin{pmatrix} " + rows + r" \end{pmatrix}", lo)

    doc.add_heading("Г.3. Подобранная геометрия железа", level=2)
    br = M["branches"]
    lam = r",\; ".join(f"{M['len_mm'][i] / 1000:.4f}" for i in range(6))
    sig = r",\; ".join(f"{M['S_cm2'][i] * 1e-4:.5f}" for i in range(6))
    doc.add_paragraph("Длины железных участков ветвей, метры:")
    _formula(doc, r"\Lambda = \mathrm{diag}\left( " + lam + r" \right)", lo)
    doc.add_paragraph("Сечения железа ветвей, квадратные метры:")
    _formula(doc, r"\Sigma = \mathrm{diag}\left( " + sig + r" \right)", lo)
    doc.add_paragraph(
        "Это единственные подбираемые числа модели: четырнадцать на шесть ветвей и ярмо, плюс два "
        "неизмеримых элемента матрицы проводимостей.")

    doc.add_heading("Г.4. Матрица вкладов ветвей в поле", level=2)
    doc.add_paragraph(
        "Строки — составляющие поля b1, a1, b2, a2; столбцы — ветви в том же порядке. Первые две "
        "строки в Тл·м на вебер, вторые две в Тл на вебер:")
    G = M["G"]
    grows = r" \\ ".join(" & ".join(f"{G[i][j]:+.2f}" for i in range(6)) for j in range(4))
    _formula(doc, r"G = \begin{pmatrix} " + grows + r" \end{pmatrix}", lo)

    doc.add_heading("Г.5. Расчётная формула", level=2)
    doc.add_paragraph(
        "Задача: обмотки c1 и c2 возбуждены согласованно, по 1808 ампер-витков каждая, остальные "
        "выключены. Вектор токов даёт вектор напоров по ветвям:")
    _formula(doc, r"T\,\mathbf{I} = \left( 1808 \;\; 0 \;\; 1808 \;\; 0 \;\; 0 \;\; 0 \;\; 0 \right)^{\mathsf{T}}", lo)
    doc.add_paragraph(
        "Подставляем всё в модель. Единственная величина, не заменённая числами, — закон "
        "намагничивания:")
    _formula(doc,
             r"\Phi = 10^{-7} \cdot Q_{0} \left( \left( 1808 \;\; 0 \;\; 1808 \;\; 0 \;\; 0 \;\; 0 \;\; 0 \right)^{\mathsf{T}}"
             r" - \Lambda\, H\!\left( \Sigma^{-1} \Phi \right) \right), \qquad \mathbf{y} = G\,\Phi", lo)
    _vars(doc, [
        ("Φ", "вектор 6, вебер", "потоки ветвей — искомое"),
        ("Q_0", "матрица 6 × 7, безразмерная", "проводимости из Г.2 без общего множителя"),
        ("Λ", "диагональная матрица 6 × 6, метр", "длины железа из Г.3"),
        ("Σ", "диагональная матрица 6 × 6, м²", "сечения железа из Г.3"),
        ("H", "функция, Тл → А/м", "кривая намагничивания, поэлементно"),
        ("G", "матрица 4 × 6, ед. коэф./Вб", "вклады ветвей из Г.4"),
        ("y", "вектор 4, единицы коэффициентов", "составляющие b1, a1, b2, a2 — ответ"),
    ])

    doc.add_heading("Г.6. Что перебирать", level=2)
    doc.add_paragraph(
        "Система неявная: искомый вектор стоит и слева, и справа. Минимизируется норма невязки:")
    _formula(doc,
             r"\mathrm{Loss}(\Phi) = \left\| \Phi - Q\left( T\mathbf{I} - \Lambda\,H(\Sigma^{-1}\Phi) \right) \right\|"
             r" \to \min_{\Phi}", lo)
    doc.add_paragraph(
        "Перебирается вектор из шести чисел — потоки ветвей. Начальное приближение даёт линейный "
        "расчёт, то есть та же формула с выключенным законом намагничивания. Найденный вектор "
        "умножается на матрицу вкладов и даёт ответ.")
    _table(doc, ["Ветвь", "Поток модели, мкВб", "Поток измеренный, мкВб"], [
        ("T+", "1094", "1100"), ("T0", "−1057", "−980"), ("T−", "1008", "1002"),
        ("B+", "205", "225"), ("B0", "123", "137"), ("B−", "275", "300"),
    ], caption="Потоки ветвей при c1 = c2 = 1808 ампер-витков")
    doc.add_paragraph(
        "Для сравнения: линейный расчёт без учёта железа дал бы в средней ветви −1067 микровебера. "
        "Измерено −980, то есть средняя ветвь приняла на восемь процентов меньше, чем должна была — "
        "она и есть узкое место цепи. Модель показывает −1057, то есть воспроизводит эффект лишь "
        "частично; это худшая из тринадцати проверочных комбинаций.")
    doc.add_paragraph(
        "Отклоняющая способность по этим потокам: модель даёт −2.13 миллитесла·метра, измерено −2.85. "
        "Разница составляет 1.7 процента от полного размаха этой величины по всем режимам корректора, "
        "но 26 процентов от значения в данном конкретном режиме. Обе цифры приведены намеренно: они "
        "описывают одно и то же расхождение и показывают, почему выбор метрики нужно оговаривать.")
    doc.add_paragraph(
        "В среднем по всем тринадцати комбинациям, ни одна из которых в настройке не участвовала, "
        "ошибка составляет 1.6 процента от размаха при медиане 7 процентов от собственного значения.")


def _ch_glossary(doc, lo, F):
    doc.add_heading("Приложение Д. Сводный список обозначений", level=1)
    _table(doc, ["Обозначение", "Тип", "Величина", "Единицы"], [
        ("n, m, k, u", "скаляры", "число обмоток, ветвей, составляющих поля, ветвей без обмоток", "штук"),
        ("I_w", "скаляр", "ток в обмотке w", "ампер"),
        ("N_w", "скаляр", "число витков обмотки w", "штук"),
        ("F", "вектор длины m", "напоры обмоток, разнесённые по ветвям", "ампер-виток"),
        ("Φ", "вектор длины m", "магнитные потоки ветвей", "вебер"),
        ("S_b", "скаляр", "площадь сечения железа ветви b", "м²"),
        ("ℓ_b", "скаляр", "длина железного участка ветви b", "м"),
        ("R_b", "скаляр", "сопротивление воздушной части пути ветви b", "(А·виток)/Вб"),
        ("D", "вектор длины m", "расход напора на железе ветвей", "ампер-виток"),
        ("U_b", "скаляр", "полный расход напора на ветви b", "ампер-виток"),
        ("H(B)", "функция Тл → А/м", "кривая намагничивания материала", "А/м"),
        ("Q", "матрица m × m", "проводимости ветвей", "Вб/(А·виток)"),
        ("R", "матрица m × m", "обратная к Q", "(А·виток)/Вб"),
        ("J", "диагональная матрица m × m", "крутизна сопротивления железа", "(А·виток)/Вб"),
        ("W", "скаляр", "энергия магнитного поля", "джоуль"),
        ("E(Φ)", "скаляр", "минимизируемая величина, определяющая решение", "джоуль"),
        ("G_b(r)", "вектор-функция из 3 составляющих", "картина поля от единичного потока ветви b", "Тл/Вб"),
        ("C_n", "комплексное число", "коэффициент составляющей поля n", "Тл/м^(n−1)"),
        ("g", "матрица k × m", "вклады потоков ветвей в составляющие поля", "ед. коэф. на вебер"),
        ("K", "матрица k × n", "матрица связи «ток → коэффициент»", "ед. коэф. на ампер"),
        ("L", "матрица n × n", "индуктивности обмоток", "генри"),
        ("Λ", "матрица n × n", "вклад рассеяния в индуктивности", "генри"),
        ("α", "матрица n × m", "доля потока ветви, сцепленная с обмоткой", "безразмерна"),
        ("QP", "скаляр", "жёсткость пучка: заряд, делённый на импульс", "1/(м·Тл)"),
        ("ζ", "комплексное число", "положение частицы в поперечном сечении", "метр"),
        ("P", "скаляр", "число подбираемых параметров модели", "штук"),
    ], caption="Все обозначения, встречающиеся в документе")


# ─────────────────────────────── сборка ────────────────────────────────

def build(path=None):
    from docx import Document
    import latex2omml as lo
    path = path or OUT_DEFAULT
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    _tcap.n = 0
    _fcap.n = 0
    _SEEN.clear()
    F = _figures()
    M = json.load(open(MODEL_JSON, encoding="utf-8"))
    doc = Document()
    doc.add_heading("Модель насыщения магнита: от токов обмоток к действию на пучок", level=0)
    doc.add_paragraph(
        "Пошаговое описание расчёта. Каждая формула сопровождается расшифровкой всех входящих величин с "
        "указанием единиц измерения и типа (скаляр, вектор, матрица); специальные термины раскрываются "
        "при первом появлении. В приложениях — три разобранных примера с измеренными числами, по которым "
        "расчёт можно повторить самостоятельно.")
    for ch in (_ch_intro, _ch_scheme, _ch_steps, _ch_params, _ch_calib,
               _ch_accuracy, _ch_material, _ch_matrices, _ch_limits,
               _ch_material_data, _ch_example_qp, _ch_example_dq):
        ch(doc, lo, F)
    _ch_example_corr(doc, lo, F, M)
    _ch_glossary(doc, lo, F)
    _body_format(doc)
    try:
        doc.save(path)
    except PermissionError:
        path = os.path.splitext(path)[0] + "_new.docx"
        doc.save(path)
        print(f"  (файл заблокирован Word — записан {os.path.basename(path)})")
    if not lo.available():
        print("  ВНИМАНИЕ: XSL Office не найден — формулы вставлены текстом, а не как объекты Word.")
    print(f"отчёт -> {path}")
    return path


def main():
    ap = argparse.ArgumentParser(description="Word-отчёт: модель насыщения магнита")
    ap.add_argument("--out", default=OUT_DEFAULT, help="путь к .docx")
    ap.add_argument("--to", choices=("odt", "pdf", "doc", "rtf", "html"),
                    help="дополнительно сконвертировать готовый отчёт (нужна LibreOffice)")
    a = ap.parse_args()
    path = build(a.out)
    if a.to:
        import office_convert
        office_convert.convert(path, a.to)


if __name__ == "__main__":
    main()
