# -*- coding: utf-8 -*-
"""Вкладка «Визуализация»: трёхмерный просмотр опыта и два вида HTML-отчётов.

Интерактивное окно PyVista открывается отдельным процессом (`visualize.py`) — тем же, что и из
командной строки, поэтому сцена, камера и наложение киковой модели ведут себя ровно так же.
Внутрь оболочки встроен только off-screen снимок (`viz_shot.py`) — как быстрый предпросмотр.

Два HTML-отчёта отвечают на разные вопросы:
  • **3D-сцена** (`viz_html.py`) — как пучок идёт по машине в пространстве;
  • **анализ траектории** (`traj_analyzer.py`) — что с ним происходит: x(S), y(S), спектр колебаний,
    тюны, огибающая и затухание.
"""
from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QPixmap
from PySide6.QtWidgets import (QCheckBox, QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox,
                               QHBoxLayout, QLabel, QLineEdit, QPushButton, QSpinBox, QVBoxLayout,
                               QWidget)

from . import paths, projects
from .runner import py_job
from .tabbase import TabBase
from .widgets import SectionLabel, open_in_os, plural

SHOT = paths.RESULTS / "gui_shot.png"


class VizTab(TabBase):
    def __init__(self, runner, log, parent=None):
        super().__init__(runner, log, parent)

        self.exp = QComboBox()
        self.exp.setMinimumWidth(320)
        self.btn_reload = QPushButton("Обновить список")
        self.btn_reload.clicked.connect(self.refresh)
        row_exp = QHBoxLayout()
        row_exp.addWidget(QLabel("Опыт:"))
        row_exp.addWidget(self.exp, 1)
        row_exp.addWidget(self.btn_reload)

        # --- интерактивный просмотр -------------------------------------------
        self.btn_view = QPushButton("Открыть 3D-окно (PyVista)")
        self.btn_live = QPushButton("Живой просмотр")
        self.btn_live.setToolTip("live_visualizer.py — следит за sp_elements.csv и перерисовывает "
                                 "сцену после правки токов")
        self.btn_export = QPushButton("Траектория кикового ядра → 3D")
        self.btn_export.setToolTip("matrix_traj_export.py inject — прогнать пучок КИКОВЫМ ядром и записать\n"
                                   "его траекторию в формате визуализатора (trajectories_exp<N>.bin).\n"
                                   "Нужен, когда 3D-картинку хочется от кикового ядра: длинные прогоны\n"
                                   "пишут только локальные .npz, а диффурный .bin может быть от старого счёта.")
        self.btn_view.clicked.connect(self._view)
        self.btn_live.clicked.connect(self._live)
        self.btn_export.clicked.connect(self._export_kick)
        self.traj_info = QLabel("")
        self.traj_info.setStyleSheet("color:#666;")
        self.traj_info.setWordWrap(True)

        box_view = QGroupBox("Интерактивный просмотр")
        v = QVBoxLayout(box_view)
        r = QHBoxLayout()
        r.addWidget(self.btn_view)
        r.addWidget(self.btn_live)
        r.addStretch(1)
        r2 = QHBoxLayout()
        r2.addWidget(self.btn_export)
        r2.addStretch(1)
        v.addLayout(r)
        v.addLayout(r2)
        v.addWidget(self.traj_info)

        # --- снимок ------------------------------------------------------------
        self.mag = QSpinBox()
        self.mag.setRange(0, 100000)
        self.mag.setSpecialValueText("— весь кадр —")
        self.view_mode = QComboBox()
        for t, hint in (("top", "сверху (виден изгиб)"), ("side", "сбоку"), ("persp", "перспектива")):
            self.view_mode.addItem(f"{t} — {hint}", t)
        self.span = QDoubleSpinBox()
        self.span.setRange(0.1, 1000.0)
        self.span.setValue(8.0)
        self.span.setSuffix(" м")
        self.btn_shot = QPushButton("Сделать снимок")
        self.btn_shot.clicked.connect(self._shot)
        self.preview = QLabel("снимок появится здесь")
        self.preview.setAlignment(Qt.AlignCenter)
        self.preview.setMinimumHeight(220)
        self.preview.setStyleSheet("border:1px solid #ccc; color:#888;")

        box_shot = QGroupBox("Снимок сцены (off-screen)")
        f = QFormLayout()
        f.addRow("Магнит (индекс):", self.mag)
        f.addRow("Ракурс:", self.view_mode)
        f.addRow("Ширина кадра:", self.span)
        sb = QHBoxLayout()
        sb.addWidget(self.btn_shot)
        sb.addStretch(1)
        vs = QVBoxLayout(box_shot)
        vs.addLayout(f)
        vs.addLayout(sb)
        vs.addWidget(self.preview, 1)

        # --- HTML-отчёты --------------------------------------------------------
        self.html_kind = QComboBox()
        self.html_kind.addItem("лёгкая: магниты блоками (для полного кольца)", "glyph")
        self.html_kind.addItem("детальная: настоящие STL-модели", "detailed")
        self.html_kind.addItem("обе версии сразу", "both")
        self.decimate = QDoubleSpinBox()
        self.decimate.setRange(-1.0, 1.0)
        self.decimate.setSingleStep(0.1)
        self.decimate.setValue(-1.0)
        self.decimate.setToolTip("Прореживание STL детальной версии: -1 — авто (~300 тыс. треугольников), "
                                 "0 — полный STL (для полного кольца это сотни МБ)")
        self.btn_html3d = QPushButton("Собрать HTML: 3D-сцена")
        self.btn_html3d.clicked.connect(self._html_scene)
        self.name3d = QLineEdit()                            # своё имя файла; пусто — имя по умолчанию
        self.name3d.setToolTip("Имя файла отчёта. Пусто — имя по умолчанию (показано серым).\n"
                               "Расширение .html добавится само; вариант «обе версии» допишет "
                               "_light и _detailed.")

        self.ck_csv = QCheckBox("дополнительно сводка в CSV")
        self.btn_html_traj = QPushButton("Собрать HTML: анализ траектории")
        self.btn_html_traj.setToolTip(
            "traj_analyzer.py — вид отчёта выбирается САМ по записи трека:\n"
            "• запись «локальные координаты» → подробные x(S), y(S), спектр колебаний, тюны;\n"
            "• запись «длинный турн-by-turn» → x(оборот), огибающая амплитуды, фазовый портрет\n"
            "  (тысячи оборотов, график прорежен).\n"
            "Что получится сейчас — написано под кнопками просмотра.")
        self.btn_html_traj.clicked.connect(self._html_traj)
        self.name_traj = QLineEdit()
        self.name_traj.setToolTip("Имя файла отчёта. Пусто — имя по умолчанию (показано серым). "
                                  "Расширение .html добавится само.")
        # Отчёты ложатся в РАЗНЫЕ папки (3D-сцена — в results/html, анализ — в results), поэтому
        # «открыть» у каждого своё: открываем сам файл, а если его ещё нет — его папку.
        self.btn_open3d = QPushButton("Открыть")
        self.btn_open3d.setToolTip("открыть готовый файл 3D-сцены (или её папку, если файла ещё нет)")
        self.btn_open3d.clicked.connect(lambda: self._open_report("3d"))
        self.btn_open_traj = QPushButton("Открыть")
        self.btn_open_traj.setToolTip("открыть готовый отчёт анализа (или его папку, если файла ещё нет)")
        self.btn_open_traj.clicked.connect(lambda: self._open_report("traj"))
        self.html_info = QLabel("")
        self.html_info.setStyleSheet("color:#666;")
        self.html_info.setWordWrap(True)

        box_html = QGroupBox("HTML-отчёты опыта")
        vh = QVBoxLayout(box_html)
        h1 = QHBoxLayout()
        h1.addWidget(QLabel("3D-сцена:"))
        h1.addWidget(self.html_kind, 1)
        h1.addWidget(QLabel("прореживание:"))
        h1.addWidget(self.decimate)
        h1a = QHBoxLayout()
        h1a.addWidget(QLabel("Имя файла:"))
        h1a.addWidget(self.name3d, 1)
        h1a.addWidget(self.btn_html3d)
        h1a.addWidget(self.btn_open3d)
        h2 = QHBoxLayout()
        h2.addWidget(QLabel("Анализ траектории:"))
        h2.addWidget(self.ck_csv)
        h2.addStretch(1)
        h2a = QHBoxLayout()
        h2a.addWidget(QLabel("Имя файла:"))
        h2a.addWidget(self.name_traj, 1)
        h2a.addWidget(self.btn_html_traj)
        h2a.addWidget(self.btn_open_traj)
        for r in (h1, h1a, h2, h2a):
            vh.addLayout(r)
        vh.addWidget(self.html_info)

        self.guard(self.btn_view, self.btn_live, self.btn_export, self.btn_shot,
                   self.btn_html3d, self.btn_html_traj)

        cols = QHBoxLayout()
        left = QVBoxLayout()
        left.addWidget(box_view)
        left.addWidget(box_html)
        left.addStretch(1)
        lw = QWidget()
        lw.setLayout(left)
        cols.addWidget(lw, 1)
        cols.addWidget(box_shot, 1)

        lay = QVBoxLayout(self)
        lay.addWidget(SectionLabel("Визуализация опытов"))
        lay.addLayout(row_exp)
        lay.addLayout(cols, 1)

    # ---- обновление ------------------------------------------------------
    def refresh(self) -> None:
        cur = self.exp.currentData()
        self.exp.clear()
        for e in projects.experiments():
            self.exp.addItem(f"опыт {e['exp_id']} — решётка {e.get('lattice_id', '?')}, "
                             f"пучок {e.get('group_id', '?')}", int(e["exp_id"]))
        target = self.project.exp if self.project else cur
        if target is not None:
            i = self.exp.findData(target)
            if i >= 0:
                self.exp.setCurrentIndex(i)
        self._update_info()

    def _exp(self) -> int | None:
        d = self.exp.currentData()
        return int(d) if d is not None else None

    @staticmethod
    def _bin_summary(path) -> str:
        """Кратко о файле глобальных траекторий: когда записан, сколько частиц и точек."""
        import datetime as dt
        import struct
        when = dt.datetime.fromtimestamp(path.stat().st_mtime).strftime("%d.%m %H:%M")
        try:
            with open(path, "rb") as f:
                data = f.read()
            counts, off = [], 0
            while off + 4 <= len(data):
                n = struct.unpack_from("<I", data, off)[0]
                off += 4 + 24 * n
                if n == 0 or off > len(data):
                    break
                counts.append(n)
        except OSError:
            return f"есть (от {when})"
        if not counts:
            return f"есть (от {when}), но частиц не разобрать"
        return f"есть (от {when}): частиц {len(counts)}, точек до {max(counts):,}".replace(",", " ")

    @staticmethod
    def _npz_summary(path) -> str:
        """Что лежит в локальных координатах и какой отчёт из этого получится.

        `traj_analyzer` сам выбирает вид отчёта по полю `mode` в файле: длинный турн-by-turn даёт
        огибающую и фазовый портрет (график по оборотам, прорежен), поэлементная запись — подробные
        x(S), y(S) и спектр. Кнопка одна, поэтому вид ожидаемого отчёта показываем заранее."""
        import datetime as dt
        when = dt.datetime.fromtimestamp(path.stat().st_mtime).strftime("%d.%m %H:%M")
        try:
            import numpy as np
            with np.load(path, allow_pickle=True) as d:
                mode = str(d["mode"]) if "mode" in d.files else ""
                turns = int(d["n_turns"]) if "n_turns" in d.files else 0
        except Exception:
            return f"есть (от {when})"
        w = plural(turns, "оборот", "оборота", "оборотов")
        if mode == "long":
            return (f"есть (от {when}): длинный режим, {turns} {w} → "
                    f"отчёт с огибающей амплитуды и фазовым портретом (график по оборотам)")
        return (f"есть (от {when}): запись по элементам, {turns} {w} → "
                f"подробный отчёт x(S), y(S), спектр колебаний и тюны")

    def _update_info(self) -> None:
        e = self._exp()
        if e is None:
            return
        glob = paths.RESULTS / f"trajectories_exp{e}.bin"
        loc = paths.RESULTS / f"trajectories_local_exp{e}.npz"
        parts = [f"глобальные траектории (3D): "
                 + (self._bin_summary(glob) if glob.exists() else "НЕТ — нужен трекинг в режиме «глобальные»")]
        parts.append("локальные координаты (анализ): "
                     + (self._npz_summary(loc) if loc.exists() else "нет — нужен режим «локальные» или «длинный»"))
        self.traj_info.setText("\n".join(parts))
        # имена по умолчанию — подсказкой в полях, чтобы было видно, КАК ляжет файл
        self.name3d.setPlaceholderText(self._default_name3d(e))
        self.name_traj.setPlaceholderText(f"traj_analysis_exp{e}.html")
        self.html_info.setText("\n".join(self._report_line(k) for k in ("3d", "traj")))

    def _default_name3d(self, exp: int) -> str:
        """Имя 3D-сцены по умолчанию (как его составляет viz_html)."""
        return (f"exp{exp}_3d_light.html + …_detailed.html" if self.html_kind.currentData() == "both"
                else f"exp{exp}_3d.html")

    @staticmethod
    def _out_path(name: str, folder) -> Path:
        """Введённое имя → путь файла отчёта (расширение и папка добавляются сами)."""
        n = name.strip().replace("/", "_").replace("\\", "_")
        if not n.lower().endswith(".html"):
            n += ".html"
        return folder / n

    def _report_path(self, kind: str) -> Path | None:
        """Куда ляжет отчёт: своё имя из поля либо имя по умолчанию.

        Отчёты пишутся в РАЗНЫЕ папки: 3D-сцена — `results/html`, анализ траектории — `results`.
        Путь показываем целиком, иначе готовый файл ищут не там, где он лежит."""
        e = self._exp()
        if e is None:
            return None
        if kind == "3d":
            name = self.name3d.text().strip()
            return self._out_path(name, paths.HTML_DIR) if name else paths.HTML_DIR / f"exp{e}_3d.html"
        name = self.name_traj.text().strip()
        return (self._out_path(name, paths.RESULTS) if name
                else paths.RESULTS / f"traj_analysis_exp{e}.html")

    def _report_line(self, kind: str) -> str:
        """Строка состояния отчёта: полный путь + есть ли файл и когда собран."""
        p = self._report_path(kind)
        if p is None:
            return ""
        title = "3D-сцена" if kind == "3d" else "анализ траектории"
        try:
            rel = p.relative_to(paths.ROOT)
        except ValueError:
            rel = p
        if p.exists():
            import datetime as dt
            st = p.stat()
            mark = (f"есть, {st.st_size / 1e6:.1f} МБ, "
                    f"{dt.datetime.fromtimestamp(st.st_mtime).strftime('%d.%m %H:%M')}")
        else:
            mark = "ещё не собран"
        return f"{title} → {rel}   ({mark})"

    def _open_report(self, kind: str) -> None:
        """Открыть готовый отчёт; нет файла — открыть его папку и сказать об этом."""
        p = self._report_path(kind)
        if p is None:
            return
        if p.exists():
            self._open(p)
            return
        self.say(f"файла ещё нет: {p} — открываю папку")
        self._open(p.parent)

    def _open(self, p) -> None:
        err = open_in_os(p)
        if err:
            self.say(err, True)

    # ---- команды ---------------------------------------------------------
    def _view(self) -> None:
        e = self._exp()
        if e is None:
            return
        if not (paths.RESULTS / f"trajectories_exp{e}.bin").exists():
            self.say(f"нет trajectories_exp{e}.bin — сначала прогоните трекинг "
                     f"(режим записи «глобальные траектории»)", True)
            return
        self.runner.run(py_job(f"3D-окно опыта {e}", "visualize.py", e))

    def _live(self) -> None:
        e = self._exp()
        if e is not None:
            self.runner.run(py_job(f"живой просмотр опыта {e}", "live_visualizer.py", e))

    def _export_kick(self) -> None:
        """Записать в формат визуализатора траекторию КИКОВОГО ядра (перезаписывает .bin этого опыта)."""
        e = self._exp()
        if e is None:
            return
        lat, grp = projects.experiment_config(e)
        if lat is None or grp is None:
            self.say(f"у опыта {e} не заданы установка и пучок в experiments.csv", True)
            return
        self.say(f"трасса кикового ядра заменит trajectories_exp{e}.bin (прежний диффурный расчёт)")
        self.runner.run(py_job(f"траектория кикового ядра, опыт {e} (решётка {lat}, пучок {grp})",
                               "matrix_traj_export.py", "inject", e, lat, grp,
                               on_done=lambda c: self._update_info()))

    def _shot(self) -> None:
        e = self._exp()
        if e is None:
            return
        args = ["--exp", e, "--view", self.view_mode.currentData(), "--span", self.span.value(),
                "--out", str(SHOT)]
        if self.mag.value() > 0:
            args += ["--mag", self.mag.value()]
        self.runner.run(py_job(f"снимок опыта {e}", "viz_shot.py", *args, on_done=self._show_shot))

    def _show_shot(self, code: int) -> None:
        if code != 0 or not SHOT.exists():
            self.preview.setText("снимок не получился — см. журнал")
            return
        pix = QPixmap(str(SHOT))
        if pix.isNull():
            self.preview.setText("снимок не читается")
            return
        self.preview.setPixmap(pix.scaled(self.preview.width(), self.preview.height(),
                                          Qt.KeepAspectRatio, Qt.SmoothTransformation))

    def _html_scene(self) -> None:
        e = self._exp()
        if e is None:
            return
        kind = self.html_kind.currentData()
        args = [e]
        if kind == "glyph":
            args.append("--glyph")
        elif kind == "both":
            args.append("--both")
        if self.decimate.value() >= 0:
            args += ["--decimate", self.decimate.value()]
        if self.name3d.text().strip():
            args += ["--out", self._out_path(self.name3d.text(), paths.HTML_DIR)]
        self.runner.run(py_job(f"HTML 3D-сцена опыта {e}", "viz_html.py", *args,
                               on_done=lambda c: self._update_info()))

    def _html_traj(self) -> None:
        e = self._exp()
        if e is None:
            return
        if not (paths.RESULTS / f"trajectories_local_exp{e}.npz").exists():
            self.say(f"нет trajectories_local_exp{e}.npz — прогоните трекинг киковым ядром "
                     f"в режиме «локальные координаты» или «длинный»", True)
            return
        args = [e, "--no-open"] + (["--csv"] if self.ck_csv.isChecked() else [])
        if self.name_traj.text().strip():
            args += ["--out", str(self._out_path(self.name_traj.text(), paths.RESULTS))]
        self.runner.run(py_job(f"HTML анализ траектории опыта {e}", "traj_analyzer.py", *args,
                               on_done=lambda c: self._update_info()))
