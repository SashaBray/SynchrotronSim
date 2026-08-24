# -*- coding: utf-8 -*-
"""Запуск инструментов проекта дочерними процессами с потоковым выводом в лог.

Одна очередь на всё приложение: команды выполняются последовательно (конвейер мастера, «матрицы →
отчёт»), интерфейс не блокируется, вывод приходит построчно по мере появления. Занятость очереди
транслируется сигналом `busy_changed`, по нему вкладки гасят свои кнопки — параллельный запуск двух
трекингов испортил бы общие таблицы БД.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from PySide6.QtCore import QObject, QProcess, Signal

from . import paths


@dataclass
class Job:
    """Одна команда: argv, человеческое имя для лога и необязательный колбэк по завершении."""
    argv: list[str]
    title: str
    cwd: str = str(paths.ROOT)
    on_done: object = None                      # callable(exit_code) | None
    show_command: bool = True
    _tag: dict = field(default_factory=dict)    # произвольные данные вызывающей стороны


def py_job(title: str, script: str, *args, **kw) -> Job:
    """Job для python-скрипта проекта: `py -3.14 <script> <args…>`."""
    argv = [paths.python_exe(), str(paths.ROOT / script)] + [str(a) for a in args if a is not None]
    return Job(argv=argv, title=title, **kw)


def exe_job(title: str, stem: str, *args, **kw) -> Job:
    """Job для C++ бинаря из build/."""
    argv = [str(paths.exe_path(stem))] + [str(a) for a in args if a is not None]
    return Job(argv=argv, title=title, **kw)


class Runner(QObject):
    """Последовательная очередь дочерних процессов."""

    line = Signal(str, bool)          # (текст, это_stderr)
    job_started = Signal(str)         # title
    job_finished = Signal(str, int)   # title, код возврата
    busy_changed = Signal(bool)
    queue_finished = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)
        self._queue: list[Job] = []
        self._cur: Job | None = None
        self._proc: QProcess | None = None
        self._stopping = False

    # ---- состояние -------------------------------------------------------
    @property
    def busy(self) -> bool:
        return self._cur is not None

    # ---- постановка задач ------------------------------------------------
    def run(self, *jobs: Job) -> None:
        """Поставить задачи в очередь (и запустить, если простаивает)."""
        self._queue.extend(jobs)
        if not self.busy:
            self._next()

    def stop(self) -> None:
        """Прервать текущий процесс и очистить очередь."""
        self._queue.clear()
        if self._proc is not None:
            self._stopping = True
            self.line.emit("[остановка по требованию пользователя]", True)
            self._proc.kill()

    # ---- внутреннее ------------------------------------------------------
    def _next(self) -> None:
        if not self._queue:
            self._cur = None
            self.busy_changed.emit(False)
            self.queue_finished.emit()
            return
        job = self._queue.pop(0)
        self._cur = job
        self._stopping = False
        self.busy_changed.emit(True)
        self.job_started.emit(job.title)

        proc = QProcess(self)
        proc.setProgram(job.argv[0])
        proc.setArguments(job.argv[1:])
        proc.setWorkingDirectory(job.cwd)
        env = proc.processEnvironment()
        for k, v in paths.child_env().items():
            env.insert(k, v)
        proc.setProcessEnvironment(env)
        proc.readyReadStandardOutput.connect(lambda: self._read(proc, False))
        proc.readyReadStandardError.connect(lambda: self._read(proc, True))
        proc.finished.connect(lambda code, _st: self._finished(code))
        proc.errorOccurred.connect(self._error)
        self._proc = proc

        if job.show_command:
            self.line.emit("$ " + " ".join(self._quote(a) for a in job.argv), False)
        proc.start()

    @staticmethod
    def _quote(a: str) -> str:
        return f'"{a}"' if " " in a else a

    def _read(self, proc: QProcess, err: bool) -> None:
        raw = (proc.readAllStandardError() if err else proc.readAllStandardOutput()).data()
        text = raw.decode("utf-8", errors="replace")
        for ln in text.splitlines():
            self.line.emit(ln, err)

    def _error(self, err) -> None:
        if err == QProcess.FailedToStart:
            self.line.emit(f"[не удалось запустить: {self._cur.argv[0] if self._cur else '?'}]", True)

    def _finished(self, code: int) -> None:
        job, self._cur, self._proc = self._cur, None, None
        if job is not None:
            if self._stopping:
                self.line.emit(f"[{job.title}: прервано]", True)
            else:
                self.line.emit(f"[{job.title}: код возврата {code}]", code != 0)
            self.job_finished.emit(job.title, code)
            if callable(job.on_done):
                try:
                    job.on_done(code)
                except Exception as e:                       # колбэк не должен ронять очередь
                    self.line.emit(f"[ошибка обработчика {job.title}: {type(e).__name__}: {e}]", True)
        if self._stopping:                                    # прервали — остаток очереди уже очищен
            self._stopping = False
        self._next()
