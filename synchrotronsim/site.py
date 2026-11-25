# -*- coding: utf-8 -*-
"""Площадка: где хранилище и где вычислитель (docs/DATA_PROTOCOL.md, этап 3, раздел 7).

    from synchrotronsim import Site
    lab = Site.load("lab")                     # профиль: адрес сервера хранилища и вычислителя
    print(lab.check())                         # связь, задержка, версии протокола, сборка ядра на сервере
    exp = lab.store.get("experiment", "fofb_demo")
    run = Simulation(exp, on=lab).run()        # счёт на вычислителе площадки (этап 3b)

`Site.local()` (по умолчанию) — всё в этом процессе, как раньше. Профили — `~/.synchrotronsim/sites.json`
(адреса, без секретов; копию можно положить в хранилище: `publish`); токены — `~/.synchrotronsim/credentials.json`
или переменная окружения `SS_TOKEN_<ИМЯ>` / `SS_TOKEN`; в хранилище и в репозиторий токены не попадают.
"""
from __future__ import annotations

import json
import os
import statistics
import time
from dataclasses import dataclass, field, asdict
from typing import Any, Dict, List, Optional

HOME = os.path.join(os.path.expanduser("~"), ".synchrotronsim")
SITES_FILE = os.path.join(HOME, "sites.json")
CRED_FILE = os.path.join(HOME, "credentials.json")


def _read_json(path: str) -> Dict[str, Any]:
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def _write_json(path: str, d: Dict[str, Any], private: bool = False) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as f:
        json.dump(d, f, ensure_ascii=False, indent=1, sort_keys=True)
    if private:
        try:
            os.chmod(tmp, 0o600)
        except OSError:
            pass
    os.replace(tmp, path)


@dataclass
class Check:
    """Итог проверки одной службы площадки."""
    service: str
    url: str
    ok: bool
    latency_ms: Optional[float] = None
    info: Dict[str, Any] = field(default_factory=dict)
    problems: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)


@dataclass
class SiteReport:
    site: str
    checks: List[Check]

    @property
    def ok(self) -> bool:
        return all(c.ok for c in self.checks)

    def __str__(self):
        lines = [f"площадка «{self.site}»: {'всё в порядке' if self.ok else 'ЕСТЬ ПРОБЛЕМЫ'}"]
        for c in self.checks:
            head = f"  {c.service:12s} {c.url or '(в этом процессе)'}  " + ("OK" if c.ok else "НЕТ")
            if c.latency_ms is not None:
                head += f", задержка {c.latency_ms:.1f} мс"
            lines.append(head)
            i = c.info
            if i.get("host"):
                lines.append(f"      сервер {i.get('host')} · {i.get('os')} · Python {i.get('python')} · ядер {i.get('cpus')}")
            for n in c.notes:
                lines.append(f"      {n}")
            for p in c.problems:
                lines.append(f"      ! {p}")
        return "\n".join(lines)


def _local_engine():
    from .run import provenance
    from .db import PROJECT_ROOT
    return ((provenance(PROJECT_ROOT) or {}).get("engine") or {}).get("extension")


def _compare_engine(c: Check, remote_engine) -> None:
    loc = None
    try:
        loc = _local_engine()
    except Exception:
        pass
    if not remote_engine:
        c.notes.append("ядро на сервере: расширение не собрано (счёт там невозможен)")
        return
    if loc and remote_engine.get("sha16") == loc.get("sha16"):
        c.notes.append(f"ядро: та же сборка, что здесь ({remote_engine.get('build_info')}) — результаты бит-в-бит")
    else:
        c.notes.append(f"ядро на сервере: {remote_engine.get('build_info')}; здесь: "
                       f"{(loc or {}).get('build_info', 'нет')} — другая сборка: бит-в-бит с этим ПК не "
                       f"гарантирован (сверять с допуском; паспорт прогона пишет, где считалось)")


@dataclass
class Site:
    """Профиль площадки: адрес сервера хранилища и вычислителя ("" — в этом процессе)."""
    name: str = "local"
    store_url: str = ""
    compute_url: str = ""
    note: str = ""

    def __post_init__(self):
        self._store = None
        self._compute = None

    # ---- профили
    @classmethod
    def local(cls) -> "Site":
        return cls("local")

    @classmethod
    def load(cls, name: str, path: Optional[str] = None) -> "Site":
        if name in ("", "local", None):
            return cls.local()
        d = _read_json(path or SITES_FILE).get(name)
        if d is None:
            raise KeyError(f"нет площадки «{name}» в {path or SITES_FILE} — добавьте: ss_server.py site add {name} …")
        return cls.from_dict(dict(d, name=name))

    @classmethod
    def from_dict(cls, d: Dict[str, Any]) -> "Site":
        return cls(name=d.get("name", "local"), store_url=d.get("store_url", ""), compute_url=d.get("compute_url", ""),
                   note=d.get("note", ""))

    def to_dict(self) -> Dict[str, Any]:
        return {k: v for k, v in asdict(self).items() if not k.startswith("_")}

    def save(self, path: Optional[str] = None) -> str:
        p = path or SITES_FILE
        d = _read_json(p)
        d[self.name] = {k: v for k, v in self.to_dict().items() if k != "name"}
        _write_json(p, d)
        return p

    @staticmethod
    def list_sites(path: Optional[str] = None) -> Dict[str, Dict[str, Any]]:
        return _read_json(path or SITES_FILE)

    @staticmethod
    def remove(name: str, path: Optional[str] = None) -> None:
        p = path or SITES_FILE
        d = _read_json(p)
        d.pop(name, None)
        _write_json(p, d)

    def publish(self, store) -> str:
        """Положить профиль (без токена) в хранилище — документ вида «site»."""
        return store.put_json("site", self.name, self.to_dict(), store.head("site", self.name), "профиль площадки")

    @classmethod
    def from_store(cls, store, name: str) -> "Site":
        return cls.from_dict(store.json_content("site", name))

    # ---- токены (только локально)
    def token(self) -> Optional[str]:
        env = os.environ.get(f"SS_TOKEN_{self.name.upper().replace('-', '_')}") or os.environ.get("SS_TOKEN")
        return env or _read_json(CRED_FILE).get(self.name) or None

    @staticmethod
    def set_token(name: str, token: str) -> str:
        d = _read_json(CRED_FILE)
        d[name] = token
        _write_json(CRED_FILE, d, private=True)
        return CRED_FILE

    # ---- клиенты
    @property
    def store(self):
        """Хранилище площадки: локальное (open_store) или клиент сервера с рабочей копией этого проекта."""
        if self._store is None:
            if not self.store_url:
                from .store import open_store
                self._store = open_store()
            else:
                from .remote import RemoteStore
                from .db import PROJECT_ROOT
                wc = PROJECT_ROOT if os.path.isdir(os.path.join(PROJECT_ROOT, "workspace", "tables")) else None
                self._store = RemoteStore(self.store_url, self.token(), working=wc, site_name=self.name)
        return self._store

    @property
    def compute(self):
        """Вычислитель площадки (клиент сервера вычислений) или None — счёт в этом процессе."""
        if self._compute is None and self.compute_url:
            from .compute import RemoteCompute
            self._compute = RemoteCompute(self.compute_url, self.token(), site=self)
        return self._compute

    # ---- проверка связи
    def check(self, pings: int = 3) -> SiteReport:
        checks = []
        if not self.store_url:
            c = Check("хранилище", "", True)
            try:
                st = self.store
                c.notes.append(f"локальное: {st!r}")
            except Exception as e:                              # noqa: BLE001
                c.ok = False; c.problems.append(f"{type(e).__name__}: {e}")
            checks.append(c)
        else:
            checks.append(self._check_store(pings))
        if self.compute_url:
            checks.append(self._check_compute(pings))
        else:
            checks.append(Check("вычислитель", "", True, notes=["счёт в этом процессе"]))
        return SiteReport(self.name, checks)

    def _check_store(self, pings: int) -> Check:
        from .remote import ping, RemoteStore
        from .wire import PROTOCOL
        c = Check("хранилище", self.store_url, False)
        try:
            p = ping(self.store_url)
        except Exception as e:                                  # noqa: BLE001
            c.problems.append(f"не отвечает: {type(e).__name__}: {e}")
            return c
        if p.get("service") != "store":
            c.problems.append(f"по адресу не сервер хранилища ({p})")
            return c
        if p.get("protocol") != PROTOCOL:
            c.problems.append(f"протокол сервера {p.get('protocol')} ≠ {PROTOCOL} (обновите код на одной из сторон)")
        lat = []
        for _ in range(max(1, pings)):
            t0 = time.perf_counter(); ping(self.store_url); lat.append((time.perf_counter() - t0) * 1e3)
        c.latency_ms = statistics.median(lat)
        rs = RemoteStore(self.store_url, self.token(), site_name=self.name)
        try:
            t0 = time.time(); info = rs.info(); t1 = time.time()
        except Exception as e:                                  # noqa: BLE001
            c.problems.append(f"нет доступа: {e}" if "токен" in str(e) else f"{type(e).__name__}: {e}")
            return c
        c.info = info
        skew = float(info.get("time", t0)) - 0.5 * (t0 + t1)
        if abs(skew) > 2.0:
            c.problems.append(f"часы сервера расходятся с этим ПК на {skew:+.1f} с (даты ревизий)")
        c.notes.append(f"хранилище {info.get('store')}, машин {info.get('machines', '?')}")
        code = (info.get("code") or {})
        if code:
            c.notes.append(f"код на сервере: {str(code.get('sha') or '?')[:10]} ({code.get('branch', '?')})"
                           + (" — есть незафиксированные правки" if code.get("dirty") else ""))
        _compare_engine(c, info.get("engine"))
        c.ok = not c.problems
        return c

    def _check_compute(self, pings: int) -> Check:
        c = Check("вычислитель", self.compute_url, False)
        try:
            from .compute import RemoteCompute
            rc = RemoteCompute(self.compute_url, self.token(), site=self)
            lat = []
            for _ in range(max(1, pings)):
                t0 = time.perf_counter(); rc.ping(); lat.append((time.perf_counter() - t0) * 1e3)
            c.latency_ms = statistics.median(lat)
            t0 = time.time(); info = rc.info(); t1 = time.time()
            c.info = info
            skew = float(info.get("time", t0)) - 0.5 * (t0 + t1)
            if abs(skew) > 2.0:
                c.problems.append(f"часы вычислителя расходятся с этим ПК на {skew:+.1f} с")
            st = info.get("store_url") or "своё локальное"
            c.notes.append(f"хранилище вычислителя: {st}; заданий в очереди {info.get('queue', 0)}")
            if self.store_url and info.get("store_url") and info["store_url"].rstrip("/") != self.store_url.rstrip("/"):
                c.notes.append(f"вычислитель смотрит в {info['store_url']}, а площадка — в {self.store_url}: ссылки на "
                               f"ревизии должны быть доступны вычислителю")
            _compare_engine(c, info.get("engine"))
        except ImportError:
            c.problems.append("клиент вычислителя ещё не подключён (этап 3b)")
            return c
        except Exception as e:                                  # noqa: BLE001
            c.problems.append(f"не отвечает: {type(e).__name__}: {e}")
            return c
        c.ok = not c.problems
        return c
