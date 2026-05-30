"""
Live 3D visualizer with auto-refresh on database changes.

Открывает PyVista-окно с решёткой и траекторией. Главный цикл периодически
поллит mtime/size файлов sp_elements.csv и hardware_library.csv. При любом
изменении:
  1. assemble_lattice_configs.py
  2. SynchrotronTracker.exe <exp_id>
  3. Перерисовка траектории (камера сохраняется).

Также можно ВРУЧНУЮ нажать 'r' или 'u' для немедленного refresh
(переопределяет встроенный PyVista "reset camera").

Запуск:
    python live_visualizer.py [exp_id]

Архитектура простая (один поток): plotter.show(interactive_update=True) ->
return immediately -> main loop вызывает plotter.update() (обработка событий,
рендер) + проверяет файлы + при изменении делает sim + refresh.
"""
import os
import sys
import json
import time
import subprocess
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import pyvista as pv

from visualize import build_combined_mesh, load_trajectories


SCRIPT_DIR = Path(__file__).resolve().parent
WS         = SCRIPT_DIR / "workspace"
EXE        = SCRIPT_DIR / "build" / "SynchrotronTracker.exe"
ASM        = SCRIPT_DIR / "assemble_lattice_configs.py"
SP_PATH    = WS / "tables" / "sp_elements.csv"
HW_PATH    = WS / "tables" / "hardware_library.csv"

# Все настройки -- из workspace/configs/visualizer.json
VIZ_CONFIG_PATH = WS / "configs" / "visualizer.json"
VIZ_CFG = {}
if VIZ_CONFIG_PATH.exists():
    try:
        with open(VIZ_CONFIG_PATH, encoding="utf-8") as f:
            VIZ_CFG = json.load(f)
    except Exception as e:
        print(f"  (warning: cannot read {VIZ_CONFIG_PATH}: {e})")

TRAJ_COLORS = VIZ_CFG.get("particle_colors",
                          ["magenta", "cyan", "lime", "yellow", "red"])
_lbls = VIZ_CFG.get("particle_labels",
                    ["P1 (ref)", "P2 (+Y)", "P3 (-Y)", "P4 (+Z)", "P5 (-Z)"])
TRAJ_LABELS = {i: lbl for i, lbl in enumerate(_lbls)}

_LV = VIZ_CFG.get("live_visualizer", {})
POLL_INTERVAL   = float(_LV.get("poll_interval_seconds", 0.5))
LOOP_SLEEP      = float(_LV.get("loop_sleep_seconds",    0.05))
HEARTBEAT_EVERY = int(_LV.get("heartbeat_every_ticks",    60))
MANUAL_KEYS     = _LV.get("manual_refresh_keys", ["r", "u"])


def stat_key(path):
    """(mtime, size) -- robust file-change signature."""
    try:
        s = path.stat()
        return (s.st_mtime, s.st_size)
    except FileNotFoundError:
        return (0.0, 0)


class LiveVisualizer:
    def __init__(self, exp_id):
        self.exp_id   = exp_id
        self.traj_bin = WS / "results" / f"trajectories_exp{exp_id}.bin"

        self.last_sp_key = stat_key(SP_PATH)
        self.last_hw_key = stat_key(HW_PATH)

        self.traj_actors = []
        self.plotter = None
        self._manual_refresh_requested = False

    # =========================================================================
    # Build lattice scene (once at startup)
    # =========================================================================
    def build_lattice_scene(self):
        ws = str(WS) + os.sep
        exp_df     = pd.read_csv(ws + "tables/experiments.csv",      sep=";", skipinitialspace=True)
        lattice_df = pd.read_csv(ws + "tables/lattice_configs.csv",  sep=";", skipinitialspace=True)
        hw_df      = pd.read_csv(ws + "tables/hardware_library.csv", sep=";", skipinitialspace=True)

        exp_row = exp_df[exp_df["exp_id"] == self.exp_id]
        if exp_row.empty:
            raise RuntimeError(f"Experiment {self.exp_id} not found")
        lat_id = int(exp_row["lattice_id"].iloc[0])
        current_lat = lattice_df[lattice_df["lattice_id"] == lat_id]
        print(f"  lattice {lat_id}: {len(current_lat)} instances")

        hw_lookup = {}
        for row in hw_df.itertuples():
            dev_id = str(getattr(row, "device_type_id", "")).strip()
            if not dev_id or dev_id == "nan": continue
            stl_name  = str(getattr(row, "stl_path", "")).strip()
            is_active = str(getattr(row, "Active", "")).strip().lower() in ("true", "1", "yes")
            hw_lookup[dev_id] = (stl_name, is_active)

        stl_cache = {}
        for stl_name, _ in hw_lookup.values():
            if stl_name in ("-", "", "nan") or stl_name in stl_cache: continue
            path = WS / "geometry" / stl_name
            if path.exists():
                m = pv.read(str(path)); m.scale(0.001, inplace=True)
                stl_cache[stl_name] = m
            else:
                stl_cache[stl_name] = None

        dev_ids = current_lat["device_type_id"].astype(str).str.strip().values
        yaws    = current_lat["yaw"].astype(float).values
        pitches = current_lat["pitch"].astype(float).values
        rolls   = current_lat["roll"].astype(float).values
        xs      = current_lat["x"].astype(float).values
        ys      = current_lat["y"].astype(float).values
        zs      = current_lat["z"].astype(float).values

        grouped = defaultdict(list)
        for i in range(len(dev_ids)):
            dev = dev_ids[i]
            if dev not in hw_lookup: continue
            stl_name, _ = hw_lookup[dev]
            if stl_cache.get(stl_name) is None: continue
            grouped[dev].append((yaws[i], pitches[i], rolls[i], xs[i], ys[i], zs[i]))

        for dev, tlist in grouped.items():
            base = stl_cache[hw_lookup[dev][0]]
            big = build_combined_mesh(base, np.asarray(tlist, dtype=np.float64))
            is_active = hw_lookup[dev][1]
            op = VIZ_CFG.get("lattice_opacity", {})
            opacity = float(op.get("EM", 0.25) if is_active else op.get("PM", 0.4))
            self.plotter.add_mesh(
                big, color="silver",
                opacity=opacity,
                show_edges=True,
                label=f"{dev} ({'EM' if is_active else 'PM'}, x{len(tlist)})",
            )

    # =========================================================================
    # Trajectory refresh (UI only -- быстро)
    # =========================================================================
    def refresh_trajectory_ui(self):
        for a in self.traj_actors:
            try: self.plotter.remove_actor(a, render=False)
            except Exception: pass
        self.traj_actors.clear()

        trajs = load_trajectories(str(self.traj_bin))
        if not trajs:
            self.plotter.render()
            return

        for i, points in enumerate(trajs):
            if len(points) < 2: continue
            poly = pv.PolyData(points)
            cells = np.full((1, len(points) + 1), len(points), dtype=np.int_)
            cells[0, 1:] = np.arange(len(points))
            poly.lines = cells
            actor = self.plotter.add_mesh(
                poly,
                color=TRAJ_COLORS[i % len(TRAJ_COLORS)],
                line_width=4,
                label=TRAJ_LABELS.get(i, f"P{i+1}"),
                render=False,
            )
            self.traj_actors.append(actor)
        self.plotter.render()

    # =========================================================================
    # Simulation run (assemble + tracker)
    # =========================================================================
    def run_simulation(self):
        t0 = time.time()
        r = subprocess.run([sys.executable, str(ASM)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  ASSEMBLE FAILED: {r.stderr.strip()[:200]}")
            return False
        r = subprocess.run([str(EXE), str(self.exp_id)], capture_output=True, text=True)
        if r.returncode != 0:
            print(f"  TRACKER FAILED: {r.stderr.strip()[:200]}")
            return False
        print(f"  sim done in {time.time()-t0:.2f}s")
        return True

    def request_manual_refresh(self):
        """Called by 'r'/'u' key. Просто ставим флаг, основной цикл подхватит."""
        self._manual_refresh_requested = True
        print("\n[manual] refresh requested via key press")

    # =========================================================================
    # Главный цикл
    # =========================================================================
    def run(self):
        print(f"\n{'='*60}\nLIVE VISUALIZER (exp {self.exp_id})\n{'='*60}")
        print(f"  Watching: {SP_PATH.name}, {HW_PATH.name}")
        print(f"  Poll interval: {POLL_INTERVAL}s")
        print(f"  Hotkeys: 'r' / 'u' = manual refresh, 'q' = quit")

        self.plotter = pv.Plotter(title=f"Live Synchrotron (exp {self.exp_id})")
        self.plotter.set_background(VIZ_CFG.get("background", "white"))

        print("\nBuilding lattice scene ...")
        t0 = time.time()
        self.build_lattice_scene()
        print(f"  ({time.time()-t0:.2f}s)")

        self.refresh_trajectory_ui()

        self.plotter.show_bounds(
            grid="back", location="outer", all_edges=False,
            xtitle="X (m)", ytitle="Y (m)", ztitle="Z (m)",
            font_size=10, color="black",
        )
        self.plotter.add_axes()
        try: self.plotter.add_legend(bcolor=None)
        except Exception: pass

        cam = VIZ_CFG.get("camera", {})
        self.plotter.camera_position = [
            tuple(cam.get("position", [7.0, 0.5, 15.0])),
            tuple(cam.get("focus",    [7.0, 0.5,  0.0])),
            tuple(cam.get("up",       [0.0, 1.0,  0.0])),
        ]
        # Параллельная (ортографическая) проекция -- по умолчанию ВКЛ.
        if cam.get("parallel_projection", True):
            self.plotter.enable_parallel_projection()

        # 'r' в PyVista по умолчанию = "reset camera". Очищаем и вешаем своё.
        for key in MANUAL_KEYS:
            try: self.plotter.clear_events_for_key(key)
            except Exception: pass
            self.plotter.add_key_event(key, self.request_manual_refresh)

        # Открываем окно НЕблокирующе
        self.plotter.show(interactive_update=True, auto_close=False)

        print("\n--- 3D window opened. Save CSV in editor to refresh. ---\n")

        # Основной цикл -- управляется ЭТИМ потоком, никаких таймеров VTK
        last_poll = time.time()
        tick = 0
        try:
            while not self._is_window_closed():
                # 1) Обработать события окна (мышь/клавиатура) и перерисовать
                try:
                    self.plotter.update()
                except Exception:
                    break

                # 2) Heartbeat
                tick += 1
                if tick == 1:
                    print(f"  [main-loop] running (heartbeat every ~3s)")
                elif tick % HEARTBEAT_EVERY == 0:
                    print(f"  [heartbeat tick {tick}, ~{tick*LOOP_SLEEP:.0f}s] watching files...")

                # 3) Проверка файлов раз в POLL_INTERVAL
                now = time.time()
                need_refresh = False
                source = ""
                if now - last_poll >= POLL_INTERVAL:
                    last_poll = now
                    sp_key = stat_key(SP_PATH)
                    hw_key = stat_key(HW_PATH)
                    if sp_key != self.last_sp_key:
                        source = f"sp_elements (mt {self.last_sp_key[0]:.1f}->{sp_key[0]:.1f}, sz {self.last_sp_key[1]}->{sp_key[1]})"
                        self.last_sp_key = sp_key
                        need_refresh = True
                    if hw_key != self.last_hw_key:
                        s2 = f"hardware_library (mt {self.last_hw_key[0]:.1f}->{hw_key[0]:.1f}, sz {self.last_hw_key[1]}->{hw_key[1]})"
                        source = source + " + " + s2 if source else s2
                        self.last_hw_key = hw_key
                        need_refresh = True

                # 4) Manual refresh request
                if self._manual_refresh_requested:
                    self._manual_refresh_requested = False
                    need_refresh = True
                    if not source:
                        source = "manual key press"
                        # При manual также обновим last_keys, чтобы watcher не дублировал
                        self.last_sp_key = stat_key(SP_PATH)
                        self.last_hw_key = stat_key(HW_PATH)

                # 5) Если нужно -- запускаем симуляцию + рисуем
                if need_refresh:
                    print(f"\n[refresh] {source}")
                    if self.run_simulation():
                        self.refresh_trajectory_ui()
                        print("  trajectory redrawn")

                # 6) Короткий sleep
                time.sleep(LOOP_SLEEP)
        finally:
            try: self.plotter.close()
            except Exception: pass
            print("\n[live_visualizer] closed")

    def _is_window_closed(self):
        """Various heuristics to detect plotter window closure."""
        if getattr(self.plotter, "_closed", False):
            return True
        try:
            if self.plotter.render_window is None:
                return True
        except Exception:
            return True
        try:
            # VTK interactor reports termination
            if self.plotter.iren.interactor.GetDone():
                return True
        except Exception:
            pass
        return False


def main():
    exp_id = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    LiveVisualizer(exp_id).run()


if __name__ == "__main__":
    main()
