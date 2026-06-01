# Synchrotron Particle Tracker

Relativistic charged-particle tracking through a synchrotron magnet lattice.
Uses **real Ansys field maps**, a **symplectic Boris pusher** for integration,
and a **stochastic optimizer** for tuning magnet currents and per-instance
calibration coefficients.

Cross-platform: builds and runs on Windows (MSVC) and Linux/macOS (gcc/clang).

GitHub: <https://github.com/SashaBray/SynchrotronSim>

---

## Contents

1. [Quick start](#quick-start)
2. [Overview](#overview)
3. [Repository layout](#repository-layout)
4. [Physical model](#physical-model)
5. [Database tables](#database-tables)
6. [Configuration files](#configuration-files)
7. [Pipeline and experiments](#pipeline-and-experiments)
8. [Optimizer loss function](#optimizer-loss-function)
9. [Build and run](#build-and-run)
10. [Diagnostic scripts](#diagnostic-scripts)
11. [Cross-platform notes](#cross-platform-notes)
12. [Performance](#performance)
13. [Troubleshooting](#troubleshooting)

---

## Quick start

Assuming the project is already built and field maps are in place:

```bash
# 1. Run experiment 4 (SP3 exit focusing) and view 3D trajectories
build/SynchrotronTracker 4         # Windows: build\SynchrotronTracker.exe 4
python visualize.py 4

# 2. Check trajectory stability (envelope + FFT)
python analyze_trajectory_stability.py 4

# 3. Optimize currents (~80 min for 10000 iterations on exp 4)
python prepare_optim_config.py workspace/optim_inputs/sp3_drift_2sp.json
build/SynchrotronOptimizer workspace/optim_inputs/sp3_drift_2sp.optimconfig

# 4. Inspect convergence and re-render with new state
python viz_optim_history.py workspace/results/sp3_drift_2sp.history.json
python assemble_lattice_configs.py
build/SynchrotronTracker 4
python visualize.py 4
```

After editing **any** DB table by hand (sp_elements, particle_groups, LTE files)
re-run `python assemble_lattice_configs.py` so the simulator picks up the
changes. If LTE has changed, run `python lte_to_lattice.py` first.

---

## Overview

**What is simulated:** 3D trajectories of charged particles in the synchrotron
magnet ring. Each magnet is represented by a precomputed field on a Cartesian
grid (Ansys output, converted to `.bin`), to which a geometric transform
(position, orientation) and calibration multipliers (`correct_coeff`,
`field_direction`, current scaling) are applied at runtime.

**Two executables:**

| Binary | Purpose |
|---|---|
| `SynchrotronTracker` | forward tracking of N particles through a chosen lattice |
| `SynchrotronOptimizer` | stochastic hill-climber with periodic kicks; finds magnet currents and calibrations minimizing a composite loss |

**Python tooling** (top level, all cross-platform):

| Script | Purpose |
|---|---|
| `lte_to_lattice.py` | LTE → `sp_elements.csv` + `lattice_layout.csv` |
| `assemble_lattice_configs.py` | sp_elements + layout → `lattice_configs.csv` (+ sub-lattices) |
| `prepare_optim_config.py` | optim JSON → flat `.optimconfig` |
| `apply_optim_history.py` | history.json → CSV (state rollback / checkpoint restore) |
| `visualize.py` | static 3D viewer (PyVista) |
| `live_visualizer.py` | live-reload viewer (auto-tracker on DB change) |
| `viz_optim_history.py` | matplotlib plots of loss / parameter trajectory |
| `analyze_trajectory_stability.py` | envelope analysis + FFT betatron tunes |
| `find_orbit_deviation_peaks.py` | locate magnets near orbit deviation maxima |
| `scan_field_maps.py` | detect / patch dropouts in `.bin` field maps |
| `_common.py` | shared helpers (`tracker_exe`, `project_root`, matplotlib backend) |

---

## Repository layout

```text
SynchrotronSim/
├── src/                                # C++ source
│   ├── main.cpp                        # SynchrotronTracker entry
│   ├── optimizer.cpp                   # SynchrotronOptimizer entry
│   ├── Integrator.h                    # Boris pusher (symplectic relativistic)
│   ├── LatticeManager.h                # magnets, transforms, spatial bucket index
│   ├── FieldMap.h                      # trilinear interp on binary maps
│   ├── FinishPlane.h                   # sign-change finish detection
│   └── AnsysConverter.h                # .fld -> .bin converter
│
├── workspace/                          # data root
│   ├── configs/                        # JSON configuration (hand-edited)
│   │   ├── lte_mapping.json            # LTE type -> hardware_library mapping
│   │   ├── sub_lattices.json           # sub-lattice definitions (102/103/104)
│   │   └── visualizer.json             # camera, colors, opacity for PyVista
│   ├── lte/                            # Elegant LTE sources
│   │   ├── SYLA_desc.lte               # device descriptions
│   │   └── SYLA_line.lte               # element sequence
│   ├── tables/                         # database (`;`-separated CSV)
│   │   ├── hardware_library.csv        # device_type_id -> field map + calibration
│   │   ├── sp_elements.csv             # SP templates (generated, then edited)
│   │   ├── lattice_layout.csv          # full ring unrolling (generated)
│   │   ├── lattice_configs.csv         # final simulator input (generated)
│   │   ├── experiments.csv             # exp_id; lat; group; dt; max_turns; stride; p1..p3
│   │   └── particle_groups.csv         # 5 particles per group (centre + 4 satellites)
│   ├── optim_inputs/                   # optimizer setups (one JSON per experiment)
│   │   ├── sp1_endpoint.json           # exp 2 (legacy: end of SP1)
│   │   ├── sp3_drift_2sp.json          # exp 4: SP3 + 2 SP free drift
│   │   ├── sp10_drift_to_sp11.json     # exp 5: SP1-10 + drift to SP11
│   │   └── junction_tune.json          # exp 6: junction tuning (SP20)
│   ├── fields_raw/                     # raw .fld files from Ansys (large, gitignored)
│   ├── fields_bin/                     # compiled .bin field maps (generated)
│   ├── geometry/                       # STL models (large, gitignored)
│   └── results/                        # outputs: trajectories_exp*.bin, *.history.json
│
├── docs/                               # documentation
│   ├── algorithms.md                   # algorithm reference (RD 50-34.698-90)
│   ├── algorithms.docx                 # Russian GOST-formatted version
│   └── build_docx.js                   # docx generator
│
├── CMakeLists.txt                      # cmake config (MSVC / gcc / clang)
├── _common.py                          # shared Python helpers (cross-platform)
├── README.md                           # this file
├── README.txt                          # plain-text mirror
└── .gitignore
```

---

## Physical model

### Integrator: relativistic Boris pusher

Equations of motion:

```
dx/dt = v
dp/dt = q · (E + v × B)
p     = γ · m · v,   γ = √(1 + p²/(m²c²))
```

The Boris algorithm splits one step into **drift–kick–drift**:

1. half-drift: `x ← x + (p/(γ·m))·dt/2`
2. sample fields: `B, E ← lattice.getTotalFields(x_{n+½})`
3. half E-kick: `p⁻ = p + q·E·dt/2`
4. magnetic rotation: `t = q·B·dt/(2·γ⁻·m)`, `p⁺ = p⁻ + (p⁻ + p⁻×t) × 2t/(1+|t|²)`
5. half E-kick: `p_{n+1} = p⁺ + q·E·dt/2`
6. half-drift: `x ← x + (p/(γ·m))·dt/2`

Properties:

- **Symplectic** (preserves phase-space volume).
- **Exactly conserves `|p|²`** for purely magnetic fields (the rotation step
  is mathematically orthogonal).
- **One field eval per step** (RK4 needs four).

In the current lattice all devices are magnetic, so `E ≡ 0` and the
E-kick steps reduce to no-ops. The code keeps the general form anyway.

### Field at a query point

For one magnet, the field is computed via **trilinear interpolation** in its
local coordinate system:

1. world → local: `transform.transformPointInv(p)`
2. trilinear interpolation on the grid (`FieldMap::getInterpolatedField`)
3. optional YZ-mirror: when `useSymmetry == true` and `local.x < 0`, sample at
   `(-local.x, local.y, local.z)` and flip `F.x`
4. local → world: `transform.rotateVector(F)`
5. scale by `correctCoeff · fieldDirection · currentScale`

The lattice sums contributions over a **2D bucket grid in xy** (the ring is
flat in z ≈ 0). Buckets of size `max(1 m, 2 × max(worldRadius))` guarantee
that each magnet registers in at most 4 buckets. At each particle step the
inner loop ranges over only the ~5 magnets in the bucket holding the current
position; drift regions hit an empty bucket and the inner loop is empty.

A bounding-sphere prescreen rejects out-of-range elements with 7 flops before
the full local-frame transform.

---

## Database tables

All CSV files in `workspace/tables/` use `;` as separator (so European
decimal commas like `1,00E-11` in `dt` are unambiguous).

### `hardware_library.csv` (8 columns, hand-edited)

| Column | Description |
|---|---|
| `device_type_id` | unique magnet-type ID (`QP_L_1A`, `DQ_I1_1A`, …) |
| `Active` | `True` for an active electromagnet, `False` for a permanent magnet |
| `map_raw_path` | filename of the Ansys `.fld` source (under `workspace/fields_raw/`) |
| `map_bin_path` | filename of the compiled `.bin` (under `workspace/fields_bin/`) |
| `stl_path` | STL model for the visualizer (`-` if none) |
| `nominal_arg` | nominal current the field map was computed at (`-` for PM) |
| `field_direction` | sign multiplier `±1` |
| `correct_coeff` | device-level field calibration |

### `sp_elements.csv` (templates, hand-edited)

Defines per-superperiod (SP) magnet content, indexed by `(sp_type,
sp_element_id)`.

| Column | Description |
|---|---|
| `sp_type` | `standard` or `junction` |
| `sp_element_id` | 1..N inside the SP |
| `device_type_id` | reference into `hardware_library.csv` |
| `arg_val` | current in Amperes, or `-` for permanent magnets |
| `correct_coeff` | per-instance field multiplier (default 1.0) |

### `lattice_layout.csv` (generated)

Full ring unrolling — position and orientation per element.

| Column | Description |
|---|---|
| `global_id` | 1..2762 |
| `sp_number` | 1..39 superperiod number |
| `sp_type` | `standard` / `junction` |
| `sp_element_id` | reference into `sp_elements` |
| `x, y, z` | element position (m) |
| `yaw, pitch, roll` | orientation (rad) |

### `lattice_configs.csv` (generated)

Final input to the simulator. The assembler joins `sp_elements` and
`lattice_layout` (substituting `arg_val` and `correct_coeff`) and appends
copies for each sub-lattice defined in `sub_lattices.json`.

### `experiments.csv` (hand-edited)

| Column | Description |
|---|---|
| `exp_id` | unique experiment ID |
| `lattice_id` | 101 (full ring), 102+ (sub-lattices) |
| `group_id` | reference to `particle_groups.csv` |
| `dt` | integration timestep (s) |
| `max_turns` | finish-plane `skipSteps` (calls to ignore before first check) |
| `traj_stride` | save 1 of every N trajectory points |
| `p1x..p3z` | three points defining the finish plane (normal via cross product) |

### `particle_groups.csv` (hand-edited)

| Column | Description |
|---|---|
| `group_id` | group ID (1, 2, 3, …) |
| `part_id` | 2..6 (5 particles per group) |
| `m`, `q` | mass (kg), charge (C) |
| `x, y, z` | initial position (m) |
| `vx, vy, vz` | initial velocity (m/s) |

Convention: 5 particles per group — 1 reference on-axis + 4 satellites at
±1 mm in local y and z, used to estimate betatron oscillations.

---

## Configuration files

`workspace/configs/`:

### `lte_mapping.json`
- `lattice_id` of the main ring (101)
- LTE type → `device_type_id` mapping
- regular expressions for DL / JL element families
- classification thresholds (lengths, K1)

### `sub_lattices.json`
List of sub-lattices automatically generated by
`assemble_lattice_configs.py`. Each entry is
`{lattice_id, max_gid, description}`. Current defaults: 102 (SP1-3),
103 (SP1-10), 104 (SP1-21).

### `visualizer.json`
Camera (position, focus, up), parallel-projection flag, particle colours and
labels (5 entries), magnet opacity (EM and PM), live-visualizer poll
interval, hotkey bindings.

### `workspace/optim_inputs/*.json`
One JSON per optimization scenario. Knobs:
- `experiment_id`, `lattice_max_gid`, `n_iterations`
- `kick_every`, `kick_scale`
- `loss.magnet_sp_numbers`, `loss.exclude_device_types`,
  `loss.extra_target_magnets[]`
- `loss.direction_loss.reference`, `loss.direction_loss.weight`
- `loss.envelope_loss.weight`
- `parameter_groups_per_position[]` — what to mutate and how

---

## Pipeline and experiments

```
SYLA_desc.lte + SYLA_line.lte
    │   (lte_to_lattice.py)
    │   └─ DL / JL straightening (5 slices on chord)
    ▼
sp_elements.csv  +  lattice_layout.csv
    │   (assemble_lattice_configs.py)
    ▼
lattice_configs.csv  (includes sub-lattices 102 / 103 / 104)
    │
    ├─→ SynchrotronTracker N  →  trajectories_expN.bin
    │       │   (visualize.py N)
    │       ▼
    │   PyVista 3D window
    │
    └─→ SynchrotronOptimizer optim.optimconfig
            │   reads sp_elements.csv & hardware_library.csv,
            │   mutates them in memory, writes them back at the end,
            │   then re-runs the tracker once with the best state
            ▼
        sp_elements.csv  (updated with best arg_val)
        hardware_library.csv  (updated with best correct_coeff)
        results/{name}.history.json  (for viz_optim_history / apply_optim_history)
```

### Experiment list (`experiments.csv`)

| exp | Lattice | Purpose |
|---|---|---|
| 1 | 101 | Sanity test: plane at `x=5`, particles cross SP1 |
| 2 | 101 | Half ring — past SP20 (yaw ≈ 180°, `skipSteps = 100k`) |
| 3 | 101 | **Full ring** — past SP1 after 360° (`skipSteps = 200k`) |
| 4 | 102 (SP1-3) | Focusing at SP3 exit + 2 SP of free drift |
| 5 | 103 (SP1-10) | Long path SP1-SP10 + drift to SP11 |
| 6 | 104 (SP1-21) | Junction tuning (SP20), particles start at SP18 first |

---

## Optimizer loss function

Composite loss with **three additive components**:

### 1. Position loss (geometry)
```
Σ_(magnet m, particle i)  w_m · d²_{m,i}
```
- `d_{m,i}` — minimum distance from particle *i*'s trajectory to magnet *m*
  centre.
- Computed **incrementally** in the tracking loop — no trajectory storage.
- Magnets filtered by `magnet_sp_numbers` and `exclude_device_types`.
- Virtual `TARGET` points may be added with custom weights.

### 2. Direction loss (parallelism)
```
w_dir · Σ_i sin²(θ_i)
```
- `θ_i` — angle between particle *i* final momentum and a reference
  direction (orbit tangent at the finish plane).
- `sin²θ ≈ θ²` for small angles.

### 3. Envelope loss (betatron amplitude growth)
```
w_env · Σ_i max(0, A_i_final − A_i_init)²
```
- `A_i(t) = distance(particle_pos, design_polyline)` — distance to the
  design orbit (element centres + prelude / postlude extrapolated along the
  tangent).
- The design polyline is rebuilt at the start of every run.
- Only **growth** is penalized — betatron damping carries zero penalty.

### Optimization JSON config (minimal example)

```json
{
  "experiment_id": 4,
  "n_iterations": 10000,
  "lattice_max_gid": 207,
  "kick_every": 225,
  "kick_scale": 0.25,
  "loss": {
    "magnet_sp_numbers": [1, 2, 3],
    "exclude_device_types": ["DP_L_PM", "DP_sh_PM", "DQ_I1_1A", "..."],
    "extra_target_magnets": [
      {"x": 127.0, "y": 42.0, "z": 0.0, "weight": 20.0}
    ],
    "direction_loss": {
      "reference": [0.891, 0.454, 0.0],
      "weight": 1000.0
    },
    "envelope_loss": { "weight": 1e6 }
  },
  "parameter_groups_per_position": [
    {
      "sp_types": ["standard"],
      "exclude_device_types": ["DP_L_PM", "DP_sh_PM"],
      "value_column": "arg_val",
      "step_abs": 0.1
    }
  ]
}
```

### Tuning knobs

| Parameter | Meaning |
|---|---|
| `kick_every` | every N iterations, all groups are perturbed at once |
| `kick_scale` | kick amplitude in A (sign drawn at random) |
| `step_abs` | normal-iteration step for `arg_val` (one random group per iter) |
| `step_percent` | step for `correct_coeff` as a fraction of current value |
| `lattice_max_gid` | restrict tracking to first N gids of the layout |

---

## Build and run

### Prerequisites

| Platform | Toolchain |
|---|---|
| Windows | MSVC 2019+ (Build Tools), CMake 3.10+, Python 3.9+ |
| Linux   | gcc 9+ or clang 10+, CMake 3.10+, OpenMP, Python 3.9+ |
| macOS   | clang with `libomp` (e.g. `brew install libomp`), CMake, Python |

### C++ build (both platforms)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release -j
```

Produces:
- `build/SynchrotronTracker[.exe]`
- `build/SynchrotronOptimizer[.exe]`

### Python environment

```bash
python -m venv venv
# Windows:
venv\Scripts\activate
# Linux / macOS:
source venv/bin/activate

pip install pandas numpy matplotlib pyvista
```

### Full pipeline

```bash
# 1. (if LTE changed) rebuild the DB tables
python lte_to_lattice.py
python assemble_lattice_configs.py

# 2. (once) convert Ansys .fld -> .bin
build/SynchrotronTracker --convert         # Windows: build\SynchrotronTracker.exe --convert

# 3. run experiment 4
build/SynchrotronTracker 4
python visualize.py 4

# 4. prepare an optim config and run optimization
python prepare_optim_config.py workspace/optim_inputs/sp3_drift_2sp.json
build/SynchrotronOptimizer workspace/optim_inputs/sp3_drift_2sp.optimconfig

# 5. after optimization — rebuild lattice & re-render
python assemble_lattice_configs.py
build/SynchrotronTracker 4
python visualize.py 4

# 6. inspect convergence and stability
python viz_optim_history.py workspace/results/sp3_drift_2sp.history.json
python analyze_trajectory_stability.py 4
```

Python scripts find the project root automatically (`_common.project_root`)
and the C++ binaries find `workspace/` automatically (`findWorkspace()` in
`src/main.cpp`), so the commands work from any working directory.

---

## Diagnostic scripts

### `analyze_trajectory_stability.py [exp_id]`
Stability analysis of all 5 particles relative to the design orbit:
- per-particle envelope `|particle − orbit|(t)`
- linear and exponential fits of amplitude growth, qualified by R²
- **FFT** to extract betatron tune frequencies
- per-particle verdict: STABLE / OSCILLATING / GROWING / DRIFTED
- special case for P1 (`A_init ≈ 0` → any drift is DRIFTED)

### `find_orbit_deviation_peaks.py [exp_id]`
Finds local maxima of `|P1 − design_polyline|` and identifies the nearest
lattice element (gid, type, sp_number). Useful for spotting miscalibrated
magnets.

### `scan_field_maps.py [name] [--fix]`
Search for "dropouts" in field maps — grid nodes with `|B| ≈ 0` surrounded
by significant field. With `--fix`, patches by averaging 6 neighbours and
writes a `.bak` backup.

### `viz_optim_history.py history.json`
Plots loss, `mean_mm`, `max_mm` over iterations and parameter evolution.

### `live_visualizer.py [exp_id]`
PyVista window with auto-reload: polls `sp_elements.csv` and
`hardware_library.csv`; on change, re-runs `assemble_lattice_configs.py`
+ `SynchrotronTracker` and refreshes the view (camera preserved).

---

## Cross-platform notes

The project runs on Windows, Linux and macOS. A few platform-specific
details:

- **Executable suffix.** `_common.exe_name(stem)` returns `stem + ".exe"` on
  Windows and `stem` elsewhere. All Python scripts use this — never hardcode
  `.exe`.
- **Workspace discovery.** Both C++ binaries call `findWorkspace(argv[0])`,
  which walks up from the binary, `argv[0]` parent and the cwd looking for a
  directory containing `workspace/`. Python scripts use
  `_common.project_root()` for the same purpose (walks up to a directory
  with `CMakeLists.txt` or `README.md`).
- **Matplotlib backend.** `_common.set_matplotlib_backend("TkAgg")` tries
  Tk first and silently falls back to Qt5Agg or Agg, so headless Linux
  builds without Tk still work.
- **Path separators.** Inside scripts everything uses `pathlib.Path` or
  forward slashes — both are accepted by Windows.
- **OpenMP on macOS.** `libomp` from Homebrew is needed (Apple clang ships
  without OpenMP). CMake should pick it up automatically once installed.

---

## Performance

Tracking time, 5 particles, exp 4 (SP1–3 + drift, ~290k steps):

| Version | Tracking | Speedup |
|---|---|---|
| Original RK4 + linear lattice scan | 3424 ms | 1× |
| Boris pusher | 740 ms | 4× |
| + Spatial bucket index | 70 ms | **47×** |
| Optimizer with incremental loss | ~110 ms/iter | **~30×** on optimizer |

Key optimizations still in the code:

- 2D bucket grid — per step ~5 magnets checked instead of 2762
- Bounding-sphere prescreening — early exit for out-of-range elements
- Boris pusher — 1 field eval / step (vs 4 in RK4)
- Incremental loss — `minDist²[i][m]` updated in the tracking loop, no
  trajectory storage
- Composite loss with weights — TARGETs pin final convergence

---

## Troubleshooting

### `lattice 104 has no rows in lattice_configs.csv`
Run `python assemble_lattice_configs.py`. The optimizer can rewrite only
its own lattice; sub-lattices are auto-regenerated by the assembler.

### `ERROR: baseline tracking failed`
At least one of the 5 particles did not reach the finish plane within the
50 M step cap. Likely causes:

- current `arg_val` values give an unstable orbit on the long path
  (e.g. lattice tuned for SP1-10 diverges on SP1-21)
- `max_turns` (`skipSteps`) too small — the plane is crossed multiple
  times along the way, more skips needed
- run `build/SynchrotronTracker N` directly to see where particles get lost

### `Access is denied` when launching the .exe (Windows)
Windows Defender may quarantine a freshly built `.exe`. Add the `build/`
folder to Defender exclusions.

### Particles start inside the first SP lens
Check `particle_groups.csv` — initial position should be ~2.94 m before the
first magnet along the orbit tangent so there is a prelude drift. For
`group_id = 1` this is `(0, 0, 0)` moving +X (SP1.first is at `x ≈ 2.94`).
For `group_id = 3` (start at SP18) — `(80.16, 333.89, 0)` moving along
`(cos 153°, sin 153°, 0)`.

### Optimization "gets stuck"
Tune the `step_abs` vs `kick_scale` balance:

- `kick_scale > step_abs` — aggressive kicks escape local minima
- `kick_scale ≈ step_abs` — conservative kicks, precise convergence
- smaller `kick_every` — kicks more often (but more progress rollback)

For final polishing, lower `step_abs` from 1.0 to 0.1 A — slower but more
accurate.

### Analyzer shows "everything growing"
Run a **fresh tracker** after optimization:

```bash
python assemble_lattice_configs.py
build/SynchrotronTracker N
python analyze_trajectory_stability.py N
```

Otherwise the analyzer is looking at a stale `trajectories_expN.bin`.

### `find_orbit_deviation_peaks.py` shows "humps" near a specific magnet
That's a **miscalibration signature** — the same magnet type across SPs
gives systematic deviations. Tweak its `correct_coeff` in
`hardware_library.csv` or include it in the optimizer groups.

---

## License and contact

Developed for SILA Synchrotron Light Source research.
GitHub: <https://github.com/SashaBray/SynchrotronSim>
