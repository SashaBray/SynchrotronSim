# ⚛️ Synchrotron Particle Tracker

Relativistic charged-particle tracking through a synchrotron magnet lattice.
Uses **real field maps** (computed in Ansys), a **symplectic Boris pusher** for
integration, and a **stochastic optimizer** for tuning magnet currents.

---

## 📖 Contents
1. [🚀 Quick Start](#quick-start)
2. [Overview](#overview)
3. [Repository Layout](#repository-layout)
4. [Physical Model](#physical-model)
5. [Database Tables](#database-tables)
6. [Configuration Files](#configuration-files)
7. [Pipeline and Experiments](#pipeline-and-experiments)
8. [Optimizer Loss Function](#optimizer-loss-function)
9. [Build and Run](#build-and-run)
10. [Diagnostic Scripts](#diagnostic-scripts)
11. [Performance](#performance)
12. [Troubleshooting](#troubleshooting)

---

## 🚀 Quick Start <a name="quick-start"></a>

**Everything is already built and configured? Start here.**

```cmd
:: 1. Run exp 4 (focusing at SP3 exit) -- quick sanity test
build\SynchrotronTracker.exe 4
python visualize.py 4

:: 2. Look at trajectory stability
python analyze_trajectory_stability.py 4

:: 3. Run the optimizer (~80 min for 10000 iters on exp 4)
build\SynchrotronOptimizer.exe workspace\optim_inputs\sp3_drift_2sp.optimconfig

:: 4. Inspect convergence
python viz_optim_history.py workspace\results\sp3_drift_2sp.history.json
python assemble_lattice_configs.py
build\SynchrotronTracker.exe 4
python visualize.py 4
```

**If you modify anything in the DB** (sp_elements, particle_groups, LTE) — re-run:
```cmd
python assemble_lattice_configs.py    :: rebuild lattice_configs.csv
```
(if LTE changed — first `python lte_to_lattice.py`)

**If you modify optim_inputs/*.json** — regenerate:
```cmd
python prepare_optim_config.py workspace\optim_inputs\sp3_drift_2sp.json
```

---

## Overview

**What is simulated:** 3D trajectories of electrons (or other charged particles)
in a synchrotron magnet ring. Each magnet is represented by a **field on a
Cartesian grid** (computed in Ansys, converted to `.bin`), to which a geometric
transform (position/orientation) and calibration multipliers
(`correct_coeff`, `field_direction`, `currentScale`) are applied.

**What's available:**
- **`SynchrotronTracker.exe`** — forward tracking of 5 particles through a
  given lattice
- **`SynchrotronOptimizer.exe`** — stochastic hill-climber with kicks; finds
  magnet currents and calibrations that minimize a composite loss
- Python tooling: LTE-driven lattice generation, visualization (PyVista),
  stability analysis (FFT, envelope), optimizer history

---

## Repository Layout

```text
SynchrotronSim/
├── src/                                # C++ source code
│   ├── main.cpp                        # SynchrotronTracker
│   ├── optimizer.cpp                   # SynchrotronOptimizer
│   ├── Integrator.h                    # Boris pusher (sympl. relativistic)
│   ├── LatticeManager.h                # magnets, transforms, spatial bucket index
│   ├── FieldMap.h                      # trilinear interp on binary maps
│   ├── FinishPlane.h                   # sign-change finish detection
│   └── AnsysConverter.h                # .fld -> .bin converter
│
├── workspace/                          # data root
│   ├── configs/                        # ★ JSON configs (external, hand-edited)
│   │   ├── lte_mapping.json            # LTE types -> hardware_library mapping
│   │   ├── sub_lattices.json           # sub-lattice definitions (102/103/104)
│   │   └── visualizer.json             # camera, colors, opacity
│   ├── lte/                            # Elegant LTE sources
│   │   ├── SYLA_desc.lte               # device definitions
│   │   └── SYLA_line.lte               # element sequence
│   ├── tables/                         # DB (`;`-separated CSV)
│   │   ├── hardware_library.csv        # device_type_id -> field map + calibrations
│   │   ├── sp_elements.csv             # SP templates (standard/junction) (generated)
│   │   ├── lattice_layout.csv          # ring unrolling, geometry (generated)
│   │   ├── lattice_configs.csv         # simulator input (generated)
│   │   ├── experiments.csv             # exp_id;lat;group;dt;max_turns;traj_stride;p1..p3
│   │   └── particle_groups.csv         # 5 particles per group (center + 4 satellites)
│   ├── optim_inputs/                   # optimizer setups
│   │   ├── sp1_endpoint.json           # exp 2 (legacy: end of SP1)
│   │   ├── sp3_drift_2sp.json          # exp 4: SP3 + 2 SP free drift
│   │   ├── sp10_drift_to_sp11.json     # exp 5: SP1..10 + drift to SP11
│   │   └── junction_tune.json          # exp 6: junction tuning (SP20)
│   ├── fields_raw/                     # raw .fld from Ansys
│   ├── fields_bin/                     # binary field maps (.bin, generated)
│   ├── geometry/                       # STL device models
│   └── results/                        # output: trajectories_exp*.bin, history.json
│
├── Python scripts (top-level):
│   ├── lte_to_lattice.py               # LTE -> sp_elements + lattice_layout
│   ├── assemble_lattice_configs.py     # sp_elements + layout -> lattice_configs
│   │                                   # + sub-lattices (102/103/104)
│   ├── prepare_optim_config.py         # JSON -> .optimconfig (flat format)
│   ├── apply_optim_history.py          # history.json -> CSV (state rollback)
│   ├── visualize.py                    # 3D viewer (PyVista)
│   ├── live_visualizer.py              # live-reload viewer
│   ├── viz_optim_history.py            # plot loss curves, params
│   ├── analyze_trajectory_stability.py # envelope analysis, FFT betatron freqs
│   ├── find_orbit_deviation_peaks.py   # local max(|particle-orbit|)
│   └── scan_field_maps.py              # field map sanity check
│
├── CMakeLists.txt                      # cmake config
├── README.md                           # this file (markdown)
├── README.txt                          # plain-text variant
└── .gitignore
```

---

## Physical Model

### Integrator: relativistic Boris pusher

Equations of motion:
```
dx/dt = v
dp/dt = q·(E + v × B)
```
where `p = γ·m·v` is the relativistic momentum, `γ = √(1 + p²/(m²c²))`.

**Boris algorithm (drift-kick-drift splitting):**
1. half-drift: `x ← x + (p/(γ·m))·dt/2`
2. field at `x_{n+1/2}`: `B, E ← getTotalFields(x)`
3. half E-kick: `p⁻ = p + q·E·dt/2`
4. magnetic rotation: `t = q·B·dt/(2·γ⁻·m)`, `p⁺ = p⁻ + (p⁻ + p⁻×t) × 2t/(1+|t|²)`
5. half E-kick: `p_{n+1} = p⁺ + q·E·dt/2`
6. half-drift: `x ← x + (p/(γ·m))·dt/2`

**Properties:**
- **Symplectic** (preserves phase-space volume)
- **Conserves |p|² exactly** for pure magnetic fields
- **1 field eval / step** (vs 4 in RK4)

### Field at a point

Magnet field is computed via **trilinear interpolation** in the local coordinate system:
```
1. world → local (transform.transformPointInv)
2. trilinear(map.data, local_pos)
3. local → world (transform.rotateVector)
4. ×= correct_coeff × field_direction × currentScale
```

Sum of fields from all magnets near the position is computed via a
**spatial bucket grid** (2D grid in xy, bucket size = 2 × max(world_radius)).
Magnets out of range are skipped entirely.

---

## Database Tables

### 1. hardware_library.csv

| Column | Description |
|---|---|
| `device_type_id` | Unique magnet type ID (QP_L_1A, DQ_I1_1A, …) |
| `Active` | True/False — active electromagnet or permanent magnet |
| `field_type` | MAGNETIC/ELECTRIC |
| `map_raw_path` | path to Ansys `.fld` file |
| `map_bin_path` | path to compiled `.bin` |
| `stl_path` | path to STL model (for visualization) |
| `use_symmetry` | True if field map is symmetric across the YZ plane |
| `nominal_arg` | nominal current at which the field map was built |
| `field_direction` | +1 / -1 to flip direction |
| `correct_coeff` | calibration field multiplier |

### 2. sp_elements.csv

Superperiod templates. For each `(sp_type, sp_element_id)` pair defines the device
type and its value.

| Column | Description |
|---|---|
| `sp_type` | "standard" or "junction" |
| `sp_element_id` | 1..N (order within SP) |
| `device_type_id` | reference to hardware_library |
| `arg_val` | current [A], or "-" for permanent magnets |
| `correct_coeff` | per-instance field multiplier |

### 3. lattice_layout.csv (generated)

Full ring unrolling. Per-element position and orientation.

| Column | Description |
|---|---|
| `global_id` | 1..2762 |
| `sp_number` | 1..39 (superperiod number) |
| `sp_type` | standard / junction |
| `sp_element_id` | reference to sp_elements |
| `x, y, z` | element position (m) |
| `yaw, pitch, roll` | orientation (rad) |

### 4. lattice_configs.csv (generated)

Final simulator input table. Joins layout with sp_elements (substituting
arg_val, correct_coeff) and adds sub-lattice copies from `sub_lattices.json`.

### 5. experiments.csv

| Column | Description |
|---|---|
| `exp_id` | unique experiment ID |
| `lattice_id` | which lattice to use (101 = full, 102+ = sub) |
| `group_id` | which particle group (from particle_groups.csv) |
| `dt` | integration timestep [s] |
| `max_turns` | finish-plane skipSteps (how many steps to skip before checking) |
| `traj_stride` | save 1 point per N (1 = all) |
| `p1..p3 (x,y,z)` | three finish-plane points (normal via cross product) |

### 6. particle_groups.csv

| Column | Description |
|---|---|
| `group_id` | group ID (1, 2, 3, …) |
| `part_id` | 2..6 (5 particles per group) |
| `m, q` | mass [kg], charge [C] |
| `x, y, z` | initial position [m] |
| `vx, vy, vz` | initial velocity [m/s] |

**Group convention:** typically 5 particles — 1 reference (on-axis) + 4 satellites
(±1 mm in local y and z) for estimating betatron oscillations.

---

## Configuration Files

Folder `workspace/configs/`:

### lte_mapping.json
- `lattice_id` of the main ring (101)
- LTE-type → hardware_library mapping
- Regex for DL/JL families
- Classification thresholds (lengths, K1)

### sub_lattices.json
- List of shortened lattice copies for experiments
- Each entry: `lattice_id`, `max_gid`, `description`
- Automatically rebuilt by `assemble_lattice_configs.py`

### visualizer.json
- Camera (position, focus, up)
- `parallel_projection`: true/false
- Colors and labels for 5 particles
- Magnet opacity (EM/PM)
- Live-visualizer parameters (poll interval, hotkeys)

---

## Pipeline and Experiments

Standard path from LTE to simulation:
```text
SYLA_desc.lte + SYLA_line.lte
    │   (lte_to_lattice.py)
    │   └─ DL/JL straightening (5 slices on chord)
    ▼
sp_elements.csv + lattice_layout.csv
    │   (assemble_lattice_configs.py)
    ▼
lattice_configs.csv  (includes sub-lattices 102/103/104)
    │
    ├─→ SynchrotronTracker.exe N  → trajectories_expN.bin
    │       │   (visualize.py N)
    │       ▼
    │   PyVista 3D window with trajectories
    │
    └─→ SynchrotronOptimizer.exe optim.optimconfig
            │   ┌─ reads sp_elements.csv, mutates, writes back
            │   └─ rebuilds lattice 104 in lattice_configs.csv
            ▼
        sp_elements.csv (new best-values)
        history.json (for viz_optim_history.py and apply_optim_history.py)
```

### Experiment list (experiments.csv)

| exp | Lattice | Purpose |
|---|---|---|
| 1 | 101 | Sanity test: plane at x=5, particles cross SP1 |
| 2 | 101 | Half ring — past SP20 (yaw=180°, skipSteps=100k) |
| 3 | 101 | **Full ring** — past SP1 after 360° (skipSteps=200k) |
| 4 | 102 (SP1-3) | Focusing at SP3 exit + 2 SP of free drift |
| 5 | 103 (SP1-10) | Long path SP1-SP10 + drift to SP11 |
| 6 | 104 (SP1-21) | Junction tuning (SP20) with particles starting from SP18 |

---

## Optimizer Loss Function

Composite loss with **three components**:

### 1. Position loss (geometry)
```
Σ_(magnet, particle)  w_m · d²_m,i
```
- `d_m,i` — min distance from particle `i`'s trajectory to magnet `m`'s center
- Computed **incrementally** in the tracking loop (no trajectory storage)
- Magnets filtered by `magnet_sp_numbers` and `exclude_device_types`
- Additional `TARGET` points (virtual magnets) with own `weight`

### 2. Direction loss (parallelism)
```
w_dir · Σ_particle  sin²(θ_i)
```
- `θ_i` — angle between final particle momentum and reference direction
- `sin²θ` ≈ θ² for small angles (quadratic scale)
- Goal: beam exits **parallel** to a given direction (orbit tangent)

### 3. Envelope loss (betatron amplitude growth)
```
w_env · Σ_particle  max(0, A_final - A_init)²
```
- `A(t) = dist(particle_pos, design_polyline)` — distance from the orbit
- design polyline is rebuilt **each run** from lattice element positions +
  prelude/postlude tangent extrapolation
- Only **growth** is penalized (damping is a useful effect)
- Includes the reference particle P1 (drift away from orbit is penalized too)

### Optimization JSON config structure

```json
{
  "experiment_id": 4,
  "n_iterations": 10000,
  "lattice_max_gid": 207,
  "kick_every": 225,
  "kick_scale": 0.25,
  "loss": {
    "magnet_sp_numbers": [1, 2, 3],
    "exclude_device_types": ["DP_L_PM", "DP_sh_PM", "DQ_I1_1A", ...],
    "extra_target_magnets": [
      {"x": ..., "y": ..., "z": ..., "weight": 20.0}
    ],
    "direction_loss": {
      "reference": [0.891, 0.454, 0.0],
      "weight": 1000.0
    },
    "envelope_loss": {
      "weight": 1e6
    }
  },
  "parameter_groups_per_position": [
    {
      "sp_types": ["standard"],
      "exclude_device_types": [...],
      "value_column": "arg_val",
      "step_abs": 0.1
    }
  ]
}
```

### Tuning parameters

| Parameter | Meaning |
|---|---|
| `kick_every` | Every N iterations all groups are perturbed at once |
| `kick_scale` | Kick amplitude in A (× sign of a random number) |
| `step_abs` | Normal-iteration step (one group per iter) for arg_val |
| `step_percent` | Step for correct_coeff (% of current value) |
| `lattice_max_gid` | Tracking limited to first N gids of the layout |

---

## Build and Run

### C++ build

```bash
cmake -B build
cmake --build build --config Release
```
Two binaries are produced:
- `build/SynchrotronTracker.exe` — tracker (reads lattice_configs)
- `build/SynchrotronOptimizer.exe` — optimizer (reads .optimconfig)

### Python environment

```bash
python -m venv venv
venv\Scripts\activate          # Windows
pip install pandas numpy matplotlib pyvista
```

### Full pipeline

```cmd
:: 1. Build DB tables from LTE (if LTE changed)
python lte_to_lattice.py
python assemble_lattice_configs.py

:: 2. (once) convert field maps from Ansys
build\SynchrotronTracker.exe --convert

:: 3. Run experiment 4
build\SynchrotronTracker.exe 4
python visualize.py 4

:: 4. Prepare optim config and run optimization
python prepare_optim_config.py workspace\optim_inputs\sp3_drift_2sp.json
build\SynchrotronOptimizer.exe workspace\optim_inputs\sp3_drift_2sp.optimconfig

:: 5. After optimization: rebuild lattice_configs, re-render trajectories
python assemble_lattice_configs.py
build\SynchrotronTracker.exe 4
python visualize.py 4

:: 6. Inspect convergence and stability
python viz_optim_history.py workspace\results\sp3_drift_2sp.history.json
python analyze_trajectory_stability.py 4
```

---

## Diagnostic Scripts

### `analyze_trajectory_stability.py [exp_id]`
Stability analysis of 5 particles relative to the design orbit:
- Envelope `|particle - orbit|` per particle
- Linear and exponential fits of amplitude growth
- **FFT** to find betatron frequencies (tune)
- Verdict: STABLE / OSCILLATING / GROWING / DRIFTED
- R² qualification (low R² ⇒ gamma is unreliable)
- Special case for P1 (a_init≈0 ⇒ any drift is DRIFTED)

### `find_orbit_deviation_peaks.py [exp_id]`
Finds local maxima of `|P1 - design_polyline|`. For each peak, identifies
the **nearest lattice element** (gid, type, sp_number). Useful for spotting
miscalibrated magnets.

### `scan_field_maps.py [name] [--fix]`
Search for "dropouts" in field maps — grid nodes with |B|≈0 surrounded by
significant field. With `--fix`, patches by averaging 6 neighbors and writes a
`.bak`.

### `viz_optim_history.py history.json`
Plots loss, mean_mm, max_mm over iterations + parameter evolution.

---

## Performance

After optimizations:

| Version | Tracking time (5 particles × SP3) | Speedup |
|---|---|---|
| Original RK4 + linear lattice scan | 3424 ms | 1× |
| Boris pusher | 740 ms | 4× |
| + Spatial bucket index | 70 ms | **47×** |
| + Incremental loss | ~110 ms/iter (incl. loss) | **~30×** on optimizer |

**Key optimizations:**
- **Spatial 2D bucket grid** — per step ~5 magnets checked instead of 2762
- **Bounding-sphere prescreening** — early exit for out-of-range elements
- **Boris pusher** — 1 field eval/step vs 4 in RK4
- **Incremental loss** — `minDist²[i][m]` is updated in the tracking loop, no
  trajectory storage
- **Composite loss with weights** — TARGET points with large weight pin
  final convergence

---

## Troubleshooting

### `lattice 104 has no rows in lattice_configs.csv`
Run `python assemble_lattice_configs.py` — the optimizer rewrites only its own
lattice, sub-lattices may be wiped. The script auto-regenerates them.

### `ERROR: baseline tracking failed`
At least one of 5 particles didn't reach the finish plane within the 50M-step
cap. Possible causes:
- Current `arg_val` in `sp_elements.csv` give an unstable orbit on long paths
  (e.g., lattice tuned for SP1-10 but tracking through SP1-21 diverges)
- `max_turns` (skipSteps) too small — the plane is crossed multiple times along
  the way, you need to skip extra crossings
- Run `build\SynchrotronTracker.exe N` to inspect actual trajectories and see
  where particles get lost

### `Access is denied` when launching .exe
Windows Defender may "eat" a freshly built .exe. Add the `build/` folder to
Defender exclusions.

### Particles start right inside the first SP lens
Check `particle_groups.csv` — position should be **~2.94 m before** the first
magnet along the orbit tangent (so there's a prelude drift). For group_id=1
this is `(0, 0, 0)` moving +X (SP1.first at x=2.94). For group_id=3 (start at
SP18) — `(80.16, 333.89, 0)` moving along (cos153°, sin153°, 0).

### Optimization "gets stuck"
Tune the `step_abs` vs `kick_scale` balance:
- **kick_scale > step_abs** → aggressive kicks can jump out of local minima
- **kick_scale ≈ step_abs** → conservative kicks, precise convergence
- **kick_every** smaller → kicks more often (but more progress rollback)

For final polishing, lower `step_abs` (e.g., from 1.0 to 0.1 A) — slower but
more accurate.

### Analyzer shows "everything growing"
Run a **fresh tracker** after optimization:
```cmd
python assemble_lattice_configs.py
build\SynchrotronTracker.exe N
python analyze_trajectory_stability.py N
```
Otherwise the analyzer is looking at stale `trajectories_expN.bin`.

### `find_orbit_deviation_peaks.py` shows "humps" around a specific magnet
This is a **miscalibration signature**: the same magnet across different SPs
gives systematic deviations. You may need to tweak its `correct_coeff` in
`hardware_library.csv` or include it in the optimizer groups.

---

## License and Contact

This project is being developed for SILA Synchrotron Light Source research.
GitHub: https://github.com/SashaBray/SynchrotronSim
