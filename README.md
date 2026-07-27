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
8. [Optimizer methods](#optimizer-methods)
9. [Optimizer loss function](#optimizer-loss-function)
10. [Matrix model (alternative engine)](#matrix-model-alternative-engine)
11. [Build and run](#build-and-run)
12. [Diagnostic scripts](#diagnostic-scripts)
13. [Cross-platform notes](#cross-platform-notes)
14. [Performance](#performance)
15. [Troubleshooting](#troubleshooting)

---

## Quick start

Assuming the project is already built and field maps are in place. The DB
currently contains lattice 101 (full ring) and optional sub-lattices for
short experiments. Optimization is launched via the `opt` wrapper which
resolves an experiment by its `experiment_id`:

```bash
# 1. Run experiment 1 (lattice 101, full ring) and view 3D trajectories
build/SynchrotronTracker 1         # Windows: build\SynchrotronTracker.exe 1
python visualize.py 1

# 2. Check trajectory stability (envelope + FFT)
python analyze_trajectory_stability.py 1

# 3. Optimize currents — choose mode
opt 4 --method=global              # Global mode (whole config); search = line_search | monte_carlo
opt 4 --method=windowed            # Windowed mode (sliding window over magnets)
opt 4 --method=gradient            # Gradient descent (FD + Adam/normalized GD + line search)

# Cross-mode orchestration (alternates windowed and global across cycles)
opt_combo 4

# Per-element fast overrides without editing JSON
opt 4 --method=pe --W=5 --iters=200 --passes=3

# 4. Inspect convergence and re-render with new state
python viz_optim_history.py workspace/results/<config>.history.json
python viz_mc_evolution.py        # parameter evolution (latest run autodetected)
python assemble_lattice_configs.py
build/SynchrotronTracker 1
python visualize.py 1
```

After editing **any** DB table by hand (sp_elements, particle_groups, LTE files)
re-run `python assemble_lattice_configs.py` so the simulator picks up the
changes. If LTE has changed, run `python lte_to_lattice.py` first, then
`python lte_to_design_orbit.py` to regenerate the design-orbit segments DB
(used by the per-element and envelope losses).

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
| `continuation_tuner.py` | **multi-superperiod driver** (track ring → find loss point → polish window → repeat); see [Continuation tuner](#continuation-tuner-multi-superperiod-periodicity) |
| `manual_kick.py` | one-shot random kick of all tunable currents/coeffs (±amp%), with ОДЗ clamps |
| `apply_optim_history.py` | history.json → CSV (state rollback / checkpoint restore) |
| `visualize.py` | static 3D viewer (PyVista); overlays kick model if `trajectories_expN_kick.bin` present |
| `live_visualizer.py` | live-reload viewer (auto-tracker on DB change) |
| `element_compare.py` | compare **diffur** (Boris/maps) vs **kick** (`matrix_tracker`, `SS_SECTOR`) per element → Word report + datasets + 3D viz configs; commands `all` / `combo` / `viz` / `viz-all` |
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
│   │   ├── sub_lattices.json           # sub-lattice definitions (extra IDs derived from 101)
│   │   └── visualizer.json             # camera, colors, opacity for PyVista
│   ├── lte/                            # Elegant LTE sources
│   │   ├── SYLA_desc.lte               # device descriptions
│   │   └── SYLA_line.lte               # element sequence
│   ├── tables/                         # database (`;`-separated CSV)
│   │   ├── hardware_library.csv        # device_type_id -> field map + calibration
│   │   ├── sp_elements.csv             # SP templates (generated, then edited)
│   │   ├── lattice_layout.csv          # full ring unrolling (generated)
│   │   ├── lattice_configs.csv         # final simulator input (generated)
│   │   ├── design_orbit_segments.csv   # per-sp_type list of line/arc segments
│   │   │                               # (generated by lte_to_design_orbit.py)
│   │   ├── experiments.csv             # exp_id; lat; group; dt; max_turns; stride; p1..p3
│   │   ├── particle_groups.csv         # 5 particles per group (centre + 4 satellites)
│   │   └── backups/                    # manual snapshots of sp_elements.csv
│   ├── optim_inputs/                   # one JSON per optimization experiment
│   ├── fields_raw/                     # raw .fld files from Ansys (large, gitignored)
│   ├── fields_bin/                     # compiled .bin field maps (generated)
│   ├── geometry/                       # STL models (large, gitignored)
│   └── results/                        # outputs:
│       ├── trajectories_exp*.bin       # tracker output per experiment
│       ├── design_orbit_lat{N}.csv     # assembled reference orbit (per lattice)
│       ├── magnet_endpoints_lat{N}.csv # per-magnet (entry,exit) on the reference orbit
│       ├── *.history.json              # optimizer iteration log
│       ├── *.history.png               # convergence plot
│       └── combined_history_exp{N}/    # per-cycle windowed+global histories + merged plot
│
├── docs/                               # documentation
│   ├── algorithms.md                   # algorithm reference (RD 50-34.698-90)
│   ├── algorithms.docx                 # Russian GOST-formatted version
│   └── build_docx.js                   # docx generator
│
├── opt.bat / opt_combo.bat             # Windows CLI wrappers for optimize.py / optimize_combined.py
├── optimize.py                         # exp_id → JSON → flat config → optimizer.exe
├── optimize_combined.py                # alternates windowed and global across cycles, builds combined plot
├── lte_to_lattice.py                   # .lte → lattice_layout.csv + sp_elements.csv
├── lte_to_design_orbit.py              # .lte → design_orbit_segments.csv
├── design_orbit.py                     # in-memory orbit assembly + magnet endpoints
├── prepare_optim_config.py             # JSON → flat .optimconfig + design_orbit/endpoints CSV
├── assemble_lattice_configs.py         # build lattice_configs.csv from layout + sp_elements
├── visualize.py                        # PyVista 3D viewer (lattice + trajectories + design orbit + kick overlay)
├── element_compare.py                  # diffur vs kick per element (report/datasets/3D viz configs)
├── analyze_trajectory_stability.py     # FFT + envelope + amplitude growth diagnostics
├── viz_optim_history.py                # convergence curves from history.json
├── viz_mc_evolution.py                 # per-parameter evolution heatmap
├── CMakeLists.txt                      # cmake config (MSVC / gcc / clang)
├── _common.py                          # shared Python helpers (cross-platform)
├── README.md                           # this file
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

### `hardware_library.csv` (10 columns, hand-edited)

| Column | Description |
|---|---|
| `device_type_id` | unique magnet-type ID (`QP_L_1A`, `DQ_I1_1A`, …) |
| `Active` | `True` for an active electromagnet, `False` for a permanent magnet |
| `field_type` | `B` (magnetic, all current devices) / `E` (electric — placeholder for future electrostatic elements) |
| `map_raw_path` | filename of the Ansys `.fld` source (under `workspace/fields_raw/`) |
| `map_bin_path` | filename of the compiled `.bin` (under `workspace/fields_bin/`) |
| `stl_path` | STL model for the visualizer (`-` if none) |
| `use_symmetry` | `1` ↔ field map is antisymmetric in x across the YZ plane (storage halved); `0` ↔ no mirror |
| `nominal_arg` | nominal current the field map was computed at (`-` for PM) |
| `field_direction` | sign multiplier `±1` |
| `correct_coeff` | device-level field calibration |

`field_type` and `use_symmetry` are read at load time and stored on each
`MagneticElement`. Today all entries are `B` / `1`, but the routing in
`MagneticElement::addContribution` and the mirror branch in
`FieldMap::getInterpolatedField` handle the other cases — adding an `E`
row or an asymmetric magnet requires no code change.

### `sp_elements.csv` (templates, hand-edited)

Defines per-superperiod (SP) magnet content, indexed by `(sp_type,
sp_element_id)`.

| Column | Description |
|---|---|
| `sp_type` | `standard` or `junction` |
| `sp_element_id` | 1..N inside the SP |
| `device_type_id` | reference into `hardware_library.csv` |
| `arg_val` | current in Amperes; **empty or `-`** for permanent magnets (`Active=False`: no current — field comes entirely from `correct_coeff`). Both empty and `-` parse to 0 and the element stays in the lattice (a PM with `currentScale=1`). |
| `correct_coeff` | per-instance field multiplier (default 1.0); the **only** tunable handle for PM dipoles |

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
`{lattice_id, max_gid, description}`. The list is hand-edited per project;
delete entries to drop unused sub-lattices, add a new entry to create one.
After editing run `assemble_lattice_configs.py` to refresh
`lattice_configs.csv`.

### `visualizer.json`
Camera (position, focus, up), parallel-projection flag, particle colours and
labels (5 entries), magnet opacity (EM and PM), live-visualizer poll
interval, hotkey bindings. Also `design_polyline` block (color, line_width,
points_per_segment) for the reference orbit overlay.

### `workspace/optim_inputs/*.json`
One JSON per optimization scenario. Top-level knobs (besides the loss block
described above):
- `experiment_id`, `lattice_max_gid`, `n_iterations`, `random_seed`
- `kick_every`, `kick_scale`, `revert_after_kicks`
- `optimization_mode` (`global` | `windowed` | `gradient`; legacy `optimization_method`/`montecarlo`/`per_element` still accepted)
- `gradient`, `per_element`, `combined`, `dynamic_loss` blocks
- `parameter_groups_per_position[]` — what to mutate and how

The companion `.optimconfig` file is a flat key/value text consumed by the
C++ binary. It is regenerated from the JSON every time `opt` (or
`prepare_optim_config.py`) runs.

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

Current DB contents (after recent cleanup — full ring + one drift demo):

| exp | Lattice | Purpose |
|---|---|---|
| 1 | 101 | Sanity test: plane at `x=5`, particles cross SP1 |
| 2 | 101 | Half ring — past SP20 (yaw ≈ 180°, `skipSteps = 100k`) |
| 3 | 101 | **Full ring** — past SP1 after 360° (`skipSteps = 200k`) |
| 4 | 105 (SP1) | SP1 + 50 m empty drift; finish plane at the end of the drift. Tests how well SP1 alone collimates the beam over a long free-flight. |

Adding a new experiment:
1. Append a row to `workspace/tables/experiments.csv`. The columns p1..p3
   are three points on the finish plane (use a `±u_perp + z=±1` pattern).
2. If you reuse an existing `lattice_id`, no further DB change needed.
   For a new sub-lattice add an entry to `workspace/configs/sub_lattices.json`
   (`{lattice_id, max_gid, description}`) and re-run
   `python assemble_lattice_configs.py`.
3. Create an optimization config in `workspace/optim_inputs/<name>.json`
   referencing the new `experiment_id`. Launch with `opt <id> --method=...`.

---

## Optimizer methods

The optimizer supports **three modes** (selectable via
`--method=global|windowed|gradient`, synonym `--mode=`, or JSON field
`optimization_mode`), plus a Python orchestrator that alternates them.
Within the `global` and `windowed` modes the per-knob **search method** is
`line_search` (successive approximation, default), `monte_carlo` (random
kicks), or `svd` (response-matrix Gauss-Newton, global mode only) — set via
`global_search_mode` / `per_element.search_mode`. Legacy tokens
`mc`/`montecarlo`→`global`, `per_element`/`pe`→`windowed` are still accepted.

> **Note (experimental).** The `svd` method and the `line_search` step
> variants below were added in the 2026-06 continuation-tuner work. They are
> **not proven to beat the plain coarse-to-fine `line_search`** on this
> problem (per-parameter loss wells are sharp, which the zoom handles well);
> treat them as exploratory levers. See `research_log/2026-06-25_*`.

#### `line_search` step variants (`PER_ELEM_LS_*` / continuation flags)

The default `line_search` does a deterministic coarse-to-fine zoom (`K` probes
× `rounds`, halving the interval around the best). These flags change the
per-knob step rule:

| Flag (continuation_tuner) | JSON `per_element` | Effect |
|---|---|---|
| `--ls-first-improvement` | `line_search_first_improve` | accept the **first** value that beats the incumbent, skip the rest of the zoom (cheaper, coarser). Tends to rail parameters. |
| `--ls-spline` | `line_search_spline` | fit a **cubic spline** through K≥5 probes, take its minimum as a candidate (verified by sim). Overshoots on sharp V-wells (can predict negative loss). |
| `--ls-parabola` | `line_search_parabola` | fit a **weighted least-squares parabola** through K≥5 probes, take the vertex (guarded: must be a minimum, in-range, else keep best probe). Saner than spline; still loses to zoom on sharp wells. |
| `--ls-full-range` | `line_search_local_current/coeff = 0` | search/approximate across the **whole ОДЗ** (±100 A / ±3) instead of the ±10 A local neighborhood (probe count auto-raised to 11). |

The spline/parabola candidate is **always re-evaluated by simulation** and the
best *evaluated* point is kept, so a bad fit never corrupts the result — it
only wastes a sim.

### SVD / response-matrix (global mode, `--search svd`)

Classic accelerator orbit correction. Builds a finite-difference response
matrix `J = ∂r/∂param` (residual `r` = per-(particle,magnet) closest-approach,
so `‖r‖²` = position loss), then takes a damped **Gauss-Newton** step via a
**truncated SVD** of `JᵀJ` (drops near-degenerate corrector combinations), with
a trust region and backtracking. Optional **periodicity term**
(`--svd-angle-weight`, JSON `loss.svd_angle_weight`) adds the centroid's
exit-momentum transverse component (→0 ⇒ beam exits the superperiod parallel
to the orbit ⇒ closed periodic orbit) to both the residual and the composite
loss. Strong for a **big correction from a wrecked state** (e.g. max
2611→403 mm); near convergence the FD Jacobian is noise-limited (the
closest-approach metric has ~0.2 mm discretization jitter) and it stalls.

### Monte-Carlo search method — random kicks (global mode)

This is a *search method* (random current changes), not a mode; it runs inside
the **global** mode when `global_search_mode: monte_carlo` (alias `random`).
The default search method is `line_search` (successive approximation).

- One random parameter group is perturbed each iteration; new state accepted
  if loss decreases (greedy hill-climb).
- Every `kick_every` iterations: **all** groups perturbed at once
  (`±kick_scale`) to escape local minima.
- **Dynamic loss reweighting** (optional): the three composite-loss
  components (position / direction / envelope) rotate their weights by
  `w_i(t) = 1 + A·cos(2πt/T − 2πi/3)` so the optimizer's "focus" shifts
  between them and helps it escape attractors. See `dynamic_loss` JSON
  block; defaults to enabled with `amplitude=0.7` and auto-period
  `T = 5×N_params`.
- Strengths: simple, robust, escapes locals via kicks + dynamic reweighting.
- Weaknesses: slow convergence for many parameters.

### Gradient descent (FD + line search)

- Numerical gradient via **central** or **forward** finite differences.
- Update rule (selectable):
  - **`normalized_gd`** (default): `Δθ_j = −lr · g_j / max|g_k|` — step
    along strict −∇L direction, scaled so the dominant component moves by
    `±lr`. Strong on early descent.
  - **`adam`**: classic Adam with momentum + bias-correction. Adaptive per
    parameter, but on small iteration budgets the sign-based first step
    can hurt.
- **Backtracking line search**: tries `α = 1.0, 0.5, 0.25, ...` until the
  step improves loss; up to `max_halvings` tries before rejecting.
- Strengths: deterministic, monotonic descent, fewer iterations.
- Weaknesses: gets stuck in local minima (no kicks); FD-grad becomes noisy
  near a minimum.

### Windowed mode (sliding window)

The **windowed** mode (`--method=windowed`; JSON block still named
`per_element`) slides a window of `window_size` magnets along the lattice and
tunes only the knobs whose magnets fall in the window. Within each window the
search method is `line_search` (default) or `monte_carlo` (`per_element.search_mode`).

- **Window**: contiguous range of magnets by global_id; slides with `stride`.
- **Particle start**: ideal launch on the design orbit at window entry, or —
  with `start_from_sp_start: true` — particles are *propagated from the
  global start through the prefix*, giving the window a realistic
  (accumulated-error) input state. The prefix simulation is cached as a
  snapshot and reused across `iters_per_window` iterations of the same
  window (recomputed per iteration when `accumulate_groups=true`).
- **Per-window loss** (`per_element.loss_mode`):
  - `endpoint` (default): `‖P1−T‖² + 0.25·Σ‖Pᵢ−T‖² + dirW·(1−cos θ)² +
    pathW·Σd²(P1, polyline)`
  - `path_only`: `Σ d²(P1, design_polyline)` — pure deviation from the
    reference orbit, sampled every `path_sample_stride` pusher steps.
  - `composite`: same composite loss as the global mode, but `magnetPos` is
    filtered by gids inside the window. Extra TARGETs (gid = −1) always
    included.
  - `composite_full`: full simulation from global start to the window's
    finish plane, no lattice mask. Loss summed over **all** lattice
    magnets with gid ≤ gid_last — sees the whole trajectory leading up
    to the window. Most expensive.
- **`accumulate_groups`**: when true, the active-knob set *grows* across
  windows within a pass — each new window can tune its own knobs **plus**
  every knob touched by earlier windows. Resets between passes. Combined
  with `iters_per_param > 0`, `m_iter = iters_per_param × N_active_groups`
  scales budget with knob count.
- **Multi-SP groups** (per-position group covering several SPs) are
  registered against every gid they touch; tuning the knob in one window
  changes its value in *all* SPs (one-knob-per-type behavior).
- **CLI overrides** (without JSON edits):
  - `--W=N` — window size
  - `--stride=N` — window stride
  - `--iters=N` — iterations per window (ignored when `iters_per_param>0`)
  - `--passes=N` — number of full passes over the lattice

### Combined (windowed ⇄ global orchestrator)

`opt_combo <exp>` runs `optimize_combined.py`, which alternates short
windowed bursts with short global bursts:

```
for cycle in 1..cycles:
    opt.exe --method=windowed  (passes = pe_passes_per_cycle)
    opt.exe --method=global --loss=mc_loss_mode  (n_iter = mc_iters_per_cycle)
```

Between phases `sp_elements.csv` is updated on disk and the next phase
picks it up. After all cycles the per-phase histories are merged into a
single `combined_history.json` and rendered as `combined_history.png`
under `workspace/results/combined_history_exp{id}/`. Configured in the
`combined` JSON block.

### Continuation tuner (multi-superperiod periodicity)

`continuation_tuner.py` drives the **shared-current** tuning of the full 10-SP
ring (lattice 107). Because the currents of one superperiod are reused by all
ten, the task is **periodicity**: make one SP map the beam back onto itself so
the next SP is a copy. The driver creeps a "frontier" (рубеж) around the ring:

```
loop:
  1. track the full ring (exp 10) → trajectories_exp10.bin
  2. find_loss_gid: рубеж = last magnet before the beam first exceeds the
     threshold (per-magnet closest-approach metric — the SAME metric the
     optimizer minimizes), excluding DP_L_PM
  3. finish = round(fraction × рубеж), clamped to the --max-sp cap
  4. polish the window [start..finish] (global line_search) → write sp_elements
  5. re-track, accept the new рубеж as-is (dynamic; rollback allowed)
```

Key flags:

| Flag | Meaning |
|---|---|
| `--cycles N` `--budget B` | cycles, and per-polish simulation budget |
| `--fraction F` | `finish = round(F × рубеж)`; `F>1` puts the finish past the loss point |
| `--max-sp X` | cap the рубеж at X superperiods (e.g. `1.0` pins the window to SP1, the stable choice — see below) |
| `--ramp-stages N` | **progressive frontier**: ramp `finish` from `finish/N` to `finish` over N stages (window grows gradually; acts as a kick substitute) |
| `--search` | `line_search` (default) \| `monte_carlo` \| `svd` |
| `--svd-angle-weight W` | periodicity term (centroid exit ∥ orbit); works in `svd` (residual) and `line_search` (loss) |
| `--ls-first-improvement` / `--ls-spline` / `--ls-parabola` / `--ls-full-range` | line_search step variants (above) |
| `--parallelism-weight` / `--envelope-weight` | loss component weights |
| `--allow-rollback` | keep the final (possibly-worse) state instead of best |

```bash
# stable recipe: pin the window to SP1, push periodicity, original zoom
python continuation_tuner.py --cycles 4 --budget 1200 --fraction 1.0 \
    --max-sp 1.0 --search line_search --svd-angle-weight 3
```

> **Status / known wall.** Extending the window **across the SP1→SP2 boundary**
> (gid 69) reliably crashes: the shared currents get retuned for the diverged
> SP2 part of the window and break SP1, so the рубеж collapses. Pinning to SP1
> (`--max-sp 1.0`) is stable but creeps slowly (~+2 gid/cycle, plateau ≈ 1.6 SP
> with the honest 0.5 mm beam). The barrier is multi-SP periodicity of the
> shared currents, not local SP1 quality. See
> `research_log/2026-06-25_continuation_sp_boundary_wall.md`.

### When to use what

| Scenario | Best choice |
|---|---|
| Fresh start, lots of zeroed currents | **windowed** with `start_from_sp_start: true` — establishes a working baseline along the lattice |
| Polishing a well-tuned state | **global** + `monte_carlo` — kicks + dynamic loss escape locals |
| Geometry is fine but focusing is off | **global** (envelope is global, windowed cannot fix it locally) |
| Get all three components down simultaneously | **opt_combo** — windowed drives geometry, global drives envelope |
| Quick basin descent | **gradient** — clean monotonic step until FD noise dominates |

---

## Optimizer loss function

The loss has **three additive components**. A `--loss=composite|orbit|envelope`
flag (or JSON `loss.mode`) chooses which are active:

| Mode | Active components | When to use |
|---|---|---|
| `composite` (default) | position + direction + envelope | Final polishing |
| `orbit` (rough) | position only | Coarse initial setup (fast convergence to design orbit) |
| `envelope` | envelope only | Suppress betatron amplitude growth when orbit is already correct |

The three components:

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
  **true reference orbit** parsed from `.lte` (not the magnet-iron centres).
  The polyline is built from `workspace/tables/design_orbit_segments.csv`
  (a database of line/arc segments per `sp_type`, derived from
  `SYLA_desc.lte` + `SYLA_line.lte`) and saved per lattice into
  `workspace/results/design_orbit_lat{N}.csv`. Regenerate with
  `python lte_to_design_orbit.py` after `.lte` changes.
- Only **growth** is penalized — betatron damping carries zero penalty.

### Reference points for position loss

The `loss.position_reference` field (string) selects what magnet positions
are used in the position-loss sum:

| Value | Source | Notes |
|---|---|---|
| `iron` (default) | `lattice_layout.csv` (iron centres) | For non-bending magnets (KQUAD, KSEXT, KOCT) iron centre is on the orbit anyway. |
| `orbit` | `magnet_endpoints_lat{N}.csv` midpoints | Uses the (entry+exit)/2 point on the true reference orbit. Matters only when bending magnets enter the loss — for DL/JL slices the iron-centre and orbit-mid can differ by up to ~17 mm. |

### Dynamic loss reweighting (MC only)

When `dynamic_loss.enabled: true`, the three component weights rotate
through a sine pattern with equally spaced phases:

```
phase_i = 2π · i / 3   (i = 0:pos, 1:dir, 2:env)
w_i(t)  = 1 + A · cos(2π·t/T − phase_i)
L_dyn   = w_pos(t)·pos + w_dir(t)·dirWeight·dir + w_env(t)·envWeight·env
```

- `amplitude` (A) ∈ [0, 1]. 0 disables (static composite). 1 makes the
  weight oscillate from 0 to 2× base.
- `period_iterations` (T): if 0, auto = `5 × N_params`.
- The MC `accept` step uses `L_dyn` (so each iteration is judged on the
  *current* component balance); `bestLoss` / `loss` written to history use
  the **static** composite — keeps the progress curve readable.

### Optimization JSON config (minimal example)

```json
{
  "experiment_id": 4,

  "optimization_mode": "global",          // "global" | "windowed" | "gradient"
  "gradient": {                           // active only when mode == "gradient"
    "fd_method": "central",               // "central" (2M evals) or "forward" (M+1)
    "update_rule": "normalized_gd",       // "normalized_gd" or "adam"
    "n_iterations": 100,
    "learning_rate": 0.1,
    "fd_step_abs": 0.05,
    "fd_step_percent": 0.005,
    "max_halvings": 6
  },

  "n_iterations": 10000,                  // used by MC mode
  "lattice_max_gid": 207,
  "kick_every": 225,
  "kick_scale": 0.25,

  "loss": {
    "mode": "composite",                  // "composite" | "orbit" | "envelope"
    "magnet_sp_numbers": [1, 2, 3],
    "position_reference": "iron",         // "iron" | "orbit"
    "exclude_device_types": ["DP_L_PM", "DQ_I1_1A", "..."],
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
      "exclude_device_types": [
        "DP_L_PM", "DP_sh_PM",            // permanent dipoles (kept via correct_coeff)
        "SP_I2_1A", "DQ_I3_1A",           // correction windings — frozen
        "corrector_I1_1A", "corrector_I2_1A"   // correctors — frozen
      ],
      "value_column": "arg_val",
      "step_abs": 0.1
    }
  ],

  "per_element": {                        // windowed-mode params (block name kept for compat)
    "window_size": 3,
    "window_stride": 1,
    "iters_per_window": 100,              // ignored when iters_per_param > 0
    "passes": 1,
    "start_from_sp_start": true,          // propagate from global start
    "accumulate_groups": false,            // groups grow within a pass when true
    "iters_per_param": 0,                  // > 0 => m_iter = this × N_active_groups
    "loss_mode": "endpoint",              // "endpoint" | "path_only" | "composite" | "composite_full"
    "path_weight": 0.0,
    "path_sample_stride": 10
  },

  "combined": {                            // windowed ⇄ global orchestrator (opt_combo)
    "cycles": 5,
    "pe_passes_per_cycle": 2,
    "mc_iters_per_cycle": 1000,
    "mc_loss_mode": "composite"
  },

  "dynamic_loss": {                        // sinusoidal weight rotation (MC only)
    "enabled": true,
    "amplitude": 0.7,
    "period_iterations": 0                 // 0 = auto (5 × N_params)
  }
}
```

### Tuning knobs

| Parameter | Meaning |
|---|---|
| `kick_every` | every N iterations, all groups are perturbed at once (MC mode). 0 disables. |
| `kick_scale` | kick amplitude in A (sign drawn at random) |
| `step_abs` | normal-iteration step for `arg_val` (one random group per iter) |
| `step_percent` | step for `correct_coeff` as a fraction of current value |
| `lattice_max_gid` | restrict tracking to first N gids of the layout |
| `gradient.learning_rate` | base step amplitude in gradient mode |
| `gradient.fd_step_abs` | FD perturbation for `arg_val` (default 0.01) |
| `gradient.fd_step_percent` | FD perturbation for `correct_coeff` (default 0.001) |
| `gradient.max_halvings` | line-search halving cap (typical 5–6) |
| `per_element.window_size` | magnets per sliding window |
| `per_element.iters_per_window` | hill-climb iterations per window (fixed mode) |
| `per_element.iters_per_param` | when > 0, replaces above by `iters_per_param × N_active_groups` |
| `per_element.passes` | number of full passes over the lattice |
| `per_element.start_from_sp_start` | propagate beam from (0,0,0) through prefix on every window |
| `per_element.accumulate_groups` | knobs from earlier windows stay tunable in later windows of the same pass |
| `per_element.loss_mode` | per-window loss formulation (`endpoint` / `path_only` / `composite` / `composite_full`) |
| `dynamic_loss.amplitude` | strength of sinusoidal reweighting (0 = static) |
| `dynamic_loss.period_iterations` | full rotation period; 0 = auto `5 × N_params` |
| `combined.cycles` | how many windowed ⇄ global alternations |

### Frozen device types (excluded from optimization)

By convention, the following are excluded from `parameter_groups_per_position`
and their `arg_val` is kept at 0 (handled outside the optimizer):

| Device type | Role | Why frozen |
|---|---|---|
| `corrector_I1_1A`, `corrector_I2_1A` | Steering correctors | Closed-orbit correction is a separate problem |
| `SP_I2_1A` | I2 winding of sextupoles | Auxiliary correction coil |
| `DQ_I3_1A` | I3 winding of dipole-quadrupoles | Auxiliary correction coil |

---

## Matrix model (alternative engine)

Alongside the field-map differential core (SynchrotronSim, a differential tracker with a symplectic Boris integrator), the project has a second, **matrix-model**
tracking engine — a fast symplectic *drift–kick–drift* tracker that reads the **same
database** but represents each element by physical multipole coefficients instead of a
3-D field map. It is ~64× faster than the Python reference and is used for lattice tuning
and whole-ring studies. Full mathematics: **`docs/` matrix-model algorithms document**.

**Canonical state vector** `(x, pₓ, y, p_y, t, δ)`: `pₓ,p_y` are transverse momenta
normalized to the reference momentum `p₀`, `δ = Δp/p₀`, `t` is path length (RF-ready),
`s` (arc length) is the independent variable. `p_full = p₀(1+δ)` is conserved inside
magnets; dispersion and chromaticity emerge from the δ-aware drift.

**Element model** — a rectangular box: centre from `lattice_configs (x,y,z)`, axis along
`yaw`, ends at centre ± `L_phys/2` (metal length, geometry), field strength integrated
over `L_eff` (`∫B·dl`). Co-located windings (SP, DQ, corrector) are merged into one magnet
(superposition). Kicks use the reference rigidity `q/p₀`: `θ₀=(q/p₀)·B·L_eff`, `k₁=(q/p₀)G`,
`k₂=(q/p₀)S`, `k₃=(q/p₀)O`. **Empty gaps** between elements are *modified drifts* — geometric
propagation from one exit plane to the next entry plane of **arbitrary orientation**
(non-parallel), with a momentum-direction check.

### Files

| File | Role |
|---|---|
| `matrix_tracker.py` / `src/MatrixTracker.h` | Tracker (Python reference / C++ header-only port, numerically identical) |
| `src/matrix_tracker_main.cpp` → `MatrixTracker.exe` | C++ tracker CLI + timing |
| `matrix_optimizer.py` / `src/matrix_optimizer_main.cpp` → `MatrixOptimizer.exe` | Successive-approximation tuner (survival loss, line search) |
| `lte_warmstart.py` | Design-K warm start from `.lte` (preserves the `K·L` integral) |
| `warmstart_ring.py` | Write design-K `mat_*` to the whole ring (incl. junction) |
| `build_knobmap.py` | Knob map: shared standard-SP knobs + unique junction knobs |
| `apply_knobmap.py` | Write tuned knobs back to `lattice_configs.mat_*` |
| `beam_injection.py` | Inject a `particle_groups` beam from its own plane; one full turn |
| `matrix_traj_export.py` | Export matrix trajectories to `visualize.py` format |
| `build_closed_ring.py` | Build a closed 40-superperiod ring (lattice 140) |

DB additions (`build_matrix_db.py`): `eff_length`/`phys_length` in `hardware_library`,
`mat_B/G/S/O` in `lattice_configs`, and tables `element_transforms`, `coupling_matrix`,
`element_connectivity`, `matrix_knobs`, `matrix_tune`.

### Tuning pipeline (full ring with injector/junction)

```bash
py -3.14 warmstart_ring.py 101          # design-K warm start (standard SP + junction)
py -3.14 build_knobmap.py 101           # 98 knobs: 32 shared standard + 66 junction
MatrixOptimizer.exe 101 60 0.03         # tune (survival of a ±1 mm beam), ~25 s
py -3.14 apply_knobmap.py 101           # tuned coefficients -> DB
py -3.14 matrix_traj_export.py 42 101   # export trajectory (or inject a beam, below)
py -3.14 visualize.py 42                # 3-D view
```

Result: the **full ring (38 standard superperiods + 1 junction = 360°)** transports a
±1 mm working beam through all 2082 magnets, junction included, with no losses.

### Beam injection from a plane

Beams live in `particle_groups.csv` (same table as the differential model). On selection the
particles are checked to be **coplanar** (⊥ their velocity); that plane is embedded into the
matrix ring, its origin set where the reference orbit (line between the nearest elements'
exit/entry planes) crosses it, direction taken from the velocity, and the ring is closed
so that tracking runs one full turn and stops back at the plane.

```bash
py -3.14 beam_injection.py 101 1                     # inject group 1, report full turn
py -3.14 matrix_traj_export.py inject 43 101 1       # export injected trajectory
py -3.14 visualize.py 43
```

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
- `build/SynchrotronTracker[.exe]` — field-map differential tracker (SynchrotronSim; symplectic Boris integrator)
- `build/SynchrotronOptimizer[.exe]` — field-map optimizer
- `build/MatrixTracker[.exe]` — matrix-model tracker (see [Matrix model](#matrix-model-alternative-engine))
- `build/MatrixOptimizer[.exe]` — matrix-model tuner

The two matrix binaries are header-only (`src/MatrixTracker.h`, `src/Workspace.h`) and use
only the C++17 standard library + OpenMP.

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
python lte_to_lattice.py                # .lte -> sp_elements.csv + lattice_layout.csv
python lte_to_design_orbit.py           # .lte -> design_orbit_segments.csv
python assemble_lattice_configs.py      # produces lattice_configs.csv

# 2. (once) convert Ansys .fld -> .bin
build/SynchrotronTracker --convert         # Windows: build\SynchrotronTracker.exe --convert

# 3. run experiment 4 (SP1 + 50 m drift demo)
build/SynchrotronTracker 4
python visualize.py 4

# 4. run optimization (opt resolves <exp_id> via experiments.csv)
opt 4 --method=global                   # global mode; search = line_search | monte_carlo
opt 4 --method=windowed                 # windowed (sliding-window) mode
opt 4 --method=gradient                 # gradient descent
opt_combo 4                             # alternates windowed and global per opt_combo cycle
opt 4 --method=pe --W=5 --passes=3      # CLI overrides without editing JSON

# 5. after optimization — rebuild lattice & re-render
python assemble_lattice_configs.py
build/SynchrotronTracker 4
python visualize.py 4

# 6. inspect convergence and stability
python viz_optim_history.py workspace/results/sp1_long_drift.history.json
python viz_mc_evolution.py        # parameter evolution from latest run
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
Plots loss, `mean_mm`, `max_mm` over iterations and **final** parameter values
(bar chart).

### `viz_mc_evolution.py [history.json] [--top N]`
Focused on **how loss fell and which parameters changed during optimization**:
- Loss vs iteration with kicks marked as vertical bars and BEST marks
- **Evolution of top-N most active parameters as line plots**
  (shows arg_val changes over the course of the run, not just final values)
- Histogram of accepted step sizes `|new - old|`
- Bar chart of activity per group (number of ACCEPT'ed proposals)
- Sidebar with run summary

If no argument given, picks the most recent `*.history.json` in
`workspace/results/` automatically.

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

Tracking time, 5 particles, a representative multi-SP run (~290k pusher steps):

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

**Rule of thumb: bigger kicks but rarer.** A large `kick_scale` with a long
`kick_every` (e.g. `kick_scale=1.0`, `kick_every=500`) explores the parameter
space more effectively than small frequent kicks. Between kicks the
hill-climber has enough budget to **settle** into the new basin; if kicks come
too often, the climber spends most of its time rolling back the previous
disturbance and never converges. Small frequent kicks tend to behave like
noise: lots of REJECTs, little progress. A single bold kick every few hundred
iterations is what actually moves the solution between basins.

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
