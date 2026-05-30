SYNCHROTRON PARTICLE TRACKER
============================

Relativistic tracking of charged particles through a synchrotron magnet
lattice. Uses real Ansys field maps, a symplectic Boris pusher for
integration, and a stochastic optimizer for tuning magnet currents and
calibration coefficients.

GitHub: https://github.com/SashaBray/SynchrotronSim


-------------------------------------------------------------------------------
REPOSITORY CONTENTS
-------------------------------------------------------------------------------

Source code:
    src/                            C++ sources
        main.cpp                    SynchrotronTracker
        optimizer.cpp               SynchrotronOptimizer
        Integrator.h                Boris pusher (relativistic, symplectic)
        LatticeManager.h            spatial bucket index, magnets, transforms
        FieldMap.h                  trilinear interpolation
        FinishPlane.h               sign-change finish detection
        AnsysConverter.h            .fld -> .bin converter
    CMakeLists.txt                  C++ build

Python tools (top-level):
    lte_to_lattice.py               LTE -> sp_elements + lattice_layout
    assemble_lattice_configs.py     sp_elements + layout -> lattice_configs
    prepare_optim_config.py         JSON -> .optimconfig for the optimizer
    apply_optim_history.py          history.json -> CSV (state rollback)
    visualize.py                    3D trajectory viewer (PyVista)
    live_visualizer.py              live-reload viewer
    viz_optim_history.py            loss/params plots from the optimizer
    analyze_trajectory_stability.py envelope analysis + FFT betatron tunes
    find_orbit_deviation_peaks.py   local max(|particle-orbit|)
    scan_field_maps.py              search/repair field map dropouts

Configurations (workspace/):
    configs/lte_mapping.json        LTE types -> hardware_library
    configs/sub_lattices.json       definitions of lattices 102/103/104
    configs/visualizer.json         camera, colors for PyVista
    lte/SYLA_desc.lte               device descriptions (Elegant LTE)
    lte/SYLA_line.lte               element sequence
    tables/hardware_library.csv     magnet calibration + paths to field maps
    tables/experiments.csv          experiment definitions
    tables/particle_groups.csv      particle initial conditions
    optim_inputs/*.json             optimizer configs for various experiments


-------------------------------------------------------------------------------
NOT INCLUDED IN THE REPOSITORY (created locally)
-------------------------------------------------------------------------------

Generated from LTE on first run:
    workspace/tables/sp_elements.csv      SP templates (standard/junction)
    workspace/tables/lattice_layout.csv   ring unrolling
    workspace/tables/lattice_configs.csv  final DB for the simulator

Large binaries / build artifacts:
    build/                                C++ objects, exe
    workspace/fields_bin/                 binary field maps
    workspace/fields_raw/*.fld            raw Ansys maps (can be gigabytes)
    workspace/results/                    trajectories, optim_history.json, .png
    workspace/geometry/*.stl              STL magnet models (for viz)
    venv/                                 Python virtual environment

NOTE: since field maps (.bin/.fld) are NOT stored in the repository,
running a simulation requires obtaining them separately (from Ansys
calculations or from the project author).


-------------------------------------------------------------------------------
PREPARING TO RUN
-------------------------------------------------------------------------------

1. C++ BUILD

    cmake -B build
    cmake --build build --config Release

   Produces:
    build/SynchrotronTracker.exe         particle tracker
    build/SynchrotronOptimizer.exe       parameter optimizer

2. PYTHON ENVIRONMENT

    python -m venv venv
    venv\Scripts\activate                 (Windows)
    source venv/bin/activate              (Linux/Mac)
    pip install pandas numpy matplotlib pyvista

3. FIELD MAPS

   Put raw Ansys maps into workspace/fields_raw/ (filenames must match
   the map_raw_path column in workspace/tables/hardware_library.csv).

   Then convert them to binary format:

    build\SynchrotronTracker.exe --convert

   This creates workspace/fields_bin/*.bin.

4. GENERATE DB FROM LTE

   The database (sp_elements.csv and lattice_layout.csv) is fully
   generated from LTE files. Run:

    python lte_to_lattice.py
    python assemble_lattice_configs.py

   lte_to_lattice.py:
     - reads workspace/lte/SYLA_desc.lte and SYLA_line.lte
     - generates workspace/tables/sp_elements.csv and lattice_layout.csv
     - preserves existing arg_val on subsequent runs

   assemble_lattice_configs.py:
     - assembles workspace/tables/lattice_configs.csv
     - additionally creates sub-lattices 102/103/104 (list from
       workspace/configs/sub_lattices.json)


-------------------------------------------------------------------------------
RUNNING AN EXPERIMENT
-------------------------------------------------------------------------------

The list of experiments lives in workspace/tables/experiments.csv. The
preconfigured experiments are:

    exp 1   sanity test (plane x=5, lattice 101)
    exp 2   past SP20 -- half ring (lattice 101)
    exp 3   full-ring diagnostics (lattice 101)
    exp 4   SP3 + 2 SP of free drift (lattice 102)
    exp 5   SP1-10 + drift to SP11 (lattice 103)
    exp 6   junction tuning, SP18-21, group 3 (lattice 104)

Run the tracker for experiment N:

    build\SynchrotronTracker.exe N

Creates workspace/results/trajectories_expN.bin (5 trajectories).

3D viewer:

    python visualize.py N

Stability analysis (envelope + FFT betatron tunes):

    python analyze_trajectory_stability.py N

Identify magnets near which the trajectory deviates the most:

    python find_orbit_deviation_peaks.py N


-------------------------------------------------------------------------------
OPTIMIZATION
-------------------------------------------------------------------------------

1. Prepare an optim-config from JSON:

    python prepare_optim_config.py workspace\optim_inputs\sp3_drift_2sp.json

   Creates workspace/optim_inputs/sp3_drift_2sp.optimconfig (flat format).

2. Run the optimization:

    build\SynchrotronOptimizer.exe workspace\optim_inputs\sp3_drift_2sp.optimconfig

   Parameters (kick_every, kick_scale, step_abs, number of iterations,
   loss-component weights) are set in the JSON config.

3. Inspect convergence curves:

    python viz_optim_history.py workspace\results\sp3_drift_2sp.history.json

4. Apply the resulting state (if the optimizer was interrupted):

    python apply_optim_history.py workspace\results\sp3_drift_2sp.history.json

   (on completion the optimizer itself writes best values into
   sp_elements.csv; apply_optim_history is only needed for rolling back
   or restoring from a checkpoint)


-------------------------------------------------------------------------------
COMPOSITE OPTIMIZER LOSS
-------------------------------------------------------------------------------

Loss = position_loss + direction_loss + envelope_loss

1. Position loss
   sum over (magnet, particle) pairs: w_m * d^2,
   where d is the minimum distance from the trajectory to the magnet
   center. Computed incrementally in the tracking loop, no trajectory
   storage. Magnets are filtered by magnet_sp_numbers and
   exclude_device_types. Additionally: TARGET points (virtual magnets)
   with their own weight.

2. Direction loss
   w_dir * sum over particles sin^2(theta_i),
   where theta_i is the angle between the particle's final momentum and
   the reference direction (tangent of the reference orbit at the finish
   point).

3. Envelope loss
   w_env * sum over particles max(0, A_final - A_init)^2,
   where A(t) = distance to the design polyline (centers of lattice
   elements + prelude/postlude along the tangent).
   Only growth is penalized (damping is a useful effect).


-------------------------------------------------------------------------------
PHYSICAL MODEL (BRIEF)
-------------------------------------------------------------------------------

Equations of motion:
    dx/dt = v
    dp/dt = q * (E + v x B)
where p = gamma * m * v is the relativistic momentum.

The integrator is a Boris pusher (drift-kick-drift splitting). It is
symplectic and exactly preserves |p|^2 for pure magnetic fields. 1 field
eval / step (vs 4 in RK4).

The magnet field is computed via trilinear interpolation in the local
coordinate system, then rotated back into world coordinates. The sum of
fields from magnets near the position is computed via a 2D spatial
bucket grid + bounding-sphere prescreening (~5 magnets per step instead
of 2762).


-------------------------------------------------------------------------------
PERFORMANCE
-------------------------------------------------------------------------------

Tracking time (5 particles, exp 4 = SP1-3 + drift, ~290k steps):

    Version                                 Tracking     Speedup
    Original RK4 + linear lattice scan      3424 ms      1x
    Boris pusher                             740 ms      4x
    + Spatial bucket index                    70 ms     47x
    + Incremental loss (optimizer)        ~110 ms/iter  ~30x

Key optimizations:
  * 2D spatial bucket grid -- ~5 magnets per step instead of 2762
  * Bounding-sphere prescreening -- early exit for far magnets
  * Boris pusher -- 1 field eval / step vs 4 in RK4
  * Incremental loss -- minDist^2[i][m] updated in the tracking loop,
    without storing the trajectory
  * Composite loss with weights -- TARGETs pin final convergence


-------------------------------------------------------------------------------
TROUBLESHOOTING
-------------------------------------------------------------------------------

"lattice 104 has no rows in lattice_configs.csv"
    Run python assemble_lattice_configs.py -- the script automatically
    regenerates sub-lattices from the config.

"ERROR: baseline tracking failed"
    At least one particle fails to reach the finish plane within 50M
    steps. Possible causes:
    - Current arg_val gives an unstable orbit on the long path
    - max_turns (skipSteps) is too small -- the plane is crossed
      multiple times along the way, extra crossings need to be skipped
    Run SynchrotronTracker.exe N to see where the beam is lost.

"Access is denied" when launching .exe (Windows)
    Defender may delete a freshly built .exe. Add the build/ folder to
    Windows Defender exclusions.

Particles start right inside the first SP lens
    Check particle_groups.csv -- position should be ~2.94 m before the
    first magnet along the orbit tangent (prelude drift).

Optimization "gets stuck"
    Play with the balance of step_abs and kick_scale in the JSON config:
    - kick_scale > step_abs: aggressive kicks jump out of local minima
    - kick_scale ~ step_abs: conservative kicks, precise convergence
    For final polishing, reduce step_abs (from 1.0 to 0.1 A).

Analyzer shows stale data
    After optimization, re-run the tracker:
        python assemble_lattice_configs.py
        build\SynchrotronTracker.exe N
        python analyze_trajectory_stability.py N
