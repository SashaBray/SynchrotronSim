SYNCHROTRON PARTICLE TRACKER
============================

Relativistic tracking of charged particles through a synchrotron magnet
lattice. Uses real Ansys field maps, a symplectic Boris pusher for
integration, and a stochastic optimizer for tuning magnet currents and
per-instance calibration coefficients.

Cross-platform: builds and runs on Windows (MSVC) and Linux/macOS
(gcc/clang).

GitHub: https://github.com/SashaBray/SynchrotronSim


-------------------------------------------------------------------------------
REPOSITORY CONTENTS
-------------------------------------------------------------------------------

Source code:
    src/                            C++ sources
        main.cpp                    SynchrotronTracker entry
        optimizer.cpp               SynchrotronOptimizer entry
        Integrator.h                Boris pusher (relativistic, symplectic)
        LatticeManager.h            spatial bucket index, magnets, transforms
        FieldMap.h                  trilinear interpolation
        FinishPlane.h               sign-change finish detection
        AnsysConverter.h            .fld -> .bin converter
    CMakeLists.txt                  cmake config (MSVC / gcc / clang)

Python tools (top-level, all cross-platform):
    _common.py                      shared helpers (exe paths, project root,
                                    matplotlib backend selection)
    lte_to_lattice.py               LTE -> sp_elements + lattice_layout
    assemble_lattice_configs.py     sp_elements + layout -> lattice_configs
    prepare_optim_config.py         JSON -> .optimconfig
    apply_optim_history.py          history.json -> CSV (state rollback)
    visualize.py                    static 3D viewer (PyVista)
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
                                    (10 columns: device_type_id, Active,
                                    field_type, map_raw_path, map_bin_path,
                                    stl_path, use_symmetry, nominal_arg,
                                    field_direction, correct_coeff)
                                    -- field_type ("B"/"E") and use_symmetry
                                    are placeholders for future electric /
                                    asymmetric devices; current rows are
                                    all B / 1, the code routes by them.
    tables/experiments.csv          experiment definitions
    tables/particle_groups.csv      particle initial conditions
    optim_inputs/*.json             optimizer configs per experiment

Documentation:
    README.md                       full documentation (markdown)
    README.txt                      this file (plain text mirror)
    docs/algorithms.md              algorithm reference (RD 50-34.698-90)
    docs/algorithms.docx            Russian GOST-formatted version
    docs/build_docx.js              docx generator script


-------------------------------------------------------------------------------
NOT INCLUDED IN THE REPOSITORY (created locally)
-------------------------------------------------------------------------------

Generated from LTE on first run:
    workspace/tables/sp_elements.csv      SP templates (standard/junction)
    workspace/tables/lattice_layout.csv   ring unrolling
    workspace/tables/lattice_configs.csv  final DB for the simulator

Large binaries / build artifacts:
    build/                                C++ objects, executables
    workspace/fields_bin/                 binary field maps
    workspace/fields_raw/*.fld            raw Ansys maps (can be gigabytes)
    workspace/results/                    trajectories, history.json, .png
    workspace/geometry/*.stl              STL magnet models (for viz)
    venv/                                 Python virtual environment

NOTE: since field maps (.bin/.fld) are NOT stored in the repository,
running a simulation requires obtaining them separately (from Ansys
calculations or from the project author).


-------------------------------------------------------------------------------
PREPARING TO RUN
-------------------------------------------------------------------------------

1. C++ BUILD

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release -j

   Produces (.exe suffix only on Windows):
    build/SynchrotronTracker      particle tracker
    build/SynchrotronOptimizer    parameter optimizer

   Toolchain:
    Windows   MSVC 2019+ (Build Tools), CMake 3.10+
    Linux     gcc 9+ or clang 10+, CMake, OpenMP
    macOS     clang with libomp (brew install libomp)

2. PYTHON ENVIRONMENT

    python -m venv venv
    venv\Scripts\activate                 (Windows)
    source venv/bin/activate              (Linux/macOS)
    pip install pandas numpy matplotlib pyvista

3. FIELD MAPS

   Put raw Ansys maps into workspace/fields_raw/ (filenames must match
   the map_raw_path column in workspace/tables/hardware_library.csv).

   Then convert them to binary format:

    build/SynchrotronTracker --convert        (Linux/macOS)
    build\SynchrotronTracker.exe --convert    (Windows)

   This creates workspace/fields_bin/*.bin.

4. GENERATE DB FROM LTE

   sp_elements.csv and lattice_layout.csv are fully generated from LTE
   files:

    python lte_to_lattice.py
    python assemble_lattice_configs.py

   lte_to_lattice.py:
     - reads workspace/lte/SYLA_desc.lte and SYLA_line.lte
     - generates workspace/tables/sp_elements.csv and lattice_layout.csv
     - preserves existing arg_val on subsequent runs

   assemble_lattice_configs.py:
     - assembles workspace/tables/lattice_configs.csv
     - additionally creates sub-lattices 102/103/104 (listed in
       workspace/configs/sub_lattices.json)


-------------------------------------------------------------------------------
RUNNING AN EXPERIMENT
-------------------------------------------------------------------------------

Experiments live in workspace/tables/experiments.csv. Preconfigured:

    exp 1   sanity test (plane x=5, lattice 101)
    exp 2   past SP20 -- half ring (lattice 101)
    exp 3   full-ring diagnostics (lattice 101)
    exp 4   SP3 + 2 SP of free drift (lattice 102)
    exp 5   SP1-10 + drift to SP11 (lattice 103)
    exp 6   junction tuning, SP18-21, group 3 (lattice 104)

Run the tracker for experiment N:

    build/SynchrotronTracker N            (Linux/macOS)
    build\SynchrotronTracker.exe N        (Windows)

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

    python prepare_optim_config.py workspace/optim_inputs/sp3_drift_2sp.json

   Creates workspace/optim_inputs/sp3_drift_2sp.optimconfig (flat format).

2. Run the optimization:

    build/SynchrotronOptimizer workspace/optim_inputs/sp3_drift_2sp.optimconfig

   Knobs (kick_every, kick_scale, step_abs, n_iterations, loss-component
   weights) are set in the JSON config.

3. Inspect convergence curves:

    python viz_optim_history.py workspace/results/sp3_drift_2sp.history.json

4. Apply the resulting state (only if the optimizer was interrupted):

    python apply_optim_history.py workspace/results/sp3_drift_2sp.history.json

   (on completion the optimizer itself writes best values back into
   sp_elements.csv and hardware_library.csv; apply_optim_history is
   only needed for rolling back or restoring from a checkpoint)


-------------------------------------------------------------------------------
COMPOSITE OPTIMIZER LOSS
-------------------------------------------------------------------------------

Loss = position_loss + direction_loss + envelope_loss

1. Position loss
   Sum over (magnet, particle) pairs: w_m * d^2, where d is the minimum
   distance from the trajectory to the magnet centre. Computed
   incrementally in the tracking loop, no trajectory storage. Magnets
   filtered by magnet_sp_numbers and exclude_device_types. Additionally:
   TARGET points (virtual magnets) with custom weight.

2. Direction loss
   w_dir * sum over particles sin^2(theta_i), where theta_i is the
   angle between the particle's final momentum and a reference
   direction (orbit tangent at the finish point).

3. Envelope loss
   w_env * sum over particles max(0, A_final - A_init)^2, where
   A(t) = distance to the design polyline (centres of lattice elements
   + prelude/postlude along the tangent). Only growth is penalized
   (betatron damping is a useful effect).


-------------------------------------------------------------------------------
PHYSICAL MODEL (BRIEF)
-------------------------------------------------------------------------------

Equations of motion:
    dx/dt = v
    dp/dt = q * (E + v x B)
where p = gamma * m * v is the relativistic momentum.

The integrator is a Boris pusher (drift-kick-drift splitting). It is
symplectic and exactly preserves |p|^2 for pure magnetic fields. 1 field
eval / step (vs 4 in RK4). In the current lattice all devices are
magnetic, so E = 0 and the E-kick steps reduce to no-ops.

The magnet field is computed via trilinear interpolation in the local
coordinate system, then rotated back into world coordinates. The sum of
fields from magnets near the position is computed via a 2D spatial
bucket grid + bounding-sphere prescreening (~5 magnets per step instead
of 2762).


-------------------------------------------------------------------------------
CROSS-PLATFORM NOTES
-------------------------------------------------------------------------------

  * Executable suffix.  _common.exe_name(stem) returns stem + ".exe" on
    Windows and stem on Linux/macOS. Python scripts use this helper --
    do not hardcode .exe.
  * Workspace discovery.  Both C++ binaries call findWorkspace(argv[0]),
    which walks up from the binary, argv[0] parent and the cwd looking
    for a directory containing workspace/. Python scripts use
    _common.project_root() (walks up to a directory with CMakeLists.txt
    or README.md).
  * Matplotlib backend.  _common.set_matplotlib_backend("TkAgg") falls
    back to Qt5Agg or Agg if Tk is unavailable, so headless Linux works.
  * Paths.  Python uses pathlib.Path / forward slashes uniformly;
    Windows accepts both.
  * OpenMP on macOS.  Install libomp via Homebrew; Apple clang ships
    without OpenMP.


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
  * Incremental loss -- minDist^2[i][m] updated in the tracking loop
  * Composite loss with weights -- TARGETs pin final convergence


-------------------------------------------------------------------------------
TROUBLESHOOTING
-------------------------------------------------------------------------------

"lattice 104 has no rows in lattice_configs.csv"
    Run python assemble_lattice_configs.py -- the script automatically
    regenerates sub-lattices.

"ERROR: baseline tracking failed"
    At least one particle fails to reach the finish plane within 50M
    steps. Likely causes:
    - arg_val gives an unstable orbit on the long path
    - max_turns (skipSteps) too small -- the plane is crossed multiple
      times along the way, extra crossings need to be skipped
    Run SynchrotronTracker N directly to see where the beam is lost.

"Access is denied" when launching the .exe (Windows)
    Defender may delete a freshly built .exe. Add the build/ folder to
    Windows Defender exclusions.

Particles start right inside the first SP lens
    Check particle_groups.csv -- position should be ~2.94 m before the
    first magnet along the orbit tangent (prelude drift).

Optimization "gets stuck"
    Play with the balance of step_abs and kick_scale in the JSON config:
    - kick_scale > step_abs: aggressive kicks jump out of local minima
    - kick_scale ~ step_abs: conservative kicks, precise convergence
    For final polishing, reduce step_abs (1.0 -> 0.1 A).

Analyzer shows stale data
    After optimization, re-run the tracker:
        python assemble_lattice_configs.py
        build/SynchrotronTracker N
        python analyze_trajectory_stability.py N
