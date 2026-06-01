# Algorithm Description

**Project:** Synchrotron Particle Tracker (SynchrotronSim)
**Document type:** Algorithm description
**Standard:** RD 50-34.698-90 (methodological guidelines to GOST 34),
section "Description of algorithm"
**Revision:** 2026-06-01

---

## Foreword

This document describes the computational models implemented in the
`SynchrotronTracker` and `SynchrotronOptimizer` programs. Each algorithm
is documented per RD 50-34.698-90 in five subsections:

1. **Purpose and characteristics** — what task is solved, type of method,
   accuracy/stability properties.
2. **Input information** — quantities consumed, sources, units, ranges.
3. **Solution results** — quantities produced, where they are stored.
4. **Mathematical description** — equations and formulas implemented.
5. **Solution algorithm** — step-by-step procedure executed in code.

References to source files use the convention `path:line` (e.g.,
`src/Integrator.h:111`). All physical quantities are in SI units.

---

## Contents

1. [Boris pusher — relativistic symplectic integrator](#1-boris-pusher)
2. [Trilinear field interpolation with YZ symmetry](#2-trilinear-field-interpolation)
3. [Spatial bucket grid for magnet lookup](#3-spatial-bucket-grid)
4. [Magnet contribution with bounding-sphere prescreening](#4-magnet-contribution)
5. [Finish-plane sign-change detection](#5-finish-plane-detection)
6. [Distance from point to design polyline](#6-distance-to-polyline)
7. [Composite loss function](#7-composite-loss-function)
8. [Stochastic optimizer (hill-climber with kicks)](#8-stochastic-optimizer)
9. [LTE → lattice generation](#9-lte-to-lattice)
10. [Sub-lattice assembly](#10-sub-lattice-assembly)
11. [Envelope analysis and betatron tune via FFT](#11-envelope-and-fft)

---

## 1. Boris pusher

**Source:** `src/Integrator.h:53` (`Integrator::borisStep`).

### 1.1 Purpose and characteristics

Single-step relativistic integrator of the equations of motion for a
charged particle in arbitrary electromagnetic field. The Boris pusher
uses drift–kick–drift splitting and is:

- **Symplectic** — preserves phase-space volume.
- **Exactly conserving |p|²** for purely magnetic fields (the rotation
  step is mathematically orthogonal).
- **Time-symmetric (DKD)**, formal accuracy O(dt²).
- **Cheap** — exactly one field evaluation per step (vs four in RK4).

References: Boris, *Proc. Fourth Conf. on Numerical Simulation of
Plasmas*, 1970; Birdsall & Langdon, *Plasma Physics via Computer
Simulation*, 1985, §4-3.

### 1.2 Input information

| Symbol | Meaning | Units | Source |
|---|---|---|---|
| `s.pos` | particle position at step start | m | `ParticleState` |
| `s.mom` | particle momentum `p = γ·m·v` | kg·m/s | `ParticleState` |
| `m0` | rest mass | kg | particle_groups.csv |
| `q` | charge | C | particle_groups.csv |
| `lattice` | accelerator (provides `getTotalFields`) | — | lattice_configs.csv |
| `dt` | timestep | s | experiments.csv |

### 1.3 Solution results

| Symbol | Meaning | Where written |
|---|---|---|
| `s.pos` | position after one step | in-place update |
| `s.mom` | momentum after one step | in-place update |

### 1.4 Mathematical description

Equations of motion:

```
dx/dt = v
dp/dt = q · (E + v × B)
p = γ · m · v,   γ = √(1 + p²/(m²c²))
```

In the current lattice all devices are magnetic, so `E ≡ 0`. The
E-kick steps below reduce to no-ops; the general form is retained for
forward compatibility.

Boris splits one step into seven half-operations:

```
1) half drift:    x_{n+½} = x_n + (p_n / (γ_n · m)) · dt/2
2) sample fields: B, E ← getTotalFields(x_{n+½})
3) half E-kick:   p⁻ = p_n + q · E · dt/2
4) γ⁻ = √(1 + |p⁻|² / (m²c²))
5) magnetic rotation:
       t = q · B · dt / (2 · γ⁻ · m)
       s = 2 · t / (1 + |t|²)
       p' = p⁻ + (p⁻ × t)
       p⁺ = p⁻ + (p' × s)       (|p⁺| = |p⁻|, exact)
6) half E-kick:   p_{n+1} = p⁺ + q · E · dt/2
7) half drift:    x_{n+1} = x_{n+½} + (p_{n+1} / (γ_{n+1} · m)) · dt/2
```

### 1.5 Solution algorithm

```
INPUT  : s = (pos, mom), m0, q, lattice, dt
OUTPUT : s updated to step n+1

1.  γ ← sqrt(1 + |mom|² / (m0²·c²))
2.  pos ← pos + mom · (dt/2) / (m0·γ)             ; half-drift
3.  (B, E) ← lattice.getTotalFields(pos)
4.  p_minus ← mom + E · q · dt/2                   ; half E-kick
5.  γ_minus ← sqrt(1 + |p_minus|² / (m0²·c²))
6.  t ← B · q · dt / (2·m0·γ_minus)
7.  s ← 2·t / (1 + |t|²)
8.  p_prime ← p_minus + cross(p_minus, t)
9.  p_plus  ← p_minus + cross(p_prime,  s)         ; rotation
10. mom ← p_plus + E · q · dt/2                    ; half E-kick
11. γ ← sqrt(1 + |mom|² / (m0²·c²))
12. pos ← pos + mom · (dt/2) / (m0·γ)              ; half-drift
```

Earlier revisions kept an analytic-drift helper (`Integrator::driftStep`)
and a predicate (`AcceleratorConfig::isDriftRegion`) for switching to
`x += v · dt` on segments without magnets. Both were removed in the
2026-06 refactor — the spatial bucket index already collapses drift-region
cost to a single empty-bucket lookup, so the separate code path provided
no measurable benefit.

---

## 2. Trilinear field interpolation

**Source:** `src/FieldMap.h:41` (`FieldMap::getInterpolatedField`).

### 2.1 Purpose and characteristics

Evaluation of vector field B (or E) at an arbitrary point inside a
magnet's local coordinate system, given a precomputed field on a regular
Cartesian grid (Ansys output). Trilinear interpolation gives C⁰
continuity and is exact at grid nodes. Optional mirroring across the YZ
plane (when `useSymmetry == true`) halves the storage of fields that are
known to be antisymmetric in x.

### 2.2 Input information

| Symbol | Meaning | Units | Source |
|---|---|---|---|
| `p` | query point in local frame | m | caller |
| `useSymmetry` | enable YZ-mirroring | bool | hardware_library.csv |
| `header.{minX..maxZ}` | grid bbox | m | binary header |
| `header.{stepX, stepY, stepZ}` | grid spacing | m | binary header |
| `header.{nx, ny, nz}` | grid size | int | binary header |
| `data[i·ny·nz + j·nz + k]` | field at node (i,j,k) | T (or V/m) | .bin |

### 2.3 Solution results

| Symbol | Meaning | Type |
|---|---|---|
| return value | interpolated field at `p` | `FieldVector` |

Returns `{0,0,0}` for points outside the grid bbox.

### 2.4 Mathematical description

Let `(i, j, k)` be the lower grid indices and `(tx, ty, tz)` the
fractional offsets in [0, 1]. Trilinear interpolation:

```
F(p) =  v000·(1-tx)(1-ty)(1-tz) + v100·tx·(1-ty)(1-tz)
      + v010·(1-tx)·ty·(1-tz)   + v110·tx·ty·(1-tz)
      + v001·(1-tx)(1-ty)·tz    + v101·tx·(1-ty)·tz
      + v011·(1-tx)·ty·tz       + v111·tx·ty·tz
```

YZ-symmetry mirror rule (when `p.x < 0` and `useSymmetry`):
sample at `(-p.x, p.y, p.z)`, then flip `F.x ← -F.x` (antisymmetric Bx,
symmetric By, Bz).

### 2.5 Solution algorithm

```
INPUT : p, useSymmetry, grid (header + data)
OUTPUT: field vector at p

1.  mirrored ← false
2.  if useSymmetry and p.x < 0:
       p.x ← -p.x; mirrored ← true
3.  if p outside [minX..maxZ]: return {0,0,0}
4.  fx ← (p.x - minX)/stepX;   i ← floor(fx);   tx ← fx - i
    fy ← (p.y - minY)/stepY;   j ← floor(fy);   ty ← fy - j
    fz ← (p.z - minZ)/stepZ;   k ← floor(fz);   tz ← fz - k
5.  clamp i, j, k to [0, n-2] (handle p == max boundary)
6.  fetch 8 corner values v000..v111 from data[]
7.  trilerp by tx, ty, tz, component-wise (see §2.4)
8.  if mirrored: F.x ← -F.x
9.  return F
```

---

## 3. Spatial bucket grid

**Source:** `src/LatticeManager.h:179` (`AcceleratorConfig::buildSpatialIndex`),
`src/LatticeManager.h:235` (`getTotalFields`).

### 3.1 Purpose and characteristics

Accelerates the per-step lookup of magnets that can contribute to the
field at the current particle position. Without indexing, the simulator
would iterate over all ~2762 elements per step. The bucket grid reduces
this to ~5 candidates per step; in drift regions the bucket is empty
and the per-step cost collapses to a single bucket lookup with no field
evaluations.

The index is built once after `addElement(...)` finishes, then used
read-only throughout tracking and optimization.

### 3.2 Input information

| Symbol | Meaning | Source |
|---|---|---|
| `elements[i].worldCenter` | sphere center | `computeBounds()` |
| `elements[i].worldRadius` | sphere radius | `computeBounds()` |

### 3.3 Solution results

| Symbol | Meaning |
|---|---|
| `bucketSize` | edge length of one bucket (m) |
| `bxMin, byMin` | lower index of grid in world units / bucketSize |
| `bxN, byN` | grid dimensions |
| `bucketGrid[bx·byN + by]` | list of element indices intersecting bucket |

The query `getTotalFields(pos)` returns `(B, E)` summed only over the
bucket containing `pos`.

### 3.4 Mathematical description

Projection of bounding spheres onto the X-Y plane (the synchrotron ring
lies in z ≈ 0):

```
bucketSize = max(1.0, 2 · max_i worldRadius_i)
```

This guarantees every element's projected disk fits inside a 2×2 group
of buckets, so it is registered at most 4 times.

Bucket coordinates of a query point:
```
bx = floor(x / bucketSize) - bxMin
by = floor(y / bucketSize) - byMin
```

### 3.5 Solution algorithm

**Build phase:**
```
1. maxR ← max_i elements[i].worldRadius
2. bucketSize ← max(1.0, 2 · maxR)
3. xLo, xHi ← min/max over elements of (worldCenter.x ± worldRadius)
   yLo, yHi ← analogous for y
4. bxMin ← floor(xLo/bucketSize);  byMin ← floor(yLo/bucketSize)
   bxN   ← floor(xHi/bucketSize) - bxMin + 1
   byN   ← floor(yHi/bucketSize) - byMin + 1
5. bucketGrid ← array of empty vectors of size bxN·byN
6. for each element i:
       bxLo..bxHi ← bucket span on x (worldCenter.x ± worldRadius)
       byLo..byHi ← bucket span on y
       for bx in [bxLo..bxHi], by in [byLo..byHi]:
           bucketGrid[bx·byN + by].push_back(i)
```

**Query phase (per particle step):**
```
1. bx ← floor(pos.x/bucketSize) - bxMin
   by ← floor(pos.y/bucketSize) - byMin
2. if (bx, by) outside grid: B, E ← 0; return
3. for each index i in bucketGrid[bx·byN + by]:
       elements[i].addContribution(pos, B, E)
```

Drift regions hit an empty bucket and the inner loop is empty, so no
field evaluations occur.

---

## 4. Magnet contribution

**Source:** `src/LatticeManager.h:114` (`MagneticElement::addContribution`).

### 4.1 Purpose and characteristics

Adds the contribution of a single magnet to the running field sum
`(B, E)` at a global query point. Performs a fast bounding-sphere
prescreen (7 flops) before the full local-frame transform and trilinear
interpolation.

### 4.2 Input information

| Symbol | Meaning | Source |
|---|---|---|
| `globalPos` | world-frame query point | tracking loop |
| `worldCenter, worldRadius²` | bounding sphere | `computeBounds()` |
| `transform` | local↔world (rotation + translation) | lattice_configs.csv |
| `map` | field grid | hardware_library.csv |
| `useSymmetry` | YZ symmetry flag (read from `use_symmetry` column) | hardware_library.csv |
| `correctCoeff` | calibration factor | hardware_library + sp_elements |
| `fieldDirection` | sign flip (±1) | hardware_library.csv |
| `currentScale` | current/nominal_current (active electromagnets) | computed at load |
| `isActive` | true ↔ apply currentScale | hardware_library.csv |
| `type` | `MAGNETIC` / `ELECTRIC` (read from `field_type` column: `B`/`M`/`MAGNETIC` → MAGNETIC; `E`/`ELECTRIC` → ELECTRIC; missing → MAGNETIC) | hardware_library.csv |

`addContribution` routes by `type` — `MAGNETIC` accumulates into
`totalB`, `ELECTRIC` into `totalE`. All current devices are magnetic, so
`totalE` stays at zero in practice and the E-kick steps of the Boris
pusher are no-ops, but the routing is live and supports adding
electrostatic elements without code changes.

### 4.3 Solution results

In-place addition into `totalB` (if `type == MAGNETIC`) or `totalE`
(if `type == ELECTRIC`). No return value.

### 4.4 Mathematical description

```
if |globalPos − worldCenter|² > worldRadius²:        return            (prescreen)
localPos = R⁻¹ · (globalPos − T)                                       (world → local)
localF   = trilinear(map, localPos, useSymmetry)
if localF == 0:                                       return            (sparse grid)
globalF  = R · localF                                                   (local → world)
scale    = correctCoeff · fieldDirection · (currentScale if isActive)
totalB or totalE += scale · globalF                                     (by type)
```

### 4.5 Solution algorithm

```
1. dx,dy,dz ← globalPos - worldCenter
   if dx²+dy²+dz² > worldRadius²:  return
2. localPos ← transform.transformPointInv(globalPos)
3. localF   ← map.getInterpolatedField(localPos, useSymmetry)
4. if localF == 0:  return
5. globalF  ← transform.rotateVector(localF)
6. scale    ← correctCoeff · fieldDirection
   if isActive:  scale *= currentScale
7. if type == MAGNETIC:  totalB += scale · globalF
   else                  totalE += scale · globalF
```

---

## 5. Finish-plane detection

**Source:** `src/FinishPlane.h:51` (`FinishPlane::hasCrossed`); setup at
`src/FinishPlane.h:30`.

### 5.1 Purpose and characteristics

Detects the moment a particle crosses an oriented plane defined by three
points. Designed for per-particle (race-free under OpenMP) state. Allows
skipping the first `skipSteps` calls to bypass earlier intersections of
the same plane with the orbit (e.g., when a half-ring plane also passes
near the start).

### 5.2 Input information

| Symbol | Meaning | Source |
|---|---|---|
| `a, b, c` | three points defining the plane | experiments.csv |
| `startPos` | initial particle position | particle_groups.csv |
| `skip` | number of initial calls to bypass | experiments.csv `max_turns` |
| `pos` | current particle position | tracking loop |

### 5.3 Solution results

`hasCrossed(pos)` returns true on the first call after the signed
distance changes sign (relative to `lastD`). Used by the tracker to stop
integration and write the trajectory.

### 5.4 Mathematical description

Plane normal:
```
n̂ = (b − a) × (c − a) / |(b − a) × (c − a)|
```
Sign convention is fixed so that `startPos` has negative signed distance:
```
n̂ ← −n̂ if (startPos − a) · n̂ > 0
```
Signed distance from arbitrary `pos`:
```
d(pos) = (pos − a) · n̂
```
Crossing condition between successive samples:
```
(lastD <  0 ∧ d ≥ 0) ∨ (lastD ≥ 0 ∧ d < 0)
```

### 5.5 Solution algorithm

**Setup (once per particle):**
```
1. u ← b − a;  v ← c − a;  n ← u × v
2. normalize n;  if degenerate, set n ← (1, 0, 0)
3. if (startPos − a) · n > 0:  n ← −n
4. callCount, lastD ← 0, 0;  initialized ← false
```

**Per-step query:**
```
1. callCount ← callCount + 1
2. if callCount ≤ skipSteps:  return false
3. d ← (pos − a) · n
4. if not initialized:  lastD ← d;  initialized ← true;  return false
5. crossed ← (lastD < 0 ∧ d ≥ 0) ∨ (lastD ≥ 0 ∧ d < 0)
6. lastD ← d
7. return crossed
```

---

## 6. Distance to polyline

**Source:** `src/optimizer.cpp:66` (`distanceToPolyline`).

### 6.1 Purpose and characteristics

Computes the minimum Euclidean distance from a point to a polyline made
of line segments. Used as the *envelope amplitude* in the envelope loss
(§7) and by `analyze_trajectory_stability.py` to compute the distance
from a particle trajectory to the design orbit.

### 6.2 Input information

| Symbol | Meaning |
|---|---|
| `p` | query point |
| `poly` | ordered list of polyline vertices (N ≥ 2) |

### 6.3 Solution results

Scalar distance ≥ 0.

### 6.4 Mathematical description

For each segment `[a, b]`, project `p` onto the line, clamp to `[0, 1]`:

```
t = clamp( ((p − a) · (b − a)) / |b − a|², 0, 1 )
c = a + t · (b − a)
d_i² = |p − c|²
```
Result:
```
dist(p, poly) = √(min_i d_i²)
```

### 6.5 Solution algorithm

```
INPUT : p, poly[0..N-1]
OUTPUT: dist ≥ 0

1. if N < 2:  return 0
2. best² ← +∞
3. for i in [0..N-2]:
     a ← poly[i];  b ← poly[i+1]
     u ← b − a;  len² ← u·u
     t ← clamp( (p − a)·u / len², 0, 1 )   if len² > 1e-30  else  0
     c ← a + t·u
     d² ← (p − c)·(p − c)
     if d² < best²:  best² ← d²
4. return √best²
```

---

## 7. Composite loss function

**Source:** `src/optimizer.cpp`, function `runTrackingAndLoss` (search by
that name).

### 7.1 Purpose and characteristics

Scalar objective minimized by the optimizer. Three additive components
target different physical defects:

| Component | Role |
|---|---|
| Position loss | particle must pass through magnet centers |
| Direction loss | beam must exit parallel to a reference direction |
| Envelope loss | betatron amplitude must not grow |

Position loss is updated *incrementally* in the tracking loop, with no
trajectory storage. Direction and envelope losses are evaluated only at
the final step of each particle.

### 7.2 Input information

| Symbol | Meaning | Source |
|---|---|---|
| trajectories (implicit) | the 5 particle paths | tracking |
| `magnet_sp_numbers` | which SPs participate in loss | optim JSON |
| `exclude_device_types` | exclude PM dipoles, DQ, … | optim JSON |
| `extra_target_magnets` | virtual targets `(x,y,z,weight)` | optim JSON |
| `direction_loss.reference` | unit vector r̂ | optim JSON |
| `direction_loss.weight` | w_dir | optim JSON |
| `envelope_loss.weight` | w_env | optim JSON |
| `design_polyline` | element centers + prelude/postlude | rebuilt each run |

### 7.3 Solution results

| Symbol | Meaning |
|---|---|
| `loss` | scalar to minimize |
| `info.{sum_d², mean_mm, max_mm, …}` | diagnostics for history.json |

### 7.4 Mathematical description

```
position_loss  = Σ_(i particle, m magnet)  w_m · min_t |x_i(t) − c_m|²
direction_loss = w_dir · Σ_i sin²(θ_i),
                 cos θ_i = (p_i_final · r̂) / |p_i_final|
envelope_loss  = w_env · Σ_i max(0, A_i_final − A_i_init)²,
                 A_i(t) = distanceToPolyline(x_i(t), design_polyline)

loss = position_loss + direction_loss + envelope_loss
```

Only **growth** of the envelope is penalized; betatron damping carries
zero penalty.

### 7.5 Solution algorithm

```
INPUT : config (weights, magnet filter, targets), lattice, particles
OUTPUT: loss, info

A) Build static state (once per call):
   1. magnets_idx ← filter elements by (sp_number ∈ S) ∧
                                       (device_type ∉ exclude)
   2. virtual_targets ← list of (c, w) from extra_target_magnets
   3. design_poly ← centers of elements 1..N + prelude
                    (extrapolate first segment by −2.94 m)
                    + postlude (extrapolate last segment by +length)

B) Per particle i:
   minDist²[m] ← +∞ for each m in magnets_idx ∪ virtual_targets
   A_init     ← distanceToPolyline(x_init,  design_poly)
   for each step until finish plane:
       Boris step
       for each m in candidates from bucket:
           d² ← |x − c_m|²
           if d² < minDist²[m]:  minDist²[m] ← d²
   A_final ← distanceToPolyline(x_final, design_poly)
   θ_i      ← angle(p_final, r̂)

C) Aggregate:
   position_loss  = Σ_(i,m) w_m · minDist²[i, m]
   direction_loss = w_dir · Σ_i sin²(θ_i)
   envelope_loss  = w_env · Σ_i max(0, A_final − A_init)²
   loss           = position_loss + direction_loss + envelope_loss

D) Fill info: sum_d², mean and max of √minDist² in mm, per-particle break.
```

If any particle fails to reach the finish plane within the 50 M-step
cap, the function returns *failure* and the optimizer rejects that
candidate state.

---

## 8. Stochastic optimizer

**Source:** `src/optimizer.cpp`, `main()` and helpers.

### 8.1 Purpose and characteristics

Single-threaded greedy hill-climber with **periodic full-vector kicks**
to escape local minima. The optimizer mutates magnet currents (`arg_val`)
and per-instance calibration multipliers (`correct_coeff`) declared in
`parameter_groups_per_position` of the optim JSON.

Properties:
- **Asymmetric step:** at each iteration one random group is perturbed
  by `±step_abs` (current) or `±step_percent · current` (coefficient).
- **Kicks:** every `kick_every` iterations *all* groups are kicked at
  once by `±kick_scale`. Sign chosen by current loss gradient sign of
  a Gaussian draw.
- **Greedy acceptance:** accept if new loss < best, else revert.
- **Checkpointing:** dump full history every `CHECKPOINT_EVERY` iters.

### 8.2 Input information

| Symbol | Meaning | Source |
|---|---|---|
| `experiment_id` | which experiment to load | optim JSON |
| `lattice_max_gid` | track only first N gids | optim JSON |
| `n_iterations` | total iterations | optim JSON |
| `kick_every` | period of full-vector kicks | optim JSON |
| `kick_scale` | kick amplitude (A or fraction) | optim JSON |
| `parameter_groups_per_position` | which (sp_type, sp_eid) and column | optim JSON |
| `step_abs` / `step_percent` | per-group step | optim JSON |
| `loss.*` | weights and filters (see §7.2) | optim JSON |

### 8.3 Solution results

| Output | Format | Path |
|---|---|---|
| Best `arg_val` per (sp_type, sp_eid) | CSV row update | `sp_elements.csv` |
| Best `correct_coeff` per device_type | CSV row update | `hardware_library.csv` |
| Iteration history (loss, params) | JSON | `results/{name}.history.json` |
| Final trajectories of best state | binary | `results/trajectories_exp{id}.bin` |

### 8.4 Mathematical description

Let `θ ∈ ℝ^M` be the vector of mutable parameters, `L(θ)` the loss
function from §7, and `g ∈ {1..G}` an iteration's randomly-chosen group.

```
Regular iter:  θ'_g = θ_g + ε_g · sign(rand() − 0.5),
               θ'_k = θ_k  for k ≠ g
               where ε_g = step_abs       (if group's column = arg_val)
                  or ε_g = step_percent · θ_g  (if column = correct_coeff)

Kick iter:     θ'   = θ + δ, where  δ_g ~ Uniform(−kick_scale,
                                                  +kick_scale) for each g

Acceptance:    θ  ← θ'  if L(θ') < L(θ_best);  L(θ_best) ← L(θ')
                else θ unchanged
```

### 8.5 Solution algorithm

```
1. Load configuration:
     - read JSON, build OptimConfig
     - read sp_elements, hardware_library, lattice_configs, particle_groups
     - filter lattice rows by lattice_id and lattice_max_gid
2. Build accelerator:
     - load field maps into a shared cache (one .bin per device_type)
     - for each row: create MagneticElement, computeBounds()
     - acc.buildSpatialIndex()
3. Build parameter list:
     - for each entry in parameter_groups_per_position:
         find matching (sp_type, sp_eid) -> collect element indices
         store reference to mutable field (arg_val or correct_coeff)
4. Initial evaluation:
     L_best ← runTrackingAndLoss(θ)
     if failure: abort
5. Loop  i = 1..n_iterations:
     if i % kick_every == 0:
         for each group g:  θ_g ← θ_g + Uniform(−kick_scale, +kick_scale)
     else:
         pick random group g
         θ_g_old ← θ_g
         θ_g     ← θ_g + step · sign(rand − 0.5)
     L_new ← runTrackingAndLoss(θ)
     if L_new is failure or L_new ≥ L_best:
         revert this iteration's mutations
     else:
         L_best ← L_new
         append (i, L_best, snapshot θ) to history
     if i % CHECKPOINT_EVERY == 0:
         write history.json (atomic via temp file)
6. Final phase:
     - write best θ back to sp_elements.csv and hardware_library.csv
     - rebuild lattice_configs.csv (preserve sub-lattices)
     - run one final tracking pass, dump trajectories_exp{id}.bin
     - write final history.json
```

---

## 9. LTE to lattice generation

**Source:** `lte_to_lattice.py`.

### 9.1 Purpose and characteristics

Parses Elegant-format LTE files (`SYLA_desc.lte`, `SYLA_line.lte`),
classifies elements per `lte_mapping.json`, lays out the ring along its
chord, and writes the editable templates `sp_elements.csv` and the
geometric unrolling `lattice_layout.csv`.

Special features:
- **DL/JL straightening** — long bend magnets are split into 5 slices
  placed on the chord rather than along the curved arc.
- **arg_val preservation** — existing `arg_val` values in
  `sp_elements.csv` are kept on re-runs, so the user does not lose
  optimization results.

### 9.2 Input information

| Path | Content |
|---|---|
| `workspace/lte/SYLA_desc.lte` | device definitions (type, length, K1, angles) |
| `workspace/lte/SYLA_line.lte` | element sequence around the ring |
| `workspace/configs/lte_mapping.json` | LTE type → device_type_id, regex, thresholds |

### 9.3 Solution results

| Path | Content |
|---|---|
| `workspace/tables/sp_elements.csv` | templates `(sp_type, sp_element_id, device_type_id, arg_val, correct_coeff)` |
| `workspace/tables/lattice_layout.csv` | unrolling `(global_id, sp_number, sp_type, sp_element_id, x, y, z, yaw, pitch, roll)` |

### 9.4 Mathematical description

For a bend with length `L` and bend angle `α`, the chord length is
`Lc = 2 · (L/α) · sin(α/2)`. Five slices are placed at fractions
`0, 0.25, 0.5, 0.75, 1.0` of the chord, each with the local yaw rotation
that matches the chord direction at that point.

Element placement walks the ring as a sequence of segments. For each
segment of length `ℓ` and yaw turn `Δψ`:
```
yaw'  = yaw + Δψ
x'    = x + ℓ · cos(yaw + Δψ/2)
y'    = y + ℓ · sin(yaw + Δψ/2)
```

### 9.5 Solution algorithm

```
1. Read JSON config -> mapping, regex patterns, classification thresholds
2. Parse SYLA_desc.lte:
     for each line "NAME: TYPE, key1=val1, key2=val2, ..."
         classify TYPE via config (DL/JL/SBM/DP/DQ/QP/SQ/COR/MARK/...)
         store (NAME, classified_type, length, angle, K1, ...)
3. Parse SYLA_line.lte:
     resolve nested LINE expansions -> flat sequence of (NAME, repeat)
4. Classify each element into (sp_type, sp_element_id):
     standard SP repeats 1..N (typical superperiod)
     junction SP at SP18-21 (geometry differs)
5. Walk the ring:
     for each element in flat sequence:
         if DL/JL: emit 5 slices on the chord
         else:     emit 1 slice
         update (x, y, yaw) per §9.4
6. Load existing sp_elements.csv if present:
     preserve arg_val for matching (sp_type, sp_eid) keys
7. Write sp_elements.csv (templates only, one row per unique key)
8. Write lattice_layout.csv (full unrolling)
```

---

## 10. Sub-lattice assembly

**Source:** `assemble_lattice_configs.py`.

### 10.1 Purpose and characteristics

Joins the editable template (`sp_elements.csv`) with the geometric
unrolling (`lattice_layout.csv`) to produce the final simulator input
`lattice_configs.csv`. Also generates **shortened lattice copies** (IDs
102, 103, 104, …) from `sub_lattices.json`, which contain only the
first `max_gid` rows of the master lattice. This lets experiments run on
a sub-ring without re-walking the LTE.

### 10.2 Input information

| Path | Content |
|---|---|
| `workspace/tables/sp_elements.csv` | templates (manually editable) |
| `workspace/tables/lattice_layout.csv` | geometry (from `lte_to_lattice.py`) |
| `workspace/configs/lte_mapping.json` | `lattice_id` of the master lattice |
| `workspace/configs/sub_lattices.json` | list of `{lattice_id, max_gid, description}` |

### 10.3 Solution results

| Path | Content |
|---|---|
| `workspace/tables/lattice_configs.csv` | one row per (lattice, instance), columns `lattice_id, instance_id, device_type_id, x, y, z, yaw, pitch, roll, arg_val, correct_coeff` |

### 10.4 Mathematical description

No arithmetic — pure relational join:
```
output =  { (LATTICE_ID, layout.global_id, sp.device_type_id,
             layout.x..roll, sp.arg_val, sp.correct_coeff)
          | layout ∈ lattice_layout,
            sp ∈ sp_elements where (sp.sp_type, sp.sp_element_id) =
                                    (layout.sp_type, layout.sp_element_id) }
        ∪ { copy of master rows with lattice_id ← sub.lattice_id and
            instance_id ≤ sub.max_gid
          | sub ∈ sub_lattices }
```

### 10.5 Solution algorithm

```
1. sp_dict ← load_sp_elements()              ; key (sp_type, sp_eid)
2. output_rows ← []
3. for each row r in lattice_layout.csv:
     key ← (r.sp_type, r.sp_element_id)
     if key not in sp_dict:  warn and skip
     dev, arg, coeff ← sp_dict[key]
     append row {LATTICE_ID, r.global_id, dev, r.x..roll, arg, coeff}
4. for each (sub_id, max_gid) in SUB_LATTICES:
     for r in output_rows where r.instance_id ≤ max_gid:
         append copy with lattice_id ← sub_id
5. Write CSV with retry-on-PermissionError (Windows file lock).
```

---

## 11. Envelope analysis and FFT betatron tune

**Source:** `analyze_trajectory_stability.py`.

### 11.1 Purpose and characteristics

Post-processing of a trajectory binary to evaluate **dynamic stability**:
- Build the envelope `A_i(t) = dist(x_i(t), design_polyline)` for each
  of the 5 particles.
- Fit linear and exponential growth models to `A_i(t)`, qualified by R².
- Run FFT of `A_i(t)` to detect betatron oscillation frequencies (tune).
- Emit a verdict per particle: STABLE / OSCILLATING / GROWING / DRIFTED.

### 11.2 Input information

| Path | Content |
|---|---|
| `workspace/results/trajectories_exp{id}.bin` | 5 trajectories |
| `workspace/tables/lattice_configs.csv` | element centers → design polyline |
| `workspace/tables/experiments.csv` | which lattice was used |

### 11.3 Solution results

Per particle:
- Envelope plot `A(t)` (matplotlib).
- Linear-fit slope `m_lin` and intercept `A₀_lin`.
- Exponential-fit growth rate `γ_exp` and amplitude `A₀_exp` with R².
- FFT spectrum and dominant frequency `f_betatron`.
- Verdict label.

### 11.4 Mathematical description

Envelope:
```
A_i(t_k) = distanceToPolyline( x_i(t_k), design_poly )
```
Linear fit:
```
A(t) ≈ m · t + A₀                  (least squares)
```
Exponential fit (on log of positive part):
```
log A(t) ≈ γ · t + log A₀          ⇒  A(t) ≈ A₀ · exp(γ · t)
R²_exp = 1 − SS_res / SS_tot
```
FFT, sample step `dt_sample = t_k − t_{k−1}`:
```
Â(f) = FFT(A(t)),
f_betatron = argmax_(f > f_min) |Â(f)|     (DC bin masked)
```

### 11.5 Solution algorithm

```
1. Load trajectories (5 lists of (Ni, 3) points).
2. Determine lattice_id from experiments.csv, load lattice_configs,
   build design_poly (centers + prelude/postlude).
3. For each particle i:
     a. A_i ← [distanceToPolyline(p, design_poly) for p in traj_i]
     b. t_i ← uniform step (dt × traj_stride × stride_inside_bin)
     c. linear fit:  m, A₀ ← polyfit(t_i, A_i, deg=1)
     d. exponential fit (only points with A > A_floor):
            log A = γ·t + log A₀;  compute R²
     e. FFT(A_i);  mask DC bin;  pick max as f_betatron
     f. classify:
            DRIFTED      if A_init ≈ 0 and m_lin > δ
            GROWING      if γ > γ_threshold and R² ≥ R²_min
            OSCILLATING  if dominant peak in FFT and γ ≤ γ_threshold
            STABLE       otherwise
4. Render plots, print verdict table.
```

---

## Appendix A. Configuration files

| File | Consumed by | Purpose |
|---|---|---|
| `workspace/configs/lte_mapping.json` | `lte_to_lattice.py`, `assemble_lattice_configs.py` | LTE type mapping + `lattice_id` |
| `workspace/configs/sub_lattices.json` | `assemble_lattice_configs.py`, `src/optimizer.cpp` | sub-lattice definitions |
| `workspace/configs/visualizer.json` | `visualize.py`, `live_visualizer.py` | camera, colors, opacity, projection |
| `workspace/optim_inputs/*.json` | `prepare_optim_config.py` → `SynchrotronOptimizer[.exe]` | optimizer setup per experiment |
| `_common.py` | all top-level Python scripts | cross-platform helpers (executable suffix, project root, matplotlib backend) |

## Appendix B. Data files

| File | Format | Description |
|---|---|---|
| `workspace/fields_bin/*.bin` | binary `GridHeader + FieldVector[nx·ny·nz]` | per-magnet field map |
| `workspace/tables/hardware_library.csv` | `;`-separated CSV, 10 columns | per `device_type` calibration |
| `workspace/tables/sp_elements.csv` | `;`-separated CSV | superperiod template (editable) |
| `workspace/tables/lattice_layout.csv` | `;`-separated CSV | ring unrolling (generated) |
| `workspace/tables/lattice_configs.csv` | `;`-separated CSV | simulator input (generated) |
| `workspace/tables/experiments.csv` | `;`-separated CSV | experiment definitions |
| `workspace/tables/particle_groups.csv` | `;`-separated CSV | initial conditions |
| `workspace/results/trajectories_exp{id}.bin` | binary | output of `SynchrotronTracker[.exe]` |
| `workspace/results/{name}.history.json` | JSON | output of `SynchrotronOptimizer[.exe]` |

Columns of `hardware_library.csv` (in file order):
`device_type_id`, `Active`, `field_type`, `map_raw_path`, `map_bin_path`,
`stl_path`, `use_symmetry`, `nominal_arg`, `field_direction`,
`correct_coeff`. The two flags `field_type` (`B`/`E`) and
`use_symmetry` (`1`/`0`) are read into `MagneticElement::type` and
`MagneticElement::useSymmetry` at load time. All current rows are
`B` / `1`, but the routing in `addContribution` and the mirror branch
in `FieldMap::getInterpolatedField` honour both flags, so introducing
an electrostatic device or an asymmetric field map needs only a CSV
edit.
