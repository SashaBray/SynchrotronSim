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

**Lorentz equations of motion** for a relativistic charged particle:

$$
\frac{d\vec{x}}{dt} = \vec{v}, \qquad
\frac{d\vec{p}}{dt} = q\,\bigl(\vec{E} + \vec{v}\times\vec{B}\bigr)
\tag{1.1}
$$

Variables:
- $\vec{x}(t)$ — particle position in the global coordinate system, in metres.
- $\vec{v}(t)$ — particle velocity in the global coordinate system, in metres per second.
- $\vec{p}(t)$ — **relativistic momentum**, in kg·m/s; the quantity actually
  integrated by the code. It coincides with the classical $m\vec{v}$ only
  at low speeds; see definition below.
- $q$ — particle charge, in coulombs. For electrons $q = -1.602\!\times\!10^{-19}$ C.
- $m$ — particle rest mass, in kilograms. For electrons $m = 9.109\!\times\!10^{-31}$ kg.
- $\vec{E}, \vec{B}$ — electric and magnetic field vectors evaluated at the
  particle's current position, in V/m and tesla respectively.
- $t$ — physical (lab-frame) time, in seconds.

The cross product $\vec{v}\times\vec{B}$ is the **magnetic Lorentz force per
unit charge**: a magnetic field bends the trajectory but does no work, so
energy is conserved exactly.

**Relativistic momentum and Lorentz factor.** Momentum and velocity are
linked by

$$
\vec{p} = \gamma\,m\,\vec{v}, \qquad
\gamma = \frac{1}{\sqrt{1 - |\vec{v}|^2/c^2}} = \sqrt{1 + \frac{|\vec{p}|^2}{m^2 c^2}}
\tag{1.2}
$$

Variables:
- $\gamma \ge 1$ — **Lorentz factor**, a dimensionless quantity that measures
  how relativistic the motion is ($\gamma = 1$ at rest, $\gamma \to \infty$
  as $|\vec{v}|\to c$).
- $c = 299\,792\,458$ m/s — speed of light (exact, by definition of the metre).
- $|\vec{v}|^2 = v_x^2 + v_y^2 + v_z^2$ — squared speed (scalar).
- $|\vec{p}|^2 = p_x^2 + p_y^2 + p_z^2$ — squared momentum (scalar).

The second form (in terms of $|\vec{p}|$) is what the code evaluates each
step: $\gamma$ is recovered directly from $\vec{p}$ without ever needing
$\vec{v}$.

In the current lattice all devices are magnetic, so $\vec{E}\equiv\vec{0}$.
The E-kick steps below leave the momentum unchanged; the general form is retained
in case future electrostatic elements are added.

**Boris splitting.** One step $\Delta t$ is applied to the state
$(\vec{x}_n, \vec{p}_n) \mapsto (\vec{x}_{n+1}, \vec{p}_{n+1})$ as seven
half-operations:

*Step 1 — half-drift.* Position advances by half a step using the current
momentum:

$$
\vec{x}_{n+1/2} = \vec{x}_n + \frac{\vec{p}_n}{\gamma_n\,m}\,\frac{\Delta t}{2}
\tag{1.3}
$$

- $\Delta t$ — integration timestep, in seconds (`experiments.csv.dt`).
- $\gamma_n = \sqrt{1 + |\vec{p}_n|^2/(m^2 c^2)}$ — Lorentz factor at step $n$.

*Step 2 — sample fields* at the midpoint position $\vec{x}_{n+1/2}$:

$$
(\vec{B},\,\vec{E}) \leftarrow \mathrm{lattice.getTotalFields}(\vec{x}_{n+1/2})
\tag{1.4}
$$

This is the **only** field evaluation per step. RK4 would need four.

*Step 3 — first half of the electric kick:*

$$
\vec{p}^{-} = \vec{p}_n + \frac{q\,\vec{E}\,\Delta t}{2}
\tag{1.5}
$$

- $\vec{p}^{-}$ — momentum after the first E-kick, before the magnetic
  rotation. When $\vec{E}=\vec{0}$ this simplifies to $\vec{p}^{-}=\vec{p}_n$.

*Step 4 — re-evaluate $\gamma$* using the post-E-kick momentum:

$$
\gamma^{-} = \sqrt{1 + \frac{|\vec{p}^{-}|^2}{m^2 c^2}}
\tag{1.6}
$$

This $\gamma^{-}$ is the "instantaneous" relativistic factor used by the
magnetic rotation that follows.

*Step 5 — magnetic rotation.* Two auxiliary vectors are formed:

$$
\vec{t} = \frac{q\,\vec{B}\,\Delta t}{2\,\gamma^{-}\,m}, \qquad
\vec{s} = \frac{2\,\vec{t}}{1 + |\vec{t}|^2}
\tag{1.7}
$$

- $\vec{t}$ — half-angle rotation vector; its magnitude equals $\tan(\theta/2)$,
  where $\theta = qB\,\Delta t / (\gamma^{-} m)$ is the Larmor rotation angle.
- $\vec{s}$ — rescaled rotation vector chosen so that the two-step update
  below preserves the momentum magnitude exactly.

The momentum is then rotated in two stages:

$$
\vec{p}' = \vec{p}^{-} + \vec{p}^{-} \times \vec{t}, \qquad
\vec{p}^{+} = \vec{p}^{-} + \vec{p}' \times \vec{s}
\tag{1.8}
$$

- $\vec{p}'$ — intermediate momentum after the first half-rotation.
- $\vec{p}^{+}$ — momentum after the full magnetic rotation; it satisfies
  $|\vec{p}^{+}| = |\vec{p}^{-}|$ as an exact algebraic identity (this is
  what makes the Boris pusher conserve energy in a pure magnetic field).

*Step 6 — second half of the electric kick:*

$$
\vec{p}_{n+1} = \vec{p}^{+} + \frac{q\,\vec{E}\,\Delta t}{2}
\tag{1.9}
$$

This and Step 3 together implement a centred E-acceleration: each one
contributes half of the impulse $q\vec{E}\,\Delta t$.

*Step 7 — second half-drift:*

$$
\vec{x}_{n+1} = \vec{x}_{n+1/2} + \frac{\vec{p}_{n+1}}{\gamma_{n+1}\,m}\,\frac{\Delta t}{2}
\tag{1.10}
$$

- $\gamma_{n+1} = \sqrt{1 + |\vec{p}_{n+1}|^2/(m^2 c^2)}$ — Lorentz factor
  at the new step.

After all seven sub-operations the state has been advanced by exactly one
$\Delta t$ with second-order accuracy in time.

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

The field is stored at the nodes of a regular Cartesian grid. A query
point $\vec{p} = (p_x, p_y, p_z)$ in local coordinates lies inside a
unit cell of that grid; trilinear interpolation blends the eight
corner values.

**Fractional cell coordinates.** First locate the cell containing the
point:

$$
i = \left\lfloor \frac{p_x - x_{\min}}{\Delta x} \right\rfloor, \quad
t_x = \frac{p_x - x_{\min}}{\Delta x} - i
\tag{2.1}
$$

(analogous definitions for $j, t_y$ and $k, t_z$).

- $i,\,j,\,k$ — integer indices of the **lower** corner of the cell
  (zero-based; the upper corner is $(i\!+\!1, j\!+\!1, k\!+\!1)$).
- $t_x, t_y, t_z \in [0, 1]$ — fractional offsets inside the cell;
  $t_\alpha = 0$ means the point sits on the lower face along axis
  $\alpha$, $t_\alpha = 1$ on the upper face.
- $x_{\min},\,y_{\min},\,z_{\min}$ — coordinates of the grid corner with
  the smallest indices, in metres (from the binary header).
- $\Delta x,\,\Delta y,\,\Delta z$ — grid spacing along each axis, in
  metres (uniform per axis but independent across axes).

**Trilinear interpolation.** Let $\vec{v}_{abc}$ denote the stored field
vector at corner $(i\!+\!a,\,j\!+\!b,\,k\!+\!c)$ for
$a, b, c \in \{0, 1\}$. The interpolated field is

$$
\vec{F}(\vec{p}) = \sum_{a, b, c \in \{0, 1\}}
       w_{abc}(t_x, t_y, t_z)\;\vec{v}_{abc}
\tag{2.2}
$$

with the corner weight

$$
w_{abc}(t_x, t_y, t_z) =
   \bigl[(1 - t_x)^{1-a}\,t_x^{a}\bigr]\,
   \bigl[(1 - t_y)^{1-b}\,t_y^{b}\bigr]\,
   \bigl[(1 - t_z)^{1-c}\,t_z^{c}\bigr]
\tag{2.3}
$$

- $\vec{v}_{abc}$ — vector field at one of the eight corners (in tesla
  for $\vec{B}$, in V/m for $\vec{E}$).
- $w_{abc}$ — non-negative weight that decreases linearly with the
  distance from the query point to that corner. The eight weights
  always sum to one, so the result is a true weighted average.

At the corners ($t_\alpha \in \{0, 1\}$ for all $\alpha$) one weight is 1
and the rest are 0, so $\vec{F}$ equals the stored corner value exactly
— the interpolant is **node-exact**.

**YZ-symmetry mirror rule.** When the field map is known to be
antisymmetric in $x$ (i.e. only half the volume is stored,
$x \in [0, x_{\max}]$), a query at negative $x$ is folded into the stored
half:

$$
\vec{F}(p_x, p_y, p_z) \;=\;
\begin{cases}
   \vec{F}_{\text{raw}}(p_x, p_y, p_z) & \text{if } p_x \ge 0 \\
   \mathrm{diag}(-1, +1, +1)\;\vec{F}_{\text{raw}}(-p_x, p_y, p_z) & \text{if } p_x < 0
\end{cases}
\tag{2.4}
$$

- $\vec{F}_{\text{raw}}$ — interpolation performed on the stored
  ($p_x \ge 0$) half of the grid.
- $\mathrm{diag}(-1, +1, +1)$ — diagonal matrix that flips the
  $x$-component and leaves the other two untouched. This encodes the
  physical antisymmetry of $B_x$ across the YZ plane (whereas
  $B_y, B_z$ are symmetric).

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
| `bxMin, byMin` | lower index of grid in global units / bucketSize |
| `bxN, byN` | grid dimensions |
| `bucketGrid[bx·byN + by]` | list of element indices intersecting bucket |

The query `getTotalFields(pos)` returns `(B, E)` summed only over the
bucket containing `pos`.

### 3.4 Mathematical description

Every magnet has a precomputed **bounding sphere** with centre
$\vec{c}_i$ (in the global coordinate system) and radius $r_i$ that fully contains
its field map. The bucket index works in the $xy$-plane only, because
the ring is flat in $z \approx 0$.

**Bucket size.** Let $r_{\max}$ be the largest bounding-sphere radius
across all elements. The grid step is chosen as

$$
h = \max\bigl(1\,\mathrm{m},\;2\,r_{\max}\bigr)
\tag{3.1}
$$

- $h$ — edge length of one square bucket, in metres.
- $r_{\max} = \max_i r_i$ — largest bounding-sphere radius across all
  magnets, in metres. The factor of 2 guarantees that the
  $xy$-projection of any bounding circle fits inside a $2 \times 2$
  group of buckets, so each magnet is registered in at most 4 buckets.
- The lower bound of 1 m prevents pathologically small buckets when all
  magnets are tiny — the bucket array would otherwise become huge.

**Bucket coordinates.** A global-frame point $(x, y, z)$ is mapped to bucket
indices $(b_x, b_y)$ by

$$
b_x = \left\lfloor \frac{x}{h} \right\rfloor - b_x^{\min}, \qquad
b_y = \left\lfloor \frac{y}{h} \right\rfloor - b_y^{\min}
\tag{3.2}
$$

- $\lfloor\cdot\rfloor$ — floor function (round down to the nearest
  integer).
- $b_x^{\min}, b_y^{\min}$ — global integer offsets chosen at build time
  so that all element-occupied buckets have non-negative indices. They
  are computed from the lowest $x$ and $y$ coordinates reached by any
  bounding sphere.

A magnet $i$ is registered in every bucket $(b_x, b_y)$ that intersects
the disk
$\{(x, y) : (x - c_{i,x})^2 + (y - c_{i,y})^2 \le r_i^2\}$;
by construction this is at most a 2×2 block of buckets.

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
| `globalPos` | query point in the global frame | tracking loop |
| `worldCenter, worldRadius²` | bounding sphere | `computeBounds()` |
| `transform` | local↔global (rotation + translation) | lattice_configs.csv |
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
pusher contribute nothing in practice, but the routing is live and supports adding
electrostatic elements without code changes.

### 4.3 Solution results

In-place addition into `totalB` (if `type == MAGNETIC`) or `totalE`
(if `type == ELECTRIC`). No return value.

### 4.4 Mathematical description

Each magnet has its own local coordinate frame in which the field map
was computed. The contribution of magnet $i$ at a global-frame query
point $\vec{x}$ is built in four geometric steps and a scaling.

**Step 1 — bounding-sphere prescreen.** Reject the magnet immediately if
the query point is outside its sphere of influence:

$$
|\vec{x} - \vec{c}_i|^2 > r_i^2 \;\;\Rightarrow\;\; \text{contribution is 0}
\tag{4.1}
$$

- $\vec{c}_i$ — centre of magnet $i$'s bounding sphere in the global
  coordinate system, in metres. Computed from the magnet's transform and field
  map bounding box.
- $r_i$ — bounding-sphere radius of magnet $i$, in metres.

Squaring both sides avoids a square root in the hot loop.

**Step 2 — global → local transform.** The global-frame point is mapped to the
magnet's local frame:

$$
\vec{x}_\text{loc} = R_i^{\!\top}\,(\vec{x} - \vec{T}_i)
\tag{4.2}
$$

- $\vec{T}_i$ — translation: the magnet's local-origin position in the
  global coordinate system, in metres.
- $R_i$ — rotation matrix (3×3, orthogonal, $R_i^{\!\top}R_i = I$) that
  takes a local vector to the global coordinate system. Built from the magnet's
  $(\text{yaw},\,\text{pitch},\,\text{roll})$ Euler angles in
  `lattice_configs.csv`.
- $R_i^{\!\top}$ is the **inverse** rotation (transpose of an
  orthogonal matrix equals its inverse), so this line undoes the
  magnet's placement.

**Step 3 — sample the local field** via trilinear interpolation
(§ 2.4):

$$
\vec{F}_\text{loc} = \mathrm{trilinear}(\mathrm{map}_i,\,\vec{x}_\text{loc},\,\text{useSymmetry}_i)
\tag{4.3}
$$

- $\mathrm{map}_i$ — magnet $i$'s field map (binary grid of vectors).
- $\text{useSymmetry}_i \in \{\text{true, false}\}$ — whether the YZ
  mirror rule is applied (see § 2.4).

If $\vec{F}_\text{loc} = \vec{0}$ (the query falls in a zero-padded part
of the grid) we return early without doing the rest.

**Step 4 — local → global rotation** of the field vector:

$$
\vec{F}_\text{global} = R_i\,\vec{F}_\text{loc}
\tag{4.4}
$$

The translation does **not** apply here — a magnetic field is a vector,
not a point, so only the rotation part of the transform matters.

**Step 5 — calibration scaling.** Multiply by all per-magnet
calibration factors:

$$
\text{scale}_i =
   c_i^{\text{hw}} \cdot c_i^{\text{sp}}
   \cdot d_i
   \cdot \begin{cases}
      I_i / I_i^{\text{nom}} & \text{if } i \text{ is active} \\
      1 & \text{otherwise}
   \end{cases}
\tag{4.5}
$$

- $c_i^{\text{hw}}$ — device-level calibration coefficient
  (`hardware_library.csv.correct_coeff`), dimensionless. One value per
  device type.
- $c_i^{\text{sp}}$ — per-instance calibration coefficient
  (`sp_elements.csv.correct_coeff`), dimensionless. One value per
  $(sp\_type, sp\_element\_id)$ pair. Defaults to 1.0 if blank.
- $d_i \in \{+1, -1\}$ — direction flip
  (`hardware_library.csv.field_direction`). Used when a magnet's
  physical orientation differs from the one assumed in the field map.
- $I_i$ — actual current through the magnet, in amperes
  (`sp_elements.csv.arg_val`). Only used for active electromagnets;
  permanent magnets carry no current and the factor is 1.
- $I_i^{\text{nom}}$ — nominal current at which the field map was
  computed, in amperes (`hardware_library.csv.nominal_arg`). The ratio
  $I_i / I_i^{\text{nom}}$ rescales the field linearly with current
  (valid for the unsaturated regime).

**Step 6 — accumulate into the right field sum:**

$$
\vec{B}_\text{total} \mathrel{+}= \text{scale}_i\,\vec{F}_\text{global} \quad
   \text{if } \text{type}_i = \text{MAGNETIC}
\tag{4.6}
$$
$$
\vec{E}_\text{total} \mathrel{+}= \text{scale}_i\,\vec{F}_\text{global} \quad
   \text{if } \text{type}_i = \text{ELECTRIC}
\tag{4.7}
$$

- $\vec{B}_\text{total},\,\vec{E}_\text{total}$ — running sums of
  magnetic and electric field at the query point, in tesla and V/m.
- $\text{type}_i$ — `MAGNETIC` (default, all current devices) or
  `ELECTRIC` (placeholder for future electrostatic elements), from
  `hardware_library.csv.field_type`.

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

The finish plane is given by three global-frame points $\vec{a},\,\vec{b},\,
\vec{c}$ that are assumed not to be collinear. They span the plane
itself; the algorithm computes a unit normal and uses it for an
oriented "signed distance" test.

**Plane normal.** Two edge vectors are formed and crossed:

$$
\hat{n} \;=\; \frac{(\vec{b} - \vec{a}) \times (\vec{c} - \vec{a})}{\bigl|(\vec{b} - \vec{a}) \times (\vec{c} - \vec{a})\bigr|}
\tag{5.1}
$$

- $\vec{a},\,\vec{b},\,\vec{c}$ — three plane-defining points from
  `experiments.csv`, in metres.
- $\hat{n}$ — **unit normal** to the plane (length 1). The cross
  product gives a vector orthogonal to the plane; the denominator
  normalises its length.

**Sign convention.** The orientation of $\hat{n}$ from a cross product
is arbitrary, so it is flipped to guarantee the particle starts on the
"behind the plane" side:

$$
\hat{n} \leftarrow -\hat{n} \quad \text{if } (\vec{x}_0 - \vec{a}) \cdot \hat{n} > 0
\tag{5.2}
$$

- $\vec{x}_0$ — initial particle position, in metres
  (`particle_groups.csv`).
- The dot product on the left is the signed projection of
  $\vec{x}_0 - \vec{a}$ onto $\hat{n}$. After this fix, the particle
  starts at a negative signed distance by construction.

**Signed distance** of an arbitrary point $\vec{x}$ from the plane:

$$
d(\vec{x}) = (\vec{x} - \vec{a}) \cdot \hat{n}
\tag{5.3}
$$

- $d > 0$ means $\vec{x}$ is on the same side as $\hat{n}$ points.
- $d < 0$ means $\vec{x}$ is on the opposite side.
- $|d|$ equals the Euclidean distance from $\vec{x}$ to the closest
  point on the plane.

**Crossing detection.** Between two consecutive samples $\vec{x}_{n-1}$
and $\vec{x}_n$ (with signed distances $d_{n-1}$ and $d_n$), a crossing
is declared when the sign flips:

$$
\bigl(d_{n-1} < 0 \;\wedge\; d_n \ge 0\bigr) \;\vee\;
\bigl(d_{n-1} \ge 0 \;\wedge\; d_n < 0\bigr)
\tag{5.4}
$$

The half-open form ($d_n \ge 0$ vs $d_n < 0$) avoids the degenerate
case $d_n = 0$ being missed.

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

The polyline is a sequence of line segments. The distance from a query
point to the polyline is the minimum distance to any single segment;
the per-segment problem is a standard projection with clamping.

**Per-segment projection.** Consider segment $[\vec{a},\,\vec{b}]$. The
point on the *infinite* line that is closest to query $\vec{p}$ has
parameter

$$
t^{*} = \frac{(\vec{p} - \vec{a}) \cdot (\vec{b} - \vec{a})}{|\vec{b} - \vec{a}|^2}
\tag{6.1}
$$

- $\vec{a},\,\vec{b}$ — endpoints of the segment, in metres.
- $\vec{p}$ — query point, in metres.
- $t^{*}$ — line parameter; $t^{*} = 0$ gives $\vec{a}$, $t^{*} = 1$
  gives $\vec{b}$, $t^{*} \in (0, 1)$ a point in between, and values
  outside $[0, 1]$ lie on the extension of the segment.
- $|\vec{b} - \vec{a}|^2$ — squared length of the segment (skipped if
  $\le 10^{-30}$ to handle a degenerate zero-length segment).

To restrict the closest point to the actual segment (not its extension)
we clamp $t^{*}$ to the unit interval:

$$
t = \mathrm{clip}(t^{*},\,0,\,1) = \min\bigl(1,\,\max(0,\,t^{*})\bigr)
\tag{6.2}
$$

The closest point on segment $[\vec{a}, \vec{b}]$ to $\vec{p}$ is then

$$
\vec{c} = \vec{a} + t\,(\vec{b} - \vec{a}), \qquad
d_i = |\vec{p} - \vec{c}|
\tag{6.3}
$$

- $\vec{c}$ — foot of the perpendicular projection (or the closest
  endpoint, when the perpendicular falls outside the segment).
- $d_i$ — Euclidean distance from $\vec{p}$ to segment $i$, in metres.

**Distance to the polyline** is the minimum over all $N\!-\!1$ segments:

$$
\mathrm{dist}(\vec{p},\,\Pi) = \min_{i = 0,\,\ldots,\,N-2}\,d_i
   \;=\; \sqrt{\,\min_{i}\,|\vec{p} - \vec{c}_i|^2\,}
\tag{6.4}
$$

- $\Pi = (\vec{p}_0, \vec{p}_1, \ldots, \vec{p}_{N-1})$ — the polyline
  (an ordered list of $N$ vertices).
- The square root is taken only once at the end; per-segment work
  compares squared distances, which is both faster and avoids the
  square-root branch in the hot inner loop.

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

The total objective is a sum of three terms:

$$
\mathcal{L} \;=\; \mathcal{L}_{\text{pos}} \;+\; \mathcal{L}_{\text{dir}} \;+\; \mathcal{L}_{\text{env}}
\tag{7.1}
$$

- $\mathcal{L}$ — scalar loss minimised by the optimizer (units: m²
  before weighting; the weights re-scale the individual terms to make
  them comparable).
- The three terms are described below in detail.

**Position loss** — particles should pass close to magnet centres
(and to user-specified virtual targets):

$$
\mathcal{L}_{\text{pos}}
   \;=\; \sum_{i = 1}^{N_p}\,\sum_{m \in \mathcal{M}} w_m\,\min_{t}\,\bigl|\vec{x}_i(t) - \vec{c}_m\bigr|^2
\tag{7.2}
$$

- $N_p$ — number of tracked particles (5 in this project: a reference
  particle plus 4 satellites).
- $\mathcal{M}$ — set of "scoring magnets" plus virtual targets, after
  filtering by `magnet_sp_numbers` and `exclude_device_types` from the
  optimizer JSON.
- $\vec{c}_m$ — geometric centre of magnet $m$ in the global coordinate system,
  in metres.
- $w_m$ — per-target weight (default 1; targets in
  `extra_target_magnets` carry their own larger weights, typically 20
  or 50, to "pin" final convergence).
- $\vec{x}_i(t)$ — position of particle $i$ at integration time $t$,
  in metres.
- $\min_t\,|\vec{x}_i(t) - \vec{c}_m|^2$ — the **closest squared
  distance** along the entire trajectory of particle $i$ to magnet $m$.
  Computed incrementally during tracking — only the running minimum is
  stored, never the full trajectory.

**Direction loss** — at the finish plane, particle momenta should be
parallel to a user-specified reference direction:

$$
\mathcal{L}_{\text{dir}}
   \;=\; w_{\text{dir}}\,\sum_{i = 1}^{N_p}\,\sin^2 \theta_i, \qquad
\cos\theta_i = \frac{\vec{p}_i^{\,\text{fin}}\cdot \hat{r}}{|\vec{p}_i^{\,\text{fin}}|}
\tag{7.3}
$$

- $w_{\text{dir}}$ — weight of the direction term (typically $10^{3}$);
  large because the angle squared is dimensionless and small in
  magnitude.
- $\vec{p}_i^{\,\text{fin}}$ — momentum of particle $i$ at the moment
  it crosses the finish plane, in kg·m/s.
- $\hat{r}$ — **unit reference direction** specified in the optimizer
  config (`direction_loss.reference`); usually the orbit tangent at the
  finish point.
- $\theta_i$ — angle between $\vec{p}_i^{\,\text{fin}}$ and $\hat{r}$,
  in radians.
- $\sin^2\theta_i = 1 - \cos^2\theta_i$ — quadratic in $\theta$ for
  small angles, so the term behaves like a least-squares angular error.

**Envelope loss** — betatron oscillations should not grow during the
pass through the lattice:

$$
\mathcal{L}_{\text{env}}
   \;=\; w_{\text{env}}\,\sum_{i = 1}^{N_p}\,\max\bigl(0,\;A_i^{\text{fin}} - A_i^{\text{init}}\bigr)^2
\tag{7.4}
$$

with

$$
A_i(t) = \mathrm{dist}\bigl(\vec{x}_i(t),\,\Pi_{\text{design}}\bigr)
\tag{7.5}
$$

- $w_{\text{env}}$ — envelope weight, typically $10^{6}$ (large, because
  $A$ is in metres and we want a metre-scale grow penalised harshly).
- $A_i(t)$ — **envelope amplitude** of particle $i$ — its distance from
  the design polyline at time $t$, computed by the algorithm in § 6.
- $\Pi_{\text{design}}$ — design orbit polyline (centres of lattice
  elements, plus a prelude/postlude segment along the orbit tangent).
  Rebuilt at the start of each optimization run.
- $A_i^{\text{init}}, A_i^{\text{fin}}$ — envelope amplitudes at the
  start and end of the run (in metres).
- The $\max(0, \cdot)$ wrapper means that **damping** (the beam getting
  closer to the orbit, $A_i^{\text{fin}} < A_i^{\text{init}}$) is *not*
  penalised — only growth.

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

The optimizer searches the parameter space by Monte-Carlo proposals
with greedy acceptance.

Notation:

- $\boldsymbol{\theta} \in \mathbb{R}^M$ — vector of all mutable
  parameters (magnet currents and per-instance calibration
  coefficients), of length $M$.
- $G$ — number of **groups** declared in `parameter_groups_per_position`.
  Each group $g$ may cover several parameter slots; mutating group $g$
  changes all its slots at once.
- $\mathcal{L}(\boldsymbol{\theta})$ — composite loss from § 7.
- $\boldsymbol{\theta}_{\text{best}}$ — current best parameter vector;
  $\mathcal{L}_{\text{best}} = \mathcal{L}(\boldsymbol{\theta}_{\text{best}})$.

**Regular iteration.** A single group $g$ is picked uniformly at random;
its parameters are shifted by a per-group step with a random sign:

$$
\theta'_g = \theta_g + \varepsilon_g \cdot \mathrm{sgn}(u - \tfrac{1}{2}), \qquad
\theta'_k = \theta_k \;\text{ for all } k \ne g
\tag{8.1}
$$

- $u \sim \mathrm{Uniform}(0, 1)$ — a fresh random number.
- $\mathrm{sgn}(u - \tfrac{1}{2}) \in \{-1, +1\}$ — the sign of the
  proposed shift, equally likely positive or negative.
- $\varepsilon_g$ — step size for group $g$, defined per column type:

$$
\varepsilon_g =
\begin{cases}
   \texttt{step\_abs}_g           & \text{if group $g$ targets the } \texttt{arg\_val} \text{ column} \\
   \texttt{step\_percent}_g \cdot \theta_g & \text{if group $g$ targets } \texttt{correct\_coeff}
\end{cases}
\tag{8.2}
$$

- $\texttt{step\_abs}_g$ — absolute step in amperes (typical 0.1–1.0).
- $\texttt{step\_percent}_g$ — fractional step (typical 0.05 ≡ 5 %),
  proportional to the current value so the relative perturbation stays
  meaningful as $\theta_g$ shrinks.

**Kick iteration.** Every `kick_every` iterations *all* groups are
perturbed simultaneously to escape local minima:

$$
\theta'_g = \theta_g + \delta_g, \quad
\delta_g \sim \mathrm{Uniform}(-\Delta,\,+\Delta) \;\;\text{for all } g
\tag{8.3}
$$

- $\Delta = \texttt{kick\_scale}$ — kick amplitude (amperes or fraction,
  depending on the group's column type).
- The independent uniform draws mean each group moves randomly and
  independently; the joint proposal is large enough to jump out of a
  basin that step-by-step search would never leave.

**Greedy acceptance.** Whether the iteration was regular or a kick, the
proposal is accepted only if it lowers the best-known loss:

$$
\boldsymbol{\theta},\;\mathcal{L}_{\text{best}} \;\leftarrow\;
\begin{cases}
   \boldsymbol{\theta}',\;\mathcal{L}(\boldsymbol{\theta}') & \text{if } \mathcal{L}(\boldsymbol{\theta}') < \mathcal{L}_{\text{best}} \\
   \boldsymbol{\theta},\;\mathcal{L}_{\text{best}} & \text{otherwise (revert)}
\end{cases}
\tag{8.4}
$$

This is a pure hill-climber — no temperature, no Metropolis acceptance
of worse states. The kicks are what give it a chance to leave local
minima.

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

**Chord length of a bend.** A bending magnet has an arc length $L$ and
bend angle $\alpha$, so its radius of curvature is $\rho = L/\alpha$.
The straight-line **chord** connecting the entry and exit faces of the
arc has length

$$
L_{\text{chord}} = 2\,\rho\,\sin\!\frac{\alpha}{2} = \frac{2 L}{\alpha}\,\sin\!\frac{\alpha}{2}
\tag{9.1}
$$

- $L$ — arc length of the bending magnet, in metres (from the LTE
  `L=` parameter).
- $\alpha$ — total bend angle, in radians (from the LTE `ANGLE=`
  parameter). Positive $\alpha$ means a left turn, negative a right
  turn.
- $\rho = L/\alpha$ — radius of curvature, in metres.
- $L_{\text{chord}}$ — chord length used to straighten the bend (slightly
  shorter than $L$ when $\alpha \ne 0$). For $\alpha \to 0$ the chord
  tends to $L$ as expected.

The bend is replaced by 5 thin slices uniformly distributed along the
chord at parameter values $\tau \in \{0,\,0.25,\,0.5,\,0.75,\,1.0\}$,
each oriented along the chord direction.

**Ring traversal.** The ring is laid out segment by segment. Each
LTE element of length $\ell$ contributes an in-plane turn $\Delta\psi$
(0 for a drift, $\alpha$ for a bend) and a displacement of its
geometric midpoint. The state update is

$$
\psi' = \psi + \Delta\psi
\tag{9.2}
$$
$$
x' = x + \ell\,\cos\!\Bigl(\psi + \tfrac{\Delta\psi}{2}\Bigr), \qquad
y' = y + \ell\,\sin\!\Bigl(\psi + \tfrac{\Delta\psi}{2}\Bigr)
\tag{9.3}
$$

- $\psi$ — yaw angle (rotation about the global $z$ axis), in radians.
- $\Delta\psi$ — yaw turn produced by this element, in radians (equal
  to $\alpha$ for a bending magnet, 0 for a straight element).
- $\ell$ — element length **along the chord** (so for a bend, $\ell$ is
  one-fifth of $L_{\text{chord}}$ per slice).
- $x, y$ — accumulated the global coordinate system of the current placement
  cursor, in metres.

Using the mid-segment yaw $\psi + \Delta\psi/2$ in the trigonometric
arguments — instead of the entry yaw $\psi$ or the exit yaw
$\psi + \Delta\psi$ — keeps the resulting polyline tangent to the true
arc at each end and balances the curvature error symmetrically across
the slice.

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

No arithmetic — the operation is a **relational join** between
`lattice_layout` (geometry) and `sp_elements` (currents and
calibrations), followed by a union with shortened copies of the master
lattice.

In set-builder notation:

$$
\mathrm{Output}
   = R_{\text{master}}
   \,\cup\,
   \bigcup_{s \,\in\, \mathrm{SubLattices}}\,R_s
\tag{10.1}
$$

where

$$
\begin{aligned}
R_{\text{master}} = \{\;
   (L_0,\;\ell.\text{global\_id},\;s.\text{device\_type\_id},\;\ell.\vec{q},\;s.\text{arg\_val},\;s.\text{correct\_coeff})
   \;\;\big|\;\;
   & \ell \in \mathrm{lattice\_layout}, \\
   & s \in \mathrm{sp\_elements}, \\
   & (s.\text{sp\_type},\,s.\text{sp\_element\_id}) = (\ell.\text{sp\_type},\,\ell.\text{sp\_element\_id})
\;\}
\end{aligned}
\tag{10.2}
$$

$$
R_s = \{\,(s.\text{lattice\_id},\;\text{rest of the row}) \;\big|\; r \in R_{\text{master}},\; r.\text{instance\_id} \le s.\text{max\_gid}\,\}
\tag{10.3}
$$

- $L_0 = 101$ — master lattice ID (from `lte_mapping.json`).
- $\ell$ — a row of `lattice_layout.csv`.
- $\ell.\vec{q} = (\ell.x, \ell.y, \ell.z, \ell.\text{yaw},
   \ell.\text{pitch}, \ell.\text{roll})$ — geometry tuple of the
   element.
- $s$ — in the first block, a row of `sp_elements.csv`. In the second
  block, a sub-lattice descriptor from `sub_lattices.json`.
- The join condition is on the composite key
  $(\text{sp\_type},\,\text{sp\_element\_id})$ — every layout row pulls
  its `device_type_id`, current and coefficient from the matching SP
  template.
- The sub-lattice union creates additional rows with a new
  `lattice_id` but only for instances with `instance_id ≤ max_gid` —
  i.e. a prefix of the master ring.

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

**Envelope amplitude.** For each particle $i$ and each saved trajectory
sample $k$:

$$
A_i(t_k) = \mathrm{dist}\bigl(\vec{x}_i(t_k),\,\Pi_{\text{design}}\bigr)
\tag{11.1}
$$

- $\vec{x}_i(t_k)$ — position of particle $i$ at the $k$-th saved time
  step, in metres.
- $\Pi_{\text{design}}$ — design orbit polyline (centres of lattice
  elements + prelude/postlude tangent extrapolation).
- $A_i(t_k)$ — instantaneous Euclidean distance from the trajectory to
  the design orbit, in metres. Computed by the algorithm of § 6.

**Linear fit** (a baseline drift model):

$$
A(t) \approx m\,t + A_0, \qquad
(m, A_0) = \arg\min_{m, A_0}\,\sum_{k}\bigl(A(t_k) - m\,t_k - A_0\bigr)^2
\tag{11.2}
$$

- $m$ — linear growth rate of the envelope, in metres per second.
  Positive $m$ means the beam is drifting away from the orbit.
- $A_0$ — extrapolated envelope at $t = 0$, in metres.
- The minimisation is the standard least-squares fit (closed form via
  `numpy.polyfit(..., deg=1)`).

**Exponential fit** (a growth/decay model):

$$
A(t) \approx A_0\,e^{\gamma\,t} \quad \Longleftrightarrow \quad
\ln A(t) \approx \ln A_0 + \gamma\,t
\tag{11.3}
$$

- $\gamma$ — exponential growth rate, in inverse seconds (sometimes
  called the **Lyapunov-like exponent** of the envelope). $\gamma > 0$
  means the amplitude grows; $\gamma < 0$ means it decays.
- The fit is linear in $\ln A$, so it is performed in log-space by
  least squares on the points where $A > A_{\text{floor}}$
  (to avoid $\log 0$).

The quality of the exponential fit is judged by the **coefficient of
determination**:

$$
R^2 = 1 - \frac{SS_{\text{res}}}{SS_{\text{tot}}}, \quad
SS_{\text{res}} = \sum_{k}\bigl(\ln A(t_k) - \widehat{\ln A(t_k)}\bigr)^2,
\;\;
SS_{\text{tot}} = \sum_{k}\bigl(\ln A(t_k) - \overline{\ln A}\bigr)^2
\tag{11.4}
$$

- $SS_{\text{res}}$ — residual sum of squares (fit error).
- $SS_{\text{tot}}$ — total sum of squares (variance of the data).
- $\overline{\ln A}$ — mean log-amplitude.
- $R^2 \in (-\infty,\,1]$ — closer to 1 means a better fit; $R^2 < 0$
  means the fit is worse than just using the mean.

A low $R^2$ disqualifies the estimate of $\gamma$ even if it is
numerically large.

**Betatron tune via FFT.** The discrete Fourier transform of the
envelope is computed at the uniform sample step
$\Delta t_{\text{sample}} = t_k - t_{k-1}$:

$$
\widehat{A}(f_n) = \sum_{k = 0}^{K-1} A(t_k)\,e^{-2\pi i\,f_n\,t_k}, \quad
f_n = \frac{n}{K\,\Delta t_{\text{sample}}}, \quad n = 0, 1, \ldots, K-1
\tag{11.5}
$$

- $K$ — number of envelope samples in the trajectory.
- $f_n$ — discrete frequency bin, in hertz.
- $\widehat{A}(f_n)$ — complex Fourier coefficient at frequency $f_n$.

The dominant **betatron frequency** is the bin with the largest
amplitude, excluding the DC bin which carries only the mean:

$$
f_{\text{betatron}} = \arg\max_{f_n > f_{\min}}\bigl|\widehat{A}(f_n)\bigr|
\tag{11.6}
$$

- $f_{\min}$ — small floor frequency that masks out the DC component
  (and optionally a few neighbouring bins to avoid leakage from a
  non-zero mean).
- $f_{\text{betatron}}$ — frequency of the strongest spectral peak, in
  hertz. Converting to a fractional **tune** requires dividing by the
  revolution frequency of the ring; the analyser does so when the
  lattice length is known.

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
