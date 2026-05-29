# Tests for the gFlex Landlab component

These tests verify the `gFlex` Landlab component
(`src/landlab/components/gflex/flexure.py`).  They require both
**Landlab** and **gFlex** to be installed; all tests are automatically
skipped when `gflex` is not importable (enforced via `conftest.py`).

Run from the repository root:

```
pytest tests/components/gflex/
```

---

## Test groups

### Component metadata
Checks that `gf.name`, `input_var_names`, `output_var_names`, and the
`_BC_OPTIONS` class attribute all have the expected values.

### Input validation
- Passing a non-`RasterModelGrid` raises `TypeError`.
- An unrecognised BC string raises `ValueError`.
- A `Periodic` BC on one side without the opposite side also being
  `Periodic` raises `ValueError`.

### Physical correctness

**Zero load → zero deflection**
With no applied load the solver must return exactly zero everywhere.

**Uniform load → isostatic equilibrium**
For a spatially uniform load *q* (Pa) with all-periodic boundary
conditions the biharmonic term D∇⁴w vanishes for constant w, so the
governing equation reduces to

    (ρ_m − ρ_fill) · g · w = −q

giving

    w = −q / ((ρ_m − ρ_fill) · g)

This result is independent of elastic thickness and domain size.  The
FD solver must match it to within 0.1 %.

**Sign convention**
A positive (downward) surface load must produce a negative (downward)
deflection.

**Linearity**
Doubling the load must double the deflection.

**Point load / Kelvin-function benchmark**
The deflection due to a concentrated vertical load *P* (N) at the origin
of an infinite 2-D elastic plate is (Turcotte & Schubert, 2002):

    w(r) = P α² / (2π D) · kei(r / α)

where

    D = E Te³ / (12(1 − ν²))        flexural rigidity [N·m]
    α = (D / (Δρ g))^(1/4)          2-D flexural parameter [m]
    kei(x) = Im[ e^{−iπ/4} K₀(x e^{iπ/4}) ]   Kelvin function

Because kei(x) ≤ 0 for x ≥ 0, a positive (downward) load produces a
negative w, consistent with the component's sign convention for
`lithosphere_surface__elevation_increment`.

The test uses a 100 × 100 grid at dx = 5 km (Te = 10 km, α ≈ 21 km,
domain ≈ 24 α wide) with `0Moment0Shear` boundary conditions.
Comparison points are sampled in the interior at radii 1.5–3.5 α from
the load centre — well away from the kei zero crossing (forebulge onset
at r/α ≈ 3.91) where the relative error metric is undefined, and at
least 3 α from every boundary to minimise boundary reflections.
Tolerance is 5 %.

**Grid-convergence study** (run separately, 2026-05-27):
The FD scheme is second-order accurate (O(dx²)).  At dx = 1250 m
(dx/α ≈ 0.06) relative errors are 0.02–0.1 %.  The table below shows
convergence orders measured between successive grid refinements
(dx = 10 000, 5 000, 2 500, 1 250 m):

| r/α  | convergence orders (coarse → fine) |
|------|------------------------------------|
| 0.97 | 2.20 → 2.07 → 2.01               |
| 1.46 | 2.68 → 2.18 → 2.04               |
| 1.95 | 1.07 → 1.88 → 1.97               |
| 2.92 | −0.07 → 1.77 → 1.95              |

The −0.07 at r = 2.92 α, dx = 10 km is a near-cancellation artefact
(the error changes sign between those two resolutions); it is not a real
anomaly.

### `topographic__elevation` interaction
- The component runs without a `topographic__elevation` field.
- When the field is present it is updated after each step.
- The field accumulates the change in deflection between steps, not the
  absolute deflection; so applying the same load a second time leaves
  `topographic__elevation` unchanged (the increment is zero).

### Repeated calls and load changes
- Calling `run_one_step` twice with the same load returns the same
  `lithosphere_surface__elevation_increment` both times.
- Changing the load between calls changes the output.

### Mirror boundary condition
`Mirror` is a valid gFlex BC (symmetry across the edge) that was previously
missing from the Landlab wrapper's `_BC_OPTIONS`.  One test verifies that
all four edges set to `Mirror` are accepted and produce a downward deflection
under a uniform positive load.

### Variable elastic thickness
The component accepts `elastic_thickness` as a scalar float, as the
string name of an existing node field, or as a NumPy array of shape
`grid.shape`.  All three forms must run without error.
