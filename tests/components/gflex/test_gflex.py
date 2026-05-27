"""Tests for the gFlex Landlab component.

Requires both landlab and gflex to be installed. All tests are
automatically skipped when gflex is not available (see conftest.py).

The Kelvin-function benchmark (``test_point_load_kelvin_function``) tests
the analytical 2-D point-load solution:

    w(r) = P * alpha² / (2π * D) * kei(r / alpha)

where alpha = (D / (drho * g))^(1/4) is the 2-D flexural parameter.
kei is the Kelvin function (imaginary part of K₀(r * exp(iπ/4)) * exp(-iπ/4)).
kei(x) < 0 for x ≥ 0, so a positive (downward) load produces negative w.

A grid-convergence study (E=65 GPa, Te=10 km, nu=0.25, rho_m=3300, g=9.81)
confirmed second-order (O(dx²)) convergence at all tested radii:

  r/alpha  convergence orders (coarse → fine)
  0.97     2.20 → 2.07 → 2.01
  1.46     2.68 → 2.18 → 2.04
  1.95     1.07 → 1.88 → 1.97
  2.92    -0.07 → 1.77 → 1.95

At dx=1250 m (dx/alpha ≈ 0.06) relative errors are 0.02–0.1 %; the 5 %
tolerance used here is conservative for the dx=5 km grid.
"""

import numpy as np
import pytest
from scipy.special import kei

from landlab import RasterModelGrid
from landlab.components.gflex.flexure import gFlex

# ---------------------------------------------------------------------------
# Component metadata
# ---------------------------------------------------------------------------


def test_component_name(gf):
    assert gf.name == "gFlex"


def test_input_var_names(gf):
    assert gf.input_var_names == ("surface_load__stress",)


def test_output_var_names(gf):
    assert "lithosphere_surface__elevation_increment" in set(gf.output_var_names)


def test_bc_options_class_attribute():
    assert hasattr(gFlex, "_BC_OPTIONS")
    assert "0Displacement0Slope" in gFlex._BC_OPTIONS
    assert "0Moment0Shear" in gFlex._BC_OPTIONS
    assert "0Slope0Shear" in gFlex._BC_OPTIONS
    assert "Periodic" in gFlex._BC_OPTIONS


def test_output_field_created_on_init(gf, grid):
    assert "lithosphere_surface__elevation_increment" in grid.at_node


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_bad_grid_type_raises_typeerror():
    with pytest.raises(TypeError):
        gFlex("not_a_grid")


@pytest.mark.parametrize("bad_bc", ["Free", "ZeroSlope", "open"])
def test_bad_bc_value_raises_valueerror(grid, bad_bc):
    with pytest.raises(ValueError):
        gFlex(grid, BC_W=bad_bc, quiet=True)


def test_unpaired_periodic_west_east_raises(grid):
    with pytest.raises(ValueError):
        gFlex(grid, BC_W="Periodic", BC_E="0Displacement0Slope", quiet=True)


def test_unpaired_periodic_north_south_raises(grid):
    with pytest.raises(ValueError):
        gFlex(grid, BC_N="Periodic", BC_S="0Displacement0Slope", quiet=True)


# ---------------------------------------------------------------------------
# Physical correctness
# ---------------------------------------------------------------------------


def test_zero_load_zero_deflection(grid):
    """No load should produce no deflection."""
    gf = gFlex(grid, quiet=True)
    gf.run_one_step()
    w = grid.at_node["lithosphere_surface__elevation_increment"]
    np.testing.assert_allclose(w, 0.0, atol=1e-10)


def test_uniform_load_isostatic_deflection():
    """Uniform load with all-periodic BCs should yield exact isostasy.

    For a spatially uniform load q [Pa] the biharmonic term D*nabla^4*w
    vanishes, so the FD solution reduces to:

        (rho_m - rho_fill) * g * w = -q
        w = -q / ((rho_m - rho_fill) * g)

    This holds to within floating-point precision for the FD solver.
    """
    rho_m = 3300.0
    rho_fill = 0.0
    g = 9.81
    q = 1e4  # Pa

    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("surface_load__stress", at="node")
    mg.at_node["surface_load__stress"][:] = q

    gf = gFlex(
        mg,
        rho_mantle=rho_m,
        rho_fill=rho_fill,
        g=g,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere_surface__elevation_increment"][mg.core_nodes]
    expected = -q / ((rho_m - rho_fill) * g)
    np.testing.assert_allclose(w, expected, rtol=1e-3)


def test_deflection_negative_for_positive_load():
    """A positive (downward) surface load should lower the surface."""
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("surface_load__stress", at="node")
    mg.at_node["surface_load__stress"][:] = 1e4

    gf = gFlex(
        mg,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere_surface__elevation_increment"][mg.core_nodes]
    assert np.all(w < 0.0)


def test_larger_load_larger_deflection():
    """Doubling the load should double the deflection (linearity)."""

    def run_with_load(q):
        mg = RasterModelGrid((10, 10), xy_spacing=25000.0)
        mg.add_zeros("surface_load__stress", at="node")
        mg.at_node["surface_load__stress"][:] = q
        gf = gFlex(
            mg,
            BC_W="Periodic",
            BC_E="Periodic",
            BC_N="Periodic",
            BC_S="Periodic",
            quiet=True,
        )
        gf.run_one_step()
        return mg.at_node["lithosphere_surface__elevation_increment"].copy()

    w1 = run_with_load(1e4)
    w2 = run_with_load(2e4)
    np.testing.assert_allclose(w2, 2.0 * w1, rtol=1e-6)


# ---------------------------------------------------------------------------
# topographic__elevation interaction
# ---------------------------------------------------------------------------


def test_topo_not_required(grid):
    """Component must run without a topographic__elevation field."""
    gf = gFlex(grid, quiet=True)
    gf.run_one_step()  # should not raise


def test_topo_updated_when_present(grid_with_topo):
    """topographic__elevation should be modified when the field exists."""
    grid_with_topo.at_node["surface_load__stress"][:] = 1e4
    gf = gFlex(
        grid_with_topo,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()
    topo = grid_with_topo.at_node["topographic__elevation"]
    assert not np.all(topo == 0.0)


def test_topo_tracks_cumulative_deflection(grid_with_topo):
    """Topographic elevation should accumulate deflection across calls."""
    grid_with_topo.at_node["surface_load__stress"][:] = 1e4
    gf = gFlex(
        grid_with_topo,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()
    topo_after_1 = grid_with_topo.at_node["topographic__elevation"].copy()
    w_after_1 = grid_with_topo.at_node[
        "lithosphere_surface__elevation_increment"
    ].copy()

    # Same load: topo should shift by the same increment again
    gf.run_one_step()
    topo_after_2 = grid_with_topo.at_node["topographic__elevation"].copy()
    w_after_2 = grid_with_topo.at_node["lithosphere_surface__elevation_increment"]

    np.testing.assert_allclose(w_after_2, w_after_1, rtol=1e-6)
    np.testing.assert_allclose(topo_after_2, topo_after_1, rtol=1e-6)


# ---------------------------------------------------------------------------
# Repeated calls and load changes
# ---------------------------------------------------------------------------


def test_repeated_calls_same_result(grid):
    """Calling run_one_step twice with the same load gives the same deflection."""
    grid.at_node["surface_load__stress"][:] = 1e4
    gf = gFlex(
        grid,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )
    gf.run_one_step()
    w1 = grid.at_node["lithosphere_surface__elevation_increment"].copy()

    gf.run_one_step()
    w2 = grid.at_node["lithosphere_surface__elevation_increment"].copy()

    np.testing.assert_allclose(w1, w2, rtol=1e-10)


def test_load_change_reflected(grid):
    """Updating the load field between calls changes the deflection."""
    gf = gFlex(
        grid,
        BC_W="Periodic",
        BC_E="Periodic",
        BC_N="Periodic",
        BC_S="Periodic",
        quiet=True,
    )

    grid.at_node["surface_load__stress"][:] = 1e4
    gf.run_one_step()
    w1 = grid.at_node["lithosphere_surface__elevation_increment"].copy()

    grid.at_node["surface_load__stress"][:] = 2e4
    gf.run_one_step()
    w2 = grid.at_node["lithosphere_surface__elevation_increment"].copy()

    assert not np.allclose(w1, w2)


# ---------------------------------------------------------------------------
# Variable elastic thickness
# ---------------------------------------------------------------------------


def test_scalar_te_runs(grid):
    gf = gFlex(grid, elastic_thickness=35000.0, quiet=True)
    gf.run_one_step()


def test_array_te_runs(grid):
    Te = np.full(grid.shape, 35000.0)
    gf = gFlex(grid, elastic_thickness=Te, quiet=True)
    gf.run_one_step()


def test_array_te_field_name_runs(grid):
    grid.add_full("elastic_thickness", 35000.0, at="node")
    gf = gFlex(grid, elastic_thickness="elastic_thickness", quiet=True)
    gf.run_one_step()


# ---------------------------------------------------------------------------
# Analytical benchmark — point load / Kelvin-function solution
# ---------------------------------------------------------------------------


def test_point_load_kelvin_function():
    """FD deflection matches the Kelvin-function infinite-plate solution.

    Set-up
    ------
    * 100 × 100 grid, dx = dy = 5 km  →  500 km domain
    * Te = 10 km  →  alpha ≈ 21 km  →  domain ≈ 24 alpha wide
    * Central cell loaded with q = 1e6 Pa
    * 0Moment0Shear BCs on all edges (free edges, minimal reflection)

    Comparison points are at radii 1.5–3.5 alpha from the load centre
    (below the forebulge onset at kei zero ≈ 3.91 alpha) and at least
    3 alpha from every boundary.  Tolerance is 5 %; the FD error at
    dx = 5 km is well below 1 % in this range (see module docstring).
    """
    E = 65e9
    Te = 10000.0
    nu = 0.25
    rho_m = 3300.0
    rho_fill = 0.0
    g = 9.81

    dx = dy = 5000.0
    nrows = ncols = 100

    D = E * Te**3 / (12.0 * (1.0 - nu**2))
    drho = rho_m - rho_fill
    alpha = (D / (drho * g)) ** 0.25

    mg = RasterModelGrid((nrows, ncols), xy_spacing=dx)
    mg.add_zeros("surface_load__stress", at="node")

    ci = nrows // 2
    cj = ncols // 2
    q_load = 1e6
    mg.at_node["surface_load__stress"][mg.grid_coords_to_node_id(ci, cj)] = q_load

    gf = gFlex(
        mg,
        Youngs_modulus=E,
        Poissons_ratio=nu,
        rho_mantle=rho_m,
        rho_fill=rho_fill,
        g=g,
        elastic_thickness=Te,
        BC_W="0Moment0Shear",
        BC_E="0Moment0Shear",
        BC_N="0Moment0Shear",
        BC_S="0Moment0Shear",
        quiet=True,
    )
    gf.run_one_step()

    w_grid = mg.at_node["lithosphere_surface__elevation_increment"].reshape(mg.shape)

    P = q_load * dx * dy

    def w_analytical(r):
        return P * alpha**2 / (2.0 * np.pi * D) * kei(r / alpha)

    min_r = 1.5 * alpha
    max_r = 3.5 * alpha
    min_from_boundary = 3.0 * alpha

    errors = []
    for dj in range(1, ncols - cj):
        r = dj * dx
        col_idx = cj + dj
        dist_from_east = (ncols - 1 - col_idx) * dx
        if r < min_r or r > max_r:
            continue
        if dist_from_east < min_from_boundary:
            continue
        w_fd = w_grid[ci, col_idx]
        w_ana = w_analytical(r)
        errors.append(abs(w_fd - w_ana) / abs(w_ana))

    assert len(errors) > 0, "No comparison points satisfied the distance criteria"
    assert max(errors) < 0.05, (
        f"Largest relative error {max(errors):.1%} exceeds 5 % tolerance; "
        "FD solution does not match the Kelvin-function analytical result."
    )
