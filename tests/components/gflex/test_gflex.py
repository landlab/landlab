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
    assert "load__normal_component_of_stress" in gf.input_var_names


def test_optional_te_in_optional_var_names(gf):
    assert "lithosphere__elastic_thickness" in gf.optional_var_names


def test_output_var_names(gf):
    assert "lithosphere__vertical_displacement" in set(gf.output_var_names)


def test_valid_bc_strings_from_gflex():
    import gflex as _gflex

    expected = {
        "zero_displacement_zero_slope",
        "zero_displacement_zero_moment",
        "zero_moment_zero_shear",
        "zero_slope_zero_shear",
        "no_outside_loads",
        "periodic",
        "clamped",
        "free",
        "mirror",
        "infinite",
    }
    assert expected <= _gflex.VALID_BC_STRINGS_2D


def test_output_field_created_on_init(gf, grid):
    assert "lithosphere__vertical_displacement" in grid.at_node


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_bad_grid_type_raises_typeerror():
    with pytest.raises(TypeError):
        gFlex("not_a_grid")


@pytest.mark.parametrize("bad_bc", ["Free", "ZeroSlope", "open"])
def test_bad_bc_value_raises_valueerror(grid, bad_bc):
    with pytest.raises(ValueError):
        gFlex(grid, bc_west=bad_bc, quiet=True)


def test_unpaired_periodic_west_east_raises(grid):
    with pytest.raises(ValueError):
        gFlex(
            grid,
            bc_west="periodic",
            bc_east="zero_displacement_zero_slope",
            quiet=True,
        )


def test_unpaired_periodic_north_south_raises(grid):
    with pytest.raises(ValueError):
        gFlex(
            grid,
            bc_north="periodic",
            bc_south="zero_displacement_zero_slope",
            quiet=True,
        )


def test_invalid_density_contrast_raises(grid):
    """rho_fill >= rho_mantle gives no restoring force; gFlex raises ValueError."""
    with pytest.raises(ValueError):
        gFlex(grid, rho_mantle=3300.0, rho_fill=3300.0, quiet=True)


def test_nonexistent_te_field_raises(grid):
    """Passing a nonexistent field name for elastic_thickness should raise."""
    with pytest.raises((KeyError, ValueError)):
        gFlex(grid, elastic_thickness="does_not_exist", quiet=True)


def test_method_uppercase_accepted(grid):
    """method='FD' (uppercase) should work via .lower() normalisation."""
    gf = gFlex(grid, method="FD", quiet=True)
    gf.run_one_step()


def test_dict_bc_passes_through(grid):
    """A dict BC passes Landlab string validation and reaches gFlex."""
    grid.at_node["load__normal_component_of_stress"][:] = 1e4
    gf = gFlex(
        grid,
        bc_west="zero_displacement_zero_slope",
        bc_east={"displacement": np.zeros(grid.shape[0])},
        bc_north="zero_displacement_zero_slope",
        bc_south="zero_displacement_zero_slope",
        quiet=True,
    )
    gf.run_one_step()
    w = grid.at_node["lithosphere__vertical_displacement"]
    assert np.any(w < 0.0)


# ---------------------------------------------------------------------------
# Physical correctness
# ---------------------------------------------------------------------------


def test_zero_load_zero_deflection(grid):
    """No load should produce no deflection."""
    gf = gFlex(grid, quiet=True)
    gf.run_one_step()
    w = grid.at_node["lithosphere__vertical_displacement"]
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
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = q

    gf = gFlex(
        mg,
        rho_mantle=rho_m,
        rho_fill=rho_fill,
        g=g,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere__vertical_displacement"][mg.core_nodes]
    expected = -q / ((rho_m - rho_fill) * g)
    np.testing.assert_allclose(w, expected, rtol=1e-3)


def test_deflection_negative_for_positive_load():
    """A positive (downward) surface load should lower the surface."""
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = 1e4

    gf = gFlex(
        mg,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere__vertical_displacement"][mg.core_nodes]
    assert np.all(w < 0.0)


def test_larger_load_larger_deflection():
    """Doubling the load should double the deflection (linearity)."""

    def run_with_load(q):
        mg = RasterModelGrid((10, 10), xy_spacing=25000.0)
        mg.add_zeros("load__normal_component_of_stress", at="node")
        mg.at_node["load__normal_component_of_stress"][:] = q
        gf = gFlex(
            mg,
            bc_west="periodic",
            bc_east="periodic",
            bc_north="periodic",
            bc_south="periodic",
            quiet=True,
        )
        gf.run_one_step()
        return mg.at_node["lithosphere__vertical_displacement"].copy()

    w1 = run_with_load(1e4)
    w2 = run_with_load(2e4)
    np.testing.assert_allclose(w2, 2.0 * w1, rtol=1e-6)


def test_negative_load_upward_deflection():
    """A negative (upward) surface load should raise the surface."""
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = -1e4

    gf = gFlex(
        mg,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere__vertical_displacement"][mg.core_nodes]
    assert np.all(w > 0.0)


def test_nonzero_rho_fill_isostatic_deflection():
    """Uniform load with rho_fill > 0 and periodic BCs matches the exact solution.

    The infill density (e.g. seawater at 1030 kg m⁻³) reduces the effective
    density contrast driving isostasy:

        w = -q / ((rho_m - rho_fill) * g)
    """
    rho_m = 3300.0
    rho_fill = 1030.0  # seawater
    g = 9.81
    q = 1e4

    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = q

    gf = gFlex(
        mg,
        rho_mantle=rho_m,
        rho_fill=rho_fill,
        g=g,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere__vertical_displacement"][mg.core_nodes]
    expected = -q / ((rho_m - rho_fill) * g)
    np.testing.assert_allclose(w, expected, rtol=1e-3)


def test_method_fft_runs():
    """FFT solver produces a downward deflection under a positive load."""
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = 1e4

    gf = gFlex(
        mg,
        method="fft",
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()

    w = mg.at_node["lithosphere__vertical_displacement"][mg.core_nodes]
    assert np.all(w < 0.0)


# ---------------------------------------------------------------------------
# Repeated calls and load changes
# ---------------------------------------------------------------------------


def test_repeated_calls_same_result(grid):
    """Calling run_one_step twice with the same load gives the same deflection."""
    grid.at_node["load__normal_component_of_stress"][:] = 1e4
    gf = gFlex(
        grid,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()
    w1 = grid.at_node["lithosphere__vertical_displacement"].copy()

    gf.run_one_step()
    w2 = grid.at_node["lithosphere__vertical_displacement"].copy()

    np.testing.assert_allclose(w1, w2, rtol=1e-10)


def test_load_change_reflected(grid):
    """Updating the load field between calls changes the deflection."""
    gf = gFlex(
        grid,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )

    grid.at_node["load__normal_component_of_stress"][:] = 1e4
    gf.run_one_step()
    w1 = grid.at_node["lithosphere__vertical_displacement"].copy()

    grid.at_node["load__normal_component_of_stress"][:] = 2e4
    gf.run_one_step()
    w2 = grid.at_node["lithosphere__vertical_displacement"].copy()

    assert not np.allclose(w1, w2)


# ---------------------------------------------------------------------------
# Variable elastic thickness
# ---------------------------------------------------------------------------


def test_mirror_bc_runs(grid):
    """Mirror BC should be accepted and produce a valid deflection."""
    grid.at_node["load__normal_component_of_stress"][:] = 1e4
    gf = gFlex(
        grid,
        bc_west="mirror",
        bc_east="mirror",
        bc_north="mirror",
        bc_south="mirror",
        quiet=True,
    )
    gf.run_one_step()
    w = grid.at_node["lithosphere__vertical_displacement"]
    assert np.all(w < 0.0)


def test_scalar_te_runs(grid):
    gf = gFlex(grid, elastic_thickness=35000.0, quiet=True)
    gf.run_one_step()


def test_array_te_runs(grid):
    Te = np.full(grid.shape, 35000.0)
    gf = gFlex(grid, elastic_thickness=Te, quiet=True)
    gf.run_one_step()


def test_variable_te_array_differs_from_uniform():
    """A spatially varying T_e array gives different deflection than uniform T_e.

    Left half: thin plate (Te=10 km); right half: thick plate (Te=60 km).
    Under a uniform load the thin half should deflect more than the thick half.
    """
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = 1e4

    Te_uniform = np.full(mg.shape, 35000.0)
    Te_variable = np.full(mg.shape, 35000.0)
    Te_variable[:, : mg.shape[1] // 2] = 10000.0  # thin west half
    Te_variable[:, mg.shape[1] // 2 :] = 60000.0  # thick east half

    def run(Te):
        gf = gFlex(
            mg,
            elastic_thickness=Te,
            bc_west="periodic",
            bc_east="periodic",
            bc_north="periodic",
            bc_south="periodic",
            quiet=True,
        )
        gf.run_one_step()
        return mg.at_node["lithosphere__vertical_displacement"].copy()

    w_uniform = run(Te_uniform)
    w_variable = run(Te_variable)

    assert not np.allclose(w_uniform, w_variable)


def test_array_te_field_name_runs(grid):
    grid.add_full("elastic_thickness", 35000.0, at="node")
    gf = gFlex(grid, elastic_thickness="elastic_thickness", quiet=True)
    gf.run_one_step()


def test_te_updated_from_lithosphere_field():
    """lithosphere__elastic_thickness field is read each step (BMI pattern).

    A stiffer plate (larger T_e) deflects less under the same load.
    """
    mg = RasterModelGrid((20, 20), xy_spacing=25000.0)
    mg.add_zeros("load__normal_component_of_stress", at="node")
    mg.at_node["load__normal_component_of_stress"][:] = 1e4
    te_field = mg.add_full(
        "lithosphere__elastic_thickness", 35000.0, at="node", dtype=float
    )

    gf = gFlex(
        mg,
        bc_west="periodic",
        bc_east="periodic",
        bc_north="periodic",
        bc_south="periodic",
        quiet=True,
    )
    gf.run_one_step()
    w_thin = mg.at_node["lithosphere__vertical_displacement"].copy()

    te_field[:] = 70000.0
    gf.run_one_step()
    w_thick = mg.at_node["lithosphere__vertical_displacement"].copy()

    assert np.all(np.abs(w_thick) < np.abs(w_thin))


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
    * no_outside_loads BCs (default; mimics infinite plate)

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
    mg.add_zeros("load__normal_component_of_stress", at="node")

    ci = nrows // 2
    cj = ncols // 2
    q_load = 1e6
    mg.at_node["load__normal_component_of_stress"][
        mg.grid_coords_to_node_id(ci, cj)
    ] = q_load

    gf = gFlex(
        mg,
        Youngs_modulus=E,
        Poissons_ratio=nu,
        rho_mantle=rho_m,
        rho_fill=rho_fill,
        g=g,
        elastic_thickness=Te,
        quiet=True,
    )
    gf.run_one_step()

    w_grid = mg.at_node["lithosphere__vertical_displacement"].reshape(mg.shape)

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
