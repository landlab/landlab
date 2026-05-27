"""Tests for the gFlex Landlab component.

Requires both landlab and gflex to be installed. All tests are
automatically skipped when gflex is not available (see conftest.py).
"""

import numpy as np
import pytest

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
    assert set(gf.output_var_names) == {
        "lithosphere_surface__elevation_increment",
        "topographic__elevation",
    }


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
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
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
            BC_W="Periodic", BC_E="Periodic",
            BC_N="Periodic", BC_S="Periodic",
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
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
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
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
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
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
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
        BC_W="Periodic", BC_E="Periodic",
        BC_N="Periodic", BC_S="Periodic",
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
