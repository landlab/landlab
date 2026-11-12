"""Test MultiClassLinearDiffusion"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal

from landlab import FieldError
from landlab import RasterModelGrid
from landlab.components import MultiClassLinearDiffusion

POROSITY = 0.4
RHO_SED = 2650.0


def _make_grid(n_classes=2, closed=True, peak=True):
    """3x5 grid with three core nodes; node 7 (the middle core node) is a mound."""
    grid = RasterModelGrid((3, 5), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(closed, closed, closed, closed)
    bedrock = grid.add_zeros("bedrock__elevation", at="node")

    if n_classes > 1:
        grains = np.zeros((grid.number_of_nodes, 2))
        grains[grid.core_nodes] = [500.0, 1500.0]
        if peak:
            grains[7] = [3000.0, 1000.0]
        total = grains.sum(axis=1)
    else:
        grains = np.zeros(grid.number_of_nodes)
        grains[grid.core_nodes] = 2000.0
        if peak:
            grains[7] = 4000.0
        total = grains

    grid.add_field("grains__mass", grains, at="node")
    soil = grid.add_field(
        "soil__depth", total / ((1.0 - POROSITY) * RHO_SED), at="node"
    )
    grid.add_field("topographic__elevation", bedrock + soil, at="node")
    return grid


def _diffuser(grid, kd=0.1):
    return MultiClassLinearDiffusion(
        grid, linear_diffusivity=kd, phi=POROSITY, rho_sed=RHO_SED
    )


def test_raise_kwargs_error():
    grid = _make_grid()
    with pytest.raises(TypeError):
        MultiClassLinearDiffusion(grid, bad_var=1)


def test_missing_grains_mass_field():
    grid = RasterModelGrid((3, 5))
    grid.add_zeros("bedrock__elevation", at="node")
    grid.add_zeros("soil__depth", at="node")
    grid.add_zeros("topographic__elevation", at="node")
    with pytest.raises(FieldError):
        MultiClassLinearDiffusion(grid)


def test_parameters_stored():
    grid = _make_grid()
    diffuser = MultiClassLinearDiffusion(
        grid, linear_diffusivity=0.2, phi=0.3, rho_sed=2000.0
    )
    assert diffuser._kd == 0.2
    assert diffuser._phi == 0.3
    assert diffuser._rho_sed == 2000.0
    assert diffuser._n_classes == 2


@pytest.mark.parametrize("n_classes", [1, 2])
def test_number_of_classes_detected(n_classes):
    grid = _make_grid(n_classes=n_classes)
    assert _diffuser(grid)._n_classes == n_classes


def test_docstring_example_values():
    grid = _make_grid()
    _diffuser(grid).run_one_step(dt=500.0)
    assert_array_almost_equal(
        grid.at_node["soil__depth"][grid.core_nodes], [1.613, 1.805, 1.613], decimal=3
    )


@pytest.mark.parametrize("n_classes", [1, 2])
def test_mass_conserved_closed_boundaries(n_classes):
    grid = _make_grid(n_classes=n_classes)
    mass_before = grid.at_node["grains__mass"].sum(axis=0).copy()
    _diffuser(grid).run_one_step(dt=500.0)
    assert_array_almost_equal(grid.at_node["grains__mass"].sum(axis=0), mass_before)


def test_mass_moves_down_slope_by_class():
    grid = _make_grid()
    before = grid.at_node["grains__mass"].copy()
    _diffuser(grid).run_one_step(dt=500.0)
    after = grid.at_node["grains__mass"]
    # Node 7 is rich in class 0, so it loses class 0
    # and its neighbours gain it.
    assert after[7, 0] < before[7, 0]
    assert np.all(after[[6, 8], 0] > before[[6, 8], 0])
    # Symmetric setup gives symmetric result.
    assert_array_almost_equal(after[6], after[8])


@pytest.mark.parametrize("n_classes", [1, 2])
def test_soil_depth_consistent_with_mass(n_classes):
    grid = _make_grid(n_classes=n_classes)
    _diffuser(grid).run_one_step(dt=500.0)
    mass = grid.at_node["grains__mass"].reshape(grid.number_of_nodes, -1).sum(axis=1)
    assert_array_almost_equal(
        grid.at_node["soil__depth"], mass / ((1.0 - POROSITY) * RHO_SED)
    )


def test_topography_is_bedrock_plus_soil():
    grid = _make_grid()
    grid.at_node["bedrock__elevation"][:] += 5.0
    grid.at_node["topographic__elevation"][:] += 5.0
    _diffuser(grid).run_one_step(dt=500.0)
    assert_array_almost_equal(
        grid.at_node["topographic__elevation"],
        grid.at_node["bedrock__elevation"] + grid.at_node["soil__depth"],
    )


def test_flat_surface_unchanged():
    grid = _make_grid(peak=False)
    mass_before = grid.at_node["grains__mass"].copy()
    soil_before = grid.at_node["soil__depth"].copy()
    _diffuser(grid).run_one_step(dt=500.0)
    assert_array_almost_equal(grid.at_node["grains__mass"], mass_before)
    assert_array_almost_equal(grid.at_node["soil__depth"], soil_before)


def test_no_negative_mass_or_soil():
    grid = _make_grid()
    _diffuser(grid, kd=10.0).run_one_step(dt=1.0e4)
    assert np.all(grid.at_node["grains__mass"] >= 0.0)
    assert np.all(grid.at_node["soil__depth"] >= 0.0)


def test_open_boundaries_lose_mass():
    grid = _make_grid(closed=False)
    mass_before = grid.at_node["grains__mass"][grid.core_nodes].sum()
    _diffuser(grid).run_one_step(dt=500.0)
    assert grid.at_node["grains__mass"][grid.core_nodes].sum() < mass_before


def test_substepping_independent_of_dt_split():
    grid_one = _make_grid()
    _diffuser(grid_one).run_one_step(dt=500.0)

    grid_many = _make_grid()
    diffuser = _diffuser(grid_many)
    for _ in range(5):
        diffuser.run_one_step(dt=100.0)

    assert_array_almost_equal(
        grid_one.at_node["soil__depth"], grid_many.at_node["soil__depth"], decimal=2
    )
