import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import ErosionDeposition
from landlab.components import FlowAccumulator
from landlab.components import Space


@pytest.fixture
def grid():
    grid = RasterModelGrid((10, 10), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    z = grid.add_zeros("node", "topographic__elevation")
    grid.add_zeros("node", "soil__depth")
    z += grid.x_of_node.copy() + grid.y_of_node.copy()
    z[25] -= 40
    z[35] -= 40
    z[26] -= 40
    z[36] -= 40
    z[24] -= 40
    z[34] -= 40

    return grid


# consider full combinitorics of solver, two phi, ED and Space, and (if space)
# initial soil depth of very large and zero.
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
@pytest.mark.parametrize("v_s", [1.5])
@pytest.mark.parametrize("dt", [2])
def test_mass_conserve_all_closed_ErosionDeposition(grid, solver, v_s, dt):
    z_init = grid.at_node["topographic__elevation"].copy()

    fa = FlowAccumulator(grid)
    fa.run_one_step()

    ed = ErosionDeposition(grid, solver=solver, v_s=v_s)
    ed.run_one_step(dt)

    dz = z_init - grid.at_node["topographic__elevation"]

    # For Erosion Deposition, porosity should not have any effect, because
    # the component operates in terms of bulk-equivalent sediment flux,
    # erosion, and deposition.

    assert_array_almost_equal(dz.sum(), 0.0, decimal=10)


@pytest.mark.parametrize("phi", [0.0, 0.3])
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
@pytest.mark.parametrize("H", [0, 1, 1000])
@pytest.mark.parametrize("v_s", [1.5])
@pytest.mark.parametrize("H_star", [0.1])
@pytest.mark.parametrize("dt", [2])
def test_mass_conserve_all_closed_Space(grid, H, solver, phi, v_s, H_star, dt):
    grid.at_node["soil__depth"][:] = H

    z_init = grid.at_node["topographic__elevation"].copy()

    fa = FlowAccumulator(grid)
    fa.run_one_step()

    ed = Space(grid, solver=solver, phi=phi, v_s=v_s, H_star=H_star)
    ed.run_one_step(dt)

    # in space, everything is either bedrock or sediment. check for
    # conservation.
    dH = grid.at_node["soil__depth"][:] - H

    # sediment is defined as having a porosity so all changes (up or down )
    # must be adjusted to mass.
    dH *= 1 - phi

    dBr = grid.at_node["bedrock__elevation"] - (z_init - H)
    mass_change = dH + dBr

    assert_array_almost_equal(mass_change.sum(), 0.0, decimal=10)


# Note that we can't make an equivalent test for with a depression finder yet
# because the depression finder can't handle no outlet on the grid.
# but what we can do is make an example in which there is a big sink in which
# almost all sediment is trapped. We can then assert that all sediment is
# either trapped OR that it is sent out of the one outlet node.
@pytest.fixture()
def grid2(grid):
    grid.status_at_node[1] = grid.BC_NODE_IS_FIXED_VALUE
    return grid


# consider full combinitorics of solver, two phi, depression finding or not,
#  ED and Space, and (if space) initial soil depth of very large and zero.


@pytest.mark.parametrize("depression_finder", [None, "DepressionFinderAndRouter"])
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
@pytest.mark.parametrize("v_s", [1.5])
@pytest.mark.parametrize("dt", [2])
def test_mass_conserve_with_depression_finder_ErosionDeposition(
    grid2, solver, depression_finder, v_s, dt
):
    assert grid2.status_at_node[1] == grid2.BC_NODE_IS_FIXED_VALUE

    z_init = grid2.at_node["topographic__elevation"].copy()

    if depression_finder is None:
        fa = FlowAccumulator(grid2)
    else:
        fa = FlowAccumulator(grid2, depression_finder=depression_finder, routing="D4")
    fa.run_one_step()

    ed = ErosionDeposition(grid2, solver=solver, v_s=v_s)
    ed.run_one_step(dt)

    dz = grid2.at_node["topographic__elevation"] - z_init

    # assert that the mass loss over the surface is exported through the one
    # outlet.
    net_change = dz[grid2.core_nodes].sum() + (
        ed.sediment_influx[1] * dt / grid2.cell_area_at_node[11]
    )
    assert_array_almost_equal(net_change, 0.0, decimal=10)


@pytest.mark.parametrize("depression_finder", [None, "DepressionFinderAndRouter"])
@pytest.mark.parametrize("phi", [0.0, 0.3])
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
@pytest.mark.parametrize("H", [0, 1000])
@pytest.mark.parametrize("v_s", [1.5])
@pytest.mark.parametrize("H_star", [0.1])
@pytest.mark.parametrize("dt", [2])
def test_mass_conserve_with_depression_finder_Space(
    grid2, H, solver, depression_finder, phi, v_s, H_star, dt
):
    grid2.at_node["soil__depth"][:] = H
    assert grid2.status_at_node[1] == grid2.BC_NODE_IS_FIXED_VALUE

    z_init = grid2.at_node["topographic__elevation"].copy()

    if depression_finder is None:
        fa = FlowAccumulator(grid2)
    else:
        fa = FlowAccumulator(grid2, depression_finder=depression_finder, routing="D4")
    fa.run_one_step()

    ed = Space(grid2, solver=solver, phi=phi, v_s=v_s, H_star=H_star)
    ed.run_one_step(dt)

    # see above test for notes.
    dH = grid2.at_node["soil__depth"][:] - H
    dH *= 1 - phi
    dBr = grid2.at_node["bedrock__elevation"] - (z_init - H)
    mass_change = dH + dBr

    # assert that the mass loss over the surface is exported through the one
    # outlet.
    net_change = mass_change[grid2.core_nodes].sum() + (
        ed.sediment_influx[1] * dt / grid2.cell_area_at_node[11]
    )
    assert_array_almost_equal(net_change, 0.0, decimal=10)
