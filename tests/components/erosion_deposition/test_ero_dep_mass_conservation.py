import pytest
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import ErosionDeposition, FlowAccumulator, Space


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


@pytest.mark.parametrize(
    "Component_SoilThickness", [(ErosionDeposition, 0), (Space, 0), (Space, 100)]
)
@pytest.mark.parametrize("phi", [0.0, 0.3])
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
def test_mass_conserve_all_closed(grid, Component_SoilThickness, solver, phi):
    Component, H = Component_SoilThickness
    grid.at_node["soil__depth"][:] = H

    z_init = grid.at_node["topographic__elevation"].copy()

    fa = FlowAccumulator(grid)
    fa.run_one_step()

    ed = Component(grid, solver=solver, phi=phi, v_s=1.5)
    ed.run_one_step(2)

    dz = grid.at_node["topographic__elevation"] - z_init

    # unpoof by phi where deposition occured so we can compare mass. We can do
    # this because only one timestep. (I think, but not sure, even with adaptive.)
    where_depo = dz > 0
    mass_change = dz.copy()
    mass_change[where_depo > 0] = dz[where_depo] * (1 - phi)
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


@pytest.mark.parametrize(
    "Component_SoilThickness", [(ErosionDeposition, 0), (Space, 0), (Space, 100)]
)
@pytest.mark.parametrize("phi", [0.0, 0.3])
@pytest.mark.parametrize("solver", ["basic", "adaptive"])
@pytest.mark.parametrize("depression_finder", [None, "DepressionFinderAndRouter"])
def test_mass_conserve_with_depression_finder(
    grid2, Component_SoilThickness, solver, depression_finder, phi
):
    Component, H = Component_SoilThickness
    dt = 2
    grid2.at_node["soil__depth"][:] = H
    assert grid2.status_at_node[1] == grid2.BC_NODE_IS_FIXED_VALUE

    z_init = grid2.at_node["topographic__elevation"].copy()

    if depression_finder is None:
        fa = FlowAccumulator(grid2)
    else:
        fa = FlowAccumulator(grid2, depression_finder=depression_finder, routing="D4")
    fa.run_one_step()

    ed = Component(grid2, solver=solver, phi=phi, v_s=1.5)
    ed.run_one_step(dt)

    dz = grid2.at_node["topographic__elevation"] - z_init

    # unpoof by phi where deposition occured so we can compare mass. We can do
    # this because only one timestep. (I think, but not sure, even with adaptive.)
    where_depo = dz > 0
    mass_change = dz.copy()
    mass_change[where_depo > 0] = dz[where_depo] * (1 - phi)

    # assert that the mass loss over the surface is exported through the one
    # outlet.
    net_change = mass_change[grid2.core_nodes].sum() + (
        ed._qs_in[1] * dt / grid2.cell_area_at_node[11]
    )

    assert_array_almost_equal(net_change, 0.0, decimal=10)
