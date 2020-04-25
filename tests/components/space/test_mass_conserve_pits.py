from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components import Space, FlowAccumulator


def test_mass_conserve_all_closed_no_depression_finder():

    grid = RasterModelGrid((10,10), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    z = grid.add_zeros("node", "topographic__elevation")
    H = grid.add_zeros("node", "soil__depth")
    z+=grid.x_of_node.copy() + grid.y_of_node.copy()
    z[25]-=40
    z[35]-=40
    z[26]-=40
    z[36]-=40
    z[24]-=40
    z[34]-=40

    z_init = z.copy()

    fa = FlowAccumulator(grid)
    fa.run_one_step()

    sp = Space(grid, phi=0.)
    sp.run_one_step(1.0)

    dz = z - z_init

    assert_array_almost_equal(dz.mean(), 0.)


# Note that we can't make an equivalent test for with a depression finder yet
# because the depression finder can't handle no outlet on the grid.
