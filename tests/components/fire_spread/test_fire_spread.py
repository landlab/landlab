import numpy as np

from landlab.components import FireSpread


def test_grass_fire_propagates(raster_grid_small):
    grid = raster_grid_small
    grid.add_zeros("fuel__model", at="cell", dtype=int)[:] = 1
    grid.add_zeros("fuel__moisture", at="cell")[:] = 0.08
    fs = FireSpread(grid, ignition_row=5, ignition_col=5)
    fs.run_one_step(dt=60)
    assert (
        grid.at_cell["fire__arrival_time"][grid.grid_coords_to_node_id(5, 5, "cell")]
        == 0.0
    )
    assert np.any(grid.at_cell["fire__arrival_time"] > 0)
