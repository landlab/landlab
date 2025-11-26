from landlab import RasterModelGrid
import numpy as np
from landlab.components.fire_spread import FireSpread


def test_grass_fire_propagates():
    from landlab import RasterModelGrid
    import numpy as np
    from landlab.components.fire_spread.fire_spread import FireSpread

    grid = RasterModelGrid((5, 5), xy_spacing=30.0)
    grid.add_zeros("fuel__model", at="cell", dtype=int)[:] = 1      # short grass
    grid.add_zeros("fuel__moisture", at="cell")[:] = 0.08

    fs = FireSpread(grid, ignition_row=2, ignition_col=2)
    fs.run_one_step(dt=60)

    arrival = grid.at_cell["fire__arrival_time"]
    print(arrival[arrival >= 0])          # debug
    assert arrival.max() > 0
    assert np.sum(arrival >= 0) > 1       # fire has spread to neighbours

