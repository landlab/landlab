from landlab import RasterModelGrid

<<<<<<< HEAD
import numpy as np

=======
from landlab.components.fire_spread import FireSpread

>>>>>>> 302f96ae87f3189411f0a971673f53a00efccd95

def test_grass_fire_propagates():
    import numpy as np

    from landlab import RasterModelGrid
    from landlab.components.fire_spread.fire_spread import FireSpread

    grid = RasterModelGrid((5, 5), xy_spacing=30.0)
    grid.add_zeros("fuel__model", at="cell", dtype=int)[:] = 1      # short grass
    grid.add_zeros("fuel__moisture", at="cell")[:] = 0.08

    fs = FireSpread(grid, ignition_row=2, ignition_col=2)
    fs.run_one_step(dt=60)
<<<<<<< HEAD

    arrival = grid.at_cell["fire__arrival_time"]
    print(arrival[arrival >= 0])          # debug
    assert arrival.max() > 0
    assert np.sum(arrival >= 0) > 1       # fire has spread to neighbours

test_grass_fire_propagates()
=======
    assert grid.at_cell["fire__arrival_time"].min() == 0.0
    assert grid.at_cell["fire__arrival_time"].max() > 0  # spread happened
>>>>>>> 302f96ae87f3189411f0a971673f53a00efccd95
