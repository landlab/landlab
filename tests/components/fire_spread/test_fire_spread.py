from landlab import RasterModelGrid
from landlab.components.fire_spread import FireSpread


def test_grass_fire_propagates():
    grid = RasterModelGrid((5, 5), xy_spacing=30.0)
    grid.add_zeros("fuel__model", at="cell", dtype=int)[:] = 1
    grid.add_zeros("fuel__moisture", at="cell")[:] = 0.08

    # Ignite the very centre (safe indices)
    fs = FireSpread(grid, ignition_row=2, ignition_col=2)

    fs.run_one_step(dt=60)
    assert grid.at_cell["fire__arrival_time"].min() == 0.0
    assert grid.at_cell["fire__arrival_time"].max() > 0  # spread happened
