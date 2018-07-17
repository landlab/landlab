import pytest
import numpy as np

from landlab import RasterModelGrid
from landlab import BAD_INDEX_VALUE as XX
from landlab.components.sink_fill import SinkFiller


@pytest.fixture
def sink_grid1():
    """Create a 7x7 test grid with a well defined hole in it."""
    sink_grid = RasterModelGrid((7, 7), spacing=1.)

    z = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
            [0.0, 2.0, 1.6, 1.5, 1.6, 2.0, 0.0],
            [0.0, 2.0, 1.7, 1.6, 1.7, 2.0, 0.0],
            [0.0, 2.0, 1.8, 2.0, 2.0, 2.0, 0.0],
            [0.0, 1.0, 0.6, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, -0.5, 0.0, 0.0, 0.0, 0.0],
        ]
    ).flatten()

    sink_grid.add_field("topographic__elevation", z, at="node", units="-")

    sink_grid.outlet = 30
    sink_grid.lake_code = 17
    sink_grid.lake = np.array([16, 17, 18, 23, 24, 25])

    return sink_grid


@pytest.fixture
def sink_grid2():
    """
    Create a 10x10 test grid with a well defined hole in it, from a flat
    surface.
    """
    sink_grid = RasterModelGrid((10, 10), spacing=1.)

    lake = np.array([44, 45, 46, 54, 55, 56, 64, 65, 66])

    z = np.ones(100, dtype=float)
    z[lake] = 0.

    sink_grid.add_field("topographic__elevation", z, at="node", units="-")
    sink_grid.lake = lake

    return sink_grid


@pytest.fixture
def sink_grid3():
    """
    Create a 10x10 test grid with two well defined holes in it, into an
    inclined surface.
    """
    sink_grid = RasterModelGrid((10, 10), spacing=1.)

    lake1 = np.array([34, 35, 36, 44, 45, 46, 54, 55, 56])
    lake2 = np.array([77, 78, 87, 88])
    guard_nodes = np.array([23, 33, 53, 63])
    lake = np.concatenate((lake1, lake2))

    z = np.ones(100, dtype=float)
    # add slope
    z += sink_grid.node_x
    z[guard_nodes] += 0.001
    z[lake] = 0.

    sink_grid.add_field("node", "topographic__elevation", z, units="-")
    sink_grid.lake1 = lake1
    sink_grid.lake2 = lake2

    return sink_grid
