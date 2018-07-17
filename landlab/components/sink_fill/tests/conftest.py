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
