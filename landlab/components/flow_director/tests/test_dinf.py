
import pytest
from landlab import VoronoiDelaunayGrid

from landlab.components.flow_director import flow_direction_dinf


def test_not_implemented_voroni():
    x = [0, 0.1, 0.2, 0.3, 1, 1.1, 1.2, 1.3, 2, 2.1, 2.2, 2.3]
    y = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    vmg = VoronoiDelaunayGrid(x, y)
    with pytest.raises(NotImplementedError):
        flow_direction_dinf.flow_directions_dinf(vmg)
