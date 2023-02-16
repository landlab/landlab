import pytest

from landlab import RasterModelGrid
from landlab.grid.mappers import map_upwind_node_link_max_to_node


def test_map_upwind_node_link_max_to_node():
    grid = RasterModelGrid((3, 4))

    grid.at_link["grad"] = (
        [-1.0, -2.0, -1.0]
        + [0.0, 0.0, 0.0, 0.0]
        + [-1.0, -2.0, -1.0]
        + [0.0, 0.0, 0.0, 0.0]
        + [-1.0, -2.0, -1.0]
    )

    expected = [0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0, 0.0, 1.0, 2.0, 1.0]

    assert map_upwind_node_link_max_to_node(grid, "grad") == pytest.approx(expected)

    out = grid.empty(at="node")
    rtn = map_upwind_node_link_max_to_node(grid, "grad", out=out)
    assert rtn is out
    assert out == pytest.approx(expected)
