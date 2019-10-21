import numpy as np

from landlab import RasterModelGrid
from landlab.components.overland_flow.generate_overland_flow_deAlmeida import (
    _active_links_at_node,
)


def test_active_links_at_node_scalar_interior():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert np.all(_active_links_at_node(grid, [6]) == np.array([[5, 9, 14, 10]]).T)


def test_active_links_at_node_scalar_boundary():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert np.all(_active_links_at_node(grid, [1]) == np.array([[-1, -1, 5, -1]]).T)


def test_active_node_with_array_arg():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert np.all(
        _active_links_at_node(grid, [6, 7])
        == np.array([[5, 9, 14, 10], [6, 10, 15, 11]]).T
    )


def test_active_links_at_node_with_no_args():
    grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    assert np.all(
        _active_links_at_node(grid)
        == np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    5,
                    6,
                    7,
                    -1,
                    -1,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    23,
                    24,
                    25,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    18,
                    19,
                    20,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    5,
                    6,
                    7,
                    -1,
                    -1,
                    14,
                    15,
                    16,
                    -1,
                    -1,
                    23,
                    24,
                    25,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    9,
                    10,
                    11,
                    12,
                    -1,
                    18,
                    19,
                    20,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        )
    )


def test_inactive_interiors():
    grid = RasterModelGrid((4, 5))
    nodes_to_close = ((6, 12),)
    grid.status_at_node[nodes_to_close] = grid.BC_NODE_IS_CLOSED
    np.all(
        _active_links_at_node(grid)
        == np.array(
            [
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    6,
                    7,
                    -1,
                    -1,
                    -1,
                    -1,
                    16,
                    -1,
                    -1,
                    23,
                    -1,
                    25,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    11,
                    12,
                    -1,
                    18,
                    -1,
                    -1,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    -1,
                    6,
                    7,
                    -1,
                    -1,
                    -1,
                    -1,
                    16,
                    -1,
                    -1,
                    23,
                    -1,
                    25,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
                [
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    11,
                    12,
                    -1,
                    18,
                    -1,
                    -1,
                    21,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                    -1,
                ],
            ]
        )
    )
