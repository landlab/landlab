import numpy as np
from numpy.testing import assert_array_equal

from landlab import NodeStatus
from landlab import RasterModelGrid

FG = NodeStatus.FIXED_GRADIENT
CB = NodeStatus.CLOSED
FV = NodeStatus.FIXED_VALUE


def test_fixed_link_boundaries_at_grid_edges():
    """test setting fixed link boundaries at grid edges"""

    grid = RasterModelGrid((3, 4))
    grid.at_node["topographic__elevation"] = np.zeros(grid.number_of_nodes)
    grid.at_link["topographic__slope"] = np.zeros(grid.number_of_links)
    grid.status_at_node[grid.nodes_at_bottom_edge] = grid.BC_NODE_IS_FIXED_GRADIENT

    assert_array_equal(
        grid.status_at_node, [FG, FG, FG, FG, FV, 0, 0, FV, FV, FV, FV, FV]
    )

    assert_array_equal(
        grid.status_at_link, [4, 4, 4, 4, 2, 2, 4, 0, 0, 0, 4, 0, 0, 4, 4, 4, 4]
    )


def test_nodata_fixed_links():
    """test setting nodata nodes to fixed gradient"""
    grid = RasterModelGrid((3, 4))
    grid["node"]["topographic__elevation"] = np.zeros(grid.number_of_nodes)
    grid["link"]["topographic__slope"] = np.zeros(grid.number_of_links)
    z = grid["node"]["topographic__elevation"]
    z[3:5] = -9999
    grid.set_nodata_nodes_to_fixed_gradient(z, -9999)
    assert_array_equal(
        grid.status_at_node, [FV, FV, FV, FG, FG, 0, 0, FV, FV, FV, FV, FV]
    )

    assert_array_equal(
        grid.status_at_link, [4, 4, 4, 4, 0, 0, 4, 2, 0, 0, 4, 0, 0, 4, 4, 4, 4]
    )


def test_fixed_gradient_and_value_boundary():
    """testing multiple boundary conditions with fixed links"""
    grid = RasterModelGrid((3, 4))
    grid["node"]["topographic__elevation"] = np.zeros(grid.number_of_nodes)
    grid["link"]["topographic__slope"] = np.zeros(grid.number_of_links)

    grid.set_closed_boundaries_at_grid_edges(False, True, False, True)
    grid.status_at_node[4] = FG
    grid.status_at_node[7] = FG
    grid.status_at_node[1] = 1

    assert_array_equal(
        grid.status_at_node, [CB, FV, CB, CB, FG, 0, 0, FG, CB, CB, CB, CB]
    )

    assert_array_equal(
        grid.status_at_link, [4, 4, 4, 4, 0, 4, 4, 2, 0, 2, 4, 4, 4, 4, 4, 4, 4]
    )
