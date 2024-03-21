import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab import VoronoiDelaunayGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorDINF
from landlab.components.flow_director import flow_direction_dinf


def test_not_implemented_voroni():
    x = [0, 0.1, 0.2, 0.3, 1, 1.1, 1.2, 1.3, 2, 2.1, 2.2, 2.3]
    y = [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3]
    vmg = VoronoiDelaunayGrid(x, y)
    with pytest.raises(NotImplementedError):
        flow_direction_dinf.flow_directions_dinf(vmg)


def test_D_infinity_low_closed_boundary_conditions():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = np.array(
        [[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0], [0, 0, 0, 0]],
        dtype="float64",
    )
    mg.add_field("topographic__elevation", z, at="node")
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    true_receivers = np.array(
        [
            [0, -1],
            [1, -1],
            [2, -1],
            [3, -1],
            [4, -1],
            [6, -1],
            [6, -1],
            [7, -1],
            [8, -1],
            [10, 6],
            [6, -1],
            [11, -1],
            [12, -1],
            [10, -1],
            [10, -1],
            [15, -1],
            [16, -1],
            [17, -1],
            [18, -1],
            [19, -1],
        ]
    )

    true_proportions = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.06058469, 0.93941531],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
        ]
    )
    assert_array_equal(fd._receivers, true_receivers)
    assert_array_equal(
        np.round(fd._proportions, decimals=6), np.round(true_proportions, decimals=6)
    )


def test_D_infinity_open_boundary_conditions():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = mg.x_of_node + 2.0 * mg.y_of_node
    mg.add_field("topographic__elevation", z, at="node")

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    true_receivers = np.array(
        [
            [0, -1],
            [1, -1],
            [2, -1],
            [3, -1],
            [4, -1],
            [1, 0],
            [2, 1],
            [7, -1],
            [8, -1],
            [5, 4],
            [6, 5],
            [11, -1],
            [12, -1],
            [9, 8],
            [10, 9],
            [15, -1],
            [16, -1],
            [17, -1],
            [18, -1],
            [19, -1],
        ]
    )

    true_proportions = np.array(
        [
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.40966553, 0.59033447],
            [0.40966553, 0.59033447],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.40966553, 0.59033447],
            [0.40966553, 0.59033447],
            [1.0, 0.0],
            [1.0, 0.0],
            [0.40966553, 0.59033447],
            [0.40966553, 0.59033447],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
            [1.0, 0.0],
        ]
    )
    assert_array_equal(fd._receivers, true_receivers)
    assert_array_equal(
        np.round(fd._proportions, decimals=6), np.round(true_proportions, decimals=6)
    )


def test_D_infinity_flat():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    mg.add_zeros("topographic__elevation", at="node")

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fd._receivers.shape)
    true_receivers[:, 0] = node_ids

    true_proportions = np.zeros(fd._proportions.shape)
    true_proportions[:, 0] = 1

    assert_array_equal(fd._receivers, true_receivers)
    assert_array_equal(
        np.round(fd._proportions, decimals=6), np.round(true_proportions, decimals=6)
    )


def test_D_infinity_flat_closed_lower():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[mg.core_nodes] += 1
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fd._receivers.shape)
    true_receivers[:, 0] = node_ids

    true_proportions = np.zeros(fd._proportions.shape)
    true_proportions[:, 0] = 1

    assert_array_equal(fd._receivers, true_receivers)
    assert_array_equal(
        np.round(fd._proportions, decimals=6), np.round(true_proportions, decimals=6)
    )


def test_D_infinity_flat_closed_upper():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[mg.core_nodes] -= 1
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fd = FlowDirectorDINF(mg)
    fd.run_one_step()

    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fd._receivers.shape)
    true_receivers[:, 0] = node_ids

    true_proportions = np.zeros(fd._proportions.shape)
    true_proportions[:, 0] = 1

    assert_array_equal(fd._receivers, true_receivers)
    assert_array_equal(
        np.round(fd._proportions, decimals=6), np.round(true_proportions, decimals=6)
    )


def test_D_infinity_SW_slope():
    mg = RasterModelGrid((10, 10))
    mg.add_field("topographic__elevation", mg.node_y + mg.node_x, at="node")
    fa = FlowAccumulator(mg, flow_director="FlowDirectorDINF")
    fa.run_one_step()

    # this one should all flow to the soutwest (third column of diagonal neighbors at node)
    node_ids = np.arange(mg.number_of_nodes)
    sw_diags = mg.diagonal_adjacent_nodes_at_node[:, 2]
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[:, 0] = sw_diags
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[:, 0] = 1

    assert_array_equal(true_receivers, fa.flow_director._receivers)
    assert_array_equal(true_proportions, fa.flow_director._proportions)


def test_D_infinity_WSW_slope():
    mg = RasterModelGrid((10, 10))
    mg.add_field(
        "topographic__elevation", mg.node_y * (2**0.5 - 1.0) + mg.node_x, at="node"
    )
    fa = FlowAccumulator(mg, flow_director="FlowDirectorDINF")
    fa.run_one_step()

    # this one should flow equally to west and southwest.
    node_ids = np.arange(mg.number_of_nodes)
    sw_diags = mg.diagonal_adjacent_nodes_at_node[:, 2]
    w_links = mg.adjacent_nodes_at_node[:, 2]
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[mg.core_nodes, 0] = w_links[mg.core_nodes]
    true_receivers[mg.core_nodes, 1] = sw_diags[mg.core_nodes]
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[mg.boundary_nodes, 0] = 1
    true_proportions[mg.core_nodes, 0] = 0.5
    true_proportions[mg.core_nodes, 1] = 0.5

    assert_array_equal(true_receivers, fa.flow_director._receivers)
