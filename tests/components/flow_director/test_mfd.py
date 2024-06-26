import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorMFD
from landlab.components.flow_director import flow_direction_mfd


def test_bad_argument_mfd():
    mg = RasterModelGrid((5, 5), xy_spacing=(1, 1))
    z = mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    neighbors_at_node = mg.adjacent_nodes_at_node
    links_at_node = mg.links_at_node
    active_link_dir_at_node = mg.active_link_dirs_at_node
    link_slope = np.arctan(mg.calc_grad_at_link(z))
    link_slope[links_at_node] * active_link_dir_at_node

    with pytest.raises(ValueError):
        flow_direction_mfd.flow_directions_mfd(
            z,
            neighbors_at_node,
            links_at_node,
            active_link_dir_at_node,
            link_slope,
            partition_method="foo",
        )


def test_mfd_on_flat_terrain():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    mg.add_zeros("topographic__elevation", at="node")

    fd = FlowDirectorMFD(mg)
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


def test_mfd_flat_closed_lower():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[mg.core_nodes] += 1
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fd = FlowDirectorMFD(mg)
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


def test_mfd_flat_closed_upper():
    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    z = mg.add_zeros("topographic__elevation", at="node")
    z[mg.core_nodes] -= 1
    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    fd = FlowDirectorMFD(mg)
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


def test_MFD_SW_slope():
    mg = RasterModelGrid((10, 10))
    mg.add_field("topographic__elevation", mg.node_y + mg.node_x, at="node")
    fa = FlowAccumulator(mg, flow_director="FlowDirectorMFD")
    fa.run_one_step()

    # this one should flow half to the west and half to the south
    w_links = mg.adjacent_nodes_at_node[:, 2]
    s_links = mg.adjacent_nodes_at_node[:, 3]
    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[mg.core_nodes, 2] = w_links[mg.core_nodes]
    true_receivers[mg.core_nodes, 3] = s_links[mg.core_nodes]
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[mg.boundary_nodes, 0] = 1
    true_proportions[mg.core_nodes, 2:] = 0.5

    assert_array_equal(true_receivers, fa.flow_director._receivers)
    assert_array_equal(true_proportions, fa.flow_director._proportions)


def test_MFD_SW_slope_w_diags():
    mg = RasterModelGrid((10, 10))
    mg.add_field("topographic__elevation", mg.node_y + mg.node_x, at="node")
    fa = FlowAccumulator(mg, flow_director="FlowDirectorMFD", diagonals=True)
    fa.run_one_step()

    # this one should part to south west, part to south, and part to west
    sw_diags = mg.diagonal_adjacent_nodes_at_node[:, 2]
    w_links = mg.adjacent_nodes_at_node[:, 2]
    s_links = mg.adjacent_nodes_at_node[:, 3]
    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[mg.core_nodes, 2] = w_links[mg.core_nodes]
    true_receivers[mg.core_nodes, 3] = s_links[mg.core_nodes]
    true_receivers[mg.core_nodes, 6] = sw_diags[mg.core_nodes]
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[mg.boundary_nodes, 0] = 1

    total_sum_of_slopes = 1.0 + 1.0 + (2.0 / 2.0**0.5)

    true_proportions[mg.core_nodes, 2] = 1.0 / total_sum_of_slopes
    true_proportions[mg.core_nodes, 3] = 1.0 / total_sum_of_slopes
    true_proportions[mg.core_nodes, 6] = (2.0 / 2.0**0.5) / total_sum_of_slopes

    assert_array_equal(true_receivers, fa.flow_director._receivers)
    assert_array_almost_equal(true_proportions, fa.flow_director._proportions)


# %%
def test_MFD_S_slope():
    mg = RasterModelGrid((10, 10))
    mg.add_field("topographic__elevation", mg.node_y, at="node")
    fa = FlowAccumulator(mg, flow_director="FlowDirectorMFD")
    fa.run_one_step()

    # this should flow totally to the south
    node_ids = np.arange(mg.number_of_nodes)
    s_links = mg.adjacent_nodes_at_node[:, 3]
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[mg.core_nodes, 3] = s_links[mg.core_nodes]
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[mg.boundary_nodes, 0] = 1
    true_proportions[mg.core_nodes, 3] = 1.0

    assert_array_equal(true_receivers, fa.flow_director._receivers)
    assert_array_equal(true_proportions, fa.flow_director._proportions)


def test_MFD_S_slope_w_diag():
    mg = RasterModelGrid((10, 10))
    mg.add_field("topographic__elevation", mg.node_y, at="node")
    fa = FlowAccumulator(mg, flow_director="FlowDirectorMFD", diagonals=True)
    fa.run_one_step()

    # this one should part to south west, part to south, and part to southeast
    sw_diags = mg.diagonal_adjacent_nodes_at_node[:, 2]
    se_diags = mg.diagonal_adjacent_nodes_at_node[:, 3]
    s_links = mg.adjacent_nodes_at_node[:, 3]
    node_ids = np.arange(mg.number_of_nodes)
    true_receivers = -1 * np.ones(fa.flow_director._receivers.shape)
    true_receivers[mg.core_nodes, 3] = s_links[mg.core_nodes]
    true_receivers[mg.core_nodes, 6] = sw_diags[mg.core_nodes]
    true_receivers[mg.core_nodes, 7] = se_diags[mg.core_nodes]
    true_receivers[mg.boundary_nodes, 0] = node_ids[mg.boundary_nodes]

    true_proportions = np.zeros(fa.flow_director._proportions.shape)
    true_proportions[mg.boundary_nodes, 0] = 1

    total_sum_of_slopes = 1.0 + 2.0 * (1.0 / 2.0**0.5)

    true_proportions[mg.core_nodes, 3] = 1.0 / total_sum_of_slopes
    true_proportions[mg.core_nodes, 6] = (1.0 / 2.0**0.5) / total_sum_of_slopes
    true_proportions[mg.core_nodes, 7] = (1.0 / 2.0**0.5) / total_sum_of_slopes

    assert_array_equal(true_receivers, fa.flow_director._receivers)
    assert_array_almost_equal(true_proportions, fa.flow_director._proportions)
