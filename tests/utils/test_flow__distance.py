import math

import numpy as np
import pytest
from numpy.testing import assert_almost_equal
from numpy.testing import assert_array_equal

from landlab import FieldError
from landlab import HexModelGrid
from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import FlowDirectorSteepest
from landlab.utils.flow__distance import calculate_flow__distance


def test_no_flow_receivers():
    """Test that correct error is raised when no flow recievers are
    on the grid."""

    # instantiate a model grid, do not run flow accumulation on it

    mg = RasterModelGrid((30, 70))

    # test that the flow distance utility will fail because of a ValueError

    with pytest.raises(FieldError):
        calculate_flow__distance(mg)


def test_no_upstream_array():
    """Test that correct error is raised when no flow__upstream_node_order."""

    # instantiate a model grid, do not run flow accumulation on it

    mg = RasterModelGrid((30, 70))

    # Add a field called topographic__elevation to mg

    mg.add_ones("topographic__elevation", at="node")

    # Run the FlowDirectorSteepest component

    fd = FlowDirectorSteepest(mg)
    fd.run_one_step()

    # test that the flow distance utility will fail because of a ValueError

    with pytest.raises(FieldError):
        calculate_flow__distance(mg)


def test_flow__distance_regular_grid_d8():
    """Test to demonstrate that flow__distance utility works as expected with
    regular grids"""

    # instantiate a model grid

    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))

    # instantiate an elevation array

    z = np.array(
        [[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0], [0, 0, 0, 0]],
        dtype="float64",
    )

    # add the elevation field to the grid

    mg.add_field("topographic__elevation", z, at="node")

    # instantiate the expected flow__distance array
    # considering flow directions calculated with D8 algorithm

    flow__distance_expected = np.array(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [0, math.sqrt(2), 1, 0],
            [0, 1 + math.sqrt(2), 2, 0],
            [0, 0, 0, 0],
        ],
        dtype="float64",
    )
    flow__distance_expected = np.reshape(
        flow__distance_expected, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # setting boundary conditions

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # calculating flow directions with FlowAccumulator component

    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)
    flow__distance = np.reshape(
        flow__distance, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # test that the flow distance utility works as expected

    assert_array_equal(flow__distance_expected, flow__distance)


def test_flow__distance_regular_grid_d4():
    """Test to demonstrate that flow__distance utility works as expected with
    regular grids"""

    # instantiate a model grid

    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))

    # instantiate an elevation array

    z = np.array(
        [[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0], [0, 0, 0, 0]],
        dtype="float64",
    )

    # add the elevation field to the grid

    mg.add_field("topographic__elevation", z, at="node")

    # instantiate the expected flow__distance array
    # considering flow directions calculated with D4 algorithm

    flow__distance_expected = np.array(
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 2, 1, 0], [0, 3, 2, 0], [0, 0, 0, 0]],
        dtype="float64",
    )
    flow__distance_expected = np.reshape(
        flow__distance_expected, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # setting boundary conditions

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # calculating flow directions with FlowAccumulator component

    fr = FlowAccumulator(mg, flow_director="D4")
    fr.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)
    flow__distance = np.reshape(
        flow__distance, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # test that the flow__distance utility works as expected

    assert_array_equal(flow__distance_expected, flow__distance)


def test_flow__distance_irregular_grid_d4():
    """Test to demonstrate that flow__distance utility works as expected with irregular grids"""

    # instantiate a model grid

    dx = 1.0
    hmg = HexModelGrid((5, 3), spacing=dx)

    # instantiate and add the elevation field

    hmg.add_field(
        "topographic__elevation", hmg.node_x + np.round(hmg.node_y), at="node"
    )

    # instantiate the expected flow__distance array

    flow__distance_expected = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            dx,
            0.0,
            0.0,
            dx,
            dx,
            2.0 * dx,
            0.0,
            0.0,
            2.0 * dx,
            2.0 * dx,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    # setting boundary conditions

    hmg.status_at_node[hmg.boundary_nodes] = hmg.BC_NODE_IS_CLOSED

    # calculating flow directions with FlowAccumulator component: D4 algorithm

    fr = FlowAccumulator(hmg, flow_director="D4")
    fr.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(hmg, add_to_grid=True, clobber=True)

    # test that the flow__distance utility works as expected

    assert_almost_equal(flow__distance_expected, flow__distance, decimal=10)


def test_flow__distance_raster_MFD_diagonals_true():
    """Test of flow__distance utility with a raster grid and MFD."""

    # instantiate a model grid

    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))

    # instantiate an elevation array

    z = np.array(
        [[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0], [0, 0, 0, 0]],
        dtype="float64",
    )

    # add the elevation field to the grid

    mg.add_field("topographic__elevation", z, at="node")

    # instantiate the expected flow__distance array
    # considering flow directions calculated with MFD algorithm

    flow__distance_expected = np.array(
        [
            [0, 0, 0, 0],
            [0, 1, 0, 0],
            [0, math.sqrt(2), 1, 0],
            [0, 1 + math.sqrt(2), 2, 0],
            [0, 0, 0, 0],
        ],
        dtype="float64",
    )
    flow__distance_expected = np.reshape(
        flow__distance_expected, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # setting boundary conditions

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # calculating flow directions with FlowAccumulator component

    fa = FlowAccumulator(
        mg, "topographic__elevation", flow_director="MFD", diagonals=True
    )
    fa.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)

    # test that the flow__distance utility works as expected

    assert_array_equal(flow__distance_expected, flow__distance)


def test_flow__distance_raster_MFD_diagonals_false():
    """Test of flow__distance utility with a raster grid and MFD."""

    # instantiate a model grid

    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))

    # instantiate an elevation array

    z = np.array(
        [[0, 0, 0, 0], [0, 21, 10, 0], [0, 31, 20, 0], [0, 32, 30, 0], [0, 0, 0, 0]],
        dtype="float64",
    )

    # add the elevation field to the grid

    mg.add_field("topographic__elevation", z, at="node")

    # instantiate the expected flow__distance array
    # considering flow directions calculated with MFD algorithm

    flow__distance_expected = np.array(
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 2, 1, 0], [0, 3, 2, 0], [0, 0, 0, 0]],
        dtype="float64",
    )
    flow__distance_expected = np.reshape(
        flow__distance_expected, mg.number_of_node_rows * mg.number_of_node_columns
    )

    # setting boundary conditions

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # calculating flow directions with FlowAccumulator component

    fa = FlowAccumulator(
        mg, "topographic__elevation", flow_director="MFD", diagonals=False
    )
    fa.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)

    # test that the flow__distance utility works as expected

    assert_array_equal(flow__distance_expected, flow__distance)


def test_flow__distance_raster_D_infinity():
    """Test of flow__distance utility with a raster grid and D infinity."""

    mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))

    # instantiate an elevation array

    z = mg.x_of_node + 3.0 * mg.y_of_node

    # add the elevation field to the grid

    mg.add_field("topographic__elevation", z, at="node")

    # instantiate the expected flow_length array

    flow__distance_expected = np.array(
        [
            [0, 0, 0, 0],
            [0, 0, 1, 0],
            [0, 1, 0 + math.sqrt(2.0), 0],
            [0, 2, 1 + math.sqrt(2.0), 0],
            [0, 0, 0, 0],
        ],
        dtype="float64",
    )

    # setting boundary conditions

    mg.set_closed_boundaries_at_grid_edges(
        bottom_is_closed=True,
        left_is_closed=True,
        right_is_closed=True,
        top_is_closed=True,
    )

    # calculating flow directions with FlowAccumulator component

    fa = FlowAccumulator(mg, "topographic__elevation", flow_director="DINF")
    fa.run_one_step()

    # calculating flow distance map

    flow__distance = calculate_flow__distance(
        mg, add_to_grid=True, clobber=True
    ).reshape(mg.shape)

    # test that the flow__distance utility works as expected

    assert_array_equal(flow__distance_expected, flow__distance)
