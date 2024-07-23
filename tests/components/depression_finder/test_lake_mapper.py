"""
test_lake_mapper:

Created on Sun Sep 27 09:52:50, 2015

@author: gtucker, amended dejh
"""

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import RasterModelGrid
from landlab.components import DepressionFinderAndRouter
from landlab.components import FlowAccumulator
from landlab.components.depression_finder.cfuncs import (
    find_lowest_node_on_lake_perimeter_c,
)

XX = RasterModelGrid.BAD_INDEX


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    mg.add_zeros("topographic__elevation", at="node")
    FlowAccumulator(mg, flow_director="MFD").run_one_step()

    with pytest.raises(NotImplementedError):
        DepressionFinderAndRouter(mg)


def test_lake_mapper():
    """Create a test grid and run a series of tests."""
    rmg = RasterModelGrid((8, 8))

    # Make topography into sinusoidal hills and depressions
    rmg.add_field(
        "topographic__elevation",
        np.sin(2 * np.pi * rmg.x_of_node / 8.0)
        * np.sin(2 * np.pi * rmg.y_of_node / 4.0),
        at="node",
    )

    # Set 3 sides of the grid to be closed boundaries
    rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    lm = DepressionFinderAndRouter(rmg)

    lm.map_depressions()

    assert "topographic__elevation" in rmg.at_node
    assert "depression__depth" in rmg.at_node
    assert "depression__outlet_node" in rmg.at_node

    is_pit = np.full(rmg.number_of_nodes, False)
    is_pit[lm.pit_node_ids] = True
    assert_array_equal(
        is_pit.reshape(rmg.shape),
        [
            [False, False, False, False, False, False, False, False],
            [False, False, False, False, False, False, True, False],
            [False, False, False, False, False, False, False, False],
            [False, False, True, False, False, False, False, False],
            [False, False, False, False, False, False, False, False],
            [False, False, False, False, False, False, True, False],
            [False, True, False, False, False, False, False, False],
            [False, False, False, False, False, False, False, False],
        ],
    )
    assert_array_equal(
        lm.flood_status.reshape(rmg.shape),
        [
            [0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 3, 3, 0],
            [0, 0, 0, 0, 0, 3, 3, 0],
            [0, 3, 3, 3, 3, 0, 0, 0],
            [0, 3, 3, 3, 3, 0, 0, 0],
            [0, 0, 0, 0, 0, 3, 3, 0],
            [0, 3, 0, 0, 0, 3, 3, 0],
            [0, 0, 0, 0, 0, 0, 0, 0],
        ],
    )
    assert_array_almost_equal(
        lm.depression_depth.reshape(rmg.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.71, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.71, 1.0, 0.71, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.71, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
        decimal=2,
    )
    assert rmg.at_node["depression__depth"] is lm.depression_depth
    assert_array_equal(
        lm.depression_outlet_map.reshape(rmg.shape),
        [
            [XX, XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, 5, 5, XX],
            [XX, XX, XX, XX, XX, 5, 5, XX],
            [XX, 5, 5, 5, 5, XX, XX, XX],
            [XX, 5, 5, 5, 5, XX, XX, XX],
            [XX, XX, XX, XX, XX, 5, 5, XX],
            [XX, 50, XX, XX, XX, 5, 5, XX],
            [XX, XX, XX, XX, XX, XX, XX, XX],
        ],
    )
    assert rmg.at_node["depression__outlet_node"] is lm.depression_outlet_map

    assert "topographic__elevation" in rmg.at_node
    assert "depression__depth" in rmg.at_node
    assert "depression__outlet_node" in rmg.at_node


def test_initial_routing(dans_grid3):
    """
    Test the action of fr.run_one_step() on the grid.
    """
    dans_grid3.fr.run_one_step()
    assert dans_grid3.mg.at_node["flow__receiver_node"] == approx(dans_grid3.r_old)
    assert dans_grid3.mg.at_node["drainage_area"] == approx(dans_grid3.A_old)


def test_rerouting_with_supplied_pits(dans_grid3):
    """
    Test with the output from a successful run of fr.run_one_step.
    """
    dans_grid3.fr.run_one_step()
    assert_array_equal(
        dans_grid3.mg.at_node["flow__link_to_receiver_node"], dans_grid3.links_old
    )
    dans_grid3.lf.map_depressions()
    assert_array_equal(dans_grid3.mg.at_node["flow__receiver_node"], dans_grid3.r_new)
    assert dans_grid3.mg.at_node["drainage_area"] == approx(dans_grid3.A_new)
    assert dans_grid3.mg.at_node["surface_water__discharge"] == approx(dans_grid3.A_new)
    assert dans_grid3.mg.at_node["flow__upstream_node_order"] == approx(
        dans_grid3.s_new
    )
    assert dans_grid3.mg.at_node["flow__link_to_receiver_node"] == approx(
        dans_grid3.links_new
    )


def test_changing_slopes(dans_grid3):
    """
    Test with the output from a successful run of fr.run_one_step.
    """
    grid = dans_grid3.mg

    slope_old = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
        [0.0, 2.0, 0.1, 0.0, 0.1, 2.0, 0.0],
        [0.0, 2.0, 0.14142136, 0.1, 0.14142136, 2.0, 0.0],
        [0.0, 2.0, 1.2, 1.0, 1.0, 2.0, 0.0],
        [0.0, 1.06066017, 1.1, 1.06066017, 1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
    slope_new = [
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0],
        [0.0, 2.0, 0.0, 0.0, 0.0, 2.0, 0.0],
        [0.0, 2.0, 0.0, 0.0, 0.1, 2.0, 0.0],
        [0.0, 2.0, 1.2, 1.0, 1.0, 2.0, 0.0],
        [0.0, 1.06066017, 1.1, 1.06066017, 1.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]

    dans_grid3.fr.run_one_step()
    assert_array_almost_equal(
        grid.at_node["topographic__steepest_slope"].reshape(grid.shape), slope_old
    )

    dans_grid3.lf.map_depressions()
    assert_array_almost_equal(
        grid.at_node["topographic__steepest_slope"].reshape(grid.shape), slope_new
    )


def test_filling_alone(dans_grid3):
    """
    Test the filler alone, w/o supplying information on the pits.

    Setting the the *pits* parameter to None causes the mapper to look for pits
    using its _find_pits method.
    """
    lf = DepressionFinderAndRouter(dans_grid3.mg, reroute_flow=False, pits=None)
    assert lf._user_supplied_pits is None

    lf.map_depressions()
    assert_array_equal(
        dans_grid3.mg.at_node["flow__receiver_node"],
        np.full(dans_grid3.mg.number_of_nodes, -1, dtype=int),
    )
    assert_array_equal(lf.depression_outlet_map, dans_grid3.depr_outlet_target)


def test_filling_supplied_pits(dans_grid3):
    """
    Test the filler without rereouting, but confusingly, where there *is*
    aready routing information available!
    Also tests the supply of an array for 'pits'
    """
    dans_grid3.fr.run_one_step()

    lf = DepressionFinderAndRouter(
        dans_grid3.mg, reroute_flow=False, pits="flow__sink_flag"
    )
    lf.map_depressions()
    assert_array_equal(dans_grid3.mg.at_node["flow__receiver_node"], dans_grid3.r_old)


def test_pits_as_IDs(dans_grid3):
    """
    Smoke test for passing specific IDs, not an array, to the mapper.
    """
    dans_grid3.fr.run_one_step()

    pits = np.nonzero(dans_grid3.mg.at_node["flow__sink_flag"])[0]
    lf = DepressionFinderAndRouter(dans_grid3.mg, pits=pits)
    lf.map_depressions()

    assert dans_grid3.mg.at_node["drainage_area"] == approx(dans_grid3.A_new)


def test_edge_draining():
    """
    This tests when the lake attempts to drain from an edge, where an issue
    is suspected.
    """
    # Create a 7x7 test grid with a well defined hole in it, AT THE EDGE.
    mg = RasterModelGrid((7, 7))

    mg.add_field(
        "topographic__elevation",
        [
            [-2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0],
            [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0],
            [0.0, -1.0, -1.0, 3.0, 4.0, 5.0, 6.0],
            [0.0, -1.0, -1.0, 3.0, 4.0, 5.0, 6.0],
            [0.0, -1.0, -1.0, 3.0, 4.0, 5.0, 6.0],
            [6.0, 6.0, 6.0, 6.0, 6.0, 6.0, 6.0],
            [-2.0, -2.0, -2.0, -2.0, -2.0, -2.0, -2.0],
        ],
        at="node",
        units="-",
    )

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)

    fr.run_one_step()
    lf.map_depressions()
    assert_array_almost_equal(
        mg.at_node["drainage_area"].reshape(mg.shape),
        [
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [15.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [0.0, 10.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
        ],
    )
    assert_array_almost_equal(
        lf.depression_outlet_map.reshape(mg.shape),
        [
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, 14, 14, XX, XX, XX, XX],
            [XX, 14, 14, XX, XX, XX, XX],
            [XX, 14, 14, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX],
        ],
    )


def test_degenerate_drainage():
    """
    This "hourglass" configuration should be one of the hardest to correctly
    re-route.
    """
    mg = RasterModelGrid((9, 5))
    mg.add_field(
        "topographic__elevation",
        [
            [1.0, 1.0001, 1.0002, 1.0003, 1.0004],
            [1.0, 1.0001, -1.0, 1.0003, 1.0004],
            [1.0, -1.0, -1.0, -1.0, 1.0004],
            [1.0, 1.0001, -1.0, 1.0003, 1.0004],
            [-0.2, 0.1, 0.0, 1.0003, 1.0004],
            [1.0, 1.0001, -1.0, 1.0003, 1.0004],
            [1.0, -1.0, -1.0, -1.0, 1.0004],
            [1.0, 1.0001, -1.0, 1.0003, 1.0004],
            [1.0, 1.0001, 1.0002, 1.0003, 1.0004],
        ],
        at="node",
    )

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    assert_array_almost_equal(
        mg.at_node["drainage_area"].reshape(mg.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 3.0, 1.0, 0.0],
            [0.0, 2.0, 4.0, 2.0, 0.0],
            [0.0, 1.0, 10.0, 1.0, 0.0],
            [21.0, 21.0, 1.0, 1.0, 0.0],
            [0.0, 1.0, 9.0, 1.0, 0.0],
            [0.0, 2.0, 2.0, 2.0, 0.0],
            [0.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )


def test_three_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits.
    """
    mg = RasterModelGrid((10, 10))
    z = [
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 1.0, 4.0, 5.0, 6.0, 4.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 1.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 1.0, 1.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    ]
    mg.add_field("topographic__elevation", z, at="node")

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node["flow__sink_flag"], flow_sinks_target)

    # test conservation of mass:
    assert mg.at_node["drainage_area"].reshape(mg.shape)[1:-1, 1].sum() == approx(
        (mg.shape[0] - 2) ** 2
    )

    assert_array_equal(
        mg.at_node["drainage_area"].reshape(mg.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [2.0, 2.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [26.0, 26.0, 25.0, 15.0, 11.0, 10.0, 9.0, 8.0, 1.0, 0.0],
            [2.0, 2.0, 1.0, 9.0, 2.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [2.0, 2.0, 1.0, 1.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.0],
            [20.0, 20.0, 19.0, 18.0, 17.0, 12.0, 3.0, 2.0, 1.0, 0.0],
            [2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 3.0, 2.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )

    assert_array_equal(
        lf.lake_map.reshape(mg.shape),
        [
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, 33, -1, -1, -1, 37, -1, -1],
            [-1, -1, -1, 33, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, 74, 74, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        ],
    )
    assert_array_equal(lf.lake_codes, [33, 37, 74])
    assert lf.number_of_lakes == 3
    assert lf.lake_areas == approx([2.0, 1.0, 2.0])
    assert lf.lake_volumes == approx([2.0, 2.0, 4.0])


def test_composite_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits, inset into each other.
    """
    mg = RasterModelGrid((10, 10))
    z = [
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 0.0, -2.0, 0.0, 0.0, 0.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 0.0, -10.0, 0.0, 0.0, -1.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 9.0],
        [0.0, 0.9, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0],
    ]
    mg.add_field("topographic__elevation", z, at="node")

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node["flow__sink_flag"], flow_sinks_target)

    # test conservation of mass:
    assert mg.at_node["drainage_area"].reshape(mg.shape)[1:-1, 1].sum() == approx(
        (mg.shape[0] - 2) ** 2
    )

    assert_array_equal(
        mg.at_node["drainage_area"].reshape(mg.shape),
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [8.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 4.0, 2.0, 2.0, 6.0, 4.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 6.0, 3.0, 12.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 8.0, 20.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 35.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0],
            [50.0, 50.0, 49.0, 13.0, 10.0, 8.0, 6.0, 4.0, 1.0, 0.0],
            [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ],
    )

    # the lake code map:
    assert_array_equal(
        lf.lake_map.reshape(mg.shape),
        [
            [XX, XX, XX, XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, 57, 57, 57, 57, 57, XX, XX],
            [XX, XX, XX, 57, 57, 57, 57, 57, XX, XX],
            [XX, XX, XX, 57, 57, 57, 57, 57, XX, XX],
            [XX, XX, XX, 57, 57, 57, 57, 57, XX, XX],
            [XX, XX, XX, 57, 57, 57, 57, 57, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX, XX, XX, XX],
            [XX, XX, XX, XX, XX, XX, XX, XX, XX, XX],
        ],
    )

    # test the remaining properties:
    assert_array_equal(lf.lake_outlets, [72])
    assert_array_equal(
        lf.depression_outlet_map.reshape(mg.shape),
        [
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, 72, 72, 72, 72, 72, -1, -1],
            [-1, -1, -1, 72, 72, 72, 72, 72, -1, -1],
            [-1, -1, -1, 72, 72, 72, 72, 72, -1, -1],
            [-1, -1, -1, 72, 72, 72, 72, 72, -1, -1],
            [-1, -1, -1, 72, 72, 72, 72, 72, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1, -1, -1, -1, -1],
        ],
    )
    assert lf.number_of_lakes == 1
    assert_array_equal(lf.lake_codes, [57])
    assert_array_almost_equal(lf.lake_areas, [25.0])
    assert_array_almost_equal(lf.lake_volumes, [63.0])


def test_D8_D4_fill(d4_grid):
    """
    Tests the functionality of D4 filling.
    """
    lfD8 = DepressionFinderAndRouter(
        d4_grid.mg1, pits=None, routing="D8", reroute_flow=False
    )
    lfD4 = DepressionFinderAndRouter(
        d4_grid.mg2, pits=None, routing="D4", reroute_flow=False
    )

    lfD8.map_depressions()
    lfD4.map_depressions()

    assert lfD8.number_of_lakes == 1
    assert lfD4.number_of_lakes == 3

    correct_D8_lake_map = np.empty(7 * 7, dtype=int)
    correct_D8_lake_map.fill(XX)
    correct_D8_lake_map[d4_grid.lake_nodes] = 10
    correct_D4_lake_map = correct_D8_lake_map.copy()
    correct_D4_lake_map[d4_grid.lake_nodes[5:]] = 32
    correct_D4_lake_map[d4_grid.lake_nodes[-2]] = 38
    correct_D8_depths = np.zeros(7 * 7, dtype=float)
    correct_D8_depths[d4_grid.lake_nodes] = 2.0
    correct_D4_depths = correct_D8_depths.copy()
    correct_D4_depths[d4_grid.lake_nodes[5:]] = 4.0
    correct_D4_depths[d4_grid.lake_nodes[-2]] = 3.0

    assert_array_equal(lfD8.lake_map, correct_D8_lake_map)
    assert_array_equal(lfD4.lake_map, correct_D4_lake_map)

    assert d4_grid.mg1.at_node["depression__depth"] == approx(correct_D8_depths)
    assert d4_grid.mg2.at_node["depression__depth"] == approx(correct_D4_depths)


def test_D8_D4_route(d4_grid):
    """
    Tests the functionality of D4 routing.
    """
    d4_grid.frD8.run_one_step()
    d4_grid.frD4.run_one_step()
    d4_grid.lfD8.map_depressions()
    d4_grid.lfD4.map_depressions()

    assert d4_grid.lfD8.number_of_lakes == 1
    assert d4_grid.lfD4.number_of_lakes == 3

    assert_array_equal(
        d4_grid.mg1.at_node["flow__receiver_node"].reshape((7, 7)),
        [
            [0, 1, 2, 3, 4, 5, 6],
            [7, 16, 10, 16, 10, 18, 13],
            [14, 14, 15, 16, 17, 18, 20],
            [21, 16, 16, 16, 18, 33, 27],
            [28, 28, 24, 24, 24, 32, 34],
            [35, 35, 38, 32, 32, 32, 41],
            [42, 43, 44, 45, 46, 47, 48],
        ],
    )
    assert_array_equal(
        d4_grid.mg2.at_node["flow__receiver_node"].reshape((7, 7)),
        [
            [0, 1, 2, 3, 4, 5, 6],
            [7, 7, 10, 17, 10, 11, 13],
            [14, 14, 15, 16, 17, 18, 20],
            [21, 21, 16, 17, 18, 33, 27],
            [28, 28, 29, 38, 31, 32, 34],
            [35, 35, 36, 37, 32, 33, 41],
            [42, 43, 44, 45, 46, 47, 48],
        ],
    )
    assert d4_grid.mg1.at_node["drainage_area"].reshape((7, 7))[:, 0].sum() == approx(
        d4_grid.mg2.at_node["drainage_area"].reshape((7, 7))[:, 0].sum()
    )


def test_find_lowest_node_on_lake_perimeter_c():
    """
    Ensures the key functionality of the cfunc is working.
    """
    mg = RasterModelGrid((7, 7), xy_spacing=0.5)
    mg.add_field(
        "topographic__elevation",
        [
            [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0],
            [0.005, 0.505, 1.005, 1.505, 2.005, 2.505, 3.005],
            [0.01, 0.51, 0.101, 0.151, 0.201, 2.51, 3.01],
            [0.015, 0.515, 0.1015, 0.1515, 0.2015, 2.515, 3.015],
            [0.02, 0.52, 0.102, 0.152, 0.202, 2.52, 3.02],
            [0.025, 0.525, 1.025, 1.525, 2.025, 2.525, 3.025],
            [0.03, 0.53, 1.03, 1.53, 2.03, 2.53, 3.03],
        ],
        at="node",
    )
    fr = FlowAccumulator(mg, flow_director="D8")
    fr.run_one_step()  # the flow "gets stuck" in the hole
    df = DepressionFinderAndRouter(mg)

    nodes_this_depression = mg.zeros("node", dtype=int)
    nodes_this_depression[0] = 16
    pit_count = 1

    assert find_lowest_node_on_lake_perimeter_c(
        df._neighbor_nodes_at_node,
        df.flood_status,
        mg.at_node["topographic__elevation"],
        nodes_this_depression,
        pit_count,
    ) == (23, 1)

    nodes_this_depression[1] = 8
    pit_count = 2
    assert find_lowest_node_on_lake_perimeter_c(
        df._neighbor_nodes_at_node,
        df.flood_status,
        mg.at_node["topographic__elevation"],
        nodes_this_depression,
        pit_count,
    ) == (0, 2)


def test_all_boundaries_are_closed():
    grid = RasterModelGrid((10, 10), xy_spacing=10.0)
    grid.set_closed_boundaries_at_grid_edges(True, True, True, True)
    grid.add_field(
        "topographic__elevation",
        [
            [0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0],
            [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0],
            [20.0, 30.0, 40.0, 50.0, 60.0, 30.0, 80.0, 90.0, 100.0, 110.0],
            [30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0],
            [40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0],
            [50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0],
            [60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0],
            [70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0],
            [80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0],
            [90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0],
        ],
        at="node",
    )
    fa = FlowAccumulator(
        grid, depression_finder="DepressionFinderAndRouter", routing="D4"
    )
    with pytest.raises(ValueError):
        fa.run_one_step()


def test_precision_in_cython():
    """
    This test confirms no issues with Cython precision matching of floats and
    Numpy array C++ doubles. Per alexmitchell @issue 1219 (now fixed).

    Bit of a smoke test on the whole component, as this is how the issue was
    found.
    """
    z = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
            [0.0, 90.0, 90.0, 90.0, 90.0, 0.0],
            [0.0, 90.0, 90.0, 1.0, 90.0, 0.0],
            [0.0, 90.0, 32.0, 32.000001, 90.0, 0.0],  # precision matters here
            [0.0, 90.0, 80.0, 90.0, 90.0, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        ]
    )
    mg = RasterModelGrid(z.shape, xy_spacing=100.0)

    mg.add_field("topographic__elevation", z, units="m")
    mg.set_closed_boundaries_at_grid_edges(False, False, False, False)

    flow_router = FlowAccumulator(
        mg, flow_director="FlowDirectorD8", depression_finder=DepressionFinderAndRouter
    )
    flow_router.run_one_step()
    flow_router.depression_finder.depression_outlet_map.reshape(mg.shape)

    # formerly, this concluded that 32.000001 < 32., and so terminated the
    # fill algo with only one lake node, but in fact:
    assert_array_equal(
        flow_router.depression_finder.depression_outlet_map.reshape(mg.shape),
        [
            [-1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1],
            [-1, -1, -1, 26, -1, -1],
            [-1, -1, 26, 26, -1, -1],
            [-1, -1, -1, -1, -1, -1],
            [-1, -1, -1, -1, -1, -1],
        ],
    )
