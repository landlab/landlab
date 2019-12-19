# -*- coding: utf-8 -*-
"""
test_lake_mapper:

Created on Sun Sep 27 09:52:50, 2015

@author: gtucker, amended dejh
"""
import numpy as np  # for use of np.round
import pytest
from numpy import pi, sin
from numpy.testing import assert_array_equal
from pytest import approx

from landlab import BAD_INDEX_VALUE as XX, RasterModelGrid
from landlab.components import DepressionFinderAndRouter, FlowAccumulator

NUM_GRID_ROWS = 8
NUM_GRID_COLS = 8
PERIOD_X = 8.0
PERIOD_Y = 4.0


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        DepressionFinderAndRouter(mg)


def create_test_grid():
    """
    Create a test grid and elevation field with sinusoidal depressions and
    hills.
    """
    # Create grid
    rmg = RasterModelGrid((NUM_GRID_ROWS, NUM_GRID_COLS))

    # Create topography field
    z = rmg.add_zeros("node", "topographic__elevation")

    # Make topography into sinusoidal hills and depressions
    z[:] = sin(2 * pi * rmg.node_x / PERIOD_X) * sin(2 * pi * rmg.node_y / PERIOD_Y)

    # Set 3 sides of the grid to be closed boundaries
    rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    return rmg


def check_fields1(grid):
    """
    Check to make sure the right fields have been created.
    """
    try:
        grid.at_node["topographic__elevation"]
        grid.at_node["flood_status_code"]
        grid.at_node["depression__depth"]
        grid.at_node["depression__outlet_node"]
        grid.at_node["is_pit"]
    except KeyError:
        print("Test failure in check_fields")
        raise


def check_array_values1(rmg, lm):
    """
    Check values of the various fields against known values.
    """

    assert_array_equal(
        lm.is_pit,
        [
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            True,
            False,
            False,
            True,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
            False,
        ],
    )

    assert_array_equal(
        lm.flood_status,
        [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            3,
            3,
            3,
            3,
            0,
            0,
            0,
            0,
            3,
            3,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            3,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
    )

    dd1 = np.round(lm.depression_depth * 100)
    assert_array_equal(
        dd1,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            71.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )

    dd1 = np.round(rmg.at_node["depression__depth"] * 100)
    assert_array_equal(
        dd1,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            71.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )

    assert_array_equal(
        lm.depression_outlet_map,
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            50,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ],
    )

    assert_array_equal(
        rmg.at_node["depression__outlet_node"],
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            50,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ],
    )


def check_fields2(grid):
    """
    Check to make sure the right fields have been created.
    """
    try:
        grid.at_node["topographic__elevation"]
        grid.at_node["flood_status_code"]
        grid.at_node["depression__depth"]
        grid.at_node["depression__outlet_node"]
        grid.at_node["is_pit"]
    except KeyError:
        print("Test failure in check_fields")
        raise


def check_array_values2(rmg, lm):
    """
    Check values of the various fields against known values.
    """
    F = False
    T = True
    assert_array_equal(
        lm.is_pit,
        [
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            T,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            T,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            T,
            F,
            F,
            T,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
            F,
        ],
    )

    assert_array_equal(
        lm.flood_status,
        [
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            3,
            3,
            3,
            3,
            0,
            0,
            0,
            0,
            3,
            3,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            3,
            0,
            0,
            0,
            3,
            3,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
            0,
        ],
    )

    dd1 = np.round(lm.depression_depth * 100)
    assert_array_equal(
        dd1,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            71.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )

    dd1 = np.round(rmg.at_node["depression__depth"] * 100)
    assert_array_equal(
        dd1,
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            71.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            71.0,
            100.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ],
    )

    assert_array_equal(
        lm.depression_outlet_map,
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            50,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ],
    )

    assert_array_equal(
        rmg.at_node["depression__outlet_node"],
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            50,
            XX,
            XX,
            XX,
            5,
            5,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ],
    )


def test_lake_mapper():
    """
    Create a test grid and run a series of tests.
    """
    # Make a test grid
    rmg = create_test_grid()

    # Instantiate a lake mapper
    # (Note that we don't need to send it an input file name, because our grid
    # already has a topographic__elevation field)
    lm = DepressionFinderAndRouter(rmg)

    # Run it on our test grid
    lm.map_depressions()

    # Run tests
    check_fields1(rmg)
    check_array_values1(rmg, lm)
    check_fields2(rmg)
    check_array_values2(rmg, lm)


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
    slope_old = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            0.0,
            0.0,
            2.0,
            0.1,
            0.0,
            0.1,
            2.0,
            0.0,
            0.0,
            2.0,
            0.14142136,
            0.1,
            0.14142136,
            2.0,
            0.0,
            0.0,
            2.0,
            1.2,
            1.0,
            1.0,
            2.0,
            0.0,
            0.0,
            1.06066017,
            1.1,
            1.06066017,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    slope_new = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            2.0,
            2.0,
            2.0,
            2.0,
            2.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            2.0,
            0.0,
            0.0,
            0.1,
            2.0,
            0.0,
            0.0,
            2.0,
            1.2,
            1.0,
            1.0,
            2.0,
            0.0,
            0.0,
            1.06066017,
            1.1,
            1.06066017,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    dans_grid3.fr.run_one_step()
    assert dans_grid3.mg.at_node["topographic__steepest_slope"] == approx(slope_old)
    dans_grid3.lf.map_depressions()
    assert dans_grid3.mg.at_node["topographic__steepest_slope"] == approx(slope_new)


def test_filling_alone(dans_grid3):
    """
    Test the filler alone, w/o supplying information on the pits.

    Setting the the *pits* parameter to None causes the mapper to look for pits
    using its _find_pits method.
    """
    dans_grid3.lf.map_depressions(pits=None, reroute_flow=False)
    assert_array_equal(
        dans_grid3.mg.at_node["flow__receiver_node"], XX * np.ones(49, dtype=int)
    )
    assert_array_equal(
        dans_grid3.lf.depression_outlet_map, dans_grid3.depr_outlet_target
    )


def test_filling_supplied_pits(dans_grid3):
    """
    Test the filler without rereouting, but confusingly, where there *is*
    aready routing information available!
    Also tests the supply of an array for 'pits'
    """
    dans_grid3.fr.run_one_step()
    dans_grid3.lf.map_depressions(
        pits=dans_grid3.mg.at_node["flow__sink_flag"], reroute_flow=False
    )
    assert_array_equal(dans_grid3.mg.at_node["flow__receiver_node"], dans_grid3.r_old)


def test_pits_as_IDs(dans_grid3):
    """
    Smoke test for passing specific IDs, not an array, to the mapper.
    """
    dans_grid3.fr.run_one_step()
    dans_grid3.lf.map_depressions(
        pits=np.where(dans_grid3.mg.at_node["flow__sink_flag"])[0]
    )
    assert dans_grid3.mg.at_node["drainage_area"] == approx(dans_grid3.A_new)


def test_edge_draining():
    """
    This tests when the lake attempts to drain from an edge, where an issue
    is suspected.
    """
    # Create a 7x7 test grid with a well defined hole in it, AT THE EDGE.
    mg = RasterModelGrid((7, 7))

    z = mg.node_x.copy()
    guard_sides = np.concatenate((np.arange(7, 14), np.arange(35, 42)))
    edges = np.concatenate((np.arange(7), np.arange(42, 49)))
    hole_here = np.array(([15, 16, 22, 23, 29, 30]))
    z[guard_sides] = z[13]
    z[edges] = -2.0  # force flow outwards from the tops of the guards
    z[hole_here] = -1.0

    A_new = np.array(
        [
            [
                [
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    15.0,
                    5.0,
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                    0.0,
                    0.0,
                    10.0,
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    4.0,
                    3.0,
                    2.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                    0.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    1.0,
                    0.0,
                ]
            ]
        ]
    ).flatten()

    depr_outlet_target = np.array(
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            14,
            14,
            XX,
            XX,
            XX,
            XX,
            XX,
            14,
            14,
            XX,
            XX,
            XX,
            XX,
            XX,
            14,
            14,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ]
    ).flatten()

    mg.add_field("node", "topographic__elevation", z, units="-")

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)

    fr.run_one_step()
    lf.map_depressions()
    assert mg.at_node["drainage_area"] == approx(A_new)
    assert lf.depression_outlet_map == approx(depr_outlet_target)


def test_degenerate_drainage():
    """
    This "hourglass" configuration should be one of the hardest to correctly
    re-route.
    """
    mg = RasterModelGrid((9, 5))
    z_init = mg.node_x.copy() * 0.0001 + 1.0
    lake_pits = np.array([7, 11, 12, 13, 17, 27, 31, 32, 33, 37])
    z_init[lake_pits] = -1.0
    z_init[22] = 0.0  # the common spill pt for both lakes
    z_init[21] = 0.1  # an adverse bump in the spillway
    z_init[20] = -0.2  # the spillway
    mg.add_field("node", "topographic__elevation", z_init)

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    #    correct_A = np.array([ 0.,   0.,   0.,   0.,   0.,
    #                           0.,   1.,   3.,   1.,   0.,
    #                           0.,   5.,   1.,   2.,   0.,
    #                           0.,   1.,  10.,   1.,   0.,
    #                          21.,  21.,   1.,   1.,   0.,
    #                           0.,   1.,   9.,   1.,   0.,
    #                           0.,   3.,   1.,   2.,   0.,
    #                           0.,   1.,   1.,   1.,   0.,
    #                           0.,   0.,   0.,   0.,   0.])

    correct_A = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            3.0,
            1.0,
            0.0,
            0.0,
            2.0,
            4.0,
            2.0,
            0.0,
            0.0,
            1.0,
            10.0,
            1.0,
            0.0,
            21.0,
            21.0,
            1.0,
            1.0,
            0.0,
            0.0,
            1.0,
            9.0,
            1.0,
            0.0,
            0.0,
            2.0,
            2.0,
            2.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )

    assert mg.at_node["drainage_area"] == approx(correct_A)


def test_three_pits():
    """
    A test to ensure the component correctly handles cases where there are
    multiple pits.
    """
    mg = RasterModelGrid((10, 10))
    z = mg.add_field("node", "topographic__elevation", mg.node_x.copy())
    # a sloping plane
    # np.random.seed(seed=0)
    # z += np.random.rand(100)/10000.
    # punch some holes
    z[33] = 1.0
    z[43] = 1.0
    z[37] = 4.0
    z[74:76] = 1.0
    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node["flow__sink_flag"], flow_sinks_target)

    # test conservation of mass:
    assert mg.at_node["drainage_area"].reshape((10, 10))[1:-1, 1].sum() == approx(
        8.0 ** 2
    )
    # ^all the core nodes

    # test the actual flow field:
    nA = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            8.0,
            8.0,
            7.0,
            6.0,
            5.0,
            4.0,
            3.0,
            2.0,
            1.0,
            0.0,
            2.0,
            2.0,
            1.0,
            1.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            26.0,
            26.0,
            25.0,
            15.0,
            11.0,
            10.0,
            9.0,
            8.0,
            1.0,
            0.0,
            2.0,
            2.0,
            1.0,
            9.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            2.0,
            2.0,
            1.0,
            1.0,
            5.0,
            4.0,
            3.0,
            2.0,
            1.0,
            0.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            3.0,
            2.0,
            1.0,
            0.0,
            20.0,
            20.0,
            19.0,
            18.0,
            17.0,
            12.0,
            3.0,
            2.0,
            1.0,
            0.0,
            2.0,
            2.0,
            1.0,
            1.0,
            1.0,
            1.0,
            3.0,
            2.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert_array_equal(mg.at_node["drainage_area"], nA)

    # test a couple more properties:
    lc = np.empty(100, dtype=int)
    lc.fill(XX)
    lc[33] = 33
    lc[43] = 33
    lc[37] = 37
    lc[74:76] = 74
    assert_array_equal(lf.lake_map, lc)
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
    z = mg.add_field("node", "topographic__elevation", mg.node_x.copy())
    # a sloping plane
    # np.random.seed(seed=0)
    # z += np.random.rand(100)/10000.
    # punch one big hole
    z.reshape((10, 10))[3:8, 3:8] = 0.0
    # dig a couple of inset holes
    z[57] = -1.0
    z[44] = -2.0
    z[54] = -10.0

    # make an outlet
    z[71] = 0.9

    fr = FlowAccumulator(mg, flow_director="D8")
    lf = DepressionFinderAndRouter(mg)
    fr.run_one_step()
    lf.map_depressions()

    flow_sinks_target = np.zeros(100, dtype=bool)
    flow_sinks_target[mg.boundary_nodes] = True
    # no internal sinks now:
    assert_array_equal(mg.at_node["flow__sink_flag"], flow_sinks_target)

    # test conservation of mass:
    assert mg.at_node["drainage_area"].reshape((10, 10))[1:-1, 1].sum() == approx(
        8.0 ** 2
    )
    # ^all the core nodes

    # test the actual flow field:
    #    nA = np.array([  0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,
    #                     8.,   8.,   7.,   6.,   5.,   4.,   3.,   2.,   1.,   0.,
    #                     1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,
    #                     1.,   1.,   1.,   4.,   2.,   2.,   8.,   4.,   1.,   0.,
    #                     1.,   1.,   1.,   8.,   3.,  15.,   3.,   2.,   1.,   0.,
    #                     1.,   1.,   1.,  13.,  25.,   6.,   3.,   2.,   1.,   0.,
    #                     1.,   1.,   1.,  45.,   3.,   3.,   5.,   2.,   1.,   0.,
    #                    50.,  50.,  49.,   3.,   2.,   2.,   2.,   4.,   1.,   0.,
    #                     1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   1.,   0.,
    #                     0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.,   0.])
    nA = np.array(
        [
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            8.0,
            8.0,
            7.0,
            6.0,
            5.0,
            4.0,
            3.0,
            2.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            4.0,
            2.0,
            2.0,
            6.0,
            4.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            6.0,
            3.0,
            12.0,
            3.0,
            2.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            8.0,
            20.0,
            4.0,
            3.0,
            2.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            35.0,
            5.0,
            4.0,
            3.0,
            2.0,
            1.0,
            0.0,
            50.0,
            50.0,
            49.0,
            13.0,
            10.0,
            8.0,
            6.0,
            4.0,
            1.0,
            0.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    )
    assert_array_equal(mg.at_node["drainage_area"], nA)

    # the lake code map:
    lc = np.array(
        [
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            57,
            57,
            57,
            57,
            57,
            XX,
            XX,
            XX,
            XX,
            XX,
            57,
            57,
            57,
            57,
            57,
            XX,
            XX,
            XX,
            XX,
            XX,
            57,
            57,
            57,
            57,
            57,
            XX,
            XX,
            XX,
            XX,
            XX,
            57,
            57,
            57,
            57,
            57,
            XX,
            XX,
            XX,
            XX,
            XX,
            57,
            57,
            57,
            57,
            57,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
            XX,
        ]
    )

    # test the remaining properties:
    assert lf.lake_outlets.size == 1
    assert lf.lake_outlets[0] == 72
    outlets_in_map = np.unique(lf.depression_outlet_map)
    assert outlets_in_map.size == 2
    assert outlets_in_map[1] == 72
    assert lf.number_of_lakes == 1
    assert lf.lake_codes[0] == 57
    assert_array_equal(lf.lake_map, lc)
    assert lf.lake_areas[0] == approx(25.0)
    assert lf.lake_volumes[0] == approx(63.0)


def test_D8_D4_fill(d4_grid):
    """
    Tests the functionality of D4 filling.
    """
    d4_grid.lfD8.map_depressions(pits=None, reroute_flow=False)
    d4_grid.lfD4.map_depressions(pits=None, reroute_flow=False)
    assert d4_grid.lfD8.number_of_lakes == 1
    assert d4_grid.lfD4.number_of_lakes == 3

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

    assert_array_equal(d4_grid.lfD8.lake_map, correct_D8_lake_map)
    assert_array_equal(d4_grid.lfD4.lake_map, correct_D4_lake_map)

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

    #    flow_recD8 = np.array([ 0,  1,  2,  3,  4,  5,  6,  7, 16, 10, 16, 10, 18,
    #                           13, 14, 14, 15, 16, 10, 18, 20, 21, 16, 16, 16, 18,
    #                           33, 27, 28, 28, 24, 24, 24, 32, 34, 35, 35, 38, 32,
    #                           32, 32, 41, 42, 43, 44, 45, 46, 47, 48])
    flow_recD8 = np.array(
        [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            16,
            10,
            16,
            10,
            18,
            13,
            14,
            14,
            15,
            16,
            17,
            18,
            20,
            21,
            16,
            16,
            16,
            18,
            33,
            27,
            28,
            28,
            24,
            24,
            24,
            32,
            34,
            35,
            35,
            38,
            32,
            32,
            32,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
        ]
    )
    #    flow_recD4 = np.array([ 0,  1,  2,  3,  4,  5,  6,  7,  7, 10, 17, 10, 11,
    #                           13, 14, 14, 15, 16, 17, 18, 20, 21, 21, 16, 17, 18,
    #                           33, 27, 28, 28, 29, 24, 31, 32, 34, 35, 35, 36, 37,
    #                           32, 33, 41, 42, 43, 44, 45, 46, 47, 48])
    flow_recD4 = np.array(
        [
            0,
            1,
            2,
            3,
            4,
            5,
            6,
            7,
            7,
            10,
            17,
            10,
            11,
            13,
            14,
            14,
            15,
            16,
            17,
            18,
            20,
            21,
            21,
            16,
            17,
            18,
            33,
            27,
            28,
            28,
            29,
            38,
            31,
            32,
            34,
            35,
            35,
            36,
            37,
            32,
            33,
            41,
            42,
            43,
            44,
            45,
            46,
            47,
            48,
        ]
    )
    assert_array_equal(d4_grid.mg1.at_node["flow__receiver_node"], flow_recD8)
    assert_array_equal(d4_grid.mg2.at_node["flow__receiver_node"], flow_recD4)
    assert d4_grid.mg1.at_node["drainage_area"].reshape((7, 7))[:, 0].sum() == approx(
        d4_grid.mg2.at_node["drainage_area"].reshape((7, 7))[:, 0].sum()
    )
