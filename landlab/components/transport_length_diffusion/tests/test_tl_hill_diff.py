#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 16:52:10 2017

@author: margauxmouchene
"""

import numpy as np
import pytest
from numpy.testing import assert_almost_equal

from landlab import RasterModelGrid
from landlab.components import (
    FlowAccumulator,
    FlowDirectorSteepest,
    TransportLengthHillslopeDiffuser,
)


def test_route_to_multiple_error_raised():
    mg = RasterModelGrid((10, 10))
    z = mg.add_zeros("node", "topographic__elevation")
    z += mg.x_of_node + mg.y_of_node
    fa = FlowAccumulator(mg, flow_director="MFD")
    fa.run_one_step()

    with pytest.raises(NotImplementedError):
        TransportLengthHillslopeDiffuser(mg, erodibility=1.0, slope_crit=0.5)


def test_tl_hill_diff():
    """Test cases where S>Sc, S=Sc and S<Sc"""

    # Test cases where S>Sc, S=Sc and S<Sc
    # Set up a 3x16 grid with closed boundaries and initial elevations.
    mg = RasterModelGrid((3, 12))
    z = np.array(
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
            5.0,
            1.9,
            1.9,
            1.9,
            1.9,
            1.3,
            1.3,
            1.3,
            1.3,
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
            0.0,
            0.0,
        ]
    )
    mg.add_field("node", "topographic__elevation", z)
    mg.set_closed_boundaries_at_grid_edges(True, True, True, True)

    # Parameter values for test
    k = 0.001
    Sc = 0.6

    # Instantiate flow director and tl hillslope diffuser
    fdir = FlowDirectorSteepest(mg)
    tl_diff = TransportLengthHillslopeDiffuser(mg, erodibility=k, slope_crit=Sc)

    # Run flow director
    fdir.run_one_step()

    # test slopes
    s_out = mg.at_node["topographic__steepest_slope"]
    s_test = np.array(
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
            3.1,
            0.0,
            0.0,
            0.0,
            0.6,
            0.0,
            0.0,
            0.0,
            0.3,
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
        ]
    )
    assert_almost_equal(s_out, s_test, decimal=10)

    # Run tl hillslope diffusion component
    tl_diff.run_one_step(1.0)

    # Test results
    # flux_out
    fo_test = np.array(
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
            0.025,
            0.0,
            0.0,
            0.0,
            0.0006,
            0.0,
            0.0,
            0.0,
            0.0003,
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
        ]
    )
    fo_out = mg.at_node["sediment__flux_out"]
    assert_almost_equal(fo_out, fo_test, decimal=10)

    # updated elevation
    elev_test = np.array(
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
            4.975,
            1.9,
            1.9,
            1.9,
            1.8994,
            1.3,
            1.3,
            1.3,
            1.2997,
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
            0.0,
            0.0,
        ]
    )
    elev_out = mg.at_node["topographic__elevation"]
    assert_almost_equal(elev_out, elev_test, decimal=10)

    # Run another time step because deposition and transfer were null
    # the first time
    fdir.run_one_step()
    tl_diff.run_one_step(1.0)

    # Test results
    # flux_out
    fo_test = np.array(
        [
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            2.47500000e-02,
            0.00000000e00,
            0.00000000e00,
            6.00000000e-07,
            5.99400000e-04,
            0.00000000e00,
            0.00000000e00,
            3.00000000e-07,
            2.99700000e-04,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
            0.00000000e00,
        ]
    )
    fo_out = mg.at_node["sediment__flux_out"]
    assert_almost_equal(fo_out, fo_test, decimal=10)

    # updated elevation
    elev_test = np.array(
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
            4.95025,
            1.925,
            1.9,
            1.8999994,
            1.8988006,
            1.3006,
            1.3,
            1.2999997,
            1.2994003,
            1.0003,
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
        ]
    )
    elev_out = mg.at_node["topographic__elevation"]
    assert_almost_equal(elev_out, elev_test, decimal=10)
