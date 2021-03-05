#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  5 08:42:24 2021

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from landlab import RasterModelGrid
from landlab.components import ListricKinematicExtender, Flexure


def test_subsidence_and_horiz_shift():
    """Test that elev subsides then shifts after 2 time steps."""
    grid = RasterModelGrid((3, 7), xy_spacing=2500.0)
    topo = grid.add_zeros("topographic__elevation", at="node")
    extender = ListricKinematicExtender(
        grid, extension_rate=0.01, fault_location=2500.0
    )
    extender.run_one_step(dt=125000.0)
    assert_array_almost_equal(
        topo[7:14],
        [0.0, 0.0, -1404.156819, -910.66907, -590.616478, -383.045648, -248.425118],
    )


def test_with_flexure():
    """Test integrating with flexure."""
    grid = RasterModelGrid((3, 9), xy_spacing=1000.0)


if __name__ == "__main__":
    test_subsidence_and_horiz_shift()
