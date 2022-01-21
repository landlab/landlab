#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Unit tests for KinwaveImplicitOverlandFlowModel.

Created on Sat Apr  1 10:49:33 2017

@author: gtucker
"""

import numpy as np
from numpy.testing import assert_equal

from landlab import RasterModelGrid
from landlab.components import LinearDiffusionOverlandFlowRouter


def test_steady_one_node():
    """Run to steady state with a single node"""
    grid = RasterModelGrid((3, 3), xy_spacing=2.0)
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)
    elev = grid.add_zeros("topographic__elevation", at="node")
    olflow = LinearDiffusionOverlandFlowRouter(
        grid, rain_rate=72.0 / (3600.0 * 1000.0), roughness=0.01, velocity_scale=1.0
    )
    for _ in range(18):
        olflow.run_one_step(20.0)
    assert_equal(round(grid.at_node["surface_water__depth"][4], 4), 0.0037)
