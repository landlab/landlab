#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Unit tests for KinwaveImplicitOverlandFlowModel.

Created on Sat Apr  1 10:49:33 2017

@author: gtucker
"""

import numpy as np

from landlab import RasterModelGrid
from landlab.components import KinwaveImplicitOverlandFlow


def test_initialization():
    """Test initialization with various parameters.
    """
    rg = RasterModelGrid((3, 4), xy_spacing=2.0)
    rg.add_zeros("node", "topographic__elevation")
    kw = KinwaveImplicitOverlandFlow(rg)

    # Make sure fields have been created
    for field_name in kw._var_mapping:
        if kw._var_mapping[field_name] == "node":
            assert field_name in kw.grid.at_node
        elif kw._var_mapping[field_name] == "link":
            assert field_name in kw.grid.at_link

    # Re-initialize, this time with fields already existing in the grid
    # (this triggers the "if" instead of "else" in the field setup in init)
    kw = KinwaveImplicitOverlandFlow(rg)


def test_first_iteration():
    """Test stuff that happens only on first iteration"""

    # Create a basic ramp
    rg = RasterModelGrid((10, 10), xy_spacing=(2, 2))
    rg.add_field("topographic__elevation", 0.1 * rg.node_y, at="node")

    # Create component and run it
    kw = KinwaveImplicitOverlandFlow(rg)
    kw.run_one_step(1.0)

    # Max gradient should be 0.1, and min should be zero
    assert round(np.amax(kw.grid.at_link["topographic__gradient"]), 2) == 0.1
    assert round(np.amin(kw.grid.at_link["topographic__gradient"]), 2) == 0.0
    assert round(np.amax(kw.sqrt_slope), 3) == 0.316
    assert round(np.amax(kw.grad_width_sum), 3) == 0.632
    assert round(np.amax(kw.alpha), 3) == 15.811


def test_steady_basic_ramp():
    """Run to steady state with basic ramp"""

    # Create a basic ramp
    rg = RasterModelGrid((10, 10), xy_spacing=(2, 2))
    rg.add_field("topographic__elevation", 0.1 * rg.node_y, at="node")

    # Create component and run it
    kw = KinwaveImplicitOverlandFlow(rg)
    for i in range(12):
        kw.run_one_step(1.0, runoff_rate=0.001)

    # Look at a column of nodes down the middle. The inflow from uphill should
    # be, from top to bottom: 0, 0.004, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028
    assert kw.disch_in[85] == 0.0
    assert round(kw.disch_in[75], 3) == 0.004
    assert round(kw.disch_in[65], 3) == 0.008
    assert round(kw.disch_in[55], 3) == 0.012
    assert round(kw.disch_in[45], 3) == 0.016
    assert round(kw.disch_in[35], 3) == 0.020
    assert round(kw.disch_in[25], 3) == 0.024
    assert round(kw.disch_in[15], 3) == 0.028

    # Try with passing in runoff
    kw = KinwaveImplicitOverlandFlow(rg, runoff_rate=360.0)
    kw.depth[:] = 0.0
    for i in range(22):
        kw.run_one_step(1.0)

    # Again, look at a column of nodes down the middle. The inflow from uphill
    # should now be 1/10 of the prior example.
    assert round(kw.disch_in[75], 4) == 0.0004
    assert round(kw.disch_in[65], 4) == 0.0008
    assert round(kw.disch_in[55], 4) == 0.0012
    assert round(kw.disch_in[45], 4) == 0.0016
    assert round(kw.disch_in[35], 4) == 0.0020
    assert round(kw.disch_in[25], 4) == 0.0024
    assert round(kw.disch_in[15], 4) == 0.0028

    # Try with default runoff rate of 1 mm/hr = 2.78e-7 m/s
    kw = KinwaveImplicitOverlandFlow(rg)
    assert round(kw.runoff_rate * 1.0e7, 2) == 2.78
    kw.depth[:] = 0.0
    for i in range(18):
        kw.run_one_step(10.0)

    # Look at a column of nodes down the middle. The inflow from uphill should
    # be, from top to bottom: 0, 0.004, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028
    assert kw.disch_in[85] == 0.0
    assert round(kw.disch_in[75], 7) == 0.0000011
    assert round(kw.disch_in[65], 7) == 0.0000022
    assert round(kw.disch_in[55], 7) == 0.0000033
    assert round(kw.disch_in[45], 7) == 0.0000044
    assert round(kw.disch_in[35], 7) == 0.0000055
    assert round(kw.disch_in[25], 7) == 0.0000066
    assert round(kw.disch_in[15], 7) == 0.0000077


def test_curved_surface():
    """Test flow across a curved surface."""

    # Create a grid
    rg = RasterModelGrid((10, 10), xy_spacing=(2, 2))
    rg.add_field(
        "topographic__elevation", 3.0 * rg.node_x ** 2 + rg.node_y ** 2, at="node"
    )

    # Create component and run it
    kw = KinwaveImplicitOverlandFlow(rg)
    for i in range(8):
        kw.run_one_step(1.0, runoff_rate=0.001)

    # The inflow discharge to each cell at steady state should equal the
    # runoff rate times the "inflow" drainage area, which is the total drainage
    # area minus the area of the cell itself. Here we'll test a column of core
    # nodes across the middle of the domain.
    area = rg.at_node["drainage_area"]
    runoff_rate = 0.001
    unit_area = 4.0
    for i in range(15, 95, 10):
        assert round(kw.disch_in[i], 6) == round(runoff_rate * (area[i] - unit_area), 6)


if __name__ == "__main__":
    test_initialization()
    test_first_iteration()
    test_steady_basic_ramp()
    test_curved_surface()
