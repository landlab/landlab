#!/usr/bin/env python2
"""
Unit tests for KinwaveImplicitOverlandFlowModel.

Created on Sat Apr  1 10:49:33 2017

@author: gtucker
"""

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import KinwaveImplicitOverlandFlow


def test_initialization():
    """Test initialization with various parameters."""
    rg = RasterModelGrid((3, 4), xy_spacing=2.0)
    rg.add_zeros("topographic__elevation", at="node")
    kw = KinwaveImplicitOverlandFlow(rg)

    # Make sure fields have been created
    for field_name in kw._info:
        if kw._info[field_name]["mapping"] == "node":
            assert field_name in kw.grid.at_node
        elif kw._info[field_name]["mapping"] == "link":
            assert field_name in kw.grid.at_link

    # Re-initialize, this time with fields already existing in the grid
    # (this triggers the "if" instead of "else" in the field setup in init)
    kw = KinwaveImplicitOverlandFlow(rg)


def test_zero_runoff_rate():
    grid = RasterModelGrid((5, 5))
    grid.add_field("topographic__elevation", 0.1 * grid.node_y, at="node")

    kw = KinwaveImplicitOverlandFlow(grid, runoff_rate=0.0)
    kw.run_one_step(1.0)

    assert np.all(grid.at_node["surface_water__depth"] == 0.0)


@pytest.mark.parametrize("var", ("runoff_rate", "roughness"))
@pytest.mark.parametrize("bad_value", (-1.0, [-1.0] * 25))
def test_negative_runoff_and_roughness(var, bad_value):
    grid = RasterModelGrid((5, 5))
    grid.add_field("topographic__elevation", 0.1 * grid.node_y, at="node")

    with pytest.raises(ValueError):
        KinwaveImplicitOverlandFlow(grid, **{var: bad_value})

    kw = KinwaveImplicitOverlandFlow(grid, **{var: 1.0})
    with pytest.raises(ValueError):
        setattr(kw, var, bad_value)


@pytest.mark.parametrize("var", ("runoff_rate", "roughness"))
@pytest.mark.parametrize("value", (2.0, [2.0] * 25))
def test_runoff_and_roughness_setter(var, value):
    grid = RasterModelGrid((5, 5))
    grid.add_field("topographic__elevation", 0.1 * grid.node_y, at="node")

    kw = KinwaveImplicitOverlandFlow(grid, **{var: 2.0})
    kw.run_one_step(1.0)
    expected = grid.at_node["surface_water__depth"].copy()

    grid.at_node["surface_water__depth"].fill(0.0)

    kw = KinwaveImplicitOverlandFlow(grid, **{var: 1.0})
    setattr(kw, var, value)
    kw.run_one_step(1.0)

    assert np.all(grid.at_node["surface_water__depth"] == expected)


@pytest.mark.parametrize("var", ("runoff_rate", "roughness"))
def test_runoff_rate_is_read_only(var):
    grid = RasterModelGrid((5, 5))
    grid.add_field("topographic__elevation", 0.1 * grid.node_y, at="node")

    kwargs = {var: np.full(grid.number_of_nodes, 2.0)}
    kw = KinwaveImplicitOverlandFlow(grid, **kwargs)
    assert (
        np.all(getattr(kw, var) == kwargs[var]) and getattr(kw, var) is not kwargs[var]
    )

    with pytest.raises(ValueError):
        getattr(kw, var)[:] = 1.0


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
    assert round(np.amax(kw._sqrt_slope), 3) == 0.316
    assert round(np.amax(kw._grad_width_sum), 3) == 0.632
    assert round(np.amax(kw._alpha), 3) == 15.811


def test_steady_basic_ramp():
    """Run to steady state with basic ramp"""

    # Create a basic ramp
    rg = RasterModelGrid((10, 10), xy_spacing=(2, 2))
    rg.add_field("topographic__elevation", 0.1 * rg.node_y, at="node")

    # Create component and run it
    kw = KinwaveImplicitOverlandFlow(rg, runoff_rate=0.001 * 3600000.0)
    for _ in range(12):
        kw.run_one_step(1.0)

    # Look at a column of nodes down the middle. The inflow from uphill should
    # be, from top to bottom: 0, 0.004, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028
    assert kw._disch_in[85] == 0.0
    assert round(kw._disch_in[75], 3) == 0.004
    assert round(kw._disch_in[65], 3) == 0.008
    assert round(kw._disch_in[55], 3) == 0.012
    assert round(kw._disch_in[45], 3) == 0.016
    assert round(kw._disch_in[35], 3) == 0.020
    assert round(kw._disch_in[25], 3) == 0.024
    assert round(kw._disch_in[15], 3) == 0.028

    # Try with passing in runoff
    kw = KinwaveImplicitOverlandFlow(rg, runoff_rate=360.0)
    kw.depth[:] = 0.0
    for _ in range(22):
        kw.run_one_step(1.0)

    # Again, look at a column of nodes down the middle. The inflow from uphill
    # should now be 1/10 of the prior example.
    assert round(kw._disch_in[75], 4) == 0.0004
    assert round(kw._disch_in[65], 4) == 0.0008
    assert round(kw._disch_in[55], 4) == 0.0012
    assert round(kw._disch_in[45], 4) == 0.0016
    assert round(kw._disch_in[35], 4) == 0.0020
    assert round(kw._disch_in[25], 4) == 0.0024
    assert round(kw._disch_in[15], 4) == 0.0028

    # Try with default runoff rate of 1 mm/hr
    kw = KinwaveImplicitOverlandFlow(rg)
    assert kw.runoff_rate == 1.0
    kw.depth[:] = 0.0
    for _ in range(18):
        kw.run_one_step(10.0)

    # Look at a column of nodes down the middle. The inflow from uphill should
    # be, from top to bottom: 0, 0.004, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028
    assert kw._disch_in[85] == 0.0
    assert round(kw._disch_in[75], 7) == 0.0000011
    assert round(kw._disch_in[65], 7) == 0.0000022
    assert round(kw._disch_in[55], 7) == 0.0000033
    assert round(kw._disch_in[45], 7) == 0.0000044
    assert round(kw._disch_in[35], 7) == 0.0000055
    assert round(kw._disch_in[25], 7) == 0.0000066
    assert round(kw._disch_in[15], 7) == 0.0000077


def test_curved_surface():
    """Test flow across a curved surface."""

    # Create a grid
    rg = RasterModelGrid((10, 10), xy_spacing=(2, 2))
    rg.add_field("topographic__elevation", 3.0 * rg.node_x**2 + rg.node_y**2, at="node")

    # Create component and run it
    kw = KinwaveImplicitOverlandFlow(rg, runoff_rate=0.001 * 3600000.0)
    for _ in range(8):
        kw.run_one_step(1.0)

    # The inflow discharge to each cell at steady state should equal the
    # runoff rate times the "inflow" drainage area, which is the total drainage
    # area minus the area of the cell itself. Here we'll test a column of core
    # nodes across the middle of the domain.
    area = rg.at_node["drainage_area"]
    runoff_rate = 0.001
    unit_area = 4.0
    for i in range(15, 95, 10):
        assert round(kw._disch_in[i], 6) == round(
            runoff_rate * (area[i] - unit_area), 6
        )


def test_kinwave_runoff_array():
    """
    Make sure that runoff_rate can be set with an array, and confirm that this
    returns the same result as setting with a float of the same magnitude.
    """
    # Set runoff_rate as float
    mg1 = RasterModelGrid((10, 10), xy_spacing=25)
    mg1.add_zeros("surface_water__depth", at="node")
    mg1.add_zeros("topographic__elevation", at="node")
    mg1.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r1 = 1.0
    kinwave1 = KinwaveImplicitOverlandFlow(
        mg1, runoff_rate=r1, roughness=0.03, depth_exp=5 / 3
    )

    # Set runoff_rate as array
    mg2 = RasterModelGrid((10, 10), xy_spacing=25)
    mg2.add_zeros("surface_water__depth", at="node")
    mg2.add_zeros("topographic__elevation", at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r2 = 1.0 * np.ones(100)
    kinwave2 = KinwaveImplicitOverlandFlow(
        mg2, runoff_rate=r2, roughness=0.03, depth_exp=5 / 3
    )

    kinwave1.run_one_step(100)
    kinwave2.run_one_step(100)
    np.testing.assert_equal(kinwave1.depth, kinwave2.depth)


def test_kinwave_roughness_array():
    """
    Make sure that roughness can be set with an array, and confirm that this
    returns the same result as setting with a float of the same magnitude.
    """

    mg1 = RasterModelGrid((10, 10), xy_spacing=25)
    mg1.add_zeros("surface_water__depth", at="node")
    mg1.add_zeros("topographic__elevation", at="node")
    mg1.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r1 = 0.03
    kinwave1 = KinwaveImplicitOverlandFlow(
        mg1, runoff_rate=1.0, roughness=r1, depth_exp=5 / 3
    )

    mg2 = RasterModelGrid((10, 10), xy_spacing=25)
    mg2.add_zeros("surface_water__depth", at="node")
    mg2.add_zeros("topographic__elevation", at="node")
    mg2.set_closed_boundaries_at_grid_edges(True, True, True, True)
    r2 = 0.03 * np.ones(100)
    kinwave2 = KinwaveImplicitOverlandFlow(
        mg2, runoff_rate=1.0, roughness=r2, depth_exp=5 / 3
    )

    kinwave1.run_one_step(100)
    kinwave2.run_one_step(100)
    np.testing.assert_equal(kinwave1.depth, kinwave2.depth)
