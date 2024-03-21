#!/usr/bin/env python2
"""
Unit tests for KinwaveImplicitOverlandFlowModel.

Created on Sat Apr  1 10:49:33 2017

@author: gtucker
"""

from numpy.testing import assert_allclose
from numpy.testing import assert_raises

from landlab import RasterModelGrid
from landlab.components import LinearDiffusionOverlandFlowRouter


def test_steady_one_node():
    """Run to steady state with a single node"""
    xy_spacing = 2.0
    roughness = 0.01
    velocity_scale = 1.0
    rain_rate = 72.0 / (3600.0 * 1000.0)

    grid = RasterModelGrid((3, 3), xy_spacing=xy_spacing)
    grid.set_closed_boundaries_at_grid_edges(False, True, True, True)
    grid.add_zeros("topographic__elevation", at="node")
    olflow = LinearDiffusionOverlandFlowRouter(
        grid, rain_rate=rain_rate, roughness=roughness, velocity_scale=velocity_scale
    )
    for _ in range(18):
        olflow.run_one_step(20.0)
    actual = grid.at_node["surface_water__depth"][4]
    expected = (grid.dx * grid.dy * olflow.rain_rate / olflow.vel_coef) ** 0.3

    assert_allclose(actual, expected, atol=1e-4)


def test_invalid_parameters():
    """Make sure invalid parameters raise the right exceptions."""
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    assert_raises(ValueError, LinearDiffusionOverlandFlowRouter, grid, roughness=-1.0)
    assert_raises(ValueError, LinearDiffusionOverlandFlowRouter, grid, rain_rate=-1.0)
    assert_raises(
        ValueError, LinearDiffusionOverlandFlowRouter, grid, infilt_depth_scale=0.0
    )
    assert_raises(ValueError, LinearDiffusionOverlandFlowRouter, grid, infilt_rate=-1.0)
    assert_raises(
        ValueError, LinearDiffusionOverlandFlowRouter, grid, velocity_scale=0.0
    )
    assert_raises(ValueError, LinearDiffusionOverlandFlowRouter, grid, cfl_factor=0.0)
    assert_raises(ValueError, LinearDiffusionOverlandFlowRouter, grid, cfl_factor=2.0)
    olflow = LinearDiffusionOverlandFlowRouter(grid)
    try:
        olflow.rain_rate = -1.0
    except ValueError:
        pass
    else:
        raise AssertionError("Setting negative rain_rate should raise ValueError")
