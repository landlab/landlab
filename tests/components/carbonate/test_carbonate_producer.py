#! /usr/bin/env python
"""
Unit tests for landlab.components.carbonate.carbonate_producer
"""
from numpy.testing import assert_allclose
from numpy.testing import assert_raises

from landlab import RasterModelGrid
from landlab.components import CarbonateProducer


def test_properties():
    grid = RasterModelGrid((3, 3))
    elev = grid.add_zeros("topographic__elevation", at="node")
    grid.add_field("sea_level__elevation", 0.0, at="grid")
    elev[:] = -1.0
    cp = CarbonateProducer(grid)
    assert cp.extinction_coefficient == 0.1
    cp.extinction_coefficient = 0.04
    assert_allclose(cp.calc_carbonate_production_rate()[4], 0.009998656)
    assert cp.max_carbonate_production_rate == 0.01
    cp.max_carbonate_production_rate = 0.015
    assert_allclose(cp.calc_carbonate_production_rate()[4], 0.014997984)
    assert cp.surface_light == 2000.0
    cp.surface_light = 2250.0
    assert_allclose(cp.calc_carbonate_production_rate()[4], 0.014999393)
    assert cp.saturating_light == 400.0
    cp.saturating_light = 50.0
    assert_allclose(cp.calc_carbonate_production_rate()[4], 0.015)
    assert cp.tidal_range == 0.0
    cp.tidal_range = 10.0
    assert_allclose(cp.calc_carbonate_production_rate()[4], 0.00824751)


def test_preexisting_field():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_field("sea_level__elevation", 0.0, at="grid")
    grid.add_zeros("carbonate_thickness", at="node")
    cp = CarbonateProducer(grid)
    assert cp._carbonate_thickness is grid.at_node["carbonate_thickness"]


def test_exception_handling():
    grid = RasterModelGrid((3, 3))
    grid.add_zeros("topographic__elevation", at="node")
    grid.add_field("sea_level__elevation", 0.0, at="grid")
    assert_raises(ValueError, CarbonateProducer, grid, extinction_coefficient=0.0)
    assert_raises(ValueError, CarbonateProducer, grid, max_carbonate_production_rate=-1)
    assert_raises(ValueError, CarbonateProducer, grid, surface_light=-1)
    assert_raises(ValueError, CarbonateProducer, grid, saturating_light=-1)
    assert_raises(ValueError, CarbonateProducer, grid, tidal_range=-1)
