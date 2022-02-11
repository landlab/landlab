#! /usr/bin/env python
"""
Unit tests for landlab.components.carbonate.carbonate_producer
"""
import numpy as np
from numpy.testing import assert_allclose

from landlab import RasterModelGrid
from landlab.components import CarbonateProducer


def test_properties():
    grid = RasterModelGrid((3, 3))
    elev = grid.add_zeros('topographic__elevation', at='node')
    sealevel = grid.add_field('sea_level__elevation', 0.0, at='grid')
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
