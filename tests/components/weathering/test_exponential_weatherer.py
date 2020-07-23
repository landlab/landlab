#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May  1 14:02:18 2020

@author: gtucker
"""

from landlab import RasterModelGrid
from landlab.components import ExponentialWeatherer


def test_create_weatherer_and_change_rate():

    grid = RasterModelGrid((3, 3), 1.0)
    grid.add_zeros("soil__depth", at="node")

    ew = ExponentialWeatherer(grid, soil_production__maximum_rate=0.0001)
    ew.maximum_weathering_rate = 0.0004
    assert ew.maximum_weathering_rate == 0.0004
