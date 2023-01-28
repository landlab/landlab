#!/usr/bin/env python3
"""
Created on Fri May  1 14:02:18 2020

@author: gtucker

Adjusted on Thu Jan 19 2023
@author: bcampforts
"""

import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.components import ExponentialWeatherer


def test_create_weatherer_and_change_rate():
    grid = RasterModelGrid((3, 3), 1.0)
    grid.add_zeros("soil__depth", at="node")
    ew = ExponentialWeatherer(grid, soil_production_maximum_rate=0.0001)
    ew.maximum_weathering_rate = 0.0004
    np.testing.assert_array_equal(
        ew.maximum_weathering_rate, np.broadcast_to(0.0004, grid.number_of_nodes)
    )
    with pytest.raises(ValueError):
        ExponentialWeatherer(grid, soil_production_maximum_rate=-0.0001)
    with pytest.raises(ValueError):
        ExponentialWeatherer(grid, soil_production_decay_depth=-0.0001)


def test_run_weatherer():
    grid = RasterModelGrid((5, 5), 1.0)
    grid.add_ones("soil__depth", at="node")
    sp_max = np.random.rand(25)
    ew = ExponentialWeatherer(grid, soil_production_maximum_rate=sp_max)
    ew.run_one_step()

    new_value = np.random.rand(25)
    new_value[7] = -0.1
    with pytest.raises(ValueError):
        ew.maximum_weathering_rate = new_value
    with pytest.raises(ValueError):
        ew.decay_depth = new_value

    ew.maximum_weathering_rate = 0.0004
    np.testing.assert_array_equal(
        ew.maximum_weathering_rate, np.broadcast_to(0.0004, grid.number_of_nodes)
    )
    ew.decay_depth = 0.5
    np.testing.assert_array_equal(
        ew.decay_depth, np.broadcast_to(0.5, grid.number_of_nodes)
    )

    # Test analytical value
    ew.run_one_step()
    soil_production_rate = 0.0004 * np.exp(-grid.at_node["soil__depth"] / 0.5)
    grid.at_node["soil_production__rate"]
    np.testing.assert_array_equal(
        soil_production_rate[grid.core_nodes],
        grid.at_node["soil_production__rate"][grid.core_nodes],
    )
