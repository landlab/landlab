"""
Unit tests for landlab.components.dimensionless.discharge
"""

import numpy as np
from numpy.testing import assert_array_almost_equal

from landlab import RasterModelGrid
from landlab.components.dimensionless_discharge import DimensionlessDischarge


def init_grid():
    """initialize grid for testing"""
    watershed_grid = RasterModelGrid((3, 3))
    _ = watershed_grid.add_ones("node", "flux")
    _ = watershed_grid.add_ones("node", "d50")
    _ = watershed_grid.add_ones("node", "dem_values")
    watershed_grid.at_node["dem_values"] = np.array([[1.1, 2, 3, 4, 2, 3, 4,
                                                      5, 3]])
    return watershed_grid


def test_dimensionless_discharge_final_values():
    """Testing for correct final dimensionless discharge values"""
    watershed_grid = init_grid()
    dd = DimensionlessDischarge(watershed_grid)
    dd.run_one_step()
    expected_values = [
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
        0.55372743,
    ]
    assert_array_almost_equal(
        watershed_grid.at_node["dimensionless_discharge"], expected_values
    )


def test_dimensionless_discharge_threshold_values():
    """Validate threshold values are correctly calculated"""
    watershed_grid = init_grid()
    dd = DimensionlessDischarge(watershed_grid)
    dd.run_one_step()
    expected_values = [
        11.09442633,
        12.63600546,
        14.73515969,
        11.01921767,
        11.72461201,
        12.5362998,
        10.94510988,
        10.94510988,
        10.94510988,
    ]
    assert_array_almost_equal(
        watershed_grid.at_node["dimensionless_discharge_threshold_value"],
        expected_values,
    )
