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
    watershed_grid.add_ones("surface_water__unit_discharge", at="node")
    watershed_grid.add_ones("channel_bottom_sediment_grain__d50_diameter", at="node")
    watershed_grid.add_ones("topographic__elevation", at="node")
    watershed_grid.at_node["topographic__elevation"] = np.array(
        [[1.1, 2, 3, 4, 2, 3, 4, 5, 3]]
    )
    return watershed_grid


def test_dimensionless_discharge_final_values():
    """Testing for correct final dimensionless discharge values"""
    watershed_grid = init_grid()
    # Setting the gravitation constant to 9.8 to make hand calculating
    # expeted values easier.
    dd = DimensionlessDischarge(watershed_grid, gravity=9.8)
    dd.run_one_step()
    # 0.55372743 is the expected value from the dimensionless discharge
    # (q*) equation in the Tang et al paper when surface_water__unit_discharge = 1
    # and channel_bottom_sediment_grain__d50_diameter = 1.
    # All values in dimensionless_discharge should be 0.5537274 after
    # running the un_one_step() function since the D50 and
    # surface_water__unit_discharge values are all 1s
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
    # Setting the gravitation constant to 9.8 to make hand calculating
    # expeted values easier.
    dd = DimensionlessDischarge(watershed_grid, gravity=9.8)
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
        watershed_grid.at_node["dimensionless_discharge_threshold"],
        expected_values,
    )
