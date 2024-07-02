import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components import FlowAccumulator
from landlab.components import VariableAbrasionGBE


def test_sediment_thickness():
    grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    elev = grid.add_zeros("topographic__elevation", at="node")
    elev[3:] = 10.0
    sed = grid.add_zeros("soil__depth", at="node")
    sed[3:] = 10.0
    fa = FlowAccumulator(grid)
    fa.run_one_step()

    eroder = VariableAbrasionGBE(grid, abrasion_coefficients=[0.002, 0.0002, 0.00002])

    sum_of_classes = np.sum(eroder._thickness_by_class, axis=0)

    assert_array_equal(sum_of_classes, sed)

    eroder.run_one_step(1.0)

    sum_of_classes = np.sum(eroder._thickness_by_class, axis=0)

    assert_array_equal(sum_of_classes, sed)


def test_value_error():
    grid = RasterModelGrid((3, 3), xy_spacing=1000.0)
    with pytest.raises(ValueError):
        VariableAbrasionGBE(
            grid, number_of_sediment_classes=3, abrasion_coefficients=[0, 0]
        )
