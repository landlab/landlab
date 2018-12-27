import numpy as np
import pytest

from landlab import RasterModelGrid
from landlab.utils.decorators import use_field_name_array_or_value


@use_field_name_array_or_value("cell")
def my_func(grid, vals):
    return grid.area_of_cell * vals


def test_use_field_name_array_or_value_raises_errors():
    grid = RasterModelGrid((4, 5), xy_spacing=(1, 2))
    bad_values = np.array([[0, 1, 2, 3, 4]])
    with pytest.raises(ValueError):
        my_func(grid, bad_values)
