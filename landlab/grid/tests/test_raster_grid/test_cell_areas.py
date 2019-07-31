import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_unit_grid_all_cells():
    rmg = RasterModelGrid((4, 4))
    assert_array_equal(rmg.area_of_cell, np.ones(4))


def test_unit_grid_one_cell():
    rmg = RasterModelGrid((4, 4))
    assert_array_equal(rmg.area_of_cell[0], 1.0)


def test_unit_grid_last_cell():
    rmg = RasterModelGrid((4, 4))
    assert_array_equal(rmg.area_of_cell[-1], 1.0)


def test_out_of_range():
    rmg = RasterModelGrid((4, 4))
    with pytest.raises(IndexError):
        rmg.area_of_cell[5]


def test_is_immutable():
    rmg = RasterModelGrid((4, 4))
    with pytest.raises(ValueError):
        rmg.area_of_cell[0] = 0.0


def test_all_cells_with_spacing():
    rmg = RasterModelGrid((4, 4), xy_spacing=10.0)
    assert_array_equal(rmg.area_of_cell, 100 * np.ones(4))
