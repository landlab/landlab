import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_number_of_corners():
    """Test number of corners on a raster."""
    grid = RasterModelGrid((4, 5))
    assert grid.number_of_corners == 12


def test_add_ones_at_corners():
    """Test add a field to corners."""
    grid = RasterModelGrid((4, 5))
    grid.add_ones("z", at="corner")

    assert_array_equal(
        grid.at_corner["z"],
        [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
    )


def test_add_field_at_corners():
    """Test add a field to corners."""
    grid = RasterModelGrid((4, 5))
    x = np.arange(grid.number_of_corners)
    grid.add_field("z", x, at="corner")

    assert_array_equal(grid.at_corner["z"], x)
