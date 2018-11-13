import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid


def test_unit_spacing():
    grid = RasterModelGrid((4, 5))
    assert_array_equal(grid.length_of_link, np.ones(31))


def test_non_unit_spacing():
    grid = RasterModelGrid((4, 5), spacing=(3, 4))
    assert_array_equal(
        grid.length_of_link,
        [
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
        ],
    )


def test_link_length():
    grid = RasterModelGrid((4, 5), spacing=(3, 4))
    assert_array_equal(
        grid.length_of_link,
        [
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
            3.,
            3.,
            3.,
            3.,
            3.,
            4.,
            4.,
            4.,
            4.,
        ],
    )


def test_active_link_length():
    grid = RasterModelGrid((4, 5), spacing=(3, 4))
    assert_array_equal(
        grid.length_of_link[grid.active_links],
        [3., 3., 3., 4., 4., 4., 4., 3., 3., 3., 4., 4., 4., 4., 3., 3., 3.],
    )
