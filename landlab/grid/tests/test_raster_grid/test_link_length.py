import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import with_setup

from landlab import RasterModelGrid


def setup_unit_grid():
    globals().update({
        '_RMG': RasterModelGrid(4, 5)
    })


def setup_grid():
    globals().update({
        '_RMG': RasterModelGrid((4, 5), spacing=(3, 4))
    })


@with_setup(setup_unit_grid)
def test_unit_spacing():
    lengths = _RMG._create_length_of_link()
    assert_array_equal(lengths, np.ones(31))


@with_setup(setup_grid)
def test_non_unit_spacing():
    assert_array_equal(_RMG._create_length_of_link(),
                       [4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.])


@with_setup(setup_grid)
def test_link_length():
    assert_array_equal(_RMG.length_of_link,
                       [4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3., 3., 3.,
                        4., 4., 4., 4.])


@with_setup(setup_grid)
def test_active_link_length():
    assert_array_equal(_RMG.length_of_link[_RMG.active_links],
                       [3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3.,
                        4., 4., 4., 4.,
                        3., 3., 3.])
