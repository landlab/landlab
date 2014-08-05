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
        '_RMG': RasterModelGrid(4, 5, 12.)
    })


@with_setup(setup_unit_grid)
def test_unit_spacing():
    lengths = _RMG._calculate_link_length()
    assert_array_equal(lengths, np.ones(31))


@with_setup(setup_grid)
def test_non_unit_spacing():
    assert_array_equal(_RMG._calculate_link_length(),
                       _RMG.node_spacing * np.ones(31))


@with_setup(setup_grid)
def test_link_length():
    assert_array_equal(_RMG.link_length, 12 * np.ones(31))


@with_setup(setup_grid)
def test_active_link_length():
    assert_array_equal(_RMG.active_link_length, 12 * np.ones(17))
