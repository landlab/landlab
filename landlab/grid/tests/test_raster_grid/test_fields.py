import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_raises

from landlab import RasterModelGrid


def test_add_field_at_node():
    """Add field at nodes."""
    grid = RasterModelGrid((4, 5))
    grid.add_field('z', np.arange(20), at='node')

    assert_array_equal(grid.at_node['z'], np.arange(20))


def test_add_field_without_at_keyword():
    """Test default is at nodes."""
    grid = RasterModelGrid((4, 5))
    grid.add_field('z', np.arange(20))

    assert_array_equal(grid.at_node['z'], np.arange(20))


def test_add_field_without_at():
    """Test raises error with wrong size array."""
    grid = RasterModelGrid((4, 5))
    assert_raises(ValueError, grid.add_field, 'z', np.arange(21))
