import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid import raster_funcs as rfuncs
from landlab.testing.tools import assert_array_is_int


def test_with_scalars():
    """Test scalar args."""
    rmg = RasterModelGrid((4, 5))
    id = rfuncs.find_nearest_node(rmg, (0.2, 0.6))
    assert id == 5
    assert id.ndim == 0
    assert_array_is_int(id)


def test_with_iterable():
    """Test iterable args."""
    rmg = RasterModelGrid((4, 5))
    id = rfuncs.find_nearest_node(rmg, ([0.2], [0.6]))
    assert_array_equal(id, np.array([5], dtype=int))
    assert_array_is_int(id)


def test_with_ndarray_with_length_0():
    """Test with 0d numpy arrays as args."""
    rmg = RasterModelGrid((4, 5))
    id = rfuncs.find_nearest_node(rmg, (np.array(0.2), np.array(0.6)))
    assert_array_equal(id, np.array(5, dtype=int))
    assert id.ndim == 0
    assert_array_is_int(id)


def test_with_ndarray():
    """Test with 1d numpy arrays as args."""
    rmg = RasterModelGrid((4, 5))
    coords = (np.array([0.1, 0.2]), np.array([3.4, 2.6]))
    id = rfuncs.find_nearest_node(rmg, coords)
    assert_array_equal(id, np.array([15, 15], dtype=int))
    assert isinstance(id, np.ndarray)


def test_non_unit_spacing():
    """Test with a grid of non-unit spacing."""
    rmg = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))
    id = rfuncs.find_nearest_node(rmg, (0.9, 0.2))
    assert id == 0


def test_beyond_grid():
    """Raise an error if points are outside the bounds of the grid."""
    rmg = RasterModelGrid((4, 5), xy_spacing=(2.0, 2.0))

    assert rfuncs.find_nearest_node(rmg, (-0.999, 0.2)) == 0
    with pytest.raises(ValueError):
        rfuncs.find_nearest_node(rmg, (-1.001, 0.2))

    assert rfuncs.find_nearest_node(rmg, (8.999, 0.2)) == 4
    with pytest.raises(ValueError):
        rfuncs.find_nearest_node(rmg, (9.001, 0.2))

    assert rfuncs.find_nearest_node(rmg, (0.2, -0.999)) == 0
    with pytest.raises(ValueError):
        rfuncs.find_nearest_node(rmg, (0.2, -1.001))

    assert rfuncs.find_nearest_node(rmg, (0.2, 6.999)) == 15
    with pytest.raises(ValueError):
        rfuncs.find_nearest_node(rmg, (0.2, 7.001))
