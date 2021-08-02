import numpy as np
from numpy.testing import assert_array_equal

from landlab import VoronoiDelaunayGrid
from landlab.grid import grid_funcs as gfuncs

v_x = np.array(
    [
        0.789,
        0.897,
        0.605,
        0.533,
        0.661,
        0.889,
        0.905,
        0.885,
        0.996,
        0.506,
        0.363,
        0.429,
        0.501,
        0.812,
        0.225,
        0.964,
        0.808,
        0.095,
        0.051,
        0.126,
    ]
)
v_y = np.array(
    [
        0.432,
        0.737,
        0.077,
        0.164,
        0.597,
        0.868,
        0.283,
        0.875,
        0.126,
        0.247,
        0.544,
        0.913,
        0.648,
        0.547,
        0.576,
        0.053,
        0.951,
        0.846,
        0.980,
        0.549,
    ]
)


def test_with_scalars_float():
    """Test scalar args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    id = gfuncs.find_nearest_node(vmg, (0.5, 0.6))
    assert_array_equal(id, np.array([12], dtype=int))
    assert id.ndim == 1  # not same as raster version
    assert id.dtype == np.int32 or id.dtype == np.int64


def test_with_scalars_int():
    """Test scalar args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    id = gfuncs.find_nearest_node(vmg, (0, 0))
    assert_array_equal(id, np.array([3], dtype=int))
    assert id.ndim == 1  # not same as raster version
    assert id.dtype == np.int32 or id.dtype == np.int64


def test_with_iterable():
    """Test iterable args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    id = gfuncs.find_nearest_node(vmg, ([0.5], [0.6]))
    assert_array_equal(id, np.array([12], dtype=int))
    assert id.dtype == np.int32 or id.dtype == np.int64


def test_with_ndarray_with_length_0():
    """Test with 0d numpy arrays as args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    id = gfuncs.find_nearest_node(vmg, (np.array(0.5), np.array(0.6)))
    assert_array_equal(id, np.array([12], dtype=int))
    assert id.ndim == 1  # again, differs from raster version
    assert id.dtype == np.int32 or id.dtype == np.int64


def test_with_ndarray_float():
    """Test with 1d numpy float arrays as args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    coords = (np.array([0.45, 0.5]), np.array([0.65, 0.6]))
    id = gfuncs.find_nearest_node(vmg, coords)
    assert_array_equal(id, np.array([12, 12], dtype=int))
    assert isinstance(id, np.ndarray)


def test_with_ndarray_int():
    """Test with 1d numpy int arrays as args."""
    vmg = VoronoiDelaunayGrid(v_x, v_y)
    coords = (np.array([0], dtype=int), np.array([0], dtype=int))
    id = gfuncs.find_nearest_node(vmg, coords)
    assert_array_equal(id, np.array([3], dtype=int))
    assert isinstance(id, np.ndarray)
