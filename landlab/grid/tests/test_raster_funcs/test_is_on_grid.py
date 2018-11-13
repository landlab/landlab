import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.grid import raster_funcs as rfuncs


def test_with_arrays():
    """Test with arrays as arg."""
    rmg = RasterModelGrid((4, 5), spacing=(2., 2.))

    coords = (np.array([1., -1.]), np.array([1., -1.]))
    assert_array_equal(rfuncs.is_coord_on_grid(rmg, coords), np.array([True, False]))


def test_just_inside():
    """Test with points just inside the grid."""
    rmg = RasterModelGrid((4, 5), spacing=(2., 2.))

    assert rfuncs.is_coord_on_grid(rmg, (0., 4.))
    assert rfuncs.is_coord_on_grid(rmg, (8. - 1e-12, 4.))
    assert rfuncs.is_coord_on_grid(rmg, (3., 0.))
    assert rfuncs.is_coord_on_grid(rmg, (3., 6. - 1e-12))


def test_just_outside():
    """Test with points just outside the grid."""
    rmg = RasterModelGrid((4, 5), spacing=(2., 2.))

    assert not rfuncs.is_coord_on_grid(rmg, (0. - 1e-12, 4.))
    assert not rfuncs.is_coord_on_grid(rmg, (8., 4.))
    assert not rfuncs.is_coord_on_grid(rmg, (3., 0. - 1e-12))
    assert not rfuncs.is_coord_on_grid(rmg, (3., 6.))


def test_just_x():
    """Test check if points are within the x bounds."""
    rmg = RasterModelGrid((4, 5), spacing=(2., 2.))
    assert rfuncs.is_coord_on_grid(rmg, (4., 1.e6), axes=(1,))
    assert not rfuncs.is_coord_on_grid(rmg, (-1., 1.), axes=(1,))
