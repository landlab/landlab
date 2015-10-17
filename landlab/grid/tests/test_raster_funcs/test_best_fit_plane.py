import numpy as np
from numpy.testing import assert_array_equal
from nose.tools import assert_equal, assert_raises

from landlab.grid import raster_funcs as rfuncs


def test_best_fit_with_bad_args():
    """Raise an error if input arrays are not the same length."""
    len_2 = np.array([0, 1])
    len_3 = np.array([0, 0, 1])

    assert_raises(ValueError,
                  rfuncs.calculate_slope_aspect_bfp, len_2, len_3, len_3)
    assert_raises(ValueError,
                  rfuncs.calculate_slope_aspect_bfp, len_3, len_2, len_3)
    assert_raises(ValueError,
                  rfuncs.calculate_slope_aspect_bfp, len_3, len_3, len_2)


def test_best_fit_plane_in_xy():
    """Best fit plane is the xy-plane."""
    (x, y, z) = (np.array([0., 1., 1.]),
                 np.array([0., 0., 1.]),
                 np.array([0., 0., 0.]), )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert_equal(slope, 0.)
    assert_equal(aspect, 90.)

    (x, y, z) = (np.array([0., 1., 1., 0.]),
                 np.array([0., 0., 1., 1.]),
                 np.array([0., 0., 0., 0.]), )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert_equal(slope, 0.)
    assert_equal(aspect, 90.)


def test_best_fit_plane_in_xz():
    """Best fit plane is the xz-plane."""
    (x, y, z) = (np.array([0., 1., 1.]),
                 np.array([0., 0., 0.]),
                 np.array([0., 0., 1.]), )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert_equal(slope, 90.)
    assert_equal(aspect, 0.)


def test_best_fit_plane_in_yz():
    """Best fit plane is the yz-plane."""
    (x, y, z) = (np.array([0., 0., 0.]),
                 np.array([0., 1., 1.]),
                 np.array([0., 0., 1.]), )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert_equal(slope, 90.)
    assert_equal(aspect, 90.)
