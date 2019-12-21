import numpy as np
import pytest

from landlab.grid import raster_funcs as rfuncs


def test_best_fit_with_bad_args():
    """Raise an error if input arrays are not the same length."""
    len_2 = np.array([0, 1])
    len_3 = np.array([0, 0, 1])

    with pytest.raises(ValueError):
        rfuncs.calculate_slope_aspect_bfp(len_2, len_3, len_3)
    with pytest.raises(ValueError):
        rfuncs.calculate_slope_aspect_bfp(len_3, len_2, len_3)
    with pytest.raises(ValueError):
        rfuncs.calculate_slope_aspect_bfp(len_3, len_3, len_2)


def test_best_fit_plane_in_xy():
    """Best fit plane is the xy-plane."""
    (x, y, z) = (
        np.array([0.0, 1.0, 1.0]),
        np.array([0.0, 0.0, 1.0]),
        np.array([0.0, 0.0, 0.0]),
    )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert slope == 0.0
    assert aspect == 90.0

    (x, y, z) = (
        np.array([0.0, 1.0, 1.0, 0.0]),
        np.array([0.0, 0.0, 1.0, 1.0]),
        np.array([0.0, 0.0, 0.0, 0.0]),
    )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert slope == 0.0
    assert aspect == 90.0


def test_best_fit_plane_in_xz():
    """Best fit plane is the xz-plane."""
    (x, y, z) = (
        np.array([0.0, 1.0, 1.0]),
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert slope == 90.0
    assert aspect == 0.0


def test_best_fit_plane_in_yz():
    """Best fit plane is the yz-plane."""
    (x, y, z) = (
        np.array([0.0, 0.0, 0.0]),
        np.array([0.0, 1.0, 1.0]),
        np.array([0.0, 0.0, 1.0]),
    )
    (slope, aspect) = rfuncs.calculate_slope_aspect_bfp(x, y, z)

    assert slope == 90.0
    assert aspect == 90.0
