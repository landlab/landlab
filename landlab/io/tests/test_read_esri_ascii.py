#! /usr/bin/env python
"""
Unit tests for landlab.io.esri_ascii module.
"""
import os

import numpy as np
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal
from six import StringIO

from landlab import RasterModelGrid
from landlab.io import (
    BadHeaderLineError,
    DataSizeError,
    KeyTypeError,
    KeyValueError,
    MismatchGridDataSizeError,
    MismatchGridXYLowerLeft,
    MismatchGridXYSpacing,
    MissingRequiredKeyError,
    read_asc_header,
    read_esri_ascii,
)

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_hugo_read_file_name():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR, "hugo_site.asc"))

    assert isinstance(grid, RasterModelGrid)

    assert field.size == 55 * 76
    assert field.shape == (55 * 76,)
    assert (grid.x_of_node.min(), grid.y_of_node.min()) == (0.0, 0.0)


def test_hugo_read_file_like():
    with open(os.path.join(_TEST_DATA_DIR, "hugo_site.asc")) as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert isinstance(grid, RasterModelGrid)

    assert field.size == 55 * 76
    assert field.shape == (55 * 76,)


def test_hugo_reshape():
    with open(os.path.join(_TEST_DATA_DIR, "hugo_site.asc")) as asc_file:
        (grid, field) = read_esri_ascii(asc_file, reshape=True)

    assert isinstance(grid, RasterModelGrid)

    assert field.shape == (55, 76)


def test_4x3_read_file_name():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR, "4_x_3.asc"))

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert (grid.x_of_node.min(), grid.y_of_node.min()) == (1.0, 2.0)
    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )


def test_4x3_read_file_like():
    with open(os.path.join(_TEST_DATA_DIR, "4_x_3.asc")) as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert isinstance(grid, RasterModelGrid)

    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )


def test_4x3_shape_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4.
5. 6. 7. 8.
9. 10. 11. 12.
        """
    )
    (grid, field) = read_esri_ascii(asc_file)
    assert field.size == 12

    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    (grid, field) = read_esri_ascii(asc_file)
    assert field.size == 12


def test_4x3_size_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10.
        """
    )
    with pytest.raises(DataSizeError):
        read_esri_ascii(asc_file)


def test_4x3_size_mismatch_with_halo():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10.
        """
    )
    with pytest.raises(DataSizeError):
        read_esri_ascii(asc_file, halo=1)


def test_grid_data_size_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((10, 10), xy_spacing=10.0)
    with pytest.raises(MismatchGridDataSizeError):
        read_esri_ascii(asc_file, grid=rmg)


def test_grid_dx_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((4, 3), xy_spacing=15.0)
    with pytest.raises(MismatchGridXYSpacing):
        read_esri_ascii(asc_file, grid=rmg)


def test_grid_lower_left_mismatch():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """
    )
    rmg = RasterModelGrid((4, 3), xy_spacing=10.0, xy_of_lower_left=(10, 15))
    with pytest.raises(MismatchGridXYLowerLeft):
        read_esri_ascii(asc_file, grid=rmg)


def test_header_missing_required_key():
    asc_file = StringIO(
        """
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(MissingRequiredKeyError):
        read_asc_header(asc_file)


def test_header_unknown_key():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
invalid_key   1
        """
    )
    with pytest.raises(BadHeaderLineError):
        read_asc_header(asc_file)


def test_header_missing_value():
    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize
NODATA_value  -9999
invalid_key   1
        """
    )
    with pytest.raises(BadHeaderLineError):
        read_asc_header(asc_file)


def test_header_bad_values():
    asc_file = StringIO(
        """
nrows         -4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(KeyValueError):
        read_asc_header(asc_file)


def test_header_missing_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    with pytest.raises(MissingRequiredKeyError):
        read_asc_header(asc_file)


def test_header_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    header = read_asc_header(asc_file)
    assert header["xllcenter"] == 1.0
    with pytest.raises(KeyError):
        header["xllcorner"]

    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """
    )
    header = read_asc_header(asc_file)
    assert header["xllcorner"] == 1.0
    with pytest.raises(KeyError):
        header["xllcenter"]


def test_header_missing_optional():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
        """
    )
    header = read_asc_header(asc_file)
    with pytest.raises(KeyError):
        header["nodata_value"]


def test_header_case_insensitive():
    asc_file = StringIO(
        """
nCoLs         3
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """
    )
    header = read_asc_header(asc_file)
    for key in ["ncols", "nrows", "xllcenter", "yllcorner", "cellsize", "nodata_value"]:
        assert key in header


def test_header_wrong_type():
    asc_file = StringIO(
        """
nCoLs         3.5
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """
    )
    with pytest.raises(KeyTypeError):
        read_asc_header(asc_file)


def test_name_keyword():
    (grid, field) = read_esri_ascii(
        os.path.join(_TEST_DATA_DIR, "4_x_3.asc"), name="air__temperature"
    )

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field, np.array([9.0, 10.0, 11.0, 6.0, 7.0, 8.0, 3.0, 4.0, 5.0, 0.0, 1.0, 2.0])
    )
    assert_array_almost_equal(grid.at_node["air__temperature"], field)
    assert grid.at_node["air__temperature"] is field


def test_halo_keyword():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR, "4_x_3.asc"), halo=1)

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field,
        np.array(
            [
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                9.0,
                10.0,
                11.0,
                -9999.0,
                -9999.0,
                6.0,
                7.0,
                8.0,
                -9999.0,
                -9999.0,
                3.0,
                4.0,
                5.0,
                -9999.0,
                -9999.0,
                0.0,
                1.0,
                2.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
            ]
        ),
    )


def test_halo_keyword_no_nodata_value():
    (grid, field) = read_esri_ascii(
        os.path.join(_TEST_DATA_DIR, "4_x_3_no_nodata_value.asc"), halo=1
    )

    assert isinstance(grid, RasterModelGrid)

    assert isinstance(field, np.ndarray)
    assert_array_equal(
        field,
        np.array(
            [
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                9.0,
                10.0,
                11.0,
                -9999.0,
                -9999.0,
                6.0,
                7.0,
                8.0,
                -9999.0,
                -9999.0,
                3.0,
                4.0,
                5.0,
                -9999.0,
                -9999.0,
                0.0,
                1.0,
                2.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
                -9999.0,
            ]
        ),
    )
