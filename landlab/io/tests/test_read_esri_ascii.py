#! /usr/bin/env python
"""
Unit tests for landlab.io.esri_ascii module.
"""
import os
import numpy as np
from numpy.testing import assert_array_equal, assert_array_almost_equal
from nose.tools import assert_equal, assert_true, assert_raises
try:
    from nose.tools import assert_is_instance, assert_list_equal, assert_is
except ImportError:
    from landlab.testing.tools import (assert_is_instance, assert_list_equal,
                                       assert_is)
from six import StringIO

from landlab.io import read_esri_ascii, read_asc_header
from landlab.io import (MissingRequiredKeyError, KeyTypeError, DataSizeError,
                        BadHeaderLineError, KeyValueError, 
                        MismatchGridDataSizeError)
from landlab import RasterModelGrid


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def test_hugo_read_file_name():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR,
                                                 'hugo_site.asc'))

    assert_is_instance(grid, RasterModelGrid)

    assert_equal(field.size, 55 * 76)
    assert_equal(field.shape, (55 * 76, ))


def test_hugo_read_file_like():
    with open(os.path.join(_TEST_DATA_DIR, 'hugo_site.asc')) as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert_is_instance(grid, RasterModelGrid)

    assert_equal(field.size, 55 * 76)
    assert_equal(field.shape, (55 * 76, ))


def test_hugo_reshape():
    with open(os.path.join(_TEST_DATA_DIR, 'hugo_site.asc')) as asc_file:
        (grid, field) = read_esri_ascii(asc_file, reshape=True)

    assert_is_instance(grid, RasterModelGrid)

    assert_true(field.shape, (55, 76))


def test_4x3_read_file_name():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR,
                                                 '4_x_3.asc'))

    assert_is_instance(grid, RasterModelGrid)

    assert_is_instance(field, np.ndarray)
    assert_array_equal(field,
                       np.array([9., 10., 11.,
                                 6.,  7.,  8.,
                                 3.,  4.,  5.,
                                 0.,  1.,  2.]))


def test_4x3_read_file_like():
    with open(os.path.join(_TEST_DATA_DIR, '4_x_3.asc')) as asc_file:
        (grid, field) = read_esri_ascii(asc_file)

    assert_is_instance(grid, RasterModelGrid)

    assert_array_equal(field,
                       np.array([9., 10., 11.,
                                 6.,  7.,  8.,
                                 3.,  4.,  5.,
                                 0.,  1.,  2.]))


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
        """)
    (grid, field) = read_esri_ascii(asc_file)
    assert_equal(field.size, 12)

    asc_file = StringIO(
        """
nrows         4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
1. 2. 3. 4. 5. 6. 7. 8. 9. 10. 11. 12.
        """)
    (grid, field) = read_esri_ascii(asc_file)
    assert_equal(field.size, 12)


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
        """)
    assert_raises(DataSizeError, read_esri_ascii, asc_file)
    
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
        """)
    rmg = RasterModelGrid((10,10),10.)
    assert_raises(MismatchGridDataSizeError, read_esri_ascii, asc_file, 
                  grid=rmg)    


def test_header_missing_required_key():
    asc_file = StringIO(
        """
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """)
    assert_raises(MissingRequiredKeyError, read_asc_header, asc_file)


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
        """)
    assert_raises(BadHeaderLineError, read_asc_header, asc_file)


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
        """)
    assert_raises(BadHeaderLineError, read_asc_header, asc_file)


def test_header_bad_values():
    asc_file = StringIO(
        """
nrows         -4
ncols         3
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """)
    assert_raises(KeyValueError, read_asc_header, asc_file)


def test_header_missing_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """)
    assert_raises(MissingRequiredKeyError, read_asc_header, asc_file)


def test_header_mutex_key():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """)
    header = read_asc_header(asc_file)
    assert_equal(header['xllcenter'], 1.)
    # with assert_raises(KeyError):
    #    header['xllcorner']
    assert_raises(KeyError, lambda k: header[k], 'xllcorner')

    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
        """)
    header = read_asc_header(asc_file)
    assert_equal(header['xllcorner'], 1.)
    assert_raises(KeyError, lambda k: header[k], 'xllcenter')


def test_header_missing_optional():
    asc_file = StringIO(
        """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
        """)
    header = read_asc_header(asc_file)
    assert_raises(KeyError, lambda k: header[k], 'nodata_value')


def test_header_case_insensitive():
    asc_file = StringIO(
        """
nCoLs         3
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """)
    header = read_asc_header(asc_file)
    for key in ['ncols', 'nrows', 'xllcenter', 'yllcorner', 'cellsize',
                'nodata_value']:
        assert_true(key in header)


def test_header_wrong_type():
    asc_file = StringIO(
        """
nCoLs         3.5
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
        """)
    assert_raises(KeyTypeError, read_asc_header, asc_file)


def test_name_keyword():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR,
                                                 '4_x_3.asc'),
                                    name='air__temperature')

    assert_is_instance(grid, RasterModelGrid)

    assert_is_instance(field, np.ndarray)
    assert_array_equal(field,
                       np.array([9., 10., 11.,
                                 6.,  7.,  8.,
                                 3.,  4.,  5.,
                                 0.,  1.,  2.]))
    assert_array_almost_equal(grid.at_node['air__temperature'], field)
    assert_is(grid.at_node['air__temperature'], field)
    
def test_halo_keyword():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR, \
                                                 '4_x_3.asc'), \
                                                 halo=1)
                                    
    assert_is_instance(grid, RasterModelGrid)

    assert_is_instance(field, np.ndarray)
    assert_array_equal(field,
                       np.array([-9999., -9999., -9999., -9999., -9999.,  
                                 -9999.,     9.,    10.,    11., -9999.,
                                 -9999.,     6.,     7.,     8., -9999.,
                                 -9999.,     3.,     4.,     5., -9999.,
                                 -9999.,     0.,     1.,     2., -9999.,
                                 -9999., -9999., -9999., -9999., -9999.]))
                                 
def test_halo_keyword_no_nodata_value():
    (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR, \
                                                '4_x_3_no_nodata_value.asc'), \
                                                halo=1)
                                    
    assert_is_instance(grid, RasterModelGrid)

    assert_is_instance(field, np.ndarray)
    assert_array_equal(field,
                       np.array([-9999., -9999., -9999., -9999., -9999.,  
                                 -9999.,     9.,    10.,    11., -9999.,
                                 -9999.,     6.,     7.,     8., -9999.,
                                 -9999.,     3.,     4.,     5., -9999.,
                                 -9999.,     0.,     1.,     2., -9999.,
                                 -9999., -9999., -9999., -9999., -9999.]))


if __name__ == '__main__':
    unittest.main()
