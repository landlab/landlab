#! /usr/bin/env python
"""
Unit tests for landlab.io.esri_ascii module.
"""

import unittest
import os
import numpy as np
from StringIO import StringIO

from landlab.io import read_esri_ascii, read_asc_header
from landlab.io import MissingRequiredKeyError, KeyTypeError, DataSizeError
from landlab import RasterModelGrid


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class TestReadEsriAsciiFileHugo(unittest.TestCase):
    def test_read_file_name(self):
        (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR,
                                                     'hugo_site.asc'))

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertEqual(field.size, 55 * 76)
        self.assertEqual(field.shape, (55 * 76, ))

    def test_read_file_like(self):
        with open(os.path.join(_TEST_DATA_DIR, 'hugo_site.asc')) as asc_file:
            (grid, field) = read_esri_ascii(asc_file)

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertEqual(field.size, 55 * 76)
        self.assertEqual(field.shape, (55 * 76, ))

    def test_reshape(self):
        with open(os.path.join(_TEST_DATA_DIR, 'hugo_site.asc')) as asc_file:
            (grid, field) = read_esri_ascii(asc_file, reshape=True)

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertTrue(field.shape, (55, 76))


class TestReadEsriAsciiFile4x3(unittest.TestCase):
    def test_read_file_name(self):
        (grid, field) = read_esri_ascii(os.path.join(_TEST_DATA_DIR,
                                                     '4_x_3.asc'))

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertIsInstance(field, np.ndarray)
        self.assertEqual(field.shape, (4 * 3, ))
        self.assertListEqual(list(field.flat), [0., 1., 2.,
                                                3., 4., 5.,
                                                6., 7., 8.,
                                                9., 10., 11.,
                                               ])

    def test_read_file_like(self):
        with open(os.path.join(_TEST_DATA_DIR, '4_x_3.asc')) as asc_file:
            (grid, field) = read_esri_ascii(asc_file)

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertIsInstance(field, np.ndarray)
        self.assertEqual(field.shape, (4 * 3, ))
        self.assertListEqual(list(field.flat), [0., 1., 2.,
                                                3., 4., 5.,
                                                6., 7., 8.,
                                                9., 10., 11.,
                                               ])

    def test_shape_mismatch(self):
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
        self.assertEqual(field.size, 12)

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
        self.assertEqual(field.size, 12)

    def test_size_mismatch(self):
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
        with self.assertRaises(DataSizeError):
            (grid, field) = read_esri_ascii(asc_file)


class TestReadEsriAsciiFileHeader(unittest.TestCase):
    def test_missing_required_key(self):
        asc_file = StringIO(
            """
nrows         4
xllcorner     1.
yllcorner     2.
cellsize      10.
NODATA_value  -9999
            """)
        with self.assertRaises(MissingRequiredKeyError):
            header = read_asc_header(asc_file)

    def test_missing_mutex_key(self):
        asc_file = StringIO(
            """
ncols         3
nrows         4
yllcorner     2.
cellsize      10.
NODATA_value  -9999
            """)
        with self.assertRaises(MissingRequiredKeyError):
            header = read_asc_header(asc_file)

    def test_mutex_key(self):
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
        self.assertEqual(header['xllcenter'], 1.)
        with self.assertRaises(KeyError):
            header['xllcorner']

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
        self.assertEqual(header['xllcorner'], 1.)
        with self.assertRaises(KeyError):
            header['xllcenter']

    def test_missing_optional(self):
        asc_file = StringIO(
            """
ncols         3
nrows         4
xllcenter     1.
yllcorner     2.
cellsize      10.
            """)
        header = read_asc_header(asc_file)
        with self.assertRaises(KeyError):
            header['nodata_value']

    def test_case_insensitive(self):
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
            self.assertTrue(key in header)

    def test_wrong_type(self):
        asc_file = StringIO(
            """
nCoLs         3.5
nrows         4
Xllcenter     1.
YLLCORNER     2.
CELLSIZE      10.
NODATA_value  -999
            """)
        with self.assertRaises(KeyTypeError):
            header = read_asc_header(asc_file)


if __name__ == '__main__':
    unittest.main()
