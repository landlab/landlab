#! /usr/bin/env python
"""
Unit tests for landlab.io.esri_ascii module.
"""

import unittest
import os
import numpy as np

from landlab.io import read_esri_ascii
from landlab import RasterModelGrid


class TestReadEsriAsciiFileHugo(unittest.TestCase):
    def test_read_file_name(self):
        (grid, field) = read_esri_ascii(os.path.join('data', 'hugo_site.asc'))

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertTrue(field.shape, (55, 76))

    def test_read_file_like(self):
        with open(os.path.join('data', 'hugo_site.asc')) as asc_file:
            (grid, field) = read_esri_ascii(asc_file)

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertTrue(field.shape, (55, 76))


class TestReadEsriAsciiFile4x3(unittest.TestCase):
    def test_read_file_name(self):
        (grid, field) = read_esri_ascii(os.path.join('data', '4_x_3.asc'))

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertIsInstance(field, np.ndarray)
        self.assertEqual(field.shape, (4, 3))
        self.assertListEqual(list(field.flat), [0., 1., 2.,
                                                3., 4., 5.,
                                                6., 7., 8.,
                                                9., 10., 11.,
                                               ])

    def test_read_file_like(self):
        with open(os.path.join('data', '4_x_3.asc')) as asc_file:
            (grid, field) = read_esri_ascii(asc_file)

        self.assertIsInstance(grid, RasterModelGrid)

        self.assertIsInstance(field, np.ndarray)
        self.assertEqual(field.shape, (4, 3))
        self.assertListEqual(list(field.flat), [0., 1., 2.,
                                                3., 4., 5.,
                                                6., 7., 8.,
                                                9., 10., 11.,
                                               ])


if __name__ == '__main__':
    unittest.main()
