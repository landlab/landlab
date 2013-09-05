#! /usr/bin/env python
"""
Unit tests for landlab.io.netcdf module.
"""

import unittest
import os
import numpy as np
from StringIO import StringIO

from landlab import RasterModelGrid
from landlab.io.netcdf import write_netcdf, NotRasterGridError
from landlab.io.netcdf.read import _get_raster_spacing


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


class TestWriteNetcdf(unittest.TestCase):
    def test_write(self):
        field = RasterModelGrid(4, 3)
        field.new_field_location('node', 12.)
        values = np.arange(12.)
        field.add_field('node', 'planet_surface__elevation', values)

        write_netcdf('test.nc', field)

        import netCDF4 as nc

        root = nc.Dataset('test.nc', 'r', format='NETCDF4')
        self.assertEqual(set(root.dimensions), set(['ni', 'nj', 'nt']))
        self.assertEqual(len(root.dimensions['ni']), 3)
        self.assertEqual(len(root.dimensions['nj']), 4)
        self.assertTrue(len(root.dimensions['nt']), 1)
        self.assertTrue(root.dimensions['nt'].isunlimited())

        self.assertEqual(set(root.variables),
                         set(['x', 'y', 'planet_surface__elevation']))

        self.assertListEqual(list(root.variables['x'][:].flat),
                             [0., 1., 2.,
                              0., 1., 2.,
                              0., 1., 2.,
                              0., 1., 2., ])
        self.assertListEqual(list(root.variables['y'][:].flat),
                             [0., 0., 0.,
                              1., 1., 1.,
                              2., 2., 2.,
                              3., 3., 3., ])
        self.assertListEqual(
            list(root.variables['planet_surface__elevation'][:].flat),
            range(12))
        root.close()


class TestGetRasterSpacing2D(unittest.TestCase):
    def test_unit_spacing(self):
        (x, y) = np.meshgrid(np.arange(5.), np.arange(4.))

        spacing = _get_raster_spacing((y, x))
        self.assertEqual(spacing, 1.)

    def test_non_unit_spacing(self):
        (x, y) = np.meshgrid(np.arange(5.) * 2, np.arange(4.) * 2)

        spacing = _get_raster_spacing((y, x))
        self.assertEqual(spacing, 2.)

    def test_uneven_spacing_axis_0(self):
        (x, y) = np.meshgrid(np.logspace(0., 2., num=5), np.arange(4.))

        with self.assertRaises(NotRasterGridError):
            spacing = _get_raster_spacing((y, x))

    def test_uneven_spacing_axis_1(self):
        (x, y) = np.meshgrid(np.arange(4.), np.logspace(0., 2., num=5))

        with self.assertRaises(NotRasterGridError):
            spacing = _get_raster_spacing((y, x))

    def test_switched_coords(self):
        (x, y) = np.meshgrid(np.arange(5.), np.arange(4.))

        spacing = _get_raster_spacing((x, y))
        self.assertEqual(spacing, 0.)


class TestGetRasterSpacing1D(unittest.TestCase):
    def test_unit_spacing(self):
        spacing = _get_raster_spacing((np.arange(5.), ))
        self.assertEqual(spacing, 1.)

    def test_non_unit_spacing(self):
        spacing = _get_raster_spacing((np.arange(5.) * 2, ))
        self.assertEqual(spacing, 2.)

    def test_uneven_spacing(self):
        (x, y) = np.meshgrid(np.logspace(0., 2., num=5), np.arange(4.))

        with self.assertRaises(NotRasterGridError):
            spacing = _get_raster_spacing((np.logspace(0., 2., num=5), ))


if __name__ == '__main__':
    unittest.main()
