#! /usr/bin/env python
"""
Unit tests for landlab.io.netcdf module.
"""

import os
import numpy as np
from StringIO import StringIO
from unittest import skipIf
from nose.tools import assert_equal

from landlab.io.netcdf import read_netcdf, WITH_NETCDF4


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def test_read_netcdf3_64bit():
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf3-64bit.nc'))
    assert_equal(grid.shape, (4, 3))


@skipIf(not WITH_NETCDF4, 'netCDF4 package not installed')
def test_read_netcdf4():
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf4.nc'))
    assert_equal(grid.shape, (4, 3))
