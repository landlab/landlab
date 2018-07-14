#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""

import os

import pytest

from landlab.io.netcdf import read_netcdf, WITH_NETCDF4


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), 'data')


def test_read_netcdf3_64bit():
    """Test read_netcdf for with 64-bit netcdf3 format."""
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf3-64bit.nc'))
    assert grid.shape == (4, 3)


@pytest.mark.skipif(not WITH_NETCDF4, reason='netCDF4 package not installed')
def test_read_netcdf4():
    """Test read_netcdf with netcdf4 format."""
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, 'test-netcdf4.nc'))
    assert grid.shape == (4, 3)
