#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""

import os

import pytest

from landlab.io.netcdf import WITH_NETCDF4, read_netcdf

_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

grid_mapping_keys = [
    "grid_mapping_name",
    "longitude_of_central_meridian",
    "false_easting",
    "false_northing",
    "latitude_of_projection_origin",
    "scale_factor_at_central_meridian",
    "long_name",
    "longitude_of_prime_meridian",
    "semi_major_axis",
    "inverse_flattening",
    "spatial_ref",
    "GeoTransform",
]


def test_read_netcdf3_64bit():
    """Test read_netcdf for with 64-bit netcdf3 format."""
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf3-64bit.nc"))
    assert grid.shape == (4, 3)


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_read_netcdf4():
    """Test read_netcdf with netcdf4 format."""
    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"))
    assert grid.shape == (4, 3)

    grid = read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"))
    assert grid.shape == (4, 3)


def test_read_llc():
    pass
