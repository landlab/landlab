#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""

import os

import pytest

from landlab import RasterModelGrid
from landlab.io import (
    MismatchGridDataSizeError,
    MismatchGridXYLowerLeft,
    MismatchGridXYSpacing,
)
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


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_read_netcdf4_bad_field_name():
    with pytest.raises(ValueError):
        read_netcdf(
            os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"),
            name="not_surface__elevation",
        )


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


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_bad_data_size():
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((10, 10))
    with pytest.raises(MismatchGridDataSizeError):
        read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"), grid=grid)


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_bad_dx():
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((4, 3), xy_spacing=10)
    with pytest.raises(MismatchGridXYSpacing):
        read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"), grid=grid)


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_bad_llc():
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((4, 3), xy_of_lower_left=(-1, -2))
    with pytest.raises(MismatchGridXYLowerLeft):
        read_netcdf(os.path.join(_TEST_DATA_DIR, "test-netcdf4.nc"), grid=grid)
