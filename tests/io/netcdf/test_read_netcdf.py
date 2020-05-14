#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""


import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io import (
    MismatchGridDataSizeError,
    MismatchGridXYLowerLeft,
    MismatchGridXYSpacing,
)
from landlab.io.netcdf import read_netcdf

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


def test_read_netcdf(datadir):
    grid = read_netcdf(datadir / "test-netcdf4.nc")
    assert grid.shape == (4, 3)
    assert grid.dy, grid.dx == (1.0, 1.0)
    assert list(grid.at_node.keys()) == ["surface__elevation"]
    assert_array_equal(
        grid.at_node["surface__elevation"],
        [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0],
    )


def test_read_netcdf_64bit(datadir):
    grid = read_netcdf(datadir / "test-netcdf3-64bit.nc")
    assert grid.shape == (4, 3)
    assert grid.dy, grid.dx == (1.0, 1.0)

    grid = RasterModelGrid((6, 5), xy_of_lower_left=(-1.0, -1.0))
    grid = read_netcdf(datadir / "test-netcdf4.nc", grid=grid, halo=1, nodata_value=-1)
    assert_array_equal(
        grid.at_node["surface__elevation"].reshape(grid.shape),
        [
            [-1.0, -1.0, -1.0, -1.0, -1.0],
            [-1.0, 0.0, 1.0, 2.0, -1.0],
            [-1.0, 3.0, 4.0, 5.0, -1.0],
            [-1.0, 6.0, 7.0, 8.0, -1.0],
            [-1.0, 9.0, 10.0, 11.0, -1.0],
            [-1.0, -1.0, -1.0, -1.0, -1.0],
        ],
    )


def test_read_netcdf4_bad_field_name(datadir):
    with pytest.raises(ValueError):
        read_netcdf(datadir / "test-netcdf4.nc", name="not_surface__elevation")


def test_read_netcdf3_64bit(datadir):
    """Test read_netcdf for with 64-bit netcdf3 format."""
    grid = read_netcdf(datadir / "test-netcdf3-64bit.nc")
    assert grid.shape == (4, 3)


def test_read_netcdf4(datadir):
    """Test read_netcdf with netcdf4 format."""
    grid = read_netcdf(datadir / "test-netcdf4.nc")
    assert grid.shape == (4, 3)

    grid = read_netcdf(datadir / "test-netcdf4.nc")
    assert grid.shape == (4, 3)


def test_bad_data_size(datadir):
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((10, 10))
    with pytest.raises(MismatchGridDataSizeError):
        read_netcdf(datadir / "test-netcdf4.nc", grid=grid)


def test_bad_dx(datadir):
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((4, 3), xy_spacing=10)
    with pytest.raises(MismatchGridXYSpacing):
        read_netcdf(datadir / "test-netcdf4.nc", grid=grid)


def test_bad_llc(datadir):
    """Test read_netcdf with netcdf4 format."""
    grid = RasterModelGrid((4, 3), xy_of_lower_left=(-1, -2))
    with pytest.raises(MismatchGridXYLowerLeft):
        read_netcdf(datadir / "test-netcdf4.nc", grid=grid)
