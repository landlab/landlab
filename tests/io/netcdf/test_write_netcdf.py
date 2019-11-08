#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""
import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import NotRasterGridError, write_netcdf
from landlab.io.netcdf.read import _get_raster_spacing


def test_netcdf_write_int64_field(tmpdir, format):
    """Test write_netcdf with a grid that has an int64 field."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12, dtype=np.int64), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", grid, format=format)

        with xr.open_dataset("test.nc") as actual:
            assert_array_equal(
                actual["topographic__elevation"],
                grid.at_node["topographic__elevation"].reshape((1, 4, 3)),
            )
            if format == "NETCDF4":
                assert actual["topographic__elevation"].dtype == "int64"
            else:
                assert actual["topographic__elevation"].dtype == "int32"


def test_netcdf_write_uint8_field(tmpdir, format):
    """Test write_netcdf with a grid that has a uint8 field."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12, dtype=np.uint8), at="node")

    with tmpdir.as_cwd():
        if format != "NETCDF4":
            with pytest.raises(RuntimeError):
                write_netcdf("test.nc", grid, format=format)
        else:
            write_netcdf("test.nc", grid, format=format)
            with xr.open_dataset("test.nc") as actual:
                assert_array_equal(
                    actual["topographic__elevation"],
                    grid.at_node["topographic__elevation"].reshape((1, 4, 3)),
                )
                assert actual["topographic__elevation"].dtype == "uint8"


def test_netcdf_write(tmpdir, format):
    """Test generic write_netcdf."""
    grid = RasterModelGrid((4, 3), xy_spacing=(1.0, 2.0), xy_of_lower_left=(-1.0, 2.0))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", grid, format=format)
        with xr.open_dataset("test.nc") as actual:
            assert set(actual.dims) == set(["ni", "nj", "nt"])
            assert actual.dims["ni"] == 3
            assert actual.dims["nj"] == 4
            assert actual.dims["nt"] == 1

            assert set(actual.variables) == set(["x", "y", "topographic__elevation"])

            assert_array_equal(actual["x"], grid.x_of_node.reshape((4, 3)))
            assert_array_equal(actual["y"], grid.y_of_node.reshape((4, 3)))
            assert_array_equal(
                actual["topographic__elevation"],
                grid.at_node["topographic__elevation"].reshape((1, 4, 3)),
            )


@pytest.mark.parametrize(
    "names", ("uplift_rate", ["uplift_rate"], ("uplift_rate",), {"uplift_rate"})
)
def test_netcdf_write_some_names(tmpdir, format, names):
    """Test write_netcdf using a ``str`` for the *names* keyword."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")
    grid.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", grid, names=names, format=format)
        with xr.open_dataset("test.nc") as actual:
            assert "topographic__elevation" not in actual
            assert_array_equal(
                actual["uplift_rate"], grid.at_node["uplift_rate"].reshape((1, 4, 3))
            )


@pytest.mark.parametrize(
    "names", (
        None,
        ("topographic__elevation", "uplift_rate"),
        ["topographic__elevation", "uplift_rate"],
        {"topographic__elevation", "uplift_rate"},
    ),
)
def test_netcdf_write_all_names(tmpdir, format, names):
    """Test write_netcdf using ``None`` for the *names* keyword."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")
    grid.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", grid, names=names, format=format)
        with xr.open_dataset("test.nc") as actual:
            for name in ["topographic__elevation", "uplift_rate"]:
                assert name in actual
                assert_array_equal(actual[name], grid.at_node[name].reshape((1, 4, 3)))


@pytest.mark.parametrize("names", ([], (), {}))
def test_netcdf_write_no_names(tmpdir, format, names):
    """Test write_netcdf using ``None`` for the *names* keyword."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")
    grid.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", grid, names=names, format=format)
        with xr.open_dataset("test.nc") as actual:
            for name in ["topographic__elevation", "uplift_rate"]:
                assert name not in actual


def test_2d_unit_spacing():
    """Test write_netcdf with a 2D grid with unit spacing."""
    (x, y) = np.meshgrid(np.arange(5.0), np.arange(4.0))

    spacing = _get_raster_spacing((y, x))
    assert spacing == 1.0


def test_2d_non_unit_spacing():
    """Test _get_raster_spacing with a 2D grid with non-unit spacing."""
    (x, y) = np.meshgrid(np.arange(5.0) * 2, np.arange(4.0) * 2)

    spacing = _get_raster_spacing((y, x))
    assert spacing == 2.0


def test_2d_uneven_spacing_axis_0():
    """Test _get_raster_spacing with a 2D grid with uneven spacing in y."""
    (x, y) = np.meshgrid(np.logspace(0.0, 2.0, num=5), np.arange(4.0))

    with pytest.raises(NotRasterGridError):
        _get_raster_spacing((y, x))


def test_2d_uneven_spacing_axis_1():
    """Test _get_raster_spacing with a 2D grid with uneven spacing in x."""
    (x, y) = np.meshgrid(np.arange(4.0), np.logspace(0.0, 2.0, num=5))

    with pytest.raises(NotRasterGridError):
        _get_raster_spacing((y, x))


def test_2d_switched_coords():
    """Test _get_raster_spacing with a 2D grid when the spacing is switched."""
    (x, y) = np.meshgrid(np.arange(5.0), np.arange(4.0))

    spacing = _get_raster_spacing((x, y))
    assert spacing == 0.0


def test_1d_unit_spacing():
    """Test _get_raster_spacing with a 1D grid with unit spacing."""
    spacing = _get_raster_spacing((np.arange(5.0),))
    assert spacing == 1.0


def test_1d_non_unit_spacing():
    """Test _get_raster_spacing with a 1D grid with non-unit spacing."""
    spacing = _get_raster_spacing((np.arange(5.0) * 2,))
    assert spacing == 2.0


def test_1d_uneven_spacing():
    """Test _get_raster_spacing with a 1D grid with uneven spacing in y."""
    with pytest.raises(NotRasterGridError):
        _get_raster_spacing((np.logspace(0.0, 2.0, num=5),))


def test_netcdf_write_at_cells(tmpdir, format):
    """Test write_netcdf using with cell fields"""
    grid = RasterModelGrid((4, 3), xy_spacing=(2.0, 3.0), xy_of_lower_left=(3.0, 0.5))
    grid.add_field("topographic__elevation", np.arange(grid.number_of_cells), at="cell")
    grid.add_field("uplift_rate", np.arange(grid.number_of_cells), at="cell")

    with tmpdir.as_cwd():
        write_netcdf("test-cells.nc", grid, format=format)
        with xr.open_dataset("test-cells.nc") as actual:
            for name in ["topographic__elevation", "uplift_rate"]:
                assert_array_equal(actual[name], grid.at_cell[name].reshape((1, 2, 1)))

            assert set(actual.dims) == set(["nv", "ni", "nj", "nt"])
            assert actual.dims["nv"] == 4
            assert actual.dims["ni"] == 1
            assert actual.dims["nj"] == 2
            assert actual.dims["nt"] == 1

            assert set(actual.variables) == {
                "x_bnds", "y_bnds", "topographic__elevation", "uplift_rate"
            }

            assert_array_equal(
                actual["x_bnds"],
                grid.x_of_corner[grid.corners_at_cell].reshape((2, 1, 4)),
            )
            assert_array_equal(
                actual["y_bnds"],
                grid.y_of_corner[grid.corners_at_cell].reshape((2, 1, 4)),
            )
