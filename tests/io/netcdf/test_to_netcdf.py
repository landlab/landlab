#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""
import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.io.netcdf import from_netcdf, to_netcdf


def test_netcdf_write_int64(tmpdir, format):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12, dtype=np.int64), at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)

        actual = xr.open_dataset("test.nc")
        values = actual["at_node:topographic__elevation"]
        assert_array_equal(values, grid.at_node["topographic__elevation"])
        if format == "NETCDF4":
            assert values.dtype == "int64"
        else:
            assert values.dtype == "int32"


def test_netcdf_write_uint8(tmpdir, format):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12, dtype=np.uint8), at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)

        actual = xr.open_dataset("test.nc")["at_node:topographic__elevation"]
        assert_array_equal(actual, grid.at_node["topographic__elevation"])
        assert actual.dtype == np.uint8 if format == "NETCDF4" else np.int8


@pytest.mark.parametrize("dtype", ("int32", "float32", "float64"))
def test_netcdf_write_dtype(tmpdir, format, dtype):
    """Test write_netcdf with a grid that has an uint8 field."""
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12, dtype=dtype), at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = xr.open_dataset("test.nc")["at_node:topographic__elevation"]

        assert_array_equal(actual, grid.at_node["topographic__elevation"])
        assert actual.dtype == dtype


def test_at_keyword(tmpdir, at):
    grid = RasterModelGrid((4, 3))

    name = "topographic__elevation"
    for src_at in {"node", "link", "patch", "corner", "face", "cell"}:
        grid.add_field(name, grid.ones(at=src_at) * 10.0, at=src_at)

    include = "at_{0}:*".format(at)
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format="NETCDF4", include=include)

        with xr.open_dataset("test.nc") as actual:
            actual_fields = set(
                [name for name in actual.variables if name.startswith("at_")]
            )
            nc_name = "at_{0}:{1}".format(at, name)

            assert actual_fields == set([nc_name])
            assert_array_equal(actual[nc_name], getattr(grid, "at_" + at)[name])


def test_raster_model_grid(tmpdir, format):
    grid = RasterModelGrid((4, 3), xy_spacing=(2, 5), xy_of_lower_left=(-2.0, 10.0))
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = from_netcdf("test.nc")
        assert (actual.dx, actual.dy) == (grid.dx, grid.dy)
        assert actual.xy_of_lower_left == grid.xy_of_lower_left


@pytest.mark.parametrize("orientation", ("horizontal", "vertical"))
@pytest.mark.parametrize("node_layout", ("rect", "hex"))
def test_hex_model_grid(tmpdir, format, orientation, node_layout):
    grid = HexModelGrid(
        shape=(4, 5),
        spacing=2.0,
        xy_of_lower_left=(-3, -5),
        orientation=orientation,
        node_layout=node_layout,
    )
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = from_netcdf("test.nc")

        assert actual.spacing == grid.spacing
        assert actual.xy_of_lower_left == grid.xy_of_lower_left
        assert actual.orientation == grid.orientation
        assert actual.node_layout == grid.node_layout


def test_layers(tmpdir, format):
    grid = RasterModelGrid((3, 4))
    grid.event_layers.add(10.0, water_depth=[1.0, 2.0])
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", include="at_layer*", format=format)
        actual = xr.open_dataset("test.nc")
        actual_fields = set(
            [name for name in actual.variables if name.startswith("at_")]
        )
        assert actual_fields == set(["at_layer:water_depth"])


@pytest.mark.parametrize("mode", ("w", "a"))
def test_with_and_without_time(tmpdir, format, mode):
    grid = RasterModelGrid((3, 4))
    grid.add_full("elevation", 1.0, at="node")
    with tmpdir.as_cwd():
        to_netcdf(grid, "test-without-time.nc", format=format, mode=mode)
        with xr.open_dataset("test-without-time.nc") as actual:
            assert "time" not in actual.dims
            assert "time" not in actual.variables
            assert actual["at_node:elevation"].dims == ("node",)

        to_netcdf(grid, "test-with-time.nc", format=format, time=10.0, mode=mode)
        with xr.open_dataset("test-with-time.nc") as actual:
            assert "time" in actual.dims
            assert "time" in actual.variables
            assert actual["time"] == [10.0]
            assert actual["at_node:elevation"].dims == ("time", "node")


def test_append_with_new_field(tmpdir, format):
    grid = RasterModelGrid((3, 4))
    grid.add_full("elevation", 1.0, at="node")
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        grid.add_full("temperature", 2.0, at="node")
        to_netcdf(grid, "test.nc", format=format, mode="a", time=10.0)
