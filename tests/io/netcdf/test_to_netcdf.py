#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""
import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal

from landlab import HexModelGrid, RasterModelGrid
from landlab.io.netcdf import to_netcdf


def pytest_generate_tests(metafunc):
    if "format" in metafunc.fixturenames:
        metafunc.parametrize("format", ("NETCDF4", "NETCDF4_CLASSIC", "NETCDF3_CLASSIC", "NETCDF3_64BIT"))
    elif "at" in metafunc.fixturenames:
        metafunc.parametrize("at", ("node", "link", "patch", "corner", "face", "cell"))


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
        if format != "NETCDF4":
            with pytest.raises(RuntimeError):
                to_netcdf(grid, "test.nc", format=format)
        else:
            to_netcdf(grid, "test.nc", format=format)

            assert_array_equal(
                xr.open_dataset("test.nc")["at_node:topographic__elevation"],
                grid.at_node["topographic__elevation"],
            )


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

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format="NETCDF4", at=at)

        actual = xr.open_dataset("test.nc")
        actual_fields = set(
            [name for name in actual.variables if name.startswith("at_")]
        )
        nc_name = "at_{0}:{1}".format(at, name)

        assert actual_fields == set([nc_name])
        assert_array_equal(actual[nc_name], getattr(grid, "at_" + at)[name])


@pytest.mark.parametrize(
    "names,expected", (
        ("elevs", ("at_link:elevs", "at_node:elevs")),
        ("temp", ("at_link:temp",)),
        (("elevs", "temp"), ("at_link:elevs", "at_node:elevs", "at_link:temp")),
        (None, ("at_link:elevs", "at_node:elevs", "at_link:temp")),
    ),
)
def test_names_keyword(tmpdir, names, expected):
    grid = RasterModelGrid((4, 3))

    grid.add_field("elevs", grid.ones(at="node"), at="node")
    grid.add_field("elevs", grid.zeros(at="link"), at="link")
    grid.add_field("temp", grid.ones(at="link"), at="link")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format="NETCDF4", names=names)

        actual = xr.open_dataset("test.nc")
        actual_fields = set(
            [name for name in actual.variables if name.startswith("at_")]
        )
        assert actual_fields == set(expected)


def test_raster_model_grid(tmpdir, format):
    grid = RasterModelGrid((4, 3), xy_spacing=(2, 5), xy_of_lower_left=(-2.0, 10.0))
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = xr.open_dataset("test.nc")
        assert_array_equal(actual["xy_spacing"], [2.0, 5.0])
        assert_array_equal(actual["xy_of_lower_left"], [-2.0, 10.0])
        assert actual.attrs["grid_type"] == "RasterModelGrid"


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
        actual = xr.open_dataset("test.nc")
        assert actual["spacing"] == 2.0
        assert_array_equal(actual["xy_of_lower_left"], [-3.0, -5.0])
        assert actual.attrs == {
            "grid_type": "HexModelGrid",
            "orientation": orientation,
            "node_layout": node_layout,
        }


def test_layers(tmpdir):
    grid = RasterModelGrid((3, 4))
    grid.event_layers.add(10.0, water_depth=[1.0, 2.0])
    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", with_layers=True)
        actual = xr.open_dataset("test.nc")
        actual_fields = set(
            [name for name in actual.variables if name.startswith("at_")]
        )
        assert actual_fields == set(["at_layer:water_depth"])
