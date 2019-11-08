#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""
import netCDF4 as nc
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import NotRasterGridError, write_netcdf
from landlab.io.netcdf.read import _get_raster_spacing


def test_netcdf_write_int64_field(tmpdir, format):
    """Test write_netcdf with a grid that has an int64 field."""
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.arange(12, dtype=np.int64), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format=format)

        root = nc.Dataset("test.nc", "r", format=format)

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(root.variables[name][:].flatten(), field.at_node[name])
            if format == "NETCDF4":
                assert root.variables[name][:].dtype == "int64"
            else:
                assert root.variables[name][:].dtype == "int32"

        root.close()


def test_netcdf_write_uint8_field(tmpdir, format):
    """Test write_netcdf with a grid that has an uint8 field."""
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.arange(12, dtype=np.uint8), at="node")

    with tmpdir.as_cwd():
        if format != "NETCDF4":
            with pytest.raises(RuntimeError):
                write_netcdf("test.nc", field, format=format)
        else:
            write_netcdf("test.nc", field, format=format)
            root = nc.Dataset("test.nc", "r", format=format)

            for name in ["topographic__elevation"]:
                assert name in root.variables
                assert_array_equal(
                    root.variables[name][:],
                    field.at_node[name].reshape((1, 4, 3)),
                )
                assert root.variables[name][:].dtype == "uint8"

            root.close()


def test_netcdf_write(tmpdir, format):
    """Test generic write_netcdf."""
    field = RasterModelGrid((4, 3), xy_spacing=(1.0, 2.0), xy_of_lower_left=(-1.0, 2.0))
    field.add_field("topographic__elevation", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format=format)
        root = nc.Dataset("test.nc", "r", format=format)

        assert set(root.dimensions) == set(["ni", "nj", "nt"])
        assert len(root.dimensions["ni"]) == 3
        assert len(root.dimensions["nj"]) == 4
        assert len(root.dimensions["nt"]) == 1
        assert root.dimensions["nt"].isunlimited()

        assert set(root.variables) == set(["x", "y", "topographic__elevation"])

        assert_array_equal(root.variables["x"], field.x_of_node.reshape((4, 3)))
        assert_array_equal(root.variables["y"], field.y_of_node.reshape((4, 3)))
        assert_array_equal(
            root.variables["topographic__elevation"],
            field.at_node["topographic__elevation"].reshape((1, 4, 3)),
        )

        root.close()


@pytest.mark.parametrize(
    "names", ("uplift_rate", ["uplift_rate"], ("uplift_rate",), {"uplift_rate"})
)
def test_netcdf_write_some_names(tmpdir, format, names):
    """Test write_netcdf using a ``str`` for the *names* keyword."""
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.arange(12.0), at="node")
    field.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, names=names, format=format)
        root = nc.Dataset("test.nc", "r", format=format)

        assert "topographic__elevation" not in root.variables
        assert "uplift_rate" in root.variables
        assert_array_equal(
            root.variables["uplift_rate"],
            field.at_node["uplift_rate"].reshape((1, 4, 3)),
        )

        root.close()


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
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.arange(12.0), at="node")
    field.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, names=names, format=format)
        root = nc.Dataset("test.nc", "r", format=format)

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name],
                field.at_node[name].reshape((1, 4, 3)),
            )
        root.close()


@pytest.mark.parametrize("names", ([], (), {}))
def test_netcdf_write_no_names(tmpdir, format, names):
    """Test write_netcdf using ``None`` for the *names* keyword."""
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.arange(12.0), at="node")
    field.add_field("uplift_rate", np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, names=names, format=format)
        root = nc.Dataset("test.nc", "r", format=format)

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name not in root.variables
        root.close()


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
    field = RasterModelGrid((4, 3), xy_spacing=(2.0, 3.0), xy_of_lower_left=(3.0, 0.5))
    field.add_field("topographic__elevation", np.arange(field.number_of_cells), at="cell")
    field.add_field("uplift_rate", np.arange(field.number_of_cells), at="cell")

    with tmpdir.as_cwd():
        write_netcdf("test-cells.nc", field, format=format)
        root = nc.Dataset("test-cells.nc", "r", format=format)

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name],
                field.at_cell[name].reshape((1, 2, 1)),
            )
            assert_array_equal(root.variables[name][:].flatten(), field.at_cell[name])

        assert set(root.dimensions) == set(["nv", "ni", "nj", "nt"])
        assert len(root.dimensions["nv"]) == 4
        assert len(root.dimensions["ni"]) == 1
        assert len(root.dimensions["nj"]) == 2
        assert len(root.dimensions["nt"]) == 1
        assert root.dimensions["nt"].isunlimited()

        assert set(root.variables) == set(
            ["x_bnds", "y_bnds", "topographic__elevation", "uplift_rate"]
        )

        assert_array_equal(
            root.variables["x_bnds"],
            field.x_of_corner[field.corners_at_cell].reshape((2, 1, 4)),
        )
        assert_array_equal(
            root.variables["y_bnds"],
            field.y_of_corner[field.corners_at_cell].reshape((2, 1, 4)),
        )

        root.close()
