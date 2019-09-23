#! /usr/bin/env python
"""Unit tests for landlab.io.netcdf module."""
import os

import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import WITH_NETCDF4, NotRasterGridError, write_netcdf
from landlab.io.netcdf.read import _get_raster_spacing

try:
    import netCDF4 as nc
except ImportError:
    pass


_TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")


def test_netcdf_write_int64_field_netcdf4(tmpdir):
    """Test write_netcdf with a grid that has an int64 field."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF4")

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(root.variables[name][:].flatten(), field.at_node[name])
            assert root.variables[name][:].dtype == "int64"

        root.close()


def test_netcdf_write_uint8_field_netcdf4(tmpdir):
    """Test write_netcdf with a grid that has an uint8 field."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12, dtype=np.uint8))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF4")

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(root.variables[name][:].flatten(), field.at_node[name])
            assert root.variables[name][:].dtype == "uint8"

        root.close()


def test_netcdf_write_as_netcdf3_64bit(tmpdir):
    """Test write_netcdf with output format 64-bit netcdf3."""
    from scipy.io import netcdf

    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF3_64BIT")

        f = netcdf.netcdf_file("test.nc", "r")

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in f.variables
            assert_array_equal(f.variables[name][:].flatten(), field.at_node[name])

        f.close()


def test_netcdf_write_as_netcdf3_classic(tmpdir):
    """Test write_netcdf with output format classic netcdf3."""
    from scipy.io import netcdf

    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF3_CLASSIC")

        f = netcdf.netcdf_file("test.nc", "r")

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in f.variables
            assert_array_equal(f.variables[name][:].flatten(), field.at_node[name])

        f.close()


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write(tmpdir):
    """Test generic write_netcdf."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF4")
        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        assert set(root.dimensions) == set(["ni", "nj", "nt"])
        assert len(root.dimensions["ni"]) == 3
        assert len(root.dimensions["nj"]) == 4
        assert len(root.dimensions["nt"]) == 1
        assert root.dimensions["nt"].isunlimited()

        assert set(root.variables) == set(["x", "y", "topographic__elevation"])

        assert_array_equal(
            root.variables["x"][:].flatten(),
            np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0, 0.0, 1.0, 2.0]),
        )
        assert_array_equal(
            root.variables["y"][:].flatten(),
            np.array([0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0]),
        )
        assert_array_equal(
            root.variables["topographic__elevation"][:].flatten(),
            field.at_node["topographic__elevation"],
        )

        root.close()


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write_as_netcdf4_classic(tmpdir):
    """Test write_netcdf to netcdf4 classic format."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, format="NETCDF4_CLASSIC")
        root = nc.Dataset("test.nc", "r", format="NETCDF4_CLASSIC")

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in root.variables
            assert_array_equal(root.variables[name][:].flatten(), field.at_node[name])

        root.close()


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write_names_keyword_as_list(tmpdir):
    """Test write_netcdf using a list for the *names* keyword."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf(
            "test.nc", field, names=["topographic__elevation"], format="NETCDF4"
        )
        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        assert "topographic__elevation" in root.variables
        assert "uplift_rate" not in root.variables
        assert_array_equal(
            root.variables["topographic__elevation"][:].flatten(),
            field.at_node["topographic__elevation"],
        )

        root.close()


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write_names_keyword_as_str(tmpdir):
    """Test write_netcdf using a ``str`` for the *names* keyword."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, names="uplift_rate", format="NETCDF4")
        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        assert "topographic__elevation" not in root.variables
        assert "uplift_rate" in root.variables
        assert_array_equal(
            root.variables["uplift_rate"][:].flatten(), field.at_node["uplift_rate"]
        )

        root.close()


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write_names_keyword_as_none(tmpdir):
    """Test write_netcdf using ``None`` for the *names* keyword."""
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", np.arange(12.0))

    with tmpdir.as_cwd():
        write_netcdf("test.nc", field, names=None, format="NETCDF4")
        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in root.variables
            assert_array_equal(root.variables[name][:].flatten(), field.at_node[name])

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


@pytest.mark.skipif(not WITH_NETCDF4, reason="netCDF4 package not installed")
def test_netcdf_write_at_cells(tmpdir):
    """Test write_netcdf using with cell fields"""
    field = RasterModelGrid((4, 3))
    field.add_field("cell", "topographic__elevation", np.arange(field.number_of_cells))
    field.add_field("cell", "uplift_rate", np.arange(field.number_of_cells))

    with tmpdir.as_cwd():
        write_netcdf("test-cells.nc", field, format="NETCDF4")
        root = nc.Dataset("test-cells.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation", "uplift_rate"]:
            assert name in root.variables
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
        root.close()


def test_write_llc():
    pass
