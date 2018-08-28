import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import WITH_NETCDF4, NotRasterGridError, write_raster_netcdf


def test_append_with_time(tmpdir):
    field = RasterModelGrid(4, 3)
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", field, append=False, format="NETCDF4", with_time=True
        )
        field.at_node["topographic__elevation"] *= 2
        write_raster_netcdf("test.nc", field, append=True, format="NETCDF4")

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:],
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]],
                ],
            )
            assert root.variables[name][:].dtype == "int64"
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 2

        root.close()


def test_without_time(tmpdir):
    field = RasterModelGrid(4, 3)
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", field, append=False, format="NETCDF4", with_time=False
        )

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:], [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]]
            )
            assert root.variables[name][:].dtype == "int64"
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        root.close()


def test_with_time(tmpdir):
    field = RasterModelGrid(4, 3)
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", field, append=False, format="NETCDF4", with_time=True
        )

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:], [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]]
            )
            assert root.variables[name][:].dtype == "int64"
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        root.close()
