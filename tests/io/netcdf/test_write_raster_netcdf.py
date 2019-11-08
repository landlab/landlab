import numpy as np
import xarray as xr
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import write_raster_netcdf


def test_append_with_time(tmpdir):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.ones(12, dtype=np.int64), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, append=False, format="NETCDF4", time=0)
        grid.at_node["topographic__elevation"] *= 2
        write_raster_netcdf("test.nc", grid, append=True, format="NETCDF4", time=1.0)

        with xr.open_dataset("test.nc") as actual:
            assert_array_equal(
                actual["topographic__elevation"],
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]],
                ],
            )
            assert actual["topographic__elevation"].dtype == "int64"
            assert "nt" in actual.dims
            assert actual.dims["nt"] == 2

            assert_array_equal(actual["t"], [0.0, 1.0])


def test_without_time(tmpdir):
    grid = RasterModelGrid((4, 3))
    grid.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, append=False, format="NETCDF4")

        with xr.open_dataset("test.nc") as actual:
            assert_array_equal(
                actual["topographic__elevation"],
                [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]],
            )
            assert actual["topographic__elevation"].dtype == "int64"
            assert actual.dims["nt"] == 1

        assert "t" not in actual


def test_with_time(tmpdir):
    grid = RasterModelGrid((4, 3))
    grid.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, append=False, format="NETCDF4", time=0.0)

        with xr.open_dataset("test.nc") as actual:
            assert_array_equal(
                actual["topographic__elevation"],
                [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]],
            )
            assert actual["topographic__elevation"].dtype == "int64"
            assert actual.dims["nt"] == 1

            assert_array_equal(actual["t"], [0.0])


def test_with_time_netcdf3(tmpdir):
    grid = RasterModelGrid((4, 3))
    grid.add_field("node", "topographic__elevation", 2.0 * np.arange(12.0))
    grid.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, format="NETCDF3_64BIT", time=10.0)
        with xr.open_dataset("test.nc") as actual:
            for name in ["uplift_rate", "topographic__elevation"]:
                assert_array_equal(
                    actual[name],
                    [
                        [
                            [0.0, 2.0, 4.0],
                            [6.0, 8.0, 10.0],
                            [12.0, 14.0, 16.0],
                            [18.0, 20.0, 22.0],
                        ]
                    ],
                )
                assert actual.dims["nt"] == 1
            assert_array_equal(actual["t"], [10.0])


def test_append_with_time_netcdf3(tmpdir):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.ones(12), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", grid, append=False, format="NETCDF3_64BIT", time=0
        )
        grid.at_node["topographic__elevation"] *= 2
        write_raster_netcdf(
            "test.nc", grid, append=True, format="NETCDF3_64BIT", time=1.0
        )

        with xr.open_dataset("test.nc") as actual:
            assert_array_equal(
                actual["topographic__elevation"],
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]],
                ],
            )
            assert actual.dims["nt"] == 2

        assert_array_equal(actual["t"], [0.0, 1.0])


def test_append_without_time_netcdf3(tmpdir, format):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.ones(12), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, append=False, format=format)
        grid.at_node["topographic__elevation"] *= 2
        write_raster_netcdf("test.nc", grid, append=True, format=format)

        with xr.open_dataset("test.nc") as actual:
            assert len(actual["nt"]) == 2
            assert_array_equal(
                actual["topographic__elevation"],
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]],
                ]
            )
            assert_array_equal(actual["t"], [0.0, 1.0])


def test_without_time_netcdf3(tmpdir, format):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", 2.0 * np.arange(12.0), at="node")
    grid.add_field("uplift_rate", 2.0 * np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", grid, format=format)

        with xr.open_dataset("test.nc") as actual:
            for name in grid.at_node:
                assert_array_equal(
                    actual[name],
                    grid.at_node[name].reshape((1, 4, 3))
                )
            assert "t" not in actual


def test_names_keyword(tmpdir, format):
    grid = RasterModelGrid((4, 3))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")
    grid.add_field("uplift_rate", 2.0 * np.arange(12.0), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", grid, format=format, names="uplift_rate"
        )

        with xr.open_dataset("test.nc") as actual:
            assert "topographic__elevation" not in actual
            assert_array_equal(
                actual["uplift_rate"], grid.at_node["uplift_rate"].reshape((1, 4, 3))
            )
            assert "nt" in actual.dims
            assert len(actual["nt"]) == 1
