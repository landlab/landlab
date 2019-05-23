import numpy as np
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import write_raster_netcdf

try:
    import netCDF4 as nc
except ImportError:
    pass


def test_append_with_time(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, append=False, format="NETCDF4", time=0)
        field.at_node["topographic__elevation"] *= 2
        write_raster_netcdf("test.nc", field, append=True, format="NETCDF4", time=1.0)

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

        assert "t" in root.variables
        assert_array_equal(root.variables["t"][:], [0.0, 1.0])

        root.close()


def test_without_time(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, append=False, format="NETCDF4")

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:], [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]]
            )
            assert root.variables[name][:].dtype == "int64"
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        assert "t" not in root.variables

        root.close()


def test_with_time(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.ones(12, dtype=np.int64))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, append=False, format="NETCDF4", time=0.0)

        root = nc.Dataset("test.nc", "r", format="NETCDF4")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:], [[[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]]]
            )
            assert root.variables[name][:].dtype == "int64"
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        assert "t" in root.variables
        assert_array_equal(root.variables["t"][:], [0.0])

        root.close()


def test_with_time_netcdf3(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", 2.0 * np.arange(12.0))
    field.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, format="NETCDF3_64BIT", time=10.0)

        root = nc.Dataset("test.nc", "r", format="NETCDF3_64BIT")

        for name in ["uplift_rate", "topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:],
                [
                    [
                        [0.0, 2.0, 4.0],
                        [6.0, 8.0, 10.0],
                        [12.0, 14.0, 16.0],
                        [18.0, 20.0, 22.0],
                    ]
                ],
            )
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        assert "t" in root.variables
        assert_array_equal(root.variables["t"][:], [10.0])

        root.close()


def test_append_with_time_netcdf3(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.ones(12), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", field, append=False, format="NETCDF3_64BIT", time=0
        )
        field.at_node["topographic__elevation"] *= 2
        write_raster_netcdf(
            "test.nc", field, append=True, format="NETCDF3_64BIT", time=1.0
        )

        root = nc.Dataset("test.nc", "r", format="NETCDF3_64BIT")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:],
                [
                    [[1, 1, 1], [1, 1, 1], [1, 1, 1], [1, 1, 1]],
                    [[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]],
                ],
            )
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 2

        assert "t" in root.variables
        assert_array_equal(root.variables["t"][:], [0.0, 1.0])

        root.close()


def test_append_without_time_netcdf3(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("topographic__elevation", np.ones(12), at="node")

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, append=False, format="NETCDF3_64BIT")
        field.at_node["topographic__elevation"] *= 2
        write_raster_netcdf("test.nc", field, append=True, format="NETCDF3_64BIT")

        root = nc.Dataset("test.nc", "r", format="NETCDF3_64BIT")

        for name in ["topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:], [[[2, 2, 2], [2, 2, 2], [2, 2, 2], [2, 2, 2]]]
            )
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        assert "t" not in root.variables

        root.close()


def test_without_time_netcdf3(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", 2.0 * np.arange(12.0))
    field.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_raster_netcdf("test.nc", field, format="NETCDF3_64BIT")

        root = nc.Dataset("test.nc", "r", format="NETCDF3_64BIT")

        for name in ["uplift_rate", "topographic__elevation"]:
            assert name in root.variables
            assert_array_equal(
                root.variables[name][:],
                [
                    [
                        [0.0, 2.0, 4.0],
                        [6.0, 8.0, 10.0],
                        [12.0, 14.0, 16.0],
                        [18.0, 20.0, 22.0],
                    ]
                ],
            )
            assert "nt" in root.dimensions
            assert len(root.dimensions["nt"]) == 1

        assert "t" not in root.variables

        root.close()


def test_names_keyword(tmpdir):
    field = RasterModelGrid((4, 3))
    field.add_field("node", "topographic__elevation", np.arange(12.0))
    field.add_field("node", "uplift_rate", 2.0 * np.arange(12.0))

    with tmpdir.as_cwd():
        write_raster_netcdf(
            "test.nc", field, format="NETCDF3_64BIT", names="uplift_rate"
        )

        root = nc.Dataset("test.nc", "r", format="NETCDF3_64BIT")

        assert "topographic__elevation" not in root.variables
        assert "uplift_rate" in root.variables

        assert_array_equal(
            root.variables["uplift_rate"][:],
            [
                [
                    [0.0, 2.0, 4.0],
                    [6.0, 8.0, 10.0],
                    [12.0, 14.0, 16.0],
                    [18.0, 20.0, 22.0],
                ]
            ],
        )
        assert "nt" in root.dimensions
        assert len(root.dimensions["nt"]) == 1

        root.close()
