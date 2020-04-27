import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.netcdf import from_netcdf, to_netcdf


@pytest.mark.parametrize("include", ((), [], set(), None))
def test_include_keyword_is_empty(tmpdir, format, include):
    grid = RasterModelGrid((4, 3), xy_spacing=(2, 5), xy_of_lower_left=(-2.0, 10.0))
    grid.add_ones("elev", at="node")
    grid.add_zeros("elev", at="link")
    grid.add_empty("temp", at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = from_netcdf("test.nc", include=include)
        assert len(actual.at_node) == 0
        assert len(actual.at_link) == 0


@pytest.mark.parametrize("include", ("*", ("*",), ("at_node:*", "at_link:*")))
@pytest.mark.parametrize("exclude", (None, ()))
def test_include_everything(tmpdir, format, include, exclude):
    grid = RasterModelGrid((4, 3), xy_spacing=(2, 5), xy_of_lower_left=(-2.0, 10.0))
    grid.add_ones("elev", at="node")
    grid.add_zeros("elev", at="link")
    grid.add_empty("temp", at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = from_netcdf("test.nc", include=include)
        assert set(actual.at_node) == set(["elev", "temp"])
        assert set(actual.at_link) == set(["elev"])


@pytest.mark.parametrize(
    "include,exclude", [(("*", "*")), ((None, None)), (([], None))]
)
def test_exclude_everything(tmpdir, format, include, exclude):
    grid = RasterModelGrid((4, 3), xy_spacing=(2, 5), xy_of_lower_left=(-2.0, 10.0))
    grid.add_ones("elev", at="node")
    grid.add_zeros("elev", at="link")
    grid.add_empty("temp", at="node")

    with tmpdir.as_cwd():
        to_netcdf(grid, "test.nc", format=format)
        actual = from_netcdf("test.nc", include=include, exclude=exclude)
        assert len(actual.at_node) == 0
        assert len(actual.at_link) == 0


@pytest.mark.parametrize(
    "grid_type", ["HexModelGrid", "RadialModelGrid", "RasterModelGrid"]
)
def test_from_grid(datadir, grid_type):
    grid = from_netcdf(datadir / "test-{0}.nc".format(grid_type))
    assert grid.__class__.__name__ == grid_type
    assert_array_equal(grid.at_node["elev"], 1.0)
    assert_array_equal(grid.at_node["temp"], 1.0)
    assert_array_equal(grid.at_link["elev"], 0.0)
