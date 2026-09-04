import io
import os
from itertools import permutations

import numpy as np
import pytest
import xarray as xr
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.io.sandsuet import _validate_field_location
from landlab.io.sandsuet import dump


@pytest.fixture
def grid_with_fields():
    grid = RasterModelGrid((3, 4), xy_spacing=(10.0, 20.0))
    grid.add_field("topographic__elevation", np.arange(12.0), at="node")
    grid.add_field("surface_water__depth", np.zeros(12), at="node")
    return grid


def test_dump_returns_bytes_when_no_stream(grid_with_fields):
    result = dump(grid_with_fields)
    assert isinstance(result, memoryview)


def test_dump_returns_none_when_stream_given(grid_with_fields):
    stream = io.BytesIO()
    result = dump(grid_with_fields, stream=stream)
    assert result is None


def test_dump_writes_to_stream(grid_with_fields):
    stream = io.BytesIO()
    dump(grid_with_fields, stream=stream)
    assert stream.seek(0, os.SEEK_END) > 0


def test_sandsuet_version_attribute(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    assert ds.attrs["sandsuet_version"] == "1.0.0"


@pytest.mark.parametrize("at", ("node", "corner", "patch", "cell"))
def test_field_values(at):
    grid = RasterModelGrid((4, 5), xy_spacing=(10.0, 20.0))
    actual = grid.add_ones("foo", at=at)
    ds = xr.open_dataset(dump(grid))

    shape = (len(ds["y"].values), len(ds["x"].values))

    assert_array_equal(ds["foo"], actual.reshape(shape))
    assert shape[0] * shape[1] == grid.number_of_elements(at)


def test_coordinates_present(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    assert "x" in ds.coords
    assert "y" in ds.coords


def test_coordinate_values(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    assert_array_equal(ds["x"].values, [0.0, 10.0, 20.0, 30.0])
    assert_array_equal(ds["y"].values, [0.0, 20.0, 40.0])


def test_fields_present_as_data_vars(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    assert set(ds) == {"surface_water__depth", "topographic__elevation"}


def test_field_values_match(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    expected = np.arange(12.0).reshape(grid_with_fields.shape)
    assert_array_equal(ds["topographic__elevation"].values, expected)

    expected = np.full(grid_with_fields.shape, 0.0)
    assert_array_equal(ds["surface_water__depth"].values, expected)


def test_field_dimensions(grid_with_fields):
    ds = xr.open_dataset(dump(grid_with_fields))
    assert ds["topographic__elevation"].dims == ("y", "x")


def test_include_single_field(grid_with_fields):
    ds = xr.open_dataset(
        dump(grid_with_fields, include="at_node:topographic__elevation")
    )
    assert "topographic__elevation" in ds
    assert "surface_water__depth" not in ds


def test_exclude_single_field(grid_with_fields):
    ds = xr.open_dataset(
        dump(grid_with_fields, exclude="at_node:topographic__elevation")
    )
    assert "topographic__elevation" not in ds
    assert "surface_water__depth" in ds


@pytest.mark.parametrize("at", ("link", "face"))
def test_bad_location_raises(at):
    grid = RasterModelGrid((3, 4), xy_spacing=10.0)
    grid.add_ones("some__flux", at=at)
    with pytest.raises(ValueError, match="fields must be at one, and only one, of"):
        dump(grid)


def test_empty_fields_produces_no_data_vars():
    grid = RasterModelGrid((3, 4), xy_spacing=10.0)
    assert dump(grid) is None

    grid.add_ones("foo", at="node")
    assert dump(grid, include="at_cell*") is None


@pytest.mark.parametrize("at", ("node", "corner", "patch", "cell"))
def test_validation_returns_location(at):
    assert _validate_field_location([f"at_{at}:foo", f"at_{at}:bar"]) == at


@pytest.mark.parametrize("fields", (set(), [], ()))
def test_validation_no_fields_return_none(fields):
    assert _validate_field_location(fields) is None


@pytest.mark.parametrize("at", ("link", "face", "foo"))
def test_vaidation_bad_location_raises(at):
    with pytest.raises(ValueError, match="fields must be at one, and only one, of"):
        _validate_field_location([f"at_{at}:foobar"])


@pytest.mark.parametrize(
    "here, there", permutations(("node", "corner", "patch", "cell"), r=2)
)
def test_vaidation_multiple_locations_raises(here, there):
    with pytest.raises(ValueError, match="fields must be at one, and only one, of"):
        _validate_field_location([f"at_{here}:foo", f"at_{there}:bar"])
