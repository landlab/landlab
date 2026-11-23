import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.data_record.data_record import DataRecord
from landlab.data_record.data_record import MissingValue
from landlab.data_record.data_record import infer_fill_values
from landlab.data_record.data_record import norm_fill_values
from landlab.data_record.data_record import spec_from_value


def test_data_record():
    grid = RasterModelGrid((4, 5))
    items = {
        "grid_element": [["node"], ["link"]],
        "element_id": [[1], [3]],
    }
    data_vars = {
        "mean_elevation": (["time"], [110.0]),
        "item_size": (["item_id", "time"], [[0.3], [0.4]]),
    }
    data_record = DataRecord(
        grid=grid, time=0.0, items=items, data_vars=data_vars, fill_value=None
    )
    data_record.add_record(
        time=50.0,
        new_record={
            "int-var": (["time"], [120]),
            "float-var": (["time"], [120.0]),
            "foo": (["time"], [120]),
        },
        fill_value={"foo": -2},
    )
    assert_array_equal(data_record.dataset["int-var"].values, [-1, 120])
    assert_array_equal(data_record.dataset["float-var"].values, [np.nan, 120.0])
    assert_array_equal(data_record.dataset["foo"].values, [-2, 120])

    data_record.add_item(
        time=0.0,
        new_item={
            "grid_element": [["node"], ["node"]],
            "element_id": [[2], [5]],
        },
        new_item_spec={
            "size": (["item_id", "time"], [[10], [5]]),
        },
    )

    assert_array_equal(
        data_record.dataset["size"],
        [
            [-1, -1],
            [-1, -1],
            [10, -1],
            [5, -1],
        ],
    )
    assert_array_equal(data_record.dataset["foo"], [-2, 120])


def test_missing_value():
    spec = MissingValue(-1)
    assert_array_equal(spec.is_missing([0, -2, -1, 5]), [False, False, True, False])

    spec = MissingValue(-1, is_missing=lambda x: np.asarray(x) < 0)
    assert_array_equal(spec.is_missing([0, -2, -1, 5]), [False, True, True, False])


@pytest.mark.parametrize("val", (1, 0, [2], [[-1]], np.int8(1)))
def test_spec_from_value_int(val):
    spec = spec_from_value(val)
    assert spec.fill_value == -1
    assert_array_equal(spec.is_missing([0, -1, -2, 1]), [False, True, False, False])


@pytest.mark.parametrize("val", (1.0, 0.0, np.nan, np.inf, [2.0], [[-1.0]]))
def test_spec_from_value_float(val):
    spec = spec_from_value(val)
    assert np.isnan(spec.fill_value)
    assert_array_equal(spec.is_missing([0, np.nan, -2, 1]), [False, True, False, False])


@pytest.mark.parametrize("val", (True, False))
def test_spec_from_value_bool(val):
    spec = spec_from_value(val)
    assert not spec.fill_value
    assert_array_equal(
        spec.is_missing([True, False, True, True]), [False, True, False, False]
    )


@pytest.mark.parametrize("val", ("foo", "", ["foo", "bar", "foobar"]))
def test_spec_from_value_str(val):
    spec = spec_from_value(val)
    assert spec.fill_value == "nan"
    assert_array_equal(
        spec.is_missing(["foo", "nan", "bar", "bar"]), [False, True, False, False]
    )


def test_norm_fill_values_empty():
    assert norm_fill_values({}) == {}


@pytest.mark.parametrize("fill_value", ({}, {"foo": -1}, {"bar": MissingValue(np.nan)}))
def test_norm_fill_values_makes_instances(fill_value):
    spec = norm_fill_values(fill_value)
    assert np.all(isinstance(x, MissingValue) for x in spec.values())


@pytest.mark.parametrize(
    "name, value", (("foo", -2), ("element_id", -5), ("grid_element", 12))
)
def test_norm_fill_values(name, value):
    spec = norm_fill_values({name: value})
    assert set(spec) == {name}
    assert spec[name].fill_value == value


def test_infer_fill_values():
    spec = infer_fill_values({"foo": [1.0, -2.0], "bar": -2})
    assert np.isnan(spec["foo"].fill_value)
    assert spec["bar"].fill_value == -1
