#! /usr/bin/env python
import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.field.graph_field import FieldDataset


def test_init():
    ds = FieldDataset("node")
    assert ds.size is None
    assert ds.fixed_size is True
    assert list(ds) == []
    assert len(ds) == 0


def test_init_with_size():
    ds = FieldDataset("node", size=1)
    assert ds.size == 1
    with pytest.raises(ValueError):
        FieldDataset("node", size=1.1)
    with pytest.raises(ValueError):
        FieldDataset("node", size=-1)


def test_init_fixed_size():
    assert FieldDataset("node", fixed_size=True).fixed_size is True
    assert FieldDataset("node", fixed_size=1).fixed_size is True
    assert FieldDataset("node", fixed_size=False).fixed_size is False
    assert FieldDataset("node", fixed_size=0).fixed_size is False
    assert FieldDataset("node", fixed_size="true").fixed_size is True
    assert FieldDataset("node", fixed_size="false").fixed_size is True


def test_size_setter():
    ds = FieldDataset("node", fixed_size=True)
    assert ds.size is None

    ds.size = 3
    assert ds.size == 3
    with pytest.raises(ValueError):
        ds.size = 4
    ds.size = 3

    ds.fixed_size = False
    ds.size = 4
    assert ds.size == 4

    with pytest.raises(ValueError):
        ds.size = 1.1
    with pytest.raises(ValueError):
        ds.size = -1


def test_add_field():
    ds = FieldDataset("node")
    assert ds.size is None
    ds.set_value("air_temperature", [0, 0, 0])
    assert ds.size == 3
    assert_array_equal(ds["air_temperature"], (0, 0, 0))
    assert len(ds) == 1
    assert list(ds) == ["air_temperature"]


def test_field_maintains_reference():
    ds = FieldDataset("node")
    array = np.arange(3)
    ds.set_value("air_temperature", array)
    assert ds["air_temperature"] is array


def test_setitem():
    ds = FieldDataset("node")
    assert ds.size is None
    ds["air_temperature"] = [0, 0, 0]
    assert ds.size == 3
    assert_array_equal(ds["air_temperature"], (0, 0, 0))
    assert len(ds) == 1
    assert list(ds) == ["air_temperature"]


def test_add_fields():
    ds = FieldDataset("node", fixed_size=True)
    ds.set_value("air_temperature", [0, 0, 0])
    ds.set_value("ground_temperature", [1, 1, 1])
    assert_array_equal(ds["air_temperature"], (0, 0, 0))
    assert_array_equal(ds["ground_temperature"], (1, 1, 1))
    assert len(ds) == 2
    assert sorted(list(ds)) == ["air_temperature", "ground_temperature"]


def test_add_fields_fixed_size():
    ds = FieldDataset("node", fixed_size=True)
    ds.set_value("air_temperature", [0, 0, 0])
    with pytest.raises(ValueError):
        ds.set_value("ground_temperature", [1, 1, 1, 1])


def test_add_fields_adjustable_size():
    ds = FieldDataset("node", fixed_size=False)
    ds.set_value("air_temperature", [0, 0, 0])
    assert ds.size == 3
    ds.set_value("air_temperature", [1, 1, 1, 1])
    assert ds.size == 4
    with pytest.raises(ValueError):
        ds.set_value("ground_temperature", [1, 1, 1])


def test_add_fields_fixed_size_setter():
    ds = FieldDataset("node", fixed_size=False)
    ds.set_value("air_temperature", [0, 0, 0])
    ds.set_value("air_temperature", [1, 1, 1, 1])
    ds.fixed_size = True
    with pytest.raises(ValueError):
        ds.set_value("air_temperature", [1, 1, 1])
    assert ds.size == 4
    ds.fixed_size = False
    ds.set_value("air_temperature", [1, 1, 1])
    assert ds.size == 3


def test_pop():
    ds = FieldDataset("node")
    ds.set_value("air_temperature", [0, 0, 0])
    ds.set_value("ground_temperature", [1, 1, 1])
    assert sorted(list(ds)) == ["air_temperature", "ground_temperature"]

    assert_array_equal(ds.pop("ground_temperature"), [1, 1, 1])
    assert sorted(list(ds)) == ["air_temperature"]

    assert_array_equal(ds.pop("air_temperature"), [0, 0, 0])
    assert sorted(list(ds)) == []

    with pytest.raises(KeyError):
        ds.pop("air_temperature")


def test_contains():
    ds = FieldDataset("node")
    ds.set_value("air_temperature", [0, 0, 0])
    assert "air_temperature" in ds
    assert "ground_temperature" not in ds
    assert "node" not in ds


def test_scalars():
    ds = FieldDataset("node")
    ds["air_temperature"] = 2
    assert ds["air_temperature"] == 2
    assert ds.size == 1

    ds["ground_temperature"] = [0]
    assert ds["ground_temperature"] == [0]


def test_size_is_one():
    ds = FieldDataset("grid")
    ds.size = 1
    ds["air_temperature"] = [1, 2, 3]
    assert_array_equal(ds["air_temperature"], [[1, 2, 3]])
    ds["ground_temperature"] = [0, 0, 0, 0, 0]
    assert_array_equal(ds["ground_temperature"], [[0, 0, 0, 0, 0]])
    ds["gravity"] = 9.81
    assert ds["gravity"] == pytest.approx(9.81)


def test_reshaping_arrays():
    ds = FieldDataset("links", size=3, fixed_size=True)
    ds["air_temperature"] = [1.0, 2.0, 3.0]
    assert_array_equal(ds["air_temperature"], [1.0, 2.0, 3.0])

    assert ds.size == 3
    ds["nodes"] = [[0, 0], [0, 0], [0, 0]]
    assert_array_equal(ds["nodes"], [[0, 0], [0, 0], [0, 0]])

    ds["nodes"] = [1, 2, 3, 4, 5, 6]
    assert_array_equal(ds["nodes"], [[1, 2], [3, 4], [5, 6]])

    with pytest.raises(ValueError):
        ds["nodes"] = [1, 2, 3, 4, 5, 6, 7]
