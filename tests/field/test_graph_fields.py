#! /usr/bin/env python
import numpy as np
import pytest
from hypothesis import given
from hypothesis.strategies import text
from numpy.testing import assert_array_equal

from landlab.field import FieldError
from landlab.field import GroupError
from landlab.field.graph_field import GraphFields as ModelDataFields


def test_init():
    fields = ModelDataFields()
    assert set() == fields.groups
    # assert_set_equal(set(), fields.groups)


def test_new_field_location():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    assert {"node"} == fields.groups
    # assert_set_equal(set(['node']), fields.groups)


def test_add_existing_group():
    fields = ModelDataFields()
    fields.new_field_location("node", size=12)
    with pytest.raises(ValueError):
        fields.new_field_location("node", size=24)


def test_add_multiple_groups():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.new_field_location("cell", 2)
    fields.new_field_location("face", 7)
    fields.new_field_location("link", 7)
    assert {"node", "cell", "face", "link"} == fields.groups
    # assert_set_equal(set(['node', 'cell', 'face', 'link']), fields.groups)


def test_ones():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.new_field_location("cell", 2)

    value_array = fields.ones("node")
    assert_array_equal(np.ones(12), value_array)

    value_array = fields.ones("cell")
    assert_array_equal(np.ones(2), value_array)


def test_add_ones():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.new_field_location("cell", 2)

    fields.add_ones("z", at="node")
    assert_array_equal(np.ones(12), fields["node"]["z"])
    assert_array_equal(np.ones(12), fields.field_values("node", "z"))

    fields.add_ones("z", at="cell")
    assert_array_equal(np.ones(2), fields["cell"]["z"])
    assert_array_equal(np.ones(2), fields.field_values("cell", "z"))


def test_add_ones_return_value():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.new_field_location("cell", 2)

    rtn_value = fields.add_ones("z", at="node")
    assert_array_equal(rtn_value, np.ones(12))
    assert rtn_value is fields["node"]["z"]
    assert rtn_value is fields.field_values("node", "z")

    rtn_value = fields.add_ones("z", at="cell")
    assert_array_equal(rtn_value, np.ones(2))
    assert rtn_value is fields["cell"]["z"]
    assert rtn_value is fields.field_values("cell", "z")


def test_add_existing_field_default():
    """Test default is to not replace existing field."""
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.add_empty("z", at="node")

    with pytest.raises(FieldError):
        fields.add_empty("z", at="node")
    with pytest.raises(FieldError):
        fields.add_ones("z", at="node")
    with pytest.raises(FieldError):
        fields.add_zeros("z", at="node")


def test_add_existing_field_with_clobber_off():
    """Test clobber raises an error with an existing field."""
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.add_empty("z", at="node")

    with pytest.raises(FieldError):
        fields.add_empty("z", at="node", clobber=False)
    with pytest.raises(FieldError):
        fields.add_ones("z", at="node", clobber=False)
    with pytest.raises(FieldError):
        fields.add_zeros("z", at="node", clobber=False)


def test_add_field_with_clobber_off():
    """Test clobber does not raise an error with an new field."""
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    fields.add_empty("a", at="node", clobber=False)
    assert "a" in fields["node"]

    fields.add_ones("b", at="node", clobber=False)
    assert "b" in fields["node"]

    fields.add_zeros("c", at="node", clobber=False)
    assert "c" in fields["node"]


def test_add_field_with_clobber():
    """Test adding a field with clobber on."""
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    assert fields.add_empty("a", at="node") is not fields.add_empty(
        "a", at="node", clobber=True
    )
    assert fields.add_ones("b", at="node") is not fields.add_ones(
        "b", at="node", clobber=True
    )
    assert fields.add_zeros("c", at="node") is not fields.add_zeros(
        "c", at="node", clobber=True
    )


def test_getitem():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    assert fields["node"] == {}
    with pytest.raises(GroupError):
        fields["cell"]
    with pytest.raises(KeyError):
        fields["cell"]


def test_at_attribute():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    assert fields.at_node == {}
    with pytest.raises(AttributeError):
        fields.at_cell

    fields.add_ones("z", at="node")
    assert_array_equal(np.ones(12), fields.at_node["z"])


def test_has_group():
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    assert fields.has_group("node")
    assert not fields.has_group("cell")


def test_delete_field():
    fields = ModelDataFields()
    fields.new_field_location("link", 17)

    assert fields.at_link == {}
    with pytest.raises(AttributeError):
        fields.at_node

    fields.add_zeros("vals", at="link")
    assert_array_equal(np.zeros(17), fields.at_link["vals"])

    with pytest.raises(KeyError):
        fields.delete_field("node", "vals")
    fields.delete_field("link", "vals")
    with pytest.raises(KeyError):
        fields.field_units("link", "vals")
    with pytest.raises(KeyError):
        fields.at_link["vals"]


def test_scalar_field():
    """Test adding a generic scalar field."""
    fields = ModelDataFields()
    fields.new_field_location("all_over_the_place", 1)

    assert fields.at_all_over_the_place == {}
    with pytest.raises(AttributeError):
        fields.at_cell

    fields.at_all_over_the_place["const"] = 1.0
    assert_array_equal(np.array(1.0), fields.at_all_over_the_place["const"])

    val = np.array(2.0)
    fields.at_all_over_the_place["const"] = val
    assert val is fields.at_all_over_the_place["const"]


def test_grid_field_as_array():
    """Test adding an array as a grid field."""
    fields = ModelDataFields()
    fields.new_field_location("grid", 1)

    fields.at_grid["const"] = [1.0, 2.0]
    assert_array_equal(fields.at_grid["const"], [[1.0, 2.0]])

    val = np.array([1.0, 2.0])
    fields.at_grid["const"] = val
    assert np.shares_memory(val, fields.at_grid["const"])

    val.shape = (1, 1, 2, 1)
    fields.at_grid["const"] = val
    assert_array_equal(fields.at_grid["const"], np.array([[1.0, 2.0]]))
    assert np.shares_memory(val, fields.at_grid["const"])


def test_grid_field_add_zeros_ones_empty():
    """Test creating scalar fields with add_zeros, add_empty, and add_ones."""
    fields = ModelDataFields()
    fields.new_field_location("grid", 1)

    with pytest.raises(ValueError):
        fields.add_zeros("value", at="grid")
    with pytest.raises(ValueError):
        fields.add_empty("value", at="grid")
    with pytest.raises(ValueError):
        fields.add_ones("value", at="grid")


def test_grid_field_zeros_ones_empty():
    """Test creating scalar fields with zeros, empty, and ones."""
    fields = ModelDataFields()
    fields.new_field_location("grid", 1)
    with pytest.raises(ValueError):
        fields.zeros("grid")
    with pytest.raises(ValueError):
        fields.empty("grid")
    with pytest.raises(ValueError):
        fields.ones("grid")


def test_nd_field():
    """Test creating fields that are nd in shape."""
    fields = ModelDataFields()
    fields.new_field_location("node", 12)

    fields.add_field("new_value", np.ones((12, 4, 5)), at="node")
    fields.add_field("newer_value", np.ones((12, 4)), at="node")

    with pytest.raises(ValueError):
        fields.add_field("newest_value", np.ones((13, 4, 5)), at="node")
    with pytest.raises(ValueError):
        fields.add_field("newestest_value", np.ones(13), at="node")


@given(name=text(), unit_str=text())
def test_setting_units(name, unit_str):
    fields = ModelDataFields()
    fields.new_field_location("node", 12)
    fields.add_field(name, np.empty(12), at="node", units=unit_str[::-1])
    assert fields.field_units("node", name) == unit_str[::-1]

    fields["node"].set_units(name, unit_str)
    assert fields.field_units("node", name) == unit_str
