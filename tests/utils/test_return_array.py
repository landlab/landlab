import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.field.errors import FieldError
from landlab.field.errors import GroupError
from landlab.field.graph_field import GraphFields
from landlab.utils.return_array import resolve_field
from landlab.utils.return_array import validate_field


@pytest.fixture
def fields():
    fields = GraphFields()
    fields.new_field_location("node", 4)
    fields.new_field_location("cell", 2)
    fields.add_field("elevation", np.array([0.0, 1.0, 2.0, 3.0]), at="node")
    return fields


class TestValidateField:
    def test_field_name_without_grid(self):
        assert validate_field("elevation") == "elevation"

    def test_field_name_with_grid(self, fields):
        assert validate_field("elevation", grid=fields, at="node") == "elevation"

    def test_missing_field_name_raises(self, fields):
        with pytest.raises(FieldError):
            validate_field("not_a_field", grid=fields, at="node")

    def test_grid_without_at_raises(self, fields):
        with pytest.raises(ValueError):
            validate_field("elevation", grid=fields)

    def test_at_without_grid_raises(self):
        with pytest.raises(ValueError):
            validate_field("elevation", at="node")

    def test_array_like_is_converted_to_ndarray(self, fields):
        rtn = validate_field([0, 1, 2, 3], grid=fields, at="node")
        assert_array_equal(rtn, [0, 1, 2, 3])
        assert isinstance(rtn, np.ndarray)

    def test_ndarray_is_not_copied(self, fields):
        value = np.array([0.0, 1.0, 2.0, 3.0])
        rtn = validate_field(value, grid=fields, at="node")
        assert rtn is value

    def test_wrong_first_dimension_raises(self, fields):
        with pytest.raises(ValueError):
            validate_field([0, 1, 2], grid=fields, at="node")

    def test_extra_trailing_dimensions_are_allowed(self, fields):
        value = np.zeros((4, 2))
        rtn = validate_field(value, grid=fields, at="node")
        assert rtn is value

    def test_shape_is_not_checked_without_grid(self):
        value = np.zeros((99, 3))
        rtn = validate_field(value)
        assert rtn is value

    def test_scalar_is_allowed_with_grid(self, fields):
        rtn = validate_field(4.0, grid=fields, at="node")
        assert rtn == 4.0

    def test_scalar_is_allowed_without_grid(self):
        rtn = validate_field(4.0)
        assert rtn == 4.0

    def test_scalar_with_unknown_group_raises(self, fields):
        with pytest.raises(GroupError):
            validate_field(4.0, grid=fields, at="not_a_group")


class TestResolveField:
    def test_field_name_resolves_to_field_values(self, fields):
        rtn = resolve_field("elevation", grid=fields, at="node")
        assert_array_equal(rtn, fields.field_values("elevation", at="node"))

    def test_missing_field_name_raises(self, fields):
        with pytest.raises(FieldError):
            resolve_field("not_a_field", grid=fields, at="node")

    def test_missing_group_raises(self, fields):
        with pytest.raises(GroupError):
            resolve_field("elevation", grid=fields, at="not_a_group")

    def test_array_is_returned_unchanged(self, fields):
        value = np.array([0.0, 1.0, 2.0, 3.0])
        rtn = resolve_field(value, grid=fields, at="node")
        assert rtn is value

    def test_scalar_is_broadcast_to_grid_size(self, fields):
        rtn = resolve_field(4.0, grid=fields, at="node")
        assert_array_equal(rtn, [4.0, 4.0, 4.0, 4.0])

    def test_scalar_broadcast_is_read_only(self, fields):
        rtn = resolve_field(4.0, grid=fields, at="node")
        with pytest.raises(ValueError):
            rtn[0] = 0.0

    def test_field_name_resolves_to_current_values(self, fields):
        """A field name is re-resolved on every call, picking up changes
        made to the field's values after it was validated.
        """
        field = validate_field("elevation", grid=fields, at="node")

        fields["node"]["elevation"][:] = -1.0

        assert_array_equal(
            resolve_field(field, grid=fields, at="node"), [-1.0, -1.0, -1.0, -1.0]
        )
