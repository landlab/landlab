import math

import numpy as np
import pytest
from numpy.testing import assert_array_equal
from requireit import ValidationError

from landlab.core.component_utils import iter_adaptive_time_steps
from landlab.core.component_utils import iter_time_steps
from landlab.core.component_utils import resolve_field
from landlab.core.component_utils import validate_field
from landlab.field.errors import FieldError
from landlab.field.errors import GroupError
from landlab.field.graph_field import GraphFields


@pytest.mark.parametrize(
    ("duration", "dt", "expected"),
    [
        (0.0, None, []),
        (5.0, None, [5.0]),
        (10.0, 2.5, [2.5, 2.5, 2.5, 2.5]),
        (10.0, 3.0, [2.5, 2.5, 2.5, 2.5]),
        (5.0, 10.0, [5.0]),
    ],
)
def test_iter_time_steps(duration, dt, expected):
    assert list(iter_time_steps(duration, dt=dt)) == expected


@pytest.mark.parametrize("duration", [-1.0, math.inf, math.nan])
def test_iter_time_steps_rejects_invalid_duration(duration):
    with pytest.raises(ValidationError):
        list(iter_time_steps(duration))


@pytest.mark.parametrize("dt", [-1.0, 0.0, math.inf, math.nan])
def test_iter_time_steps_rejects_invalid_dt(dt):
    with pytest.raises(ValidationError):
        list(iter_time_steps(1.0, dt=dt))


def test_iter_adaptive_time_steps_caps_final_step():
    expected = [3.0, 3.0, 3.0, 1.0]
    actual = list(iter_adaptive_time_steps(10.0, calc_dt=lambda: 3.0))
    assert actual == expected


def test_iter_adaptive_time_steps_stops_on_none():
    steps = iter([2.0, 2.0, None])

    expected = [2.0, 2.0]
    actual = list(iter_adaptive_time_steps(10.0, calc_dt=lambda: next(steps)))
    assert actual == expected


def test_iter_adaptive_time_steps_accepts_infinite_step():
    assert list(iter_adaptive_time_steps(10.0, calc_dt=lambda: math.inf)) == [10.0]


def test_iter_adaptive_time_steps_does_not_calculate_step_for_zero_duration():
    def calc_dt():
        raise AssertionError("calc_dt should not be called")

    with pytest.raises(AssertionError, match="calc_dt should not be called"):
        list(iter_adaptive_time_steps(10.0, calc_dt=calc_dt))
    assert list(iter_adaptive_time_steps(0.0, calc_dt=calc_dt)) == []


def test_iter_adaptive_time_steps_stops_within_tolerance():
    step = 1.0 - 5.0e-13

    assert list(iter_adaptive_time_steps(1.0, calc_dt=lambda: step)) == [step]


def test_iter_adaptive_time_steps_honors_zero_tolerance():
    steps = list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0 - 5.0e-13, rtol=0.0))

    assert len(steps) == 2
    assert sum(steps) == 1.0


@pytest.mark.parametrize("max_steps", (1, np.int64(1), np.int32(1)))
def test_iter_adaptive_time_steps_accepts_int(max_steps):
    actual = list(
        iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, max_steps=max_steps)
    )
    assert actual == [1.0]


def test_iter_adaptive_time_steps_raises_when_max_steps_exceeded():
    with pytest.raises(RuntimeError, match="unable to reach 2.0 in 1 substeps"):
        list(iter_adaptive_time_steps(2.0, calc_dt=lambda: 1.0, max_steps=1))


@pytest.mark.parametrize("max_steps", [0, -1, 1.5])
def test_iter_adaptive_time_steps_rejects_invalid_max_steps(max_steps):
    with pytest.raises(ValidationError):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, max_steps=max_steps))


@pytest.mark.parametrize("rtol", [-1.0, 1.0, math.inf, math.nan])
def test_iter_adaptive_time_steps_rejects_invalid_rtol(rtol):
    with pytest.raises(ValidationError):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: 1.0, rtol=rtol))


@pytest.mark.parametrize("step", [-1.0, 0.0, math.nan])
def test_iter_adaptive_time_steps_rejects_invalid_step(step):
    with pytest.raises(RuntimeError, match="step must be positive"):
        list(iter_adaptive_time_steps(1.0, calc_dt=lambda: step))


def test_iter_adaptive_time_steps_detects_step_too_small_to_advance():
    steps = iter([1.0, 1.0e-20])

    with pytest.raises(RuntimeError, match="too small relative to the elapsed time"):
        list(iter_adaptive_time_steps(2.0, calc_dt=lambda: next(steps), rtol=0.0))


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
