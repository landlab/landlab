import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.core._validate import require_between
from landlab.core._validate import require_negative
from landlab.core._validate import require_nonnegative
from landlab.core._validate import require_nonpositive
from landlab.core._validate import require_one_of
from landlab.core._validate import require_positive
from landlab.core._validate import validate_array
from landlab.core.errors import ValidationError


@pytest.mark.parametrize(
    "value, allowed",
    (
        (0, (0, 1)),
        (1, ("0", 1)),
        ("foo", ("foo", "bar", "baz")),
        ("", ("", 0, True)),
        ("b", "foobar"),
        (0, ([1], [2], 0)),
        ([1], ([1], [2])),
    ),
)
@pytest.mark.parametrize("allowed_type", (list, tuple, set))
def test_require_one_of_ok(value, allowed, allowed_type):
    try:
        allowed_type(allowed)
    except TypeError:
        pytest.skip("allowed contains unhashable items")
    assert require_one_of(value, allowed=allowed_type(allowed)) == value


@pytest.mark.parametrize(
    "value, allowed",
    (
        (-1, (0, 1)),
        (0, ("0", 1)),
        ("fu", ("foo", "bar", "baz")),
        (1, ()),
        ("bar", "foobar"),
        (0, ([0],)),
        ([0], (0, 1)),
    ),
)
@pytest.mark.parametrize("allowed_type", (list, tuple, set))
def test_require_one_of_not_ok(value, allowed, allowed_type):
    try:
        allowed_type(allowed)
    except TypeError:
        pytest.skip("allowed contains unhashable items")
    with pytest.raises(ValidationError, match="invalid value"):
        require_one_of(value, allowed=allowed_type(allowed))


@pytest.mark.parametrize("value", (1.0, (-1, 1), np.asarray([-1, 1])))
@pytest.mark.parametrize(
    "kwds",
    (
        {"a_min": -2.0},
        {"a_min": -1.0, "inclusive_min": True},
        {"a_max": 2.0},
        {"a_max": 1.0, "inclusive_max": True},
        {"a_min": -2.0, "a_max": 2.0},
        {"a_min": -1.0, "a_max": 1.0, "inclusive_min": True, "inclusive_max": True},
    ),
)
def test_require_between_valid(value, kwds):
    assert require_between(value, **kwds) is value


@pytest.mark.parametrize("value", (1.0, (-1, 1), np.asarray([-1, 1])))
@pytest.mark.parametrize(
    "kwds",
    (
        {"a_min": 2.0},
        {"a_min": 1.0, "inclusive_min": False},
        {"a_max": 0.0},
        {"a_max": 1.0, "inclusive_max": False},
        {"a_min": -1.0, "a_max": 1.0, "inclusive_min": True, "inclusive_max": False},
        {"a_min": 1.0, "a_max": 1.0, "inclusive_min": False, "inclusive_max": True},
    ),
)
def test_require_require_is_invalid(value, kwds):
    with pytest.raises(ValidationError, match="value must be"):
        require_between(value, **kwds)


@pytest.mark.parametrize("value,ok", ((-1.0, True), (0.0, False), (1.0, False)))
def test_require_negative(value, ok):
    if ok:
        assert require_negative(value) is value
    else:
        with pytest.raises(ValidationError, match="value must be <"):
            require_negative(value)


@pytest.mark.parametrize("value,ok", ((1.0, True), (0.0, False), (-1.0, False)))
def test_require_positive(value, ok):
    if ok:
        assert require_positive(value) is value
    else:
        with pytest.raises(ValidationError, match="value must be >"):
            require_positive(value)


@pytest.mark.parametrize("value,ok", ((-1.0, False), (0.0, True), (1.0, True)))
def test_require_nonnegative(value, ok):
    if ok:
        assert require_nonnegative(value) is value
    else:
        with pytest.raises(ValidationError, match="value must be >="):
            require_nonnegative(value)


@pytest.mark.parametrize("value,ok", ((-1.0, True), (0.0, True), (1.0, False)))
def test_require_nonpositive(value, ok):
    if ok:
        assert require_nonpositive(value) is value
    else:
        with pytest.raises(ValidationError, match="value must be <="):
            require_nonpositive(value)


@pytest.mark.parametrize("array", ([1.0, 2.0], [[1, 2]], [[1, 2], [3, 4]], []))
def test_validate_array_noop(array):
    x = np.asarray(array)
    actual = validate_array(x)
    assert actual is x
    assert_array_equal(actual, array)


@pytest.mark.parametrize("array", ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0], []))
def test_validate_returns_same_object_when_valid(array):
    x = np.asarray(array)
    actual = validate_array(x, dtype=x.dtype, shape=x.shape)
    assert actual is x
    assert_array_equal(actual, array)


@pytest.mark.parametrize(
    "dtype",
    [np.float32, "float32", np.dtype("float32")],
)
@pytest.mark.parametrize("array", ([1.0, 2.0, 3.0, 4.0], []))
def test_validate_dtype_accepts_multiple_specifiers(array, dtype):
    x = np.asarray(array, dtype=np.float32)
    actual = validate_array(x, dtype=dtype, shape=x.shape)
    assert actual is x
    assert_array_equal(actual, array)


@pytest.mark.parametrize(
    "array,dtype",
    (
        ([1.0, 2.0, 3.0], int),
        ([1.0, 2.0, 3.0], bool),
        ([1.0, 2.0, 3.0], complex),
        ([1, 2, 3], float),
        ([1, 2, 3], complex),
        ([1, 2, 3], bool),
    ),
)
def test_validate_dtype_mismatch_raises(array, dtype):
    with pytest.raises(ValidationError, match="incorrect type"):
        validate_array(np.asarray(array), dtype=dtype)


@pytest.mark.parametrize(
    "array,shape",
    [
        ([0, 0, 0, 0, 0, 0], (2, 3)),
        ([[0, 0, 0], [0, 0, 0]], (3, 2)),
        ([[0, 0], [0, 0], [0, 0]], (6,)),
    ],
)
def test_validate_shape_mismatch(array, shape):
    with pytest.raises(ValidationError, match="incorrect shape"):
        validate_array(np.asarray(array), shape=shape)


def test_validate_requires_writable_passes_when_writable():
    x = np.zeros(5)
    actual = validate_array(x, writable=True)
    assert actual is x


def test_validate_requires_writable_raises_when_readonly():
    x = np.arange(5)
    x.setflags(write=False)
    with pytest.raises(ValidationError):
        validate_array(x, writable=True)


@pytest.mark.parametrize("array", ([1, 2, 3], [], [[1, 2, 3], [4, 5, 6]]))
def test_validate_contiguous_requirement_passes_for_c_contiguous(array):
    x = np.asarray(array).copy(order="C")
    assert x.flags.c_contiguous
    actual = validate_array(x, contiguous=True)
    assert_array_equal(actual, array)


@pytest.mark.parametrize(
    "array",
    (
        np.arange(12).reshape((3, 4)).T,
        np.arange(12).reshape((3, 4))[:, ::2],
    ),
)
def test_validate_contiguous_requirement_raises_for_noncontiguous(array):
    with pytest.raises(ValidationError):
        validate_array(array, contiguous=True)
