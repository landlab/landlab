import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.core._validate import validate_array
from landlab.core.errors import ValidationError


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
