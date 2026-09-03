import numpy as np
import pytest
from numpy.testing import assert_array_equal
from requireit import ValidationError

from landlab.utils.jaggedarray import JaggedArray
from landlab.utils.jaggedarray import padded_row_contains
from landlab.utils.jaggedarray import unravel


@pytest.mark.parametrize("dtype0", (int, np.intp, np.long, np.longlong, float))
@pytest.mark.parametrize("dtype1", (int, np.intp, np.long, np.longlong))
def test_unravel_jaggedarray(dtype0, dtype1):
    data = np.empty(24, dtype=dtype0)
    data[:] = np.arange(24)
    actual = unravel(data, np.asarray([0, 4, 8, 12, 16, 20, 24], dtype=dtype1))
    expected = data.copy().reshape((6, 4))

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("dtype0", (int, np.intp, np.long, np.longlong, float))
@pytest.mark.parametrize("dtype1", (int, np.intp, np.long, np.longlong))
def test_unravel_jaggedarray_with_padding(dtype0, dtype1):
    data = np.empty(24, dtype=dtype0)
    data[:] = np.arange(24)

    actual = unravel(data, np.asarray([0, 2, 6, 11, 17, 21, 24], dtype=dtype1), pad=-1)
    expected = np.asarray(
        [
            [0, 1, -1, -1, -1, -1],
            [2, 3, 4, 5, -1, -1],
            [6, 7, 8, 9, 10, -1],
            [11, 12, 13, 14, 15, 16],
            [17, 18, 19, 20, -1, -1],
            [21, 22, 23, -1, -1, -1],
        ],
        dtype=dtype0,
    )

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("dtype", (int, np.intp, np.long, np.longlong))
def test_padded_row_contains(dtype):
    rows = [[-1, 2, 3], [4, -1, 6, 7], [], [8, 9]]
    array = JaggedArray(rows)
    length_of_row = np.diff(array.offset)
    padded_array = unravel(array.array, array.offset, pad=-1).astype(dtype)

    actual = padded_row_contains(padded_array, length_of_row, value=-1)
    assert_array_equal(actual, [True, True, False, False])


def test_padded_row_contains_with_out():
    rows = [[-1, 2, 3], [4, -1, 6, 7], [], [8, 9]]
    array = JaggedArray(rows)
    length_of_row = np.diff(array.offset)
    padded_array = unravel(array.array, array.offset, pad=-1)

    out = np.ones(len(rows), dtype=bool)
    actual = padded_row_contains(padded_array, length_of_row, value=-1, out=out)
    assert actual is out
    assert_array_equal(actual, [True, True, False, False])


@pytest.mark.parametrize("dtype", (float, int, np.intp, np.long, np.longlong, np.uint8))
def test_padded_row_contains_bad_out_type(dtype):
    rows = [[-1, 2, 3], [4, -1, 6, 7], [], [8, 9]]
    array = JaggedArray(rows)
    length_of_row = np.diff(array.offset)
    padded_array = unravel(array.array, array.offset, pad=-1)

    out = np.ones(len(rows), dtype=dtype)
    with pytest.raises(ValidationError, match="out must have dtype"):
        padded_row_contains(padded_array, length_of_row, value=-1, out=out)
