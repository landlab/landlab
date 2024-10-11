import numpy as np
import pytest
from numpy.testing import assert_array_equal

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
