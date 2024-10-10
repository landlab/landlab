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
