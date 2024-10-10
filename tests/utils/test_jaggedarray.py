import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.utils.jaggedarray import unravel


@pytest.mark.parametrize("dtype", (int,))
def test_unravel_jaggedarray(dtype):
    actual = unravel(
        np.arange(24, dtype=int),
        np.asarray([0, 4, 8, 12, 16, 20, 24], dtype=int),
    )

    assert_array_equal(actual, np.arange(24).reshape((6, 4)))
