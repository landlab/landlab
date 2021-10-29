from functools import partial

import numpy as np
import pytest
from numpy.testing import assert_array_equal


from landlab.io.shapefile.read_shapefile import _infer_data_type


@pytest.mark.parametrize("in_type", [int, float])
@pytest.mark.parametrize(
    "astype", [list, tuple, np.asarray, partial(np.asarray, dtype=object)]
)
def test_infer_dtype(astype, in_type):
    values = np.asarray([1, 2, 3], dtype=in_type)
    array = _infer_data_type(astype(values))
    assert array.dtype == in_type
    assert_array_equal(array, values)


@pytest.mark.parametrize("dtype", [int, float, np.int32, np.int64])
@pytest.mark.parametrize("in_type", [int, float, object])
@pytest.mark.parametrize(
    "astype", [list, tuple, np.asarray, partial(np.asarray, dtype=object)]
)
def test_dtype_keyword(astype, in_type, dtype):
    values = np.asarray([1, 2, 3], dtype=in_type)
    array = _infer_data_type(astype(values), dtype=dtype)
    assert array.dtype == dtype
    assert_array_equal(array, [1, 2, 3])
