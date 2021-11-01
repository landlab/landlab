import numpy as np
import pytest
from numpy.testing import assert_array_equal


from landlab.io.shapefile.read_shapefile import _infer_data_type


@pytest.mark.parametrize("src_type", [np.int32, np.int64, float, bool])
@pytest.mark.parametrize("as_iterable", [list, tuple, np.asarray])
def test_infer_dtype(as_iterable, src_type):
    values = np.asarray([1, 2, 3], dtype=src_type)
    array = _infer_data_type(as_iterable(values))
    assert array.dtype == src_type
    assert_array_equal(array, values)


@pytest.mark.parametrize("as_iterable", [list, tuple, np.asarray])
@pytest.mark.parametrize(
    "values,dst_type",
    [([1, 2, 3], np.int_), ([1.0, 2.0, 3.0], float), ([True, True, False], bool)],
)
def test_infer_dtype_from_object(as_iterable, values, dst_type):
    values = np.asarray(values, dtype=object)
    array = _infer_data_type(as_iterable(values))
    assert array.dtype == dst_type
    assert_array_equal(array, values)


@pytest.mark.parametrize("dst_type", [int, float, np.int32, np.int64])
@pytest.mark.parametrize("src_type", [int, float, object])
@pytest.mark.parametrize("as_iterable", [list, tuple, np.asarray])
def test_dtype_keyword(as_iterable, src_type, dst_type):
    values = np.asarray([1, 2, 3], dtype=src_type)
    array = _infer_data_type(as_iterable(values), dtype=dst_type)
    assert array.dtype == dst_type
    assert_array_equal(array, [1, 2, 3])
