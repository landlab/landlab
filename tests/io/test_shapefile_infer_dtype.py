import hypothesis.extra.numpy as hynp
import numpy as np
import pytest
from hypothesis import assume
from hypothesis import given
from hypothesis.strategies import integers
from numpy.testing import assert_array_equal

from landlab.io.shapefile import _infer_data_type


@pytest.mark.parametrize("src_type", [np.int32, np.int64, float, bool])
@pytest.mark.parametrize("as_iterable", [list, tuple, np.asarray])
def test_infer_dtype(as_iterable, src_type):
    values = np.asarray([1, 2, 3], dtype=src_type)
    array = _infer_data_type(as_iterable(values))
    assert array.dtype == src_type
    assert_array_equal(array, values)


@pytest.mark.parametrize(
    "values,dst_type",
    [
        ([1, 2, 3], np.int_),
        ([1, 2, 3.0], float),
        ([1.0, 2.0, 3.0], float),
        ([1j, 2j, 3j], complex),
        ([True, True, False], bool),
    ],
)
def test_infer_dtype_from_object(values, dst_type):
    values = np.asarray(values, dtype=object)
    array = _infer_data_type(values)
    assert array.dtype == dst_type
    assert_array_equal(array, values)


@pytest.mark.parametrize(
    "values,dst_type,expected",
    [
        ([1, 2, 3.0], float, [1.0, 2.0, 3.0]),
        ([1.0, 2.0, 3j], complex, [1 + 0j, 2 + 0j, 3j]),
        ([1, 2, 3j], complex, [1 + 0j, 2 + 0j, 3j]),
        ([True, False, 1], int, [1, 0, 1]),
        ([True, False, 1.0], float, [1.0, 0.0, 1.0]),
        ([None, 1.0, 2.0], float, [np.nan, 1.0, 2.0]),
        ([None, 1.0, 2j], complex, [np.nan * 1j, 1.0, 2j]),
    ],
)
def test_infer_dtype_from_mixed(values, dst_type, expected):
    values = np.asarray(values, dtype=object)
    array = _infer_data_type(values)
    assert array.dtype == dst_type
    assert_array_equal(array, expected)


@given(
    values=hynp.arrays(
        dtype=hynp.floating_dtypes()
        | hynp.integer_dtypes()
        | hynp.complex_number_dtypes(),
        shape=hynp.array_shapes(),
        elements=integers(0, 2**7 - 1),
    ),
    dst_type=hynp.floating_dtypes()
    | hynp.integer_dtypes()
    | hynp.complex_number_dtypes(),
)
@pytest.mark.parametrize("as_iterable", [list, tuple, np.asarray])
def test_dtype_keyword(as_iterable, values, dst_type):
    array = _infer_data_type(as_iterable(values), dtype=dst_type)
    assert array.dtype == dst_type
    assert_array_equal(array, values)


@given(
    values=hynp.arrays(
        dtype=hynp.floating_dtypes()
        | hynp.integer_dtypes()
        | hynp.complex_number_dtypes(),
        shape=hynp.array_shapes(),
        elements=integers(0, 2**7 - 1),
    ),
    dst_type=hynp.floating_dtypes()
    | hynp.integer_dtypes()
    | hynp.complex_number_dtypes(),
)
def test_object_arrays(values, dst_type):
    assume(
        not np.issubdtype(values.dtype, np.complexfloating)
        or np.issubdtype(dst_type, np.complexfloating)
    )
    values = np.asarray(values, dtype=object)
    array = _infer_data_type(values, dtype=dst_type)
    assert array.dtype == dst_type
    assert_array_equal(array, values)


@given(
    values=hynp.arrays(
        dtype=hynp.complex_number_dtypes(),
        shape=hynp.array_shapes(),
        elements=integers(0, 8),
    ),
    dst_type=hynp.floating_dtypes() | hynp.integer_dtypes(),
)
def test_from_complex(values, dst_type):
    array = _infer_data_type(values, dtype=dst_type)
    assert_array_equal(array.real, values.real)


@given(values=hynp.arrays(dtype=hynp.array_dtypes(), shape=hynp.array_shapes()))
def test_array_not_copied(values):
    array = _infer_data_type(values, dtype=values.dtype)
    assert array is values
