import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.object.ext.at_node import reorder_rows
from landlab.graph.object.ext.at_node import reorder_rows_inplace


@pytest.mark.parametrize(
    "dtype", ("int8", "int", "long", "longlong", "float", "double")
)
def test_reorder_rows_dtype(dtype):
    value_at_row = np.asarray(
        [
            [0, 1, 2, 3],
            [0, 1, 2, 3],
            [0, 1, 2, 3],
            [0, 1, 2, 3],
        ],
        dtype=dtype,
    )
    sorted_cols = np.asarray(
        [
            [1, 2, 3, 0],
            [2, 3, 0, 1],
            [3, 0, 1, 2],
            [0, 1, 2, 3],
        ]
    )
    initial_array = value_at_row.copy()
    actual = reorder_rows(value_at_row, sorted_cols)

    assert_array_equal(value_at_row, initial_array)
    assert_array_equal(
        actual,
        [
            [1, 2, 3, 0],
            [2, 3, 0, 1],
            [3, 0, 1, 2],
            [0, 1, 2, 3],
        ],
    )


@pytest.mark.parametrize("n", (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
def test_reorder_rows(n):
    value_at_row = np.tile([0, 1, 2, 3], (n, 1))
    sorted_cols = value_at_row.copy()
    for row in sorted_cols:
        np.random.shuffle(row)

    actual = reorder_rows(value_at_row, sorted_cols)

    assert_array_equal(actual, sorted_cols)


@pytest.mark.parametrize("n", (0, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048))
def test_reorder_rows_inplace(n):
    value_at_row = np.tile([0, 1, 2, 3], (n, 1))
    sorted_cols = value_at_row.copy()
    for row in sorted_cols:
        np.random.shuffle(row)

    expected = reorder_rows(value_at_row, sorted_cols)
    reorder_rows_inplace(value_at_row, sorted_cols)

    assert_array_equal(value_at_row, sorted_cols)
    assert_array_equal(value_at_row, expected)
