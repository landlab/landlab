import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph._reindexer import Reindexer
from landlab.graph._reindexer import sort_order_from_coords


@pytest.mark.parametrize("order", ([0, 1, 2, 3], [3, 2, 1, 0], []))
def test_order_length(order):
    r = Reindexer(order)
    assert len(r) == len(order)


@pytest.mark.parametrize("order", ([1, 0, 2], []))
def test_repr(order):
    r = Reindexer(order)
    assert repr(r) == f"Reindexer(size={len(order)})"


@pytest.mark.parametrize(
    "bad_order",
    [
        [0, 0, 2],  # duplicate
        [0, 1, 3],  # out of range
        [-1, 0, 1],  # negative
        [1, 2],  # missing element
    ],
)
def test_invalid_order_raises(bad_order):
    with pytest.raises(ValueError, match="^order must be a permutation"):
        Reindexer(bad_order)


@pytest.mark.parametrize(
    "order, array, expected",
    [
        ([0, 1, 2], [10, 20, 30], [10, 20, 30]),  # identity
        ([2, 0, 1], [10, 20, 30], [30, 10, 20]),  # rotation
        ([2, 1, 0], [10, 20, 30], [30, 20, 10]),  # reversal
    ],
)
def test_reorder_1d(order, array, expected):
    r = Reindexer(order)
    assert_array_equal(r.reorder(array), expected)


def test_reorder_2d():
    r = Reindexer([1, 0])
    array = np.array([[1, 2, 5], [3, 4, 0]])
    assert_array_equal(r.reorder(array), [[3, 4, 0], [1, 2, 5]])


def test_reorder_into_out():
    r = Reindexer([2, 0, 1])
    array = np.array([10, 20, 30])
    out = np.full(3, -999, dtype=int)
    result = r.reorder(array, out=out)
    assert result is out
    assert_array_equal(out, [30, 10, 20])


def test_reorder_out_wrong_shape_raises():
    r = Reindexer([1, 0])
    with pytest.raises(ValueError, match="^out must have the same shape"):
        r.reorder([10, 20], out=np.empty(3, dtype=int))


def test_reorder_out_shared_memory_raises():
    r = Reindexer([1, 0])
    array = np.array([10, 20, -1])
    with pytest.raises(ValueError, match="^out must not share memory"):
        r.reorder(array[:2], out=array[-2:])


def test_reorder_wrong_length_raises():
    r = Reindexer([1, 0])
    with pytest.raises(ValueError, match="^array must have length 2"):
        r.reorder([10, 20, 30])


@pytest.mark.parametrize(
    "values", ([], np.array([], dtype=float), np.array([], dtype=int))
)
def test_reorder_empty(values):
    r = Reindexer([])
    assert_array_equal(r.reorder(values), values)


@pytest.mark.parametrize(
    "order, array, expected",
    [
        ([0, 1, 2], [0, 1, 2], [0, 1, 2]),  # identity
        ([2, 0, 1], [0, 1, 2], [1, 2, 0]),  # rotation: old→new
        ([2, 1, 0], [0, 1, 2], [2, 1, 0]),  # reversal
    ],
)
def test_remap_1d(order, array, expected):
    r = Reindexer(order)
    assert_array_equal(r.remap(array), expected)


def test_remap_2d():
    r = Reindexer([1, 0, 2])
    array = np.array([[0, 1], [2, 0]])
    actual = r.remap(array)
    assert_array_equal(actual, [[1, 0], [2, 1]])


def test_remap_with_fill_value():
    r = Reindexer([2, 0, 1])
    array = np.array([0, -1, 1, -1, 2])
    actual = r.remap(array, fill_value=-1)
    assert_array_equal(actual, [1, -1, 2, -1, 0])


def test_remap_fill_value_preserved():
    r = Reindexer([1, 0])
    array = np.array([-1, -1])
    actual = r.remap(array, fill_value=-1)
    assert_array_equal(actual, [-1, -1])


def test_remap_into_out():
    r = Reindexer([1, 0, 2])
    array = np.array([0, 1, 2])
    out = np.full(3, -999, dtype=int)
    actual = r.remap(array, out=out)
    assert actual is out
    assert_array_equal(actual, [1, 0, 2])


def test_remap_out_shared_memory_raises():
    r = Reindexer([1, 0])
    array = np.array([0, 1])
    with pytest.raises(ValueError, match="^out must not share memory"):
        r.remap(array, out=array)


def test_remap_out_wrong_shape_raises():
    r = Reindexer([1, 0])
    array = np.array([0, 1])
    with pytest.raises(ValueError, match="^out must have the same shape"):
        r.remap(array, out=np.empty_like(array).reshape((2, -1)))


def test_remap_out_wrong_dtype_raises():
    r = Reindexer([1, 0])
    array = np.array([0, 1])
    with pytest.raises(ValueError, match="^out must have dtype"):
        r.remap(array, out=np.empty_like(array, dtype=float))


def test_remap_out_of_range_raises():
    r = Reindexer([1, 0])
    with pytest.raises(ValueError, match="^value must be <= 1"):
        r.remap([0, 5])


@pytest.mark.parametrize(
    "array", ([], np.array([], dtype=float), np.array([], dtype=int))
)
def test_remap_empty(array):
    r = Reindexer([])
    assert_array_equal(r.remap(array), array)


@pytest.mark.parametrize(
    "ids, expected_order",
    [
        ([0, 1, 2], [0, 1, 2]),  # already sorted
        ([2, 0, 1], [1, 2, 0]),  # argsort of ids
        ([3, 1, 2], [1, 2, 0]),  # general case
    ],
)
def test_from_ids_order(ids, expected_order):
    r = Reindexer.from_ids(ids)
    # reorder(ids) should give sorted ids
    assert_array_equal(r.reorder(ids), sorted(ids))


def test_from_ids_stable():
    # Ties resolved stably: first occurrence keeps lower new index
    ids = [1, 1, 0]
    r = Reindexer.from_ids(ids)
    result = r.reorder(ids)
    assert result[0] == 0  # the 0 comes first
    assert result[1] == result[2] == 1


def test_from_xy_sorts_by_y_then_x():
    # Default: primary key is last column (y), tie-break by x
    xy = np.array([[1.0, 0.0], [0.0, 1.0], [0.0, 0.0], [1.0, 1.0]])
    r = Reindexer.from_xy(xy)
    sorted_xy = r.reorder(xy)
    # Should be sorted by y first (ascending), then x
    assert_array_equal(sorted_xy[:, 1], sorted(xy[:, 1]))


def test_from_xy_axes_swap():
    xy = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    _ = Reindexer.from_xy(xy)  # primary key: col 1 (y)
    r_swapped = Reindexer.from_xy(xy, axes=(1, 0))  # primary key: col 0 (x)
    # With swapped axes, col-0 values should be non-decreasing
    sorted_xy = r_swapped.reorder(xy)
    assert all(sorted_xy[i, 0] <= sorted_xy[i + 1, 0] for i in range(len(xy) - 1))


def test_sort_order_from_coords_basic():
    xy = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    order = sort_order_from_coords(xy)
    assert_array_equal(order, [0, 2, 1, 3])


@pytest.mark.parametrize(
    "axes, expected",
    [
        ((0, 1), [0, 2, 1, 3]),  # primary key: col 1 (y)
        ((1, 0), [0, 1, 2, 3]),  # primary key: col 0 (x)
    ],
)
def test_sort_order_from_coords_axes(axes, expected):
    xy = np.array([[0.0, 0.0], [0.0, 1.0], [1.0, 0.0], [1.0, 1.0]])
    assert_array_equal(sort_order_from_coords(xy, axes=axes), expected)


def test_sort_order_from_coords_tol():
    # Two points differ by less than tol — treated as equal, sorted by other axis
    xy = np.array([[0.0, 0.0], [0.001, 0.0], [1.0, 0.0]])
    order = sort_order_from_coords(xy, tol=0.01)
    # First two x-values are within tol, so secondary key (y=0 for all) doesn't
    # distinguish them — but both come before x=1.0
    assert order[2] == 2


def test_sort_order_from_coords_not_2d_raises():
    with pytest.raises(ValueError, match="coords must be 2D"):
        sort_order_from_coords(np.array([0.0, 1.0, 2.0]))


def test_sort_order_from_coords_bad_axes_raises():
    xy = np.array([[0.0, 0.0], [1.0, 1.0]])
    with pytest.raises(ValueError, match="^axes must be a permutation of"):
        sort_order_from_coords(xy, axes=(0, 2))


def test_sort_order_from_coords_negative_tol_raises():
    xy = np.array([[0.0, 0.0], [1.0, 1.0]])
    with pytest.raises(ValueError, match="^tol must be > 0.0"):
        sort_order_from_coords(xy, tol=-0.1)
