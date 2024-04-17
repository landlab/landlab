import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab.graph.quantity.ext.of_element import count_of_children_at_parent
from landlab.graph.quantity.ext.of_element import diff_of_children_at_parent
from landlab.graph.quantity.ext.of_element import max_of_children_at_parent
from landlab.graph.quantity.ext.of_element import mean_of_children_at_parent
from landlab.graph.quantity.ext.of_element import min_of_children_at_parent


def numpy_count(elements_at_element, out):
    out[:] = np.sum(elements_at_element != -1, axis=1)


def numpy_mean(elements_at_element, values, out):
    np.mean(
        values[elements_at_element], axis=1, where=elements_at_element != -1, out=out
    )


def numpy_min(elements_at_element, values, out):
    np.min(
        values[elements_at_element],
        axis=1,
        where=elements_at_element != -1,
        out=out,
        initial=np.inf,
    )


def numpy_max(elements_at_element, values, out):
    np.max(
        values[elements_at_element],
        axis=1,
        where=elements_at_element != -1,
        out=out,
        initial=-np.inf,
    )


def numpy_diff(children_at_parent, value_at_parent, value_at_child, out):
    np.subtract(
        value_at_child[children_at_parent],
        value_at_parent.reshape((-1, 1)),
        where=children_at_parent != -1,
        out=out,
    )

    return out


@pytest.mark.parametrize("impl", ["numpy", "cython"])
@pytest.mark.benchmark(group="diff")
def test_diff_of_elements_bench(benchmark, impl):
    rng = np.random.default_rng(seed=1973)
    elements_at_element = rng.integers(-1, 1000, size=(10000, 10))

    value_at_child = np.arange(1000, dtype=float)
    value_at_parent = np.arange(10000, dtype=float)
    actual = np.full_like(elements_at_element, -999, dtype=float)

    if impl == "numpy":
        func = numpy_diff
    else:
        func = diff_of_children_at_parent

    benchmark(
        func, np.asarray(elements_at_element), value_at_parent, value_at_child, actual
    )

    expected = np.full_like(actual, -999)
    numpy_diff(
        np.asarray(elements_at_element), value_at_parent, value_at_child, expected
    )

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("impl", ["numpy", "cython"])
def test_mean_of_elements(impl):
    elements_at_element = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
    ]
    values = np.arange(np.max(elements_at_element) + 1, dtype=float)
    actual = np.empty(len(elements_at_element), dtype=float)
    expected = [1.5, 5.5, 9.5]

    if impl == "numpy":
        func = numpy_mean
    else:
        func = mean_of_children_at_parent

    func(np.asarray(elements_at_element), values, actual)

    assert_array_equal(actual, expected)


def test_min_of_elements():
    elements_at_element = [
        [0, 1, 2, 3],
        [4, 5, 6, -1],
        [8, -1, 10, -1],
        [-1, -1, -1, -1],
    ]
    values = np.arange(np.max(elements_at_element) + 1, dtype=float)
    actual = np.full(len(elements_at_element), -999, dtype=float)
    expected = [0, 4, 8, -999]

    min_of_children_at_parent(np.asarray(elements_at_element), values, actual)

    assert_array_equal(actual, expected)


def test_max_of_elements():
    elements_at_element = [
        [0, 1, 2, 3],
        [4, 5, 6, -1],
        [8, -1, 10, -1],
        [-1, -1, -1, -1],
    ]
    values = np.arange(np.max(elements_at_element) + 1, dtype=float)
    actual = np.full(len(elements_at_element), -999, dtype=float)
    expected = [3, 6, 10, -999]

    max_of_children_at_parent(np.asarray(elements_at_element), values, actual)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("impl", ["numpy", "cython"])
@pytest.mark.benchmark(group="min")
def test_min_of_elements_bench(benchmark, impl):
    rng = np.random.default_rng(seed=1973)
    elements_at_element = rng.integers(-1, 1000, size=(10000, 10))

    values = np.arange(1000, dtype=float)
    actual = np.full(len(elements_at_element), -999, dtype=float)

    if impl == "numpy":
        func = numpy_min
    else:
        func = min_of_children_at_parent

    benchmark(func, np.asarray(elements_at_element), values, actual)

    expected = np.full_like(actual, -999)
    numpy_min(np.asarray(elements_at_element), values, expected)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("impl", ["numpy", "cython"])
@pytest.mark.benchmark(group="max")
def test_max_of_elements_bench(benchmark, impl):
    rng = np.random.default_rng(seed=1973)
    elements_at_element = rng.integers(-1, 1000, size=(10000, 10))

    values = np.arange(1000, dtype=float)
    actual = np.full(len(elements_at_element), -999, dtype=float)

    if impl == "numpy":
        func = numpy_max
    else:
        func = max_of_children_at_parent

    benchmark(func, np.asarray(elements_at_element), values, actual)

    expected = np.full_like(actual, -999)
    numpy_max(np.asarray(elements_at_element), values, expected)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("impl", ["numpy", "cython"])
@pytest.mark.benchmark(group="mean")
def test_mean_of_elements_bench(benchmark, impl):
    rng = np.random.default_rng(seed=1973)
    elements_at_element = rng.integers(-1, 1000, size=(10000, 10))

    values = np.arange(1000, dtype=float)
    actual = np.empty(len(elements_at_element), dtype=float)

    if impl == "numpy":
        func = numpy_mean
    else:
        func = mean_of_children_at_parent

    benchmark(func, np.asarray(elements_at_element), values, actual)

    expected = np.empty_like(actual)
    numpy_mean(np.asarray(elements_at_element), values, expected)

    assert_array_equal(actual, expected)


def test_count_of_elements():
    elements_at_element = [
        [0, 1, 2, 3],
        [4, 5, 6, 7],
        [8, 9, 10, 11],
    ]
    actual = np.empty(len(elements_at_element), dtype=int)
    expected = [4, 4, 4]

    count_of_children_at_parent(np.asarray(elements_at_element), actual)

    assert_array_equal(actual, expected)


def test_count_of_elements_with_missing():
    elements_at_element = [
        [0, 1, 2, -1],
        [-1, -1, 6, 7],
        [8, 9, 10, 11],
        [-1, -1, -1, -1],
    ]
    actual = np.empty(len(elements_at_element), dtype=int)
    expected = [3, 2, 4, 0]

    count_of_children_at_parent(np.asarray(elements_at_element), actual)

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("impl", ["numpy", "cython"])
def test_count_of_elements_benchmark(benchmark, impl):
    elements_at_element = np.random.randint(-1, 2, size=10000).reshape((1000, 10))
    actual = np.empty(len(elements_at_element), dtype=int)

    if impl == "numpy":
        func = numpy_count
    else:
        func = count_of_children_at_parent

    benchmark(func, elements_at_element, actual)

    expected = np.sum(elements_at_element != -1, axis=1)

    assert_array_equal(actual, expected)
