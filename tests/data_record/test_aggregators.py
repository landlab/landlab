import numpy as np
import pytest
import requireit
from numpy.testing import assert_allclose
from numpy.testing import assert_array_equal
from scipy.stats import gmean

from landlab.data_record._aggregators import (
    aggregate_items_as_count as _aggregate_items_as_count,
)
from landlab.data_record._aggregators import (
    aggregate_items_as_mean as _aggregate_items_as_mean,
)
from landlab.data_record._aggregators import (
    aggregate_items_as_sum as _aggregate_items_as_sum,
)
from landlab.data_record.aggregators import aggregate_items_as_count
from landlab.data_record.aggregators import aggregate_items_as_gmean
from landlab.data_record.aggregators import aggregate_items_as_mean
from landlab.data_record.aggregators import aggregate_items_as_sum


def test_count_bench_cython():
    n_links = 1000
    n_parcels = 100000
    out = np.empty(n_links, dtype=int)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    _aggregate_items_as_count(out, link_of_parcel)

    assert_array_equal(out[0], 100000)
    assert_array_equal(out[1:], 0)


@pytest.mark.slow
def test_count_bench():
    n_links = 1000
    n_parcels = 100000
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    out = aggregate_items_as_count(link_of_parcel, size=n_links)

    assert_array_equal(out[0], 100000)
    assert_array_equal(out[1:], 0)


def test_sum_bench_cython():
    n_links = 1000
    n_parcels = 100000
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    _aggregate_items_as_sum(out, link_of_parcel, value_of_parcel)

    assert_array_equal(out[0], 100000.0)
    assert_array_equal(out[1:], 0.0)


def test_sum_bench():
    n_links = 1000
    n_parcels = 100000
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    out = aggregate_items_as_sum(link_of_parcel, value_of_parcel, size=n_links)

    assert_array_equal(out[0], 100000.0)
    assert_array_equal(out[1:], 0.0)


def test_mean_bench_cython():
    n_links = 100
    n_parcels = 100000
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    _aggregate_items_as_mean(out, link_of_parcel, value_of_parcel, weight_of_parcel)

    assert_array_equal(out[0], 1.0)
    assert_array_equal(out[1:], 0.0)


def test_mean_bench():
    n_links = 100
    n_parcels = 100000
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    out = aggregate_items_as_mean(
        link_of_parcel,
        value_of_parcel,
        weights=weight_of_parcel,
        size=n_links,
    )

    assert_array_equal(out[0], 1.0)
    assert_array_equal(out[1:], 0.0)


def test_sum():
    n_links = 10
    n_parcels = 100
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    out = aggregate_items_as_sum(link_of_parcel, value_of_parcel, size=n_links)
    assert_array_equal(out, 10.0)


def test_count():
    n_links = 10
    n_parcels = 100
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    out = aggregate_items_as_count(link_of_parcel)
    assert_array_equal(out, 10)


def test_sum_with_negative_links():
    n_parcels = 100
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.full(n_parcels, -1, dtype=int)

    out = aggregate_items_as_sum(link_of_parcel, value_of_parcel)

    assert_array_equal(out, 0.0)


def test_mean():
    n_links = 10
    n_parcels = 100
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    out = aggregate_items_as_mean(
        link_of_parcel, value_of_parcel, weights=weight_of_parcel
    )

    assert_array_equal(out, 1.0)


def test_mean_with_negative_links():
    n_parcels = 100
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.full(n_parcels, -1, dtype=int)

    out = aggregate_items_as_mean(
        link_of_parcel, value_of_parcel, weights=weight_of_parcel
    )

    assert_array_equal(out, 0.0)


def test_aggregate_items_as_gmean():
    ids = np.asarray([0, 0, 1, 1])
    values = np.asarray([1.0, 4.0, 3.0, 12.0])

    actual = aggregate_items_as_gmean(ids, values)
    expected = [gmean(values[ids == 0]), gmean(values[ids == 1])]

    assert_allclose(actual, expected)


def test_aggregate_items_as_gmean_with_weights():
    ids = np.array([0, 0, 1, 1, 1])
    values = np.array([1.0, 4.0, 1.0, 3.0, 9.0])
    weights = np.array([1.0, 1.0, 1.0, 2.0, 1.0])

    actual = aggregate_items_as_gmean(ids, values, weights=weights)
    expected = [
        gmean(values[ids == 0], weights=weights[ids == 0]),
        gmean(values[ids == 1], weights=weights[ids == 1]),
    ]

    assert_allclose(actual, expected)


def test_aggregate_items_as_gmean_with_where():
    ids = np.array([0, 0, 1, 1, 1])
    values = np.array([1.0, 4.0, 1.0, 3.0, 9.0])
    where = np.array([True, True, False, False, True])

    actual = aggregate_items_as_gmean(ids, values, where=where)
    expected = [gmean(values[(ids == 0) & where]), gmean(values[(ids == 1) & where])]

    assert_allclose(actual, expected)


def test_aggregate_items_as_gmean_ignores_negative_ids():
    ids = [0, -1, 0, 1]
    values = [2.0, 100.0, 8.0, 9.0]

    actual = aggregate_items_as_gmean(ids, values)

    assert_allclose(actual, [4.0, 9.0])


def test_aggregate_items_as_gmean_leaves_unselected_out_values_unchanged():
    ids = np.array([0, 2])
    values = np.array([4.0, 9.0])
    out = np.array([-1.0, -2.0, -3.0, -4.0])

    actual = aggregate_items_as_gmean(ids, values, out=out)

    assert actual is out
    assert_allclose(actual, [4.0, -2.0, 9.0, -4.0])


def test_aggregate_items_as_gmean_nowhere_is_noop():
    ids = [0, 1, 2]
    values = [4.0, 9.0, 16.0]
    where = [False, False, False]
    out = np.array([-1.0, -2.0, -3.0])

    actual = aggregate_items_as_gmean(ids, values, where=where, out=out)

    assert actual is out
    assert_array_equal(actual, [-1.0, -2.0, -3.0])


def test_aggregate_items_as_gmean_only_validates_selected_values():
    ids = [0, 1]
    values = [4.0, -1.0]
    where = [True, False]

    actual = aggregate_items_as_gmean(ids, values, where=where, out=np.full(2, np.nan))

    assert_allclose(actual, [4.0, np.nan], equal_nan=True)


@pytest.mark.parametrize("bad_values", [[0.0], [-1.0]])
def test_aggregate_items_as_gmean_rejects_nonpositive_selected_values(bad_values):
    with pytest.raises(requireit.ValidationError):
        aggregate_items_as_gmean([0], bad_values)


def test_aggregate_items_as_gmean_rejects_negative_selected_weights():
    with pytest.raises(requireit.ValidationError):
        aggregate_items_as_gmean([0], [1.0], weights=[-1.0])


def test_aggregate_items_as_gmean_allows_zero_weight():
    actual = aggregate_items_as_gmean([0, 0], [2.0, 8.0], weights=[0.0, 1.0])

    assert_allclose(actual, [8.0])


def test_aggregate_items_as_gmean_rejects_out_that_is_too_small():
    out = np.empty(2, dtype=float)
    with pytest.raises(requireit.ValidationError):
        aggregate_items_as_gmean([0, 2], [1.0, 4.0], out=out)


@pytest.mark.parametrize(
    "kwds",
    [
        pytest.param({"values": [1.0, 2.0]}, id="values"),
        pytest.param({"weights": [1.0, 1.0]}, id="weights"),
        pytest.param({"where": [True, False]}, id="where"),
    ],
)
def test_aggregate_items_as_gmean_rejects_shape_mismatch(kwds):
    args = {"ids": [0], "values": [1.0]}
    args.update(kwds)

    with pytest.raises(requireit.ValidationError):
        aggregate_items_as_gmean(**args)
