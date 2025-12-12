import numpy as np
import pytest
from numpy.testing import assert_array_equal

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
