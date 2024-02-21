import numpy as np
from numpy.testing import assert_array_equal

from landlab.data_record.aggregators import (
    aggregate_items_as_count,
    aggregate_items_as_mean,
    aggregate_items_as_sum,
)


def test_count_bench(benchmark):
    n_links = 1000
    n_parcels = 100000
    out = np.empty(n_links, dtype=int)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    benchmark(aggregate_items_as_count, out, n_links, link_of_parcel, n_parcels)

    assert_array_equal(out[0], 100000)
    assert_array_equal(out[1:], 0)


def test_sum_bench(benchmark):
    n_links = 1000
    n_parcels = 100000
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    benchmark(
        aggregate_items_as_sum,
        out,
        n_links,
        link_of_parcel,
        n_parcels,
        value_of_parcel,
    )

    assert_array_equal(out[0], 100000.0)
    assert_array_equal(out[1:], 0.0)


def test_mean_bench(benchmark):
    n_links = 100
    n_parcels = 100000
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    benchmark(
        aggregate_items_as_mean,
        out,
        n_links,
        link_of_parcel,
        n_parcels,
        value_of_parcel,
        weight_of_parcel,
    )

    assert_array_equal(out[0], 1.0)
    assert_array_equal(out[1:], 0.0)


def test_sum():
    n_links = 10
    n_parcels = 100
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    aggregate_items_as_sum(out, n_links, link_of_parcel, n_parcels, value_of_parcel)

    assert_array_equal(out, 10.0)


def test_count():
    n_links = 10
    n_parcels = 100
    out = np.empty(n_links, dtype=int)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    aggregate_items_as_count(out, n_links, link_of_parcel, n_parcels)

    assert_array_equal(out, 10)


def test_sum_with_negative_links():
    n_links = 10
    n_parcels = 100
    out = np.full(n_links, -999, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.full(n_parcels, -1, dtype=int)

    aggregate_items_as_sum(out, n_links, link_of_parcel, n_parcels, value_of_parcel)

    assert_array_equal(out, 0.0)


def test_mean():
    n_links = 10
    n_parcels = 100
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    aggregate_items_as_mean(
        out, n_links, link_of_parcel, n_parcels, value_of_parcel, weight_of_parcel
    )

    assert_array_equal(out, 1.0)


def test_mean_with_negative_links():
    n_links = 10
    n_parcels = 100
    out = np.full(n_links, -999, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.full(n_parcels, -1, dtype=int)

    aggregate_items_as_mean(
        out, n_links, link_of_parcel, n_parcels, value_of_parcel, weight_of_parcel
    )

    assert_array_equal(out, 0.0)
