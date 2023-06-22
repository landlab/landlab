import numpy as np
from numpy.testing import assert_array_equal

from landlab.components.network_sediment_transporter.aggregate_parcels import (
    aggregate_parcels_at_link_mean,
    aggregate_parcels_at_link_sum,
)


def test_sum_bench(benchmark):
    n_links = 1000
    n_parcels = 100000
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.zeros(n_parcels, dtype=int)

    benchmark(
        aggregate_parcels_at_link_sum,
        out,
        n_links,
        value_of_parcel,
        link_of_parcel,
        n_parcels,
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
        aggregate_parcels_at_link_mean,
        out,
        n_links,
        value_of_parcel,
        weight_of_parcel,
        link_of_parcel,
        n_parcels,
    )

    assert_array_equal(out[0], 1.0)
    assert_array_equal(out[1:], 0.0)


def test_sum():
    n_links = 10
    n_parcels = 100
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    aggregate_parcels_at_link_sum(
        out, n_links, value_of_parcel, link_of_parcel, n_parcels
    )

    assert_array_equal(out, 10)


def test_mean():
    n_links = 10
    n_parcels = 100
    out = np.empty(n_links, dtype=float)
    value_of_parcel = np.ones(n_parcels, dtype=float)
    weight_of_parcel = np.ones(n_parcels, dtype=float)
    link_of_parcel = np.arange(n_parcels, dtype=int) // n_links

    aggregate_parcels_at_link_mean(
        out, n_links, value_of_parcel, weight_of_parcel, link_of_parcel, n_parcels
    )

    assert_array_equal(out, 1)
