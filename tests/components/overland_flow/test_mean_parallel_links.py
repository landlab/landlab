import numpy as np
import pytest
from numpy.testing import assert_array_equal

from landlab import RasterModelGrid
from landlab.components.overland_flow._calc import calc_weighted_mean_of_parallel_links


def weighted_mean_reference(vals, parallels, status, where, theta):
    out = np.array(vals)
    for link in where:
        l, r = parallels[link]

        total, n_neighbors = 0.0, 0
        if l != -1 and status[l] != 4:
            n_neighbors += 1
            total += vals[l]
        if r != -1 and status[r] != 4:
            n_neighbors += 1
            total += vals[r]

        if n_neighbors > 0:
            out[link] = theta * vals[link] + (1 - theta) * total / n_neighbors
        else:
            out[link] = vals[link]

    return out


@pytest.mark.parametrize("values_dtype", (np.float32, np.float64))
@pytest.mark.parametrize("where_dtype", (np.int32, np.int64))
@pytest.mark.parametrize("theta", (0.0, 0.25, 0.5, 0.75, 1.0))
def test_weighted_mean(values_dtype, where_dtype, theta):
    grid = RasterModelGrid((3, 4))

    values_at_link = np.asarray(
        [
            *[0, 1, 3],
            *[3, 5, 9, 11],
            *[0, 1, 3],
            *[3, 5, 9, 11],
            *[0, 1, 3],
        ],
        dtype=values_dtype,
    )
    n_links = len(values_at_link)

    actual = np.empty_like(values_at_link)
    rtn = calc_weighted_mean_of_parallel_links(
        values_at_link,
        grid.parallel_links_at_link.astype(where_dtype),
        status_at_link=np.zeros(n_links, dtype=np.uint8),
        where=np.arange(n_links, dtype=where_dtype),
        theta=theta,
        out=actual,
    )
    assert rtn is actual

    expected = weighted_mean_reference(
        values_at_link,
        grid.parallel_links_at_link,
        np.zeros(n_links),
        np.arange(n_links),
        theta,
    )

    assert_array_equal(actual, expected)


@pytest.mark.parametrize("theta", (0.0, 0.25, 0.5, 0.75, 1.0))
def test_weighted_mean_missing_neighbors(theta):
    grid = RasterModelGrid((3, 4))

    values_at_link = np.asarray(
        [
            *[0, 9, 3],
            *[3, 5, 9, 11],
            *[0, 7, 3],
            *[3, 5, 9, 11],
            *[0, 3.14, 3],
        ],
        dtype=float,
    )
    status_at_link = [
        *[4, 0, 4],
        *[0, 0, 0, 0],
        *[4, 0, 4],
        *[0, 0, 0, 0],
        *[4, 0, 4],
    ]
    where = [
        *[0, 1, 0],
        *[0, 0, 0, 0],
        *[0, 1, 0],
        *[0, 0, 0, 0],
        *[0, 1, 0],
    ]
    where = np.flatnonzero(where)

    actual = values_at_link.copy()
    calc_weighted_mean_of_parallel_links(
        values_at_link,
        grid.parallel_links_at_link.astype(int),
        status_at_link=np.asarray(status_at_link, dtype=np.uint8),
        where=where,
        theta=theta,
        out=actual,
    )
    expected = weighted_mean_reference(
        values_at_link,
        grid.parallel_links_at_link,
        status_at_link,
        where,
        theta,
    )

    assert_array_equal(actual[where], [9, 7, 3.14])
    assert_array_equal(actual, expected)
