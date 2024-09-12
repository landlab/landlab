from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray

from landlab.data_record._aggregators import (
    aggregate_items_as_count as _aggregate_items_as_count,
)
from landlab.data_record._aggregators import (
    aggregate_items_as_mean as _aggregate_items_as_mean,
)
from landlab.data_record._aggregators import (
    aggregate_items_as_sum as _aggregate_items_as_sum,
)


def aggregate_items_as_sum(
    ids: ArrayLike, values: ArrayLike, size: int | None = None
) -> NDArray[np.floating]:
    """Find the sum of values associated with an id.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    values : array_like
        The value associated with the corresponding id in the `id` array.
    size : int, optional
        The size of the output array. This is useful if the `ids`
        array doesn't contain all possible ids.

    Returns
    -------
    ndarray of int
        The sum of the values at each id.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_sum
    >>> aggregate_items_as_sum([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([3., 3., 0., 3., 1., 5.])
    >>> aggregate_items_as_sum([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5], size=8)
    array([3., 3., 0., 3., 1., 5., 0., 0.])

    Negative ids are ignored.

    >>> aggregate_items_as_sum([0, -1, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([1., 3., 0., 3., 1., 5.])
    """
    values = np.asarray(values, dtype=float)
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=float)

    _aggregate_items_as_sum(out, ids, values)

    return out


def aggregate_items_as_mean(
    ids: ArrayLike,
    values: ArrayLike,
    weights: ArrayLike | None = None,
    size: int | None = None,
) -> NDArray[np.floating]:
    """Find the mean of values associated with an id.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    values : array_like
        The value associated with the corresponding id in the `id` array.
    size : int, optional
        The size of the output array. This is useful if the `ids`
        array doesn't contain all possible ids.

    Returns
    -------
    ndarray of int
        The mean of the values at each id.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_mean
    >>> aggregate_items_as_mean([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([1.5, 3. , 0. , 3. , 1. , 5. ])
    >>> aggregate_items_as_mean([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5], size=8)
    array([1.5, 3. , 0. , 3. , 1. , 5. , 0. , 0. ])

    Negative ids are ignored.

    >>> aggregate_items_as_mean([0, -1, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([1., 3., 0., 3., 1., 5.])
    """
    values = np.asarray(values)
    if weights is None:
        weights = np.ones_like(values)
    else:
        weights = np.asarray(weights, dtype=values.dtype)
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=float)

    assert len(values) == len(weights)

    _aggregate_items_as_mean(out, ids, values, weights)

    return out


def aggregate_items_as_count(
    ids: ArrayLike, size: int | None = None
) -> NDArray[np.int_]:
    """Count the number of time an id appears in an array.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    size : int, optional
        The size of the output array. This is useful if the `ids`
        array doesn't contain all possible ids.

    Returns
    -------
    ndarray of int
        The number of times each id appears.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_count
    >>> aggregate_items_as_count([1, 2, 3, 3, 1, 5])
    array([0, 2, 1, 2, 0, 1])
    >>> aggregate_items_as_count([1, 2, 3, 3, 1, 5], size=8)
    array([0, 2, 1, 2, 0, 1, 0, 0])

    Negative ids are ignored.

    >>> aggregate_items_as_count([1, 2, 3, 3, -1, 5])
    array([0, 1, 1, 2, 0, 1])
    """
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=int)

    _aggregate_items_as_count(out, ids)

    return out


def _validate_size(ids: NDArray[np.int_], size: int | None = None):
    if size is None:
        size = ids.max() + 1
    else:
        assert (
            size >= ids.max() + 1
        ), "size must be greater than or equal to the largest input id"
    return size
