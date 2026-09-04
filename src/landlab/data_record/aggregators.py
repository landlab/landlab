from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray
from requireit import require_array
from requireit import require_greater_than
from requireit import require_greater_than_or_equal
from requireit import require_length_at_least

from landlab.data_record._aggregators import (
    aggregate_items_as_count as _aggregate_items_as_count,
)
from landlab.data_record._aggregators import (
    aggregate_items_as_gmean as _aggregate_items_as_gmean,
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


def aggregate_items_as_gmean(
    ids: ArrayLike,
    values: ArrayLike,
    *,
    weights: ArrayLike | None = None,
    where: ArrayLike | None = None,
    out: NDArray[np.floating] | None = None,
) -> NDArray[np.floating]:
    """Compute a weighted geometric mean of values grouped by integer IDs.

    Group ``values`` by ``ids`` and compute the weighted geometric mean
    for each group.

    Parameters
    ----------
    ids : array_like of int, shape (n,)
        Integer group labels for each value. Negative IDs are ignored.
    values : array_like, shape (n,)
        Values to aggregate. Must be strictly positive.
    weights : array_like, shape (n,), optional
        Weights associated with each value. If not provided, all weights
        are taken to be 1.
    where : array_like of bool, shape (n,), optional
        Boolean mask indicating which items to include. If not provided,
        all items are included.
    out : ndarray of float, shape (n_groups,), optional
        Output array. If provided, results are written in-place.

    Returns
    -------
    out : ndarray of float, shape (n_groups,)
        Weighted geometric mean for each group. Elements with no valid
        items are left unchanged in ``out`` if provided, or contain
        uninitialized values if ``out`` was created internally.

    Examples
    --------
    >>> import numpy as np

    >>> ids = [0, 0, 1, 1, 1]
    >>> values = [1.0, 4.0, 1.0, 3.0, 9.0]
    >>> weights = [1.0, 1.0, 1.0, 2.0, 1.0]

    >>> aggregate_items_as_gmean(ids, values, weights=weights)
    array([2., 3.])

    Use ``where`` to filter items:

    >>> where = [True, True, False, False, True]
    >>> aggregate_items_as_gmean(ids, values, where=where)
    array([2., 9.])

    Reuse an output array:

    >>> out = np.full(2, np.nan)
    >>> where = [True, True, False, False, False]
    >>> aggregate_items_as_gmean(ids, values, out=out, where=where)
    array([ 2., nan])
    >>> out
    array([ 2., nan])
    """
    ids, values = np.asarray(ids), np.asarray(values)

    require_array(ids, dtype=np.integer, contiguous=True, shape=("n",), name="ids")

    n_items, max_id = ids.shape[0], np.max(ids)

    max_id = np.max(ids)
    if out is None:
        out = np.empty(max_id + 1, dtype=float)
    out = require_array(np.asarray(out), shape=("n_items",), name="out")
    require_length_at_least(out, max_id + 1, name="out")

    if where is None:
        where = np.full(n_items, True, dtype=np.bool_)
    where = require_array(np.asarray(where), shape=(n_items,), name="where")

    if not np.any(where):
        return out

    if weights is None:
        weights = np.ones(n_items, dtype=float)
    weights = np.asarray(weights)

    out = require_array(out, writable=True, contiguous=True, dtype=float, name="out")
    values = require_array(values, shape=(n_items,), contiguous=True, name="values")
    weights = require_array(weights, shape=(n_items,), contiguous=True, name="weights")
    where = require_array(where, contiguous=True, dtype=np.bool_, name="where")

    require_greater_than(values[where], 0.0, name="values")
    require_greater_than_or_equal(weights[where], 0.0, name="weights")

    _aggregate_items_as_gmean(out, ids, values, weights, where)

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
