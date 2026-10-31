from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray
from requireit import require_array
from requireit import require_between
from requireit import require_greater_than
from requireit import require_greater_than_or_equal
from requireit import require_length_at_least

from landlab.data_record._aggregators import AGGREGATE_ITEMS_MEMORY_ERROR
from landlab.data_record._aggregators import AGGREGATE_ITEMS_ZERO_TOTAL_WEIGHT
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


@dataclass(frozen=True)
class PreparedAggregatorInputs:
    ids: NDArray[np.integer]
    active: NDArray[np.bool_]
    out: NDArray


def aggregate_items_as_sum(
    ids: ArrayLike,
    values: ArrayLike,
    *,
    where: ArrayLike | None = None,
    out: NDArray[np.floating] | None = None,
) -> NDArray[np.floating]:
    """Find the sum of values associated with an id.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    values : array_like
        The value associated with the corresponding id in the `id` array.
    where : array_like of bool, shape (n,), optional
        Boolean mask indicating which items to include. If not provided,
        all items are included.
    out : ndarray of float, shape (n_groups,), optional
        Output array. If provided, results are written in-place.

    Returns
    -------
    ndarray of int
        The sum of the values at each id.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_sum
    >>> aggregate_items_as_sum([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([3., 3., 0., 3., 1., 5.])

    >>> out = np.full(8, -1, dtype=float)
    >>> aggregate_items_as_sum([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5], out=out)
    array([3., 3., 0., 3., 1., 5., 0., 0.])
    """
    inputs = _prepare_aggregator_inputs(ids, where=where, out=out)
    n_items = len(inputs.ids)

    values = np.ascontiguousarray(values, dtype=float)
    require_array(values, shape=(n_items,), contiguous=True, name="values")

    _aggregate_items_as_sum(inputs.out, inputs.ids, values, inputs.active)

    return inputs.out


def aggregate_items_as_mean(
    ids: ArrayLike,
    values: ArrayLike,
    *,
    weights: ArrayLike | None = None,
    where: ArrayLike | None = None,
    out: NDArray[np.floating] | None = None,
) -> NDArray[np.floating]:
    """Find the mean of values associated with an id.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    values : array_like
        The value associated with the corresponding id in the `id` array.
    where : array_like of bool, shape (n,), optional
        Boolean mask indicating which items to include. If not provided,
        all items are included.
    out : ndarray of float, shape (n_groups,), optional
        Output array. If provided, results are written in-place.

    Returns
    -------
    ndarray of int
        The mean of the values at each id.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_mean
    >>> aggregate_items_as_mean([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5])
    array([1.5, 3. , 0. , 3. , 1. , 5. ])

    >>> out = np.zeros(8, dtype=float)
    >>> aggregate_items_as_mean([0, 0, 1, 3, 4, 5], [1, 2, 3, 3, 1, 5], out=out)
    array([1.5, 3. , 0. , 3. , 1. , 5. , 0. , 0. ])
    """
    inputs = _prepare_aggregator_inputs(ids, where=where, out=out)
    n_items = len(inputs.ids)

    if not np.any(inputs.active):
        return inputs.out

    values = np.ascontiguousarray(values)

    if weights is None:
        weights = np.ones(n_items, dtype=float)
    weights = np.ascontiguousarray(weights)

    values = require_array(values, shape=(n_items,), contiguous=True, name="values")
    weights = require_array(weights, shape=(n_items,), contiguous=True, name="weights")

    require_greater_than_or_equal(weights[inputs.active], 0.0, name="weights")

    status = _aggregate_items_as_mean(
        inputs.out, inputs.ids, values, weights, inputs.active
    )

    if status == AGGREGATE_ITEMS_MEMORY_ERROR:  # pragma: no cover
        raise MemoryError("failed to allocate temporary workspace")
    if status == AGGREGATE_ITEMS_ZERO_TOTAL_WEIGHT:
        raise ValueError("weights must sum to a positive value for selected groups")

    return inputs.out


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
        Boolean mask indicating which items to include. If provided, ``out``
        must also be provided so that groups with no selected items have a
        defined output value.
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
    >>> aggregate_items_as_gmean(ids, values, where=where, out=np.empty(2, dtype=float))
    array([2., 9.])

    Reuse an output array:

    >>> out = np.full(2, np.nan)
    >>> where = [True, True, False, False, False]
    >>> aggregate_items_as_gmean(ids, values, out=out, where=where)
    array([ 2., nan])
    >>> out
    array([ 2., nan])
    """
    inputs = _prepare_aggregator_inputs(ids, where=where, out=out)
    n_items = len(inputs.ids)

    values = np.ascontiguousarray(values)

    if not np.any(inputs.active):
        return inputs.out

    if weights is None:
        weights = np.ones(n_items, dtype=float)
    weights = np.ascontiguousarray(weights)

    values = require_array(values, shape=(n_items,), contiguous=True, name="values")
    weights = require_array(weights, shape=(n_items,), contiguous=True, name="weights")

    require_greater_than(values[inputs.active], 0.0, name="values")
    require_greater_than_or_equal(weights[inputs.active], 0.0, name="weights")

    status = _aggregate_items_as_gmean(
        inputs.out, inputs.ids, values, weights, inputs.active
    )

    if status == AGGREGATE_ITEMS_MEMORY_ERROR:  # pragma: no cover
        raise MemoryError("failed to allocate temporary workspace")
    if status == AGGREGATE_ITEMS_ZERO_TOTAL_WEIGHT:
        raise ValueError("weights must sum to a positive value for selected groups")

    return inputs.out


def aggregate_items_as_count(
    ids: ArrayLike,
    *,
    where: ArrayLike | None = None,
    out: NDArray[np.floating] | None = None,
) -> NDArray[np.int_]:
    """Count the number of time an id appears in an array.

    Parameters
    ----------
    ids : array_like of int
        An array of ids.
    where : array_like of bool, shape (n,), optional
        Boolean mask indicating which items to include. If not provided,
        all items are included.
    out : ndarray of float, shape (n_groups,), optional
        Output array. If provided, results are written in-place.

    Returns
    -------
    ndarray of int
        The number of times each id appears.

    Examples
    --------
    >>> from landlab.data_record.aggregators import aggregate_items_as_count
    >>> aggregate_items_as_count([1, 2, 3, 3, 1, 5])
    array([0, 2, 1, 2, 0, 1])

    >>> aggregate_items_as_count([1, 2, 3, 3, 1, 5], out=np.empty(8, dtype=float))
    array([0, 2, 1, 2, 0, 1, 0, 0])
    """
    inputs = _prepare_aggregator_inputs(ids, where=where, out=out, dtype=int)

    _aggregate_items_as_count(inputs.out, inputs.ids, inputs.active)

    return inputs.out


def _prepare_aggregator_inputs(
    ids: ArrayLike,
    *,
    dtype=float,
    where: ArrayLike | None = None,
    out: NDArray[np.floating | np.integer] | None = None,
) -> PreparedAggregatorInputs:
    if where is not None and out is None:
        raise ValueError("'out' must be provided when 'where' is used")

    ids = np.ascontiguousarray(ids)

    require_array(ids, dtype=np.integer, contiguous=True, shape=("n",), name="ids")

    n_items = ids.shape[0]

    active = ids >= 0
    if where is not None:
        where = np.ascontiguousarray(where)
        require_array(where, shape=(n_items,), name="where")
        active &= where

    if np.any(active):
        n_groups = int(np.max(ids[active])) + 1
    else:
        n_groups = 0

    if out is None:
        out = np.empty(n_groups, dtype=dtype)
    out = require_array(np.asarray(out), shape=("n_groups",), name="out")
    require_length_at_least(out, n_groups, name="out")

    require_between(ids[active], 0, len(out), inclusive_max=False, name="ids")
    active = require_array(active, contiguous=True, dtype=np.bool_, name="where")
    out = require_array(out, writable=True, contiguous=True, dtype=dtype, name="out")

    return PreparedAggregatorInputs(ids=ids, active=active, out=out)
