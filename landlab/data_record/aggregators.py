import numpy as np
from numpy.typing import ArrayLike, NDArray

from landlab.data_record.aggregators import (
    aggregate_items_as_count as _aggregate_items_as_count,
    aggregate_items_as_mean as _aggregate_items_as_mean,
    aggregate_items_as_sum as _aggregate_items_as_sum,
)


def aggregate_items_as_sum(values: ArrayLike, ids: ArrayLike) -> NDArray[np.floating]:
    values = np.asarray(values, dtype=float)
    ids = np.asarray(ids, dtype=int)

    out = np.zeros(ids.max() + 1, dtype=float)

    assert len(values) >= len(out)

    _aggregate_items_as_sum(out, len(out), ids, len(ids), values)

    return out


def aggregate_items_as_mean(
    values: ArrayLike, ids: ArrayLike, weights: ArrayLike | None = None
) -> NDArray[np.floating]:
    values = np.asarray(values)
    if weights is None:
        weights = np.ones_like(values)
    else:
        weights = np.asarray(weights, dtype=values.dtype)
    ids = np.asarray(ids, dtype=int)

    out = np.zeros(ids.max() + 1, dtype=float)

    assert len(values) == len(weights)
    assert len(values) >= len(out)

    _aggregate_items_as_mean(out, len(out), ids, len(ids), values, weights)

    return out


def aggregate_items_as_count(ids: ArrayLike) -> NDArray[np.int_]:
    ids = np.asarray(ids, dtype=int)
    out = np.zeros(ids.max() + 1, dtype=int)

    _aggregate_items_as_count(out, len(out), ids, len(ids))

    return out
