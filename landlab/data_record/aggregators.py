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
    values = np.asarray(values, dtype=float)
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=float)

    # assert len(values) >= len(out)

    _aggregate_items_as_sum(out, len(out), ids, len(ids), values)

    return out


def aggregate_items_as_mean(
    ids: ArrayLike,
    values: ArrayLike,
    weights: ArrayLike | None = None,
    size: int | None = None,
) -> NDArray[np.floating]:
    values = np.asarray(values)
    if weights is None:
        weights = np.ones_like(values)
    else:
        weights = np.asarray(weights, dtype=values.dtype)
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=float)

    assert len(values) == len(weights)
    # assert len(values) >= len(out)

    _aggregate_items_as_mean(out, len(out), ids, len(ids), values, weights)

    return out


def aggregate_items_as_count(
    ids: ArrayLike, size: int | None = None
) -> NDArray[np.int_]:
    ids = np.asarray(ids, dtype=int)

    size = _validate_size(ids, size=size)

    out = np.empty(size, dtype=int)

    _aggregate_items_as_count(out, len(out), ids, len(ids))

    return out


def _validate_size(ids: NDArray[np.int_], size: int | None = None):
    if size is None:
        size = ids.max() + 1
    else:
        assert (
            size >= ids.max() + 1
        ), "size must be greater than or equal to the largest input id"
    return size
