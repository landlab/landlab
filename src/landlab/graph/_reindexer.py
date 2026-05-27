from __future__ import annotations

from collections.abc import Iterable

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray
from requireit import raise_as
from requireit import require_array
from requireit import require_between
from requireit import require_length
from requireit import require_positive


class Reindexer:
    """Map arrays from an old indexing scheme to a new one.

    Parameters
    ----------
    order : array_like of int
        Permutation such that ``order[new_id] == old_id``. In other words,
        this gives the old index of each item in the new ordering.
    """

    def __init__(self, order: ArrayLike) -> None:
        order = np.asarray(order)
        if order.ndim == 1 and order.size == 0:
            order = order.astype(np.intp)

        with raise_as(ValueError):
            new_to_old = require_array(
                order, shape=("n",), dtype=np.intp, allow_cast=True, name="order"
            )
        new_to_old = np.asarray(new_to_old, dtype=np.intp)

        n = new_to_old.size
        if n and (
            np.any(new_to_old < 0)
            or np.any(new_to_old >= n)
            or np.unique(new_to_old).size != n
        ):
            raise ValueError("order must be a permutation of [0, n)")

        self._new_to_old = new_to_old
        self._old_to_new = np.empty_like(self._new_to_old, dtype=np.intp)
        self._old_to_new[self._new_to_old] = np.arange(
            self._new_to_old.size, dtype=np.intp
        )

    def __len__(self) -> int:
        return self._new_to_old.size

    def __repr__(self) -> str:
        return f"Reindexer(size={len(self)})"

    @classmethod
    def from_xy(cls, xy: NDArray, axes: Iterable[int] | None = None) -> Reindexer:
        """Create a reindexer from lexicographic sorting of coordinates."""
        return cls(sort_order_from_coords(xy, axes=axes))

    @classmethod
    def from_ids(cls, ids: ArrayLike) -> Reindexer:
        """Create a reindexer that orders elements by integer ids.

        Parameters
        ----------
        ids : array_like of int
            Array of integer ids that define the desired ordering. The returned
            reindexer sorts elements in ascending order of `ids`, such that
            elements with smaller ids appear first.

        Returns
        -------
        Reindexer
            A reindexer whose ordering is given by ``np.argsort(ids)``.

        Notes
        -----
        This constructor is useful for deriving an ordering from a related
        element type. For example, if each cell corresponds to a node via
        `node_at_cell`, then ``Reindexer.from_ids(node_at_cell)`` produces
        a cell ordering consistent with the current node ordering.

        Ties in `ids` are resolved by the behavior of ``np.argsort``.
        """
        with raise_as(ValueError):
            ids = require_array(
                np.asarray(ids), shape=("n_ids",), dtype=np.integer, name="ids"
            )
        ids = np.asarray(ids, dtype=np.intp)

        return cls(np.argsort(ids, kind="stable"))

    def reorder(self, array: ArrayLike, out: NDArray | None = None) -> NDArray:
        """Reorder an array along axis 0.

        Parameters
        ----------
        array : array_like
            Array whose first axis uses the old indexing.
        out : ndarray, optional
            Output array. Must have the same shape as `array`.

        Returns
        -------
        ndarray
            Array reordered so that axis 0 uses the new indexing.
        """
        with raise_as(ValueError):
            array = require_length(np.asarray(array), length=len(self), name="array")

        if out is None:
            out = np.empty_like(array)
        elif np.shares_memory(array, out):
            raise ValueError("out must not share memory with array")

        if out.shape != array.shape:
            raise ValueError("out must have the same shape as array")
        with raise_as(ValueError):
            require_array(out, dtype=array.dtype, name="out")

        out[:] = array[self._new_to_old]
        return out

    def remap(
        self,
        array: ArrayLike,
        *,
        fill_value: int | None = None,
        out: NDArray | None = None,
    ) -> NDArray:
        """Remap integer ids from the old indexing to the new indexing.

        Parameters
        ----------
        array : array_like of int
            Integer id array expressed in the old indexing.
        fill_value : int, optional
            Sentinel value that should be left unchanged.
        out : ndarray, optional
            Output array. Must have the same shape as `array`.

        Returns
        -------
        ndarray
            Array with ids remapped into the new indexing.
        """
        with raise_as(ValueError):
            array = require_array(
                np.asarray(array, dtype=np.intp), dtype=np.integer, name="array"
            )

        if out is None:
            out = np.empty_like(array)
        elif np.shares_memory(array, out):
            raise ValueError("out must not share memory with array")

        if out.shape != array.shape:
            raise ValueError("out must have the same shape as array")
        with raise_as(ValueError):
            require_array(out, dtype=np.integer, name="out")

        if fill_value is None:
            with raise_as(ValueError):
                require_between(array, a_min=0, a_max=len(self) - 1)
            out[...] = self._old_to_new[array]
        else:
            mask = array != fill_value
            valid_ids = array[mask]
            with raise_as(ValueError):
                require_between(valid_ids, a_min=0, a_max=len(self) - 1)
            out[...] = array
            out[mask] = self._old_to_new[valid_ids]

        return out


def sort_order_from_coords(
    coords: ArrayLike,
    axes: Iterable[int] | None = None,
    tol: float | None = None,
) -> NDArray:
    """Return indices that lexicographically sort coordinates.

    Parameters
    ----------
    coords : array_like of ``shape (n, ndim)``
        Coordinates to sort.
    axes : tuple of int, optional
        Permutation of coordinate axes that defines the sort priority.
        The axes are interpreted from **least** significant to most
        **significant** key (i.e., the last axis in the tuple is the
        primary sort key). By default, axes are ``range(ndim)``, so the
        last column is the primary key.
    tol : float, optional
        Quantization bin width. Coordinates are rounded to the nearest
        multiple of `tol` before sorting, so two values are treated as
        equal only if they are within ``tol / 2`` of the same multiple.
        If ``None``, coordinates are compared exactly.

    Returns
    -------
    ndarray of int
        Indices that reorder `coords` into lexicographic order.

    Notes
    -----
    Sorting is performed using ``np.lexsort``. Because of its semantics,
    the last key has highest priority. For example, with ``axes=(1, 0)``,
    sorting is primarily by column 0, and ties are broken by column 1.

    Examples
    --------
    ::

        1 --- 3
        |     |
        0 --- 2

    >>> xy = np.asarray([[0, 0], [0, 1], [1, 0], [1, 1]], dtype=float)

    The default is to first sort by column 1, and then column 0.

    >>> order = sort_order_from_coords(xy)
    >>> order
    array([0, 2, 1, 3])
    >>> xy[order]
    array([[0., 0.],
           [1., 0.],
           [0., 1.],
           [1., 1.]])

    When specifying ``axes``, the last value gives the dimension with
    the highest priority. For example, ``axes=(1, 0)`` sorts primarily
    by column 0, and ties are broken by column 1.

    >>> order = sort_order_from_coords(xy, axes=(1, 0))
    >>> order
    array([0, 1, 2, 3])
    >>> xy[order]
    array([[0., 0.],
           [0., 1.],
           [1., 0.],
           [1., 1.]])
    """
    xy = np.asarray(coords)
    if tol is not None:
        with raise_as(ValueError):
            tol = require_positive(tol, name="tol")

    if xy.ndim != 2:
        raise ValueError("coords must be 2D")

    n_dim = xy.shape[1]
    axes = range(n_dim) if axes is None else axes

    if set(axes) != set(range(n_dim)):
        raise ValueError(f"axes must be a permutation of {tuple(range(n_dim))}")

    if tol is not None:
        xy = np.round(xy / tol).astype(np.int64)

    return np.lexsort(tuple(xy[:, axis] for axis in axes))
