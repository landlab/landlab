import numpy as np

from .ext.remap_element import map_pairs_to_values as _map_pairs_to_values
from .ext.remap_element import (
    map_rolling_pairs_to_values as _map_rolling_pairs_to_values,
)
from .ext.remap_element import pair_isin as _pair_isin


def pair_isin(src, pairs, out=None, sorter=None, sorted=False):
    """Check if integer-pairs are contained in source set.

    Parameters
    ----------
    src : ndarray of int, size *(N, 2)*
        Integer pairs that form the source set.
    pairs : ndarray of int, size *(M, 2)*
        Integer pairs to check if the are contained in the source set.
    out : ndarray of bool, size *(M,)*, optional
        Buffer to place the result. If not provided, a new array will be allocated.
    sorter : ndarray of int, size *(N,)*, optional
        Array of indices that sorts the *src*, as would be returned by *argsort*.
        If not provided, *src* is assumed to already be sorted.
    sorted : bool, optional
        Indicate if the source pairs are already sorted.

    Returns
    -------
    ndarray of bool
        Array that indicates if the pair is contained in the source set.
    """
    if not sorted and sorter is None:
        sorter = np.argsort(src[:, 0])
    if sorter is not None:
        src = src[sorter]

    result = np.empty(len(pairs), dtype=np.uint8)
    _pair_isin(np.ascontiguousarray(src), np.ascontiguousarray(pairs), result)

    if out is None:
        out = result.astype(dtype=bool, copy=False)
    else:
        out[:] = result.astype(dtype=bool, copy=False)
    return out


def map_pairs_to_values(mapping, pairs, out=None, sorter=None, sorted=False):
    """Return the values for integer pairs from a mapping.

    Parameters
    ----------
    mapping : tuple of ndarray of int
        Integer pair to value mapping as *(pairs, values)* where *pairs* is
        *ndarray* of shape *(M, 2)* and *values* an array of length *M*.
    pairs : ndarray of int of shape *(N, 2)*
        Integer pairs to get the values of.
    out : ndarray of bool, size *(N,)*, optional
        Buffer to place the result. If not provided, a new array will be allocated.
    sorter : ndarray of int, size *(M,)*, optional
        Array of indices that sorts the *src*, as would be returned by *argsort*.
        If not provided, *src* is assumed to already be sorted.
    sorted : bool, optional
        Indicate if the mapping key pairs are already sorted.

    Returns
    -------
    ndarray of int
        Array of values of the given integer pairs.

    Examples
    --------
    >>> from landlab.graph.sort.intpair import map_pairs_to_values

    >>> keys = [[0, 1], [1, 1], [2, 1], [3, 1], [4, 1]]
    >>> values = [0, 10, 20, 30, 40]
    >>> pairs = [[1, 1], [3, 1]]
    >>> map_pairs_to_values((keys, values), pairs)
    array([10, 30])
    """
    keys, values = np.asarray(mapping[0]), np.asarray(mapping[1])
    pairs = np.asarray(pairs)

    if out is None:
        out = np.empty(len(pairs), dtype=int)

    if not sorted and sorter is None:
        sorter = np.argsort(keys[:, 0])
    if sorter is not None:
        keys = keys[sorter]
        values = values[sorter]

    _map_pairs_to_values(
        np.ascontiguousarray(keys), np.ascontiguousarray(values), pairs, out
    )

    return out


def map_rolling_pairs_to_values(
    mapping, pairs, out=None, sorter=None, sorted=False, size_of_row=None
):
    """Return the values for integer pairs given as a 2D matrix of rolling
    pairs.

    Parameters
    ----------
    mapping : tuple of ndarray of int
        Integer pair to value mapping as *(pairs, values)* where *pairs* is
        *ndarray* of shape *(N, 2)* and *values* an array of length *N*.
    pairs : ndarray of int of shape *(M, L)*
        Integer pairs to get the values of.
    out : ndarray of bool, size *(M, L)*, optional
        Buffer to place the result. If not provided, a new array will be allocated.
    sorter : ndarray of int, size *(N,)*, optional
        Array of indices that sorts the *src*, as would be returned by *argsort*.
        If not provided, *src* is assumed to already be sorted.
    sorted : bool, optional
        Indicate if the mapping key pairs are already sorted.

    Returns
    -------
    ndarray of int
        Array of values of the given integer pairs.

    Examples
    --------
    >>> from landlab.graph.sort.intpair import map_rolling_pairs_to_values

    >>> keys = [[0, 1], [1, 2], [2, 3], [3, 4], [4, 0]]
    >>> values = [0, 10, 20, 30, 40]
    >>> pairs = [[0, 1, 2, 3], [0, 2, 3, 4]]
    >>> map_rolling_pairs_to_values((keys, values), pairs)
    array([[ 0, 10, 20, -1],
           [-1, 20, 30, 40]])
    """
    keys, values = np.asarray(mapping[0]), np.asarray(mapping[1])
    pairs = np.asarray(pairs)

    if out is None:
        out = np.empty_like(pairs, dtype=int)

    if size_of_row is None:
        size_of_row = np.full(len(pairs), pairs.shape[1], dtype=int)
    else:
        size_of_row = np.asarray(size_of_row)
        out[:] = -1

    if not sorted and sorter is None:
        sorter = np.argsort(keys[:, 0])
    if sorter is not None:
        keys = keys[sorter]
        values = values[sorter]

    _map_rolling_pairs_to_values(
        np.ascontiguousarray(keys),
        np.ascontiguousarray(values),
        np.ascontiguousarray(pairs),
        np.ascontiguousarray(size_of_row),
        out,
    )

    return out
