from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray


def calc_cumulative_path_length(points: ArrayLike) -> NDArray[np.float_]:
    """Calculate the cumulative path length along an n-dimensional curve.

    Parameters
    ----------
    points : array_like
        An array of shape (n_points, n_dims) representing coordinates in n-dimensional
        space.

    Returns
    -------
    ndarray
        A 1D array of cumulative distances of shape (n_points,).

    Examples
    --------
    >>> calc_cumulative_path_length([(0, 0), (0, 1), (2, 1)])
    array([0., 1., 3.])
    >>> calc_cumulative_path_length([(0, 0)])
    array([0.])
    >>> calc_cumulative_path_length([])
    array([], dtype=float64)
    """
    points = np.asarray(points)
    if points.size == 0:
        return np.array([], dtype=float)
    return np.concatenate(([0.0], calc_path_segment_lengths(points).cumsum()))


def calc_path_length(points: ArrayLike) -> float:
    """Calculate the total path length of an n-dimensional curve.

    Parameters
    ----------
    points : array_like
        An array of shape (n_points, n_dims) representing coordinates in n-dimensional
        space.

    Returns
    -------
    float
        Total path length.

    Examples
    --------
    >>> calc_path_length([(0, 0), (0, 1), (2, 1)])
    3.0
    >>> calc_path_length([])
    0.0
    >>> calc_path_length([(0, 1)])
    0.0
    """
    return calc_path_segment_lengths(points).sum()


def calc_path_segment_lengths(points: ArrayLike) -> NDArray[np.float_]:
    """Calculate the distance between each pair of consecutive points in a path.

    Parameters
    ----------
    points : array_like
        An array of shape (n_points, n_dims) representing coordinates in n-dimensional
        space.

    Returns
    -------
    ndarray of shape (n_points - 1,)
        Euclidean distance between each pair of consecutive points.

    Notes
    -----
    If `points` has fewer than two elements, the function returns an empty array.

    Examples
    --------
    >>> calc_path_segment_lengths([(0, 0), (0, 1), (2, 1)])
    array([1., 2.])
    >>> calc_path_segment_lengths([(0, 0, 0, 0), (0, 1, 0, 0)])
    array([1.])
    >>> calc_path_segment_lengths([(0.0, 0.0)])
    array([], dtype=float64)
    """
    points = np.asarray(points)
    if points.size == 0:
        return np.array([], dtype=float)
    if points.ndim != 2:
        raise ValueError(
            f"Expected a 2D array of shape (n_points, n_dims), got shape {points.shape}"
        )
    return np.linalg.norm(np.diff(points, axis=0), axis=1)
