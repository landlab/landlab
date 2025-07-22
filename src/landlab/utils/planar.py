from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray


def calc_cumulative_path_length(points: ArrayLike) -> NDArray[np.float_]:
    """Calculate the cumulative path length along a 2D curve.

    Parameters
    ----------
    points : array_like
        An array of shape (n_points, 2) representing x, y coordinates.

    Returns
    -------
    np.ndarray
        A 1D array of cumulative distances of shape (n_points,).

    Examples
    --------
    >>> calc_cumulative_path_length([(0, 0), (0, 1), (2, 1)])
    array([0., 1., 3.])
    """
    return np.linalg.norm(np.diff(points, prepend=[points[0]], axis=0), axis=1).cumsum()


def calc_path_length(points: ArrayLike) -> float:
    """Calculate the total path length of a 2D curve.

    Parameters
    ----------
    points : array_like
        An array of shape (n_points, 2) representing x, y coordinates.

    Returns
    -------
    float
        Total path length.

    Examples
    --------
    >>> calc_path_length([(0, 0), (0, 1), (2, 1)])
    3.0
    """
    return np.linalg.norm(np.diff(points, axis=0), axis=1).sum()
