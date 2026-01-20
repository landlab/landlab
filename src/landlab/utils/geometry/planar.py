from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray
from scipy.spatial import cKDTree


class NearestNodeFinder:
    def __init__(self, coords_of_node: ArrayLike, leafsize: int = 64):
        coords_of_node = np.asarray(coords_of_node, dtype=float)
        if coords_of_node.ndim != 2:
            raise ValueError(
                "`coords_of_node` must have shape (n_nodes, n_dims),"
                f" got {coords_of_node.shape}."
            )
        self._tree = cKDTree(coords_of_node, leafsize=leafsize)

    def __call__(self, coords_of_target) -> np.ndarray:
        """Find indices of the nearest point in the source for each target.

        Parameters
        ----------
        coords_of_target : (n_targets, n_dims) or (n_dims,) array-like
            Coordinates of the target points.

        Returns
        -------
        nearest_nodes : int or (n_targets,) ndarray of int
            Indices into `coords_of_node`. If a single target is provided,
            return a scalar `int`.

        Examples
        --------
        >>> find_nearest = NearestNodeFinder(
        ...     [[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 2.0]]
        ... )
        >>> find_nearest([-0.25, 0.75])
        0
        >>> find_nearest([[-0.25, 0.75], [0.75, 1.75]])
        array([0, 2])
        """
        coords_of_target = np.atleast_2d(np.asarray(coords_of_target, dtype=float))
        _, nearest_nodes = self._tree.query(coords_of_target, k=1)
        return nearest_nodes if coords_of_target.shape[0] > 1 else nearest_nodes.item()


def find_nearest_node(coords_of_node: ArrayLike, coords_of_target: ArrayLike):
    """Find the index of the nearest node(s) to given target point(s).

    Parameters
    ----------
    coords_of_node : array_like, shape (n_nodes, n_dims)
        Coordinates of the source nodes.
    coords_of_target : array_like, shape (n_targets, n_dims) or (n_dims,)
        Coordinates of the target point(s).

    Returns
    -------
    nearest_nodes : int or (n_targets,) ndarray of int
        Indices into `coords_of_node`. If a single target is provided,
        return a scalar `int`.

    Examples
    --------
    >>> find_nearest_node(
    ...     [[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 2.0]], [0.25, 0.75]
    ... )
    0
    >>> find_nearest_node(
    ...     [[0.0, 0.0], [1.0, 0.0], [1.0, 2.0], [0.0, 2.0]],
    ...     [[0.25, 0.75], [0.75, 1.75]],
    ... )
    array([0, 2])
    """
    return NearestNodeFinder(coords_of_node)(coords_of_target)


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
