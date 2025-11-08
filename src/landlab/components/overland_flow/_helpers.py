from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike
from numpy.typing import NDArray


def update_water_depths(
    h_at_node: ArrayLike,
    q_at_node: ArrayLike,
    *,
    dt: float = 1.0,
    rainfall_rate_at_node: ArrayLike | None = None,
    where: ArrayLike | bool = True,
    out: NDArray | None = None,
) -> NDArray:
    """Update water depths

    Calculate the change in water depths on nodes by finding the
    difference between inputs (rainfall) and the inputs/outputs (flux
    divergence of discharge)

    Parameters
    ----------
    h_at_node : array_like
        Current water depths.
    q_at_node : array_like
        Water discharge at each node. A positive value means
        water is leaving the node.
    dt : float, optional
        Timestep.
    rainfall_rate_at_node : array_like, optional
        Rainfall rate at each node. If None, treated as 0.0.
    where : array_like of bool, optional
        Boolean mask selecting which nodes to update. If ``None``,
        all nodes are updated (default NumPy ufunc behavior).
    out : ndarray, optional
        Destination array. If ``None``, create a new array.

    Returns
    -------
    out : ndarray
        Updated water depths.

    Examples
    --------
    >>> update_water_depths([1.0, 1.0, 11.0, 1.0, 1.0], [0.5, -3.0, 5.0, 2.0, 1.0])
    array([0.5, 4. , 6. , 0. , 0. ])
    >>> update_water_depths(
    ...     [1.0, 1.0, 11.0, 1.0, 1.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     rainfall_rate_at_node=1.0,
    ... )
    array([ 2., 2., 12., 2., 2.])
    """
    h = np.asarray(h_at_node)
    q = np.asarray(q_at_node)

    n_nodes = h.size

    if q.shape != h.shape:
        raise ValueError(f"q_at_node shape {q.shape} must match h_at_node {h.shape}")

    if rainfall_rate_at_node is None:
        r = 0.0
    else:
        r = np.asarray(rainfall_rate_at_node)
        if r.ndim != 0 and r.shape != h.shape:
            raise ValueError(
                f"rainfall_rate_at_node shape {r.shape} must be scalar or {h.shape}"
            )

    where = np.asarray(where, dtype=bool)

    if out is None:
        out = np.empty_like(h, dtype=float)
        out[...] = h

    if not (out.size == n_nodes and out.ndim == 1):
        raise ValueError(f"out array must be of length {n_nodes}, got {out.size}")

    np.subtract(r, q, out=out, where=where)
    np.multiply(out, dt, out=out, where=where)
    np.add(h, out, out=out, where=where)
    np.maximum(out, 0.0, out=out, where=where)

    return out


def calc_grad_at_link(
    values: ArrayLike,
    *,
    length_of_link: ArrayLike = 1.0,
    axis: int = -1,
    where: ArrayLike | bool = True,
    out: NDArray | None = None,
) -> NDArray:
    """Calculate gradients at raster links.

    The function computes the difference in ``values`` across each link of a
    raster grid in the canonical Landlab ordering.

    Parameters
    ----------
    values : array_like of float, shape (n_rows, n_cols)
        Node values on a regular raster grid.
    length_of_link : array_like of float or scalar, optional
        Length of each link. Can be a scalar applied to all links or an array
        of shape ``(n_links,)`` matching the output.
    where : array_like of bool, optional
        Boolean mask of shape ``(n_links,)`` selecting which link gradients
        to compute. Where ``False``, the corresponding entries in ``out`` are
        left unchanged. If omitted, all links are computed.
    out : ndarray, optional
        Optional output buffer of shape ``(n_links,)`` into which results are
        written. If provided, it must have the correct size and dtype
        compatible with ``values``.

    Returns
    -------
    out : ndarray of float
        Array of link gradients in Landlab order. The dtype is determined by
        the input.

    Examples
    --------
    >>> z = [
    ...     [0, 2, 4],
    ...     [1, 3, 5],
    ... ]
    >>> calc_grad_at_link(z, axis=0)
    array([[1., 1., 1.]])
    >>> calc_grad_at_link(z, axis=1)
    array([[2., 2.],
           [2., 2.]])
    >>> calc_grad_at_link([0.0, 1.0, 2.0, 4.0, 8.0])
    array([1., 1., 2., 4.])
    """
    values = np.asarray(values)
    length_of_link = np.asarray(length_of_link)

    if values.ndim == 0:
        raise ValueError("values must be at least 1d")

    axis = np.core.numeric.normalize_axis_index(axis, values.ndim)

    values = np.moveaxis(values, axis, -1)

    out_shape = values.shape[:-1] + (values.shape[-1] - 1,)

    if out is None:
        out = np.empty(out_shape, dtype=float)
        out_provided_by_us = True
    else:
        out = np.moveaxis(out, axis, -1)
        if out.shape != out_shape:
            raise ValueError("shape mismatch for out")
        out_provided_by_us = False

    np.subtract(
        values[..., 1:],
        values[..., :-1],
        out=out,
        where=where,
    )
    np.divide(out, length_of_link, out=out, where=where)

    if out_provided_by_us:
        return np.moveaxis(out, -1, axis)
    return out
