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


def calc_mean_of_neighbors(
    values: ArrayLike,
    *,
    weight: float = 0.5,
    axis: int = -1,
    out: NDArray | None = None,
):
    """Calculate a weighted mean of a value with its neighbors.

    The mean is calculated as,

        mean = weight * value_at_center + (1.0 - weight) * 0.5 * (value_at_left + value_at_right)

    Parameters
    ----------
    values : array_like
        Values to average.
    weight : float, optional
        Weighting value between 0 and 1. A value of zero
        means to average the two neighbors. A value of one
        means to use the center value.
    out : NDArray, optional
        Optional output buffer, of the same shape as ``values``,
        into which results are written.

    Returns
    -------
    out : ndarray of float
        Array of calculated means.

    Examples
    --------
    >>> calc_mean_of_neighbors([1, 1, 2, 3, 5], weight=1.0)
    array([1., 1., 2., 3., 5.])
    >>> calc_mean_of_neighbors([1, 1, 2, 3, 5], weight=0.0)
    array([1. , 1.5, 2. , 3.5, 5. ])
    >>> x = [
    ...     [0, 1, 3, 6],
    ...     [0, 1, 3, 6],
    ...     [0, 1, 3, 6],
    ...     [0, 1, 3, 6],
    ... ]
    >>> calc_mean_of_neighbors(x, weight=0.0, axis=1)
    array([[0. , 1.5, 3.5, 6. ],
           [0. , 1.5, 3.5, 6. ],
           [0. , 1.5, 3.5, 6. ],
           [0. , 1.5, 3.5, 6. ]])
    >>> calc_mean_of_neighbors(x, weight=0.0, axis=0)
    array([[0., 1., 3., 6.],
           [0., 1., 3., 6.],
           [0., 1., 3., 6.],
           [0., 1., 3., 6.]])
    """
    if not (0.0 <= weight <= 1.0):
        raise ValueError(f"weight must be between 0 and 1, got {weight}")

    values = np.asarray(values)
    original_shape = values.shape

    values = np.moveaxis(values, axis, -1)

    if out is None:
        out = np.empty_like(values, dtype=float)
    else:
        out = np.moveaxis(out, axis, -1)
        if out.shape != values.shape:
            raise ValueError(
                f"out must have the same shape as values {out.shape} != {values.shape}"
            )

    out[..., 1:-1] = weight * values[..., 1:-1] + (1.0 - weight) * 0.5 * (
        values[..., :-2] + values[..., 2:]
    )
    out[..., 0] = values[..., 0]
    out[..., -1] = values[..., -1]

    return np.moveaxis(out, -1, axis)


def calc_flow_height(
    z_at_node: ArrayLike,
    h_at_node: ArrayLike,
    *,
    axis: int = -1,
    where: ArrayLike | bool = True,
    out: NDArray | None = None,
):
    """Compute Bates-style flow height along one axis of an array.

    For adjacent pairs along the chosen `axis`, this computes

        max(z + h at i, z + h at i+1) - max(z at i, z at i+1)

    which corresponds to the Bates et al. (2010) shallow-water update using
    the higher water surface on the pair minus the higher bed elevation.

    The output has the same shape as the inputs except that the selected axis
    is reduced by one (because it represents link-wise values between neighbors).

    Parameters
    ----------
    z_at_node : array_like
        Bed elevation array.
    h_at_node : array_like
        Water depth array. Must have the same shape as `z_at_node`.
    axis : int, optional
        Axis along which to compute pairwise link values.
    where : bool or array_like of bool, optional
        Mask selecting elements to update (NumPy ufunc semantics). Must be
        broadcastable to the output shape. Where `False`, the corresponding
        entries in `out` are left unchanged. If `out` is not provided, those
        entries will be uninitialized.
    out : ndarray, optional
        Destination array. Must have the same shape as the result
        (i.e., input shape with the selected axis reduced by one). When
        provided, only entries where `where` is `True` are overwritten.

    Returns
    -------
    ndarray
        The result array. If `out` is provided, it is returned.

    Examples
    --------
    >>> z = [
    ...     [10.0, 11.0, 13.0],
    ...     [9.0, 10.0, 12.0],
    ... ]
    >>> h = [
    ...     [0.5, 0.2, 0.3],
    ...     [0.4, 0.1, 0.6],
    ... ]
    >>> calc_flow_height(z, h, axis=0).shape
    (1, 3)
    >>> calc_flow_height(z, h, axis=1).shape
    (2, 2)

    >>> where = [
    ...     [True, False],
    ...     [True, True],
    ... ]
    >>> out = np.zeros((2, 2))
    >>> calc_flow_height(z, h, axis=1, where=where, out=out)
    array([[0.2, 0. ],
           [0.1, 0.6]])
    """
    z = np.asarray(z_at_node)
    h = np.asarray(h_at_node)

    if z.shape != h.shape:
        raise ValueError(
            f"z_at_node and h_at_node must be the same shape, got {z_at_node.shape}"
            f" and {h_at_node.shape}."
        )

    axis = np.lib.array_utils.normalize_axis_index(axis, z.ndim)

    z = np.moveaxis(z_at_node, axis, -1)
    h = np.moveaxis(h_at_node, axis, -1)

    out_shape = z.shape[:-1] + (z.shape[-1] - 1,)

    if out is None:
        out = np.empty(out_shape, dtype=float)
    else:
        out = np.moveaxis(out, axis, -1)
        if out.shape != out_shape:
            raise ValueError("shape mismatch for out")

    z_max = np.empty_like(z[..., :-1], dtype=float)
    w_tail = np.empty_like(z_max)
    w_head = np.empty_like(z_max)

    np.add(z[..., :-1], h[..., :-1], where=where, out=w_tail)
    np.add(z[..., 1:], h[..., 1:], where=where, out=w_head)
    np.maximum(w_tail, w_head, where=where, out=out)
    np.maximum(z[..., :-1], z[..., 1:], where=where, out=z_max)
    np.subtract(out, z_max, where=where, out=out)

    return np.moveaxis(out, -1, axis)
