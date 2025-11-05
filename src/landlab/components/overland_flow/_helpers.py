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
    nodes: ArrayLike | None = None,
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
    nodes : array_like of int, optional
        Node indices to update (must be unique). If ``None``, update all nodes.
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
    q = np.asarray(q_at_node)
    h = np.asarray(h_at_node)

    if q.shape != h.shape:
        raise ValueError(f"q_at_node shape {q.shape} must match h_at_node {h.shape}")

    if rainfall_rate_at_node is None:
        r = 0.0
        is_scalar_rainfall = True
    else:
        r = np.asarray(rainfall_rate_at_node)
        is_scalar_rainfall = np.ndim(r) == 0
        if not is_scalar_rainfall and r.shape != h.shape:
            raise ValueError(
                f"rainfall_rate_at_node shape {r.shape} must be scalar or {h.shape}"
            )

    if nodes is None:
        idx = slice(None)
    else:
        idx = np.asarray(nodes, dtype=np.intp)

    if out is None:
        out = np.empty_like(h, dtype=np.result_type(h, q, r, float))
        out[...] = h
    elif nodes is None or nodes is Ellipsis:
        if out.shape != h.shape:
            raise ValueError(f"out shape {out.shape} must match {h.shape}")
        out[...] = h

    if is_scalar_rainfall:
        tmp = out[idx] + (r - q[idx]) * dt
    else:
        tmp = out[idx] + (r[idx] - q[idx]) * dt

    np.maximum(tmp, 0.0, out=tmp)
    out[idx] = tmp

    return out
