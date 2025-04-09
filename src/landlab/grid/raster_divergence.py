#! /usr/bin/env python
"""Calculate flux divergence on a raster grid."""
import numpy as np

from landlab.grid.ext.raster_divergence import _calc_flux_div_at_node
from landlab.grid.ext.raster_divergence import _calc_net_face_flux_at_cell
from landlab.utils.decorators import use_field_name_or_array


@use_field_name_or_array("link")
def calc_flux_div_at_node(grid, unit_flux, out=None):
    """Calculate divergence of link-based fluxes at nodes.

    Given a flux per unit width across each face in the grid, calculate the net
    outflux (or influx, if negative) divided by cell area, at each node (zero
    or "out" value for nodes without cells).

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    unit_flux : ndarray or field name
        Flux per unit width along links (x number of links).

    Returns
    -------
    ndarray (x number of nodes)
        Flux divergence at nodes.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_divergence import calc_flux_div_at_node

    >>> grid = RasterModelGrid((3, 4), xy_spacing=10.0)
    >>> z = grid.add_zeros("topographic__elevation", at="node")
    >>> z[5] = 50.0
    >>> z[6] = 36.0
    >>> z.reshape(grid.shape)
    array([[  0.,   0.,   0.,   0.],
           [  0.,  50.,  36.,   0.],
           [  0.,   0.,   0.,   0.]])
    >>> grads = grid.calc_grad_at_link(z)
    >>> grads[grid.horizontal_links].reshape((3, 3))
    array([[ 0. ,  0. ,  0. ],
           [ 5. , -1.4, -3.6],
           [ 0. ,  0. ,  0. ]])
    >>> grads[grid.vertical_links].reshape((2, 4))
    array([[ 0. ,  5. ,  3.6,  0. ],
           [ 0. , -5. , -3.6,  0. ]])
    >>> calc_flux_div_at_node(grid, -grads).reshape(grid.shape)
    array([[0.  ,  0.  ,  0.  ,  0.  ],
           [0.  ,  1.64,  0.94,  0.  ],
           [0.  ,  0.  ,  0.  ,  0.  ]])

    >>> grid.set_status_at_node_on_edges(right=grid.BC_NODE_IS_CLOSED)
    >>> grid.set_status_at_node_on_edges(top=grid.BC_NODE_IS_CLOSED)
    >>> unit_flux_at_links = grid.zeros(at="link")
    >>> unit_flux_at_links[grid.active_links] = -grads[grid.active_links]
    >>> calc_flux_div_at_node(grid, unit_flux_at_links).reshape(grid.shape)
    array([[0.  ,  0.  ,  0.  ,  0.  ],
           [0.  ,  1.14,  0.22,  0.  ],
           [0.  ,  0.  ,  0.  ,  0.  ]])
    >>> _ = grid.add_field("neg_grad_at_link", -grads, at="link")
    >>> calc_flux_div_at_node(grid, "neg_grad_at_link").reshape(grid.shape)
    array([[0.  ,  0.  ,  0.  ,  0.  ],
           [0.  ,  1.64,  0.94,  0.  ],
           [0.  ,  0.  ,  0.  ,  0.  ]])

    Notes
    -----
    Performs a numerical flux divergence operation on nodes.

    :meta landlab: info-node, gradient
    """
    if unit_flux.size != grid.number_of_links:
        raise ValueError(
            f"unit_flux must be num links long (expected {grid.number_of_links},"
            f" got {unit_flux.size})"
        )
    if out is None:
        out = grid.empty(at="node")
    elif out.size != grid.number_of_nodes:
        raise ValueError("output buffer length mismatch with number of nodes")

    _calc_flux_div_at_node(grid.shape, (grid.dx, grid.dy), unit_flux, out)

    return out


def calc_net_face_flux_at_cell(grid, unit_flux_at_face, out=None):
    """Calculate net face fluxes at cells.

    Given a flux per unit width across each face in the grid, calculate the net
    outflux (or influx, if negative) at each cell.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    unit_flux_at_faces : ndarray or field name
        Flux per unit width associated with faces.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray (x number of cells)
        Net flux at cells.
    """
    if len(unit_flux_at_face) != grid.number_of_faces:
        raise ValueError("Parameter unit_flux_at_face must be num faces long")

    if out is None:
        out = grid.empty(at="cell")

    _calc_net_face_flux_at_cell(
        grid.shape, (grid.dx, grid.dy), np.asarray(unit_flux_at_face), out
    )

    return out
