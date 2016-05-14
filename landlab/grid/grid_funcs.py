"""Utility functions that operate on landlab grids."""


import numpy as np
from six.moves import range


def resolve_values_on_active_links(grid, active_link_values):
    """Resolve active-link values into x and y directions.

    Takes a set of values defined on active links, and returns those values
    resolved into the x and y directions.  Two link arrays are returned:
    x, then y.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    active_link_values : ndarray
        Values on active links.

    Returns
    -------
    tuple of ndarray
        Values resolved into x-component and y-component.
    """
    link_lengths = grid.length_of_link[grid.active_links]
    return (
        np.multiply(((grid.node_x[grid._activelink_tonode] -
                      grid.node_x[grid._activelink_fromnode]) /
                     link_lengths), active_link_values),
        np.multiply(((grid.node_y[grid._activelink_tonode] -
                      grid.node_y[grid._activelink_fromnode]) /
                     link_lengths), active_link_values))


def resolve_values_on_links(grid, link_values):
    """Resolve link values into x and y directions.

    Takes a set of values defined on active links, and returns those values
    resolved into the x and y directions.  Two link arrays are returned:
    x, then y.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    link_values : ndarray
        Values on links.

    Returns
    -------
    tuple of ndarray
        Values resolved into x-component and y-component.
    """
    return (
        np.multiply(((grid.node_x[grid.node_at_link_head] -
                      grid.node_x[grid.node_at_link_tail]) /
                     grid.length_of_link), link_values),
        np.multiply(((grid.node_y[grid.node_at_link_head] -
                      grid.node_y[grid.node_at_link_tail]) /
                     grid.length_of_link), link_values))


def calculate_flux_divergence_at_nodes(grid, active_link_flux, out=None):
    """Calculate flux divergence at grid nodes.

    Same as calculate_flux_divergence_at_active_cells, but works with and
    returns a list of net unit fluxes that corresponds to all nodes, rather
    than just active cells.

    Note that we don't compute net unit fluxes at
    boundary nodes (which don't have active cells associated with them, and
    often don't have cells of any kind, because they are on the perimeter),
    but simply return zeros for these entries. The advantage is that the
    caller can work with node-based arrays instead of active-cell-based
    arrays.

    Parameters
    ----------
    grid : ModelGrid
        A ModelGrid.
    active_link_flux : ndarray
        Fluxes at active links.
    out : ndarray, optional
        Buffer to hold the result.

    Returns
    -------
    ndarray
        Net unit fluxes at nodes.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.grid_funcs import calculate_flux_divergence_at_nodes

    >>> grid = RasterModelGrid((4, 5))
    >>> link_flux = np.ones(grid.number_of_active_links, dtype=float)
    >>> flux_at_node = calculate_flux_divergence_at_nodes(grid, link_flux)
    ...     # doctest: +NORMALIZE_WHITESPACE
    >>> flux_at_node
    array([ 0.,  1.,  1.,  1.,  0.,
            1.,  0.,  0.,  0., -1.,
            1.,  0.,  0.,  0., -1.,
            0., -1., -1., -1.,  0.])
    >>> flux_at_node[grid.core_nodes]
    array([ 0., 0., 0., 0., 0., 0.])

    This is *deprecated*. Instead use ``calc_flux_div_at_node`. Notice that
    fluxes at non-core nodes are handled differently. However, these boundary
    nodes don't have "flux" anyway and so should be ignored.

    >>> grid = RasterModelGrid((4, 5))
    >>> link_flux = grid.zeros(at='link')
    >>> link_flux[grid.active_links] = 1.
    >>> flux_at_node = grid.calc_flux_div_at_node(link_flux)
    >>> flux_at_node
    array([ 0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.,
            0.,  0.,  0.,  0.,  0.])
    >>> flux_at_node[grid.core_nodes]
    array([ 0., 0., 0., 0., 0., 0.])
    """
    assert len(active_link_flux) == grid.number_of_active_links, (
        "incorrect length of active_link_flux array")

    # If needed, create net_unit_flux array
    if out is None:
        out = grid.empty(at='node')
    out.fill(0.)
    net_unit_flux = out

    assert len(net_unit_flux) == grid.number_of_nodes

    # Create a flux array one item longer than the number of active links.
    # Populate it with flux times face width (so, total flux rather than
    # unit flux). Here, face_width is an array with one entry for each
    # active link, so we are multiplying the unit flux at each link by the
    # width of its corresponding face.
    flux = np.zeros(len(active_link_flux) + 1)
    flux[:len(active_link_flux)] = active_link_flux * grid.width_of_face

    # Next, we need to add up the incoming and outgoing fluxes.
    #
    # Notes:
    #    1) because "net flux" is defined as positive outward, we add the
    #       outflux and subtract the influx
    #    2) the loop is over the number of rows in the inlink/outlink
    #       matrices. This dimension is equal to the maximum number of links
    #       attached to a node, so should be of order 6 or 7 and won't
    #       generally increase with the number of nodes in the grid.
    #
    for i in range(np.size(grid._node_active_inlink_matrix, 0)):
        net_unit_flux += flux[grid._node_active_outlink_matrix[i][:]]
        net_unit_flux -= flux[grid._node_active_inlink_matrix[i][:]]

    # Now divide by cell areas ... where there are core cells.
    node_at_active_cell = grid.node_at_cell[grid.core_cells]
    net_unit_flux[node_at_active_cell] /= grid.area_of_cell[grid.core_cells]

    return net_unit_flux
