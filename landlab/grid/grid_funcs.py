"""
Utility functions that operate on landlab grids.
------------------------------------------------

"""
import numpy as np
from .cfuncs import _get_closest_nodes


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
        np.multiply(
            (
                (
                    grid.node_x[grid.node_at_link_head]
                    - grid.node_x[grid.node_at_link_tail]
                )
                / grid.length_of_link
            ),
            link_values,
        ),
        np.multiply(
            (
                (
                    grid.node_y[grid.node_at_link_head]
                    - grid.node_y[grid.node_at_link_tail]
                )
                / grid.length_of_link
            ),
            link_values,
        ),
    )


def find_nearest_node(grid, coords):
    """Node nearest a point or array of points.

    Find the index to the node(s) nearest the given x, y coordinates.
    Coordinates are provided as numpy arrays in the *coords* tuple.

    Returns the indices of the nodes nearest the given coordinates.

    Parameters
    ----------
    grid : grid object
        Any Landlab grid object.
    coords : tuple of array-like
        Coordinates of points; (x, y). Note this ordering is not LL standard.

    Returns
    -------
    array-like
        IDs of the nearest nodes.

    Notes
    -----
    This base function is fairly slow. For a raster, use the equivalent, much
    faster, raster_funcs.find_nearest_node.

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> hmg = HexModelGrid((4, 5))
    >>> find_nearest_node(hmg, [-1., -1.])
    array([0])
    >>> find_nearest_node(hmg, (np.array([3.1, 4.2]), np.array([0.2, 1.7])))
    array([ 2, 15])

    LLCATS: NINF SUBSET
    """
    x = coords[0]
    y = coords[1]
    # Now, coerce inputs into standard form to feed the cfunc:
    if not isinstance(x, np.ndarray):
        if type(x) in (float, int):
            x = [
                x,
            ]
        x = np.array(x, dtype=float)
    else:
        if not x.dtype is np.dtype(float):
            x = x.astype(float)
    if not isinstance(y, np.ndarray):
        if type(y) in (float, int):
            y = [
                y,
            ]
        y = np.array(y, dtype=float)
    else:
        if not y.dtype is np.dtype(float):
            y = y.astype(float)
    out = np.empty_like(x, dtype=int)
    _get_closest_nodes(out, grid.node_x, grid.node_y, x, y)
    return out
