import numpy as np
from six.moves import range

from ..core.utils import make_optional_arg_into_id_array


def neighbor_active_link_at_cell(grid, inds, *args):
    """neighbor_active_link_at_cell(grid, link_ids [, cell_ids])

    Return an array of the active link ids for neighbors of *cell_id* cells.
    *link_ids* is an index into the links of a cell as measured
    clockwise starting from the south.

    If *cell_ids* is not given, return neighbors for all cells in the grid.

    Parameters
    ----------
    grid : RasterModelGrid
        Source grid.
    link_inds : array_like
        IDs of links
    cell_ids : array_like, optional
        IDs of cells for which to get links

    """
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]
    links = grid._active_links_at_node(node_ids).T

    if not isinstance(inds, np.ndarray):
        inds = np.array(inds)

    return links[range(len(cell_ids)), inds]


def neighbor_node_at_cell(grid, inds, *args):
    """ node_id_of_cell_neighbor(grid, neighbor_ids [, cell_ids])

    Return an array of the node ids for neighbors of *cell_id* cells.
    *neighbor_ids* is an index into the neighbors of a cell as measured
    clockwise starting from the south.

    If *cell_ids* is not given, return neighbors for all cells in the grid.

    Parameters
    ----------
    grid : RasterModelGrid
        Input grid.
    neighbor_ids : array_like
        IDs of the neighbor nodes.
    cell_ids : array_like, optional
        IDs of cell about which to get neighbors.

    Returns
    -------
    ndarray
        Node IDs for given neighbors of cells.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_funcs import neighbor_node_at_cell
    >>> grid = RasterModelGrid((4, 5), xy_spacing=1.0)
    >>> neighbor_node_at_cell(grid, 0, 0)
    array([1])

    Get the lower and the the upper neighbors for all the cells.

    >>> neighbor_node_at_cell(grid, 0)
    array([1, 2, 3, 6, 7, 8])
    >>> neighbor_node_at_cell(grid, 2)
    array([11, 12, 13, 16, 17, 18])

    As an alternative to the above, use fancy-indexing to get both sets of
    neighbors with one call.

    >>> neighbor_node_at_cell(grid, np.array([0, 2]), [1, 4])
    array([[ 2, 12],
           [ 7, 17]])
    """
    cell_ids = make_optional_arg_into_id_array(grid.number_of_cells, *args)
    node_ids = grid.node_at_cell[cell_ids]
    neighbors = grid.active_adjacent_nodes_at_node[node_ids]

    if not isinstance(inds, np.ndarray):
        inds = np.array(inds)

    # return neighbors[range(len(cell_ids)), 3 - inds]
    return np.take(np.take(neighbors, range(len(cell_ids)), axis=0), 3 - inds, axis=1)


def calculate_slope_aspect_bfp(xs, ys, zs):
    """Calculate slope and aspect.

    .. codeauthor:: Katy Barnhart <katherine.barnhart@colorado.edu>

    Fits a plane to the given N points with given *xs*, *ys*, and *zs* values
    using single value decomposition.

    Returns a tuple of (*slope*, *aspect*) based on the normal vector to the
    best fit plane.

    .. note::

        This function does not check if the points fall on a line, rather
        than a plane.
    """
    if not len(xs) == len(ys) == len(zs):
        raise ValueError("array must be the same length")

    # step 1: subtract the centroid from the points
    # step 2: create a 3XN matrix of the points for SVD
    # in python, the unit normal to the best fit plane is
    # given by the third column of the U matrix.
    mat = np.vstack((xs - np.mean(xs), ys - np.mean(ys), zs - np.mean(zs)))
    U, _, _ = np.linalg.svd(mat)
    normal = U[:, 2]

    # step 3: calculate the aspect
    asp = 90.0 - np.degrees(np.arctan2(normal[1], normal[0]))
    asp = asp % 360.0

    # step 4: calculate the slope
    slp = 90.0 - np.degrees(np.arcsin(normal[2]))

    return slp, asp


def find_nearest_node(rmg, coords, mode="raise"):
    """Find the node nearest a point.

    Find the index to the node nearest the given x, y coordinates.
    Coordinates are provided as numpy arrays in the *coords* tuple.
    *coords* is tuple of coordinates, one for each dimension.

    The *mode* keyword to indicate what to do if a point is outside of the
    grid. Valid choices are the same as that used with the numpy function
    `ravel_multi_index`.

    A point is considered to be outside of the grid if it is outside the
    perimeter of the grid by one half of the grid spacing.

    Parameters
    ----------
    rmg : RasterModelGrid
        The source grid.
    coords : tuple
        Coordinates of point as (x, y)
    mode : {'raise', 'wrap', 'clip'}, optional
        Action to take if a point is outside of the grid.

    Returns
    -------
    array_like :
        Indices of the nodes nearest the given coordinates.

    Examples
    --------
    Create a grid of 4 by 5 nodes with unit spacing.

    >>> import landlab
    >>> from landlab.grid.raster_funcs import find_nearest_node
    >>> rmg = landlab.RasterModelGrid((4, 5))

    The points can be either a tuple of scalars or of arrays.

    >>> find_nearest_node(rmg, (0.2, 0.6))
    5
    >>> find_nearest_node(rmg, (np.array([1.6, 3.6]), np.array([2.3, .7])))
    array([12,  9])

    The *mode* keyword indicates what to do if a point is outside of the
    grid.

    >>> find_nearest_node(rmg, (-0.6, 0.6), mode='raise')
    Traceback (most recent call last):
        ...
    ValueError: invalid entry in coordinates array
    >>> find_nearest_node(rmg, (-0.6, 0.6), mode='clip')
    5
    >>> find_nearest_node(rmg, (-0.6, 0.6), mode='wrap')
    9
    """
    if isinstance(coords[0], np.ndarray):
        return _find_nearest_node_ndarray(rmg, coords, mode=mode)
    else:
        return find_nearest_node(
            rmg, (np.array(coords[0]), np.array(coords[1])), mode=mode
        )


def _find_nearest_node_ndarray(rmg, coords, mode="raise"):
    """Find the node nearest to a point.

    Parameters
    ----------
    rmg : RasterModelGrid
        A RasterModelGrid.
    coords : tuple of float
        Coordinates of test points as *x*, then *y*.
    mode : {'raise', 'wrap', 'clip'}, optional
        What to do with out-of-bounds indices (as with
        numpy.ravel_multi_index).

    Returns
    -------
    ndarray
        Nodes that are closest to the points.

    Examples
    --------
    >>> from landlab.grid.raster_funcs import _find_nearest_node_ndarray
    >>> from landlab import RasterModelGrid
    >>> import numpy as np
    >>> grid = RasterModelGrid((4, 5))
    >>> _find_nearest_node_ndarray(grid, (.25, 1.25))
    5
    >>> _find_nearest_node_ndarray(grid, (.75, 2.25))
    11

    >>> grid = RasterModelGrid((4, 5), xy_spacing=(3, 4))
    >>> _find_nearest_node_ndarray(grid, (3.1, 4.1))
    6
    """
    column_indices = np.int_(np.around((coords[0] - rmg.node_x[0]) / rmg.dx))
    row_indices = np.int_(np.around((coords[1] - rmg.node_y[0]) / rmg.dy))

    return rmg.grid_coords_to_node_id(row_indices, column_indices, mode=mode)


def _value_is_in_bounds(value, bounds):
    """Check if a value is within bounds.

    Parameters
    ----------
    value : float or ndarray
        The test value.
    bounds : (lower, upper)
        The lower and upper bounds.

    Returns
    -------
    bool
        ``True`` if the value is within the bounds. Otherwise, ``False``.

    Examples
    --------
    >>> from landlab.grid.raster_funcs import _value_is_in_bounds
    >>> import numpy as np
    >>> _value_is_in_bounds(.5, (0, 1))
    True
    >>> _value_is_in_bounds(1, (0, 1))
    False
    >>> _value_is_in_bounds(0, (0, 1))
    True
    >>> _value_is_in_bounds(np.array((0, 1)), (0, 1))
    array([ True, False], dtype=bool)
    """
    dummy = value >= bounds[0]
    dummy &= value < bounds[1]
    return dummy


def _value_is_within_axis_bounds(rmg, value, axis):
    """Check if a value is within the bounds of a grid axis.

    Parameters
    ----------
    rmg : RasterModelGrid
        A RasterModelGrid.
    value : float
        The test value.
    axis : int
        The axis.

    Returns
    -------
    bool
        ``True`` if the value is within the axis bounds. Otherwise, ``False``.

    Examples
    --------
    >>> from landlab.grid.raster_funcs import _value_is_within_axis_bounds
    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid((4, 5))
    >>> _value_is_within_axis_bounds(rmg, 3.1, 0)
    False
    >>> _value_is_within_axis_bounds(rmg, 2.9, 0)
    True
    >>> _value_is_within_axis_bounds(rmg, 4.1, 1)
    False
    >>> _value_is_within_axis_bounds(rmg, 3.9, 1)
    True
    """
    axis_coord = rmg.node_axis_coordinates(axis)
    return _value_is_in_bounds(value, (axis_coord[0], axis_coord[-1]))


def is_coord_on_grid(rmg, coords, axes=(0, 1)):
    """Check if coordinates are contained on a grid.

    Parameters
    ----------
    rmg : RasterModelGrid
        Source grid.
    coords : tuple
        Coordinates of point as (x, y)
    axes : tuple, optional
        Check bounds only on a particular axis

    Examples
    --------
    Create a grid that ranges from x=0 to x=4, and y=0 to y=3.

    >>> from landlab import RasterModelGrid
    >>> from landlab.grid.raster_funcs import is_coord_on_grid
    >>> grid = RasterModelGrid((4, 5))
    >>> is_coord_on_grid(grid, (3.999, 2.999))
    True

    Check two points with one call. Numpy broadcasting rules apply for the
    point coordinates.

    >>> is_coord_on_grid(grid, ([3.9, 4.1], 2.9))
    array([ True, False], dtype=bool)

    >>> is_coord_on_grid(grid, ([3.9, 4.1], 2.9), axes=(0, ))
    array([ True,  True], dtype=bool)
    """
    coords = np.broadcast_arrays(*coords)

    is_in_bounds = _value_is_within_axis_bounds(rmg, coords[1 - axes[0]], axes[0])
    for axis in axes[1:]:
        is_in_bounds &= _value_is_within_axis_bounds(rmg, coords[1 - axis], axis)

    return is_in_bounds
