#! /usr/bin/env python
"""Grid element mappers that are specific to raster grids.

Mapping functions unique to raster grids
++++++++++++++++++++++++++++++++++++++++

.. autosummary::

    ~landlab.grid.raster_mappers.map_sum_of_inlinks_to_node
    ~landlab.grid.raster_mappers.map_mean_of_inlinks_to_node
    ~landlab.grid.raster_mappers.map_max_of_inlinks_to_node
    ~landlab.grid.raster_mappers.map_min_of_inlinks_to_node
    ~landlab.grid.raster_mappers.map_sum_of_outlinks_to_node
    ~landlab.grid.raster_mappers.map_mean_of_outlinks_to_node
    ~landlab.grid.raster_mappers.map_max_of_outlinks_to_node
    ~landlab.grid.raster_mappers.map_min_of_outlinks_to_node
    ~landlab.grid.raster_mappers.map_mean_of_links_to_node
    ~landlab.grid.raster_mappers.map_mean_of_horizontal_links_to_node
    ~landlab.grid.raster_mappers.map_mean_of_horizontal_active_links_to_node
    ~landlab.grid.raster_mappers.map_mean_of_vertical_links_to_node
    ~landlab.grid.raster_mappers.map_mean_of_vertical_active_links_to_node
"""

import numpy as np

from .ext.raster_mappers import (
    map_max_of_link_nodes_to_link as _map_max_of_link_nodes_to_link,
)


def _node_out_link_ids(shape):
    """Links leaving each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.raster_mappers import _node_out_link_ids
    >>> (vert, horiz) = _node_out_link_ids((3, 4))
    >>> vert
    array([[ 3,  4,  5,  6],
           [10, 11, 12, 13],
           [-1, -1, -1, -1]])
    >>> horiz
    array([[ 0,  1,  2, -1],
           [ 7,  8,  9, -1],
           [14, 15, 16, -1]])
    """
    from ..graph.structured_quad.structured_quad import StructuredQuadGraphTopology

    layout = StructuredQuadGraphTopology(shape)

    node_horizontal_link_ids = np.empty(shape, int)
    node_horizontal_link_ids[:, :-1] = layout.horizontal_links.reshape(
        (shape[0], shape[1] - 1)
    )
    node_horizontal_link_ids[:, -1] = -1

    node_vertical_link_ids = np.empty(shape, int)
    node_vertical_link_ids[:-1, :] = layout.vertical_links.reshape(
        (shape[0] - 1, shape[1])
    )
    node_vertical_link_ids[-1, :] = -1

    return node_vertical_link_ids, node_horizontal_link_ids


def _node_in_link_ids(shape):
    """Links entering each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    tuple :
        Tuple of array of link IDs as (vertical_links, horizontal_links).

    Examples
    --------
    >>> from landlab.grid.raster_mappers import _node_in_link_ids
    >>> (vert, horiz) = _node_in_link_ids((3, 4))
    >>> vert
    array([[-1, -1, -1, -1],
           [ 3,  4,  5,  6],
           [10, 11, 12, 13]])
    >>> horiz
    array([[-1,  0,  1,  2],
           [-1,  7,  8,  9],
           [-1, 14, 15, 16]])
    """
    from ..graph.structured_quad.structured_quad import StructuredQuadGraphTopology

    layout = StructuredQuadGraphTopology(shape)

    node_horizontal_link_ids = np.empty(shape, int)
    node_horizontal_link_ids[:, 1:] = layout.horizontal_links.reshape(
        (shape[0], shape[1] - 1)
    )
    node_horizontal_link_ids[:, 0] = -1

    node_vertical_link_ids = np.empty(shape, int)
    node_vertical_link_ids[1:, :] = layout.vertical_links.reshape(
        (shape[0] - 1, shape[1])
    )
    node_vertical_link_ids[0, :] = -1

    return node_vertical_link_ids, node_horizontal_link_ids


def _number_of_links_per_node(shape):
    """Number of links touching each node.

    Parameters
    ----------
    shape : tuple of int
        Shape of grid of nodes.

    Returns
    -------
    ndarray :
        Array of number of links per node.

    Examples
    --------
    >>> from landlab.grid.raster_mappers import _number_of_links_per_node
    >>> _number_of_links_per_node((3, 4))
    array([[2, 3, 3, 2],
           [3, 4, 4, 3],
           [2, 3, 3, 2]])
    """
    from ..graph.structured_quad.structured_quad import StructuredQuadGraphTopology

    layout = StructuredQuadGraphTopology(shape)

    n_links_at_node = np.full(shape[0] * shape[1], 4, int)
    n_links_at_node[layout.perimeter_nodes] = 3
    n_links_at_node[layout.corner_nodes] = 2

    return n_links_at_node.reshape(shape)


def map_max_of_link_nodes_to_link(grid, value_at_node, out=None):
    """Map the maximum of a link's nodes to the link.

    map_max_of_link_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function evaluates the value of ``var_name`` at both the "to" and
    "from" node. The maximum value of the two node values is then mapped to
    the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    value_at_node : array or field name
        Values defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_max_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> _ = grid.add_field(
    ...     "z",
    ...     [
    ...         [0, 1, 2, 3],
    ...         [7, 6, 5, 4],
    ...         [8, 9, 10, 11],
    ...     ],
    ...     at="node",
    ... )
    >>> map_max_of_link_nodes_to_link(grid, "z")
    array([ 1,  2,  3,  7,  6,  5,  4,  7,  6,  5,  8,  9, 10, 11,  9, 10, 11])

    >>> values_at_links = grid.empty(at="link", dtype=grid.at_node["z"].dtype)
    >>> rtn = map_max_of_link_nodes_to_link(grid, "z", out=values_at_links)
    >>> values_at_links
    array([ 1,  2,  3,  7,  6,  5,  4,  7,  6,  5,  8,  9, 10, 11,  9, 10, 11])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if isinstance(value_at_node, str):
        value_at_node = grid.at_node[value_at_node]

    if out is None:
        out = grid.empty(at="link", dtype=value_at_node.dtype)

    _map_max_of_link_nodes_to_link(out, value_at_node, grid.shape)

    return out


def map_sum_of_inlinks_to_node(grid, var_name, out=None):
    """Map the sum of links entering a node to the node.

    map_sum_of_inlinks_to_node takes an array *at the links* and finds the
    inlink values for each node in the grid. it sums the inlinks and returns
    values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_sum_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_sum_of_inlinks_to_node(rmg, "z")
    array([  0.,   0.,   1.,   2.,   3.,  11.,  13.,  15.,  10.,  25.,  27.,
            29.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)

    south, west = _node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = values_at_links[south] + values_at_links[west]

    return out


def map_mean_of_inlinks_to_node(grid, var_name, out=None):
    """Map the mean of links entering a node to the node.

    map_mean_of_inlinks_to_node takes an array *at the links* and finds the
    inlink values for each node in the grid. It finds the average of
    the inlinks and returns values at the nodes.

    This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_mean_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_mean_of_inlinks_to_node(rmg, "z")
    array([  0. ,   0. ,   0.5,   1. ,   1.5,   5.5,   6.5,   7.5,   5. ,
            12.5,  13.5,  14.5])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    south, west = _node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = 0.5 * (values_at_links[south] + values_at_links[west])

    return out


def map_max_of_inlinks_to_node(grid, var_name, out=None):
    """Map the maximum of links entering a node to the node.

    map_max_of_inlinks_to_node takes an array *at the links* and finds the
    inlink values for each node in the grid. it finds the maximum value at the
    the inlinks and returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_max_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_max_of_inlinks_to_node(rmg, "z")
    array([  0.,   0.,   1.,   2.,
             3.,   7.,   8.,   9.,
            10.,  14.,  15.,  16.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    south, west = _node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.maximum(values_at_links[south], values_at_links[west])

    return out


def map_min_of_inlinks_to_node(grid, var_name, out=None):
    """Map the minimum of links entering a node to the node.

    map_min_of_inlinks_to_node takes an array *at the links* and finds the
    inlink values for each node in the grid. it finds the minimum value at the
    the inlinks and returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_min_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_min_of_inlinks_to_node(rmg, "z")
    array([  0.,   0.,   0.,   0.,   0.,   4.,   5.,   6.,   0.,  11.,  12.,
            13.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    south, west = _node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.minimum(values_at_links[south], values_at_links[west])

    return out


def map_sum_of_outlinks_to_node(grid, var_name, out=None):
    """Map the sum of links leaving a node to the node.

    map_sum_of_outlinks_to_node takes an array *at the links* and finds the
    outlink values for each node in the grid. it sums the outlinks and returns
    values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_sum_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_sum_of_outlinks_to_node(rmg, "z")
    array([  3.,  5.,  7.,   6.,  17.,  19.,  21.,  13.,  14.,  15.,  16.,
             0.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    north, east = _node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = values_at_links[north] + values_at_links[east]

    return out


def map_mean_of_outlinks_to_node(grid, var_name, out=None):
    """Map the mean of links leaving a node to the node.

    map_mean_of_outlinks_to_node takes an array *at the links* and finds the
    outlink values for each node in the grid. it finds the average of
    the outlinks and returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_mean_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_mean_of_outlinks_to_node(rmg, "z")
    array([  1.5,   2.5,   3.5,   3. ,   8.5,   9.5,  10.5,   6.5,   7. ,
             7.5,   8. ,   0. ])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    north, east = _node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = 0.5 * (values_at_links[north] + values_at_links[east])

    return out


def map_max_of_outlinks_to_node(grid, var_name, out=None):
    """Map the max of links leaving a node to the node.

    map_max_of_outlinks_to_node takes an array *at the links* and finds the
    outlink values for each node in the grid. it finds the maximum value at the
    the outlinks and returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_max_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_max_of_outlinks_to_node(rmg, "z")
    array([  3.,   4.,   5.,   6.,  10.,  11.,  12.,  13.,  14.,  15.,  16.,
             0.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    north, east = _node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.maximum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_min_of_outlinks_to_node(grid, var_name, out=None):
    """Map the min of links leaving a node to the node.

    map_min_of_outlinks_to_node takes an array *at the links* and finds the
    outlink values for each node in the grid. It finds the minimum value at the
    the outlinks and returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_min_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_min_of_outlinks_to_node(rmg, "z")
    array([0.,  1.,  2.,  0.,  7.,  8.,  9.,  0.,  0.,  0.,  0.,  0.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)
    north, east = _node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.minimum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_mean_of_links_to_node(grid, var_name, out=None):
    """Map the mean of links touching a node to the node.

    map_mean_all_links_to_node takes an array *at the links* and finds the
    average of all ~existing~ link neighbor values for each node in the grid.
    it returns values at the nodes.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_mean_of_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_mean_of_links_to_node(rmg, "z")
    array([  1.5       ,   1.66666667,   2.66666667,   4.        ,
             6.66666667,   7.5       ,   8.5       ,   9.33333333,
            12.        ,  13.33333333,  14.33333333,  14.5       ])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    values_at_links = np.append(values_at_links, 0)

    north, east = _node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    south, west = _node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)

    number_of_links = _number_of_links_per_node(grid.shape)
    number_of_links = number_of_links.reshape(number_of_links.size)
    number_of_links.astype(float, copy=False)
    out[:] = (
        values_at_links[north]
        + values_at_links[east]
        + values_at_links[south]
        + values_at_links[west]
    ) / number_of_links

    return out


def map_mean_of_horizontal_links_to_node(grid, var_name, out=None):
    """Map the mean of links in the x direction touching a node to the node.

    map_mean_of_horizontal_links_to_node takes an array *at the links* and
    finds the average of all horizontal (x-direction) link neighbor values
    for each node in the grid.
    It returns an array at the nodes of the mean of these values. If a link
    is absent, it is ignored.
    Note that here a positive returned value means flux to the east, and
    a negative to the west.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_mean_of_horizontal_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_mean_of_horizontal_links_to_node(rmg, "z")
    array([  0. ,   0.5,   1.5,   2. ,   7. ,   7.5,   8.5,   9. ,  14. ,
            14.5,  15.5,  16. ])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    hoz_links = grid.links_at_node[:, [0, 2]]
    hoz_link_dirs = np.fabs(grid.link_dirs_at_node[:, [0, 2]])
    # ^retain "true" directions of links
    valid_links = values_at_links[hoz_links] * hoz_link_dirs  # invalids = 0
    num_valid_links = hoz_link_dirs.sum(axis=1)
    np.divide(valid_links.sum(axis=1), num_valid_links, out=out)
    return out


def map_mean_of_horizontal_active_links_to_node(grid, var_name, out=None):
    """Map the mean of active links in the x direction touching node to the
    node.

    map_mean_of_horizontal_active_links_to_node takes an array *at the links*
    and finds the average of all horizontal (x-direction) link neighbor values
    for each node in the grid.
    It returns an array at the nodes of the mean of these values. If a link
    is absent, it is ignored. If a node has no active links, it receives 0.
    Note that here a positive returned value means flux to the east, and
    a negative to the west.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import (
    ...     map_mean_of_horizontal_active_links_to_node,
    ... )
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", -np.arange(17, dtype=float), at="link")
    >>> rmg.status_at_node[rmg.nodes_at_left_edge] = rmg.BC_NODE_IS_CLOSED
    >>> map_mean_of_horizontal_active_links_to_node(rmg, "z")
    array([ 0. ,  0. ,  0. ,  0. ,  0. , -8. , -8.5, -9. ,  0. ,  0. ,  0. ,
            0. ])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.zeros(centering="node", dtype=float)
    else:
        out.fill(0.0)

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    hoz_links = grid.links_at_node[:, [0, 2]]
    hoz_link_dirs = np.fabs(grid.active_link_dirs_at_node[:, [0, 2]])
    # ^retain "true" directions of links; no inactives now
    valid_links = values_at_links[hoz_links] * hoz_link_dirs  # invalids = 0
    num_valid_links = hoz_link_dirs.sum(axis=1)
    good_nodes = num_valid_links != 0
    out[good_nodes] = valid_links.sum(axis=1)[good_nodes] / num_valid_links[good_nodes]
    return out


def map_mean_of_vertical_links_to_node(grid, var_name, out=None):
    """Map the mean of links in the y direction touching a node to the node.

    map_mean_of_vertical_links_to_node takes an array *at the links* and
    finds the average of all vertical (y-direction) link neighbor values
    for each node in the grid.
    It returns an array at the nodes of the mean of these values. If a link
    is absent, it is ignored.
    Note that here a positive returned value means flux to the north, and
    a negative to the south.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import map_mean_of_vertical_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", np.arange(17.0), at="link")
    >>> map_mean_of_vertical_links_to_node(rmg, "z")
    array([  3. ,   4. ,   5. ,   6. ,   6.5,   7.5,   8.5,   9.5,  10. ,
            11. ,  12. ,  13. ])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(centering="node")

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    vert_links = grid.links_at_node[:, [1, 3]]
    vert_link_dirs = np.fabs(grid.link_dirs_at_node[:, [1, 3]])
    # ^retain "true" directions of links
    valid_links = values_at_links[vert_links] * vert_link_dirs  # invalids = 0
    num_valid_links = vert_link_dirs.sum(axis=1)
    np.divide(valid_links.sum(axis=1), num_valid_links, out=out)
    return out


def map_mean_of_vertical_active_links_to_node(grid, var_name, out=None):
    """Map the mean of active links in the y direction touching node to the
    node.

    map_mean_of_vertical_active_links_to_node takes an array *at the links*
    and finds the average of all vertical (y-direction) link neighbor values
    for each node in the grid.
    It returns an array at the nodes of the mean of these values. If a link
    is absent, it is ignored. If a node has no active links, it receives 0.
    Note that here a positive returned value means flux to the north, and
    a negative to the south.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.raster_mappers import (
    ...     map_mean_of_vertical_active_links_to_node,
    ... )
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field("z", -np.arange(17, dtype=float), at="link")
    >>> rmg.status_at_node[rmg.nodes_at_bottom_edge] = rmg.BC_NODE_IS_CLOSED
    >>> map_mean_of_vertical_active_links_to_node(rmg, "z")
    array([  0.,   0.,   0.,   0.,   0., -11., -12.,   0.,   0., -11., -12.,
             0.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.zeros(centering="node", dtype=float)
    else:
        out.fill(0.0)

    if type(var_name) is str:
        values_at_links = grid.at_link[var_name]
    else:
        values_at_links = var_name
    vert_links = grid.links_at_node[:, [1, 3]]
    vert_link_dirs = np.fabs(grid.active_link_dirs_at_node[:, [1, 3]])
    # ^retain "true" directions of links; no inactives now
    valid_links = values_at_links[vert_links] * vert_link_dirs  # invalids = 0
    num_valid_links = vert_link_dirs.sum(axis=1)
    good_nodes = num_valid_links != 0
    out[good_nodes] = valid_links.sum(axis=1)[good_nodes] / num_valid_links[good_nodes]
    return out


def map_link_vector_components_to_node_raster(grid, data_at_link):
    """Map (x,y) vector components of data_at_link onto nodes.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rmg = RasterModelGrid((3, 4))
    >>> link_data = np.arange(rmg.number_of_links)
    >>> x, y = map_link_vector_components_to_node_raster(rmg, link_data)
    >>> x[5:7]
    array([7.5,  8.5])
    >>> y[5:7]
    array([7.5,  8.5])
    """
    x = grid.map_mean_of_horizontal_links_to_node(data_at_link)
    y = grid.map_mean_of_vertical_links_to_node(data_at_link)
    return x, y
