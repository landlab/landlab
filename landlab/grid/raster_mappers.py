#! /usr/bin/env python
"""Grid element mappers that are specific to raster grids."""

from __future__ import division

import numpy as np

from landlab.grid.structured_quad import links


def map_sum_of_inlinks_to_node(grid, var_name, out=None):
    """Map the sum of links entering a node to the node.

    map_inlink_sums_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it sums the inlinks and returns
    a field at the nodes with the same var_name as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_sum_of_inlinks_to_node(rmg, 'z')
    array([  0.,   8.,   9.,  10.,   0.,  12.,  14.,  16.,   4.,  19.,  21.,
            23.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)

    south, west = links._node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = values_at_links[south] + values_at_links[west]

    return out


def map_mean_of_inlinks_to_node(grid, var_name, out=None):
    """Map the mean of links entering a node to the node.

    map_inlink_average_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the average of
    the inlinks and returns a field at the nodes with the same var_name
    as the link field.

    This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_inlinks_to_node(rmg, 'z')
    array([  0. ,   4. ,   4.5,   5. ,   0. ,   6. ,   7. ,   8. ,   2. ,
             9.5,  10.5,  11.5])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = 0.5 * (values_at_links[south] + values_at_links[west])

    return out


def map_max_of_inlinks_to_node(grid, var_name, out=None):
    """Map the maximum of links entering a node to the node.

    map_max_inlink_value_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the maximum value at the
    the inlinks and returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_max_of_inlinks_to_node(rmg, 'z')
    array([  0.,   8.,   9.,  10.,   0.,  11.,  12.,  13.,   4.,  14.,  15.,
            16.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.maximum(values_at_links[south], values_at_links[west])

    return out


def map_min_of_inlinks_to_node(grid, var_name, out=None):
    """Map the minimum of links entering a node to the node.

    map_min_inlink_value_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the minimum value at the
    the inlinks and returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_min_of_inlinks_to_node(rmg, 'z')
    array([ 0.,  0.,  0.,  0.,  0.,  1.,  2.,  3.,  0.,  5.,  6.,  7.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.minimum(values_at_links[south], values_at_links[west])

    return out


def map_sum_of_outlinks_to_node(grid, var_name, out=None):
    """Map the sum of links leaving a node to the node.

    map_outlink_sums_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it sums the outlinks and returns
    a field at the nodes with the same var_name as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_sum_of_outlinks_to_node(rmg, 'z')
    array([  8.,  10.,  12.,   3.,  15.,  17.,  19.,   7.,  14.,  15.,  16.,
             0.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = values_at_links[north] + values_at_links[east]

    return out


def map_mean_of_outlinks_to_node(grid, var_name, out=None):
    """Map the mean of links leaving a node to the node.

    map_outlink_average_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the average of
    the outlinks and returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_outlinks_to_node(rmg, 'z')
    array([ 4. ,  5. ,  6. ,  1.5,  7.5,  8.5,  9.5,  3.5,  7. ,  7.5,  8. ,
            0. ])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = 0.5 * (values_at_links[north] + values_at_links[east])

    return out


def map_max_of_outlinks_to_node(grid, var_name, out=None):
    """Map the max of links leaving a node to the node.

    map_max_outlink_value_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the maximum value at the
    the outlinks and returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_max_of_outlinks_to_node(rmg, 'z')
    array([  8.,   9.,  10.,   3.,  11.,  12.,  13.,   7.,  14.,  15.,  16.,
             0.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.maximum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_min_of_outlinks_to_node(grid, var_name, out=None):
    """Map the min of links leaving a node to the node.

    map_min_outlink_value_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the minimum value at the
    the outlinks and returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_min_of_outlinks_to_node(rmg, 'z')
    array([ 0.,  1.,  2.,  0.,  4.,  5.,  6.,  0.,  0.,  0.,  0.,  0.])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.minimum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_mean_of_links_to_node(grid, var_name, out=None):
    """Map the mean of links touching a node to the node.

    map_average_all_links_to_node takes a field *at the links* and finds the
    average of all ~existing~ link neighbor values for each node in the grid.
    it returns a field at the nodes with the same var_name
    as the link field.

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
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
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_links_to_node(rmg, 'z')
    array([  4.        ,   6.        ,   7.        ,   6.5       ,
             5.        ,   7.25      ,   8.25      ,   7.66666667,
             9.        ,  11.33333333,  12.33333333,  11.5       ])
    """
    if out is None:
        out = grid.empty(centering='node')

    values_at_links = grid.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)

    north, east = links._node_out_link_ids(grid.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    south, west = links._node_in_link_ids(grid.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)

    number_of_links = links.number_of_links_per_node(grid.shape)
    number_of_links = number_of_links.reshape(number_of_links.size)
    number_of_links.astype(float, copy=False)
    out[:] = (values_at_links[north] + values_at_links[east] +
              values_at_links[south] + values_at_links[west]) / number_of_links

    return out
