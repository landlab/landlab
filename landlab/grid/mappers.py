#! /usr/bin/env python
"""Map values from one grid element to another.

Each link has a *tail* and *head* node. The *tail* nodes are located at the
start of a link, while the head nodes are located at end of a link.

Below, the numbering scheme for links in `RasterModelGrid` is illustrated
with an example of a four-row by five column grid (4x5). In this example,
each * (or X) is a node, the lines represent links, and the ^ and > symbols
indicate the direction and *head* of each link. Link heads in the
`RasterModelGrid` always point in the cardinal directions North (N) or East
(E).::

    *--27-->*--28-->*--29-->*--30-->*
    ^       ^       ^       ^       ^
   10      11      12      13      14
    |       |       |       |       |
    *--23-->*--24-->*--25-->*--26-->*
    ^       ^       ^       ^       ^
    5       6       7       8       9
    |       |       |       |       |
    *--19-->*--20-->X--21-->*--22-->*
    ^       ^       ^       ^       ^
    0       1       2       3       4
    |       |       |       |       |
    *--15-->*--16-->*--17-->*--18-->*
  
For example, node 'X' has four link-neighbors. From south and going clockwise,
these neighbors are [2, 20, 7, 21]. Both link 2 and link 20 have node 'X' as
their 'head' node, while links 7 and 21 have node 'X' as their tail node. 
"""
from __future__ import division

import numpy as np
from landlab.grid.structured_quad import links


def map_link_head_node_to_link(mg, var_name, out=None):
    """Map values from a link head nodes to links.

    Iterate over a grid and identify the node at the *head*. For each link,
    the value of *var_name* at the *head* node is mapped to the corresponding
    link. 
    
    In a RasterModelGrid, each one node has two adjacent "link heads". This
    means each node value is mapped to two corresponding links. 

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_link_head_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', np.arange(12.))
    >>> map_link_head_node_to_link(rmg, 'z')
    array([  4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,   1.,   2.,   3.,
             5.,   6.,   7.,   9.,  10.,  11.])

    >>> values_at_links = rmg.empty(centering='link')
    >>> rtn = map_link_head_node_to_link(rmg, 'z', out=values_at_links)
    >>> values_at_links
    array([  4.,   5.,   6.,   7.,   8.,   9.,  10.,  11.,   1.,   2.,   3.,
             5.,   6.,   7.,   9.,  10.,  11.])
    >>> rtn is values_at_links
    True
    """
    values_at_nodes = mg.at_node[var_name]
    if out is None:
        out = mg.empty(centering='link')
    out[:] = values_at_nodes[mg.node_at_link_head]

    return out


def map_link_tail_node_to_link(mg, var_name, out=None):
    """Map values from a link tail nodes to links.

    map_values_from_link_tail_node_to_link iterates across the grid and
    identifies the node at the "tail", or the "from" node for each link. For
    each link, the value of 'var_name' at the "from" node is mapped to the 
    corresponding link. 
    
    In a RasterModelGrid, each one node has two adjacent "link tails". This means
    each node value is mapped to two corresponding links. 

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_link_tail_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', np.arange(12.))
    >>> map_link_tail_node_to_link(rmg, 'z')
    array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   0.,   1.,   2.,
             4.,   5.,   6.,   8.,   9.,  10.])

    >>> values_at_links = rmg.empty(centering='link')
    >>> rtn = map_link_tail_node_to_link(rmg, 'z', out=values_at_links)
    >>> values_at_links
    array([  0.,   1.,   2.,   3.,   4.,   5.,   6.,   7.,   0.,   1.,   2.,
             4.,   5.,   6.,   8.,   9.,  10.])
    >>> rtn is values_at_links
    True
    """
    if out is None:
        out = mg.empty(centering='link')

    values_at_nodes = mg.at_node[var_name]
    out[:] = values_at_nodes[mg.node_at_link_tail]

    return out


def map_min_of_link_nodes_to_link(mg, var_name, out=None):
    """Map the minimum of a link's nodes to the link.

    map_link_end_node_min_value_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function evaluates the value of 'var_name' at both the "to" and "from" node.
    The minimum value of the two node values is then mapped to the link.

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_min_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', [[0, 1, 2, 3], [7, 6, 5, 4], [8, 9, 10, 11]])
    >>> map_min_of_link_nodes_to_link(rmg, 'z')
    array([  0.,   1.,   2.,   3.,   7.,   6.,   5.,   4.,   0.,   1.,   2.,
             6.,   5.,   4.,   8.,   9.,  10.])

    >>> values_at_links = rmg.empty(centering='link')
    >>> rtn = map_min_of_link_nodes_to_link(rmg, 'z', out=values_at_links)
    >>> values_at_links
    array([  0.,   1.,   2.,   3.,   7.,   6.,   5.,   4.,   0.,   1.,   2.,
             6.,   5.,   4.,   8.,   9.,  10.])
    >>> rtn is values_at_links
    True
    """
    if out is None:
        out = mg.empty(centering='link')

    values_at_nodes = mg.at_node[var_name]
    np.minimum(values_at_nodes[mg.node_at_link_head],
               values_at_nodes[mg.node_at_link_tail],
               out=out)

    return out


def map_max_of_link_nodes_to_link(mg, var_name, out=None):
    """Map the maximum of a link's nodes to the link.

    map_link_end_node_max_value_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function evaluates the value of 'var_name' at both the "to" and "from" node.
    The maximum value of the two node values is then mapped to the link.

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_max_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', [[0, 1, 2, 3], [7, 6, 5, 4], [8, 9, 10, 11]])
    >>> map_max_of_link_nodes_to_link(rmg, 'z')
    array([  7.,   6.,   5.,   4.,   8.,   9.,  10.,  11.,   1.,   2.,   3.,
             7.,   6.,   5.,   9.,  10.,  11.])

    >>> values_at_links = rmg.empty(centering='link')
    >>> rtn = map_max_of_link_nodes_to_link(rmg, 'z', out=values_at_links)
    >>> values_at_links
    array([  7.,   6.,   5.,   4.,   8.,   9.,  10.,  11.,   1.,   2.,   3.,
             7.,   6.,   5.,   9.,  10.,  11.])
    >>> rtn is values_at_links
    True
    """
    if out is None:
        out = mg.empty(centering='link')

    values_at_nodes = mg.at_node[var_name]
    np.maximum(values_at_nodes[mg.node_at_link_head],
               values_at_nodes[mg.node_at_link_tail],
               out=out)

    return out


def map_mean_of_link_nodes_to_link(mg, var_name, out=None):
    """Map the mean of a link's nodes to the link.

    map_values_from_link_end_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.  
    This function takes the sum of the two values of 'var_name' at both the "to"
    and "from" node. The average value of the two node values of 'var_name'
    is then mapped to the link.

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', np.arange(12.))
    >>> map_mean_of_link_nodes_to_link(rmg, 'z')
    array([  2. ,   3. ,   4. ,   5. ,   6. ,   7. ,   8. ,   9. ,   0.5,
             1.5,   2.5,   4.5,   5.5,   6.5,   8.5,   9.5,  10.5])

    >>> values_at_links = rmg.empty(centering='link')
    >>> rtn = map_mean_of_link_nodes_to_link(rmg, 'z', out=values_at_links)
    >>> values_at_links
    array([  2. ,   3. ,   4. ,   5. ,   6. ,   7. ,   8. ,   9. ,   0.5,
             1.5,   2.5,   4.5,   5.5,   6.5,   8.5,   9.5,  10.5])
    >>> rtn is values_at_links
    True
    """
    if out is None:
        out = mg.empty(centering='link')

    values_at_nodes = mg.at_node[var_name]
    out[:] = 0.5 * (values_at_nodes[mg.node_at_link_head] +
                    values_at_nodes[mg.node_at_link_tail])

    return out


def map_node_to_cell(mg, var_name, out=None):
    """Map values for nodes to cells.

    map_values_from_cell_node_to_cell iterates across the grid and
    identifies the all node values of 'var_name'.  
    
    This function takes node values of 'var_name' and mapes that value to the
    corresponding cell area for each node.

    Parameters
    ----------
    mg : ModelGrid
        A landlab ModelGrid.
    var_name : str
        Name of variable field defined at nodes.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at cells.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_node_to_cell
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('node', 'z', np.arange(12.))
    >>> map_node_to_cell(rmg, 'z')
    array([ 5., 6.])

    >>> values_at_cells = rmg.empty(centering='cell')
    >>> rtn = map_node_to_cell(rmg, 'z', out=values_at_cells)
    >>> values_at_cells
    array([ 5., 6.])
    >>> rtn is values_at_cells
    True
    """
    if out is None:
        out = mg.empty(centering='cell')

    values_at_nodes = mg.at_node[var_name]
    out[:] = values_at_nodes[mg.node_at_cell]

    return out


def map_sum_of_inlinks_to_node(mg, var_name, out=None):
    """Map the sum of links entering a node to the node.

    map_inlink_sums_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it sums the inlinks and returns
    a field at the nodes with the same var_name as the link field. 
    
    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_sum_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_sum_of_inlinks_to_node(rmg, 'z')
    array([  0.,   8.,   9.,  10.,   0.,  12.,  14.,  16.,   4.,  19.,  21.,
            23.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)

    south, west = links._node_in_link_ids(mg.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = values_at_links[south] + values_at_links[west]

    return out


def map_mean_of_inlinks_to_node(mg, var_name, out=None):
    """Map the mean of links entering a node to the node.

    map_inlink_average_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the average of
    the inlinks and returns a field at the nodes with the same var_name
    as the link field. 
    
    This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_mean_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_inlinks_to_node(rmg, 'z')
    array([  0. ,   4. ,   4.5,   5. ,   0. ,   6. ,   7. ,   8. ,   2. ,
             9.5,  10.5,  11.5])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(mg.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = 0.5 * (values_at_links[south] + values_at_links[west])

    return out

def map_max_of_inlinks_to_node(mg, var_name, out=None):
    """Map the maximum of links entering a node to the node.

    map_max_inlink_value_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the maximum value at the
    the inlinks and returns a field at the nodes with the same var_name
    as the link field. 
    
    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_max_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_max_of_inlinks_to_node(rmg, 'z')
    array([  0.,   8.,   9.,  10.,   0.,  11.,  12.,  13.,   4.,  14.,  15.,
            16.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(mg.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.maximum(values_at_links[south], values_at_links[west])

    return out


def map_min_of_inlinks_to_node(mg, var_name, out=None):
    """Map the minimum of links entering a node to the node.

    map_min_inlink_value_to_node takes a field *at the links* and finds the
    inlink values for each node in the grid. it finds the minimum value at the
    the inlinks and returns a field at the nodes with the same var_name
    as the link field. 
    
    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_min_of_inlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_min_of_inlinks_to_node(rmg, 'z')
    array([ 0.,  0.,  0.,  0.,  0.,  1.,  2.,  3.,  0.,  5.,  6.,  7.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    south, west = links._node_in_link_ids(mg.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)
    out[:] = np.minimum(values_at_links[south], values_at_links[west])

    return out


def map_sum_of_outlinks_to_node(mg, var_name, out=None):
    """Map the sum of links leaving a node to the node.

    map_outlink_sums_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it sums the outlinks and returns
    a field at the nodes with the same var_name as the link field. 
    
    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_sum_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_sum_of_outlinks_to_node(rmg, 'z')
    array([  8.,  10.,  12.,   3.,  15.,  17.,  19.,   7.,  14.,  15.,  16.,
             0.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(mg.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = values_at_links[north] + values_at_links[east]

    return out


def map_mean_of_outlinks_to_node(mg, var_name, out=None):
    """Map the mean of links leaving a node to the node.

    map_outlink_average_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the average of
    the outlinks and returns a field at the nodes with the same var_name
    as the link field. 

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_mean_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_outlinks_to_node(rmg, 'z')
    array([ 4. ,  5. ,  6. ,  1.5,  7.5,  8.5,  9.5,  3.5,  7. ,  7.5,  8. ,
            0. ])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(mg.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    out[:] = 0.5 * (values_at_links[north] + values_at_links[east])

    return out


def map_max_of_outlinks_to_node(mg, var_name, out=None):
    """Map the max of links leaving a node to the node.

    map_max_outlink_value_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the maximum value at the
    the outlinks and returns a field at the nodes with the same var_name
    as the link field. 

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_max_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_max_of_outlinks_to_node(rmg, 'z')
    array([  8.,   9.,  10.,   3.,  11.,  12.,  13.,   7.,  14.,  15.,  16.,
             0.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(mg.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.maximum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_min_of_outlinks_to_node(mg, var_name, out=None):
    """Map the min of links leaving a node to the node.

    map_min_outlink_value_to_node takes a field *at the links* and finds the
    outlink values for each node in the grid. it finds the minimum value at the
    the outlinks and returns a field at the nodes with the same var_name
    as the link field. 

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_min_of_outlinks_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_min_of_outlinks_to_node(rmg, 'z')
    array([ 0.,  1.,  2.,  0.,  4.,  5.,  6.,  0.,  0.,  0.,  0.,  0.])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)
    north, east = links._node_out_link_ids(mg.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    np.minimum(values_at_links[north], values_at_links[east], out=out)

    return out


def map_mean_of_links_to_node(mg, var_name, out=None):
    """Map the mean of links touching a node to the node.

    map_average_all_links_to_node takes a field *at the links* and finds the
    average of all ~existing~ link neighbor values for each node in the grid. 
    it returns a field at the nodes with the same var_name
    as the link field. 

    .. note::

        This considers all inactive links to have a value of 0.

    Parameters
    ----------
    mg : ModelGrid
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
    >>> from landlab.grid.mappers import map_mean_of_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field('link', 'z', np.arange(17.))
    >>> map_mean_of_links_to_node(rmg, 'z')
    array([  4.        ,   6.        ,   7.        ,   6.5       ,
             5.        ,   7.25      ,   8.25      ,   7.66666667,
             9.        ,  11.33333333,  12.33333333,  11.5       ])
    """
    if out is None:
        out = mg.empty(centering='node')

    values_at_links = mg.at_link[var_name]
    values_at_links = np.append(values_at_links, 0)

    north, east = links._node_out_link_ids(mg.shape)
    north, east = north.reshape(north.size), east.reshape(east.size)
    south, west = links._node_in_link_ids(mg.shape)
    south, west = south.reshape(south.size), west.reshape(west.size)

    number_of_links = links.number_of_links_per_node(mg.shape)
    number_of_links = number_of_links.reshape(number_of_links.size)
    number_of_links.astype(float, copy=False)
    out[:] = (values_at_links[north] + values_at_links[east] +
              values_at_links[south] + values_at_links[west]) / number_of_links

    return out
