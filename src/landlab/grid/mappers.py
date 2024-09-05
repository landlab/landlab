#! /usr/bin/env python
"""Map values from one grid element to another.

Grid mapping functions
+++++++++++++++++++++++

.. autosummary::

    ~landlab.grid.mappers.map_link_head_node_to_link
    ~landlab.grid.mappers.map_link_tail_node_to_link
    ~landlab.grid.mappers.map_min_of_link_nodes_to_link
    ~landlab.grid.mappers.map_max_of_link_nodes_to_link
    ~landlab.grid.mappers.map_mean_of_link_nodes_to_link
    ~landlab.grid.mappers.map_value_at_min_node_to_link
    ~landlab.grid.mappers.map_value_at_max_node_to_link
    ~landlab.grid.mappers.map_node_to_cell
    ~landlab.grid.mappers.map_min_of_node_links_to_node
    ~landlab.grid.mappers.map_max_of_node_links_to_node
    ~landlab.grid.mappers.map_upwind_node_link_max_to_node
    ~landlab.grid.mappers.map_downwind_node_link_max_to_node
    ~landlab.grid.mappers.map_upwind_node_link_mean_to_node
    ~landlab.grid.mappers.map_downwind_node_link_mean_to_node
    ~landlab.grid.mappers.map_value_at_upwind_node_link_max_to_node
    ~landlab.grid.mappers.map_value_at_downwind_node_link_max_to_node
    ~landlab.grid.mappers.map_link_vector_components_to_node
    ~landlab.grid.mappers.map_node_to_link_linear_upwind
    ~landlab.grid.mappers.map_node_to_link_lax_wendroff

Each link has a *tail* and *head* node. The *tail* nodes are located at the
start of a link, while the head nodes are located at end of a link.

Below, the numbering scheme for links in :class:`~.RasterModelGrid` is illustrated
with an example of a four-row by five column grid (4x5). In this example,
each ``*`` (or ``X``) is a node, the lines represent links, and the ``^`` and ``>`` symbols
indicate the direction and *head* of each link. Link heads in the
:class:`~.RasterModelGrid` always point in the cardinal directions North (N) or East
(E).::

    *--27-->*--28-->*--29-->*--30-->*
    ^       ^       ^       ^       ^
   22      23      24      25      26
    |       |       |       |       |
    *--18-->*--19-->*--20-->*--21-->*
    ^       ^       ^       ^       ^
    13      14      15      16     17
    |       |       |       |       |
    *---9-->*--10-->X--11-->*--12-->*
    ^       ^       ^       ^       ^
    4       5       6       7       8
    |       |       |       |       |
    *--0--->*---1-->*--2--->*---3-->*

For example, node ``X`` has four link-neighbors. From south and going clockwise,
these neighbors are ``[6, 10, 15, 11]``. Both link 6 and link 10 have node ``X`` as
their *head* node, while links 15 and 11 have node ``X`` as their *tail* node.
"""

import numpy as np


def cartesian_to_polar(x, y):
    """Return 2d polar coordinates (r, theta) equivalent to given cartesian
    coordinates (x, y).

    Examples
    --------
    >>> r, theta = cartesian_to_polar(1.0, 1.0)
    >>> int(r * 1000)
    1414
    >>> int(theta * 1000)
    785
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x)  # r, theta


def map_link_head_node_to_link(grid, var_name, out=None):
    """Map values from a link head nodes to links.

    Iterate over a grid and identify the node at the *head*. For each link,
    the value of *var_name* at the *head* node is mapped to the corresponding
    link.

    In a RasterModelGrid, each one node has two adjacent "link heads". This
    means each node value is mapped to two corresponding links.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
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
    >>> from landlab.grid.mappers import map_link_head_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["z"] = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    >>> map_link_head_node_to_link(rmg, "z")
    array([  1.,   2.,   3.,   4.,   5.,   6.,   7.,   5.,   6.,   7.,   8.,
          9.,  10.,  11.,   9.,  10.,  11.])

    >>> values_at_links = rmg.empty(at="link")
    >>> rtn = map_link_head_node_to_link(rmg, "z", out=values_at_links)
    >>> values_at_links
    array([  1.,   2.,   3.,   4.,   5.,   6.,   7.,   5.,   6.,   7.,   8.,
          9.,  10.,  11.,   9.,  10.,  11.])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    if out is None:
        out = grid.empty(at="link")
    out[:] = var_name[grid.node_at_link_head]

    return out


def map_link_tail_node_to_link(grid, var_name, out=None):
    """Map values from a link tail nodes to links.

    map_link_tail_node_to_link iterates across the grid and
    identifies the node at the "tail", or the "from" node for each link. For
    each link, the value of 'var_name' at the "from" node is mapped to the
    corresponding link.

    In a RasterModelGrid, each one node has two adjacent "link tails". This
    means each node value is mapped to two corresponding links.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
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
    >>> from landlab.grid.mappers import map_link_tail_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["z"] = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    >>> map_link_tail_node_to_link(rmg, "z")
    array([  0.,   1.,   2.,   0.,   1.,   2.,   3.,   4.,   5.,   6.,   4.,
             5.,   6.,   7.,   8.,   9.,  10.])

    >>> values_at_links = rmg.empty(at="link")
    >>> rtn = map_link_tail_node_to_link(rmg, "z", out=values_at_links)
    >>> values_at_links
    array([  0.,   1.,   2.,   0.,   1.,   2.,   3.,   4.,   5.,   6.,   4.,
             5.,   6.,   7.,   8.,   9.,  10.])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    out[:] = var_name[grid.node_at_link_tail]

    return out


def map_min_of_link_nodes_to_link(grid, var_name, out=None):
    """Map the minimum of a link's nodes to the link.

    map_min_of_link_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function evaluates the value of 'var_name' at both the "to" and
    "from" node. The minimum value of the two node values is then mapped to
    the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
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
    >>> from landlab.grid.mappers import map_min_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field(
    ...     "z",
    ...     [
    ...         [0, 1, 2, 3],
    ...         [7, 6, 5, 4],
    ...         [8, 9, 10, 11],
    ...     ],
    ...     at="node",
    ... )
    >>> map_min_of_link_nodes_to_link(rmg, "z")
    array([  0.,   1.,   2.,   0.,   1.,   2.,   3.,   6.,   5.,   4.,   7.,
             6.,   5.,   4.,   8.,   9.,  10.])

    >>> values_at_links = rmg.empty(at="link")
    >>> rtn = map_min_of_link_nodes_to_link(rmg, "z", out=values_at_links)
    >>> values_at_links
    array([  0.,   1.,   2.,   0.,   1.,   2.,   3.,   6.,   5.,   4.,   7.,
             6.,   5.,   4.,   8.,   9.,  10.])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    np.minimum(
        var_name[grid.node_at_link_head], var_name[grid.node_at_link_tail], out=out
    )

    return out


def map_max_of_link_nodes_to_link(grid, var_name, out=None):
    """Map the maximum of a link's nodes to the link.

    map_max_of_link_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function evaluates the value of 'var_name' at both the "to" and
    "from" node. The maximum value of the two node values is then mapped to
    the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
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
    >>> from landlab.grid.mappers import map_max_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field(
    ...     "z",
    ...     [
    ...         [0, 1, 2, 3],
    ...         [7, 6, 5, 4],
    ...         [8, 9, 10, 11],
    ...     ],
    ...     at="node",
    ... )
    >>> map_max_of_link_nodes_to_link(rmg, "z")
    array([  1.,   2.,   3.,   7.,   6.,   5.,   4.,   7.,   6.,   5.,   8.,
             9.,  10.,  11.,   9.,  10.,  11.])

    >>> values_at_links = rmg.empty(at="link")
    >>> rtn = map_max_of_link_nodes_to_link(rmg, "z", out=values_at_links)
    >>> values_at_links
    array([  1.,   2.,   3.,   7.,   6.,   5.,   4.,   7.,   6.,   5.,   8.,
             9.,  10.,  11.,   9.,  10.,  11.])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    np.maximum(
        var_name[grid.node_at_link_head], var_name[grid.node_at_link_tail], out=out
    )

    return out


def map_mean_of_link_nodes_to_link(grid, var_name, out=None):
    """Map the mean of a link's nodes to the link.

    map_mean_of_link_nodes_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function takes the sum of the two values of 'var_name' at both the
    "to" and "from" node. The average value of the two node values of
    'var_name' is then mapped to the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
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
    >>> from landlab.grid.mappers import map_mean_of_link_nodes_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["z"] = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])
    >>> map_mean_of_link_nodes_to_link(rmg, "z")
    array([  0.5,   1.5,   2.5,   2. ,   3. ,   4. ,   5. ,   4.5,   5.5,
             6.5,   6. ,   7. ,   8. ,   9. ,   8.5,   9.5,  10.5])

    >>> values_at_links = rmg.empty(at="link")
    >>> rtn = map_mean_of_link_nodes_to_link(rmg, "z", out=values_at_links)
    >>> values_at_links
    array([  0.5,   1.5,   2.5,   2. ,   3. ,   4. ,   5. ,   4.5,   5.5,
             6.5,   6. ,   7. ,   8. ,   9. ,   8.5,   9.5,  10.5])
    >>> rtn is values_at_links
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    out[:] = 0.5 * (var_name[grid.node_at_link_head] + var_name[grid.node_at_link_tail])

    return out


def map_value_at_min_node_to_link(grid, control_name, value_name, out=None):
    """Map the the value found in one node array to a link, based on the
    minimum value found in a second node field or array.

    map_value_at_min_node_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function evaluates the value of 'control_name' at both the "to" and
    "from" node. The value of 'value_name' at the node with the minimum value
    of the two values of 'control_name' is then mapped to the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    control_name : array or field name
        Name of field defined at nodes or a node array that dictates which end
        of the link to draw values from.
    value_name : array or field name
        Name of field defined at nodes or  node array from which values are
        drawn, based on control_name.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_value_at_min_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field(
    ...     "z",
    ...     [
    ...         [0, 1, 2, 3],
    ...         [7, 6, 5, 4],
    ...         [8, 9, 10, 11],
    ...     ],
    ...     at="node",
    ... )
    >>> _ = rmg.add_field(
    ...     "vals_to_map",
    ...     [
    ...         [0, 10, 20, 30],
    ...         [70, 60, 50, 40],
    ...         [80, 90, 100, 110],
    ...     ],
    ...     at="node",
    ... )
    >>> map_value_at_min_node_to_link(rmg, "z", "vals_to_map")
    array([   0.,   10.,   20.,    0.,   10.,   20.,   30.,   60.,   50.,
             40.,   70.,   60.,   50.,   40.,   80.,   90.,  100.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(control_name) is str:
        control_name = grid.at_node[control_name]
    if type(value_name) is str:
        value_name = grid.at_node[value_name]
    head_control = control_name[grid.node_at_link_head]
    tail_control = control_name[grid.node_at_link_tail]
    head_vals = value_name[grid.node_at_link_head]
    tail_vals = value_name[grid.node_at_link_tail]

    out[:] = np.where(tail_control < head_control, tail_vals, head_vals)
    return out


def map_value_at_max_node_to_link(grid, control_name, value_name, out=None):
    """Map the the value found in one node array to a link, based on the
    maximum value found in a second node field or array.

    map_value_at_max_node_to_link iterates across the grid and
    identifies the node values at both the "head" and "tail" of a given link.
    This function evaluates the value of 'control_name' at both the "to" and
    "from" node. The value of 'value_name' at the node with the maximum value
    of the two values of 'control_name' is then mapped to the link.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    control_name : array or field name
        Name of field defined at nodes or a node array that dictates which end
        of the link to draw values from.
    value_name : array or field name
        Name of field defined at nodes or  node array from which values are
        drawn, based on control_name.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at links.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_value_at_max_node_to_link
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> _ = rmg.add_field(
    ...     "z",
    ...     [
    ...         [0, 1, 2, 3],
    ...         [7, 6, 5, 4],
    ...         [8, 9, 10, 11],
    ...     ],
    ...     at="node",
    ... )
    >>> _ = rmg.add_field(
    ...     "vals_to_map",
    ...     [
    ...         [0, 10, 20, 30],
    ...         [70, 60, 50, 40],
    ...         [80, 90, 100, 110],
    ...     ],
    ...     at="node",
    ... )
    >>> map_value_at_max_node_to_link(rmg, "z", "vals_to_map")
    array([  10.,   20.,   30.,   70.,   60.,   50.,   40.,   70.,   60.,
             50.,   80.,   90.,  100.,  110.,   90.,  100.,  110.])

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="link")

    if type(control_name) is str:
        control_name = grid.at_node[control_name]
    if type(value_name) is str:
        value_name = grid.at_node[value_name]
    head_control = control_name[grid.node_at_link_head]
    tail_control = control_name[grid.node_at_link_tail]
    head_vals = value_name[grid.node_at_link_head]
    tail_vals = value_name[grid.node_at_link_tail]

    out[:] = np.where(tail_control > head_control, tail_vals, head_vals)
    return out


def map_node_to_cell(grid, var_name, out=None):
    """Map values for nodes to cells.

    map_node_to_cell iterates across the grid and
    identifies the all node values of 'var_name'.

    This function takes node values of 'var_name' and mapes that value to the
    corresponding cell area for each node.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at nodes.
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
    >>> _ = rmg.add_field("z", np.arange(12.0), at="node")
    >>> map_node_to_cell(rmg, "z")
    array([5.,  6.])

    >>> values_at_cells = rmg.empty(at="cell")
    >>> rtn = map_node_to_cell(rmg, "z", out=values_at_cells)
    >>> values_at_cells
    array([5.,  6.])
    >>> rtn is values_at_cells
    True

    :meta landlab: info-cell, info-node, map
    """
    if out is None:
        out = grid.empty(at="cell")

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    out[:] = var_name[grid.node_at_cell]

    return out


def map_min_of_node_links_to_node(grid, var_name, out=None):
    """Map the minimum value of a nodes' links to the node.

    map_min_of_node_links_to_node iterates across the grid and
    identifies the link values at each link connected to  a node.
    This function finds the minimum value of 'var_name' of each set
    of links, and then maps this value to the node. Note no attempt is made
    to honor the directionality of the links.

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
    >>> from landlab.grid.mappers import map_min_of_node_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = np.arange(rmg.number_of_links)
    >>> map_min_of_node_links_to_node(rmg, "grad")
    array([  0.,   0.,   1.,   2.,
             3.,   4.,   5.,   6.,
            10.,  11.,  12.,  13.])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_min_of_node_links_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes
    array([  0.,   0.,   1.,   2.,
             3.,   4.,   5.,   6.,
            10.,  11.,  12.,  13.])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    values_at_linksX = np.empty(grid.number_of_links + 1, dtype=float)
    values_at_linksX[-1] = np.finfo(dtype=float).max
    if type(var_name) is str:
        values_at_linksX[:-1] = grid.at_link[var_name]
    else:
        values_at_linksX[:-1] = var_name
    np.amin(values_at_linksX[grid.links_at_node], axis=1, out=out)

    return out


def map_max_of_node_links_to_node(grid, var_name, out=None):
    """Map the maximum value of a nodes' links to the node.

    map_max_of_node_links_to_node iterates across the grid and
    identifies the link values at each link connected to  a node.
    This function finds the maximum value of 'var_name' of each set
    of links, and then maps this value to the node. Note no attempt is made
    to honor the directionality of the links.

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
    >>> from landlab.grid.mappers import map_max_of_node_links_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = np.arange(rmg.number_of_links)
    >>> map_max_of_node_links_to_node(rmg, "grad")
    array([  3.,   4.,   5.,   6.,
            10.,  11.,  12.,  13.,
            14.,  15.,  16.,  16.])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_max_of_node_links_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes
    array([  3.,   4.,   5.,   6.,
            10.,  11.,  12.,  13.,
            14.,  15.,  16.,  16.])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    values_at_linksX = np.empty(grid.number_of_links + 1, dtype=float)
    values_at_linksX[-1] = np.finfo(dtype=float).min
    if type(var_name) is str:
        values_at_linksX[:-1] = grid.at_link[var_name]
    else:
        values_at_linksX[:-1] = var_name
    np.amax(values_at_linksX[grid.links_at_node], axis=1, out=out)

    return out


def map_upwind_node_link_max_to_node(grid, var_name, out=None):
    """Map the largest magnitude of the links bringing flux into the node to
    the node.

    map_upwind_node_link_max_to_node iterates across the grid and identifies
    the link values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links bringing flux into the
    node, then maps the maximum magnitude of 'var_name' found on these links
    onto the node. If no upwind link is found, the value will be recorded as
    zero.

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
    >>> from landlab.grid.mappers import map_upwind_node_link_max_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.1,
    ...     -1.2,
    ...     -1.3,
    ...     1.4,
    ...     1.5,
    ...     1.6,
    ...     -1.7,
    ...     -1.8,
    ...     -1.9,
    ...     2.0,
    ...     2.1,
    ...     2.2,
    ...     -2.3,
    ...     2.4,
    ...     2.5,
    ...     2.6,
    ...     -2.7,
    ... ]
    >>> map_upwind_node_link_max_to_node(rmg, "grad").reshape((3, 4))
    array([[1.4,  1.5,  1.6,  1.3],
           [2.1,  2.2,  2. ,  2.4],
           [2.5,  2.6,  2.3,  2.7]])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_upwind_node_link_max_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes.reshape((3, 4))
    array([[1.4,  1.5,  1.6,  1.3],
           [2.1,  2.2,  2. ,  2.4],
           [2.5,  2.6,  2.3,  2.7]])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    if type(var_name) is str:
        var_name = grid.at_link[var_name]
    values_at_links = var_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    np.amax(-values_at_links, axis=1, out=out)

    return out


def map_downwind_node_link_max_to_node(grid, var_name, out=None):
    """Map the largest magnitude of the links carrying flux from the node to
    the node.

    map_downwind_node_link_max_to_node iterates across the grid and identifies
    the link values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links carrying flux out of the
    node, then maps the maximum magnitude of 'var_name' found on these links
    onto the node. If no downwind link is found, the value will be recorded as
    zero.

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
    >>> from landlab.grid.mappers import map_downwind_node_link_max_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ... ]
    >>> map_downwind_node_link_max_to_node(rmg, "grad")
    array([1.,  2.,  1.,  0.,
           1.,  2.,  1.,  0.,
           1.,  2.,  1.,  0.])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_downwind_node_link_max_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes
    array([1.,  2.,  1.,  0.,
           1.,  2.,  1.,  0.,
           1.,  2.,  1.,  0.])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    if type(var_name) is str:
        var_name = grid.at_link[var_name]
    values_at_links = var_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    steepest_links_at_node = np.amax(values_at_links, axis=1)
    np.fabs(steepest_links_at_node, out=out)

    return out


def map_upwind_node_link_mean_to_node(grid, var_name, out=None):
    """Map the mean magnitude of the links bringing flux into the node to the
    node.

    map_upwind_node_link_mean_to_node iterates across the grid and identifies
    the link values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links bringing flux into the
    node, then maps the mean magnitude of 'var_name' found on these links
    onto the node. Links with zero values are not included in the means,
    and zeros are returned if no upwind links are found.

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
    >>> from landlab.grid.mappers import map_upwind_node_link_mean_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     -2.0,
    ...     -3.0,
    ...     -4.0,
    ...     -5.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     -1.0,
    ...     -2.0,
    ...     -3.0,
    ...     -4.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ... ]
    >>> map_upwind_node_link_mean_to_node(rmg, "grad")
    array([0. ,  1. ,  2. ,  1. ,
           2. ,  2. ,  3. ,  3. ,
           1. ,  1.5,  2.5,  2.5])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_upwind_node_link_mean_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes
    array([0. ,  1. ,  2. ,  1. ,
           2. ,  2. ,  3. ,  3. ,
           1. ,  1.5,  2.5,  2.5])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")
    out[:] = 0.0

    if type(var_name) is str:
        var_name = grid.at_link[var_name]
    values_at_links = var_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    vals_in_positive = -values_at_links
    vals_above_zero = vals_in_positive > 0.0
    total_vals = np.sum(vals_in_positive * vals_above_zero, axis=1)
    link_count = np.sum(vals_above_zero, axis=1)
    np.divide(total_vals, link_count, out=out, where=link_count != 0)
    out[link_count == 0] = 0.0

    return out


def map_downwind_node_link_mean_to_node(grid, var_name, out=None):
    """Map the mean magnitude of the links carrying flux out of the node to the
    node.

    map_downwind_node_link_mean_to_node iterates across the grid and identifies
    the link values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links carrying flux out of the
    node, then maps the mean magnitude of 'var_name' found on these links
    onto the node. Links with zero values are not included in the means,
    and zeros are returned if no upwind links are found.

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
    >>> from landlab.grid.mappers import map_downwind_node_link_mean_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     -2.0,
    ...     -3.0,
    ...     -4.0,
    ...     -5.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     -1.0,
    ...     -2.0,
    ...     -3.0,
    ...     -4.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ... ]
    >>> map_downwind_node_link_mean_to_node(rmg, "grad")
    array([1.5,  2.5,  2.5,  5. ,
           1. ,  2. ,  2. ,  4. ,
           1. ,  2. ,  1. ,  0. ])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_downwind_node_link_mean_to_node(rmg, "grad", out=values_at_nodes)
    >>> values_at_nodes
    array([1.5,  2.5,  2.5,  5. ,
           1. ,  2. ,  2. ,  4. ,
           1. ,  2. ,  1. ,  0. ])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")
    out[:] = 0.0

    if type(var_name) is str:
        var_name = grid.at_link[var_name]
    values_at_links = var_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    vals_in_positive = values_at_links
    vals_above_zero = vals_in_positive > 0.0
    total_vals = np.sum(vals_in_positive * vals_above_zero, axis=1)
    link_count = np.sum(vals_above_zero, axis=1)
    np.divide(total_vals, link_count, out=out, where=link_count != 0)

    return out


def map_value_at_upwind_node_link_max_to_node(grid, control_name, value_name, out=None):
    """Map the the value found in one link array to a node, based on the
    largest magnitude value of links bringing fluxes into the node, found in a
    second node array or field.

    map_upwind_node_link_max_to_node iterates across the grid and identifies
    the link control_values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links bringing flux into the
    node, then identifies the link with the maximum magnitude. The value of the
    second field 'value_name' at these links is then mapped onto the node.
    If no upwind link is found, the value will be recorded as zero.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    control_name : array or field name
        Values defined at nodes that dictate which end of the link to
        draw values from.
    value_name : array or field name
        Values defined at nodes from which values are drawn, based on
        control_name.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_value_at_upwind_node_link_max_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ... ]
    >>> rmg.at_link["vals"] = np.arange(rmg.number_of_links, dtype=float)
    >>> map_value_at_upwind_node_link_max_to_node(rmg, "grad", "vals")
    array([  0.,   0.,   1.,   2.,
             0.,   7.,   8.,   9.,
             0.,  14.,  15.,  16.])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_value_at_upwind_node_link_max_to_node(
    ...     rmg, "grad", "vals", out=values_at_nodes
    ... )
    >>> values_at_nodes
    array([  0.,   0.,   1.,   2.,
             0.,   7.,   8.,   9.,
             0.,  14.,  15.,  16.])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    if type(control_name) is str:
        control_name = grid.at_link[control_name]
    if type(value_name) is str:
        value_name = grid.at_link[value_name]
    values_at_nodes = control_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    which_link = np.argmax(-values_at_nodes, axis=1)
    invalid_links = values_at_nodes >= 0.0
    link_vals_without_invalids = value_name[grid.links_at_node]
    link_vals_without_invalids[invalid_links] = 0.0
    out[:] = link_vals_without_invalids[np.arange(grid.number_of_nodes), which_link]

    return out


def map_value_at_downwind_node_link_max_to_node(
    grid, control_name, value_name, out=None
):
    """Map the the value found in one link array to a node, based on the
    largest magnitude value of links carrying fluxes out of the node, found in
    a second node array or field.

    map_downwind_node_link_max_to_node iterates across the grid and identifies
    the link control_values at each link connected to a node. It then uses the
    link_dirs_at_node data structure to identify links carrying flux out of the
    node, then identifies the link with the maximum magnitude. The value of the
    second field 'value_name' at these links is then mapped onto the node.
    If no downwind link is found, the value will be recorded as zero.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    control_name : array or field name
        Values defined at nodes that dictate which end of the link to
        draw values from.
    value_name : array or field name
        Values defined at nodes from which values are drawn, based on
        control_name.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_value_at_downwind_node_link_max_to_node
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_link["grad"] = [
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     0.0,
    ...     -1.0,
    ...     -2.0,
    ...     -1.0,
    ... ]
    >>> rmg.at_link["vals"] = np.arange(rmg.number_of_links, dtype=float)
    >>> map_value_at_downwind_node_link_max_to_node(rmg, "grad", "vals")
    array([  0.,   1.,   2.,   0.,
             7.,   8.,   9.,   0.,
            14.,  15.,  16.,   0.])

    >>> values_at_nodes = rmg.add_empty("z", at="node")
    >>> rtn = map_value_at_downwind_node_link_max_to_node(
    ...     rmg, "grad", "vals", out=values_at_nodes
    ... )
    >>> values_at_nodes
    array([  0.,   1.,   2.,   0.,
             7.,   8.,   9.,   0.,
            14.,  15.,  16.,   0.])
    >>> rtn is values_at_nodes
    True

    :meta landlab: info-node, info-link, map
    """
    if out is None:
        out = grid.empty(at="node")

    if type(control_name) is str:
        control_name = grid.at_link[control_name]
    if type(value_name) is str:
        value_name = grid.at_link[value_name]
    values_at_nodes = control_name[grid.links_at_node] * grid.link_dirs_at_node
    # this procedure makes incoming links NEGATIVE
    which_link = np.argmax(values_at_nodes, axis=1)
    invalid_links = values_at_nodes <= 0.0
    link_vals_without_invalids = value_name[grid.links_at_node]
    link_vals_without_invalids[invalid_links] = 0.0
    out[:] = link_vals_without_invalids[np.arange(grid.number_of_nodes), which_link]

    return out


def map_mean_of_patch_nodes_to_patch(
    grid, var_name, ignore_closed_nodes=True, out=None
):
    """Map the mean value of nodes around a patch to the patch.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at nodes.
    ignore_closed_nodes : bool
        If True, do not incorporate closed nodes into calc. If all nodes are
        masked at a patch, record zero if out is None or leave the existing
        value if out.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at patches.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_mean_of_patch_nodes_to_patch
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> map_mean_of_patch_nodes_to_patch(rmg, "vals")
    array([4.5, 3.5, 2.5,
           3.5, 2.5, 1.5])

    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> rmg.status_at_node[rmg.node_x > 1.5] = rmg.BC_NODE_IS_CLOSED
    >>> ans = np.zeros(6, dtype=float)
    >>> _ = map_mean_of_patch_nodes_to_patch(rmg, "vals", out=ans)
    >>> ans
    array([4.5, 4. , 0. ,
           3.5, 3. , 0. ])

    :meta landlab: info-patch, info-node, map
    """
    if out is None:
        out = np.zeros(grid.number_of_patches, dtype=float)

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    values_at_nodes = var_name[grid.nodes_at_patch]
    if ignore_closed_nodes:
        values_at_nodes = np.ma.masked_where(
            grid.status_at_node[grid.nodes_at_patch] == grid.BC_NODE_IS_CLOSED,
            values_at_nodes,
            copy=False,
        )
        meanvals = np.mean(values_at_nodes, axis=1)
        if type(meanvals.mask) is not np.bool_:
            gooddata = np.logical_not(meanvals.mask)
            out[gooddata] = meanvals.data[gooddata]
        else:
            if not meanvals.mask:
                out[:] = meanvals.data
    else:
        np.mean(values_at_nodes, axis=1, out=out)

    return out


def map_max_of_patch_nodes_to_patch(grid, var_name, ignore_closed_nodes=True, out=None):
    """Map the maximum value of nodes around a patch to the patch.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at nodes.
    ignore_closed_nodes : bool
        If True, do not incorporate closed nodes into calc. If all nodes are
        masked at a patch, record zero if out is None or leave the existing
        value if out.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at patches.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_max_of_patch_nodes_to_patch
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> map_max_of_patch_nodes_to_patch(rmg, "vals")
    array([5., 4., 3.,
           4., 4., 3.])

    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> rmg.status_at_node[rmg.node_x > 1.5] = rmg.BC_NODE_IS_CLOSED
    >>> ans = np.zeros(6, dtype=float)
    >>> _ = map_max_of_patch_nodes_to_patch(rmg, "vals", out=ans)
    >>> ans
    array([5., 4., 0.,
           4., 4., 0.])

    :meta landlab: info-patch, info-node, map
    """
    if out is None:
        out = np.zeros(grid.number_of_patches, dtype=float)

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    values_at_nodes = var_name[grid.nodes_at_patch]
    if ignore_closed_nodes:
        values_at_nodes = np.ma.masked_where(
            grid.status_at_node[grid.nodes_at_patch] == grid.BC_NODE_IS_CLOSED,
            values_at_nodes,
            copy=False,
        )
        maxvals = values_at_nodes.max(axis=1)
        if type(maxvals.mask) is not np.bool_:
            gooddata = np.logical_not(maxvals.mask)
            out[gooddata] = maxvals.data[gooddata]
        else:
            if not maxvals.mask:
                out[:] = maxvals.data
    else:
        np.amax(values_at_nodes, axis=1, out=out)

    return out


def map_min_of_patch_nodes_to_patch(grid, var_name, ignore_closed_nodes=True, out=None):
    """Map the minimum value of nodes around a patch to the patch.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at nodes.
    ignore_closed_nodes : bool
        If True, do not incorporate closed nodes into calc. If all nodes are
        masked at a patch, record zero if out is None or leave the existing
        value if out.
    out : ndarray, optional
        Buffer to place mapped values into or `None` to create a new array.

    Returns
    -------
    ndarray
        Mapped values at patches.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_min_of_patch_nodes_to_patch
    >>> from landlab import RasterModelGrid

    >>> rmg = RasterModelGrid((3, 4))
    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> map_min_of_patch_nodes_to_patch(rmg, "vals")
    array([4., 3., 2.,
           2., 1., 0.])

    >>> rmg.at_node["vals"] = [
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [5.0, 4.0, 3.0, 2.0],
    ...     [3.0, 2.0, 1.0, 0.0],
    ... ]
    >>> rmg.status_at_node[rmg.node_x > 1.5] = rmg.BC_NODE_IS_CLOSED
    >>> ans = np.zeros(6, dtype=float)
    >>> _ = map_min_of_patch_nodes_to_patch(rmg, "vals", out=ans)
    >>> ans
    array([4., 4., 0.,
           2., 2., 0.])

    :meta landlab: info-patch, info-node, map
    """
    if out is None:
        out = np.zeros(grid.number_of_patches, dtype=float)

    if type(var_name) is str:
        var_name = grid.at_node[var_name]
    values_at_nodes = var_name[grid.nodes_at_patch]
    if ignore_closed_nodes:
        values_at_nodes = np.ma.masked_where(
            grid.status_at_node[grid.nodes_at_patch] == grid.BC_NODE_IS_CLOSED,
            values_at_nodes,
            copy=False,
        )
        minvals = values_at_nodes.min(axis=1)
        if type(minvals.mask) is not np.bool_:
            gooddata = np.logical_not(minvals.mask)
            out[gooddata] = minvals.data[gooddata]
        else:
            if not minvals.mask:
                out[:] = minvals.data
    else:
        np.amin(values_at_nodes, axis=1, out=out)

    return out


def map_link_vector_sum_to_patch(grid, var_name, ignore_inactive_links=True, out=None):
    """Map the vector sum of links around a patch to the patch.

    The resulting vector is returned as a length-2 list, with the two
    items being arrays of the x component and the y component of the resolved
    vectors at the patches, respectively.

    Parameters
    ----------
    grid : ModelGrid
        A landlab ModelGrid.
    var_name : array or field name
        Values defined at links.
    ignore_inactive_links : bool
        If True, do not incorporate inactive links into calc. If all links are
        inactive at a patch, record zero if out is None or leave the existing
        value if out.
    out : len-2 list of npatches-long arrays, optional
        Buffer to place mapped values into or ``None`` to create a new array.

    Returns
    -------
    len-2 list of arrays
        [x_component_of_link_vals_at_patch, y_component_of_link_vals_at_patch].

    Examples
    --------
    >>> import numpy as np
    >>> from landlab.grid.mappers import map_link_vector_sum_to_patch
    >>> from landlab import HexModelGrid

    >>> mg = HexModelGrid((4, 3))
    >>> interior_nodes = mg.status_at_node == mg.BC_NODE_IS_CORE
    >>> exterior_nodes = mg.status_at_node != mg.BC_NODE_IS_CORE

    Add a ring of closed nodes at the edge:

    >>> mg.status_at_node[exterior_nodes] = mg.BC_NODE_IS_CLOSED

    This gives us 5 core nodes, 7 active links, and 3 present patches

    >>> (mg.number_of_core_nodes == 5 and mg.number_of_active_links == 7)
    True
    >>> A = mg.add_ones("vals", at="link")
    >>> A.fill(9.0)  # any old values on the inactive links
    >>> A[mg.active_links] = np.array([1.0, -1.0, 1.0, -1.0, -1.0, -1.0, -1.0])

    This setup should give present patch 0 pure east, patch 1 zero (vorticity),
    and patch 2 westwards and downwards components.

    >>> xcomp, ycomp = map_link_vector_sum_to_patch(mg, "vals")
    >>> xcomp, ycomp = np.round(xcomp, decimals=5), np.round(ycomp, decimals=5)
    >>> np.allclose(xcomp[(6, 9, 10),], [2.0, 0.0, -1.0])
    True
    >>> np.allclose(ycomp[(6, 9, 10),] / np.sqrt(3.0), [0.0, 0.0, -1.0])
    True

    These are the patches with ``LinksStatus.INACTIVE`` on all three sides:

    >>> absent_patches = np.array([0, 1, 2, 4, 8, 11, 12, 15, 16, 17, 18])
    >>> np.allclose(xcomp[absent_patches], 0.0)
    True
    >>> np.allclose(ycomp[absent_patches], 0.0)
    True

    Now demonstrate the remaining functionality:

    >>> A = mg.at_link["vals"].copy()
    >>> A.fill(1.0)
    >>> _ = map_link_vector_sum_to_patch(
    ...     mg, A, ignore_inactive_links=False, out=[xcomp, ycomp]
    ... )
    >>> np.allclose(xcomp[absent_patches], 0.0)
    False
    >>> np.allclose(ycomp[absent_patches], 0.0)
    False

    :meta landlab: info-patch, info-link, map
    """
    if out is None:
        out = [
            np.zeros(grid.number_of_patches, dtype=float),
            np.zeros(grid.number_of_patches, dtype=float),
        ]
    else:
        assert len(out) == 2

    if type(var_name) is str:
        var_name = grid.at_link[var_name]
    angles_at_links = grid.angle_of_link  # CCW round tail
    hoz_cpt = np.cos(angles_at_links)
    vert_cpt = np.sin(angles_at_links)
    hoz_vals = var_name * hoz_cpt
    vert_vals = var_name * vert_cpt
    hoz_vals_at_patches = hoz_vals[grid.links_at_patch]
    vert_vals_at_patches = vert_vals[grid.links_at_patch]
    if ignore_inactive_links:
        linkmask = grid.status_at_link[grid.links_at_patch] == grid.BC_LINK_IS_INACTIVE
        hoz_vals_at_patches = np.ma.array(
            hoz_vals_at_patches, mask=linkmask, copy=False
        )
        vert_vals_at_patches = np.ma.array(
            vert_vals_at_patches, mask=linkmask, copy=False
        )
        hoz_sum = np.sum(hoz_vals_at_patches, axis=1)
        vert_sum = np.sum(vert_vals_at_patches, axis=1)
        if type(hoz_sum.mask) is not np.bool_:  # the 2 comps have same mask
            gooddata = np.logical_not(hoz_sum.mask)
            out[0][gooddata] = hoz_sum.data[gooddata]
            out[1][gooddata] = vert_sum.data[gooddata]
        else:
            if not hoz_sum.mask:
                out[0][:] = hoz_sum.data
                out[1][:] = vert_sum.data

    else:
        hoz_sum = np.sum(hoz_vals_at_patches, axis=1, out=out[0])
        vert_sum = np.sum(vert_vals_at_patches, axis=1, out=out[1])

    return out


def map_link_vector_components_to_node(grid, data_at_link):
    """Map (x,y) components of link data data_at_link onto nodes.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid, HexModelGrid

    >>> grid = RasterModelGrid((3, 4))
    >>> link_data = np.arange(grid.number_of_links)

    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vx[5:7]
    array([7.5, 8.5])

    >>> grid = HexModelGrid((3, 3))
    >>> link_data = np.zeros(grid.number_of_links) + 0.5 * 3.0**0.5
    >>> link_data[np.isclose(grid.angle_of_link, 0.0)] = 0.0
    >>> vx, vy = map_link_vector_components_to_node(grid, link_data)
    >>> vy
    array([0.,  0.,  0.,  0.,  1.,  1.,  0.,  0.,  0.,  0.])
    """
    from landlab import HexModelGrid
    from landlab import RasterModelGrid

    if isinstance(grid, HexModelGrid):
        from .hex_mappers import map_link_vector_components_to_node_hex

        return map_link_vector_components_to_node_hex(grid, data_at_link)
    elif isinstance(grid, RasterModelGrid):
        from .raster_mappers import map_link_vector_components_to_node_raster

        return map_link_vector_components_to_node_raster(grid, data_at_link)
    else:
        raise NotImplementedError("Only available for HexModelGrid")


def map_node_to_link_linear_upwind(grid, v, u, out=None):
    """Assign values to links from upstream nodes.

    Assign to each link the value `v` associated with whichever of its two
    nodes lies upstream, according to link value `u`.

    Consider a link k with tail node `t(k)` and head node `t(h)`. Nodes have
    value `v(n)`. We want to assign a value to links, `v'(k)`. The assignment is::

        v'(k) = v(t(k)) where u(k) > 0,
        v'(k) = v(h(k)) where u(k) <= 0

    As an example, consider 3x5 raster grid with the following values
    at the nodes in the central row::

        0---1---2---3---4

    Consider a uniform velocity value `u = 1` at the horizontal links.
    The mapped link values should be::

        .-0-.-1-.-2-.-3-.

    If `u < 0`, the link values should be::

        .-1-.-2-.-3-.-4-.

    Parameters
    ----------
    grid : ModelGrid
        A *Landlab* grid.
    v : (n_nodes,), ndarray
        Values at grid nodes.
    u : (n_links,) ndarray
        Values at grid links.
    out : (n_links,) ndarray, optional
        If provided, place calculated values in this array. Otherwise, create a
        new array.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 5))
    >>> grid.at_node["node_value"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 1.0, 2.0, 3.0, 4.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> v = grid.at_node["node_value"]

    >>> u = grid.add_zeros("advection_speed", at="link")
    >>> u[grid.horizontal_links] = 1.0

    Set values for the middle row of horizontal links.

    >>> val_at_link = map_node_to_link_linear_upwind(grid, v, u)
    >>> val_at_link[9:13]
    array([0.,  1.,  2.,  3.])
    >>> val_at_link = map_node_to_link_linear_upwind(grid, v, -u)
    >>> val_at_link[9:13]
    array([1.,  2.,  3.,  4.])
    """
    if out is None:
        out = np.empty_like(v, shape=(grid.number_of_links,))

    u_is_positive = u > 0.0
    out[u_is_positive] = v[grid.node_at_link_tail[u_is_positive]]
    out[~u_is_positive] = v[grid.node_at_link_head[~u_is_positive]]
    return out


def map_node_to_link_lax_wendroff(grid, v, c, out=None):
    """Assign values to links using a weighted combination of node values.

    Assign to each link a weighted combination of values `v` at nodes
    using the Lax-Wendroff method for upwind weighting.

    `c` is a scalar or link vector that gives the link-parallel signed
    Courant number. Where `c` is positive, velocity is in the direction of
    the link; where negative, velocity is in the opposite direction.

    As an example, consider 3x5 raster grid with the following values
    at the nodes in the central row::

        0---1---2---3---4

    Consider a uniform Courant value `c = +0.2` at the horizontal links.
    The mapped link values should be::

        .-0.4-.-1.4-.-2.4-.-3.4-.

    Values at links when `c = -0.2`::

        .-0.6-.-1.6-.-2.6-.-3.6-.

    Parameters
    ----------
    grid : ModelGrid
        A *Landlab* grid.
    v : (n_nodes,) ndarray
        Values at grid nodes.
    c : float or (n_links,) ndarray
        Courant number to use at links.
    out : (n_links,) ndarray, optional
        If provided, place calculated values in this array. Otherwise, create a
        new array.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid

    >>> grid = RasterModelGrid((3, 5))
    >>> grid.at_node["node_value"] = [
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 1.0, 2.0, 3.0, 4.0],
    ...     [0.0, 0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> v = grid.at_node["node_value"]

    >>> c = grid.add_zeros("courant_number", at="link")
    >>> c[grid.horizontal_links] = 0.2

    Set values for the middle row of horizontal links.

    >>> val_at_link = map_node_to_link_lax_wendroff(grid, v, c)
    >>> val_at_link[9:13]
    array([0.4,  1.4,  2.4,  3.4])
    >>> val_at_link = map_node_to_link_lax_wendroff(grid, v, -c)
    >>> val_at_link[9:13]
    array([0.6,  1.6,  2.6,  3.6])
    """
    if out is None:
        out = np.empty_like(v, shape=(grid.number_of_links,))

    out[:] = 0.5 * (
        (1 + c) * v[grid.node_at_link_tail] + (1 - c) * v[grid.node_at_link_head]
    )

    return out


def map_vectors_to_links(grid, ux, uy, out=None):
    """Map magnitude and sign of vectors with components (ux, uy) onto grid links.

    Examples
    --------
    >>> from landlab import HexModelGrid
    >>> import numpy
    >>> hmg = HexModelGrid((3, 2))
    >>> (numpy.round(10 * map_vectors_to_links(hmg, 1.0, 0.0))).astype(int)
    array([10, -5,  5, -5,  5, 10, 10,  5, -5,  5, -5, 10])
    """
    if out is None:
        out = np.zeros(grid.number_of_links)
    u, theta_u = cartesian_to_polar(ux, uy)
    theta = theta_u - grid.angle_of_link
    out[:] = u * np.cos(theta)
    return out
