"""

This file contain a series of functions that prepares and set nodes and links
characteristics that are used in boundary conditions in the river_bed_dynamics
landlab component. Examples of applications are given in each function.

.. codeauthor:: Angel Monsalve
.. codecoauthors: Sam Anderson, Nicole Gasparini, Elowyn Yager

"""

import numpy as np


def fixed_links(sed_transp__bedload_gsd_fix_link):
    """
    Search and identify links defined as having a fixed bed load GSD
    This function is used only when the bedload gsd is fixed on a link.
    This function is automatically called during the execution of the component
    but for verification purposes we test it here.
    The example is based on the main example given in river_bed_dynamics.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from . import _nodes_and_links_info as info

    >>> grid = RasterModelGrid((5, 5))
    >>> gsd_loc = [
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ...     [0, 1.0, 1.0, 1.0, 0],
    ... ]
    >>> gsd = [[128, 100], [64, 90], [32, 80], [16, 50], [8, 20], [2, 0]]
    >>> qb_fix_gsd = np.zeros((grid.number_of_links, np.array(gsd).shape[0] - 1))

    Let's impose a bedload gsd to links # 15 and 29

    >>> qb_fix_gsd[[15, 29], :] = np.array([0.15, 0.30, 0.2, 0.2, 0.15])

    Let's check which link is has an fix bedload gsd (it should be 29)

    >>> print(list(info.fixed_links(qb_fix_gsd)))
    [15, 29]

    """
    # Gives the links id for which the bed load gsd is fixed
    return np.where(np.any(sed_transp__bedload_gsd_fix_link > 0, axis=1))[0]


def outlet_nodes(grid):
    """
    Search and identify the node upstream the outlet to apply boundary
    conditions

    This function is automatically called during the execution of the component
    but for verification purposes we test it here.
    The example is based on the main example given in river_bed_dynamics.
    We will explore the horizontal link id of the outlet in different cases

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from . import _nodes_and_links_info as info

    >>> grid = RasterModelGrid((5, 5))

    Case1: In this topography the outlet is at the left edge

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.07, 1.08, 1.09, 1.09, 1.09],
    ...     [1.06, 1.07, 1.08, 1.09, 1.09],
    ...     [1.00, 1.03, 1.07, 1.08, 1.09],
    ...     [1.06, 1.07, 1.08, 1.09, 1.09],
    ...     [1.07, 1.08, 1.09, 1.09, 1.09],
    ... ]
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])
    >>> (
    ...     out_id,
    ...     upstream_out_id,
    ...     outlet_links,
    ...     closed_nodes,
    ...     boundary_links,
    ... ) = info.outlet_nodes(grid)

    The nodes where the flow will exit is:

    >>> out_id[0]  # [0] is used to avoid displaying dtype
    10

    All links connected to the outlet and the node upstream the outlet are:

    >>> outlet_links
    array([13, 14, 18, 19, 22, 23])

    Case 2: In this topography the outlet is at the top edge

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.09, 1.09, 1.09, 1.09, 1.09],
    ...     [1.09, 1.09, 1.08, 1.09, 1.09],
    ...     [1.09, 1.08, 1.07, 1.08, 1.09],
    ...     [1.08, 1.07, 1.03, 1.07, 1.08],
    ...     [1.07, 1.06, 1.00, 1.06, 1.07],
    ... ]
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])
    >>> (
    ...     out_id,
    ...     upstream_out_id,
    ...     outlet_links,
    ...     closed_nodes,
    ...     boundary_links,
    ... ) = info.outlet_nodes(grid)

    >>> out_id[0]
    22

    >>> outlet_links
    array([24, 28, 29, 33, 37, 38])

    Case 3: In this topography the outlet is at the right edge

    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.09, 1.09, 1.09, 1.08, 1.07],
    ...     [1.09, 1.09, 1.08, 1.07, 1.06],
    ...     [1.09, 1.08, 1.07, 1.03, 1.00],
    ...     [1.09, 1.09, 1.08, 1.07, 1.06],
    ...     [1.09, 1.09, 1.09, 1.08, 1.07],
    ... ]
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])
    >>> (
    ...     out_id,
    ...     upstream_out_id,
    ...     outlet_links,
    ...     closed_nodes,
    ...     boundary_links,
    ... ) = info.outlet_nodes(grid)

    >>> out_id[0]
    14

    >>> outlet_links
    array([16, 17, 20, 21, 25, 26])

    """

    # Gives the ID of the outlet node
    (out_id,) = np.where(grid.status_at_node == 1)
    upstream_out_id = grid.active_adjacent_nodes_at_node[out_id][
        np.where(grid.active_adjacent_nodes_at_node[out_id] != -1)
    ]
    (closed_id,) = np.where(grid.status_at_node == 4)

    outlet_nodes = np.sort(np.concatenate((out_id, upstream_out_id)))
    outlet_links = np.unique(grid.links_at_node[outlet_nodes].ravel())
    outlet_links = outlet_links[outlet_links >= 0]

    if closed_id.size > 0:
        closed_nodes = np.sort(np.array(closed_id)).flatten()
    else:
        closed_nodes = []

    closed_links = np.unique(grid.links_at_node[closed_nodes].ravel())
    closed_links = closed_links[closed_links >= 0]

    links_at_border_cells = links_at_border(grid)
    boundary_links = np.sort(
        np.unique(np.hstack([closed_links, links_at_border_cells]))
    )

    return out_id, upstream_out_id, outlet_links, closed_nodes, boundary_links


def links_at_border(grid):
    """
    Search and identify links connected to edge nodes.
    This function is automatically called during the execution of the component
    but for verification purposes we test it here.
    The example is based on the main example given in river_bed_dynamics.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from . import _nodes_and_links_info as info

    >>> grid = RasterModelGrid((5, 5))
    >>> grid.at_node["topographic__elevation"] = [
    ...     [1.07, 1.08, 1.09, 1.09, 1.09],
    ...     [1.06, 1.07, 1.08, 1.09, 1.09],
    ...     [1.00, 1.03, 1.07, 1.08, 1.09],
    ...     [1.06, 1.07, 1.08, 1.09, 1.09],
    ...     [1.07, 1.08, 1.09, 1.09, 1.09],
    ... ]
    >>> grid.set_watershed_boundary_condition(grid.at_node["topographic__elevation"])

    >>> links_at_border_cells = info.links_at_border(grid)
    >>> links_at_border_cells
    array([ 0,  3,  4,  5,  6,  7,  8,  9, 12, 18, 21, 27, 30, 31, 32, 33, 34,
           35, 36, 39])

    This is a list of all links that could lead the flow directly to a potential outlet
    For example, link 1 is not listed because is horizontal and the outlet through node
    1 or 2 will be vertical
    """
    links_at_border_cells = np.sort(
        np.hstack(
            (
                grid.links_at_node[:, 2][grid.nodes_at_right_edge],
                grid.links_at_node[:, 3][grid.nodes_at_top_edge],
                grid.links_at_node[:, 0][grid.nodes_at_left_edge],
                grid.links_at_node[:, 1][grid.nodes_at_bottom_edge],
            )
        )
    )
    return links_at_border_cells
