#! /usr/bin/env python
"""Functions to work with watersheds of model grids."""

from landlab import FieldError
import numpy as np


def get_watershed_mask(grid, outlet_id):
    """
    Get the watershed of an outlet returned as a boolean array.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    outlet_id : integer
        The id of the outlet node.

    Returns
    -------
    watershed_mask : boolean ndarray
        True elements of this array correspond to nodes with flow that is
        received by the outlet. The length of the array is equal to the grid
        number of nodes.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils import get_watershed_mask

    >>> rmg = RasterModelGrid((7, 7), 1)
    >>> z = np.array([
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
    ...     -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
    ...     -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
    ...     -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
    ...     -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
    ...     -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])
    >>> rmg.at_node['topographic__elevation'] = z

    Only the bottom boundary is set to open.
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> rmg.set_fixed_value_boundaries_at_grid_edges(False, False, False, True)

    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()

    >>> get_watershed_mask(rmg, 2)
    array([False, False,  True, False, False, False, False, False, False,
            True, False, False, False, False, False,  True,  True,  True,
            True,  True, False, False,  True,  True,  True,  True,  True,
           False, False,  True,  True,  True,  True,  True, False, False,
            True,  True,  True,  True,  True, False, False, False, False,
           False, False, False, False], dtype=bool)
    """
    if 'flow__receiver_node' not in grid.at_node:
        raise FieldError("A 'flow__receiver_node' field is required at the "
                         "nodes of the input grid.")

    grid_nodes = grid.nodes.flatten()
    receiver_at_node = grid.at_node['flow__receiver_node']

    # Prepare output.
    watershed_mask = np.zeros(grid.number_of_nodes, dtype=bool)

    for node in grid_nodes:
        # Follow flow path of each node.
        receiver_node = receiver_at_node[node]
        outlet_not_found = True

        while outlet_not_found:
            node_flows_to_outlet = any([receiver_node == outlet_id, 
                                        receiver_at_node[receiver_node] ==
                                        outlet_id])
            node_is_outlet = node == outlet_id

            if node_flows_to_outlet or node_is_outlet:
                watershed_mask[node] = True
                outlet_not_found = False

            else:
                receiver_node = receiver_at_node[receiver_node]
       
                if receiver_node == receiver_at_node[receiver_node]:
                    # Receiver_node is a pit.
                    outlet_not_found = False        

    return watershed_mask


def get_watershed_nodes(grid, outlet_id):
    """
    Get the watershed of an outlet returned as a list of nodes.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    outlet_id : integer
        The id of the outlet node.

    Returns
    -------
    watershed_nodes : integer ndarray
        The ids of the nodes that flow to the node with the id, outlet_id.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils import get_watershed_nodes

    >>> rmg = RasterModelGrid((7, 7), 1)
    >>> z = np.array([
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
    ...     -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
    ...     -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
    ...     -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
    ...     -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
    ...     -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])

    >>> rmg.at_node['topographic__elevation'] = z
    >>> rmg.set_watershed_boundary_condition_outlet_id(2, z,
    ...                                                nodata_value=-9999.)

    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()

    Get the nodes of two watersheds.
    >>> mainstem_watershed_nodes = get_watershed_nodes(rmg, 2)
    >>> tributary_watershed_nodes = get_watershed_nodes(rmg, 24)

    Given the watershed boundary conditions, the number of mainstem watershed
    nodes should be equal to the number of core nodes plus the outlet node.
    >>> len(mainstem_watershed_nodes) == rmg.number_of_core_nodes + 1
    True
    """
    ws_mask = get_watershed_mask(grid, outlet_id)
    ws_nodes = np.where(ws_mask)[0]

    return ws_nodes


def get_watershed_masks_with_area_threshold(grid, critical_area):
    """
    Get masks of all of the watersheds with a minimum drainage area size.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    critical_area : integer or float
        The minimum drainage area of the watersheds to identify.

    Returns
    -------
    watershed_masks : integer ndarray
        The length of the array is equal to the grid number of nodes. Values of
        this array are the watershed ids. The value of a watershed id is the
        node id of the watershed outlet. Nodes with a value of -1 have only
        downstream nodes with drainage areas below `critical_area`.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils import get_watershed_masks_with_area_threshold

    Create a grid with a node spacing of 200 meter.
    >>> rmg = RasterModelGrid((7, 7), 200)
    >>> z = np.array([
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
    ...     -9999.,    26.,     0.,    26.,    30.,    34., -9999.,
    ...     -9999.,    28.,     1.,    28.,     5.,    32., -9999.,
    ...     -9999.,    30.,     3.,    30.,    10.,    34., -9999.,
    ...     -9999.,    32.,    11.,    32.,    15.,    38., -9999.,
    ...     -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])
    >>> rmg.at_node['topographic__elevation'] = z
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)

    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()

    Get the masks of watersheds greater than or equal to 80,000 square-meters.
    >>> critical_area = 80000
    >>> mask = get_watershed_masks_with_area_threshold(rmg, critical_area)

    # Verify that all mask null nodes have a drainage area below critical area.
    >>> null_nodes = np.where(mask == -1)[0]
    >>> A = rmg.at_node['drainage_area'][null_nodes]
    >>> below_critical_area_nodes = A < critical_area
    >>> np.all(below_critical_area_nodes)
    True
    """

    nodes = grid.nodes.flatten()
    core_nodes = grid.node_is_core(nodes)
    boundary_nodes = grid.node_is_boundary(nodes)

    # Find outlets along grid boundary.
    A = grid.at_node['drainage_area']
    stream_nodes = A >= critical_area
    outlet_nodes = np.all([stream_nodes, boundary_nodes], 0)
    outlet_nodes = np.where(outlet_nodes)[0]

    # Get watersheds of boundary outlets.
    watershed_masks = np.repeat(-1, grid.number_of_nodes)
    for outlet in outlet_nodes:
        _assign_outlet_id_to_watershed(grid, outlet, watershed_masks)

    # Get internal watersheds.
    node_is_unresolved = np.all([A >= critical_area, watershed_masks == -1,
                                 core_nodes], 0)
    while np.any(node_is_unresolved):
        outlet = np.all([A == max(A[node_is_unresolved]), node_is_unresolved],
                         0)

        # Get the outlet with the greatest drainage area of the unresolved
        # nodes.
        outlet = np.where(outlet)[0][0]
        _assign_outlet_id_to_watershed(grid, outlet, watershed_masks)
        node_is_unresolved = np.all([A >= critical_area, watershed_masks == -1,
                                     core_nodes], 0)

    return watershed_masks


def _assign_outlet_id_to_watershed(grid, outlet, masks):
    ws_mask = get_watershed_mask(grid, outlet)
    ws_number_of_nodes = len(np.where(ws_mask)[0])
    masks[ws_mask] = np.repeat(outlet, ws_number_of_nodes)


def get_watershed_outlet(grid, source_node_id):
    """
    Get the downstream-most node (the outlet) of the source node.

    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    source_node_id : integer
        The id of the node in which to identify its outlet.

    Returns
    -------
    outlet_node : integer
        The id of the node that is the downstream-most node (the outlet) of the
        source node.

    Examples
    --------

    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils import get_watershed_outlet

    >>> rmg = RasterModelGrid((7, 7), 1)
    >>> z = np.array([
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
    ...     -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
    ...     -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
    ...     -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
    ...     -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
    ...     -9999.,    34.,    32.,    34.,    36.,    40., -9999.,
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])

    >>> rmg.at_node['topographic__elevation'] = z
    >>> imposed_outlet = 2
    >>> rmg.set_watershed_boundary_condition_outlet_id(imposed_outlet, z,
    ...                                                nodata_value=-9999.)

    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()

    Get the grid watershed outlet.
    >>> determined_outlet = get_watershed_outlet(rmg, 40)
    >>> determined_outlet == imposed_outlet
    True
    """
    if 'flow__receiver_node' not in grid.at_node:
        raise FieldError("A 'flow__receiver_node' field is required at the "
                         "nodes of the input grid.")

    receiver_at_node = grid.at_node['flow__receiver_node']
    receiver_node = receiver_at_node[source_node_id]
    outlet_not_found = True

    while outlet_not_found:
        node_is_outlet = receiver_node == source_node_id
        node_is_boundary = grid.node_is_boundary(receiver_node)
        node_is_pit = receiver_node == receiver_at_node[receiver_node]

        if node_is_outlet or node_is_boundary or node_is_pit:
            outlet_not_found = False
            outlet_node = receiver_node
        else:
            receiver_node = receiver_at_node[receiver_node]

    return outlet_node
