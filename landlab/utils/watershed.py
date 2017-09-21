#! /usr/bin/env python
"""Functions to identify watersheds of model grid nodes."""

from landlab import FieldError
import numpy as np


def get_watershed_array(grid, outlet_id):
    """
    Gets the watershed of an outlet returned as a boolean array.
    
    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    outlet_id : integer
        The id of the outlet node.
    
    Returns
    --------
    watershed_array : boolean ndarray
        True elements of this array correspond to nodes with flow that is
        received by the outlet. Array length is equal to the grid number of
        nodes.
        
    Examples
    --------
    
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils.watershed import get_watershed_ids
    
    >>> rmg = RasterModelGrid((5, 5), 1)
    >>> z = np.array([-9999., -9999., -9999., -9999., -9999.,
    ...               -9999.,    28.,     1.,    28., -9999.,
    ...               -9999.,    30.,     5.,    30., -9999.,
    ...               -9999.,    34.,    32.,    34., -9999.,
    ...               -9999., -9999., -9999., -9999., -9999.])
    >>> rmg.at_node['topographic__elevation'] = z
    
    Only the bottom boundary is set to open.
    >>> rmg.set_closed_boundaries_at_grid_edges(True, True, True, False)
    >>> rmg.set_fixed_value_boundaries_at_grid_edges(False, False, False, True)
     
    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()
    
    >>> get_watershed_array(rmg, 2)
    array([False, False,  True, False, False, False, False,  True, False,
           False, False,  True,  True,  True, False, False,  True,  True,
            True, False, False, False, False, False, False], dtype=bool)
    """
    if 'flow__receiver_node' not in grid.at_node:
        raise FieldError("This method requires a 'flow__receiver_node' "
                         "field at the nodes of the input grid.")
    
    grid_nodes = grid.nodes.flatten()
    receiver_at_node = grid.at_node['flow__receiver_node']

    # Prepare output.
    watershed_array = np.zeros(len(grid_nodes), dtype=bool)

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
                watershed_array[node] = True
                outlet_not_found = False
                
            else:
                receiver_node = receiver_at_node[receiver_node]
                
                if receiver_node == receiver_at_node[receiver_node]:
                    # Receiver_node is a pit.
                    outlet_not_found = False

    return watershed_array

def get_watershed_ids(grid, outlet_id):
    """
    Gets the watershed of an outlet returned as a list of node ids.
    
    Parameters
    ----------
    grid : RasterModelGrid
        A landlab RasterModelGrid.
    outlet_id : integer
        The id of the outlet node.
    
    Returns
    --------
    watershed_ids : integer ndarray
        The ids of the nodes that flow to the node with the id, outlet_id.
        
    Examples
    --------
    
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowRouter
    >>> from landlab.utils.watershed import get_watershed_ids

    >>> rmg = RasterModelGrid((7, 7), 1)
    >>> z = np.array([
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.,
    ...     -9999.,    26.,     0.,    30.,    32.,    34., -9999.,
    ...     -9999.,    28.,     1.,    25.,    28.,    32., -9999.,
    ...     -9999.,    30.,     3.,     3.,    11.,    34., -9999.,
    ...     -9999.,    32.,    11.,    25.,    18.,    38., -9999.,
    ...     -9999.,    34.,    32.,    38.,    36.,    40., -9999.,
    ...     -9999., -9999., -9999., -9999., -9999., -9999., -9999.])
    
    >>> rmg.at_node['topographic__elevation'] = z
    >>> rmg.set_watershed_boundary_condition_outlet_id(2, z,
    ...                                                nodata_value=-9999.)
     
    Route flow.
    >>> fr = FlowRouter(rmg)
    >>> fr.run_one_step()
    
    Get ids of two watersheds.
    >>> mainstem_watershed_ids = get_watershed_ids(rmg, 2)
    >>> tributary_watershed_ids = get_watershed_ids(rmg, 24)
    
    Given the watershed boundary conditions, the number of mainstem watershed
    nodes should be equal to the number of core nodes plus the outlet node.
    >>> len(mainstem_watershed_ids) == rmg.number_of_core_nodes + 1
    True
    
    The tributary recieves flow only from a portion of the grid nodes.
    >>> tributary_watershed_ids
    array([12, 18, 19, 24, 25, 26, 31, 32, 33, 39, 40])
    """
    ws_array = get_watershed_array(grid, outlet_id)
    ws_ids = np.where(ws_array)[0]
    
    return ws_ids
