import numpy as np

from landlab import FieldError

from cfuncs import _get_watershed_mask


def get_watershed_masks_with_area_threshold(grid, critical_area):
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
    grid_nodes = grid.nodes.flatten()

    unassigned_nodes = grid_nodes[masks == -1]
    flow_receivers = grid.at_node['flow__receiver_node']

    ws_mask = np.zeros(len(grid_nodes), dtype=int)

    _get_watershed_mask(unassigned_nodes, flow_receivers, np.float(outlet),
                        ws_mask)

    ws_mask = np.array(ws_mask, dtype=bool)

    ws_number_of_nodes = len(np.where(ws_mask)[0])
    masks[ws_mask] = np.repeat(outlet, ws_number_of_nodes)
