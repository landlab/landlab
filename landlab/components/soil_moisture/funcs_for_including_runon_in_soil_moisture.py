# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 11:48:13 2017

@author: Sai Nudurupati & Erkan Istanbulluoglu
"""

import numpy as np
from landlab.components import FlowRouter


def get_ordered_cells_for_soil_moisture(grid, outlet_id=None):
    """
    Runs Landlab's FlowRouter and DepressionFinderAndRouter to
    route flow. Also orders the cells in the descending order of
    channel length (upstream cell order).

    Parameters
    ----------
    grid: grid object
        RasterModelGrid
    outlet_id: int (Optional)
        Outlet id to be set

    Returns
    -------
    ordered_cells: np.array(dtype=int)
        cells ordered in descending order of channel length
    grid: grid object
        updated RasterModelGrid
    """

    if outlet_id is None:
        outlet_id = np.argmin(grid.at_node['topographic__elevation'])
    grid.set_watershed_boundary_condition_outlet_id(
        outlet_id,
        grid.at_node['topographic__elevation'],
        nodata_value=-9999., )
    flw_r = FlowRouter(grid)
    flw_r.run_one_step()
    r = grid.at_node['flow__receiver_node'][grid.node_at_core_cell]
    R = np.zeros(grid.number_of_nodes, dtype=int)
    R[grid.node_at_core_cell] = r
    channel_length = np.zeros(grid.number_of_nodes, dtype=int)
    # Compute channel lengths for each node in the wtrshd (node_at_core_cell)
    for node in grid.node_at_core_cell:
        node_c = node.copy()
        while R[node_c] != node_c:
            channel_length[node] += 1
            node_c = R[node_c]
    grid.at_node['channel_length'] = channel_length
    # Sorting nodes in the ascending order of channel length
    # NOTE: length of ordered_nodes = grid.number_of_core_cells
    ordered_nodes = grid.node_at_core_cell[
        np.argsort(channel_length[grid.node_at_core_cell])]
    # Sorting nodes in the descending order of channel length
    ordered_nodes = ordered_nodes[::-1]
    dd = 1  # switch 2 for while loop
    count_loops = 0  # No. of loops while runs
    while dd:
        dd = 0
        count_loops += 1
        sorted_order = list(ordered_nodes)
        alr_counted_ = []
        for node_ in sorted_order:
            donors = []
            donors = list(grid.node_at_core_cell[np.where(r == node_)[0]])
            if len(donors) != 0:
                for k in range(0, len(donors)):
                    if donors[k] not in alr_counted_:
                        sorted_order.insert(donors[k], sorted_order.pop(sorted_order.index(node_)))
                        dd = 1
            alr_counted_.append(node_)
        ordered_nodes = np.array(sorted_order)
    ordered_cells = grid.cell_at_node[ordered_nodes]
    return ordered_cells, grid
