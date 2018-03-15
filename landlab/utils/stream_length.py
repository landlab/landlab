#! /usr/bin/env python
"""Functions to calculate stream distance."""

from landlab.grid import RasterModelGrid
from landlab import FieldError
import numpy as np


def calculate_stream_length(grid, add_to_grid=False, noclobber=True):
    """
    Calculate the along stream length from node to outlet.

    This utility calculates the along stream distance based on the results of
    running flow accumulation on the grid. It will use the connectivity
    used by the FlowAccumulator (e.g. D4 or D8). It requires that a route to
    one method was used.

    Parameters
    ----------
    grid : ModelGrid
    add_to_grid : boolean, optional
        Flag to indicate if the stream length field should be added to the 
        grid. Default is False. The field name used is ``stream__length``.
    noclobber : boolean, optional
        Flag to indicate if adding the field to the grid should clobber an 
        existing field with the same name. Default is True. 
        
    Examples
    ________
    
    
    """
    # check that flow__reciever nodes exists
    if 'flow__receiver_node' not in grid.at_node:
        raise FieldError("A 'flow__receiver_node' field is required at the "
                         "nodes of the input grid.")
    if 'flow__upstream_node_order' not in grid.at_node:
        raise FieldError("A 'flow__upstream_node_order' field is required at the "
                         "nodes of the input grid.")

    # get an array of grid nodes
    grid_nodes = grid.nodes.flatten()

    # get the reciever nodes
    flow__receiver_node = grid.at_node['flow__receiver_node']

    if grid_nodes.size != flow__receiver_node.size:
        raise ValueError('A route to one method must be used with the calculate '
                         'stream distance method.')

    # get the upstream node order
    flow__upstream_node_order = grid.at_node['flow__upstream_node_order']

    # get downstream flow link lengths, result depends on type of grid.
    if isinstance(grid, RasterModelGrid):
        flow_link_lengths = grid.length_of_d8[grid.at_node['flow__link_to_receiver_node']]
    else:
        flow_link_lengths = grid.length_of_link[grid.at_node['flow__link_to_receiver_node']]

    # create an array that representes the outlet lengths.
    stream__length = np.ones(grid.nodes.size)

    # iterate through the flow__upstream_node_order, this will already have
    # identified the locations of the outlet nodes and have
    for node in flow__upstream_node_order:

        # if not an outlet
        if flow__receiver_node[node] is not node:

            # get the stream length of the downstream node
            downstream_stream_length = stream__length[flow__receiver_node[node]]

            # get the stream segment length from this node to its downstream
            # neigbor
            stream_increment_length = flow_link_lengths[node]

            # set the total stream length of this node
            stream__length[node] = downstream_stream_length + stream_increment_length

    # store on the grid
    if add_to_grid:
        grid.add_field('node', 'stream__length', stream__length, noclobber=noclobber)

    return stream__length