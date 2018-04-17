#! /usr/bin/env python
"""Functions to calculate stream distance."""
import numpy as np

from landlab import BAD_INDEX_VALUE, RasterModelGrid, FieldError


def calculate_stream_length(grid, add_to_grid=False, noclobber=True):
    """
    Calculate the along stream length from node to outlet.

    This utility calculates the along stream distance based on the results of
    running flow accumulation on the grid. It will use the connectivity
    used by the FlowAccumulator (e.g. D4, D8, Dinf).

    Parameters
    ----------
    grid : ModelGrid
    add_to_grid : boolean, optional
        Flag to indicate if the stream length field should be added to the
        grid. Default is False. The field name used is ``stream__length``.
    noclobber : boolean, optional
        Flag to indicate if adding the field to the grid should not clobber an
        existing field with the same name. Default is True.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.stream_length import calculate_stream_length
    >>> mg = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> elev = np.array([0.,  0.,  0., 0.,
        ...              0., 21., 10., 0.,
        ...              0., 31., 20., 0.,
        ...              0., 32., 30., 0.,
        ...              0.,  0.,  0., 0.])
    >>> mg.add_field('node','topographic__elevation', z)
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    >>> fr = FlowAccumulator(mg, flow_director = 'D8')
    >>> fr.run_one_step()
    >>> flow_length = calculate_stream_length(mg, add_to_grid=True, noclobber=False)
    >>> #Set to zero the flow length values of boundary nodes and watershed's outlet
    >>> flow_length[mg.boundary_nodes] = 0
    >>> outlet_id = 6
    >>> flow_length[outlet_id] = 0
    >>> mg.at_node['flow_length']
    >>> array([  0,  0,  0,  0,
                 0,  1,  0,  0,
                 0,  1.414,  1, 0,
                 0, 2.414, 2, 0,
                 0, 0, 0, 0])
        
        Now, let's change to D4 the flow_director method, that does not 
        consider diagonal links bewtween nodes.
        
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.stream_length import calculate_stream_length
    >>> mg = RasterModelGrid((5, 4), spacing=(1, 1))
    >>> elev = np.array([0.,  0.,  0., 0.,
        ...              0., 21., 10., 0.,
        ...              0., 31., 20., 0.,
        ...              0., 32., 30., 0.,
        ...              0.,  0.,  0., 0.])
    >>> mg.add_field('node','topographic__elevation', z)
    >>> mg.set_closed_boundaries_at_grid_edges(bottom_is_closed=True, left_is_closed=True, right_is_closed=True, top_is_closed=True)
    >>> fr = FlowAccumulator(mg, flow_director = 'D4')
    >>> fr.run_one_step()
    >>> flow_length = calculate_stream_length(mg, add_to_grid=True, noclobber=False)
    >>> #Set to zero the flow length values of boundary nodes and watershed's outlet
    >>> flow_length[mg.boundary_nodes] = 0
    >>> outlet_id = 6
    >>> flow_length[outlet_id] = 0
    >>> mg.at_node['flow_length']
    >>> array([  0,  0,  0,  0,
                 0,  1,  0,  0,
                 0,  2,  1, 0,
                 0, 3, 2, 0,
                 0, 0, 0, 0])
                 
        The stream_length utility can also work on irregular grids. For the example we 
        will use a Hexagonal Model Grid, a special type of Voroni Grid that has 
        regularly spaced hexagonal cells. We will also set the dx spacing such 
        that each cell has an area of one.  
        
    >>> from landlab import HexModelGrid, CLOSED_BOUNDARY
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.stream_length import calculate_stream_length
    >>> dx=(2./(3.**0.5))**0.5
    >>> mg = HexModelGrid(5,3, dx)
    >>> mg.add_field('topographic__elevation', mg.node_x + np.round(mg.node_y), at = 'node')
    >>> hmg.status_at_node[0] = CLOSED_BOUNDARY
    >>> fr = FlowAccumulator(hmg, flow_director = 'D4')
    >>> fr.run_one_step()
    >>> flow_length = calculate_stream_length(hmg, add_to_grid=True, noclobber=False)
    >>> #Even in this case, boundary and outlet nodes should have a flow length equal to zero
    >>> flow_length[hmg.boundary_nodes] = 0
    >>> outlet_id = 4
    >>> flow_length[outlet_id] = 0
    
        
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

    # get the reciever nodes, depending on if this is to-one, or to-multiple,
    # we'll need to get a different at-node field.
    if 'flow__receiver_nodes' in grid.at_node:
        to_one = False
        flow__receiver_node = grid.at_node['flow__receiver_nodes']
    else:
        to_one = True
        flow__receiver_node = grid.at_node['flow__receiver_node']

    # get the upstream node order
    flow__upstream_node_order = grid.at_node['flow__upstream_node_order']

    # get downstream flow link lengths, result depends on type of grid.
    if isinstance(grid, RasterModelGrid):
        if to_one:
            flow_link_lengths = grid.length_of_d8[grid.at_node['flow__link_to_receiver_node']]
        else:
            flow_link_lengths = grid.length_of_d8[grid.at_node['flow__links_to_receiver_nodes']]
    else:
        if to_one:
            flow_link_lengths = grid.length_of_link[grid.at_node['flow__link_to_receiver_node']]
        else:
            flow_link_lengths = grid.length_of_link[grid.at_node['flow__links_to_receiver_nodes']]

    # create an array that representes the outlet lengths.
    stream__length = np.zeros(grid.nodes.size)

    # iterate through the flow__upstream_node_order, this will already have
    # identified the locations of the outlet nodes and have
    for node in flow__upstream_node_order:

        # get flow recievers
        reciever = flow__receiver_node[node]

        # assess if this is a to one (D8/D4) or to multiple (Dinf, MFD)
        # flow directing method.
        if to_one:
            potential_outlet = reciever
        else:
            # if this is an outlet, the first element of the recievers will be
            # the nodes ID.
            potential_outlet = reciever[0]

        # assess if this is an outlet or not.
        if potential_outlet == node:
            not_outlet = False
        else:
            not_outlet = True

        # if not an outlet
        if not_outlet:

            # deal with the two cases of route to one and route to multiple.
            if to_one:
                # get the stream length of the downstream node
                downstream_stream_length = stream__length[reciever]

                # get the stream segment length from this node to its downstream
                # neigbor
                stream_increment_length = flow_link_lengths[node]

            else:
                # non-existant links are coded with -1
                useable_recievers = np.where(reciever != BAD_INDEX_VALUE)[0]

                # we will have the stream flow to the downstream node with the
                # shortest distance to the outlet.
                # in the event of a tie, we will choose the shorter link length.

                # get the stream lengths of the downstream nodes
                potential_downstream_stream_lengths = stream__length[flow__receiver_node[node]][useable_recievers]

                # get the stream segment lengths from this node to its downstream
                # neighbor
                potential_stream_increment_lengths = flow_link_lengths[node][useable_recievers]

                # get the lowest downstream stream length.
                downstream_stream_length = np.min(potential_downstream_stream_lengths)

                # determine which of the stream increments flowed to this
                # downstream neighbor.
                which_link = np.where(potential_downstream_stream_lengths == downstream_stream_length)[0]

                # and choose the smallest of these links.
                stream_increment_length = np.min(potential_stream_increment_lengths[which_link])

            # set the total stream length of this node
            stream__length[node] = downstream_stream_length + stream_increment_length

    # store on the grid
    if add_to_grid:
        grid.add_field('node', 'stream__length', stream__length, noclobber=noclobber)

    return stream__length
