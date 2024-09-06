#! /usr/bin/env python
"""Functions to calculate flow distance."""
import numpy as np

from landlab import FieldError
from landlab import RasterModelGrid


def calculate_flow__distance(grid, add_to_grid=False, clobber=False):
    """Calculate the along flow distance from node to outlet.

    This utility calculates the along flow distance based on the results of
    running flow accumulation on the grid. It will use the connectivity
    used by the FlowAccumulator (e.g. D4, D8, Dinf).

    Parameters
    ----------
    grid : ModelGrid
    add_to_grid : boolean, optional
        Flag to indicate if the stream length field should be added to the
        grid. Default is False. The field name used is ``flow__distance``.
    clobber : boolean, optional
        Flag to indicate if adding the field to the grid should not clobber an
        existing field with the same name. Default is False.

    Returns
    -------
    flow__distance : float ndarray
        The distance that has to be covered from an imaginary flow, located in
        each node of the grid, to reach the watershed's outlet.

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.flow__distance import calculate_flow__distance
    >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    >>> mg.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 21.0, 10.0, 0.0],
    ...     [0.0, 31.0, 20.0, 0.0],
    ...     [0.0, 32.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=True,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> fr.run_one_step()
    >>> flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)
    >>> mg.at_node["flow__distance"]
    array([0.        ,  0.        ,  0.        ,  0.        ,
           0.        ,  1.        ,  0.        ,  0.        ,
           0.        ,  1.41421356,  1.        ,  0.        ,
           0.        ,  2.41421356,  2.        ,  0.        ,
           0.        ,  0.        ,  0.        ,  0.        ])

    Now, let's change to D4 the flow_director method, which does not
    consider diagonal links bewtween nodes.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.flow__distance import calculate_flow__distance
    >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    >>> mg.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 21.0, 10.0, 0.0],
    ...     [0.0, 31.0, 20.0, 0.0],
    ...     [0.0, 32.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=True,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> fr = FlowAccumulator(mg, flow_director="D4")
    >>> fr.run_one_step()
    >>> flow__distance = calculate_flow__distance(mg, add_to_grid=True, clobber=True)
    >>> mg.at_node["flow__distance"]
    array([0.,  0.,  0.,  0.,
           0.,  1.,  0.,  0.,
           0.,  2.,  1.,  0.,
           0.,  3.,  2.,  0.,
           0.,  0.,  0.,  0.])

    The flow__distance utility can also work on irregular grids. For the example we
    will use a Hexagonal Model Grid, a special type of Voroni Grid that has
    regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.flow__distance import calculate_flow__distance
    >>> dx = 1
    >>> hmg = HexModelGrid((5, 3), spacing=dx)
    >>> _ = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x + np.round(hmg.node_y),
    ...     at="node",
    ... )
    >>> hmg.status_at_node[hmg.boundary_nodes] = hmg.BC_NODE_IS_CLOSED
    >>> hmg.status_at_node[0] = hmg.BC_NODE_IS_FIXED_VALUE
    >>> fr = FlowAccumulator(hmg, flow_director="D4")
    >>> fr.run_one_step()
    >>> flow__distance = calculate_flow__distance(hmg, add_to_grid=True, clobber=True)
    >>> hmg.at_node["flow__distance"]
    array([0.,  0.,  0.,
           0.,  1.,  2.,  0.,
           0.,  2.,  2.,  3.,  0.,
           0.,  3.,  3.,  0.,
           0.,  0.,  0.])
    """
    # check that flow__receiver nodes exists
    if "flow__receiver_node" not in grid.at_node:
        raise FieldError(
            "A 'flow__receiver_node' field is required at the "
            "nodes of the input grid."
        )
    if "flow__upstream_node_order" not in grid.at_node:
        raise FieldError(
            "A 'flow__upstream_node_order' field is required at the "
            "nodes of the input grid."
        )

    # get the reciever nodes, depending on if this is to-one, or to-multiple,
    # we'll need to get a different at-node field.
    if grid.at_node["flow__receiver_node"].size != grid.size("node"):
        to_one = False
    else:
        to_one = True
    flow__receiver_node = grid.at_node["flow__receiver_node"]

    # get the upstream node order
    flow__upstream_node_order = grid.at_node["flow__upstream_node_order"]

    # get downstream flow link lengths, result depends on type of grid.
    if isinstance(grid, RasterModelGrid):
        flow_link_lengths = grid.length_of_d8[
            grid.at_node["flow__link_to_receiver_node"]
        ]
    else:
        flow_link_lengths = grid.length_of_link[
            grid.at_node["flow__link_to_receiver_node"]
        ]

    # create an array that representes the outlet lengths.
    flow__distance = np.zeros(grid.nodes.size)

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
                downstream_stream_length = flow__distance[reciever]

                # get the stream segment length from this node to its downstream
                # neigbor
                stream_increment_length = flow_link_lengths[node]

            else:
                # non-existant links are coded with -1
                useable_receivers = np.where(reciever != grid.BAD_INDEX)[0]

                # we will have the stream flow to the downstream node with the
                # shortest distance to the outlet.
                # in the event of a tie, we will choose the shorter link length.

                # get the flow distances of the downstream nodes
                potential_downstream_stream_lengths = flow__distance[
                    flow__receiver_node[node]
                ][useable_receivers]

                # get the stream segment lengths from this node to its downstream
                # neighbor
                potential_stream_increment_lengths = flow_link_lengths[node][
                    useable_receivers
                ]

                # get the lowest downstream stream length.
                downstream_stream_length = np.min(potential_downstream_stream_lengths)

                # determine which of the stream increments flowed to this
                # downstream neighbor.
                which_link = np.where(
                    potential_downstream_stream_lengths == downstream_stream_length
                )[0]

                # and choose the smallest of these links.
                stream_increment_length = np.min(
                    potential_stream_increment_lengths[which_link]
                )

            # set the total stream length of this node
            flow__distance[node] = downstream_stream_length + stream_increment_length

    # store on the grid
    if add_to_grid:
        grid.add_field("flow__distance", flow__distance, at="node", clobber=clobber)

    return flow__distance
