#! /usr/bin/env python
"""Functions to calculate flow distance from divide."""
import numpy as np

from landlab import FieldError
from landlab import RasterModelGrid


def calculate_distance_to_divide(
    grid, longest_path=True, add_to_grid=False, clobber=False
):
    """Calculate the along flow distance from drainage divide to point.

    This utility calculates the along flow distance based on the results of
    running flow accumulation on the grid. It will use the connectivity
    used by the FlowAccumulator (e.g. D4, D8, Dinf).

    Parameters
    ----------
    grid : ModelGrid
    longest_path : bool, optional
        Take the longest (or shortest) path to a drainage divide. Default is
        true.
    add_to_grid : boolean, optional
        Flag to indicate if the stream length field should be added to the
        grid. Default is False. The field name used is ``distance_to_divide``.
    clobber : boolean, optional
        Flag to indicate if adding the field to the grid should not clobber an
        existing field with the same name. Default is False.

    Returns
    -------
    distance_to_divide : float ndarray
        The distance that has to be covered from an imaginary flow, located in
        each node of the grid, to reach the watershed's outlet.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.distance_to_divide import calculate_distance_to_divide
    >>> mg = RasterModelGrid((5, 4))
    >>> mg.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 10.0, 10.0, 0.0],
    ...     [0.0, 20.0, 20.0, 0.0],
    ...     [0.0, 30.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=False,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> fr.run_one_step()
    >>> distance_to_divide = calculate_distance_to_divide(
    ...     mg,
    ...     add_to_grid=True,
    ...     clobber=True,
    ... )
    >>> mg.at_node["distance_to_divide"]
    array([0.,  3.,  3.,  0.,
           0.,  2.,  2.,  0.,
           0.,  1.,  1.,  0.,
           0.,  0.,  0.,  0.,
           0.,  0.,  0.,  0.])

    Now, let's change to MFD the flow_director method, which routes flow to
    multiple nodes.

    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.distance_to_divide import calculate_distance_to_divide
    >>> mg = RasterModelGrid((5, 4), xy_spacing=(1, 1))
    >>> mg.at_node["topographic__elevation"] = [
    ...     [0.0, 0.0, 0.0, 0.0],
    ...     [0.0, 10.0, 10.0, 0.0],
    ...     [0.0, 20.0, 20.0, 0.0],
    ...     [0.0, 30.0, 30.0, 0.0],
    ...     [0.0, 0.0, 0.0, 0.0],
    ... ]
    >>> mg.set_closed_boundaries_at_grid_edges(
    ...     bottom_is_closed=False,
    ...     left_is_closed=True,
    ...     right_is_closed=True,
    ...     top_is_closed=True,
    ... )
    >>> fr = FlowAccumulator(mg, flow_director="MFD")
    >>> fr.run_one_step()
    >>> distance_to_divide = calculate_distance_to_divide(
    ...     mg,
    ...     add_to_grid=True,
    ...     clobber=True,
    ... )
    >>> mg.at_node["distance_to_divide"]
    array([0.,  3.,  3.,  0.,
           0.,  2.,  2.,  0.,
           0.,  1.,  1.,  0.,
           0.,  0.,  0.,  0.,
           0.,  0.,  0.,  0.])

    The distance_to_divide utility can also work on irregular grids. For the
    example we will use a Hexagonal Model Grid, a special type of Voroni Grid
    that has regularly spaced hexagonal cells.

    >>> from landlab import HexModelGrid
    >>> from landlab.components import FlowAccumulator
    >>> from landlab.utils.distance_to_divide import calculate_distance_to_divide
    >>> dx = 1
    >>> hmg = HexModelGrid((5, 3), dx)
    >>> _ = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x + np.round(hmg.node_y),
    ...     at="node",
    ... )
    >>> hmg.status_at_node[hmg.boundary_nodes] = hmg.BC_NODE_IS_CLOSED
    >>> hmg.status_at_node[0] = hmg.BC_NODE_IS_FIXED_VALUE
    >>> fr = FlowAccumulator(hmg, flow_director="D4")
    >>> fr.run_one_step()
    >>> distance_to_divide = calculate_distance_to_divide(
    ...     hmg,
    ...     add_to_grid=True,
    ...     clobber=True,
    ... )
    >>> hmg.at_node["distance_to_divide"]
    array([3.,  0.,  0.,
        0.,  2.,  1.,  0.,
      0.,  1.,  1.,  0.,  0.,
        0.,   0.,  0.,  0.,
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

    if "drainage_area" not in grid.at_node:
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
    drainage_area = grid.at_node["drainage_area"]

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

    # create an array that representes the distance to the divide.
    distance_to_divide = np.zeros(grid.nodes.size)

    if not longest_path:
        distance_to_divide[:] = 2 * grid.size("node") * np.max(flow_link_lengths)

    # iterate through the flow__upstream_node_order backwards.
    for node in reversed(flow__upstream_node_order):
        # if drainage are is equal to node cell area, set distance to zeros
        # this should handle the drainage divide cells as boundary cells have
        # their area set to zero.
        if drainage_area[node] == grid.cell_area_at_node[node]:
            distance_to_divide[node] = 0

        # get flow recievers
        reciever = flow__receiver_node[node]

        if to_one:
            # if not processing an outlet node.
            if reciever != node:
                if longest_path:
                    cond = (
                        distance_to_divide[reciever]
                        < distance_to_divide[node] + flow_link_lengths[node]
                    )
                else:
                    cond = (
                        distance_to_divide[reciever]
                        > distance_to_divide[node] + flow_link_lengths[node]
                    )

                if cond:
                    distance_to_divide[reciever] = (
                        distance_to_divide[node] + flow_link_lengths[node]
                    )

        else:
            # non-existant links are coded with -1
            useable_receivers = np.where(reciever != grid.BAD_INDEX)[0]

            for idx in range(len(useable_receivers)):
                r = reciever[useable_receivers][idx]
                fll = flow_link_lengths[node][useable_receivers][idx]

                # if not processing an outlet node.
                if r != node:
                    if longest_path:
                        cond = distance_to_divide[r] < distance_to_divide[node] + fll
                    else:
                        cond = distance_to_divide[r] > distance_to_divide[node] + fll

                    if cond:
                        distance_to_divide[r] = distance_to_divide[node] + fll

    # store on the grid
    if add_to_grid:
        grid.add_field(
            "distance_to_divide", distance_to_divide, at="node", clobber=clobber
        )

    return distance_to_divide
