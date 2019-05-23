#!/usr/env/python

"""
lossy_flow_accumulator.py: Component to accumulate flow and calc drainage area,
while permitting gain or loss of discharge during flow.

DEJH, late 2018
"""

import sys

from landlab.components.flow_accum import (
    FlowAccumulator,
    flow_accum_bw,
    flow_accum_to_n,
)

if sys.version_info[0] >= 3:
    from inspect import signature


class LossyFlowAccumulator(FlowAccumulator):

    """Component to calculate drainage area and accumulate flow, while
    permitting dynamic loss or gain of flow downstream.

    This component is closely related to the FlowAccumulator, in that
    this is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge. However,
    this component additionally requires the passing of a function that
    describes how discharge is lost or gained downstream,
    f(Qw, nodeID, linkID, grid). See the Examples below to see how this works
    in practice.

    Optionally, spatially variable runoff can be set either by the model grid
    field 'water__unit_flux_in' or the input variable *runoff_rate**.

    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished
    together. Note that the DepressionFinderAndRouter is not particularly
    intelligent when running on lossy streams, and in particular, it will
    reroute flow around pits even when they are in fact not filled due to loss.

    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even
    if the  gradients from them point inwards to the main body of the grid.
    This is because under Landlab definitions, perimeter nodes lack cells, so
    cannot accumulate any discharge.

    LossyFlowAccumulator stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*
        -  Node array of discharges: *'surface_water__discharge'*
        -  Node array of discharge loss in transit (vol/sec). This is the
            total loss across all of the downstream links:
            *'surface_water__discharge_loss'*
        -  Node array containing downstream-to-upstream ordered list of node
           IDs: *'flow__upstream_node_order'*
        -  Node array of all but the first element of the delta data structure:
            *flow__data_structure_delta*. The first element is always zero.
        -  At-grid D data structure: *flow__data_structure_D*

    The FlowDirector component will add additional ModelGrid fields; see the
    `FlowAccumulator component <https://landlab.readthedocs.io/en/latest/landlab.components.flow_accum.html>`_
    for full details. These are:

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
            there is no receiver: *'flow__receiver_node'*
        -  Node array of flow proportions: *'flow__receiver_proportions'*
        -  Node array of links carrying flow:  *'flow__link_to_receiver_node'*
        -  Node array of downhill slopes from each receiver:
            *'topographic__steepest_slope'*
        -  Boolean node array of all local lows: *'flow__sink_flag'*

    The primary method of this class is :func:`run_one_step`.

    `run_one_step` takes the optional argument update_flow_director (default is
    True) that determines if the flow_director is re-run before flow is
    accumulated.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    flow_director : string, class, instance of class.
        A string of method or class name (e.g. 'D8' or 'FlowDirectorD8'), an
        uninstantiated FlowDirector class, or an instance of a FlowDirector
        class. This sets the method used to calculate flow directions.
        Default is 'FlowDirectorSteepest'
    runoff_rate : field name, array, or float, optional (m/time)
        If provided, sets the runoff rate and will be assigned to the grid
        field 'water__unit_flux_in'. If a spatially and and temporally variable
        runoff rate is desired, pass this field name and update the field
        through model run time. If both the field and argument are present at
        the time of initialization, runoff_rate will *overwrite* the field. If
        neither are set, defaults to spatially constant unit input.
    depression_finder : string, class, instance of class, optional
        A string of class name (e.g., 'DepressionFinderAndRouter'), an
        uninstantiated DepressionFinder class, or an instance of a
        DepressionFinder class.
        This sets the method for depression finding.
    loss_function : Python function, optional
        A function of the form f(Qw, [node_ID, [linkID, [grid]]]), where Qw is
        the discharge at a node, node_ID the ID of the node at which the loss
        is to be calculated, linkID is the ID of the link down which the
        outflow drains (or a d8 ID if the routing is d8), and grid is a Landlab
        ModelGrid. The function then returns the new discharge at the node
        after the function is applied.
        Note that if a linkID is needed, a nodeID must also be specified, even
        if only as a dummy parameter; similarly, if a grid is to be passed, all
        of the preceding parameters must be specified. Both nodeID and linkID
        are required to permit spatially variable losses, and also losses
        dependent on flow path geometry (e.g., flow length). The grid is passed
        to allow fields or grid properties describing values across the grid
        to be accessed for the loss calculation (see examples).
        This function expects (float, [int, [int, [ModelGrid]]]), and
        return a single float, the new discharge value. This behavior is
        verified during component instantiation.
    **kwargs : optional
        Any additional parameters to pass to a FlowDirector or
        DepressionFinderAndRouter instance (e.g., partion_method for
        FlowDirectorMFD). This will have no effect if an instantiated
        component is passed using the flow_director or depression_finder
        keywords.

    Examples
    --------
    These examples pertain only to the LossyFlowAccumulator. See the main
    FlowAccumulator documentation for more generic and comprehensive examples.

    First, a very simple example. Here's a 50% loss of discharge every time
    flow moves along a node:

    >>> import numpy as np
    >>> from landlab import RasterModelGrid, HexModelGrid
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab.components import DepressionFinderAndRouter

    >>> mg = RasterModelGrid((3, 5), xy_spacing=(2, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    >>> z = mg.add_field('topographic__elevation',
    ...                  mg.node_x + mg.node_y,
    ...                  at='node')

    >>> def mylossfunction(qw):
    ...     return 0.5 * qw

    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest,
    ...                           routing='D4', loss_function=mylossfunction)
    >>> fa.run_one_step()

    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[ 0.,  0.,  0.,  0.,  0.],
           [ 6.,  6.,  4.,  2.,  0.],
           [ 0.,  0.,  0.,  0.,  0.]])
    >>> mg.at_node['surface_water__discharge'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 1.75,  3.5 ,  3.  ,  2.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])
    >>> mg.at_node['surface_water__discharge_loss'].reshape(mg.shape)
    array([[ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ],
           [ 0.  ,  1.75,  1.5 ,  1.  ,  0.  ],
           [ 0.  ,  0.  ,  0.  ,  0.  ,  0.  ]])

    Here we use a spatially distributed field to derive loss terms, and also
    use a filled, non-raster grid.

    >>> dx=(2./(3.**0.5))**0.5  # area to be 100.
    >>> hmg = HexModelGrid(5,3, dx, xy_of_lower_left=(-1.0745, 0.))
    >>> z = hmg.add_field('topographic__elevation',
    ...                   hmg.node_x**2 + np.round(hmg.node_y)**2,
    ...                   at = 'node')
    >>> z[9] = -10.  # poke a hole
    >>> lossy = hmg.add_zeros('node', 'mylossterm', dtype=float)
    >>> lossy[14] = 1.  # suppress all flow from node 14

    Without loss looks like this:

    >>> fa = LossyFlowAccumulator(hmg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest,
    ...                           depression_finder=DepressionFinderAndRouter)
    >>> fa.run_one_step()
    >>> hmg.at_node['flow__receiver_node'] # doctest: +NORMALIZE_WHITESPACE
    array([ 0,  1,  2,
            3,  0,  9,  6,
            7,  9,  4,  9, 11,
           12,  9,  9, 15,
           16, 17, 18])
    >>> hmg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 7.,  0.,  0.,
            0.,  7.,  1.,  0.,
            0.,  1.,  6.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])
    >>> hmg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
    array([ 7.,  0.,  0.,
            0.,  7.,  1.,  0.,
            0.,  1.,  6.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])

    With loss looks like this:

    >>> def mylossfunction2(Qw, nodeID, linkID, grid):
    ...     return (1. - grid.at_node['mylossterm'][nodeID]) * Qw
    >>> fa = LossyFlowAccumulator(hmg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorSteepest,
    ...                           depression_finder=DepressionFinderAndRouter,
    ...                           loss_function=mylossfunction2)
    >>> fa.run_one_step()
    >>> hmg.at_node['drainage_area'] # doctest: +NORMALIZE_WHITESPACE
    array([ 7.,  0.,  0.,
            0.,  7.,  1.,  0.,
            0.,  1.,  6.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])
    >>> hmg.at_node['surface_water__discharge'] # doctest: +NORMALIZE_WHITESPACE
    array([ 6.,  0.,  0.,
            0.,  6.,  1.,  0.,
            0.,  1.,  5.,  1.,  0.,
            0., 1.,  1.,  0.,
            0.,  0.,  0.])
    >>> np.allclose(hmg.at_node['surface_water__discharge_loss'],
    ...             lossy*hmg.at_node['surface_water__discharge'])
    True

    (Loss is only happening from the node, 14, that we set it to happen at.)


    Finally, note we can use the linkIDs to create flow-length-dependent
    effects:

    >>> from landlab.components import FlowDirectorMFD
    >>> mg = RasterModelGrid((4, 6), xy_spacing=(1, 2))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    >>> z = mg.add_field('node', 'topographic__elevation', 2.*mg.node_x)
    >>> z[9] = 8.
    >>> z[16] = 6.5  # force the first node sideways

    >>> L = mg.add_zeros('node', 'spatialloss')
    >>> mg.at_node['spatialloss'][9] = 1.
    >>> mg.at_node['spatialloss'][13] = 1.
    >>> def fancyloss(Qw, nodeID, linkID, grid):
    ...     # now a true transmission loss:
    ...     Lt = (1. - 1./grid.length_of_link[linkID]**2)
    ...     Lsp = grid.at_node['spatialloss'][nodeID]
    ...     return Qw * (1. - Lt) * (1. - Lsp)

    >>> fa = LossyFlowAccumulator(mg, 'topographic__elevation',
    ...                           flow_director=FlowDirectorMFD,
    ...                           loss_function=fancyloss)
    >>> fa.run_one_step()

    >>> mg.at_node['drainage_area'].reshape(mg.shape)
    array([[  0. ,   0. ,   0. ,   0. ,   0. ,   0. ],
           [  5.6,   5.6,   3.6,   2. ,   2. ,   0. ],
           [ 10.4,  10.4,   8.4,   6.4,   4. ,   0. ],
           [  0. ,   0. ,   0. ,   0. ,   0. ,   0. ]])
    >>> mg.at_node['surface_water__discharge'].reshape(mg.shape)
    array([[ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ],
           [ 4. ,  4. ,  2. ,  2. ,  2. ,  0. ],
           [ 0. ,  8.5,  6.5,  4.5,  2.5,  0. ],
           [ 0. ,  0. ,  0. ,  0. ,  0. ,  0. ]])
    """

    _name = "LossyFlowAccumulator"

    _input_var_names = ("topographic__elevation", "water__unit_flux_in")

    _output_var_names = (
        "drainage_area",
        "surface_water__discharge",
        "surface_water__discharge_loss",
        "flow__upstream_node_order",
        "flow__nodes_not_in_stack",
        "flow__data_structure_delta",
        "flow__data_structure_D",
    )

    _var_units = {
        "topographic__elevation": "m",
        "flow__receiver_node": "m",
        "water__unit_flux_in": "m/s",
        "drainage_area": "m**2",
        "surface_water__discharge": "m**3/s",
        "surface_water__discharge_loss": "m**3/s",
        "flow__upstream_node_order": "-",
        "flow__data_structure_delta": "-",
        "flow__data_structure_D": "-",
        "flow__nodes_not_in_stack": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "flow__receiver_node": "node",
        "water__unit_flux_in": "node",
        "drainage_area": "node",
        "surface_water__discharge": "node",
        "surface_water__discharge_loss": "node",
        "flow__upstream_node_order": "node",
        "flow__nodes_not_in_stack": "grid",
        "flow__data_structure_delta": "node",
        "flow__data_structure_D": "grid",
    }
    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "flow__receiver_node": "Node array of receivers (node that "
        + "receives flow from current node)",
        "drainage_area": "Upstream accumulated surface area contributing "
        + "to the node's discharge",
        "surface_water__discharge": "Discharge of water through each node",
        "surface_water__discharge_loss": "Total volume of water per second "
        + "lost during all flow out of the node",
        "water__unit_flux_in": "External volume water per area per time "
        + "input to each node (e.g., rainfall rate)",
        "flow__upstream_node_order": "Node array containing "
        + "downstream-to-upstream ordered list of node IDs",
        "flow__data_structure_delta": "Node array containing the "
        + "elements delta[1:] of the data structure 'delta' used for "
        + "construction of the downstream-to-upstream node array",
        "flow__data_structure_D": "Grid array containing the data structure "
        + "D used for construction of the downstream-to-upstream node array",
        "flow__nodes_not_in_stack": "Boolean value indicating if there "
        + "are any nodes that have not yet been added to the stack stored "
        + "in flow__upstream_node_order.",
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_director="FlowDirectorSteepest",
        runoff_rate=None,
        depression_finder=None,
        loss_function=None,
        **kwargs
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and
        compatability for the flow_director and depression_finder
        keyword arguments, tests the argument of runoff_rate, and
        initializes new fields.
        """
        super(LossyFlowAccumulator, self).__init__(
            grid,
            surface=surface,
            flow_director=flow_director,
            runoff_rate=runoff_rate,
            depression_finder=depression_finder,
            **kwargs
        )

        if loss_function is not None:
            if sys.version_info[0] >= 3:
                sig = signature(loss_function)
                num_params = len(sig.parameters)
            else:  # Python 2
                num_params = loss_function.func_code.co_argcount
            # save the func for loss, and do a quick test on its inputs:
            if num_params == 1:
                # check the func takes a single value and turns it into a new
                # single value:
                if not isinstance(loss_function(1.0), float):
                    raise TypeError(
                        "The loss_function should take a float, and return " "a float."
                    )
                # now, for logical consistency in our calls to
                # find_drainage_area_and_discharge, wrap the func so it has two
                # arguments:

                def lossfunc(Qw, dummyn, dummyl, dummygrid):
                    return float(loss_function(Qw))

                self._lossfunc = lossfunc

            elif num_params == 2:
                # check the func takes a single value and turns it into a new
                # single value:
                if not isinstance(loss_function(1.0, 0), float):
                    raise TypeError(
                        "The loss_function should take (float, int), and "
                        "return a float."
                    )
                # now, for logical consistency in our calls to
                # find_drainage_area_and_discharge, wrap the func so it has two
                # arguments:

                def lossfunc(Qw, nodeID, dummyl, dummygrid):
                    return float(loss_function(Qw, nodeID))

                self._lossfunc = lossfunc

            elif num_params == 3:
                # check the func takes (float, int) and turns it into a new
                # single value:
                if not isinstance(loss_function(1.0, 0, 0), float):
                    raise TypeError(
                        "The loss_function should take (float, int, int), "
                        "and return a float."
                    )

                def lossfunc(Qw, nodeID, linkID, dummygrid):
                    return float(loss_function(Qw, nodeID, linkID))

                self._lossfunc = lossfunc

            elif num_params == 4:
                # this time, the test is too hard to implement cleanly so just
                self._lossfunc = loss_function
            else:
                raise ValueError(
                    "The loss_function must have only a single argument, "
                    "which should be the discharge at a node; a pair of "
                    "arguments, which should be the discharge at a node and "
                    "the node ID; or three arguments, which should be the "
                    "discharge at a node, the node ID, and the link along "
                    "which that discharge will flow."
                )
        else:
            # make a dummy
            def lossfunc(Qw, dummyn, dummyl, dummygrid):
                return float(Qw)

            self._lossfunc = lossfunc

        # add the new loss discharge field if necessary:
        if "surface_water__discharge_loss" not in grid.at_node:
            self.grid.add_zeros(
                "node", "surface_water__discharge_loss", dtype=float, noclobber=False
            )

    def _accumulate_A_Q_to_one(self, s, r):
        """Accumulate area and discharge for a route-to-one scheme."""
        link = self._grid.at_node["flow__link_to_receiver_node"]
        a, q = flow_accum_bw.find_drainage_area_and_discharge_lossy(
            s,
            r,
            link,
            self._lossfunc,
            self._grid,
            self.node_cell_area,
            self._grid.at_node["water__unit_flux_in"],
        )
        return a, q

    def _accumulate_A_Q_to_n(self, s, r, p):
        """Accumulate area and discharge for a route-to-one scheme."""
        link = self._grid.at_node["flow__link_to_receiver_node"]
        a, q = flow_accum_to_n.find_drainage_area_and_discharge_to_n_lossy(
            s,
            r,
            link,
            p,
            self._lossfunc,
            self._grid,
            self.node_cell_area,
            self._grid.at_node["water__unit_flux_in"],
        )
        return a, q


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
