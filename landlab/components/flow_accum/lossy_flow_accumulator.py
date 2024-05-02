"""Accumulate flow and calc drainage area, while permitting gain or loss
of discharge during flow.

DEJH, late 2018
"""

import sys

from landlab.components.flow_accum import FlowAccumulator
from landlab.components.flow_accum import flow_accum_bw
from landlab.components.flow_accum import flow_accum_to_n

if sys.version_info[0] >= 3:
    from inspect import signature


class LossyFlowAccumulator(FlowAccumulator):
    """Component to calculate drainage area and accumulate flow, while
    permitting dynamic loss or gain of flow downstream.

    This component is closely related to the :class:`.FlowAccumulator`,
    in that this is accomplished by first finding flow directions by a user-specified
    method and then calculating the drainage area and discharge. However,
    this component additionally requires the passing of a function that
    describes how discharge is lost or gained downstream::

        f(Qw, nodeID, linkID, grid)

    See the Examples below to see how this works in practice.

    Optionally, spatially variable runoff can be set either by the model grid
    field ``"water__unit_flux_in"`` or the input variable ``runoff_rate``.

    Optionally a depression finding component can be specified and flow
    directing, depression finding, and flow routing can all be accomplished
    together. Note that the :class:`.DepressionFinderAndRouter`
    is not particularly intelligent when running on lossy streams, and in particular,
    it will reroute flow around pits even when they are in fact not filled due to loss.

    .. note::

        The perimeter nodes *NEVER* contribute to the accumulating flux, even
        if the  gradients from them point inwards to the main body of the grid.
        This is because under Landlab definitions, perimeter nodes lack cells, so
        cannot accumulate any discharge.

    :class:`~.LossyFlowAccumulator` stores as :class:`~.ModelGrid` fields:

    -  Node array of drainage areas: ``"drainage_area"``
    -  Node array of discharges: ``"surface_water__discharge"``
    -  Node array of discharge loss in transit (vol / sec). This is the
       total loss across all of the downstream links:
       ``"surface_water__discharge_loss"``
    -  Node array containing downstream-to-upstream ordered list of node
       IDs: ``"flow__upstream_node_order"``
    -  Node array of all but the first element of the delta data structure:
       ``"flow__data_structure_delta"``. The first element is always zero.

    The :class:`FlowDirector` component will add additional
    :class:`~.ModelGrid` fields; see the
    :class:`~.FlowAccumulator` for full details. These are:

    -  Node array of receivers (nodes that receive flow), or *ITS OWN ID* if
       there is no receiver: ``"flow__receiver_node"``.
    -  Node array of flow proportions: ``"flow__receiver_proportions"``.
    -  Node array of links carrying flow:  ``"flow__link_to_receiver_node"``.
    -  Node array of downhill slopes from each receiver:
       ``"topographic__steepest_slope"``.
    -  Boolean node array of all local lows: ``"flow__sink_flag"``.

    The primary method of this class is
    :meth:`~landlab.components.LossyFlowAccumulator.run_one_step`.

    Examples
    --------
    These examples pertain only to the :class:`~.LossyFlowAccumulator`.
    See the main :class:`~.FlowAccumulator` documentation for more
    generic and comprehensive examples.

    First, a very simple example. Here's a 50% loss of discharge every time
    flow moves along a node:

    >>> import numpy as np
    >>> from landlab import RasterModelGrid, HexModelGrid
    >>> from landlab.components import FlowDirectorSteepest
    >>> from landlab.components import DepressionFinderAndRouter

    >>> mg = RasterModelGrid((3, 5), xy_spacing=(2, 1))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    >>> z = mg.add_field("topographic__elevation", mg.node_x + mg.node_y, at="node")

    >>> def mylossfunction(qw):
    ...     return 0.5 * qw
    ...

    >>> fa = LossyFlowAccumulator(
    ...     mg,
    ...     "topographic__elevation",
    ...     flow_director=FlowDirectorSteepest,
    ...     loss_function=mylossfunction,
    ... )
    >>> fa.run_one_step()

    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[0., 0., 0., 0., 0.],
           [6., 6., 4., 2., 0.],
           [0., 0., 0., 0., 0.]])
    >>> mg.at_node["surface_water__discharge"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  ],
           [1.75, 3.5 , 3.  , 2.  , 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  ]])
    >>> mg.at_node["surface_water__discharge_loss"].reshape(mg.shape)
    array([[0.  , 0.  , 0.  , 0.  , 0.  ],
           [0.  , 1.75, 1.5 , 1.  , 0.  ],
           [0.  , 0.  , 0.  , 0.  , 0.  ]])

    Here we use a spatially distributed field to derive loss terms, and also
    use a filled, non-raster grid.

    >>> dx = (2.0 / (3.0**0.5)) ** 0.5  # area to be 100.0
    >>> hmg = HexModelGrid((5, 3), spacing=dx, xy_of_lower_left=(-1.0745, 0.0))
    >>> z = hmg.add_field(
    ...     "topographic__elevation",
    ...     hmg.node_x**2 + np.round(hmg.node_y) ** 2,
    ...     at="node",
    ... )
    >>> z[9] = -10.0  # poke a hole
    >>> lossy = hmg.add_zeros("mylossterm", dtype=float, at="node")
    >>> lossy[14] = 1.0  # suppress all flow from node 14

    Without loss looks like this:

    >>> fa = LossyFlowAccumulator(
    ...     hmg,
    ...     "topographic__elevation",
    ...     flow_director=FlowDirectorSteepest,
    ...     depression_finder=DepressionFinderAndRouter,
    ... )
    >>> fa.run_one_step()
    >>> hmg.at_node["flow__receiver_node"]
    array([ 0,  1,  2,
            3,  0,  9,  6,
            7,  9,  4,  9, 11,
           12,  9,  9, 15,
           16, 17, 18])
    >>> np.round(hmg.at_node["drainage_area"])
    array([7., 0., 0.,
           0., 7., 1., 0.,
           0., 1., 6., 1., 0.,
           0., 1., 1., 0.,
           0., 0., 0.])
    >>> np.round(hmg.at_node["surface_water__discharge"])
    array([7., 0., 0.,
           0., 7., 1., 0.,
           0., 1., 6., 1., 0.,
           0., 1., 1., 0.,
           0., 0., 0.])

    With loss looks like this:

    >>> def mylossfunction2(Qw, nodeID, linkID, grid):
    ...     return (1.0 - grid.at_node["mylossterm"][nodeID]) * Qw
    ...
    >>> fa = LossyFlowAccumulator(
    ...     hmg,
    ...     "topographic__elevation",
    ...     flow_director=FlowDirectorSteepest,
    ...     depression_finder=DepressionFinderAndRouter,
    ...     loss_function=mylossfunction2,
    ... )
    >>> fa.run_one_step()
    >>> np.round(hmg.at_node["drainage_area"])
    array([7., 0., 0.,
           0., 7., 1., 0.,
           0., 1., 6., 1., 0.,
           0., 1., 1., 0.,
           0., 0., 0.])
    >>> np.round(hmg.at_node["surface_water__discharge"])
    array([6., 0., 0.,
           0., 6., 1., 0.,
           0., 1., 5., 1., 0.,
           0., 1., 1., 0.,
           0., 0., 0.])
    >>> np.allclose(
    ...     hmg.at_node["surface_water__discharge_loss"],
    ...     lossy * hmg.at_node["surface_water__discharge"],
    ... )
    True

    (Loss is only happening from the node, 14, that we set it to happen at.)

    Finally, note we can use the *linkIDs* to create flow-length-dependent
    effects:

    >>> from landlab.components import FlowDirectorMFD
    >>> mg = RasterModelGrid((4, 6), xy_spacing=(1, 2))
    >>> mg.set_closed_boundaries_at_grid_edges(True, True, False, True)
    >>> z = mg.add_field("topographic__elevation", 2.0 * mg.node_x, at="node")
    >>> z[9] = 8.0
    >>> z[16] = 6.5  # force the first node sideways

    >>> L = mg.add_zeros("spatialloss", at="node")
    >>> mg.at_node["spatialloss"][9] = 1.0
    >>> mg.at_node["spatialloss"][13] = 1.0
    >>> def fancyloss(Qw, nodeID, linkID, grid):
    ...     # now a true transmission loss:
    ...     Lt = 1.0 - 1.0 / grid.length_of_link[linkID] ** 2
    ...     Lsp = grid.at_node["spatialloss"][nodeID]
    ...     return Qw * (1.0 - Lt) * (1.0 - Lsp)
    ...

    >>> fa = LossyFlowAccumulator(
    ...     mg,
    ...     "topographic__elevation",
    ...     flow_director=FlowDirectorMFD,
    ...     loss_function=fancyloss,
    ... )
    >>> fa.run_one_step()

    >>> mg.at_node["drainage_area"].reshape(mg.shape)
    array([[ 0. ,  0. , 0. ,  0. ,  0. ,  0. ],
           [ 5.6,  5.6, 3.6,  2. ,  2. ,  0. ],
           [10.4, 10.4, 8.4,  6.4,  4. ,  0. ],
           [ 0. ,  0. , 0. ,  0. ,  0. ,  0. ]])
    >>> mg.at_node["surface_water__discharge"].reshape(mg.shape)
    array([[0. , 0. , 0. , 0. , 0. , 0. ],
           [4. , 4. , 2. , 2. , 2. , 0. ],
           [0. , 8.5, 6.5, 4.5, 2.5, 0. ],
           [0. , 0. , 0. , 0. , 0. , 0. ]])

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    None Listed

    **Additional References**

    Braun, J., Willett, S. (2013). A very efficient O(n), implicit and parallel
    method to solve the stream power equation governing fluvial incision and
    landscape evolution. Geomorphology  180-181(C), 170-179.
    https://dx.doi.org/10.1016/j.geomorph.2012.10.008
    """

    _name = "LossyFlowAccumulator"

    _info = {
        "drainage_area": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flow__data_structure_delta": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Node array containing the elements delta[1:] of the data "
                "structure 'delta' used for construction of the "
                "downstream-to-upstream node array"
            ),
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
        },
        "surface_water__discharge_loss": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Total volume of water per second lost during all flow out of the node",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
        "water__unit_flux_in": {
            "dtype": float,
            "intent": "in",
            "optional": True,
            "units": "m/s",
            "mapping": "node",
            "doc": (
                "External volume water per area per time input to each node "
                "(e.g., rainfall rate)"
            ),
        },
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_director="FlowDirectorSteepest",
        runoff_rate=None,
        depression_finder=None,
        loss_function=None,
        **kwargs,
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and
        compatability for the flow_director and depression_finder
        keyword arguments, tests the argument of runoff_rate, and
        initializes new fields.

        Parameters
        ----------
        grid : ModelGrid
            A Landlab grid.
        surface : str, int or array of int
            The surface to direct flow across.  Can be a field name at node or
            an array of length *node*.
        flow_director : str, class or instance of class.
            A string of method or class name (e.g. ``'D8'`` or ``'FlowDirectorD8'``), an
            uninstantiated FlowDirector class, or an instance of a
            :class:`FlowDirector` class. This sets the method
            used to calculate flow directions.
        runoff_rate : field name, array, or float, optional (m/time)
            If provided, sets the runoff rate and will be assigned to the grid
            field ``'water__unit_flux_in'``. If a spatially and and temporally variable
            runoff rate is desired, pass this field name and update the field
            through model run time. If both the field and argument are present at
            the time of initialization, runoff_rate will *overwrite* the field. If
            neither are set, defaults to spatially constant unit input.
        depression_finder : str, class, instance of class, optional
            A string of class name (e.g., ``'DepressionFinderAndRouter'``), an
            uninstantiated :class:`~.DepressionFinder` class,
            or an instance of a :class:`~.DepressionFinder` class.
            This sets the method for depression finding.
        loss_function : func, optional
            A function of the form ``f(Qw, [node_ID, [linkID, [grid]]])``, where Qw is
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
        """

        # add the new loss discharge field if necessary:
        if "surface_water__discharge_loss" not in grid.at_node:
            grid.add_zeros(
                "node", "surface_water__discharge_loss", dtype=float, clobber=True
            )

        super().__init__(
            grid,
            surface=surface,
            flow_director=flow_director,
            runoff_rate=runoff_rate,
            depression_finder=depression_finder,
            **kwargs,
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

    def _accumulate_A_Q_to_one(self, s, r):
        """Accumulate area and discharge for a route-to-one scheme."""
        link = self._grid.at_node["flow__link_to_receiver_node"]
        a, q = flow_accum_bw.find_drainage_area_and_discharge_lossy(
            s,
            r,
            link,
            self._lossfunc,
            self._grid,
            self._node_cell_area,
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
            self._node_cell_area,
            self._grid.at_node["water__unit_flux_in"],
        )
        return a, q
