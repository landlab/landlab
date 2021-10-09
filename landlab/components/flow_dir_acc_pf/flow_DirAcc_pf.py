#!/usr/env/python

"""
flow_DirAcc_Pf.py: Component to fill or breach a DEM, accumulate flow and calculate drainage area using the priority flood algorithm.

FlowDirAccPf is a wrapper of the RichDEM package: https://richdem.readthedocs.io/en/latest/flow_metrics.html

The componet merges a filling/breaching algorithm, a flow director as well as a flow accumulator.
Moreover, the component supports the definitation of two flow accumulator fields associated to the same grid.
This prevents the user from updating the filling/breaching algorithms in between calculation of flow accumulater one and two.

@author: benjaminCampforts
"""

import copy as cp
import os
import sys
import numpy as np
import numpy.matlib as npm
import richdem as rd
from landlab import Component, FieldError, RasterModelGrid
from landlab.utils.return_array import return_array_at_node
from .cfuncs import _FA_D8, _FD_D8
from .suppress_stdout_stderr import suppress_stdout_stderr
# Codes for depression status
_UNFLOODED = 0
_PIT = 1
_CURRENT_LAKE = 2
_FLOODED = 3
# Flow metrics resulting in single flow
PSINGLE_FMs = ["D8", "D4", "Rho8", "Rho4"]
# Flow metrics resulting in multiple flow
PMULTIPLE_FMs = ["Quinn", "Freeman", "Holmgren", "Dinf"]


class FlowDirAccPf(Component):

    """Component to accumulate flow and calculate drainage area based RICHDEM software package.
    See also: https://richdem.readthedocs.io/en/latest/


    NOTE: The perimeter nodes  NEVER contribute to the accumulating flux, even
    if the  gradients from them point inwards to the main body of the grid.
    This is because under Landlab definitions, perimeter nodes lack cells, so
    cannot accumulate any discharge.

    FlowAccumulatorPf stores as ModelGrid fields:

        -  Node array of drainage areas: *'drainage_area'*

        -  Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or _FLOODED) *'flood_status_code'*

        -  Node array of discharges: *'surface_water__discharge'*

        -  Distance to receiver: *'Distance to receiver'*

        -  External volume water per area per time input to each node *'water__unit_flux_in'*

        -  Node array containing downstream-to-upstream ordered list of node
            IDs: *'flow__upstream_node_order'*

        -  Node array of receivers (nodes that receive flow), or ITS OWN ID if
            there is no receiver: *'flow__receiver_node'*. This array is 2D for
            RouteToMany methods and has the shape

            (n-nodes x max number of receivers).
        -  Node array of flow proportions: *'flow__receiver_proportions'*. This
            array is 2D, for RouteToMany methods and has the shape
            (n-nodes x max number of receivers).

        -  Node array of downhill slopes from each receiver:
            *'topographic__steepest_slope'* This array is 2D for RouteToMany
            methods and has the shape (n-nodes x max number of receivers).


        -  Node array of links carrying flow
            *'flow__link_to_receiver_node'*

        - Node array of proportion of flow sent to each receiver '*flow__receiver_proportions*'

        - Depression free land surface topographic elevation, at closed borders,
            value equals -1 *'deprFree_elevation'*

        # The following fields are required when an additional
        # hillslope flowrouting scheme is required, can be completed
        # with flow acc and discharge if required


        - Node array containing downstream-to-upstream ordered list of node IDs
            *'hill_flow__upstream_node_order'*

        - Node array of receivers (node that receives flow from current node)
            *'hill_flow__receiver_node'*
        - The steepest *downhill* slope *'hill_topographic__steepest_slope'*

        - Node array of proportion of flow sent to each receiver
            *'hill_flow__receiver_proportions'*


    The primary method of this class is :func:`run_one_step`.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    surface : field name at node or array of length node
        The surface to direct flow across.
    flow_metric: Various options possible:
            D8 (O’Callaghan and Mark, 1984) {default}
            Rho8 (Fairfield and Leymarie, 1991)
            Quinn (1991)
            Freeman (1991)
            Holmgren (1994)
            D∞ (Tarboton, 1997)
            for details and comparison, see https://richdem.readthedocs.io/en/latest/flow_metrics.html


    runoff_rate : field name, array, or float, optional (m/time)
        If provided, sets the runoff rate and will be assigned to the grid field
        'water__unit_flux_in'. If a spatially and and temporally variable runoff
        rate is desired, pass this field name and update the field through model
        run time. If both the field and argument are present at the time of
        initialization, runoff_rate will *overwrite* the field. If neither are
        set, defaults to spatially constant unit input.
        Both a runoff_rate array and the 'water__unit_flux_in' field are
        permitted to contain negative values, in which case they mimic
        transmission losses rather than e.g. rain inputs.

    updateFlowDepressions {True} Build-in depression handler. Can be through filling or breaching (see below). Default is True

    updateHillDepressions {False} Only needed if DEM needs to be filled seperately for second (hill flow) flow accumulator.
    Default behavior is not to execute a seperate filling procedure in between the first and the second flow accumualtor.



    depressionHandler : string, "fill" or "breach" indicating wheather to use a
        Depression-Filling or breaching algorithm to process depressions
            - fill {default}: Depression-Filling
                Depression-filling is often used to fill in all the depressions
                in a DEM to the level of their lowest outlet or spill-point.
                See also: https://richdem.readthedocs.io/en/latest/depression_filling.html
            - breach {default}: Complete Breaching
                Depression-breaching is used to dig channels from the pit cells
                of a DEM to the nearest cells (in priority-flood sense) outside
                of the depression containing the pit. This resolves the depression
                as all cells in the depression now have a drainage path to the
                edge of the DEM.
                See also: https://richdem.readthedocs.io/en/latest/depression_breaching.html

    exponent:  Some methods require an exponent (see flow_metric) Default {1}
    epsilon: boolean, default {true}. If True, an epsilon gradient is imposed
        to all flat regions. This ensures that there is always a local gradient.
    accumulateFlow : if True flow directions and acummulations will be calcualted.
        Set to False when only interested in flow directions
    accumulateFlowHill: if True flow directions and acummulations will be calcualted
        for second FD component (Hill). Set to False when only interested in flow
        directions.
    seperate_Hill_Flow: boolean {False}
        For some applications (e.g. HyLands) both single and
        multiple flow direction and accumulation is required.
        By calculating them in the same component, we can optimize procedures
        invovled with filling and breaching of DEMs
    update_HillFlow_instanteneous: boolean {True}
        update seperate hillslope director and accumulator simultaneously on update
        Set to if other operations have to be performed in between updating the
        principle flow properties and the hillslope properties
    hill_flow_metric:
            various options possible:
            D8 (O’Callaghan and Mark, 1984)
            D4 (O’Callaghan and Mark, 1984)
            Rho8 (Fairfield and Leymarie, 1991)
            Rho4 (Fairfield and Leymarie, 1991)
            Quinn (1991) {default}
            Freeman (1991)
            Holmgren (1994)
            D∞ (Tarboton, 1997)
            for details and comparison, see https://richdem.readthedocs.io/en/latest/flow_metrics.html
    hill_exponent Some methods require an exponent (see flow_metric) {1}
    suppress_out:
        suppress verbose of priority flood algorithm




    Examples
    --------


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion cell digital elevation models on desktops or clusters. Environmental Modelling & Software 92, 202–212. doi: 10.1016/j.envsoft.2017.02.022

    **Additional References**

    https://richdem.readthedocs.io/en/latest/

    """

    _name = "FlowDirAccPf"

    _unit_agnostic = True

    _info = {
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "drainage_area": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flood_status_code": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or _FLOODED).",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m**3/s",
            "mapping": "node",
            "doc": "Volumetric discharge of surface water",
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
            "doc": "External volume water per area per time input to each node (e.g., rainfall rate)",
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "SQUARED_length_adjacent": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Length to adjacent nodes, squared (calcualted in advance to save time during calculation",
        },
        "flow__receiver_proportions": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "deprFree_elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Filled land surface topographic elevation, at closed borders, value equals -1!",
        },
        # The following fields are required when an additional
        # hillslope flowrouting scheme is required, can be completed
        # with flow acc and discharge if required
        "hill_drainage_area": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "hill_surface_water__discharge": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "hill_flow__upstream_node_order": {
            "dtype": int,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "hill_flow__receiver_node": {
            "dtype": int,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "hill_topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "hill_flow__receiver_proportions": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
    }

    def __init__(
        self,
        grid,
        surface="topographic__elevation",
        flow_metric="D8",
        runoff_rate=None,
        updateFlowDepressions=True,
        depressionHandler="fill",
        exponent=1,
        epsilon=True,
        accumulateFlow=True,
        accumulateFlowHill=False,
        seperate_Hill_Flow=False,
        updateHillDepressions=False,
        update_HillFlow_instanteneous=True,
        hill_flow_metric="Quinn",
        hill_exponent=1,
        suppress_out=False,
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests imput types and
        compatability for the flow_metric and depression_finder
        keyword arguments, tests the argument of runoff, and
        initializes new fields.
        """
        super(FlowDirAccPf, self).__init__(grid)
        # Keep a local reference to the grid

        self._suppress_out = suppress_out

        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self._test_water_inputs(grid, runoff_rate)

        # Grid type testing
        if not isinstance(self._grid, RasterModelGrid):
            raise FieldError(
                "Flow Accumulator Priority flood only works with regular raster grids, use default Landlab flow accumulator instead"
            )

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.0
        self._node_cell_area = node_cell_area
        self._runoff_rate = runoff_rate

        if (flow_metric in PSINGLE_FMs) or (flow_metric in PMULTIPLE_FMs):
            self._flow_metric = flow_metric
        else:
            raise ValueError(
                [
                    "flow metric should be one of these single flow directors : "
                    + str(PSINGLE_FMs)
                    + " or multiple flow directors"
                    + str(PMULTIPLE_FMs)
                ]
            )
        if (hill_flow_metric in PSINGLE_FMs) or (hill_flow_metric in PMULTIPLE_FMs):
            self._hill_flow_metric = hill_flow_metric
        else:
            raise ValueError(
                [
                    "flow metric should be one of these single flow directors : "
                    + str(PSINGLE_FMs)
                    + " or multiple flow directors"
                    + str(PMULTIPLE_FMs)
                ]
            )

        if depressionHandler == "fill" or depressionHandler == "breach":
            self._depressionHandler = depressionHandler
        else:
            raise ValueError("depressionHandler should be fill or breach")

        self._depressionHandler = depressionHandler
        self._exponent = exponent
        self._epsilon = epsilon
        self._seperate_Hill_Flow = seperate_Hill_Flow
        self._update_HillFlow_instanteneous = update_HillFlow_instanteneous

        self._updateFlowDepressions = updateFlowDepressions
        self._updateHillDepressions = updateHillDepressions

        self._hill_exponent = hill_exponent

        if self._seperate_Hill_Flow:
            # Adjust dict
            self._info["hill_drainage_area"]["optional"] = False
            self._info["hill_surface_water__discharge"]["optional"] = False
            self._info["hill_flow__upstream_node_order"]["optional"] = False
            self._info["hill_flow__receiver_node"]["optional"] = False
            self._info["hill_topographic__steepest_slope"]["optional"] = False
            self._info["hill_flow__receiver_proportions"]["optional"] = False

        self._accumulateFlow = accumulateFlow
        self._accumulateFlowHill = accumulateFlowHill

        if not self._accumulateFlow:
            self._info["drainage_area"]["optional"] = True
            self._info["surface_water__discharge"]["optional"] = True

        if not self._accumulateFlowHill:
            self._info["hill_drainage_area"]["optional"] = True
            self._info["hill_surface_water__discharge"]["optional"] = True

        self.initialize_output_fields()

        # Make aliases
        if self._accumulateFlow:
            self._drainage_area = self.grid.at_node["drainage_area"]
            self._discharges = self.grid.at_node["surface_water__discharge"]
        self._sort = self.grid.at_node["flow__upstream_node_order"]
        # if multiple flow algorithm is made, the dimensions of the slope and receiver fields change (8 colums for all neightbors)
        if flow_metric in PMULTIPLE_FMs:
            self.grid.at_node["topographic__steepest_slope"] = np.zeros(
                (self.grid.number_of_nodes, 8)
            )
            self.grid.at_node["flow__receiver_node"] = np.zeros(
                (self.grid.number_of_nodes, 8), dtype=int
            )
            self.grid.at_node["flow__receiver_proportions"] = np.zeros(
                (self.grid.number_of_nodes, 8)
            )
            self.grid.at_node["flow__link_to_receiver_node"] = np.zeros(
                (self.grid.number_of_nodes, 8)
            )
        self._slope = self.grid.at_node["topographic__steepest_slope"]
        self._rcvs = self.grid.at_node["flow__receiver_node"]
        self._prps = self.grid.at_node["flow__receiver_proportions"]
        self._recvr_link = self.grid.at_node["flow__link_to_receiver_node"]

        if self._seperate_Hill_Flow:
            if self._accumulateFlowHill:
                self._hill_drainage_area = self.grid.at_node["hill_drainage_area"]
                self._hill_discharges = self.grid.at_node[
                    "hill_surface_water__discharge"
                ]
            if hill_flow_metric in PMULTIPLE_FMs:
                self.grid.at_node["hill_topographic__steepest_slope"] = np.zeros(
                    (self.grid.number_of_nodes, 8)
                )
                self.grid.at_node["hill_flow__receiver_node"] = np.zeros(
                    (self.grid.number_of_nodes, 8), dtype=int
                )
                self.grid.at_node["hill_flow__receiver_proportions"] = np.zeros(
                    (self.grid.number_of_nodes, 8)
                )
            self._hill_slope = self.grid.at_node["hill_topographic__steepest_slope"]
            self._hill_rcvs = self.grid.at_node["hill_flow__receiver_node"]
            self._hill_prps = self.grid.at_node["hill_flow__receiver_proportions"]

        # Create properties specific to RichDEM
        if self._suppress_out:
            with suppress_stdout_stderr():
                sys.stdout = open(os.devnull, "w")
                sys.stderr = open(os.devnull, "w")

                self._create_properties()

                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        else:
            self._create_properties()

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid["node"]["drainage_area"]

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid["node"]["surface_water__discharge"]

    def _create_properties(self):
        self._deprFreeDEM = cp.deepcopy(
            rd.rdarray(
                self.grid.at_node["topographic__elevation"].reshape(self.grid.shape),
                no_data=-9999,
            )
        )

        # Calculate SQUARED length adjacent
        self.grid.at_node["SQUARED_length_adjacent"] = np.concatenate(
            (
                np.ones((self.grid.number_of_nodes, 4)),
                2 * np.ones((self.grid.number_of_nodes, 4)),
            ),
            axis=1,
        )

        self._closed = np.zeros(self._grid.number_of_nodes)
        self._closed[self._grid.status_at_node == 4] = 1
        self._closed = rd.rdarray(self._closed.reshape(self._grid.shape), no_data=-9999)

    def _test_water_inputs(self, grid, runoff_rate):
        """Test inputs for runoff_rate and water__unit_flux_in."""
        if "water__unit_flux_in" not in grid.at_node:
            if runoff_rate is None:
                # assume that if runoff rate is not supplied, that the value
                # should be set to one everywhere.
                grid.add_ones("water__unit_flux_in", at="node", dtype=float)
            else:
                runoff_rate = return_array_at_node(grid, runoff_rate)
                grid.at_node["water__unit_flux_in"] = runoff_rate
        else:
            if runoff_rate is not None:
                print(
                    "FlowAccumulator found both the field "
                    + "'water__unit_flux_in' and a provided float or "
                    + "array for the runoff_rate argument. THE FIELD IS "
                    + "BEING OVERWRITTEN WITH THE SUPPLIED RUNOFF_RATE!"
                )
                runoff_rate = return_array_at_node(grid, runoff_rate)
                grid.at_node["water__unit_flux_in"] = runoff_rate

    def calculateFlowDir_Acc(self, hillFlow=False, updateDepressions=True):
        if hillFlow:
            flowMetric = self._hill_flow_metric
        else:
            flowMetric = self._flow_metric

        # 1: Remove depressions
        if updateDepressions:
            self.removeDepressions()

        # 2: Flow directions and accumulation
        # D8 flow accumulation in richDEM seems not to differntiatie between caridal and diagonal cells,
        #   so we provide an alternative D8 implementation strategy
        if flowMetric == "D8":
            self.FA_DA_D8(hillFlow=hillFlow)
        else:

            # Calculate flow direction (proportion) and accumulation using RichDEM
            # self._deprFreeDEM[self._closed == 1] = 1e6
            props_Pf = rd.FlowProportions(
                dem=self._deprFreeDEM, method=flowMetric, exponent=self._exponent
            )
            # self._deprFreeDEM[self._closed == 1] = -1

            # Calculate flow accumulation using RichDEM
            if hillFlow:
                if self._accumulateFlowHill:
                    self.accumulate_flow_RD(props_Pf, hillFlow=hillFlow)
            else:
                if self._accumulateFlow:
                    self.accumulate_flow_RD(props_Pf, hillFlow=hillFlow)

            # Convert flow proportions to landlab structure
            props_Pf = props_Pf.reshape(
                props_Pf.shape[0] * props_Pf.shape[1], props_Pf.shape[2]
            )
            props_Pf_col0 = props_Pf[:, 0]
            props_Pf = np.column_stack(
                (
                    props_Pf[:, 5],
                    props_Pf[:, 7],
                    props_Pf[:, 1],
                    props_Pf[:, 3],
                    props_Pf[:, 6],
                    props_Pf[:, 8],
                    props_Pf[:, 2],
                    props_Pf[:, 4],
                )
            )
            rcvrs = np.concatenate(
                (
                    self._grid.adjacent_nodes_at_node,
                    self._grid.diagonal_adjacent_nodes_at_node,
                ),
                axis=1,
            )
            rcvrs[props_Pf <= 0] = -1
            val = np.arange(0, props_Pf.shape[0])
            rcvrs[props_Pf_col0 == -1, 0] = val[props_Pf_col0 == -1]

            # Links
            recvr_link = np.array(self._grid.d8s_at_node)
            recvr_link[props_Pf <= 0] = -1

            # make sure sum of proportions equals 1
            if flowMetric in PSINGLE_FMs:
                # update slope; make pointer to slope location first, in order to not break connection to the slope memory location previously established
                slope_temp = np.divide(
                    np.subtract(
                        npm.repmat(
                            self.grid.at_node["topographic__elevation"].reshape(
                                self.grid.number_of_nodes, 1
                            ),
                            1,
                            8,
                        ),
                        self.grid.at_node["topographic__elevation"][rcvrs],
                    ),
                    self.grid.dx
                    * np.sqrt(self.grid.at_node["SQUARED_length_adjacent"]),
                )
                slope_temp[rcvrs == -1] = 0

            else:
                props_Pf[props_Pf_col0 == -1, 0] = 1
                props_Pf = props_Pf.astype(np.float64)  # should be float64
                # Now, make sure sum is 1 in 64 bits
                props_Pf[props_Pf == -1] = 0
                rc64_temp = props_Pf / npm.repmat(
                    np.reshape(props_Pf.sum(axis=1), [props_Pf.shape[0], 1]), 1, 8
                )
                props_Pf[props_Pf[:, 0] != 1, :] = rc64_temp[props_Pf[:, 0] != 1, :]
                props_Pf[props_Pf == 0] = -1

            if hillFlow:
                self._hill_prps[:] = props_Pf
                if flowMetric in PSINGLE_FMs:
                    idx = np.argmax(rcvrs, axis=1)
                    self._hill_rcvs[:] = rcvrs[np.arange(rcvrs.shape[0]), idx]
                    self._hill_slope[:] = slope_temp[np.arange(rcvrs.shape[0]), idx]
                else:
                    self._hill_rcvs[:] = rcvrs
                    # update slope; make pointer to slope location first, in order to not break connection to the slope memory location previously established
                    self._hill_slope[:] = np.divide(
                        np.subtract(
                            npm.repmat(
                                self.grid.at_node["topographic__elevation"].reshape(
                                    self.grid.number_of_nodes, 1
                                ),
                                1,
                                8,
                            ),
                            self.grid.at_node["topographic__elevation"][rcvrs],
                        ),
                        self.grid.dx
                        * np.sqrt(self.grid.at_node["SQUARED_length_adjacent"]),
                    )
                    self._hill_slope[rcvrs == -1] = 0

            else:

                if flowMetric in PSINGLE_FMs:
                    idx = np.argmax(rcvrs, axis=1)
                    self._prps[:] = props_Pf[np.arange(rcvrs.shape[0]), idx]
                    self._rcvs[:] = rcvrs[np.arange(rcvrs.shape[0]), idx]
                    self._slope[:] = slope_temp[np.arange(rcvrs.shape[0]), idx]
                    self._recvr_link[:] = recvr_link[np.arange(rcvrs.shape[0]), idx]
                else:
                    self._prps[:] = props_Pf
                    self._rcvs[:] = rcvrs
                    # update slope; make pointer to slope location first, in order to not break connection to the slope memory location previously established
                    self._slope[:] = np.divide(
                        np.subtract(
                            npm.repmat(
                                self.grid.at_node["topographic__elevation"].reshape(
                                    self.grid.number_of_nodes, 1
                                ),
                                1,
                                8,
                            ),
                            self.grid.at_node["topographic__elevation"][rcvrs],
                        ),
                        self.grid.dx
                        * np.sqrt(self.grid.at_node["SQUARED_length_adjacent"]),
                    )
                    self._slope[rcvrs == -1] = 0

                    self._recvr_link[:] = recvr_link

    def FA_DA_D8(self, hillFlow=False):

        # self._deprFreeDEM[self._closed == 1] = 1e6

        # calcualte D8 in C++
        r = self.grid.number_of_node_rows
        c = self.grid.number_of_node_columns
        dx = self.grid.dx
        sq2 = np.sqrt(2)
        activeCells = np.array(self._grid.status_at_node != 4 + 0, dtype=np.int64)
        receivers = np.array(self.grid.status_at_node, dtype=np.int64)
        distance_receiver = np.zeros((receivers.shape), dtype=float)
        cores = self.grid.core_nodes
        # Make boundaries to save time with conditionals in c loops
        receivers[np.nonzero(self._grid.status_at_node)] = -1
        steepest_slope = np.zeros((receivers.shape), dtype=float)
        el = self._deprFreeDEM.reshape(self.grid.number_of_nodes)
        el_ori = self.grid.at_node["topographic__elevation"]
        dist = np.multiply([1, 1, 1, 1, sq2, sq2, sq2, sq2], dx)
        ngb = np.zeros((8,), dtype=np.int64)
        el_d = np.zeros((8,), dtype=float)

        # Links
        adj_link = np.array(self._grid.d8s_at_node, dtype=np.int64)
        recvr_link = np.zeros((receivers.shape), dtype=np.int64) - 1

        _FD_D8(
            receivers,
            distance_receiver,
            steepest_slope,
            np.array(el),
            el_ori,
            dist,
            ngb,
            cores,
            activeCells,
            el_d,
            r,
            c,
            dx,
            sq2,
            adj_link,
            recvr_link,
        )

        # Calcualte flow acc
        do_FA = False
        if hillFlow:
            if self._accumulateFlowHill:
                do_FA = True
                a = self._hill_drainage_area
                q = self._hill_discharges
        else:
            if self._accumulateFlow:
                do_FA = True
                a = self._drainage_area
                q = self._discharges

        if do_FA:
            if any(self.grid.at_node["water__unit_flux_in"] != 1):
                wg_q = (
                    self.grid.at_node["water__unit_flux_in"]
                    * self.grid.dx
                    * self.grid.dx
                )
                # Only core nodes (status == 0) need to receive a weight
                wg_q[np.nonzero(self._grid.status_at_node)] = 0
                dis = wg_q
            else:
                dis = (
                    np.zeros(self.grid.number_of_nodes, dtype=np.int64)
                    + self._node_cell_area
                )

            da = (
                np.zeros(self.grid.number_of_nodes, dtype=np.int64)
                + self._node_cell_area
            )
            stack = np.int64(self._sort)
            nb = np.int64(len(da))
            _FA_D8(nb, da, dis, stack, receivers)

            a[:] = da
            q[:] = dis

        # Closed nodes flow to themselves
        val = np.arange(0, receivers.shape[0])
        receivers[receivers == -1] = val[receivers == -1]

        # Restore depression free DEM
        # self._deprFreeDEM[self._closed == 1] = -1

        if hillFlow:
            self._hill_prps[:] = np.ones_like(receivers)
            self._hill_rcvs[:] = receivers
            self._hill_slope[:] = steepest_slope
        else:
            self._prps[:] = np.ones_like(receivers)
            self._rcvs[:] = receivers
            self._slope[:] = steepest_slope
            self._recvr_link[:] = recvr_link

    def removeDepressions(self):
        self._deprFreeDEM = cp.deepcopy(
            rd.rdarray(
                self.grid.at_node["topographic__elevation"].reshape(self.grid.shape),
                no_data=-9999,
            )
        )
        if self._depressionHandler == "fill":
            rd.FillDepressions(self._deprFreeDEM, self._epsilon, in_place=True)
        elif self._depressionHandler == "breach":
            rd.BreachDepressions(self._deprFreeDEM, in_place=True)
        self._sort[:] = np.argsort(
            np.array(self._deprFreeDEM.reshape(self.grid.number_of_nodes))
        )
        self.grid.at_node["deprFree_elevation"] = self._deprFreeDEM

    def accumulate_flow_RD(self, props_Pf, hillFlow=False):
        if not hillFlow:
            a = self._drainage_area
            q = self._discharges
        else:
            a = self._hill_drainage_area
            q = self._hill_discharges

        # Create weight for flow accum: both open (status ==1) and closed nodes (status ==4) will have zero weight
        wg = np.ones(self.grid.number_of_nodes) * self.grid.dx * self.grid.dx
        # Only core nodes (status == 0) need to receive a weight
        wg[np.nonzero(self._grid.status_at_node)] = 0
        wg = rd.rdarray(wg.reshape(self.grid.shape), no_data=-9999)

        a[:] = np.array(
            rd.FlowAccumFromProps(props=props_Pf, weights=wg).reshape(
                self.grid.number_of_nodes
            )
        )

        if any(self.grid.at_node["water__unit_flux_in"] != 1):
            wg = self.grid.at_node["water__unit_flux_in"] * self.grid.dx * self.grid.dx
            # Only core nodes (status == 0) need to receive a weight
            wg[np.nonzero(self._grid.status_at_node)] = 0
            wg = rd.rdarray(wg.reshape(self.grid.shape), no_data=-9999)
            q_pf = rd.FlowAccumFromProps(props=self._props_Pf, weights=wg)
            q[:] = np.array(q_pf.reshape(self.grid.number_of_nodes))
        else:
            q[:] = self._drainage_area

    def update_Hill_FDFA(self, updateDepressions=False):
        if not self._seperate_Hill_Flow:
            raise ValueError(
                "If hillslope properties are updated, the seperate_Hill_Flow property of the FlowDirAccPf class should be True upon initialisation"
            )
        if self._suppress_out:
            with suppress_stdout_stderr():
                sys.stdout = open(os.devnull, "w")
                sys.stderr = open(os.devnull, "w")
                self.calculateFlowDir_Acc(
                    hillFlow=True, updateDepressions=updateDepressions
                )
                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        else:
            self.calculateFlowDir_Acc(
                hillFlow=True, updateDepressions=updateDepressions
            )

    def run_one_step(self):
        if self._suppress_out:
            with suppress_stdout_stderr():
                sys.stdout = open(os.devnull, "w")
                sys.stderr = open(os.devnull, "w")
                self.calculateFlowDir_Acc(
                    hillFlow=False, updateDepressions=self._updateFlowDepressions
                )
                if self._seperate_Hill_Flow:
                    if self._update_HillFlow_instanteneous:
                        self.calculateFlowDir_Acc(
                            hillFlow=True, updateDepressions=self._updateHillDepressions
                        )

                sys.stdout = sys.__stdout__
                sys.stderr = sys.__stderr__
        else:
            self.calculateFlowDir_Acc(
                hillFlow=False, updateDepressions=self._updateFlowDepressions
            )
            if self._seperate_Hill_Flow:
                if self._update_HillFlow_instanteneous:
                    self.calculateFlowDir_Acc(
                        hillFlow=True, updateDepressions=self._updateHillDepressions
                    )


if __name__ == "__main__":  # pragma: no cover
    import doctest

    doctest.testmod()
