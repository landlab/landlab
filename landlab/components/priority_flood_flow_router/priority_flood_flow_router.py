"""Fill or breach a DEM, accumulate flow and calculate drainage area using
the priority flood algorithm.

PriorityFloodFlowRouter is a wrapper of the RichDEM package:
https://richdem.readthedocs.io/en/latest/flow_metrics.html

The component merges a filling/breaching algorithm, a flow director as well
as a flow accumulator.  Moreover, the component supports the definition of
two flow accumulator fields associated to the same grid.  This prevents the
user from updating the filling/breaching algorithms in between calculation
of flow accumulator one and two.

@author: benjaminCampforts
"""

import copy as cp
from functools import partial

import numpy as np
import numpy.matlib as npm

from landlab import Component
from landlab import FieldError
from landlab import RasterModelGrid
from landlab.grid.nodestatus import NodeStatus
from landlab.utils.return_array import return_array_at_node

from ...utils.suppress_output import suppress_output
from .cfuncs import _D8_FlowAcc
from .cfuncs import _D8_flowDir

# Codes for depression status
_UNFLOODED = 0
_PIT = 1
_CURRENT_LAKE = 2
_FLOODED = 3
# Flow metrics resulting in single flow
PSINGLE_FMs = ["D8", "D4", "Rho8", "Rho4"]
# Flow metrics resulting in multiple flow
PMULTIPLE_FMs = ["Quinn", "Freeman", "Holmgren", "Dinf"]


class PriorityFloodFlowRouter(Component):
    """Component to accumulate flow and calculate drainage area based RICHDEM software package.

    See also: https://richdem.readthedocs.io/en/latest/


    .. note::

        The perimeter nodes  NEVER contribute to the accumulating flux, even
        if the  gradients from them point inwards to the main body of the grid.
        This is because under Landlab definitions, perimeter nodes lack cells, so
        cannot accumulate any discharge.

    *FlowAccumulatorPf* stores as *ModelGrid* fields:

    - *'drainage_area'*: Node array of drainage areas
    - *'flood_status_code'*: Map of flood status (_PIT, _CURRENT_LAKE, _UNFLOODED, or _FLOODED).
    - *'surface_water__discharge'*: Node array of discharges.
    - *'Distance to receiver'*: Distance to receiver
    - *'water__unit_flux_in'*: External volume water per area per time input to each node.
    - *'flow__upstream_node_order'*: Node array containing downstream-to-upstream ordered
      list of node IDs.
    - *'flow__receiver_node'*: Node array of receivers (nodes that receive flow),
      or ITS OWN ID if there is no receiver. This array is 2D for *RouteToMany*
      methods and has the shape *(n-nodes x max number of receivers)*.
    - *'flow__receiver_proportions'*: Node array of flow proportions. This
      array is 2D, for *RouteToMany* methods and has the shape
      *(n-nodes x max number of receivers)*.
    - *'topographic__steepest_slope'*: Node array of downhill slopes from each receiver.
      This array is 2D for *RouteToMany* methods and has the shape
      *(n-nodes x max number of receivers)*.
    - *'flow__link_to_receiver_node'*: Node array of links carrying flow.
    - *'flow__receiver_proportion's*: Node array of proportion of flow sent to each receiver.
    - *'depression_free_elevation'*: Depression free land surface topographic
      elevation, at closed borders, value equals -1.

    The following fields are required when an additional hillslope flowrouting
    scheme is required, can be completed with flow acc and discharge if required:

    - *'hill_flow__upstream_node_order'*: Node array containing downstream-to-upstream
      ordered list of node IDs
    - *'hill_flow__receiver_node'*: Node array of receivers (node that receives flow
      from current node)
    - *'hill_topographic__steepest_slope'*: The steepest *downhill* slope.
    - *'hill_flow__receiver_proportions'*: Node array of proportion of flow sent to each
      receiver

    The primary method of this class is :func:`run_one_step`.

    Parameters
    ----------
    grid : ModelGrid
        A Landlab grid.
    surface : str or array_like, optional
        The surface to direct flow across. An at-node field name or an array
        of length *n_node*.
    flow_metric : str, optional
        String has to be one of 'D8' (O’Callaghan and Mark, 1984), 'Rho8'
        (Fairfield and Leymarie, 1991), 'Quinn' (1991), 'Freeman' (1991),
        'Holmgren' (1994), 'Dinf' (Tarboton, 1997). For details and comparison,
        see https://richdem.readthedocs.io/en/latest/flow_metrics.html
    runoff_rate : str, array_like, or float, optional
        If provided, sets the runoff rate (m / time) and will be assigned to the grid field
        'water__unit_flux_in'. If a spatially and and temporally variable runoff
        rate is desired, pass this field name and update the field through model
        run time. If both the field and argument are present at the time of
        initialization, runoff_rate will *overwrite* the field. If neither are
        set, defaults to spatially constant unit input.
        Both a runoff_rate array and the 'water__unit_flux_in' field are
        permitted to contain negative values, in which case they mimic
        transmission losses rather than e.g. rain inputs.
    update_flow_depressions : bool, optional
        Build-in depression handler. Can be through filling or breaching (see below).
    update_hill_depressions : bool, optional
        Only needed if DEM needs to be filled separately for second (hill flow)
        flow accumulator.  Default behavior is not to execute a separate filling
        procedure in between the first and the second flow accumulator.
    depression_handler : str, optional
        Must be one of 'fill or 'breach'.
        Depression-Filling or breaching algorithm to process depressions

        - 'fill': Depression-Filling.
          Depression-filling is often used to fill in all the depressions
          in a DEM to the level of their lowest outlet or spill-point.
          See also: https://richdem.readthedocs.io/en/latest/depression_filling.html
        - 'breach': Complete Breaching.
          Depression-breaching is used to dig channels from the pit cells
          of a DEM to the nearest cells (in priority-flood sense) outside
          of the depression containing the pit. This resolves the depression
          as all cells in the depression now have a drainage path to the
          edge of the DEM.
          See also: https://richdem.readthedocs.io/en/latest/depression_breaching.html
    exponent : float, optional
        Some methods require an exponent (see flow_metric) Default {1}
    epsilon : bool, optional
        If ``True``, an epsilon gradient is imposed to all flat regions. This ensures
        that there is always a local gradient.
    accumulate_flow : bool, optional
        If ``True`` flow directions and accumulations will be calculated.
        Set to ``False`` when only interested in flow directions
    accumulate_flow_hill : bool, optional
        If ``True`` flow directions and accumulations will be calculated
        for second FD component (Hill). Set to ``False`` when only interested in flow
        directions.
    separate_hill_flow : bool, optional
        For some applications (e.g. *HyLands*) both single and
        multiple flow direction and accumulation is required.
        By calculating them in the same component, we can optimize procedures
        involved with filling and breaching of DEMs
    update_hill_flow_instantaneous : bool, optional
        Update separate hillslope director and accumulator simultaneously on update.
        Set if other operations have to be performed in between updating the
        principle flow properties and the hillslope properties.
    hill_flow_metric : str, optional
        Must be one 'D8' (O’Callaghan and Mark, 1984),'D4' (O’Callaghan and Mark, 1984),
        'Rho8' (Fairfield and Leymarie, 1991), 'Rho4' (Fairfield and Leymarie, 1991),
        'Quinn' (1991) {default},'Freeman' (1991), 'Holmgren' (1994),
        'Dinf' (Tarboton, 1997).
        For details and comparison, see
        https://richdem.readthedocs.io/en/latest/flow_metrics.html
    hill_exponent : float, optional
        Some methods require an exponent (see flow_metric)
    suppress_out : bool, optional
        Suppress verbose of priority flood algorithm


    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Barnes, R., 2017. Parallel non-divergent flow accumulation for trillion
    cell digital elevation models on desktops or clusters. Environmental
    Modelling & Software 92, 202–212. doi: 10.1016/j.envsoft.2017.02.022

    **Additional References**

    https://richdem.readthedocs.io/en/latest/

    """

    _name = "PriorityFloodFlowRouter"

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
            "doc": (
                "External volume water per area per time input to each node "
                "(e.g., rainfall rate)"
            ),
        },
        "topographic__steepest_slope": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "The steepest *downhill* slope",
        },
        "squared_length_adjacent": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": (
                "Length to adjacent nodes, squared (calcualted in advance to "
                "save time during calculation"
            ),
        },
        "flow__receiver_proportions": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of proportion of flow sent to each receiver.",
        },
        "depression_free_elevation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "Filled land surface topographic elevation, at closed borders, "
                "value equals -1!"
            ),
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
        update_flow_depressions=True,
        depression_handler="fill",
        exponent=1,
        epsilon=True,
        accumulate_flow=True,
        accumulate_flow_hill=False,
        separate_hill_flow=False,
        update_hill_depressions=False,
        update_hill_flow_instantaneous=True,
        hill_flow_metric="Quinn",
        hill_exponent=1,
        suppress_out=True,
    ):
        """Initialize the FlowAccumulator component.

        Saves the grid, tests grid type, tests input types and
        compatibility for the flow_metric and depression_finder
        keyword arguments, tests the argument of runoff, and
        initializes new fields.
        """
        self._richdem = self.load_richdem()

        super().__init__(grid)
        # Keep a local reference to the grid

        self._suppress_output = partial(
            suppress_output, out=suppress_out, err=suppress_out
        )

        # STEP 1: Testing of input values, supplied either in function call or
        # as part of the grid.
        self._validate_water_inputs(grid, runoff_rate)

        # Grid type testing
        if not isinstance(self._grid, RasterModelGrid):
            raise FieldError(
                "Flow Accumulator Priority flood only works with regular raster grids, "
                "use default Landlab flow accumulator instead"
            )

        # save elevations as class properites.
        self._surface = surface
        self._surface_values = return_array_at_node(grid, surface)

        node_cell_area = self._grid.cell_area_at_node.copy()
        node_cell_area[self._grid.closed_boundary_nodes] = 0.0
        self._node_cell_area = node_cell_area
        self._runoff_rate = runoff_rate

        if (flow_metric in PSINGLE_FMs) or (flow_metric in PMULTIPLE_FMs):
            self._flow_metric = flow_metric
        else:
            raise ValueError(
                "flow metric should be one of these single flow directors: "
                f"{', '.join(repr(x) for x in PSINGLE_FMs)} or multiple flow directors: "
                f"{', '.join(repr(x) for x in PMULTIPLE_FMs)}"
            )
        if (hill_flow_metric in PSINGLE_FMs) or (hill_flow_metric in PMULTIPLE_FMs):
            self._hill_flow_metric = hill_flow_metric
        else:
            raise ValueError(
                "flow metric should be one of these single flow directors:"
                f"{', '.join(repr(x) for x in PSINGLE_FMs)} or multiple flow directors: "
                f"{', '.join(repr(x) for x in PMULTIPLE_FMs)}"
            )

        if depression_handler in ("fill", "breach"):
            self._depression_handler = depression_handler
        else:
            raise ValueError(
                "depression_handler should be one of 'fill' or 'breach'"
                f" (got {depression_handler!r})"
            )

        self._epsilon = epsilon
        self._exponent = exponent
        self._separate_hill_flow = separate_hill_flow
        self._update_hill_flow_instantaneous = update_hill_flow_instantaneous

        self._update_flow_depressions = update_flow_depressions
        self._update_hill_depressions = update_hill_depressions

        self._hill_exponent = hill_exponent

        if self._separate_hill_flow:
            # Adjust dict
            self._info["hill_drainage_area"]["optional"] = False
            self._info["hill_surface_water__discharge"]["optional"] = False
            self._info["hill_flow__upstream_node_order"]["optional"] = False
            self._info["hill_flow__receiver_node"]["optional"] = False
            self._info["hill_topographic__steepest_slope"]["optional"] = False
            self._info["hill_flow__receiver_proportions"]["optional"] = False
        else:
            self._info["hill_drainage_area"]["optional"] = True
            self._info["hill_surface_water__discharge"]["optional"] = True
            self._info["hill_flow__upstream_node_order"]["optional"] = True
            self._info["hill_flow__receiver_node"]["optional"] = True
            self._info["hill_topographic__steepest_slope"]["optional"] = True
            self._info["hill_flow__receiver_proportions"]["optional"] = True

        self._accumulate_flow = accumulate_flow
        self._accumulate_flow_hill = accumulate_flow_hill

        if not self._accumulate_flow:
            self._info["drainage_area"]["optional"] = True
            self._info["surface_water__discharge"]["optional"] = True
        else:
            self._info["drainage_area"]["optional"] = False
            self._info["surface_water__discharge"]["optional"] = False

        if not self._accumulate_flow_hill:
            self._info["hill_drainage_area"]["optional"] = True
            self._info["hill_surface_water__discharge"]["optional"] = True
        else:
            self._info["hill_drainage_area"]["optional"] = False
            self._info["hill_surface_water__discharge"]["optional"] = False

        self.initialize_output_fields()

        # Make aliases
        if self._accumulate_flow:
            self._drainage_area = self.grid.at_node["drainage_area"]
            self._discharges = self.grid.at_node["surface_water__discharge"]
        self._sort = self.grid.at_node["flow__upstream_node_order"]
        # if multiple flow algorithm is made, the dimensions of the slope
        # and receiver fields change (8 colums for all neightbors)
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

        if self._separate_hill_flow:
            if self._accumulate_flow_hill:
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
        self._create_richdem_properties()

    @staticmethod
    def load_richdem():
        try:
            import _richdem  # noqa: F401
            import richdem
        except ModuleNotFoundError as exc:
            raise ModuleNotFoundError(
                "PriorityFloodFlowRouter requires richdem but richdem is not installed. "
                "You can install richdem either from source "
                "(https://github.com/r-barnes/richdem), or through conda "
                "(conda install richdem -c conda-forge) or pip (pip install richdem)."
            ) from exc
        return richdem

    @property
    def surface_values(self):
        """Values of the surface over which flow is directed."""
        return self._surface_values

    def _changed_surface(self):
        """Check if the surface values have changed.

        If the surface values are stored as a field, it is important to
        check if they have changed since the component was instantiated.
        """
        if isinstance(self._surface, str):
            self._surface_values = return_array_at_node(self._grid, self._surface)

    @property
    def node_drainage_area(self):
        """Return the drainage area."""
        return self._grid["node"]["drainage_area"]

    @property
    def node_water_discharge(self):
        """Return the surface water discharge."""
        return self._grid["node"]["surface_water__discharge"]

    def _create_richdem_properties(self):
        self._depression_free_dem = cp.deepcopy(
            self._richdem.rdarray(
                self._surface_values.reshape(self.grid.shape),
                no_data=-9999,
            )
        )
        self._depression_free_dem.geotransform = [0, 1, 0, 0, 0, -1]

        # Calculate SQUARED length adjacent
        self.grid.at_node["squared_length_adjacent"] = np.concatenate(
            (
                np.ones((self.grid.number_of_nodes, 4)),
                2 * np.ones((self.grid.number_of_nodes, 4)),
            ),
            axis=1,
        )

        self._closed = np.zeros(self._grid.number_of_nodes, dtype=np.uint8)
        self._closed[self._grid.status_at_node == NodeStatus.CLOSED] = 1
        self._closed = self._richdem.rdarray(
            self._closed.reshape(self._grid.shape), no_data=-9999
        )
        self._closed.geotransform = [0, 1, 0, 0, 0, -1]

    def _validate_water_inputs(self, grid, runoff_rate):
        """Test inputs for runoff_rate and water__unit_flux_in."""

        if "water__unit_flux_in" not in grid.at_node:
            grid.add_empty("water__unit_flux_in", at="node")
            runoff_rate = 1.0 if runoff_rate is None else runoff_rate
            grid.at_node["water__unit_flux_in"][:] = runoff_rate
        elif runoff_rate is not None:
            grid.at_node["water__unit_flux_in"][:] = runoff_rate

    def calc_flow_dir_acc(self, hill_flow=False, update_depressions=True):
        """Calculate flow direction and accumulation using the richdem package"""
        if hill_flow:
            flow_metric = self._hill_flow_metric
        else:
            flow_metric = self._flow_metric

        # 1: Remove depressions
        if update_depressions:
            self.remove_depressions(flow_metric=flow_metric)
        if flow_metric == "D8":
            self._FlowAcc_D8(hill_flow=hill_flow)
        else:
            closed_boundary_values = self._depression_free_dem[self._closed == 1]
            self._depression_free_dem[self._closed == 1] = np.inf
            # Calculate flow direction (proportion) and accumulation using RichDEM
            with self._suppress_output():
                props_Pf = self._richdem.FlowProportions(
                    dem=self._depression_free_dem,
                    method=flow_metric,
                    exponent=self._exponent,
                )
            self._depression_free_dem[self._closed == 1] = closed_boundary_values

            # Calculate flow accumulation using RichDEM
            if (hill_flow and self._accumulate_flow_hill) or (
                not hill_flow and self._accumulate_flow
            ):
                self._accumulate_flow_RD(props_Pf, hill_flow=hill_flow)

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

            slope_temp = (
                self._surface_values.reshape(-1, 1) - self._surface_values[rcvrs]
            ) / (self.grid.dx * np.sqrt(self.grid.at_node["squared_length_adjacent"]))

            if flow_metric in PSINGLE_FMs:
                slope_temp[rcvrs == -1] = 0
            else:
                props_Pf[props_Pf_col0 == -1, 0] = 1
                props_Pf = props_Pf.astype(np.float64)  # should be float64
                # Now, make sure sum is 1 in 64 bits
                props_Pf[props_Pf == -1] = 0
                proportion_matrix = npm.repmat(
                    np.reshape(props_Pf.sum(axis=1), [props_Pf.shape[0], 1]), 1, 8
                )
                rc64_temp = np.where(
                    proportion_matrix == 0, props_Pf, props_Pf / proportion_matrix
                )
                props_Pf[props_Pf[:, 0] != 1, :] = rc64_temp[props_Pf[:, 0] != 1, :]
                props_Pf[props_Pf == 0] = -1

            if hill_flow:
                if flow_metric in PSINGLE_FMs:
                    ij_at_max = range(len(rcvrs)), np.argmax(rcvrs, axis=1)
                    self._hill_prps[:] = props_Pf[ij_at_max]
                    self._hill_rcvs[:] = rcvrs[ij_at_max]
                    self._hill_slope[:] = slope_temp[ij_at_max]
                else:
                    self._hill_prps[:] = props_Pf
                    self._hill_rcvs[:] = rcvrs
                    self._hill_slope[:] = slope_temp
                    self._hill_slope[rcvrs == -1] = 0

            else:
                if flow_metric in PSINGLE_FMs:
                    ij_at_max = range(len(rcvrs)), np.argmax(rcvrs, axis=1)
                    self._prps[:] = props_Pf[ij_at_max]
                    self._rcvs[:] = rcvrs[ij_at_max]
                    self._slope[:] = slope_temp[ij_at_max]
                    self._recvr_link[:] = recvr_link[ij_at_max]
                else:
                    self._prps[:] = props_Pf
                    self._rcvs[:] = rcvrs
                    self._slope[:] = slope_temp
                    self._slope[rcvrs == -1] = 0
                    self._recvr_link[:] = recvr_link

    def _FlowAcc_D8(self, hill_flow=False):
        """
        Function to calcualte flow accumulation using the D8 flow algorithm.

        Parameters
        ----------
        hill_flow : Boolean, optional
            Defines which instance of flow accumulation is updated.
            If FALSE, the first, default instance is updated.
            If TRUE, the second, hillslope, instance is updated.
            The default is False.

        Returns
        -------
        None.

        """
        c = self.grid.number_of_node_columns
        dx = self.grid.dx
        activeCells = np.array(
            self._grid.status_at_node != NodeStatus.CLOSED + 0, dtype=int
        )
        receivers = np.array(self.grid.status_at_node, dtype=int)
        distance_receiver = np.zeros((receivers.shape), dtype=float)
        cores = self.grid.core_nodes
        activeCores = cores[activeCells[cores] == 1]
        # Make boundaries to save time with conditionals in c loops
        receivers[np.nonzero(self._grid.status_at_node)] = -1
        steepest_slope = np.zeros((receivers.shape), dtype=float)
        el_dep_free = self._depression_free_dem.reshape(self.grid.number_of_nodes)
        el_ori = self._surface_values
        dist = np.multiply(
            [1, 1, 1, 1, np.sqrt(2), np.sqrt(2), np.sqrt(2), np.sqrt(2)], dx
        )
        ngb = np.zeros((8,), dtype=int)
        el_d = np.zeros((8,), dtype=float)

        # Links
        adj_link = np.array(self._grid.d8s_at_node, dtype=int)
        recvr_link = np.zeros((receivers.shape), dtype=int) - 1

        _D8_flowDir(
            receivers,
            distance_receiver,
            steepest_slope,
            np.array(el_dep_free),
            el_ori,
            dist,
            ngb,
            activeCores,
            activeCells,
            el_d,
            c,
            dx,
            adj_link,
            recvr_link,
        )

        # Calcualte flow acc
        do_FA = False
        if hill_flow:
            if self._accumulate_flow_hill:
                do_FA = True
                a = self._hill_drainage_area
                q = self._hill_discharges
        else:
            if self._accumulate_flow:
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
                wg_q[np.nonzero(self._grid.status_at_node)] = NodeStatus.CORE
                dis = wg_q

            else:
                dis = np.full(self.grid.number_of_nodes, self._node_cell_area)

            da = np.array(self._node_cell_area)
            stack_flip = np.flip(self._sort)
            # Filter out donors giving to receivers being -1
            stack_flip = stack_flip[receivers[stack_flip] != -1]

            _D8_FlowAcc(da, dis, stack_flip, receivers)

            a[:] = da
            q[:] = dis

        # Closed nodes flow to themselves
        val = np.arange(0, receivers.shape[0])
        receivers[receivers == -1] = val[receivers == -1]

        # Restore depression free DEM
        # self._depression_free_dem[self._closed == 1] = -1

        if hill_flow:
            self._hill_prps[:] = 1
            self._hill_rcvs[:] = receivers
            self._hill_slope[:] = steepest_slope
        else:
            self._prps[:] = 1
            self._rcvs[:] = receivers
            self._slope[:] = steepest_slope
            self._recvr_link[:] = recvr_link

    def remove_depressions(self, flow_metric="D8"):
        self._depression_free_dem = cp.deepcopy(
            self._richdem.rdarray(
                self._surface_values.reshape(self.grid.shape),
                no_data=-9999,
            )
        )
        self._depression_free_dem.geotransform = [0, 1, 0, 0, 0, -1]
        closed_boundary_values = self._depression_free_dem[self._closed == 1]
        self._depression_free_dem[self._closed == 1] = np.inf

        if flow_metric in ("D4", "Rho4"):
            topology = "D4"
        else:
            topology = "D8"
        with self._suppress_output():
            if self._depression_handler == "fill":
                self._richdem.FillDepressions(
                    self._depression_free_dem,
                    epsilon=self._epsilon,
                    in_place=True,
                    topology=topology,
                )

            elif self._depression_handler == "breach":
                self._richdem.BreachDepressions(
                    self._depression_free_dem, in_place=True, topology=topology
                )

        self._sort[:] = np.argsort(
            np.array(self._depression_free_dem.reshape(self.grid.number_of_nodes))
        )

        self._depression_free_dem[self._closed == 1] = closed_boundary_values
        self.grid.at_node["depression_free_elevation"] = self._depression_free_dem

        self.grid.at_node["flood_status_code"] = np.where(
            self.grid.at_node["depression_free_elevation"]
            == self.grid.at_node["topographic__elevation"],
            0,
            3,
        )

    def _accumulate_flow_RD(self, props_Pf, hill_flow=False):
        """
        Function to accumualte flow using the richdem package

        Parameters
        ----------
        props_Pf : float
            flow proportions calcualte with the RichDEM package using the
            FlowProportions function
        hill_flow : Boolean, optional
            Defines which instance of flow accumulation is updated.
            If FALSE, the first, default instance is updated.
            If TRUE, the second, hillslope, instance is updated.
            The default is False.

        Returns
        -------
        None.

        """
        if not hill_flow:
            a = self._drainage_area
            q = self._discharges
        else:
            a = self._hill_drainage_area
            q = self._hill_discharges

        # Create weight for flow accum: both open (status ==1) and closed
        # nodes (status ==4) will have zero weight
        wg = np.full(self.grid.number_of_nodes, self.grid.dx**2)

        # Only core nodes (status == 0) need to receive a weight
        wg[self._grid.status_at_node != NodeStatus.CORE] = 0
        wg = self._richdem.rdarray(
            wg.reshape(self.grid.shape),
            no_data=-9999,
        )
        wg.geotransform = [0, 1, 0, 0, 0, -1]

        with self._suppress_output():
            a[:] = np.array(
                self._richdem.FlowAccumFromProps(props=props_Pf, weights=wg).reshape(
                    self.grid.number_of_nodes
                )
            )

        if any(self.grid.at_node["water__unit_flux_in"] != 1):
            wg = self.grid.at_node["water__unit_flux_in"] * self.grid.dx * self.grid.dx
            # Only core nodes (status == 0) need to receive a weight
            wg[self._grid.status_at_node != NodeStatus.CORE] = 0
            wg = self._richdem.rdarray(wg.reshape(self.grid.shape), no_data=-9999)
            wg.geotransform = [0, 1, 0, 0, 0, -1]
            with self._suppress_output():
                q_pf = self._richdem.FlowAccumFromProps(props=props_Pf, weights=wg)
            q[:] = np.array(q_pf.reshape(self.grid.number_of_nodes))
        else:
            q[:] = self._drainage_area

    def update_hill_fdfa(self, update_depressions=False):
        if not self._separate_hill_flow:
            raise ValueError(
                "If hillslope properties are updated, the separate_hill_flow "
                "property of the PriorityFloodFlowRouter class should be "
                "True upon initialisation"
            )
        self.calc_flow_dir_acc(hill_flow=True, update_depressions=update_depressions)

    def run_one_step(self):
        self.calc_flow_dir_acc(
            hill_flow=False, update_depressions=self._update_flow_depressions
        )
        if self._separate_hill_flow and self._update_hill_flow_instantaneous:
            self.calc_flow_dir_acc(
                hill_flow=True, update_depressions=self._update_hill_depressions
            )
