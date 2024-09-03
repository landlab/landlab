#! /usr/env/python
"""Fastscape stream power erosion."""

# This module attempts to "component-ify" GT's Fastscape stream
# power erosion.
# Created DEJH, March 2014.


import numpy as np

from landlab import Component
from landlab import RasterModelGrid
from landlab.utils.return_array import return_array_at_node

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import brent_method_erode_fixed_threshold
from .cfuncs import brent_method_erode_variable_threshold


class FastscapeEroder(Component):
    r"""Fastscape stream power erosion.

    This class uses the Braun-Willett Fastscape approach to calculate the
    amount of erosion at each node in a grid, following a stream power
    framework. This should allow it to be stable against larger timesteps
    than an explicit stream power scheme.

    Note that although this scheme is nominally implicit, and will reach a
    numerically-correct solution under topographic steady state regardless of
    timestep length, the accuracy of transient solutions is *not* timestep
    independent (see Braun & Willett 2013, Appendix B for further details).
    Although the scheme remains significantly more robust and permits longer
    timesteps than a traditional explicit solver under such conditions, it
    is still possible to create numerical instability through use of too long
    a timestep while using this component. The user is cautioned to check their
    implementation is behaving stably before fully trusting it.

    Stream power erosion is implemented as:

    .. math::

        E = K  A ^ m  S ^ n -
               \textit{threshold_sp}

    if :math:`K A ^ m S ^ n > \textit{threshold_sp}`, and:

    .. math:: E = 0,

    if :math:`K A^m S^n <= \textit{threshold_sp}`.

    This module assumes you have already run
    :func:`landlab.components.flow_accum.flow_accumulator.FlowAccumulator.run_one_step`
    in the same timestep. It looks for 'flow__upstream_node_order',
    'flow__link_to_receiver_node', 'drainage_area', 'flow__receiver_node', and
    'topographic__elevation' at the nodes in the grid. 'drainage_area' should
    be in area upstream, not volume (i.e., set runoff_rate=1.0 when calling
    FlowAccumulator.run_one_step).

    The primary method of this class is :func:`run_one_step`.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, FastscapeEroder

    >>> grid = RasterModelGrid((5, 5), xy_spacing=10.0)
    >>> z = [
    ...     [7.0, 7.0, 7.0, 7.0, 7.0],
    ...     [7.0, 5.0, 3.2, 6.0, 7.0],
    ...     [7.0, 2.0, 3.0, 5.0, 7.0],
    ...     [7.0, 1.0, 1.9, 4.0, 7.0],
    ...     [7.0, 0.0, 7.0, 7.0, 7.0],
    ... ]
    >>> z = grid.add_field("topographic__elevation", z, at="node")
    >>> fr = FlowAccumulator(grid, flow_director="D8")
    >>> sp = FastscapeEroder(grid, K_sp=1.0)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=1.0)
    >>> z
    array([7.        , 7.        , 7.        , 7.        , 7.        ,
           7.        , 2.92996598, 2.02996598, 4.01498299, 7.        ,
           7.        , 0.85993197, 1.87743897, 3.28268321, 7.        ,
           7.        , 0.28989795, 0.85403051, 2.42701526, 7.        ,
           7.        , 0.        , 7.        , 7.        , 7.        ])

    >>> grid = RasterModelGrid((3, 7), xy_spacing=1.0)
    >>> z = np.array(grid.node_x**2.0)
    >>> z = grid.add_field("topographic__elevation", z, at="node")
    >>> grid.status_at_node[grid.nodes_at_left_edge] = grid.BC_NODE_IS_FIXED_VALUE
    >>> grid.status_at_node[grid.nodes_at_top_edge] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[grid.nodes_at_bottom_edge] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[grid.nodes_at_right_edge] = grid.BC_NODE_IS_CLOSED
    >>> fr = FlowAccumulator(grid, flow_director="D8")
    >>> sp = FastscapeEroder(grid, K_sp=0.1, m_sp=0.0, n_sp=2.0, threshold_sp=2.0)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=10.0)
    >>> z.reshape(grid.shape)[1, :]
    array([ 0.        ,  1.        ,  4.        ,  8.52493781,
           13.29039716, 18.44367965, 36.        ])

    >>> grid = RasterModelGrid((3, 7), xy_spacing=1.0)
    >>> z = np.array(grid.node_x**2.0)
    >>> z = grid.add_field("topographic__elevation", z, at="node")
    >>> grid.status_at_node[grid.nodes_at_left_edge] = grid.BC_NODE_IS_FIXED_VALUE
    >>> grid.status_at_node[grid.nodes_at_top_edge] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[grid.nodes_at_bottom_edge] = grid.BC_NODE_IS_CLOSED
    >>> grid.status_at_node[grid.nodes_at_right_edge] = grid.BC_NODE_IS_CLOSED
    >>> cell_area = 1.0
    >>> fr = FlowAccumulator(grid, flow_director="D8", runoff_rate=2.0)
    >>> grid.at_node["water__unit_flux_in"]
    array([2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
           2., 2., 2., 2., 2., 2., 2., 2.])
    >>> K_field = grid.ones(at="node")  # K can be a field
    >>> sp = FastscapeEroder(
    ...     grid,
    ...     K_sp=K_field,
    ...     m_sp=1.0,
    ...     n_sp=0.6,
    ...     threshold_sp=grid.node_x,
    ...     discharge_field="surface_water__discharge",
    ... )
    >>> fr.run_one_step()
    >>> sp.run_one_step(1.0)
    >>> z.reshape(grid.shape)[1, :]
    array([ 0.        ,  0.0647484 ,  0.58634455,  2.67253503,
            8.49212152, 20.92606987, 36.        ])

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

    _name = "FastscapeEroder"

    _unit_agnostic = True

    _info = {
        "drainage_area": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m**2",
            "mapping": "node",
            "doc": "Upstream accumulated surface area contributing to the node's discharge",
        },
        "flow__link_to_receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "ID of link downstream of each node, which carries the discharge",
        },
        "flow__receiver_node": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array of receivers (node that receives flow from current node)",
        },
        "flow__upstream_node_order": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "Node array containing downstream-to-upstream ordered list of node IDs",
        },
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(
        self,
        grid,
        K_sp=0.001,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=0.0,
        discharge_field="drainage_area",
        erode_flooded_nodes=True,
    ):
        """Initialize the Fastscape stream power component. Note: a timestep,
        dt, can no longer be supplied to this component through the input file.
        It must instead be passed directly to the run method.

        Parameters
        ----------
        grid : ModelGrid
            A grid.
        K_sp : float, array, or field name
            K in the stream power equation (units vary with other parameters).
        m_sp : float, optional
            m in the stream power equation (power on drainage area).
        n_sp : float, optional
            n in the stream power equation (power on slope).
        threshold_sp : float, array, or field name
            Erosion threshold in the stream power equation.
        discharge_field : float, field name, or array, optional
            Discharge [L^2/T]. The default is to use the grid field
            'drainage_area'. To use custom spatially/temporally varying
            rainfall, use 'water__unit_flux_in' to specify water input to the
            FlowAccumulator and use "surface_water__discharge" for this
            keyword argument.
        erode_flooded_nodes : bool (optional)
            Whether erosion occurs in flooded nodes identified by a
            depression/lake mapper (e.g., DepressionFinderAndRouter). When set
            to false, the field *flood_status_code* must be present on the grid
            (this is created by the DepressionFinderAndRouter). Default True.
        """
        super().__init__(grid)

        if "flow__receiver_node" in grid.at_node and grid.at_node[
            "flow__receiver_node"
        ].size != grid.size("node"):
            raise NotImplementedError(
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that FastscapeEroder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )

        if not erode_flooded_nodes and "flood_status_code" not in self._grid.at_node:
            raise ValueError(
                "In order to not erode flooded nodes another component "
                "must create the field *flood_status_code*. You want to "
                "run a lake mapper/depression finder."
            )

        self._erode_flooded_nodes = erode_flooded_nodes

        # use setter for K defined below
        self.K = K_sp

        self._m = float(m_sp)
        self._n = float(n_sp)

        if isinstance(threshold_sp, (float, int)):
            self._thresholds = float(threshold_sp)
        else:
            self._thresholds = return_array_at_node(grid, threshold_sp)

        self._A = return_array_at_node(grid, discharge_field)

        # make storage variables
        self._A_to_the_m = grid.zeros(at="node")
        self._alpha = grid.empty(at="node")

    @property
    def K(self):
        """Erodibility (units depend on m_sp)."""
        return self._K

    @K.setter
    def K(self, new_val):
        self._K = return_array_at_node(self._grid, new_val)

    def run_one_step(self, dt):
        """Erode for a single time step.

        This method implements the stream power erosion across one time
        interval, dt, following the Braun-Willett (2013) implicit Fastscape
        algorithm.

        This follows Landlab standardized component design, and supercedes the
        old driving method :func:`erode`.

        Parameters
        ----------
        dt : float
            Time-step size
        """
        if not self._erode_flooded_nodes:
            flood_status = self._grid.at_node["flood_status_code"]
            flooded_nodes = np.nonzero(flood_status == _FLOODED)[0]
        else:
            flooded_nodes = []

        upstream_order_IDs = self._grid.at_node["flow__upstream_node_order"]
        flow_receivers = self._grid["node"]["flow__receiver_node"]
        z = self._grid.at_node["topographic__elevation"]

        defined_flow_receivers = np.not_equal(
            self._grid.at_node["flow__link_to_receiver_node"], self._grid.BAD_INDEX
        )

        if isinstance(self._grid, RasterModelGrid):
            flow_link_lengths = self._grid.length_of_d8[
                self._grid.at_node["flow__link_to_receiver_node"][
                    defined_flow_receivers
                ]
            ]
        else:
            flow_link_lengths = self._grid.length_of_link[
                self._grid.at_node["flow__link_to_receiver_node"][
                    defined_flow_receivers
                ]
            ]

        np.power(self._A, self._m, out=self._A_to_the_m)
        self._alpha[defined_flow_receivers] = (
            self._K[defined_flow_receivers]
            * dt
            * self._A_to_the_m[defined_flow_receivers]
            / (flow_link_lengths**self._n)
        )

        # Handle flooded nodes, if any (no erosion there)
        if len(flooded_nodes) > 0:
            self._alpha[flooded_nodes] = 0.0
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            self._alpha[reversed_flow] = 0.0

        threshsdt = self._thresholds * dt

        # solve using Brent's Method in Cython for Speed
        if isinstance(self._thresholds, float):
            brent_method_erode_fixed_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, self._alpha, self._n, z
            )
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, self._alpha, self._n, z
            )
