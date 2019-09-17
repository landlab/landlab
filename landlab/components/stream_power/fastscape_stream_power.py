#! /usr/env/python
"""Fastscape stream power erosion."""

# This module attempts to "component-ify" GT's Fastscape stream
# power erosion.
# Created DEJH, March 2014.


import numpy as np

from landlab import BAD_INDEX_VALUE as UNDEFINED_INDEX, Component, RasterModelGrid

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import (
    brent_method_erode_fixed_threshold,
    brent_method_erode_variable_threshold,
)


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

        E = K (\textit{rainfall_intensity} \, A) ^ m  S ^ n -
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
    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> from landlab.components import FlowAccumulator, FastscapeEroder

    >>> grid = RasterModelGrid((5, 5), xy_spacing=10.)
    >>> z = np.array([7.,  7.,  7.,  7.,  7.,
    ...               7.,  5., 3.2,  6.,  7.,
    ...               7.,  2.,  3.,  5.,  7.,
    ...               7.,  1., 1.9,  4.,  7.,
    ...               7.,  0.,  7.,  7.,  7.])
    >>> z = grid.add_field('topographic__elevation', z, at='node')
    >>> fr = FlowAccumulator(grid, flow_director='D8')
    >>> sp = FastscapeEroder(grid, K_sp=1.)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=1.)
    >>> z  # doctest: +NORMALIZE_WHITESPACE
    array([ 7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
            7.        ,  2.92996598,  2.02996598,  4.01498299,  7.        ,
            7.        ,  0.85993197,  1.87743897,  3.28268321,  7.        ,
            7.        ,  0.28989795,  0.85403051,  2.42701526,  7.        ,
            7.        ,  0.        ,  7.        ,  7.        ,  7.        ])

    >>> grid = RasterModelGrid((3, 7), xy_spacing=1.)
    >>> z = np.array(grid.node_x ** 2.)
    >>> z = grid.add_field('topographic__elevation', z, at='node')
    >>> grid.status_at_node[grid.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> fr = FlowAccumulator(grid, flow_director='D8')
    >>> sp = FastscapeEroder(grid, K_sp=0.1, m_sp=0., n_sp=2.,
    ...                      threshold_sp=2.)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=10.)
    >>> z.reshape(grid.shape)[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   1.        ,   4.        ,   8.52493781,
            13.29039716,  18.44367965,  36.        ])

    >>> grid = RasterModelGrid((3, 7), xy_spacing=1.)
    >>> z = np.array(grid.node_x ** 2.)
    >>> z = grid.add_field('topographic__elevation', z, at='node')
    >>> grid.status_at_node[grid.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> grid.status_at_node[grid.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> fr = FlowAccumulator(grid, flow_director='D8')
    >>> K_field = grid.ones(at='node') # K can be a field
    >>> sp = FastscapeEroder(grid, K_sp=K_field, m_sp=1., n_sp=0.6,
    ...                      threshold_sp=grid.node_x,
    ...                      rainfall_intensity=2.)
    >>> fr.run_one_step()
    >>> sp.run_one_step(1.)
    >>> z.reshape(grid.shape)[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   0.0647484 ,   0.58634455,   2.67253503,
             8.49212152,  20.92606987,  36.        ])
    """

    _name = "FastscapeEroder"

    _input_var_names = set(
        (
            "topographic__elevation",
            "drainage_area",
            "flow__link_to_receiver_node",
            "flow__upstream_node_order",
            "flow__receiver_node",
        )
    )

    _output_var_names = set(("topographic__elevation",))

    _var_units = {
        "topographic__elevation": "m",
        "drainage_area": "m**2",
        "flow__link_to_receiver_node": "-",
        "flow__upstream_node_order": "-",
        "flow__receiver_node": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "drainage_area": "node",
        "flow__link_to_receiver_node": "node",
        "flow__upstream_node_order": "node",
        "flow__receiver_node": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "drainage_area": "Upstream accumulated surface area contributing to the node's "
        "discharge",
        "flow__link_to_receiver_node": "ID of link downstream of each node, which carries the discharge",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
    }

    def __init__(
        self,
        grid,
        K_sp=None,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=0.0,
        rainfall_intensity=1.0,
        discharge_field="drainage_area",
        erode_flooded_nodes=True,
    ):
        """
        Initialize the Fastscape stream power component. Note: a timestep,
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
        rainfall intensity : float, array, or field name; optional
            Modifying factor on drainage area to convert it to a true water
            volume flux in (m/time). i.e., E = K * (r_i*A)**m * S**n
        discharge_field : string; optional
            Name of field to use for discharge proxy. Defaults to 'drainage_area',
            which means the component will expect the driver or another component
            to have created and populated a 'drainage_area' field. To use a
            different field, such as 'surface_water__discharge', give its name in
            this argument.
        erode_flooded_nodes : bool (optional)
            Whether erosion occurs in flooded nodes identified by a
            depression/lake mapper (e.g., DepressionFinderAndRouter). When set
            to false, the field *flood_status_code* must be present on the grid
            (this is created by the DepressionFinderAndRouter). Default True.
        """
        super(FastscapeEroder, self).__init__(grid)

        if "flow__receiver_node" in grid.at_node:
            if grid.at_node["flow__receiver_node"].size != grid.size("node"):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The landlab development team has not "
                    "verified that FastscapeEroder is compatible with "
                    "route-to-multiple methods. Please open a GitHub Issue "
                    "to start this process."
                )
                raise NotImplementedError(msg)

        if not erode_flooded_nodes:
            if "flood_status_code" not in self._grid.at_node:
                msg = (
                    "In order to not erode flooded nodes another component "
                    "must create the field *flood_status_code*. You want to "
                    "run a lake mapper/depression finder."
                )
                raise ValueError(msg)

        self._erode_flooded_nodes = erode_flooded_nodes

        self._K = K_sp  # overwritten below in special cases
        self._m = float(m_sp)
        self._n = float(n_sp)
        if isinstance(threshold_sp, (float, int)):
            self._thresholds = float(threshold_sp)
        else:
            if isinstance(threshold_sp, str):
                self._thresholds = self._grid.at_node[threshold_sp]
            else:
                self._thresholds = threshold_sp
            assert self._thresholds.size == self._grid.number_of_nodes

        # make storage variables
        self._A_to_the_m = grid.zeros(at="node")
        self._alpha = grid.empty(at="node")

        if self._K is None:
            raise ValueError(
                "K_sp must be set as a float, node array, or "
                + "field name. It was None."
            )

        # now handle the inputs that could be float, array or field name:
        # some support here for old-style inputs
        if isinstance(K_sp, str):
            if K_sp == "array":
                self._K = None
            else:
                self._K = self._grid.at_node[K_sp]
        elif isinstance(K_sp, (float, int)):
            self._K = float(K_sp)
        else:
            self._K = np.asarray(K_sp, dtype=float)
            if len(self._K) != self._grid.number_of_nodes:
                raise TypeError("Supplied value of K_sp is not n_nodes long")

        if isinstance(rainfall_intensity, str):
            raise ValueError(
                "This component can no longer handle "
                + "spatially variable runoff directly. Use "
                + "FlowAccumulator with specified "
                + "water__unit_flux_in, or use StreamPowerEroder"
                + "component instead of FastscapeEroder."
            )
            if rainfall_intensity == "array":
                self._r_i = None
            else:
                self._r_i = self._grid.at_node[rainfall_intensity]
        elif isinstance(rainfall_intensity, (float, int)):  # a float
            self._r_i = float(rainfall_intensity)
        elif len(rainfall_intensity) == self._grid.number_of_nodes:
            raise ValueError(
                "This component can no longer handle "
                "spatially variable runoff directly. Use "
                "FlowAccumulator with specified "
                "water__unit_flux_in, or use StreamPowerEroder"
                "component instead of FastscapeEroder."
            )
            self._r_i = np.array(rainfall_intensity)
        else:
            raise TypeError("Supplied type of rainfall_intensity was " "not recognised")

        # Handle option for area vs discharge
        self._discharge_field = discharge_field
        self._verify_output_fields()

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
            self._grid.at_node["flow__link_to_receiver_node"], UNDEFINED_INDEX
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
        # make arrays from input the right size
        if isinstance(self._K, np.ndarray):
            K_here = self._K[defined_flow_receivers]
        else:
            K_here = self._K

        r_i_here = self._r_i

        if self._K is None:  # "old style" setting of array
            assert K_if_used is not None
            self._K = K_if_used

        n = float(self._n)

        np.power(
            self._grid["node"][self._discharge_field], self._m, out=self._A_to_the_m
        )
        self._alpha[defined_flow_receivers] = (
            r_i_here ** self._m
            * K_here
            * dt
            * self._A_to_the_m[defined_flow_receivers]
            / (flow_link_lengths ** self._n)
        )

        alpha = self._alpha

        # Handle flooded nodes, if any (no erosion there)
        if len(flooded_nodes) > 0:
            alpha[flooded_nodes] = 0.0
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            alpha[reversed_flow] = 0.0

        threshsdt = self._thresholds * dt

        # solve using Brent's Method in Cython for Speed
        if isinstance(self._thresholds, float):
            brent_method_erode_fixed_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z
            )
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z
            )
