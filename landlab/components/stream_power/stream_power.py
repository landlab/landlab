import numpy as np

from landlab import Component
from landlab import MissingKeyError
from landlab.utils.return_array import return_array_at_node

from ..depression_finder.lake_mapper import _FLOODED
from .cfuncs import brent_method_erode_fixed_threshold
from .cfuncs import brent_method_erode_variable_threshold


class StreamPowerEroder(Component):
    """Erode where channels are.

    Implemented as:

    .. math::
        E = K A^m S^n - sp_{crit},

    and if :math:`E < 0`, :math:`E = 0`.

    If ``channel_width_field`` is declared and ``True``, the module instead implements:

    .. math::
        E = K A^m S^n / W - sp_{crit}

    Note that although the Braun-Willett (2013) scheme that underlies this
    component is nominally implicit, and will reach a numerically-correct
    solution under topographic steady state regardless of timestep length, the
    accuracy of transient solutions is *not* timestep independent (see
    Braun & Willett 2013, Appendix B for further details).
    Although the scheme remains significantly more robust and permits longer
    timesteps than a traditional explicit solver under such conditions, it
    is still possible to create numerical instability through use of too long
    a timestep while using this component. The user is cautioned to check their
    implementation is behaving stably before fully trusting it.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import FlowAccumulator, StreamPowerEroder

    >>> mg = RasterModelGrid((5, 5), xy_spacing=10.0)
    >>> z = np.array(
    ...     [
    ...         [7.0, 7.0, 7.0, 7.0, 7.0],
    ...         [7.0, 5.0, 3.2, 6.0, 7.0],
    ...         [7.0, 2.0, 3.0, 5.0, 7.0],
    ...         [7.0, 1.0, 1.9, 4.0, 7.0],
    ...         [7.0, 0.0, 7.0, 7.0, 7.0],
    ...     ]
    ... )
    >>> z = mg.add_field("topographic__elevation", z, at="node")
    >>> fr = FlowAccumulator(mg, flow_director="D8")
    >>> sp = StreamPowerEroder(mg, K_sp=1.0)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=1.0)
    >>> z
    array([7.        , 7.        , 7.        , 7.        , 7.        ,
           7.        , 2.92996598, 2.02996598, 4.01498299, 7.        ,
           7.        , 0.85993197, 1.87743897, 3.28268321, 7.        ,
           7.        , 0.28989795, 0.85403051, 2.42701526, 7.        ,
           7.        , 0.        , 7.        , 7.        , 7.        ])

    >>> mg2 = RasterModelGrid((3, 7))
    >>> z = np.array(mg2.node_x**2.0)
    >>> z = mg2.add_field("topographic__elevation", z, at="node")
    >>> mg2.status_at_node[mg2.nodes_at_left_edge] = mg2.BC_NODE_IS_FIXED_VALUE
    >>> mg2.status_at_node[mg2.nodes_at_top_edge] = mg2.BC_NODE_IS_CLOSED
    >>> mg2.status_at_node[mg2.nodes_at_bottom_edge] = mg2.BC_NODE_IS_CLOSED
    >>> mg2.status_at_node[mg2.nodes_at_right_edge] = mg2.BC_NODE_IS_CLOSED
    >>> fr2 = FlowAccumulator(mg2, flow_director="D8")
    >>> sp2 = StreamPowerEroder(mg2, K_sp=0.1, m_sp=0.0, n_sp=2.0, threshold_sp=2.0)
    >>> fr2.run_one_step()
    >>> sp2.run_one_step(dt=10.0)
    >>> z.reshape((3, 7))[1, :]
    array([ 0.        ,  1.        ,  4.        ,  8.52493781,
           13.29039716, 18.44367965, 36.        ])

    >>> mg3 = RasterModelGrid((5, 5), xy_spacing=2.0)
    >>> z = mg.node_x / 100.0
    >>> z = mg3.add_field("topographic__elevation", z, at="node")
    >>> mg3.status_at_node[mg3.nodes_at_left_edge] = mg3.BC_NODE_IS_FIXED_VALUE
    >>> mg3.status_at_node[mg3.nodes_at_top_edge] = mg3.BC_NODE_IS_CLOSED
    >>> mg3.status_at_node[mg3.nodes_at_bottom_edge] = mg3.BC_NODE_IS_CLOSED
    >>> mg3.status_at_node[mg3.nodes_at_right_edge] = mg3.BC_NODE_IS_CLOSED
    >>> mg3.at_node["water__unit_flux_in"] = mg3.node_y
    >>> fr3 = FlowAccumulator(mg3, flow_director="D8")
    >>> sp3 = StreamPowerEroder(
    ...     mg3,
    ...     K_sp=1.0,
    ...     sp_type="Unit",
    ...     a_sp=1.0,
    ...     b_sp=0.5,
    ...     c_sp=1.0,
    ...     discharge_field="surface_water__discharge",
    ... )
    >>> fr3.run_one_step()
    >>> sp3.run_one_step(1.0)
    >>> z
    array([0.        , 0.1       , 0.2       , 0.3       , 0.4       ,
           0.        , 0.02898979, 0.0859932 , 0.17463772, 0.4       ,
           0.        , 0.02240092, 0.06879049, 0.14586033, 0.4       ,
           0.        , 0.01907436, 0.05960337, 0.12929386, 0.4       ,
           0.        , 0.1       , 0.2       , 0.3       , 0.4       ])

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

    _name = "StreamPowerEroder"

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
        threshold_sp=0.0,
        sp_type="set_mn",
        m_sp=0.5,
        n_sp=1.0,
        a_sp=None,
        b_sp=None,
        c_sp=None,
        channel_width_field=1.0,
        discharge_field="drainage_area",
        erode_flooded_nodes=True,
    ):
        """Initialize the StreamPowerEroder.

        Parameters
        ----------
        grid : ModelGrid
            A grid.
        K_sp : float, array, or field name
            K in the stream power equation (units vary with other parameters).
        threshold_sp : positive float, array, or field name, optional
            The threshold stream power, below which no erosion occurs. This
            threshold is assumed to be in "stream power" units, i.e., if
            sp_type is 'Shear_stress', the value should be tau**a.
        sp_type : {'set_mn', 'Total', 'Unit', 'Shear_stress'}
            Controls how the law is implemented. If 'set_mn', use the supplied
            values of m_sp and n_sp. Else, component will derive values of m and n
            from supplied values of a_sp, b_sp, and c_sp, following Whipple and
            Tucker:

            *  If ``'Total'``, ``m = a * c``, ``n = a``.
            *  If ``'Unit'``, ``m = a * c *(1 - b)``, ``n = a``.
            *  If ``'Shear_stress'``, ``m = 2 * a * c * (1 - b) / 3``,
               ``n = 2 * a / 3``.

        m_sp : float, optional
            m in the stream power equation (power on drainage area). Overridden if
            a_sp, b_sp, and c_sp are supplied.
        n_sp : float, optional, ~ 0.5<n_sp<4.
            n in the stream power equation (power on slope). Overridden if
            a_sp, b_sp, and c_sp are supplied.
        a_sp : float, optional
            The power on the SP/shear term to get the erosion rate; the "erosional
            process" term. Only used if sp_type is not 'set_mn'.
        b_sp : float, optional
            The power on discharge to get width; the "hydraulic geometry" term.
            Only used if sp_type in ('Unit', 'Shear_stress').
        c_sp : float, optional
            The power on area to get discharge; the "basin hydology" term. Only
            used if sp_type is not 'set_mn'.
        channel_width_field : None, float, array, or field name, optional
            If not None, component will look for node-centered data describing
            channel width or if an array, will take the array as the channel
            widths. It will use the widths to implement incision ~ stream power
            per unit width. If sp_type is 'set_mn', follows the equation given
            above. If sp_type in ('Unit', 'Shear_stress'), the width value will
            be implemented directly. W has no effect if sp_type is 'Total'.
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
                "verified that StreamPowerEroder is compatible with "
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

        self._A = return_array_at_node(grid, discharge_field)
        self._elevs = return_array_at_node(grid, "topographic__elevation")
        self._sp_crit = return_array_at_node(grid, threshold_sp)

        # use setter for K defined below
        self.K = K_sp

        assert np.all(self._sp_crit >= 0.0)

        if discharge_field == "drainage_area":
            self._use_Q = False
        else:
            self._use_Q = True

        if channel_width_field is None:
            self._use_W = False
        else:
            self._use_W = True
            self._W = return_array_at_node(grid, channel_width_field)

        if np.any(threshold_sp != 0.0):
            self._set_threshold = True
            # ^flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
        else:
            self._set_threshold = False

        self._type = sp_type
        if sp_type == "set_mn":
            assert (float(m_sp) >= 0.0) and (
                float(n_sp) >= 0.0
            ), "m and n must be positive"
            self._m = float(m_sp)
            self._n = float(n_sp)
            assert (
                (a_sp is None) and (b_sp is None) and (c_sp is None)
            ), "If sp_type is 'set_mn', do not pass values for a, b, or c!"
        else:
            assert sp_type in ("Total", "Unit", "Shear_stress"), (
                "sp_type not recognised. It must be 'set_mn', 'Total', "
                + "'Unit', or 'Shear_stress'."
            )
            assert (
                m_sp == 0.5 and n_sp == 1.0
            ), "Do not set m and n if sp_type is not 'set_mn'!"
            assert float(a_sp) >= 0.0, "a must be positive"
            self._a = float(a_sp)
            if b_sp is not None:
                assert float(b_sp) >= 0.0, "b must be positive"
                self._b = float(b_sp)
            else:
                assert self._use_W, "b was not set"
                self._b = 0.0
            if c_sp is not None:
                assert float(c_sp) >= 0.0, "c must be positive"
                self._c = float(c_sp)
            else:
                assert self._use_Q, "c was not set"
                self._c = 1.0
            if self._type == "Total":
                self._n = self._a
                self._m = self._a * self._c  # ==_a if use_Q
            elif self._type == "Unit":
                self._n = self._a
                self._m = self._a * self._c * (1.0 - self._b)
                # ^ ==_a iff use_Q&use_W etc
            elif self._type == "Shear_stress":
                self._m = 2.0 * self._a * self._c * (1.0 - self._b) / 3.0
                self._n = 2.0 * self._a / 3.0
            else:
                raise MissingKeyError(
                    "Not enough information was provided on the exponents to use!"
                )

        # m and n will always be set, but care needs to be taken to include Q
        # and W directly if appropriate

        self._stream_power_erosion = self._grid.zeros(centering="node")
        self._alpha = self._grid.zeros("node")

    @property
    def K(self):
        """Erodibility (units depend on m_sp)."""
        return self._K

    @K.setter
    def K(self, new_val):
        self._K = return_array_at_node(self._grid, new_val)

    def run_one_step(self, dt):
        """A simple, explicit implementation of a stream power algorithm.

        If you are routing across flooded depressions in your flow routing
        scheme, be sure to set *erode_flooded_nodes* flag in the instantiation
        of the component to ensure erosion cannot occur in the lake. Erosion
        is always zero if the gradient is adverse, but can still procede as
        usual on the entry into the depression unless this flag is set.

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

        upstream_order_IDs = self._grid["node"]["flow__upstream_node_order"]

        defined_flow_receivers = np.not_equal(
            self._grid["node"]["flow__link_to_receiver_node"], self._grid.BAD_INDEX
        )

        try:
            length_of_link = self._grid.length_of_d8
        except AttributeError:
            length_of_link = self._grid.length_of_link

        flow_link_lengths = length_of_link[
            self._grid.at_node["flow__link_to_receiver_node"][defined_flow_receivers]
        ]
        flow_receivers = self._grid["node"]["flow__receiver_node"]

        # Operate the main function:
        if self._use_W:
            self._alpha[defined_flow_receivers] = (
                self._K[defined_flow_receivers]
                * dt
                * self._A[defined_flow_receivers] ** self._m
                / self._W[defined_flow_receivers]
                / (flow_link_lengths**self._n)
            )

        else:
            self._alpha[defined_flow_receivers] = (
                self._K[defined_flow_receivers]
                * dt
                * self._A[defined_flow_receivers] ** self._m
                / (flow_link_lengths**self._n)
            )

        # Handle flooded nodes, if any (no erosion there)
        if flooded_nodes is not None:
            self._alpha[flooded_nodes] = 0.0

        reversed_flow = self._elevs < self._elevs[flow_receivers]
        # this check necessary if flow has been routed across
        # depressions
        self._alpha[reversed_flow] = 0.0

        threshdt = self._sp_crit * dt

        # solve using Brent's Method in Cython for Speed
        if isinstance(threshdt, float):
            brent_method_erode_fixed_threshold(
                upstream_order_IDs,
                flow_receivers,
                threshdt,
                self._alpha,
                self._n,
                self._elevs,
            )
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs,
                flow_receivers,
                threshdt,
                self._alpha,
                self._n,
                self._elevs,
            )
