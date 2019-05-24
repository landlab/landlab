# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy as np

from landlab import BAD_INDEX_VALUE as UNDEFINED_INDEX, Component
from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError
from landlab.utils.decorators import use_file_name_or_kwds

from .cfuncs import (
    brent_method_erode_fixed_threshold,
    brent_method_erode_variable_threshold,
)


class StreamPowerEroder(Component):
    """Erode where channels are.

    Implemented as:

    .. math::
        E = K A^m S^n - sp_{crit},

    and if :math:`E < 0`, :math:`E = 0`.

    If ``use_W`` is declared and ``True``, the module instead implements:

    .. math::
        E = K A^m S^n / W - sp_{crit}

    DEJH Sept 2013, major modifications Sept 14 and May 16. This component
    now wraps Fastscape-style functionality under the hood.

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

    NB: If you want spatially or temporally variable runoff, pass the
    runoff values at each pixel to the flow router using the input argument
    *use_Q*.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> from landlab.components import FlowAccumulator, StreamPowerEroder
    >>> mg = RasterModelGrid((5, 5), xy_spacing=10.)
    >>> z = np.array([7.,  7.,  7.,  7.,  7.,
    ...               7.,  5., 3.2,  6.,  7.,
    ...               7.,  2.,  3.,  5.,  7.,
    ...               7.,  1., 1.9,  4.,  7.,
    ...               7.,  0.,  7.,  7.,  7.])
    >>> z = mg.add_field('node', 'topographic__elevation', z)
    >>> fr = FlowAccumulator(mg, flow_director='D8')
    >>> sp = StreamPowerEroder(mg, K_sp=1.)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=1.)
    >>> z  # doctest: +NORMALIZE_WHITESPACE
    array([ 7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
            7.        ,  2.92996598,  2.02996598,  4.01498299,  7.        ,
            7.        ,  0.85993197,  1.87743897,  3.28268321,  7.        ,
            7.        ,  0.28989795,  0.85403051,  2.42701526,  7.        ,
            7.        ,  0.        ,  7.        ,  7.        ,  7.        ])

    >>> mg2 = RasterModelGrid((3, 7))
    >>> z = np.array(mg2.node_x**2.)
    >>> z = mg2.add_field('node', 'topographic__elevation', z)
    >>> mg2.status_at_node[mg2.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> fr2 = FlowAccumulator(mg2, flow_director='D8')
    >>> sp2 = StreamPowerEroder(mg2, K_sp=0.1, m_sp=0., n_sp=2.,
    ...                         threshold_sp=2.)
    >>> fr2.run_one_step()
    >>> sp2.run_one_step(dt=10.)
    >>> z.reshape((3, 7))[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   1.        ,   4.        ,   8.52493781,
            13.29039716,  18.44367965,  36.        ])

    >>> mg3 = RasterModelGrid((5, 5), xy_spacing=2.)
    >>> z = mg.node_x/100.
    >>> z = mg3.add_field('node', 'topographic__elevation', z)
    >>> mg3.status_at_node[mg3.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> mg3.at_node['water__unit_flux_in'] = mg3.node_y
    >>> fr3 = FlowAccumulator(mg3, flow_director='D8')
    >>> Q = mg3.at_node['surface_water__discharge']
    >>> sp3 = StreamPowerEroder(mg3, K_sp=1., sp_type='Unit', a_sp=1.,
    ...                         b_sp=0.5, c_sp=1., use_Q=Q)
    >>> fr3.run_one_step()
    >>> sp3.run_one_step(1.)
    >>> z
    array([ 0.        ,  0.1       ,  0.2       ,  0.3       ,  0.4       ,
            0.        ,  0.02898979,  0.0859932 ,  0.17463772,  0.4       ,
            0.        ,  0.02240092,  0.06879049,  0.14586033,  0.4       ,
            0.        ,  0.01907436,  0.05960337,  0.12929386,  0.4       ,
            0.        ,  0.1       ,  0.2       ,  0.3       ,  0.4       ])
    """

    _name = "StreamPowerEroder"

    _input_var_names = (
        "topographic__elevation",
        "flow__link_to_receiver_node",
        "drainage_area",
        "flow__receiver_node",
        "flow__upstream_node_order",
        "topographic__steepest_slope",
    )

    _output_var_names = ("topographic__elevation",)

    _var_units = {
        "topographic__elevation": "m",
        "drainage_area": "m**2",
        "flow__link_to_receiver_node": "-",
        "flow__receiver_node": "-",
        "flow__upstream_node_order": "-",
        "topographic__steepest_slope": "-",
    }

    _var_mapping = {
        "topographic__elevation": "node",
        "drainage_area": "node",
        "flow__link_to_receiver_node": "node",
        "flow__receiver_node": "node",
        "flow__upstream_node_order": "node",
        "topographic__steepest_slope": "node",
    }

    _var_doc = {
        "topographic__elevation": "Land surface topographic elevation",
        "drainage_area": "Upstream accumulated surface area contributing to the node's "
        "discharge",
        "flow__link_to_receiver_node": "ID of link downstream of each node, which carries the discharge",
        "flow__receiver_node": "Node array of receivers (node that receives flow from current "
        "node)",
        "flow__upstream_node_order": "Node array containing downstream-to-upstream ordered list of "
        "node IDs",
        "topographic__steepest_slope": "Node array of steepest *downhill* slopes",
    }

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        K_sp=None,
        threshold_sp=0.0,
        sp_type="set_mn",
        m_sp=0.5,
        n_sp=1.0,
        a_sp=None,
        b_sp=None,
        c_sp=None,
        use_W=None,
        use_Q=None,
        **kwds
    ):
        """Initialize the StreamPowerEroder

        Parameters
        ----------
        grid : ModelGrid
            A grid.
        K_sp : float, array, or field name
            K in the stream power equation (units vary with other parameters).
        threshold_sp : positive float, optional
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
        use_W : None, array, or field name, optional
            If not None, component will look for node-centered data describing
            channel width in grid.at_node[use_W] or if an array, will take the
            array as the channel widths. It will use the widths to implement
            incision ~ stream power per unit width. If sp_type is 'set_mn',
            follows the equation given above. If sp_type in ('Unit',
            'Shear_stress'), the width value will be implemented directly. W has no
            effect if sp_type is 'Total'.
        use_Q : None, array, or field name, optional
            If not None, the equation becomes E=K*Q**m*S**n. Effectively sets c=1
            in Wh&T's 1999 derivation, if you are setting m and n through a, b,
            and c.
        """
        if "flow__receiver_node" in grid.at_node:
            if grid.at_node["flow__receiver_node"].size != grid.size("node"):
                msg = (
                    "A route-to-multiple flow director has been "
                    "run on this grid. The landlab development team has not "
                    "verified that StreamPowerEroder is compatible with "
                    "route-to-multiple methods. Please open a GitHub Issue "
                    "to start this process."
                )
                raise NotImplementedError(msg)

        if type(use_Q) is str and use_Q == "water__discharge":
            use_Q = "surface_water__discharge"
        self._grid = grid

        self.use_K = False  # grandfathered in; only if K_sp == 'array'
        if type(K_sp) is np.ndarray:
            self._K_unit_time = K_sp
        else:
            try:
                self._K_unit_time = self.grid.zeros("node", dtype=float)
                self._K_unit_time.fill(K_sp)
            except ValueError:  # could not cast => was a str
                if K_sp == "array":
                    self.use_K = True
                else:
                    self._K_unit_time = grid.at_node[K_sp]

        assert np.all(threshold_sp >= 0.0)
        # for now, enforce threshold as a float
        assert type(threshold_sp) in (float, int)
        try:
            self.sp_crit = float(threshold_sp)
        except TypeError:
            try:
                self.sp_crit = self.grid.at_node[threshold_sp]
            except TypeError:  # was an array
                self.sp_crit = threshold_sp
                assert self.sp_crit.size == self.grid.number_of_nodes
        if np.any(threshold_sp != 0.0):
            self.set_threshold = True
            # ^flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
        else:
            self.set_threshold = False
        try:
            self.tstep = kwds["dt"]
        except KeyError:
            self.tstep = None
            # retained for back compatibility; undocumented functionality
        if type(use_W) is bool:  # again for back-compatibility
            self.use_W = use_W
            self._W = None
        elif use_W is None:
            self.use_W = False
            self._W = None
        else:
            self.use_W = True
            try:
                self._W = self.grid.at_node[use_W]
            except (FieldError, TypeError):
                assert use_W.size == self._grid.number_of_nodes
                self._W = use_W
        if type(use_Q) is bool:
            self.use_Q = use_Q
            self._Q = None
        elif use_Q is None:
            self.use_Q = False
            self._Q = None
        else:
            self.use_Q = True
            try:
                self._Q = self.grid.at_node[use_Q]
            except (FieldError, TypeError):
                assert use_Q.size == self._grid.number_of_nodes
                self._Q = use_Q
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
                assert self.use_W, "b was not set"
                self._b = 0.0
            if c_sp is not None:
                assert float(c_sp) >= 0.0, "c must be positive"
                self._c = float(c_sp)
            else:
                assert self.use_Q, "c was not set"
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
                    "Not enough information was provided " + "on the exponents to use!"
                )
        # m and n will always be set, but care needs to be taken to include Q
        # and W directly if appropriate

        self.stream_power_erosion = grid.zeros(centering="node")
        self.alpha = self.grid.zeros("node")

    def erode(
        self,
        grid,
        dt,
        elevs="topographic__elevation",
        drainage_areas="drainage_area",
        flow_receiver="flow__receiver_node",
        order_upstream="flow__upstream_node_order",
        slopes_at_nodes="topographic__steepest_slope",
        link_mapping="flow__link_to_receiver_node",
        link_slopes=None,
        slopes_from_elevs=None,
        W_if_used=None,
        Q_if_used=None,
        K_if_used=None,
        flooded_nodes=None,
    ):
        """
        .. note:: deprecated
            This run method is now DEPRECATED. Use the fully standardized
            method :func:`run_one_step` instead.

        A simple, explicit implementation of a stream power algorithm.

        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        dt : float
            Component time step.

        elevs : str or ndarray, optional
            Elevations on the grid, either a field string or nnodes-long array.

        drainage_areas: str or ndarray, optional
            Tells the component where to look for the drainage area values.
            Change to another string to override which grid field the
            component looks at, or pass a nnodes-long array of drainage
            areas values directly instead.

        flow_receiver, order_upstream : str or ndarray, optional
            The downstream node to which each node flows and the ordering of
            the nodes in the network starting at the outlet, respectively,
            are both necessary as inputs to allow stability testing.

            If you already have slopes defined at nodes on the grid, pass them
            to the component with *slopes_at_nodes*. The same syntax is
            expected: string gives a name in the grid fields, an array gives
            values direct.

            Alternatively, set *link_slopes* (and *link_mapping*) if this
            data
            is only available at links. 'topographic__derivative_of_elevation'
            is the default field name for link slopes. Override this name by
            setting the variable as the appropriate string, or override use of
            grid fields altogether by passing an array. *link_mapping*
            controls how the component maps these link values onto the arrays.
            We assume there is always a 1:1 mapping (pass the values already
            projected onto the nodes using slopes_at_nodes if not). Other
            components, e.g., flow_routing.route_flow_dn, may provide the
            necessary outputs to make the mapping easier: e.g., just pass
            'flow__link_to_receiver_node' from that module (the default name).
            If the component cannot find an existing mapping through this
            parameter, it will derive one on the fly, at considerable cost of
            speed (see on-screen reports).

        slopes_from_elevs : str, optional
            Allows the module to create gradients internally
            from elevations rather than have them provided. Set to True to
            force the component to look for the data in the location specified
            by elevs. Using this option is
            considerably slower than any of the alternatives, as it also has to
            calculate the link_mapping from stratch each time.

            In both these cases, at present the mapping is to use the maximum
            slope of *any* link attached to the node as the representative
            node slope. This is primarily for speed, but may be a good idea
            to modify later.

        W_if_used, Q_if_used : str or ndarray, optional
            Must be provided if you set *use_W* and *use_Q* respectively in
            the component initialization. They can be either field names or
            nnodes arrays as in the other cases.

            If you are routing across flooded depressions in your flow routing
            scheme, be sure to set *flooded_nodes* with a boolean array or
            array of IDs to ensure erosion cannot occur in the lake. Erosion
            is always zero if the gradient is adverse, but can still procede as
            usual on the entry into the depression unless *flooded_nodes* is
            set.

            NB: If you want spatially or temporally variable runoff, pass the
            runoff values at each pixel to the flow router, then pass
            discharges at each node using *Q_if_used* to this component.

        Returns
        -------
        tuple
            Tuple of (*grid*, *modified_elevs*, *stream_power_erosion*);
            modifies grid elevation fields to reflect updates. Note the value
            stream_power_erosion is not an excess stream power; any specified
            erosion threshold is not incorporated into it.
        """
        if grid.at_node["flow__receiver_node"].size != grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that StreamPowerEroder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)
        if type(order_upstream) is str:
            upstream_order_IDs = grid.at_node[order_upstream]
        else:
            upstream_order_IDs = self._grid["node"][order_upstream]

        defined_flow_receivers = np.not_equal(
            self._grid["node"][link_mapping], UNDEFINED_INDEX
        )

        try:
            length_of_link = self._grid.length_of_d8
        except AttributeError:
            length_of_link = self._grid.length_of_link

        flow_link_lengths = length_of_link[
            self._grid.at_node[link_mapping][defined_flow_receivers]
        ]
        flow_receivers = self.grid["node"][flow_receiver]

        if W_if_used is not None:
            assert self.use_W, (
                "Widths were provided, but you didn't set "
                + "the use_W flag in your input file! "
                + "Aborting..."
            )
            assert self._W is None, (
                "Do not pass W to the run method "
                + "if you also set them at initialization!"
            )

        if Q_if_used is not None:
            assert self.use_Q, (
                "Discharges were provided, but you didn't "
                + "set the use_Q flag in your input file! "
                + "Aborting..."
            )
            assert self._Q is None, (
                "Do not pass Q to the run method "
                + "if you also set them at initialization!"
            )

        if K_if_used is not None:
            assert self.use_K, (
                "An array of erodabilities was provided, "
                + "but you didn't set K_sp to 'array' in your "
                + "input file! Aborting..."
            )
            try:
                _K_unit_time = grid.at_node[K_if_used]
            except TypeError:
                _K_unit_time = K_if_used
        else:
            # little move to save a bit of memory management time...
            if flooded_nodes is not None:
                _K_unit_time = self._K_unit_time.copy()
            else:
                _K_unit_time = self._K_unit_time

        if type(elevs) is str:
            z = grid.at_node[elevs]
        else:
            z = elevs

        if type(drainage_areas) is str:
            A = grid.at_node[drainage_areas]
        else:
            A = drainage_areas

        # Disable incision in flooded nodes, as appropriate
        if flooded_nodes is not None:
            _K_unit_time[flooded_nodes] = 0.0

        # Operate the main function:
        if self.use_W is False and self.use_Q is False:  # normal case
            self.alpha[defined_flow_receivers] = (
                _K_unit_time[defined_flow_receivers]
                * dt
                * A[defined_flow_receivers] ** self._m
                / (flow_link_lengths ** self._n)
            )
            # Handle flooded nodes, if any (no erosion there)
            if flooded_nodes is not None:
                self.alpha[flooded_nodes] = 0.0
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across
            # depressions
            self.alpha[reversed_flow] = 0.0

            threshdt = self.sp_crit * dt

        elif self.use_W:
            if self._W is None:
                try:
                    W = grid.at_node[W_if_used]
                except TypeError:
                    W = W_if_used
            else:
                W = self._W
            if self.use_Q:  # use both Q and W direct
                if self._Q is None:
                    try:
                        Q_direct = grid.at_node[Q_if_used]
                    except TypeError:
                        Q_direct = Q_if_used
                else:
                    Q_direct = self._Q
                self.alpha[defined_flow_receivers] = (
                    _K_unit_time[defined_flow_receivers]
                    * dt
                    * Q_direct[defined_flow_receivers] ** self._m
                    / W[defined_flow_receivers]
                    / (flow_link_lengths ** self._n)
                )
                # Handle flooded nodes, if any (no erosion there)
                if flooded_nodes is not None:
                    self.alpha[flooded_nodes] = 0.0
                reversed_flow = z < z[flow_receivers]
                # this check necessary if flow has been routed across
                # depressions
                self.alpha[reversed_flow] = 0.0

                threshdt = self.sp_crit * dt

            else:  # just W to be used
                self.alpha[defined_flow_receivers] = (
                    _K_unit_time[defined_flow_receivers]
                    * dt
                    * A[defined_flow_receivers] ** self._m
                    / W[defined_flow_receivers]
                    / flow_link_lengths ** self._n
                )
                # Handle flooded nodes, if any (no erosion there)
                if flooded_nodes is not None:
                    self.alpha[flooded_nodes] = 0.0
                reversed_flow = z < z[flow_receivers]
                # this check necessary if flow has been routed across
                # depressions
                self.alpha[reversed_flow] = 0.0

                threshdt = self.sp_crit * dt

        else:  # just use_Q
            if self._Q is None:
                try:
                    Q_direct = grid.at_node[Q_if_used]
                except TypeError:
                    assert type(Q_if_used) in (np.ndarray, list)
                    Q_direct = Q_if_used
            else:
                Q_direct = self._Q
            self.alpha[defined_flow_receivers] = (
                _K_unit_time[defined_flow_receivers]
                * dt
                * Q_direct[defined_flow_receivers] ** self._m
                / flow_link_lengths ** self._n
            )
            # Handle flooded nodes, if any (no erosion there)
            if flooded_nodes is not None:
                self.alpha[flooded_nodes] = 0.0
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across
            # depressions
            self.alpha[reversed_flow] = 0.0

            threshdt = self.sp_crit * dt

        # solve using Brent's Method in Cython for Speed
        if isinstance(threshdt, float):
            brent_method_erode_fixed_threshold(
                upstream_order_IDs, flow_receivers, threshdt, self.alpha, self._n, z
            )
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs, flow_receivers, threshdt, self.alpha, self._n, z
            )

        return grid, z, self.stream_power_erosion

    def run_one_step(self, dt, flooded_nodes=None, **kwds):
        """
        A simple, explicit implementation of a stream power algorithm.

        This component now looks exclusively for the field
        'topographic__steepest_slope' at each node to determine the local
        slope (previoiusly it was possible to map values from links explicitly
        within the component, but this functionality is now deprecated).

        If you are routing across flooded depressions in your flow routing
        scheme, be sure to set *flooded_nodes* with a boolean array or array
        of IDs to ensure erosion cannot occur in the lake. Erosion
        is always zero if the gradient is adverse, but can still procede as
        usual on the entry into the depression unless *flooded_nodes* is set.

        Parameters
        ----------
        dt : float
            Time-step size
        flooded_nodes : ndarray of int (optional)
            IDs of nodes that are flooded and should have no erosion. If not
            provided but flow has still been routed across depressions, erosion
            may still occur beneath the apparent water level (though will
            always still be positive).
        """
        self.erode(grid=self._grid, dt=dt, flooded_nodes=flooded_nodes)
