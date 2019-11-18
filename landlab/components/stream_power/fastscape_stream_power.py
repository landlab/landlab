#! /usr/env/python
"""Fastscape stream power erosion."""

# This module attempts to "component-ify" GT's Fastscape stream
# power erosion.
# Created DEJH, March 2014.
from __future__ import print_function

import numpy as np
from six import string_types

from landlab import BAD_INDEX_VALUE as UNDEFINED_INDEX, Component, RasterModelGrid
from landlab.utils.decorators import use_file_name_or_kwds

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
    >>> previous_z = z.copy()
    >>> sp.run_one_step(1., rainfall_intensity_if_used=0.)
    >>> np.allclose(z, previous_z)
    True
    """

    _name = "FastscapeEroder"

    _input_var_names = (
        "topographic__elevation",
        "drainage_area",
        "flow__link_to_receiver_node",
        "flow__upstream_node_order",
        "flow__receiver_node",
    )

    _output_var_names = ("topographic__elevation",)

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

    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        K_sp=None,
        m_sp=0.5,
        n_sp=1.0,
        threshold_sp=0.0,
        rainfall_intensity=1.0,
        discharge_name="drainage_area",
        **kwds
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
        discharge_name : string; optional
            Name of field to use for discharge proxy. Defaults to 'drainage_area',
            which means the component will expect the driver or another component
            to have created and populated a 'drainage_area' field. To use a
            different field, such as 'surface_water__discharge', give its name in
            this argument.
        """
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

        self._grid = grid

        self.K = K_sp  # overwritten below in special cases
        self.m = float(m_sp)
        self.n = float(n_sp)
        if isinstance(threshold_sp, (float, int)):
            self.thresholds = float(threshold_sp)
        else:
            if isinstance(threshold_sp, string_types):
                self.thresholds = self.grid.at_node[threshold_sp]
            else:
                self.thresholds = threshold_sp
            assert self.thresholds.size == self.grid.number_of_nodes

        # make storage variables
        self.A_to_the_m = grid.zeros(at="node")
        self.alpha = grid.empty(at="node")

        if self.K is None:
            raise ValueError(
                "K_sp must be set as a float, node array, or "
                + "field name. It was None."
            )

        # now handle the inputs that could be float, array or field name:
        # some support here for old-style inputs
        if isinstance(K_sp, string_types):
            if K_sp == "array":
                self.K = None
            else:
                self.K = self._grid.at_node[K_sp]
        elif isinstance(K_sp, (float, int)):
            self.K = float(K_sp)
        else:
            self.K = np.asarray(K_sp, dtype=float)
            if len(self.K) != self.grid.number_of_nodes:
                raise TypeError("Supplied value of K_sp is not n_nodes long")

        if isinstance(rainfall_intensity, string_types):
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
        elif len(rainfall_intensity) == self.grid.number_of_nodes:
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

        # We now forbid changing of the field name
        if "value_field" in kwds.keys():
            raise ValueError(
                "This component can no longer support variable"
                'field names. Use "topographic__elevation".'
            )

        # Handle option for area vs discharge
        self.discharge_name = discharge_name

    def erode(
        self,
        grid_in,
        dt=None,
        K_if_used=None,
        flooded_nodes=None,
        rainfall_intensity_if_used=None,
    ):
        """Erode using stream power erosion.

        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power
        scheme.

        This driving method for this component is now superceded by the new,
        standardized wrapper :func:`run_one_step`, but is retained for
        back compatibility.

        Set *K_if_used* as a field name or nnodes-long array if you set
        *K_sp* as *"array"* during initialization.

        It returns the grid, in which it will have modified the value of
        *value_field*, as specified in component initialization.

        Parameters
        ----------
        grid_in : a grid
            This is a dummy argument maintained for component back-
            compatibility. It is superceded by the copy of the grid passed
            during initialization.
        dt : float
            Time-step size. If you are calling the deprecated function
            :func:`gear_timestep`, that method will supercede any value
            supplied here.
        K_if_used : array (optional)
            Set this to an array if you set K_sp to 'array' in your input file.
        flooded_nodes : ndarray of int (optional)
            IDs of nodes that are flooded and should have no erosion. If not
            provided but flow has still been routed across depressions, erosion
            may still occur beneath the apparent water level (though will
            always still be positive).
        rainfall_intensity_if_used : float or None (optional)
            Supply to drive this component with a time-varying spatially
            constant rainfall.

        Returns
        -------
        grid
            A reference to the grid.
        """
        if self._grid.at_node["flow__receiver_node"].size != self._grid.size("node"):
            msg = (
                "A route-to-multiple flow director has been "
                "run on this grid. The landlab development team has not "
                "verified that FastscapeEroder is compatible with "
                "route-to-multiple methods. Please open a GitHub Issue "
                "to start this process."
            )
            raise NotImplementedError(msg)

        upstream_order_IDs = self._grid.at_node["flow__upstream_node_order"]
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
        if isinstance(self.K, np.ndarray):
            K_here = self.K[defined_flow_receivers]
        else:
            K_here = self.K
        if rainfall_intensity_if_used is not None:
            assert type(rainfall_intensity_if_used) in (float, np.float64, int)
            r_i_here = float(rainfall_intensity_if_used)
        else:
            r_i_here = self._r_i

        if dt is None:
            dt = self.dt
        assert dt is not None, (
            "Fastscape component could not find a dt to "
            + "use. Pass dt to the run_one_step() method."
        )

        if self.K is None:  # "old style" setting of array
            assert K_if_used is not None
            self.K = K_if_used

        n = float(self.n)

        np.power(self._grid["node"][self.discharge_name], self.m, out=self.A_to_the_m)
        self.alpha[defined_flow_receivers] = (
            r_i_here ** self.m
            * K_here
            * dt
            * self.A_to_the_m[defined_flow_receivers]
            / (flow_link_lengths ** self.n)
        )

        flow_receivers = self._grid["node"]["flow__receiver_node"]
        alpha = self.alpha

        # Handle flooded nodes, if any (no erosion there)
        if flooded_nodes is not None:
            alpha[flooded_nodes] = 0.0
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            alpha[reversed_flow] = 0.0

        threshsdt = self.thresholds * dt

        # solve using Brent's Method in Cython for Speed
        if isinstance(self.thresholds, float):
            brent_method_erode_fixed_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z
            )
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z
            )

        return self._grid

    def run_one_step(
        self, dt, flooded_nodes=None, rainfall_intensity_if_used=None, **kwds
    ):
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
        flooded_nodes : ndarray of int (optional)
            IDs of nodes that are flooded and should have no erosion. If not
            provided but flow has still been routed across depressions, erosion
            may still occur beneath the apparent water level (though will
            always still be positive).
        rainfall_intensity_if_used : float or None (optional)
            Supply to drive this component with a time-varying spatially
            constant rainfall.
        """
        self.erode(
            grid_in=self._grid,
            dt=dt,
            flooded_nodes=flooded_nodes,
            rainfall_intensity_if_used=rainfall_intensity_if_used,
        )
