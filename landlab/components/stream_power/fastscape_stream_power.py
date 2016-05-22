#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""
from __future__ import print_function

import numpy
import warnings
from landlab import ModelParameterDictionary, Component
from landlab.core.model_parameter_dictionary import MissingKeyError, \
    ParameterValueError
from landlab.utils.decorators import use_file_name_or_kwds
from landlab.field.scalar_data_fields import FieldError
from scipy.optimize import newton, fsolve

UNDEFINED_INDEX = -1


class FastscapeEroder(Component):
    '''
    This class uses the Braun-Willett Fastscape approach to calculate the
    amount of erosion at each node in a grid, following a stream power
    framework. This should allow it to be stable against larger timesteps
    than an explicit stream power scheme.

    Stream power erosion is implemented as::

        E = K * (rainfall_intensity*A)**m * S**n - threshold_sp,

    if K * A**m * S**n > threshold_sp, and::

        E = 0,

    if K * A**m * S**n <= threshold_sp.

    This module assumes you have already run
    :func:`landlab.components.flow_routing.route_flow_dn.FlowRouter.route_flow`
    in the same timestep. It looks for 'flow__upstream_node_order',
    'flow__link_to_receiver_node', 'drainage_area', 'flow__receiver_node', and
    'topographic__elevation' at the nodes in the grid. 'drainage_area' should
    be in area upstream, not volume (i.e., set runoff_rate=1.0 when calling
    FlowRouter.route_flow).

    The primary method of this class is :func:`run_one_step`.

    Construction::

        FastscapeEroder(grid, K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=0.,
                        rainfall_intensity=1.)

    Parameters
    ----------
    grid : ModelGrid
        A grid.
    K_sp : float, array, or field name
        K in the stream power equation (units vary with other parameters).
    m_sp : float, optional
        m in the stream power equation (power on drainage area).
    n_sp : float, optional, ~ 0.5<n_sp<4.
        n in the stream power equation (power on slope).
        Performance will be VERY degraded if n < 1.
    threshold_sp : float, array, or field name
        The threshold stream power.
    rainfall intensity : float, array, or field name; optional
        Modifying factor on drainage area to convert it to a true water
        volume flux in (m/time). i.e., E = K * (r_i*A)**m * S**n

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> from landlab.components import FlowRouter
    >>> mg = RasterModelGrid((5, 5), 10.)
    >>> z = np.array([7.,  7.,  7.,  7.,  7.,
    ...               7.,  5., 3.2,  6.,  7.,
    ...               7.,  2.,  3.,  5.,  7.,
    ...               7.,  1., 1.9,  4.,  7.,
    ...               7.,  0.,  7.,  7.,  7.])
    >>> z = mg.add_field('node', 'topographic__elevation', z)
    >>> fr = FlowRouter(mg)
    >>> sp = FastscapeEroder(mg, K_sp=1.)
    >>> fr.run_one_step()
    >>> sp.run_one_step(dt=1.)
    >>> z  # doctest: +NORMALIZE_WHITESPACE
    array([ 7.        ,  7.        ,  7.        ,  7.        ,  7.        ,
            7.        ,  2.92996598,  2.02996598,  4.01498299,  7.        ,
            7.        ,  0.85993197,  1.87743897,  3.28268321,  7.        ,
            7.        ,  0.28989795,  0.85403051,  2.42701526,  7.        ,
            7.        ,  0.        ,  7.        ,  7.        ,  7.        ])

    >>> mg2 = RasterModelGrid((3, 7), 1.)
    >>> z = np.array(mg2.node_x**2.)
    >>> z = mg2.add_field('node', 'topographic__elevation', z)
    >>> mg2.status_at_node[mg2.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> mg2.status_at_node[mg2.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> fr2 = FlowRouter(mg2)
    >>> sp2 = FastscapeEroder(mg2, K_sp=0.1, m_sp=0., n_sp=2.,
    ...                       threshold_sp=2.)
    >>> fr2.run_one_step()
    >>> sp2.run_one_step(dt=10.)
    >>> z.reshape((3, 7))[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   1.        ,   4.        ,   8.52493781,
            13.29039716,  18.44367965,  36.        ])

    >>> mg3 = RasterModelGrid((3, 7), 1.)
    >>> z = np.array(mg3.node_x**2.)
    >>> z = mg3.add_field('node', 'topographic__elevation', z)
    >>> mg3.status_at_node[mg3.nodes_at_left_edge] = FIXED_VALUE_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_top_edge] = CLOSED_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_bottom_edge] = CLOSED_BOUNDARY
    >>> mg3.status_at_node[mg3.nodes_at_right_edge] = CLOSED_BOUNDARY
    >>> fr3 = FlowRouter(mg3)
    >>> K_field = mg3.ones('node')  # K can be a field
    >>> sp3 = FastscapeEroder(mg3, K_sp=K_field, m_sp=1., n_sp=0.6,
    ...                       threshold_sp=mg3.node_x,
    ...                       rainfall_intensity=z)
    >>> fr3.run_one_step()
    >>> sp3.run_one_step(1.)
    >>> z.reshape((3, 7))[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   0.        ,   0.18508544,   0.42869679,
             0.85354205,   2.05732094,  36.        ])
    '''

    _name = 'FastscapeEroder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'flow__link_to_receiver_node',
        'flow__upstream_node_order',
        'flow__receiver_node',
    )

    _output_var_names = (
        'topographic__elevation',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'drainage_area': 'm**2',
        'flow__link_to_receiver_node': '-',
        'flow__upstream_node_order': '-',
        'flow__receiver_node': '-',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'drainage_area': 'node',
        'flow__link_to_receiver_node': 'node',
        'flow__upstream_node_order': 'node',
        'flow__receiver_node': 'node',
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'flow__link_to_receiver_node':
            'ID of link downstream of each node, which carries the discharge',
        'flow__upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow__receiver_node':
            'Node array of receivers (node that receives flow from current '
            'node)',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, K_sp=None, m_sp=0.5, n_sp=1., threshold_sp=0.,
                 rainfall_intensity=1., **kwds):
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
        """
        self._grid = grid

        self.K = K_sp  # overwritten below in special cases
        self.m = float(m_sp)
        self.n = float(n_sp)
        if type(threshold_sp) in (float, int):
            self.thresholds = float(threshold_sp)
        else:
            if type(threshold_sp) is str:
                self.thresholds = self.grid.at_node[threshold_sp]
            else:
                self.thresholds = threshold_sp
            assert self.thresholds.size == self.grid.number_of_nodes

        # make storage variables
        self.A_to_the_m = grid.zeros(at='node')
        self.alpha = grid.empty(at='node')
        self.alpha_by_flow_link_lengthtothenless1 = numpy.empty_like(
                                                        self.alpha)

        self.grid._diagonal_links_at_node  # calc number of diagonal links

        if self.K is None:
            raise ValueError('K_sp must be set as a float, node array, or ' +
                             'field name. It was None.')

        # now handle the inputs that could be float, array or field name:
        # some support here for old-style inputs
        if type(K_sp) is str:
            if K_sp == 'array':
                self.K = None
            else:
                self.K = self._grid.at_node[K_sp]
        elif type(K_sp) in (float, int):  # a float
            self.K = float(K_sp)
        elif len(K_sp) == self.grid.number_of_nodes:
            self.K = numpy.array(K_sp)
        else:
            raise TypeError('Supplied type of K_sp ' +
                            'was not recognised, or array was ' +
                            'not nnodes long!')

        if type(rainfall_intensity) is str:
            if rainfall_intensity == 'array':
                self.r_i = None
            else:
                self.r_i = self._grid.at_node[rainfall_intensity]
        elif type(rainfall_intensity) in (float, int):  # a float
            self.r_i = float(rainfall_intensity)
        elif len(rainfall_intensity) == self.grid.number_of_nodes:
            self.r_i = numpy.array(rainfall_intensity)
        else:
            raise TypeError('Supplied type of rainfall_' +
                            'intensity was not recognised, or array was ' +
                            'not nnodes long!')

        # We now forbid changing of the field name
        if 'value_field' in kwds.keys():
            raise ValueError('This component can no longer support variable' +
                             'field names. Use "topographic__elevation".')

    def erode(self, grid_in, dt=None, K_if_used=None, flooded_nodes=None):
        """
        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power
        scheme.

        This driving method for this component is now superceded by the new,
        standardized wrapper :func:`run_one_step`, but is retained for
        back compatibility.

        Set 'K_if_used' as a field name or nnodes-long array if you set K_sp as
        'array' during initialization.

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

        Returns
        -------
        grid
            A reference to the grid.
        """
        upstream_order_IDs = self._grid['node']['flow__upstream_node_order']
        z = self._grid['node']['topographic__elevation']
        defined_flow_receivers = numpy.not_equal(self._grid['node'][
            'flow__link_to_receiver_node'], UNDEFINED_INDEX)
        flow_link_lengths = self._grid._length_of_link_with_diagonals[
            self._grid['node']['flow__link_to_receiver_node'][
                defined_flow_receivers]]

        # make arrays from input the right size
        if type(self.K) is numpy.ndarray:
            K_here = self.K[defined_flow_receivers]
        else:
            K_here = self.K
        if type(self.r_i) is numpy.ndarray:
            r_i_here = self.r_i[defined_flow_receivers]
        else:
            r_i_here = self.r_i

        if dt is None:
            dt = self.dt
        assert dt is not None, ('Fastscape component could not find a dt to ' +
                                'use. Pass dt to the run_one_step() method.')

        if self.K is None:  # "old style" setting of array
            assert K_if_used is not None
            self.K = K_if_used

        numpy.power(self._grid['node']['drainage_area'], self.m,
                    out=self.A_to_the_m)
        self.alpha[defined_flow_receivers] = r_i_here**self.m * K_here * dt * \
            self.A_to_the_m[defined_flow_receivers] / flow_link_lengths

        flow_receivers = self._grid['node']['flow__receiver_node']
        n_nodes = upstream_order_IDs.size
        alpha = self.alpha

        # Handle flooded nodes, if any (no erosion there)
        if flooded_nodes is not None:
            alpha[flooded_nodes] = 0.
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            alpha[reversed_flow] = 0.

        self.alpha_by_flow_link_lengthtothenless1[
            defined_flow_receivers] = (alpha[defined_flow_receivers] /
                                       flow_link_lengths**(self.n - 1.))
        alpha_divided = self.alpha_by_flow_link_lengthtothenless1
        n = float(self.n)
        threshdt = self.thresholds * dt
        if type(self.thresholds) is float:
            from .cfuncs import erode_with_link_alpha_fixthresh
            erode_with_link_alpha_fixthresh(upstream_order_IDs, flow_receivers,
                                            threshdt, alpha_divided, n, z)
        else:
            from .cfuncs import erode_with_link_alpha_varthresh
            erode_with_link_alpha_varthresh(upstream_order_IDs, flow_receivers,
                                            threshdt, alpha_divided, n, z)
            # # This replicates the cython for testing:
            # for i in range(upstream_order_IDs.size):
            #     src_id = upstream_order_IDs[i]
            #     dst_id = flow_receivers[src_id]
            #     thresh = threshdt[i]
            #     if src_id != dst_id:
            #         next_z = z[src_id]
            #         prev_z = 0.
            #         while True:
            #         #for j in xrange(2):
            #             z_diff = next_z - z[dst_id]
            #             f = alpha_divided[src_id] * pow(z_diff, n - 1.)
            #             # if z_diff -> 0, pow -> nan (in reality, inf)
            #             # print (f, prev_z, next_z, z_diff, z[dst_id])
            #             next_z = next_z - ((next_z - z[src_id] + (
            #                 f*z_diff - thresh).clip(0.)) / (1. + n * f))
            #             if next_z < z[dst_id]:
            #                 next_z = z[dst_id] + 1.e-15
            #                 # ^maintain connectivity
            #             if next_z != 0.:
            #                 if (numpy.fabs((next_z - prev_z)/next_z) <
            #                     1.48e-08) or (n == 1.):
            #                     break
            #             else:
            #                 break
            #             prev_z = next_z
            #         if next_z < z[src_id]:
            #             z[src_id] = next_z

        return self._grid

    def run_one_step(self, dt, flooded_nodes=None, **kwds):
        """
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
        """
        self.erode(grid_in=self._grid, dt=dt, flooded_nodes=flooded_nodes)
