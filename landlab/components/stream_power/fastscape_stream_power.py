#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""
from __future__ import print_function

import numpy
from landlab import Component

from landlab.utils.decorators import use_file_name_or_kwds

from landlab import BAD_INDEX_VALUE as UNDEFINED_INDEX

from .cfuncs import (brent_method_erode_fixed_threshold,
                     brent_method_erode_variable_threshold)


class FastscapeEroder(Component):
    '''
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
    rainfall_intensity : float; optional
        Modifying factor on drainage area to convert it to a true water
        volume flux in (m/time). i.e., E = K * (r_i*A)**m * S**n. For a time
        varying rainfall intensity, pass rainfall_intensity_if_used to
        `run_one_step`. For a spatially variable rainfall, use the
        StreamPowerEroder component.

    Examples
    --------
    >>> import numpy as np
    >>> from landlab import RasterModelGrid
    >>> from landlab import CLOSED_BOUNDARY, FIXED_VALUE_BOUNDARY
    >>> from landlab.components import FlowRouter, FastscapeEroder
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
    ...                       rainfall_intensity=2.)
    >>> fr3.run_one_step()
    >>> sp3.run_one_step(1.)
    >>> z.reshape((3, 7))[1, :]  # doctest: +NORMALIZE_WHITESPACE
    array([  0.        ,   0.0647484 ,   0.58634455,   2.67253503,
             8.49212152,  20.92606987,  36.        ])
    >>> previous_z = z.copy()
    >>> sp3.run_one_step(1., rainfall_intensity_if_used=0.)
    >>> np.allclose(z, previous_z)
    True
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

        try:
            self.grid._diagonal_links_at_node  # calc number of diagonal links
        except AttributeError:
            pass  # was not a raster

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
        elif (type(K_sp) is numpy.ndarray and
              len(K_sp) == self.grid.number_of_nodes):
            self.K = K_sp
        else:
            raise TypeError('Supplied type of K_sp ' +
                            'was not recognised, or array was ' +
                            'not nnodes long!')

        if type(rainfall_intensity) is str:
            raise ValueError('This component can no longer handle ' +
                             'spatially variable rainfall. Use ' +
                             'StreamPowerEroder.')
            if rainfall_intensity == 'array':
                self._r_i = None
            else:
                self._r_i = self._grid.at_node[rainfall_intensity]
        elif type(rainfall_intensity) in (float, int):  # a float
            self._r_i = float(rainfall_intensity)
        elif len(rainfall_intensity) == self.grid.number_of_nodes:
            raise ValueError('This component can no longer handle ' +
                             'spatially variable rainfall. Use ' +
                             'StreamPowerEroder.')
            self._r_i = numpy.array(rainfall_intensity)
        else:
            raise TypeError('Supplied type of rainfall_' +
                            'intensity was not recognised!')

        # We now forbid changing of the field name
        if 'value_field' in kwds.keys():
            raise ValueError('This component can no longer support variable' +
                             'field names. Use "topographic__elevation".')

    def erode(self, grid_in, dt=None, K_if_used=None, flooded_nodes=None,
              rainfall_intensity_if_used=None):
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
        rainfall_intensity_if_used : float or None (optional)
            Supply to drive this component with a time-varying spatially
            constant rainfall.

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
        if rainfall_intensity_if_used is not None:
            assert type(rainfall_intensity_if_used) in (float, numpy.float64,
                                                        int)
            r_i_here = float(rainfall_intensity_if_used)
        else:
            r_i_here = self._r_i

        if dt is None:
            dt = self.dt
        assert dt is not None, ('Fastscape component could not find a dt to ' +
                                'use. Pass dt to the run_one_step() method.')

        if self.K is None:  # "old style" setting of array
            assert K_if_used is not None
            self.K = K_if_used

        n = float(self.n)
        numpy.power(self._grid['node']['drainage_area'], self.m,
                    out=self.A_to_the_m)
        self.alpha[defined_flow_receivers] = (
            r_i_here**self.m * K_here * dt * self.A_to_the_m[
                defined_flow_receivers] / (flow_link_lengths**self.n))

        flow_receivers = self._grid['node']['flow__receiver_node']
        alpha = self.alpha

        # Handle flooded nodes, if any (no erosion there)
        if flooded_nodes is not None:
            alpha[flooded_nodes] = 0.
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            alpha[reversed_flow] = 0.

        threshsdt = self.thresholds * dt

        # solve using Brent's Method in Cython for Speed
        if type(self.thresholds) == float:
            brent_method_erode_fixed_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z)
        else:
            brent_method_erode_variable_threshold(
                upstream_order_IDs, flow_receivers, threshsdt, alpha, n, z)

        return self._grid

    def run_one_step(self, dt, flooded_nodes=None,
                     rainfall_intensity_if_used=None, **kwds):
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
        rainfall_intensity_if_used : float or None (optional)
            Supply to drive this component with a time-varying spatially
            constant rainfall.
        """
        self.erode(grid_in=self._grid, dt=dt, flooded_nodes=flooded_nodes,
                   rainfall_intensity_if_used=rainfall_intensity_if_used)
