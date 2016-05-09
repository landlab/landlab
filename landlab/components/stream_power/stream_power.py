# -*- coding: utf-8 -*-
from __future__ import print_function

import numpy as np
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY, Component

from landlab.core.model_parameter_dictionary import MissingKeyError, \
    ParameterValueError
from landlab.field.scalar_data_fields import FieldError
from landlab.grid.base import BAD_INDEX_VALUE
from landlab.utils.decorators import use_file_name_or_kwds


class StreamPowerEroder(Component):
    """Erode where channels are.

    Implemented as:

    .. math::
        E = K A^m S^n - sp_{crit},

    and if :math:`E < 0`, :math:`E = 0`.

    If ``use_W`` is declared and ``True``, the module instead implements:

    .. math::
        E = K A^m S^n / W - sp_{crit}

    DEJH Sept 2013, major modifications Sept 14.

    NB: If you want spatially or temporally variable runoff, pass the
    runoff values at each pixel to the flow router using the input argument
    *use_Q*.

    Construction::

        StreamPowerEroder(grid, K_sp=None, threshold_sp=0., sp_type='set_mn',
                          m_sp=0.5, n_sp=1., a_sp=None, b_sp=None, c_sp=None,
                          use_W=None, use_Q=None):

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

    _name = 'StreamPowerEroder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'flow_receiver',
        'upstream_node_order',
        'topographic__steepest_slope'
    )

    _output_var_names = (
        'topographic__elevation',
        'stream_power_erosion'
    )

    _var_units = {
        'topographic__elevation': 'm',
        'drainage_area': 'm**2',
        'flow_receiver': '-',
        'upstream_node_order': '-',
        'topographic__steepest_slope': '-',
        'stream_power_erosion': 'variable'
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'drainage_area': 'node',
        'flow_receiver': 'node',
        'upstream_node_order': 'node',
        'topographic__steepest_slope': 'node',
        'stream_power_erosion': 'node'
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'flow_receiver':
            'Node array of receivers (node that receives flow from current '
            'node)',
        'upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'topographic__steepest_slope':
            'Node array of steepest *downhill* slopes',
        'stream_power_erosion':
            ('The value dt*K*A**m*S**n. Note the incorporation of time, and ' +
             'that any threshold is NOT included in this value.')
    }

    @use_file_name_or_kwds
    def __init__(self, grid, K_sp=None, threshold_sp=0., sp_type='set_mn',
                 m_sp=0.5, n_sp=1., a_sp=None, b_sp=None, c_sp=None,
                 use_W=None, use_Q=None, **kwds):
        self._grid = grid
        self.fraction_gradient_change = 1.
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links+1)
        # ^needs to be filled with values in execution
        self.count_active_links = np.zeros_like(
            self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1

        active_nodes = grid.status_at_node != 4
        self._K_unit_time = np.empty(active_nodes.sum(), dtype=float)
        self.use_K = False  # grandfathered in; only if K_sp == 'array'
        if type(K_sp) is np.ndarray:
            self._K_unit_time[:] = K_sp[active_nodes]
        else:
            try:
                self._K_unit_time.fill(K_sp)
            except ValueError:  # could not cast => was a str
                if K_sp == 'array':
                    self.use_K = True
                else:
                    self._K_unit_time = grid.at_node[K_sp]

        assert float(threshold_sp) >= 0.
        self.sp_crit = float(threshold_sp)
        if threshold_sp != 0.:
            self.set_threshold = True
            # ^flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
        else:
            self.set_threshold = False
        try:
            self.tstep = kwds['dt']
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
        if sp_type is 'set_mn':
            assert (float(m_sp) >= 0.) and (float(n_sp) >= 0.), \
                "m and n must be positive"
            self._m = float(m_sp)
            self._n = float(n_sp)
            assert ((a_sp is None) and (b_sp is None) and (c_sp is None)), (
                "If sp_type is 'set_mn', do not pass values for a, b, or c!")
        else:
            assert sp_type in ('Total', 'Unit', 'Shear_stress'), (
                "sp_type not recognised. It must be 'set_mn', 'Total', " +
                "'Unit', or 'Shear_stress'.")
            assert (m_sp == 0.5 and n_sp == 1.), \
                "Do not set m and n if sp_type is not 'set_mn'!"
            assert float(a_sp) >= 0., "a must be positive"
            self._a = float(a_sp)
            if b_sp is not None:
                assert float(b_sp) >= 0., "b must be positive"
                self._b = float(b_sp)
            else:
                assert self.use_W, "b was not set"
                self._b = 0.
            if c_sp is not None:
                assert float(c_sp) >= 0., "c must be positive"
                self._c = float(c_sp)
            else:
                assert self.use_Q, "c was not set"
                self._c = 1.
            if self._type == 'Total':
                self._n = self._a
                self._m = self._a*self._c  # ==_a if use_Q
            elif self._type == 'Unit':
                self._n = self._a
                self._m = self._a*self._c*(1.-self._b)
                # ^ ==_a iff use_Q&use_W etc
            elif self._type == 'Shear_stress':
                self._m = 2.*self._a*self._c*(1.-self._b)/3.
                self._n = 2.*self._a/3.
            else:
                raise MissingKeyError('Not enough information was provided ' +
                                      'on the exponents to use!')
        # m and n will always be set, but care needs to be taken to include Q
        # and W directly if appropriate

        self.stream_power_erosion = grid.zeros(centering='node')
        grid.add_zeros('stream_power_erosion', at='node')

    def initialize(self, grid, params_file):
        r"""
        .. note:: deprecated
            Use the ``__init__`` method directly.

        params_file is the name of the text file containing the parameters
        needed for this stream power component.

        Module erodes where channels are, implemented as

        .. math::
            E = K A^m S^n - sp_{crit}

        and if ``E < 0``, ``E = 0``.

        If ``use_W`` is declared and ``True``, the module instead implements:

        .. math::
            E = K A^m S^n / W - sp_{crit}

        **Parameters for input file**:

        Obligatory:

        *  K_sp : positive float, the prefactor. This is defined per unit
           time, not per tstep. Type the string 'array' to cause the
           component's erode method to look for an array of values of K
           (see documentation for 'erode').

        Alternatives:

        *either*:

        *  ``m_sp``: positive float, the power on ``A`` *and*
        *  ``n_sp``: positive float, the power on ``S``

        *or*

        *  ``sp_type``: str. Must be one of ``'Total'``, ``'Unit'``, or
           ``'Shear_stress'``.

        and (following Whipple & Tucker 1999)

        *  ``a_sp``: +ve float. The power on the SP/shear term to get the
           erosion rate.
        *  ``b_sp``: +ve float. The power on discharge to get width,
           "hydraulic
           geometry". Unnecessary if sp_type='Total'.
        *  ``c_sp``: +ve float. The power on area to get discharge, "basin
           hydology".

        If:

        *  ``'Total'``, :math:`m = a c, n = a`.
        *  ``'Unit'``, :math:`m = a c (1 - b), n = a`.
        *  ``'Shear_stress'``, :math:`m =2 a c (1 - b) / 3, n = 2 a / 3`.


        Options:

        *  ``threshold_sp``: +ve float; the threshold sp_crit. Defaults to 0.
           This threshold is assumed to be in "stream power" units, i.e.,
           if 'Shear_stress', the value should be tau**a.
        *  ``dt``: +ve float. If set, this is the fixed timestep for this
           component. Can be overridden easily as a parameter in erode().
           If not set (default), this parameter MUST be set in erode().
        *  ``use_W``: ``bool``; if True, component will look for node-centered
           data describing channel width in ``grid.at_node['channel_width']``,
           and use it to implement incision ~ stream power per unit width.
           Defaults to ``False``. If you set sp_m and sp_n, follows the
           equation given above. If you set sp_type, it will be ignored if
           'Total', but used directly if you want 'Unit' or
           'Shear_stress'.
        *  ``use_Q``: ``Bool``. If true, the equation becomes

           .. math::
               E = K Q^m S^n

           Effectively sets c=1 in Wh&T's 1999 derivation, if you are
           setting m and n through a, b, and c.
        """
        self._grid = grid
        self.fraction_gradient_change = 1.
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links+1)
        # ^needs to be filled with values in execution
        self.count_active_links = np.zeros_like(
            self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1
        inputs = ModelParameterDictionary(params_file)
        try:
            self._K_unit_time = np.full((grid.status_at_node != 4).sum(),
                                        inputs.read_float('K_sp'))
        except ParameterValueError:  # it was a string
            self.use_K = True
        else:
            self.use_K = False
        try:
            self.sp_crit = inputs.read_float('threshold_sp')
            self.set_threshold = True
            # ^flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
            # print("Found a threshold to use: ", self.sp_crit)
        except MissingKeyError:
            self.sp_crit = 0.
            self.set_threshold = False
        try:
            self.tstep = inputs.read_float('dt')
        except MissingKeyError:
            pass
        try:
            self.use_W = inputs.read_bool('use_W')
        except MissingKeyError:
            self.use_W = False
        try:
            self.use_Q = inputs.read_bool('use_Q')
        except MissingKeyError:
            self.use_Q = False
        try:
            self._m = inputs.read_float('m_sp')
        except MissingKeyError:
            self._type = inputs.read_string('sp_type')
            self._a = inputs.read_float('a_sp')
            try:
                self._b = inputs.read_float('b_sp')
            except MissingKeyError:
                if self.use_W:
                    self._b = 0.
                else:
                    raise NameError('b was not set')
            try:
                self._c = inputs.read_float('c_sp')
            except MissingKeyError:
                if self.use_Q:
                    self._c = 1.
                else:
                    raise NameError('c was not set')
            if self._type == 'Total':
                self._n = self._a
                self._m = self._a*self._c  # ==_a if use_Q
            elif self._type == 'Unit':
                self._n = self._a
                self._m = self._a*self._c*(1.-self._b)
                # ^ ==_a iff use_Q&use_W etc
            elif self._type == 'Shear_stress':
                self._m = 2.*self._a*self._c*(1.-self._b)/3.
                self._n = 2.*self._a/3.
            else:
                raise MissingKeyError('Not enough information was provided ' +
                                      'on the exponents to use!')
        else:
            self._n = inputs.read_float('n_sp')
        # m and n will always be set, but care needs to be taken to include Q
        # and W directly if appropriate

        self.stream_power_erosion = grid.zeros(centering='node')

    def erode(self, grid, dt, node_elevs='topographic__elevation',
              node_drainage_areas='drainage_area',
              flow_receiver='flow_receiver',
              node_order_upstream='upstream_node_order',
              slopes_at_nodes='topographic__steepest_slope',
              link_node_mapping='links_to_flow_receiver',
              link_slopes=None, slopes_from_elevs=None,
              W_if_used=None, Q_if_used=None, K_if_used=None,
              flooded_nodes=None):
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

        node_elevs : str or ndarray, optional
            Elevations on the grid, either a field string or nnodes-long array.

        node_drainage_areas: str or ndarray, optional
            Tells the component where to look for the drainage area values.
            Change to another string to override which grid field the
            component looks at, or pass a nnodes-long array of drainage
            areas values directly instead.

        flow_receiver, node_order_upstream : str or ndarray, optional
            The downstream node to which each node flows and the ordering of
            the nodes in the network starting at the outlet, respectively,
            are both necessary as inputs to allow stability testing.

            If you already have slopes defined at nodes on the grid, pass them
            to the component with *slopes_at_nodes*. The same syntax is
            expected: string gives a name in the grid fields, an array gives
            values direct.

            Alternatively, set *link_slopes* (and *link_node_mapping*) if this
            data
            is only available at links. 'topographic__derivative_of_elevation'
            is the default field name for link slopes. Override this name by
            setting the variable as the appropriate string, or override use of
            grid fields altogether by passing an array. *link_node_mapping*
            controls how the component maps these link values onto the arrays.
            We assume there is always a 1:1 mapping (pass the values already
            projected onto the nodes using slopes_at_nodes if not). Other
            components, e.g., flow_routing.route_flow_dn, may provide the
            necessary outputs to make the mapping easier: e.g., just pass
            'links_to_flow_receiver' from that module (the default name). If
            the component cannot find an existing mapping through this
            parameter, it will derive one on the fly, at considerable cost of
            speed (see on-screen reports).

        slopes_from_elevs : str, optional
            Allows the module to create gradients internally
            from elevations rather than have them provided. Set to True to
            force the component to look for the data in the location specified
            by node_elevs. Using this option is
            considerably slower than any of the alternatives, as it also has to
            calculate the link_node_mapping from stratch each time.

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
            modifies grid elevation fields to reflect updates; creates and
            maintains ``grid.at_node['stream_power_erosion']``. Note the value
            stream_power_erosion is not an excess stream power; any specified
            erosion threshold is not incorporated into it.
        """
        active_nodes = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]

        if W_if_used is not None:
            assert self.use_W, ("Widths were provided, but you didn't set " +
                                "the use_W flag in your input file! " +
                                "Aborting...")
            assert self._W is None, ("Do not pass W to the run method " +
                                     "if you also set them at initialization!")

        if Q_if_used is not None:
            assert self.use_Q, ("Discharges were provided, but you didn't " +
                                "set the use_Q flag in your input file! " +
                                "Aborting...")
            assert self._Q is None, ("Do not pass Q to the run method " +
                                     "if you also set them at initialization!")

        if K_if_used is not None:
            assert self.use_K, ("An array of erodabilities was provided, " +
                                "but you didn't set K_sp to 'array' in your " +
                                "input file! Aborting...")
            try:
                self._K_unit_time = grid.at_node[K_if_used][active_nodes]
            except TypeError:
                self._K_unit_time = K_if_used[active_nodes]

        if type(node_elevs) is str:
            node_z = grid.at_node[node_elevs]
        else:
            node_z = node_elevs

        # Perform check on whether we use grid or direct fed data:
        try:
            self.slopes = grid.at_node[slopes_at_nodes]
        except TypeError:
            if type(slopes_at_nodes) is np.ndarray:
                self.slopes = slopes_at_nodes
            else:
                raise TypeError('slopes_at_nodes input not recognised')
        except FieldError:
            if slopes_from_elevs is True:
                S_links = (node_z[grid.node_at_link_tail] -
                           node_z[grid.node_at_link_head])/grid.length_of_link
            else:
                if link_slopes:
                    if type(link_slopes) is str:
                        S_links = grid.at_link[link_slopes]
                    else:
                        S_links = link_slopes
                else:
                    try:
                        S_links = grid.at_link['topographic__derivative_of_' +
                                               'elevation']
                    except FieldError:
                        # ...for back compatibility
                        S_links = grid.at_link['planet_surface__derivative_' +
                                               'of_elevation']

            # put the slopes onto the nodes
            try:
                self.slopes = S_links[grid.at_node[link_node_mapping]]
            except TypeError:
                try:
                    self.slopes = S_links[link_node_mapping]
                except IndexError:
                    # need to do the mapping on the fly.
                    # we're going to use the max slope (i.e., + or -) of *all*
                    # adjacent nodes.
                    # This isn't ideal. It should probably just be the outs...
                    # i.e., np.max(self.link_S_with_trailing_blank[grid.node_
                    # outlinks] AND -self.link_S_with_trailing_blank[grid.
                    # node_inlinks])
                    self.link_S_with_trailing_blank[:-1] = S_links
                    self.slopes = np.amax(np.fabs(
                        self.link_S_with_trailing_blank[grid.links_at_node.T]),
                        axis=0)

        if type(node_drainage_areas) is str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        if type(flow_receiver) is str:
            flow_receiver = grid.at_node[flow_receiver]

        if type(node_order_upstream) is str:
            node_order_upstream = grid.at_node[node_order_upstream]

        # Disable incision in flooded nodes, as appropriate
        if flooded_nodes is not None:
            self._K_unit_time[flooded_nodes[active_nodes]] = 0.

        # Operate the main function:
        if self.use_W is False and self.use_Q is False:  # normal case
            stream_power_active_nodes = (self._K_unit_time * dt *
                                         node_A[active_nodes]**self._m *
                                         self.slopes[active_nodes]**self._n)
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
                stream_power_active_nodes = (self._K_unit_time * dt *
                                             Q_direct[active_nodes]**self._m *
                                             self.slopes[active_nodes] **
                                             self._n / W[active_nodes])
            else:  # just W to be used
                stream_power_active_nodes = (self._K_unit_time * dt *
                                             node_A[active_nodes]**self._m *
                                             self.slopes[active_nodes] **
                                             self._n / W[active_nodes])
        else:  # just use_Q
            if self._Q is None:
                try:
                    Q_direct = grid.at_node[Q_if_used]
                except TypeError:
                    assert type(Q_if_used) in (np.ndarray, list)
                    Q_direct = Q_if_used
            else:
                Q_direct = self._Q
            stream_power_active_nodes = (self._K_unit_time * dt *
                                         Q_direct[active_nodes]**self._m *
                                         self.slopes[active_nodes]**self._n)

        # Note that we save "stream_power_erosion" incorporating both K and a.
        # Most definitions would need this value /K then **(1/a) to give actual
        # stream power (unit, total, whatever), and it does not yet include the
        # threshold
        self.stream_power_erosion[active_nodes] = stream_power_active_nodes
        grid.at_node['stream_power_erosion'][:] = self.stream_power_erosion
        erosion_increment = (self.stream_power_erosion - self.sp_crit).clip(0.)

        # this prevents any node from incising below any node downstream of it
        # we have to go in upstream order in case our rate is so big we impinge
        # on baselevels > 1 node away

        elev_dstr = node_z[flow_receiver]
        # ^we substract erosion_increment[flow_receiver] in the loop, as it
        # can update

        method = 'cython'
        if method == 'cython':
            from .cfuncs import erode_avoiding_pits

            erode_avoiding_pits(node_order_upstream, flow_receiver, node_z,
                                erosion_increment)
        else:
            for i in node_order_upstream:
                elev_this_node_before = node_z[i]
                elev_this_node_after = (elev_this_node_before -
                                        erosion_increment[i])
                elev_dstr_node_after = (elev_dstr[i] -
                                        erosion_increment[flow_receiver[i]])
                if elev_this_node_after < elev_dstr_node_after:
                    erosion_increment[i] = (elev_this_node_before -
                                            elev_dstr_node_after)*0.999999
                # ^we add a tiny elevation excess to prevent the module from
                # ever totally severing its own flow paths
        # clip the erosion increments one more time to remove regatives
        # introduced by any pit filling algorithms or the above procedure:
        node_z -= erosion_increment.clip(0.)

        self._grid = grid

        return grid, node_z, self.stream_power_erosion

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
