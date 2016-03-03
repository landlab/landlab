#! /usr/env/python

"""
This module attempts to "component-ify" GT's Fastscape stream power erosion.
Created DEJH, March 2014.
"""
from __future__ import print_function

import numpy
from landlab import ModelParameterDictionary, Component
from landlab.core.model_parameter_dictionary import MissingKeyError, \
    ParameterValueError
from landlab.utils.decorators import use_file_name_or_kwds
from landlab.field.scalar_data_fields import FieldError
from scipy.optimize import newton, fsolve

UNDEFINED_INDEX = numpy.iinfo(numpy.int32).max


class FastscapeEroder(Component):
    '''
    This class uses the Braun-Willett Fastscape approach to calculate the
    amount of erosion at each node in a grid, following a stream power
    framework.

    On initialization, it takes *grid*, a reference to a ModelGrid, and
    *input_stream*, a string giving the filename (and optionally, path) of the
    required input file.

    It needs to be supplied with the key variables:

        *K_sp*

        *m_sp*

    ...which it will draw from the supplied input file. *n_sp*  can be any
    value ~ 0.5<n_sp<4., but note that performance will be EXTREMELY degraded
    if n<1.

    If you want to supply a spatial variation in K, set K_sp to the string
    'array', and pass a field name or array to the erode method's K_if_used
    argument.

    *dt*, *rainfall_intensity*, and *value_field* are optional variables.

    *dt* is a fixed timestep, and *rainfall_intensity* is a parameter which
    modulates K_sp (by a product, r_i**m_sp) to reflect the direct influence of
    rainfall intensity on erosivity. *value_field* is a string giving the name
    of the field containing the elevation data in the grid. It defaults to
    'topographic__elevation' if not supplied.

    This module assumes you have already run
    :func:`landlab.components.flow_routing.route_flow_dn.FlowRouter.route_flow`
    in the same timestep. It looks for 'upstream_node_order',
    'links_to_flow_receiver', 'drainage_area', 'flow_receiver', and
    'topographic__elevation' at the nodes in the grid. 'drainage_area' should
    be in area upstream, not volume (i.e., set runoff_rate=1.0 when calling
    FlowRouter.route_flow).

    The primary method of this class is :func:`erode`.
    '''

    _name = 'FastscapeEroder'

    _input_var_names = (
        'topographic__elevation',
        'drainage_area',
        'links_to_flow_receiver',
        'upstream_node_order',
        'flow_receiver',
    )

    _output_var_names = (
        'topographic__elevation',
    )

    _var_units = {
        'topographic__elevation': 'm',
        'drainage_area': 'm**2',
        'links_to_flow_receiver': '-',
        'upstream_node_order': '-',
        'flow_receiver': '-',
    }

    _var_mapping = {
        'topographic__elevation': 'node',
        'drainage_area': 'node',
        'links_to_flow_receiver': 'node',
        'upstream_node_order': 'node',
        'flow_receiver': 'node',
    }

    _var_doc = {
        'topographic__elevation': 'Land surface topographic elevation',
        'drainage_area':
            "Upstream accumulated surface area contributing to the node's "
            "discharge",
        'links_to_flow_receiver':
            'ID of link downstream of each node, which carries the discharge',
        'upstream_node_order':
            'Node array containing downstream-to-upstream ordered list of '
            'node IDs',
        'flow_receiver':
            'Node array of receivers (node that receives flow from current '
            'node)',
    }

    @use_file_name_or_kwds
    def __init__(self, grid, K_sp=None, m_sp=0.5, n_sp=1.,
                 rainfall_intensity=1., **kwds):
        """Initialize the Fastscape stream power component.

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

        self.K = K_sp
        self.m = m_sp
        self.n = n_sp
        self.r_i = rainfall_intensity
        self.use_K = True
        self.use_ri = True  # these ones overwritten below
        self.dt = None  # this is a dummy to allow the old-style component
        # style to still work

        # make storage variables
        self.A_to_the_m = grid.zeros(at='node')
        self.alpha = grid.empty(at='node')
        self.alpha_by_flow_link_lengthtothenless1 = numpy.empty_like(
                                                        self.alpha)

        self._grid.diagonal_links_at_node()  # calc number of diagonal links

        if self.n != 1.:
            self.nonlinear_flag = True
            if self.n < 1.:
                print("***WARNING: With n<1 performance of the Fastscape" +
                      " algorithm is slow!***")
        else:
            self.nonlinear_flag = False

        if self.K is None:
            raise ValueError('K_sp must be set as a float, node array, or ' +
                             'field name. It was None.')

        def func_for_newton(x, last_step_elev, receiver_elev,
                            alpha_by_flow_link_lengthtothenless1, n):
            y = (x - last_step_elev + alpha_by_flow_link_lengthtothenless1 *
                 (x - receiver_elev)**n)
            return y

        def func_for_newton_diff(x, last_step_elev, receiver_elev,
                                 alpha_by_flow_link_lengthtothenless1, n):
            y = (1. + n * alpha_by_flow_link_lengthtothenless1 *
                 (x - receiver_elev)**(n - 1.))
            return y

        self.func_for_newton = func_for_newton
        self.func_for_newton_diff = func_for_newton_diff

        # now handle the inputs that could be float, array or field name:
        input_to_property = {k_sp: self.K, rainfall_intensity: self.r_i}
        input_to_flag = {k_sp: self.use_K, rainfall_intensity: self.use_ri}
        for input_param in input_to_property.keys():
            if type(input_param) is str:
                if input_param is 'array':
                    input_to_property[input_param] = None
                else:
                    input_to_property[input_param] = self._grid.at_node[
                        input_param]
            elif type(input_param) in (float, int):  # a float
                input_to_flag[input_param] = False
            elif len(input_param) == self._grid.number_of_nodes:
                pass  # array of right length
            else:
                raise TypeError('Supplied type of either K or rainfall ' +
                                'intensity was not recognised, or array was ' +
                                'not nnodes long!')

        # We now forbid changing of the field name
        if 'value_field' in kwds.keys():
            raise ValueError('This component can no longer support variable' +
                             'field names. Use "topographic__elevation".')

    # this should now be redundant, but retained for back compatibility
    def gear_timestep(self, dt_in, rainfall_intensity_in=None):
        self.dt = dt_in
        if rainfall_intensity_in is not None:
            self.r_i = rainfall_intensity_in
        return self.dt, self.r_i

    def erode(self, grid_in, dt=None, K_if_used=None, flooded_nodes=None):
        """
        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power
        scheme.

        This driving method for this component is now superceded by the new,
        standardized wrapper :func:`run_one_timestep`, but is retained for
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
        dt : float (optional)
            Time-step size. If you call :func:`gear_timestep`, that method will
            supercede any value supplied here.
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
        upstream_order_IDs = self._grid['node']['upstream_node_order']
        z = self._grid['node']['topographic__elevation']
        defined_flow_receivers = numpy.not_equal(self._grid['node'][
            'links_to_flow_receiver'], UNDEFINED_INDEX)
        flow_link_lengths = self._grid.link_length[self._grid['node'][
            'links_to_flow_receiver'][defined_flow_receivers]]

        # make arrays from input the right size
        if self.use_K:
            K_here = self.K[defined_flow_receivers]
        else:
            K_here = self.K
        if self.use_ri:
            r_i_here = self.r_i[defined_flow_receivers]
        else:
            r_i_here = self.r_i

        if self.dt is not None:
            dt = self.dt
        else:
            assert dt is not None
            dt = dt_in

        if self.K is None:  # "old style" setting of array
            assert K_if_used is not None
            self.K = K_if_used

        numpy.power(self._grid['node']['drainage_area'], self.m,
                    out=self.A_to_the_m)
        self.alpha[defined_flow_receivers] = r_i_here**self.m * K_here * dt * \
            self.A_to_the_m[defined_flow_receivers] / flow_link_lengths

        flow_receivers = self._grid['node']['flow_receiver']
        n_nodes = upstream_order_IDs.size
        alpha = self.alpha

        # Handle flooded nodes, if any (no erosion there)
        if flooded_nodes is not None:
            alpha[flooded_nodes] = 0.
        else:
            reversed_flow = z < z[flow_receivers]
            # this check necessary if flow has been routed across depressions
            alpha[reversed_flow] = 0.

        method = 'cython'

        if self.nonlinear_flag is False:  # n==1
            if method == 'cython':
                from .cfuncs import erode_with_alpha
                erode_with_alpha(upstream_order_IDs, flow_receivers, alpha, z)
            else:
                for i in upstream_order_IDs:
                    j = flow_receivers[i]
                    if i != j:
                        z[i] = (z[i] + alpha[i] * z[j])/(1.0 + alpha[i])
        else:  # general, nonlinear n case
            print('Non-linear')
            self.alpha_by_flow_link_lengthtothenless1[
                defined_flow_receivers] = (alpha[defined_flow_receivers] /
                                           flow_link_lengths**(self.n - 1.))
            alpha_divided = self.alpha_by_flow_link_lengthtothenless1
            n = float(self.n)
            if method == 'cython':
                if n < 1.:
                    # this is SLOOOOOOOOOOOW...
                    for i in upstream_order_IDs:
                        j = flow_receivers[i]
                        func_for_newton = self.func_for_newton
                        func_for_newton_diff = self.func_for_newton_diff
                        if i != j:
                            z[i] = fsolve(func_for_newton, z[i],
                                          args=(z[i], z[j], alpha_divided[i],
                                                n))
                else:
                    from .cfuncs import erode_with_link_alpha
                    erode_with_link_alpha(upstream_order_IDs, flow_receivers,
                                          alpha_divided, n, z)
            else:
                for i in upstream_order_IDs:
                    j = flow_receivers[i]
                    func_for_newton = self.func_for_newton
                    func_for_newton_diff = self.func_for_newton_diff
                    if i != j:
                        if n >= 1.:
                            z[i] = newton(func_for_newton, z[i],
                                          fprime=func_for_newton_diff,
                                          args=(z[i], z[j], alpha_divided[i],
                                                n), maxiter=10)
                        else:
                            z[i] = fsolve(func_for_newton, z[i],
                                          args=(z[i], z[j], alpha_divided[i],
                                                n))

        return self._grid

    def run_one_timestep(self, dt, flooded_nodes=None):
        """
        This method implements the stream power erosion, following the Braun-
        Willett (2013) implicit Fastscape algorithm. This should allow it to
        be stable against larger timesteps than an explicit stream power
        scheme.

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
