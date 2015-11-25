from __future__ import print_function

import numpy as np
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY

from landlab.core.model_parameter_dictionary import MissingKeyError, ParameterValueError
from landlab.field.scalar_data_fields import FieldError
from landlab.grid.base import BAD_INDEX_VALUE

class StreamPowerEroder(object):
    """
    This component is now verified stable for simple m,n specified, followup-to-
    Fastscape-flow-routing cases. Threshold appears stable.
    The more exciting cases (e.g., specifying a,b,c; forcing with W or Q) are
    untested, but should run.
    There is as yet no explicit stabilization check on the timestep. If your
    run destabilizes, try reducing dt.
    See, e.g., ./examples/simple_sp_driver.py

    DEJH Sept 2013, major modifications Sept 14.
    This component *should* run on any grid, but untested.
    """

    def __init__(self, grid, params):
        self.initialize(grid, params)

#This draws attention to a potential problem. It will be easy to have modules update z, but because noone "owns" the data, to forget to also update dz/dx...
#How about a built in grid utility that updates "derived" data (i.e., using only grid fns, e.g., slope, curvature) at the end of any given tstep loop?
#Or an explicit flagging system for all variables in the modelfield indicating if they have been updated this timestep. (Currently implemented)
#Or wipe the existance of any derived grid data at the end of a timestep entirely, so modules find they don't have it next timestep.

    def initialize(self, grid, params_file):
        '''
        params_file is the name of the text file containing the parameters
        needed for this stream power component.

        Module erodes where channels are, implemented as

        E = K * A**m * S**n - sp_crit,

        and if E<0, E=0.

        If 'use_W' is declared and True, the module instead implements:

        E = K * A**m * S**n / W - sp_crit

        ***Parameters for input file***
        OBLIGATORY:
            K_sp -> positive float, the prefactor. This is defined per unit
                time, not per tstep. Type the string 'array' to cause the
                component's erode method to look for an array of values of K
                (see documentation for 'erode').
        ALTERNATIVES:
        *either*
            m_sp -> positive float, the power on A
        and
            n_sp -> positive float, the power on S
        *or*
            sp_type -> String. Must be one of 'Total', 'Unit', or 'Shear_stress'.
        and (following Whipple & Tucker 1999)
            a_sp -> +ve float. The power on the SP/shear term to get the erosion
                rate.
            b_sp -> +ve float. The power on discharge to get width, "hydraulic
                geometry". Unnecessary if sp_type='Total'.
            c_sp -> +ve float. The power on area to get discharge, "basin
                hydology".
            ... If 'Total', m=a*c, n=a.
            ... If 'Unit', m=a*c*(1-b), n=a.
            ... If 'Shear_stress', m=2*a*c*(1-b)/3, n = 2*a/3.
        OPTIONS:
            threshold_sp -> +ve float; the threshold sp_crit. Defaults to 0.
                This threshold is assumed to be in "stream power" units, i.e.,
                if 'Shear_stress', the value should be tau**a.
            dt -> +ve float. If set, this is the fixed timestep for this
                component. Can be overridden easily as a parameter in erode().
                If not set (default), this parameter MUST be set in erode().
            use_W -> Bool; if True, component will look for node-centered data
                describing channel width in grid.at_node['channel_width'], and
                use it to implement incision ~ stream power per unit width.
                Defaults to False. If you set sp_m and sp_n, follows the
                equation given above. If you set sp_type, it will be ignored if
                'Total', but used directly if you want 'Unit' or 'Shear_stress'.
            use_Q -> Bool. If true, the equation becomes E=K*Q**m*S**n.
                Effectively sets c=1 in Wh&T's 1999 derivation, if you are
                setting m and n through a, b, and c.

        '''
        self.grid = grid
        self.fraction_gradient_change = 1.
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links+1) #needs to be filled with values in execution
        self.count_active_links = np.zeros_like(self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1
        inputs = ModelParameterDictionary(params_file)
        try:
            self._K_unit_time = inputs.read_float('K_sp')
        except ParameterValueError: #it was a string
            self.use_K = True
        else:
            self.use_K = False
        try:
            self.sp_crit = inputs.read_float('threshold_sp')
            self.set_threshold = True #flag for sed_flux_dep_incision to see if the threshold was manually set.
            print("Found a threshold to use: ", self.sp_crit)
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
                self._m = self._a*self._c #==_a if use_Q
            elif self._type == 'Unit':
                self._n = self._a
                self._m = self._a*self._c*(1.-self._b) #==_a iff use_Q&use_W etc
            elif self._type == 'Shear_stress':
                self._m = 2.*self._a*self._c*(1.-self._b)/3.
                self._n = 2.*self._a/3.
            else:
                raise MissingKeyError('Not enough information was provided on the exponents to use!')
        else:
            self._n = inputs.read_float('n_sp')
        #m and n will always be set, but care needs to be taken to include Q and W directly if appropriate

        self.stream_power_erosion = grid.zeros(centering='node')
        ##Flags for self-building of derived data:
        #self.made_link_gradients = False
        ##This will us the MPD once finalized
        ##Now perform checks for existance of needed data items:
        #try:
        #    _ = self.grid.at_link['planet_surface__derivative_of_elevation']
        #except FieldError:
        #    self.made_link_gradients = True


    def erode(self, grid, dt, node_elevs='topographic__elevation',
            node_drainage_areas='drainage_area',
            flow_receiver='flow_receiver',
            node_order_upstream='upstream_node_order',
            slopes_at_nodes='topographic__steepest_slope',
            link_node_mapping='links_to_flow_receiver',
            link_slopes=None, slopes_from_elevs=None,
            W_if_used=None, Q_if_used=None, K_if_used=None):
        """
        A simple, explicit implementation of a stream power algorithm.

        *grid* & *dt* are the grid object and timestep (float) respectively.

        *node_elevs* is the elevations on the grid, either a field string or
        nnodes-long array.

        *Node_drainage_areas* tells the component where to look for the drainage
        area values. Change to another string to override which grid field the
        component looks at, or pass a nnodes-long array of drainage areas values
        directly instead.

        *flow_receiver* and *node_order_upstream*, the downstream node to which
        each node flows and the ordering of the nodes in the network starting
        at the outlet, respectively,
        are both necessary as inputs to allow stability testing.

        If you already have slopes defined at nodes on the grid, pass them to
        the component with *slopes_at_nodes*. The same syntax is expected:
        string gives a name in the grid fields, an array gives values direct.

        Alternatively, set *link_slopes* (and *link_node_mapping*) if this data
        is only available at links. 'planet_surface__derivative_of_elevation'
        is the default field name for link slopes. Override this name by
        setting the variable as the appropriate string, or override use of
        grid fields altogether by passing an array. *link_node_mapping* controls
        how the component maps these link values onto the arrays. We assume
        there is always a 1:1 mapping (pass the values already projected onto
        the nodes using slopes_at_nodes if not). Other components, e.g.,
        flow_routing.route_flow_dn, may provide the necessary outputs to make
        the mapping easier: e.g., just pass 'links_to_flow_reciever' from that
        module (the default name). If the component cannot find an existing
        mapping through this parameter, it will derive one on the fly, at
        considerable cost of speed (see on-screen reports).

        *slopes_from_elevs* allows the module to create gradients internally
        from elevations rather than have them provided. Set to True to force
        the component to look for the data in the location specified by
        node_elevs. Using this option is
        considerably slower than any of the alternatives, as it also has to
        calculate the link_node_mapping from stratch each time.

        In both these cases, at present the mapping is to use the maximum
        slope of --any-- link attached to the node as the representative node
        slope. This is primarily for speed, but may be a good idea to modify
        later.

        *W_if_used* and *Q_if_used* must be provided if you set use_W and use_Q
        respectively in the component initialization. They can be either field
        names or nnodes arrays as in the other cases.

        NB: If you want spatially or temporally variable runoff, pass the
        runoff values at each pixel to the flow router, then pass discharges
        at each node using *Q_if_used* to this component.

        RETURNS (grid, modified_elevs, stream_power_erosion); modifies grid elevation
        fields to reflect updates; creates and maintains
        grid.at_node['stream_power_erosion']. Note the value stream_power_erosion
        is not an excess stream power; any specified erosion threshold is not
        incorporated into it.
        """
        active_nodes = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]

        if W_if_used!=None:
            assert self.use_W, "Widths were provided, but you didn't set the use_W flag in your input file! Aborting..."

        if Q_if_used!=None:
            assert self.use_Q, "Discharges were provided, but you didn't set the use_Q flag in your input file! Aborting..."

        if K_if_used!=None:
            assert self.use_K, "An array of erodabilities was provided, but you didn't set K_sp to 'array' in your input file! Aborting..."
            try:
                self._K_unit_time = grid.at_node[K_if_used][active_nodes]
            except TypeError:
                self._K_unit_time = K_if_used[active_nodes]

        if type(node_elevs)==str:
            node_z = grid.at_node[node_elevs]
        else:
            node_z = node_elevs

        #Perform check on whether we use grid or direct fed data:
        try:
            self.slopes = grid.at_node[slopes_at_nodes]
        except TypeError:
            if type(slopes_at_nodes)==np.ndarray:
                self.slopes = slopes_at_nodes
            else:
                raise TypeError('slopes_at_nodes input not recognised')
        except FieldError:
            if slopes_from_elevs==True:
                S_links = (node_z[grid.node_at_link_tail]-node_z[grid.node_at_link_head])/grid.link_length
            else:
                if link_slopes:
                    if type(link_slopes)==str:
                        S_links = grid.at_link[link_slopes]
                    else:
                        S_links = link_slopes
                else:
                    S_links = grid.at_link['planet_surface__derivative_of_elevation']

            #put the slopes onto the nodes
            try:
                self.slopes = S_links[grid.at_node[link_node_mapping]]
            except TypeError:
                try:
                    self.slopes = S_links[link_node_mapping]
                except IndexError:
                    #need to do the mapping on the fly.
                    #we're going to use the max slope (i.e., + or -) of *all* adjacent nodes.
                    #This isn't ideal. It should probably just be the outs...
                    #i.e., np.max(self.link_S_with_trailing_blank[grid.node_outlinks] AND -self.link_S_with_trailing_blank[grid.node_inlinks])
                    self.link_S_with_trailing_blank[:-1] = S_links
                    self.slopes = np.amax(np.fabs(self.link_S_with_trailing_blank[grid.node_links]),axis=0)

        if type(node_drainage_areas)==str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        if type(flow_receiver)==str:
            flow_receiver = grid.at_node[flow_receiver]

        if type(node_order_upstream)==str:
            node_order_upstream = grid.at_node[node_order_upstream]

        #Operate the main function:
        if self.use_W==False and self.use_Q==False: #normal case
            stream_power_active_nodes = self._K_unit_time * dt * node_A[active_nodes]**self._m * self.slopes[active_nodes]**self._n
        elif self.use_W:
            try:
                W = grid.at_node[W_if_used]
            except TypeError:
                W = W_if_used
            if self.use_Q: #use both Q and W direct
                try:
                    Q_direct = grid.at_node[Q_if_used]
                except TypeError:
                    Q_direct = Q_if_used
                stream_power_active_nodes = self._K_unit_time * dt * Q_direct[active_nodes]**self._m * self.slopes[active_nodes]**self._n / W
            else: #just W to be used
                stream_power_active_nodes = self._K_unit_time * dt * node_A[active_nodes]**self._m * self.slopes[active_nodes]**self._n / W
        else: #just use_Q
            try:
                Q_direct = grid.at_node[Q_if_used]
            except TypeError:
                assert type(Q_if_used) in (np.ndarray, list)
                Q_direct = Q_if_used
            stream_power_active_nodes = self._K_unit_time * dt * Q_direct[active_nodes]**self._m * self.slopes[active_nodes]**self._n

        #Note that we save "stream_power_erosion" incorporating both K and a. Most definitions would need this value /K then **(1/a) to give actual stream power (unit, total, whatever), and it does not yet include the threshold
        self.stream_power_erosion[active_nodes] = stream_power_active_nodes
        grid.at_node['stream_power_erosion'] = self.stream_power_erosion
        #print "max stream power: ", self.stream_power_erosion.max()
        erosion_increment = (self.stream_power_erosion - self.sp_crit).clip(0.)

        #this prevents any node from incising below any node downstream of it
        #we have to go in upstream order in case our rate is so big we impinge on baselevels > 1 node away

        elev_dstr = node_z[flow_receiver]# we substract erosion_increment[flow_receiver] in the loop, as it can update

        method = 'cython'
        if method == 'cython':
            from .cfuncs import erode_avoiding_pits

            erode_avoiding_pits(node_order_upstream, flow_receiver, node_z,
                                erosion_increment)
        else:
            for i in node_order_upstream:
                elev_this_node_before = node_z[i]
                elev_this_node_after = elev_this_node_before - erosion_increment[i]
                elev_dstr_node_after = elev_dstr[i] - erosion_increment[flow_receiver[i]]
                if elev_this_node_after<elev_dstr_node_after:
                    erosion_increment[i] = (elev_this_node_before - elev_dstr_node_after)*0.999999 #we add a tiny elevation excess to prevent the module from ever totally severing its own flow paths

        node_z -= erosion_increment


        self.grid = grid

        return grid, node_z, self.stream_power_erosion
