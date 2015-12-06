from __future__ import print_function

import numpy as np
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY
from time import sleep
from landlab.utils import structured_grid as sgrid
import pylab

from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError


class TransportLimitedEroder(object):
    """
    This component implements transport limited erosion for a landscape in which
    flow directions are fully convergent. i.e., all nodes in the landscape have
    a single, uniquely defined downstream node.

    The module can in principle take multiple transport laws, but at the moment
    only Meyer-Peter Muller (MPM) is implemented.

    There is as yet no explicit stabilization check on the timestep. If your
    run destabilizes, try reducing dt.
    See, e.g., ./examples/simple_sp_driver.py

    Assumes grid does not deform or change during run.

    Note it is vital to ensure all your units match. t is assumed to be in
    years. Any length value is assumed to be in meters (including both dx
    and the uplift rate...!)

    DEJH Sept 14.
    Currently only runs on a raster grid
    Not well tested
    """

    def __init__(self, grid, params):
        self.initialize(grid, params)

# This draws attention to a potential problem. It will be easy to have modules update z, but because noone "owns" the data, to forget to also update dz/dx...
# How about a built in grid utility that updates "derived" data (i.e., using only grid fns, e.g., slope, curvature) at the end of any given tstep loop?
# Or an explicit flagging system for all variables in the modelfield indicating if they have been updated this timestep. (Currently implemented)
# Or wipe the existance of any derived grid data at the end of a timestep
# entirely, so modules find they don't have it next timestep.

    def initialize(self, grid, params_file):
        '''
        params_file is the name of the text file containing the parameters
        needed for this stream power component.

        ***Parameters for input file***
        OBLIGATORY:
            * Qc -> String. Controls how to set the carrying capacity.
                Either 'MPM', or a string giving the name of the model field
                where capacity values are stored on nodes.
                At the moment, only 'MPM' is permitted as a way to set the
                capacity automatically, but expansion would be trivial.
                If 'from_array', the module will attempt to set the capacity
                Note capacities must be specified as volume flux.
            *

            ...Then, assuming you set Qc=='MPM':
            * b_sp, c_sp -> Floats. These are the powers on discharge and
                drainage area in the equations used to control channel width and
                basin hydrology, respectively:
                        W = k_w * Q**b_sp
                        Q = k_Q * A**c_sp
                These parameters are used to constrain flow depth, and may be
                omitted if use_W or use_Q are set.
            *k_Q, k_w, mannings_n -> floats. These are the prefactors on the
                basin hydrology and channel width-discharge relations, and n
                from the Manning's equation, respectively. These are
                needed to allow calculation of shear stresses and hence carrying
                capacities from the local slope and drainage area alone.
                The equation for depth used to derive shear stress and hence
                carrying capacity contains a prefactor:
                    mannings_n*(k_Q**(1-b)/K_w)**0.6
                (so shear = fluid_density*g*depth_equation_prefactor*A**(0.6*c*(1-b)*S**0.7 !)
                Don't know what to set these values to? k_w=0.002, k_Q=1.e-9,
                mannings_n=0.03 give vaguely plausible numbers (e.g., for a
                drainage area ~400km2, like Boulder Creek at Boulder,
                => depth~2.5m, width~35m, shear stress ~O(1000Pa)).
            *Dchar -> float.  The characteristic grain diameter in meters
                (==D50 in most cases) used to calculate Shields numbers
                in the channel. If you want to define Dchar values at each node,
                don't set, and use the Dchar_if_used argument in erode()
                instead.

        OPTIONS:
            *rock_density -> in kg/m3 (defaults to 2700)
            *sediment_density -> in kg/m3 (defaults to 2700)
            *fluid_density -> in most cases water density, in kg/m3 (defaults to 1000)
            *g -> acceleration due to gravity, in m/s**2 (defaults to 9.81)

            *threshold_shields -> +ve float; the threshold taustar_crit.
                Defaults to 0.047, or if 'slope_sensitive_threshold' is set True,
                becomes a weak function of local slope following Lamb et al
                (2008):
                    threshold_shields=0.15*S**0.25
            *slope_sensitive_threshold -> bool, defaults to 'False'.
                If true, threshold_shields is set according to the Lamb
                equation. An exception will be raised if threshold_shields is
                also set.
            *Parker_epsilon -> float, defaults to 0.4. This is Parker's (1978)
                epsilon, which is used in the relation
                    tau - tauc = tau * (epsilon/(epsilon+1))
                The 0.4 default is appropriate for coarse (gravelly) channels.
                The value approaches infinity as the river banks become more
                cohesive.
            *dt -> +ve float. If set, this is the fixed timestep for this
                component. Can be overridden easily as a parameter in erode().
                If not set (default), this parameter MUST be set in erode().
            *use_W -> Bool; if True, component will look for node-centered data
                describing channel width in grid.at_node['channel_width'], and
                use it to implement incision ~ stream power per unit width.
                Defaults to False.
            *use_Q -> Bool. Overrides the basin hydrology relation, using an
                local water discharge value assumed already calculated and
                stored in grid.at_node['discharge'].
            *C_MPM -> float. Defaults to 1. Allows tuning of the MPM prefactor,
                which is calculated as
                    Qc = 8.*C_MPM*(taustar - taustarcrit)**1.5
                In almost all cases, tuning depth_equation_prefactor' is
                preferred to tuning this parameter.
            *return_capacity -> bool (default False). NOT YET IMPLEMENTED.
                If True, this component
                will save the calculated capacity in the field
                'fluvial_sediment_transport_capacity'. (Requires some additional
                math, so is suppressed for speed by default).

        '''
        self.grid = grid
        # needs to be filled with values in execution
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links + 1)
        self.count_active_links = np.zeros_like(
            self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1
        inputs = ModelParameterDictionary(params_file)
        try:
            self.g = inputs.read_float('g')
        except MissingKeyError:
            self.g = 9.81
        try:
            self.rock_density = inputs.read_float('rock_density')
        except MissingKeyError:
            self.rock_density = 2700.
        try:
            self.sed_density = inputs.read_float('sediment_density')
        except MissingKeyError:
            self.sed_density = 2700.
        try:
            self.fluid_density = inputs.read_float('fluid_density')
        except MissingKeyError:
            self.fluid_density = 1000.
        self.rho_g = self.fluid_density * self.g

        try:
            self.Qc = inputs.read_string('Qc')
        except MissingKeyError:
            raise MissingKeyError("Qc must be 'MPM' or a grid field name!")
        else:
            if self.Qc == 'MPM':
                self.calc_cap_flag = True
            else:
                self.calc_cap_flag = False

        try:
            self.lamb_flag = inputs.read_bool('slope_sensitive_threshold')
        except:
            self.lamb_flag = False
        try:
            self.shields_crit = inputs.read_float('threshold_shields')
            # flag for sed_flux_dep_incision to see if the threshold was
            # manually set.
            self.set_threshold = True
            print("Found a threshold to use: ", self.shields_crit)
            assert self.lamb_flag == False
        except MissingKeyError:
            if not self.lamb_flag:
                self.shields_crit = 0.047
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
            self.return_capacity = inputs.read_bool('return_capacity')
        except MissingKeyError:
            self.return_capacity = False

        try:
            self._b = inputs.read_float('b_sp')
        except MissingKeyError:
            if self.use_W:
                self._b = 0.
            else:
                if self.calc_cap_flag:
                    raise NameError('b was not set')
        try:
            self._c = inputs.read_float('c_sp')
        except MissingKeyError:
            if self.use_Q:
                self._c = 1.
            else:
                if self.calc_cap_flag:
                    raise NameError('c was not set')
        try:
            self.Dchar_in = inputs.read_float('Dchar')
        except MissingKeyError:
            pass

        # assume Manning's equation to set the power on A for shear stress:
        self.shear_area_power = 0.6 * self._c * (1. - self._b)

        self.k_Q = inputs.read_float('k_Q')
        self.k_w = inputs.read_float('k_w')
        mannings_n = inputs.read_float('mannings_n')
        if mannings_n < 0. or mannings_n > 0.2:
            print("***STOP. LOOK. THINK. You appear to have set Manning's n outside it's typical range. Did you mean it? Proceeding...***")
            sleep(2)
        self.depth_prefactor = self.rho_g * mannings_n * \
            (self.k_Q**(1. - self._b) / self.k_w)**0.6
        # Note the depth_prefactor we store already holds rho*g
        try:
            epsilon = inputs.read_float('Parker_epsilon')
        except MissingKeyError:
            epsilon = 0.4

        try:
            self.C_MPM = inputs.read_float('C_MPM')
        except MissingKeyError:
            self.C_MPM = 1.
        try:
            self.shields_prefactor = 1. / \
                ((self.sed_density - self.fluid_density) * self.g * self.Dchar_in)
            self.MPM_prefactor = 8. * self.C_MPM * \
                np.sqrt(self.relative_weight * self.Dchar_in *
                        self.Dchar_in * self.Dchar_in)
            self.MPM_prefactor_alt = 4. * \
                self.g**(-2. / 3.) / self.excess_SG / \
                self.fluid_density / self.sed_density
        except AttributeError:
            # have to set these manually as needed
            self.shields_prefactor_noD = 1. / \
                ((self.sed_density - self.fluid_density) * self.g)
        self.diffusivity_prefactor = 8. * np.sqrt(8. * self.g) / (self.sed_density / self.fluid_density - 1.) * (epsilon / (
            epsilon + 1.))**1.5 * mannings_n**(5. / 6.) * self.k_w**-0.9 * self.k_Q**(0.9 * (1. - self._b))  # ...this is multiplied by A**c(1-0.1*(1-b))
        # we consciously skip out a factor of S**0.05-->1. in the diffusion prefactor, to avoid delinearizing the diffusion. Only a possible problem at tiny S (20% error @S==0.01; 37% error @S==10**-4)
        # we could include this as a static adjustment in the actual looping code (i.e., just multiply by S**0.05, and don't work with it as part of the problem)
        # in reality, Manning's n changes downstream too, so... whatever
        self.diffusivity_power_on_A = 0.9 * self._c * \
            (1. - self._b)  # i.e., q/D**(1/6)

        self.cell_areas = np.empty(grid.number_of_nodes)
        self.cell_areas.fill(np.mean(grid.area_of_cell))
        self.cell_areas[grid.node_at_cell] = grid.area_of_cell
        self.dx2 = grid.dx ** 2
        self.dy2 = grid.dy ** 2
        self.bad_neighbor_mask = np.equal(
            grid.get_active_neighbors_at_node(bad_index=-1), -1)

    def erode(self, grid, dt, node_drainage_areas='drainage_area',
              node_elevs='topographic__elevation',
              W_if_used=None, Q_if_used=None,
              Dchar_if_used=None, io=None):
        """
        This method calculates the fluvial sediment transport capacity at all
        nodes, then erodes or deposits sediment as required.

        *grid* & *dt* are the grid object and timestep (float) respectively.

        *Node_drainage_areas* tells the component where to look for the drainage
        area values. Change to another string to override which grid field the
        component looks at, or pass a nnodes-long array of drainage areas values
        directly instead.

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
        the component to look for the data in grid.at_node['topographic__elevation'];
        set to 'name_of_field' to override this name, or pass an nnode-array
        to use those values as elevations instead. Using this option is
        considerably slower than any of the alternatives, as it also has to
        calculate the link_node_mapping from stratch each time.

        In both these cases, at present the mapping is to use the maximum
        slope of --any-- link attached to the node as the representative node
        slope. This is primarily for speed, but may be a good idea to modify
        later.

        *W_if_used* and *Q_if_used* must be provided if you set use_W and use_Q
        respectively in the component initialization. They can be either field
        names or nnodes arrays as in the other cases.

        *Dchar_if_used* must be set as a grid field string or nnoodes-long array
        if 'Dchar' as a float was not provided in the input file. (If it was,
        this will be overridden).

        RETURNS XXXXXX
        """
        dx = grid.dx
        dy = grid.dy
        dx2 = self.dx2
        dy2 = self.dy2
        nrows = grid.number_of_node_rows
        ncols = grid.number_of_node_columns

        try:
            self.Dchar = self.Dchar_in
        except AttributeError:
            try:
                self.Dchar = grid.at_node[Dchar_if_used]
            except FieldError:
                assert type(Dchar_if_used) == np.ndarray
                self.Dchar = Dchar_if_used

        if type(node_elevs) == str:
            node_z = grid.at_node[node_elevs]
        else:
            node_z = node_elevs
        node_z_asgrid = node_z.view().reshape((nrows, ncols))

        if type(node_drainage_areas) == str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        # this still "feels" the influence of closed nodes. Need to
        # get rid of them.

        all_nodes_diffusivity = self.diffusivity_prefactor * \
            node_A**self.diffusivity_power_on_A
        # ALT
        neighbor_nodes = grid.get_active_neighbors_at_node(bad_index=-1)
        # the -1 lets us get *some* value for all nodes, which we then mask:
        neighbor_diffusivities = np.ma.array(
            all_nodes_diffusivity[neighbor_nodes], mask=self.bad_neighbor_mask)
        # pylab.figure(1)
        # pylab.imshow(neighbor_diffusivities[:,3].reshape((nrows,ncols)))
        # pylab.colorbar()
        #order is E,N,W,S
        mean_diffusivities_byspacing = neighbor_diffusivities + \
            all_nodes_diffusivity.reshape(
                (grid.number_of_nodes, 1))  # not yet _byspacing...
        mean_diffusivities_byspacing[:, [0, 2]] /= (2. * dx)
        mean_diffusivities_byspacing[:, [1, 3]] /= (2. * dy)  # now complete
        # for in the loop
        rate_of_z_change_store = np.ma.empty(
            (mean_diffusivities_byspacing.shape[0], 2), dtype=float)
        # set up the Von Neumann stability criterion:
        np.sum(mean_diffusivities_byspacing[
               :, [0, 2]] / dx, axis=1, out=rate_of_z_change_store[:, 0])
        np.sum(mean_diffusivities_byspacing[
               :, [1, 3]] / dy, axis=1, out=rate_of_z_change_store[:, 1])
        max_sum_of_Ds = np.amax(np.sum(rate_of_z_change_store, axis=1))

        # adjust the tstep for Von Neumann stability here:
        delta_t_internal = 1. / max_sum_of_Ds
        num_reps_internal = int(dt // delta_t_internal)
        # print 'num reps: ', num_reps_internal
        dt_excess = dt % delta_t_internal

        for reps in xrange(num_reps_internal + 1):
            # check if last loop
            if reps == num_reps_internal:
                delta_t_internal = dt_excess

            node_gradients = grid.calculate_gradient_along_node_links(node_z)
            # pylab.figure(2)
            # pylab.imshow(node_gradients[:,3].reshape((nrows,ncols)))
            # pylab.colorbar()
            # pylab.show()
            # this method returns (nnodes,4), & masked values where links are inactive
            # ordering is E,N,W,S, same as nieghbors and hence Ds
            # the masks shoud be propagating forward still
            D_slope_product = node_gradients * mean_diffusivities_byspacing
            np.sum(D_slope_product[:, :2], axis=1,
                   out=rate_of_z_change_store[:, 0])
            np.sum(-D_slope_product[:, 2:], axis=1,
                   out=rate_of_z_change_store[:, 1])
            # this use of sum is necessary to preserve the right masking
            rate_of_z_change = np.sum(rate_of_z_change_store, axis=1)
            node_z[grid.core_nodes] += delta_t_internal * \
                rate_of_z_change[grid.core_nodes]

        ##########
        # all_nodes_diffusivity_asgrid = all_nodes_diffusivity.view().reshape((nrows,ncols)) #view prevents a copy from being made accidentally
        #interior_nodes_diffusivity = all_nodes_diffusivity_asgrid[1:-1,1:-1]
        # the "subsets" are only defined on the interior nodes (non-perimeter)
        #D_Ebydx2 = (all_nodes_diffusivity_asgrid[2:,1:-1]+interior_nodes_diffusivity)/2./dx2
        #D_Wbydx2 = (all_nodes_diffusivity_asgrid[:-2,1:-1]+interior_nodes_diffusivity)/2./dx2
        # D_Nbydy2 = (all_nodes_diffusivity_asgrid[1:-1,:-2]+interior_nodes_diffusivity)/2./dy2 #remember, the grid is upside down by this reshape
        #D_Sbydy2 = (all_nodes_diffusivity_asgrid[1:-1,2:]+interior_nodes_diffusivity)/2./dy2
        #sum_of_Ds = D_Ebydx2 + D_Wbydx2 + D_Nbydy2 + D_Sbydy2
        #
        # adjust the tstep for Von Neumann stability here:
        #delta_t_internal = 1./np.amax(sum_of_Ds)
        #num_reps_internal = int(dt//delta_t_internal)
        # print 'num reps: ', num_reps_internal
        #dt_excess = dt%delta_t_internal
        # we won't be adjusting the flow paths during step subdivision, so this fortunately stays linear, and we only do it once...
        #one_less_delta_t_times_sumofDs = 1. - delta_t_internal*sum_of_Ds
        #
        # for reps in xrange(num_reps_internal+1):
        #    #check if last loop
        #    if reps==num_reps_internal:
        #        delta_t_internal = dt_excess
        #
        #    z_ij = node_z_asgrid[1:-1,1:-1]
        #    z_iplus1j = node_z_asgrid[2:,1:-1]
        #    z_iminus1j = node_z_asgrid[:-2,1:-1]
        #    z_ijplus1 = node_z_asgrid[1:-1,:-2]
        #    z_ijminus1 = node_z_asgrid[1:-1,2:]
        #
        #    #new_z_internal = one_less_delta_t_times_sumofDs * z_ij + delta_t_internal*(D_Ebydx2*z_iplus1j
        #    #                                                                          +D_Wbydx2*z_iminus1j
        #    #                                                                          +D_Nbydy2*z_ijplus1
        #    #                                                                          +D_Sbydy2*z_ijminus1)
        #    new_z_internal = z_ij + delta_t_internal*(D_Ebydx2*(z_iplus1j-z_ij)-D_Wbydx2*(z_ij-z_iminus1j)+D_Nbydy2*(z_ijplus1-z_ij)-D_Sbydy2*(z_ij-z_ijminus1))
        #    print 'mean_diff', np.mean(np.fabs(node_z_asgrid[1:-1,1:-1]-new_z_internal))
        #    print 'max_diff', np.max(np.fabs(node_z_asgrid[1:-1,1:-1]-new_z_internal))
        #    #we don't add uplift as we do it manually outside the component
        #    node_z_asgrid[1:-1,1:-1] = new_z_internal
        #    #repin the BCs, in case there are any interior boundaries
        #    #grid.at_node[node_elevs][sgrid.interior_nodes((nrows,ncols))] = node_z_asgrid.ravel()
        #    #grid.update_boundary_nodes()
        #    #...remember, because everything has hopefully been passed by reference, the
        #    #fields should already have updated to reflect changes to z here
        #
        # to see if this is actually necessary
        ##grid.at_node[node_elevs][sgrid.interior_nodes((nrows,ncols))] = node_z_asgrid.ravel()
        self.grid = grid

        active_nodes = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]
        if io:
            try:
                io[active_nodes] += node_z_asgrid.ravel()[active_nodes]
            except TypeError:
                if type(io) == str:
                    elev_name = io
            else:
                return grid, io

        else:
            elev_name = 'topographic__elevation'

        return grid, grid.at_node[elev_name], all_nodes_diffusivity
