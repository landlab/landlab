from __future__ import print_function

import numpy as np
from time import sleep
from landlab import ModelParameterDictionary, CLOSED_BOUNDARY

from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError
from landlab.grid.base import BAD_INDEX_VALUE


class SedDepEroder(object):
    """
    This class implements sediment flux dependent fluvial incision. It is built
    on the back of the simpler stream power class, stream_power.py, also in this
    component, and follows its limitations - we require single flow directions,
    provided to the class. See the docstrings of stream_power.py for more detail
    on required initialization and operating parameters.

    Assumes grid does not deform or change during run.
    """

    def __init__(self, grid, params):
        self.initialize(grid, params)

    def initialize(self, grid, params_file):
        """
        This module implements sediment flux dependent channel incision following:

        E = f(Qs, Qc) * stream_power - sp_crit,

        where stream_power is the stream power (often ==K*A**m*S**n) provided by the
        stream_power.py component. Note that under this incision paradigm, sp_crit
        is assumed to be controlled exclusively by sediment mobility, i.e., it is
        not a function of bedrock resistance. If you want it to represent a bedrock
        resistance term, be sure to set Dchar if you use the MPM transport capacity
        relation, and do not use the flag 'slope_sensitive_threshold'.

        The component currently assumes that the threshold on bed incision is
        controlled by the threshold of motion of its sediment cover. This means
        there is assumed interplay between the supplied Shields number,
        characteristic grain size, and shear stress threshold.

        This calculation has a tendency to be slow, and can easily result in
        numerical instabilities. These instabilities are suppressed by retaining a
        memory of what the sediment flux was in the last time step, and weighting
        the next timestep by that value. XXXmore detail needed. Possibilities:
            1. weight by last timestep/2timesteps (what about early ones?)
            2. do it iteratively; only do incision one the sed flux you are using
            stabilises (so the previous iter "seed" becomes less important)

        Parameters needed in the initialization file follow those for
        stream_power.py. However, we now require additional input terms for the
        f(Qs,Qc) term:
        REQUIRED:
            *Set the stream power terms using a, b, and c NOT m, n.
            *...remember, any threshold set is set as tau**a, not just tau.
            *k_Q, k_w, mannings_n -> floats. These are the prefactors on the
                basin hydrology and channel width-discharge relations, and n
                from the Manning's equation, respectively. These are
                needed to allow calculation of shear stresses and hence carrying
                capacities from the local slope and drainage area alone.
                Don't know what to set these values to? k_w=2.5, k_Q=2.5e-7,
                mannings_n=0.05 give vaguely plausible numbers with b=0.5,
                c = 1.(e.g., for a drainage area ~350km2, like Boulder Creek
                at Boulder, => depth~1.3m, width~23m, shear stress ~O(200Pa)
                for an "annual-ish" flood). [If you want to continue playing
                with calibration, the ?50yr return time 2013 floods produced
                depths ~2.3m with Q~200m3/s]
            *sed_dependency_type -> 'None', 'linear_decline', 'parabolic',
                'almost_parabolic', 'generalized_humped'. For definitions, see Gasparini
                et al., 2006; Hobley et al., 2011.
            *Qc -> This input controls the sediment capacity used by the component.
                It can calculate sediment carrying capacity for itself if this
                parameter is a string 'MPM', which will cause the component to use a
                slightly modified version of the Meyer-Peter Muller equation (again, see
                Hobley et al., 2011). Alternatively, it can be another string denoting
                the grid field name in which a precalculated capacity is stored.

        Depending on which options are specified above, further parameters may be
        required:
            *If sed_dependency_type=='generalized_humped', need the shape parameters
            used by Hobley et al:
                kappa_hump
                nu_hump
                phi_hump
                c_hump
                Note the onus is on the user to ensure that these parameters result
                in a viable shape, i.e., one where the maximum is 1 and there is
                indeed a hump in the profile. If these parameters are NOT specified,
                they will default to the form of the curve for Leh valley as found
                in Hobley et al 2011: nu=1.13; phi=4.24; c=0.00181; kappa=13.683.

            *If Qc=='MPM', these parameters may optionally be provided:
                Dchar -> characteristic grain size (i.e., D50) on the bed, in m.
                C_MPM -> the prefactor in the MPM relation. Defaults to 1, as in the
                    relation sensu stricto, but can be modified to "tune" the
                    equations to a known point where sediment deposition begins.
                    In cases where k_Q and k_w are not known from real data, it
                    is recommended these parameters be tuned in preference to C.

            *...if Dchar is NOT provided, the component will attempt to set (and will
                report) an appropriate characteristic grain size, such that it is
                consistent both with the threshold provided *and* a critical Shields
                number of 0.05. (If you really, really want to, you can override this
                critical Shields number too; use parameter *threshold_Shields*).

        OPTIONAL:
            *rock_density -> in kg/m3 (defaults to 2700)
            *sediment_density -> in kg/m3 (defaults to 2700)
            *fluid_density -> in most cases water density, in kg/m3 (defaults to 1000)
            *g -> acceleration due to gravity, in m/s**2 (defaults to 9.81)

            *threshold_shear_stress -> a float.
                If not provided, may be overridden by the following parameters.
                If it is not, defaults to 0.
            *slope_sensitive_threshold -> a boolean, defaults to FALSE.
                In steep mountain environments, the critical Shields number for particle
                motion appears to be weakly sensitive to the local slope, as
                taustar_c=0.15*S**0.25 (Lamb et al, 2008). If this flag is set to TRUE,
                the critical threshold in the landscape is allowed to become slope
                sensitive as well, in order to be consistent with this equation.
                This modification was used by Hobley et al., 2011.

            *set_threshold_from_Dchar -> a boolean, defaults to FALSE.
                Use this flag to force an appropriate threshold value from a
                provided Dchar. i.e., this is the inverse of the procedure that is
                used to find Dchar if it isn't provided. No threshold can be
                specified in the parameter file, and Dchar must be specified.

            *return_stream_properties -> bool (default False).
                If True, this component will save the calculations for
                'channel_width', 'channel_depth', and 'channel_discharge' in
                those grid fields. (Requires some
                additional math, so is suppressed for speed by default).
        """
        #this is the fraction we allow any given slope in the grid to evolve by in one go (suppresses numerical instabilities)
        self.fraction_gradient_change = 0.25
        self.pseudoimplicit_repeats = 5
        self.grid = grid
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links+1) #needs to be filled with values in execution
        self.count_active_links = np.zeros_like(self.link_S_with_trailing_blank, dtype=int)
        self.count_active_links[:-1] = 1
        inputs = ModelParameterDictionary(params_file)
        try:
            self.thresh = inputs.read_float('threshold_shear_stress')
            self.set_threshold = True #flag for sed_flux_dep_incision to see if the threshold was manually set.
            print("Found a shear stress threshold to use: ", self.thresh)
        except MissingKeyError:
            print("Found no incision threshold to use.")
            self.thresh = 0.
            self.set_threshold = False
        try:
            self._a = inputs.read_float('a_sp')
        except:
            print("a not supplied. Setting power on shear stress to 1...")
            self._a = 1.
        try:
            self._b = inputs.read_float('b_sp')
        except MissingKeyError:
            #if self.use_W:
            #    self._b = 0.
            #else:
            raise NameError('b was not set')
        try:
            self._c = inputs.read_float('c_sp')
        except MissingKeyError:
            #if self.use_Q:
            #    self._c = 1.
            #else:
            raise NameError('c was not set') #we need to restore this functionality later
        #'To use the sed flux dependent model, you must set a,b,c not m,n. Try a=1,b=0.5,c=1...?'

        self._K_unit_time = inputs.read_float('K_sp')/31557600. #...because we work with dt in seconds
        #set gravity
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
        self.relative_weight = (self.sed_density-self.fluid_density)/self.fluid_density*self.g #to accelerate MPM calcs
        self.rho_g = self.fluid_density*self.g

        self.k_Q = inputs.read_float('k_Q')
        self.k_w = inputs.read_float('k_w')
        mannings_n = inputs.read_float('mannings_n')
        self.mannings_n = mannings_n
        if mannings_n<0. or mannings_n>0.2:
            print("***STOP. LOOK. THINK. You appear to have set Manning's n outside its typical range. Did you mean it? Proceeding...***")
            sleep(2)

        self.diffusivity_power_on_A = 0.9*self._c*(1.-self._b) #i.e., q/D**(1/6)

        self.type = inputs.read_string('sed_dependency_type')
        try:
            self.Qc = inputs.read_string('Qc')
        except MissingKeyError:
            self.Qc = None
        else:
            if self.Qc=='MPM':
                self.calc_cap_flag = True
            else:
                self.calc_cap_flag = False
        try:
            self.override_threshold = inputs.read_bool('set_threshold_from_Dchar')
        except MissingKeyError:
            self.override_threshold = False
        try:
            self.shields_crit = inputs.read_float('threshold_Shields')
        except MissingKeyError:
            self.shields_crit = 0.05
        try:
            self.return_ch_props = inputs.read_bool('return_stream_properties')
        except MissingKeyError:
            self.return_ch_props = False

        #now conditional inputs
        if self.type == 'generalized_humped':
            try:
                self.kappa = inputs.read_float('kappa_hump')
                self.nu = inputs.read_float('nu_hump')
                self.phi = inputs.read_float('phi_hump')
                self.c = inputs.read_float('c_hump')
            except MissingKeyError:
                self.kappa = 13.683
                self.nu = 1.13
                self.phi = 4.24
                self.c = 0.00181
                print('Adopting inbuilt parameters for the humped function form...')

        try:
            self.lamb_flag = inputs.read_bool('slope_sensitive_threshold')
            #this is going to be a nightmare to implement...
        except:
            self.lamb_flag = False

        if self.Qc == 'MPM':
            try:
                self.Dchar_in = inputs.read_float('Dchar')
            except MissingKeyError:
                assert self.thresh > 0., "Can't automatically set characteristic grain size if threshold is 0 or unset!"
                #remember the threshold getting set is already tau**a
                self.Dchar_in = self.thresh/self.g/(self.sed_density-self.fluid_density)/self.shields_crit
                print('Setting characteristic grain size from the Shields criterion...')
                print('Characteristic grain size is: ', self.Dchar_in)
            try:
                self.C_MPM = inputs.read_float('C_MPM')
            except MissingKeyError:
                self.C_MPM = 1.

        if self.override_threshold:
            print("Overriding any supplied threshold...")
            try:
                self.thresh = self.shields_crit*self.g*(self.sed_density-self.fluid_density)*self.Dchar_in
            except AttributeError:
                self.thresh = self.shields_crit*self.g*(self.sed_density-self.fluid_density)*inputs.read_float('Dchar')
            print("Threshold derived from grain size and Shields number is: ", self.thresh)

        self.cell_areas = np.empty(grid.number_of_nodes)
        self.cell_areas.fill(np.mean(grid.area_of_cell))
        self.cell_areas[grid.node_at_cell] = grid.area_of_cell
        #new 11/12/14
        self.point6onelessb = 0.6*(1.-self._b)
        self.shear_stress_prefactor = self.fluid_density*self.g*(self.mannings_n/self.k_w)**0.6

        if self.set_threshold is False or self.override_threshold:
            try:
                self.shields_prefactor_to_shear = (self.sed_density-self.fluid_density)*self.g*self.Dchar_in
            except AttributeError: #no Dchar
                self.shields_prefactor_to_shear_noDchar = (self.sed_density-self.fluid_density)*self.g

        twothirds = 2./3.
        self.Qs_prefactor = 4.*self.C_MPM**twothirds*self.fluid_density**twothirds/(self.sed_density-self.fluid_density)**twothirds*self.g**(twothirds/2.)*mannings_n**0.6*self.k_w**(1./15.)*self.k_Q**(0.6+self._b/15.)/self.sed_density**twothirds
        self.Qs_thresh_prefactor = 4.*(self.C_MPM*self.k_w*self.k_Q**self._b/self.fluid_density**0.5/(self.sed_density-self.fluid_density)/self.g/self.sed_density)**twothirds
        #both these are divided by sed density to give a vol flux
        self.Qs_power_onA = self._c*(0.6+self._b/15.)
        self.Qs_power_onAthresh = twothirds*self._b*self._c




    def get_sed_flux_function(self, rel_sed_flux):
        if self.type == 'generalized_humped':
            "Returns K*f(qs,qc)"
            sed_flux_fn = self.kappa*(rel_sed_flux**self.nu + self.c)*np.exp(-self.phi*rel_sed_flux)
        elif self.type == 'linear_decline':
            sed_flux_fn = (1.-rel_sed_flux)
        elif self.type == 'parabolic':
            raise MissingKeyError('Pure parabolic (where intersect at zero flux is exactly zero) is currently not supported, sorry. Try almost_parabolic instead?')
            sed_flux_fn = 1. - 4.*(rel_sed_flux-0.5)**2.
        elif self.type == 'almost_parabolic':
            sed_flux_fn = np.where(rel_sed_flux>0.1, 1. - 4.*(rel_sed_flux-0.5)**2., 2.6*rel_sed_flux+0.1)
        elif self.type == 'None':
            sed_flux_fn = 1.
        else:
            raise MissingKeyError('Provided sed flux sensitivity type in input file was not recognised!')
        return sed_flux_fn


    def get_sed_flux_function_pseudoimplicit(self, sed_in, trans_cap_vol_out, prefactor_for_volume, prefactor_for_dz):
        rel_sed_flux_in = sed_in/trans_cap_vol_out
        rel_sed_flux = rel_sed_flux_in
        if self.type == 'generalized_humped':
            "Returns K*f(qs,qc)"
            def sed_flux_fn_gen(rel_sed_flux_in):
                return self.kappa*(rel_sed_flux_in**self.nu + self.c)*np.exp(-self.phi*rel_sed_flux_in)
        elif self.type == 'linear_decline':
            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1.-rel_sed_flux_in
        elif self.type == 'parabolic':
            raise MissingKeyError('Pure parabolic (where intersect at zero flux is exactly zero) is currently not supported, sorry. Try almost_parabolic instead?')
            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1. - 4.*(rel_sed_flux_in-0.5)**2.
        elif self.type == 'almost_parabolic':
            def sed_flux_fn_gen(rel_sed_flux_in):
                return np.where(rel_sed_flux_in>0.1, 1. - 4.*(rel_sed_flux_in-0.5)**2., 2.6*rel_sed_flux_in+0.1)
        elif self.type == 'None':
            def sed_flux_fn_gen(rel_sed_flux_in):
                return 1.
        else:
            raise MissingKeyError('Provided sed flux sensitivity type in input file was not recognised!')

        for i in xrange(self.pseudoimplicit_repeats):
            sed_flux_fn = sed_flux_fn_gen(rel_sed_flux)
            sed_vol_added = prefactor_for_volume*sed_flux_fn
            rel_sed_flux = rel_sed_flux_in + sed_vol_added/trans_cap_vol_out
            #print rel_sed_flux
            if rel_sed_flux>=1.:
                rel_sed_flux = 1.
                break
            if rel_sed_flux<0.:
                rel_sed_flux = 0.
                break
        last_sed_flux_fn = sed_flux_fn
        sed_flux_fn = sed_flux_fn_gen(rel_sed_flux)
        #this error could alternatively be used to break the loop
        error_in_sed_flux_fn = sed_flux_fn-last_sed_flux_fn
        dz = prefactor_for_dz*sed_flux_fn
        sed_flux_out = rel_sed_flux*trans_cap_vol_out
        return dz, sed_flux_out, rel_sed_flux, error_in_sed_flux_fn


    def erode(self, grid, dt=None, node_elevs='topographic__elevation',
                node_drainage_areas='drainage_area',
                node_receiving_flow='flow_receiver',
                node_order_upstream='upstream_node_order',
                node_slope='topographic__steepest_slope',
                steepest_link='links_to_flow_receiver',
                runoff_rate_if_used=None,
                #W_if_used=None, Q_if_used=None,
                stability_condition='loose',
                Dchar_if_used=None, io=None):

        """
        Note this method must be passed both 'receiver' and 'upstream_order',
        either as strings for field access, or nnode- long arrays of the
        relevant IDs. These are most easily used as the
        outputs from route_flow_dn.
        Note that you must supply either slopes_at_nodes or link_slopes.
        """
        if runoff_rate_if_used != None:
            runoff_rate = runoff_rate_if_used
            assert type(runoff_rate) in (int, float, np.ndarray)
        else:
            runoff_rate = 1.

        if dt==None:
            dt = self.tstep
        try:
            self.Dchar=self.Dchar_in
        except AttributeError:
            try:
                self.Dchar=grid.at_node[Dchar_if_used]
            except FieldError:
                assert type(Dchar_if_used)==np.ndarray
                self.Dchar=Dchar_if_used
            if not self.set_threshold:
                assert self.override_threshold, "You need to confirm to the module you intend it to internally calculate a shear stress threshold, with set_threshold_from_Dchar in the input file."
                #we need to adjust the thresholds for the Shields number & gs dynamically:
                variable_thresh = self.shields_crit*self.g*(self.sed_density-self.fluid_density)*self.Dchar
        else:
            assert Dchar_if_used is None, "Trouble ahead... you can't provide Dchar both in the input file and as an array!"

        if type(node_elevs)==str:
            node_z = grid.at_node[node_elevs]
        else:
            node_z = node_elevs

        if type(node_drainage_areas)==str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        if type(node_receiving_flow)==str:
            flow_receiver = grid.at_node[node_receiving_flow]
        else:
            flow_receiver = node_receiving_flow

        #new V3:
        if type(node_order_upstream)==str:
            s_in = grid.at_node[node_order_upstream]
        else:
            s_in = node_order_upstream

        if type(node_slope)==str:
            node_S = grid.at_node[node_slope]
        else:
            node_S = node_slope

        if self.lamb_flag:
            variable_shields_crit = 0.15*node_S**0.25
            try:
                variable_thresh = variable_shields_crit*self.shields_prefactor_to_shear
            except AttributeError:
                variable_thresh = variable_shields_crit*self.shields_prefactor_to_shear_noDchar*self.Dchar


        if type(steepest_link)==str:
            link_length = np.empty(grid.number_of_nodes,dtype=float)
            link_length.fill(np.nan)
            draining_nodes = np.not_equal(grid.at_node[steepest_link], BAD_INDEX_VALUE)
            core_draining_nodes = np.intersect1d(np.where(draining_nodes)[0], grid.core_nodes, assume_unique=True)
            link_length[core_draining_nodes] = grid.link_length[grid.at_node[steepest_link][core_draining_nodes]]
            #link_length=grid.dx
        else:
            link_length = grid.link_length[steepest_link]

        node_Q = self.k_Q*runoff_rate*node_A**self._c
        shear_stress_prefactor_timesAparts = self.shear_stress_prefactor*node_Q**self.point6onelessb
        try:
            transport_capacities_thresh = self.thresh*self.Qs_thresh_prefactor*runoff_rate**(0.66667*self._b)*node_A**self.Qs_power_onAthresh
        except AttributeError:
            transport_capacities_thresh = variable_thresh*self.Qs_thresh_prefactor*runoff_rate**(0.66667*self._b)*node_A**self.Qs_power_onAthresh

        transport_capacity_prefactor_withA = self.Qs_prefactor*runoff_rate**(0.6+self._b/15.)*node_A**self.Qs_power_onA

        internal_t = 0.
        break_flag = False
        dt_secs = dt*31557600.
        counter = 0
        rel_sed_flux = np.empty_like(node_Q)
        #excess_vol_overhead = 0.

        while 1: #use the break flag, to improve computational efficiency for runs which are very stable
            #we assume the drainage structure is forbidden to change during the whole dt
            #print "loop..."
            #note slopes will be *negative* at pits
            #track how many loops we perform:
            counter += 1
            #print counter
            downward_slopes = node_S.clip(0.)
            #positive_slopes = np.greater(downward_slopes, 0.)
            slopes_tothe07 = downward_slopes**0.7
            transport_capacities_S = transport_capacity_prefactor_withA*slopes_tothe07
            trp_diff = (transport_capacities_S - transport_capacities_thresh).clip(0.)
            transport_capacities = np.sqrt(trp_diff*trp_diff*trp_diff)
            shear_stress = shear_stress_prefactor_timesAparts*slopes_tothe07
            shear_tothe_a = shear_stress**self._a

            dt_this_step = dt_secs-internal_t #timestep adjustment is made AFTER the dz calc
            node_vol_capacities = transport_capacities*dt_this_step

            sed_into_node = np.zeros(grid.number_of_nodes, dtype=float)
            dz = np.zeros(grid.number_of_nodes, dtype=float)
            len_s_in = s_in.size
            cell_areas = self.cell_areas

            for i in s_in[::-1]: #work downstream
                try:
                    cell_area = cell_areas[i]
                except TypeError: #it's a float, not an array
                    cell_area = cell_areas
                sed_flux_into_this_node = sed_into_node[i]
                node_capacity = transport_capacities[i] #we work in volume flux, not volume per se here
                node_vol_capacity = node_vol_capacities[i]
                if sed_flux_into_this_node < node_vol_capacity: #note incision is forbidden at capacity
#                    sed_flux_ratio = sed_flux_into_this_node/node_capacity
#                    fqsqc=self.get_sed_flux_function(sed_flux_ratio)
#                    try:
#                        thresh = variable_thresh
#                    except: #it doesn't exist
#                        thresh = self.thresh
#                    dz_here = self._K_unit_time*dt_this_step*fqsqc*(shear_tothe_a[i]-thresh).clip(0.) #let's define down as +ve
#                    vol_pass_attempted = dz_here*cell_area + sed_flux_into_this_node
#                    if vol_pass_attempted > node_vol_capacity:
#                        #it's vital we don't allow the final node to erode more than it can remove from the node!!
#                        #=>must modify dz_here, not just drop that sed
#                        excess_volume = vol_pass_attempted-node_vol_capacity
#                        dz_reduction = excess_volume/cell_area
#                        dz_here -= dz_reduction
#                        #add a "sneak" to stop the transition point developing into a sed deposition shock:
#                        if dz_here < 0.:
#                            #excess_vol_overhead += (dz_reduction-dz_here)*cell_area
#                            #node_vol_capacities[flow_receiver[i]] += excess_vol_overhead
#                            dz_here = 0. #this sed packet gets to be "excess overhead" -> it gets freely transported downstream.
#                            #...but problems will arise if the system terminates before a boundary cell! (Mass leak)
#                            #also, upstream node is still incising freely; it can't feel this.
#                            #slopes can still reverse, if the initial gradient was small
#                            #=>no node may incise more than the next node down in the previous tstep?
#                            #could use this to set dt_internal...
#                            #do the "incising" part as normal, then pause once depo occurs,...
#                            #then reduce dt_int so it can't dig further than the next node anywhere...
#                            #then RERUN, don't scale dz
#                        vol_pass = node_vol_capacity
#                    else:
#                        vol_pass = vol_pass_attempted
##implementing the pseudoimplicit method instead:
                    try:
                        thresh = variable_thresh
                    except: #it doesn't exist
                        thresh = self.thresh
                    dz_prefactor = self._K_unit_time*dt_this_step*(shear_tothe_a[i]-thresh).clip(0.)
                    vol_prefactor = dz_prefactor*cell_area
                    dz_here, sed_flux_out, rel_sed_flux_here, error_in_sed_flux = self.get_sed_flux_function_pseudoimplicit(sed_flux_into_this_node, node_vol_capacity, vol_prefactor, dz_prefactor)
                    #note now dz_here may never create more sed than the out can transport...
                    assert sed_flux_out <= node_vol_capacity, 'failed at node '+str(s_in.size-i)+' with rel sed flux '+str(sed_flux_out/node_capacity)
                    rel_sed_flux[i] = rel_sed_flux_here
                    vol_pass = sed_flux_out#*dt_this_step
                else:
                    rel_sed_flux[i] = 1.
                    vol_dropped = sed_flux_into_this_node - node_vol_capacity
                    dz_here = -vol_dropped/cell_area
                    vol_pass = node_vol_capacity

                dz[i] -= dz_here
                sed_into_node[flow_receiver[i]] += vol_pass

#            #perform the fractional-slope-change stability analysis
#            elev_diff = node_z - node_z[flow_receiver]
#            delta_dz = dz - dz[flow_receiver] #remember, here dz is DOWN (unlike the TL case)
#            ###delta_dz = np.fabs(delta_dz)
#            ##excess_fraction = delta_dz/elev_diff
#            ##most_flattened_nodes = np.argmax(np.fabs(excess_fraction[grid.core_nodes]))
#            node_flattening = self.fraction_gradient_change*elev_diff - delta_dz #note the condition is that gradient may not change by >X%, not must be >0
#            #note all these things are zero for a pit node
#            most_flattened_nodes = np.argmin(node_flattening[grid.core_nodes])
#            most_flattened_nodes = np.take(grid.core_nodes, most_flattened_nodes) #get it back to node number, not core_node number
#            ##most_flattened_val = np.take(excess_fraction, most_flattened_nodes)
#            ##abs_most_flattened_val = np.fabs(most_flattened_val)
#            most_flattened_val = np.take(node_flattening, most_flattened_nodes)
#            print 'most flattened val: ', most_flattened_val
#            if most_flattened_val>=0.:
#            ##if abs_most_flattened_val<self.fraction_gradient_change:
#                break_flag = True #all nodes are stable
#            else: # a fraction > the critical fraction
#            ###need to think about if we can assume dz and dt behave linearly in these cases...
#            #first impression is that it's OK. Still assume everything works linearly *within any given tstep*
#                most_extreme_elev_diff = np.take(elev_diff, most_flattened_nodes)
#                most_extreme_delta_dz = np.take(delta_dz, most_flattened_nodes)
#                print 'elev_diff: ', most_extreme_elev_diff
#                print 'delta dz: ', most_extreme_delta_dz
#                print '***'
#                dt_fraction = 0.5*self.fraction_gradient_change*most_extreme_elev_diff/most_extreme_delta_dz
#                ##dt_fraction = self.fraction_gradient_change/abs_most_flattened_val
###persistent failure here happens when elev diff manages to reverse itself, which should be forbidden
###this might be showing that a linear model *isn't* adequate
###how about doing all the stability *before* the analysis? Once we have slope, capacity & fldir for all nodes,
###worst case scenario is that all the potential sed gets dumped in this node
###(better case - deduct the transport capacity out). ...isn't this just the TL stability condition??
#
#                #print 'dt_fraction: ', dt_fraction
#                #correct those elevs
#                dz *= dt_fraction
#                dt_this_step *= dt_fraction
            break_flag = True

            node_z[grid.core_nodes] += dz[grid.core_nodes]

            if break_flag:
                break
            #do we need to reroute the flow/recalc the slopes here? -> NO, slope is such a minor component of Diff we'll be OK
            #BUT could be important not for the stability, but for the actual calc. So YES.
            node_S = np.zeros_like(node_S)
            #print link_length[core_draining_nodes]
            node_S[core_draining_nodes] = (node_z-node_z[flow_receiver])[core_draining_nodes]/link_length[core_draining_nodes]
            internal_t += dt_this_step #still in seconds, remember

        self.grid=grid

        active_nodes = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]
        if io:
            try:
                io[active_nodes] = node_z[active_nodes]
            except TypeError:
                if type(io)==str:
                    elev_name = io
            else:
                return grid, io

        else:
            elev_name = node_elevs

        if self.return_ch_props:
            #add the channel property field entries,
            #'channel_width', 'channel_depth', and 'channel_discharge'
            W = self.k_w*node_Q**self._b
            H = shear_stress/self.rho_g/node_S #...sneaky!
            grid.at_node['channel_width'] = W
            grid.at_node['channel_depth'] = H

        grid.at_node['channel_discharge'] = node_Q
        grid.at_node['channel_bed_shear_stress'] = shear_stress
        grid.at_node['fluvial_sediment_transport_capacity'] = transport_capacities
        grid.at_node['fluvial_sediment_flux_into_node'] = sed_into_node
        grid.at_node['relative_sediment_flux'] = rel_sed_flux
        #elevs set automatically to the name used in the function call.
        self.iterations_in_dt = counter

        return grid, grid.at_node[elev_name]
