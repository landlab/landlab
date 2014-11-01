import numpy as np
from time import sleep
from scipy import weave
from landlab import ModelParameterDictionary

from landlab.core.model_parameter_dictionary import MissingKeyError
from landlab.field.scalar_data_fields import FieldError
from scipy.weave.build_tools import CompileError

from landlab.components.stream_power.stream_power import StreamPowerEroder

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
                The equation for depth used to derive shear stress and hence
                carrying capacity contains a prefactor:
                    mannings_n*(k_Q**(1-b)/K_w)**0.6
                If tuning this equation to adjust sediment capacities, the 
                parameter 'depth_equation_prefactor', equal to this value, can be
                set instead. 
                (so shear = fluid_density*g*depth_equation_prefactor*A**(0.6*c*(1-b)*S**0.7 !)
                It's also recommended to set 'sp_type' to 'Shear_stress'.
            *sed_dependency_type -> 'None', 'linear_decline', 'parabolic', 
                'almost_parabolic', 'generalized_humped'. For definitions, see Gasparini
                et al., 2006; Hobley et al., 2011.
            *Qc -> This input controls the sediment capacity used by the component.
                It can either calculate sediment carrying capacity for itself if this
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
                number of 0.06. (If you really, really want to, you can override this
                critical Shields number too; use parameter 'critical_Shields').
        
        OPTIONAL:
            *rock_density -> in kg/m3 (defaults to 2700)
            *sediment_density -> in kg/m3 (defaults to 2700)
            *fluid_density -> in most cases water density, in kg/m3 (defaults to 1000)
            *g -> acceleration due to gravity, in m/s**2 (defaults to 9.81)
            
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
        """
        
        inputs = ModelParameterDictionary(params_file)
        #create a initialization of stream power, which will provide the stream power value to modify:
        self.simple_sp = StreamPowerEroder(grid, params_file)
        self.simple_sp.no_erode = True #suppress the ability of the module to do any erosion
        self.threshold = self.simple_sp.sp_crit
        try:
            self._a = self.simple_sp._a
            self._b = self.simple_sp._b
            self._c = self.simple_sp._c
            self.shear_area_power = 0.6*self._c*(1.-self._b)
        except:
            raise MissingKeyError('To use the sed flux dependent model, you must set a,b,c not m,n. Try a=1,b=0.5,c=1...?')
        
        self._K_unit_time = inputs.read_float('K_sp')
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
        
        try:
            k_Q = inputs.read_float('k_Q')
        except MissingKeyError:
            self.depth_prefactor = self.rho_g*inputs.read_float('depth_equation_prefactor')
        else:
            k_w = inputs.read_float('k_w')
            mannings_n = inputs.read_float('mannings_n')
            if mannings_n<0. or mannings_n>0.2:
                print "***STOP. LOOK. THINK. You appear to have set Manning's n outside it's typical range. Did you mean it? Proceeding...***"
                sleep(2)
            self.depth_prefactor = self.rho_g*mannings_n*(k_Q**(1.-self._b)/k_w)**0.6
            ##Note the depth_prefactor we store already holds rho*g
        
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
            override_threshold = inputs.read_bool('set_threshold_from_Dchar')
        except MissingKeyError:
            override_threshold = False
        try:
            self.shields_crit = inputs.read_float('critical_Shields')
        except MissingKeyError:
            self.shields_crit = 0.06
            
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
                print 'Adopting inbuilt parameters for the humped function form...'
                
                def get_sed_flux_function(self, rel_sed_flux):
                    """
                    Returns K*f(qs,qc)
                    """
                    sed_flux_fn = self.kappa*(rel_sed_flux**self.nu + self.c)*np.exp(-self.phi*rel_sed_flux)
                    return sed_flux_fn
        
        elif self.type == 'linear_decline':
            def get_sed_flux_function(self, rel_sed_flux):
                sed_flux_fn = (1.-rel_sed_flux)
                return sed_flux_fn
                
        elif self.type == 'parabolic':
            def get_sed_flux_function(self, rel_sed_flux):
                sed_flux_fn = 1. - 4.*(rel_sed_flux-0.5)**2.
                return sed_flux_fn
        
        elif self.type == 'almost_parabolic':
            def get_sed_flux_function(self, rel_sed_flux):
                sed_flux_fn = np.where(rel_sed_flux>0.1, 1. - 4.*(rel_sed_flux-0.5)**2., 2.6*rel_sed_flux+0.1)
                return sed_flux_fn
        
        elif self.type == 'None':
            def get_sed_flux_function(self, rel_sed_flux):
                sed_flux_fn = 1.
                return sed_flux_fn
        
        else:
            raise MissingKeyError('Provided sed flux sensitivity type in input file was not recognised!')
                
        try:
            self.S_sensitive_thresh = inputs.read_bool('slope_sensitive_threshold')
            #this is going to be a nightmare to implement...
        except:
            self.S_sensitive_thresh = False
            
        if self.Qc == 'MPM':
            try:
                self.Dchar = inputs.read_float('Dchar')
            except MissingKeyError:
                assert self.threshold > 0., "Can't automatically set characteristic grain size if threshold is 0 or unset!"
                #remember the threshold getting set is already tau**a
                self.Dchar = self.threshold**(1./self._a)/self.g/(self.sed_density-self.fluid_density)/self.shields_crit
                print 'Setting characteristic grain size from the Shields criterion...'
                print 'Characteristic grain size is: ', self.Dchar
            try:
                self.eight_C_MPM = 8.*inputs.read_float('C_MPM')
            except MissingKeyError:
                self.eight_C_MPM = 8.
            self.shields_prefactor = 1./((self.sed_sensity-self.fluid_density)*self.g*self.Dchar)
            self.MPM_prefactor = self.eight_C_MPM*np.sqrt(self.relative_weight*self.Dchar*self.Dchar*self.Dchar)
        
        if override_threshold:
            assert self.simple_sp.set_threshold==False, 'Threshold cannot be manually set if you wish it to be generated from Dchar!'
            try:
                self.threshold = (self.shields_crit*self.g*(self.sed_density-self.fluid_density)*self.Dchar)**self._a
            except AttributeError:
                self.threshold = (self.shields_crit*self.g*(self.sed_density-self.fluid_density)*inputs.read_float('Dchar'))**self._a
            self.simple_sp.sp_crit = self.threshold
        
        def sed_capacity_equation(self, grid, shields_stress, slopes_at_nodes=None):
            """
            slopes_at_nodes must be provided as nnodes array if slope_sensitive_threshold.
            This is a volume flux.
            At the moment, this only supports MPM for sed flux capacity calculation,
            but is set up to allow easy additions later if needed.
            """
            if self.S_sensitive_thresh:
                thresh = 0.15*slopes_at_nodes**0.25
            else:
                thresh = self.shields_crit
            if self.Qc == 'MPM':
                capacity = self.MPM_prefactor*(shields_stress-thresh)**1.5
            elif type(self.Qc)==str:
                try:
                    capacity = grid.at_node[self.Qc]
                except:
                    raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided input...')
            elif type(self.Qc)==np.ndarray:
                capacity=self.Qc
            else:
                raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided inputs...')
            return capacity
        
        self.cell_areas = np.empty(grid.number_of_nodes)
        self.cell_areas.fill(np.mean(grid.cell_areas))
        self.cell_areas[grid.cell_nodes] = grid.cell_areas
        
        ###No, we won't do this. ID the first step with try: self.past_sed_flux; except AttributeError: ...
        ##self.past_sed_flux = np.empty(grid.number_of_nodes) #this is where we're going to store the previous sed flux...
        ##self.first_iter = True #this will let us identify the "startup" condition, which we need to allow to become stable...
        
    def weave_iter_sed_flux(self, grid, incision_iter, r, s):
        '''
        Needs to be passed at node values for 'flow_receiver' and 
        'upstream_ID_order' (i.e., r and s). These are most likely the output
        values from route_flow_dn.
        fqsqc is used to force the result; it isn't dynamically adjusted here.
        '''
        code =  """
                int donor;
                int rcvr;
                for (int i = num_pts-1; i > -1; i--) {
                    donor = s[i];
                    rcvr = r[donor];
                    if (donor != rcvr) {
                        sed_flux[rcvr] += sed_flux[donor];
                    }
                }
                """
        #others we need
        #K = self._K_unit_time
        #a = self.a
        num_pts = len(s)
        #cell_areas = self.cell_areas #nnodes long
        sed_flux = np.empty(num_pts, dtype=np.float64)
        sed_flux[:] = incision_iter
        weave.inline(code, ['num_pts', 's', 'r', 'sed_flux'])
        return sed_flux
        
    def erode(self, grid, dt, node_drainage_areas='planet_surface__drainage_area', 
                slopes_at_nodes=None, link_slopes=None, link_node_mapping='links_to_flow_reciever', 
                receiver='flow_reciever', upstream_order='upstream_ID_order', 
                slopes_from_elevs=None, W_if_used=None, Q_if_used=None, io=None):
        """
        Note this method must be passed both 'receiver' and 'upstream_order',
        either as strings for field access, or nnode- long arrays of the 
        relevant IDs. These are most easily used as the
        outputs from route_flow_dn.
        """
        try:
            r_in = grid.at_node[receiver]
        except FieldError:
            r_in = receiver
            assert r_in.size == grid.number_of_nodes
        try:
            s_in = grid.at_node[upstream_order]
        except FieldError:
            s_in = upstream_order
            assert s_in.size == grid.number_of_nodes
        #get the stream power part of the equation from the simple module:
        #_,_,simple_stream_power = self.simple_sp.erode(grid, dt, node_drainage_areas, slopes_at_nodes, link_slopes, link_node_mapping, slopes_from_elevs, W_if_used, Q_if_used, io)
        #slopes = self.simple_sp.slopes
        ######slopes needs to come across from simple_sp
        #this stuff deriving self.slopes is a direct port from simple_sp
        if not slopes_at_nodes:
            if slopes_from_elevs:
                if slopes_from_elevs == True:
                    node_z = grid.at_node['planet_surface__elevation']
                elif type(slopes_from_elevs) == str:
                    node_z = grid.at_node[slopes_from_elevs]
                else:
                    node_z = slopes_from_elevs
                S_links = (node_z[grid.node_index_at_link_head]-node_z[grid.node_index_at_link_tail])/grid.link_length
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
            except FieldError:
                try:
                    self.slopes = S_links[link_node_mapping]
                except IndexError:
                    #need to do the mapping on the fly.
                    #we're going to use the max slope (i.e., + or -) of *all* adjacent nodes.
                    #This isn't ideal. It should probably just be the outs...
                    #i.e., np.max(self.link_S_with_trailing_blank[grid.node_outlinks] AND -self.link_S_with_trailing_blank[grid.node_inlinks])
                    self.link_S_with_trailing_blank[:-1] = S_links
                    self.slopes = np.max(np.fabs(self.link_S_with_trailing_blank[grid.node_links]))
        else:
            try:
                self.slopes = grid.at_node[slopes_at_nodes]
            except FieldError:
                self.slopes = slopes_at_nodes
        slopes = self.slopes
        
        if type(node_drainage_areas)==str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        #######need to handle cells where Qs > Qc!!!

        #calc shear stress
        shear_stress = self.depth_prefactor*node_A**self.shear_area_power*slopes**0.7
        #calc Shields
        shields_stress = self.shields_prefactor*shear_stress
        capacity = self.sed_capacity_equation(grid, shields_stress, slopes)
        stability_ratio = 10. #initialize high, as if the exceptions trip, we will want the second loop to run first time
        try:
            #this assumes both the sed flux is defined in the last step, and that the fractional difference between the sed flux functions is <0.05 everywhere
            #if not, we bail and use pseudoimplicit convergence (slooooooooow)
            #note the 0.05 is basically introducing some slop into the solutions (<5%), which we trade off for speed
            #improve accuracy by turning down the timestep
            fqsqc = self.get_sed_flux_function(self.last_sed_flux/capacity)
            incision_iter = self._K_unit_time*dt*fqsqc*(shear_stress**self.a-self.threshold).clip(0.)*self.cell_areas #a volume
            sed_flux_iter = self.weave_iter_sed_flux(grid, incision_iter, r_in, s_in)
            self.last_sed_flux = (self.last_sed_flux, sed_flux_iter)/2.#moving window to pre-emptively smooth numerical oscillation
            fqsqc2 = self.get_sed_flux_function(self.last_sed_flux/capacity)
            if not np.allclose(fqsqc, fqsqc2, rtol=0.05):
                stability_ratio = np.amax(np.fabs(fqsqc/fqsqc2 - 1.)) #zero is good. >0.05 is bad.
                raise AssertionError
            incision_iter = self._K_unit_time*dt*fqsqc2*(shear_stress**self.a-self.threshold).clip(0.)*self.cell_areas #a volume
            sed_flux_iter = self.weave_iter_sed_flux(grid, incision_iter, r_in, s_in)
            self.last_sed_flux = (self.last_sed_flux, sed_flux_iter)/2.
        except (AttributeError, AssertionError):
            print "Pseudoimplicit solution this loop..."
            try:
                fqsqc = self.get_sed_flux_function(self.last_sed_flux/capacity)
            except AttributeError:
                fqsqc = self.get_sed_flux_function(0.1) #arbitrary uniform starting relative sed flux to get us going
            incision_iter = self._K_unit_time*dt*fqsqc*(shear_stress**self.a-self.threshold).clip(0.)*self.cell_areas #a volume
            sed_flux_iter = self.weave_iter_sed_flux(grid, incision_iter, r_in, s_in)
            try:
                self.last_sed_flux = (self.last_sed_flux, sed_flux_iter)/2.#moving window to pre-emptively smooth numerical oscillation
            except AttributeError:
                self.last_sed_flux = sed_flux_iter
            counter=0
            while stability_ratio > 0.05:
                fqsqc2 = self.get_sed_flux_function(self.last_sed_flux/capacity)
                incision_iter = self._K_unit_time*dt*fqsqc2*(shear_stress**self.a-self.threshold).clip(0.)*self.cell_areas #a volume
                sed_flux_iter = self.weave_iter_sed_flux(grid, incision_iter, r_in, s_in)
                stability_ratio = np.amax(np.fabs(fqsqc/fqsqc2 - 1.))
                fqsqc = fqsqc2[:]
                if counter>100:
                    print "Bailing pseudoimplicit loop after 100 iterations. stability ratio at stop is ", stability_ratio
                    break
                counter+=1
        