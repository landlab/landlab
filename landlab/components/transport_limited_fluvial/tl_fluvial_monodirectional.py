import numpy as np
from landlab import ModelParameterDictionary
from time import sleep
from scipy import weave
from scipy.weave.build_tools import CompileError

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
            
        '''
        self.grid = grid
        self.link_S_with_trailing_blank = np.zeros(grid.number_of_links+1) #needs to be filled with values in execution
        self.count_active_links = np.zeros_like(self.link_S_with_trailing_blank, dtype=int)
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
        self.relative_weight = (self.sed_density-self.fluid_density)/self.fluid_density*self.g #to accelerate MPM calcs
        self.excess_SG = (self.sed_density-self.fluid_density)/self.fluid_density
        self.rho_g = self.fluid_density*self.g     
        
        try:
            self.Qc = inputs.read_string('Qc')
        except MissingKeyError:
            raise MissingKeyError("Qc must be 'MPM' or a grid field name!")
        else:
            if self.Qc=='MPM':
                self.calc_cap_flag = True
            else:
                self.calc_cap_flag = False
                
        try:
            self.lamb_flag = inputs.read_bool('slope_sensitive_threshold')
        except:
            self.lamb_flag = False
        try:
            self.shields_crit = inputs.read_float('threshold_shields')
            self.set_threshold = True #flag for sed_flux_dep_incision to see if the threshold was manually set.
            print "Found a threshold to use: ", self.shields_crit
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
            
        #assume Manning's equation to set the power on A for shear stress:
        self.shear_area_power = 0.6*self._c*(1.-self._b)
        
        try:
            self.k_Q = inputs.read_float('k_Q')
        except MissingKeyError:
            self.depth_prefactor = self.rho_g*inputs.read_float('depth_equation_prefactor')
        else:
            self.k_w = inputs.read_float('k_w')
            mannings_n = inputs.read_float('mannings_n')
            if mannings_n<0. or mannings_n>0.2:
                print "***STOP. LOOK. THINK. You appear to have set Manning's n outside it's typical range. Did you mean it? Proceeding...***"
                sleep(2)
            self.depth_prefactor = self.rho_g*mannings_n*(self.k_Q**(1.-self._b)/self.k_w)**0.6
            ##Note the depth_prefactor we store already holds rho*g   

        try:
            self.C_MPM = inputs.read_float('C_MPM')
        except MissingKeyError:
            self.C_MPM = 1.
        try:
            self.shields_prefactor = 1./((self.sed_density-self.fluid_density)*self.g*self.Dchar_in)
            self.MPM_prefactor = 8.*self.C_MPM*np.sqrt(self.relative_weight*self.Dchar_in*self.Dchar_in*self.Dchar_in)
            self.MPM_prefactor_alt = 4.*self.g**(-2./3.)/self.excess_SG/self.fluid_density/self.sed_density
        except AttributeError:
            #have to set these manually as needed
            self.shields_prefactor_noD = 1./((self.sed_sensity-self.fluid_density)*self.g)

        self.cell_areas = np.empty(grid.number_of_nodes)
        self.cell_areas.fill(np.mean(grid.cell_areas))
        self.cell_areas[grid.cell_node] = grid.cell_areas
        
        
    def sed_capacity_equation(self, grid, shields_stress, slopes_at_nodes=None, areas_at_node=None):
        """
        slopes_at_nodes must be provided as nnodes array if slope_sensitive_threshold.
        This is a volume flux.
        At the moment, this only supports MPM for sed flux capacity calculation,
        but is set up to allow easy additions later if needed.
        """
        if self.lamb_flag:
            thresh = 0.15*slopes_at_nodes**0.25
        else:
            thresh = self.shields_crit
        if self.Qc == 'MPM':
            print "in capacity:"
            print 'prefactor ', self.MPM_prefactor
            print 'max Shields ', np.amax(shields_stress[self.grid.core_nodes])
            print 'thresh ', np.amax(thresh)
            try:
                capacity = self.MPM_prefactor*((shields_stress-thresh).clip(0.))**1.5
            except AttributeError:
                capacity = 8.*self.C_MPM*np.sqrt(self.relative_weight*self.Dchar*self.Dchar*self.Dchar)*((shields_stress-thresh).clip(0.))**1.5
        elif type(self.Qc)==str:
            try:
                capacity = grid.at_node[self.Qc]
            except:
                raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided input...')
        elif type(self.Qc)==np.ndarray:
            capacity=self.Qc
        else:
            raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided inputs...')
        capacity *= self.k_w*areas_at_node**self._b*31557600.
        return capacity #returned capacity is volume flux/yr
    
    def sed_capacity_equation_alt_method(self, grid, shear_stress, slopes_at_nodes=None, areas_at_node=None):
        """
        slopes_at_nodes must be provided as nnodes array if slope_sensitive_threshold.
        This is a volume flux.
        At the moment, this only supports MPM for sed flux capacity calculation,
        but is set up to allow easy additions later if needed.
        """
        if self.lamb_flag:
            thresh = 0.15*slopes_at_nodes**0.25
        else:
            thresh = self.shields_crit
        if self.Qc == 'MPM':
            capacity = self.MPM_prefactor_alt*((shear_stress-self.excess_SG*self.Dchar*self.rho_g*thresh).clip(0.))
        elif type(self.Qc)==str:
            try:
                capacity = grid.at_node[self.Qc]
            except:
                raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided input...')
        elif type(self.Qc)==np.ndarray:
            capacity=self.Qc
        else:
            raise TypeError('sed_flux_dep_incision does not know how to set sed transport capacity from the provided inputs...')
        capacity *= self.k_w*areas_at_node**self._b*31557600.
        return capacity #returned capacity is volume flux/yr
    
    
    def weave_iter_sed_trp_balance(self, grid, capacity, r, s, dt):
        '''
        Needs to be passed at node values for 'flow_receiver' and 
        'upstream_ID_order' (i.e., r and s). These are most likely the output
        values from route_flow_dn.
        Returns dz/dt, the rate of elevation change.
        '''
        code =  """
                int donor;
                int rcvr;
                for (int i = num_pts-1; i > -1; i--) {
                    donor = s[i];
                    rcvr = r[donor];
                    dz[donor] = (sed_into_cell[donor] - capacity[donor]*dt)/cell_areas[donor];
                    if (donor != rcvr) {
                        sed_into_cell[rcvr] += capacity[donor];
                    }
                }
                """
        num_pts = len(s)
        cell_areas = self.cell_areas #nnodes long
        sed_into_cell = np.zeros(capacity.size, dtype=np.float64) #this is, more accurately, a sed *flux* in
        dz = np.zeros_like(sed_into_cell)
        #print type(cell_areas), type(s), type(r), type(capacity), type(sed_into_cell)
        try:
            raise CompileError #force the python loop for debug
            weave.inline(code, ['num_pts', 'cell_areas', 's', 'r', 'capacity', 'sed_into_cell', 'dz', 'dt'])
        except CompileError:
            for donor in s[::-1]:
                rcvr = r[donor]
                #by the time we get to donor, no more sed will be added; it's "done":
    ########where should the next line be? Can an internally drained node erode down/deposit?
                dz[donor] = (sed_into_cell[donor]-capacity[rcvr]*dt)/cell_areas[donor] #if capacity is greater than sed_in, it will erode (-ve value here)
                #we look downstream to ensure continuity with the BCs at grid margins
                #print "dz: ", dz[donor]
                #capacity's worth of sed leaves donor & goes to rcvr
                if donor != rcvr:
                    #print "saved", donor, rcvr
                    #print "dz: ", dz[donor]
                    sed_into_cell[rcvr] += capacity[donor]
        #print dz
        return dz

        
    def erode(self, grid, dt, node_drainage_areas='drainage_area', 
                slopes_at_nodes=None, link_slopes=None, link_node_mapping='links_to_flow_receiver', 
                receiver='flow_receiver', upstream_order='upstream_ID_order', 
                slopes_from_elevs=None, W_if_used=None, Q_if_used=None,
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
        the component to look for the data in grid.at_node['topographic_elevation'];
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
        try:
            self.Dchar=self.Dchar_in
        except AttributeError:
            try:
                self.Dchar=grid.at_node[Dchar_if_used]
            except FieldError:
                assert type(Dchar_if_used)==np.ndarray
                self.Dchar=Dchar_if_used
            
        #Perform check on whether we use grid or direct fed data:
        if slopes_at_nodes==None:
            if slopes_from_elevs:
                if slopes_from_elevs == True:
                    node_z = grid.at_node['topographic_elevation']
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
                    #i.e., np.amax(self.link_S_with_trailing_blank[grid.node_outlinks] AND -self.link_S_with_trailing_blank[grid.node_inlinks])
                    self.link_S_with_trailing_blank[:-1] = S_links
                    self.slopes = np.amax(np.fabs(self.link_S_with_trailing_blank[grid.node_links]),axis=0)
        else:
            try:
                self.slopes = grid.at_node[slopes_at_nodes]
            except FieldError:
                self.slopes = slopes_at_nodes
        
        #print 'Max slope ', np.amax(self.slopes[self.grid.core_nodes])
        #print 'Mean slope ', np.mean(self.slopes[self.grid.core_nodes])
        if type(node_drainage_areas)==str:
            node_A = grid.at_node[node_drainage_areas]
        else:
            node_A = node_drainage_areas

        #print 'Area ', np.amax(node_A)
        #calc Shields
######need to add ability to do this with W or Q scaling specified
        shear_stress = self.depth_prefactor*node_A**self.shear_area_power*self.slopes**0.7
        try:
            shields_stress = self.shields_prefactor*shear_stress
        except AttributeError:
            shields_stress = self.shields_prefactor_noD/self.Dchar*shear_stress
        #print 'max_shear ', np.amax(shear_stress[self.grid.core_nodes])
        #capacity = self.sed_capacity_equation(grid, shields_stress, self.slopes, node_A)
        capacity = self.sed_capacity_equation_alt_method(grid, shear_stress, self.slopes, node_A)
        #print capacity
        #need to strip closed nodes from s:
        #...but in fact we use the mask on capacity, not a BC check, as some open nodes can still have no slopes defined (so are functionally closed) - i.e., the corners
        try:
            s_in = np.delete(s_in, np.where(np.in1d(s_in,np.where(capacity.mask),assume_unique=True))) #yuck! #this isn't very memory efficient...
        except AttributeError: #provided S array might not be masked...
            s_in = np.delete(s_in, np.where(self.grid.is_boundary(s_in, boundary_flag=4)))
        #print "r: ", r_in
        #print "s: ", s_in
        dz = self.weave_iter_sed_trp_balance(grid, capacity, r_in, s_in, dt)
        #print 'max inc rate ', np.amax(np.fabs(dz[self.grid.core_nodes]))
        active_nodes = grid.get_active_cell_node_ids()
        if io:
            try:
                io[active_nodes] += dz[active_nodes]
            except TypeError:
                if type(io)==str:
                    elev_name = io
            else:
                return grid, io, capacity
            
        else:
            elev_name = 'topographic_elevation'

        grid.at_node[elev_name][active_nodes] += dz[active_nodes]
        grid.at_node['sediment_flux_capacity'] = capacity
        
        return grid, grid.at_node[elev_name], capacity
        
