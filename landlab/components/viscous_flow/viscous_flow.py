#! /usr/env/python
"""

Component that models Poiseuille flow

Created March 2016 KAK


"""




from __future__ import print_function

import numpy as np
from six.moves import range

from landlab import ModelParameterDictionary, Component, FieldError
from landlab import create_and_initialize_grid
from landlab.core.model_parameter_dictionary import MissingKeyError


#_VERSION = 'make_all_data'
#_VERSION = 'explicit'



class ViscousFlowModel(Component):
    """
    This component implements Poiseuille flow of a field in the supplied
    ModelGrid.

    This components requires the following parameters be set in the input file,
    *input_stream*, set in the component initialization:

        'linear_diffusivity', the diffusivity to use

    Optional inputs are:

        'uplift_rate', if you want this component to include the uplift
            internally
        'dt', the model timestep (assumed constant)
        'values_to_diffuse', a string giving the name of the grid field
        containing the data to diffuse.

    Supply *dt* to the diffuser through the diffuse() argument.
    This allows you to set a dynamic timestep for this class.
    If 'values_to_diffuse' is not provided, defaults to
    'topographic__elevation'.

    No particular units are necessary where they are not specified, as long as
    all units are internally consistent.

    The component takes *grid*, the ModelGrid object, and (optionally)
    *current_time* and *input_stream*. If *current_time* is not set, it defaults
    to 0.0. If *input_stream* is not set in instantiation of the class,
    :func:`initialize` with *input_stream* as in input must be called instead.
    *Input_stream* is the filename of (& optionally, path to) the parameter
    file.

    At the moment, this diffuser can only work with constant diffusivity.
    Spatially variable diffusivity hopefully coming soon.

    The primary method of this class is :func:`diffuse`.
    
    ==========================
    
    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> rg = RasterModelGrid(3, 4)
    >>> vfm = ViscousFlowModel(rg)
    """

    _name = 'ViscousFlowModel'

    _input_var_names = set(['salt__thickness',
                            'topographic__elevation',
                            'salt__bottom'
                            ])

    _output_var_names = set(['topographic__elevation',
                             'salt__discharge',
                             'pressure__gradient',
                             'salt__velocity',
                             'oveburden__thickness'
                             ])

    _var_units = {'topographic__elevation' : 'm',
                  'salt__thickness' : 'm',
                  'salt__bottom' : 'm',
                  'salt__discharge' : 'm^2/s',
                  'pressure__gradient' : 'Pa/m',
                  'salt__velocity' : 'm/s',
                  'overburden__thickness' : 'm'
                  }

    _var_mapping = {'topographic__elevation' : 'node',
                    'salt__thickness' : 'node',
                    'salt__bottom' : 'node',
                    'overburden__thickness' : 'node',
                    'salt__discharge' : 'link',
                    'pressure__gradient' : 'link',
                    'salt__velocity' : 'link',
                    }

    _var_doc = {'topographic__elevation' : 'Land surface topographic elevation; can be overwritten in initialization',
                 'salt__thickness' : 'Vertical thickness of subsurface salt',
                 'salt__bottom' : 'elevation of the base of the salt layer',
                 'salt__discharge' : 'salt flux calculated using Poiseiulle flow',
                 'pressure__gradient' : 'Pressure gradient in the middle of salt layer',
                 'salt__velocity' : 'salt velocity in the x direction',
                 'overburden__thickness' : 'rock thickness above the salt, elevation minus top salt'
                  }

    def __init__(self, grid, viscosity=1e18, rho_overburden=2400.0, 
                 rho_salt=2200.0, slope_salt=0, g=9.81, 
                 initial_salt_thickness=300.0, salt_bottom=-500.0):
        self._grid = grid
        self.viscosity = viscosity
        self.rho_overburden = rho_overburden
        self.rho_salt = rho_salt
        self.slope_salt = slope_salt*np.pi/180
        self.g = g
        
        
        for name in self._var_mapping:
            if self._var_mapping[name] == 'node':
                if name not in self.grid.at_node:
                    self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._var_mapping:
            if self._var_mapping[name] == 'link':
                if name not in self.grid.at_link:
                    self.grid.add_zeros('link', name, units=self._var_units[name])

        self.h = self.grid.at_node['salt__thickness']
        self.h[:] = initial_salt_thickness
        self.z = self.grid.at_node['topographic__elevation']
        self.Q = self.grid.at_link['salt__discharge']
        self.dPdx = self.grid.at_link['pressure__gradient']
        self.vel = self.grid.at_link['salt__velocity']
        self.salt_bottom = self.grid.at_node['salt__bottom']
        self.salt_bottom[:] = salt_bottom
        self.T = self.grid.at_node['overburden__thickness']
        self.T[:] = self.z - (self.salt_bottom + self.h)
        #print(self.grid.at_node)
        #print(self.grid.at_link)


    def input_timestep(self, timestep_in):
        """
        Allows the user to set a dynamic (evolving) timestep manually as part of
        a run loop.
        """
        self.timestep_in = timestep_in
        self.tstep_ratio = timestep_in/self.dt


    def run_one_step_explicit(self, mg, z, g, qs, dqsds, dzdt, delt):

        # Take the smaller of delt or built-in time-step size self.dt
        dt = min(self.dt, delt)

        # Calculate the gradients and sediment fluxes
        g = mg.calculate_gradients_at_active_links(z)
        qs = -self.kd*g

        # Calculate the net deposition/erosion rate at each node
        dqsds = mg.calculate_flux_divergence_at_nodes(qs)

        # Calculate the total rate of elevation change
        dzdt = self.uplift_rate - dqsds

        # Update the elevations
        z[self.interior_cells] = z[self.interior_cells] \
                                 + dzdt[self.interior_cells] * dt

        # Update current time and return it
        self.current_time += dt

        return z, g, qs, dqsds, dzdt


    def run_one_step_internal(self, delt):

        # Take the smaller of delt or built-in time-step size self.dt
        dt = min(self.dt, delt)

        # Calculate the gradients and sediment fluxes
        self.g = self._grid.calculate_gradients_at_active_links(self.z)
        self.qs = -self.kd*self.g

        # Calculate the net deposition/erosion rate at each node
        self.dqsds = self._grid.calculate_flux_divergence_at_nodes(self.qs)

        # Calculate the total rate of elevation change
        dzdt = self.uplift_rate - self.dqsds

        # Update the elevations
        self.z[self.interior_cells] += dzdt[self.interior_cells] * dt

        # Update current time and return it
        self.current_time += dt
        return self.current_time

    def flow(self, dt):
        """
        This is the primary method of the class. Call it to perform an iteration
        of the model. Takes *dt*, the current timestep.

        The modelgrid must contain the field to diffuse, which defaults to
        'topographic__elevation'. This can be overridden with the
        values_to_diffuse property in the input file.

        See the class docstring for a list of the other properties necessary
        in the input file for this component to run.

        To improve stability, this component can incorporate uplift into its
        internal mechanism. To use this, set *internal_uplift* to True, and . If you only have one module that requires this, do not add
        uplift manually in your loop; this method will include uplift
        automatically. If more than one of your components has this requirement,
        set *num_uplift_implicit_comps* to the total number of components that
        do.

        You can suppress this behaviour by setting *internal_uplift* to False.

        """
        # Take the smaller of delt or built-in time-step size self.dt
        # self.tstep_ratio = dt/self.dt
        # repeats = int(self.tstep_ratio//1.)
        # extra_time = self.tstep_ratio-repeats
        # z = self._grid.at_node[self.values_to_diffuse]

        # core_nodes = self._grid.node_at_core_cell
#
#        for i in range (1):
#            P_base = self.rho_overburden*self.T*self.g + self.rho_salt*self.h*self.g
#            dPbase_dx = self._grid.calculate_gradients_at_active_links(P_base) # function will change 
#            #dPbase_dx = self._grid.calc_grad_at_link(P_base) # function will change: HERE IS NEW VER!
#            h_link_iso = self._grid.map_value_at_max_node_to_link(P_base, 'salt__thickness')
#            #Q_iso[self._grid.active_links] = -(1/(12*self.viscosity))*(dPbase_dx[self._grid.active_links])*h_link_iso[self._grid.active_links]**3
#            self.Q[self._grid.active_links] = -(1/(12*self.viscosity))*dPbase_dx*h_link_iso[self._grid.active_links]**3
#            #dhdt_iso = -self._grid.calc_flux_div_at_node(Q_iso) # FOR NEW VER!
#            dhdt_iso = -self._grid.calculate_flux_divergence_at_nodes(self.Q[self._grid.active_links])
#            deltah_iso = dhdt_iso*dt
#            self.h[self._grid.core_nodes] = self.h[self._grid.core_nodes] + deltah_iso[self._grid.core_nodes]
#        return self._grid
        
        for i in range(1):
            salt_top = self.salt_bottom + self.h
            self.T[:] = self.z - salt_top
            # print('T=')
            # print(self.T)
            # print(self._grid.at_node['overburden__thickness'])

#            # Isostatic equilibrium
#            P_base = self.rho_overburden*self.g*self.T + self.rho_salt*self.g*self.h
#            P_iso = np.mean(P_base)
#            self.h[:] = (P_iso - self.rho_overburden*self.g*self.T)/(self.rho_salt*self.g)  
#            self.h[self.h < 0] = 0
#            print('h=')
#            print(self.h)
#            # print(self._grid.at_node['salt__thickness'])

            # Calculate pressure, pressure gradient, flux, dhdt, velocities
            # update h, T
            P = self.rho_overburden*self.T*self.g + self.rho_salt*(self.h/2)*self.g
            # print('P=')
            # print(P)
            # include salt in this?
            self.dPdx[self._grid.active_links] = self._grid.calculate_gradients_at_active_links(P)
            # print('dpdx=')
            # print(self.dPdx)
            # h_link=self._grid.map_mean_of_link_nodes_to_link('salt__thickness')
            h_link = self._grid.map_value_at_max_node_to_link(P, 'salt__thickness')
            # print('hlink=')
            # print(h_link)
            self.Q[self._grid.active_links] = -(1/(12*self.viscosity))*(self.dPdx[self._grid.active_links]*h_link[self._grid.active_links]**3)
            # print('Q=')
            # print(self.Q)            
            uavg = self.Q/h_link  # average horizontal velocity
            dhdt = -self._grid.calculate_flux_divergence_at_nodes(self.Q[self._grid.active_links])
            deltah = dhdt*dt
            # need to look at dt from earlier in code
            
            # print('h=')
            # print(self.h)
            # print(self._grid.at_node['salt__thickness'])
            # print('z=')
            # print(self.z)
            # print(self._grid.at_node['topographic__elevation'])
            
            self.h[self._grid.core_nodes] = self.h[self._grid.core_nodes] + deltah[self._grid.core_nodes]
            self.z[self._grid.core_nodes] = self.z[self._grid.core_nodes] + deltah[self._grid.core_nodes]

        #return the grid
        return self._grid


    def run_until_explicit(self, mg, t, z, g, qs, dqsds, dzdt):

        while self.current_time < t:
            remaining_time = t - self.current_time
            z, g, qs, dqsds, dzdt = self.run_one_step_explicit(mg, z, g, qs, dqsds, dzdt, remaining_time)

        return z, g, qs, dqsds, dzdt


    def run_until_internal(self, t):

        while self.current_time < t:
            remaining_time = t - self.current_time
            self.run_one_step_internal(remaining_time)


    def run_until(self, t):  # this is just a temporary duplicate

        while self.current_time < t:
            remaining_time = t - self.current_time
            self.run_one_step_internal(remaining_time)


    def get_time_step(self):
        """
        Returns time-step size.
        """
        return self.dt

    @property
    def time_step(self):
        """
        Returns time-step size (as a property).
        """
        return self.dt


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    