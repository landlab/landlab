#! /usr/env/python
"""

2D numerical model of shallow-water flow over topography read from a DEM, using
the Bates et al. (2010) algorithm for storage-cell inundation modeling.

1st stab at componentizing this routine, following GT, by DEJH, Oct 2013.

"""

from landlab import ModelParameterDictionary, CLOSED_BOUNDARY
import numpy as np
import six


class SurfaceFlowTransport(object):

    def __init__(self, grid, input_stream):
        #grid here is a true model field. i.e., we should be able to do grid.at_node['topographic__elevation']
        #input_stream is a text file, entered in format './my_file.txt'
        self.grid = grid

        inputs = ModelParameterDictionary(input_stream)
        self.n = inputs.read_float('n')              # roughness coefficient (Manning's n)
        self.g = inputs.read_float('g')               # gravitational acceleration (m/s2)
        self.alpha = inputs.read_float('alpha')       # time-step factor (ND; from Bates et al., 2010)
        self.tau_crit = inputs.read_float('tau_crit') # critical shear stress, pascals
        self.mpm = inputs.read_float('mpm')          # sed trans coefficient
        self.erode_start_time = inputs.read_float('erode_start_time') # allows an offset between when flow starts and when we start to allow erosion


        #test the necessary fields are all already present:
        try:
            self.z = grid.at_node['topographic__elevation']
        except:
            six.print_('elevations not found in grid!')
        try:
            self.h = grid.at_node['planet_surface__water_depth']
        except:
            six.print_('initial water depths not found in grid!')

        #build the internally necessary params:
        self.rhog = 9810.          # water unit weight, kg/m2s2 (N/m3)
        self.q = grid.create_active_link_array_zeros()       # unit discharge (m2/s)
        self.dhdt = grid.create_node_array_zeros()           # rate of water-depth change
        self.tau = grid.create_active_link_array_zeros()     # shear stress (Pa)
        self.qs = grid.create_active_link_array_zeros()      # sediment flux (m2/s)
        self.dqsds = grid.create_node_array_zeros()
        self.dzdt = self.dhdt
        self.dzaccum = grid.create_node_array_zeros()
        self.zm = grid.create_node_array_zeros()
        self.zm[:] = self.z[:]

    def set_and_return_dynamic_timestep(self):
        # Calculate time-step size for this iteration (Bates et al., eq 14)
        self.dtmax = self.alpha*min(self.grid.dx, self.grid.dy)/np.sqrt(self.g*np.amax(self.h))
        return self.dtmax

    def set_timestep(self, timestep_in):
        """
        This fn allows overriding of the inbuilt dynamic timestepping, if, e.g.,
        another component requires a shorter timestep.
        Function assumes you have already called self.
        set_and_return_dynamic_timestep() to establish what the min should be.
        """
        if timestep_in <= self.dtmax:
            self.dtmax = timestep_in
        else:
            raise RuntimeError('Attempting to manually set an unstable timestep! Abort! Abort!')

    def transport_sed(self, elapsed_time):
        #load data for speed and clarity:
        grid = self.grid
        z = self.z
        h = self.h
        g = self.g
        n = self.n
        tau_crit = self.tau_crit
        q = self.q
        qs = self.qs
        tau = self.tau
        rhog = self.rhog
        alpha = self.alpha
        mpm = self.mpm
        zm = self.zm
        dtmax = self.dtmax
        erode_start_time = self.erode_start_time
        ten_thirds = 10./3.

        interior_cells = np.where(grid.status_at_node != CLOSED_BOUNDARY)[0]


        # Calculate the effective flow depth at active links. Bates et al. 2010
        # recommend using the difference between the highest water-surface
        # and the highest bed elevation between each pair of cells.
        zmax = grid.max_of_link_end_node_values(z)
        w = h+z   # water-surface height
        wmax = grid.max_of_link_end_node_values(w)
        hflow = wmax - zmax

        # Calculate water-surface slopes
        water_surface_slope = grid.calculate_gradients_at_active_links(w)

        # Calculate the unit discharges (Bates et al., eq 11)
        q = (q-g*hflow*dtmax*water_surface_slope)/ \
            (1.+g*hflow*dtmax*n*n*abs(q)/(hflow**ten_thirds))

        # Calculate shear stress and sediment flux
        tau = -rhog*hflow*water_surface_slope
        tauex = abs(tau)-tau_crit
        tauex[np.where(tauex<0.)] = 0.
        qs = np.sign(tau)*mpm*pow(tauex,1.5)

        # Calculate water-flux divergence at nodes
        dqds = grid.calculate_flux_divergence_at_nodes(q)
        dqsds = grid.calculate_flux_divergence_at_nodes(qs)

        # Calculate rate of change of water depth
        dhdt = -dqds
        dzdt = -dqsds

        # Second time-step limiter (experimental): make sure you don't allow
        # water-depth to go negative
        excess_time = 0.
        if np.amin(dhdt) < 0.:
            shallowing_locations = np.where(dhdt<0.)
            time_to_drain = -h[shallowing_locations]/dhdt[shallowing_locations]
            dtmax2 = alpha*np.amin(time_to_drain)
            min_timestep_ratio = int(dtmax//dtmax2)
            if min_timestep_ratio: #any value other than 0, i.e., dtmax2<dtmax
                excess_time = dtmax/dtmax2 - min_timestep_ratio
                dt = dtmax2
            else:
                dt =  dtmax
        else:
            min_timestep_ratio = 0
            dt = dtmax

        #perform a loop if we had to subdivide the tstep above:
        for i in xrange(min_timestep_ratio+1):
            # Update the water-depth field
            h[interior_cells] += dhdt[interior_cells]*dt
            if elapsed_time >= erode_start_time:
                #z[interior_cells] = z[interior_cells] + dzdt[interior_cells]*dt
                zm[interior_cells] += dzdt[interior_cells]*dt
                #dzaccum[interior_cells] += dzdt[interior_cells]*dt
        if excess_time:
            h[interior_cells] += dhdt[interior_cells]*dt*excess_time
            if elapsed_time >= erode_start_time:
                zm[interior_cells] += dzdt[interior_cells]*dt*excess_time

        #Now we need to embed the results back into the object, and the grid
        self.z = z
        self.h = h
        self.q = q
        self.qs = qs
        self.tau = tau
        self.zm = zm
        self.grid['node']['topographic__elevation'] = zm
        self.grid['node']['planet_surface__water_depth'] = h

        return self.grid

