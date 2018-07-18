#! /usr/bin/env python
import numpy as np
from scipy.interpolate import interp1d
from landlab import RasterModelGrid, Component


class Fluvial(Component):

    _name = 'Sand Percent Calculator'

    _time_units = 'y'

    _input_var_names = ()

    _output_var_names = (
        'delta_sediment_sand__volume_fraction',
    )

    _var_units = {
        'delta_sediment_sand__volume_fraction': '%',
    }

    _var_mapping = {
        'delta_sediment_sand__volume_fraction': 'grid',
    }

    _var_doc = {
        'delta_sediment_sand__volume_fraction': 'delta sand fraction',
    }

    def __init__(self, grid, sand_frac, start=0., **kwds):
        """Generate percent sand/mud for fluvial section.
        
        Parameters
        ----------
        grid: ModelGrid
            A landlab grid.
        sand_frac: str
            Name of csv-formatted sea-level file.

        """
        super(FluvialSand, self).__init__(grid, **kwds)
        
        # fixed parameters
        sand_grain = .001      # grain size = 1 mm 
        alpha = 10.            # ratio of channel depth to channel belt thickness  */
        beta = .1              # beta*h is flow depth of flood, beta = .1 to .5   */
        lambdap = .30
        flood_period = 10.     #  recurrence time of floods ~1-10 y  */
        basin_width = 5000.    #  Basin width or river spacing of 20 km */
        basin_length = 500000. #length for downstream increase in diffusion */

        #Upstream boundary conditions  */
        mud_vol = sediment_load/(1.-sand_frac)
        sand_vol = sediment_load
        qs = 10.*sqrt(9.8*(sand_density/1000.-1.))*(sand_grain ** 1.5);
            # m^2/s  units */

        # upstream diffusivity is set by equilibrium slope */
        diffusion = sediment_load/plain_slope
        qw = diffusion/0.61
        conc_mud[0] = mud_vol/qw

        channel_width = sand_vol*basin_width/qs/31536000.
        
        x = grid.x_of_node.reshape(grid.shape)
        land = x[1] < init_shore
        slope = np.gradient(z[1, land]) / dx
        channel_depth[land] = (sand_density-1000.)/1000.*sand_grain/slope[land]

        # Type of channelization */
        if channel_width/channel_depth >75.:
            epsilon = 0.8;  # braided 0.3-0.5  */
        else:
            epsilon = 0.125; # meandering  0.1-0.15  */
        width_cb = channel_width/epsilon
        
        #Original: r_cb = (model.new_height[i]-model.height[i]+model.d_sl);        
        r_cb = dz = self.grid.at_node['bedrock_surface__increment_of_elevation']
        # original: r_b = model.thickness[i];
        r_b = self.grid.at_node['bedrock_surface__increment_of_elevation']

        for i  in range (1, land):

            # all rates  per timestep */
            # channelbelt deposition  */
            if r_cb[i] < 0.: r_cb[i] = 0.

            # floodplain deposition  */
            r_fp[i] = beta*channel_depth/flood_period*conc_mud[i]*dt*1000.
            if r_fp > r_cb: r_fp = r_cb

            #Find avulsion rate and sand density   */
            
            if dz[i] > 0.: 

                bigN = alpha*(r_cb[i] - r_fp)/r_b[i];
                if bigN > 1.: r_cb *=  bigN;  
            # rate is bigger because of avulsions */
            
                if r_cb <= 0. :
                    r_cb = 0.
                    self.grid.percent_sand[i] = 1.

                else :
                    bigN = alpha*(r_cb - r_fp)/r_b;
                    self.grid.percent_sand[i] = 1.- (1-width_cb/basin_width)*exp(-1.*width_cb/basin_width*bigN);
            else :
                self.grid.percent_sand[i] = 0; #NULL;*/

            # adjust parameters for next downstream point */
            if dz[i] > 0.:
                #printf ("%d %f %f %f %f %f %f %f %f %f %f %f\n",i, model.fc[i], bigN, model.thickness[i], r_cb,r_fp, r_b, 
                #epsilon, conc_mud[i], qw, channel_width, channel_depth);
                sand_vol -= self.grid.percent_sand[i]* self.grid.spacing*(dz[i]*dz[i+1])/2/dt 
                mud_vol  -= (1.-self.grid.percent_sand[i])*self.grid.spacing*(dz[i]*dz[i+1])/2/dt
            diffusion = sediment_load/slope[i]*(1.+i*self.grid.spacing/basin_length) #question is i correct?
            qw = diffusion/0.61;
            conc_mud[i+1] = mud_vol/qw;
            channel_depth = (sand_density-1000.)*sand_grain/slope[1]
            if channel_width/channel_depth >75.:
                epsilon = 0.8;  # braided 0.3-0.5  */
            else :
                epsilon = 0.125; # meandering  0.1-0.15  */
            width_cb = channel_width/epsilon
