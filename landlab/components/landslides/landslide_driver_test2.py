# ##!/usr/bin/env python
"""
THIS IS TEST CODE!!!

Purpose is to gather input data, run landslide component, and visualize/store
data and outputs.  Input data is provide by the user and consists of elevation
from a DEM to provide topographic traits such as slope and contributing area.
User also supplies soil characteristics from a soil survey, land cover, or
other source, including transmissivity, cohesion, internal angle of friction,
density, and thickness.

@author: rstrauch - Univerity of Washington
Created on Thu Aug 20 16:47:11 2015
Last edit Aug 28, 2015
"""

# %% Import libraries and components
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
# from landlab import RasterModelGrid
from landlab.plot.imshow import imshow_node_grid
from landlab.io import read_esri_ascii
from statsmodels.distributions.empirical_distribution import ECDF
from collections import defaultdict
from landslide_component_test import Factor_of_Safety
# from landslide_component_v01 import Factor_of_Safety

# %% Instantiate classes and objects
# Values for TESTING
number_of_interations = 2500
contributing_area = 400                                            # [m]
slope = 30                                                   # [degrees]
soil_transmissivity__mean = 10                                # [m2/day]
combined_cohesion__mean = 17000                              # [Pa or kg/m-s2]
soil_internal_angle_friction__mean = 34                     # [degrees]
soil_density	= 2000	                                            # density of soil [kg/m3]
soil_thickness__mean = 1.2                                     # [m]

indir_template = 'H:\Big_Data\VIC_Fluxes\NOCA_fluxes\\test_fluxes\{}'           # Location of inputs
outdir_template = 'H:\Big_Data\VIC_Fluxes\NOCA_fluxes\\test_fluxes\output\{}' # where you want your output stored
# static_fields = ['Site_traits']                                                 # Location of static traits, like DEM
model_simulations = ['data_sum']                                                # Location of hydrology model simulations
d = defaultdict(list)

# plotting parameters
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['legend.fontsize'] = 15
# %% load or create the data in fields and Set boundary conditions

# load in DEM from ArcGIS
(mg, z) = read_esri_ascii('boston_elev.txt', name='topographic__elevation')
mg.at_node.keys()     # loads DEM grid with elevation
mg.set_nodata_nodes_to_closed(mg.at_node['topographic__elevation'], -9999)      # set boundary conditions closed where not data
# calculate slope from DEM
#slp = mg.calculate_slope_aspect_at_nodes_Burrough(ids=None, vals='topographic__elevation')
#mg.add_field('node', 'slope', slp)                                       # unitless...get theta by np.arctan(theta)
## Load contributing area from text file...unless we can use route_flow component???
#(mg1,ca) = read_esri_ascii('boston_ca.txt')
#mg.add_field('node', 'contributing_area', ca, units='m')
## load transmissivity from text file
#(mg1,T) = read_esri_ascii('boston_T.txt')
#mg.add_field('node', 'soil_transmissivity__mean', T, units='m2/day')
## load cohesion from text file
#(mg1,C) = read_esri_ascii('boston_C.txt')
#mg.add_field('node', 'combined_cohesion__mean', C, units='Pa')
## Load internal angle of friction from text file
#(mg1,phi) = read_esri_ascii('boston_phi.txt')
#mg.add_field('node', 'soil_internal_angle_friction__mean', phi, units='degrees')
## set soil density value and assign to all nodes
#rho = mg.add_ones('node', 'soil_density')
#mg.at_node['soil_density'][:] = 2000
## load soil thickness from text file
#(mg1,hs) = read_esri_ascii('boston_hs.txt')
#mg.add_field('node', 'soil_thickness__mean', hs, units='m')


#%% Run Factor-of-Safety model in landslide component

# should this next part be in the component???  but you do this for each 
def write_files_data(df, outdir, hydro_data_file):                             # writes results data to file
    out_path = os.path.join(outdir, hydro_data_file)
    np.savetxt(out_path,df,fmt=('%8.3f','%8.3f','%8.3f','%8.3f'))
    
# Run component for each hydrology model simulation of recharge
for simulation in model_simulations:
    indir = indir_template.format(simulation)
    outdir = outdir_template.format(simulation)
    for hydro_data_file in os.listdir(indir):                                  # Data file from hydrology simulation
        lat = hydro_data_file.split('_')[1]
        lon = hydro_data_file.split('_')[2]
        in_path = os.path.join(indir, hydro_data_file)
        daily_recharge = np.genfromtxt(in_path, usecols=3)                     # put annual max daily recharge into np array
        print daily_recharge        
        number_of_data = len(daily_recharge)        
        print number_of_data
#        x_vals = daily_recharge
#        print len(x_vals)
#        x_vals.sort()
#        Fx=ECDF(x_vals)
        Fx=ECDF(daily_recharge)             # probability associated with daily recharge values
        FS = Factor_of_Safety(number_of_interations, contributing_area, slope,
            soil_transmissivity__mean, combined_cohesion__mean, 
            soil_internal_angle_friction__mean, soil_density, soil_thickness__mean,
            daily_recharge, Fx(daily_recharge), number_of_data)   # parameters & data passed to class
        data_extract = np.transpose([FS.daily_recharge, FS.WI__mean, FS.FS__mean,
                FS.Probability_of_failure])
        df = pd.DataFrame(data_extract)
# load Probability of failure to grid
print FS.WI__mean
print FS.FS__mean
print FS.Probability_of_failure
#mg.add_field('node', 'Probability_of_failure', FS.Probability_of_failure)

#        write_files_data(df, outdir, hydro_data_file) 

# calculate contributing area
# import from ArcGIS...flow direction and CA


#%% Finalize and output

plt.figure('Elevations from the DEM')  # new fig, with a title
imshow_node_grid(mg, 'topographic__elevation', cmap='terrain', grid_units=['30m','30m'])
#plt.figure('Probability of Failure')
#imshow_node_grid(mg, 'Probability_of_failure', cmap='OrRd')
#%% run the main() function

#if __name__ == "__main__":
#    main()

#%% check the following
# units in degrees or radians everywhere for slopes/angles?