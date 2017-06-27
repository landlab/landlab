#! /usr/env/python
"""

GSA_PitFilled_SC.py

This is a driver file to run sample storms over the D4 pitfilled Spring Creek
DEM.

"""


from __future__ import print_function

from landlab.components import OverlandFlow, DetachmentLtdErosion, SoilInfiltrationGreenAmpt
from landlab.components.detachment_ltd_erosion.generate_erosion_by_depth_slope import DepthSlopeProductErosion
from landlab.utils.depth_dependent_roughness import depth_dependent_mannings_n
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.grid.raster_mappers import map_max_of_inlinks_to_node , map_max_of_outlinks_to_node
import time
import numpy as np
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid
from matplotlib import cm
import matplotlib as mpl
def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

   4         cmap: colormap instance, eg. cm.jet.
   5         N: number of colors.
   6
   7     Example
   8         x = resize(arange(100), (5,100))
   9         djet = cmap_discretize(cm.jet, 5)
  10         imshow(x, cmap=djet)
  11     """
    if type(cmap) == str:
        cmap = cm.get_cmap(cmap)
    colors_i = np.concatenate((np.linspace(0, 1., N), (0.,0.,0.,0.)))
    colors_rgba = cmap(colors_i)
    indices = np.linspace(0, 1., N+1)
    cdict = {}
    for ki,key in enumerate(('red','green','blue')):
        cdict[key] = [ (indices[i], colors_rgba[i-1,ki], colors_rgba[i,ki]) for i in xrange(N+1) ]
     # Return colormap object.
    return mpl.colors.LinearSegmentedColormap(cmap.name + "_%d"%N, cdict, 1024)
start_time = time.time()


DAs = [19405, 39016]
links= [110756] # 117782

filename_input = 'SC_30m_Mar17.asc'#SpringCreek_PitsFilled_30m.asc'

(rmg, z) = read_esri_ascii(filename_input)

elapsed_time = 0.0

## Model Run Time in seconds
model_run_time = 28800.0#43200.0

## Lists for saving data
discharge_at_outlet = []
gage_hydrographs = []
midstream_node = []
midstream_incision =[]
hydrograph_time_sec = []
hydrograph_time_hrs = []
outlet_incisograph = []

drainage_area = 'Pitfilled_SC_DrainageArea.txt'
drainage_area = np.loadtxt(drainage_area)
unique_da = np.unique(drainage_area)


## Setting initial fields...
rmg['node']['topographic__elevation'] = z
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)
rmg['node']['total_incision_depth'] = np.zeros(rmg.number_of_nodes)
rmg['node']['soil_water_infiltration__depth'] = np.zeros(rmg.number_of_nodes)
rmg['node']['soil_water_infiltration__depth'] += 10**-3
rmg.set_watershed_boundary_condition(z, -9999)
outlet_id = np.where(rmg.status_at_node == 1)[0]
link_to_sample = 110756
total_incision_depth = np.zeros(rmg.number_of_nodes)
depth_at_outlet=[]

storm_flag = 'MeanStorm'
save_flag = 'False'

if storm_flag == 'LargeStorm':
    rainfall_mmhr = 80.
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 3600.

if storm_flag == 'MeanStorm':
    rainfall_mmhr = 22.0
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 3600.
    
rmg['link']['mannings_n'] = np.ones(rmg.number_of_links)
rmg['link']['mannings_n'] *= 0.055
conductivity1 = 2.7778*(10**-7)
lowest_3_conductivity = 0.00000083333# 3 mm/hr
low_6_conductivity = 0.000001666666667 # 6 mm/hr
mid_9_conductivity = 0.0000025 # 9 mm/hr
high_conductivity =  0.0000033333 # 12 mm/hr

#mid_conductivity = (1.94*10**-6)#0.00001944444444444
#high_conductivity = 0.00002222222222

low_soil_moisture = 0.05
default_soil_moisture=0.15
mid_soil_moisture = 0.2
high_soil_moisture = 0.35

of = OverlandFlow(rmg, steep_slopes=True, alpha=0.15, mannings_n='mannings_n', h_init=0.00001)
dle = DepthSlopeProductErosion(rmg, k_e = (4.0 * (10**-9)))#, tau_crit=0.1)#.1)#232#704)
siga = SoilInfiltrationGreenAmpt(rmg, surface_water_minimum_depth=of.h_init,
                                 hydraulic_conductivity=mid_9_conductivity, 
                                 initial_soil_moisture_content=0.20, 
                                 soil_bulk_density=1700., rock_density=2700., 
                                 soil_type ='sand',    
                                 volume_fraction_coarse_fragments=0.6)
peak_I = np.zeros(rmg.number_of_nodes)
filenameS_template = r'/Users/Jordan/Desktop/Movies/WSS_movie_22mmhr_K6_SM35/{0:01d}.png'
filenameT_template = r'/Users/Jordan/Desktop/Movies/Tau_movie_22mmhr_K6_SM35/{0:01d}.png'
y= 0
outlet_tau_list = []
upstream_tau_list=[]
#print(rainfall_mmhr, 'Infil K = ', siga._Ks ,'n = 0.05', 'soil moisture =', round((0.37) - siga._Md, 2))
## Running the overland flow component.
while elapsed_time < model_run_time:
    
    of.dt = of.calc_time_step()
    
    #if (elapsed_time < storm_duration) and ((elapsed_time + of.dt) > storm_duration):
       # of.dt = storm_duration - (elapsed_time - 0.1)
    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.

    if elapsed_time < (storm_duration):

        of.rainfall_intensity =  rainfall_ms

    ## Then the elapsed time exceeds the storm duration, rainfall ceases.
    else:

        of.rainfall_intensity = 0.0
        
    depth_dependent_mannings_n(rmg, min_mannings_n=0.055, index_flow_depth=0.05)
    
    rmg['link']['mannings_n'] = rmg.map_min_of_link_nodes_to_link('mannings_n')

    of.overland_flow(dt=of.dt)

    node_slope = ((of.water_surface_slope[rmg.links_at_node] *
                                rmg.active_link_dirs_at_node))
    
    incision_Q = np.abs(of.q * rmg.dx)[rmg.links_at_node]

    rmg['node']['surface_water__discharge'] = (incision_Q[np.arange(len(
                                node_slope)), np.argmax(node_slope, axis=1)])

    node_slope = node_slope.max(axis=1)
   
    rmg['node']['water_surface__slope'] = node_slope

    siga.run_one_step(of.dt)
    dle.erode(of.dt, slope = 'water_surface__slope')
    
    
   # dle.erode(of.dt, slope = 'water_surface__slope')
    peak_I = np.maximum(np.abs(dle.E), peak_I)
    z[outlet_id] -= np.min(dle.dz[rmg.neighbors_at_node[outlet_id]])
    outlet_tau_list.append(dle.tau[55744])
    upstream_tau_list.append(dle.tau[DAs])
    hydrograph_time_sec.append(elapsed_time)
    depth_at_outlet.append(of.h[55744])
    discharge_at_outlet.append(rmg.at_node['surface_water__discharge'][55744])

    midstream_node.append(rmg.at_node['surface_water__discharge'][DAs])
    midstream_incision.append(np.abs(dle.E[DAs]))
    outlet_incisograph.append(dle.E[55744])
    total_incision_depth += (np.abs(dle.dz))

    y+=1
    elapsed_time += of.dt

sq_dx = rmg.dx ** 2
volume_at_cells = total_incision_depth * sq_dx
print('total incision volume: ', sum(volume_at_cells))

rmg['node']['total_incision_depth'] = total_incision_depth

calc_water_mass = round(np.abs((np.trapz(hydrograph_time_sec, (np.abs(
                    discharge_at_outlet) * rmg.dx)))), 2)
theoretical_water_mass = round(((max(drainage_area)) *
                    (rainfall_ms) * storm_duration), 2)
percent_error = round(((np.abs(calc_water_mass) - theoretical_water_mass) /
                    theoretical_water_mass * 100), 2)

print('Percent error: ', percent_error)
endtime = time.time()
print('\n', 'Overall model run time: ', round(endtime - start_time, 2), ' seconds, for alpha = ', of.alpha)


plt.figure(1)
#plt.title(rainfall_mmhr, '  ', siga._Ks, '   ', round(0.35 - siga._Md, 2))
hr = [x / 3600. for x in hydrograph_time_sec]
plt.plot(hr, np.abs(discharge_at_outlet)*rmg.dx)

plt.title('Outlet Hydrograph')
plt.xlim(0, 8)


print(max(np.abs(discharge_at_outlet)*rmg.dx))
#print(max(np.abs(discharge_at_outlet)*rmg.dx))
q_midstream = np.array(midstream_node)
i_midstream = np.array(midstream_incision)
#Ks = str(int(siga._Ks / (2.7778*(10**-7))) + 1)
#SM = str(int((round((0.37) - siga._Md, 2))*100.))
#string_for_save = str(int(rainfall_mmhr)) + 'mmhr_K' + Ks + '_SM' + SM 

rmg.at_node['peak_I'] = peak_I

tau_flag = 'TauCrit' # 'TauCrit'

#write_esri_ascii('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_max_erosion_rate_grid.asc', rmg, 'peak_I')
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_outlet_incisograph.txt', outlet_incisograph)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_midstream_incisograph.txt', i_midstream[:,1])
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_upstream_incisograph.txt', i_midstream[:,0])
#write_esri_ascii('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_total_eroded_depth_grid.asc', rmg, 'total_incision_depth')
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_total_sed_yield.txt', volume_at_cells)
upstream_tau_list = np.array(upstream_tau_list)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_tau_graph.txt', outlet_tau_list)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I80mmhr/' + tau_flag + '/Q' + string_for_save + '_tau_time.txt', hr)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_tau_graph_upstream.txt', upstream_tau_list[:,0])
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/IncisionResults/I22mmhr/' + tau_flag + '/Q' + string_for_save + '_tau_graph_midstream.txt', upstream_tau_list[:,1])








#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I10mmhr_1hr/K1mmhr_5cm/Q10mmhr_60m_ddn5cm_SIGA_35.txt', np.abs(discharge_at_outlet)*rmg.dx)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I10mmhr_1hr/K1mmhr_5cm/T10mmhr_60m_ddn5cm_SIGA_35.txt', hr)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I10mmhr_1hr/K1mmhr_5cm/Q10mmhr_60m_ddn5cm_SIGA_35_upstream.txt', q_midstream[:,0])
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I10mmhr_1hr/K1mmhr_5cm/Q10mmhr_60m_ddn5cm_SIGA_35_midstream.txt', q_midstream[:,1])
###
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I80mmhr_1hr/K1mmhr_5cm/Q22mmhr_60m_ddn5cm_NOSIGA.txt', np.abs(discharge_at_outlet)*rmg.dx)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I80mmhr_1hr/K1mmhr_5cm/T22mmhr_60m_ddn5cm_NOSIGA.txt', hr)
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I80mmhr_1hr/K1mmhr_5cm/Q22mmhr_60m_ddn5cm_NOSIGA_upstream.txt', q_midstream[:,0])
#np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/I80mmhr_1hr/K1mmhr_5cm/Q22mmhr_60m_ddn5cm_NOSIGA_midstream.txt', q_midstream[:,1])