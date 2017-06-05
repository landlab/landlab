#! /usr/env/python
"""

GSA_PitFilled_SC.py

This is a driver file to run sample storms over the D4 pitfilled Spring Creek
DEM.

"""
from __future__ import print_function

from landlab.components import OverlandFlow, DetachmentLtdErosion, SoilInfiltrationGreenAmpt
#from landlab.components.detachment_ltd_erosion.generate_detachment_ltd_erosion import DetachmentLtdErosion
from landlab.utils.depth_dependent_roughness import depth_dependent_mannings_n
from landlab.io import read_esri_ascii, write_esri_ascii
from landlab.grid.raster_mappers import map_max_of_inlinks_to_node , map_max_of_outlinks_to_node
import time
import numpy as np
from matplotlib import pyplot as plt
from landlab.plot import imshow_grid

start_time = time.time()


DAs = [53161, 7913, 46890, 46367, 24091, 53977, 57293, 17503, 35487]
links= [110756] # 117782

filename_input = 'SpringCreek_PitsFilled_30m.asc'

(rmg, z) = read_esri_ascii(filename_input)

elapsed_time = 0.0

## Model Run Time in seconds
model_run_time = 43200.0

## Lists for saving data
discharge_at_outlet = []
gage_hydrographs = []
#discharge_at_outlet2 = []
#discharge_at_gage =  []
hydrograph_time_sec = []
hydrograph_time_hrs = []
other_hydrographs = []

drainage_area = 'Pitfilled_SC_DrainageArea.txt'
drainage_area = np.loadtxt(drainage_area)
unique_da = np.unique(drainage_area)


## Setting initial fields...
rmg['node']['topographic__elevation'] = z
rmg['node']['water_surface__slope'] = np.zeros(rmg.number_of_nodes)
rmg['node']['total_incision_depth'] = np.zeros(rmg.number_of_nodes)
rmg['node']['soil_water_infiltration__depth'] = np.ones(rmg.number_of_nodes)
rmg['node']['soil_water_infiltration__depth'] *= 0.2
rmg.set_watershed_boundary_condition(z, -9999)

link_to_sample = 110756
#link_to_sample_2 = 110484
#gage_link = 117782


storm_flag = 'ThirdStorm'
save_flag = 'False'

if storm_flag == 'FirstStorm':
    rainfall_mmhr = 15.0
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 7200.
    filenameH_template = r'/Users/Jordan/Desktop/VideoImages/15mmhr/Depth/'
    filenameQ_template = r'/Users/Jordan/Desktop/VideoImages/15mmhr/Discharge/'
    filenameErate_template = r'/Users/Jordan/Desktop/VideoImages/15mmhr/ErosionRate/'
    filenameE_template = r'/Users/Jordan/Desktop/VideoImages/15mmhr/Erosion/'
#if storm_flag == 'SecondStorm':
#    rainfall_mmhr = 20.0 # * 2.)
#    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
#    storm_duration = 10800.
if storm_flag == 'ThirdStorm':
    rainfall_mmhr = 25.0
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 3600.
    filenameH_template = r'/Users/Jordan/Desktop/VideoImages/25mmhr/Depth/'
    filenameQ_template = r'/Users/Jordan/Desktop/VideoImages/25mmhr/Discharge/'
    filenameErate_template = r'/Users/Jordan/Desktop/VideoImages/25mmhr/ErosionRate/'
    filenameE_template = r'/Users/Jordan/Desktop/VideoImages/25mmhr/Erosion/'
#if storm_flag == 'FourthStorm':
#    rainfall_mmhr = 30.0
#    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
#    storm_duration = 10800.
if storm_flag == 'FifthStorm':
    rainfall_mmhr = 35.0
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 3600.
    filenameH_template = '/Users/Jordan/Desktop/VideoImages/35mmhr/Depth/'
    filenameQ_template = '/Users/Jordan/Desktop/VideoImages/35mmhr/Discharge/'
    filenameErate_template = r'/Users/Jordan/Desktop/VideoImages/35mmhr/ErosionRate/'
    filenameE_template = r'/Users/Jordan/Desktop/VideoImages/35mmhr/Erosion/'
if storm_flag == 'UngagedStorm':
#if storm_flag == 'UngagedStorm':
#if storm_flag == 'UngagedStorm':
    rainfall_mmhr =50.0
    rainfall_ms = rainfall_mmhr * 2.77778 * (10**-7)
    storm_duration = 10800.
#K_erodibility = placeholder_value

of = OverlandFlow(rmg, steep_slopes=True, alpha=0.2, mannings_n=0.055)
dle = DetachmentLtdErosion(rmg, K_sp = (2.5 * (10**-6)))
#siga = SoilInfiltrationGreenAmpt(rmg,
#                                 hydraulic_conductivity=1.94*(10**-5), 
#                                 soil_bulk_density=1700., rock_density=2700., 
#                                 initial_soil_moisture_content=1.,
#                                 surface_water_minimum_depth = 10**-8,
#                                 soil_type ='sand',    
#                                 volume_fraction_coarse_fragments=0.6)

#of.dt = 1
print(rainfall_mmhr, ' test')
old_q = of.q
peak_q = np.zeros(rmg.number_of_links)
    dle = DepthSlopeProductErosion()
#old_incision_rate = dle.I
#peak_incision_rate = np.zeros(rmg.number_of_nodes)
total_incision_depth = np.zeros(rmg.number_of_nodes)
i=0

u = []
v=[]
## Running the overland flow component.
while elapsed_time < model_run_time:
    print("start loop: ", max(rmg['node']['surface_water__depth']))
    of.dt = of.calc_time_step()
    if (elapsed_time < storm_duration) and ((elapsed_time + of.dt) > storm_duration):
        print('in loop')
        of.dt = storm_duration - (elapsed_time - 0.1)
    ## The storm starts when the model starts. While the elapsed time is less
    ## than the storm duration, we add water to the system as rainfall.
    if elapsed_time < (storm_duration):

        of.rainfall_intensity =  rainfall_ms

    ## Then the elapsed time exceeds the storm duration, rainfall ceases.
    else:

        of.rainfall_intensity = 0.0

    ## Generating overland flow based on the deAlmeida solution.
    of.overland_flow()# = 1.)#of.dt)
    print("after OF:", max(of.h))
    node_slope = (of.water_surface_slope[rmg.links_at_node] * 
                                          rmg.link_dirs_at_node)

    incision_Q = np.abs(of.q * rmg.dx)[rmg.links_at_node]

    rmg['node']['surface_water__discharge'] = (incision_Q[np.arange(len(
                                node_slope)), np.argmax(node_slope, axis=1)])
    

    node_slope = node_slope.max(axis=1)
  #  np.savetxt('/Users/Jordan/Desktop/pre_depths.txt', rmg['node']['surface_water__depth'])
    rmg['node']['water_surface__slope'] = node_slope

    
    #rmg['node']['total_I'] += np.abs(dle.I) * of.dt
    #u = rmg['node']['surface_water__depth']
    #siga.run_one_step(of.dt)
   #np.savetxt('/Users/Jordan/Desktop/post_depths.txt', rmg['node']['surface_water__depth'])
   # dle.erode(rmg, slope = 'water_surface__slope')
   # print('DLE: ' , dle.I[6009], 'OF.H: ', of.h[6009], 'Q:', rmg.at_node['surface_water__discharge'][6009])
   # peak_q = np.maximum(np.abs(of.q)*rmg.dx, peak_q)
   # old_q[:] = np.abs(of.q)
    #v = rmg['node']['surface_water__depth']
    #print("after infil:", max(of.h))
    hydrograph_time_sec.append(elapsed_time)
    discharge_at_outlet.append(of.q[link_to_sample])
    #gage_hydrographs.append((np.abs(of.q)*rmg.dx)[links])
    #other_hydrographs.append(rmg.at_node['surface_water__discharge'][DAs])
    #total_incision_depth += (np.abs(dle.I) * of.dt)

    i+=1
    print(elapsed_time, of.dt)
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

average_peak_discharge = []
maximum_peak_discharge =[]
rmg.add_field('link','peak_discharges', peak_q)
in_peaks = map_max_of_inlinks_to_node(rmg, 'peak_discharges')
out_peaks = map_max_of_outlinks_to_node(rmg, 'peak_discharges')
max_peaks = np.maximum(in_peaks, out_peaks)
average_incision_depth=[]
maximumum_incision_depth=[]
for each in unique_da:

    locations  = np.where(drainage_area == each)

    q_at_nodes = max_peaks[locations]
    I_at_nodes = total_incision_depth[locations]
    average_peak_discharge.append(np.average(q_at_nodes))
    maximum_peak_discharge.append(np.amax(q_at_nodes, axis=0))

    average_incision_depth.append(np.average(I_at_nodes))
    maximumum_incision_depth.append(np.amax(I_at_nodes, axis=0))


#
#if save_flag == 'True':
#    if storm_flag == 'FirstStorm':
#        print('saving to first storm...')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/Hydrographs.txt', other_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/GageOutletHydrographs.txt', gage_hydrographs)
#
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/StormTime_hr.txt', hydrograph_time_hrs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/StormTime_sec.txt', hydrograph_time_sec)
#
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/AvgPeakDischarge.txt', average_peak_discharge)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/MaxPeakDischarge.txt', maximum_peak_discharge)
#
#        write_esri_ascii('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/total_incision.asc', rmg, 'total_incision_depth')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/AvgIncisionDepth.txt', average_incision_depth)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/15mmhr/MaxIncisionDepth.txt', maximumum_incision_depth)
##    if storm_flag == 'SecondStorm':
##        print('saving to second storm...')
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/Hydrographs.txt', other_hydrographs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/GageOutletHydrographs.txt', gage_hydrographs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/StormTime_hr.txt', hydrograph_time_hrs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/StormTime_sec.txt', hydrograph_time_sec)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/AvgPeakDischarge.txt', average_peak_discharge)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/20mmhr/MaxPeakDischarge.txt', maximum_peak_discharge)
#    if storm_flag == 'ThirdStorm':
#        print('saving to third storm...')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/Hydrographs.txt', other_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/GageOutletHydrographs.txt', gage_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/StormTime_hr.txt', hydrograph_time_hrs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/StormTime_sec.txt', hydrograph_time_sec)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/AvgPeakDischarge.txt', average_peak_discharge)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/MaxPeakDischarge.txt', maximum_peak_discharge)
#        write_esri_ascii('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/total_incision.asc', rmg, 'total_incision_depth')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/AvgIncisionDepth.txt', average_incision_depth)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/25mmhr/MaxIncisionDepth.txt', maximumum_incision_depth)
##    if storm_flag == 'FourthStorm':
##        print('saving to fourth storm...')
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/Hydrographs.txt', other_hydrographs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/GageOutletHydrographs.txt', gage_hydrographs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/StormTime_hr.txt', hydrograph_time_hrs)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/StormTime_sec.txt', hydrograph_time_sec)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/AvgPeakDischarge.txt', average_peak_discharge)
##        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/30mmhr/MaxPeakDischarge.txt', maximum_peak_discharge)
#    if storm_flag == 'FifthStorm':
#        print('saving to fifth storm...')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/Hydrographs.txt', other_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/GageOutletHydrographs.txt', gage_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/StormTime_hr.txt', hydrograph_time_hrs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/StormTime_sec.txt', hydrograph_time_sec)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/AvgPeakDischarge.txt', average_peak_discharge)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/MaxPeakDischarge.txt', maximum_peak_discharge)
#        write_esri_ascii('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/total_incision.asc', rmg, 'total_incision_depth')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/AvgIncisionDepth.txt', average_incision_depth)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/35mmhr/MaxIncisionDepth.txt', maximumum_incision_depth)
#    if storm_flag == 'UngagedStorm':
#        print('saving to ungaged storm...')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/Hydrographs.txt', other_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/GageOutletHydrographs.txt', gage_hydrographs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/StormTime_hr.txt', hydrograph_time_hrs)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/StormTime_sec.txt', hydrograph_time_sec)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/AvgPeakDischarge.txt', average_peak_discharge)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/MaxPeakDischarge.txt', maximum_peak_discharge)
#        write_esri_ascii('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/total_incision.asc', rmg, 'total_incision_depth')
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/AvgIncisionDepth.txt', average_incision_depth)
#        np.savetxt('/Users/Jordan/Desktop/GSA 2016/Pitfilled/UngagedStorm/MaxIncisionDepth.txt', maximumum_incision_depth)
#else:
#    pass

endtime = time.time()
print('\n', 'Overall model run time: ', round(endtime - start_time, 2), ' seconds, for alpha = ', of.alpha)
storm_dur_min = storm_duration / 60.
hr = [x / 3600. for x in hydrograph_time_sec]
np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/Q25_InitMoisture1_60min.txt', np.abs(discharge_at_outlet)*rmg.dx)
np.savetxt('/Users/Jordan/Desktop/Chapter 2 - Buffalo Creek/InfiltrationTest/T25_InitMoisture1_60min.txt', hr)
#plt.figure(1)
#plt.plot(hydrograph_time_hrs, other_hydrographs)
#plt.xlim(0, 5)
#plt.ylim(0, 20)
#plt.title('Upstream Hydrographs')
#plt.savefig('/Users/Jordan/Desktop/Upstream_' + str(rainfall_mmhr) + '_' + str(storm_dur_min) +'min.png')
plt.figure(2)
plt.plot(hr, np.abs(discharge_at_outlet)*rmg.dx)
plt.title('Outlet Hydrograph')
plt.xlim(0, 8)
plt.ylim(0, 250)
#plt.savefig('/Users/Jordan/Desktop/Outlet_' + str(rainfall_mmhr) + '_' + str(storm_dur_min) +'min.png')
















    #dle.erode(of.dt, slope='water_surface__slope')

    #peak_incision_rate = np.maximum(np.abs(dle.I), old_incision_rate)

    #old_incision_rate[:]  = np.abs(dle.I)
    #rmg['node']['total_incision_depth'] += np.abs(dle.I) * of.dt

    ## Append time and discharge to their lists to save data and for plotting.
    #rmg['node']['topographic__elevation'] += uplift*of.dt





#filename_input0 = '/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/ASCII_Files/Storm_event_0_.asc'
#filename_input500 = '/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/ASCII_Files/Storm_event_500_.asc'

### Now the ASCII is read, assuming that it is standard ESRI format.
#
#rmg.at_node['total_incision'] = overall_I
#write_esri_ascii('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Total_Incision.asc', rmg, 'total_incision')
#
#total_incision_upstream = overall_I[DAs]
#np.savetxt('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Total_Incision_upstream.txt', total_incision_upstream)
#rmg.at_node['total_incision'] = overall_I
#write_esri_ascii('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Total_Incision.asc', rmg, 'total_incision')
#average_peak_incision = []
#maximum_peak_incision =[]

#average_incision = []
#maximum_incision =[]
#
#for each in unique_da:
#
#    locations  = np.where(drainage_area == each)
#
#    #q_at_nodes = max_peaks[locations]
#    I_at_nodes = overall_I[locations]
#
#
#    average_incision.append(np.average(I_at_nodes))
#    maximum_incision.append(np.amax(I_at_nodes, axis=0))
#
#total_incision_upstream = overall_I[DAs]
#np.savetxt('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Total_Incision_upstream.txt', total_incision_upstream)
#
#np.savetxt('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Average_Incision_DA.txt', average_incision)
#np.savetxt('/Users/Jordan/Desktop/Paper1_Runs/Square/StandardStorm/TotalIncision/Maximum_Incision_DA.txt',maximum_incision)