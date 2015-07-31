
# Brigid Lynch
# Import required libraries
from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.stream_power import StreamPowerEroder
from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab.components.diffusion.diffusion import LinearDiffuser #the two different diffusion formulations
import os, random #to manipulate path names
import time #timekeeping
import numpy as np #numerical operations
import pylab #to plot results
from landlab import RasterModelGrid #import grid class
from landlab import ModelParameterDictionary #handles input from input file
from landlab.plot.imshow import imshow_node_grid
from scipy.io import netcdf
from landlab.io.netcdf import write_netcdf



os.chdir('/Users/blynch/Research/Working_Directories/landlab/driver_bml')

#get needed properties to build grid:
input_file = 'driver6_params.txt'


#initialize an object that will supply the parameters:
inputs = ModelParameterDictionary(input_file) 
x = inputs.read_int('x') 
y = inputs.read_int('y') 
dx = inputs.read_float('dx') 
dt = inputs.read_float('dt') 
run_time = inputs.read_float('run_time') 
uplift = inputs.read_float('uplift_rate') 
#init_elev = inputs.read_float('init_elev') 
initial_slope = inputs.read_float('initial_slope') 
nt = int(run_time//dt) #divide and truncate (number of timesteps) diffusion_driver

#Instatiate the grid object, 
mg = RasterModelGrid(x, y, dx) #sp

#Create fields
os.chdir('/Users/blynch/Research/Working_Directories/wrf_hydro/Andes') #Location of topo file
t=netcdf.netcdf_file('Fulldom_hires_netcdf_file.nc')
topography = t.variables['TOPOGRAPHY']
t = []
for i in range(0,1990):
    for j in range(0,1990):
        t.append(topography[i][j])
topographic__elevation = np.array(t,dtype=float) #making input file into an array so it can be read by landlab

topography = mg.add_field('node', 'topographic__elevation', topographic__elevation) #create the field


os.chdir('/Users/blynch/Research/Working_Directories/landlab/driver_bml')


#instantiate the components:
fr = FlowRouter(mg)
sp = StreamPowerEroder(mg, input_file)
diffuse = PerronNLDiffuse(mg, input_file)
lin_diffuse = LinearDiffuser(grid=mg, input_stream=input_file)


#set boundary conditions
mg.set_fixed_value_boundaries_at_grid_edges(True, True, True, True) 



print( 'Running ...' )
time_on = time.time()



#perform the loop

uplifted_nodes = mg.get_core_nodes() #perform block uplift of interior of grid, leave boundary nodes at original elevations
    
#stream power loop

elapsed_time = 0. #total time in simulation
while elapsed_time < run_time:
    print elapsed_time
    if elapsed_time+dt>run_time:
        print "Short step!"
        dt = run_time - elapsed_time
    
    os.chdir('/Users/blynch/Research/Working_Directories/wrf_hydro/Andes/discharge_single') #Location of discharge maps
    #randnum=random.randint(0,40)
    #discharge_filename=('streamflow' + str(randnum) +'.nc') #match discharge_filename to file names in discharge map directory
    #print (discharge_filename) 

    f=netcdf.netcdf_file('streamflow2.nc')
    discharge = f.variables['flow']
    f = []
    for i in range(0,1990):
        for j in range(0,1990):
            f.append(discharge[j][i])
    water__volume_flux = np.array(f,dtype=float)  #making input file into an array so it can be read by landlab

    water__volume_flux = mg.add_field('node', 'water__volume_flux', water__volume_flux) #create the field 
    

    os.chdir('/Users/blynch/Research/Working_Directories/landlab')
    mg = fr.route_flow()     
    #erode
    mg,_,_ = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes= 'topographic__steepest_slope', Q_if_used='water__volume_flux')
   
    #add uplift
    mg.at_node['topographic__elevation'][mg.core_nodes] +=uplift*dt
    elapsed_time += dt
    
    #diffusion driver loop
    for i in xrange(nt): # loop nt times
    #this line performs functionality of component,swap between these two
        mg = lin_diffuse.diffuse(dt) #linear diffusion
        #mg = diffuse.diffuse(mg, i*dt) #nonlinear diffusion
        #mg = fr.route_flow()
        #mg = sp.erode(mg, dt, node_drainage_areas='drainage_area', slopes_at_nodes= 'topographic__steepest_slope', Q_if_used='water__volume_flux')
        
    ##plot N-S cross section from this stage in the run onto figure 1. 
    #pylab.figure(2)
    #elev_r = mg.node_vector_to_raster(mg['node']['topographic__elevation']) #turn the 1-D array of elevation values 
    #    #into a spatially accurate 2-D gridded format, for plotting
    #drainage_area = mg.node_vector_to_raster(mg['node']['drainage_area'])
    #im = pylab.plot(mg.dx*np.arange(x), elev_r[:,int(y//2)])
    #im = pylab.plot(mg.dx*np.arange(x), drainage_area[:,int(y//2)])
    #    #square brackets: subset of nodes to use
    #    #this type of data extraction from a larger data structure is called slicing/ fancy indexing 
    print mg.at_node['water__volume_flux_in']
    print mg.at_node['water__volume_flux']
        
time_off = time.time()
print 'Elapsed time: ', time_off-time_on

#Finalize and plot
#elev = mg['node']['topographic__elevation']
#evel_r = mg.node_vector_to_raster(elev) 
#drain_a = mg['node']['drainage_area']
#drainage_area = mg.node_vector_to_raster(drain_a)

#Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(2)
pylab.close()
pylab.figure(3)
pylab.close()
pylab.figure(4)
pylab.close()

#Plot topography
pylab.figure(1)
im = imshow_node_grid(mg, 'topographic__elevation', cmap='terrain') #display a colored image

#pylab.figure(2)
#im = pylab.plot(dx*np.arange(x), elev_r[:,int(y//2)]) #display a colored image
#pylab.title('Vertical Cross Section')

pylab.figure(3)
im = imshow_node_grid(mg, 'drainage_area')
pylab.show()

pylab.figure(4)
im = imshow_node_grid(mg, 'water__volume_flux', cmap='PuBu')
pylab.show()

pylab.figure(5)
im = imshow_node_grid(mg, 'water__volume_flux_in')
pylab.show()

os.chdir('/Users/blynch/Research/Working_Directories/landlab/driver_bml')

write_netcdf('topography_out.nc', mg, names='topographic__elevation') 
write_netcdf('discharge_out.nc', mg, names='water__volume_flux')
write_netcdf('drainage_area_out.nc', mg, names='drainage_area')





print('Done. ')
