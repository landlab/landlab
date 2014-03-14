from landlab.components.flow_routing.route_flow_dn import FlowRouter
from landlab.components.stream_power.fastscape_stream_power import SPEroder
from landlab.components.nonlinear_diffusion.Perron_nl_diffuse import PerronNLDiffuse
from landlab.components.diffusion.diffusion import DiffusionComponent
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
from landlab import ModelParameterDictionary

from landlab import RasterModelGrid
import numpy as np
import pylab

#get the needed properties to build the grid:
input_file = './germany_test_params.txt'
inputs = ModelParameterDictionary(input_file)
nrows = inputs.read_int('nrows')
ncols = inputs.read_int('ncols')
dx = inputs.read_float('dx')
leftmost_elev = inputs.read_float('leftmost_elevation')
initial_slope = inputs.read_float('initial_slope')
uplift_rate = inputs.read_float('uplift_rate')
runtime = inputs.read_float('total_time')
try:
    dt = inputs.read_float('dt')
    nt = int(runtime//dt)
    uplift_per_step = uplift_rate * dt
except:
    dynamic_timesteps=True
    pd = PrecipitationDistribution()
    pd.initialize(input_file)
else:
    dynamic_timesteps = False

mg = RasterModelGrid(nrows, ncols, dx)
mg.set_inactive_boundaries(False, True, False, True)

#create the elevation field in the grid
mg.create_node_array_zeros('planet_surface__elevation')
z = mg.create_node_array_zeros() + leftmost_elev
z += initial_slope*np.amax(mg.node_y) - initial_slope*mg.node_y
mg['node'][ 'planet_surface__elevation'] = z + np.random.rand(len(z))/100000.

# Display a message
print 'Running ...' 

#instantiate the components:
fr = FlowRouter(mg)
sp = SPEroder(mg, input_file)
diffuse = PerronNLDiffuse(mg, input_file)
lin_diffuse = DiffusionComponent(grid=mg)
lin_diffuse.initialize(input_file)

#perform the loops:
if not dynamic_timesteps:
    for i in xrange(nt):
        mg['node']['planet_surface__elevation'][mg.get_interior_nodes()] += uplift_per_step
        mg = fr.route_flow(grid=mg)
        mg = sp.erode(mg)
        mg = diffuse.diffuse(mg, i*dt)
        #mg = lin_diffuse.diffuse(mg, dt)
        print 'Completed loop ', i
else:
    elapsed_time = 0.
    breaker = True
    track_storm_durations = []
    track_interstorm_durations = []
    track_storm_intensities = []
    while breaker:
        #A storm occurs
        dt = pd.storm_duration
        if elapsed_time+dt > runtime:
            dt = runtime-elapsed_time
            breaker = False
        uplift_per_step = uplift_rate*dt
        mg['node']['planet_surface__elevation'] += uplift_per_step
        mg = fr.route_flow(grid=mg)
        sp.gear_timestep(dt, pd.intensity)
        mg = sp.erode(mg)
        diffuse.gear_timestep(dt)
        mg = diffuse.diffuse(mg, elapsed_time)
        track_storm_durations.append(dt)
        track_storm_intensities.append(pd.intensity)
        elapsed_time += dt
        print 'storm duration ', pd.storm_duration
        print 'storm intensity ', pd.intensity
        #an interval between storms occurs - diffusion is permitted still
        dt = pd.interstorm_duration
        if elapsed_time+dt > runtime:
            dt = runtime-elapsed_time
            breaker = False
        if elapsed_time > runtime:
            track_interstorm_durations.append(0.)
            break
        uplift_per_step = uplift_rate*dt
        mg['node']['planet_surface__elevation'] += uplift_per_step
        diffuse.gear_timestep(dt)
        mg = diffuse.diffuse(mg, elapsed_time)
        elapsed_time += dt
        print 'interstorm duration ', pd.interstorm_duration
        print 'total time ', elapsed_time
        track_interstorm_durations.append(dt)
        pd.update()

print np.amin(mg['node']['drainage_area'])

#Finalize and plot
#elev = mg['node']['planet_surface__elevation']
elev = fr.node_water_discharge
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
pylab.figure(1)
pylab.close()
pylab.figure(1)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')

elev = mg['node']['planet_surface__elevation']
#elev = fr.node_water_discharge
elev_r = mg.node_vector_to_raster(elev)
# Clear previous plots
pylab.figure(2)
im = pylab.imshow(elev_r, cmap=pylab.cm.RdBu)  # display a colored image
pylab.colorbar(im)
pylab.title('Topography')

drainage_areas = mg['node']['drainage_area'][mg.get_interior_nodes()]
steepest_slopes = mg['node']['steepest_slope'][mg.get_interior_nodes()]
pylab.figure(3)
pylab.loglog(drainage_areas, steepest_slopes, 'x')
pylab.xlabel('Upstream drainage area, m^2')
pylab.ylabel('Maximum slope')

print('Done.')

pylab.show()