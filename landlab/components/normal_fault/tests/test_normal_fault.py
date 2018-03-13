from landlab import RasterModelGrid
from landlab.components import NormalFault, FastscapeEroder, FlowAccumulator
grid = RasterModelGrid((6, 6), spacing=10)


z = grid.add_zeros('node', 'topographic__elevation')

param_dict = {'faulted_surface': 'topographic__elevation',
              'active_uplift_throw_rate' : 0.5,
              'fault_dip_angle': 90.0,
              'active_uplift_start_time': 10, 
              'active_uplift_end_time': 20,
              'fault_trace_dict': {'y1': 0,
                                   'x1': 0, 
                                   'y2': 30, 
                                   'x2': 60},
             'include_boundaries': True}

nf = NormalFault(grid, params=param_dict)
fr = FlowAccumulator(grid)
fs = FastscapeEroder(grid, K_sp=0.01)

print(nf.faulted_nodes.reshape(grid.shape))

print(z.reshape(grid.shape))

dt = 1.0
for i in range(30):
    nf.run_one_step(dt)
    fr.run_one_step()
    fs.run_one_step(dt)

print(z.reshape(grid.shape))

#%%

from landlab import RasterModelGrid
from landlab.components import NormalFault, FastscapeEroder, FlowAccumulator
grid = RasterModelGrid((30, 30), spacing=10)

from landlab.plot import imshow_grid

z = grid.add_zeros('node', 'topographic__elevation')

param_dict = {'faulted_surface': 'topographic__elevation',
              'active_uplift_throw_rate' : 0.001,
              'fault_dip_angle': 90.0,
              'active_uplift_start_time': 10, 
              'fault_trace_dict': {'y1': 0,
                                   'x1': 0, 
                                   'y2': 30, 
                                   'x2': 60},
             'include_boundaries': False}

nf = NormalFault(grid, params=param_dict)
fr = FlowAccumulator(grid)
fs = FastscapeEroder(grid, K_sp=0.01)


dt = 100.0
for i in range(300):
    nf.run_one_step(dt)
    fr.run_one_step()
    fs.run_one_step(dt)
    
imshow_grid(grid, 'topographic__elevation')
