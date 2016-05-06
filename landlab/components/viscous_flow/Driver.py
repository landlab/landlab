# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 13:11:16 2016

@author: katherinekravitz
"""


from viscous_flow import ViscousFlowModel
from landlab import RasterModelGrid
from landlab.io.netcdf import write_netcdf

rg = RasterModelGrid(5, 5, 25.0)
rg.set_closed_boundaries_at_grid_edges(True, True, True, True)

vfm = ViscousFlowModel(rg, viscosity=1.0e10, rho_overburden=2600.0, 
                       rho_salt=2200.0, slope_salt=0.0, g=9.8)
             
dt = 1.0
run_duration = 10.0
output_interval = 2.0
output_filename = 'test_salt_output'

nt = int(run_duration / dt)
next_output = output_interval
output_slice = 0


vfm._grid.at_node['topographic__elevation'][12] = 10.0
# vfm._grid.at_node['topographic__elevation'][540] = 1050.0


def write_output(filename, time_slice, grid):
    outname = filename + str(time_slice).zfill(4) + '.nc'
    write_netcdf(outname, grid)
    

write_output(output_filename, output_slice, rg)
output_slice += 1

for i in range(nt):
    
    vfm.flow(dt)
    
    if i*dt >= next_output:
        write_output(output_filename, output_slice, rg)
        output_slice += 1
        next_output += output_interval
        
    
    

