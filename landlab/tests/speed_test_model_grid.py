"""
simple script to run speed tests of various functions in model grid
"""

from landlab import RasterModelGrid
import time
import numpy

mg = RasterModelGrid(20, 30, 1.0)

nt = 1000

s = mg.create_node_dvector()
g = numpy.zeros(mg.num_active_links)

start_time = time.time()

for i in range(nt):
    
    g = mg.calculate_gradients_at_active_links(s, g)
    
time1 = time.time()

for i in range(nt):
    
    g = mg.calculate_gradients_at_active_links_slow(s, g)
    
time2 = time.time()

for i in range(nt):
    
    divg = mg.calculate_flux_divergence_at_nodes(g)
    
time3 = time.time()

for i in range(nt):
    
    divg = mg.calculate_flux_divergence_at_nodes_slow(g)

time4 = time.time()
  
print('Elapsed time with fast gradient algo: '+str(time1-start_time))
print('Elapsed time with slow gradient algo: '+str(time2-time1))
print('Elapsed time with fast node-divergence algo: '+str(time3-time2))
print('Elapsed time with slow node-divergence algo: '+str(time4-time3))

   