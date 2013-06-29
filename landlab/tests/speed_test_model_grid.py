"""
simple script to run speed tests of various functions in model grid
"""

from landlab import RasterModelGrid
import time
import numpy

mg = RasterModelGrid(20, 30, 1.0)

nt = 3000

s = mg.create_node_dvector()
g = numpy.zeros(mg.num_active_links)

start_time = time.time()

for i in range(nt):
    
    g = mg.calculate_gradients_at_active_links(s, g)
    
print('Elapsed time with fast algo: '+str(time.time()-start_time))

start_time = time.time()

for i in range(nt):
    
    g = mg.calculate_gradients_at_active_links_slow(s, g)
    
print('Elapsed time with slow algo: '+str(time.time()-start_time))

    