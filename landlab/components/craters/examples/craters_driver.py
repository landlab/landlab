import landlab.components.craters as craters
from landlab import RasterModelGrid
import numpy
from pylab import show, imshow, colorbar, plot
from copy import copy

craters=reload(craters)

#craters_data = craters.one_crater_then_degrade()

class data(object):
    '''
        This is where all the whole-grid data lives, as arrays over the various elements of the grid.
        '''
    #Data goes here!!!
    def __init__(self, grid):
        self.elev = grid.create_node_dvector() #some data
        self.flag_already_in_the_list = grid.create_node_dvector()
        self.craters_over_max_radius_not_plotted = grid.create_node_dvector()
        self.impact_sequence = []

#def dig_one_crater_driver(nr, nc, dx, rel_x, rel_y, radius=numpy.nan, angle=numpy.nan):
#    '''
#        This is an ad-hoc script to dig one crater.
#        '''
#    
#    #Setup the grid & impactor object
#    mg = RasterModelGrid()
#    mg.initialize(nr, nc, dx)
#    vectors = data(mg)
#    vectors.elev[:] = 0.
#    cr = craters.impactor()
#
#    #Dig the crater
#    cr.excavate_a_crater_optimized(mg, vectors, forced_radius=radius, forced_angle=angle, forced_pos=[rel_x, rel_y])
#
#    #Finalize
#    elev_raster = mg.node_vector_to_raster(vectors.elev, flip_vertically=True)
#    return elev_raster


#elev_raster = dig_one_crater_driver(120, 120, 0.025, 0.5, 0.5, 0.75, 90.)
#imshow(elev_raster)
#colorbar()
#show()

##This code snippet to produce a landscape w. 500000 impacts, but where the first half of the impacts are 10x bigger on average.
#mg = RasterModelGrid()
#mg.initialize(1200, 1200, 0.0025)
#vectors = data(mg)
#vectors.elev[:] = 0.
#
#dict_of_vectors_50m_min, dict_of_vectors_5m_min = craters.ten_times_reduction(mg, vectors)
#
##This code snippet produces a smaller (800x800) grid where we let the basic crater routine run for a long time (1M impacts). The aim is to see if we get a steady state landscape.

#mg2 = RasterModelGrid()
#dict_of_vectors2_5m_min = {}
#mg2.initialize(600, 600, 0.0025)
#vectors2 = data(mg2)
#vectors2.elev[:] = 0.

#for i in xrange(20, 50):
#    print 'LOOP NUMBER ', i
#    mg2, vectors2, profile, xsec = craters.dig_some_craters(mg2, vectors2, nt_in=10000)
#    dict_of_vectors2_5m_min[i] = copy(vectors2)

#This code snippet to produce a landscape w. 200000 impacts, but where the first quarter of the impacts are 3x bigger on average. The output interval is only 5k impacts this time.
mg3 = RasterModelGrid()
mg3.initialize(1200, 1200, 0.0025)
vectors = data(mg3)
vectors.elev[:] = 0.

dict1_1200_50k_15m_min, dict1_1200_150k_5m_min = craters.step_reduce_size(mg3, vectors, loops=[10,30], interval=5000, min_radius_in=[0.015,0.005])

vectors.elev[:] = 0.

dict2_1200_150k_10m_min, dict2_1200_150k_5m_min = craters.step_reduce_size(mg3, vectors, loops=[30,30], interval=5000, min_radius_in=[0.01,0.005])


