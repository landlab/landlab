# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 10:20:19 2018

Playing with Margaux's DataRecord Tutorial


@author: pfeif
"""
import numpy as np
from landlab import RasterModelGrid
from landlab.data_record import DataRecord

from landlab import imshow_grid
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, subplot, xlabel, ylabel, title, legend, figure


# %%
grid_3 = RasterModelGrid((5, 5), (2, 2))
initial_boulder_sizes_3 = np.array([[10], [4], [8], [3], [5]])

boulders_3 = {'grid_element' : 'node',
            'element_id' : np.array([[6], [11], [12], [17], [12]])}

dr_3 = DataRecord(grid_3,
                  time=[0.],
                  items=boulders_3,
                  data_vars={'boulder_size' : (['item_id', 'time'], initial_boulder_sizes_3)}, 
                  attrs={'boulder_size' : 'm'})

for t in range(5):

    dr_3.add_record(time=np.array([t]))

    dr_3.ffill_grid_element_and_id() 

f = dr_3['time']== [0.0]

f.to_dataframe()
np.shape(f)

s = dr_3.calc_aggregate_value(func=np.sum,
                              data_variable='boulder_size', 
                              at='node', 
                              filter_array=f)



# %%

grid_3 = RasterModelGrid((5, 5), (2, 2))

initial_boulder_sizes_3 = np.array([[10], [4], [8], [3], [5]])
boulder_lithologies = np.array(['sandstone', 'granite', 'sandstone', 'sandstone', 'limestone']) #same as above, already run

boulders_3 = {'grid_element' : 'node',
            'element_id' : np.array([[6], [11], [12], [17], [12]])}

dr_3 = DataRecord(grid_3,
                  time=[0.],
                  items=boulders_3,
                  data_vars={'boulder_size' : (['item_id', 'time'], initial_boulder_sizes_3),
                             'boulder_litho': (['item_id'], boulder_lithologies)}, 
                  attrs={'boulder_size' : 'm'})


dt = 100
total_time = 1000

time_index = 1

for t in range(dt, total_time, dt):

    # create a new time coordinate:
    dr_3.add_record(time=np.array([t]))

    # this propagates grid_element and element_id values forward in time (instead of the 'nan' default filling):
    dr_3.ffill_grid_element_and_id() 

#    for i in range(0, dr_3.number_of_items):
#        # value of block erodibility:
#        if dr_3['boulder_litho'].values[i] == 'limestone':
#            k_b = 10**-5
#        elif dr_3['boulder_litho'].values[i] == 'sandstone':
#            k_b = 3*10**-6
#        elif dr_3['boulder_litho'].values[i] == 'granite':
#            k_b = 3*10**-7
#        else:
#            print('Unknown boulder lithology')
#
#        dr_3['boulder_size'].values[i, time_index] = dr_3['boulder_size'].values[i, time_index-1] - k_b * dr_3['boulder_size'].values[i, time_index-1] * dt
#
#    time_index += 1

# %%
grid = RasterModelGrid((5, 5), (2, 2))
initial_thing_size = np.array([[10], [4], [8], [3], [5]])

thing = {'grid_element' : 'node',
            'element_id' : np.array([[6], [11], [12], [17], [12]])}

dr_3 = DataRecord(grid_3,
                  time=[0.],
                  items=thing,
                  data_vars={'thing_size' : (['item_id', 'time'], initial_thing_size)}, 
                  attrs={'thing_size' : 'm'})

for t in range(5):
    dr_3.add_record(time=np.array([t]))
    dr_3.ffill_grid_element_and_id() 

# OPTION 1 - yields ValueError: The truth value of an array with more than one element is ambiguous. Use a.any() or a.all()
f = dr_3['time']== [0.0] and dr_3['item_id']<= 4

# OPTION 2 - works
f = dr_3['thing_size']<= 4
f[:,0:-1]= False # workaround to only look at current timestep
        
        

f.to_dataframe()
print('Shape of f = ', np.shape(f))

s = dr_3.calc_aggregate_value(func=np.sum,
                              data_variable='thing_size', 
                              at='node', 
                              filter_array=f)


