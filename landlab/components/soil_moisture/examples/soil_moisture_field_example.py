#########################################################
##  
##    Example for soil_moisture_field.py
##
<<<<<<< HEAD
##    Sai Nudurupati   15May2014
=======
##    Sai Nudurupati and Erkan Istanbulluoglu - 15May2014
>>>>>>> FETCH_HEAD
##
#########################################################
from landlab import RasterModelGrid
from landlab.components.soil_moisture.soil_moisture_field import SoilMoisture
<<<<<<< HEAD
from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
import numpy as np
import matplotlib.pyplot as plt
=======
#from landlab.components.uniform_precip.generate_uniform_precip import PrecipitationDistribution
import numpy as np
#import matplotlib.pyplot as plt
>>>>>>> FETCH_HEAD

grid = RasterModelGrid( 10, 10, 1. )
grid.add_zeros('cell','VegetationCover',units = 'None')
grid.add_zeros('cell','LiveLeafAreaIndex',units = 'None')
grid.add_zeros('cell','PotentialEvapotranspiration',units = 'mm')
grid['cell']['InitialSaturationFraction'] = np.random.rand(grid.number_of_cells)
grid['cell']['VegetationsCover'] = np.random.rand(grid.number_of_cells)
grid['cell']['LiveLeafAreaIndex'] = np.ones(grid.number_of_cells)*3
grid['cell']['PotentialEvapotranspiration'] = np.ones(grid.number_of_cells)*6
current_time = 0.56
SM = SoilMoisture( grid )
#PD = PrecipitationDistribution()
#PD.update()
#precip = PD.get_storm_depth()
#Tb = PD.get_interstorm_event_duration()
#Tr = PD.get_precipitation_event_duration()
current_time = SM.update( current_time )