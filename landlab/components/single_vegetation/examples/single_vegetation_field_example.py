#########################################################
##  
##    Example for single_vegetation_field.py
##
##    Sai Nudurupati and Erkan Istanbulluoglu - 16May2014
##
#########################################################
from landlab import RasterModelGrid
from landlab.components.soil_moisture.soil_moisture_field import SoilMoisture
from landlab.components.single_vegetation.single_vegetation_field import Vegetation
import numpy as np
from landlab.plot.imshow
#import matplotlib.pyplot as plt

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
VEG = Vegetation( grid )
current_time = SM.update( current_time )
VEG.update()
