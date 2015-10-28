# Sai Nudurupati and Erkan Istanbulluoglu- 16May2014 :
# Example to use potential_evapotranspiration_field.py

#import landlab
from landlab import RasterModelGrid
from landlab.components.radiation.radiation_field import Radiation
from landlab.components.pet.potential_evapotranspiration_field import (
    PotentialEvapotranspiration)
import numpy as np
import matplotlib.pyplot as plt
from landlab.plot.imshow import imshow_grid

grid = RasterModelGrid( 100, 100, 20. )
elevation = np.random.rand(grid.number_of_nodes) * 1000
grid.add_zeros('node','Elevation',units = 'm')
grid['node']['Elevation'] = elevation
rad = Radiation( grid )
PET = PotentialEvapotranspiration( grid )
current_time = 0.56
rad.update( current_time )
PET.update( ConstantPotentialEvapotranspiration = 10.0 )

plt.figure(0)
imshow_grid(grid,'RadiationFactor', values_at = 'cell',
            grid_units = ('m','m'))

plt.figure(1)
imshow_grid(grid,'PotentialEvapotranspiration', values_at = 'cell',
            grid_units = ('m','m'))
plt.savefig('PET_test')
plt.show()
