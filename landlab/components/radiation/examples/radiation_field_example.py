# Sai Nudurupati and Erkan Istanbulluoglu - 14May2014 : 
# Example to use radiation_field.py

import landlab
from landlab import RasterModelGrid
from landlab.components.radiation.radiation_field import Radiation
import numpy as np
import matplotlib.pyplot as plt
from landlab.plot.imshow import imshow_field

grid = RasterModelGrid( 100, 100, 20. )
elevation = np.random.rand(grid.number_of_nodes) * 1000
grid.add_zeros('node','Elevation',units = 'm')
grid['node']['Elevation'] = elevation
rad = Radiation( grid )
current_time = 0.56
rad.update( current_time )

plt.figure(0)
imshow_field(grid,'TotalShortWaveRadiation',
                values_at = 'cell', grid_units = ('m','m'))
plt.show()