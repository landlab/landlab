# Sai Nudurupati and Erkan Istanbulluoglu - 14May2014 :
# Example to use radiation_field.py
"""
    Article author: Sai S. Nudurupati (saisiddu@uw.edu) and Erkan Istanbulluoglu
    Date: 22 May 2014

"""
"""
    Radiation component calculates total incoming shortwave radiation and
    relative radiation factor (ratio of total radiation incident on the surface
    with respect to the total radiation incident on flat surface).
    This example demonstrates a simple use case for radiation_field.py.
    In this example, a random elevation field of 100m X 100m is created
    with a cell area of 1 m^2 each. A raster grid with this elevation field
    is created. Total incident short wave radiation and radiation factor
    on a given Julian day (at noon) is calculated and plotted.
"""

"""
    import landlab's raster grid library 'RasterModelGrid'
"""
from landlab import RasterModelGrid
"""
    import 'Radiation' class from 'components' library under radiation package.
"""
from landlab.components.radiation.radiation_field import Radiation
"""
    import landlab's plotting function that has the capability of plotting
    'fields' stored on the grid.
"""
from landlab.plot.imshow import imshow_grid
"""
    import Numpy library. We will use 'random' module from numpy for this
    tutorial.
"""
import numpy as np
"""
    import pyplot module from matplotlib. This is a handy plotting library.
    For easier access (and typing), call 'matplotlib.pyplot' as 'plt'
"""
import matplotlib.pyplot as plt

"""
    RasterModelGrid module creates a raster grid of size defined by its first
    two arguments with a spatial resolution defined by its third optional
    argument. Hence the grid has 100 X 100 nodes (a total of 10000 nodes) with
    spacing between two nodes being 1.0 unit.
"""
grid = RasterModelGrid(100,100, 1.)
"""
    Create a random elevation field. 'np.random.rand' function returns an
    array of random numbers that follow a uniform distribution with range(0,1)
    and with size defined by its' argument.
"""
elevation = np.random.rand(grid.number_of_nodes) * 1000
"""
    Creat a nodal field called 'Elevation' on the grid with units in metres and
    populate it with zeros.
"""
grid.add_zeros('node','Elevation',units = 'm')
"""
    This 'Elevation' field stored on the grid can be accessed as following:
"""
grid['node']['Elevation'] = elevation
"""
    Instantiate an object for 'Radiation' Class. This instantiation associates
    the object 'rad' with the capabilities of the class 'Radiation'. This
    initiation requires an input of a grid. Creation of the object
    automatically associates this grid to the object 'rad'.
"""
rad = Radiation( grid )
"""
   Set random time for our example. Time is in years.
"""
current_time = 0.56
"""
    'Radiation' class has an update function (like CSDMS BMI component). This
    function takes current_time as an input argument and calculates Total Short
    Wave Radiation incident on each cell at noon of the julian day represented
    by current_time input. It also calculates the Radiation Factor which
    represents the ratio of total short wave radiation incident on a grid cell
    and flat surface. Hence Radiation Factor will be 1.0 if the grid cell has a
    slope of 0.0 . Therefore, two cellular fields are created on the grid
    (which means that the two arrays of length equivalent to number of cells
    that the grid has are created and stored in conjunction with the grid.
    Whenever this grid is transferred, these two cellular fields go with them.
"""
rad.update( current_time )
"""
    Create a figure window available from pyplot library. This allows separating
    figures while plotting multiple figures.
"""
plt.figure(0)
"""
    Plot the cellular field 'TotalShortWaveRadiation' that is available on the
    grid. imshow_grid is a Landlab plotting tool that reads the input of
    grid, the variable name that needs to be plotted, and type of field (whether
    it is 'cell' or 'node', etc... It also reads optional inputs (keywords),
    grid_units (units of grid X and Y axes , e.g. 'm'). For more options, please refer
    documentation for landlab.plot.imshow.
"""
imshow_grid(grid,'TotalShortWaveRadiation', values_at = 'cell',
            grid_units = ('m','m'))

"""
    The plot created can be saved using the function 'savefig' available in
    pyplot library. This file will be saved in the current directory that your
    python shell is in. You can know the current directory by using 'pwd'
    command in the shell. (This documentation works if Enthought Canopy is used.
    It might work for other packages too but is not tested by the author.)
"""
plt.savefig('Radiation')
"""
    Plot another figure.
"""
plt.figure(1)
"""
    Using Landlab plotting tool, lets plot another variable 'RadiationFactor'.

"""
imshow_grid(grid,'RadiationFactor', values_at = 'cell', grid_units = ('m','m'))
"""
    Save this figure.
"""
plt.savefig('RadiationFactor')
"""
    Figure windows generated by pyplot library do not pop up by default. We have
    to use function 'show' availabe in pyplot library to allow these figures to
    pop up.
"""
plt.show()
"""
    Please note that this is a simple implementation of radiation_field
    component intended to familiarize its use. Please refer to the documentation
    of landlab/components/radiation/radiation_field for more information to
    use this component.
"""
