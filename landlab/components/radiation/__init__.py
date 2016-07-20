"""Landlab component that simulates relative incidence shortwave radiation
on sloped surface.

Landlab component that computes 1D and 2D total incident shortwave
radiation. This code also computes relative incidence shortwave radiation
compared to a flat surface. Ref: Bras, Rafael L. Hydrology: an introduction
to hydrologic science. Addison Wesley Publishing Company, 1990.

.. codeauthor:: Sai Nudurupati & Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import Radiation

Create a grid on which to calculate incident shortwave radiation

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> Radiation.input_var_names
('topographic__elevation',)

Check the units for the fields.

>>> Radiation.var_units('topographic__elevation')
'm'

Create the input fields.

>>> grid['node']['topographic__elevation'] = np.array([
...       0., 0., 0., 0.,
...       1., 1., 1., 1.,
...       2., 2., 2., 2.,
...       3., 4., 4., 3.,
...       4., 4., 4., 4.])

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Radiation.var_help('topographic__elevation')
name: topographic__elevation
description:
  elevation of the ground surface relative to some datum
units: m
at: node
intent: in

Check the output variable names

>>> sorted(Radiation.output_var_names) # doctest: +NORMALIZE_WHITESPACE
['radiation__incoming_shortwave_flux',
 'radiation__net_shortwave_flux',
 'radiation__ratio_to_flat_surface']

Instantiate the 'Radiation' component to work on this grid, and run it.

>>> rad = Radiation(grid)

Run the *update* method to update output variables with current time
>>> current_time = 0.5
>>> rad.update(current_time)

>>> rad.grid.at_cell['radiation__ratio_to_flat_surface']
array([ 0.38488566,  0.38488566,  0.33309785,  0.33309785,  0.37381705,
        0.37381705])

>>> rad.grid.at_cell['radiation__incoming_shortwave_flux']
array([ 398.33664988,  398.33664988,  344.73895668,  344.73895668,
        386.88120966,  386.88120966])
"""
from .radiation import Radiation


__all__ = ['Radiation', ]
