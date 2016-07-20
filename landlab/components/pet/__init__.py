"""Landlab component that simulates potential evapotranspiration rate.

Potential Evapotranspiration Component calculates spatially distributed
potential evapotranspiration based on input radiation factor (spatial
distribution of incoming radiation) using chosen method such as constant
or Priestley Taylor. Ref: ASCE-EWRI Task Committee Report Jan 2005.

.. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.pet import PotentialEvapotranspiration

Create a grid on which to calculate potential evapotranspiration rate.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> PotentialEvapotranspiration.input_var_names
('radiation__ratio_to_flat_surface',)

Check the units for the fields.

>>> PotentialEvapotranspiration.var_units('radiation__ratio_to_flat_surface')
'None'

Create the input fields.

>>> grid['cell']['radiation__ratio_to_flat_surface'] = np.array([
...       0.38488566, 0.38488566,
...       0.33309785, 0.33309785,
...       0.37381705, 0.37381705])

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> PotentialEvapotranspiration.var_help('radiation__ratio_to_flat_surface')
name: radiation__ratio_to_flat_surface
description:
  ratio of total incident shortwave radiation on sloped surface
  to flat surface
units: None
at: cell
intent: in

Check the output variable names

>>> sorted(PotentialEvapotranspiration.output_var_names)
['radiation__incoming_shortwave_flux',
 'radiation__net_flux',
 'radiation__net_longwave_flux',
 'radiation__net_shortwave_flux',
 'surface__potential_evapotranspiration_rate']

Instantiate the 'PotentialEvapotranspiration' component to work on this grid,
and run it.

>>> PET = PotentialEvapotranspiration(grid, method='PriestleyTaylor')

Run the *update* method to update output variables with current time
>>> current_time = 0.5
>>> PET.update(current_time)

>>> PET.grid.at_cell['radiation__incoming_shortwave_flux']
array([ 33.09968448,  33.09968448,  28.64599771,  28.64599771,
        32.14779789,  32.14779789])

>>> PET.grid.at_cell['radiation__net_flux']
array([ 13.9764353 ,  13.9764353 ,  12.09585347,  12.09585347,
        13.57449849,  13.57449849])

>>> PET.grid.at_cell['radiation__net_shortwave_flux']
array([ 13.23987379,  13.23987379,  11.45839908,  11.45839908,
        12.85911915,  12.85911915])

>>> PET.grid.at_cell['surface__potential_evapotranspiration_rate']
array([ 0.25488065,  0.25488065,  0.22058551,  0.22058551,  0.24755075,
        0.24755075])
"""
from .potential_evapotranspiration_field import PotentialEvapotranspiration


__all__ = ['PotentialEvapotranspiration', ]
