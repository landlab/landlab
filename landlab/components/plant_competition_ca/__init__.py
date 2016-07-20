""" Landlab component that simulates inter-species plant competition.

This code is based on Cellular Automata Tree Grass Shrub Simulator (CATGraSS).
It simulates spatial competition of multiple plant functional types through
establishment and mortality. In the current code, tree, grass and
shrubs are used.

Ref: Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
ecohydrological role of aspect controlled radiation on tree grass shrub
coexistence in a semiarid climate. Water Resources Research,
49(5), 2872-2895.

.. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import VegCA

Create a grid on which to calculate soil moisture saturation fraction.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> sorted(VegCA.input_var_names)
['vegetation__cumulative_water_stress', 'vegetation__plant_functional_type']

Check the units for the fields.

>>> VegCA.var_units('vegetation__cumulative_water_stress')
'None'

Create the input fields.

>>> grid['cell']['vegetation__plant_functional_type']= np.arange(
...        0, grid.number_of_cells, dtype=int)

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> VegCA.var_help('vegetation__cumulative_water_stress')
name: vegetation__cumulative_water_stress
description:
  cumulative vegetation__water_stress over the growing season
units: None
at: cell
intent: in

>>> grid['cell']['vegetation__cumulative_water_stress'] = np.ones(
...        grid.number_of_cells)

Instantiate the 'VegCA' component to work on this grid,
and run it.

>>> ca_veg = VegCA(grid)

Run the *update* method to update output variables with current time

>>> ca_veg.update()

Check the output variable names

>>> sorted(VegCA.output_var_names)
['plant__age', 'plant__live_index']

>>> np.allclose(grid['cell']['plant__age'],0)
False
"""
from landlab.components.plant_competition_ca.plant_competition_ca import VegCA


__all__ = ['VegCA']
