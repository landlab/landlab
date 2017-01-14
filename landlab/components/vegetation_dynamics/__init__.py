"""Landlab component that implements vegetation dynamics model.

Landlab component that simulates net primary productivity,
above ground biomass, and leaf area index at each cell based on
inputs of root-zone average soil moisture.

Ref: Zhou, X., Istanbulluoglu, E., & Vivoni, E. R. (2013). Modeling the
ecohydrological role of aspect controlled radiation on tree grass shrub
coexistence in a semiarid climate. Water Resources Research,
49(5), 2872-2895.

.. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components import Vegetation

Create a grid on which to simulate vegetation dynamics.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> sorted(Vegetation.input_var_names)  # doctest: +NORMALIZE_WHITESPACE
['surface__evapotranspiration',
 'surface__potential_evapotranspiration_30day_mean',
 'surface__potential_evapotranspiration_rate',
 'vegetation__plant_functional_type',
 'vegetation__water_stress']

Check the units for an input field

>>> Vegetation.var_units('surface__evapotranspiration')
'mm'

Create the input fields.

>>> grid['cell']['vegetation__plant_functional_type']= np.zeros(
...        grid.number_of_cells, dtype=int)

>>> grid['cell']['surface__evapotranspiration'] = (
...        0.2 * np.ones(grid.number_of_cells))

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> Vegetation.var_help('surface__potential_evapotranspiration_30day_mean')
name: surface__potential_evapotranspiration_30day_mean
description:
  30 day mean of surface__potential_evapotranspiration
units: mm
at: cell
intent: in

>>> grid['cell']['surface__potential_evapotranspiration_rate']= np.array([
...        0.25547770, 0.25547770, 0.22110221,
...        0.22110221, 0.24813062, 0.24813062])

>>> grid['cell']['surface__potential_evapotranspiration_30day_mean']= (
...        np.array([0.25547770, 0.25547770, 0.22110221,
...                  0.22110221, 0.24813062, 0.24813062]))

>>> grid['cell']['vegetation__water_stress'] = (
...        0.01 * np.ones(grid.number_of_cells))

Instantiate the 'SoilMoisture' component to work on this grid,
and run it.

>>> Veg = Vegetation(grid)

Run the *update* method to update output variables

>>> Veg.update()

Check the output variable names

>>> sorted(Vegetation.output_var_names)  # doctest: +NORMALIZE_WHITESPACE
['vegetation__cover_fraction',
 'vegetation__dead_biomass',
 'vegetation__dead_leaf_area_index',
 'vegetation__live_biomass',
 'vegetation__live_leaf_area_index']

>>> grid['cell']['vegetation__cover_fraction']
array([ 0.77686984,  0.77686984,  0.77686984,  0.77686984,  0.77686984,
        0.77686984])

>>> grid['cell']['vegetation__live_leaf_area_index']
array([ 0.47371568,  0.47371568,  0.47371568,  0.47371568,  0.47371568,
        0.47371568])
"""
from .vegetation_dynamics import Vegetation


__all__ = ['Vegetation', ]
