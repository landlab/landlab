""" Landlab component that simulates root-zone average soil moisture
saturation fraction.

This component uses a single soil moisture layer and models soil moisture
loss through transpiration by plants, evaporation by bare soil, and leakage.
The solution of water balance is based on Laio et. al 2001. The component
requires fields of initial soil moisture, rainfall input (if any), time
to the next storm and potential transpiration.

Ref: Laio, F., Porporato, A., Ridolfi, L., & Rodriguez-Iturbe, I. (2001).
Plants in water-controlled ecosystems: active role in hydrologic processes
and response to water stress: II. Probabilistic soil moisture dynamics.
Advances in Water Resources, 24(7), 707-723.

.. codeauthor:: Sai Nudurupati and Erkan Istanbulluoglu

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.soil_moisture import SoilMoisture

Create a grid on which to calculate soil moisture saturation fraction.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> sorted(SoilMoisture.input_var_names)   # doctest: +NORMALIZE_WHITESPACE
['rainfall__daily_depth',
 'soil_moisture__initial_saturation_fraction',
 'surface__potential_evapotranspiration_rate',
 'vegetation__cover_fraction',
 'vegetation__live_leaf_area_index',
 'vegetation__plant_functional_type']

Check the units for the fields.

>>> SoilMoisture.var_units('surface__potential_evapotranspiration_rate')
'mm'

Create the input fields.

>>> grid['cell']['surface__potential_evapotranspiration_rate']= np.array([
...        0.25547770, 0.25547770, 0.22110221,
...        0.22110221, 0.24813062, 0.24813062])

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> SoilMoisture.var_help('vegetation__cover_fraction')
name: vegetation__cover_fraction
description:
  fraction of land covered by vegetation
units: None
at: cell
intent: in

>>> grid['cell']['soil_moisture__initial_saturation_fraction'] = (
...        0.75 * np.ones(grid.number_of_cells))

>>> grid['cell']['vegetation__plant_functional_type']= (
...        np.zeros(grid.number_of_cells, dtype=int))

>>> grid['cell']['vegetation__live_leaf_area_index']= (
...        2. * np.ones(grid.number_of_cells))

>>> grid['cell']['vegetation__cover_fraction']= (
...        np.ones(grid.number_of_cells))

Instantiate the 'SoilMoisture' component to work on this grid,
and run it.

>>> SM = SoilMoisture(grid)

Run the *update* method to update output variables with current time

>>> current_time = 0.5

>>> grid['cell']['rainfall__daily_depth'] = 25. * np.ones(grid.number_of_cells)

>>> current_time = SM.update(current_time)

Check the output variable names

>>> sorted(SoilMoisture.output_var_names)   # doctest: +NORMALIZE_WHITESPACE
['soil_moisture__root_zone_leakage',
 'soil_moisture__saturation_fraction',
 'surface__evapotranspiration',
 'surface__runoff',
 'vegetation__water_stress']

>>> grid['cell']['soil_moisture__root_zone_leakage']
array([ 29.97001453,  29.97001453,  29.97001453,  29.97001453,
        29.97001453,  29.97001453])

>>> SM.grid.at_cell['soil_moisture__saturation_fraction']
array([ 0.70372004,  0.70372004,  0.70372004,  0.70372004,  0.70372004,
        0.70372004])

>>> grid['cell']['surface__evapotranspiration']
array([ 0.0001,  0.0001,  0.0001,  0.0001,  0.0001,  0.0001])

>>> grid['cell']['surface__runoff']
array([ 0.,  0.,  0.,  0.,  0.,  0.])
"""

from .soil_moisture_dynamics import SoilMoisture
from .infiltrate_soil_green_ampt import SoilInfiltrationGreenAmpt


__all__ = ['SoilMoisture',
           'SoilInfiltrationGreenAmpt', ]
