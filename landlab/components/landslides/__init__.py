""" Landlab component that simulates relative wetness, mean factor-of-safety,
and probability of failure.

Relative wetness and factor-of-safety are based on the infinite slope
stability model driven by topographic and soils inputs and recharge provided
by user in a "landslide_driver" file. In addition, the component simulates
the probability of failure for each node based on Monte Carlo simulations
of the factor-of-safety as the number of simulations with factor-of-safety
<= 1.0 divided by the number of simulations.

Modified to be more generic in the inputs for greater usability as well
as accomodate functionality with new release of Landlab version 1.

.. codauthor:: R.Strauch & E.Istanbulluoglu - University of Washington
Created on Thu Aug 20, 2015
Last edit July 12, 2016

Examples
----------
>>> from landlab import RasterModelGrid
>>> from landlab.components.landslides import LandslideProbability
>>> import numpy as np

Create a grid on which to calculate landslide probability.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

Check the number of core nodes.

>>> grid.number_of_core_nodes
6

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> sorted(LandslideProbability.input_var_names)  # doctest: +NORMALIZE_WHITESPACE
['soil__density',
 'soil__internal_friction_angle',
 'soil__maximum_total_cohesion',
 'soil__minimum_total_cohesion',
 'soil__mode_total_cohesion',
 'soil__thickness',
 'soil__transmissivity',
 'topographic__slope',
 'topographic__specific_contributing_area']

Check the units for the fields.

>>> LandslideProbability.var_units('topographic__specific_contributing_area')
'm'

Create an input field.

>>> grid['node']['topographic__slope'] = np.random.rand(grid.number_of_nodes)

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> LandslideProbability.var_help('soil__transmissivity')  # doctest: +NORMALIZE_WHITESPACE
name: soil__transmissivity
description:
  mode rate of water transmitted through a unit width of
saturated soil
units: m2/day
at: node
intent: in

Additional required fields for component.

>>> scatter_dat = np.random.random_integers(1, 10, grid.number_of_nodes)
>>> grid['node']['topographic__specific_contributing_area'] = np.sort(
...      np.random.random_integers(30, 900, grid.number_of_nodes))
>>> grid['node']['soil__transmissivity'] = np.sort(
...      np.random.random_integers(5, 20, grid.number_of_nodes), -1)
>>> grid['node']['soil__mode_total_cohesion'] = np.sort(
...      np.random.random_integers(30, 900, grid.number_of_nodes))
>>> grid['node']['soil__minimum_total_cohesion'] = (
...      grid.at_node['soil__mode_total_cohesion'] - scatter_dat)
>>> grid['node']['soil__maximum_total_cohesion'] = (
...      grid.at_node['soil__mode_total_cohesion'] + scatter_dat)
>>> grid['node']['soil__internal_friction_angle'] = np.sort(
...      np.random.random_integers(26, 40, grid.number_of_nodes))
>>> grid['node']['soil__thickness'] = np.sort(
...      np.random.random_integers(1, 10, grid.number_of_nodes))
>>> grid['node']['soil__density'] = (2000. * np.ones(grid.number_of_nodes))

Instantiate the 'LandslideProbability' component to work on this grid,
and run it.

>>> LS_prob = LandslideProbability(grid)
>>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
True

Run the *calculate_landslide_probability* method to update output
variables with grid

>>> LS_prob.calculate_landslide_probability()

Check the output variable names.

>>> sorted(LS_prob.output_var_names) # doctest: +NORMALIZE_WHITESPACE
['landslide__mean_factor_of_safety',
 'landslide__probability_of_failure',
 'soil__mean_relative_wetness']

Check the output from the component, including array at one node.

>>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
False
>>> core_nodes = LS_prob.grid.core_nodes
>>> (isinstance(LS_prob.landslide__factor_of_safety_histogram[
...      core_nodes[0]], np.ndarray) == True)
True
"""
from .landslide import LandslideProbability


__all__ = ['LandslideProbability', ]
