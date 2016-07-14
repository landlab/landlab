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

>>> sorted(LandslideProbability.input_var_names) \
                # doctest: +NORMALIZE_WHITESPACE
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

>>> LandslideProbability.var_help('soil__transmissivity') \
                        # doctest: +NORMALIZE_WHITESPACE
name: soil__transmissivity
description:
  mode rate of water transmitted through a unit width of
saturated soil
units: m2/day
at: node
intent: in

Additional required fields for component.

>>> scatter_dat = np.random.random_integers(1, 10, grid.number_of_nodes)
>>> grid['node']['topographic__specific_contributing_area'] = \
         np.sort(np.random.random_integers(30, 900, grid.number_of_nodes))
>>> grid['node']['soil__transmissivity'] = \
         np.sort(np.random.random_integers(5, 20, grid.number_of_nodes), -1)
>>> grid['node']['soil__mode_total_cohesion'] = \
         np.sort(np.random.random_integers(30, 900, grid.number_of_nodes))
>>> grid['node']['soil__minimum_total_cohesion'] = \
         grid.at_node['soil__mode_total_cohesion'] - scatter_dat
>>> grid['node']['soil__maximum_total_cohesion'] = \
         grid.at_node['soil__mode_total_cohesion'] + scatter_dat
>>> grid['node']['soil__internal_friction_angle'] = \
         np.sort(np.random.random_integers(26, 40, grid.number_of_nodes))
>>> grid['node']['soil__thickness'] = \
         np.sort(np.random.random_integers(1, 10, grid.number_of_nodes))
>>> grid['node']['soil__density'] = \
         2000. * np.ones(grid.number_of_nodes)

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
>>> isinstance(LS_prob.landslide__factor_of_safety_histogram[ \
    core_nodes[0]], np.ndarray) == True
True
"""

# %% Import Libraries
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np


# %% Instantiate Object


class LandslideProbability(Component):
    """
    Landlab component designed to calculate a probability of failure at
    each grid node based on the infinite slope stability model
    stability index (Factor of Safety).

    The driving force for failure is provided by the user in the form of
    groundwater recharge, simply user provided minimum and maximum annual
    peak values of recharge. The model uses topographic and soils
    characteristics provided as input in the landslide_driver.

    A LandslideProbability calcuation function provides the user with the
    mean soil relative wetness, mean factor-of-safety, and probabilty
    of failure at each node.

    Construction::
        LandslideProbability(grid, number_of_simulations"=250,
        rechare_minimum=5., groundwater__recharge_maximum=120.)

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    number_of_simulations: float, optional
        Number of simulations to run Monte Carlo.
    groundwater__recharge_minimum: float, optional
        User provided minimum annual maximum recharge\
        recharge (mm/day).
    groundwater__recharge_maximum: float, optional
        User provided maximum annual maximum recharge\
        recharge (mm/day).

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.landslides import LandslideProbability
    >>> import numpy as np

    >>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    >>> LS_prob = LandslideProbability(grid)
    >>> LS_prob.name
    'Landslide Probability'
    >>> sorted(LandslideProbability.input_var_names) \
              # doctest: +NORMALIZE_WHITESPACE
    ['soil__density',
     'soil__internal_friction_angle',
     'soil__maximum_total_cohesion',
     'soil__minimum_total_cohesion',
     'soil__mode_total_cohesion',
     'soil__thickness',
     'soil__transmissivity',
     'topographic__slope',
     'topographic__specific_contributing_area']
    >>> sorted(LS_prob.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['landslide__mean_factor_of_safety',
     'landslide__probability_of_failure',
     'soil__mean_relative_wetness']
    >>> sorted(LS_prob.units) # doctest: +NORMALIZE_WHITESPACE
    [('landslide__mean_factor_of_safety', 'None'),
     ('landslide__probability_of_failure', 'None'),
     ('soil__density', 'kg/m3'),
     ('soil__internal_friction_angle', 'degrees'),
     ('soil__maximum_total_cohesion', 'Pa or kg/m-s2'),
     ('soil__mean_relative_wetness', 'None'),
     ('soil__minimum_total_cohesion', 'Pa or kg/m-s2'),
     ('soil__mode_total_cohesion', 'Pa or kg/m-s2'),
     ('soil__thickness', 'm'),
     ('soil__transmissivity', 'm2/day'),
     ('topographic__slope', 'tan theta'),
     ('topographic__specific_contributing_area', 'm')]

    >>> LS_prob.grid.number_of_node_rows
    5
    >>> LS_prob.grid.number_of_node_columns
    4
    >>> LS_prob.grid is grid
    True

    >>> grid['node']['topographic__slope'] = \
        np.random.rand(grid.number_of_nodes)
    >>> scatter_dat = np.random.random_integers(1, 10, grid.number_of_nodes)
    >>> grid['node']['topographic__specific_contributing_area'] = \
             np.sort(np.random.random_integers(30, 900, grid.number_of_nodes))
    >>> grid['node']['soil__transmissivity'] = \
             np.sort(np.random.random_integers(5, 20, grid.number_of_nodes),-1)
    >>> grid['node']['soil__mode_total_cohesion'] = \
             np.sort(np.random.random_integers(30, 900, grid.number_of_nodes))
    >>> grid['node']['soil__minimum_total_cohesion'] = \
             grid.at_node['soil__mode_total_cohesion'] - scatter_dat
    >>> grid['node']['soil__maximum_total_cohesion'] = \
             grid.at_node['soil__mode_total_cohesion'] + scatter_dat
    >>> grid['node']['soil__internal_friction_angle'] = \
             np.sort(np.random.random_integers(26, 40, grid.number_of_nodes))
    >>> grid['node']['soil__thickness']= \
             np.sort(np.random.random_integers(1, 10, grid.number_of_nodes))
    >>> grid['node']['soil__density'] = \
             2000. * np.ones(grid.number_of_nodes)

    >>> LS_prob = LandslideProbability(grid)
    >>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
    True
    >>> LS_prob.calculate_landslide_probability()
    >>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
    False
    >>> core_nodes = LS_prob.grid.core_nodes
    >>> isinstance(LS_prob.landslide__factor_of_safety_histogram[ \
        core_nodes[0]], np.ndarray) == True
    True
    """

# component name
    _name = 'Landslide Probability'
    __version__ = '1.0'
# component requires these values to do its calculation, get from driver
    _input_var_names = (
        'topographic__specific_contributing_area',
        'topographic__slope',
        'soil__transmissivity',
        'soil__mode_total_cohesion',
        'soil__minimum_total_cohesion',
        'soil__maximum_total_cohesion',
        'soil__internal_friction_angle',
        'soil__density',
        'soil__thickness',
        )

#  component creates these output values
    _output_var_names = (
        'soil__mean_relative_wetness',
        'landslide__mean_factor_of_safety',
        'landslide__probability_of_failure',
        )

# units for each parameter and output
    _var_units = {
        'topographic__specific_contributing_area': 'm',
        'topographic__slope': 'tan theta',
        'soil__transmissivity': 'm2/day',
        'soil__mode_total_cohesion': 'Pa or kg/m-s2',
        'soil__minimum_total_cohesion': 'Pa or kg/m-s2',
        'soil__maximum_total_cohesion': 'Pa or kg/m-s2',
        'soil__internal_friction_angle': 'degrees',
        'soil__density': 'kg/m3',
        'soil__thickness': 'm',
        'soil__mean_relative_wetness': 'None',
        'landslide__mean_factor_of_safety': 'None',
        'landslide__probability_of_failure': 'None',
        }

# grid centering of each field and variable
    _var_mapping = {
        'topographic__specific_contributing_area': 'node',
        'topographic__slope': 'node',
        'soil__transmissivity': 'node',
        'soil__mode_total_cohesion': 'node',
        'soil__minimum_total_cohesion': 'node',
        'soil__maximum_total_cohesion': 'node',
        'soil__internal_friction_angle': 'node',
        'soil__density': 'node',
        'soil__thickness': 'node',
        'soil__mean_relative_wetness': 'node',
        'landslide__mean_factor_of_safety': 'node',
        'landslide__probability_of_failure': 'node',
        }

# short description of each field
    _var_doc = {
        'topographic__specific_contributing_area':
            ('specific contributing (upslope area/cell face )' +
             ' that drains to node'),
        'topographic__slope':
        'slope of surface at node represented by tan theta',
        'soil__transmissivity':
            ('mode rate of water transmitted' +
             ' through a unit width of saturated soil'),
        'soil__mode_total_cohesion':
        'mode of combined root and soil cohesion at node',
        'soil__minimum_total_cohesion':
        'minimum of combined root and soil cohesion at node',
        'soil__maximum_total_cohesion':
        'maximum of combined root and soil cohesion at node',
        'soil__internal_friction_angle':
            ('critical angle just before failure' +
             ' due to friction between particles'),
        'soil__density': 'wet bulk density of soil',
        'soil__thickness': 'soil depth to restrictive layer',
        'soil__mean_relative_wetness':
            ('Indicator of soil wetness;' +
             ' relative depth perched water table' +
             ' within the soil layer'),
        'landslide__mean_factor_of_safety':
            ('(FS) dimensionless index of stability' +
             ' based on infinite slope stabiliity model'),
        'landslide__probability_of_failure':
            ('number of times FS is <1 out of number of' +
             ' interations user selected'),
        }

# Run Component
    @use_file_name_or_kwds
    def __init__(self, grid, number_of_simulations=250.,
                 groundwater__recharge_minimum=20.,
                 groundwater__recharge_maximum=120., **kwds):

        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        number_of_simulations: int, optional
            number of simulations to run Monte Carlo (None)
        groundwater__recharge_minimum: float, optional
            Minimum annual maximum recharge (mm/d)
        groundwater__recharge_maximum: float, optional
            Maximum annual maximum rechage (mm/d)
        g: float, optional
            acceleration due to gravity (m/sec^2)
        """

        # Store grid and parameters and do unit conversions
        self._grid = grid
        self.n = number_of_simulations
        self.recharge_min = groundwater__recharge_minimum/1000.0  # mm->m
        self.recharge_max = groundwater__recharge_maximum/1000.0
        self.g = 9.81

        super(LandslideProbability, self).__init__(grid)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._nodal_values = self.grid['node']

        # Raise an error if somehow someone is using this weird functionality
        if self._grid is None:
            raise ValueError('You must now provide an existing grid!')

    def calculate_factor_of_safety(self, i):

        """
        Method calculates factor-of-safety stability index by using
        node specific parameters, creating distributions of these parameters,
        and calculating the index by sampling these distributions 'n' times.

        The index is calculated from the 'infinite slope stabilty
        factor-of-safety equation' in the format of Pack RT, Tarboton DG,
        and Goodwin CN (1998)The SINMAP approach to terrain stability mapping.

        Parameters
        ----------
        i: int
            index of core node ID.
        """

        # generate distributions to sample from to provide input parameters
        # currently triangle distribution using mode, min, & max
        self.a = self.grid['node'][
            'topographic__specific_contributing_area'][i]
        self.theta = self.grid['node']['topographic__slope'][i]
        self.Tmode = self.grid['node']['soil__transmissivity'][i]
        self.Cmode = self.grid['node']['soil__mode_total_cohesion'][i]
        self.Cmin = self.grid['node']['soil__minimum_total_cohesion'][i]
        self.Cmax = self.grid['node']['soil__maximum_total_cohesion'][i]
        self.phi_mode = self.grid['node']['soil__internal_friction_angle'][i]
        self.rho = self.grid['node']['soil__density'][i]
        self.hs_mode = self.grid['node']['soil__thickness'][i]

        # Transmissivity (T)
        Tmin = self.Tmode-(0.3*self.Tmode)
        Tmax = self.Tmode+(0.3*self.Tmode)
        self.T = np.random.triangular(Tmin, self.Tmode, Tmax, size=self.n)
        # Cohesion
        # if provide fields of min and max C, uncomment 2 lines below
        #    Cmin = Cmode-0.3*self.Cmode
        #    Cmax = Cmode+0.3*self.Cmode
        self.C = np.random.triangular(self.Cmin, self.Cmode,
                                      self.Cmax, size=self.n)
        # phi - internal angle of friction provided in degrees
        phi_min = self.phi_mode-0.18*self.phi_mode
        phi_max = self.phi_mode+0.32*self.phi_mode
        self.phi = np.random.triangular(phi_min, self.phi_mode,
                                        phi_max, size=self.n)
        # soil thickness
        hs_min = self.hs_mode-0.3*self.hs_mode
        hs_max = self.hs_mode+0.3*self.hs_mode
        self.hs = np.random.triangular(hs_min, self.hs_mode,
                                       hs_max, size=self.n)
        # recharge distribution
        self.Re = np.random.uniform(self.recharge_min,
                                    self.recharge_max, size=self.n)
        # calculate Factor of Safety for n number of times
        # calculate components of FS equation
        self.C_dim = self.C/(self.hs*self.rho*self.g)  # demensionless cohesion
        self.Rel_wetness = ((self.Re)/self.T)*(self.a/np.sin(
            np.arctan(self.theta)))                       # relative wetness
        np.place(self.Rel_wetness, self.Rel_wetness > 1, 1.0)
        # maximum Rel_wetness = 1.0
        self.soil__mean_relative_wetness = np.mean(self.Rel_wetness)
        self.Y = np.tan(np.radians(self.phi))*(1 - (self.Rel_wetness*0.5))
        # convert from degrees; 0.5 = water to soil density ratio
        # calculate Factor-of-safety
        self.FS = (self.C_dim/np.sin(np.arctan(self.theta))) + (
            np.cos(np.arctan(self.theta)) *
            (self.Y/np.sin(np.arctan(self.theta))))
        self.FS_store = np.array(self.FS)        # array of factor of safety
        self.FS_distribution = self.FS_store
        self.landslide__mean_factor_of_safety = np.mean(self.FS)
        count = 0
        for val in self.FS:                   # find how many FS values <= 1
            if val <= 1.0:
                count = count + 1
        self.FS_L1 = float(count)     # number with unstable FS values (<=1)
        # probability: No. unstable values/total No. of values (n)
        self.landslide__probability_of_failure = self.FS_L1/self.n

    def calculate_landslide_probability(self, **kwds):

        """
        Method creates arrays for output variables then loops through all
        the core nodes to run the method 'calculate_factor_of_safety.'
        Some output variables are assigned as fields to nodes. One output
        parameter is an factor-of-safety distribution at each node.

        Parameters
        ----------
        self.landslide__factor_of_safety_histogram: numpy.ndarray([
            self.grid.number_of_nodes, self.n], dtype=float)
            This is an output - distribution of factor-of-safety from
            Monte Carlo simulations (units='None')
        """

        # Create arrays for data with -9999 as default to store output
        self.mean_Relative_Wetness = -9999*np.ones(self.grid.number_of_nodes,
                                                   dtype='float')
        self.mean_FS = -9999*np.ones(self.grid.number_of_nodes, dtype='float')
        self.prob_fail = -9999*np.ones(
            self.grid.number_of_nodes, dtype='float')
        self.landslide__factor_of_safety_histogram = -9999*np.ones(
            [self.grid.number_of_nodes, self.n], dtype='float')
        # Run factor of safety Monte Carlo for all core nodes in domain
        # i refers to each core node id
        for i in self.grid.core_nodes:
            self.calculate_factor_of_safety(i)
            # Populate storage arrays with calculated values
            self.mean_Relative_Wetness[i] = self.soil__mean_relative_wetness
            self.mean_FS[i] = self.landslide__mean_factor_of_safety
            self.prob_fail[i] = self.landslide__probability_of_failure
            self.landslide__factor_of_safety_histogram[i] = \
                self.FS_distribution
            # stores FS values from last loop (node)
        # replace unrealistic values in arrays
        self.mean_Relative_Wetness[
            self.mean_Relative_Wetness < 0.] = 0.  # so can't be negative
        self.mean_FS[self.mean_FS < 0.] = 0.       # can't be negative
        self.mean_FS[self.mean_FS == np.inf] = 0.  # to deal with NaN in data
        self.prob_fail[self.prob_fail < 0.] = 0.   # can't be negative
        # assign output fields to nodes
        self.grid['node']['soil__mean_relative_wetness'] =\
            self.mean_Relative_Wetness
        self.grid['node']['landslide__mean_factor_of_safety'] = self.mean_FS
        self.grid['node']['landslide__probability_of_failure'] = self.prob_fail
