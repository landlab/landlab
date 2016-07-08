""" Landlab component that simulates relative wetness, mean factor-of-safety,
and probability of failure.

Relative wetness and factor-of-safety are based on the infinite slope
stability model driven by topographic and soils inputs and recharge provided
by user in a "landslide_driver" file. In addition, the component simulates
the probability of failure for each node based on Monte Carlo simulations
of the factor-of-safety as the number of simulations <= 1.0 divided by
the number of simulations.

Modified to be more generic in the inputs for greater usability as well
as accomodate functionality with new release of Landlab version 1.

.. codauthor:: R.Strauch & E.Istanbulluoglu - University of Washington
Created on Thu Aug 20, 2015
Last edit July 7, 2016

Examples
----------
>>> from landlab import RasterModelGrid
>>> from landlab.components import LandslideProbability
>>> import numpy as np

Create a grid on which to calculate landslide probability.

>>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))

The grid will need some input data. To check the names of the fields
that provide the input to this component, use the *input_var_names*
class property.

>>> LandslideProbabiity.input_var_names
('number_of_simulations', 'contributing_area',
 'slope', 'soil_transmissivity__mode', 'combined_cohesion__mode',
 'combined_cohesion__minimum', 'combined_cohesion__maximum',
 'soil_internal_angle_friction__mode', 'soil_density',
 'soil_thickness__mode', 'recharge__minimum', 'recharge__maximum')

Check the units for the fields.

>>> LandslideProbability.var_units('recharge__minimum')
    'mm/day'

Create the input fields.

>>> grid['node']['slope'] = np.random.rand(5,4)

If you are not sure about one of the input or output variables, you can
get help for specific variables.

>>> LandslideProbability.var_help('soil_transmissivity__mode')
name: soil_transmissivity__mode
description:
  mode rate of water transmitted through a unit width of saturated soil
units: m2/day
at: node
intent: in

Additional required fields for component.

>>> scatter_dat = np.random.random_integers(1, 10, (5,4))
>>> grid['node']['contributing_area']= \
         np.sort(np.random.random_integers(30, 900, (5,4)))
>>> grid['node']['soil_transmissivity']= \
         np.sort(np.random.random_integers(5, 20, (5,4)),-1)
>>> grid['node']['combined_cohesion__mode']= \
         np.sort(np.random.random_integers(30, 900, (5,4)))
>>> grid['node']['combined_cohesion__minimum']= \
         grid.at_node['combined_cohesion__mode'] - scatter_dat
>>> grid['node']['combined_cohesion__maximum']= \
         grid.at_node['combined_cohesion__mode'] + scatter_dat
>>> grid['node']['soil_internal_angle_friction__mode']= \
         np.sort(np.random.random_integers(26, 40, (5,4)))
>>> grid['node']['soil_thickness__mode']= \
         np.sort(np.random.random_integers(1, 10, (5,4)))
>>> grid['node']['soil_density']= \
         2000. * np.ones(grid.number_of_nodes)

Instantiate the 'LandslideProbability' component to work on this grid,
and run it.

>>> LS_prob = LandslideProbability(grid)

Check the output variable names.

>>> sorted(LS_prob.output_var_names) # doctest: +NORMALIZE_WHITESPACE
['Relative_Wetness__mean', 'Factor_of_Safety__mean',
 'Probability_of_failure', 'Factor_of_Safety__distribution']

Check the output from the component.

>>> grid['node']['Relative_Wetness__mean']
>>> grid['node']['Factor_of_Safety__mean']
>>> grid['node']['Probability_of_failure]
>>> grid['node']['Factor_of_Safety__distribution'][1]
"""

# %% Import Libraries
from landlab import Component
from ...utils.decorators import use_file_name_or_kwds
import numpy as np

_VALID_METHODS = set(['Grid', 'Multi'])


def assert_method_is_valid(method):
    if method not in _VALID_METHODS:
        raise ValueError('%s: Invalid method name' % method)

# %% Instantiate Object


class LandslideProbability(Component):
    """
    Landlab component designed to calculate a probability of failure at
    each grid nodes based on the infinite slope stability model
    stability index (Factor of Safety).

    The driving force for failure is provided by the user in the form of
    groundwater recharge, simply user provide minimum and maximum annual
    peak values. The model uses topographic and soils characteristics
    provided as input in the landslide_driver.

    A LandslideProbability calcuation function provides the user with the
    mean soil relative wetness, mean factor-of-safety, and probabilty
    of failure at each node.

    Construction:: # NEED do I need the "method = 'Grid'???
        LandslideProbability(grid, method='Grid', number_of_simulations"=250,
        rechare_minimum=5., recharge__maximum=120.)

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    method: {'Grid'}, optional
        Currently, only default is available.
    number_of_simulations: float, optional
        Number of simulations to run Monte Carlo.
    recharge__minimum: float, optional
        User provided minimum annual maximum recharge\
        recharge (mm/day).
    recharge__maximum: float, optional
        User provided maximum annual maximum recharge\
        recharge (mm/day).

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import LandslideProbability
    >>> import numpy as np

    >>> grid = RasterModelGrid((5, 4), spacing=(0.2, 0.2))
    >>> LS_prob = LandslideProbability(grid)
    >>> LS_prob.name
    >>> LS_prob.input_var_names
    ('contributing_area','slope', 'soil_transmissivity__mode',
     'combined_cohesion__mode', 'combined_cohesion__minimum',
     'combined_cohesion__maximum', 'soil_internal_angle_friction__mode',
     'soil_density', 'soil_thickness__mode')
    >>> sorted(LS_prob.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['Relative_Wetness__mean', 'Factor_of_Safety__mean',
        'Probability_of_failure', 'Factor_of_Safety__distribution']
    >>> sorted(LS_prob.units) # doctest: +NORMALIZE_WHITESPACE
    [('contributing_area', 'm'),
     ('slope', 'tan theta'),
     ('soil_transmissivity__mode', 'm2/day'),
     ('combined_cohesion__mode', 'Pa or kg/m-s2'),
     ('combined_cohesion__minimum', 'Pa or kg/m-s2'),
     ('combined_cohesion__maximum', 'Pa or kg/m-s2'),
     ('soil_internal_angle_friction__mode', 'degrees'),
     ('soil_density', 'kg/m3'),
     ('soil_thickness__mode', 'm'),
     ('Relative_Wetness__mean', 'None'),
     ('Factor-of-Safety__mean', 'None'),
     ('Probability_of_failure', 'None'),
     ('Factor_of_Safety__distribution', 'None')]

    >>> LS_prob.grid.number_of_node_rows
    5
    >>> LS_prob.grid.number_of_node_columns
    4
    >>> LS_prob.grid is grid
    True
    >>> np.allclose(grid.at_node['slope'], 0.)
    True
# NEED this last input, might raise a fielderror(name)
# because I haven't filled the field yet?

    >>> grid['node']['slope']= np.random.rand(5,4)
    >>> scatter_dat = np.random.random_integers(1, 10, (5,4))
    >>> grid['node']['contributing_area']= \
             np.sort(np.random.random_integers(30, 900, (5,4)))
    >>> grid['node']['soil_transmissivity']= \
             np.sort(np.random.random_integers(5, 20, (5,4)),-1)
    >>> grid['node']['combined_cohesion__mode']= \
             np.sort(np.random.random_integers(30, 900, (5,4)))
    >>> grid['node']['combined_cohesion__minimum']= \
             grid.at_node['combined_cohesion__mode'] - scatter_dat
    >>> grid['node']['combined_cohesion__maximum']= \
             grid.at_node['combined_cohesion__mode'] + scatter_dat
    >>> grid['node']['soil_internal_angle_friction__mode']= \
             np.sort(np.random.random_integers(26, 40, (5,4)))
    >>> grid['node']['soil_thickness__mode']= \
             np.sort(np.random.random_integers(1, 10, (5,4)))
    >>> grid['node']['soil_density']= \
             2000. * np.ones(grid.number_of_nodes)

    >>> LS_prob = LandslideProbability(grid)
    >>> np.all(grid.at_node['Probability_of_failure'] == 0.)
    False
    """

# component name
    _name = 'Landslide Probability'
    __version__ = '1.0'
# component requires these values to do its calculation, get from driver
    _input_var_names = (
        'contributing_area',
        'slope',
        'soil_transmissivity__mode',
        'combined_cohesion__mode',
        'combined_cohesion__minimum',
        'combined_cohesion__maximum',
        'soil_internal_angle_friction__mode',
        'soil_density',
        'soil_thickness__mode',
        )

#  component creates these output values
    _output_var_names = (
        'Relative_Wetness__mean',
        'Factor_of_Safety__mean',
        'Probability_of_failure',
        'Factor_of_Safety__distribution',
        )

# units for each parameter and output
    _var_units = {
        'contributing_area': 'm',
        'slope': 'tan theta',
        'soil_transmissivity__mode': 'm2/day',
        'combined_cohesion__mode': 'Pa or kg/m-s2',
        'combined_cohesion__minimum': 'Pa or kg/m-s2',
        'combined_cohesion__maximum': 'Pa or kg/m-s2',
        'soil_internal_angle_friction__mode': 'degrees',
        'soil_density': 'kg/m3',
        'soil_thickness__mode': 'm',
        'Relative_Wetness__mean': 'None',
        'Factor-of-Safety__mean': 'None',
        'Probability_of_failure': 'None',
        'Factor_of_Safety__distribution': 'None',
        }

# grid centering of each field
    _var_mapping = {
        'contributing_area': 'node',
        'slope': 'node',
        'soil_transmissivity__mode': 'node',
        'combined_cohesion__mode': 'node',
        'combined_cohesion__minimum': 'node',
        'combined_cohesion__maximum': 'node',
        'soil_internal_angle_friction__mode': 'node',
        'soil_density': 'node',
        'soil_thickness__mode': 'node',
        'Relative_Wetness__mean': 'node',
        'Factor_of_Safety__mean': 'node',
        'Probability_of_failure': 'node',
        'Factor_of_Safety__distribution': 'node',
        }

# short description of each field
    _var_doc = {
        'contributing_area': 'specific contributing area\
        (upslope area/cell face length) that drains to node',
        'slope': 'slope of surface at node represented by tan theta',
        'soil_transmissivity__mode': 'mode rate of water transmitted\
        through a unit width of saturated soil',
        'combined_cohesion__mode': 'mode combined root and soil\
        cohesion at node',
        'combined_cohesion__minimum': 'minimum combined root and soil\
        cohesion at node',
        'combined_cohesion__maximum': 'maximum combined root and soil\
        cohesion at node',
        'soil_internal_angle_friction__mode': 'critical angle just before\
        failure due to friction between particles',
        'soil_density': 'wet bulk density of soil',
        'soil_thickness__mode': 'soil depth to restrictive layer',
        'Relative_Wetness__mean': 'Indicator of soil wetness; relative depth\
        perched water table within the soil layer',
        'Factor_of_Safety__mean': '(FS) dimensionless index of stability\
        based on infinite slope stabiliity model',
        'Probability_of_failure': 'number of times FS is <1 out of number of\
        interations user selected',
        'Factor_of_Safety__distribution': 'distribution of factor of safety\
        from Monte Carlo simulations',
        }

# Run Component
    @use_file_name_or_kwds
    def __init__(self, grid, method='Grid', number_of_simulations=1000.,
                 recharge__minimum=20., recharge__maximum=120., **kwds):
        self.n = number_of_simulations
        self.a = self.grid['node']['contributing_area']
        self.theta = self.grid['node']['slope ']
        self.Tmode = self.grid['node']['soil_transmissivity__mode']
        self.Cmode = self.grid['node']['combined_cohesion__mode']
        self.Cmin = self.grid['node']['combined_cohesion__minimum']
        self.Cmax = self.grid['node']['combined_cohesion__maximum']
        self.phi_mode = self.grid['node']['soil_internal_angle_friction__mode']
        self.rho = self.grid['node']['soil_density']
        self.hs_mode = self.grid['node']['soil_thickness__mode']
        self.recharge_min = recharge__minimum
        self.recharge_max = recharge__maximum
        self.g = 9.81
        self.update()

        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        method : {'Grid'}, optional
            Currently, only default is available.
        number_of_simulations: int, optional
            number of simulations to run Monte Carlo (None)
        recharge__minimum: float, optional
            Minimum annual maximum recharge (mm/d).
        recharge__maximum: float, optional
            Maximum annual maximum rechage (mm/d).
        """
        self._method = kwds.pop('method', 'Grid')  # NEED???

        assert_method_is_valid(self._method)

        super(LandslideProbability, self).__init__(grid)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._nodal_values = self.grid['node']

    def calculate_factor_of_safety(self):
        # generate distributions to sample from to provide input parameters
        # currently triangle distribution using mode, min, & max
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
        self.Rel_wetness = ((self.Re/1000.0)/self.T)*(self.a/np.sin(
            np.arctan(self.theta)))                       # relative wetness
        np.place(self.Rel_wetness, self.Rel_wetness > 1, 1.0)
        # maximum Rel_wetness = 1.0
        self.Relative_Wetness__mean = np.mean(self.Rel_wetness)
        self.Y = np.tan(np.radians(self.phi))*(1 - (self.Rel_wetness*0.5))
        # convert from degrees; 0.5 = water to soil density ratio
        # calculate Factor-of-safety
        self.FS = (self.C_dim/np.sin(np.arctan(self.theta))) + (
            np.cos(np.arctan(self.theta)) *
            (self.Y/np.sin(np.arctan(self.theta))))
        self.FS_store = np.array(self.FS)        # array of factor of safety
        self.Factor_of_Safety__mean = np.mean(self.FS)
        count = 0
        for val in self.FS:                   # find how many FS values <= 1
            if val <= 1.0:
                count = count + 1
        self.FS_L1 = float(count)     # number with unstable FS values (<=1)
        # probability: No. unstable values/total No. of values (n)
        self.Probability_of_failure = self.FS_L1/self.n
        self.Factor_of_Safety__distribution = self.FS_store

    def update(self, **kwds):

        # Create arrays for data with -9999 as default to store output
        self.mean_Relative_Wetness = -9999*np.ones(self.grid.number_of_nodes,
                                                   dtype='float')
        self.mean_FS = -9999*np.ones(self.grid.number_of_nodes, dtype='float')
        self.prob_fail = -9999*np.ones(
            self.grid.number_of_nodes, dtype='float')
        self.FS_dist = -9999*np.ones(
            [self.grid.number_of_nodes, self.n], dtype='float')
        # Run factor of safety Monte Carlo for all core nodes in domain
        # i refers to each core node id
        for i in self.grid.core_nodes:
            self.calculate_factor_of_safety(
                self.n, self.grid['node']['contributing_area'][i],
                self.grid['node']['slope'][i],
                self.grid['node']['soil_transmissivity__mode'][i],
                self.grid['node']['combined_cohesion__mode'][i],
                self.grid['node']['combined_cohesion__minimum'][i],
                self.grid['node']['combined_cohesion__maximum'][i],
                self.grid['node']['soil_internal_angle_friction__mode'][i],
                self.grid['node']['soil_density'][i],
                self.grid['node']['soil_thickness__mode'][i],
                self.Re)   # parameters & data passed to FS class
            # Populate storage arrays with calculated values
            self.mean_Relative_Wetness[i] = self.Relative_Wetness__mean
            self.mean_FS[i] = self.Factor_of_Safety__mean
            self.prob_fail[i] = self.Probability_of_failure
            self.FS_dist[i] = self.Factor_of_Safety__distribution
            # stores FS values from last loop (node)
        # replace unrealistic values in arrays
        self.mean_Relative_Wetness[
            self.mean_Relative_Wetness < 0.] = 0.  # so can't be negative
        self.mean_FS[self.mean_FS < 0.] = 0.       # can't be negative
        self.mean_FS[self.mean_FS == np.inf] = 0.  # to deal with NaN in data
        self.prob_fail[self.prob_fail < 0.] = 0.   # can't be negative
        self.FS_dist[self.FS_dist < 0.] = 0.     # can't be negative
