#!/usr/env/python

"""
landslide_probability.py:

Landlab component that simulates relative wetness and probability of failure.

Relative wetness and factor-of-safety are based on the infinite slope
stability model driven by topographic and soils inputs and recharge provided
by user in a "landslide_driver" file. In addition, the component simulates
the probability of failure for each node based on Monte Carlo simulations
of the factor-of-safety as the number of simulations with factor-of-safety
<= 1.0 divided by the number of simulations.

Modified to be more generic in the inputs for greater usability as well
as accomodate functionality with new release of Landlab version 1.

.. codauthor:: R.Strauch, E.Istanbulluoglu, & S.Nudurupati
University of Washington
Created on Thu Aug 20, 2015
Last edit April 04, 2017
"""

# %% Import Libraries
from landlab import Component
from landlab.utils.decorators import use_file_name_or_kwds
import numpy as np
from scipy import interpolate
from statsmodels.distributions.empirical_distribution import ECDF
import copy

# %% Instantiate Object


class LandslideProbability(Component):
    """
    Landlab component designed to calculate a probability of failure at
    each grid node based on the infinite slope stability model
    stability index (Factor of Safety).

    The driving force for failure is provided by the user in the form of
    groundwater recharge; 4 options for providing recharge are included.
    The model uses topographic and soils characteristics provided as input
    by the user in the landslide_driver.

    A LandslideProbability calcuation function provides the user with the
    mean soil relative wetness and probabilty of failure at each node.

    Construction::
        LandslideProbability(grid, number_of_iterations=250,
        recharge_minimum=5., groundwater__recharge_maximum=120.)

    Parameters
    ----------
    grid: RasterModelGrid
        A grid.
    number_of_iterations: float, optional
        Number of iterations to run Monte Carlo.
    groundwater__recharge_distribution: str, optional
        single word indicating recharge distribution, either 'uniform',
        'lognormal', 'lognormal_spatial,' or 'data_driven_spatial' (None)
    groundwater__recharge_min_value: float, optional
        minium groundwater recharge for 'uniform' (mm/d)
    groundwater__recharge_max_value: float, optional
        maximum groundwater recharge for 'uniform' (mm/d)
    groundwater__recharge_mean: float, optional
        mean grounwater recharge for 'lognormal' (mm/d)
    groundwater__recharge_standard_deviation: float, optional
        standard deviation of grounwater recharge for 'lognormal' (mm/d)
    groundwater__recharge_HSD_inputs: list, optional
        list of 3 dictionaries in order - HSD_dict {Hydroligic Source
        Domain (HSD) keys: recharge numpy array values}, {node IDs keys:
        list of HSD_Id values}, HSD_fractions {node IDS keys: list of
        HSD fractions values} (none)

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
    
    >>> grid.at_node['topographic__slope'] = np.random.rand(grid.number_of_nodes)
    
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
    
    >>> scatter_dat = np.random.randint(1, 10, grid.number_of_nodes)
    >>> grid.at_node['topographic__specific_contributing_area'] = np.sort(
    ...      np.random.randint(30, 900, grid.number_of_nodes))
    >>> grid.at_node['soil__transmissivity'] = np.sort(
    ...      np.random.randint(5, 20, grid.number_of_nodes), -1)
    >>> grid.at_node['soil__mode_total_cohesion'] = np.sort(
    ...      np.random.randint(30, 900, grid.number_of_nodes))
    >>> grid.at_node['soil__minimum_total_cohesion'] = (
    ...      grid.at_node['soil__mode_total_cohesion'] - scatter_dat)
    >>> grid.at_node['soil__maximum_total_cohesion'] = (
    ...      grid.at_node['soil__mode_total_cohesion'] + scatter_dat)
    >>> grid.at_node['soil__internal_friction_angle'] = np.sort(
    ...      np.random.randint(26, 40, grid.number_of_nodes))
    >>> grid.at_node['soil__thickness'] = np.sort(
    ...      np.random.randint(1, 10, grid.number_of_nodes))
    >>> grid.at_node['soil__density'] = (2000. * np.ones(grid.number_of_nodes))
    
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
    ['landslide__probability_of_failure', 'soil__mean_relative_wetness']
    
    Check the output from the component, including array at one node.
    
    >>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
    False
    >>> core_nodes = LS_prob.grid.core_nodes
    >>> (isinstance(LS_prob.landslide__factor_of_safety_distribution[
    ...      core_nodes[0]], np.ndarray) == True)
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
        'landslide__probability_of_failure':
            ('number of times FS is <=1 out of number of' +
             ' iterations user selected'),
        }

# Run Component
    @use_file_name_or_kwds
    def __init__(self, grid, number_of_iterations=250,
                 groundwater__recharge_distribution='uniform',
                 groundwater__recharge_min_value=20.,
                 groundwater__recharge_max_value=120.,
                 groundwater__recharge_mean=None,
                 groundwater__recharge_standard_deviation=None,
                 groundwater__recharge_HSD_inputs=[],
                 seed=0, **kwds):

        """
        Parameters
        ----------
        grid: RasterModelGrid
            A grid.
        number_of_iterations: int, optional
            number of iterations to run Monte Carlo simulation (None)
        groundwater__recharge_distribution: str, optional
            single word indicating recharge distribution, either 'uniform',
            'lognormal', 'lognormal_spatial,' or 'data_driven_spatial' (None)
        groundwater__recharge_min_value: float, optional
            minium groundwater recharge for 'uniform' (mm/d)
        groundwater__recharge_max_value: float, optional
            maximum groundwater recharge for 'uniform' (mm/d)
        groundwater__recharge_mean: float, optional
            mean grounwater recharge for 'lognormal' (mm/d)
        groundwater__recharge_standard_deviation: float, optional
            standard deviation of grounwater recharge for 'lognormal' (mm/d)
        groundwater__recharge_HSD_inputs: list, optional
            list of 3 dictionaries in order - HSD_dict {Hydroligic Source
            Domain (HSD) keys: recharge numpy array values}, {node IDs keys:
            list of HSD_Id values}, HSD_fractions {node IDS keys: list of
            HSD fractions values} (none)
        g: float, optional
            acceleration due to gravity (m/sec^2)
        seed: int, optional
            seed for random number generation. if seed is assigned any value
            other than the default value of zero, it will create different
            sequence. to create a certain sequence repititively, use the same
            value as input for seed.
        """

        # Initialize seeded random number generation        
        self.seed_generator(seed)

        # Store grid and parameters and do unit conversions
        self._grid = grid
        self.n = int(number_of_iterations)
        self.g = 9.81
        self.groundwater__recharge_distribution = (
            groundwater__recharge_distribution)
        # Following code will deal with the input distribution and associated
        # parameters
        # Uniform distribution
        if self.groundwater__recharge_distribution == 'uniform':
            self.recharge_min = groundwater__recharge_min_value
            self.recharge_max = groundwater__recharge_max_value
            self.Re = np.random.uniform(self.recharge_min, self.recharge_max,
                                        size=self.n)
            self.Re /= 1000. # Convert mm to m
        # Lognormal Distribution - Uniform in space
        elif self.groundwater__recharge_distribution == 'lognormal_uniform':
            assert (groundwater__recharge_mean != None), (
                'Input mean of the distribution!')
            assert (groundwater__recharge_standard_deviation != None), (
                'Input standard deviation of the distribution!')
            self.recharge_mean = groundwater__recharge_mean
            self.recharge_stdev = groundwater__recharge_standard_deviation
            self.mu_lognormal = np.log((self.recharge_mean**2)/np.sqrt(
                self.recharge_stdev**2 + self.recharge_mean**2))
            self.sigma_lognormal = np.sqrt(np.log((self.recharge_stdev**2)/(
                self.recharge_mean**2)+1))
            self.Re = np.random.lognormal(mean=self.mu_lognormal,
                                          sigma=self.sigma_lognormal,
                                          size=self.n)
            self.Re /= 1000. # Convert mm to m
        # Lognormal Distribution - Variable in space                                  
        elif self.groundwater__recharge_distribution == 'lognormal_spatial':
            assert (groundwater__recharge_mean.shape[0] == (
                self.grid.number_of_core_nodes)), (
                'Input array should be of the length of grid.number_of_nodes!')
            assert (groundwater__recharge_standard_deviation.shape[0] == (
                self.grid.number_of_core_nodes)), (
                'Input array should be of the length of grid.number_of_nodes!')
            self.recharge_mean = groundwater__recharge_mean
            self.recharge_stdev = groundwater__recharge_standard_deviation
        # Custom HSD inputs - Hydrologic Source Domain -> Model Domain
        elif self.groundwater__recharge_distribution == 'data_driven_spatial':
            self.HSD_dict = groundwater__recharge_HSD_inputs[0]
            self.HSD_id_dict = groundwater__recharge_HSD_inputs[1]
            self.fract_dict = groundwater__recharge_HSD_inputs[2]
            self._interpolate_HSD_dict()
            
        super(LandslideProbability, self).__init__(grid)

        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros('node', name, units=self._var_units[name])

        self._nodal_values = self.grid.at_node

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
        self.a = self.grid.at_node[
            'topographic__specific_contributing_area'][i]
        self.theta = self.grid.at_node['topographic__slope'][i]
        self.Tmode = self.grid.at_node['soil__transmissivity'][i]
        self.Cmode = self.grid.at_node['soil__mode_total_cohesion'][i]
        self.Cmin = self.grid.at_node['soil__minimum_total_cohesion'][i]
        self.Cmax = self.grid.at_node['soil__maximum_total_cohesion'][i]
        self.phi_mode = self.grid.at_node['soil__internal_friction_angle'][i]
        self.rho = self.grid.at_node['soil__density'][i]
        self.hs_mode = self.grid.at_node['soil__thickness'][i]

        # recharge distribution based on distribution type
        if self.groundwater__recharge_distribution == 'data_driven_spatial':
            self._calculate_HSD_recharge(i)
            self.Re /= 1000.0  # mm->m
        elif self.groundwater__recharge_distribution == 'lognormal_spatial':
            mu_lognormal = np.log((self.recharge_mean[i]**2)/np.sqrt(
                self.recharge_stdev[i]**2 + self.recharge_mean[i]**2))
            sigma_lognormal = np.sqrt(np.log((self.recharge_stdev[i]**2)/(
                self.recharge_mean[i]**2)+1))
            self.Re = np.random.lognormal(mean=mu_lognormal,
                                          sigma=sigma_lognormal,
                                          size=self.n)
            self.Re /= 1000. # Convert mm to m

        # Transmissivity (T)
        Tmin = self.Tmode-(0.3*self.Tmode)
        Tmax = self.Tmode+(0.1*self.Tmode)
        self.T = np.random.triangular(Tmin, self.Tmode, Tmax, size=self.n)
        # Cohesion
        # if provide fields of min and max C, uncomment 2 lines below
        #    Cmin = self.Cmode-0.3*self.Cmode
        #    Cmax = self.Cmode+0.3*self.Cmode
        self.C = np.random.triangular(self.Cmin, self.Cmode,
                                      self.Cmax, size=self.n)
        # phi - internal angle of friction provided in degrees
        phi_min = self.phi_mode-0.18*self.phi_mode
        phi_max = self.phi_mode+0.32*self.phi_mode
        self.phi = np.random.triangular(phi_min, self.phi_mode,
                                        phi_max, size=self.n)
        # soil thickness
        hs_min = min(0.005, self.hs_mode-0.3*self.hs_mode)
        hs_max = self.hs_mode+0.1*self.hs_mode
        self.hs = np.random.triangular(hs_min, self.hs_mode,
                                       hs_max, size=self.n)
        self.hs[self.hs <= 0.] = 0.0001

        # calculate Factor of Safety for n number of times
        # calculate components of FS equation
        self.C_dim = self.C/(self.hs*self.rho*self.g)  # dimensionless cohesion
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
        self.landslide__factor_of_safety_distribution: numpy.ndarray([
            self.grid.number_of_nodes, self.n], dtype=float)
            This is an output - distribution of factor-of-safety from
            Monte Carlo simulation (units='None')
        """

        # Create arrays for data with -9999 as default to store output
        self.mean_Relative_Wetness = -9999*np.ones(self.grid.number_of_nodes,
                                                   dtype='float')
        self.prob_fail = -9999*np.ones(
            self.grid.number_of_nodes, dtype='float')
        self.landslide__factor_of_safety_distribution = -9999*np.ones(
            [self.grid.number_of_nodes, self.n], dtype='float')
        # Run factor of safety Monte Carlo for all core nodes in domain
        # i refers to each core node id
        for i in self.grid.core_nodes:
            self.calculate_factor_of_safety(i)
            # Populate storage arrays with calculated values
            self.mean_Relative_Wetness[i] = self.soil__mean_relative_wetness
            self.prob_fail[i] = self.landslide__probability_of_failure
            self.landslide__factor_of_safety_distribution[i] = (
                self.FS_distribution)
            # stores FS values from last loop (node)
        # replace unrealistic values in arrays
        self.mean_Relative_Wetness[
            self.mean_Relative_Wetness < 0.] = 0.  # so can't be negative
        self.prob_fail[self.prob_fail < 0.] = 0.   # can't be negative
        # assign output fields to nodes
        self.grid.at_node['soil__mean_relative_wetness'] = (
            self.mean_Relative_Wetness)
        self.grid.at_node['landslide__probability_of_failure'] = self.prob_fail


    def seed_generator(self, seed=0):
        """Seed the random-number generator. This method will create the same
        sequence again by re-seeding with the same value (default value is
        zero). To create a sequence other than the default, assign non-zero
        value for seed.
        """
        np.random.seed(seed)


    def _interpolate_HSD_dict(self):
        HSD_dict = copy.deepcopy(self.HSD_dict)
        # First generate interpolated Re for each HSD grid
        Yrand = np.sort(np.random.rand(self.n))
        # n random numbers (0 to 1) in a column
        for vkey in HSD_dict.keys():
            if isinstance(HSD_dict[vkey], int):
                continue       # loop back up if value is integer, not array
            Re_temp = HSD_dict[vkey]	 # an array of years Re for 1 HSD grid
            Fx = ECDF(Re_temp)  # instantiate to get probabilities with Re
            Fx_ = Fx(Re_temp)    # probability array associated with Re data
            # interpolate function based on recharge data & probability
            f = interpolate.interp1d(Fx_, Re_temp, bounds_error=False,
                                     fill_value=min(Re_temp))
            # array of Re interpolated from Yrand probabilities (n count)
            Re_interpolated = f(Yrand)
            # replace values in HSD_dict with interpolated Re
            HSD_dict[vkey] = Re_interpolated

        self.interpolated_HSD_dict = HSD_dict


    def _calculate_HSD_recharge(self, i):
        store_Re = np.zeros(self.n)
        HSD_id_list = self.HSD_id_dict[i]
        fract_list = self.fract_dict[i]
        for j in range(0, len(HSD_id_list)):
            Re_temp = self.interpolated_HSD_dict[HSD_id_list[j]]
            fract_temp = fract_list[j]
            Re_adj = (Re_temp*fract_temp)
            store_Re = np.vstack((store_Re, np.array(Re_adj)))
        self.Re = np.sum(store_Re, 0)
