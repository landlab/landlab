#!/usr/env/python
"""Landlab component that simulates landslide probability of failure as well
as mean relative wetness and probability of saturation.

Relative wetness and factor-of-safety are based on the infinite slope
stability model driven by topographic and soils inputs and recharge provided
by user as inputs to the component. For each node, component simulates mean
relative wetness as well as the probability of saturation based on Monte Carlo
simulation of relative wetness where the probability is the number of
iterations with relative wetness >= 1.0 divided by the number of iterations.
Probability of failure for each node is also simulated in the Monte Carlo
simulation as the number of iterations with factor-of-safety <= 1.0
divided by the number of iterations.

.. codeauthor:: R.Strauch, E.Istanbulluoglu, & S.S.Nudurupati
University of Washington

Ref 1: Strauch et. al. 2017, 'A hydro-climatological approach to predicting
regional landslide probability using Landlab, Earth Surface Dynamics, In prep.

Ref 2: 'The Landlab LandslideProbability Component User Manual' @
https://github.com/RondaStrauch/pub_strauch_etal_esurf/blob/master/LandslideComponentUsersManual.docx

Created on Thu Aug 20, 2015
Last edit June 7, 2017
"""

import copy

import numpy as np
import scipy.constants
from scipy import interpolate
from statsmodels.distributions.empirical_distribution import ECDF

from landlab import Component
from landlab.utils.decorators import use_file_name_or_kwds


class LandslideProbability(Component):
    """Landslide probability component using the infinite slope stability
    model.

    Landlab component designed to calculate probability of failure at
    each grid node based on the infinite slope stability model
    stability index (Factor of Safety).

    The driving force for failure is provided by the user in the form of
    groundwater recharge; four options for providing recharge are supported.
    The model uses topographic and soil characteristics provided as input
    by the user.

    The main method of the LandslideProbability class is
    calculate_landslide_probability(), which calculates the mean soil relative
    wetness, probability of soil saturation, and probability of failure at
    each node based on a Monte Carlo simulation.

    **Usage:**

    Option 1 - Uniform recharge::

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__recharge_distribution='uniform',
                             groundwater__recharge_min_value=5.,
                             groundwater__recharge_max_value=121.)

    Option 2 - Lognormal recharge::

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__recharge_distribution='lognormal',
                             groundwater__recharge_mean=30.,
                             groundwater__recharge_standard_deviation=0.25)

    Option 3 - Lognormal_spatial recharge::

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__recharge_distribution='lognormal_spatial',
                             groundwater__recharge_mean=np.random.randint(20, 120, grid_size),
                             groundwater__recharge_standard_deviation=np.random.rand(grid_size))

    Option 4 - Data_driven_spatial recharge::

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__recharge_distribution='data_driven_spatial',
                             groundwater__recharge_HSD_inputs=[HSD_dict, HSD_id_dict,
                             fract_dict])

    Examples
    ----------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components.landslides import LandslideProbability
    >>> import numpy as np

    Create a grid on which to calculate landslide probability.

    >>> grid = RasterModelGrid((5, 4), xy_spacing=(0.2, 0.2))

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
     'soil__saturated_hydraulic_conductivity',
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
      mode rate of water transmitted through a unit width of saturated
      soil - either provided or calculated with Ksat and soil depth
    units: m2/day
    at: node
    intent: in

    Additional required fields for component.

    >>> scatter_dat = np.random.randint(1, 10, grid.number_of_nodes)
    >>> grid.at_node['topographic__specific_contributing_area'] = np.sort(
    ...      np.random.randint(30, 900, grid.number_of_nodes))
    >>> grid.at_node['soil__transmissivity'] = np.sort(
    ...      np.random.randint(5, 20, grid.number_of_nodes), -1)
    >>> grid.at_node['soil__saturated_hydraulic_conductivity'] = np.sort(
    ...      np.random.randint(2, 10, grid.number_of_nodes), -1)
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

    >>> ls_prob = LandslideProbability(grid)
    >>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
    True

    Run the *calculate_landslide_probability* method to update output
    variables with grid

    >>> ls_prob.calculate_landslide_probability()

    Check the output variable names.

    >>> sorted(ls_prob.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['landslide__probability_of_failure',
     'soil__mean_relative_wetness',
     'soil__probability_of_saturation']

    Check the output from the component, including array at one node.

    >>> np.allclose(grid.at_node['landslide__probability_of_failure'], 0.)
    False
    >>> core_nodes = ls_prob.grid.core_nodes
    """

    # component name
    _name = "Landslide Probability"
    __version__ = "1.0"
    # component requires these values to do its calculation, get from driver
    _input_var_names = (
        "topographic__specific_contributing_area",
        "topographic__slope",
        "soil__transmissivity",
        "soil__saturated_hydraulic_conductivity",
        "soil__mode_total_cohesion",
        "soil__minimum_total_cohesion",
        "soil__maximum_total_cohesion",
        "soil__internal_friction_angle",
        "soil__density",
        "soil__thickness",
    )

    #  component creates these output values
    _output_var_names = (
        "soil__mean_relative_wetness",
        "landslide__probability_of_failure",
        "soil__probability_of_saturation",
    )

    # units for each parameter and output
    _var_units = {
        "topographic__specific_contributing_area": "m",
        "topographic__slope": "tan theta",
        "soil__transmissivity": "m2/day",
        "soil__saturated_hydraulic_conductivity": "m/day",
        "soil__mode_total_cohesion": "Pa or kg/m-s2",
        "soil__minimum_total_cohesion": "Pa or kg/m-s2",
        "soil__maximum_total_cohesion": "Pa or kg/m-s2",
        "soil__internal_friction_angle": "degrees",
        "soil__density": "kg/m3",
        "soil__thickness": "m",
        "soil__mean_relative_wetness": "None",
        "landslide__probability_of_failure": "None",
        "soil__probability_of_saturation": "None",
    }

    # grid centering of each field and variable
    _var_mapping = {
        "topographic__specific_contributing_area": "node",
        "topographic__slope": "node",
        "soil__transmissivity": "node",
        "soil__saturated_hydraulic_conductivity": "node",
        "soil__mode_total_cohesion": "node",
        "soil__minimum_total_cohesion": "node",
        "soil__maximum_total_cohesion": "node",
        "soil__internal_friction_angle": "node",
        "soil__density": "node",
        "soil__thickness": "node",
        "soil__mean_relative_wetness": "node",
        "landslide__probability_of_failure": "node",
        "soil__probability_of_saturation": "node",
    }

    # short description of each field
    _var_doc = {
        "topographic__specific_contributing_area": (
            "specific contributing (upslope area/cell face )" + " that drains to node"
        ),
        "topographic__slope": "slope of surface at node represented by tan theta",
        "soil__transmissivity": (
            "mode rate of water transmitted"
            + " through a unit width of saturated soil - "
            + "either provided or calculated with Ksat "
            + "and soil depth"
        ),
        "soil__saturated_hydraulic_conductivity": (
            "mode rate of water transmitted"
            + " through soil - provided if transmissivity "
            + "is NOT provided to calculate tranmissivity "
            + " with soil depth"
        ),
        "soil__mode_total_cohesion": "mode of combined root and soil cohesion at node",
        "soil__minimum_total_cohesion": "minimum of combined root and soil cohesion at node",
        "soil__maximum_total_cohesion": "maximum of combined root and soil cohesion at node",
        "soil__internal_friction_angle": (
            "critical angle just before failure" + " due to friction between particles"
        ),
        "soil__density": "wet bulk density of soil",
        "soil__thickness": "soil depth to restrictive layer",
        "soil__mean_relative_wetness": (
            "Indicator of soil wetness;"
            + " relative depth perched water table"
            + " within the soil layer"
        ),
        "landslide__probability_of_failure": (
            "number of times FS is <=1 out of number of" + " iterations user selected"
        ),
        "soil__probability_of_saturation": (
            "number of times relative wetness is >=1 out of"
            + " number of iterations user selected"
        ),
    }

    # Run Component
    @use_file_name_or_kwds
    def __init__(
        self,
        grid,
        number_of_iterations=250,
        groundwater__recharge_distribution="uniform",
        groundwater__recharge_min_value=20.0,
        groundwater__recharge_max_value=120.0,
        groundwater__recharge_mean=None,
        groundwater__recharge_standard_deviation=None,
        groundwater__recharge_HSD_inputs=[],
        seed=0,
        **kwds
    ):
        """
        Parameters
        ----------
        grid: RasterModelGrid
            A raster grid.
        number_of_iterations: int, optional
            Number of iterations to run Monte Carlo simulation (default=250).
        groundwater__recharge_distribution: str, optional
            single word indicating recharge distribution, either 'uniform',
            'lognormal', 'lognormal_spatial,' or 'data_driven_spatial'.
             (default='uniform')
        groundwater__recharge_min_value: float, optional (mm/d)
            minium groundwater recharge for 'uniform' (default=20.)
        groundwater__recharge_max_value: float, optional (mm/d)
            maximum groundwater recharge for 'uniform' (default=120.)
        groundwater__recharge_mean: float, optional (mm/d)
            mean grounwater recharge for 'lognormal'
            and 'lognormal_spatial' (default=None)
        groundwater__recharge_standard_deviation: float, optional (mm/d)
            standard deviation of grounwater recharge for 'lognormal'
            and 'lognormal_spatial' (default=None)
        groundwater__recharge_HSD_inputs: list, optional
            list of 3 dictionaries in order (default=[]) - HSD_dict
            {Hydrologic Source Domain (HSD) keys: recharge numpy array values},
            {node IDs keys: list of HSD_Id values}, HSD_fractions {node IDS
            keys: list of HSD fractions values} (none)
            Note: this input method is a very specific one, and to use this method,
            one has to refer Ref 1 & Ref 2 mentioned above, as this set of
            inputs require rigorous pre-processing of data.
        g: float, optional (m/sec^2)
            acceleration due to gravity.
        seed: int, optional
            seed for random number generation. if seed is assigned any value
            other than the default value of zero, it will create different
            sequence. To create a certain sequence repititively, use the same
            value as input for seed.
        """
        # Initialize seeded random number generation
        self._seed_generator(seed)

        super(LandslideProbability, self).__init__(grid)

        # Store grid and parameters and do unit conversions
        self.n = int(number_of_iterations)
        self._g = kwds.get("g", scipy.constants.g)
        self.groundwater__recharge_distribution = groundwater__recharge_distribution
        # Following code will deal with the input distribution and associated
        # parameters
        # Uniform distribution
        if self.groundwater__recharge_distribution == "uniform":
            self._recharge_min = groundwater__recharge_min_value
            self._recharge_max = groundwater__recharge_max_value
            self._Re = np.random.uniform(
                self._recharge_min, self._recharge_max, size=self.n
            )
            self._Re /= 1000.0  # Convert mm to m
        # Lognormal Distribution - Uniform in space
        elif self.groundwater__recharge_distribution == "lognormal":
            assert (
                groundwater__recharge_mean is not None
            ), "Input mean of the distribution!"
            assert (
                groundwater__recharge_standard_deviation is not None
            ), "Input standard deviation of the distribution!"
            self._recharge_mean = groundwater__recharge_mean
            self._recharge_stdev = groundwater__recharge_standard_deviation
            self._mu_lognormal = np.log(
                (self._recharge_mean ** 2)
                / np.sqrt(self._recharge_stdev ** 2 + self._recharge_mean ** 2)
            )
            self._sigma_lognormal = np.sqrt(
                np.log((self._recharge_stdev ** 2) / (self._recharge_mean ** 2) + 1)
            )
            self._Re = np.random.lognormal(
                self._mu_lognormal, self._sigma_lognormal, self.n
            )
            self._Re /= 1000.0  # Convert mm to m
        # Lognormal Distribution - Variable in space
        elif self.groundwater__recharge_distribution == "lognormal_spatial":
            assert groundwater__recharge_mean.shape[0] == (
                self.grid.number_of_nodes
            ), "Input array should be of the length of grid.number_of_nodes!"
            assert groundwater__recharge_standard_deviation.shape[0] == (
                self.grid.number_of_nodes
            ), "Input array should be of the length of grid.number_of_nodes!"
            self._recharge_mean = groundwater__recharge_mean
            self._recharge_stdev = groundwater__recharge_standard_deviation
        # Custom HSD inputs - Hydrologic Source Domain -> Model Domain
        elif self.groundwater__recharge_distribution == "data_driven_spatial":
            self._HSD_dict = groundwater__recharge_HSD_inputs[0]
            self._HSD_id_dict = groundwater__recharge_HSD_inputs[1]
            self._fract_dict = groundwater__recharge_HSD_inputs[2]
            self._interpolate_HSD_dict()

        # Check if all input fields are initialized
        for name in self._input_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros(
                    name, at=self._var_mapping[name], units=self._var_units[name]
                )

        # Check if all output fields are initialized
        for name in self._output_var_names:
            if name not in self.grid.at_node:
                self.grid.add_zeros(
                    name, at=self._var_mapping[name], units=self._var_units[name]
                )

        # Create a switch to imply whether Ksat is provided.
        if np.all(self.grid.at_node["soil__saturated_hydraulic_conductivity"] == 0):
            self.Ksat_provided = 0  # False
        else:
            self.Ksat_provided = 1  # True

        self._nodal_values = self.grid.at_node

    def calculate_factor_of_safety(self, i):
        """Method to calculate factor of safety.

        Method calculates factor-of-safety stability index by using
        node specific parameters, creating distributions of these parameters,
        and calculating the index by sampling these distributions 'n' times.

        The index is calculated from the 'infinite slope stabilty
        factor-of-safety equation' in the format of Pack RT, Tarboton DG,
        and Goodwin CN (1998),The SINMAP approach to terrain stability mapping.

        Parameters
        ----------
        i: int
            index of core node ID.
        """
        # generate distributions to sample from to provide input parameters
        # currently triangle distribution using mode, min, & max
        self._a = np.float32(
            self.grid.at_node["topographic__specific_contributing_area"][i]
        )
        self._theta = np.float32(self.grid.at_node["topographic__slope"][i])
        self._Tmode = np.float32(self.grid.at_node["soil__transmissivity"][i])
        self._Ksatmode = np.float32(
            self.grid.at_node["soil__saturated_hydraulic_conductivity"][i]
        )
        self._Cmode = np.float32(self.grid.at_node["soil__mode_total_cohesion"][i])
        self._Cmin = np.float32(self.grid.at_node["soil__minimum_total_cohesion"][i])
        self._Cmax = np.float32(self.grid.at_node["soil__maximum_total_cohesion"][i])
        self._phi_mode = np.float32(
            self.grid.at_node["soil__internal_friction_angle"][i]
        )
        self._rho = np.float32(self.grid.at_node["soil__density"][i])
        self._hs_mode = np.float32(self.grid.at_node["soil__thickness"][i])

        # recharge distribution based on distribution type
        if self.groundwater__recharge_distribution == "data_driven_spatial":
            self._calculate_HSD_recharge(i)
            self._Re /= 1000.0  # mm->m
        elif self.groundwater__recharge_distribution == "lognormal_spatial":
            mu_lognormal = np.log(
                (self._recharge_mean[i] ** 2)
                / np.sqrt(self._recharge_stdev[i] ** 2 + self._recharge_mean[i] ** 2)
            )
            sigma_lognormal = np.sqrt(
                np.log(
                    (self._recharge_stdev[i] ** 2) / (self._recharge_mean[i] ** 2) + 1
                )
            )
            self._Re = np.random.lognormal(mu_lognormal, sigma_lognormal, self.n)
            self._Re /= 1000.0  # Convert mm to m

        # Cohesion
        # if don't provide fields of min and max C, uncomment 2 lines below
        #    Cmin = self._Cmode-0.3*self._Cmode
        #    Cmax = self._Cmode+0.3*self._Cmode
        self._C = np.random.triangular(self._Cmin, self._Cmode, self._Cmax, size=self.n)

        # phi - internal angle of friction provided in degrees
        phi_min = self._phi_mode - 0.18 * self._phi_mode
        phi_max = self._phi_mode + 0.32 * self._phi_mode
        self._phi = np.random.triangular(phi_min, self._phi_mode, phi_max, size=self.n)
        # soil thickness
        # hs_min = min(0.005, self._hs_mode-0.3*self._hs_mode) # Alternative
        hs_min = self._hs_mode - 0.3 * self._hs_mode
        hs_max = self._hs_mode + 0.1 * self._hs_mode
        self._hs = np.random.triangular(hs_min, self._hs_mode, hs_max, size=self.n)
        self._hs[self._hs <= 0.0] = 0.005
        if self.Ksat_provided:
            # Hydraulic conductivity (Ksat)
            Ksatmin = self._Ksatmode - (0.3 * self._Ksatmode)
            Ksatmax = self._Ksatmode + (0.1 * self._Ksatmode)
            self._Ksat = np.random.triangular(
                Ksatmin, self._Ksatmode, Ksatmax, size=self.n
            )
            self._T = self._Ksat * self._hs
        else:
            # Transmissivity (T)
            Tmin = self._Tmode - (0.3 * self._Tmode)
            Tmax = self._Tmode + (0.1 * self._Tmode)
            self._T = np.random.triangular(Tmin, self._Tmode, Tmax, size=self.n)

        # calculate Factor of Safety for n number of times
        # calculate components of FS equation
        self._C_dim = self._C / (
            self._hs * self._rho * self._g
        )  # dimensionless cohesion
        self._rel_wetness = ((self._Re) / self._T) * (
            self._a / np.sin(np.arctan(self._theta))
        )  # relative wetness
        # calculate probability of saturation
        countr = 0
        for val in self._rel_wetness:  # find how many RW values >= 1
            if val >= 1.0:
                countr = countr + 1  # number with RW values (>=1)
        # probability: No. high RW values/total No. of values (n)
        self._soil__probability_of_saturation = np.float32(countr) / self.n
        # Maximum Rel_wetness = 1.0
        np.place(self._rel_wetness, self._rel_wetness > 1, 1.0)
        self._soil__mean_relative_wetness = np.mean(self._rel_wetness)
        Y = np.tan(np.radians(self._phi)) * (1 - (self._rel_wetness * 0.5))
        # convert from degrees; 0.5 = water to soil density ratio
        # calculate Factor-of-safety
        self._FS = (self._C_dim / np.sin(np.arctan(self._theta))) + (
            np.cos(np.arctan(self._theta)) * (Y / np.sin(np.arctan(self._theta)))
        )
        count = 0
        for val in self._FS:  # find how many FS values <= 1
            if val <= 1.0:
                count = count + 1  # number with unstable FS values (<=1)
        # probability: No. unstable values/total No. of values (n)
        self._landslide__probability_of_failure = np.float32(count) / self.n

    def calculate_landslide_probability(self, **kwds):
        """Main method of Landslide Probability class.

        Method creates arrays for output variables then loops through all
        the core nodes to run the method 'calculate_factor_of_safety.'
        Output parameters probability of failure, mean relative wetness,
        and probability of saturation are assigned as fields to nodes.
        """
        # Create arrays for data with -9999 as default to store output
        self.mean_Relative_Wetness = np.full(self.grid.number_of_nodes, -9999.0)
        self.prob_fail = np.full(self.grid.number_of_nodes, -9999.0)
        self.prob_sat = np.full(self.grid.number_of_nodes, -9999.0)
        # Run factor of safety Monte Carlo for all core nodes in domain
        # i refers to each core node id
        for i in self.grid.core_nodes:
            self.calculate_factor_of_safety(i)
            # Populate storage arrays with calculated values
            self.mean_Relative_Wetness[i] = self._soil__mean_relative_wetness
            self.prob_fail[i] = self._landslide__probability_of_failure
            self.prob_sat[i] = self._soil__probability_of_saturation
        # Values can't be negative
        self.mean_Relative_Wetness[self.mean_Relative_Wetness < 0.0] = 0.0
        self.prob_fail[self.prob_fail < 0.0] = 0.0
        # assign output fields to nodes
        self.grid.at_node["soil__mean_relative_wetness"] = self.mean_Relative_Wetness
        self.grid.at_node["landslide__probability_of_failure"] = self.prob_fail
        self.grid.at_node["soil__probability_of_saturation"] = self.prob_sat

    def _seed_generator(self, seed=0):
        """Method to initiate random seed.

        Seed the random-number generator. This method will create the same
        sequence again by re-seeding with the same value (default value is
        zero). To create a sequence other than the default, assign non-zero
        value for seed.
        """
        np.random.seed(seed)

    def _interpolate_HSD_dict(self):
        """Method to extrapolate input data.

        This method uses a non-parametric approach to expand the input
        recharge array to the length of number of iterations. Output is
        a new dictionary of interpolated recharge for each HSD id.
        """
        HSD_dict = copy.deepcopy(self._HSD_dict)
        # First generate interpolated Re for each HSD grid
        Yrand = np.sort(np.random.rand(self.n))
        # n random numbers (0 to 1) in a column
        for vkey in HSD_dict.keys():
            if isinstance(HSD_dict[vkey], int):
                continue  # loop back up if value is integer (e.g. -9999)
            Re_temp = HSD_dict[vkey]  # an array of annual Re for 1 HSD grid
            Fx = ECDF(Re_temp)  # instantiate to get probabilities with Re
            Fx_ = Fx(Re_temp)  # probability array associated with Re data
            # interpolate function based on recharge data & probability
            f = interpolate.interp1d(
                Fx_, Re_temp, bounds_error=False, fill_value=min(Re_temp)
            )
            # array of Re interpolated from Yrand probabilities (n count)
            Re_interpolated = f(Yrand)
            # replace values in HSD_dict with interpolated Re
            HSD_dict[vkey] = Re_interpolated

        self._interpolated_HSD_dict = HSD_dict

    def _calculate_HSD_recharge(self, i):
        """Method to calculate recharge based on upstream fractions.

        This method calculates the resultant recharge at node i of the
        model domain, using recharge of contributing HSD ids and the areal
        fractions of upstream contributing HSD ids. Output is a numpy array
        of recharge at node i.
        """
        store_Re = np.zeros(self.n)
        HSD_id_list = self._HSD_id_dict[i]
        fract_list = self._fract_dict[i]
        for j in range(0, len(HSD_id_list)):
            Re_temp = self._interpolated_HSD_dict[HSD_id_list[j]]
            fract_temp = fract_list[j]
            Re_adj = Re_temp * fract_temp
            store_Re = np.vstack((store_Re, np.array(Re_adj)))
        self._Re = np.sum(store_Re, 0)
