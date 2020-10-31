#!/usr/env/python
"""Landlab component that simulates landslide probability of failure as well as
mean relative wetness and probability of saturation.

Relative wetness and factor-of-safety are based on the infinite slope
stability model driven by topographic and soils inputs and recharge provided
by user as inputs to the component. For each node, component simulates mean
relative wetness as well as the probability of saturation based on Monte Carlo
simulation of relative wetness where the probability is the number of
iterations with relative wetness >= 1.0 divided by the number of iterations.
Probability of failure for each node is also simulated in the Monte Carlo
simulation as the number of iterations with factor-of-safety <= 1.0
divided by the number of iterations.

.. codeauthor:: R.Strauch, C.Bandaragoda, E.Istanbulluoglu, & S.S.Nudurupati

University of Washington

Ref 1: Strauch et. al. 2017, 'A hydro-climatological approach to predicting
regional landslide probability using Landlab, Earth Surface Dynamics, In prep.

Ref 2: 'The Landlab LandslideProbability Component User Manual' @
https://github.com/RondaStrauch/pub_strauch_etal_esurf/blob/master/LandslideComponentUsersManual.pdf

Created on Thu Aug 20, 2015
Last feature edits June 7, 2017
Landlab v2 and PEP8 updates by mcflugen and kbarnhart 2018-2019
Component Changes: April 28, 2020
Add depth to water table as an optional hydrology input (alternative to recharge)

"""

import copy

import numpy as np
import pandas as pd
import scipy.constants
from scipy import interpolate
from statsmodels.distributions.empirical_distribution import ECDF

from landlab import Component


class LandslideProbability(Component):
    """Landslide probability component using the infinite slope stability
    model.

    Landlab component designed to calculate probability of failure at
    each grid node based on the infinite slope stability model
    stability index (Factor of Safety).

    The driving force for failure is provided by the user in the form of
    groundwater recharge OR depth to groundwater. Four options for providing
    recharge plus four options of depth to groundwater are supported.
    The model uses topographic and soil characteristics provided as input
    by the user.

    The main method of the LandslideProbability class is
    `calculate_landslide_probability()`, which calculates the mean soil
    relative wetness, mean water table depth, mean recharge, probability of soil saturation,
    and probability of failure at each node based on a Monte Carlo simulation.

    **Usage:**

    Option 1 - Uniform recharge

    .. code-block:: python

        LandslideProbability(grid,
                             number_of_iterations=250,
                             groundwater__recharge_distribution='uniform',
                             groundwater__recharge_min_value=5.,
                             groundwater__recharge_max_value=120.)

    Option 2 - Lognormal recharge

    .. code-block:: python

        LandslideProbability(grid,
                             number_of_iterations=250,
                             groundwater__recharge_distribution='lognormal',
                             groundwater__recharge_mean=30.,
                             groundwater__recharge_standard_deviation=0.25)

    Option 3 - Lognormal_spatial recharge

    .. code-block:: python

        LandslideProbability(grid,
                             number_of_iterations=250,
                             groundwater__recharge_distribution='lognormal_spatial',
                             groundwater__recharge_mean=np.random.randint(20, 120, grid_size),
                             groundwater__recharge_standard_deviation=np.random.rand(grid_size))

    Option 4 - Data_driven_spatial recharge

    .. code-block:: python

        LandslideProbability(grid,
                             number_of_iterations=250,
                             groundwater__recharge_distribution='data_driven_spatial',
                             groundwater__recharge_HSD_inputs=[HSD_dict,
                                                               HSD_id_dict,
                                                               fract_dict])
    Option 5 - Uniform depth to groundwater

    .. code-block:: python

        LandslideProbability(grid,
                             number_of_iterations=250,
                             groundwater__depth_distribution='uniform',
                             groundwater__depth_min_value=0.01,
                             groundwater__depth_max_value=1.5)

    Option 6 - Lognormal depth to groundwater

    .. code-block:: python

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__depth_distribution='lognormal',
                             groundwater__depth_mean=0.5.,
                             groundwater__depth_standard_deviation=0.1)

    Option 7 - Lognormal_spatial depth to groundwater

    .. code-block:: python

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__depth_distribution='lognormal_spatial',
                             groundwater__depth_mean=np.random.randint(0, 2, grid_size),
                             groundwater__depth_standard_deviation=np.random.rand(grid_size))

    Option 8 - Data_driven_spatial depth to groundwater

    .. code-block:: python

        LandslideProbability(grid, number_of_iterations=250,
                             groundwater__depth_distribution='data_driven_spatial',
                             groundwater__depth_HSD_inputs={HSD_dict})

    Examples
    --------
    >>> from landlab import RasterModelGrid
    >>> from landlab.components import LandslideProbability
    >>> import landlab.grid.landslide_unitgrid as landslide_unitgrid
    >>> import numpy as np

    Create a grid on which to calculate landslide probability.

    >>> grid = landslide_unitgrid.build_grid_unitarea()

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

    >>> LandslideProbability.var_help('soil__mean_watertable_depth')
    name: soil__mean_watertable_depth
    description:
        Mean depth to water table from surface to perched water table within
        the soil layer
        units: m
        unit agnostic: False
        at: node
        intent: out

    >>> LandslideProbability.var_help('soil__mean_recharge')
    name: soil__mean_recharge
    description:
        Mean recharge to the soil layer
        units: mm/day
        unit agnostic: False
        at: node
        intent: out


    Instantiate a new grid with default parameters for each model instance


    >>> grid_r1 = landslide_unitgrid.build_grid_unitarea()

    >>> grid_d1 = landslide_unitgrid.build_grid_unitarea()

    Instantiate the 'LandslideProbability' component to work on each grid,
    and run it.
    >>> LS_prob1_r = LandslideProbability(grid_r1,groundwater__recharge_distribution='uniform')

    >>> LS_prob1_d = LandslideProbability(grid_d1, groundwater__depth_distribution='uniform')

    Run the *calculate_landslide_probability* method to update output
    variables with grid

    >>> LS_prob1_r.calculate_landslide_probability()
    >>> LS_prob1_d.calculate_landslide_probability()

    Check the output variable names.

    >>> sorted(LS_prob1_d.output_var_names) # doctest: +NORMALIZE_WHITESPACE
    ['landslide__probability_of_failure', 'soil__mean_recharge', 'soil__mean_relative_wetness', 'soil__mean_watertable_depth', 'soil__probability_of_saturation']

    References
    ----------
    **Required Software Citation(s) Specific to this Component**

    Strauch, R., Istanbulluoglu, E., Nudurupati, S., Bandaragoda, C.,
    Gasparini, N., Tucker, G. (2018). A hydroclimatological approach to
    predicting regional landslide probability using Landlab Earth Surface
    Dynamics  6(1), 49-75. https://dx.doi.org/10.5194/esurf-6-49-2018

    **Additional References**

    Publication Pending

    """

    # component name
    _name = "Landslide Probability"

    _unit_agnostic = False

    __version__ = "2.0"

    _cite_as = """
    @article{strauch2018hydroclimatological,
      author = {Strauch, Ronda and Istanbulluoglu, Erkan and Nudurupati,
      Sai Siddhartha and Bandaragoda, Christina and Gasparini, Nicole M and
      Tucker, Gregory E},
      title = {{A hydroclimatological approach to predicting regional landslide
      probability using Landlab}},
      issn = {2196-6311},
      doi = {10.5194/esurf-6-49-2018},
      pages = {49--75},
      number = {1},
      volume = {6},
      journal = {Earth Surface Dynamics},
      year = {2018}
    }
    """
    _info = {
        "landslide__probability_of_failure": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "node",
            "doc": "number of times FS is <=1 out of number of iterations user selected",
        },
        "soil__density": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "kg/m3",
            "mapping": "node",
            "doc": "wet bulk density of soil",
        },
        "soil__internal_friction_angle": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "degrees",
            "mapping": "node",
            "doc": "critical angle just before failure due to friction between particles",
        },
        "soil__maximum_total_cohesion": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "Pa or kg/m-s2",
            "mapping": "node",
            "doc": "maximum of combined root and soil cohesion at node",
        },
        "soil__mean_relative_wetness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "node",
            "doc": "Indicator of soil wetness; relative depth perched water table within the soil layer",
        },
        "soil__mean_watertable_depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Mean depth to water table from surface to perched water table within the soil layer",
        },
        "soil__mean_recharge": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "mm/day",
            "mapping": "node",
            "doc": "Mean recharge to the soil layer",
        },
        "soil__minimum_total_cohesion": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "Pa or kg/m-s2",
            "mapping": "node",
            "doc": "minimum of combined root and soil cohesion at node",
        },
        "soil__mode_total_cohesion": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "Pa or kg/m-s2",
            "mapping": "node",
            "doc": "mode of combined root and soil cohesion at node",
        },
        "soil__probability_of_saturation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "None",
            "mapping": "node",
            "doc": "number of times relative wetness is >=1 out of number of iterations user selected",
        },
        "soil__saturated_hydraulic_conductivity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m/day",
            "mapping": "node",
            "doc": "mode rate of water transmitted through soil. If transmissivity is NOT provided, the component calculates transmissivity using Ksat and soil depth",
        },
        "soil__thickness": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "soil depth to restrictive layer",
        },
        "soil__transmissivity": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "m2/day",
            "mapping": "node",
            "doc": "mode rate of water transmitted through a unit width of saturated soil - either provided or calculated with Ksat and soil depth",
        },
        "topographic__slope": {
            "dtype": float,
            "intent": "in",
            "optional": False,
            "units": "tan theta",
            "mapping": "node",
            "doc": "gradient of the ground surface",
        },
        "topographic__specific_contributing_area": {
            "dtype": int,
            "intent": "in",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "specific contributing (upslope area/cell face ) that drains to node",
        },
    }

    def __init__(
        self,
        grid,
        number_of_iterations=250,
        g=scipy.constants.g,
        groundwater__recharge_distribution=None,
        groundwater__recharge_min_value=20.0,
        groundwater__recharge_max_value=120.0,
        groundwater__recharge_mean=None,
        groundwater__recharge_standard_deviation=None,
        groundwater__recharge_HSD_inputs=[],
        groundwater__depth_distribution=None,
        groundwater__depth_min_value=0.01,
        groundwater__depth_max_value=3.0,
        groundwater__depth_mean=None,
        groundwater__depth_standard_deviation=None,
        groundwater__depth_HSD_inputs={},
        seed=0,
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
            minimum groundwater recharge for 'uniform' (default=20.0)
        groundwater__recharge_max_value: float, optional (mm/d)
            maximum groundwater recharge for 'uniform' (default=120.0)
        groundwater__recharge_mean: float, optional (mm/d)
            mean groundwater recharge for 'lognormal'
            and 'lognormal_spatial' (default=None)
        groundwater__recharge_standard_deviation: float, optional (mm/d)
            standard deviation of groundwater recharge for 'lognormal'
            and 'lognormal_spatial' (default=None)
        groundwater__recharge_HSD_inputs: list, optional
            list of 3 dictionaries in order (default=[]) - HSD_dict
            {Hydrologic Source Domain (HSD) keys: recharge numpy array values},
            {node IDs keys: list of HSD_Id values}, HSD_fractions {node IDS
            keys: list of HSD fractions values} (none)
            Note: this input method is a very specific one, and to use this method,
            one has to refer Ref 1 & Ref 2 mentioned above, as this set of
            inputs require rigorous pre-processing of data.
        groundwater__depth_distribution: str, optional
            single word indicating depth to water table distribution, either
            'uniform', 'lognormal', 'lognormal_spatial,' or
            'data_driven_spatial'.
             (default=None)
        groundwater__depth_min_value: float, optional (m)
            minimum groundwater depth to water table for 'uniform' (default=0.01)
        groundwater__depth_max_value: float, optional (m)
            maximum groundwater depth for 'uniform' (default=2.)
        groundwater__depth_mean: float, optional (m)
            mean groundwater depth to water table for 'lognormal'
            and 'lognormal_spatial' (default=None)
        groundwater__depth_standard_deviation: float, optional (m)
            standard deviation of groundwater depth to water table for
            'lognormal' and 'lognormal_spatial' (default=None)
        groundwater__depth_HSD_inputs: dictionary, optional
            one dictionary (default={}) - HSD_dict
            {Hydrologic Source Domain (HSD) Grid equal size to Landlab
            node ID keys: groundwater depth distribution numpy array values}
        g: float, optional (m/sec^2)
            acceleration due to gravity.
        seed: int, optional
            seed for random number generation. if seed is assigned any value
            other than the default value of zero, it will create different
            sequence. To create a certain sequence repetitively, use the same
            value as input for seed.
        """
        # Initialize seeded random number generation
        self._seed_generator(seed)

        super().__init__(grid)

        # Store parameters and do unit conversions
        self._n = int(number_of_iterations)
        self._g = g
        self._groundwater__recharge_distribution = groundwater__recharge_distribution
        self._groundwater__depth_distribution = groundwater__depth_distribution

        # Following code will deal with the input distribution and associated
        # parameters for recharge hydrologic forcing
        # Recharge Uniform distribution
        if self._groundwater__recharge_distribution == "uniform":
            self._recharge_min = groundwater__recharge_min_value
            self._recharge_max = groundwater__recharge_max_value
            self._Re = np.random.uniform(
                self._recharge_min, self._recharge_max, size=self._n
            )
            self._Re /= 1000.0  # Convert mm to m
        # Recharge Lognormal Distribution - Uniform in space
        elif self._groundwater__recharge_distribution == "lognormal":
            assert (
                groundwater__recharge_mean is not None
            ), "Input groundwater__recharge_mean and try again!"
            assert (
                groundwater__recharge_standard_deviation is not None
            ), "Input missing standard deviation of the distribution!"
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
                self._mu_lognormal, self._sigma_lognormal, self._n
            )
            self._Re /= 1000.0  # Convert mm to m
        # Recharge Lognormal Distribution - Variable in space
        elif self._groundwater__recharge_distribution == "lognormal_spatial":
            assert groundwater__recharge_mean.shape[0] == (
                self._grid.number_of_nodes
            ), "Input array should be of the length of grid.number_of_nodes!"
            assert groundwater__recharge_standard_deviation.shape[0] == (
                self._grid.number_of_nodes
            ), "Input array should be of the length of grid.number_of_nodes!"
            self._recharge_mean = groundwater__recharge_mean
            self._recharge_stdev = groundwater__recharge_standard_deviation
        # Recharge Custom HSD inputs - Hydrologic Source Domain -> Model Domain
        elif self._groundwater__recharge_distribution == "data_driven_spatial":
            self._HSD_dict = groundwater__recharge_HSD_inputs[0]
            self._HSD_id_dict = groundwater__recharge_HSD_inputs[1]
            self._fract_dict = groundwater__recharge_HSD_inputs[2]
            self._interpolate_HSD_dict()

        # Following code will deal with the input distribution and associated
        # parameters for depth to water table hydrologic forcing
        # Depth to water table - Uniform distribution
        if self._groundwater__depth_distribution == "uniform":
            self._depth_min = groundwater__depth_min_value
            self._depth_max = groundwater__depth_max_value
            self._De = np.random.uniform(self._depth_min, self._depth_max, size=self._n)
        # Depth to water table - Lognormal Distribution - Uniform in space
        elif self._groundwater__depth_distribution == "lognormal":
            assert (
                groundwater__depth_mean is not None
            ), "Input mean of the distribution!"
            assert (
                groundwater__depth_standard_deviation is not None
            ), "Input standard deviation of the distribution!"
            self._depth_mean = groundwater__depth_mean
            self._depth_stdev = groundwater__depth_standard_deviation
            self._mu_lognormal = np.log(
                (self._depth_mean ** 2)
                / np.sqrt(self._depth_stdev ** 2 + self._depth_mean ** 2)
            )
            self._sigma_lognormal = np.sqrt(
                np.log((self._depth_stdev ** 2) / (self._depth_mean ** 2) + 1)
            )
            self._De = np.random.lognormal(
                self._mu_lognormal, self._sigma_lognormal, self._n
            )

        # Depth to water table - Lognormal Distribution - Variable in space
        elif self._groundwater__depth_distribution == "lognormal_spatial":

            # assert groundwater__depth_mean.shape[0] == (
            #    self._grid.number_of_nodes
            # ), "Input array should be of the length of grid.number_of_nodes!"
            # assert (groundwater__depth_standard_deviation.shape[0] == (
            #    self._grid.number_of_nodes
            # ), "Input array should be of the length of grid.number_of_nodes!"
            self._depth_mean = groundwater__depth_mean
            self._depth_stdev = groundwater__depth_standard_deviation

        # Depth to water table - Hydrologic Source Domain -> Model Domain
        elif self._groundwater__depth_distribution == "data_driven_spatial":
            self._HSD_dict = groundwater__depth_HSD_inputs

        # Check if all output fields are initialized
        self.initialize_output_fields()

        self._nodal_values = self._grid.at_node

    def _calculate_factor_of_safety(self, i):
        """Method to calculate factor of safety.

        Method calculates factor-of-safety stability index by using
        node specific parameters, creating distributions of these parameters,
        and calculating the index by sampling these distributions 'n' times.

        The index is calculated from the 'infinite slope stability
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
            self._grid.at_node["topographic__specific_contributing_area"][i]
        )
        self._theta = np.float32(self._grid.at_node["topographic__slope"][i])
        self._Tmode = np.float32(self._grid.at_node["soil__transmissivity"][i])
        self._Ksatmode = np.float32(
            self._grid.at_node["soil__saturated_hydraulic_conductivity"][i]
        )
        self._Cmode = np.int32(self._grid.at_node["soil__mode_total_cohesion"][i])
        self._Cmin = np.int32(self._grid.at_node["soil__minimum_total_cohesion"][i])
        self._Cmax = np.int32(self._grid.at_node["soil__maximum_total_cohesion"][i])
        self._phi_mode = np.float32(
            self._grid.at_node["soil__internal_friction_angle"][i]
        )
        self._rho = np.int32(self._grid.at_node["soil__density"][i])
        self._hs_mode = np.float32(self._grid.at_node["soil__thickness"][i])

        # Create a switch to imply whether Recharge or Depth to Water Table Forcing
        if self._groundwater__depth_distribution is not None:

            # Depth to water table distribution based on distribution type
            if self._groundwater__depth_distribution == "data_driven_spatial":

                self._calculate_HSD_groundwater_depth(i)

            elif self._groundwater__depth_distribution == "lognormal_spatial":
                mu_lognormal = np.log(
                    (self._depth_mean[i] ** 2)
                    / np.sqrt(self._depth_stdev[i] ** 2 + self._depth_mean[i] ** 2)
                )
                sigma_lognormal = np.sqrt(
                    np.log((self._depth_stdev[i] ** 2) / (self._depth_mean[i] ** 2) + 1)
                )
                self._De = np.random.lognormal(mu_lognormal, sigma_lognormal, self._n)

            self._hw_dist = self._hs_mode - self._De
            self._hw_dist[
                np.where(self._hw_dist < 0)
            ] = 0  # no water in soil column when De input > soil thickness

            # output mean depth of water
            self._soil__mean_watertable_depth = np.mean(self._De)

        if self._groundwater__recharge_distribution is not None:
            # recharge distribution based on distribution type
            if self._groundwater__recharge_distribution == "data_driven_spatial":
                self._calculate_HSD_recharge(i)
                self._Re /= 1000.0  # mm->m

            elif self._groundwater__recharge_distribution == "lognormal_spatial":
                mu_lognormal = np.log(
                    (self._recharge_mean[i] ** 2)
                    / np.sqrt(
                        self._recharge_stdev[i] ** 2 + self._recharge_mean[i] ** 2
                    )
                )
                sigma_lognormal = np.sqrt(
                    np.log(
                        (self._recharge_stdev[i] ** 2) / (self._recharge_mean[i] ** 2)
                        + 1
                    )
                )
                self._Re = np.random.lognormal(mu_lognormal, sigma_lognormal, self._n)
                self._Re /= 1000.0  # Convert mm to m

            # output mean recharge
            self._soil__mean_recharge = np.mean(self._Re)

        # Cohesion
        if np.all(self._grid.at_node["soil__minimum_total_cohesion"]) is not None:
            self._C = np.random.triangular(
                self._Cmin, self._Cmode, self._Cmax, size=self._n
            )
        else:
            Cmin = self._Cmode - 0.3 * self._Cmode
            Cmax = self._Cmode + 0.3 * self._Cmode
            self._C = np.random.triangular(
                self._Cmin, self._Cmode, self._Cmax, size=self._n
            )

        # phi - internal angle of friction provided in degrees
        phi_min = self._phi_mode - 0.18 * self._phi_mode
        phi_max = self._phi_mode + 0.32 * self._phi_mode
        self._phi = np.random.triangular(phi_min, self._phi_mode, phi_max, size=self._n)
        # soil thickness
        # hs_min = min(0.005, self._hs_mode-0.3*self._hs_mode) # Alternative
        hs_min = self._hs_mode - 0.3 * self._hs_mode
        hs_max = self._hs_mode + 0.1 * self._hs_mode
        self._hs = np.random.triangular(hs_min, self._hs_mode, hs_max, size=self._n)
        self._hs[self._hs <= 0.0] = 0.005

        # if Ksat provided (if T is on grid, it's not used in the calculation):
        # if np.all(self._grid.at_node['soil__saturated_hydraulic_conductivity']) is not None:
        #    # Hydraulic conductivity (Ksat)
        #    Ksatmin = self._Ksatmode - (0.3 * self._Ksatmode)
        #    Ksatmax = self._Ksatmode + (0.1 * self._Ksatmode)
        #    self._Ksat = np.random.triangular(
        #        Ksatmin, self._Ksatmode, Ksatmax, size=self._n
        #    )
        self._T = self._Ksatmode * self._hs
        # else:
        # Transmissivity (T); Ksat not provided
        Tmin = self._Tmode - (0.3 * self._Tmode)
        Tmax = self._Tmode + (0.1 * self._Tmode)
        self._T = np.random.triangular(Tmin, self._Tmode, Tmax, size=self._n)
        self._Ksat = self._T / self._hs

        # calculate Factor of Safety for n number of times
        # calculate components of FS equation
        # dimensionless cohesion
        self._C_dim = self._C / (self._hs * self._rho * self._g)

        # relative wetness
        sat_threshold = 0.01  # numerical approximation to accomodate precision of 'saturated depth to water'
        # value for saturated depth that is a not-negative not-zero value; RW = 1  at this depth.

        if self._groundwater__recharge_distribution is not None:
            self._rel_wetness = ((self._Re) / self._T) * (
                self._a / np.sin(np.arctan(self._theta))
            )

        elif self._groundwater__depth_distribution is not None:

            if self._groundwater__depth_distribution == "data_driven_spatial":

                self._rel_wetness = (self._interp_hw_dist) / (
                    self._hs_mode - sat_threshold
                )

            else:
                self._rel_wetness = (self._hs_mode - self._De) / (
                    self._hs_mode - sat_threshold
                )

        # calculate probability of saturation
        countr = 0
        for val in self._rel_wetness:  # find how many RW values >= 1
            if val >= 1.0:
                countr = countr + 1  # number with RW values (>=1)

        # probability: No. high RW values/total No. of values (n)

        self._soil__probability_of_saturation = np.float32(countr) / self._n

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
        self._landslide__probability_of_failure = np.float32(count) / self._n

    def calculate_landslide_probability(self):
        """Main method of Landslide Probability class.

        Method creates arrays for output variables then loops through
        all the core nodes to run the method
        'calculate_factor_of_safety.' Output parameters probability of
        failure, mean relative wetness, mean water table depth, and probability of saturation
        are assigned as fields to nodes.
        """
        # Create arrays for data with -9999 as default to store output
        self._mean_Relative_Wetness = np.full(self._grid.number_of_nodes, -9999.0)
        self._prob_fail = np.full(self._grid.number_of_nodes, -9999.0)
        self._prob_sat = np.full(self._grid.number_of_nodes, -9999.0)
        if self._groundwater__depth_distribution is not None:
            self._mean_watertable_depth = np.full(self._grid.number_of_nodes, -9999.0)
        else:
            self._mean_watertable_depth = np.full(self._grid.number_of_nodes, np.NaN)
        if self._groundwater__recharge_distribution is not None:
            self._mean_recharge = np.full(self._grid.number_of_nodes, -9999.0)
        else:
            self._mean_recharge = np.full(self._grid.number_of_nodes, np.NaN)

        # Run factor of safety Monte Carlo for all core nodes in domain
        # i refers to each core node id
        for i in self._grid.core_nodes:
            self._calculate_factor_of_safety(i)
            # Populate storage arrays with calculated values
            self._mean_Relative_Wetness[i] = self._soil__mean_relative_wetness
            self._prob_fail[i] = self._landslide__probability_of_failure
            self._prob_sat[i] = self._soil__probability_of_saturation
            if self._groundwater__depth_distribution is not None:
                self._mean_watertable_depth[i] = self._soil__mean_watertable_depth
            if self._groundwater__recharge_distribution is not None:
                self._mean_recharge[i] = self._soil__mean_recharge

        # Values can't be negative
        self._mean_Relative_Wetness[self._mean_Relative_Wetness < 0.0] = 0.0
        self._prob_fail[self._prob_fail < 0.0] = 0.0
        self._prob_sat[self._prob_sat < 0.0] = 0.0
        if self._groundwater__depth_distribution is not None:
            self._mean_watertable_depth[self._mean_watertable_depth < 0.0] = 0.0
        if self._groundwater__recharge_distribution is not None:
            self._mean_recharge[self._mean_recharge < 0.0] = 0.0

        # assign output fields to nodes
        self._grid.at_node["soil__mean_relative_wetness"] = self._mean_Relative_Wetness
        self._grid.at_node["landslide__probability_of_failure"] = self._prob_fail
        self._grid.at_node["soil__probability_of_saturation"] = self._prob_sat
        self._grid.at_node["soil__mean_watertable_depth"] = self._mean_watertable_depth
        self._grid.at_node["soil__mean_recharge"] = self._mean_recharge

    def _seed_generator(self, seed=0):
        """Method to initiate random seed.

        Seed the random-number generator. This method will create the
        same sequence again by re-seeding with the same value (default
        value is zero). To create a sequence other than the default,
        assign non-zero value for seed.
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
        Yrand = np.sort(np.random.rand(self._n))
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

    def _interpolate_HSD_array(self):
        """Method to extrapolate input data.

        This method uses a non-parametric approach to expand the input
        depth to water table array to the length of number of iterations. Output is
        a new array of interpolated values for each HSD id.
        """
        Yrand = np.sort(np.random.rand(self._n))
        Fx = ECDF(self._hw_dist)
        Fx_ = Fx(self._hw_dist)
        f = interpolate.interp1d(
            Fx_, self._hw_dist, bounds_error=False, fill_value=min(self._hw_dist)
        )
        # array of hw_dist  (height of water for RW calculations) interpolated from Yrand probabilities (n count)

        self._interp_hw_dist = f(Yrand)

    def _calculate_HSD_recharge(self, i):
        """Method to calculate recharge based on upstream fractions.

        This method calculates the resultant recharge at node i of the
        model domain, using recharge of contributing HSD ids and the
        areal fractions of upstream contributing HSD ids. Output is a
        numpy array of recharge at node i.
        """
        store_Re = np.zeros(self._n)
        HSD_id_list = self._HSD_id_dict[i]
        fract_list = self._fract_dict[i]
        for j in range(0, len(HSD_id_list)):
            Re_temp = self._interpolated_HSD_dict[HSD_id_list[j]]
            fract_temp = fract_list[j]
            Re_adj = Re_temp * fract_temp
            store_Re = np.vstack((store_Re, np.array(Re_adj)))
        self._Re = np.sum(store_Re, 0)

    def _calculate_HSD_groundwater_depth(self, i):
        """Method to calculate groundwater depth.

        This method calculates the resultant groundwater depth at node i of the
        model domain, using depth to water table. Output is a numpy array
        of depth to groundwater at node i.
        """
        self._hw_dist = self._hs_mode - self._HSD_dict[i]
        self._hw_dist[
            np.where(self._hw_dist < 0)
        ] = 0  # no water in soil column when De input > soil thickness

        # interpolate
        self._interpolate_HSD_array()
        self._De = self._HSD_dict[i]

    def _build_unitgrid(rowcol, edgelength):
        """ Customize UNIT TEST for bigger grids. Testing the main method 'calculate_landslide_probability()' with
        'uniform' method. Useful for testing run time and designing model experiment limits before GIS work.

        Input parameters:  dimensions of grid
        """

        grid = RasterModelGrid(rowcol, xy_spacing=edgelength)

        ls_prob = landslide_pars_ongrid(
            grid, unit_default_value
        )  # see default values given in dictionary above

        return ls_prob
