"""Landlab component that simulates fragmentation of soil grains through time.

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.soil_grading import SoilGrading

Create a grid on which to simulate soil grading.
>>> mg = RasterModelGrid((4, 5))

Create topographic_elevation and bedrock__elevation fields
>>> _ = mg.add_zeros("topographic__elevation", at="node")
>>> _ = mg.add_zeros("bedrock__elevation", at="node")

Initialise the SoilGrading component
>>> soil_grading = SoilGrading(mg)

Soil grading component created a soil layer:
>>> mg.at_node["soil__depth"][mg.core_nodes]
array([2., 2., 2., 2., 2., 2.])

Run one step of soil fragmentation
>>> soil_grading.run_one_step()
>>> mg.at_node["median_size__weight"][mg.core_nodes]
array([0.0001, 0.0001, 0.0001, 0.0001, 0.0001, 0.0001])
"""

import random
import warnings

import numpy as np

from landlab import Component


class SoilGrading(Component):
    """Simulate fragmentation of soil grains through time.

    Landlab component that simulates grading of soil particles through time
    based on mARM (Cohen et al., 2009, 2010) approach.

    The fragmentation process is controlled by a weathering transition matrix which
    defines the relative mass change in each soil grain size class (grading class)
    as a result of the fracturing of particles in the weathering mechanism.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    Cohen, S., Willgoose, G., & Hancock, G. (2009). The mARM spatially distributed soil
    evolution model: A computationally efficient modeling framework and analysis of
    hillslope soil surface organization. Journal of Geophysical Research: Earth Surface,
    114(F3).

    Cohen, S., Willgoose, G., & Hancock, G. (2010). The mARM3D spatially distributed
    soil evolution model: Three‐dimensional model framework and analysis of hillslope
    and landform responses. Journal of Geophysical Research: Earth Surface, 115(F4).
    """

    _name = "SoilGrading"
    _unit_agnostic = True
    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Topographic elevation at node",
        },
        "bedrock__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": True,
            "units": "m",
            "mapping": "node",
            "doc": "Bedrock elevation at node",
        },
        "soil__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Soil depth at node",
        },
        "grains__weight": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "kg/m^2",
            "mapping": "node",
            "doc": "Weight per unit area of grains in each size class stored at node",
        },
        "median_size__weight": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": (
                "The median grain size in each node based on soil grains"
                " weight distribution"
            ),
        },
        "grains_classes__size": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "The size of each grain size class",
        },
        "bed_grains__proportions": {
            "dtype": float,
            "intent": "out",
            "optional": True,
            "units": "[-]",
            "mapping": "node",
            "doc": "Proportional weight of each grain size class in the bed layer",
        },
    }

    def __init__(
        self,
        grid,
        meansizes=(0.0001, 0.0005, 0.002),
        limits=None,
        grains_weight=None,
        fragmentation_pattern=(0, 1),
        A_factor=1,
        initial_median_size=None,
        soil_density=2650,
        phi=0.5,
        initial_total_soil_weight=2650,
        std=None,
        CV=0.6,
        is_bedrock_distribution_flag=False,
        seed=2024,
    ):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        meansizes : float (m)
            The mean grain size in each size class.
            Can be provided as list or tuple.
        limits : float (m)
            2D array with the limits of each size class
        fragmentation_pattern : float (-)
            A list of floats describes how much mass transfer from a parent grain to
            daughters. The list must be in a size >=2 and < number of size classes,
            while the first element is percentage that remains in the parent grain.
            The other elements are the proportion of parent grain that is weathered for
            each daughter size class. The sum of fragmentation pattern list should
            be <=1. Default value set to [0, 1] which means that the entire mass of
            parent grain is transfer to the next (smaller) size class.
        A_factor : float (-. 0-1), optional
            Factor that control the fragmentation rate.
        initial_median_size : float (m), optional
            The initial median grain size in the soil.
        grains_weight : float (kg/m^2), optional
            The weight per unit area of each grain size class.
            Can be provided as list with number of elements at the size of number of classes
            or as a 2D array at the size of n_nodes x n_classes
        soil_density : float (kg/m^3), optional
            Density of the soil particles.
        phi : float (-. 0-1), optional
            Soil porosity.
        initial_total_soil_weight : float (kg/m^2), optional
            The total initial soil weight per unit area (taking into account all grain size classes).
        std: float, optional
            The standard deviation of grain size distribution.
        CV: float, optional
            Coefficient of variance of the grain size distribution.
        is_bedrock_distribution_flag: bool, optional
            A flag to indicate if the grain size distribution is generated for soil or
            bedrock layer.
        seed: float, optional
            Provide seed to set grain size distribution.
            If not provided, seed is set to 2024.
        """

        super().__init__(grid)
        self._fragmentation_pattern = fragmentation_pattern
        self._A_factor = A_factor
        self._soil_density = float(soil_density)
        self._phi = phi
        self._seed = seed
        self._CV = CV
        self._is_bedrock_distribution_flag = is_bedrock_distribution_flag
        random.seed(seed)
        
        # Get number of classes
        self._get_n_classes(meansizes=meansizes)
        
        # Create a 2D array for meansizes at node
        self._meansizes = self._create_2D_array_for_input_var(meansizes, "meansizes")

        # Set grading limits
        self.set_grading_limits(limits=limits)

        # Note: Landlabs' init_out_field procedure will not work
        # for the 'grains__weight' and 'grains_classes__size' fields
        # because the shape of these fields is: n_nodes x n_grain_sizes.
        if not grid.has_field("median_size__weight", at="node"):
            print("creating MSW")
            grid.add_zeros("median_size__weight", at="node")
        if not grid.has_field("grains_classes__size", at="node"):
            print("creating GCS")
            grid.at_node["grains_classes__size"] = np.ones(
                (grid.number_of_nodes, self._n_sizes)
            )
        if not grid.has_field("bed_grains__proportions", at="node"):
            print("creating BGP")
            grid.at_node["bed_grains__proportions"] = np.ones(
                (grid.number_of_nodes, self._n_sizes)
            )

        # Create fields for soil depth, topographic elevation and bedrock elevation
        if not grid.has_field("soil__depth"):
#MOVETHIS            warnings.warn(
#                "Soil depth is rewrite due to inconsistent with grains__weight",
#                stacklevel=2,
#            )
            print("creating SD")
            grid.add_zeros("soil__depth", at="node", clobber=True)
    
        if "topographic__elevation" not in grid.at_node:
            grid.add_zeros("topographic__elevation", at="node", clobber=True)

        if "bedrock__elevation" not in grid.at_node:
            grid.add_zeros("bedrock__elevation", at="node", clobber=True)

        # Update the weight in each size class.
        # In case grains_weight not provided, the weights will be spread around
        # the initial_median_size assuming normal distribution
        if not grid.has_field("grains__weight", at="node"):
            self._grid.at_node["grains__weight"] = np.zeros(
                (grid.number_of_nodes, self._n_sizes)
            )
            if grains_weight is None:
                if initial_median_size is None:
                    self._initial_median_size = self._meansizes[
                        self._grid.core_nodes[0], int(self._n_sizes / 2)
                    ]
                else:
                    self._get_initial_median_size(initial_median_size=initial_median_size)
                if std is None:
                    std = self._initial_median_size * self._CV
                self._std = std
                self._initial_total_soil_weight = initial_total_soil_weight
                self.generate_weight_distribution()
            else:
                self._grains_weight = self._create_2D_array_for_input_var(
                    grains_weight, "grains__weight"
                )

            # Update mass
            self._update_mass(self._grains_weight)
        else:
            print("SG was given a preexisting GW")

        # TODO: do we need to update soil, topo, or bedrock if grains__weight is provided as a field?

        # Update bed grains proportions
        self.update_bed_grains_proportions()

        # Check if the fragmentation pattern provided is valid
        self.check_fragmentation_pattern()

        # Store meansizes at grains_classes__size field and verify
        # that the number of classes match the number of classes at the grain__weight field
        if np.ndim(self._grid.at_node["grains_classes__size"])==1:
            self._grid.at_node["grains_classes__size"] *= self._meansizes[:,0]
        else:
            self._grid.at_node["grains_classes__size"] *= self._meansizes
        self._check_match_weights_n_classes()

        # Transition matrix
        self.create_transition_mat()

        # Get the median size
        self.update_median_grain_size()

    @property
    def A(self):
        """Transition matrix."""
        return self._A

    @property
    def A_factor(self):
        """Fragmentation rate"""
        return self._A_factor

    @staticmethod
    def create_grading(
        grading_name,
        grain_max_size,
        power_of=1 / 3,
        n_sizes=10,
        precent_of_volume_in_spread=10,
    ):
        """
        Parameters
        ----------
        grading_name : str
             Grading name should be provided in the from pX-AAA-BBB-CCC-DDD where X is
             the total number of daughter particles the grading fractions have, AAAA
             is the percentage of volume that remains in the parent size fraction
             after fragmentation, BBB and CCC (and so on) are the proportion of
             daughter particles in the next (smaller) size classes.
         grain_max_size : float (m)
             The maximal grain size represented in the grading distribution
         power_of : float (-)
             The power relation between parent and dauther grain volume
         n_sizes : int (-)
             The number of classes represented in the distribution
         precent_of_volume_in_spread : float (-)
             The percentage of mass from parent spread between the daughter grains
              in case of 'spread' flag is found in the grading name.
        """

        # Make sure that grading_name is valid
        grading_name_split = grading_name.split("-")
        is_p = grading_name_split[0][0] == "p"
        is_size = np.size(grading_name_split) > 2
        is_N = len(grading_name_split[0]) == 2 and grading_name_split[0][1].isnumeric()
        is_numeric = True
        sum_values = 0
        for s in grading_name.split("-")[1:]:
            if s != "spread":
                is_numeric *= s.isnumeric()
                if is_numeric:
                    sum_values += int(s)
        is_smaller_than_100 = sum_values <= 100 or "spread" in grading_name
        if is_p * is_size * is_numeric * is_smaller_than_100 * is_N is False:
            raise ValueError("grading name provided not valid")

        # Create grading distribution based on the grading name
        N = int(grading_name_split[0][1])
        parent = np.copy(grain_max_size)
        limits = np.zeros((n_sizes, 2))
        limits[0, 1] = grain_max_size
        for cnt in range(0, n_sizes):
            limits[cnt, 1] = parent
            parent = parent * (1 / N) ** (power_of)
            limits[cnt, 0] = parent
        limits = limits[::-1]  # sort in ascend order
        meansizes = np.asarray((limits[:, 0] + limits[:, 1]) / 2)

        # Create fragmentation vector describing the fragmentation pattern
        precents = np.array(
            [
                float(s)
                for s in grading_name.split("-")
                if s.replace(".", "", 1).isdigit()
            ]
        )
        if "spread" in grading_name:
            mass_precent_to_spread = 100 - precents[0] - np.sum(precents[1:])
            precents_to_add = (
                np.ones((1, n_sizes - len(precents))) * precent_of_volume_in_spread
            )
            if np.sum(precents_to_add) > mass_precent_to_spread:
                precents_to_add *= mass_precent_to_spread / np.sum(precents_to_add)
            precents = np.append(precents, precents_to_add)
        fragmentation_pattern = precents / 100

        return meansizes, limits, fragmentation_pattern

    def check_fragmentation_pattern(self):
        """
        This procedure verifies that the fragmentation pattern provided is valid
        based on the expected fragmentation format and the number of classes
        """
        if (
            len(self._fragmentation_pattern) < 2
            or len(self._fragmentation_pattern) > len(self._meansizes)
            or np.sum(self._fragmentation_pattern) > 1.0
        ):
            raise ValueError("fragmentation pattern provided not valid")

    def create_transition_mat(self):
        """
        This procedure creates a transition matrix that control the
        weathering / fragmentation pattern. The matrix represent the
        proportion of mass that weathered from each size class to other classes
        """

        self._A = np.zeros((self._n_sizes, self._n_sizes))
        self._A_factor = np.ones_like(self._A) * self._A_factor

        self._A[0, 0] = -round(
            1
            - self._fragmentation_pattern[0]
            - np.sum(self._fragmentation_pattern[1:]),
            10,
        )
        for i in range(1, self._n_sizes):

            self._A[i, i] = -(1 - (self._fragmentation_pattern[0]))
            cnti = i - 1  # rows
            cnt = 1
            while cnti >= 0 and cnt <= (len(self._fragmentation_pattern) - 1):
                if cnti == 0 and cnt <= (len(self._fragmentation_pattern) - 2):
                    self._A[cnti, i] = (1 - self._fragmentation_pattern[0]) - np.sum(
                        self._fragmentation_pattern[1:cnt]
                    )
                else:
                    self._A[cnti, i] = self._fragmentation_pattern[cnt]
                cnt += 1
                cnti -= 1

    def set_grading_limits(self, limits=None):
        """
        This procedure verifies that the array provided and describe grain size limits is valid.
        If not limits array provided, a limit array will be created based on meansizes.
        """

        if limits is None:
            lowers = (
                self._meansizes[self._grid.core_nodes[0], :-1]
                + self._meansizes[self._grid.core_nodes[0], 1:]
            ) * 0.5
            lowers = np.insert(lowers, 0, 0.0)
            uppers = np.concatenate((lowers[1:], [np.inf]))
            limits = np.empty((self._n_sizes, 2))
            limits[:, 0] = lowers
            limits[:, 1] = uppers

        elif isinstance(limits, list):
            limits = np.array(limits)

        if np.shape(limits)[0] != self._n_sizes or np.shape(limits)[1] != 2:
            raise ValueError("limits array must be in shape of n_sizes x 2")
        elif (
            np.all(limits[:, 1] > limits[:, 0])
            * np.all(np.diff(limits[:, 1]) > 0)
            * np.all(np.diff(limits[:, 0]) > 0)
        ) is False:
            raise ValueError("limits array must in ascending order")
        else:
            self._limits = limits

    def generate_weight_distribution(
        self, median_size=None, is_bedrock_distribution_flag=False
    ):
        """
        This procedure split the total soil weight between all size classes
        assuming normal distribution around the median size class. Different
        grain size distribution can be assigned to the initial soil layer and
        the bedrock layer
        """

        if not is_bedrock_distribution_flag:
            median_size = self._initial_median_size
            total_soil_weight = self._initial_total_soil_weight
            grains_weight__distribution = self._generate_normal_distribution(
                median_size=median_size, total_soil_weight=total_soil_weight
            )
            grains_weight__distribution = self._create_2D_array_for_input_var(
                grains_weight__distribution
            )
            self._update_mass(grains_weight__distribution)
        else:
            total_bedrock_weight = 10000
            # SoilGrading assumes that bedrock thickness is unlimited
            # and only the bedrock elevation is tracked.
            # The variable total_bedrock_weight is set to a large enough number just for
            # generate grain size distribution without relation to bedrock thickness.
            if median_size is None:
                grains_weight__distribution = self.g_state0
            else:
                if self._std is None:
                    self._std = self._CV * median_size
                grains_weight__distribution = self._generate_normal_distribution(
                    median_size=median_size,
                    total_soil_weight=total_bedrock_weight,
                )
            proportions = np.divide(
                grains_weight__distribution,
                np.sum(grains_weight__distribution),
                where=grains_weight__distribution > 0,
            )
            self.update_bed_grains_proportions(proportions=proportions)

    def _update_mass(self, grains_weight__distribution):
        """
        This procedure update the weight per unit area of each grain class. Based on the total weight
        at node, the soil depth field is updated. Then, topography field is also updated.
        """

        self.g_state0 = grains_weight__distribution
        if np.ndim(self._grid.at_node["grains__weight"])>1:
            self._grid.at_node["grains__weight"][self._grid.core_nodes, :] = (
                grains_weight__distribution[self._grid.core_nodes, :]
            )
            layer_depth = np.sum(
                self._grid.at_node["grains__weight"][self._grid.core_nodes], 1
            ) / (self._soil_density * (1 - self._phi))

        else:
            self._grid.at_node["grains__weight"][self._grid.core_nodes] = (
                grains_weight__distribution[self._grid.core_nodes, 0]
            )
            layer_depth = (self._grid.at_node["grains__weight"][self._grid.core_nodes] 
                           / (self._soil_density * (1 - self._phi)))

        self._grid.at_node["soil__depth"][self._grid.core_nodes] += layer_depth
        self._grid.at_node["topographic__elevation"] = (
            self._grid.at_node["soil__depth"] + self._grid.at_node["bedrock__elevation"]
        )

    def _generate_normal_distribution(
        self, median_size=None, total_soil_weight=None, std=None
    ):
        """
        This procedure spread mass between all grains classes assuming normal distribution
        centered at the median size class
        """
        if median_size is None:
            median_size = self._initial_median_size
        if std is None:
            self._std = self._CV * median_size
        if total_soil_weight is None:
            total_soil_weight = self._initial_total_soil_weight

        lower = np.min(self._limits)
        upper = np.max(self._limits)

        values = []
        if median_size < lower:
            grains_weight__distribution = np.zeros_like(self._meansizes)
            grains_weight__distribution[:, 0] = total_soil_weight
            warnings.warn(
                "Median size requested is smaller than the smallest mean size"
                "in the distribution",
                stacklevel=2,
            )

        elif median_size > upper:
            grains_weight__distribution = np.zeros_like(self._meansizes)
            grains_weight__distribution[:, -1] = total_soil_weight
            warnings.warn(
                "Median size requested is larger than the largest mean size"
                "in the distribution",
                stacklevel=2,
            )

        else:
            while np.size(values) < total_soil_weight:
                sample = random.gauss(median_size, self._std)
                if sample >= lower and sample <= upper:
                    values.append(sample)

            grains_weight__distribution = np.histogram(
                values, np.append(self._limits[:, 0], np.max(self._limits))
            )

        return grains_weight__distribution[0]

    def update_median_grain_size(self):
        """
        The median grain size at each node is defined as the size of the class closest
        to the median based on the weight in each size class
        """
        if np.ndim(self._grid.at_node["grains__weight"])>1:
            cumsum_gs = np.cumsum(self._grid.at_node["grains__weight"], axis=1)
            sum_gs = np.sum(self._grid.at_node["grains__weight"], axis=1)
            self._grid.at_node["median_size__weight"][sum_gs <= 0] = 0
            sum_gs_exp = np.expand_dims(sum_gs, -1)
    
            fraction_from_total = np.divide(
                cumsum_gs,
                sum_gs_exp,
                out=np.zeros_like(cumsum_gs),
                where=sum_gs_exp != 0,
            )
            fraction_from_total[fraction_from_total < 0.5] = np.inf
            median_val_indx = np.argmin(
                fraction_from_total - 0.5,
                axis=1,
            )
    
            self._grid.at_node["median_size__weight"][self._grid.core_nodes] = (
                self._meansizes[
                    self._grid.core_nodes, median_val_indx[self._grid.core_nodes]
                ]
            )
        else:
            self._grid.at_node["median_size__weight"][self._grid.core_nodes]= self._meansizes[self._grid.core_nodes,0]
            
    def run_one_step(self, A_factor=None):
        """
        The run_one_step procedure transform mass from parent grain size classes to
        daughters based on the fragmentation pattern.
        """
        if np.any(A_factor is None):
            A_factor = self._A_factor

        temp_g_weight = np.moveaxis(
            np.dot(
                self._A * A_factor,
                np.swapaxes(
                    np.reshape(
                        self._grid.at_node["grains__weight"][self._grid.nodes],
                        (self._grid.shape[0], self._grid.shape[1], self._n_sizes),
                    ),
                    1,
                    2,
                ),
            ),
            0,
            -1,
        )
        
        if self._n_sizes==1:
            self._grid.at_node["grains__weight"] = np.reshape(temp_g_weight, 
                                                              (self._grid.shape[0] * 
                                                               self._grid.shape[1], self._n_sizes))[:,0]
        else:
            self._grid.at_node["grains__weight"] += np.reshape(
                temp_g_weight, (self._grid.shape[0] * self._grid.shape[1], self._n_sizes)
            )
        self.update_median_grain_size()

    def _create_2D_array_for_input_var(self, input_var, var_name="None"):
        """ ""
        This procedure create a 2D array with dimensions of n_nodes x n_classes for a various
        input types.
        """

        if np.ndim(input_var) == 2:
            input_var_array = input_var
        elif isinstance(input_var, int) or isinstance(input_var, float):
            input_var_array = np.zeros(
                (np.size(self._grid.nodes.flatten()), self._n_sizes)
            )
            input_var_array[:, :] = input_var
        elif np.ndim(input_var) <= 1:
            if isinstance(input_var, list):
                input_var = np.array(input_var)
            elif isinstance(input_var, tuple):
                input_var = np.array(list(input_var))
            input_var_array = (
                np.ones((np.size(self._grid.nodes.flatten()), np.size(input_var)))
                * input_var[np.newaxis, :]
            )
        else:
            raise ValueError(f"{var_name} array format is invalid")

        return input_var_array

    def update_bed_grains_proportions(self, proportions=None):
        """ ""
        This procedure set the weight proportions of grain classes in the bed layer.
        By default, the proportion in the bed layer will set to the initial proportion of the soil layer.
        """

        if proportions is None:
            proportions = np.divide(
                self.g_state0,
                np.sum(self.g_state0, 1)[:, np.newaxis],
                where=self.g_state0 > 0,
            )
        else:
            proportions = self._create_2D_array_for_input_var(
                proportions, "bed_grains_proportions"
            )

        try:
            if np.ndim(self._grid.at_node["bed_grains__proportions"])==1:
                self._grid.at_node["bed_grains__proportions"][:] = proportions[:,0]
            else:
                self._grid.at_node["bed_grains__proportions"][:] = proportions

        except:
            raise ValueError(
                "Proportions array must be in shape of n_nodes x n_classes"
            )

    def _check_match_weights_n_classes(self):
        """
        This procedure verifies that the number of classes in grains__weight
        field and in grains_classes__size, match each other.
        """
        if np.ndim(self._grid.at_node["grains__weight"])==1:
            if np.ndim(self._grid.at_node["grains_classes__size"])!=1:
                raise ValueError(
                    "Grain weights provided do not match the number of classes"
                )

        elif (
            np.shape(self._grid.at_node["grains__weight"])[1]
            != np.shape(self._grid.at_node["grains_classes__size"])[1]
        ):
            raise ValueError(
                "Grain weights provided do not match the number of classes"
            )

    def _get_initial_median_size(self, initial_median_size):

        if isinstance(initial_median_size, int):
            self._initial_median_size = np.float(initial_median_size)
        elif isinstance(initial_median_size, float):
            self._initial_median_size = initial_median_size
        else:
            raise ValueError(
                "Initial median size must be float or integer. \n "
                "For setting initial spatial-diffrences in median grain size, \n"
                "grains_weight input parameter should be changed"
            )
        
    def _get_n_classes(self, meansizes):
        if np.ndim(meansizes) == 2:
            self._n_sizes = np.shape(input_var_array)[1]
        elif isinstance(meansizes, int) or isinstance(meansizes, float):
            self._n_sizes = 1
        elif np.ndim(meansizes) <= 1:
            self._n_sizes = np.shape(meansizes)[0]
        else:
            raise ValueError(f"meansizes format is invalid")


    def update_mass_based_on_outsource_dz(self,
                                erosion='landslide__erosion',
                                deposition='landslide__deposition',
                                proportions='bed_grains__proportions',
                                ):


        """Update the sediment mass according to information on elevation change from an external source.
        By default, this procedure is built to work with the output fields from the BedrockLandslider component
        describing elevation change from landslides.

        Parameters
        ----------
        erosion: array (float)
            Erosion (dz) at node
        deposition : array (float)
            Deposition (dz) at node
        proportions : array (float)
            Proportional weight of each grain class in the bed layer

        """

        # Make sure the inputs format is valid
        (erosion,
        deposition,
        proportions) = self._check_outsource_inputs(erosion,
                                                     deposition,
                                                     proportions)


        # Just a pointer
        grains_weight = self._grid.at_node['grains__weight']
        if np.ndim(grains_weight)==1:
            grains_weight = grains_weight[:,np.newaxis]

        # Get the total eroded/desposited mass in kg/m2
        erosion_mass = (erosion * self._soil_density *
                                        (1 - self._phi))
        deposition_mass = (deposition * self._soil_density * (1 - self._phi))

        erosion_mass = erosion_mass[:,np.newaxis]
        deposition_mass = deposition_mass[:,np.newaxis]


        # Store the mass of eroded bedrock.
        # The mass of soil will be removed later.
        bedrock_out_mass_per_class = np.sum(proportions[:,:] * erosion_mass,0)

        # Store the mass of eroded soil.
        # Here the mass will be added later.
        soil_out_mass_per_class = np.zeros_like(bedrock_out_mass_per_class)


        # Operate only over nodes with action
        non_zero_erosion_indices = np.where(erosion > 0)[0]
        non_zero_deposition_indices = np.where(deposition > 0)[0]

        if np.any(non_zero_erosion_indices):

            # Get the fraction of each grain class at node.
            a = np.sum(grains_weight[non_zero_erosion_indices, :], axis=1)[:,np.newaxis]
            b = grains_weight[non_zero_erosion_indices, :]
            grains_fractions = np.divide(b, a, where= a!= 0)

            # Partitioning the eroded soil mass across grain classes
            soil_erosion_mass = erosion_mass[non_zero_erosion_indices,:]
            soil_erosion_mass_per_class = grains_fractions * soil_erosion_mass

            # Avoid negative mass
            soil_erosion_mass_per_class = np.min((grains_weight[non_zero_erosion_indices],
                                        soil_erosion_mass_per_class),axis=0)


            # Update grains_weight field according to removed soil
            grains_weight[non_zero_erosion_indices, :] -= soil_erosion_mass_per_class

            # Update and partitioning the removed mass across despoisted grain classes
            # Remove the soil mass from the bedrock vector and add it to the soil vector
            bedrock_out_mass_per_class -= np.sum(soil_erosion_mass_per_class,0)
            soil_out_mass_per_class += np.sum(soil_erosion_mass_per_class,0)



        # Now we will collect all the removed mass and assume
        # it mixed fully before deposition.
        tot_out_mass_per_class = bedrock_out_mass_per_class+soil_out_mass_per_class

        # Get the fraction of each sediment class for deposition
        tot_deposition_mass=np.sum(tot_out_mass_per_class)
        depoistion_ratios_per_class = np.divide(tot_out_mass_per_class,
                                                tot_deposition_mass,
                                                where=tot_out_mass_per_class!=0)

        # Partitioning the deposited mass based on the dz input
        if np.any(non_zero_deposition_indices):
            grains_weight[non_zero_deposition_indices, :] +=(
                    deposition_mass[non_zero_deposition_indices] * depoistion_ratios_per_class)

        
    def _test_input_outsource_dz(self, var):

        """Verify inputs dimensions.
        """


        if isinstance(var, str):
            try:
                var = self._grid.at_node[var]
            except:
                raise ValueError(f"{var} field not exists")

        elif np.shape(var) != np.shape(self._grid.nodes.flatten()):
            raise ValueError(f"Input dimension should match the number of nodes")
        return var

    def _check_outsource_inputs(self,
                                erosion,
                                deposition,
                                proportions,
                                ):


        """Verify outsource inputs format and dimensions. 
        An error will be raised if an unexpected format is given
        """

        erosion = self._test_input_outsource_dz(erosion)
        deposition = self._test_input_outsource_dz(deposition)
        
        if isinstance(proportions, str):
            try:
                proportions = self._grid.at_node[proportions]
            except:
                raise ValueError(f"{var} field not exists")

        elif np.shape(proportions) != np.shape(self._grid.at_node['bed_grains__proportions']):
            raise ValueError(f"Proportions input dimensions should match number of nodes x number of classes")
        
        if np.ndim(proportions)==1:
            proportions = proportions[:,np.newaxis]
        return erosion, deposition, proportions
