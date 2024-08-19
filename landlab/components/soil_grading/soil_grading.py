"""Landlab component that simulates fragmentation of soil grains through time.

Examples
--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.soil_grading import SoilGrading

Create a grid on which to simulate soil grading.
>>> mg = RasterModelGrid((4, 5))

Create topographic_elevation and bedrock__elevation fields
>>> mg.add_zeros("topographic__elevation", at="node")
>>> mg.add_zeros("bedrock__elevation", at="node")

Initialise the SoilGrading component
>>> soil_grading = SoilGrading(mg)

Soil grading component created a soil layer:
>>> mg.at_node["soil__depth"][mg.core_nodes]
array([2., 2., 2., 2., 2., 2.])
       
Run one step of soil fragmentation
>>> soil_grading.run_one_step()
>>> mg.at_node["median_size__weight"][mg.core_nodes]
array([0.0028249, 0.0028249, 0.0028249, 0.0028249, 0.0028249, 0.0028249])

"""

import random
import warnings

import numpy as np

from landlab import Component


class SoilGrading(Component):
    """Simulate fragmentation of soil grains through time .

    Landlab component that simulates grading of soil particles through time
    based on mARM (Cohen et al., 2009, 2010) approach.

    The fragmentation process is controlled by weathering transition matrix which defines
    the relative mass change in each soil grain size class (grading class) as a result of
    the fracturing of particles in the weathering mechanism.

    The primary method of this class is :func:`run_one_step`.

    References
    ----------
    Cohen, S., Willgoose, G., & Hancock, G. (2009). The mARM spatially distributed soil
    evolution model: A computationally efficient modeling framework and analysis of hillslope
    soil surface organization. Journal of Geophysical Research: Earth Surface, 114(F3).

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
            "units": "kg",
            "mapping": "node",
            "doc": "Weight of grains in each size class stored at node",
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
            "units": "m",
            "mapping": "node",
            "doc": "Proportional weight of each grain size class in bedrock",
        },
    }

    def __init__(
        self,
        grid,
        grading_name="p2-0-100",
        n_of_grainsize_classes=10,
        meansizes=None,
        alpha=1,
        A_factor=0.001,
        soil_density=2650,
        phi=0.5,
        grain_max_size=0.02,
        power_of=1 / 3,
        initial_median_size=0.002,
        initial_total_soil_weight=2650,
        std=None,
        CV=0.6,
        is_bedrock_distribution_flag=False,
        precent_of_volume_in_spread=10,
        seed=2024,
    ):
        """
        Parameters
        ----------
        grid : RasterModelGrid
            A grid.
        grading_name : string, optional
            The name of the fragmentation pattern.
        n_of_grainsize_classes : int or float, optional
            The number of grain size classes in the soil.
        meansizes : float (m), optional
            The mean size of grain size in each clas.
        alpha : float (-. 0-1), optional
            The fraction of parent grain size that weathered to daughters.
        A_factor : float (-. 0-1), optional
            Factor that control the fragmentation rate.
        soil_density : float (kg/m^3), optional
            Density of the soil particles.
        phi : float (-. 0-1), optional
            Soil porosity.
        grain_max_size : float (m), optional
            The maximal grain size in the distribution.
        power_of : float, optional
            A parameter that controls the geometry relation between the grain sizes in the soil.
        initial_median_size : float (m), optional
            The initial median grain size in the soil.
        initial_total_soil_weight : float (kg), optional
            The total initial soil weight (taking into account all grain size classes).
        std: float, optional
            The standard deviation of grain size distribution.
        CV: float, optional
            Coefficient of variance of the grain size distribution.
        is_bedrock_distribution_flag: bool, optional
            A flag to indicate if the grain size distribution is generated for soil or
            bedrock layer.
        precent_of_volume_in_spread: float, optional
            The precent of volume transferred from parent to daughter in case of 'spread'
             grading pattern
        seed: float, optional
            Provide seed to set grain size distribution.
            If not provided, seed is set to 2024.
        """

        super().__init__(grid)

        self._grading_name = grading_name
        self._meansizes = meansizes
        self._alpha = alpha
        self._A_factor = A_factor

        if self._meansizes is not None:
            self._n_sizes = np.size(self._meansizes)
            self._grain_max_size = np.max(self._meansizes)
            self._input_sizes_flag = True

        else:
            self._input_sizes_flag = False
            self._n_sizes = n_of_grainsize_classes
            self._grain_max_size = grain_max_size

        self._N = int(self._grading_name.split("-")[0][1:])
        self._soil_density = float(soil_density)
        self._phi = phi
        self._power_of = power_of
        self._initial_median_size = initial_median_size
        self._CV = CV
        self._initial_total_soil_weight = initial_total_soil_weight
        if std is None:
            std = self._initial_median_size * self._CV
        self._std = std
        self._is_bedrock_distribution_flag = is_bedrock_distribution_flag
        self._precent_of_volume_in_spread = precent_of_volume_in_spread
        self._seed = seed
        random.seed(seed)

        # Note: Landlabs' init_out_field procedure will not work
        # for the 'grains__weight' and 'grains_classes__size' fields
        # because the shape of these fields is: n_nodes x n_grain_sizes.
        grid.add_field(
            "median_size__weight",
            np.zeros(grid.shape[0] * grid.shape[1]),
            at="node",
            dtype=float,
        )
        grid.add_field(
            "grains_classes__size",
            np.ones((int(grid.shape[0] * grid.shape[1]), self._n_sizes)),
            at="node",
            dtype=float,
        )
        grid.add_field(
            "grains__weight",
            np.zeros((int(grid.shape[0] * grid.shape[1]), self._n_sizes)),
            at="node",
            dtype=float,
        )
        grid.add_field(
            "bed_grains__proportions",
            np.ones((int(grid.shape[0] * grid.shape[1]), self._n_sizes)),
            at="node",
            dtype=float,
        )

        if grid.has_field("soil__depth"):
            warnings.warn(
                "Soil depth is rewrite due to inconsistent with grains__weight",
                stacklevel=2,
            )
        grid.add_zeros("soil__depth", at="node", clobber=True)

        if "topographic__elevation" not in grid.at_node:
            grid.add_field(
                "topographic__elevation",
                np.zeros_like(self._grid.nodes.flatten(), dtype=float),
                at="node",
                dtype=float,
            )

        if "bedrock__elevation" not in grid.at_node:
            grid.add_field(
                "bedrock__elevation",
                np.zeros_like(self._grid.nodes.flatten(), dtype=float),
                at="node",
                dtype=float,
            )

        self.create_transition_mat()
        self.set_grading_classes()
        self.create_dist()
        self.update_median__size()

    @property
    def A(self):
        """Transition matrix."""
        return self._A

    @property
    def A_factor(self):
        """Fragmentation rate"""
        return self._A_factor

    @property
    def grading_name(self):
        """The name of fragmentation pattern"""
        return self._grading_name

    def create_transition_mat(self):
        """
        This procedure creates a transition matrix that control the
        weathering / fragmentation pattern. The matrix represent the
        proportion of mass that weathered from each size class to other classes
        """

        self._A = np.zeros((self._n_sizes, self._n_sizes))
        precents = np.array(
            [
                float(s)
                for s in self._grading_name.split("-")
                if s.replace(".", "", 1).isdigit()
            ]
        )
        if "spread" in self._grading_name:
            precents_to_add = (
                np.ones((1, self._n_sizes)) * self._precent_of_volume_in_spread
            )
            precents = np.append(precents, precents_to_add)
        alphas_fractios = precents / 100
        self._A_factor = np.ones_like(self._A) * self._A_factor

        for i in range(self._n_sizes):

            if i == 0:
                self._A[i, i] = 0
            elif i == self._n_sizes:
                self._A[i, i] = -(self._alpha - (self._alpha * alphas_fractios[0]))
            else:
                self._A[i, i] = -(self._alpha - (self._alpha * alphas_fractios[0]))
                cnti = i - 1  # rows,
                cnt = 1
                while cnti >= 0 and cnt <= (len(alphas_fractios) - 1):
                    self._A[cnti, i] = self._alpha * alphas_fractios[cnt]
                    if cnti == 0 and cnt <= (len(alphas_fractios) - 1):
                        self._A[cnti, i] = (1 - alphas_fractios[0]) - np.sum(
                            alphas_fractios[1:cnt]
                        )
                        cnt += 1
                        cnti -= 1
                    cnt += 1
                    cnti -= 1

    def set_grading_classes(
        self,
    ):

        input_sizes_flag = self._input_sizes_flag

        def lower_limit_of(maxsize):
            lower_limit = maxsize * (1 / self._N) ** (power_of)
            return lower_limit

        upperlimits = []
        lowerlimits = []
        num_of_size_classes_plusone = self._n_sizes + 1

        if input_sizes_flag:
            meansizes = self._meansizes
            self._meansizes = np.array(meansizes)
            self._upperlims = np.array(meansizes)
            self._lowerlims = np.insert(np.array(meansizes), 0, 0)
        else:
            maxsize = self._grain_max_size
            power_of = self._power_of
            for _ in range(num_of_size_classes_plusone - 1):
                upperlimits.append(maxsize)
                maxsize = lower_limit_of(maxsize)
                lowerlimits.append(maxsize)

            self._upperlims = np.sort(upperlimits)
            self._lowerlims = np.sort(lowerlimits)
            self._meansizes = (
                np.array(self._upperlims) + np.array(self._lowerlims)
            ) / 2

        self._grid.at_node["grains_classes__size"][self._grid.nodes] *= self._meansizes

    def create_dist(
        self, median_size=None, std=None, is_bedrock_distribution_flag=False
    ):
        """
        This procedure split the total soil weight between all size classes
        assuming normal distribution around the median size class. Different
        grain size distribution can be assigned to the initial soil layer and
        the bedrock layer
        """

        if not is_bedrock_distribution_flag:

            is_bedrock_distribution_flag = self._is_bedrock_distribution_flag
            median_size = self._initial_median_size
            std = self._std
            total_soil_weight = self._initial_total_soil_weight
            grains_weight__distribution = self._generate_normal_distribution(
                median_size=median_size, std=std, total_soil_weight=total_soil_weight
            )

        else:
            total_bedrock_weight = 10000  # SoilGrading assumes that bedrock thickness is unlimited
                                          # and only the bedrock elevation is tracked.
                                          # The variable total_bedrock_weight is set to a large enough number just for
                                          # generate grain size distribution without relation to bedrock thickness.
            if median_size is None:
                grains_weight__distribution = self.g_state0
            else:
                median_size = median_size
                if std is None:
                    std = self._CV * median_size
                else:
                    std = std
                grains_weight__distribution = self._generate_normal_distribution(
                    median_size=median_size,
                    std=std,
                    total_soil_weight=total_bedrock_weight,
                )

        if not is_bedrock_distribution_flag:
            self.g_state = (
                np.ones(
                    (
                        int(np.shape(self._grid)[0]),
                        int(np.shape(self._grid)[1]),
                        int(len(grains_weight__distribution)),
                    )
                )
                * grains_weight__distribution
            )

            self.g_state0 = grains_weight__distribution
            self._grid.at_node["grains__weight"][self._grid.core_nodes, :] = 1
            self._grid.at_node["grains__weight"] *= grains_weight__distribution
            layer_depth = np.sum(self.g_state0) / (
                self._soil_density * self._grid.dx * self._grid.dx
            )
            layer_depth /= 1 - self._phi

            self._grid.at_node["soil__depth"][self._grid.core_nodes] += layer_depth
            self._grid.at_node["topographic__elevation"] = (
                self._grid.at_node["soil__depth"]
                + self._grid.at_node["bedrock__elevation"]
            )

            self.g_state_bedrock = grains_weight__distribution
            self._grid.at_node["bed_grains__proportions"][self._grid.core_nodes] = 1
            self._grid.at_node["bed_grains__proportions"][
                self._grid.core_nodes
            ] *= np.divide(self.g_state_bedrock, np.sum(self.g_state_bedrock))

        else:
            self.g_state_bedrock = grains_weight__distribution
            self._grid.at_node["bed_grains__proportions"][self._grid.core_nodes] = 1
            self._grid.at_node["bed_grains__proportions"][
                self._grid.core_nodes
            ] *= np.divide(self.g_state_bedrock, np.sum(self.g_state_bedrock))

    def _generate_normal_distribution(
        self, median_size=None, std=None, total_soil_weight=None
    ):

        if median_size is None:
            median_size = self._initial_median_size
        if std is None:
            std = self._CV * median_size
        if total_soil_weight is None:
            total_soil_weight = self._initial_total_soil_weight

        lower = np.min(self._lowerlims)
        upper = np.max(self._upperlims)

        values = []
        if median_size < lower:
            grains_weight__distribution = np.zeros_like(self._upperlims)
            grains_weight__distribution[0] = total_soil_weight
            warnings.warn(
                "Median size provided is smaller than the distribution lower bound",
                stacklevel=2,
            )

        elif median_size > upper:
            grains_weight__distribution = np.zeros_like(self._upperlims)
            grains_weight__distribution[-1] = total_soil_weight
            warnings.warn(
                "Median size provided is larger than the distribution upper bound",
                stacklevel=2,
            )

        else:
            while np.size(values) < total_soil_weight:
                sample = random.gauss(median_size, std)
                if sample >= lower and sample <= upper:
                    values.append(sample)

            grains_weight__distribution = np.histogram(
                values, np.insert(self._upperlims, 0, 0)
            )[0]

        return grains_weight__distribution

    def update_median__size(self):
        """
        The median grain size at each node is defined as the size of the class closest
        to the median based on the weight in each size class
        """

        cumsum_gs = np.cumsum(self._grid.at_node["grains__weight"], axis=1)
        sum_gs = np.sum(self._grid.at_node["grains__weight"], axis=1)
        self._grid.at_node["median_size__weight"][sum_gs <= 0] = 0
        sum_gs_exp = np.expand_dims(sum_gs, -1)

        median_val_indx = np.argmin(
            np.abs(
                np.divide(
                    cumsum_gs,
                    sum_gs_exp,
                    out=np.zeros_like(cumsum_gs),
                    where=sum_gs_exp != 0,
                )
                - 0.5
            ),
            axis=1,
        )

        self._grid.at_node["median_size__weight"][self._grid.core_nodes] = (
            self._meansizes[median_val_indx[self._grid.core_nodes]]
        )

    def run_one_step(self, A_factor=None):

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

        self._grid.at_node["grains__weight"] += np.reshape(
            temp_g_weight, (self._grid.shape[0] * self._grid.shape[1], self._n_sizes)
        )
        self.update_median__size()
