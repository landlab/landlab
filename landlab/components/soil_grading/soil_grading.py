
"""Landlab component that simulates fragmentation of soil grains through time.

--------
>>> import numpy as np
>>> from landlab import RasterModelGrid
>>> from landlab.components.soil_grading import SoilGrading

Create a grid on which to simulate soil grading.
>>> grid = RasterModelGrid((4, 5))

Create topographic_elevation and bedrock__elevation fields
>>> mg.add_zeros("topographic__elevation", at="node")
>>> mg.add_zeros("bedrock__elevation", at="node")

Initialise the SoilGrading component
>>> soil_grading = SoilGrading(mg)

Soil grading component created a soil layer:
>>>  mg.at_node['soil__depth']
>>>  array([ 1.66666667,  1.66666667,  1.66666667,  1.66666667,  1.66666667,
        1.66666667,  1.66666667,  1.66666667,  1.66666667,  1.66666667,
        1.66666667,  1.66666667,  1.66666667,  1.66666667,  1.66666667,
        1.66666667,  1.66666667,  1.66666667,  1.66666667,  1.66666667])

Run one step of soil fragmentation
>>> soil_grading.run_one_step()
>>> mg.at_node['median_size__weight']
>>> array([ 0.0028249,  0.0028249,  0.0028249,  0.0028249,  0.0028249,
        0.0028249,  0.0028249,  0.0028249,  0.0028249,  0.0028249,
        0.0028249,  0.0028249,  0.0028249,  0.0028249,  0.0028249,
        0.0028249,  0.0028249,  0.0028249,  0.0028249,  0.0028249])
        
"""
import numpy as np
from landlab import Component
import scipy.stats as stats
import random
import time
from landlab.grid.nodestatus import NodeStatus
import warnings


class SoilGrading(Component):
    """Simulate fragmentation of soil grains through time .

        Landlab component that simulates grading of soil particles through time
        based on mARM (Cohen et al., 2009, 2010) approach.

        The fragmentation process is controlled by weathering transition matrix which defines
        the relative mass change in each soil grain size class (grading class) as a result of the fracturing
        of particles in the weathering mechanism.

        The primary method of this class is :func:`run_one_step`.

        References
        ----------
        Cohen, S., Willgoose, G., & Hancock, G. (2009). The mARM spatially distributed soil evolution model:
        A computationally efficient modeling framework and analysis of hillslope soil surface organization.
        Journal of Geophysical Research: Earth Surface, 114(F3).

        Cohen, S., Willgoose, G., & Hancock, G. (2010). The mARM3D spatially distributed soil evolution model:
        Three‐dimensional model framework and analysis of hillslope and landform responses.
        Journal of Geophysical Research: Earth Surface, 115(F4).
        """


    _name = 'SoilGrading'
    _unit_agnostic = True
    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional":True,
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
            "doc": "Weight of grains in each size fraction stored at node",
        },
        "median_size__weight": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "The median grain size in each node based on grains weight distribution",
        },
        "grains_fractions__size": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "The mean size of grain size fractions",
        },
        "bed_grains__proportions":
            {
                "dtype": float,
                "intent": "out",
                "optional": True,
                "units": "m",
                "mapping": "node",
                "doc": "Proportion of each grain size fraction in the bedrock layer",
            },
        }


    def __init__(self,
                 grid,
                 grading_name = 'p2-0-100',         # Fragmentation pattern (string)
                 n_of_grainsize_classes = 10,       # Number of size classes
                 meansizes = None,
                 alpha = 1,                         # Fragmentation rate
                 A_factor = 0.001,
                 soil_density = 2650,               # Soil density [kg/m3]
                 phi = 0.4,                         # Soil porosity [-]
                 grain_max_size = 0.02,
                 power_of = 1 / 3,
                 initial_median_size = 0.002,       # Initial median grain size [m]
                 initial_total_soil_weight = 2650, 
                 std = None,
                 CV = 0.6,
                 is_bedrock_distribution_flag = False,
                 precent_of_volume_in_spread = 10,
    ):
        """
                Parameters
                ----------
                grid : RasterModelGrid
                    A grid.
                grading_name : string
                    The name of the fragmentation pattern.
                n_of_grainsize_classes : int or float 
                    The number of grain size classes in the soil.
                meansizes : float (m)
                    The mean size of grain size in each clas.
                alpha : float (-. 0-1)
                    The fraction of parent grain size that weathered to daughters.
                A_factor : float (-. 0-1)
                    Factor that control the fragmentation rate.
                soil_density : float (kg/m^3)
                    Density of the soil particles.
                phi : float (-. 0-1)
                    Soil porosity.
                grain_max_size : float (m)
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
                    A flag to indicate if the grain size distribution is generated for soil or bedrock layer.
                precent_of_volume_in_spread: float, optional
                    The precent of volume transferred from parent to daughter in case of 'spread' grading pattern



                """
        super(SoilGrading, self).__init__(grid)

        self._grading_name = grading_name
        self._meansizes  = meansizes
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

        self._N = int(self._grading_name.split('-')[0][1:])
        self._soil_density = float(soil_density)
        self._phi = phi
        self._power_of = power_of
        self._initial_median_size = initial_median_size
        self._CV = CV
        self._initial_total_soil_weight = initial_total_soil_weight
        if std == None:
            std = self._initial_median_size * self._CV
        self._std = std
        self._is_bedrock_distribution_flag  = is_bedrock_distribution_flag
        self._precent_of_volume_in_spread = precent_of_volume_in_spread

        # Create out fields.
        # Note: Landlabs' init_out_field procedure will not work for the 'grains__weight' and 'grains_fractions__size' fields
        # because the shape of these fields is: n_nodes x n_grain_sizes.
        grid.add_field('median_size__weight',np.zeros((grid.shape[0], grid.shape[1])), at='node', dtype=float)
        grid.add_field("grains_fractions__size", np.ones((grid.shape[0], grid.shape[1], self._n_sizes )), at="node", dtype=float)
        grid.add_field("grains__weight", np.ones((grid.shape[0], grid.shape[1], self._n_sizes)), at="node", dtype=float)
        grid.add_field("bed_grains__proportions", np.ones((self._grid.shape[0], self._grid.shape[1], self._n_sizes)), at="node", dtype=float)
        
        try:
            grid.add_field("soil__depth", np.zeros((grid.shape[0], grid.shape[1])), at="node", dtype=float)
        except KeyError as exc:
            raise ValueError(f"Soil field already exists") from exc


        if "topographic__elevation" not in grid.at_node:
            grid.add_field("topographic__elevation",
                           np.zeros_like(self._grid.nodes.flatten()), at="node", dtype=float)

        if "bedrock__elevation" not in grid.at_node:
            grid.add_field("bedrock__elevation",
                           np.zeros_like(self._grid.nodes.flatten()), at="node", dtype=float)
    
        
        self.create_transition_mat()
        self.set_grading_classes()
        self.create_dist()

    def create_transition_mat(self):

        # A matrix is this is the volume/weight that weathered from each size fraction in each step
        # A matrix controls the fragmentation pattern
        # A_factor matrix is the fragmentation factor / rate

        self._A = np.zeros((self._n_sizes, self._n_sizes))
        precents = np.array([float(s) for s in self._grading_name.split('-') if s.replace('.', '', 1).isdigit()])
        if 'spread' in self._grading_name:
            precents_to_add = np.ones((1, self._n_sizes)) * self._precent_of_volume_in_spread
            precents = np.append(precents, precents_to_add)
        alphas_fractios = precents / 100
        self._A_factor = np.ones_like(self._A) * self._A_factor

        for i in range(self._n_sizes):

            if i == 0:
                self._A[i, i] = 0
            elif i == self._n_sizes:
                self._A[i, i] = -(self._alpha - (
                            self._alpha * alphas_fractios[0])
                                 )
            else:
                self._A[i, i] = -(self._alpha - (
                        self._alpha * alphas_fractios[0])
                                 )
                cnti = i - 1 # rows,
                cnt = 1
                while cnti >= 0 and cnt <= (len(alphas_fractios) - 1):
                    self._A[cnti, i] = (self._alpha * alphas_fractios[cnt])
                    if cnti == 0 and cnt <= (len(alphas_fractios) - 1):
                        self._A[cnti, i] = (1 - alphas_fractios[0]) - np.sum(alphas_fractios[1:cnt])
                        cnt += 1
                        cnti -= 1
                    cnt += 1
                    cnti -= 1

    def set_grading_classes(self, ):

        input_sizes_flag = self._input_sizes_flag
        # The grain-size classes could be a-priori set based on
        # a given geometery relations and known fragmentation pattern

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
            self._lowerlims = np.insert(np.array(meansizes),0,0)
        else:
            maxsize = self._grain_max_size
            power_of = self._power_of
            for _ in range(num_of_size_classes_plusone-1):
                upperlimits.append(maxsize)
                maxsize = lower_limit_of(maxsize)
                lowerlimits.append(maxsize)

            self._upperlims = np.sort(upperlimits)
            self._lowerlims = np.sort(lowerlimits)
            self._meansizes = ( np.array(self._upperlims) + np.array(self._lowerlims) )/ 2

        self.grid.at_node["grains_fractions__size"] *= self._meansizes
        self.update_median__size()

    def create_dist(self, median_size = None, std = None,
                    is_bedrock_distribution_flag = False):

        if not is_bedrock_distribution_flag:

            is_bedrock_distribution_flag = self._is_bedrock_distribution_flag
            median_size = self._initial_median_size
            std = self._std
            total_soil_weight = self._initial_total_soil_weight
            grains_weight__distribution = self._generate_normal_distribution(median_size = median_size,
                                         std=std,
                                         total_soil_weight=total_soil_weight
                                         )

        else:

            total_bedrock_weight = 10000 # There is no meaning for this number
            if median_size == None:
                grains_weight__distribution = self.g_state0
            else:
                median_size = median_size 
                if std == None:
                    std = self._CV * median_size
                else:
                    std =std
                grains_weight__distribution = self._generate_normal_distribution(median_size = median_size,
                                                                           std=std,
                                                                           total_soil_weight = total_bedrock_weight
                                                                           )
                    

        if not is_bedrock_distribution_flag:
            self.g_state = np.ones(
                (int(np.shape(self.grid)[0]), int(np.shape(self.grid)[1]), int(len(grains_weight__distribution)))) * grains_weight__distribution

            self.g_state0 = grains_weight__distribution
            self.grid.at_node['grains__weight'] *= grains_weight__distribution
            layer_depth = np.sum(self.g_state0) / (self._soil_density * self.grid.dx * self.grid.dx )  # meter
            layer_depth /= (1-self._phi) # Porosity correction

            self.grid.at_node['soil__depth'] +=  layer_depth
            self.grid.at_node['topographic__elevation'] = self.grid.at_node['soil__depth'] + self.grid.at_node['bedrock__elevation']

            self.g_state_bedrock = grains_weight__distribution
            self.grid.at_node["bed_grains__proportions"][:] = 1
            self.grid.at_node["bed_grains__proportions"] *= np.divide(self.g_state_bedrock, np.sum(self.g_state_bedrock))

        else:
            self.g_state_bedrock = grains_weight__distribution
            self.grid.at_node["bed_grains__proportions"][:] = 1
            self.grid.at_node["bed_grains__proportions"] *= np.divide(self.g_state_bedrock, np.sum(self.g_state_bedrock))


    def _generate_normal_distribution(self, median_size = None,
                                     std=None,
                                     total_soil_weight=None
                                     ):

        if median_size == None:
            median_size = self._initial_median_size
        if std == None:
            std = self._CV * median_size
        if total_soil_weight == None:
            total_soil_weight = self._initial_total_soil_weight

        lower = np.min(self._lowerlims)
        upper = np.max(self._upperlims)

        values = []
        if median_size < lower:
            grains_weight__distribution = np.zeros_like(self._upperlims)
            grains_weight__distribution[0] = total_soil_weight
            warnings.warn("Median size provided is smaller than the distribution lower bound")

        elif median_size > upper:
            grains_weight__distribution = np.zeros_like(self._upperlims)
            grains_weight__distribution[-1] = total_soil_weight
            warnings.warn("Median size provided is larger than the distribution upper bound")

        else:
            while np.size(values) < total_soil_weight:
                sample = random.gauss(median_size, std)
                if sample >= lower and sample <= upper:
                    values.append(sample)

            grains_weight__distribution = np.histogram(values, np.insert(self._upperlims, 0, 0))[0]

        return grains_weight__distribution

    def update_median__size(self):

        # Median size based on weight distribution
        cumsum_gs = np.cumsum(self.grid.at_node['grains__weight'], axis=1)
        sum_gs = np.sum(self.grid.at_node['grains__weight'], axis=1)
        self.grid.at_node['median_size__weight'][sum_gs <= 0] = 0
        sum_gs_exp = np.expand_dims(sum_gs, -1)

        median_val_indx = np.argmax(np.where(
                np.divide(
                    cumsum_gs,
                    sum_gs_exp,
                    out=np.zeros_like(cumsum_gs),
                    where=sum_gs_exp != 0)>=0.5, 1, 0)
            ,axis=1
        )

        self.grid.at_node['median_size__weight'] = self._meansizes[median_val_indx[:]]

    def run_one_step(self, A_factor=None):
        # Run one step of fragmentation
        # 'update_median__size' procedure is called after

        if np.any(A_factor == None):
            A_factor = self._A_factor

        temp_g_weight = np.moveaxis(
            np.dot(
                self._A * A_factor, np.swapaxes(
                    np.reshape(self.grid.at_node['grains__weight'], (self.grid.shape[0],
                                                        self.grid.shape[1],
                                                        self._n_sizes)),
                    1,2)),
            0, -1)
        
        
        self.grid.at_node['grains__weight'] += np.reshape(temp_g_weight,
                                                         (self.grid.shape[0] * self.grid.shape[1],
                                                          self._n_sizes))
        self.update_median__size()